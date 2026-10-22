#MOUTHWASH of Gerard D. & Stephens M., Biostatistics (2020): https://doi.org/10.1093/biostatistics/kxy029
#after the rotation, the rows of interest give β̃ⱼ = Y₂ⱼ / R₂₂ ~ N(βⱼ + Γⱼ'α, λ sⱼ²) with Γ the factor loadings
#estimated on the residual rows; α (confounders, including the compositional offset if offset_factor = true),
#the ash prior g on β and the variance inflation λ are estimated jointly by maximising the marginal likelihood
#with control genes (βⱼ = 0) α is their generalised least squares fit (RUV-4), which stays valid when most genes
#change (a global shift of transcription); without controls α is identified by the sparsity of β through g

"""
    factor_analysis(Y, q; offset_factor = true, iter = 5)

heteroscedastic factor analysis Y ≈ L F, residual variances σ² (columns), by alternating weighted PCA
with offset_factor = true the first row of F is fixed to ones (per-sample shift common to all genes)
"""
function factor_analysis(Y::AbstractMatrix{<:Real}, q::Integer; offset_factor::Bool = true, iter::Integer = 5)
    m, p = size(Y)
    nfix = Int(offset_factor)
    @assert m > q + nfix "not enough residual degrees of freedom"
    σ2 = vec(var(Y; dims = 1))
    F = zeros(q + nfix, p)
    for _ in 1:iter
        W = 1 ./ σ2
        l0 = offset_factor ? (Y * W) ./ sum(W) : zeros(m) #weighted mean shift of each row
        R = Y .- l0
        D = svd(R .* sqrt.(W)')
        F1 = D.Vt[1:q, :] ./ sqrt.(W)'
        res = R - (D.U[:, 1:q] .* D.S[1:q]') * D.Vt[1:q, :] ./ sqrt.(W)'
        σ2 = max.(vec(sum(abs2, res; dims = 1)) ./ (m - q - nfix), 1.0e-8)
        F = offset_factor ? [ones(1, p); F1] : F1
    end
    return (F = F, σ2 = σ2)
end

#maximise Q(λ) = Σⱼₖ wⱼₖ log N(rⱼ; 0, λ sⱼ² + σₖ²) over log λ ∈ [-5, 5] by golden-section search
function update_λ(r, s2, σ2, W; lo = -5.0, hi = 5.0, tol = 1.0e-6)
    Q(logλ) = sum(W[j, k] == 0 ? 0.0 : W[j, k] * lognormpdf(r[j], exp(logλ) * s2[j] + σ2[k]) for j in eachindex(r), k in eachindex(σ2))
    g = (sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    Qc, Qd = Q(c), Q(d)
    while b - a > tol
        if Qc > Qd
            b, d, Qd = d, c, Qc
            c = b - g * (b - a)
            Qc = Q(c)
        else
            a, c, Qc = c, d, Qd
            d = a + g * (b - a)
            Qd = Q(d)
        end
    end
    return exp((a + b) / 2)
end

#joint estimation of g, α and λ by block coordinate ascent:
#π by mix-SQP, α by weighted least squares (MM step; fixed by the controls when there are any), λ by a
#one-dimensional search
function mouthwash_fit(β̃, s, Γ, controls; nullweight = 10, estimate_λ = true, tol = 1.0e-7, maxiter = 500)
    p = length(β̃)
    s2 = s .^ 2
    free = setdiff(1:p, controls)
    #α from the controls only (weights 1 / sⱼ², λ cancels), or a starting value from all genes without controls
    if isempty(controls)
        α = Γ' \ β̃
    else
        ΓC = Γ[:, controls]
        α = Symmetric(ΓC * (ΓC' ./ s2[controls])) \ (ΓC * (β̃[controls] ./ s2[controls]))
    end
    λ = 1.0
    r = β̃ - Γ' * α
    grid = ash_grid(r[free], s[free])
    σ2 = grid .^ 2
    K = length(σ2)
    prior = [nullweight; ones(K - 1)]
    Lik = zeros(p, K)
    logscale = zeros(p)
    W = zeros(p, K)
    π = fill(1 / K, K)
    obj = -Inf
    for _ in 1:maxiter
        mixture_lik!(Lik, logscale, r, s2, σ2, λ)
        π, _ = mixture_weights(Lik, logscale, free, prior; π0 = π)
        #observed penalised log-likelihood (controls only use the point mass)
        objnew = mixture_objective(π, Lik, logscale, free, prior) + sum(j -> log(Lik[j, 1]) + logscale[j], controls; init = 0.0)
        converged = objnew - obj < tol * (abs(objnew) + tol)
        obj = objnew
        converged && break
        #posterior component probabilities
        for j in free
            w = view(Lik, j, :) .* π
            W[j, :] .= w ./ sum(w)
        end
        for j in controls
            W[j, :] .= 0.0
            W[j, 1] = 1.0
        end
        #α (sparsity only): weighted least squares with weights cⱼ = Σₖ wⱼₖ / (λ sⱼ² + σₖ²)
        if isempty(controls)
            c = [sum(W[j, k] / (λ * s2[j] + σ2[k]) for k in 1:K) for j in 1:p]
            α = Symmetric(Γ * (c .* Γ')) \ (Γ * (c .* β̃))
            r = β̃ - Γ' * α
        end
        estimate_λ && (λ = update_λ(r, s2, σ2, W))
    end
    post = ash_posterior(r, λ .* s2, σ2, π)
    for j in controls #controls are fixed to zero
        post.mean[j] = post.sd[j] = 0.0
        post.lfsr[j] = post.lfdr[j] = 1.0
    end
    e = evalues(r, λ .* s2, grid, π)
    e[controls] .= 0.0
    return merge((α = α, λ = λ, g = (σ = grid, π = π), loglik = obj, evalue = e), post)
end

struct MouthwashFit
    β̃::Vector{Float64} #unadjusted estimates (rotated rows of interest / R₂₂)
    s::Vector{Float64} #their standard errors (before variance inflation)
    α::Vector{Float64} #confounder coefficients (first = compositional offset if offset_factor = true)
    λ::Float64 #variance inflation
    g::NamedTuple #ash prior (σ grid and weights π)
    loglik::Float64
    mean::Vector{Float64} #posterior mean of β
    sd::Vector{Float64}
    lfsr::Vector{Float64}
    lfdr::Vector{Float64}
    svalue::Vector{Float64}
    qvalue::Vector{Float64}
    evalue::Vector{Float64} #compound e-values (use with ebh for FDR control under dependence)
    controls::Vector{Int}
end

"""
    mouthwash(Y, X, q; controls = Int[], cov_of_interest = size(X, 2), offset_factor = true, nullweight = 10, estimate_λ = true, shrink_var = :vash)

MOUTHWASH on Gaussian data Y (samples x genes) with design X and q latent factors
controls are optional negative controls; one covariate of interest (default the last column of X)
"""
function mouthwash(
        Y::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}, q::Integer;
        controls::AbstractVector{<:Integer} = Int[], cov_of_interest::Integer = size(X, 2),
        offset_factor::Bool = true, nullweight::Real = 10, estimate_λ::Bool = true, shrink_var::Symbol = :vash
    )
    rot = rotate(Y, X, [cov_of_interest])
    fa = factor_analysis(rot.Yres, q; offset_factor = offset_factor)
    R22 = rot.R22[1, 1]
    β̃ = vec(rot.Yint) ./ R22
    d = size(rot.Yres, 1) - q - Int(offset_factor) #residual degrees of freedom of the variances
    s = sqrt.(shrink_variances(fa.σ2, d, shrink_var)) ./ abs(R22)
    fit = mouthwash_fit(β̃, s, fa.F, collect(controls); nullweight = nullweight, estimate_λ = estimate_λ)
    return MouthwashFit(
        β̃, s, fit.α, fit.λ, fit.g, fit.loglik, fit.mean, fit.sd, fit.lfsr, fit.lfdr,
        cumulative_rate(fit.lfsr), cumulative_rate(fit.lfdr), fit.evalue, collect(controls),
    )
end
