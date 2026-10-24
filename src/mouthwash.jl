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

#maximise Q(λ) = Σⱼₖ wⱼₖ log N(rⱼ; mₖ, λ sⱼ² + σₖ²) (or the t / mixture of t densities of the noise model N) over
#log λ ∈ [-5, 5] by golden-section search
function update_λ(r, N::NoiseModel, σ2, W; means = zeros(length(σ2)), lo = -5.0, hi = 5.0, tol = 1.0e-6)
    Q(logλ) = sum(
        W[j, k] == 0 ? 0.0 : W[j, k] * lognoise(N, j, r[j] - means[k], σ2[k], exp(logλ))
            for j in eachindex(r), k in eachindex(σ2)
    )
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
#with alpha_uncertainty = true the variance of each corrected estimate rⱼ = β̃ⱼ - Γⱼ'α̂ includes the uncertainty of α̂
#(as in the CATE variance of Wang et al.): λ (sⱼ² + Γⱼ' V Γⱼ), V = (Γ_C diag(1 / s_C²) Γ_C')⁻¹ with controls; without
#controls, V is the inverse weighted Gram matrix of the α step at convergence, and the prior is refitted
#Γcov (one covariance matrix per gene, or empty) is the uncertainty of the estimated loadings: the correction Γⱼ'α
#then has the error variance lⱼ(α) = α'Γcovⱼα, which is added to the noise of gene j (inside λ, as the uncertainty
#of α̂); it vanishes without confounding (α ≈ 0) and, as Γcovⱼ scales with the latent variance σⱼ² and not with the
#Poisson noise, weighs most on highly expressed genes, which a single multiplicative λ cannot cover
#the noise model of β̃ (default: a t with df degrees of freedom and squared scale sⱼ²) can be a mixture of t
#(vash); s is then the effective standard error used for the α weights and the grid
#with controls the non-null part of the prior may be centred at an estimated mode instead of zero (see ash, mode)
function mouthwash_fit(
        β̃, s, Γ, controls; nullweight = 10, estimate_λ = true, alpha_uncertainty = true, df = Inf,
        noise::NoiseModel = NoiseModel(s .^ 2, df), mode::Symbol = :auto, mode_threshold::Real = 0.05,
        Γcov::AbstractVector = Matrix{Float64}[], tol = 1.0e-7, maxiter = 500
    )
    p = length(β̃)
    s2base = s .^ 2
    loading_var(α) = isempty(Γcov) ? zeros(p) : [dot(α, Γcov[j] * α) for j in 1:p]
    #normal equations of α with uncertain loadings: Σⱼ cⱼ E[ΓⱼΓⱼ'] = Σⱼ cⱼ (ΓⱼΓⱼ' + Γcovⱼ) (errors in variables); along
    #a direction that the loadings do not determine (a factor nearly switched off) Γcov dominates and α stays finite,
    #where Γ diag(c) Γ' alone would be singular
    gram(genes, c) = Symmetric(Γ[:, genes] * (Γ[:, genes]' .* c) .+ (isempty(Γcov) ? 0.0 : sum(c[i] .* Γcov[j] for (i, j) in enumerate(genes))))
    free = setdiff(1:p, controls)
    #α from the controls only (weights 1 / (sⱼ² + lⱼ(α)), λ cancels; fixed-point iteration until α settles, as lⱼ
    #depends on α), or a starting value from all genes without controls
    if isempty(controls)
        α = gram(1:p, 1 ./ s2base) \ (Γ * (β̃ ./ s2base))
    else
        ΓC = Γ[:, controls]
        α = zeros(size(Γ, 1))
        G = Symmetric(ΓC * ΓC')
        for _ in 1:(isempty(Γcov) ? 1 : 50)
            wC = 1 ./ (s2base[controls] .+ loading_var(α)[controls])
            G = gram(controls, wC)
            αnew = G \ (ΓC * (β̃[controls] .* wC))
            done = maximum(abs, αnew .- α; init = 0.0) < 1.0e-10 * (1 + maximum(abs, αnew; init = 0.0))
            α = αnew
            done && break
        end
    end
    l = loading_var(α)
    s2 = s2base .+ l
    N = add_noise(noise, l)
    if !isempty(controls)
        if alpha_uncertainty #Var(rⱼ) = λ (sⱼ² + lⱼ + Γⱼ' G⁻¹ Γⱼ) for the non-control genes
            a = zeros(p)
            a[free] .= vec(sum(Γ[:, free] .* (G \ Γ[:, free]); dims = 1))
            s2 = s2 .+ a
            N = add_noise(N, a)
        end
    end
    λ = 1.0
    r = β̃ - Γ' * α
    #grid from the standard errors that include the uncertainty of α̂ (the variances the likelihood uses)
    grid = ash_grid(r[free], sqrt.(s2[free]))
    #block coordinate ascent for the prior components (σ2, means)
    function fitloop(σ2, means, λ, α, r)
        K = length(σ2)
        prior = [nullweight; ones(K - 1)]
        Lik = zeros(p, K)
        logscale = zeros(p)
        W = zeros(p, K)
        π = fill(1 / K, K)
        obj = -Inf
        converged = false
        for _ in 1:maxiter
            mixture_lik!(Lik, logscale, r, N, σ2, λ; means = means)
            π, _ = mixture_weights(Lik, logscale, free, prior; π0 = π)
            #observed penalised log-likelihood (controls only use the point mass)
            objnew = mixture_objective(π, Lik, logscale, free, prior) + sum(j -> log(Lik[j, 1]) + logscale[j], controls; init = 0.0)
            #absolute change: with a t likelihood the α and λ steps are not exact ascent steps, and a decrease of
            #the objective must not pass for convergence
            converged = abs(objnew - obj) < tol * (abs(objnew) + tol)
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
            #α (sparsity only, components centred at zero): weighted least squares with weights
            #cⱼ = Σₖ wⱼₖ / (λ sⱼ² + σₖ²)
            if isempty(controls)
                c = [sum(W[j, k] / (λ * s2[j] + σ2[k]) for k in 1:K) for j in 1:p]
                α = gram(1:p, c) \ (Γ * (c .* β̃))
                r = β̃ - Γ' * α
                if !isempty(Γcov) #the loading term follows α (l, s2 and N of the enclosing function: the steps after
                    #the loop use their last values)
                    l = loading_var(α)
                    s2 = s2base .+ l
                    N = add_noise(noise, l)
                end
            end
            estimate_λ && (λ = update_λ(r, N, σ2, W; means = means))
        end
        return (; π, obj, W, λ, α, r, Lik, logscale, prior, converged)
    end
    σ2, means = mode_components(grid, 0.0)
    f = fitloop(σ2, means, λ, α, r)
    #mode of the non-null effects (with controls only: without them a shift shared by the genes cannot be told
    #apart from the compositional offset), chosen on the non-control genes at the fitted inflation
    μ, gain = isempty(controls) ? (0.0, NaN) : choose_mode(mode, f.r, N, f.λ, free, grid, nullweight, mode_threshold)
    if μ != 0
        σ2, means = mode_components(grid, μ)
        f = fitloop(σ2, means, f.λ, f.α, f.r)
    end
    (; π, obj, W, λ, α, r, Lik, logscale, prior, converged) = f
    K = length(σ2)
    Nv = scale_noise(N, λ)
    if alpha_uncertainty && isempty(controls)
        #Var(α̂) from the weighted least squares of the last α step, then refit the prior with the corrected variances
        c = [sum(W[j, k] / (λ * s2[j] + σ2[k]) for k in 1:K) for j in 1:p]
        Vα = inv(gram(1:p, c))
        Nv = add_noise(Nv, vec(sum(Γ .* (Vα * Γ); dims = 1)))
        mixture_lik!(Lik, logscale, r, Nv, σ2, 1.0; means = means)
        π, obj = mixture_weights(Lik, logscale, free, prior; π0 = π)
    end
    post = ash_posterior(r, Nv, σ2, π; means = means)
    for j in controls #controls are fixed to zero
        post.mean[j] = post.sd[j] = 0.0
        post.lfsr[j] = post.lfdr[j] = 1.0
    end
    e = evalues(r, Nv, sqrt.(σ2), π; means = means)
    e[controls] .= 0.0
    return merge((α = α, λ = λ, g = (σ = sqrt.(σ2), μ = means, π = π, mode = μ, mode_gain = gain), loglik = obj, evalue = e, converged = converged), post)
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
    converged::Bool #relative change of the objective below tol before maxiter
end

#noise model of the estimates for a likelihood option: :tmix the mixture of t of the vash posterior (exact for the
#moderated estimates), :t a single t with the averaged degrees of freedom, :normal the normal
function noise_model(likelihood::Symbol, s2, df, mixw, mixv, mixdf)
    likelihood === :tmix && return NoiseModel(mixw, mixv, mixdf)
    likelihood === :t && return NoiseModel(s2, df)
    likelihood === :normal && return NoiseModel(s2, Inf)
    throw(ArgumentError("unknown likelihood $likelihood (use :tmix, :t or :normal)"))
end

"""
    mouthwash(Y, X, q; controls = Int[], cov_of_interest = size(X, 2), offset_factor = true, nullweight = 10, estimate_λ = true, shrink_var = :vash, var_trend = false, likelihood = :t, mode = :auto)

MOUTHWASH on Gaussian data Y (samples x genes) with design X and q latent factors
controls are optional negative controls; one covariate of interest (default the last column of X)
with var_trend = true the variances are shrunk around their trend in the mean of each column of Y
"""
function mouthwash(
        Y::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}, q::Integer;
        controls::AbstractVector{<:Integer} = Int[], cov_of_interest::Integer = size(X, 2),
        offset_factor::Bool = true, nullweight::Real = 10, estimate_λ::Bool = true, shrink_var::Symbol = :vash,
        alpha_uncertainty::Bool = true, likelihood::Symbol = :t, var_trend::Bool = false, mode::Symbol = :auto,
        mode_threshold::Real = 0.05
    )
    rot = rotate(Y, X, [cov_of_interest])
    fa = factor_analysis(rot.Yres, q; offset_factor = offset_factor)
    R22 = rot.R22[1, 1]
    β̃ = vec(rot.Yint) ./ R22
    d = size(rot.Yres, 1) - q - Int(offset_factor) #residual degrees of freedom of the variances
    shrunk = shrink_variances(fa.σ2, d, shrink_var; covariate = var_trend ? vec(mean(Y; dims = 1)) : nothing)
    s = sqrt.(shrunk.s2) ./ abs(R22)
    noise = noise_model(likelihood, s .^ 2, shrunk.df, shrunk.mix.w, shrunk.mix.s2 ./ R22^2, shrunk.mix.df)
    fit = mouthwash_fit(
        β̃, s, fa.F, collect(controls); nullweight = nullweight, estimate_λ = estimate_λ,
        alpha_uncertainty = alpha_uncertainty, noise = noise,
        mode = mode, mode_threshold = mode_threshold
    )
    return MouthwashFit(
        β̃, s, fit.α, fit.λ, fit.g, fit.loglik, fit.mean, fit.sd, fit.lfsr, fit.lfdr,
        cumulative_rate(fit.lfsr), cumulative_rate(fit.lfdr), fit.evalue, collect(controls), fit.converged,
    )
end
