#variational Poisson-lognormal model (Chiquet J., Mariadassou M. & Robin S., Ann. Appl. Stat. (2018): https://doi.org/10.1214/18-AOAS1177)
#with a design X, latent factors orthogonal to X and an optional offset factor (loadings fixed to ones)
#    yᵢⱼ ~ Poisson(exp(oᵢ + Zᵢⱼ)),   Zᵢⱼ = xᵢ'Bⱼ + wᵢ'Γⱼ + εᵢⱼ,   εᵢⱼ ~ N(0, σⱼ²),   X'W = 0
#the variational distribution is q(Zᵢⱼ) = N(Mᵢⱼ, Sᵢⱼ²); in the M-step the rotation by the QR decomposition of X
#makes B the exact fit of the X rows of Q'M and (W, Γ) a weighted PCA of the residual rows

#W₀(exp(t)), the Lambert W of an exponential, by Newton iterations on w + log(w) = t (no overflow for large t)
function lambertw_exp(t::Real)
    w = t > 1 ? t - log(t) : exp(t)
    for _ in 1:50
        Δ = (w + log(w) - t) / (1 + 1 / w)
        w = max(w - Δ, w / 10)
        abs(Δ) <= 1.0e-12 * w && break
    end
    return w
end

#variational E-step for one cell: maximise y M - exp(o + M + S² / 2) - ((M - μ)² + S²) / (2σ²) + log(S²) / 2
#M has a closed form through the Lambert W function, S² is a contractive fixed point
function estep_cell(y, o, μ, σ2, M, S2; inner = 3)
    for _ in 1:inner
        M = μ + σ2 * y - lambertw_exp(log(σ2) + o + μ + σ2 * y + S2 / 2)
        S2 = 1 / (exp(o + M + S2 / 2) + 1 / σ2)
    end
    return M, S2
end

#evidence lower bound (without the constant -Σ log(yᵢⱼ!))
function pln_elbo(Y, o, M, S2, μ, σ2)
    elbo = 0.0
    for j in axes(Y, 2), i in axes(Y, 1)
        elbo += Y[i, j] * (o[i] + M[i, j]) - exp(o[i] + M[i, j] + S2[i, j] / 2) -
            0.5 * (log(σ2[j]) + ((M[i, j] - μ[i, j])^2 + S2[i, j]) / σ2[j]) + 0.5 * (log(S2[i, j]) + 1)
    end
    return elbo
end

#M-step: fitted means μ = X B + W Γ (W orthogonal to X) and the residual variances σ²
function pln_mstep(M, S2, QR, k, q, nfix, σ2, shrink_var)
    n, p = size(M)
    Mrot = QR.Q' * M
    Mres = Mrot[(k + 1):n, :]
    #weighted PCA of the residual rows (columns weighted by 1 / σ²), offset factor first
    w = 1 ./ σ2
    l0 = nfix > 0 ? (Mres * w) ./ sum(w) : zeros(n - k)
    R = Mres .- l0
    D = svd(R .* sqrt.(w)')
    Lres = D.U[:, 1:q] .* D.S[1:q]'
    Γ = D.Vt[1:q, :] ./ sqrt.(w)'
    fitted = [Mrot[1:k, :]; l0 .+ Lres * Γ]
    μ = QR.Q * fitted
    #residual variances with n - k - q - nfix degrees of freedom, shrunk across genes
    d = n - k - q - nfix
    σ2 = max.(vec(sum(abs2, M - μ; dims = 1) + sum(S2; dims = 1)) ./ d, 1.0e-8)
    shrunk = shrink_variances(σ2, d, shrink_var)
    B = UpperTriangular(QR.R) \ Mrot[1:k, :]
    return μ, shrunk.s2, B, nfix > 0 ? [ones(1, p); Γ] : Γ, shrunk.df
end

struct PLNFit
    M::Matrix{Float64} #variational means of the latent log-abundances (offsets removed)
    S2::Matrix{Float64} #variational variances
    μ::Matrix{Float64} #fitted means X B + W Γ
    B::Matrix{Float64} #coefficients of X (k x p), confounded by W
    Γ::Matrix{Float64} #loadings (first row = ones with the offset factor)
    σ2::Vector{Float64} #lognormal overdispersion of each gene
    offsets::Vector{Float64}
    df::Vector{Float64} #degrees of freedom of the overdispersions after shrinkage
    elbo::Float64
    iterations::Int
end

"""
    pln(Y, X, q; offsets = log library sizes, offset_factor = true, maxiter = 500, tol = 1e-7, shrink_var = :vash)

fit the PLN model to counts Y (samples x genes) with design X and q latent factors
offsets default to the log library sizes
"""
function pln(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        offsets::AbstractVector{<:Real} = vec(log.(sum(Y; dims = 2))), offset_factor::Bool = true,
        maxiter::Integer = 500, tol::Real = 1.0e-7, shrink_var::Symbol = :vash
    )
    n, p = size(Y)
    k = size(X, 2)
    nfix = Int(offset_factor)
    @assert n - k > q + nfix "not enough residual degrees of freedom"
    o = collect(float.(offsets))
    QR = qr(X)
    #initialisation from the log counts
    M = log.(Y .+ 0.5) .- o
    S2 = 1 ./ (Y .+ 1.0)
    σ2 = vec(var(M; dims = 1))
    μ, σ2, B, Γ, df = pln_mstep(M, S2, QR, k, q, nfix, σ2, shrink_var)
    elbo = -Inf
    iter = 0
    for outer in 1:maxiter
        iter = outer
        Threads.@threads for j in 1:p
            for i in 1:n
                M[i, j], S2[i, j] = estep_cell(Y[i, j], o[i], μ[i, j], σ2[j], M[i, j], S2[i, j])
            end
        end
        μ, σ2, B, Γ, df = pln_mstep(M, S2, QR, k, q, nfix, σ2, shrink_var)
        elbonew = pln_elbo(Y, o, M, S2, μ, σ2)
        converged = abs(elbonew - elbo) < tol * abs(elbonew)
        elbo = elbonew
        converged && break
    end
    return PLNFit(M, S2, μ, B, Γ, σ2, o, df, elbo, iter)
end

"""
    pln_mouthwash(Y, X, q; controls = Int[], cov_of_interest = size(X, 2), offset_factor = true, nullweight = 10, estimate_λ = true, kwargs...)

empirical Bayes branch: PLN fit (stage 1) followed by MOUTHWASH (stage 2)
β̃ⱼ = q'Mⱼ / R₂₂ with q the column of Q for the covariate of interest; its standard error uses the latent-scale
variance of each cell vᵢⱼ = σⱼ² + 1 / λᵢⱼ (lognormal overdispersion + Poisson noise, λᵢⱼ the expected count),
the variances of M itself are too small in the residual rows that are shrunk toward the fitted factors
"""
function pln_mouthwash(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        controls::AbstractVector{<:Integer} = Int[], cov_of_interest::Integer = size(X, 2),
        offset_factor::Bool = true, nullweight::Real = 10, estimate_λ::Bool = true, alpha_uncertainty::Bool = true,
        likelihood::Symbol = :t, kwargs...
    )
    fit = pln(Y, X, q; offset_factor = offset_factor, kwargs...)
    k = size(X, 2)
    F = qr(X[:, [setdiff(1:k, cov_of_interest); cov_of_interest]])
    qint = (F.Q * Matrix(I, size(X, 1), k))[:, k]
    R22 = F.R[k, k]
    β̃ = vec(qint' * fit.M) ./ R22
    λcount = exp.(fit.offsets .+ fit.μ .+ fit.σ2' ./ 2)
    s = sqrt.(vec(sum(qint .^ 2 .* (fit.σ2' .+ 1 ./ λcount); dims = 1))) ./ abs(R22)
    mw = mouthwash_fit(
        β̃, s, fit.Γ, collect(controls); nullweight = nullweight, estimate_λ = estimate_λ,
        alpha_uncertainty = alpha_uncertainty, df = likelihood === :t ? fit.df : Inf
    )
    return (
        pln = fit,
        mouthwash = MouthwashFit(
            β̃, s, mw.α, mw.λ, mw.g, mw.loglik, mw.mean, mw.sd, mw.lfsr, mw.lfdr,
            cumulative_rate(mw.lfsr), cumulative_rate(mw.lfdr), mw.evalue, collect(controls),
        ),
    )
end
