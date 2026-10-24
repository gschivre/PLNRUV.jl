#the Poisson-lognormal layer shared by plnpca and pln_mouthwash (Chiquet J., Mariadassou M. & Robin S., Ann. Appl.
#Stat. (2018): https://doi.org/10.1214/18-AOAS1177): the variational step of one count, the fit result, and the
#empirical Bayes branch (PLN fit followed by MOUTHWASH)
#    yᵢⱼ ~ Poisson(exp(oᵢ + Zᵢⱼ)),   Zᵢⱼ = xᵢ'Bⱼ + wᵢ'Γⱼ + εᵢⱼ,   εᵢⱼ ~ N(0, σⱼ²),   X'W = 0,   q(Zᵢⱼ) = N(Mᵢⱼ, Sᵢⱼ²)

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

#variational E-step for one count: maximise y M - exp(o + M + S² / 2) - ((M - μ)² + S²) / (2σ²) + log(S²) / 2
#M has a closed form through the Lambert W function, S² is a contractive fixed point
function estep_count(y, o, μ, σ2, M, S2; inner = 3)
    for _ in 1:inner
        M = μ + σ2 * y - lambertw_exp(log(σ2) + o + μ + σ2 * y + S2 / 2)
        S2 = 1 / (exp(o + M + S2 / 2) + 1 / σ2)
    end
    return M, S2
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
    σ2mix::NamedTuple #posterior mixture of each overdispersion (weights w, squared scales s2, df; genes x components)
    Γcov::Vector{Matrix{Float64}} #posterior covariance of each gene's column of Γ
    elbo::Float64
    iterations::Int
    converged::Bool #change of the ELBO per count below tol before maxiter
end

"""
    pln_mouthwash(Y, X, q; controls = Int[], cov_of_interest = size(X, 2), offset_factor = true, nullweight = 10, estimate_λ = true, likelihood = :t, mode = :auto, kwargs...)

empirical Bayes branch: PLN fit (stage 1) followed by MOUTHWASH (stage 2)
bⱼ = u'Mⱼ / Rₖₖ with u the column of Q for the covariate of interest; its standard error uses the latent-scale
variance of each count σⱼ² + 1 / ȳᵢⱼ (lognormal variance + Poisson noise, ȳᵢⱼ the expected count),
the variances of M itself are too small in the residual rows that are shrunk toward the fitted factors
with controls, mode = :auto centres the non-null part of the effect prior at an estimated mode when a shift shared by
many genes (a global change of transcription) raises the held-out likelihood (see ash); the null stays at zero
the PLN layer is fitted by plnpca (kwargs go to it); the posterior covariance of its loadings is passed to MOUTHWASH,
which adds the error of the confounding correction, α'Γcovⱼα, to the noise of each gene
"""
function pln_mouthwash(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        controls::AbstractVector{<:Integer} = Int[], cov_of_interest::Integer = size(X, 2),
        offset_factor::Bool = true, nullweight::Real = 10, estimate_λ::Bool = true, alpha_uncertainty::Bool = true,
        likelihood::Symbol = :t, mode::Symbol = :auto, mode_threshold::Real = 0.05, kwargs...
    )
    fit = plnpca(Y, X, q; offset_factor = offset_factor, kwargs...)
    k = size(X, 2)
    F = qr(X[:, [setdiff(1:k, cov_of_interest); cov_of_interest]])
    qint = (F.Q * Matrix(I, size(X, 1), k))[:, k]
    R22 = F.R[k, k]
    β̃ = vec(qint' * fit.M) ./ R22
    λcount = exp.(fit.offsets .+ fit.μ .+ fit.σ2' ./ 2)
    s = sqrt.(vec(sum(qint .^ 2 .* (fit.σ2' .+ 1 ./ λcount); dims = 1))) ./ abs(R22)
    #per variance component m of the overdispersion: Σᵢ qᵢ² (σ̃²ⱼₘ + 1 / λᵢⱼ) / R₂₂² (Σᵢ qᵢ² = 1)
    poisson = vec(sum(qint .^ 2 ./ λcount; dims = 1))
    mixv = (fit.σ2mix.s2 .+ poisson) ./ R22^2
    noise = noise_model(likelihood, s .^ 2, fit.df, fit.σ2mix.w, mixv, fit.σ2mix.df)
    mw = mouthwash_fit(
        β̃, s, fit.Γ, collect(controls); nullweight = nullweight, estimate_λ = estimate_λ,
        alpha_uncertainty = alpha_uncertainty, noise = noise, Γcov = fit.Γcov,
        mode = mode, mode_threshold = mode_threshold
    )
    return (
        pln = fit,
        mouthwash = MouthwashFit(
            β̃, s, mw.α, mw.λ, mw.g, mw.loglik, mw.mean, mw.sd, mw.lfsr, mw.lfdr,
            cumulative_rate(mw.lfsr), cumulative_rate(mw.lfdr), mw.evalue, collect(controls), mw.converged,
        ),
    )
end
