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
#by exact coordinate steps, so that each one increases it: M has a closed form through the Lambert W function, S² is
#the root of its own equation (see s2_given_m); a single fixed-point step S² <- 1 / (exp(o + M + S² / 2) + 1 / σ²)
#instead oscillates with a growing amplitude for a zero count with a large σ², which can drive a very sparse gene to
#an ever larger variance
function estep_count(y, o, μ, σ2, M, S2; inner = 3)
    for _ in 1:inner
        M = μ + σ2 * y - lambertw_exp(log(σ2) + o + μ + σ2 * y + S2 / 2)
        S2 = s2_given_m(o + M, σ2, S2)
    end
    return M, S2
end

#S² maximising -exp(c + S² / 2) - S² / (2σ²) + log(S²) / 2 (concave): the root of 1 / S² = exp(c + S² / 2) + 1 / σ²,
#written φ(S²) = c + S² / 2 + log S² - log(1 - S² / σ²) = 0 (increasing, without the exponential, whose convexity
#makes Newton steps crawl); the root lies between 1 / (exp(c + σ² / 2) + 1 / σ²) and 1 / (exp(c) + 1 / σ²): Newton
#steps from S20 kept inside this bracket, bisection (geometric while the bracket is wide) when a step leaves it
function s2_given_m(c, σ2, S20)
    lo = 1 / (exp(c + σ2 / 2) + 1 / σ2)
    hi = 1 / (exp(c) + 1 / σ2)
    S = clamp(S20, lo, hi)
    for _ in 1:100
        φ = c + S / 2 + log(S) - log1p(-S / σ2)
        φ < 0 ? (lo = S) : (hi = S)
        Snew = S - φ / (1 / 2 + 1 / S + 1 / (σ2 - S))
        lo < Snew < hi || (Snew = hi > 4lo ? sqrt(lo * hi) : (lo + hi) / 2)
        done = abs(Snew - S) <= 1.0e-12 * S
        S = Snew
        done && break
    end
    return S
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
    pln_mouthwash(Y, X, q; controls = Int[], contrast = last column of X, offset_factor = true, nullweight = 10, estimate_λ = true, likelihood = :t, mode = :auto, kwargs...)

empirical Bayes branch: PLN fit (stage 1) followed by MOUTHWASH (stage 2), for counts Y (samples x genes), design X
and q latent factors
contrast: the effect tested, c'B, a vector with one weight per column of X (default: the last column; or
cov_of_interest = the index of a column); a matrix with one contrast per column gives one result per contrast from a
single fit (see test_contrast to add more later, and contrast_vector to build contrasts from column names)
genes whose counts do not determine every coefficient (a group, a level or a time point without any count) are left
out of the fit and reported as not tested (mouthwash.separated, estimates NaN)
the raw estimate of gene j is a'Mⱼ with a = X (X'X)⁻¹ c; its standard error uses the latent-scale variance of each
count σⱼ² + 1 / ȳᵢⱼ (lognormal variance + Poisson noise, ȳᵢⱼ the expected count); with controls, mode = :auto
centres the non-null part of the effect prior at an estimated mode when a shift shared by many genes (a global
change of transcription) raises the held-out likelihood (see ash); the null stays at zero
the PLN layer is fitted by plnpca (kwargs go to it); the posterior covariance of its loadings is passed to MOUTHWASH,
which adds the error of the confounding correction, α'Γcovⱼα, to the noise of each gene
returns (pln = the PLN fit of the tested genes, mouthwash = a MouthwashFit, or a vector of them, over all genes,
tested = the indices of the tested genes, ngenes)
"""
function pln_mouthwash(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        controls::AbstractVector{<:Integer} = Int[], cov_of_interest::Integer = size(X, 2), contrast = nothing,
        offset_factor::Bool = true, nullweight::Real = 10, estimate_λ::Bool = true, alpha_uncertainty::Bool = true,
        likelihood::Symbol = :t, mode::Symbol = :auto, mode_threshold::Real = 0.05, kwargs...
    )
    #genes whose observed samples do not determine every coefficient (a group, a level or a time point without any
    #count): their effects are set by nothing in the data, so they are left out and reported as not tested
    sep = separated_genes(Y, X)
    keep = findall(!, sep)
    fit = plnpca(Y[:, keep], X, q; offset_factor = offset_factor, kwargs...)
    res = (pln = fit, mouthwash = nothing, tested = keep, ngenes = size(Y, 2))
    if contrast === nothing
        contrast = zeros(size(X, 2))
        contrast[cov_of_interest] = 1
    end
    opts = (; controls, nullweight, estimate_λ, alpha_uncertainty, likelihood, mode, mode_threshold)
    mw = contrast isa AbstractVector ? test_contrast(res, X, contrast; opts...) :
        [test_contrast(res, X, c; opts...) for c in eachcol(contrast)]
    return (pln = fit, mouthwash = mw, tested = keep, ngenes = size(Y, 2))
end

"""
    test_contrast(res, X, c; controls = Int[], nullweight = 10, estimate_λ = true, likelihood = :t, mode = :auto)

MOUTHWASH for the contrast c'B (a vector with one weight per column of X) on the PLN fit of a result of pln_mouthwash
(res.pln, on the genes res.tested): the fit does not depend on the contrast, so further contrasts cost only the
MOUTHWASH step; the raw estimate of gene j is a'Mⱼ with a = X (X'X)⁻¹ c, and the confounders (with the compositional
offset) are estimated anew for each contrast; returns a MouthwashFit over all genes (the separated genes untested)
"""
function test_contrast(
        res, X::AbstractMatrix{<:Real}, c::AbstractVector{<:Real}; controls::AbstractVector{<:Integer} = Int[], kwargs...
    )
    length(c) == size(X, 2) || throw(DimensionMismatch("the contrast has $(length(c)) weights for $(size(X, 2)) columns of X"))
    keep = res.tested
    mw = contrast_test(res.pln, X, c; controls = findall(in(Set(controls)), keep), kwargs...)
    return expand_genes(mw, keep, res.ngenes, collect(controls))
end

"""
    contrast_vector(names, spec)

the contrast vector of a design whose columns are named names: spec is a column name (its coefficient), or pairs
name => weight, e.g. ["time: 6h" => 1, "time: 2h" => -1] for the difference between two time points
"""
contrast_vector(names::AbstractVector{<:AbstractString}, spec::AbstractString) = contrast_vector(names, [spec => 1])
function contrast_vector(names::AbstractVector{<:AbstractString}, spec::AbstractVector{<:Pair})
    c = zeros(length(names))
    for (name, w) in spec
        i = findfirst(==(string(name)), names)
        i === nothing && throw(ArgumentError("no column named \"$name\" (the columns are: $(join(names, ", ")))"))
        c[i] += w
    end
    return c
end

"""
    design_matrix(formula, data)

the design matrix X and the names of its columns for a formula such as @formula(0 ~ genotype * time) and a table
of samples (a DataFrame, a NamedTuple of vectors or any Tables.jl table), with the intercept unless the formula
removes it (0 ~ 0 + ...); available after `using StatsModels`, as are pln_mouthwash(Y, formula, data, q; contrast)
and choose_q(Y, formula, data), where contrast may also be given by column names (see contrast_vector)
"""
function design_matrix end

#genes whose samples with a positive count leave some coefficient of the design undetermined (X restricted to them
#has rank below k); the analogue of complete separation in logistic regression: no finite estimate exists
function separated_genes(Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real})
    k = size(X, 2)
    return [rank(X[view(Y, :, j) .> 0, :]) < k for j in axes(Y, 2)]
end

#MOUTHWASH for the estimate a'Mⱼ of a linear combination of the coefficients, a = X (X'X)⁻¹ c: for the column of
#interest (c = eₖ) a is the column of Q over Rₖₖ; its standard error uses the latent-scale variance of each count,
#σⱼ² + 1 / ȳᵢⱼ (the variances of M itself are too small in the residual rows, pulled towards the fitted factors)
function contrast_test(
        fit::PLNFit, X::AbstractMatrix{<:Real}, c::AbstractVector{<:Real}; controls = Int[], nullweight = 10,
        estimate_λ = true, alpha_uncertainty = true, likelihood = :t, mode = :auto, mode_threshold = 0.05
    )
    X = float(X) #an integer X'X would be solved in BigFloat, and every later step with it
    a = X * (Symmetric(X' * X) \ c)
    β̃ = vec(a' * fit.M)
    λcount = exp.(fit.offsets .+ fit.μ .+ fit.σ2' ./ 2)
    poisson = vec(sum(a .^ 2 ./ λcount; dims = 1))
    s = sqrt.(sum(abs2, a) .* fit.σ2 .+ poisson)
    #per variance component m of the overdispersion: Σᵢ aᵢ² σ̃²ⱼₘ + Σᵢ aᵢ² / λᵢⱼ
    mixv = sum(abs2, a) .* fit.σ2mix.s2 .+ poisson
    noise = noise_model(likelihood, s .^ 2, fit.df, fit.σ2mix.w, mixv, fit.σ2mix.df)
    mw = mouthwash_fit(
        β̃, s, fit.Γ, collect(controls); nullweight = nullweight, estimate_λ = estimate_λ,
        alpha_uncertainty = alpha_uncertainty, noise = noise, Γcov = fit.Γcov, mode = mode, mode_threshold = mode_threshold
    )
    return MouthwashFit(
        β̃, s, mw.α, mw.λ, mw.g, mw.loglik, mw.mean, mw.sd, mw.lfsr, mw.lfdr, cumulative_rate(mw.lfsr),
        cumulative_rate(mw.lfdr), mw.evalue, collect(controls), Int[], mw.converged,
    )
end

#results of the tested genes (keep) placed among all p genes: the separated ones get no estimate (NaN) and an e-value
#of 0; controls are given in the numbering of all genes
function expand_genes(mw::MouthwashFit, keep, p, controls)
    full(v, x0) = (out = fill(x0, p); out[keep] .= v; out)
    sep = setdiff(1:p, keep)
    return MouthwashFit(
        full(mw.β̃, NaN), full(mw.s, NaN), mw.α, mw.λ, mw.g, mw.loglik, full(mw.mean, NaN), full(mw.sd, NaN),
        full(mw.lfsr, NaN), full(mw.lfdr, NaN), full(mw.svalue, NaN), full(mw.qvalue, NaN), full(mw.evalue, 0.0),
        setdiff(controls, sep), sep, mw.converged,
    )
end
