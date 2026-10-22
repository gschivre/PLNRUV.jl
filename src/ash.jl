#adaptive shrinkage of Stephens M., Biostatistics (2017): https://doi.org/10.1093/biostatistics/kxw041
#normal means problem rⱼ ~ N(βⱼ, sⱼ²) with the prior g = Σₖ πₖ N(0, σₖ²), σ₁ = 0 being the point mass at zero
#the weights π are estimated by maximum a posteriori with a Dirichlet(nullweight, 1, ..., 1) prior (mix-SQP)

#log density of N(0, v) at x
lognormpdf(x::Real, v::Real) = -0.5 * (log(2π) + log(v) + x^2 / v)

#normalising constant of the Student t density with ν degrees of freedom (ν = Inf: the normal)
tconst(ν::Real) = isinf(ν) ? -0.5 * log(2π) : loggamma((ν + 1) / 2) - loggamma(ν / 2) - 0.5 * log(ν * π)

#log density at x of a t with ν degrees of freedom and squared scale v (normal when ν = Inf); c = tconst(ν)
logtpdf(x::Real, v::Real, ν::Real, c::Real) = isinf(ν) ? lognormpdf(x, v) : c - 0.5 * log(v) - (ν + 1) / 2 * log1p(x^2 / (ν * v))

#degrees of freedom of each estimate: a scalar (Inf for the normal likelihood) or one value per gene
gene_df(ν::Real, p) = fill(float(ν), p)
gene_df(ν::AbstractVector, p) = collect(float.(ν))

#noise model of the estimates: rⱼ - βⱼ ~ Σₘ wⱼₘ t(νⱼₘ) with squared scale vⱼₘ (normal when νⱼₘ = Inf)
#one column is a single t or normal; several columns are the exact likelihood when the variances were shrunk by a
#mixture prior (vash): given the variance component m, the moderated estimate follows a t with 2aₘ + d degrees of
#freedom (Lemma 1 of Lu M. & Stephens M., arXiv:1901.10679 (2019)), and wⱼₘ is the posterior weight of m
struct NoiseModel
    w::Matrix{Float64}
    v::Matrix{Float64}
    ν::Matrix{Float64}
    c::Matrix{Float64} #normalising constants tconst(ν)
    function NoiseModel(w::AbstractMatrix, v::AbstractMatrix, ν::AbstractMatrix)
        @assert size(w) == size(v) == size(ν)
        return new(float.(w), float.(v), float.(ν), tconst.(ν))
    end
end
NoiseModel(v::AbstractVector, ν = Inf) =
    NoiseModel(ones(length(v), 1), reshape(collect(float.(v)), :, 1), reshape(gene_df(ν, length(v)), :, 1))
Base.length(N::NoiseModel) = size(N.v, 1)

#degrees of freedom of a t(ν) with squared scale v convolved with an independent N(0, (λ - 1) v), λ > 1:
#the variance inflation of MOUTHWASH is additive noise that the variance estimate does not see (loadings,
#sample-specific variances, ...), not a stretch of the t; matching the excess kurtosis, 6 / (ν - 4) scaled by the
#squared share 1 / λ² of the t part, gives a t with 4 + (ν - 4) λ² degrees of freedom and squared scale λ v
#(λ ≤ 1 deflates the estimated variances, the degrees of freedom are kept)
inflated_df(ν::Real, λ::Real) = λ > 1 && isfinite(ν) && ν > 4 ? 4 + (ν - 4) * λ^2 : float(ν)

#the noise model with every squared scale inflated by λ (with the degrees of freedom of the additive inflation),
#or increased by a vector a of normal variances (degrees of freedom kept: a is small next to the scales)
scale_noise(N::NoiseModel, λ::Real) = NoiseModel(N.w, λ .* N.v, inflated_df.(N.ν, λ))
add_noise(N::NoiseModel, a::AbstractVector) = NoiseModel(N.w, N.v .+ a, N.ν)

#log density of gene j at x when the prior component of variance σ2 is convolved with the noise (the t
#approximation of the convolution), with the noise inflated by λ (see inflated_df)
function lognoise(N::NoiseModel, j::Integer, x::Real, σ2::Real, λ::Real = 1.0)
    function term(k)
        N.w[j, k] == 0 && return -Inf
        ν = inflated_df(N.ν[j, k], λ)
        c = ν == N.ν[j, k] ? N.c[j, k] : tconst(ν)
        return log(N.w[j, k]) + logtpdf(x, λ * N.v[j, k] + σ2, ν, c)
    end
    M = size(N.w, 2)
    M == 1 && return term(1) - log(N.w[j, 1])
    m = maximum(term, 1:M)
    return m + log(sum(k -> exp(term(k) - m), 1:M))
end

#grid of prior standard deviations: from min(s) / 10 to 2 * sqrt(max(r² - s²)), multiplicative step mult
function ash_grid(r::AbstractVector{<:Real}, s::AbstractVector{<:Real}; mult::Real = sqrt(2))
    σmin = minimum(s) / 10
    σmax = 2 * sqrt(max(maximum(r .^ 2 .- s .^ 2), 0.0))
    σmax = max(σmax, 8 * σmin)
    K = ceil(Int, log(σmax / σmin) / log(mult)) + 1
    return [0.0; σmin .* mult .^ (0:(K - 1))]
end

#component likelihoods Lⱼₖ = N(rⱼ; mₖ, λ sⱼ² + σₖ²) scaled by their row maximum (returned in logscale), mₖ the
#component means (0 for ash, see mode_components); with ν (degrees of freedom of each sⱼ) the normal is replaced by
#a t with scale √(λ sⱼ² + σₖ²), the usual approximation of the convolution of a t likelihood with a normal mixture
#(exact for the point masses); with a NoiseModel N, Lⱼₖ = Σₘ wⱼₘ t(rⱼ - mₖ; νⱼₘ, λ vⱼₘ + σₖ²)
mixture_lik!(Lik, logscale, r, s2, σ2, λ; ν = Inf, means = zeros(length(σ2))) =
    mixture_lik!(Lik, logscale, r, NoiseModel(s2, ν), σ2, λ; means = means)
function mixture_lik!(Lik, logscale, r, N::NoiseModel, σ2, λ; means = zeros(length(σ2)))
    for j in eachindex(r)
        m = -Inf
        for k in eachindex(σ2)
            Lik[j, k] = lognoise(N, j, r[j] - means[k], σ2[k], λ)
            m = max(m, Lik[j, k])
        end
        logscale[j] = m
        for k in eachindex(σ2)
            Lik[j, k] = exp(Lik[j, k] - m)
        end
    end
    return Lik
end

#penalised log-likelihood of the mixture weights (rows in free only)
function mixture_objective(π, Lik, logscale, free, prior)
    den = Lik * π
    obj = sum(k -> (prior[k] - 1) * log(max(π[k], 1.0e-300)), eachindex(π))
    for j in free
        obj += log(max(den[j], 1.0e-300)) + logscale[j]
    end
    return obj
end

#non-negative quadratic program min ½ y'Hy - b'y, y ≥ 0 (H positive definite), by the active-set method of
#Lawson & Hanson (1974): add the variable with the largest negative gradient, solve on the passive set, step back if needed
function nnqp(H, b; maxiter = 10 * length(b))
    #symmetric diagonal scaling (unit diagonal) keeps y ≥ 0 and removes most of the ill-conditioning
    D = 1 ./ sqrt.(diag(H))
    return D .* nnqp_scaled(D .* H .* D', D .* b, maxiter)
end
function nnqp_scaled(H, b, maxiter)
    K = length(b)
    tol = 1.0e-10 * max(maximum(abs, b), 1.0)
    passive = falses(K)
    y = zeros(K)
    w = copy(b) #minus the gradient at y
    for _ in 1:maxiter
        candidates = findall(k -> !passive[k] && w[k] > tol, 1:K)
        isempty(candidates) && break
        passive[candidates[argmax(w[candidates])]] = true
        while true
            P = findall(passive)
            z = zeros(K)
            z[P] .= Symmetric(H[P, P]) \ b[P]
            if all(>(0), z[P])
                y = z
                break
            end
            #move towards z until a passive variable reaches zero, and release it
            neg = filter(k -> z[k] <= 0, P)
            α = minimum(y[k] / (y[k] - z[k]) for k in neg)
            y .+= α .* (z .- y)
            for k in P
                if y[k] <= tol
                    y[k] = 0.0
                    passive[k] = false
                end
            end
        end
        w = b .- H * y
    end
    return y
end

#mix-SQP of Kim Y., Carbonetto P., Stephens M. & Anitescu M., J. Comput. Graph. Stat. (2020): https://doi.org/10.1080/10618600.2019.1689985
#maximise F(x) = Σⱼ log(Lⱼ'x) + Σₖ (priorₖ - 1) log(xₖ) - W Σₖ xₖ over x ≥ 0 (W = number of rows + Σₖ (priorₖ - 1)),
#whose maximiser lies on the simplex; each step solves the quadratic approximation as a non-negative QP, then a line search
function mixture_weights(Lik, logscale, free, prior; π0 = fill(1 / size(Lik, 2), size(Lik, 2)), tol = 1.0e-8, maxiter = 100)
    L = Lik[collect(free), :]
    K = size(L, 2)
    a = prior .- 1.0
    W = size(L, 1) + sum(a)
    F(x) = sum(log, max.(L * x, 1.0e-300)) + sum(k -> a[k] > 0 ? a[k] * log(x[k]) : 0.0, 1:K) - W * sum(x)
    x = 0.99 .* π0 ./ sum(π0) .+ 0.01 / K #strictly positive start
    Fx = F(x)
    for _ in 1:maxiter
        u = max.(L * x, 1.0e-300)
        #gradient of F and minus its Hessian (the Dirichlet terms only where aₖ > 0: other xₖ may reach zero)
        barrier = [a[k] > 0 ? a[k] / x[k] : 0.0 for k in 1:K]
        g = L' * (1 ./ u) .+ barrier .- W
        H = L' * (L ./ u .^ 2)
        H[diagind(H)] .+= barrier ./ max.(x, 1.0e-300) .+ 1.0e-8
        d = nnqp(H, H * x .+ g) .- x
        ascent = dot(g, d) #predicted ascent: zero at a KKT point
        ascent <= tol * W && break
        #step capped so the components with a Dirichlet weight stay positive, then Armijo backtracking
        t = 1.0
        for k in 1:K
            a[k] > 0 && d[k] < 0 && (t = min(t, 0.99 * x[k] / -d[k]))
        end
        Fnew = F(x .+ t .* d)
        while Fnew < Fx + 1.0e-4 * t * ascent && t > 1.0e-12
            t /= 2
            Fnew = F(x .+ t .* d)
        end
        t <= 1.0e-12 && break
        x = x .+ t .* d
        Fx = Fnew
    end
    π = x ./ sum(x)
    return π, mixture_objective(π, Lik, logscale, free, prior)
end

#posterior summaries of βⱼ given rⱼ ~ N(βⱼ, vⱼ) and the prior Σₖ πₖ N(mₖ, σₖ²)
#(with a t likelihood the component weights use the t density and the conditional posteriors the normal ones;
#with a NoiseModel the weights run over the pairs of prior and noise components); the lfsr counts the point mass at
#zero on both sides, a point mass at mₖ ≠ 0 on the side of its sign
ash_posterior(r, v, σ2, π; ν = Inf, means = zeros(length(σ2))) = ash_posterior(r, NoiseModel(v, ν), σ2, π; means = means)
function ash_posterior(r, N::NoiseModel, σ2, π; means = zeros(length(σ2)))
    p = length(r)
    K, M = length(π), size(N.w, 2)
    out = (mean = zeros(p), sd = zeros(p), lfsr = zeros(p), lfdr = zeros(p))
    w = zeros(K, M)
    for j in 1:p
        for k in 1:K, l in 1:M
            w[k, l] = π[k] == 0 || N.w[j, l] == 0 ? -Inf :
                log(π[k]) + log(N.w[j, l]) + logtpdf(r[j] - means[k], N.v[j, l] + σ2[k], N.ν[j, l], N.c[j, l])
        end
        w .= exp.(w .- maximum(w))
        w ./= sum(w)
        m1 = m2 = pneg = ppos = pzero = 0.0
        for k in 1:K, l in 1:M
            w[k, l] == 0 && continue
            if σ2[k] == 0
                mk = means[k]
                mk == 0 ? (pzero += w[k, l]) : mk > 0 ? (ppos += w[k, l]) : (pneg += w[k, l])
                m1 += w[k, l] * mk
                m2 += w[k, l] * mk^2
            else
                v = N.v[j, l]
                vk = σ2[k] * v / (σ2[k] + v)
                mk = (σ2[k] * r[j] + v * means[k]) / (σ2[k] + v)
                m1 += w[k, l] * mk
                m2 += w[k, l] * (vk + mk^2)
                pk = cdf(Normal(), -mk / sqrt(vk))
                pneg += w[k, l] * pk
                ppos += w[k, l] * (1 - pk)
            end
        end
        out.mean[j] = m1
        out.sd[j] = sqrt(max(m2 - m1^2, 0.0))
        out.lfsr[j] = min(pneg, ppos) + pzero
        out.lfdr[j] = pzero
    end
    return out
end

#q-values (from lfdr) or s-values (from lfsr): running mean of the sorted local rates
function cumulative_rate(x::AbstractVector{<:Real})
    o = sortperm(x)
    out = similar(x, Float64)
    acc = 0.0
    for (i, j) in enumerate(o)
        acc += x[j]
        out[j] = acc / i
    end
    return out
end

#prior components (variances, means): ash's point mass and scale mixture at zero, or with a mode μ ≠ 0 the point
#mass at zero (the null) followed by the grid of scales centred at μ (its first scale, 0, a point mass at μ):
#g = π₀ δ₀ + Σₖ πₖ N(μ, σₖ²), unimodal around μ apart from the null, as ash with mode = "estimate" but with the null
#kept at zero (where the controls put it)
mode_components(grid, μ) = μ == 0 ? (grid .^ 2, zeros(length(grid))) : ([0.0; grid .^ 2], [0.0; fill(float(μ), length(grid))])

#mixture weights fitted on the rows `rows` for the components (σ2, means), with the component likelihoods of all rows
function prior_fit(r, N::NoiseModel, λ, rows, σ2, means, nullweight)
    Lik = zeros(length(r), length(σ2))
    logscale = zeros(length(r))
    mixture_lik!(Lik, logscale, r, N, σ2, λ; means = means)
    π, obj = mixture_weights(Lik, logscale, rows, [nullweight; ones(length(σ2) - 1)])
    return π, obj, Lik, logscale
end

#mode of the non-null part: golden-section search of the marginal likelihood between the 1% and 99% quantiles of r
function estimate_mode(r, N::NoiseModel, λ, rows, grid, nullweight)
    obj(μ) = prior_fit(r, N, λ, rows, mode_components(grid, μ)..., nullweight)[2]
    g = (sqrt(5) - 1) / 2
    a, b = quantile(r[rows], [0.01, 0.99])
    tol = 1.0e-3 * (b - a)
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = obj(c), obj(d)
    while b - a > tol
        if fc > fd
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = obj(c)
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = obj(d)
        end
    end
    return (a + b) / 2
end

#gate of the mode: held-out log-likelihood gain per gene of the prior with an estimated mode over the prior at zero,
#by two-fold cross-validation over the genes (both priors fitted on one half, scored on the other); a mode fitted to
#noise or to residual confounding barely improves held-out genes (≤ 0.02 per gene on thinned GTEx, where pure nulls
#give ± 0.01), a global shift of transcription does (≥ 0.24 for a log2 shift of 1 in 60 % of the genes, 0.06–0.14
#for 0.5 in 40 %)
function mode_gain(r, N::NoiseModel, λ, rows, grid, nullweight; rng::AbstractRNG = Xoshiro(1))
    perm = shuffle(rng, collect(rows))
    h = length(perm) ÷ 2
    gain = 0.0
    for (tr, te) in ((perm[1:h], perm[(h + 1):end]), (perm[(h + 1):end], perm[1:h]))
        π0, _, L0, ls0 = prior_fit(r, N, λ, tr, mode_components(grid, 0.0)..., nullweight)
        μ = estimate_mode(r, N, λ, tr, grid, nullweight)
        π1, _, L1, ls1 = prior_fit(r, N, λ, tr, mode_components(grid, μ)..., nullweight)
        for j in te
            gain += log(max(dot(view(L1, j, :), π1), 1.0e-300)) + ls1[j] - log(max(dot(view(L0, j, :), π0), 1.0e-300)) - ls0[j]
        end
    end
    return gain / length(perm)
end

#prior mode of the non-null effects: :zero (ash), :estimate, or :auto (estimated when its held-out gain per gene
#exceeds threshold, see mode_gain); return the mode and the gain (NaN when not computed)
function choose_mode(mode::Symbol, r, N::NoiseModel, λ, rows, grid, nullweight, threshold)
    mode === :zero && return 0.0, NaN
    mode in (:auto, :estimate) || throw(ArgumentError("unknown mode $mode (use :auto, :zero or :estimate)"))
    gain = mode === :auto ? mode_gain(r, N, λ, rows, grid, nullweight) : NaN
    mode === :auto && gain <= threshold && return 0.0, gain
    return estimate_mode(r, N, λ, rows, grid, nullweight), gain
end

struct AshFit
    σ::Vector{Float64} #prior standard deviations (σ[1] = 0 is the point mass at zero)
    μ::Vector{Float64} #prior means (all zero unless a mode was estimated)
    π::Vector{Float64} #prior weights
    loglik::Float64 #penalised log-likelihood at the optimum
    mean::Vector{Float64}
    sd::Vector{Float64}
    lfsr::Vector{Float64}
    lfdr::Vector{Float64}
    svalue::Vector{Float64}
    qvalue::Vector{Float64}
end

"""
    ash(r, s; nullweight = 10, grid = ash_grid(r, s), df = Inf, mode = :auto, mode_threshold = 0.05)

ash for estimates r with standard errors s, nullweight is the Dirichlet weight of the point mass at zero
the non-null part of the prior is unimodal at zero (mode = :zero) or around an estimated mode μ, the null staying at
zero (mode = :estimate); mode = :auto estimates μ only when it raises the held-out log-likelihood by more than
mode_threshold per estimate (a shift shared by many effects, such as a global change of transcription)
"""
function ash(
        r::AbstractVector{<:Real}, s::AbstractVector{<:Real}; nullweight::Real = 10, grid = ash_grid(r, s), df = Inf,
        mode::Symbol = :auto, mode_threshold::Real = 0.05
    )
    @assert length(r) == length(s) && all(>(0), s)
    N = NoiseModel(s .^ 2, df)
    μ, _ = choose_mode(mode, r, N, 1.0, eachindex(r), grid, nullweight, mode_threshold)
    σ2, means = mode_components(grid, μ)
    π, obj, _, _ = prior_fit(r, N, 1.0, eachindex(r), σ2, means, nullweight)
    post = ash_posterior(r, N, σ2, π; means = means)
    return AshFit(sqrt.(σ2), means, π, obj, post.mean, post.sd, post.lfsr, post.lfdr, cumulative_rate(post.lfsr), cumulative_rate(post.lfdr))
end

"""
    evalues(r, v, σ, π)

compound e-values of Ramdas A. & Wang R., Hypothesis testing with e-values (2024), Theorem 9.21 / Ignatiadis N.,
Wang R. & Ramdas A., arXiv:2409.19812 (2024): for rⱼ ~ N(βⱼ, vⱼ) and the prior Σₖ πₖ N(0, σₖ²) (σ₁ = 0 the point
mass), the likelihood ratio of the non-null mixture to the null, eⱼ = qⱼ(rⱼ) / N(rⱼ; 0, vⱼ) with
qⱼ = Σₖ₌₂ πₖ N(0, vⱼ + σₖ²) / (1 - π₁), equivalently eⱼ = π₁ / (1 - π₁) × (1 - lfdrⱼ) / lfdrⱼ;
with an estimated prior they are asymptotically valid compound e-values
"""
evalues(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, σ::AbstractVector{<:Real}, π::AbstractVector{<:Real}; ν = Inf, means = zeros(length(σ))) =
    evalues(r, NoiseModel(v, ν), σ, π; means = means)
function evalues(r::AbstractVector{<:Real}, N::NoiseModel, σ::AbstractVector{<:Real}, π::AbstractVector{<:Real}; means = zeros(length(σ)))
    π0 = π[1]
    e = zeros(length(r))
    for j in eachindex(r)
        l0 = lognoise(N, j, r[j], 0.0)
        num = sum(π[k] * exp(lognoise(N, j, r[j] - means[k], σ[k]^2) - l0) for k in 2:length(π); init = 0.0)
        e[j] = num / max(1 - π0, 1.0e-300)
    end
    return e
end

"""
    ebh(e; level = 0.05)

e-BH procedure of Wang R. & Ramdas A., J. R. Stat. Soc. B (2022): reject the k* hypotheses with the largest
e-values, k* = max{k : k e₍ₖ₎ / K ≥ 1 / level}; controls the FDR under arbitrary dependence for compound e-values.
Return a Boolean vector of discoveries
"""
function ebh(e::AbstractVector{<:Real}; level::Real = 0.05)
    K = length(e)
    o = sortperm(e; rev = true)
    kstar = 0
    for k in 1:K
        k * e[o[k]] / K >= 1 / level && (kstar = k)
    end
    disc = falses(K)
    disc[o[1:kstar]] .= true
    return disc
end
