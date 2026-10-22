#empirical Bayes shrinkage of variances: the sample variances follow s²ⱼ ~ σ²ⱼ χ²_d / d

#inverse of the trigamma function by Newton iterations (Smyth G.K., Stat. Appl. Genet. Mol. Biol. (2004))
function trigamma_inverse(x::Real)
    x > 1.0e7 && return 1 / sqrt(x)
    x < 1.0e-6 && return 1 / x
    y = 0.5 + 1 / x
    for _ in 1:50
        tri = trigamma(y)
        Δ = tri * (1 - tri / x) / polygamma(2, y)
        y += Δ
        -Δ / y < 1.0e-8 && break
    end
    return y
end

"""
    squeeze_var(s2, d)

limma moderation of Smyth G.K., Stat. Appl. Genet. Mol. Biol. (2004): https://doi.org/10.2202/1544-6115.1027
prior σ²ⱼ ~ scaled inverse χ²(d₀, s₀²) estimated by the method of moments on log s²
return the posterior variances (d₀ s₀² + d s²ⱼ) / (d₀ + d) and the prior (d₀, s₀²)
"""
function squeeze_var(s2::AbstractVector{<:Real}, d::Real)
    p = length(s2)
    e = log.(s2) .- digamma(d / 2) .+ log(d / 2)
    ē = mean(e)
    v = sum(abs2, e .- ē) / (p - 1) - trigamma(d / 2)
    if v <= 0 #no evidence of variance heterogeneity: complete pooling
        s02 = exp(ē + digamma(d / 2) - log(d / 2))
        return (s2 = fill(s02, p), d0 = Inf, s02 = s02, df = fill(Inf, p))
    end
    d0 = 2 * trigamma_inverse(v)
    s02 = exp(ē + digamma(d0 / 2) - log(d0 / 2))
    return (s2 = (d0 * s02 .+ d .* s2) ./ (d0 + d), d0 = d0, s02 = s02, df = fill(d0 + d, p))
end

"""
    vash(s2, d; shapes = 2.0 .^ (-2:0.5:12))

vash of Lu M. & Stephens M., Bioinformatics (2016): https://doi.org/10.1093/bioinformatics/btw483
prior on the variances: mixture of inverse gammas σ²ⱼ ~ Σₖ πₖ InvGamma(aₖ, aₖ c) with a common scale c
and a grid of shapes aₖ (limma is the one-component case), so that s²ⱼ / c ~ Σₖ πₖ F(d, 2aₖ)
π by mix-SQP (as in ash), c by golden-section search; return the posterior 1 / E[1 / σ²ⱼ]
and the posterior mixture over the shapes with nonzero weight (mix: weights w, squared scales s2 and degrees of
freedom df, genes x components), from which the exact likelihood of a moderated estimate is a mixture of t
"""
function vash(s2::AbstractVector{<:Real}, d::Real; shapes = 2.0 .^ (-2:0.5:12))
    p = length(s2)
    K = length(shapes)
    Lik = zeros(p, K)
    logscale = zeros(p)
    prior = ones(K)
    #log density of s² / c ~ F(d, d₂) with d₂ = 2aₖ, including the Jacobian 1 / c:
    #(d / 2) log(d / d₂) - logB(d / 2, d₂ / 2) + (d / 2 - 1) log(s² / c) - ((d + d₂) / 2) log(1 + d s² / (c d₂)) - log(c)
    d2 = 2 .* shapes
    const_k = (d / 2) .* log.(d ./ d2) .- logbeta.(d / 2, d2 ./ 2)
    logs2 = log.(s2)
    πwarm = fill(1 / K, K) #warm start across the evaluations of the search
    function fit_weights(logc)
        c = exp(logc)
        for j in 1:p
            common = (d / 2 - 1) * (logs2[j] - logc) - logc
            m = -Inf
            for k in 1:K
                Lik[j, k] = const_k[k] + common - ((d + d2[k]) / 2) * log1p(d * s2[j] / (c * d2[k]))
                m = max(m, Lik[j, k])
            end
            logscale[j] = m
            for k in 1:K
                Lik[j, k] = exp(Lik[j, k] - m)
            end
        end
        π, obj = mixture_weights(Lik, logscale, 1:p, prior; π0 = πwarm)
        πwarm .= π
        return π, obj
    end
    #golden-section search of the common scale between the 1% and 99% quantiles of the sample variances
    g = (sqrt(5) - 1) / 2
    a, b = log.(quantile(s2, [0.01, 0.99]))
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = fit_weights(x1)[2], fit_weights(x2)[2]
    while b - a > 1.0e-4
        if f1 > f2
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a)
            f1 = fit_weights(x1)[2]
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a)
            f2 = fit_weights(x2)[2]
        end
    end
    logc = (a + b) / 2
    π, _ = fit_weights(logc)
    c = exp(logc)
    #posterior of the precisions: Σₖ wⱼₖ Gamma(aₖ + d / 2, rate aₖ c + d s²ⱼ / 2)
    post = zeros(p)
    df = zeros(p) #effective degrees of freedom: the predictive is a mixture of t with 2aₖ + d degrees of freedom
    active = findall(>(1.0e-8), π)
    mix = (w = zeros(p, length(active)), s2 = zeros(p, length(active)), df = zeros(p, length(active)))
    for j in 1:p
        w = view(Lik, j, :) .* π
        w ./= sum(w)
        post[j] = 1 / sum(w[k] * (shapes[k] + d / 2) / (shapes[k] * c + d * s2[j] / 2) for k in 1:K)
        df[j] = sum(w[k] * (2 * shapes[k] + d) for k in 1:K)
        for (l, k) in enumerate(active)
            mix.w[j, l] = w[k]
            mix.s2[j, l] = (shapes[k] * c + d * s2[j] / 2) / (shapes[k] + d / 2)
            mix.df[j, l] = 2 * shapes[k] + d
        end
        mix.w[j, :] ./= sum(view(mix.w, j, :))
    end
    return (s2 = post, π = π, c = c, shapes = collect(shapes), df = df, mix = mix)
end

"""
    variance_trend(s2, x)

trend of the variances in a covariate x (e.g. the mean log-abundance of each gene, as in limma-trend):
a cubic polynomial fitted to log s² by least squares, evaluated at x clamped to its 1% and 99% quantiles
(no extrapolation of the cubic in the tails); return the trend on the variance scale
the constant bias of E[log s²] is not removed, it is absorbed by the scale of the prior (s₀² or c)
"""
function variance_trend(s2::AbstractVector{<:Real}, x::AbstractVector{<:Real})
    lo, hi = quantile(x, [0.01, 0.99])
    xc = clamp.(x, lo, hi)
    z = (xc .- mean(xc)) ./ max(std(xc), eps()) #standardised for a well-conditioned design
    A = [ones(length(z)) z z .^ 2 z .^ 3]
    return exp.(A * (A \ log.(s2)))
end

#dispatch on the shrinkage method (:vash, :limma or :none): the shrunk variances, their degrees of freedom and the
#posterior mixture (one component for limma and none, the active shapes for vash)
#with a covariate the variances are shrunk around their trend in it: method applied to s² / trend, then rescaled
#(an exchangeable prior pulls the overdispersions of low-count genes down and those of high-count genes up)
function shrink_variances(s2::AbstractVector{<:Real}, d::Real, method::Symbol; covariate = nothing)
    single(s2, df) = (s2 = s2, df = df, mix = (w = ones(length(s2), 1), s2 = reshape(s2, :, 1), df = reshape(df, :, 1)))
    method === :none && return single(collect(float.(s2)), fill(float(d), length(s2)))
    t = covariate === nothing ? ones(length(s2)) : variance_trend(s2, covariate)
    if method === :limma
        lim = squeeze_var(s2 ./ t, d)
        return single(lim.s2 .* t, lim.df)
    elseif method === :vash
        va = vash(s2 ./ t, d)
        return (s2 = va.s2 .* t, df = va.df, mix = (w = va.mix.w, s2 = va.mix.s2 .* t, df = va.mix.df))
    end
    throw(ArgumentError("unknown variance shrinkage method $method (use :vash, :limma or :none)"))
end
