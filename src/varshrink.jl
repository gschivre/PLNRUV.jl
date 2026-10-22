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
        return (s2 = fill(s02, p), d0 = Inf, s02 = s02)
    end
    d0 = 2 * trigamma_inverse(v)
    s02 = exp(ē + digamma(d0 / 2) - log(d0 / 2))
    return (s2 = (d0 * s02 .+ d .* s2) ./ (d0 + d), d0 = d0, s02 = s02)
end

"""
    vash(s2, d; shapes = 2.0 .^ (-2:0.5:12))

vash of Lu M. & Stephens M., Bioinformatics (2016): https://doi.org/10.1093/bioinformatics/btw483
prior on the variances: mixture of inverse gammas σ²ⱼ ~ Σₖ πₖ InvGamma(aₖ, aₖ c) with a common scale c
and a grid of shapes aₖ (limma is the one-component case), so that s²ⱼ / c ~ Σₖ πₖ F(d, 2aₖ)
π by mix-SQP (as in ash), c by golden-section search; return the posterior 1 / E[1 / σ²ⱼ]
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
    for j in 1:p
        w = view(Lik, j, :) .* π
        w ./= sum(w)
        post[j] = 1 / sum(w[k] * (shapes[k] + d / 2) / (shapes[k] * c + d * s2[j] / 2) for k in 1:K)
    end
    return (s2 = post, π = π, c = c, shapes = collect(shapes))
end

#dispatch on the shrinkage method (:vash, :limma or :none)
function shrink_variances(s2::AbstractVector{<:Real}, d::Real, method::Symbol)
    method === :none && return collect(float.(s2))
    method === :limma && return squeeze_var(s2, d).s2
    method === :vash && return vash(s2, d).s2
    throw(ArgumentError("unknown variance shrinkage method $method (use :vash, :limma or :none)"))
end
