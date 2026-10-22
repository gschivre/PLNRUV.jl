#convergence diagnostics of Vehtari A., Gelman A., Simpson D., Carpenter B. & Bürkner P.-C., Bayesian Anal. (2021):
#https://doi.org/10.1214/20-BA1221 (rank-normalised split-R̂ and bulk effective sample size)
#the draws of one quantity are a matrix with one column per chain

#split each chain in two halves (detects chains that drift)
function split_chains(x::AbstractMatrix)
    h = size(x, 1) ÷ 2
    return hcat(x[1:h, :], x[(end - h + 1):end, :])
end

#normal scores of the pooled ranks (Blom's offset)
function rank_normalize(x::AbstractMatrix)
    N = length(x)
    z = similar(x, Float64)
    z[sortperm(vec(x))] .= quantile.(Normal(), ((1:N) .- 3 / 8) ./ (N + 1 / 4))
    return z
end

#within-chain variance and the pooled estimate of the posterior variance
function chain_variances(x::AbstractMatrix)
    n = size(x, 1)
    W = mean(var(x; dims = 1))
    B = n * var(vec(mean(x; dims = 1)))
    return W, (n - 1) / n * W + B / n
end

#classical R̂ on split chains
function rhat_basic(x::AbstractMatrix)
    W, varplus = chain_variances(x)
    return W > 0 ? sqrt(varplus / W) : NaN
end

"""
    rhat(draws)

rank-normalised split-R̂ (maximum of the bulk and folded versions) of draws with one column per chain;
values above 1.01 indicate that the chains have not mixed
"""
function rhat(x::AbstractMatrix)
    s = split_chains(x)
    return max(rhat_basic(rank_normalize(s)), rhat_basic(rank_normalize(abs.(s .- median(s)))))
end

#autocovariance at lag t of one chain (denominator n)
function autocov(x::AbstractVector, t::Integer)
    μ = mean(x)
    n = length(x)
    return sum((x[i] - μ) * (x[i + t] - μ) for i in 1:(n - t); init = 0.0) / n
end

"""
    ess(draws)

bulk effective sample size of draws with one column per chain: rank-normalised split chains,
autocorrelations combined across chains, and Geyer's initial monotone sequence
"""
function ess(x::AbstractMatrix)
    z = rank_normalize(split_chains(x))
    n, m = size(z)
    W, varplus = chain_variances(z)
    varplus > 0 || return NaN
    ρ(t) = 1 - (W - mean(autocov(view(z, :, c), t) for c in 1:m)) / varplus
    τ = -1.0
    previous = Inf
    for k in 0:((n - 2) ÷ 2)
        pair = ρ(2k) + ρ(2k + 1)
        pair <= 0 && break
        pair = min(pair, previous) #monotone sequence
        τ += 2pair
        previous = pair
    end
    return n * m / max(τ, 1 / log10(n * m))
end
