#number of latent factors by count splitting (Neufeld A., Gao L. L., Popp J., Battle A. & Witten D., Biostatistics
#(2024): https://doi.org/10.1093/biostatistics/kxac047): each count is split in two, the PLN model is fitted on the
#first part for each q and scored by the predictive log-likelihood of the second
#binomial thinning, Y₁ ~ Binomial(Y, ε) (the default), gives two parts independent given the latent values Zᵢⱼ, but
#they share the latent overdispersion: a fit that learns it is rewarded on the second part; with plnpca (random
#scores and loadings, factors that the data do not support switched off) this did not lead to too many factors on
#simulations, where it finds the true number in one pass
#beta-binomial thinning (data thinning, Neufeld A., Dharamshi A., Gao L. L. & Witten D., J. Mach. Learn. Res. (2024)),
#Y₁ | Y ~ BetaBinomial(Y, εr, (1 - ε)r), splits the latent noise too: for a count that is Poisson given a rate of mean
#m and variance v, Cov(Y₁, Y₂) = ε(1 - ε)(v r - m²) / (r + 1), zero at r = m² / v = 1 / CV², which is
#1 / (exp(σⱼ²) - 1) for the lognormal rate of gene j whatever its mean (binomial thinning is r → ∞); the parts are
#then uncorrelated (independent for a gamma rate), and an overestimated σⱼ² (r too small) makes them negatively
#correlated, which can only favour fewer factors; on GTEx it keeps fewer factors than the binomial split
#σⱼ² is taken from a fit at the current choice of q, starting from q = 0 (σ² overestimated: conservative) and
#repeated until the choice is stable; when the latent noise cannot be identified (very low counts) the choice stays
#at q = 0, whose σ² absorbs the factors' variance: the conservative side

#Gauss–Hermite nodes and weights for ∫ f(x) exp(-x²) dx (Golub–Welsch: eigen-decomposition of the Jacobi matrix)
function gausshermite(m::Integer)
    J = SymTridiagonal(zeros(m), [sqrt(k / 2) for k in 1:(m - 1)])
    E = eigen(J)
    return E.values, sqrt(π) .* E.vectors[1, :] .^ 2
end

#predictive log-likelihood of the test counts given a PLN fit on the training counts: each count is
#Poisson(exp(oᵢ + Zᵢⱼ)) with Zᵢⱼ ~ N(μᵢⱼ, σⱼ²) a fresh latent value (Gauss–Hermite quadrature)
function pln_predictive(Ytest, otest, μ, σ2; nodes = 20)
    x, w = gausshermite(nodes)
    lw = log.(w ./ sqrt(π))
    ll = 0.0
    for j in axes(Ytest, 2), i in axes(Ytest, 1)
        y = Ytest[i, j]
        term(k) = (η = otest[i] + μ[i, j] + sqrt(2 * σ2[j]) * x[k]; lw[k] + y * η - exp(η))
        m = maximum(term, eachindex(x))
        ll += m + log(sum(k -> exp(term(k) - m), eachindex(x))) - loggamma(y + 1)
    end
    return ll
end

#split of the counts: Y₁ | Y ~ BetaBinomial(Y, εrⱼ, (1 - ε)rⱼ) per gene (Binomial(Y, ε) when rⱼ is infinite)
function split_counts(rng, Y, r, ε)
    Y1 = similar(Y)
    for j in axes(Y, 2), i in axes(Y, 1)
        y = Y[i, j]
        Y1[i, j] = y == 0 ? 0 : isfinite(r[j]) ? rand(rng, BetaBinomial(y, ε * r[j], (1 - ε) * r[j])) : rand(rng, Binomial(y, ε))
    end
    return Y1, Y .- Y1
end

"""
    choose_q(Y, X; qs = 0:6, ε = 0.5, nsplit = 2, split = :binomial, maxrounds = 5, offset_factor = true, rng, kwargs...)

number of latent factors of the PLN model by count splitting: Y is split in two (binomial thinning; split =
:betabinomial for beta-binomial thinning with rⱼ = 1 / (exp(σⱼ²) - 1), which splits the latent noise as well as the
Poisson noise), the model is fitted on the first part for each q in qs and scored by the predictive log-likelihood
of the second part (averaged over nsplit splits); with the beta-binomial split σⱼ² comes from a fit at the current
choice, from q = 0, repeated until the choice is stable (at most maxrounds); kwargs are passed to plnpca
the automatic relevance determination of plnpca does not replace this choice: with a large q it keeps most factors
(real data have more weak structure than noise), and each one costs residual degrees of freedom; but above the
number of factors the data support, the chosen bound qmax no longer matters (the extra factors are switched off),
and the scores of such bounds differ only by noise: q is the number of factors kept by the fit at qmax
return q (factors kept), the chosen bound qmax, the scores of each bound in qs (last round) and the choices of
the successive rounds
"""
function choose_q(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}; qs = 0:6, ε::Real = 0.5, nsplit::Integer = 2,
        split::Symbol = :binomial, maxrounds::Integer = 5, offset_factor::Bool = true,
        rng::AbstractRNG = Random.default_rng(), kwargs...
    )
    @assert 0 < ε < 1
    @assert split in (:betabinomial, :binomial) "split must be :betabinomial or :binomial"
    n, k = size(X)
    qs = [q for q in qs if n - k > q + Int(offset_factor)]
    function scores(r)
        score = zeros(length(qs))
        for _ in 1:nsplit
            Ytrain, Ytest = split_counts(rng, Y, r, ε)
            #drop the genes with no training count in any sample (no information to fit them)
            keep = vec(sum(Ytrain; dims = 1)) .> 0
            Ytr, Yte = Ytrain[:, keep], Ytest[:, keep]
            otr = vec(log.(sum(Ytr; dims = 2)))
            ote = vec(log.(sum(Yte; dims = 2)))
            #bounds in increasing order, each fit started from the previous one with factors added from its residuals
            #(the same optimum as a fit from scratch, in a fraction of the sweeps)
            prev = nothing
            for l in sortperm(qs)
                q = qs[l]
                init = nothing
                if prev !== nothing && q > size(prev.state.g, 1)
                    init = (prev.state, prev.priors)
                    while size(init[1].g, 1) < q
                        init = plnpca_add_factor(init[1], init[2], float.(X))
                    end
                end
                f = plnpca_run(Ytr, X, q; offsets = otr, offset_factor = offset_factor, init = init, kwargs...)
                score[l] += pln_predictive(Yte, ote, latent_mean(f.state, float.(X)), f.state.b ./ f.state.a) / nsplit
                prev = f
            end
        end
        return score
    end
    if split === :binomial
        score = scores(fill(Inf, size(Y, 2)))
        qmax = qs[argmax(score)]
        return (q = kept(plnpca(Y, X, qmax; offset_factor = offset_factor, kwargs...), offset_factor), qmax = qmax, qs = qs, score = score, path = [qmax])
    end
    q, score, path = first(qs), Float64[], Int[]
    pilot = nothing
    for _ in 1:maxrounds
        #the σⱼ² of the split from a fit at the current choice
        pilot = plnpca(Y, X, q; offset_factor = offset_factor, kwargs...)
        score = scores(1 ./ expm1.(pilot.σ2))
        qnew = qs[argmax(score)]
        push!(path, qnew)
        qnew == q && break
        q = qnew
        pilot = nothing
    end
    #factors kept at the chosen bound (the pilot of the last round when the choice was stable)
    pilot === nothing && (pilot = plnpca(Y, X, q; offset_factor = offset_factor, kwargs...))
    return (q = kept(pilot, offset_factor), qmax = q, qs = qs, score = score, path = path)
end

#number of hidden factors kept by a fit (the rows of Γ after the offset factor)
kept(fit, offset_factor) = size(fit.Γ, 1) - Int(offset_factor)
