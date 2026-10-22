#number of latent factors by count splitting (Neufeld A., Gao L. L., Popp J., Battle A. & Witten D., Biostatistics
#(2024): https://doi.org/10.1093/biostatistics/kxac047): each count is split by binomial thinning,
#Y_train ~ Binomial(Y, ε) and Y_test = Y - Y_train, which are independent Poisson given the latent values
#Zᵢⱼ; the PLN model is fitted on Y_train for each q and scored by the predictive log-likelihood of Y_test
#the latent overdispersion εᵢⱼ is shared by the two parts, so an overfitted low-rank mean is still partly
#rewarded on Y_test: the criterion scores the mean μ only (the per-cell latent values M are not used) and
#integrates a fresh latent value, which removes most but not all of this optimism

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

"""
    choose_q(Y, X; qs = 0:6, ε = 0.5, nsplit = 2, offset_factor = true, rng, kwargs...)

number of latent factors of the PLN model by count splitting: Y is split by binomial thinning with probability ε,
the model is fitted on the first part for each q in qs and scored by the predictive log-likelihood of the second part
(averaged over nsplit splits); kwargs are passed to pln
return the chosen q and the scores of each q
"""
function choose_q(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}; qs = 0:6, ε::Real = 0.5, nsplit::Integer = 2,
        offset_factor::Bool = true, rng::AbstractRNG = Random.default_rng(), kwargs...
    )
    @assert 0 < ε < 1
    n, k = size(X)
    qs = [q for q in qs if n - k > q + Int(offset_factor)]
    score = zeros(length(qs))
    for _ in 1:nsplit
        Ytrain = [rand(rng, Binomial(y, ε)) for y in Y]
        Ytest = Y .- Ytrain
        #drop the genes with no training count in any sample (no information to fit them)
        keep = vec(sum(Ytrain; dims = 1)) .> 0
        Ytr, Yte = Ytrain[:, keep], Ytest[:, keep]
        otr = vec(log.(sum(Ytr; dims = 2)))
        ote = vec(log.(sum(Yte; dims = 2)))
        for (l, q) in enumerate(qs)
            fit = pln(Ytr, X, q; offsets = otr, offset_factor = offset_factor, kwargs...)
            score[l] += pln_predictive(Yte, ote, fit.μ, fit.σ2) / nsplit
        end
    end
    return (q = qs[argmax(score)], qs = qs, score = score)
end
