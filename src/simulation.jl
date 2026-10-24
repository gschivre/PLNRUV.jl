#data-based RNA-seq simulation by binomial thinning of Gerard D., BMC Bioinformatics (2020): https://doi.org/10.1186/s12859-020-3450-9
#thinning a mixed-Poisson count y ~ Poisson(λ) with probability q gives Poisson(q * λ),
#so the Poisson-lognormal model stays exact and the true log2 fold changes are known

"""
    surrogate_factors(counts, nfactors; pseudocount = 0.5)

surrogate variables of a count matrix (genes in rows, samples in columns) from the PCA of the doubly centred log counts
return a nsamples x nfactors matrix whose columns have mean 0 and mean square 1
"""
function surrogate_factors(counts::AbstractMatrix{<:Integer}, nfactors::Integer; pseudocount::Real = 0.5)
    n = size(counts, 2)
    @assert 0 <= nfactors < n "the number of factors must be smaller than the number of samples"
    L = log.(counts .+ pseudocount)
    L .-= mean(L; dims = 1) #centred log-ratio of each sample (removes sequencing depth)
    L .-= mean(L; dims = 2) #centre each gene so the factors are orthogonal to the intercept
    return svd(L).V[:, 1:nfactors] .* sqrt(n)
end

"""
    assign_groups(rng, U, target_cor; group_prop = 0.5)

assign samples to two groups so that the group indicator has a target correlation with the columns of U
a latent w = U * r + noise is thresholded at its quantile, and for a Gaussian w
cor(1{w > t}, U[:, k]) = r[k] * φ(t) / sqrt(p * (1 - p)) with t the (1 - p) quantile
"""
function assign_groups(rng::AbstractRNG, U::AbstractMatrix, target_cor::AbstractVector; group_prop::Real = 0.5)
    n = size(U, 1)
    @assert size(U, 2) == length(target_cor)
    t = quantile(Normal(), 1 - group_prop)
    r = target_cor .* (sqrt(group_prop * (1 - group_prop)) / pdf(Normal(), t))
    @assert sum(abs2, r) < 1 "target correlations too large for this group proportion"
    w = U * r .+ sqrt(1 - sum(abs2, r)) .* randn(rng, n)
    x = zeros(Int, n)
    x[partialsortperm(w, 1:round(Int, group_prop * n); rev = true)] .= 1
    return x
end

"""
    rerandomise_groups(rng, U, target_cor; group_prop = 0.5, margin = 0.15, maxtries = 1000)

assign samples to two groups at random, and keep the assignment only if its realised correlations with the columns
of U are within margin of the target: maxₗ ||cor(x, Uₗ)| - |target_corₗ|| ≤ margin (the sign of a factor is
arbitrary); rerandomisation (Morgan & Rubin 2012): first maxtries uniform random assignments (the assignment is then
uniform among those that meet the target), then, when none met it, maxtries proposals of assign_groups, which aims
at the target on average (with many samples a uniform split almost never reaches a strong target: the chance
correlations have a standard deviation of about 1 / √n)
return (x, tries), or nothing when no assignment met it in 2 maxtries draws (with few samples the possible
correlations are few, and a target may be out of reach)
"""
function rerandomise_groups(
        rng::AbstractRNG, U::AbstractMatrix, target_cor::AbstractVector; group_prop::Real = 0.5, margin::Real = 0.15,
        maxtries::Integer = 1000
    )
    n = size(U, 1)
    @assert size(U, 2) == length(target_cor)
    n1 = round(Int, group_prop * n)
    meets(x) = all(abs(abs(cor(x, u)) - abs(c)) <= margin for (u, c) in zip(eachcol(U), target_cor))
    for t in 1:maxtries
        x = zeros(Int, n)
        x[randperm(rng, n)[1:n1]] .= 1
        meets(x) && return (x, t)
    end
    #proposals from the latent threshold; a target that it cannot reach on average (e.g. one beyond the largest
    #correlation a binary split can have) gives up here
    sum(abs2, target_cor .* (sqrt(group_prop * (1 - group_prop)) / pdf(Normal(), quantile(Normal(), 1 - group_prop)))) < 1 ||
        return nothing
    for t in 1:maxtries
        x = assign_groups(rng, U, target_cor; group_prop = group_prop)
        meets(x) && return (x, maxtries + t)
    end
    return nothing
end

"""
    thin_counts([rng,] counts, X, B)

thin a count matrix (genes in rows, samples in columns) given a design X (samples x d, no intercept)
and log2 coefficients B (genes x d): gene j in sample i is thinned with probability 2^(s[i, j] - max_i(s[i, j])), s = X * B'
"""
function thin_counts(rng::AbstractRNG, counts::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real})
    p, n = size(counts)
    @assert size(X, 1) == n && size(B, 1) == p && size(X, 2) == size(B, 2)
    S = B * X'
    out = similar(counts, Int)
    for j in 1:p
        smax = maximum(view(S, j, :))
        for i in 1:n
            q = exp2(S[j, i] - smax)
            out[j, i] = q == 1.0 ? counts[j, i] : rand(rng, Binomial(counts[j, i], q))
        end
    end
    return out
end
thin_counts(counts::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}) =
    thin_counts(Random.default_rng(), counts, X, B)

"""
    thin_2group([rng,] counts; ngenes, nsamp, nctrl = 0, prop_null = 0.9, signal = Normal(0, 0.8), ctrl_signal = Dirac(0), group_prop = 0.5, target_cor = Float64[])

simulate a two-group differential expression data set from real counts (genes in rows, samples in columns)
ngenes genes and nsamp samples are drawn at random, nctrl genes are negative controls,
a proportion prop_null of the other genes are null and the log2 fold changes of the rest are drawn from signal
ctrl_signal allows controls that are only approximately null (e.g. Normal(0, 0.1))
target_cor sets the correlation between the group indicator and the first surrogate variables (confounding):
with assignment = :rerandomise (default) the realised correlations are within margin of it (rerandomise_groups; when
no assignment of a sample draw reaches it, new samples are drawn, at most maxdraws times), with :threshold they
match it only on average (assign_groups); an empty target_cor gives random groups; the realised correlations with
nrecord surrogate variables (at least length(target_cor)) are returned in all cases; the assignment uses its own
random stream, so data sets that differ only by the target do not share their assignment draws
return the thinned counts with samples in rows and genes in columns (the orientation used by the models)
"""
function thin_2group(
        rng::AbstractRNG, counts::AbstractMatrix{<:Integer};
        ngenes::Integer = size(counts, 1), nsamp::Integer = size(counts, 2), nctrl::Integer = 0,
        prop_null::Real = 0.9, signal::UnivariateDistribution = Normal(0.0, 0.8),
        ctrl_signal::UnivariateDistribution = Dirac(0.0),
        group_prop::Real = 0.5, target_cor::AbstractVector{<:Real} = Float64[], assignment::Symbol = :rerandomise,
        margin::Real = 0.15, maxtries::Integer = 1000, maxdraws::Integer = 20, nrecord::Integer = 2
    )
    p, n = size(counts)
    @assert ngenes <= p && nsamp <= n
    @assert 0 <= nctrl <= ngenes
    @assert assignment in (:rerandomise, :threshold)
    genes = randperm(rng, p)[1:ngenes]
    arng = Xoshiro(rand(rng, UInt64) ⊻ hash(target_cor))
    nf = min(max(nrecord, length(target_cor)), nsamp - 1)
    n1 = round(Int, group_prop * nsamp)

    #design: random groups, or groups whose realised correlation with the surrogate variables of the subsample meets
    #the target; a sample draw where the target is out of reach is replaced by another one (redraws counts them)
    local samples, sub, U, x
    tries, redraws = 1, 0
    while true
        samples = randperm(rng, n)[1:nsamp]
        sub = counts[genes, samples]
        U = surrogate_factors(sub, nf)
        if isempty(target_cor)
            x = shuffle(arng, [ones(Int, n1); zeros(Int, nsamp - n1)])
            break
        elseif assignment === :threshold
            x = assign_groups(arng, U[:, 1:length(target_cor)], target_cor; group_prop = group_prop)
            break
        end
        a = rerandomise_groups(arng, U[:, 1:length(target_cor)], target_cor; group_prop, margin, maxtries)
        if a !== nothing
            x, tries = a
            break
        end
        redraws += 1
        redraws > maxdraws && error("no group assignment within $margin of the target $target_cor in $maxdraws sample draws")
    end
    realized_cor = [cor(x, u) for u in eachcol(U)]

    #log2 fold changes: controls first, then non-null, then null genes, in a random gene order
    nsignal = round(Int, (ngenes - nctrl) * (1 - prop_null))
    β = [rand(rng, ctrl_signal, nctrl); rand(rng, signal, nsignal); zeros(ngenes - nctrl - nsignal)]
    role = [fill(:control, nctrl); fill(:signal, nsignal); fill(:null, ngenes - nctrl - nsignal)]
    perm = randperm(rng, ngenes)
    β = β[perm]
    role = role[perm]

    out = thin_counts(rng, sub, reshape(x, :, 1), reshape(β, :, 1))

    return (
        counts = permutedims(out), X = [ones(Int, nsamp) x], β = β,
        controls = findall(==(:control), role), signal = findall(==(:signal), role), null = findall(==(:null), role),
        genes = genes, samples = samples, factors = U, realized_cor = realized_cor,
        tries = tries, redraws = redraws,
    )
end
thin_2group(counts::AbstractMatrix{<:Integer}; kwargs...) = thin_2group(Random.default_rng(), counts; kwargs...)

"""
    filter_genes(counts; min_mean = 0.0)

index of the genes (rows) that are observed at least once and with a mean count ≥ min_mean
"""
function filter_genes(counts::AbstractMatrix{<:Integer}; min_mean::Real = 0.0)
    return [i for (i, c) in enumerate(eachrow(counts)) if any(>(0), c) && mean(c) >= min_mean]
end
