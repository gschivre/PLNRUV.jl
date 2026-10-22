#cross-check of the binomial thinning simulator against the R package seqgendiff
#run from the benchmark environment (julia --project=benchmark benchmark/seqgendiff_check.jl <gtex counts file>)
using PLNRUV, RCall, Random, Statistics

R"suppressPackageStartupMessages(library(seqgendiff))"

#1. thinning probabilities for a fixed design and fixed log2 coefficients
#with constant large counts the thinned / original ratio is the thinning probability up to O(1e-3)
function check_thin_counts(; p = 300, n = 20, d = 2, depth = 10^7, seed = 1)
    rng = Xoshiro(seed)
    counts = fill(depth, p, n)
    X = [repeat([0, 1], n ÷ 2) randn(rng, n, d - 1)]
    B = 1.5 .* randn(rng, p, d)
    ours = thin_counts(rng, counts, X, B) ./ depth
    R"set.seed($seed)"
    theirs = rcopy(R"thin_diff(mat = $counts, design_fixed = $X, coef_fixed = $B)$mat") ./ depth
    return maximum(abs, ours - theirs)
end

#2. realized |correlation| between the group indicator and the surrogate variables
#seqgendiff's surrogate variables (log2 counts without per-sample centring, so they include sequencing depth)
#are given to our assign_groups so both designs are compared on the same factors
#(PLNRUV itself centres each sample: confounding comes from biological / unwanted factors, not depth)
function check_confounding(counts; ngenes = 2_000, nsamp = 100, target_cor = [0.6, 0.2], nrep = 30, seed = 2)
    rng = Xoshiro(seed)
    R"set.seed($seed)"
    ours = zeros(nrep, length(target_cor))
    theirs = zeros(nrep, length(target_cor))
    for r in 1:nrep
        sub = counts[randperm(rng, size(counts, 1))[1:ngenes], randperm(rng, size(counts, 2))[1:nsamp]]
        out = R"thin_2group(mat = $sub, prop_null = 1, group_prop = 0.5, corvec = $target_cor)"
        x = rcopy(R"$out$designmat[, 1]")
        sv = rcopy(R"$out$sv")
        sv ./= sqrt.(mean(abs2, sv .- mean(sv; dims = 1); dims = 1)) #unit variance as assign_groups expects
        ours[r, :] .= abs.([cor(assign_groups(rng, sv, target_cor), u) for u in eachcol(sv)])
        theirs[r, :] .= abs.([cor(x, u) for u in eachcol(sv)])
    end
    return (
        target = target_cor, ours = vec(mean(ours; dims = 1)), ours_sd = vec(std(ours; dims = 1)),
        seqgendiff = vec(mean(theirs; dims = 1)), seqgendiff_sd = vec(std(theirs; dims = 1)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    gtex = read_gct(ARGS[1])
    M = gtex.counts[filter_genes(gtex.counts; min_mean = 10), :]
    println("max |q_PLNRUV - q_seqgendiff| = ", check_thin_counts())
    println(check_confounding(M))
end
