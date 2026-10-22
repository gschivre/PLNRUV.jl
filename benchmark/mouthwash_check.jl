#cross-check of the Gaussian MOUTHWASH against vicar::mouthwash (Gerard & Stephens)
#run from the benchmark environment (julia --project=benchmark benchmark/mouthwash_check.jl)
using PLNRUV, RCall, Random, Distributions, Statistics

R"suppressPackageStartupMessages(library(vicar))"

#Gaussian log-abundance with confounders correlated with the group and an optional compositional shift
function gauss_sim(rng; n = 40, p = 1_000, c = 100, q = 2, shift = 0.0)
    X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
    β = [randn(rng, 1, p); [zeros(1, c) (rand(rng, 1, p - c) .< 0.3) .* randn(rng, 1, p - c) .* 1.5]]
    Z = X * [0.0 0.0; 1.2 -0.8] + randn(rng, n, q)
    A = X * β + Z * randn(rng, q, p) + randn(rng, n, p) .* sqrt.(rand(rng, Gamma(2, 0.1), p))'
    A .+= shift .* X[:, 2] .+ 0.3 .* randn(rng, n)
    return (Y = A, X = X, β = β[2, :], controls = 1:c)
end

#RMSE, bias, discoveries at svalue < 0.05 and their false discovery proportion (null or wrong sign)
function eval_fit(m, sv, β, idx)
    disc = idx[sv[idx] .< 0.05]
    return (
        rmse = sqrt(mean(abs2, m[idx] - β[idx])), bias = mean(m[idx] - β[idx]), ndisc = length(disc),
        fdp = isempty(disc) ? 0.0 : mean((β[disc] .== 0) .| (sign.(β[disc]) .!= sign.(m[disc]))),
    )
end

#vicar settings matched to PLNRUV: normal likelihood, normal mixture, no limma variance shrinkage
function vicar_mouthwash(Y, X, q)
    out = R"""mouthwash(Y = $Y, X = $X, k = $q, cov_of_interest = 2, include_intercept = FALSE,
    likelihood = "normal", mixing_dist = "normal", limmashrink = FALSE, verbose = FALSE)"""
    return (mean = rcopy(R"$out$result$PosteriorMean"), svalue = rcopy(R"$out$result$svalue"))
end

#1. same model (q = 3 free factors, no offset factor): both implementations should agree
#2. PLNRUV with the offset factor vs vicar, with and without a compositional shift correlated with the group
function check_mouthwash(; nrep = 3, seed = 300)
    rows = []
    for r in 1:nrep, shift in (0.0, 1.0)
        d = gauss_sim(Xoshiro(seed + r); shift = shift)
        idx = setdiff(1:length(d.β), d.controls)
        same_ours = mouthwash(d.Y, d.X, 3; offset_factor = false)
        same_vicar = vicar_mouthwash(d.Y, d.X, 3)
        ours = mouthwash(d.Y, d.X, 2)
        vicar2 = vicar_mouthwash(d.Y, d.X, 2)
        push!(
            rows, (
                rep = r, shift = shift,
                same_model = (ours = eval_fit(same_ours.mean, same_ours.svalue, d.β, idx), vicar = eval_fit(same_vicar.mean, same_vicar.svalue, d.β, idx)),
                offset_vs_vicar = (ours = eval_fit(ours.mean, ours.svalue, d.β, idx), vicar = eval_fit(vicar2.mean, vicar2.svalue, d.β, idx)),
            ),
        )
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(println, check_mouthwash())
end
