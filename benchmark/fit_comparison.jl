#the PLN fit (plnpca) compared on the same data sets with two variance priors and two likelihoods, all with negative
#controls, on binomially thinned GTEx data (scenarios and seeds of the grid, signal 0.8); one .jls per data set in
#results/fit_comparison/ (truth and each fit's log2 fold changes and s-values); the saved data sets also hold the
#outputs of the older fit (pln, fixed scores and loadings, removed from the package after this comparison)
#run from the benchmark environment, after include("grid_benchmark.jl") and loading the GTEx counts M:
#    include("benchmark/fit_comparison.jl"); run_fit_comparison(M)
#figures: benchmark/figures/fit_comparison_figures.jl
using Serialization

const CMP_FITS = [
    ("plnpca_trend", (; var_prior = :trend)),
    ("plnpca_common", (; var_prior = :common)),
    ("plnpca_trend_normal", (; var_prior = :trend, likelihood = :normal)),
]
const CMP_NPG = [2, 3, 5, 10]
const CMP_SCENARIOS = ["balanced", "confounded", "global_increase"]

#nrep data sets per cell; a balanced design out of reach (NA) is replaced by the next seed, at most maxseeds seeds;
#a data set saved earlier only gets the fits it lacks
function run_fit_comparison(
        M; nrep = 5, maxseeds = 30, ngenes = 2_000, nctrl = 200, seed = 1,
        dir = joinpath(@__DIR__, "results", "fit_comparison")
    )
    mkpath(dir)
    scs = Dict(sc.name => sc for sc in grid_scenarios(0.8))
    for g in CMP_NPG, name in CMP_SCENARIOS
        sc = scs[name]
        found = 0
        for r in 1:maxseeds
            found >= nrep && break
            STOP[] && return nothing
            key = "$(name)_n$(g)_r$(r)"
            path = joinpath(dir, key * ".jls")
            saved = isfile(path) ? deserialize(path) : nothing
            if saved !== nothing && all(haskey(saved.fits, f) for (f, _) in CMP_FITS)
                found += 1
                continue
            end
            sim = simulate(seed + 100g + r, M, sc; ngenes = ngenes, nsamp = 2g, nctrl = nctrl)
            if sim === nothing
                logline("fit comparison $key: NA (target out of reach)")
                continue
            end
            saved === nothing || @assert saved.β == sim.β "data set $key differs from the saved one"
            q = nfactors(2g)
            fits = saved === nothing ? Dict{String, Any}() : saved.fits
            for (fname, kw) in CMP_FITS
                haskey(fits, fname) && continue
                t = @elapsed out = try
                    fit = pln_mouthwash(sim.counts, sim.X, q; controls = sim.controls, kw...)
                    (lfc = fit.mouthwash.mean ./ log(2), s = fit.mouthwash.svalue, converged = fit.pln.converged)
                catch err
                    err isa InterruptException && rethrow()
                    (lfc = fill(NaN, ngenes), s = fill(NaN, ngenes), converged = false, error = sprint(showerror, err))
                end
                fits[fname] = merge(out, (time = t,))
            end
            serialize(
                path, (
                    β = sim.β, controls = sim.controls, meancount = vec(mean(sim.counts; dims = 1)), q = q,
                    realized_cor = sim.realized_cor, fits = fits,
                )
            )
            found += 1
            logline("fit comparison $key")
        end
    end
    logline("fit comparison done")
    return nothing
end
