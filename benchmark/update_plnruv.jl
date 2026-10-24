#refresh the Julia rows (PLNRUV and the oracles) of the 20-sample benchmarks after a change of the package, keeping the
#saved outputs of the R methods (results/<name>_datasets/*.jls); the data sets are simulated again from their seeds
#and checked against the saved truth
#run from the benchmark environment: julia -t auto --project=benchmark benchmark/update_plnruv.jl <gtex counts file>
include(joinpath(@__DIR__, "grid_benchmark.jl"))

is_julia_method(m) = startswith(m, "PLNRUV_") || startswith(m, "oracle_")

function refresh_scenarios(M, scenarios, file; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1)
    dir = replace(file, ".csv" => "_datasets")
    names = Set(sc.name for sc in scenarios)
    lines = readlines(file)
    rows = filter(l -> split(l, ",")[1] ∉ names, lines[2:end])
    cfile = replace(file, ".csv" => "_curves.csv")
    clines = readlines(cfile)
    crows = filter(l -> split(l, ",")[1] ∉ names, clines[2:end])
    for sc in scenarios, r in 1:nrep
        key = "$(sc.name)_r$(r)"
        saved = deserialize(joinpath(dir, key * ".jls"))
        sim = thin_2group(
            Xoshiro(seed + r), M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = sc.prop_null,
            signal = sc.signal, ctrl_signal = sc.ctrl, target_cor = sc.target_cor
        )
        @assert sim.β == saved.β && sim.controls == saved.controls "data set $key differs from the saved one"
        res = Dict(m => (v.lfc, v.q) for (m, v) in saved.methods if !is_julia_method(m))
        merge!(res, safely(julia_methods, sim.counts, sim.X, sim.controls, q))
        merge!(res, safely(oracle_methods, sim, sc))
        idx = setdiff(1:ngenes, sim.controls)
        meancount = vec(mean(sim.counts; dims = 1))
        for (method, (lfc, qv)) in sort(collect(res); by = first)
            x = metrics(lfc, qv, sim.β, idx)
            push!(rows, join((sc.name, r, method, x.ndisc, x.fdp, x.power, x.rmse), ","))
            append!(crows, curve_rows(join((sc.name, r, method), ","), lfc, qv, sim.β, meancount, idx))
        end
        save_dataset(dir, key, sim, res)
        logline("refresh $(basename(file)) $key")
    end
    write(file, join([lines[1]; rows], "\n") * "\n")
    write(cfile, join([clines[1]; crows], "\n") * "\n")
    return nothing
end

function refresh_all(M)
    dir = joinpath(@__DIR__, "results")
    refresh_scenarios(M, MAIN_SCENARIOS, joinpath(dir, "gtex_benchmark.csv"))
    refresh_scenarios(M, CTRL_SCENARIOS, joinpath(dir, "gtex_controls.csv"))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    gtex = read_gct(ARGS[1])
    refresh_all(gtex.counts[filter_genes(gtex.counts; min_mean = 10), :])
end
