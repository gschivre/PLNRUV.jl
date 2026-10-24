#refresh the Julia rows (PLNRUV and the oracles) of the 20-sample benchmarks after a change of the package, keeping the
#saved outputs of the R methods (results/<name>_datasets/*.jls); the data sets are simulated again from their seeds
#and checked against the saved truth
#run from the benchmark environment: julia -t auto --project=benchmark benchmark/update_plnruv.jl <gtex counts file>
include(joinpath(@__DIR__, "grid_benchmark.jl"))

is_julia_method(m) = startswith(m, "PLNRUV_") || startswith(m, "oracle_")

#the surrogate variables of a saved data set, recomputed from its genes and samples (deterministic): data sets saved
#before thin_2group returned them for every scenario lack them for random groups
function saved_factors(M, saved)
    U = surrogate_factors(M[saved.genes, saved.samples], 2)
    size(saved.factors, 2) > 0 && @assert U[:, 1:size(saved.factors, 2)] ≈ saved.factors "surrogate variables differ from the saved ones"
    return U
end

function refresh_scenarios(M, scenarios, file; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1)
    dir = replace(file, ".csv" => "_datasets")
    names = Set(sc.name for sc in scenarios)
    lines = readlines(file)
    rows = filter(l -> split(l, ",")[1] ∉ names, lines[2:end])
    cfile = replace(file, ".csv" => "_curves.csv")
    clines = readlines(cfile)
    crows = filter(l -> split(l, ",")[1] ∉ names, clines[2:end])
    for sc in scenarios, r in 1:nrep
        STOP[] && return nothing
        key = "$(sc.name)_r$(r)"
        path = joinpath(dir, key * ".jls")
        isfile(path) || continue #NA: target confounding out of reach
        saved = deserialize(path)
        #the Julia methods run on the saved data set, the one the R methods saw; the simulation from the seed is only a
        #check (the random streams of Random and Distributions may change between versions)
        @assert hasproperty(saved, :counts) "data set $key was saved without its counts: run the benchmark again"
        sim = (
            counts = saved.counts, X = saved.X, β = saved.β, controls = saved.controls, null = saved.null,
            signal = saved.signal, factors = saved_factors(M, saved), genes = saved.genes, samples = saved.samples,
            realized_cor = saved.realized_cor, tries = saved.tries, redraws = saved.redraws,
        )
        again = simulate(seed + r, M, sc; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl)
        if again === nothing || again.counts != sim.counts || again.X != sim.X
            logline("refresh $key: the simulation from its seed differs from the saved data set (the saved one is used)")
        end
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

#the Julia rows (PLNRUV and the oracles) of the grid, on its saved data sets (results/grid_benchmark_datasets): the
#previous Julia rows are dropped and recomputed, the R rows kept; one data set at a time, appended
function refresh_grid(M, file = joinpath(@__DIR__, "results", "grid_benchmark.csv"))
    dir = replace(file, ".csv" => "_datasets")
    cfile = replace(file, ".csv" => "_curves.csv")
    keep(l) = !is_julia_method(split(l, ",")[6])
    lines, clines = readlines(file), readlines(cfile)
    write(file, join([lines[1]; filter(keep, lines[2:end])], "\n") * "\n")
    write(cfile, join([clines[1]; filter(l -> !is_julia_method(split(l, ",")[5]), clines[2:end])], "\n") * "\n")
    scs = Dict((sc.name, s) => sc for s in SIGNAL for sc in grid_scenarios(s))
    for f in sort(readdir(dir))
        STOP[] && return nothing
        m = match(r"^(.*)_n(\d+)_s([\d.]+)_r(\d+)\.jls$", f)
        name, g, s, r = m[1], parse(Int, m[2]), parse(Float64, m[3]), parse(Int, m[4])
        saved = deserialize(joinpath(dir, f))
        sim = (
            counts = saved.counts, X = saved.X, β = saved.β, controls = saved.controls, null = saved.null,
            signal = saved.signal, factors = saved_factors(M, saved), genes = saved.genes, samples = saved.samples,
            realized_cor = saved.realized_cor, tries = saved.tries, redraws = saved.redraws,
        )
        sc = scs[(name, s)]
        q = nfactors(2g)
        res = Dict(k => (v.lfc, v.q) for (k, v) in saved.methods if !is_julia_method(k))
        merge!(res, safely(julia_methods, sim.counts, sim.X, sim.controls, q))
        merge!(res, safely(oracle_methods, sim, sc))
        idx = setdiff(1:length(sim.β), sim.controls)
        meancount = vec(mean(sim.counts; dims = 1))
        rows, crows = String[], String[]
        for (method, (lfc, qv)) in res
            is_julia_method(method) || continue
            x = metrics(lfc, qv, sim.β, idx)
            push!(rows, join((name, g, s, r, q, method, x.ndisc, x.fdp, x.power, x.rmse), ","))
            append!(crows, curve_rows(join((name, g, s, r, method), ","), lfc, qv, sim.β, meancount, idx))
        end
        open(io -> foreach(l -> println(io, l), rows), file, "a")
        open(io -> foreach(l -> println(io, l), crows), cfile, "a")
        save_dataset(dir, replace(f, ".jls" => ""), sim, res)
        logline("refresh grid $f")
    end
    return nothing
end

#PLNRUV with q chosen by choose_q (binomial split) on a subset of the saved grid data sets: rows
#PLNRUV_EB_controls_chooseq and PLNRUV_EB_sparsity_chooseq, with the chosen q in the q column; run after refresh_grid
#(which drops all Julia rows)
function chooseq_grid(
        file = joinpath(@__DIR__, "results", "grid_benchmark.csv"); npergroup = [5, 10],
        scenarios = ["randomised", "confounded"], signal = 0.8, nrep = 10
    )
    dir = replace(file, ".csv" => "_datasets")
    cfile = replace(file, ".csv" => "_curves.csv")
    for g in npergroup, name in scenarios, r in 1:nrep
        STOP[] && return nothing
        key = "$(name)_n$(g)_s$(signal)_r$(r)"
        path = joinpath(dir, key * ".jls")
        isfile(path) || continue
        saved = deserialize(path)
        t = @elapsed begin
            c = choose_q(saved.counts, saved.X; qs = 0:min(5, 2g - 4), rng = Xoshiro(r))
            fits = Dict(
                "PLNRUV_EB_controls_chooseq" => pln_mouthwash(saved.counts, saved.X, c.q; controls = saved.controls).mouthwash,
                "PLNRUV_EB_sparsity_chooseq" => pln_mouthwash(saved.counts, saved.X, c.q).mouthwash,
            )
        end
        idx = setdiff(1:length(saved.β), saved.controls)
        rows, crows = String[], String[]
        methods = Dict(saved.methods)
        for (m, mw) in fits
            lfc, qv = mw.mean ./ log(2), mw.svalue
            x = metrics(lfc, qv, saved.β, idx)
            push!(rows, join((name, g, signal, r, c.q, m, x.ndisc, x.fdp, x.power, x.rmse), ","))
            append!(crows, curve_rows(join((name, g, signal, r, m), ","), lfc, qv, saved.β, saved.meancount, idx))
            methods[m] = (lfc = lfc, q = qv)
        end
        open(io -> foreach(l -> println(io, l), rows), file, "a")
        open(io -> foreach(l -> println(io, l), crows), cfile, "a")
        serialize(path, merge(saved, (methods = methods,)))
        logline("choose_q grid $key: q = $(c.q) (bound $(c.qmax)), $(round(t; digits = 0)) s")
    end
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
