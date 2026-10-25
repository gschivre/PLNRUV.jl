#FDP and power as the number of replicates and the signal strength vary, on binomially thinned GTEx data, with two
#oracles: known confounders and normalisation (oracle_design), and in addition the true effect prior (oracle_prior)
#run from the benchmark environment: julia -t auto --project=benchmark benchmark/grid_benchmark.jl <gtex counts file>
#runs the 20-sample benchmark (exact and misspecified controls) and the grid; for each, in benchmark/results/:
#<name>.csv (summary at level 0.05), <name>_curves.csv (calibration and ranked curves by expression tercile) and
#<name>_datasets/ (one .jls per simulated data set: truth and every method's output)
include(joinpath(@__DIR__, "gtex_benchmark.jl"))
using LinearAlgebra, Serialization

#nominal levels of the calibration / FDP-power curves, and observed FDP targets of the ranked curves
const LEVELS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]
const FDP_TARGETS = [0.01, 0.05, 0.1]

#genes of idx by tercile of their mean count (low, mid, high), and all of them
function strata(meancount, idx)
    t = quantile(meancount[idx], [1 / 3, 2 / 3])
    return [
        "all" => idx, "low" => idx[meancount[idx] .<= t[1]],
        "mid" => idx[t[1] .< meancount[idx] .<= t[2]], "high" => idx[meancount[idx] .> t[2]],
    ]
end

#curve rows of one method: at each nominal level (kind "nominal": genes with q < level) and at each observed FDP
#target when genes are ranked by q (kind "ranked": the largest prefix whose FDP is at most the target);
#a discovery is false when the gene is null or its sign is wrong
function curve_rows(prefix, lfc, qv, β, meancount, idx)
    rows = String[]
    for (name, genes) in strata(meancount, idx)
        nsig = max(count(!=(0), β[genes]), 1)
        wrong = (β[genes] .== 0) .| (sign.(β[genes]) .!= sign.(lfc[genes]))
        for level in LEVELS
            d = qv[genes] .< level
            nd = count(d)
            fdp = nd == 0 ? 0.0 : count(d .& wrong) / nd
            push!(rows, join((prefix, name, "nominal", level, nd, fdp, count(d .& .!wrong) / nsig), ","))
        end
        order = sortperm(qv[genes])
        fdps = cumsum(wrong[order]) ./ (1:length(order))
        for t in FDP_TARGETS
            k = something(findlast(<=(t), fdps), 0)
            push!(rows, join((prefix, name, "ranked", t, k, k == 0 ? 0.0 : fdps[k], count(.!wrong[order[1:k]]) / nsig), ","))
        end
    end
    return rows
end

#set STOP[] = true to stop a background run after the current data set
const STOP = Ref(false)

#progress line in results/run.log (background runs have no eval to report to)
function logline(msg)
    open(io -> println(io, Libc.strftime("%Y-%m-%d %H:%M:%S", time()), " ", msg), joinpath(@__DIR__, "results", "run.log"), "a")
    return nothing
end

const CURVE_HEADER = "stratum,kind,level,ndisc,fdp,power"

#everything a figure might need later, as plain types (no package structs, so the file stays readable when
#PLNRUV changes): the simulated data set itself (counts, design, factors, truth), so that later methods run on
#exactly the data the saved methods saw, and each method's log2 fold changes and adjusted p- or s-values
function save_dataset(dir, key, sim, res)
    mkpath(dir)
    data = (
        julia = string(VERSION), β = sim.β, controls = sim.controls, null = sim.null, signal = sim.signal,
        meancount = vec(mean(sim.counts; dims = 1)), counts = sim.counts, X = sim.X, factors = sim.factors,
        genes = sim.genes, samples = sim.samples, realized_cor = sim.realized_cor, tries = sim.tries,
        redraws = sim.redraws, methods = Dict(m => (lfc = l, q = v) for (m, (l, v)) in res),
    )
    serialize(joinpath(dir, key * ".jls"), data)
    return nothing
end

#the simulated data set of a scenario, or nothing when its target confounding is out of reach (too few samples,
#e.g. no balanced 2 + 2 design: the three possible splits span the residual space of the intercept)
function simulate(seed, M, sc; ngenes, nsamp, nctrl)
    return try
        thin_2group(
            Xoshiro(seed), M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = sc.prop_null,
            signal = sc.signal, ctrl_signal = get(sc, :ctrl, Dirac(0.0)), target_cor = sc.target_cor
        )
    catch err
        err isa ErrorException && occursin("no group assignment", err.msg) || rethrow()
        nothing
    end
end

#one line of <name>_design.csv: realised correlations with the first two surrogate variables, assignment draws
#(tries) and sample redraws of the data set, or NA when the target was out of reach
function design_row(file, prefix, sim)
    dfile = replace(file, ".csv" => "_design.csv")
    isfile(dfile) || write(
        dfile, (count(==(','), prefix) == 1 ? "scenario,rep" : "scenario,npergroup,signal,rep") *
            ",rho1,rho2,tries,redraws\n"
    )
    ρ = sim === nothing ? ["NA", "NA"] : [length(sim.realized_cor) >= l ? string(round(sim.realized_cor[l]; digits = 3)) : "NA" for l in 1:2]
    tail = sim === nothing ? "NA,NA" : "$(sim.tries),$(sim.redraws)"
    open(io -> println(io, join((prefix, ρ..., tail), ",")), dfile, "a")
    return nothing
end

const NPERGROUP = [2, 3, 5, 10, 15]
const SIGNAL = [0.4, 0.8, 1.6]

#scenarios at signal strength s: the effect sd (and the mean shift of the global scenario) scale with s, s = 0.8 being
#the scenarios of gtex_benchmark.jl
grid_scenarios(s) = [
    (name = "randomised", prop_null = 0.8, signal = Normal(0.0, s), target_cor = Float64[]),
    (name = "balanced", prop_null = 0.8, signal = Normal(0.0, s), target_cor = [0.0, 0.0]),
    (name = "confounded", prop_null = 0.8, signal = Normal(0.0, s), target_cor = [0.6, 0.2]),
    (name = "global_increase", prop_null = 0.4, signal = Normal(1.25s, 0.625s), target_cor = [0.6, 0.2]),
    (name = "global_randomised", prop_null = 0.4, signal = Normal(1.25s, 0.625s), target_cor = Float64[]),
]

#PLNRUV and the oracles use the number of factors chosen by choose_q (one seed per replicate); the R methods use the
#fixed nfactors(n)
pln_q(counts, X, r) = choose_q(counts, X; qs = 0:5, rng = Xoshiro(r)).q

is_julia_method(m) = startswith(m, "PLNRUV_") || startswith(m, "oracle_")

#latent factors that the residual degrees of freedom allow: the compositional offset takes one of the n - 2
nfactors(n) = clamp(n - 4, 0, 2)

#each method separately, so that one method failing on a small design (e.g. cate or vicar with q = 0) only loses
#its own rows
function safely(f, args...)
    return try
        f(args...)
    catch err
        @warn "method failed" f err
        Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    end
end

#DESeq2 median-of-ratios size factors on the truly null genes (untouched by the thinning): the oracle normalisation
oracle_offsets(counts, nulls) = log.(size_factors(permutedims(counts), nulls))

#oracle_design: PLNRUV with the true normalisation (DESeq2 size factors on the truly null genes, no compositional
#offset to estimate) and its own latent factors, q as PLNRUV (choose_q); the surrogate variables of the scenarios are
#not given as covariates: they are principal components of the same samples, and with few samples regressing on them
#leaves the quietest directions to estimate the variances (anti-conservative); oracle_prior: the same estimates with
#the true effect prior in place of the fitted ash prior
function oracle_methods(sim, sc, q)
    #only the null genes: the controls are untouched only when they are exact
    fit = pln_mouthwash(
        sim.counts, float.(sim.X), q; controls = sim.controls, offsets = oracle_offsets(sim.counts, sim.null),
        offset_factor = false
    )
    mw = fit.mouthwash
    out = Dict("oracle_design" => untested(mw.mean ./ log(2), mw.svalue))
    #true prior in natural log units; the controls are part of the null spike
    P = size(sim.counts, 2)
    #the PLN fit covers the tested genes only (keep); the computation below is on them, then expanded
    keep = fit.tested
    p = length(keep)
    nsig = length(sim.signal)
    d = sc.signal
    π = [1 - nsig / P, nsig / P]
    σ2 = [0.0, (std(d) * log(2))^2]
    means = [0.0, mean(d) * log(2)]
    #the corrected estimates and their noise as in mouthwash_fit (default t likelihood): control-based α, its
    #uncertainty added to the variances, then the fitted inflation
    #loading uncertainty: its error variance α'Cov(Γⱼ)α, and E[ΓΓ'] in the normal equations of α
    Γ, s2, ctrl, Γcov = fit.pln.Γ, mw.s[keep] .^ 2, findall(in(Set(sim.controls)), keep), fit.pln.Γcov
    a = zeros(p)
    l = size(Γ, 1) > 0 ? [dot(mw.α, Γcov[j] * mw.α) for j in 1:p] : zeros(p)
    if size(Γ, 1) > 0 #factors switched off by the fit carry no loadings
        w = 1 ./ (s2[ctrl] .+ l[ctrl])
        G = Symmetric(Γ[:, ctrl] * (Γ[:, ctrl]' .* w) .+ sum(w[i] .* Γcov[j] for (i, j) in enumerate(ctrl)))
        free = setdiff(1:p, ctrl)
        a[free] .= vec(sum(Γ[:, free] .* (G \ Γ[:, free]); dims = 1))
    end
    r = size(Γ, 1) > 0 ? mw.β̃[keep] - Γ' * mw.α : mw.β̃[keep]
    N = PLNRUV.scale_noise(PLNRUV.NoiseModel(s2 .+ l .+ a, fit.pln.df), mw.λ)
    post = PLNRUV.ash_posterior(r, N, σ2, π; means = means)
    lfc, sv = zeros(P), ones(P)
    lfc[keep] .= post.mean ./ log(2)
    sv[keep] .= PLNRUV.cumulative_rate(post.lfsr)
    out["oracle_prior"] = (lfc, sv)
    return out
end

#set R_ONLY[] = true to run only the R methods (the PLNRUV and oracle rows are added later on the saved data sets)
const R_ONLY = Ref(false)

function grid_methods(sim, sc, q, qpln)
    c, X, ctrl = sim.counts, sim.X, sim.controls
    out = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    merge!(out, safely(r_methods, c, X, ctrl, q))
    merge!(out, safely(paper_method, c, X, ctrl, q))
    R_ONLY[] && return out
    merge!(out, safely(julia_methods, c, X, ctrl, qpln))
    merge!(out, safely(oracle_methods, sim, sc, qpln))
    return out
end

function run_grid(
        M; nrep = 10, ngenes = 2_000, nctrl = 200, seed = 1, outfile = nothing, done = Set{Tuple}(),
        npergroup = NPERGROUP, signals = SIGNAL, scenarios = nothing
    )
    header = "scenario,npergroup,signal,rep,q,method,ndisc,fdp,power,rmse"
    if outfile !== nothing && !isfile(outfile)
        write(outfile, header * "\n")
    end
    for g in npergroup, s in signals, sc in grid_scenarios(s), r in 1:nrep
        scenarios === nothing || sc.name in scenarios || continue
        (sc.name, g, s, r) in done && continue
        STOP[] && return nothing
        #the same seeds across signal strengths and scenarios: the same samples and genes, only the effects differ
        sim = simulate(seed + 100g + r, M, sc; ngenes = ngenes, nsamp = 2g, nctrl = nctrl)
        outfile === nothing || design_row(outfile, join((sc.name, g, s, r), ","), sim)
        if sim === nothing
            logline("grid $(sc.name) n/group = $g signal = $s rep $r: NA (target out of reach)")
            continue
        end
        q = nfactors(2g)
        qpln = R_ONLY[] ? q : pln_q(sim.counts, sim.X, r)
        idx = setdiff(1:ngenes, sim.controls)
        res = grid_methods(sim, sc, q, qpln)
        meancount = vec(mean(sim.counts; dims = 1))
        rows, crows = String[], String[]
        for (method, (lfc, qv)) in res
            x = metrics(lfc, qv, sim.β, idx)
            push!(rows, join((sc.name, g, s, r, is_julia_method(method) ? qpln : q, method, x.ndisc, x.fdp, x.power, x.rmse), ","))
            append!(crows, curve_rows(join((sc.name, g, s, r, method), ","), lfc, qv, sim.β, meancount, idx))
        end
        #appended after each data set, so an interrupted run can resume (pass the finished cells as done)
        if outfile !== nothing
            open(io -> foreach(l -> println(io, l), rows), outfile, "a")
            cfile = replace(outfile, ".csv" => "_curves.csv")
            isfile(cfile) || write(cfile, "scenario,npergroup,signal,rep,method," * CURVE_HEADER * "\n")
            open(io -> foreach(l -> println(io, l), crows), cfile, "a")
            save_dataset(replace(outfile, ".csv" => "_datasets"), "$(sc.name)_n$(g)_s$(s)_r$(r)", sim, res)
        end
        logline("grid $(sc.name) n/group = $g signal = $s rep $r")
    end
    return nothing
end

#misspecified negative controls, 20 samples: noisy (each control has a small effect of its own) or biased (the
#controls shift together with the group, as a spike-in whose amount follows the condition)
const CTRL_SCENARIOS = [
    (name = "confounded_noisy_ctrl", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = [0.6, 0.2], ctrl = Normal(0.0, 0.2)),
    (name = "confounded_biased_ctrl", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = [0.6, 0.2], ctrl = Normal(0.15, 0.05)),
    (name = "global_noisy_ctrl", prop_null = 0.4, signal = Normal(1.0, 0.5), target_cor = [0.6, 0.2], ctrl = Normal(0.0, 0.2)),
    (name = "global_biased_ctrl", prop_null = 0.4, signal = Normal(1.0, 0.5), target_cor = [0.6, 0.2], ctrl = Normal(0.15, 0.05)),
]

#the scenarios of gtex_benchmark.jl, with exact controls
const MAIN_SCENARIOS = [merge(sc, (ctrl = Dirac(0.0),)) for sc in SCENARIOS]

#20-sample benchmark on a list of scenarios (MAIN_SCENARIOS: gtex_benchmark.csv; CTRL_SCENARIOS: gtex_controls.csv),
#all methods; same seeds as gtex_benchmark.jl, so a misspecified-control data set
#differs from its exact-control twin only by the control effects; rows in the format of gtex_benchmark.csv, appended
#after each data set (an interrupted run resumes)
function run_scenarios(M, scenarios, file; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1)
    isfile(file) || write(file, "scenario,rep,method,ndisc,fdp,power,rmse\n")
    done = Set((String(f[1]), parse(Int, f[2])) for f in (split(l, ",") for l in Iterators.drop(eachline(file), 1)))
    for sc in scenarios, r in 1:nrep
        (sc.name, r) in done && continue
        STOP[] && return nothing
        sim = simulate(seed + r, M, sc; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl)
        design_row(file, join((sc.name, r), ","), sim)
        if sim === nothing
            logline("$(basename(file)) $(sc.name) rep $r: NA (target out of reach)")
            continue
        end
        idx = setdiff(1:ngenes, sim.controls)
        res = grid_methods(sim, sc, q, R_ONLY[] ? q : pln_q(sim.counts, sim.X, r))
        meancount = vec(mean(sim.counts; dims = 1))
        rows, crows = String[], String[]
        for (method, (lfc, qv)) in res
            x = metrics(lfc, qv, sim.β, idx)
            push!(rows, join((sc.name, r, method, x.ndisc, x.fdp, x.power, x.rmse), ","))
            append!(crows, curve_rows(join((sc.name, r, method), ","), lfc, qv, sim.β, meancount, idx))
        end
        open(io -> foreach(l -> println(io, l), rows), file, "a")
        cfile = replace(file, ".csv" => "_curves.csv")
        isfile(cfile) || write(cfile, "scenario,rep,method," * CURVE_HEADER * "\n")
        open(io -> foreach(l -> println(io, l), crows), cfile, "a")
        save_dataset(replace(file, ".csv" => "_datasets"), "$(sc.name)_r$(r)", sim, res)
        logline("$(basename(file)) $(sc.name) rep $r")
    end
    return nothing
end

#all the GTEx benchmarks in sequence as a background task of the session (the eval returns at once); each step
#resumes from its files, progress and errors go to results/run.log, STOP[] = true stops after the current data set
function launch(M)
    dir = joinpath(@__DIR__, "results")
    STOP[] = false
    return @async try
        logline("start")
        run_scenarios(M, MAIN_SCENARIOS, joinpath(dir, "gtex_benchmark.csv"))
        run_scenarios(M, CTRL_SCENARIOS, joinpath(dir, "gtex_controls.csv"))
        grid = joinpath(dir, "grid_benchmark.csv")
        run_grid(M; outfile = grid, done = finished_cells(grid))
        logline(STOP[] ? "stopped" : "done")
    catch err
        logline("ERROR " * sprint(showerror, err, catch_backtrace())[1:min(end, 3000)])
    end
end

#the (scenario, n per group, signal, rep) cells already in a results file
function finished_cells(file)
    isfile(file) || return Set{Tuple}()
    cells = Set{Tuple}()
    for l in Iterators.drop(eachline(file), 1)
        f = split(l, ",")
        push!(cells, (String(f[1]), parse(Int, f[2]), parse(Float64, f[3]), parse(Int, f[4])))
    end
    return cells
end

if abspath(PROGRAM_FILE) == @__FILE__
    gtex = read_gct(ARGS[1])
    M = gtex.counts[filter_genes(gtex.counts; min_mean = 10), :]
    file = joinpath(@__DIR__, "results", "grid_benchmark.csv")
    mkpath(dirname(file))
    run_scenarios(M, MAIN_SCENARIOS, joinpath(@__DIR__, "results", "gtex_benchmark.csv"))
    run_scenarios(M, CTRL_SCENARIOS, joinpath(@__DIR__, "results", "gtex_controls.csv"))
    run_grid(M; outfile = file, done = finished_cells(file))
end
