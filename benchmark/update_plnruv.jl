#regenerate the PLNRUV rows of benchmark/results/gtex_benchmark.csv after a change of the package, keeping the rows of
#the R methods (same simulated data: same seeds and scenarios as gtex_benchmark.jl); no RCall needed
#run: julia -t auto --project=. benchmark/update_plnruv.jl <gtex counts file> [bayes = false]
using PLNRUV, Random, Distributions, Statistics

const SCENARIOS = [
    (name = "balanced", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = Float64[]),
    (name = "confounded", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = [0.6, 0.2]),
    (name = "global_increase", prop_null = 0.4, signal = Normal(1.0, 0.5), target_cor = [0.6, 0.2]),
]

#false discovery proportion (null or wrong sign) and power at level 0.05, RMSE of the log2 fold changes
function metrics(lfc, q, β, idx; level = 0.05)
    disc = idx[q[idx] .< level]
    signal = idx[β[idx] .!= 0]
    wrong = (β[disc] .== 0) .| (sign.(β[disc]) .!= sign.(lfc[disc]))
    return (
        ndisc = length(disc), fdp = isempty(disc) ? 0.0 : mean(wrong),
        power = sum(.!wrong) / max(length(signal), 1), rmse = sqrt(mean(abs2, lfc[idx] - β[idx])),
    )
end

#the PLNRUV methods of gtex_benchmark.jl (log2 fold changes and s-values)
function plnruv_methods(counts, X, controls, q; bayes = true)
    out = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    p = size(counts, 2)
    eb = pln_mouthwash(counts, X, q; controls = controls).mouthwash
    out["PLNRUV_EB_controls"] = (eb.mean ./ log(2), eb.svalue)
    ebf = pln_mouthwash(counts, X, q).mouthwash
    out["PLNRUV_EB_sparsity"] = (ebf.mean ./ log(2), ebf.svalue)
    g = mouthwash(log_ratio(counts, controls), X, q; controls = controls)
    out["PLNRUV_gaussian_logratio"] = (g.mean ./ log(2), g.svalue)
    if bayes
        s = summarize(pln_ruvb(counts, X, controls, q; burnin = 1_000, nsample = 250, thin = 4))
        m, sv = zeros(p), ones(p)
        est = setdiff(1:p, controls)
        m[est] .= s.ash_mean[1, est] ./ log(2)
        sv[est] .= s.ash_svalue[1, est]
        out["PLNRUV_bayes"] = (m, sv)
    end
    return out
end

function update_plnruv(M, file; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1, bayes = false)
    lines = readlines(file)
    header, old = lines[1], lines[2:end]
    redo = bayes ? ("PLNRUV_",) : ("PLNRUV_EB", "PLNRUV_gaussian")
    kept = filter(l -> !any(startswith(split(l, ",")[3], m) for m in redo), old)
    rows = String[]
    for sc in SCENARIOS, r in 1:nrep
        sim = thin_2group(
            Xoshiro(seed + r), M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = sc.prop_null,
            signal = sc.signal, target_cor = sc.target_cor
        )
        idx = setdiff(1:ngenes, sim.controls)
        for (method, (lfc, qv)) in plnruv_methods(sim.counts, sim.X, sim.controls, q; bayes = bayes)
            x = metrics(lfc, qv, sim.β, idx)
            push!(rows, join((sc.name, r, method, x.ndisc, x.fdp, x.power, x.rmse), ","))
        end
        println(sc.name, " ", r)
    end
    write(file, join([header; kept; rows], "\n") * "\n")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    gtex = read_gct(ARGS[1])
    M = gtex.counts[filter_genes(gtex.counts; min_mean = 10), :]
    bayes = length(ARGS) > 1 ? parse(Bool, ARGS[2]) : false
    update_plnruv(M, joinpath(@__DIR__, "results", "gtex_benchmark.csv"); bayes = bayes)
end
