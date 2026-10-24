#figure of the benchmark by number of samples per group (benchmark/results/grid_benchmark.csv): FDP, power and RMSE,
#one column per scenario, all methods, effect sd 0.8; random groups at 2 to 10 samples per group, the scenarios with a
#target confounding from 5 per group (with 2 or 3 a target is met only on selected sets of samples)
#run: julia --project=benchmark/figures benchmark/figures/sample_size_figures.jl [output directory]
#writes sample_size.png and strong_signal.png (default: docs/src/figures)
using CairoMakie, Statistics

const SROOT = joinpath(@__DIR__, "..", "..")
const SCSV = joinpath(SROOT, "benchmark", "results", "grid_benchmark.csv")

#methods: key, label, colour, line style, marker; PLNRUV in blue, the oracles in black
const SMETHODS = [
    ("PLNRUV_EB_controls", "PLNRUV (controls)", "#2a78d6", :solid, :circle),
    ("PLNRUV_EB_sparsity", "PLNRUV (sparsity)", "#2a78d6", :dash, :circle),
    ("PLNRUV_gaussian_logratio", "PLNRUV Gaussian log-ratios", "#2a78d6", :dot, :diamond),
    ("oracle_design", "Oracle: true normalisation", "#0b0b0b", :solid, :rect),
    ("oracle_prior", "Oracle: + true effect prior", "#0b0b0b", :dash, :rect),
    ("DESeq2", "DESeq2", "#eb6834", :solid, :utriangle),
    ("DESeq2_controls", "DESeq2 (control size factors)", "#eb6834", :dash, :utriangle),
    ("edgeR", "edgeR", "#1baf7a", :solid, :dtriangle),
    ("limma_voom", "limma-voom", "#eda100", :solid, :ltriangle),
    ("cate_nc", "CATE (negative controls)", "#6250d6", :solid, :rtriangle),
    ("vicar_mouthwash", "MOUTHWASH (vicar)", "#008300", :solid, :pentagon),
    ("paper_RNARx_CATE", "Spike-in CATE + DESeq2 (2025)", "#e34948", :solid, :hexagon),
]
const SSCEN = [
    ("randomised", "Randomised"), ("balanced", "Balanced"), ("confounded", "Confounded"), ("global_increase", "Global increase"),
]

function srows(file)
    lines = readlines(file)
    h = split(lines[1], ",")
    return [Dict(zip(h, split(l, ","))) for l in lines[2:end]]
end

#mean and standard error of a metric per n per group, for one scenario, signal and method
function scurve(rows, sc, s, m, metric)
    out = NTuple{4, Float64}[]
    for g in (2, 3, 5, 10)
        x = [
            parse(Float64, r[metric]) for r in rows if r["scenario"] == sc && parse(Float64, r["signal"]) == s &&
                r["method"] == m && parse(Int, r["npergroup"]) == g
        ]
        isempty(x) && continue
        push!(out, (g, mean(x), length(x) > 1 ? std(x) / sqrt(length(x)) : 0.0, length(x)))
    end
    return out
end

#FDP, power and RMSE against the number of samples per group, one column per scenario; each scenario is shown at the
#sizes where it exists (a target confounding is met without selecting samples only from 5 per group)
const FIGCOLS = [
    ("randomised", "Random groups"), ("balanced", "Balanced"), ("confounded", "Confounded"),
    ("global_randomised", "Global increase,\nrandom groups"), ("global_increase", "Global increase,\nconfounded"),
]

function sample_size_figure(rows, s; columns = FIGCOLS, ns = [2, 3, 5, 10], title)
    fig = Figure(size = (300 * length(columns) + 120, 1020), fontsize = 14)
    metrics_ = (("fdp", "FDP"), ("power", "Power"), ("rmse", "RMSE of log₂FC"))
    for (j, (sc, slabel)) in enumerate(columns), (i, (metric, mlabel)) in enumerate(metrics_)
        ax = Axis(
            fig[i, j]; title = i == 1 ? slabel : "", xlabel = i == length(metrics_) ? "Samples per group" : "",
            ylabel = j == 1 ? mlabel : "", xscale = log2, xticks = (ns, string.(ns)), xgridcolor = (:gray, 0.15),
            ygridcolor = (:gray, 0.15), limits = ((minimum(ns) / 1.15, maximum(ns) * 1.15), (0, nothing))
        )
        metric == "fdp" && hlines!(ax, [0.05]; color = (:gray, 0.6), linestyle = :dot, linewidth = 1)
        for (key, label, col, ls, mk) in SMETHODS
            c = filter(x -> x[1] in ns, scurve(rows, sc, s, key, metric))
            isempty(c) && continue
            xs, ms, ses = first.(c), getindex.(c, 2), getindex.(c, 3)
            errorbars!(ax, xs, ms, ses; color = (col, 0.5), whiskerwidth = 4)
            scatterlines!(ax, xs, ms; color = col, linestyle = ls, linewidth = 2, marker = mk, markersize = 9)
        end
    end
    elems = [[LineElement(color = c, linestyle = ls, linewidth = 2), MarkerElement(color = c, marker = mk, markersize = 9)] for (_, _, c, ls, mk) in SMETHODS]
    Legend(fig[length(metrics_) + 1, 1:length(columns)], elems, [l for (_, l, _, _, _) in SMETHODS]; orientation = :horizontal, nbanks = length(columns) >= 4 ? 3 : 6, framevisible = false)
    Label(fig[0, 1:length(columns)], title; fontsize = 16)
    return fig
end

function sample_size_figures(out = joinpath(SROOT, "docs", "src", "figures"))
    rows = srows(SCSV)
    f = sample_size_figure(
        rows, 0.8; title = "Effect sd 0.8: FDP and power at 0.05 (s-value or adjusted p-value), RMSE of the log₂ fold changes; mean ± s.e. over the data sets; dotted line: 0.05"
    )
    save(joinpath(out, "sample_size.png"), f; px_per_unit = 2)
    #stronger effects with 2 and 3 samples per group, where the effects of sd 0.8 are mostly too small to be found
    f2 = sample_size_figure(
        rows, 1.6; columns = [FIGCOLS[1], FIGCOLS[4]], ns = [2, 3],
        title = "Effect sd 1.6, 2 and 3 samples per group: FDP, power and RMSE; mean ± s.e. over the data sets"
    )
    save(joinpath(out, "strong_signal.png"), f2; px_per_unit = 2)
    return f, f2
end
