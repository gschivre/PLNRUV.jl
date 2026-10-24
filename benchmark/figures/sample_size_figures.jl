#figure of the benchmark by number of samples per group (benchmark/results/grid_benchmark.csv): FDP and power, one
#column per scenario, all methods; by default the scenarios whose confounding is set exactly (balanced, confounded,
#global increase) at 5 and 10 samples per group: with 2 or 3 per group a target confounding is reached only on
#selected sets of samples
#run: julia --project=benchmark/figures benchmark/figures/sample_size_figures.jl [output directory]
#writes sample_size.png (default: docs/src/figures)
using CairoMakie, Statistics

const SROOT = joinpath(@__DIR__, "..", "..")
const SCSV = joinpath(SROOT, "benchmark", "results", "grid_benchmark.csv")

#methods: key, label, colour, line style, marker; PLNRUV in blue, the oracles in black
const SMETHODS = [
    ("PLNRUV_EB_controls", "PLNRUV (controls)", "#2a78d6", :solid, :circle),
    ("PLNRUV_EB_sparsity", "PLNRUV (sparsity)", "#2a78d6", :dash, :circle),
    ("PLNRUV_gaussian_logratio", "PLNRUV Gaussian log-ratios", "#2a78d6", :dot, :diamond),
    ("PLNRUV_EB_controls_chooseq", "PLNRUV (controls), q by choose_q", "#2a78d6", nothing, :star5),
    ("oracle_design", "Oracle: known factors + normalisation", "#0b0b0b", :solid, :rect),
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

function sample_size_figure(rows, s, ns; title, scenarios = SSCEN)
    fig = Figure(size = (350 * length(scenarios) + 100, 760), fontsize = 14)
    for (j, (sc, slabel)) in enumerate(scenarios), (i, (metric, mlabel)) in enumerate((("fdp", "FDP"), ("power", "Power")))
        ax = Axis(
            fig[i, j]; title = i == 1 ? slabel : "", xlabel = i == 2 ? "Samples per group" : "", ylabel = j == 1 ? mlabel : "",
            xscale = log2, xticks = (ns, string.(ns)), xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15),
            limits = ((minimum(ns) / 1.15, maximum(ns) * 1.15), metric == "fdp" ? (0, nothing) : (0, nothing))
        )
        metric == "fdp" && hlines!(ax, [0.05]; color = (:gray, 0.6), linestyle = :dot, linewidth = 1)
        for (key, label, col, ls, mk) in SMETHODS
            c = filter(x -> x[1] in ns, scurve(rows, sc, s, key, metric))
            isempty(c) && continue
            xs, ms, ses = first.(c), getindex.(c, 2), getindex.(c, 3)
            errorbars!(ax, xs, ms, ses; color = (col, 0.5), whiskerwidth = 4)
            if ls === nothing
                scatter!(ax, xs, ms; color = col, marker = mk, markersize = 14, strokecolor = :white, strokewidth = 1)
            else
                scatterlines!(ax, xs, ms; color = col, linestyle = ls, linewidth = 2, marker = mk, markersize = 9)
            end
        end
    end
    elems = [
        ls === nothing ? MarkerElement(color = c, marker = mk, markersize = 14) :
            [LineElement(color = c, linestyle = ls, linewidth = 2), MarkerElement(color = c, marker = mk, markersize = 9)]
            for (_, _, c, ls, mk) in SMETHODS
    ]
    Legend(fig[3, 1:length(scenarios)], elems, [l for (_, l, _, _, _) in SMETHODS]; orientation = :horizontal, nbanks = 5, framevisible = false)
    Label(fig[0, 1:length(scenarios)], title; fontsize = 16)
    return fig
end

function sample_size_figures(out = joinpath(SROOT, "docs", "src", "figures"))
    rows = srows(SCSV)
    f = sample_size_figure(
        rows, 0.8, [5, 10]; scenarios = SSCEN[2:end],
        title = "FDP and power at 0.05 (s-value or adjusted p-value), mean ± s.e. over 10 data sets; dotted line: 0.05"
    )
    save(joinpath(out, "sample_size.png"), f; px_per_unit = 2)
    return f
end

if abspath(PROGRAM_FILE) == @__FILE__
    sample_size_figures(isempty(ARGS) ? joinpath(SROOT, "docs", "src", "figures") : ARGS[1])
end
