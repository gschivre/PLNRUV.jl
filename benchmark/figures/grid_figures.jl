#figure of the replicates x signal grid (benchmark/results/grid_benchmark.csv): absolute FDP and power against the
#number of replicates per group, one panel per scenario and signal strength; the main methods in colour, the others
#as thin grey lines; cells with fewer than minrep replicates are left out (partial runs)
#run: julia --project=benchmark/figures benchmark/figures/grid_figures.jl
#writes docs/src/figures/grid.png
using CairoMakie, Statistics

const GROOT = joinpath(@__DIR__, "..", "..")
const GCSV = joinpath(GROOT, "benchmark", "results", "grid_benchmark.csv")

#highlighted methods: key, label, colour, line style
const HIGHLIGHT = [
    ("oracle_design", "Oracle: known confounders + normalisation", colorant"#0b0b0b", :solid),
    ("oracle_prior", "Oracle: + true effect prior", colorant"#0b0b0b", :dash),
    ("PLNRUV_EB_controls", "PLNRUV EB (controls)", colorant"#2a78d6", :solid),
    ("PLNRUV_EB_sparsity", "PLNRUV EB (sparsity)", colorant"#2a78d6", :dash),
    ("PLNRUV_gaussian_logratio", "PLNRUV Gaussian log-ratios", colorant"#2a78d6", :dot),
    ("PLNRUV_PLN_EB_controls", "PLNRUV EB, older fit (controls)", colorant"#8fb8e8", :solid),
    ("paper_RNARx_CATE", "Spike-in CATE + DESeq2 (2025)", colorant"#eb6834", :solid),
    ("DESeq2", "DESeq2", colorant"#52514e", :solid),
]
const GSCEN = [
    ("balanced", "Balanced"), ("confounded", "Confounded"), ("global_increase", "Global increase"),
]
const GGREY = colorant"#c9c8c3"
const GINK2 = colorant"#52514e"
const GSURF = colorant"#fcfcfb"

function grid_rows(file)
    lines = readlines(file)
    h = split(lines[1], ",")
    return [Dict(zip(h, split(l, ","))) for l in lines[2:end]]
end

#mean and SE of a metric per n per group, for one scenario, signal and method (cells with at least minrep reps)
function curve(rows, sc, s, m, metric; minrep = 10)
    ns = sort(unique(parse(Int, r["npergroup"]) for r in rows))
    out = NTuple{3, Float64}[]
    for g in ns
        x = [
            parse(Float64, r[metric]) for r in rows if r["scenario"] == sc && parse(Float64, r["signal"]) == s &&
                r["method"] == m && parse(Int, r["npergroup"]) == g
        ]
        length(x) >= minrep && push!(out, (g, mean(x), std(x) / sqrt(length(x))))
    end
    return out
end

function grid_figure(rows; minrep = 10)
    signals = sort(unique(parse(Float64, r["signal"]) for r in rows))
    methods = unique(r["method"] for r in rows)
    others = setdiff(methods, first.(HIGHLIGHT))
    ns = sort(unique(parse(Int, r["npergroup"]) for r in rows))
    fig = Figure(; size = (1250, 1500), backgroundcolor = GSURF, fontsize = 13)
    for (b, (metric, mlabel)) in enumerate((("fdp", "FDP"), ("power", "Power")))
        for (i, (sc, sctitle)) in enumerate(GSCEN), (j, s) in enumerate(signals)
            ax = Axis(
                fig[(b - 1) * 3 + i, j]; backgroundcolor = GSURF, xscale = log10, xticks = (ns, string.(ns)),
                title = b == 1 && i == 1 ? "effect sd = $s" : "", ylabel = j == 1 ? "$sctitle\n$mlabel" : "",
                xlabel = b == 2 && i == 3 ? "replicates per group" : "", ylabelcolor = GINK2, xlabelcolor = GINK2,
                topspinevisible = false, rightspinevisible = false, xgridvisible = false,
            )
            metric == "fdp" && hlines!(ax, [0.05]; color = GINK2, linestyle = :dash, linewidth = 1)
            for m in others
                c = curve(rows, sc, s, m, metric; minrep = minrep)
                isempty(c) || lines!(ax, first.(c), [x[2] for x in c]; color = GGREY, linewidth = 1)
            end
            for (m, _, col, ls) in HIGHLIGHT
                c = curve(rows, sc, s, m, metric; minrep = minrep)
                isempty(c) && continue
                g, μ, se = first.(c), [x[2] for x in c], [x[3] for x in c]
                rangebars!(ax, g, μ .- se, μ .+ se; color = col, whiskerwidth = 4)
                scatterlines!(ax, g, μ; color = col, linestyle = ls, markersize = 7, linewidth = 2)
            end
            ylims!(ax, 0, metric == "fdp" ? 1 : 1)
        end
    end
    Legend(
        fig[7, 1:3],
        [[LineElement(; color = col, linestyle = ls, linewidth = 2) for (_, _, col, ls) in HIGHLIGHT]; LineElement(; color = GGREY)],
        [[h[2] for h in HIGHLIGHT]; "other methods (edgeR, limma-voom, CATE, vicar, DESeq2 control size factors)"];
        orientation = :horizontal, nbanks = 3, framevisible = false,
    )
    Label(
        fig[0, 1:3],
        "Binomially thinned GTEx muscle data, 2 000 genes (200 controls); mean ± 1 SE over $minrep replicates. Rows 1-3: false discovery proportion (dashed: 5 %); rows 4-6: power.";
        color = GINK2,
    )
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    save(joinpath(GROOT, "docs", "src", "figures", "grid.png"), grid_figure(grid_rows(GCSV)); px_per_unit = 2)
end
