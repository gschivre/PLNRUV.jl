#figures of the GTEx benchmark (benchmark/results/gtex_benchmark.csv)
#run: julia --project=benchmark/figures benchmark/figures/make_figures.jl
#writes docs/src/figures/benchmark.png (and .svg)
using CairoMakie, Statistics

const ROOT = joinpath(@__DIR__, "..", "..")
const CSV = joinpath(ROOT, "benchmark", "results", "gtex_benchmark.csv")
const OUT = joinpath(ROOT, "docs", "src", "figures")

#display names, in plotting order (top to bottom), and their family
const METHODS = [
    ("PLNRUV_bayes", "PLNRUV fully Bayesian", :plnruv),
    ("PLNRUV_EB_controls", "PLNRUV EB (controls)", :plnruv),
    ("PLNRUV_EB_sparsity", "PLNRUV EB (sparsity)", :plnruv),
    ("PLNRUV_gaussian_logratio", "PLNRUV Gaussian log-ratios", :plnruv),
    ("paper_RNARx_CATE", "Spike-in CATE + DESeq2 (2025)", :paper),
    ("DESeq2_controls", "DESeq2 (control size factors)", :reference),
    ("DESeq2", "DESeq2", :reference),
    ("edgeR", "edgeR", :reference),
    ("limma_voom", "limma-voom", :reference),
    ("cate_nc", "CATE (negative controls)", :reference),
    ("vicar_mouthwash", "MOUTHWASH (vicar)", :reference),
]
const SCENARIOS = [
    ("balanced", "Balanced\n20 % DE, random groups"),
    ("confounded", "Confounded\ngroups correlated with 2 hidden factors"),
    ("global_increase", "Global increase\n60 % of genes up, confounded"),
]
const METRICS = [(:fdp, "False discovery\nproportion"), (:power, "Power"), (:rmse, "RMSE of log₂FC")]

#reference palette (light): categorical slots 1-2 and the secondary text ink for the references
const COLORS = Dict(:plnruv => colorant"#2a78d6", :paper => colorant"#eb6834", :reference => colorant"#8a8984")
const INK = colorant"#0b0b0b"
const INK2 = colorant"#52514e"
const GRID = colorant"#e6e5e1"
const SURFACE = colorant"#fcfcfb"

function read_rows(file)
    lines = readlines(file)
    header = split(lines[1], ",")
    return [Dict(zip(header, split(l, ","))) for l in lines[2:end]]
end

#mean and standard error over replicates of one metric, for a scenario and a method
function stat(rows, scenario, method, metric)
    x = [parse(Float64, r[String(metric)]) for r in rows if r["scenario"] == scenario && r["method"] == method]
    isempty(x) && return (NaN, NaN, 0)
    return (mean(x), length(x) > 1 ? std(x) / sqrt(length(x)) : 0.0, length(x))
end

function benchmark_figure(rows)
    methods = [m for m in METHODS if any(r -> r["method"] == m[1], rows)]
    nm = length(methods)
    fig = Figure(; size = (1250, 260 + 230 * length(METRICS)), backgroundcolor = SURFACE, fontsize = 13)
    for (j, (sc, sctitle)) in enumerate(SCENARIOS), (i, (metric, mlabel)) in enumerate(METRICS)
        ax = Axis(
            fig[i, j];
            backgroundcolor = SURFACE, xgridcolor = GRID, ygridvisible = false,
            yticks = (1:nm, [m[2] for m in methods]), yreversed = true,
            yticklabelsvisible = j == 1, yticklabelcolor = INK, xticklabelcolor = INK2,
            leftspinevisible = false, rightspinevisible = false, topspinevisible = false, bottomspinecolor = INK2,
            xtickcolor = INK2, yticksvisible = false,
            title = i == 1 ? sctitle : "", titlecolor = INK, titlesize = 14, titlegap = 10,
            xlabel = j == 2 ? replace(mlabel, "\n" => " ") : "", xlabelcolor = INK2,
            ylabel = j == 1 ? mlabel : "", ylabelcolor = INK2,
        )
        metric === :fdp && vlines!(ax, [0.05]; color = INK2, linestyle = :dash, linewidth = 1.5)
        for (k, (key, _, fam)) in enumerate(methods)
            μ, se, n = stat(rows, sc, key, metric)
            isnan(μ) && continue
            rangebars!(ax, [k], [μ - se], [μ + se]; direction = :x, color = COLORS[fam], linewidth = 2, whiskerwidth = 6)
            scatter!(ax, [μ], [k]; color = COLORS[fam], markersize = 10, strokecolor = SURFACE, strokewidth = 2)
        end
        metric === :power && xlims!(ax, -0.02, 1.0)
        (metric === :fdp || metric === :rmse) && xlims!(ax; low = 0)
    end
    Legend(
        fig[length(METRICS) + 1, 1:length(SCENARIOS)],
        [MarkerElement(; color = COLORS[f], marker = :circle, markersize = 10) for f in (:plnruv, :paper, :reference)],
        ["PLNRUV.jl", "Spike-in CATE pipeline of Schivre et al. (2025)", "Reference methods (R)"];
        orientation = :horizontal, framevisible = false, labelcolor = INK,
    )
    Label(
        fig[0, 1:length(SCENARIOS)],
        "Binomially thinned GTEx muscle data: 2 000 genes (200 negative controls), 20 samples; mean ± 1 SE over replicates. Dashed line: nominal 5 % FDR.";
        color = INK2, fontsize = 13,
    )
    rowgap!(fig.layout, 12)
    return fig
end

rows = read_rows(CSV)
mkpath(OUT)
fig = benchmark_figure(rows)
save(joinpath(OUT, "benchmark.png"), fig; px_per_unit = 2)
save(joinpath(OUT, "benchmark.svg"), fig)
