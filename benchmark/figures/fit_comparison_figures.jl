#figures of the fit comparison (benchmark/fit_comparison.jl): ROC curves and RMSE of the log2 fold changes of the older
#PLN fit and the default fit (two variance priors), by scenario and number of samples per group
#run: julia --project=benchmark/figures benchmark/figures/fit_comparison_figures.jl [output directory]
#writes fit_comparison_roc.png and fit_comparison_rmse.png
using CairoMakie, Statistics, Serialization

const DIR = joinpath(@__DIR__, "..", "results", "fit_comparison")
const FITS = [
    ("pln", "Older fit (pln, removed)", "#eb6834", :solid),
    ("plnpca_trend", "Default fit (plnpca), trend prior", "#2a78d6", :solid),
    ("plnpca_common", "Default fit (plnpca), common prior", "#1baf7a", :dash),
    ("plnpca_trend_normal", "Default fit (plnpca), trend prior, normal likelihood", "#6250d6", :dashdot),
]
const NPG = [2, 3, 5, 10]
const SCEN = [
    ("balanced", "Balanced"), ("confounded", "Confounded"), ("global_increase", "Global increase"),
]
const FPR = collect(range(0, 1; length = 201))

#ROC of one data set: genes ranked by s-value, positives the changed genes (controls left out); TPR at the grid FPR
#(step interpolation), and the (FPR, TPR) of the discoveries at s < 0.05
function roc(d, fit)
    idx = setdiff(eachindex(d.β), d.controls)
    s, pos = fit.s[idx], d.β[idx] .!= 0
    any(isnan, s) && return nothing
    o = sortperm(s)
    tp, fp = cumsum(pos[o]) ./ max(count(pos), 1), cumsum(.!pos[o]) ./ max(count(.!pos), 1)
    tpr = [
        let k = findlast(<=(f), fp)
            k === nothing ? 0.0 : tp[k]
        end for f in FPR
    ]
    disc = s .< 0.05
    op = (count(disc .& .!pos) / max(count(.!pos), 1), count(disc .& pos) / max(count(pos), 1))
    return (tpr = tpr, op = op)
end

rmse(d, fit) = (idx = setdiff(eachindex(d.β), d.controls); sqrt(mean(abs2, fit.lfc[idx] .- d.β[idx])))

function load_cells()
    cells = Dict{Tuple{String, Int}, Vector{Any}}()
    isdir(DIR) || return cells
    for f in readdir(DIR; join = true)
        endswith(f, ".jls") || continue
        m = match(r"^(.*)_n(\d+)_r\d+\.jls$", basename(f))
        push!(get!(cells, (m[1], parse(Int, m[2])), Any[]), deserialize(f))
    end
    return cells
end

function roc_figure(cells, out)
    fig = Figure(size = (1300, 1000), fontsize = 14)
    for (i, (sname, slabel)) in enumerate(SCEN), (j, g) in enumerate(NPG)
        ds = get(cells, (sname, g), Any[])
        ax = Axis(
            fig[i, j]; title = "$slabel, $g per group (n = $(length(ds)))", xlabel = i == length(SCEN) ? "False positive rate" : "",
            ylabel = j == 1 ? "True positive rate" : "", limits = (0, 1, 0, 1), aspect = 1,
            xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15)
        )
        lines!(ax, [0, 1], [0, 1]; color = (:gray, 0.5), linewidth = 1, linestyle = :dot)
        if isempty(ds)
            text!(ax, 0.5, 0.5; text = "not reachable", align = (:center, :center), color = :gray40)
            continue
        end
        for (key, label, col, ls) in FITS
            rs = filter(!isnothing, [roc(d, d.fits[key]) for d in ds if haskey(d.fits, key)])
            isempty(rs) && continue
            lines!(ax, FPR, vec(mean(reduce(hcat, [r.tpr for r in rs]); dims = 2)); color = col, linewidth = 2, linestyle = ls, label = label)
            scatter!(ax, [mean(r.op[1] for r in rs)], [mean(r.op[2] for r in rs)]; color = col, markersize = 10, strokecolor = :white, strokewidth = 2)
        end
    end
    Legend(
        fig[length(SCEN) + 1, 1:length(NPG)], [LineElement(color = c, linestyle = ls, linewidth = 2) for (_, _, c, ls) in FITS],
        [l for (_, l, _, _) in FITS]; orientation = :horizontal, framevisible = false
    )
    Label(fig[0, 1:length(NPG)], "ROC curves (mean over data sets); dot: discoveries at s-value < 0.05"; fontsize = 16)
    save(joinpath(out, "fit_comparison_roc.png"), fig; px_per_unit = 2)
    return fig
end

function rmse_figure(cells, out)
    fig = Figure(size = (1300, 420), fontsize = 14)
    for (i, (sname, slabel)) in enumerate(SCEN)
        ax = Axis(
            fig[1, i]; title = slabel, xlabel = "Samples per group", ylabel = i == 1 ? "RMSE of log2 fold changes" : "",
            xticks = NPG, xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15)
        )
        for (k, (key, label, col, ls)) in enumerate(FITS)
            xs, ms, ses = Float64[], Float64[], Float64[]
            for g in NPG
                v = [rmse(d, d.fits[key]) for d in get(cells, (sname, g), Any[]) if haskey(d.fits, key)]
                v = filter(isfinite, v)
                isempty(v) && continue
                push!(xs, g * (1 + 0.03 * (k - (length(FITS) + 1) / 2))) #small offsets so the error bars do not overlap
                push!(ms, mean(v))
                push!(ses, length(v) > 1 ? std(v) / sqrt(length(v)) : 0.0)
            end
            isempty(xs) && continue
            errorbars!(ax, xs, ms, ses; color = col, whiskerwidth = 6)
            scatterlines!(ax, xs, ms; color = col, linestyle = ls, linewidth = 2, markersize = 9, label = label)
        end
        ax.xscale = log2
    end
    Legend(
        fig[2, 1:length(SCEN)], [LineElement(color = c, linestyle = ls, linewidth = 2) for (_, _, c, ls) in FITS],
        [l for (_, l, _, _) in FITS]; orientation = :horizontal, framevisible = false
    )
    Label(fig[0, 1:length(SCEN)], "RMSE of the log2 fold changes (mean ± s.e. over data sets; controls left out)"; fontsize = 16)
    save(joinpath(out, "fit_comparison_rmse.png"), fig; px_per_unit = 2)
    return fig
end

const NOMINAL = collect(range(0.0025, 0.2; length = 80))

#achieved FDP of the discoveries at s-value < level (null genes and wrong signs are false; 0 without discoveries, so the
#mean over data sets estimates the FDR)
function achieved_fdp(d, fit)
    idx = setdiff(eachindex(d.β), d.controls)
    s, β, lfc = fit.s[idx], d.β[idx], fit.lfc[idx]
    any(isnan, s) && return nothing
    wrong = (β .== 0) .| (sign.(β) .!= sign.(lfc))
    return [
        let disc = s .< l
            count(disc) == 0 ? 0.0 : count(disc .& wrong) / count(disc)
        end for l in NOMINAL
    ]
end

function fdp_figure(cells, out)
    fig = Figure(size = (1300, 1000), fontsize = 14)
    for (i, (sname, slabel)) in enumerate(SCEN), (j, g) in enumerate(NPG)
        ds = get(cells, (sname, g), Any[])
        curves = Dict(key => filter(!isnothing, [achieved_fdp(d, d.fits[key]) for d in ds if haskey(d.fits, key)]) for (key, _, _, _) in FITS)
        top = maximum((maximum(vec(mean(reduce(hcat, c); dims = 2))) for c in values(curves) if !isempty(c)); init = 0.0)
        ax = Axis(
            fig[i, j]; title = "$slabel, $g per group (n = $(length(ds)))",
            xlabel = i == length(SCEN) ? "Nominal level (s-value cutoff)" : "", ylabel = j == 1 ? "Achieved FDP" : "",
            limits = (0, 0.2, 0, max(0.25, 1.05 * top)), xticks = 0:0.05:0.2,
            xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15)
        )
        lines!(ax, [0, 0.2], [0, 0.2]; color = (:gray, 0.6), linewidth = 1, linestyle = :dot)
        vlines!(ax, [0.05]; color = (:gray, 0.4), linewidth = 1)
        for (key, label, col, ls) in FITS
            isempty(curves[key]) && continue
            lines!(ax, NOMINAL, vec(mean(reduce(hcat, curves[key]); dims = 2)); color = col, linewidth = 2, linestyle = ls, label = label)
        end
    end
    Legend(
        fig[length(SCEN) + 1, 1:length(NPG)], [LineElement(color = c, linestyle = ls, linewidth = 2) for (_, _, c, ls) in FITS],
        [l for (_, l, _, _) in FITS]; orientation = :horizontal, framevisible = false
    )
    Label(fig[0, 1:length(NPG)], "Achieved FDP against the nominal level (mean over data sets; dotted: target, grey line: 0.05)"; fontsize = 16)
    save(joinpath(out, "fit_comparison_fdp.png"), fig; px_per_unit = 2)
    return fig
end

function fit_comparison_figures(out = joinpath(@__DIR__, "..", "results"))
    cells = load_cells()
    return roc_figure(cells, out), rmse_figure(cells, out), fdp_figure(cells, out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_comparison_figures(isempty(ARGS) ? joinpath(@__DIR__, "..", "results") : ARGS[1])
end
