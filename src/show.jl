#display of the fits in the REPL: log fold changes are shown on the log2 scale

#format a probability, with a floor for very small values
fmt_p(x) = isnan(x) ? "     —" : x < 1.0e-4 ? " <1e-4" : lpad(string(round(x; digits = 4)), 6)
fmt_x(x, d = 2) = isnan(x) ? "—" : string(round(x; digits = d))

#table of the top genes: one row per gene, columns given as (header, width, formatter)
function print_table(io, header, rows)
    println(io, "  ", join((lpad(h, w) for (h, w) in header), "  "))
    for r in rows
        println(io, "  ", join((lpad(c, w) for (c, (_, w)) in zip(r, header)), "  "))
    end
    return nothing
end

Base.show(io::IO, fit::MouthwashFit) = print(io, "MouthwashFit(", length(fit.mean), " genes)")

function Base.show(io::IO, ::MIME"text/plain", fit::MouthwashFit; ntop::Integer = 10, level::Real = 0.05)
    p = length(fit.mean)
    est = setdiff(1:p, fit.controls, fit.separated)
    printstyled(io, "Empirical Bayes fit (MOUTHWASH)"; bold = true)
    println(io)
    fit.converged || printstyled(io, "  not converged: increase maxiter or check the data\n"; color = :yellow)
    println(io, "  ", p, " genes (", length(fit.controls), " controls", isempty(fit.controls) ? ": confounders identified by sparsity" : "", ")")
    isempty(fit.separated) || println(io, "  ", length(fit.separated), " genes not tested: no counts in some group, level or time point, so some effect cannot be estimated")
    println(io, "  Confounders   α = [", join(fmt_x.(fit.α, 3), ", "), "] (first: compositional offset, if fitted) · variance inflation λ = ", fmt_x(fit.λ, 2))
    centre = fit.g.mode == 0 ? "" : string(", centred at the mode ", fmt_x(fit.g.mode, 3), " (shared shift of the non-null genes)")
    println(io, "  Effect prior  π₀ = ", fmt_x(fit.g.π[1], 3), " (proportion of null genes) over a grid of ", length(fit.g.σ), " scales", centre)
    disc = est[fit.svalue[est] .< level]
    up = count(j -> fit.mean[j] > 0, disc)
    println(io, "  Discoveries   ", length(disc), " genes with s-value < ", level, " (", up, " up, ", length(disc) - up, " down)")
    top = est[sortperm(fit.lfsr[est])][1:min(ntop, length(est))]
    println(io)
    println(io, "  Top ", length(top), " genes (log₂ fold changes; raw = confounded estimate before adjustment)")
    header = [("gene", 6), ("raw", 7), ("log₂FC", 7), ("sd", 6), ("lfsr", 6), ("s-value", 7)]
    rows = [
        (
            string(j), fmt_x(fit.β̃[j] / log(2)), fmt_x(fit.mean[j] / log(2)), fmt_x(fit.sd[j] / log(2)),
            fmt_p(fit.lfsr[j]), fmt_p(fit.svalue[j]),
        ) for j in top
    ]
    print_table(io, header, rows)
    return nothing
end
