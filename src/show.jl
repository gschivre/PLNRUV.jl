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

Base.show(io::IO, fit::RUVBFit) = print(
    io, "RUVBFit(", size(fit.β, 2), " genes, ", size(fit.β, 4), " chains × ", size(fit.β, 3), " draws)"
)

function Base.show(io::IO, ::MIME"text/plain", fit::RUVBFit; ntop::Integer = 10, level::Real = 0.05)
    k₂, p, S, C = size(fit.β)
    s = summarize(fit)
    est = setdiff(1:p, fit.controls)
    branch = isempty(fit.acceptance) ? "Gaussian scale" : "Poisson-lognormal"
    printstyled(io, "Fully Bayesian RUVB fit"; bold = true)
    println(io, " (", branch, ")")
    println(io, "  ", C, " chains × ", S, " draws · ", p, " genes (", length(fit.controls), " controls) · ", k₂, " covariate(s) of interest")
    #convergence
    R = filter(!isnan, vec(s.rhat[:, est]))
    E = filter(!isnan, vec(s.ess[:, est]))
    nbad = count(>(1.01), R)
    print(io, "  Convergence   R̂ ≤ 1.01 for ", length(R) - nbad, "/", length(R), " coefficients (max ", fmt_x(maximum(R), 3), ")")
    println(io, " · bulk ESS min ", round(Int, minimum(E)), ", median ", round(Int, median(E)))
    if !isempty(fit.acceptance)
        lo, hi = extrema(fit.acceptance)
        println(io, "                latent Metropolis-Hastings acceptance ", fmt_x(lo, 3), lo ≈ hi ? "" : "–" * fmt_x(hi, 3))
    end
    if nbad > 0
        printstyled(
            io, "  ⚠ ", nbad, " coefficient(s) have R̂ > 1.01: the chains have not mixed, run longer chains (burnin, nsample)\n";
            color = :yellow
        )
    end
    for i in 1:k₂
        sv = s.ash_svalue[i, est]
        disc = est[sv .< level]
        up = count(j -> s.ash_mean[i, j] > 0, disc)
        k₂ > 1 && println(io, "  Covariate ", fit.cov_of_interest[i], ":")
        println(io, "  Discoveries   ", length(disc), " genes with s-value < ", level, " (", up, " up, ", length(disc) - up, " down)")
        top = est[sortperm(s.ash_lfsr[i, est])][1:min(ntop, length(est))]
        println(io)
        println(io, "  Top ", length(top), " genes (log₂ fold changes; shrunk = after ash, use it with the s-value)")
        header = [("gene", 6), ("log₂FC", 7), ("95% interval", 17), ("shrunk", 7), ("lfsr", 6), ("s-value", 7), ("R̂", 6), ("ESS", 5)]
        rows = [
            (
                string(j), fmt_x(s.mean[i, j] / log(2)),
                "[" * fmt_x(s.lower[i, j] / log(2)) * ", " * fmt_x(s.upper[i, j] / log(2)) * "]",
                fmt_x(s.ash_mean[i, j] / log(2)), fmt_p(s.ash_lfsr[i, j]), fmt_p(s.ash_svalue[i, j]),
                fmt_x(s.rhat[i, j], 3), string(round(Int, s.ess[i, j])),
            ) for j in top
        ]
        print_table(io, header, rows)
    end
    return nothing
end

Base.show(io::IO, fit::MouthwashFit) = print(io, "MouthwashFit(", length(fit.mean), " genes)")

function Base.show(io::IO, ::MIME"text/plain", fit::MouthwashFit; ntop::Integer = 10, level::Real = 0.05)
    p = length(fit.mean)
    est = setdiff(1:p, fit.controls)
    printstyled(io, "Empirical Bayes fit (MOUTHWASH)"; bold = true)
    println(io)
    println(io, "  ", p, " genes (", length(fit.controls), " controls", isempty(fit.controls) ? ": confounders identified by sparsity" : "", ")")
    println(io, "  Confounders   α = [", join(fmt_x.(fit.α, 3), ", "), "] (first: compositional offset, if fitted) · variance inflation λ = ", fmt_x(fit.λ, 2))
    println(io, "  Effect prior  π₀ = ", fmt_x(fit.g.π[1], 3), " (proportion of null genes) over a grid of ", length(fit.g.σ), " scales")
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
