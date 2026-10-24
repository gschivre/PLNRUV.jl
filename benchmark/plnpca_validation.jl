#validation of the experimental plnpca fit against pln (decision E of design/plnpca.md): small-sample calibration on
#PLN-exact simulations and on thinned GTEx, both fits with q chosen by choose_q, on the same data sets
#run from the benchmark environment, after loading the GTEx counts M (genes x samples, see gct.jl):
#    include("benchmark/plnpca_validation.jl"); run_smalln(M)
#rows are appended to results/plnpca_smalln.csv after each data set (an interrupted run resumes)
using PLNRUV, Random, Distributions, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "gct.jl"))

const FITS = [:pln, :plnpca]
const NPG = [2, 3, 5]
const SMALLN_SCENARIOS = [
    (name = "balanced", target_cor = Float64[]),
    (name = "confounded", target_cor = [0.6, 0.2]),
]

#PLN-exact two-group data with GTEx-like gene means and library sizes: latent variances that fall with expression
#(about 0.45 below 40 reads, 0.06 above 2500, log slope -0.5), two hidden factors (loadings sd 0.4) correlated with
#the group as in thin_2group, 20 % of the non-control genes changed (log2 fold changes N(0, 0.8))
function sim_exact(rng, M, nsamp, sc; ngenes = 2_000, nctrl = 200, prop_null = 0.8)
    genes = sort(randperm(rng, size(M, 1))[1:ngenes])
    libs = vec(sum(M; dims = 1))
    logmean = log.(vec(mean(M[genes, :] ./ libs' .* median(libs); dims = 2)) .+ 0.5)
    σ2 = clamp.(0.45 .* exp.(-0.5 .* (logmean .- log(40))), 0.06, 0.45)
    n = nsamp
    U = randn(rng, n, 2)
    grp = assign_groups(rng, U, [sc.target_cor; zeros(2 - length(sc.target_cor))])
    X = [ones(n) grp]
    β = zeros(ngenes)
    nz = (nctrl + 1):ngenes
    β[nz] .= (rand(rng, length(nz)) .> prop_null) .* rand(rng, Normal(0, 0.8), length(nz))
    Γ = 0.4 .* randn(rng, 2, ngenes)
    o = log.(rand(rng, libs, n) ./ median(libs))
    Z = logmean' .+ grp .* (β' .* log(2)) .+ U * Γ .+ randn(rng, n, ngenes) .* sqrt.(σ2)'
    Y = [rand(rng, Poisson(exp(o[i] + Z[i, j]))) for i in 1:n, j in 1:ngenes]
    return (counts = Y, X = X, β = β, controls = collect(1:nctrl))
end

function sim_gtex(rng, M, nsamp, sc; ngenes = 2_000, nctrl = 200)
    return thin_2group(
        rng, M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = 0.8, signal = Normal(0.0, 0.8),
        target_cor = sc.target_cor
    )
end

#FDP (null or wrong sign) and power at s-value < 0.05, log2 fold changes
function smalln_metrics(mean_lfc, svalue, β, controls)
    idx = setdiff(eachindex(β), controls)
    d = idx[svalue[idx] .< 0.05]
    wrong = (β[d] .== 0) .| (sign.(β[d]) .!= sign.(mean_lfc[d]))
    return (ndisc = length(d), fdp = isempty(d) ? 0.0 : mean(wrong), power = sum(.!wrong) / max(count(!=(0), β[idx]), 1))
end

#one fit: q by choose_q then pln_mouthwash with controls, both with the same fit; a failure is recorded, not thrown
function smalln_fit(sim, fit)
    t = @elapsed res = try
        cq = choose_q(sim.counts, sim.X; qs = 0:max(0, min(4, size(sim.X, 1) - 4)), fit = fit, rng = Xoshiro(1))
        mw = pln_mouthwash(sim.counts, sim.X, cq.q; controls = sim.controls, fit = fit)
        (q = cq.q, mw = mw, error = "")
    catch err
        err isa InterruptException && rethrow()
        (q = -1, mw = nothing, error = first(sprint(showerror, err), 80))
    end
    res.mw === nothing && return (q = -1, ndisc = 0, fdp = NaN, power = NaN, converged = false, time = t, error = res.error)
    x = smalln_metrics(res.mw.mouthwash.mean ./ log(2), res.mw.mouthwash.svalue, sim.β, sim.controls)
    return (q = res.q, x..., converged = res.mw.pln.converged && res.mw.mouthwash.converged, time = t, error = "")
end

function run_smalln(M; nrep = 10, seed = 1, file = joinpath(@__DIR__, "results", "plnpca_smalln.csv"))
    header = "data,scenario,npergroup,rep,fit,q,ndisc,fdp,power,converged,time,error"
    isfile(file) || write(file, header * "\n")
    done = Set(Tuple(split(l, ",")[1:5]) for l in Iterators.drop(eachline(file), 1))
    for data in ("exact", "gtex"), sc in SMALLN_SCENARIOS, g in NPG, r in 1:nrep
        all((data, sc.name, string(g), string(r), string(f)) in done for f in FITS) && continue
        rng = Xoshiro(seed + 100g + r)
        sim = data == "exact" ? sim_exact(rng, M, 2g, sc) : sim_gtex(rng, M, 2g, sc)
        for f in FITS
            (data, sc.name, string(g), string(r), string(f)) in done && continue
            x = smalln_fit(sim, f)
            row = join((data, sc.name, g, r, f, x.q, x.ndisc, x.fdp, x.power, x.converged, round(x.time; digits = 1), replace(x.error, "," => ";")), ",")
            open(io -> println(io, row), file, "a")
        end
    end
    return nothing
end
