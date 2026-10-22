#benchmark of PLNRUV against R methods on binomially thinned GTEx data
#run from the benchmark environment: julia -t auto --project=benchmark benchmark/gtex_benchmark.jl <gtex counts file> [nrep]
#results are written to benchmark/results/ (one CSV row per method, scenario and replicate)
using PLNRUV, RCall, Random, Distributions, Statistics

R"""
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(limma); library(cate); library(vicar)})
"""

#each R method returns log2 fold changes and adjusted p-values (BH) or s-values for the group effect
function r_methods(counts, X, controls, q)
    Yt = permutedims(counts) #genes x samples for Bioconductor
    grp = X[:, 2]
    ctrl = collect(controls)
    out = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    R"""
    grp <- factor($grp)
    dds <- DESeqDataSetFromMatrix($Yt, data.frame(grp = grp), ~ grp)
    res <- results(DESeq(dds, quiet = TRUE), name = "grp_1_vs_0")
    deseq <- cbind(res$log2FoldChange, res$padj)
    dds <- estimateSizeFactors(dds, controlGenes = $ctrl)
    res <- results(DESeq(dds, quiet = TRUE), name = "grp_1_vs_0")
    deseq_ctrl <- cbind(res$log2FoldChange, res$padj)
    y <- calcNormFactors(DGEList($Yt))
    design <- model.matrix(~ grp)
    y <- estimateDisp(y, design)
    tt <- topTags(glmQLFTest(glmQLFit(y, design)), n = Inf, sort.by = "none")$table
    edger <- cbind(tt$logFC, tt$FDR)
    v <- voom(y, design)
    tt <- topTable(eBayes(lmFit(v, design)), coef = 2, n = Inf, sort.by = "none")
    voomr <- cbind(tt$logFC, tt$adj.P.Val)
    logcpm <- t(cpm(y, log = TRUE))
    ct <- cate(~ grp | 1, data.frame(grp = grp), Y = logcpm, r = $q, adj.method = "nc", nc = $ctrl, calibrate = FALSE)
    cater <- cbind(as.numeric(ct$beta), p.adjust(as.numeric(ct$beta.p.value), "BH"))
    mw <- mouthwash(Y = logcpm, X = design, k = $q, cov_of_interest = 2, include_intercept = FALSE, verbose = FALSE)
    vicarr <- cbind(mw$result$PosteriorMean, mw$result$svalue)
    """
    for (name, sym) in [
            ("DESeq2", "deseq"), ("DESeq2_controls", "deseq_ctrl"), ("edgeR", "edger"), ("limma_voom", "voomr"),
            ("cate_nc", "cater"), ("vicar_mouthwash", "vicarr"),
        ]
        m = rcopy(reval(sym))
        out[name] = (m[:, 1], replace(m[:, 2], NaN => 1.0))
    end
    return out
end

#PLNRUV methods: log2 fold changes and s-values
function julia_methods(counts, X, controls, q; bayes = true)
    out = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    p = size(counts, 2)
    eb = pln_mouthwash(counts, X, q; controls = controls).mouthwash
    out["PLNRUV_EB_controls"] = (eb.mean ./ log(2), eb.svalue)
    ebf = pln_mouthwash(counts, X, q).mouthwash
    out["PLNRUV_EB_sparsity"] = (ebf.mean ./ log(2), ebf.svalue)
    g = mouthwash(log_ratio(counts, controls), X, q; controls = controls)
    out["PLNRUV_gaussian_logratio"] = (g.mean ./ log(2), g.svalue)
    if bayes
        s = summarize(pln_ruvb(counts, X, controls, q; burnin = 1_000, nsample = 250, thin = 4).fit)
        m, sv = zeros(p), ones(p)
        est = setdiff(1:p, controls)
        m[est] .= s.ash_mean[1, est] ./ log(2)
        sv[est] .= s.ash_svalue[1, est]
        out["PLNRUV_bayes"] = (m, sv)
    end
    return out
end

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

const SCENARIOS = [
    (name = "balanced", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = Float64[]),
    (name = "confounded", prop_null = 0.8, signal = Normal(0.0, 0.8), target_cor = [0.6, 0.2]),
    (name = "global_increase", prop_null = 0.4, signal = Normal(1.0, 0.5), target_cor = [0.6, 0.2]),
]

function run_benchmark(M; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1, bayes = true, outfile = nothing)
    rows = []
    for sc in SCENARIOS, r in 1:nrep
        rng = Xoshiro(seed + r)
        sim = thin_2group(
            rng, M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = sc.prop_null,
            signal = sc.signal, target_cor = sc.target_cor
        )
        idx = setdiff(1:ngenes, sim.controls)
        res = merge(r_methods(sim.counts, sim.X, sim.controls, q), julia_methods(sim.counts, sim.X, sim.controls, q; bayes = bayes))
        for (method, (lfc, qv)) in res
            push!(rows, (scenario = sc.name, rep = r, method = method, metrics(lfc, qv, sim.β, idx)...))
        end
        if outfile !== nothing
            open(outfile, "w") do io
                println(io, "scenario,rep,method,ndisc,fdp,power,rmse")
                for x in rows
                    println(io, join((x.scenario, x.rep, x.method, x.ndisc, x.fdp, x.power, x.rmse), ","))
                end
            end
        end
    end
    return rows
end

#mean of each metric by scenario and method
function summarize_benchmark(rows)
    keys_ = unique((x.scenario, x.method) for x in rows)
    return [
        (
            scenario = s, method = m,
            fdp = mean(x.fdp for x in rows if x.scenario == s && x.method == m),
            power = mean(x.power for x in rows if x.scenario == s && x.method == m),
            rmse = mean(x.rmse for x in rows if x.scenario == s && x.method == m),
        ) for (s, m) in keys_
    ]
end

if abspath(PROGRAM_FILE) == @__FILE__
    gtex = read_gct(ARGS[1])
    M = gtex.counts[filter_genes(gtex.counts; min_mean = 10), :]
    nrep = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 10
    mkpath(joinpath(@__DIR__, "results"))
    rows = run_benchmark(M; nrep = nrep, outfile = joinpath(@__DIR__, "results", "gtex_benchmark.csv"))
    foreach(println, summarize_benchmark(rows))
end
