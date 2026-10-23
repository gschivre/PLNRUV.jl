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
    #each method in its own tryCatch: one failing on a small design (e.g. cate or vicar with q = 0) gives NULL
    R"""
    grp <- factor($grp)
    design <- model.matrix(~ grp)
    safe <- function(expr) tryCatch(expr, error = function(e) NULL)
    deseq <- safe({
        dds <- DESeqDataSetFromMatrix($Yt, data.frame(grp = grp), ~ grp)
        res <- results(DESeq(dds, quiet = TRUE), name = "grp_1_vs_0")
        cbind(res$log2FoldChange, res$padj)
    })
    deseq_ctrl <- safe({
        dds <- estimateSizeFactors(DESeqDataSetFromMatrix($Yt, data.frame(grp = grp), ~ grp), controlGenes = $ctrl)
        res <- results(DESeq(dds, quiet = TRUE), name = "grp_1_vs_0")
        cbind(res$log2FoldChange, res$padj)
    })
    y <- estimateDisp(calcNormFactors(DGEList($Yt)), design)
    edger <- safe({
        tt <- topTags(glmQLFTest(glmQLFit(y, design)), n = Inf, sort.by = "none")$table
        cbind(tt$logFC, tt$FDR)
    })
    voomr <- safe({
        tt <- topTable(eBayes(lmFit(voom(y, design), design)), coef = 2, n = Inf, sort.by = "none")
        cbind(tt$logFC, tt$adj.P.Val)
    })
    logcpm <- t(cpm(y, log = TRUE))
    cater <- safe({
        ct <- cate(~ grp | 1, data.frame(grp = grp), Y = logcpm, r = $q, adj.method = "nc", nc = $ctrl, calibrate = FALSE)
        cbind(as.numeric(ct$beta), p.adjust(as.numeric(ct$beta.p.value), "BH"))
    })
    vicarr <- safe({
        #vicar defaults: normal likelihood, limma shrinkage
        mw <- mouthwash(Y = logcpm, X = design, k = $q, cov_of_interest = 2, include_intercept = FALSE, verbose = FALSE)
        cbind(mw$result$PosteriorMean, mw$result$svalue)
    })
    """
    for (name, sym) in [
            ("DESeq2", "deseq"), ("DESeq2_controls", "deseq_ctrl"), ("edgeR", "edger"), ("limma_voom", "voomr"),
            ("cate_nc", "cater"), ("vicar_mouthwash", "vicarr"),
        ]
        rcopy(reval("is.null($sym)")) && continue
        m = rcopy(reval(sym))
        #genes the method cannot test (NA in R: all-zero counts, outliers) get no effect and are never called
        lfc = Float64[ismissing(x) || isnan(x) ? 0.0 : x for x in m[:, 1]]
        padj = Float64[ismissing(x) || isnan(x) ? 1.0 : x for x in m[:, 2]]
        out[name] = (lfc, padj)
    end
    return out
end

#DESeq2 median-of-ratios size factors computed on a set of reference genes (rows of d: genes, columns: samples)
function size_factors(d, ref)
    ratios = [d[g, :] ./ exp(mean(log.(d[g, :]))) for g in ref if minimum(d[g, :]) > 0]
    sf = [median(r[i] for r in ratios) for i in axes(d, 2)]
    return sf ./ exp(mean(log.(sf)))
end

#the spike-in/reference CATE pipeline of Schivre et al. (Nat. Commun. 2025): size factors on the controls,
#CATE (ESA, negative controls) on the log2 counts, corrected counts, size factors recomputed on the corrected
#controls, then DESeq2 with these size factors
function paper_method(counts, X, controls, q)
    Yt = permutedims(counts)
    ctrl = collect(controls)
    sf = size_factors(Yt, ctrl)
    logcts = log2.(counts .+ 0.5) .- log2.(sf)
    grp = X[:, 2]
    #without factor (q = 0, very small designs) there is nothing to correct: DESeq2 with control size factors
    corrected = Yt
    if q > 0
        R"""
        ct <- cate::cate(~ grp | 1, data.frame(grp = factor($grp)), Y = $logcts, r = $q, fa.method = "esa",
                         adj.method = "nc", nc = $ctrl, nc.var.correction = FALSE, calibrate = FALSE)
        """
        ZΓ = rcopy(R"ct$Z %*% t(ct$Gamma)")
        corrected = round.(Int, max.(0.0, permutedims(2 .^ (logcts .- ZΓ) .* sf .- 0.5)))
    end
    sfc = size_factors(corrected, ctrl)
    R"""
    dds <- DESeq2::DESeqDataSetFromMatrix($corrected, data.frame(grp = factor($grp)), ~ grp)
    DESeq2::sizeFactors(dds) <- $sfc
    res <- DESeq2::results(DESeq2::DESeq(dds, quiet = TRUE), name = "grp_1_vs_0")
    paper <- cbind(res$log2FoldChange, res$padj)
    """
    m = rcopy(R"paper")
    lfc = Float64[ismissing(x) || isnan(x) ? 0.0 : x for x in m[:, 1]]
    padj = Float64[ismissing(x) || isnan(x) ? 1.0 : x for x in m[:, 2]]
    return Dict("paper_RNARx_CATE" => (lfc, padj))
end

#run one extra method on the same simulated data sets as run_benchmark (same seeds) and append to the CSV
function run_extra(M, method; nrep = 10, ngenes = 2_000, nctrl = 200, nsamp = 20, q = 2, seed = 1, outfile = nothing)
    rows = []
    for sc in SCENARIOS, r in 1:nrep
        rng = Xoshiro(seed + r)
        sim = thin_2group(
            rng, M; ngenes = ngenes, nsamp = nsamp, nctrl = nctrl, prop_null = sc.prop_null,
            signal = sc.signal, target_cor = sc.target_cor
        )
        idx = setdiff(1:ngenes, sim.controls)
        for (name, (lfc, qv)) in method(sim.counts, sim.X, sim.controls, q)
            push!(rows, (scenario = sc.name, rep = r, method = name, metrics(lfc, qv, sim.β, idx)...))
        end
    end
    if outfile !== nothing
        open(outfile, "a") do io
            for x in rows
                println(io, join((x.scenario, x.rep, x.method, x.ndisc, x.fdp, x.power, x.rmse), ","))
            end
        end
    end
    return rows
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
        s = summarize(pln_ruvb(counts, X, controls, q; burnin = 1_000, nsample = 250, thin = 4))
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
