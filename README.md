# PLNRUV.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://gschivre.github.io/PLNRUV.jl/dev/)
[![CI](https://github.com/gschivre/PLNRUV.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/gschivre/PLNRUV.jl/actions/workflows/CI.yml)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Code Style: Runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

*Differential expression from RNA-seq counts, adjusted for hidden confounders and for the compositional nature of
sequencing.*

Sequencing measures the *share* of each gene among the reads, not its amount per cell, so a global change of
transcription is invisible to standard normalisation. Hidden confounders add their own bias. PLNRUV.jl treats both
with one model:

1. **Poisson-lognormal counts**, fitted by variational empirical Bayes: no pseudo-count, no prior to choose.
2. **Hidden confounders as latent factors**, separated from the effects as in CATE and MOUTHWASH.
3. **The compositional shift as one more factor**, with the same loading for every gene.

The confounders and the shift are identified by **negative control genes** or by **sparsity** (most genes do not
change). A global change is measured relative to the controls, and only as well as the controls are truly unchanged.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/gschivre/PLNRUV.jl")
```

## Quick start

```julia
using PLNRUV, StatsModels

# counts: samples x genes (integers); samples: a table with one row per sample; controls: control genes (optional)
f = @formula(0 ~ genotype * time)
q = choose_q(counts, f, samples).q                              # number of hidden factors
fit = pln_mouthwash(counts, f, samples, q; controls,
                    contrast = ["time: 6h" => 1, "time: 2h" => -1])   # any contrast; fit.coefnames lists the names
res = fit.mouthwash
res.mean     # posterior mean log fold change (natural log)
res.svalue   # report the genes with svalue < 0.05
res.α[1]     # the compositional shift, relative to the controls
```

Without StatsModels, give a design matrix `X` instead of the formula; the last column is tested by default.

## Status

Version 0.1, under development. Benchmarked on experiments built from GTEx skeletal muscle by binomial thinning,
against DESeq2, edgeR, limma-voom, CATE and MOUTHWASH ([validation](docs/src/validation.md)). With valid negative
controls PLNRUV keeps the false discovery rate near its nominal level from 5 samples per group, also when the groups
are confounded or most genes change; with 2 or 3 per group it is as close to it as an oracle. Controls that shift with
the condition bias every estimate ([limits](docs/src/limits.md)). A fully Bayesian version is kept on the
[`fully-bayesian`](https://github.com/gschivre/PLNRUV.jl/tree/fully-bayesian) branch.

The [manual](https://gschivre.github.io/PLNRUV.jl/dev/) explains the method step by step, from the biology to the
equations, with all references. The GTEx data used by the benchmark are not distributed with the package; download
them from the [GTEx portal](https://gtexportal.org).

## Authorship and AI assistance

PLNRUV.jl is developed by Geoffrey Schivre. The research question is his: adapting confounder adjustment to the
compositional nature of sequencing data, after observing a global increase of transcription in light-grown
*Arabidopsis* cotyledons (Schivre et al. 2025, *Nat. Commun.* 16, 11228). His earlier Julia scripts were the starting
point, and the modelling directions are his.

Claude (Opus 5.5, Anthropic), used as a coding and analysis agent through Claude Code, contributed under his
direction: the formulation of the compositional shift as a latent factor, the implementation, the tests, the
benchmarks and the manual. The author chose the questions and the methods, and is responsible for the package and its
use.

## Licence

MIT, see [`LICENSE`](LICENSE).
