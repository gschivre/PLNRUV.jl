# PLNRUV.jl

*Differential expression from RNA-seq counts, adjusted for hidden confounders and for the compositional nature of
sequencing.*

Standard RNA-seq normalisation assumes that the total amount of mRNA per cell does not change between conditions. That
assumption fails in many systems: hypertranscription in stem cells and cancers, global shutdown in quiescence, changes
of cell size and ploidy, a seedling exposed to light. Sequencing measures how RNA is *shared* among genes, not how
much there is, so a global change is invisible. Hidden confounders (batch, sample quality, cell-type mix) add their
own bias. PLNRUV.jl treats both problems with one model:

1. **Poisson-lognormal counts**: independent Poisson counts on a latent Gaussian layer. Zeros and low counts need no
   pseudo-count.
2. **Hidden confounders as latent factors**, separated from the effects of interest by the rotation of Wang et al.
   (CATE), and identified by negative control genes or by the sparsity of the effects (RUV-4, RUVB, MOUTHWASH of
   Gerard & Stephens).
3. **The compositional offset as one more factor**, whose loadings are fixed to one: a shift common to all genes of a
   sample. Controls or sparsity estimate it together with the other confounders.

Two branches share the model:

- an **empirical Bayes** branch (`pln_mouthwash`): variational PLN fit, then MOUTHWASH with an ash prior on the
  effects and vash/limma shrinkage of the variances. There is no prior to choose. It reports posterior means, local
  false sign rates and s-values;
- a **fully Bayesian** branch (`ruvb`): a Gibbs sampler that treats the effects as missing data and propagates the
  uncertainty of the confounders. It currently works on a Gaussian scale (log-ratios).

The package is pure Julia. Comparisons with R packages live in the separate `benchmark/` environment.

## Status

Early development. The components are validated on simulations, but not yet benchmarked on real data. The interface
may change.

## Installation

The package is not registered. Install it from its repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/gschivre/PLNRUV.jl")
```

## Quick start

```julia
using PLNRUV

# counts: samples x genes (integers); X: design with an intercept, the covariate of interest last;
# controls: indices of negative control genes (optional)
fit = pln_mouthwash(counts, X, 2; controls = controls)   # 2 latent factors + the compositional offset
res = fit.mouthwash
res.mean     # posterior mean log fold change (natural log)
res.lfsr     # local false sign rate
res.svalue   # declare the genes with svalue < 0.05

fit_free = pln_mouthwash(counts, X, 2)                   # without controls: sparsity identifies the confounders
```

## Simulating realistic data

`thin_2group` creates two-group experiments with known effects from real counts, by binomial thinning (Gerard 2020).
It can add negative controls and groups correlated with hidden factors. The validation uses GTEx v8 gene read counts
(skeletal muscle). **The GTEx data are not distributed with this package.** Download them from the
[GTEx portal](https://gtexportal.org) under its terms of use, and read them with `read_gct`.

## Learn more

The manual explains the method step by step, from the biology to the equations. It is written in
[`docs/src`](docs/src):

- [the problem](docs/src/problem.md): global transcription changes, compositions and confounders;
- [the model](docs/src/model.md), [separating the confounders](docs/src/confounders.md),
  [empirical Bayes](docs/src/empirical_bayes.md), [fitting the PLN layer](docs/src/pln.md) and
  [the fully Bayesian branch](docs/src/bayes.md);
- [simulation by binomial thinning](docs/src/simulation.md) and [validation](docs/src/validation.md);
- [assumptions and limits](docs/src/limits.md).

Build the HTML manual with `julia --project=docs docs/make.jl`.

## Authorship and AI assistance

PLNRUV.jl is developed by Geoffrey Schivre. The research question is the author's: adapting the confounder adjustment
of Wang et al. to the compositional nature of sequencing data. It arose from the observation of a global increase of
transcription in light-grown cotyledons. The author's earlier Julia scripts (a compositional CATE, Bayesian factor
models, a binomial-thinning simulator and an ash-type shrinkage) were the starting point. The modelling directions
were also the author's: a Poisson-lognormal model instead of a multinomial one, empirical Bayes and fully Bayesian
branches inspired by the Stephens lab, variance shrinkage with vash, and validation by binomial thinning of GTEx data.

Claude (Opus 5.5, Anthropic), an AI model used as a coding and analysis agent through Claude Code, contributed under
the author's direction:

- the review of the earlier scripts and the development plan;
- the formulation of the compositional offset as a latent factor with loadings fixed to one, and the latent-scale
  standard errors of the PLN branch;
- the implementation of the package (simulation, RUVB sampler, ash with mix-SQP, MOUTHWASH, variance shrinkage,
  variational PLN), its tests, its cross-checks against seqgendiff and vicar, and the manual.

The author chose the questions and the modelling directions, and is responsible for the package and for its use in
publications.

## Licence

MIT, see [`LICENSE`](LICENSE).
