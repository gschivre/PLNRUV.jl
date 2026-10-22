# PLNRUV.jl

*Differential expression from RNA-seq counts, adjusted for hidden confounders and for the compositional nature of
sequencing.*

An RNA-seq experiment compares the expression of thousands of genes between conditions. Two things make the comparison
harder than it looks:

- **Hidden confounders.** Batch, the time a sample spent on the bench, the cell-type mix of a tissue: they shift the
  expression of many genes at once. When they are correlated with the condition, they create false differences, or
  hide true ones.
- **Composition.** A sequencer reads a fixed budget of molecules from each library. The reads tell us how the RNA of a
  sample is *shared* among genes, not how much RNA there is. If a few genes are induced massively, every other gene
  loses reads, and looks repressed. And when the whole transcriptome rises or falls, as in stem cells, cancers,
  quiescence or a seedling exposed to light, standard normalisation hides the change.

This package handles both with a single device. It models the counts with a Poisson-lognormal (PLN) distribution. It
treats the hidden confounders as latent factors, following Wang et al.'s CATE and Gerard & Stephens' RUVB and
MOUTHWASH. And it shows that the compositional offset is just one more latent factor, whose loadings are known: it
shifts every gene of a sample by the same amount. Negative control genes, or the assumption that most genes do not
change, identify all of these factors at once.

## The story in six steps

1. [**The problem**](problem.md). Counts, compositions, and confounders. Why the change of a gene is not identified
   from the reads alone, and what restores it.
2. [**The model**](model.md). Independent Poisson counts on top of a latent Gaussian layer: a regression on the design,
   latent factors, and gene-specific overdispersion. The compositional offset becomes a factor with loadings fixed to
   one.
3. [**Separating the confounders from the effects**](confounders.md). The rotation of Wang et al. splits the samples
   into a part that carries the effects and a part that only sees the confounders. Controls or sparsity then pin down
   the confounders.
4. [**Empirical Bayes: ash, MOUTHWASH and variance shrinkage**](empirical_bayes.md). The prior on the effects, the
   confounders and the variances are all estimated from the data. No prior to choose.
5. [**Fitting the Poisson-lognormal layer**](pln.md). A variational EM whose E-step has a closed form (Lambert W),
   and whose M-step is the Gaussian machinery of step 3.
6. [**The fully Bayesian branch**](bayes.md). RUVB: a Gibbs sampler that treats the effects as missing data and
   propagates every source of uncertainty.

Then: [**simulating data by binomial thinning**](simulation.md), [**validation**](validation.md),
[**assumptions and limits**](limits.md), [notation and references](notation.md), and the [API](api.md).

## How to read this manual

The main text is written for a biologist who is comfortable with probabilities. It explains every idea in words, then
states it as an equation. Each modelling choice is labelled:

- **assumption**: a choice about the biology or the measurement, which could be wrong;
- **result**: a mathematical consequence of the assumptions, shown in the text.

The derivations are in folded boxes like this one:

!!! details "Derivation"
    The step-by-step argument behind a result, for readers who want to check it. The main text never depends on it.

Each page ends with the functions that implement it.

## Quick start

```julia
using PLNRUV, Random, Distributions

# real counts (genes x samples) from GTEx, thinned to create a two-group experiment with known effects
gtex = read_gct("gene_reads_muscle_skeletal.gct")
M = gtex.counts[filter_genes(gtex.counts; min_mean = 10), :]
sim = thin_2group(Xoshiro(1), M; ngenes = 5_000, nsamp = 20, nctrl = 500, prop_null = 0.8,
    signal = Normal(0, 0.8), target_cor = [0.6, 0.2])   # the group is correlated with two hidden factors

# empirical Bayes branch: PLN fit with 2 latent factors + the compositional offset, then MOUTHWASH
fit = pln_mouthwash(sim.counts, sim.X, 2; controls = sim.controls)
res = fit.mouthwash
res.mean     # posterior mean of each log fold change (natural log; divide by log(2) for log2)
res.lfsr     # local false sign rate: probability that the sign of the effect is wrong
res.svalue   # s-value: declare the genes with svalue < 0.05 to control the false sign rate at 5 %

# without control genes, the sparsity of the effects identifies the confounders
fit_free = pln_mouthwash(sim.counts, sim.X, 2)

# fully Bayesian branch on a Gaussian scale (log-ratios), with the posterior of every coefficient
post = summarize(ruvb(log_ratio(sim.counts), sim.X, sim.controls, 2))
post.mean, post.lower, post.upper, post.lfsr
```
