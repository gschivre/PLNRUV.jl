# 6. Simulation by binomial thinning

A method for real data should be tested on data that look real: real correlations between genes, real
overdispersion, real zeros and real hidden factors. Parametric simulations get some of these wrong. Binomial thinning
(Gerard 2020) starts from a real count matrix and adds known effects to it.

## 6.1 Thinning adds an exact effect

Take two groups of real samples, with ``x_i = 1`` in group 1 and ``x_i = 0`` in group 0, and a log₂ fold change
``\beta_j``. Keep each read of gene ``j`` in sample ``i`` independently with probability

```math
\rho_{ij} = 2^{\,x_i\beta_j - \max_l x_l\beta_j} .
```

**Result** — if the original count is Poisson given any random mean, the thinned count is Poisson with that mean
times ``\rho_{ij}``. So the expected count of group 1 is multiplied by exactly ``2^{\beta_j}`` relative to group 0.
Everything else in the data is kept: correlations, factors, overdispersion.

Thinning the changed genes also lowers the library sizes of the samples where they were thinned. This is the
compositional effect of [§1.1](problem.md), created naturally.

## 6.2 Controls and the truth

A set of genes is drawn as negative controls, with ``\beta_j = 0``. To test robustness, controls can also carry small
effects (`ctrl_signal`). Among the other genes, a proportion `prop_null` is null. The others get effects drawn from
`signal`. Because thinning acts on the expected counts, the true effect of each gene relative to the controls is
``\beta_j``.

## 6.3 Groups correlated with hidden factors

To create confounding, the groups are not assigned at random:

1. The factors ``U`` of the real samples are estimated by a PCA of their log counts. Each sample is centred first,
   which removes sequencing depth; then each gene is centred.
2. A score ``\zeta_i = U_i^\top h + \sqrt{1 - \lVert h\rVert^2}\,\epsilon_i``, with ``\epsilon_i \sim \mathcal N(0,
   1)``, is drawn for each sample. The samples with the largest scores form group 1.

**Result** — if a fraction ``f`` of the samples is in group 1, the correlation between the group and factor ``l`` is
``h_l\,\phi(z_f)/\sqrt{f(1-f)}``, where ``\phi`` is the standard normal density and ``z_f`` its ``1 - f`` quantile.
The weights ``h`` are set to reach the target correlations (`target_cor`).

!!! details "Derivation"
    ``U_l`` and ``\zeta`` are jointly normal with unit variances and covariance ``h_l``. Stein's lemma gives
    ``\mathbb E[U_l \mathbf 1\{\zeta > z_f\}] = h_l\,\phi(z_f)``. The group indicator has variance ``f(1-f)``.

The samples are centred on purpose. Sequencing depth is technical, and every method corrects for it. The confounding
that matters comes from the biological and unwanted variation of [§1.4](problem.md). seqgendiff does not centre the
samples, so its first factor mostly follows depth.

## Functions

- [`thin_2group`](@ref): the complete two-group simulation.
- [`thin_counts`](@ref): thinning for any design and coefficient matrix.
- [`surrogate_factors`](@ref), [`assign_groups`](@ref): the confounded group assignment.
- [`filter_genes`](@ref): keep the genes observed often enough.

Reading the GTEx count files is part of the benchmark scripts (`benchmark/gct.jl`), not of the package.
