# 6. Simulating data by binomial thinning

A method for real data should be tested on data that look real: with the real correlations between genes, the real
overdispersion, the real zeros and the real hidden factors. Parametric simulations get all of these wrong in some way.
Binomial thinning (Gerard 2020) starts from a real count matrix, here GTEx skeletal muscle (803 samples), and adds
known effects to it.

## 6.1 Thinning adds an effect exactly

Take a gene in a set of real samples, and a log fold change ``\beta_j`` (in log2 units) between two groups. Keep each
read of sample ``i`` independently with probability

```math
\pi_{ij} = 2^{\,x_i\beta_j - \max_l x_l\beta_j} ,
```

where ``x_i \in \{0, 1\}`` is the group. The thinned count is ``\text{Binomial}(y_{ij}, \pi_{ij})``.

**Result** — if the original count is Poisson with any random mean ``\Lambda_{ij}`` (a Poisson-lognormal count, for
instance), the thinned count is Poisson with mean ``\pi_{ij}\Lambda_{ij}``. So the expected count of group 1 is
multiplied by exactly ``2^{\beta_j}`` relative to group 0, and everything else in the data (correlations, factors,
overdispersion) is kept.

Thinning the differentially expressed genes reduces the library sizes of the samples where they were thinned. This
is the compositional effect of [§1.4](problem.md), created naturally: the simulated data are compositional in the
same way as real data.

## 6.2 Negative controls and the truth

A set of ``n_c`` genes is drawn as negative controls, with ``\beta_j = 0`` exactly. They can also be made
approximately null (`ctrl_signal = Normal(0, 0.1)`), to test robustness. Among the other genes, a proportion
`prop_null` is null, and the effects of the rest are drawn from `signal`. Because thinning acts on the expected
counts, the true effect relative to the controls is ``\beta_j`` itself.

## 6.3 Confounding: groups correlated with hidden factors

In real data the group is often correlated with hidden factors. To reproduce this, the groups are not assigned at
random:

1. the hidden factors ``U`` of the subsample are estimated by the PCA of its doubly centred log counts: each sample
   is centred (a centred log-ratio, which removes sequencing depth), then each gene;
2. a latent score ``w = U r + \sqrt{1 - \lVert r\rVert^2}\,\epsilon`` is drawn, and the samples with the largest
   ``w`` form group 1.

**Result** — for a Gaussian ``w`` and a group proportion ``p``, the correlation between the group indicator and
factor ``k`` is ``r_k\,\varphi(t)/\sqrt{p(1-p)}``, with ``t`` the ``1-p`` quantile of the standard normal. So
``r`` is chosen to hit a target correlation (`target_cor`).

!!! details "Derivation"
    For jointly Gaussian ``(U_k, w)`` with unit variances and covariance ``r_k``, Stein's lemma gives
    ``\mathbb E[U_k \mathbf 1\{w > t\}] = r_k\,\varphi(t)``. The indicator has variance ``p(1-p)`` and ``U_k`` has
    variance one, which gives the correlation.

The factors are centred per sample on purpose. Sequencing depth is technical, and every method normalises for it.
The confounding that matters comes from biological or unwanted variation. (seqgendiff does not centre the samples,
so its leading factor largely tracks depth.)

## Functions

- [`read_gct`](@ref), [`filter_genes`](@ref): read GTEx count files and filter genes.
- [`thin_2group`](@ref): the complete two-group simulation.
- [`thin_counts`](@ref): thinning for an arbitrary design and coefficient matrix.
- [`surrogate_factors`](@ref), [`assign_groups`](@ref): the confounded group assignment.
