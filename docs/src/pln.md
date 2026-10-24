# 5. Fitting the PLN layer

The latent log-abundances ``Z_{ij}`` of [§2](model.md) are not observed. The likelihood of a count is an integral
with no closed form:

```math
p(y_{ij} \mid \mu_{ij}, \sigma_j^2) = \int \text{Poisson}\big(y_{ij};\, e^{o_i + z}\big)\,
\mathcal N(z;\, \mu_{ij}, \sigma_j^2)\, dz .
```

The package maximises instead a lower bound of the log-likelihood, the ELBO, by variational inference (Chiquet et al.
2018). The fit, [`plnpca`](@ref), also treats the factors, the loadings, the coefficients and the variances as
random, with priors estimated from all genes. This is variational empirical Bayes.

## 5.1 The model as a hierarchy

[`plnpca`](@ref) completes the model of [§2](model.md) with priors:

- **Factors.** The scores of the ``q`` hidden factors are random, ``\mathcal N(0, 1)``, and live in the residual space
  of the design (``X^\top W = 0``). Their loadings are learned only from the variation that the design cannot explain,
  as in the rotation of [§3](confounders.md). A factor therefore cannot absorb the effects of interest.
- **Loadings.** ``\Gamma_{lj} \sim \mathcal N(0, \tau_l^2)``, with one variance per factor. A factor that the data do
  not support gets ``\tau_l^2 \to 0`` and is switched off (automatic relevance determination).
- **Coefficients.** ``B_j`` has a weak normal prior. It keeps the estimates finite for genes that are zero in a whole
  group.
- **Variances.** ``\sigma_j^2 \sim \text{InvGamma}(a_0,\, a_0 \sigma_0^2(x_j))``, the limma prior with ``d_0 = 2a_0``
  ([§4.5](empirical_bayes.md)), whose scale ``\sigma_0^2`` is a smooth function of the mean log-abundance ``x_j`` of
  gene ``j`` (§5.4).

The offset factor ``w_{i0}`` stays a free parameter, also in the residual space of the design. The hyperparameters
``\tau_l^2``, ``a_0`` and ``\sigma_0^2(\cdot)`` are estimated from all genes at every iteration.

## 5.2 The variational approximation

- **Approximation** — the posterior is replaced by a product of simple distributions: a normal
  ``\mathcal N(M_{ij}, S^2_{ij})`` for each latent value, normals for the scores, and for each gene a joint
  normal–inverse-gamma distribution of its coefficients, loadings and variance.

Only the latent values ``Z_{ij}`` see the Poisson counts. Given them, the rest of the model is a linear Gaussian
factor model. So every update has a closed form, and each one increases the ELBO, except in the first iterations on
low counts (§5.4). The fit stops when the ELBO changes
by less than a tolerance per count; the field `converged` reports whether this happened before `maxiter`.

## 5.3 The latent values: a closed form

For one count ``y_{ij}``, the ELBO is concave in ``M_{ij}``. **Result** — its maximum has a closed form:

```math
M_{ij} = \mu_{ij} + \sigma_j^2 y_{ij} - \mathcal W_0\Big(\sigma_j^2\, e^{\,o_i + \mu_{ij} + \sigma_j^2 y_{ij} +
S^2_{ij}/2}\Big),
```

where ``\mathcal W_0`` is the Lambert W function (the solution of ``w e^w = x``), ``\mu_{ij}`` the current mean and
``\sigma_j^2`` the current value of ``1 / \mathbb E[1/\sigma_j^2]``. The variance is the solution of ``S^2_{ij} = 1 /
(e^{\,o_i + M_{ij} + S^2_{ij}/2} + 1/\sigma_j^2)``, found by fixed-point iteration.

!!! details "Derivation"
    The derivative in ``M`` is zero when ``y - e^{a + M} - (M - \mu)/\sigma^2 = 0``, with ``a = o + S^2/2``. Let
    ``u = \mu + \sigma^2 y - M``. The equation becomes ``u\,e^{u} = \sigma^2 e^{a + \mu + \sigma^2 y}``, so ``u =
    \mathcal W_0(\sigma^2 e^{a + \mu + \sigma^2 y})``. For large counts the argument overflows. The package therefore
    computes ``\mathcal W_0(e^t)`` directly, by Newton's method on ``w + \log w = t``. The derivative in ``S^2`` is
    zero when ``1/S^2 = e^{a + M} + 1/\sigma^2``.

## 5.4 The variances

**Degrees of freedom.** Gene ``j`` has ``n`` counts, one per sample, but its own coefficients and loadings use some of them. Their
number is the trace of the gene's hat matrix, ``t_j \approx k + q``. The posterior of the variance counts only the
rest:

```math
q(\sigma_j^2) = \text{InvGamma}\Big(a_0 + \tfrac{n - t_j}{2},\;\; a_0 \sigma_0^2(x_j) + \tfrac{1}{2}\,\text{RSS}_j\Big),
```

where ``\text{RSS}_j`` is the expected residual sum of squares of the latent values. Its mean is the restricted
(REML) estimate, as in limma ([§4.5](empirical_bayes.md)).

**Why ``n - t_j`` and not ``n``.** A plain product of independent distributions (mean field) for the coefficients
and the variance gives the shape ``a_0 + n/2``. It treats every count as a free observation of ``\sigma_j^2``. The
mean is the same, but the posterior is too narrow. With few samples this matters: with 4 samples, 2 of the 4
degrees of freedom go to the coefficients. The estimate of ``a_0`` then grows far too large. All variances are pulled
to one value, the fit converges very slowly, and low-count genes get false discoveries. The joint
normal–inverse-gamma distribution of §5.2 avoids this.

**The trend in expression.** On real data the variances depend on expression. On GTEx muscle they follow a U-shape:
high for weakly expressed genes, lowest in the upper middle, and higher again for the most expressed genes
([§7.3](validation.md)). A curve that can only decrease, such as the dispersion trend of DESeq2, does not follow it.
The scale ``\sigma_0^2`` is therefore a local fit: at each point of a grid over ``x``, a local-linear likelihood over the 20 %
of genes nearest in expression (`span = 0.2`), with tricube weights. It is held flat beyond the 1 % and 99 % quantiles
of ``x``. The covariate ``x_j`` is the mean of the initial log counts of gene ``j``, fixed during the fit.

`var_prior = :common` gives a single scale, as limma without trend. It pulls the variances of weakly and highly
expressed genes towards the middle.

**Degrees of freedom passed on.** The shrunk variance is ``\sigma_j^2 = 1/\mathbb E[1/\sigma_j^2]``, with ``2a_0 + n -
t_j`` degrees of freedom. The ``t`` likelihood of MOUTHWASH uses them ([§4.4](empirical_bayes.md)).

**Counts without information.** A count of 0 or 1 says little about its latent value: its variational variance
``S^2_{ij}`` stays close to ``\sigma_j^2``, and the count only returns the current value of ``\sigma_j^2``. Such a count
is not an observation of the variance. The shape above therefore also leaves out ``\sum_i S^2_{ij}/\sigma_j^2``, the
share of the counts that the data leave unresolved: ``t_j`` becomes ``u_j = t_j + \sum_i S^2_{ij}/\sigma_j^2``, at
most ``n``. Without this, low counts make the variances look far more certain than they are, ``a_0`` grows without
bound, and the fit drifts.

The share of uninformative counts depends on ``\sigma_j^2`` and ``S^2_{ij}``, which are far from each other at the
start. With low counts the ELBO can therefore decrease slightly in the first iterations; it then increases.

**Weak factors.** Scores and loadings are updated in turn, and for a weak factor a change of one is almost undone by
the other: alternating updates then take hundreds of small steps. After each iteration, a parameter-expanded step
(PX-VB, Qi & Jaakkola 2007) finds the best rescaling and rotation of scores and loadings that leaves their product,
and so the fit, unchanged. It cannot decrease the ELBO, and brings such fits from hundreds of iterations to a few
dozen.

**Slow convergence.** Very sparse counts, such as single-cell data, make the fit slow. When all genes have exactly
the same variance, ``a_0`` grows without bound, which is the right answer, and the trend has nothing to fit: the fit
may then stop at `maxiter` with `converged = false`. Real data are far from this case.

## 5.5 From the PLN fit to MOUTHWASH

The raw estimate of each gene is the rotated row of interest of ``M``: ``b_j = u^\top M_{\cdot j} / R_{kk}``, where
``u`` is the column ``k`` of ``Q``.

Its standard error uses the model's variance of each count on the latent scale: the biological variance plus the
Poisson noise, ``\sigma_j^2 + 1/\bar y_{ij}``:

```math
s_j = \frac{1}{\lvert R_{kk}\rvert}\sqrt{\sum_i u_i^2\, \big(\sigma_j^2 + 1/\bar y_{ij}\big)} .
```

The spread of the residual rows of ``M`` cannot be used instead: the fit pulls them towards the factors, so their
spread underestimates the noise.

**Uncertainty of the loadings.** The fit also returns the posterior covariance of each gene's loadings, `Γcov`.
MOUTHWASH uses it: the correction ``\Gamma_j^\top\alpha`` has the error variance ``\alpha^\top
\text{Cov}(\Gamma_j)\,\alpha``, which is added to the noise of the gene ([§4.3](empirical_bayes.md)). The estimates
``b_j``, their standard errors, the loadings and their covariances then go into MOUTHWASH.

## 5.6 Without the PLN layer: Gaussian log-ratios

[`log_ratio`](@ref) followed by the Gaussian [`mouthwash`](@ref) is a lighter alternative. The counts are turned into
log-ratios, ``\log(y_{ij} + 0.5)`` minus the mean over a set of reference genes (all genes, or the controls), and
analysed as Gaussian data: rotation, factor analysis, variance shrinkage by vash or limma, and MOUTHWASH. The offset
factor plays the same role as in the PLN fit.

It is fast and works well when counts are high. Its limits:

- **Low counts.** The pseudocount and the logarithm distort low counts, and the Poisson noise is not modelled: a
  gene with 3 reads is treated as precisely as one with 3 000.
- **Fixed loadings.** The factor analysis gives point estimates of the loadings. Their uncertainty (§5.5) is not
  passed to MOUTHWASH, so confounded designs with few samples can be slightly anti-conservative.
- **A common variance scale.** The variance trend is off by default (`var_trend = false`), and when on it is a cubic in
  the mean, not the local fit of §5.4.

## Functions

- [`plnpca`](@ref): the PLN fit.
- [`pln_mouthwash`](@ref): the PLN fit followed by MOUTHWASH.
- [`log_ratio`](@ref), [`mouthwash`](@ref): the Gaussian branch.
