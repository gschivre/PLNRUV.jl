# 5. Fitting the PLN layer

The latent log-abundances ``Z_{ij}`` of [§2](model.md) are not observed. The likelihood of a count is an integral
with no closed form:

```math
p(y_{ij} \mid \mu_{ij}, \sigma_j^2) = \int \text{Poisson}\big(y_{ij};\, e^{o_i + z}\big)\,
\mathcal N(z;\, \mu_{ij}, \sigma_j^2)\, dz .
```

The package maximises instead a lower bound of the log-likelihood, the ELBO, by variational inference (Chiquet et al.
2018). The default fit, [`plnpca`](@ref), also treats the factors, the loadings, the coefficients and the variances
as random, with priors estimated from all genes. This is variational empirical Bayes.

## 5.1 The model as a hierarchy

[`plnpca`](@ref) completes the model of [§2](model.md) with priors:

- **Factors.** The scores of the ``q`` hidden factors are random, ``\mathcal N(0, 1)``, and live in the residual space
  of the design (``X^\top W = 0``). Their loadings are learned only from the variation that the design cannot explain,
  as in the rotation of [§3](confounders.md). A factor therefore cannot absorb the effects of interest.
- **Loadings.** ``\Gamma_{lj} \sim \mathcal N(0, \tau_l^2)``, with one variance per factor. A factor that the data do
  not support gets ``\tau_l^2 \to 0`` and is switched off (automatic relevance determination). Too many factors then
  cost little.
- **Coefficients.** ``B_j`` has a weak normal prior. It keeps the estimates finite for genes that are zero in a whole
  group.
- **Variances.** ``\sigma_j^2 \sim \text{InvGamma}(a_0, b_0)``, the prior of limma ([§4.5](empirical_bayes.md)).

The offset factor ``w_{i0}`` stays a free parameter, also in the residual space of the design. The hyperparameters
``\tau_l^2``, ``a_0`` and ``b_0`` are estimated from all genes at every iteration.

## 5.2 The variational approximation

- **Approximation** — the posterior is replaced by a product of simple distributions: a normal
  ``\mathcal N(M_{ij}, S^2_{ij})`` for each latent value, normals for the scores, the loadings and the coefficients,
  and an inverse gamma for each variance.

Only the latent values ``Z_{ij}`` see the Poisson counts. Given them, the rest of the model is a linear Gaussian
factor model. So every update has a closed form, and each one increases the ELBO. The fit stops when the ELBO changes
by less than a tolerance per cell; the field `converged` reports whether this happened before `maxiter`.

## 5.3 The latent values: a closed form

For one cell ``(i, j)``, the ELBO is concave in ``M_{ij}``. **Result** — its maximum has a closed form:

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

The update of ``\sigma_j^2`` uses the expected squared residual of each cell. This expectation includes the
uncertainty of the coefficients, the scores and the loadings. **Result** — in the linear Gaussian model, this gives
the restricted (REML) estimate: the residual sum of squares divided by about ``n - k - q`` rather than ``n``. Without
it, the variances would be too small with few samples.

The estimate passed on is ``\sigma_j^2 = 1/\mathbb E[1/\sigma_j^2]``. Its degrees of freedom are ``2a_0 + n - k - q -
1``, as for a variance shrunk by limma with ``d_0 = 2a_0`` prior degrees of freedom.

## 5.5 From the PLN fit to MOUTHWASH

The raw estimate of each gene is the rotated row of interest of ``M``: ``b_j = u^\top M_{\cdot j} / R_{kk}``, where
``u`` is the column ``k`` of ``Q``.

Its standard error uses the model's variance of each cell on the latent scale: the biological variance plus the
Poisson noise, ``\sigma_j^2 + 1/\bar y_{ij}``:

```math
s_j = \frac{1}{\lvert R_{kk}\rvert}\sqrt{\sum_i u_i^2\, \big(\sigma_j^2 + 1/\bar y_{ij}\big)} .
```

The spread of the residual rows of ``M`` cannot be used instead: the fit pulls them towards the factors, so their
spread underestimates the noise. The estimates ``b_j``, their standard errors and the loadings then go into
MOUTHWASH ([§4.3](empirical_bayes.md)).

## 5.6 The older fit

`fit = :pln` selects the earlier fit, [`pln`](@ref). It treats the scores and loadings as fixed parameters, estimated
by a weighted PCA of the residual rows of ``M``, and shrinks the variances with vash or limma at each iteration. With
two or three samples per group it can diverge, and its false discovery rate is far above the nominal level on real
data. The comparison is in [§7.5](validation.md).

**Slow convergence.** Two situations make the fit slow:

- very sparse counts, such as single-cell data, where each cell carries little information;
- low counts together with variances that are nearly the same for all genes, where the prior pools them completely.

The variances then approach their value slowly from above, and the fit may stop at `maxiter` with `converged =
false`. The variances are then overestimated: the error is conservative. On bulk RNA-seq data the variances differ
between genes, and the fit converges in a few tens of iterations.

## Functions

- [`plnpca`](@ref): the default fit of the PLN layer.
- [`pln`](@ref): the older fit.
- [`pln_mouthwash`](@ref): the PLN fit followed by MOUTHWASH.
