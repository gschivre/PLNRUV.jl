# 5. Fitting the PLN layer

The latent log-abundances ``Z_{ij}`` of [§2](model.md) are not observed. The likelihood of a count is an integral
with no closed form:

```math
p(y_{ij} \mid \mu_{ij}, \sigma_j^2) = \int \text{Poisson}\big(y_{ij};\, e^{o_i + z}\big)\,
\mathcal N(z;\, \mu_{ij}, \sigma_j^2)\, dz .
```

Following Chiquet et al. (2018), the package maximises a lower bound of the log-likelihood instead, the ELBO, by
variational EM.

## 5.1 The variational approximation

- **Approximation** — the posterior of each latent value is replaced by an independent normal distribution,
  ``\mathcal N(M_{ij}, S^2_{ij})``.

With a normal approximation, every term of the ELBO is explicit, because ``\mathbb E[e^{Z}] = e^{M + S^2/2}``:

```math
\text{ELBO} = \sum_{i,j}\Big[ y_{ij}(o_i + M_{ij}) - e^{\,o_i + M_{ij} + S^2_{ij}/2}
- \frac{(M_{ij} - \mu_{ij})^2 + S^2_{ij}}{2\sigma_j^2} - \tfrac12\log\sigma_j^2
+ \tfrac12\log S^2_{ij}\Big] + \text{const}.
```

The algorithm alternates an E-step (update ``M`` and ``S^2``) and an M-step (update ``\mu`` and ``\sigma^2``).

## 5.2 E-step: a closed form

For one cell ``(i, j)``, the ELBO is concave in ``M_{ij}``. **Result** — its maximum has a closed form:

```math
M_{ij} = \mu_{ij} + \sigma_j^2 y_{ij} - \mathcal W_0\Big(\sigma_j^2\, e^{\,o_i + \mu_{ij} + \sigma_j^2 y_{ij} +
S^2_{ij}/2}\Big),
```

where ``\mathcal W_0`` is the Lambert W function (the solution of ``w e^w = x``). The variance is the solution of
``S^2_{ij} = 1 / (e^{\,o_i + M_{ij} + S^2_{ij}/2} + 1/\sigma_j^2)``, found by fixed-point iteration. Three
alternations of the two updates are enough.

!!! details "Derivation"
    The derivative in ``M`` is zero when ``y - e^{a + M} - (M - \mu)/\sigma^2 = 0``, with ``a = o + S^2/2``. Let
    ``u = \mu + \sigma^2 y - M``. The equation becomes ``u\,e^{u} = \sigma^2 e^{a + \mu + \sigma^2 y}``, so ``u =
    \mathcal W_0(\sigma^2 e^{a + \mu + \sigma^2 y})``. For large counts the argument overflows. The package therefore
    computes ``\mathcal W_0(e^t)`` directly, by Newton's method on ``w + \log w = t``. The derivative in ``S^2`` is
    zero when ``1/S^2 = e^{a + M} + 1/\sigma^2``.

The E-step runs in parallel over genes.

## 5.3 M-step: the Gaussian machinery

With ``M`` and ``S^2`` fixed, the ELBO depends on the model parameters only through

```math
-\sum_{i,j} \frac{(M_{ij} - \mu_{ij})^2 + S^2_{ij}}{2\sigma_j^2} - \frac n2 \sum_j \log\sigma_j^2 ,
\qquad \mu_{ij} = x_i^\top B_j + w_i^\top \Gamma_j .
```

This is a Gaussian factor regression of the matrix ``M``, solved with the rotation of [§3](confounders.md):

- the first ``k`` rows of ``Q^\top M`` give ``B`` exactly;
- a weighted PCA of the residual rows gives the factors ``W`` and the loadings ``\Gamma``, offset factor first.

The variance of gene ``j`` is its residual sum of squares, including the variational variances, divided by an
effective number of degrees of freedom:

```math
\hat\sigma_j^2 = \frac{\sum_i \big[(M_{ij} - \mu_{ij})^2 + S^2_{ij}\big]}{n - (n - d)\,\bar\kappa_j},
\qquad \bar\kappa_j = \frac 1n \sum_i \Big(1 - \frac{S^2_{ij}}{\sigma_j^2}\Big).
```

``\kappa_{ij}`` is the information that cell ``(i, j)`` carries about its latent value: 1 when ``Z_{ij}`` is known
exactly, 0 when the cell carries no information. The divisor is ``d`` for well-measured genes, as in a regression. It
is ``n`` for genes with almost no reads, whose cells cannot pay for the fitted parameters. The variances are then
shrunk across genes ([§4.5](empirical_bayes.md)).

**Convergence.** The E-step increases the ELBO. The M-step with variance shrinkage is not an exact maximisation, so
the ELBO can decrease slightly between iterations. The fit stops when the ELBO changes by less than a relative
tolerance, in either direction. The field `converged` of the result reports whether this happened before `maxiter`.

## 5.4 From the PLN fit to MOUTHWASH

The raw estimate of each gene is the rotated row of interest of ``M``: ``b_j = u^\top M_{\cdot j} / R_{kk}``, where
``u`` is the column ``k`` of ``Q``.

Its standard error cannot come from the residual rows of ``M``. The E-step pulls them towards the fitted factors, so
their spread underestimates the noise. The standard error uses instead the model's variance of each cell on the
latent scale. It is the biological variance plus the Poisson noise, ``\sigma_j^2 + 1/\bar y_{ij}``:

```math
s_j = \frac{1}{\lvert R_{kk}\rvert}\sqrt{\sum_i u_i^2\, \big(\sigma_j^2 + 1/\bar y_{ij}\big)} .
```

The estimates ``b_j``, their standard errors ``s_j`` and the loadings ``\Gamma`` then go into MOUTHWASH
([§4.3](empirical_bayes.md)). The check of these standard errors is in [§7.3](validation.md).

## Functions

- [`pln`](@ref): variational EM for the PLN layer.
- [`pln_mouthwash`](@ref): the PLN fit followed by MOUTHWASH.
