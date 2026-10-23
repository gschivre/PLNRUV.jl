# 5. Fitting the Poisson-lognormal layer

The latent log abundances ``Z_{ij}`` of [§2](model.md) are not observed. The likelihood of the counts requires the
integral

```math
p(y_{ij} \mid \mu_{ij}, \sigma_j^2) = \int \text{Poisson}\big(y_{ij};\, e^{o_i + z}\big)\,
\mathcal N(z;\, \mu_{ij}, \sigma_j^2)\, dz ,
```

which has no closed form. Following Chiquet et al. (2018), the package maximises instead a lower bound of the
log-likelihood, the ELBO, by variational EM.

## 5.1 The variational approximation

- **approximation** — the posterior of each latent value is replaced by an independent Gaussian:
  ``q(Z_{ij}) = \mathcal N(M_{ij}, S^2_{ij})``.

With a Gaussian ``q``, every term of the ELBO is explicit, because ``\mathbb E_q[e^{Z}] = e^{M + S^2/2}``:

```math
\text{ELBO} = \sum_{i,j}\Big[ y_{ij}(o_i + M_{ij}) - e^{\,o_i + M_{ij} + S^2_{ij}/2}
- \frac{(M_{ij} - \mu_{ij})^2 + S^2_{ij}}{2\sigma_j^2} - \tfrac12\log\sigma_j^2
+ \tfrac12\log S^2_{ij}\Big] + \text{const}.
```

The algorithm alternates two steps, each of which increases the ELBO.

## 5.2 E-step: a closed form with the Lambert W function

For one cell, the ELBO is concave in ``M`` and in ``S^2``. **Result** — its maximum in ``M`` has a closed form:

```math
M_{ij} = \mu_{ij} + \sigma_j^2 y_{ij} - W_0\!\Big(\sigma_j^2\, e^{\,o_i + \mu_{ij} + \sigma_j^2 y_{ij} +
S^2_{ij}/2}\Big),
```

where ``W_0`` is the principal branch of the Lambert W function (``W e^W = x``). The variance solves the fixed point
``S^2_{ij} = 1 / (e^{\,o_i + M_{ij} + S^2_{ij}/2} + 1/\sigma_j^2)``, a contraction. Three alternations of the two
updates are enough.

!!! details "Derivation"
    Setting the derivative in ``M`` to zero: ``y - e^{a + M} - (M - \mu)/\sigma^2 = 0`` with ``a = o + S^2/2``. Write
    ``u = \mu + \sigma^2 y - M``. The equation becomes ``e^{a + \mu + \sigma^2 y - u} = u/\sigma^2``, that is
    ``u\,e^{u} = \sigma^2 e^{a + \mu + \sigma^2 y}``, so ``u = W_0(\sigma^2 e^{a + \mu + \sigma^2 y})``. For large
    counts the argument overflows, so ``W_0(e^t)`` is computed directly by Newton's method on ``w + \log w = t``,
    from ``w = t - \log t``. The derivative in ``S^2`` gives ``-\tfrac12 e^{a+M} - 1/(2\sigma^2) + 1/(2S^2) = 0``.

The E-step is a loop over the ``n \times p`` cells, run in parallel over genes.

## 5.3 M-step: the Gaussian machinery

With ``q`` fixed, the ELBO depends on the model parameters only through

```math
-\sum_{i,j} \frac{(M_{ij} - x_i^\top B_j - w_i^\top\Gamma_j)^2 + S^2_{ij}}{2\sigma_j^2} - \frac n2 \sum_j
\log\sigma_j^2 .
```

This is a Gaussian factor regression of the matrix ``M`` on the design, solved with the rotation of
[§3](confounders.md). The ``k`` design rows of ``Q^\top M`` give ``B = R^{-1}(Q^\top M)_{1:k}`` exactly. A weighted
PCA of the ``n - k`` residual rows gives the factors ``W`` (orthogonal to ``X``) and the loadings ``\Gamma``, with the
offset factor first. The overdispersions are the residual sums of squares, plus the variational variances:

```math
\sigma_j^2 = \frac{1}{n - k - q - 1}\sum_i \Big[(M_{ij} - \mu_{ij})^2 + S^2_{ij}\Big],
```

shrunk across genes by vash or limma ([§4.3](empirical_bayes.md)).

## 5.4 From the PLN fit to MOUTHWASH

The confounded estimate of each gene is the rotated row of interest of ``M``: ``\hat{\tilde\beta}_j = q^\top
M_{\cdot j} / R_{22}``, with ``q`` the column of ``Q`` for the covariate of interest. Its standard error needs care.

**Result** — the residual rows of ``M`` cannot provide it: the E-step pulls them towards the fitted factors, so
their spread underestimates the noise. On simulations this gave a false discovery proportion of 25 % at a nominal 5 %.

Instead, the standard error uses the model's own variance of each cell on the latent scale: the lognormal
overdispersion plus the Poisson noise, ``v_{ij} = \sigma_j^2 + 1/\lambda_{ij}``, where ``\lambda_{ij} = e^{o_i +
\mu_{ij} + \sigma_j^2/2}`` is the expected count:

```math
s_j = \frac{1}{\lvert R_{22}\rvert}\sqrt{\sum_i q_i^2\, v_{ij}} .
```

Then ``\hat{\tilde\beta}_j``, ``s_j`` and the PLN loadings ``\Gamma`` go into MOUTHWASH ([§4.2](empirical_bayes.md)).
Without confounding, this standard error needs no inflation (``\hat\lambda \approx 1``) and 95 % intervals cover
95–96 % of the true effects ([§7.3](validation.md)).

## Functions

- [`pln`](@ref): variational EM for the PLN model with design, factors and offset factor.
- [`pln_mouthwash`](@ref): PLN fit followed by MOUTHWASH.
