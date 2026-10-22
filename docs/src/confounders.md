# 3. Separating the confounders from the effects

This page works on a Gaussian scale: a matrix ``Y`` (``n`` samples ``\times`` ``p`` genes) of log abundances, or of
the latent means of the PLN fit ([§5](pln.md)). For a gene ``j``:

```math
Y_{\cdot j} = X \beta_j + Z \Gamma_j + E_{\cdot j},
\qquad Z = X\alpha + W,
\qquad E_{\cdot j} \sim \mathcal N(0, \sigma_j^2 I_n),
```

with ``Z`` the ``n \times r`` matrix of hidden factors (the offset factor included, with ``\Gamma_{1j} = 1``).

## 3.1 The rotation

The design is split as ``X = [X_1\; X_2]``: nuisance covariates ``X_1`` (``k_1`` columns, with the intercept) and the
covariate of interest ``X_2`` (``k_2`` columns, usually one). Write the QR decomposition ``X = QR``, with ``Q`` an
``n \times n`` orthogonal matrix and ``R`` upper triangular, and rotate the data: ``\tilde Y = Q^\top Y``. Because
``Q`` is orthogonal, the noise ``Q^\top E_{\cdot j}`` is still ``\mathcal N(0, \sigma_j^2 I_n)``. The rows of
``\tilde Y`` split into three blocks (Wang et al. 2017):

| rows | number | content |
|---|---|---|
| nuisance | ``k_1`` | ``R_{11}``(nuisance coefficients) ``+ \dots``: carries no information on ``\beta_2``, dropped |
| interest | ``k_2`` | ``R_{22}\,(\beta_{2j} + \alpha_2^\top \Gamma_j) + u^\top\Gamma_j + `` noise |
| residual | ``n - k`` | ``\tilde W \Gamma_j + `` noise, with no trace of ``X`` |

**Result** — the residual rows see the confounders but not the effects. A factor analysis on them estimates the
loadings ``\Gamma_j`` and the variances ``\sigma_j^2`` without any risk of absorbing the effects of interest. This is
how Wang et al. avoid the main pitfall of factor-based corrections (factors that "explain away" the treatment).

!!! details "Derivation"
    ``Q^\top X = R``, which is zero below row ``k``. The rows ``k+1, \dots, n`` of ``Q^\top Y_{\cdot j}`` are
    therefore ``Q_{\text{res}}^\top W \Gamma_j + Q_{\text{res}}^\top E_{\cdot j}``: the design has disappeared. Row
    ``k_1 + 1`` (for ``k_2 = 1``) is ``R_{22}(\beta_{2j} + \alpha_2^\top\Gamma_j) + q^\top W \Gamma_j + q^\top
    E_{\cdot j}``, with ``q`` the corresponding column of ``Q``; the nuisance coefficients do not appear because
    ``R`` is upper triangular. The term ``u^\top \Gamma_j = q^\top W \Gamma_j`` is the realisation of the factors in
    that row. It is the same vector ``u`` for every gene, so dividing by ``R_{22}``, it is absorbed into the
    confounder coefficient: ``\alpha \leftarrow \alpha_2 + u / R_{22}``. The ``\alpha`` estimated below is this
    data-specific value, which is what removes the confounding in this data set.

Dividing the row of interest by ``R_{22}`` gives, for each gene, a confounded estimate and its standard error:

```math
\hat{\tilde\beta}_j = \tilde Y_{\text{int}, j} / R_{22}
\;\sim\; \mathcal N\!\big(\beta_j + \alpha^\top\Gamma_j,\; s_j^2\big),
\qquad s_j = \sigma_j / \lvert R_{22}\rvert .
```

## 3.2 Factor analysis on the residual rows

The residual block ``Y_{\text{res}}`` (``(n - k) \times p``) is fitted by a heteroscedastic factor model,
``Y_{\text{res}} \approx L\,\Gamma``, with gene-specific residual variances ``\sigma_j^2``. The fit alternates two
steps, in the spirit of the early-stopping alternation (ESA) of Owen & Wang (2016):

1. with the variances fixed, the rank-``q`` approximation that minimises ``\sum_j \lVert Y_{\text{res},j} - L
   \Gamma_j\rVert^2 / \sigma_j^2`` is a singular value decomposition of ``Y_{\text{res}}\,\text{diag}(\sigma)^{-1}``;
   the offset factor (loadings one) is fitted first, as a weighted mean of each row;
2. with the fit fixed, ``\sigma_j^2`` is the residual sum of squares divided by its ``n - k - q - 1`` degrees of
   freedom.

**Number of factors.** [`choose_q`](@ref) chooses ``q`` by count splitting (Neufeld et al. 2024): every count is
split by binomial thinning, ``Y^{\text{train}} \sim \text{Binomial}(Y, \varepsilon)``, ``Y^{\text{test}} = Y -
Y^{\text{train}}``, two Poisson matrices that are independent given the latent values. The PLN model is fitted on
``Y^{\text{train}}`` for each ``q`` and scored by the predictive log-likelihood of ``Y^{\text{test}}``, where each
count is ``\text{Poisson}(e^{o_i + Z_{ij}})`` with a fresh ``Z_{ij} \sim \mathcal N(\mu_{ij}, \sigma_j^2)``
(Gauss–Hermite quadrature). The two parts share the latent overdispersion, so a low-rank mean that has learned it is
still partly rewarded; scoring the mean ``\mu`` rather than the per-cell latent values limits this optimism. On
PLN simulations with two factors it recovers ``q = 2``. On thinned GTEx muscle it keeps adding factors up to the
largest ``q`` tried (5–6 of 6), consistent with the residual structure left by ``q = 2``, although with 20 samples
more factors make the error rates worse ([§8](validation.md)).

## 3.3 Identifying the confounders

The confounded estimates follow a regression with unknown coefficients ``\alpha`` and known regressors
``\Gamma_j``: ``\hat{\tilde\beta}_j = \beta_j + \Gamma_j^\top\alpha + \text{noise}``. Two assumptions identify
``\alpha``:

- **negative controls** (RUV-4, Gagnon-Bartsch et al. 2013): for ``j \in C``, ``\beta_j = 0``, so a regression of
  ``\hat{\tilde\beta}_C`` on ``\Gamma_C`` estimates ``\alpha``. It needs at least ``q + 1`` controls (one per factor,
  the offset factor included).
- **sparsity** (CATE-RR of Wang et al., MOUTHWASH of Gerard & Stephens 2020): most ``\beta_j`` are zero, so ``\alpha``
  is the value that makes the residuals ``\hat{\tilde\beta}_j - \Gamma_j^\top\alpha`` as concentrated at zero as
  possible. MOUTHWASH makes this precise with a likelihood, [§4.2](empirical_bayes.md).

**Result** — the compositional shift is estimated like any other confounder: it is the coefficient ``\alpha_1`` of
the offset factor, whose regressor is ``\Gamma_{1j} = 1``. With controls, it is the average confounded effect of the
controls, once the other factors are removed. Without controls, it is the mode of the effects, as in LinDA (Zhou et
al. 2022) or ANCOM-BC, but estimated jointly with the other confounders.

## 3.4 Why a fixed offset factor rather than one more free factor

A free factor could, in principle, learn loadings close to one. But its loadings are estimated from the residual
rows, where the compositional shift may be weak or absent: the shift correlated with the design lives in the row of
interest. With the loadings fixed, nothing has to be learned, and the shift is captured exactly. On simulations with
a shift correlated with the group ([§8.3](validation.md)), a free extra factor leaves an RMSE of 0.15–0.26 and a false
discovery proportion of 6–15 %; the fixed offset factor gives 0.10–0.12 and 3–7 %.

## Functions

- [`factor_analysis`](@ref): heteroscedastic factor analysis with the optional offset factor.
- [`mouthwash`](@ref): rotation, factor analysis and the joint estimation of [§4](empirical_bayes.md).
