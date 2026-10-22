# 2. The model

## 2.1 Counts on a latent Gaussian layer

The model has two layers. The bottom layer is the measurement. The top layer is the biology, on a log scale where
effects add up.

- **assumption** — given the latent log abundance ``Z_{ij}``, the counts are independent Poisson variables:

```math
y_{ij} \mid Z_{ij} \;\sim\; \text{Poisson}\big(e^{\,o_i + Z_{ij}}\big),
```

  where ``o_i`` is a known offset, by default the log library size ``\log \sum_j y_{ij}``.

- **assumption** — the latent log abundances follow a regression on the design, latent factors and independent
  Gaussian noise:

```math
Z_{ij} \;=\; x_i^\top B_j \;+\; w_i^\top \Gamma_j \;+\; \varepsilon_{ij},
\qquad \varepsilon_{ij} \sim \mathcal N(0, \sigma_j^2).
```

Here ``x_i`` is the row of the design matrix ``X`` (``n \times k``, with an intercept), ``B_j`` the coefficients of
gene ``j``, ``w_i`` the ``q`` latent factors of sample ``i``, ``\Gamma_j`` their loadings on gene ``j``, and
``\sigma_j^2`` the *overdispersion* of gene ``j``.

This is the Poisson-lognormal (PLN) model of Aitchison & Ho (1989), with a low-rank covariance as in Chiquet et al.
(2018). Integrating ``Z`` out, a count has mean and variance

```math
\mathbb E[y_{ij}] = \lambda_{ij} = e^{o_i + \mu_{ij} + \sigma_j^2/2},
\qquad
\operatorname{Var}[y_{ij}] = \lambda_{ij} + \lambda_{ij}^2\,\big(e^{\sigma_j^2} - 1\big),
```

with ``\mu_{ij} = x_i^\top B_j + w_i^\top \Gamma_j``. This is the quadratic mean–variance relationship of the negative
binomial. The Poisson part dominates for small counts, and the lognormal part (biological variability) dominates
for large counts. Zeros need no special treatment and no pseudo-count.

!!! details "Why Poisson-lognormal and not negative binomial"
    Both have the same mean–variance relationship. The lognormal layer has two advantages here. Effects, factors and
    noise all add up on the same log scale, so the latent layer is an ordinary Gaussian factor model, and every
    method for confounders (CATE, RUV, MOUTHWASH) applies to it. And it is closed under binomial thinning: thinning
    a Poisson(``\lambda``) count with probability ``\pi`` gives a Poisson(``\pi\lambda``) count. The simulations of
    [§7](simulation.md), which thin real counts, therefore stay exactly inside the model.

## 2.2 The confounders and the offset factor

The latent factors ``w_i`` play the role of the hidden confounders of [§1.5](problem.md), with the part correlated
with the design already absorbed in ``B_j``:

```math
B_j \;=\; \beta_j + \alpha^\top \Gamma_j \qquad\text{(for the covariates of interest)} .
```

So ``B_j`` is the *confounded* effect ``\tilde\beta_j``, and the factors ``w_i`` are, by construction, uncorrelated
with the design (``X^\top W = 0``, see [§3](confounders.md)).

- **assumption** (optional, on by default) — the first factor has loadings fixed to one: ``\Gamma_{1j} = 1`` for
  every gene.

This is the *offset factor* of [§1.6](problem.md). Its value ``w_{i1}`` is a per-sample shift, common to all genes.
It absorbs the difference between the true offset (with the unknown total RNA content) and the library size used in
``o_i``. Its part correlated with the design is the compositional shift ``\delta = \alpha_1``.

## 2.3 What is identified

- **result** — without further information, ``\beta_j`` and ``\alpha`` are not separately identified: only
  ``B_j = \beta_j + \alpha^\top\Gamma_j`` is.
- **assumption** — either a set ``C`` of negative control genes has ``\beta_j = 0`` for the covariate of interest,
  or the distribution of the ``\beta_j`` has a point mass at zero (sparsity).

With either assumption, ``\alpha`` and therefore every ``\beta_j`` are identified, including the compositional
shift. The two options use the same code: a control gene is a gene whose prior is the point mass at zero
([§4.2](empirical_bayes.md)).

## 2.4 Two ways to fit it

The package offers two branches, which share the model and differ in how they handle uncertainty:

| | empirical Bayes ([§4](empirical_bayes.md), [§5](pln.md)) | fully Bayesian ([§6](bayes.md)) |
|---|---|---|
| counts | PLN, variational EM | PLN, Gibbs sampler with Metropolis–Hastings latent updates |
| confounders ``\alpha``, loadings ``\Gamma`` | point estimates | integrated out by MCMC |
| prior on the effects | estimated from the data (ash) | flat (RUVB), then ash on the posterior |
| variances | shrunk across genes (vash or limma) | half-``t`` prior, sampled |
| output | posterior mean, lfsr, s-value | full posterior samples |

The empirical Bayes branch is fast and has no prior to choose. The fully Bayesian branch propagates the uncertainty
of the confounders, which matters when the factors are hard to estimate (few samples, low counts).

## Functions

- [`pln`](@ref): fit the PLN model with a design, latent factors and the offset factor.
- [`pln_mouthwash`](@ref): the complete empirical Bayes branch.
