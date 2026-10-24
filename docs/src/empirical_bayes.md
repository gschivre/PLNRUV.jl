# 4. Empirical Bayes

A Bayesian answer needs a prior on the effects. A prior chosen by hand is hard to defend. With thousands of genes,
the distribution of the effects can instead be *estimated* from the data, and then used as the prior of each gene.
This is empirical Bayes. The package estimates three things this way:

- the distribution of the effects (ash, §4.1–4.2);
- the confounders, together with that distribution (MOUTHWASH, §4.3);
- the distribution of the gene variances (limma or vash, §4.5).

## 4.1 Adaptive shrinkage (ash)

Start from estimates ``b_j`` with standard errors ``s_j``: ``b_j \sim \mathcal N(\beta_j, s_j^2)``.

- **Assumption** — the effects come from a distribution ``g`` that is unimodal and centred at zero. It is written as
  a mixture of normals that includes a point mass at zero (Stephens 2017):

```math
g(\beta) = \sum_{l=0}^{L} \pi_l\, \mathcal N(\beta;\, 0,\, \tau_l^2),\qquad \tau_0 = 0 .
```

The standard deviations ``\tau_l`` form a fixed, dense grid, from ``\min_j s_j / 10`` to ``2\sqrt{\max_j (b_j^2 -
s_j^2)}``. Only the weights ``\pi_l`` are estimated. So ``g`` can take almost any unimodal symmetric shape: spiky,
heavy-tailed, or all its mass at zero.

The weights maximise the marginal likelihood of the estimates, plus a penalty ``9 \log \pi_0`` that favours the
point mass at zero (a Dirichlet(10, 1, …, 1) prior). The penalty makes the error rates conservative:

```math
\hat\pi = \arg\max_{\pi}\; \sum_j \log \sum_l \pi_l\, \mathcal N\big(b_j;\, 0,\, s_j^2 + \tau_l^2\big)
\;+\; 9 \log \pi_0 .
```

This problem is concave. It is solved by mix-SQP (Kim et al. 2020), which converges in a few steps where EM can need
thousands.

!!! details "Posterior of each effect"
    Given ``\hat\pi``, the posterior of ``\beta_j`` is again a mixture. Component ``l`` has weight ``\propto
    \hat\pi_l\,\mathcal N(b_j; 0, s_j^2 + \tau_l^2)``. For ``l \ge 1`` it is a normal with mean ``b_j\, \tau_l^2 /
    (\tau_l^2 + s_j^2)`` and variance ``\tau_l^2 s_j^2 / (\tau_l^2 + s_j^2)``. Component 0 is the point mass at zero.

**Outputs.**

- The **posterior mean** of each effect. A noisy estimate is pulled towards zero more than a precise one.
- The **local false sign rate** (lfsr): the probability that the sign of the posterior mean is wrong,
  ``\text{lfsr}_j = \min\{P(\beta_j \ge 0 \mid \text{data}),\, P(\beta_j \le 0 \mid \text{data})\}``.
- The **s-value**: the average lfsr of all genes at least as significant. Reporting the genes with s-value below
  0.05 keeps the expected proportion of wrong signs among them at 5 %.
- The local false discovery rate ``P(\beta_j = 0 \mid \text{data})`` and its q-value.

The lfsr is preferred to the local false discovery rate. A discovery with the wrong sign is an error too. And the
lfsr does not depend much on the exact mass at zero, which is hard to estimate.

## 4.2 A prior for global changes

The unimodal prior of §4.1 assumes that most genes do not change. Suppose instead that the controls are valid and
that most other genes go up together. The effects are then a point mass at zero (the unchanged genes) plus a bump
around a common shift ``\bar\beta`` ([figure of §1.3](problem.md)). A prior centred at zero can reach the bump only
with wide components. These also cover small values, so they take mass from the point mass at zero. The null genes
then look less null than they are, and false discoveries increase.

With `mode = :auto`, the non-null part of the prior is centred at an estimated ``\bar\beta``. The point mass stays at
zero, where the controls put it:

```math
g(\beta) = \pi_0\, \delta_0 + \sum_{l=1}^{L} \pi_l\, \mathcal N(\beta;\, \bar\beta,\, \tau_l^2) .
```

ash in R has a similar option, `mode = "estimate"`, but it moves the point mass to ``\bar\beta`` too. That makes the
bulk of changing genes the reference, which is what the controls are meant to avoid.

``\bar\beta`` maximises the marginal likelihood. It is used only when it improves the fit of genes it was not fitted
on. The genes are split in two halves. Both priors are fitted on one half and scored on the other. The shifted prior
is kept only if it gains more than `mode_threshold = 0.05` in log-likelihood per gene. Otherwise the prior is the
unimodal one of §4.1.

**This option needs controls.** Without controls, a shift shared by most genes cannot be told apart from the
compositional offset ``\alpha_0``. Sparsity then sets the most common change to zero, by definition.

## 4.3 MOUTHWASH: confounders and effects together

After the rotation of [§3](confounders.md), each gene has a raw estimate ``b_j \sim \mathcal N(\beta_j +
\Gamma_j^\top\alpha,\; \lambda s_j^2)``. Gerard & Stephens (2020) estimate ``\alpha`` inside the ash likelihood:

```math
(\hat\pi, \hat\alpha, \hat\lambda) = \arg\max\;
\sum_{j \notin C}\log\sum_l \pi_l\,\mathcal N\big(b_j - \Gamma_j^\top\alpha;\; 0,\; \lambda s_j^2 + \tau_l^2\big)
+ \sum_{j \in C}\log\mathcal N\big(b_j - \Gamma_j^\top\alpha;\; 0,\; \lambda s_j^2\big)
+ 9\log\pi_0 .
```

The second sum is over the controls, whose prior is the point mass at zero. Without controls, only the point mass of
``g`` identifies ``\alpha`` (sparsity). With controls, ``\alpha`` is fitted on the controls only, by weighted least
squares. This stays valid when most other genes change.

The **variance inflation** ``\lambda`` absorbs what the standard errors miss, mainly the error in the estimated
loadings ``\Gamma_j``.

The maximisation alternates three steps:

1. ``\pi`` by mix-SQP, with ``\alpha`` and ``\lambda`` fixed;
2. ``\alpha`` by weighted least squares (without controls);
3. ``\lambda`` by a one-dimensional search.

It stops when the objective changes by less than a relative tolerance, in either direction. The field `converged`
of the result reports whether this happened before the iteration limit.

The effects are then shrunk by ash, applied to the corrected estimates ``b_j - \Gamma_j^\top\hat\alpha``.

**Uncertainty of ``\hat\alpha``.** The corrected estimates share the error of ``\hat\alpha``. Its variance is added
to each gene: ``s_j^2 + \Gamma_j^\top V_\alpha \Gamma_j``, where ``V_\alpha`` is the covariance of ``\hat\alpha``
(as in Wang et al. 2017). The grid of ``\tau_l`` is built from these standard errors.

## 4.4 The noise model

**The ``t`` likelihood.** Each ``s_j`` rests on a variance estimated from ``d`` degrees of freedom (§4.5). Given the
prior on the variances, the shrunk estimate follows a ``t`` distribution with ``d_0 + d`` degrees of freedom (Lemma 1
of Lu & Stephens 2019). With `likelihood = :t`, the default, the normal densities above are replaced by these ``t``
densities.

**Inflation adds noise.** The inflation ``\lambda`` stands for noise that the variance estimate does not see. Scaling
the ``t`` by ``\lambda`` would stretch its tails over this extra noise. Instead, the package adds normal noise of
variance ``(\lambda - 1) s_j^2`` to the ``t``. The sum is approximated by a ``t`` with squared scale ``\lambda
s_j^2`` and the same excess kurtosis, which gives the degrees of freedom

```math
\nu_\lambda = 4 + (\nu - 4)\,\lambda^2 \qquad (\lambda > 1;\ \text{otherwise } \nu_\lambda = \nu).
```

!!! details "Derivation"
    A ``t_\nu`` has excess kurtosis ``6 / (\nu - 4)``. Adding independent normal noise does not add kurtosis, but
    dilutes it: the ``t`` part carries a share ``1/\lambda`` of the variance, so the excess kurtosis of the sum is
    ``6 / ((\nu - 4)\lambda^2)``. Setting this equal to ``6 / (\nu_\lambda - 4)`` gives the formula.

**Mixture of ``t``.** With vash (§4.5) the prior on the variances is a mixture. The exact likelihood is then a mixture
of ``t``, one per variance component ``m``, with ``2a_m + d`` degrees of freedom. `likelihood = :tmix` uses this
mixture. `likelihood = :t` uses a single ``t`` with averaged degrees of freedom. `likelihood = :normal` ignores the
uncertainty of the variances.

## 4.5 Variance shrinkage

Each variance ``\hat\sigma_j^2`` is estimated from few degrees of freedom: ``d = n - k - q - 1``, about six for ten
samples. A gene with a small sample variance by chance gets a large, spurious significance. Pooling the variances
across genes prevents this.

- **Assumption** — ``\hat\sigma_j^2 \mid \sigma_j^2 \sim \sigma_j^2\,\chi^2_d / d``.

**limma** (Smyth 2004) puts one scaled inverse-``\chi^2`` prior on the variances, ``\sigma_j^2 \sim d_0 \sigma_0^2 /
\chi^2_{d_0}``. It estimates ``d_0`` and ``\sigma_0^2`` from the moments of ``\log \hat\sigma_j^2``. The shrunk
variance is

```math
\tilde\sigma_j^2 = \frac{d_0 \sigma_0^2 + d\, \hat\sigma_j^2}{d_0 + d} .
```

**vash** (Lu & Stephens 2016) does for the variances what ash does for the effects. Its prior is a mixture of inverse
gammas, with a common scale ``c`` and a grid of shapes ``a_m``:

```math
\sigma_j^2 \sim \sum_m \omega_m\, \text{InvGamma}(a_m,\; a_m c) .
```

limma is the case of one component. A mixture can fit variances that are bimodal or heavy-tailed. The weights
``\omega_m`` are estimated by mix-SQP, and ``c`` by a one-dimensional search. The shrunk variance is ``\tilde\sigma_j^2
= 1 / \mathbb E[\sigma_j^{-2} \mid \hat\sigma_j^2]``.

The option `shrink_var = :vash` (default), `:limma` or `:none` applies wherever a variance is estimated: in the
Gaussian MOUTHWASH and in each M-step of the PLN fit.

**Variances that depend on expression.** Both priors pull every variance towards the same scale. But RNA-seq
overdispersions decrease with expression: this is the mean–dispersion trend fitted by DESeq2 and edgeR. A common
scale then shrinks the variances of low-count genes too much and those of high-count genes too little. With
`var_trend = true`, the variances are shrunk around a cubic trend of ``\log \hat\sigma_j^2`` in the mean
log-abundance of each gene (as in limma-trend, Law et al. 2014). This option is off by default: on real data it made
calibration worse ([§7.5](validation.md)).

## Functions

- [`ash`](@ref): adaptive shrinkage of estimates with standard errors.
- [`mouthwash`](@ref): rotation, factor analysis and the joint estimation of ``g``, ``\alpha`` and ``\lambda``.
- [`squeeze_var`](@ref), [`vash`](@ref): variance shrinkage.
- [`evalues`](@ref), [`ebh`](@ref): compound e-values and the e-BH procedure.
