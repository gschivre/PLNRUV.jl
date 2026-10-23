# 4. Empirical Bayes: ash, MOUTHWASH and variance shrinkage

A Bayesian answer needs a prior on the effects. Biologists are rightly wary of priors chosen by hand. Empirical Bayes
removes the choice: with thousands of genes, the distribution of the effects can be *estimated* from the data, and
then used as the prior of each gene. Three quantities are estimated this way: the distribution of the effects
(ash), the confounders (MOUTHWASH), and the distribution of the variances (vash or limma).

## 4.1 Adaptive shrinkage (ash)

Start from estimates ``\hat\beta_j`` with standard errors ``s_j``: ``\hat\beta_j \sim \mathcal N(\beta_j, s_j^2)``.

- **assumption** — the effects are drawn from a unimodal distribution centred at zero, written as a mixture of
  normals with a point mass at zero (Stephens 2017):

```math
g(\beta) \;=\; \sum_{k=0}^{K} \pi_k\, \mathcal N(\beta;\, 0,\, \sigma_k^2),\qquad \sigma_0 = 0 .
```

The grid ``\sigma_1 < \dots < \sigma_K`` is fixed and dense: from ``\min_j s_j / 10`` to ``2\sqrt{\max_j
(\hat\beta_j^2 - s_j^2)}``, in steps of ``\sqrt 2``. Only the weights ``\pi`` are estimated, so ``g`` can take almost
any unimodal symmetric shape: spiky, heavy-tailed, or all mass at zero.

The weights maximise the marginal likelihood, with a Dirichlet(10, 1, …, 1) penalty that favours the null component
(the "nullweight"; it makes the false sign rate conservative):

```math
\hat\pi = \arg\max_{\pi \in \Delta}\; \sum_j \log \sum_k \pi_k\, \mathcal N\!\big(\hat\beta_j;\, 0,\, s_j^2 +
\sigma_k^2\big) \;+\; 9 \log \pi_0 .
```

This is a concave problem over the simplex. It is solved by mix-SQP (Kim et al. 2020): sequential quadratic
programming, where each step solves a small non-negative quadratic program by an active-set method. It converges in
about ten steps, where EM can need thousands when the mixture components overlap.

!!! details "Posterior of each effect"
    Given ``\hat\pi``, the posterior of ``\beta_j`` is again a mixture. Component ``k`` has weight
    ``w_{jk} \propto \hat\pi_k\,\mathcal N(\hat\beta_j; 0, s_j^2 + \sigma_k^2)`` and, for ``k \ge 1``, is
    ``\mathcal N(m_{jk}, v_{jk})`` with ``v_{jk} = \sigma_k^2 s_j^2/(\sigma_k^2 + s_j^2)`` and
    ``m_{jk} = \hat\beta_j\, \sigma_k^2/(\sigma_k^2 + s_j^2)``; component 0 is the point mass at zero. The posterior
    mean is ``\sum_k w_{jk} m_{jk}``.

**Outputs.**

- the *posterior mean*, a shrunken estimate: a noisy estimate is pulled towards zero more than a precise one;
- the **local false sign rate** (lfsr), the probability that we get the sign of ``\beta_j`` wrong if we report the
  sign of its posterior mean:

```math
\text{lfsr}_j = \min\big\{P(\beta_j \ge 0 \mid \text{data}),\; P(\beta_j \le 0 \mid \text{data})\big\} ;
```

- the *s-value*, the average lfsr of all genes at least as significant. Declaring the genes with s-value below 0.05
  controls the expected proportion of wrong signs among the discoveries at 5 %;
- the local false discovery rate ``P(\beta_j = 0 \mid \text{data})`` and its q-value.

The lfsr is preferred to the lfdr: a "discovery" of the wrong sign is an error too, and the lfsr does not depend on
the exact mass at zero, which is hard to estimate.

**A global change of transcription breaks unimodality at zero.** When most genes change together, as in a global
increase of transcription measured against controls, the effects are a point mass at zero (the unchanged genes)
plus a bump away from it, around a shared shift ``\mu``. A prior unimodal at zero can only reach the bump with wide
components centred at zero, which also cover small values: they take mass from the null, ``\pi_0`` is
underestimated, and the lfsr of the null genes is too small. On thinned GTEx with 60 % of the genes shifted by a
log2 fold change of about 1, ash finds ``\hat\pi_0 \approx 0.22`` instead of 0.40 and the FDP at ``s < 0.05`` is
about 9 %. (The CATE-based pipeline of Schivre et al. 2025 meets the same issue.) With `mode = :auto` the non-null
part of the prior is centred at an estimated mode while the null stays at zero, where the controls put it:

```math
g(\beta) \;=\; \pi_0\, \delta_0 \;+\; \sum_{k=0}^{K} \pi_k\, \mathcal N(\beta;\, \mu,\, \sigma_k^2) .
```

This is ash's `mode = "estimate"` with one difference: ash moves the point mass to the mode too, which turns the bulk
of changing genes into the reference, the assumption the controls are there to avoid. ``\mu`` maximises the
marginal likelihood (golden-section search). It is used only when it improves the fit of genes it was not fitted
on: the genes are split in two halves, both priors are fitted on one half and scored on the other, and the mode is
kept when the held-out log-likelihood gain exceeds `mode_threshold = 0.05` per gene. Otherwise the prior is plain
ash. On thinned GTEx the gain is at most 0.003 per gene without a shared shift and at most 0.016 with residual
confounding (pure nulls, the controls, give ± 0.01), against 0.24–0.33 for a shift of 1 in 60 % of the genes and
0.06–0.14 for a shift of 0.5 in 40 %. A shift too small to pass the gate falls back to ash. In MOUTHWASH the
mode is estimated only with controls: without them, a shift shared by the genes cannot be told apart from the
compositional offset.

## 4.2 MOUTHWASH: confounders and effects together

After the rotation of [§3](confounders.md), each gene gives a confounded estimate
``\hat{\tilde\beta}_j \sim \mathcal N(\beta_j + \Gamma_j^\top\alpha,\; \lambda s_j^2)``. Gerard & Stephens (2020)
estimate ``\alpha`` *inside* the ash likelihood:

```math
(\hat\pi, \hat\alpha, \hat\lambda) = \arg\max\;
\sum_{j \notin C}\log\sum_k \pi_k\,\mathcal N\!\big(\hat{\tilde\beta}_j - \Gamma_j^\top\alpha;\; 0,\;
\lambda s_j^2 + \sigma_k^2\big)
\;+\; \sum_{j \in C}\log\mathcal N\!\big(\hat{\tilde\beta}_j - \Gamma_j^\top\alpha;\; 0,\; \lambda s_j^2\big)
\;+\; 9\log\pi_0 .
```

The second sum is over the control genes, whose prior is the point mass at zero. With no controls, sparsity alone
(the point mass of ``g``) identifies ``\alpha``. The **variance inflation** ``\lambda`` absorbs what the standard
errors miss: mainly the uncertainty in the estimated loadings ``\Gamma_j`` and in ``\alpha``.

The maximisation is a block coordinate ascent; each block increases the objective:

1. ``\pi`` by mix-SQP, with ``\alpha`` and ``\lambda`` fixed;
2. ``\alpha`` by weighted least squares, the weights being the posterior precisions ``c_j = \sum_k w_{jk} / (\lambda
   s_j^2 + \sigma_k^2)`` (a minorise–maximise step);
3. ``\lambda`` by a one-dimensional golden-section search.

The effects are then shrunk by ash on the corrected estimates ``\hat{\tilde\beta}_j - \Gamma_j^\top\hat\alpha``, with
standard errors ``\sqrt{\hat\lambda}\, s_j``.

**Inflation and the ``t`` likelihood.** With `likelihood = :t` (the default) the noise of each estimate is a ``t``
whose degrees of freedom are those of the shrunk variance: given the prior of the variances, the moderated estimate
follows exactly a ``t_{d_0 + d}`` (Lemma 1 of Lu & Stephens 2019). The inflation ``\lambda`` is a different source of
noise: the variance estimate does not see it (loadings, sample-specific variances, …), and on GTEx its residuals have
no excess kurtosis. Multiplying the ``t`` by ``\lambda`` would stretch the ``t`` tails over it. It is added instead,
``t_\nu(0, v) * \mathcal N(0, (\lambda - 1) v)``, approximated by a ``t`` of squared scale ``\lambda v`` whose excess
kurtosis matches: the ``t`` share ``1/\lambda`` of the variance scales ``6 / (\nu - 4)`` by ``1/\lambda^2``, so

```math
\nu_\lambda = 4 + (\nu - 4)\,\lambda^2 \qquad (\lambda > 1;\ \nu_\lambda = \nu \text{ otherwise}).
```

On thinned GTEx (10 replicates) this brings the confounded scenario from 2.2 % to 5.0 % FDP at ``s < 0.05`` with
15 % more power, and leaves the balanced one at 4.5 %.

With vash the predictive is a mixture of ``t``: given the shape ``a_m`` the estimate follows a ``t_{2a_m + d}`` with
the posterior scale of that component. `likelihood = :tmix` uses this mixture,
``\sum_m w_{jm}\, t_{2a_m + d}(\beta_j, \tilde s_{jm})``, instead of a single ``t`` with the averaged degrees of freedom
(`:t`). The two differ when vash keeps a component of nearly fixed variance (``a_m`` large) next to a heavy-tailed
one: averaging the degrees of freedom then gives an almost normal likelihood.

## 4.3 Variance shrinkage

Each ``s_j`` rests on a variance estimated from few degrees of freedom: ``d = n - k - q - 1`` residual rows, about
six for a 10-sample experiment. A gene that happens to have a small sample variance gets a large, spurious
significance. Pooling the variances across genes fixes this.

- **assumption** — the sample variances follow ``s_j^2 \mid \sigma_j^2 \sim \sigma_j^2\,\chi^2_d / d``.

**limma** (Smyth 2004) puts one scaled inverse-``\chi^2`` prior on the variances, ``\sigma_j^2 \sim d_0 s_0^2 /
\chi^2_{d_0}``. It estimates ``d_0`` and ``s_0^2`` from the moments of ``\log s_j^2``, and returns the posterior
estimates

```math
\tilde s_j^2 = \frac{d_0 s_0^2 + d\, s_j^2}{d_0 + d} .
```

**vash** (Lu & Stephens 2016) does for the variances what ash does for the effects. The prior is a mixture of inverse
gammas with a common scale ``c`` and a grid of shapes ``a_k``:

```math
\sigma_j^2 \sim \sum_k \pi_k\, \text{InvGamma}(a_k,\; a_k c),
\qquad\text{so that}\qquad
s_j^2 / c \;\sim\; \sum_k \pi_k\, F(d,\, 2a_k) .
```

limma is the one-component case. The mixture can fit variances that are bimodal or heavy-tailed, as overdispersions
often are. The weights are estimated by mix-SQP, as in ash, and ``c`` by a golden-section search. The posterior
variance is ``\tilde s_j^2 = 1 / \mathbb E[\sigma_j^{-2} \mid s_j^2]``, a mixture of gamma means.

Both are available everywhere a variance is estimated (`shrink_var = :vash`, `:limma` or `:none`): in the Gaussian
MOUTHWASH, and inside each M-step of the PLN fit, where the shrunk overdispersions also stabilise the E-step of
low-count genes. vash is the default. On bimodal variances it beats limma, and on the low-count PLN simulations of
[§7.3](validation.md) it gives the best-calibrated false sign rate.

Both priors are exchangeable: every variance is pulled toward the same scale. RNA-seq overdispersions are not
exchangeable, they decrease with expression (the mean–dispersion trend that DESeq2 and edgeR fit). On thinned GTEx
muscle the PLN overdispersion falls from about 0.45 for genes below 40 reads to about 0.06 above 2500, a slope of about
``-0.5`` on the log scale, so it is not extra Poisson noise (which would give ``-1``). A common scale then shrinks the
overdispersions of low-count genes down and those of high-count genes up, and the standard errors are too small at low
counts and too large at high counts, a bias that the single inflation ``\lambda`` of MOUTHWASH cannot remove. With
`var_trend = true` the variances are shrunk around a cubic trend of ``\log s_j^2`` in the mean log-abundance of each
gene (limma-trend, Law et al. 2014): vash or limma is applied to ``s_j^2 / t_j`` and the result is multiplied by
``t_j``.

## Functions

- [`ash`](@ref): adaptive shrinkage of estimates with standard errors.
- [`mouthwash`](@ref): rotation, factor analysis, and the joint estimation of ``g``, ``\alpha`` and ``\lambda``.
- [`squeeze_var`](@ref), [`vash`](@ref): variance shrinkage.
