# 6. The fully Bayesian branch (RUVB)

MOUTHWASH plugs in point estimates of the loadings ``\Gamma`` and of the confounders ``\alpha``. When they are poorly
estimated (few samples, weak factors, low counts), that uncertainty is lost, and the single inflation factor
``\lambda`` can only partly make up for it. The fully Bayesian branch keeps it. It follows RUVB of Gerard & Stephens
(2021), on a Gaussian scale such as the log-ratios of [`log_ratio`](@ref).

## 6.1 The effects as missing data

After the rotation of [§3](confounders.md), keep the ``k_2`` rows of interest and the ``n - k`` residual rows
(``m = k_2 + n - k`` rows in all), and fit one factor model to all of them, ``\tilde Y = L F + E``.

- **result** — if the non-control genes had no effect, their rows of interest would follow the same factor model as
  the controls. So treat those cells as **missing**: the factor model fitted to everything else predicts what they
  would be without an effect.
- The effect is the difference between what is observed and what is predicted:

```math
\beta_{2j} = R_{22}^{-1}\big(\tilde Y^{\text{obs}}_{\text{int},j} - \tilde Y^{\text{pred}}_{\text{int},j}\big).
```

Each draw of the predicted cells from their posterior predictive distribution gives a draw of ``\beta_{2j}``. The
posterior of the effects thus includes the uncertainty of the loadings, the factors and the variances. The control
genes are what make the prediction possible: their rows of interest are observed, and they anchor the factor values
``L_{\text{int}}`` (at least ``q + 1`` controls are needed, the offset factor included).

## 6.2 Priors

The priors are the half-``t`` priors of Huang & Wand (2013), written as gamma mixtures so that every conditional
distribution is conjugate:

```math
F_{\cdot j} \sim \mathcal N(0, \sigma_j^2 I_q),\quad
\sigma_j^{-2} \sim \text{Gamma}(\nu/2,\ \text{rate } \nu\phi),\quad
\phi \sim \text{Gamma}(1/2,\ \text{rate } 1/a_2),
```

```math
L_{\cdot l} \sim \mathcal N(0, \psi_l I_m),\quad
\psi_l^{-1} \sim \text{Gamma}(\tau/2,\ \text{rate } \tau\eta_l),\quad
\eta_l \sim \text{Gamma}(1/2,\ \text{rate } 1/b_2).
```

The residual standard deviations and the scales of the factors thus have half-``t`` priors with ``\nu`` (resp.
``\tau``) degrees of freedom and scale ``\sqrt{a_2}`` (resp. ``\sqrt{b_2}``). These are weakly informative. With
thousands of genes, the shared hyperparameters ``\phi`` and ``\eta_l`` are learned from the data, and the results are
insensitive to ``a_2``, ``b_2``. With the offset factor, the first row of ``F`` is fixed to ones.

## 6.3 The Gibbs sampler

One sweep updates, in turn:

1. the hyperparameters ``\eta_l``, ``\psi_l``, ``\phi`` (gamma conditionals);
2. the precisions ``\sigma_j^{-2}`` (gamma, with ``m + q`` degrees of freedom from the residuals and the factors);
3. the rows of ``L``: ``\mathcal N(P^{-1}F\Sigma^{-1}\tilde y_i,\; P^{-1})`` with ``P = F\Sigma^{-1}F^\top +
   \text{diag}(\psi^{-1})``;
4. the free rows of ``F``: ``\mathcal N(G^{-1}L_1^\top r_j,\; \sigma_j^2 G^{-1})`` with ``G = L_1^\top L_1 + I``;
5. the missing cells, from their posterior predictive ``\mathcal N(L_i^\top F_{\cdot j}, \sigma_j^2)``, followed by the
   effect draw of §6.1.

Chains run in parallel, each with its own random number generator and an overdispersed start (truncated SVD of the
residual rows plus noise). [`summarize`](@ref) returns the posterior mean, sd, credible interval and lfsr of each
effect.

## 6.4 Status

The Gaussian RUVB is complete and validated ([§8.2](validation.md)). The PLN layer for this branch is the next step:
the latent ``Z_{ij}`` are added to the Gibbs sampler, updated by an independence Metropolis–Hastings step whose
proposal is centred at the Lambert W mode of [§5.2](pln.md). The empirical Bayes fit will provide the defaults of the
hyperparameters and the starting values.

## Functions

- [`ruvb`](@ref): the Gibbs sampler, with the optional offset factor.
- [`summarize`](@ref), [`lfsr`](@ref): posterior summaries.
- [`log_ratio`](@ref): log-ratios to a set of reference genes, as Gaussian input.
