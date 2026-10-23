# 8. Assumptions and limits

| assumption | where | what would break it, and how to check |
|---|---|---|
| counts are Poisson given the latent abundance | [§2.1](model.md) | technical overdispersion beyond lognormal (PCR duplicates); multi-mapping reads give non-integer counts, which must be rounded |
| latent noise is Gaussian with one variance per gene | [§2.1](model.md) | outlier samples; heavy-tailed residuals |
| the hidden factors act linearly on the log scale, with ``q`` factors | [§2.2](model.md) | choose ``q`` with `choose_q` ([§3.2](confounders.md)); too few factors leave confounding, too many cost power and, with few samples, fit the noise |
| negative controls are truly null | [§2.3](model.md) | a shift of the controls moves every estimate; test with `ctrl_signal` in the simulator |
| or: most effects are zero and the effect distribution is unimodal at zero | [§4.1](empirical_bayes.md) | global changes of transcription (most genes up); then controls (spike-ins) are required |
| one covariate of interest at a time | [§3.1](confounders.md) | for several contrasts, run once per covariate; a joint (mash) prior is planned |
| the variational posterior is Gaussian and independent across cells | [§5.1](pln.md) | its uncertainty is not used for the standard errors, which come from the model ([§5.4](pln.md)) |
| the estimated loadings ``\Gamma`` and confounders ``\alpha`` are plugged in (empirical Bayes branch) | [§4.2](empirical_bayes.md) | slight excess of false discoveries with weak factors and low counts, and with very few residual degrees of freedom (2–3 samples per group); let `choose_q` choose fewer factors |

Some questions cannot be answered from the reads alone, whatever the model:

- a change common to all genes (total transcription) is invisible without spike-ins or controls ([§1.4](problem.md));
- a hidden factor perfectly aligned with the design cannot be told apart from the effect: the controls, or sparsity,
  decide how it is split, and they can be wrong.
