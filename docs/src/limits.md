# 8. Assumptions and limits

## 8.1 The assumptions

| assumption | where | what breaks it | what to do |
|---|---|---|---|
| counts are Poisson given the latent abundance | [§2.1](model.md) | technical noise beyond the lognormal (PCR duplicates); non-integer counts | round the counts |
| one Gaussian variance per gene | [§2.1](model.md) | outlier samples; samples with different overall variability | inspect the samples; see §8.3 |
| hidden factors act linearly on the log scale | [§1.4](problem.md) | too few factors leave confounding; too many fit the noise | choose ``q`` with `choose_q` ([§3.3](confounders.md)) |
| hidden factors also vary within groups | [§3.1](confounders.md) | a factor that only differs between groups cannot be learned from the residual rows | it is then inseparable from the effect: see §8.2 |
| controls are truly null | [§2.3](model.md) | controls that shift with the condition, even slightly | compare with the sparsity branch; test with `ctrl_signal` ([§7.4](validation.md)) |
| or: most effects are zero | [§4.1](empirical_bayes.md) | a global change of transcription | use controls |
| a global change multiplies all genes equally | [§1.5](problem.md) | a change that favours some genes | read the gene effects as relative to the common part |
| one contrast at a time | [§5.6](pln.md) | a question about several coefficients at once ("any change over time") | test the contrasts one by one |
| every coefficient estimable | [§5.6](pln.md) | a gene with no count in a group, level or time point | it is reported as not tested |
| confounders are plugged in | [§4.3](empirical_bayes.md) | weak factors, very few samples; loadings plugged in by the Gaussian log-ratio branch ([§5.7](pln.md)) | see §8.4 and §8.5 |

## 8.2 What the reads cannot tell

Some questions have no answer from the reads alone, whatever the model:

- **A change common to all genes** is invisible without controls ([§1.3](problem.md)).
- **A hidden factor aligned with the condition** cannot be told apart from the effect. The controls, or sparsity,
  decide how the apparent effect is split, and they can be wrong.
- **A factor caused by the condition** (a change of the cell-type mix, for example) is part of the effect. Removing it
  changes the question ([§1.4](problem.md)).

## 8.3 Interpreting a global change

The estimate ``\hat\alpha_0`` and the effects ``\beta_j`` are always relative to a reference. What the reference is
depends on the assumption used.

**With controls.** ``\hat\alpha_0`` is the change of the bulk of genes relative to the controls. It is a change of the
RNA content per cell only if all of these hold:

1. **The controls do not change per cell.** Spike-ins must be added in proportion to the number of cells, not to the
   amount of RNA. Endogenous reference genes must be stable in this system, which needs independent evidence.
2. **The controls go through the same measurement as the other genes.** Different extraction or capture efficiency
   adds a shift between the controls and the genes.
3. **The controls are not shifted together.** Small independent errors average out; a shared shift does not
   ([§7.4](validation.md)). It moves ``\hat\alpha_0`` and every gene with it.
4. **The change is the same for all genes.** Otherwise ``\hat\alpha_0`` is only its common part ([§1.5](problem.md)).

**Without controls.** ``\hat\alpha_0`` puts the most common change at zero. The effects are relative to the bulk of
genes. This branch cannot detect a global change, and should not be used to claim its absence.

A global change should therefore be reported together with its identifying assumption, for example: "relative to
spike-ins added per cell, the transcripts increase by …". The spike-in pipeline of Schivre et al. (2025)
and an independent measurement, such as RNA per cell or imaging, are useful checks.

## 8.4 Few samples

The benchmark covers 2 to 10 samples per group ([§7.4](validation.md)). With valid controls, PLNRUV is as close to
the nominal level as an oracle that knows the true normalisation. What few samples cost:

- **Two samples per group.** No hidden factor can be fitted: the design and the offset leave one residual degree of
  freedom. Confounding is then not corrected, and a confounded design gives false discoveries. With random groups even
  the oracle is somewhat above the nominal level.
- **Three samples per group.** Each gene has one or two residual degrees of freedom. Only strong effects can be
  found.
- **Chance confounding.** With 2 or 3 samples per group, a random split is often almost aligned with a main source of
  variation of the samples. Look at the main factors of your samples before trusting the results.
- **The Gaussian branch** is not calibrated with 2 or 3 samples per group. Use the PLN branch there.
- **Ten samples per group or more.** The ``t`` likelihood is conservative: the achieved false discovery proportion is
  below the nominal level, and some power is lost. The normal likelihood (`likelihood = :normal`) is calibrated there,
  but not with fewer samples.

Real samples also differ in overall variability, which the model does not describe. With few samples, treat the
discoveries as candidates to confirm.

## 8.5 The Gaussian branch

The Gaussian [`mouthwash`](@ref) on log-ratios ([§5.7](pln.md)) does not model the Poisson noise, distorts low counts
through the pseudocount, treats the loadings as exact, and tests one column of the design (no contrasts yet). Prefer the PLN fit at low depth and in confounded designs
with few samples.
