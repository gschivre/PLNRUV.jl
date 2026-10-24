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
| one covariate of interest | [§3.1](confounders.md) | several contrasts | run once per contrast |
| loadings and confounders are plugged in | [§4.3](empirical_bayes.md) | weak factors, low counts, very few samples | see §8.3 |

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

With fewer than five samples per group, the false discovery rate is not controlled ([§7.5](validation.md)). The
variance of each gene then rests on one or two degrees of freedom. Real samples also differ in overall variability,
which the model does not describe. With such designs, treat the discoveries as candidates to confirm.
