# 1. The problem

## 1.1 When the whole transcriptome moves

Most analyses of RNA-seq data rest on a silent assumption: between the conditions compared, the total amount of
messenger RNA per cell stays the same, and only its distribution among genes changes. Normalisation by library size,
TMM or DESeq2's median of ratios all build this assumption in (Evans et al. 2018). It is often false.

- **Hypertranscription.** Stem and progenitor cells, regenerating tissues and early embryos raise their global
  transcriptional output (Percharde et al. 2017; Kim et al. 2023), and so do many aggressive cancers (Zatzman et al.
  2022). The oncogene c-Myc acts as a general amplifier of the genes already expressed (Nie et al. 2012). In that
  setting, Lovén et al. (2012) showed that standard normalisation hides the amplification, and can even reverse its
  apparent direction.
- **Global shutdown.** Quiescent yeast cells reduce their transcriptome several-fold (Marguerat et al. 2012), through
  genome-wide repressive mechanisms (McKnight et al. 2015; Swygert et al. 2019).
- **Cell size and ploidy.** Transcript content scales with cell size and with endoreduplication, which differ between
  tissues, developmental stages and genotypes (Bourdon et al. 2012; Pirrello et al. 2018; Robinson et al. 2018; Coate &
  Doyle 2015).
- **Environmental switches in plants.** The amount of RNA per cell follows the time of day (Laosuntisuk et al. 2024).
  When a dark-grown seedling is exposed to light, RNA Polymerase II activity increases, and the cotyledons roughly
  double their transcript levels within a day. Standard RNA-seq analyses of the same samples detected no such
  increase (Schivre et al. 2025).

The problem has been known since the microarray era (van de Peppel et al. 2003). Its root is that sequencing measures
proportions, not amounts: the information about the absolute content is lost when the sample is collected and
normalised, because sequencing counts are compositional (Quinn et al. 2018). Spike-ins added in a fixed proportion to
the cells (whole-cell spike-ins; Taruttis et al. 2017) can recover it, but they are rarely available.

A global change of transcription is also a confounder of a particular kind. Like a batch effect, it acts on many
genes at once. Unlike a batch effect, it acts on *all* of them, by the same factor. This page states the problem in
those terms. The rest of the manual shows that handling it as a confounder, with the tools developed for hidden
confounders, resolves it whenever there are negative controls (spike-ins, stable genes) or when most genes do not
change.

## 1.2 What a count measures

An RNA-seq library is a sample of the RNA molecules of a tissue. After sequencing, each read is assigned to a gene,
and the data are a table of counts ``y_{ij}``: the number of reads of gene ``j`` in sample ``i``, for ``n`` samples
and ``p`` genes (typically ``n`` between 6 and a few hundred, ``p`` around 20 000).

Write ``a_{ij}`` for the *absolute* abundance of gene ``j`` in sample ``i`` (molecules per cell, say). The expected
count is not ``a_{ij}``. The sequencer reads a budget of ``N_i`` molecules from library ``i``, and that budget is
shared among the genes in proportion to their abundance:

```math
\mathbb{E}[y_{ij}] \;=\; N_i \, \frac{a_{ij}}{\sum_{l} a_{il}} .
```

Taking logs, the expected log count is the log abundance plus a term that depends on the sample only:

```math
\log \mathbb{E}[y_{ij}] \;=\; \underbrace{\log N_i - \log \textstyle\sum_l a_{il}}_{o_i} \;+\; \log a_{ij} .
```

The sequencing depth ``N_i`` is known (it is the library size, up to filtering). The total RNA content
``\sum_l a_{il}`` is not. So the **offset** ``o_i`` is unknown, and it is different in every sample.

## 1.3 Are the data compositional?

Two views of sequencing data are common:

- the **multinomial** view: the library size ``N_i`` is fixed, and the reads are a multinomial draw over the genes;
- the **Poisson** view: each gene's count is an independent Poisson variable, and ``N_i`` is random.

**Result** — the two views give the same inference, as long as each sample has its own free offset. This is the
multinomial–Poisson transform (Baker 1994). The library size is in practice random: a sequencing run is shared
between projects, and reads are filtered after sequencing. But even if it were fixed, nothing would change.

!!! details "Derivation"
    Let ``y_{ij} \sim \text{Poisson}(\lambda_{ij})`` independently, with ``\lambda_{ij} = e^{o_i + \eta_{ij}}``. The
    total ``N_i = \sum_j y_{ij}`` is Poisson with mean ``\Lambda_i = e^{o_i} \sum_j e^{\eta_{ij}}``, and given
    ``N_i`` the counts are multinomial with probabilities ``e^{\eta_{ij}} / \sum_l e^{\eta_{il}}``, which do not
    depend on ``o_i``. The Poisson likelihood therefore factorises as
    ``\text{Poisson}(N_i; \Lambda_i) \times \text{Multinomial}(y_i; N_i, \text{softmax}(\eta_i))``. The first factor
    involves ``o_i`` through ``\Lambda_i``; with ``o_i`` free, maximising over it makes ``\Lambda_i = N_i`` and leaves
    the multinomial likelihood of ``\eta``, unchanged.

So the counts are not compositional *in their distribution*. They are compositional *in their mean*: the unknown
total RNA content sits in the offset ``o_i``. That is the part that matters.

## 1.4 The consequence: fold changes are only relative

Suppose a treatment multiplies the abundance of gene ``j`` by ``e^{\beta_j}``. Write the treatment as a covariate
``x_i`` (1 for treated samples, 0 for controls). Then

```math
\log \mathbb{E}[y_{ij}] = o_i + \log a^{(0)}_{j} + x_i \beta_j .
```

Now add any constant ``\delta`` to every ``\beta_j`` and subtract ``x_i\delta`` from every offset ``o_i``. The
expected counts do not change. **Result** — from the reads alone, the fold changes are identified only up to a
common shift: ``\beta_j`` and ``\beta_j + \delta`` fit the data equally well.

This is the "total mRNA" problem of §1.1. A global shift, like a general increase of transcription, is invisible. Worse, a strong change of a few highly expressed genes moves ``\sum_l a_{il}``, and therefore ``o_i``,
and every other gene appears to move in the opposite direction. Normalisation methods (TMM, DESeq2's median of
ratios, log-ratio transforms) all make a choice of ``\delta``, usually by assuming that most genes do not change.

Two kinds of information fix ``\delta``:

- **negative controls**: genes known not to change (``\beta_j = 0`` for ``j`` in a control set ``C``), for example
  spike-ins or stable reference genes;
- **sparsity**: the assumption that most genes do not change, so the distribution of the ``\beta_j`` has a mode at
  zero.

The two are not equivalent. Sparsity fails exactly in the situations of §1.1: when transcription rises globally, most
genes go up together. The mode of their changes is then the global change itself, and sparsity would set it to zero.
To measure a global change, controls are needed. Sparsity is enough when the question is about the changes of genes
relative to the bulk.

## 1.5 Hidden confounders

The second difficulty has the same shape. Unmeasured variables, such as batch, RNA quality or the cell-type
composition of a tissue, change the expression of many genes at once. Write them as a vector ``z_i`` of ``r``
hidden factors per sample, acting on gene ``j`` through loadings ``\gamma_j``:

```math
\log a_{ij} = \mu_j + x_i \beta_j + z_i^\top \gamma_j + \varepsilon_{ij} .
```

If the factors are independent of the treatment, they only add noise. If they are correlated with it (treated
samples processed on another day, say), part of ``z_i^\top\gamma_j`` looks like a treatment effect. Wang et al.
(2017) write this correlation as a regression of the factors on the design, ``z_i = \alpha^\top x_i + w_i``, so that
the apparent effect of the treatment on gene ``j`` is

```math
\tilde\beta_j = \beta_j + \alpha^\top \gamma_j .
```

**Result** — ``\beta_j`` is not identified either: for any ``\alpha``, the true effects can be traded against the
confounding. The same two kinds of information resolve it. Negative controls have ``\tilde\beta_j = \alpha^\top
\gamma_j``, which determines ``\alpha``. Without controls, sparsity does the same: the right ``\alpha`` is the one
that makes most ``\beta_j = \tilde\beta_j - \alpha^\top\gamma_j`` equal to zero.

## 1.6 One problem, not two

Compare §1.4 and §1.5. The compositional shift adds ``x_i\delta`` to every gene: it is a hidden factor whose effect
on gene ``j`` is ``\delta \times 1``. In other words, **the compositional offset is a confounder whose loading is the
same for every gene**. Its loading is known, and only its value, ``\delta``, has to be estimated, together with the
other confounders.

This is the idea of the package. A latent factor with loadings fixed to one (the *offset factor*) is added to the
hidden factors. Controls or sparsity then identify the compositional shift and the confounders at once, with the
same machinery. [§3](confounders.md) explains how, and [§8](validation.md) shows what happens when the offset factor
is left out.
