# Notation and references

## Notation

Each symbol has one meaning in the whole manual. Logarithms are natural; log₂ is used only to report simulated effects.

**Indices and data**

| symbol | meaning |
|---|---|
| ``i = 1, \dots, n`` | samples |
| ``j = 1, \dots, p`` | genes |
| ``y_{ij}`` | read count of gene ``j`` in sample ``i`` |
| ``N_i``, ``o_i`` | library size of sample ``i`` and the offset ``o_i = \log N_i`` |
| ``X``, ``x_i`` | design matrix (``n \times k``, intercept first, covariate of interest last) and its row ``i`` |
| ``C`` | the set of negative control genes |

**Model**

| symbol | meaning |
|---|---|
| ``Z_{ij}`` | latent log-abundance of gene ``j`` in sample ``i`` |
| ``\mu_{ij}`` | mean of ``Z_{ij}`` given the design and the factors |
| ``\sigma_j^2`` | latent variance (overdispersion) of gene ``j`` |
| ``\bar y_{ij}`` | expected count, ``e^{o_i + \mu_{ij} + \sigma_j^2/2}`` |
| ``B_j`` | coefficients of gene ``j`` on the design (``k`` values) |
| ``\beta_j`` | effect of interest of gene ``j``: the true change, the target of inference |
| ``q`` | number of free hidden factors |
| ``w_i`` | factors of sample ``i``: the offset factor ``w_{i0}`` (a shift common to all genes), then the ``q`` hidden factors |
| ``\Gamma_j`` | loadings of gene ``j`` on the factors (``q + 1`` values); ``\Gamma_{0j} = 1`` for the offset factor |
| ``\alpha``, ``\alpha_0`` | association of the factors with the covariate of interest (``q + 1`` values); ``\alpha_0`` is the compositional shift |

**Estimation**

| symbol | meaning |
|---|---|
| ``Q``, ``R`` | QR decomposition of the design, ``X = QR`` |
| ``u``, ``R_{kk}`` | the column of ``Q`` and the diagonal entry of ``R`` for the covariate of interest |
| ``b_j``, ``s_j`` | raw (confounded) estimate of gene ``j`` and its standard error |
| ``d`` | residual degrees of freedom, ``n - k - q - 1`` |
| ``\hat\sigma_j^2``, ``\tilde\sigma_j^2`` | estimated and shrunk variances |
| ``d_0``, ``\sigma_0^2`` | degrees of freedom and scale of the limma prior on the variances |
| ``c``, ``a_m``, ``\omega_m`` | scale, shapes and weights of the vash prior on the variances |
| ``g`` | prior distribution of the effects ``\beta_j`` (ash) |
| ``\pi_l``, ``\tau_l`` | weights and standard deviations of its components; ``\tau_0 = 0`` is the point mass at zero |
| ``\bar\beta`` | centre of the non-null part of ``g`` (zero unless a global change is detected) |
| ``\lambda`` | variance inflation of MOUTHWASH |
| ``\nu`` | degrees of freedom of the ``t`` likelihood |
| ``M_{ij}``, ``S^2_{ij}`` | variational mean and variance of ``Z_{ij}`` |

**Simulation and model choice**

| symbol | meaning |
|---|---|
| ``\rho_{ij}`` | probability of keeping a read in binomial thinning |
| ``U`` | factors estimated from the real counts, used to build confounded groups |
| ``\theta``, ``r_j`` | fraction of the counts in the training part, and precision of the beta-binomial split |

**Error rates**

| term | meaning |
|---|---|
| lfsr | local false sign rate: the probability that the reported sign of ``\beta_j`` is wrong |
| s-value | the average lfsr of all genes at least as significant |
| FDP | false discovery proportion: wrong discoveries (null genes or wrong signs) among all discoveries |

## References

- Aitchison J. & Ho C. H. (1989) The multivariate Poisson-log normal distribution. *Biometrika* 76, 643–653.
- Baker S. G. (1994) The multinomial–Poisson transformation. *J. R. Stat. Soc. D* 43, 495–504.
- Bourdon M. et al. (2012) Evidence for karyoplasmic homeostasis during endoreduplication and a ploidy-dependent
  increase in gene transcription during tomato fruit growth. *Development* 139, 3817–3826.
- Chiquet J., Mariadassou M. & Robin S. (2018) Variational inference for probabilistic Poisson PCA. *Ann. Appl. Stat.*
  12, 2674–2698.
- Coate J. E. & Doyle J. J. (2015) Variation in transcriptome size: are we getting the message? *Chromosoma* 124,
  27–43.
- Evans C., Hardin J. & Stoebel D. M. (2018) Selecting between-sample RNA-Seq normalization methods from the
  perspective of their assumptions. *Brief. Bioinform.* 19, 776–792.
- Ferreira P. G. et al. (2018) The effects of death and post-mortem cold ischemia on human tissue transcriptomes.
  *Nat. Commun.* 9, 490.
- Gagnon-Bartsch J. A., Jacob L. & Speed T. P. (2013) Removing unwanted variation from high dimensional data with
  negative controls. Technical report, UC Berkeley.
- Gallego Romero I., Pai A. A., Tung J. & Gilad Y. (2014) RNA-seq: impact of RNA degradation on transcript
  quantification. *BMC Biol.* 12, 42.
- Gerard D. (2020) Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics* 21, 206.
- Gerard D. & Stephens M. (2020) Empirical Bayes shrinkage and false discovery rate estimation, allowing for unwanted
  variation. *Biostatistics* 21, 15–32.
- Gerard D. & Stephens M. (2021) Unifying and generalizing methods for removing unwanted variation based on negative
  controls. *Stat. Sin.* 31, 1145–1166.
- Kim Y., Carbonetto P., Stephens M. & Anitescu M. (2020) A fast algorithm for maximum likelihood estimation of
  mixture proportions using sequential quadratic programming. *J. Comput. Graph. Stat.* 29, 261–273.
- Kim Y.-K. et al. (2023) Absolute scaling of single-cell transcriptomes identifies pervasive hypertranscription in
  adult stem and progenitor cells. *Cell Rep.* 42, 111978.
- Laosuntisuk K. et al. (2024) A normalization method that controls for total RNA abundance affects the identification
  of differentially expressed genes, revealing bias toward morning-expressed responses. *Plant J.* 118, 1241–1257.
- Law C. W., Chen Y., Shi W. & Smyth G. K. (2014) voom: precision weights unlock linear model analysis tools for
  RNA-seq read counts. *Genome Biol.* 15, R29.
- Lawson C. L. & Hanson R. J. (1974) *Solving Least Squares Problems*. Prentice-Hall.
- Leek J. T. et al. (2010) Tackling the widespread and critical impact of batch effects in high-throughput data.
  *Nat. Rev. Genet.* 11, 733–739.
- Lovén J. et al. (2012) Revisiting global gene expression analysis. *Cell* 151, 476–482.
- Lu M. & Stephens M. (2016) Variance adaptive shrinkage (vash): flexible empirical Bayesian estimation of variances.
  *Bioinformatics* 32, 3428–3434.
- Lu M. & Stephens M. (2019) Empirical Bayes estimation of normal means, accounting for uncertainty in estimated
  standard errors. arXiv:1901.10679.
- Marguerat S. et al. (2012) Quantitative analysis of fission yeast transcriptomes and proteomes in proliferating and
  quiescent cells. *Cell* 151, 671–683.
- McKnight J. N. et al. (2015) Global promoter targeting of a conserved lysine deacetylase for transcriptional
  shutoff during quiescence entry. *Mol. Cell* 59, 732–743.
- Neufeld A., Dharamshi A., Gao L. L. & Witten D. (2024) Data thinning for convolution-closed distributions.
  *J. Mach. Learn. Res.* 25(57), 1–35.
- Neufeld A., Gao L. L., Popp J., Battle A. & Witten D. (2024) Inference after latent variable estimation for
  single-cell RNA sequencing data. *Biostatistics* 25, 270–287.
- Nie Z. et al. (2012) c-Myc is a universal amplifier of expressed genes in lymphocytes and embryonic stem cells.
  *Cell* 151, 68–79.
- Owen A. B. & Wang J. (2016) Bi-cross-validation for factor analysis. *Stat. Sci.* 31, 119–139.
- Percharde M., Bulut-Karslioglu A. & Ramalho-Santos M. (2017) Hypertranscription in development, stem cells, and
  regeneration. *Dev. Cell* 40, 9–21.
- Pirrello J. et al. (2018) Transcriptome profiling of sorted endoreduplicated nuclei from tomato fruits: how the
  global shift in expression ascribed to DNA ploidy influences RNA-Seq data normalization and interpretation.
  *Plant J.* 93, 387–398.
- Quinn T. P. et al. (2018) Understanding sequencing data as compositions: an outlook and review. *Bioinformatics* 34,
  2870–2878.
- Robinson D. O. et al. (2018) Ploidy and size at multiple scales in the Arabidopsis sepal. *Plant Cell* 30,
  2308–2329.
- Schivre G. et al. (2025) Genome-scale transcriptome augmentation during *Arabidopsis thaliana* photomorphogenesis.
  *Nat. Commun.* 16, 11228.
- Shen-Orr S. S. & Gaujoux R. (2013) Computational deconvolution: extracting cell type-specific information from
  heterogeneous samples. *Curr. Opin. Immunol.* 25, 571–578.
- Smyth G. K. (2004) Linear models and empirical Bayes methods for assessing differential expression in microarray
  experiments. *Stat. Appl. Genet. Mol. Biol.* 3, 3.
- Stegle O., Parts L., Piipari M., Winn J. & Durbin R. (2012) Using probabilistic estimation of expression residuals
  (PEER) to obtain increased power and interpretability of gene expression analyses. *Nat. Protoc.* 7, 500–507.
- Stephens M. (2017) False discovery rates: a new deal. *Biostatistics* 18, 275–294.
- Swygert S. G. et al. (2019) Condensin-dependent chromatin compaction represses transcription globally during
  quiescence. *Mol. Cell* 73, 533–546.
- Taruttis F. et al. (2017) External calibration with *Drosophila* whole-cell spike-ins delivers absolute mRNA fold
  changes from human RNA-Seq and qPCR data. *Biotechniques* 62, 53–61.
- van de Peppel J. et al. (2003) Monitoring global messenger RNA changes in externally controlled microarray
  experiments. *EMBO Rep.* 4, 387–393.
- Wang J., Zhao Q., Hastie T. & Owen A. B. (2017) Confounder adjustment in multiple hypothesis testing. *Ann. Stat.*
  45, 1863–1894.
- Zatzman M. et al. (2022) Widespread hypertranscription in aggressive human cancers. *Sci. Adv.* 8, eabn0238.
- Zhou H., He K., Chen J. & Zhang X. (2022) LinDA: linear models for differential abundance analysis of microbiome
  compositional data. *Genome Biol.* 23, 95.
