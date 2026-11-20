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
| ``X``, ``x_i`` | design matrix (``n \times k``, intercept first) and its row ``i`` |
| ``c`` | the contrast tested, ``c^\top B_j`` (default: the last coefficient) |
| ``a`` | weights of the raw estimate, ``a = X (X^\top X)^{-1} c`` |
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
| ``\alpha``, ``\alpha_0`` | association of the factors with the contrast tested (``q + 1`` values); ``\alpha_0`` is the compositional shift |

**Estimation**

| symbol | meaning |
|---|---|
| ``Q``, ``R`` | QR decomposition of the design, ``X = QR`` |
| ``u``, ``R_{kk}`` | the last column of ``Q`` and the last diagonal entry of ``R``; ``a = u / R_{kk}`` for the last coefficient |
| ``b_j``, ``s_j`` | raw (confounded) estimate of gene ``j`` and its standard error |
| ``d`` | residual degrees of freedom, ``n - k - q - 1`` |
| ``\hat\sigma_j^2``, ``\tilde\sigma_j^2`` | estimated and shrunk variances |
| ``d_0``, ``\sigma_0^2`` | degrees of freedom and scale of the limma prior on the variances; in the PLN fit ``d_0 = 2a_0`` and the scale is a function ``\sigma_0^2(x_j)`` of expression |
| ``x_j`` | mean log-abundance of gene ``j``, the covariate of the variance trend |
| ``t_j`` | degrees of freedom used by the coefficients and loadings of gene ``j`` (``\approx k + q``) |
| ``\text{Cov}(\Gamma_j)`` | posterior covariance of the loadings of gene ``j`` |
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

```@bibliography
*
```
