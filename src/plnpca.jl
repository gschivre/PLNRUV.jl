#Poisson-lognormal factor model fitted by variational empirical Bayes, the PLN layer of pln_mouthwash and choose_q
#    yᵢⱼ | Zᵢⱼ ~ Poisson(exp(oᵢ + Zᵢⱼ)),   Zᵢⱼ | ... ~ N(xᵢ'Bⱼ + ℓᵢ + uᵢ'Γⱼ, σⱼ²)
#    U = Q₂V,   Vᵣ ~ N(0, I_q),   Γⱼₗ ~ N(0, τₗ²),   Bⱼ ~ N(b̄, τB² I),   σⱼ² ~ InvGamma(a0, a0 h(xⱼ))
#the scores U (rows uᵢ) live in the residual space of X, spanned by the n - k orthonormal columns of Q₂, as the
#residual rows of the CATE rotation: the loadings are learned from the variation that the design cannot explain, and
#the factors cannot absorb the effects of interest
#variational distributions, mean field: q(Zᵢⱼ) = N(Mᵢⱼ, S2ᵢⱼ), q(Vᵣ) = N(m̃ᵣ, S) (so E[uᵢ] = mᵢ with m = Q₂m̃ and
#Cov(uᵢ) = hᵢS, hᵢ = (Q₂Q₂')ᵢᵢ), q(Γⱼ) = N(gⱼ, Gⱼ), q(Bⱼ) = N(βⱼ, Vⱼ), q(σⱼ²) = InvGamma(aⱼ, bⱼ); the offset ℓ is a
#free parameter in the residual space of X (its component along X belongs to B). Given q(Z), the rest is a linear
#Gaussian factor model: every update is conjugate and exact, so the ELBO increases at every step. The expected squared residual seen by σⱼ² includes the uncertainty of Bⱼ, Γⱼ and uᵢ,
#which gives the REML-type divisor (n - k - q, in the linear Gaussian model) instead of n
#hyperparameters (τₗ², a0, h) are estimated from all genes at every sweep (empirical Bayes; one prior variance per
#factor: automatic relevance determination); b̄ and τB are fixed (b̄ the median of the initial coefficients)
#h is the prior scale of the variances as a smooth function of the mean log-abundance xⱼ of each gene (fixed at its
#initial value): genes are shrunk towards genes of similar expression; on GTEx the variances follow a U-shaped trend
#in expression, which no monotone parametric curve (such as the dispersion trend of DESeq2) follows, so h is a local
#fit (var_prior = :common gives a single scale, as limma without trend)

mutable struct PLNPCAState
    M::Matrix{Float64} #n x p
    S2::Matrix{Float64} #n x p
    β::Matrix{Float64} #k x p
    V::Vector{Matrix{Float64}} #p, k x k
    g::Matrix{Float64} #q x p
    G::Vector{Matrix{Float64}} #p, q x q
    m::Matrix{Float64} #n x q, in the residual space of X
    S::Matrix{Float64} #q x q, covariance of each row of V
    h::Vector{Float64} #n, diagonal of Q₂Q₂' (Cov(uᵢ) = hᵢ S)
    ℓ::Vector{Float64} #n
    a::Vector{Float64} #p
    b::Vector{Float64} #p
end

mutable struct PLNPCAPriors
    b̄::Vector{Float64}
    τB::Float64
    τ2::Vector{Float64}
    a0::Float64
    b0::Vector{Float64} #a0 h(xⱼ), the rate of the prior of each σⱼ²
    x::Vector{Float64} #covariate of the variance trend (mean log-abundance), empty for a common prior
    span::Float64
end

invσ2(st) = st.a ./ st.b
logσ2(st) = log.(st.b) .- digamma.(st.a)

#prior mean of the latent values, xᵢ'βⱼ + ℓᵢ + mᵢ'gⱼ (n x p)
latent_mean(st, X) = X * st.β .+ st.ℓ .+ st.m * st.g

#E[(Zᵢⱼ - xᵢ'Bⱼ - ℓᵢ - uᵢ'Γⱼ)²] under q (n x p): squared residual of the means, plus S2ᵢⱼ, xᵢ'Vⱼxᵢ, mᵢ'Gⱼmᵢ and
#hᵢ (gⱼ'Sgⱼ + tr(GⱼS))
function expected_sq(st, X)
    E = (st.M .- latent_mean(st, X)) .^ 2 .+ st.S2
    n, p = size(st.M)
    for j in 1:p
        V = st.V[j]
        for i in 1:n
            E[i, j] += dot(view(X, i, :), V * view(X, i, :))
        end
        if size(st.g, 1) > 0
            gj = view(st.g, :, j)
            G = st.G[j]
            c = dot(gj, st.S * gj) + tr(G * st.S)
            for i in 1:n
                E[i, j] += st.h[i] * c + dot(view(st.m, i, :), G * view(st.m, i, :))
            end
        end
    end
    return E
end

#KL(N(μ, Σ) ‖ N(μ0, D)) with D diagonal (vector d)
kl_gauss(μ, Σ, μ0, d) = 0.5 * (sum(diag(Σ) ./ d) + sum((μ .- μ0) .^ 2 ./ d) - length(μ) + sum(log.(d)) - logdet(Σ))

function plnpca_elbo(st, Y, X, o, pr::PLNPCAPriors)
    n, p = size(Y)
    q = size(st.g, 1)
    k = size(X, 2)
    w, lσ = invσ2(st), logσ2(st)
    E = expected_sq(st, X)
    elbo = 0.0
    for j in 1:p, i in 1:n
        elbo += Y[i, j] * (o[i] + st.M[i, j]) - exp(o[i] + st.M[i, j] + st.S2[i, j] / 2) +
            0.5 * (log(st.S2[i, j]) + 1) - 0.5 * (lσ[j] + w[j] * E[i, j])
    end
    #normal-inverse-gamma family (see plnpca_variances!): uⱼ of the degrees of freedom of gene j carry no information
    #on σⱼ², -(n - uⱼ) / 2 E[log σⱼ²] - (wⱼ Σᵢ Eᵢⱼ - uⱼ) / 2 in place of the mean-field term
    elbo += 0.5 * sum(unresolved_dofs(st, X, w) .* (lσ .+ 1))
    if q > 0
        #KL of the n - k rows of V: Σᵣ ‖m̃ᵣ‖² = ‖m‖² because Q₂ has orthonormal columns
        elbo -= 0.5 * (sum(st.h) * (tr(st.S) - q - logdet(st.S)) + sum(abs2, st.m))
        elbo -= sum(kl_gauss(view(st.g, :, j), st.G[j], zeros(q), pr.τ2) for j in 1:p)
    end
    elbo -= sum(kl_gauss(view(st.β, :, j), st.V[j], pr.b̄, fill(pr.τB^2, k)) for j in 1:p)
    for j in 1:p #E[log p(σⱼ²)] + entropy of q(σⱼ²)
        elbo += pr.a0 * log(pr.b0[j]) - loggamma(pr.a0) - (pr.a0 + 1) * lσ[j] - pr.b0[j] * w[j] +
            st.a[j] + log(st.b[j]) + loggamma(st.a[j]) - (1 + st.a[j]) * digamma(st.a[j])
    end
    return elbo
end

#latent values of the counts: the Lambert-W E-step of the PLN fit, prior N(xᵢ'βⱼ + ℓᵢ + mᵢ'gⱼ, 1 / E[1 / σⱼ²])
function plnpca_cells!(st, Y, X, o)
    μ = latent_mean(st, X)
    σ2 = 1 ./ invσ2(st)
    Threads.@threads for j in axes(Y, 2)
        for i in axes(Y, 1)
            st.M[i, j], st.S2[i, j] = estep_count(Y[i, j], o[i], μ[i, j], σ2[j], st.M[i, j], st.S2[i, j])
        end
    end
    return st
end

#coefficients: q(Bⱼ) Gaussian, precision wⱼ X'X + I / τB², mean from the residual Mⱼ - ℓ - m gⱼ
function plnpca_coefficients!(st, X, pr)
    w = invσ2(st)
    XtX = X' * X
    k = size(X, 2)
    R = st.M .- st.ℓ .- st.m * st.g
    for j in axes(st.M, 2)
        V = inv(Symmetric(w[j] .* XtX + I(k) ./ pr.τB^2))
        st.V[j] = Matrix(V)
        st.β[:, j] = V * (w[j] .* (X' * view(R, :, j)) .+ pr.b̄ ./ pr.τB^2)
    end
    return st
end

#offset: weighted least squares of the residual on the samples, projected on the residual space of X (the weights
#Σⱼ wⱼ are the same for every sample, so the projection is exact)
function plnpca_offset!(st, X, Q2)
    size(Q2, 2) == 0 && return st
    w = invσ2(st)
    R = st.M .- X * st.β .- st.m * st.g
    ℓ = (R * w) ./ sum(w)
    st.ℓ .= Q2 * (Q2' * ℓ)
    return st
end

#scores: q(Vᵣ) Gaussian with the common covariance S = (I + Σⱼ wⱼ (gⱼgⱼ' + Gⱼ))⁻¹, and m = Q₂m̃ the projection on
#the residual space of X of the unconstrained update
function plnpca_scores!(st, X, Q2)
    q = size(st.g, 1)
    q == 0 && return st
    w = invσ2(st)
    P = Matrix{Float64}(I, q, q)
    for j in axes(st.g, 2)
        P .+= w[j] .* (st.g[:, j] * st.g[:, j]' .+ st.G[j])
    end
    st.S = Matrix(inv(Symmetric(P)))
    R = st.M .- X * st.β .- st.ℓ
    st.m .= Q2 * (Q2' * ((R .* w') * st.g' * st.S))
    return st
end

#loadings: q(Γⱼ) Gaussian, precision wⱼ Σᵢ (mᵢmᵢ' + S) + diag(1 / τₗ²)
function plnpca_loadings!(st, X, pr)
    q = size(st.g, 1)
    q == 0 && return st
    w = invσ2(st)
    A = st.m' * st.m .+ sum(st.h) .* st.S
    R = st.M .- X * st.β .- st.ℓ
    MR = st.m' * R #q x p
    for j in axes(st.g, 2)
        G = inv(Symmetric(w[j] .* A + Diagonal(1 ./ pr.τ2)))
        st.G[j] = Matrix(G)
        st.g[:, j] = G * (w[j] .* view(MR, :, j))
    end
    return st
end

#variances: q(σⱼ²) = InvGamma(a0 + (n - uⱼ) / 2, b0 + (Σᵢ E[(Zᵢⱼ - ...)²] - uⱼ / wⱼ) / 2), with uⱼ the degrees of
#freedom of gene j that carry no information on σⱼ²:
#- tⱼ = wⱼ (tr(X'X Vⱼ) + tr(A Gⱼ)), those used by its coefficients and loadings (the trace of its hat matrix);
#- wⱼ Σᵢ S2ᵢⱼ, those of its counts that leave their latent value uncertain: S2ᵢⱼ wⱼ = 1 / (1 + ȳᵢⱼσⱼ²) is 0 for a
#  count whose latent value is known and 1 for a count without information (S2ᵢⱼ = σⱼ²: it only returns σⱼ²)
#their part of the expected squared residuals is proportional to σⱼ², and moves from the rate to the shape: this is
#q(σⱼ²) of the normal-inverse-gamma family q(σⱼ²) q(Bⱼ, Γⱼ, Zⱼ | σⱼ²); the mean field gives the shape a0 + n / 2,
#with the same mean, but takes every count as a free observation of σⱼ², which with few samples or low counts makes
#q(σⱼ²) too concentrated and biases a0 upwards (towards complete pooling)
function plnpca_variances!(st, X, pr)
    n = size(st.M, 1)
    w = invσ2(st)
    E = vec(sum(expected_sq(st, X); dims = 1))
    u = unresolved_dofs(st, X, w)
    st.a .= pr.a0 .+ (n .- u) ./ 2
    st.b .= pr.b0 .+ (E .- u ./ w) ./ 2
    return st
end

#uⱼ = tⱼ + wⱼ Σᵢ S2ᵢⱼ (see plnpca_variances!), at most n (the ELBO term treats uⱼ as fixed within a sweep: with
#low counts the ELBO can decrease slightly in the first sweeps, while S2 and σ² are far from each other, and then
#increases monotonically): the two parts overlap when the counts carry little
#information (coefficients fitted to uninformative counts do not use information either), and a gene with uⱼ = n
#carries none on σⱼ², whose q is then the prior
unresolved_dofs(st, X, w) = min.(hat_traces(st, X, w) .+ w .* vec(sum(st.S2; dims = 1)), size(st.M, 1))

#tⱼ = wⱼ (tr(X'X Vⱼ) + tr(A Gⱼ)), with A = m'm + Σᵢ hᵢ S as in the loading update
function hat_traces(st, X, w)
    XtX = X' * X
    q = size(st.g, 1)
    A = q > 0 ? st.m' * st.m .+ sum(st.h) .* st.S : zeros(0, 0)
    return [w[j] * (dot(XtX, st.V[j]) + (q > 0 ? dot(A, st.G[j]) : 0.0)) for j in eachindex(w)]
end

const A0MAX = 1.0e3

#hyperparameters from all genes: τₗ² = meanⱼ (gⱼₗ² + Gⱼₗₗ); (a0, h) maximising Σⱼ E[log p(σⱼ²)], with
#E[log p(σⱼ²)] = a0 log(a0 hⱼ) - log Γ(a0) - (a0 + 1) E[log σⱼ²] - a0 hⱼ E[1 / σⱼ²]:
#h maximises Σⱼ (log hⱼ - hⱼ wⱼ) whatever a0 (local likelihood, see variance_scale_trend; a common h is 1 / mean w),
#then log a0 - ψ(a0) = mean E[log σⱼ²] - mean log hⱼ + mean hⱼwⱼ - 1 (Newton)
function plnpca_hyper!(pr, st)
    for l in eachindex(pr.τ2)
        pr.τ2[l] = max(mean(st.g[l, j]^2 + st.G[j][l, l] for j in axes(st.g, 2)), 1.0e-10)
    end
    w, lσ = invσ2(st), logσ2(st)
    h = isempty(pr.x) ? fill(1 / mean(w), length(w)) : variance_scale_trend(pr.x, w, pr.span)
    c = mean(lσ) - mean(log.(h)) + mean(h .* w) - 1
    #a0 is capped at 1000 (2000 prior degrees of freedom: complete pooling for any sample size); when all genes share
    #one variance around the trend, the estimate grows without bound
    a0 = c <= 1.0e-8 ? A0MAX : clamp(0.5 / c, 1.0e-3, A0MAX)
    for _ in 1:50
        f = log(a0) - digamma(a0) - c
        a0new = clamp(a0 - f / (1 / a0 - trigamma(a0)), a0 / 10, A0MAX)
        done = abs(a0new - a0) < 1.0e-10 * a0
        a0 = a0new
        done && break
    end
    pr.a0 = a0
    pr.b0 .= a0 .* h
    return pr
end

#prior scale h(x) of the variances: at each of ngrid points x₀ between the 1 % and 99 % quantiles of x, the
#local-linear likelihood maxθ Σₖ Kₖ (θ₀ + θ₁dₖ - exp(θ₀ + θ₁dₖ) wₖ), dₖ = xₖ - x₀, over the span · p genes nearest to
#x₀ with tricube weights K (concave, Newton); log h is interpolated linearly between the grid points and held flat
#outside them (no extrapolation); span 0.2 is the bandwidth that predicted best the variances of 200 GTEx samples
#from those of 4 to 10 others (the cross-validated error is flat between 0.05 and 0.2: the largest is the smoothest)
function variance_scale_trend(x, w, span; ngrid = 50)
    p = length(x)
    o = sortperm(x)
    xs, ws = x[o], w[o]
    lo, hi = quantile(xs, [0.01, 0.99])
    grid = hi > lo ? collect(range(lo, hi; length = ngrid)) : [lo]
    k = min(p, max(200, ceil(Int, span * p)))
    θ = zeros(length(grid))
    for (gi, x0) in enumerate(grid)
        #the k genes nearest to x0
        i = clamp(searchsortedfirst(xs, x0), 1, p)
        a, b = i, i
        while b - a + 1 < k
            if b == p || (a > 1 && x0 - xs[a - 1] <= xs[b + 1] - x0)
                a -= 1
            else
                b += 1
            end
        end
        d = xs[a:b] .- x0
        r = maximum(abs, d) * (1 + 1.0e-8) + eps()
        K = (1 .- (abs.(d) ./ r) .^ 3) .^ 3
        wl = ws[a:b]
        t0, t1 = -log(sum(K .* wl) / sum(K)), 0.0
        for _ in 1:50
            e = exp.(t0 .+ t1 .* d) .* wl .* K
            g0, g1 = sum(K) - sum(e), sum(K .* d) - sum(e .* d)
            H00, H01, H11 = sum(e), sum(e .* d), sum(e .* d .^ 2)
            D = H00 * H11 - H01^2
            D <= 1.0e-12 * H00 * H11 && (t0 += g0 / H00; break)
            s0, s1 = (H11 * g0 - H01 * g1) / D, (H00 * g1 - H01 * g0) / D
            t0 += s0
            t1 += s1
            abs(s0) + abs(s1) * (r + 1) < 1.0e-10 && break
        end
        θ[gi] = t0
    end
    length(grid) == 1 && return fill(exp(θ[1]), p)
    return [exp(interp_flat(grid, θ, xj)) for xj in x]
end

#linear interpolation, constant beyond the ends
function interp_flat(xg, yg, x)
    x <= xg[1] && return yg[1]
    x >= xg[end] && return yg[end]
    i = searchsortedlast(xg, x)
    u = (x - xg[i]) / (xg[i + 1] - xg[i])
    return (1 - u) * yg[i] + u * yg[i + 1]
end

#remove the factors whose prior variance τₗ² fell below threshold (switched off by ARD); true if any was removed
function plnpca_prune!(st, pr, threshold)
    keep = findall(>=(threshold), pr.τ2)
    length(keep) == length(pr.τ2) && return false
    st.g = st.g[keep, :]
    st.G = [G[keep, keep] for G in st.G]
    st.m = st.m[:, keep]
    st.S = st.S[keep, keep]
    pr.τ2 = pr.τ2[keep]
    return true
end

#one sweep of the coordinate ascent
function plnpca_sweep!(st, pr, Y, X, o, Q2, offset_factor)
    plnpca_cells!(st, Y, X, o)
    plnpca_coefficients!(st, X, pr)
    offset_factor && plnpca_offset!(st, X, Q2)
    plnpca_scores!(st, X, Q2)
    plnpca_loadings!(st, X, pr)
    plnpca_variances!(st, X, pr)
    plnpca_hyper!(pr, st)
    plnpca_expand!(st, pr)
    return st
end

#parameter-expanded step (PX-VB, Qi & Jaakkola 2007, after the PX-EM of Liu, Rubin & Wu 1998): the scores and the
#loadings enter the likelihood only through products that a transformation U -> U R, Γ -> R⁻¹Γ leaves unchanged, so
#R can maximise the prior terms alone, in closed form: it whitens the scores (Σᵣ E[vᵣvᵣ'] = N I, N = n - k rows) and
#makes Σⱼ E[ΓⱼΓⱼ'] diagonal (the ARD prior is diagonal), R = √N A^(-1/2) W with A = Σᵣ E[vᵣvᵣ'] and W the
#eigenvectors of A^(1/2) Σⱼ E[ΓⱼΓⱼ'] A^(1/2); then τ² at its optimum; R = I is one candidate, so the ELBO cannot
#decrease; alternating updates of scores and loadings crawl along the directions where one compensates the other
#(weak factors: hundreds of sweeps), the step jumps along them
function plnpca_expand!(st, pr)
    q = size(st.g, 1)
    q == 0 && return st
    N = sum(st.h)
    E = eigen(Symmetric(st.m' * st.m .+ N .* st.S))
    Ah = E.vectors * Diagonal(sqrt.(E.values)) * E.vectors'
    W = eigen(Symmetric(Ah * Symmetric(st.g * st.g' .+ sum(st.G)) * Ah)).vectors
    R = sqrt(N) .* (E.vectors * Diagonal(1 ./ sqrt.(E.values)) * E.vectors') * W
    Ri = inv(R)
    st.m = st.m * R
    st.S = Matrix(Symmetric(R' * st.S * R))
    st.g = Ri * st.g
    st.G = [Matrix(Symmetric(Ri * G * Ri')) for G in st.G]
    for l in 1:q
        pr.τ2[l] = max(mean(st.g[l, j]^2 + st.G[j][l, l] for j in axes(st.g, 2)), 1.0e-10)
    end
    return st
end

#one more factor for a warm start: scores and loadings from the leading singular vector of the residuals of the
#current fit (in the residual space of X), the other parts of the state kept
function plnpca_add_factor(st, pr, X)
    n, p = size(st.M)
    k = size(X, 2)
    Q2 = (qr(X).Q * Matrix(I, n, n))[:, (k + 1):n]
    R = Q2 * (Q2' * (st.M .- X * st.β .- st.ℓ .- st.m * st.g))
    D = svd(R)
    mnew = D.U[:, 1] .* sqrt(n)
    gnew = (D.S[1] .* D.Vt[1, :]) ./ sqrt(n)
    q = size(st.g, 1)
    st2, pr2 = deepcopy(st), deepcopy(pr)
    st2.m = [st.m mnew]
    st2.g = [st.g; gnew']
    st2.S = [st.S zeros(q, 1); zeros(1, q) 1 / n]
    st2.G = [[G zeros(q, 1); zeros(1, q) 0.0] for G in st.G]
    pr2.τ2 = [pr.τ2; mean(abs2, gnew)]
    return st2, pr2
end

"""
    plnpca(Y, X, q; offsets = log library sizes, offset_factor = true, τB = 5, var_prior = :trend, span = 0.2,
           maxiter = 1000, tol = 1e-7, prune = 1e-3)

Poisson-lognormal factor model fitted by variational empirical Bayes: counts Y (samples x genes), design X and q
latent factors, with random factor scores in the residual space of X, random loadings, random coefficients (whose
uncertainty gives the overdispersions a REML-type correction) and an inverse-gamma prior on the overdispersions,
whose scale follows the mean log-abundance of the genes (var_prior = :trend, a local fit over a fraction span of the
genes; :common for a single scale, which shrinks the variances of weakly and of highly expressed genes towards the
middle); every update is exact and the ELBO increases at every sweep; the loadings have one prior variance per
factor (automatic relevance determination) and the hyperparameters are estimated from all genes; a factor whose
prior variance falls below prune times the median variance is removed (q is an upper bound); maxiter counts sweeps,
and the fit stops when the ELBO changes by less than tol per count (gene in a sample)
return a PLNFit (σ2 = 1 / E[1 / σⱼ²], df = 2a0 + n - tⱼ as for a variance shrunk by limma, tⱼ ≈ k + q the degrees
of freedom used by the coefficients and loadings of gene j, Γcov the posterior covariance of the loadings); the fit
of pln_mouthwash and choose_q
on very sparse counts, or low counts with nearly equal variances (complete pooling), the variances approach their value
slowly from above and the fit may stop at maxiter (converged = false), with overestimated, conservative variances
"""
function plnpca(Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer; offset_factor::Bool = true, kwargs...)
    r = plnpca_run(Y, X, q; offset_factor = offset_factor, kwargs...)
    st, pr = r.state, r.priors
    n, p = size(Y)
    k = size(X, 2)
    σ2 = st.b ./ st.a
    df = 2 .* st.a
    #the factors switched off by ARD were removed during the fit: they carry no confounding, and their coefficients α
    #in MOUTHWASH would not be identified
    Γ = offset_factor ? [ones(1, p); st.g] : st.g
    σ2mix = (w = ones(p, 1), s2 = reshape(σ2, :, 1), df = reshape(df, :, 1))
    #posterior covariance of each column of Γ (the loadings of the offset factor are known)
    Γcov = offset_factor ? [cat(zeros(1, 1), G; dims = (1, 2)) for G in st.G] : copy(st.G)
    return PLNFit(st.M, st.S2, latent_mean(st, X), st.β, Γ, σ2, r.offsets, df, σ2mix, Γcov, r.elbo, r.sweeps, r.converged)
end

#the fit itself: the variational state and the hyperparameters (prior variance of each factor in priors.τ2)
function plnpca_run(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        offsets::AbstractVector{<:Real} = vec(log.(sum(Y; dims = 2))), offset_factor::Bool = true, τB::Real = 5.0,
        var_prior::Symbol = :trend, span::Real = 0.2, maxiter::Integer = 1000, tol::Real = 1.0e-7, prune::Real = 1.0e-3,
        init = nothing
    )
    var_prior in (:trend, :common) || throw(ArgumentError("unknown var_prior $var_prior (use :trend or :common)"))
    n, p = size(Y)
    k = size(X, 2)
    nfix = Int(offset_factor)
    @assert n - k > q + nfix "not enough residual degrees of freedom"
    X = collect(float.(X))
    o = collect(float.(offsets))
    F = qr(X)
    Q2 = (F.Q * Matrix(I, n, n))[:, (k + 1):n] #orthonormal basis of the residual space of X
    if init !== nothing #warm start from the state and priors of a previous fit (choose_q)
        st, pr = deepcopy(init[1]), deepcopy(init[2])
        return plnpca_iterate!(st, pr, Y, X, o, Q2, offset_factor, maxiter, tol, prune)
    end
    #initialisation from the log counts: least squares, then an SVD of the residuals (in the residual space of X)
    M = log.(Y .+ 0.5) .- o
    β = F \ M
    R = M .- X * β
    ℓ = offset_factor ? Q2 * (Q2' * vec(mean(R; dims = 2))) : zeros(n)
    R .-= ℓ
    D = svd(R)
    m = q > 0 ? D.U[:, 1:q] .* sqrt(n) : zeros(n, 0)
    g = q > 0 ? (D.S[1:q] .* D.Vt[1:q, :]) ./ sqrt(n) : zeros(0, p)
    s2 = max.(vec(var(R .- m * g; dims = 1)), 1.0e-3)
    st = PLNPCAState(
        M, 1 ./ (Y .+ 1.0), β, [zeros(k, k) for _ in 1:p], g, [zeros(q, q) for _ in 1:p], m,
        Matrix{Float64}(I, q, q) ./ n, vec(sum(abs2, Q2; dims = 2)), ℓ, fill(n / 2 + 1, p), s2 .* (n / 2)
    )
    x = var_prior === :trend ? vec(mean(M; dims = 1)) : Float64[] #mean log-abundance, the covariate of the trend
    pr = PLNPCAPriors(
        vec(median(β; dims = 2)), τB, q > 0 ? vec(mean(abs2, g; dims = 2)) : Float64[], 1.0, fill(median(s2), p), x,
        float(span)
    )
    plnpca_hyper!(pr, st)
    return plnpca_iterate!(st, pr, Y, X, o, Q2, offset_factor, maxiter, tol, prune)
end

#the sweeps of plnpca_run, from a given state
function plnpca_iterate!(st, pr, Y, X, o, Q2, offset_factor, maxiter, tol, prune)
    n, p = size(Y)
    elbo = plnpca_elbo(st, Y, X, o, pr)
    sweeps = 0
    converged = false
    while sweeps < maxiter
        plnpca_sweep!(st, pr, Y, X, o, Q2, offset_factor)
        sweeps += 1
        elbonew = plnpca_elbo(st, Y, X, o, pr)
        #a factor that ARD switches off has its prior variance decrease only geometrically, which keeps the ELBO
        #moving: it is removed once negligible, and the fit goes on with the others (the ELBO jumps at a removal)
        if plnpca_prune!(st, pr, prune * median(st.b ./ st.a))
            elbo = plnpca_elbo(st, Y, X, o, pr)
            continue
        end
        #change per count: the ELBO itself is dominated by constant terms (Σ yᵢⱼoᵢ), so a change relative to it would
        #stop the fit long before the parameters settle
        converged = abs(elbonew - elbo) < tol * n * p
        elbo = elbonew
        converged && break
    end
    return (state = st, priors = pr, offsets = o, elbo = elbo, sweeps = sweeps, converged = converged)
end
