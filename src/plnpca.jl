#Poisson-lognormal factor model fitted by variational empirical Bayes (design/plnpca.md)
#    yᵢⱼ | Zᵢⱼ ~ Poisson(exp(oᵢ + Zᵢⱼ)),   Zᵢⱼ | ... ~ N(xᵢ'Bⱼ + ℓᵢ + uᵢ'Γⱼ, σⱼ²)
#    U = Q₂V,   Vᵣ ~ N(0, I_q),   Γⱼₗ ~ N(0, τₗ²),   Bⱼ ~ N(b̄, τB² I),   σⱼ² ~ InvGamma(a0, b0)
#the scores U (rows uᵢ) live in the residual space of X, spanned by the n - k orthonormal columns of Q₂, as the
#residual rows of the CATE rotation: the loadings are learned from the variation that the design cannot explain, and
#the factors cannot absorb the effects of interest
#variational distributions, mean field: q(Zᵢⱼ) = N(Mᵢⱼ, S2ᵢⱼ), q(Vᵣ) = N(m̃ᵣ, S) (so E[uᵢ] = mᵢ with m = Q₂m̃ and
#Cov(uᵢ) = hᵢS, hᵢ = (Q₂Q₂')ᵢᵢ), q(Γⱼ) = N(gⱼ, Gⱼ), q(Bⱼ) = N(βⱼ, Vⱼ), q(σⱼ²) = InvGamma(aⱼ, bⱼ); the offset ℓ is a
#free parameter in the residual space of X (its component along X belongs to B). Given q(Z), the rest is a linear
#Gaussian factor model: every update is conjugate and exact, so the ELBO increases at every step. The expected squared residual seen by σⱼ² includes the uncertainty of Bⱼ, Γⱼ and uᵢ,
#which gives the REML-type divisor (n - k - q, in the linear Gaussian model) instead of n
#hyperparameters (τₗ², a0, b0) are estimated from all genes at every sweep (empirical Bayes; one prior variance per
#factor: automatic relevance determination); b̄ and τB are fixed (b̄ the median of the initial coefficients)

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
    b0::Float64
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
    if q > 0
        #KL of the n - k rows of V: Σᵣ ‖m̃ᵣ‖² = ‖m‖² because Q₂ has orthonormal columns
        elbo -= 0.5 * (sum(st.h) * (tr(st.S) - q - logdet(st.S)) + sum(abs2, st.m))
        elbo -= sum(kl_gauss(view(st.g, :, j), st.G[j], zeros(q), pr.τ2) for j in 1:p)
    end
    elbo -= sum(kl_gauss(view(st.β, :, j), st.V[j], pr.b̄, fill(pr.τB^2, k)) for j in 1:p)
    for j in 1:p #E[log p(σⱼ²)] + entropy of q(σⱼ²)
        elbo += pr.a0 * log(pr.b0) - loggamma(pr.a0) - (pr.a0 + 1) * lσ[j] - pr.b0 * w[j] +
            st.a[j] + log(st.b[j]) + loggamma(st.a[j]) - (1 + st.a[j]) * digamma(st.a[j])
    end
    return elbo
end

#cells: the Lambert-W E-step of the PLN fit, prior N(xᵢ'βⱼ + ℓᵢ + mᵢ'gⱼ, 1 / E[1 / σⱼ²])
function plnpca_cells!(st, Y, X, o)
    μ = latent_mean(st, X)
    σ2 = 1 ./ invσ2(st)
    Threads.@threads for j in axes(Y, 2)
        for i in axes(Y, 1)
            st.M[i, j], st.S2[i, j] = estep_cell(Y[i, j], o[i], μ[i, j], σ2[j], st.M[i, j], st.S2[i, j])
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

#variances: q(σⱼ²) = InvGamma(a0 + n / 2, b0 + Σᵢ E[(Zᵢⱼ - ...)²] / 2)
function plnpca_variances!(st, X, pr)
    n = size(st.M, 1)
    E = expected_sq(st, X)
    st.a .= pr.a0 + n / 2
    st.b .= pr.b0 .+ vec(sum(E; dims = 1)) ./ 2
    return st
end

#hyperparameters from all genes: τₗ² = meanⱼ (gⱼₗ² + Gⱼₗₗ); (a0, b0) maximising Σⱼ E[log p(σⱼ²)]:
#b0 = p a0 / Σⱼ E[1 / σⱼ²] and log a0 - ψ(a0) = mean E[log σⱼ²] + log mean E[1 / σⱼ²] (Newton)
function plnpca_hyper!(pr, st)
    for l in eachindex(pr.τ2)
        pr.τ2[l] = max(mean(st.g[l, j]^2 + st.G[j][l, l] for j in axes(st.g, 2)), 1.0e-10)
    end
    w, lσ = invσ2(st), logσ2(st)
    c = mean(lσ) + log(mean(w))
    a0 = c <= 1.0e-8 ? 1.0e4 : clamp(0.5 / c, 1.0e-3, 1.0e4)
    for _ in 1:50
        f = log(a0) - digamma(a0) - c
        a0new = clamp(a0 - f / (1 / a0 - trigamma(a0)), a0 / 10, 1.0e4)
        done = abs(a0new - a0) < 1.0e-10 * a0
        a0 = a0new
        done && break
    end
    pr.a0 = a0
    pr.b0 = length(w) * a0 / sum(w)
    return pr
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
    return st
end

"""
    plnpca(Y, X, q; offsets = log library sizes, offset_factor = true, τB = 5, maxiter = 1000, tol = 1e-8)

experimental: Poisson-lognormal factor model fitted by variational empirical Bayes (design/plnpca.md), the model of
pln with random factor scores in the residual space of X, random loadings, random coefficients (whose uncertainty
gives the overdispersions a REML-type correction) and an inverse-gamma prior on the overdispersions; every update is
exact and the ELBO increases at every sweep; the loadings have one prior variance per factor (automatic relevance
determination) and the hyperparameters are estimated from all genes; maxiter counts sweeps, and the fit stops when
the relative change of the ELBO falls below tol
return a PLNFit (σ2 = 1 / E[1 / σⱼ²], df = 2a0 + n - k - q - 1 as for a variance shrunk by limma), usable in place
of pln by pln_mouthwash and choose_q (fit = :plnpca)
"""
function plnpca(Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer; offset_factor::Bool = true, kwargs...)
    r = plnpca_run(Y, X, q; offset_factor = offset_factor, kwargs...)
    st, pr = r.state, r.priors
    n, p = size(Y)
    k = size(X, 2)
    σ2 = st.b ./ st.a
    df = fill(2 * pr.a0 + n - k - q - Int(offset_factor), p)
    Γ = offset_factor ? [ones(1, p); st.g] : st.g
    σ2mix = (w = ones(p, 1), s2 = reshape(σ2, :, 1), df = reshape(df, :, 1))
    return PLNFit(st.M, st.S2, latent_mean(st, X), st.β, Γ, σ2, r.offsets, df, σ2mix, r.elbo, r.sweeps, r.converged)
end

#the fit itself: the variational state and the hyperparameters (prior variance of each factor in priors.τ2)
function plnpca_run(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        offsets::AbstractVector{<:Real} = vec(log.(sum(Y; dims = 2))), offset_factor::Bool = true, τB::Real = 5.0,
        maxiter::Integer = 1000, tol::Real = 1.0e-8
    )
    n, p = size(Y)
    k = size(X, 2)
    nfix = Int(offset_factor)
    @assert n - k > q + nfix "not enough residual degrees of freedom"
    X = collect(float.(X))
    o = collect(float.(offsets))
    F = qr(X)
    Q2 = (F.Q * Matrix(I, n, n))[:, (k + 1):n] #orthonormal basis of the residual space of X
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
    pr = PLNPCAPriors(vec(median(β; dims = 2)), τB, q > 0 ? vec(mean(abs2, g; dims = 2)) : Float64[], 1.0, median(s2))
    plnpca_hyper!(pr, st)
    elbo = plnpca_elbo(st, Y, X, o, pr)
    sweeps = 0
    converged = false
    while sweeps < maxiter
        plnpca_sweep!(st, pr, Y, X, o, Q2, offset_factor)
        sweeps += 1
        elbonew = plnpca_elbo(st, Y, X, o, pr)
        converged = abs(elbonew - elbo) < tol * abs(elbonew)
        elbo = elbonew
        converged && break
    end
    return (state = st, priors = pr, offsets = o, elbo = elbo, sweeps = sweeps, converged = converged)
end

#the PLN layer used by pln_mouthwash and choose_q
function fit_latent(method::Symbol, Y, X, q; kwargs...)
    method === :pln && return pln(Y, X, q; kwargs...)
    method === :plnpca && return plnpca(Y, X, q; kwargs...)
    throw(ArgumentError("unknown fit $method (use :pln or :plnpca)"))
end
