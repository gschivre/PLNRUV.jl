#Poisson-lognormal factor model fitted by variational empirical Bayes (design/plnpca.md)
#    yᵢⱼ | Zᵢⱼ ~ Poisson(exp(oᵢ + Zᵢⱼ)),   Zᵢⱼ | ... ~ N(xᵢ'Bⱼ + ℓᵢ + uᵢ'Γⱼ, σⱼ²)
#    uᵢ ~ N(0, I_q),   Γⱼₗ ~ N(0, τₗ²),   Bⱼ ~ N(b̄, τB² I),   σⱼ² ~ InvGamma(a0, b0)
#variational distributions, mean field: q(Zᵢⱼ) = N(Mᵢⱼ, S2ᵢⱼ), q(uᵢ) = N(mᵢ, S), q(Γⱼ) = N(gⱼ, Gⱼ), q(Bⱼ) = N(βⱼ, Vⱼ),
#q(σⱼ²) = InvGamma(aⱼ, bⱼ); the offset ℓ is a free parameter in the residual space of X (its component along X
#belongs to B). Given q(Z), the rest is a linear Gaussian factor model: every update is conjugate and exact, so the
#ELBO increases at every step. The expected squared residual seen by σⱼ² includes the uncertainty of Bⱼ, Γⱼ and uᵢ,
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
    m::Matrix{Float64} #n x q
    S::Matrix{Float64} #q x q, common to the samples
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

#E[(Zᵢⱼ - xᵢ'Bⱼ - ℓᵢ - uᵢ'Γⱼ)²] under q (n x p): squared residual of the means, plus S2ᵢⱼ, xᵢ'Vⱼxᵢ, gⱼ'Sgⱼ, mᵢ'Gⱼmᵢ
#and tr(GⱼS)
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
                E[i, j] += c + dot(view(st.m, i, :), G * view(st.m, i, :))
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
        elbo -= sum(kl_gauss(view(st.m, i, :), st.S, zeros(q), ones(q)) for i in 1:n)
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

#scores: q(uᵢ) Gaussian with the common covariance S = (I + Σⱼ wⱼ (gⱼgⱼ' + Gⱼ))⁻¹
function plnpca_scores!(st, X)
    q = size(st.g, 1)
    q == 0 && return st
    w = invσ2(st)
    P = Matrix{Float64}(I, q, q)
    for j in axes(st.g, 2)
        P .+= w[j] .* (st.g[:, j] * st.g[:, j]' .+ st.G[j])
    end
    st.S = Matrix(inv(Symmetric(P)))
    R = st.M .- X * st.β .- st.ℓ
    st.m .= (R .* w') * st.g' * st.S
    return st
end

#loadings: q(Γⱼ) Gaussian, precision wⱼ Σᵢ (mᵢmᵢ' + S) + diag(1 / τₗ²)
function plnpca_loadings!(st, X, pr)
    q = size(st.g, 1)
    q == 0 && return st
    w = invσ2(st)
    n = size(st.M, 1)
    A = st.m' * st.m .+ n .* st.S
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
function plnpca_sweep!(st, pr, Y, X, o, Q2)
    plnpca_cells!(st, Y, X, o)
    plnpca_coefficients!(st, X, pr)
    plnpca_offset!(st, X, Q2)
    plnpca_scores!(st, X)
    plnpca_loadings!(st, X, pr)
    plnpca_variances!(st, X, pr)
    plnpca_hyper!(pr, st)
    return st
end

#the primary parameters, extrapolated by SQUAREM (positive ones on the log scale); the covariances Vⱼ, Gⱼ, S, the
#shapes aⱼ and the hyperparameters are recomputed from them by the next sweep
plnpca_pack(st) = [vec(st.M); vec(log.(st.S2)); vec(st.β); vec(st.g); vec(st.m); st.ℓ; log.(st.b)]
function plnpca_unpack!(st, θ)
    i = 0
    for A in (st.M, st.S2, st.β, st.g, st.m)
        len = length(A)
        A .= reshape(θ[(i + 1):(i + len)], size(A))
        i += len
    end
    st.S2 .= exp.(st.S2)
    st.ℓ .= θ[(i + 1):(i + length(st.ℓ))]
    i += length(st.ℓ)
    st.b .= exp.(θ[(i + 1):end])
    return st
end

#one SQUAREM cycle (Varadhan R. & Roland C., Scand. J. Stat. (2008), scheme S3): two sweeps θ₁ = F(θ₀), θ₂ = F(θ₁),
#step α = -‖r‖ / ‖v‖ within [-stepmax, -1] (r = θ₁ - θ₀, v = θ₂ - 2θ₁ + θ₀), extrapolation θ₀ - 2αr + α²v and one
#stabilising sweep from it, kept when its ELBO is at least that of θ₂ (every sweep is an ascent of the exact ELBO, so
#this safeguard keeps the objective increasing); stepmax starts at 1 and is multiplied by 4 when the bound is used
function plnpca_squarem!(st, pr, Y, X, o, Q2, stepmax)
    θ0 = plnpca_pack(st)
    st1, pr1 = deepcopy(st), deepcopy(pr)
    plnpca_sweep!(st1, pr1, Y, X, o, Q2)
    θ1 = plnpca_pack(st1)
    st2, pr2 = deepcopy(st1), deepcopy(pr1)
    plnpca_sweep!(st2, pr2, Y, X, o, Q2)
    θ2 = plnpca_pack(st2)
    e2 = plnpca_elbo(st2, Y, X, o, pr2)
    r, v = θ1 .- θ0, θ2 .- 2 .* θ1 .+ θ0
    nv = norm(v)
    α = nv == 0 ? -1.0 : clamp(-norm(r) / nv, -stepmax[], -1.0)
    atbound = α == -stepmax[]
    nsweeps = 2
    if α < -1
        stx, prx = deepcopy(st2), deepcopy(pr2)
        plnpca_unpack!(stx, θ0 .- 2α .* r .+ α^2 .* v)
        ex = try
            plnpca_sweep!(stx, prx, Y, X, o, Q2)
            nsweeps += 1
            plnpca_elbo(stx, Y, X, o, prx)
        catch err
            err isa InterruptException && rethrow()
            -Inf
        end
        if isfinite(ex) && ex >= e2
            atbound && (stepmax[] *= 4)
            return stx, prx, ex, nsweeps
        end
    else
        atbound && (stepmax[] *= 4)
    end
    return st2, pr2, e2, nsweeps
end

"""
    plnpca(Y, X, q; offsets = log library sizes, offset_factor = true, τB = 5, maxiter = 1000, tol = 1e-8, accelerate = true)

experimental: Poisson-lognormal factor model fitted by variational empirical Bayes (design/plnpca.md), the same model
as pln with random factor scores and loadings, random coefficients (whose uncertainty gives the overdispersions a
REML-type correction) and an inverse-gamma prior on the overdispersions; every update is exact and the ELBO increases
at every sweep; the loadings have one prior variance per factor (automatic relevance determination) and the
hyperparameters are estimated from all genes; accelerate = true extrapolates the sweeps by SQUAREM; maxiter counts
sweeps, and the fit stops when the relative change of the ELBO falls below tol
return the state, the priors, the ELBO path, the number of sweeps and whether it converged
"""
function plnpca(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, q::Integer;
        offsets::AbstractVector{<:Real} = vec(log.(sum(Y; dims = 2))), offset_factor::Bool = true, τB::Real = 5.0,
        maxiter::Integer = 1000, tol::Real = 1.0e-8, accelerate::Bool = true
    )
    n, p = size(Y)
    k = size(X, 2)
    @assert n - k > q + Int(offset_factor) "not enough residual degrees of freedom"
    X = collect(float.(X))
    o = collect(float.(offsets))
    F = qr(X)
    Q2 = offset_factor ? (F.Q * Matrix(I, n, n))[:, (k + 1):n] : zeros(n, 0)
    #initialisation from the log counts: least squares, then an SVD of the residuals
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
        M, 1 ./ (Y .+ 1.0), β, [zeros(k, k) for _ in 1:p], g, [zeros(q, q) for _ in 1:p], m, Matrix{Float64}(I, q, q) ./ n,
        ℓ, fill(n / 2 + 1, p), s2 .* (n / 2)
    )
    pr = PLNPCAPriors(vec(median(β; dims = 2)), τB, q > 0 ? vec(mean(abs2, g; dims = 2)) : Float64[], 1.0, median(s2))
    plnpca_hyper!(pr, st)
    elbo = [plnpca_elbo(st, Y, X, o, pr)]
    sweeps = 0
    converged = false
    stepmax = Ref(1.0)
    while sweeps < maxiter
        if accelerate
            st, pr, e, ns = plnpca_squarem!(st, pr, Y, X, o, Q2, stepmax)
        else
            plnpca_sweep!(st, pr, Y, X, o, Q2)
            e, ns = plnpca_elbo(st, Y, X, o, pr), 1
        end
        sweeps += ns
        push!(elbo, e)
        if abs(elbo[end] - elbo[end - 1]) < tol * abs(elbo[end])
            converged = true
            break
        end
    end
    return (state = st, priors = pr, elbo = elbo, offsets = o, sweeps = sweeps, converged = converged)
end
