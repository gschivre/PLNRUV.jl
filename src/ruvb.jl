#RUVB of Gerard D. & Stephens M., Stat. Sin. (2021): https://doi.org/10.5705/ss.202018.0345
#the data Y (samples x genes, Gaussian scale e.g. log-ratios) are rotated by the QR decomposition of X = [X₁ X₂]
#(X₂ the covariates of interest); the rows of X₁ are dropped, a Bayesian factor analysis is fitted on the remaining rows
#where the X₂ rows of the non-control genes are treated as missing, and β₂ = R₂₂⁻¹ (Y_obs - Y_imputed)
#
#factor model on the m = k₂ + n - k kept rows: Ỹ = L * F + E, E ~ MatrixNormal(0, I, Σ), Σ = diag(σ²)
#priors (half-t of Huang & Wand (2013) through gamma mixtures):
#    F[:, j] ~ N(0, σ²ⱼ I)              σ⁻²ⱼ ~ Gamma(ν / 2, rate ν ϕ)       ϕ ~ Gamma(1 / 2, rate 1 / a2)
#    L[:, l] ~ N(0, ψₗ I)               ψ⁻¹ₗ ~ Gamma(τ / 2, rate τ ηₗ)      ηₗ ~ Gamma(1 / 2, rate 1 / b2)
#with offset_factor = true the first row of F is fixed to ones: a factor that shifts all genes of a sample
#by the same amount (compositional / normalisation offset), its X₂ loadings are identified by the controls

struct RUVBFit
    β::Array{Float64, 4} #k₂ x p x nsample x nchain, controls are NaN (fixed to 0 by assumption)
    controls::Vector{Int}
    cov_of_interest::Vector{Int}
    acceptance::Vector{Float64} #acceptance rate of the latent updates per chain (empty on a Gaussian scale)
end

#Gibbs sampler state of one chain
mutable struct RUVBState
    L::Matrix{Float64} #m x qa
    F::Matrix{Float64} #qa x p, first row fixed to ones if there is an offset factor
    σ2inv::Vector{Float64}
    ψinv::Vector{Float64}
    η::Vector{Float64}
    ϕ::Float64
    Yt::Matrix{Float64} #rotated data with the missing block imputed
end

#initialise a chain: truncated SVD of the residual rows, loadings of the X₂ rows by least squares on the controls
function ruvb_init(rng::AbstractRNG, Ỹ, k₂, q, controls, noncontrols, nfix)
    m, p = size(Ỹ)
    Yres = Ỹ[(k₂ + 1):m, :]
    mres = m - k₂
    l0 = nfix > 0 ? vec(mean(Yres; dims = 2)) : zeros(mres)
    D = svd(Yres .- l0)
    F1 = Diagonal(D.S[1:q]) * D.Vt[1:q, :] ./ sqrt(mres)
    F1 .+= 0.1 * std(F1) .* randn(rng, q, p) #overdispersed starting points across chains
    F = nfix > 0 ? [ones(1, p); F1] : F1
    Lres = nfix > 0 ? [l0 D.U[:, 1:q] .* sqrt(mres)] : D.U[:, 1:q] .* sqrt(mres)
    Lint = Ỹ[1:k₂, controls] / F[:, controls]
    L = [Lint; Lres]
    Yt = copy(Ỹ)
    Yt[1:k₂, noncontrols] .= Lint * F[:, noncontrols]
    σ2inv = 1 ./ max.(vec(mean(abs2, view(Yt - L * F, (k₂ + 1):m, :); dims = 1)), 1.0e-8)
    ψinv = 1 ./ max.(vec(mean(abs2, L; dims = 1)), 1.0e-8)
    return RUVBState(L, F, σ2inv, ψinv, ones(size(L, 2)), 1.0, Yt)
end

#one Gibbs sweep
function ruvb_step!(rng::AbstractRNG, s::RUVBState, k₂, noncontrols, nfix, a2, b2, ν, τ)
    m, p = size(s.Yt)
    qa = size(s.L, 2)
    q = qa - nfix
    free = (nfix + 1):qa

    #hyperparameters of the half-t priors
    for l in 1:qa
        s.η[l] = rand(rng, Gamma((τ + 1) / 2, inv(1 / b2 + τ * s.ψinv[l])))
        s.ψinv[l] = rand(rng, Gamma((τ + m) / 2, inv(τ * s.η[l] + sum(abs2, view(s.L, :, l)) / 2)))
    end
    s.ϕ = rand(rng, Gamma((p * ν + 1) / 2, inv(1 / a2 + ν * sum(s.σ2inv))))

    #residual variances (the free rows of F share the scale σ²ⱼ)
    E = s.Yt - s.L * s.F
    for j in 1:p
        ssq = sum(abs2, view(E, :, j)) + sum(abs2, view(s.F, free, j))
        s.σ2inv[j] = rand(rng, Gamma((ν + m + q) / 2, inv(ν * s.ϕ + ssq / 2)))
    end

    #loadings: rows of L ~ N(P⁻¹ F Σ⁻¹ ỹᵢ, P⁻¹) with P = F Σ⁻¹ F' + diag(ψ⁻¹)
    FS = s.F .* s.σ2inv'
    C = cholesky(Symmetric(FS * s.F' + Diagonal(s.ψinv)))
    s.L .= (C \ (FS * s.Yt'))' .+ randn(rng, m, qa) / C.U'

    #factors: F[free, j] ~ N(G⁻¹ L₁' (ỹⱼ - L₀ F₀ⱼ), σ²ⱼ G⁻¹) with G = L₁' L₁ + I
    L1 = view(s.L, :, free)
    R = nfix > 0 ? s.Yt - view(s.L, :, 1:nfix) * view(s.F, 1:nfix, :) : s.Yt
    G = cholesky(Symmetric(L1' * L1 + I))
    s.F[free, :] .= (G \ (L1' * R)) .+ (G.U \ randn(rng, q, p)) ./ sqrt.(s.σ2inv')

    #impute the X₂ rows of the non-control genes from the posterior predictive
    for j in noncontrols
        σ = 1 / sqrt(s.σ2inv[j])
        for i in 1:k₂
            s.Yt[i, j] = dot(view(s.L, i, :), view(s.F, :, j)) + σ * randn(rng)
        end
    end
    return s
end

"""
    ruvb(Y, X, controls, q; cov_of_interest = [size(X, 2)], offset_factor = true, a2 = 1, b2 = 1, ν = 2, τ = 2, burnin = 2500, nsample = 500, thin = 10, nchain = 4, rng)

fit RUVB with q latent factors (plus the offset factor if offset_factor = true)
controls are the indices of the negative control genes, cov_of_interest the columns of X whose effects are estimated
"""
function ruvb(
        Y::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}, controls::AbstractVector{<:Integer}, q::Integer;
        cov_of_interest::AbstractVector{<:Integer} = [size(X, 2)], offset_factor::Bool = true,
        a2::Real = 1.0, b2::Real = 1.0, ν::Real = 2, τ::Real = 2,
        burnin::Integer = 2_500, nsample::Integer = 500, thin::Integer = 10, nchain::Integer = 4,
        rng::AbstractRNG = Random.default_rng()
    )
    n, p = size(Y)
    k = size(X, 2)
    @assert size(X, 1) == n
    nfix = Int(offset_factor)
    @assert length(controls) >= q + nfix "need at least q + offset_factor control genes"
    k₂ = length(cov_of_interest)
    @assert n - k >= q + 1 "not enough residual degrees of freedom for q factors"

    #rotate and drop the nuisance rows
    rot = rotate(Y, X, cov_of_interest)
    Ỹ = [rot.Yint; rot.Yres]
    R22inv = inv(rot.R22)
    noncontrols = setdiff(1:p, controls)
    Yobs = Ỹ[1:k₂, noncontrols]

    β = fill(NaN, k₂, p, nsample, nchain)
    seeds = rand(rng, UInt64, nchain)
    tasks = map(1:nchain) do c
        Threads.@spawn begin
            crng = Xoshiro(seeds[c])
            s = ruvb_init(crng, Ỹ, k₂, q, controls, noncontrols, nfix)
            for _ in 1:burnin
                ruvb_step!(crng, s, k₂, noncontrols, nfix, a2, b2, ν, τ)
            end
            for t in 1:nsample
                for _ in 1:thin
                    ruvb_step!(crng, s, k₂, noncontrols, nfix, a2, b2, ν, τ)
                end
                β[:, noncontrols, t, c] .= R22inv * (Yobs - view(s.Yt, 1:k₂, noncontrols))
            end
        end
    end
    foreach(wait, tasks)
    return RUVBFit(β, collect(controls), collect(cov_of_interest), Float64[])
end

"""
    lfsr(x)

local false sign rate of a vector of posterior samples
"""
lfsr(x::AbstractVector{<:Real}) = min(mean(>=(0), x), mean(<=(0), x))

"""
    summarize(fit::RUVBFit; level = 0.95, shrink = true)

posterior summaries of the coefficients of interest (k₂ x p matrices, NaN for the controls), with the convergence
of each coefficient (rank-normalised split-R̂ and bulk ESS)
the flat prior of RUVB gives posterior means, sds, credible intervals and lfsr that behave like unshrunk estimates;
with shrink = true they are also passed to ash, as in Gerard & Stephens: ash_mean, ash_lfsr and ash_svalue
are the shrunk estimates and the error rates to use for discoveries
"""
function summarize(fit::RUVBFit; level::Real = 0.95, shrink::Bool = true)
    k₂, p, S, C = size(fit.β)
    draws = reshape(fit.β, k₂, p, S * C)
    α = (1 - level) / 2
    stat(f) = [isnan(draws[i, j, 1]) ? NaN : f(view(draws, i, j, :)) for i in 1:k₂, j in 1:p]
    #convergence: rank-normalised split-R̂ and bulk ESS from the draws of each chain
    conv(f) = [isnan(fit.β[i, j, 1, 1]) ? NaN : f(view(fit.β, i, j, :, :)) for i in 1:k₂, j in 1:p]
    out = (
        mean = stat(mean), sd = stat(std), lower = stat(x -> quantile(x, α)), upper = stat(x -> quantile(x, 1 - α)),
        lfsr = stat(lfsr), rhat = conv(rhat), ess = conv(ess),
    )
    shrink || return out
    ash_mean, ash_lfsr, ash_svalue = fill(NaN, k₂, p), fill(NaN, k₂, p), fill(NaN, k₂, p)
    est = setdiff(1:p, fit.controls)
    for i in 1:k₂
        a = ash(out.mean[i, est], out.sd[i, est])
        ash_mean[i, est] .= a.mean
        ash_lfsr[i, est] .= a.lfsr
        ash_svalue[i, est] .= a.svalue
    end
    return merge(out, (ash_mean = ash_mean, ash_lfsr = ash_lfsr, ash_svalue = ash_svalue))
end
