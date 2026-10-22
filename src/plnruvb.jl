#fully Bayesian branch: RUVB (Gerard & Stephens) on the latent layer of the Poisson-lognormal model
#    yᵢⱼ ~ Poisson(exp(oᵢ + Aᵢⱼ)),   Q'A = rotated latent log-abundances following the RUVB factor model
#one sweep: (1) the RUVB step on the rotated latent matrix (factors, variances, imputation of the missing block),
#(2) the mean of A given the RUVB state, rotated back, (3) each Aᵢⱼ by an independence Metropolis-Hastings step
#whose Gaussian proposal is centred at the mode of its full conditional (Lambert W) with the local curvature

#full conditional of one latent value: log π(a) = y a - exp(o + a) - (a - μ)² / (2σ²)
latent_logdensity(a, y, o, μ, σ2) = y * a - exp(o + a) - (a - μ)^2 / (2σ2)

#independence Metropolis-Hastings update of every latent value; return the number of accepted moves
function mh_latent!(rng::AbstractRNG, A, Y, o, μ, σ2; scale = 1.2)
    accepted = 0
    for j in axes(A, 2), i in axes(A, 1)
        y, oi, m, s2 = Y[i, j], o[i], μ[i, j], σ2[j]
        mode = m + s2 * y - lambertw_exp(log(s2) + oi + m + s2 * y)
        sd = sqrt(scale / (exp(oi + mode) + 1 / s2))
        prop = mode + sd * randn(rng)
        logr = latent_logdensity(prop, y, oi, m, s2) - latent_logdensity(A[i, j], y, oi, m, s2) +
            ((prop - mode)^2 - (A[i, j] - mode)^2) / (2sd^2)
        if log(rand(rng)) < logr
            A[i, j] = prop
            accepted += 1
        end
    end
    return accepted
end

"""
    pln_ruvb(Y, X, controls, q; cov_of_interest = [size(X, 2)], offset_factor = true, offsets = log library sizes, a2 = 1, b2 = 1, ν = 2, τ = 2, burnin = 2500, nsample = 500, thin = 10, nchain = 4, rng)

fully Bayesian RUVB on counts Y (samples x genes) through the Poisson-lognormal latent layer
the chains start from the variational PLN fit (empirical Bayes branch); return an RUVBFit (see summarize) that
also records the acceptance rate of the latent Metropolis-Hastings updates
"""
function pln_ruvb(
        Y::AbstractMatrix{<:Integer}, X::AbstractMatrix{<:Real}, controls::AbstractVector{<:Integer}, q::Integer;
        cov_of_interest::AbstractVector{<:Integer} = [size(X, 2)], offset_factor::Bool = true,
        offsets::AbstractVector{<:Real} = vec(log.(sum(Y; dims = 2))),
        a2::Real = 1.0, b2::Real = 1.0, ν::Real = 2, τ::Real = 2,
        burnin::Integer = 2_500, nsample::Integer = 500, thin::Integer = 10, nchain::Integer = 4,
        rng::AbstractRNG = Random.default_rng()
    )
    n, p = size(Y)
    k = size(X, 2)
    nfix = Int(offset_factor)
    @assert length(controls) >= q + nfix "need at least q + offset_factor control genes"
    @assert n - k >= q + 1 "not enough residual degrees of freedom for q factors"
    k₂ = length(cov_of_interest)
    k₁ = k - k₂
    o = collect(float.(offsets))
    noncontrols = setdiff(1:p, controls)
    nuisance = setdiff(1:k, cov_of_interest)
    QR = qr(X[:, [nuisance; cov_of_interest]])
    R22inv = inv(UpperTriangular(QR.R[(k₁ + 1):k, (k₁ + 1):k]))
    kept = [(k₁ + 1):k; (k + 1):n] #rows of Q'A used by RUVB (the nuisance rows have a flat prior)

    #starting latent values from the variational fit
    A0 = pln(Y, X, q; offset_factor = offset_factor, offsets = o).M

    β = fill(NaN, k₂, p, nsample, nchain)
    acceptance = zeros(nchain)
    seeds = rand(rng, UInt64, nchain)
    tasks = map(1:nchain) do c
        Threads.@spawn begin
            crng = Xoshiro(seeds[c])
            A = copy(A0)
            Arot = QR.Q' * A
            s = ruvb_init(crng, Arot[kept, :], k₂, q, controls, noncontrols, nfix)
            M̃ = zeros(n, p)
            accepted = 0
            for it in 1:(burnin + nsample * thin)
                #(1) RUVB step on the rotated latent values (the missing block keeps its imputed values)
                ruvb_step!(crng, s, k₂, noncontrols, nfix, a2, b2, ν, τ)
                σ2 = 1 ./ s.σ2inv
                #(2) rotated mean: flat-prior draw for the nuisance rows, R₂₂β + L F for the rows of interest,
                #L F for the residual rows; R₂₂β = observed - imputed for the non-control genes
                LF = s.L * s.F
                for j in 1:p, i in 1:k₁
                    M̃[i, j] = Arot[i, j] + sqrt(σ2[j]) * randn(crng)
                end
                M̃[kept, :] .= LF
                M̃[(k₁ + 1):k, noncontrols] .+= Arot[(k₁ + 1):k, noncontrols] .- s.Yt[1:k₂, noncontrols]
                if it > burnin && (it - burnin) % thin == 0
                    t = (it - burnin) ÷ thin
                    β[:, noncontrols, t, c] .= R22inv * (Arot[(k₁ + 1):k, noncontrols] .- s.Yt[1:k₂, noncontrols])
                end
                #(3) latent values given their mean, then the observed part of the RUVB data
                accepted += mh_latent!(crng, A, Y, o, QR.Q * M̃, σ2)
                Arot .= QR.Q' * A
                s.Yt[1:k₂, controls] .= Arot[(k₁ + 1):k, controls]
                s.Yt[(k₂ + 1):end, :] .= Arot[(k + 1):n, :]
            end
            acceptance[c] = accepted / ((burnin + nsample * thin) * n * p)
        end
    end
    foreach(wait, tasks)
    return RUVBFit(β, collect(controls), collect(cov_of_interest), acceptance)
end
