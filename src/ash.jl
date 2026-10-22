#adaptive shrinkage of Stephens M., Biostatistics (2017): https://doi.org/10.1093/biostatistics/kxw041
#normal means problem rⱼ ~ N(βⱼ, sⱼ²) with the prior g = Σₖ πₖ N(0, σₖ²), σ₁ = 0 being the point mass at zero
#the weights π are estimated by maximum a posteriori with a Dirichlet(nullweight, 1, ..., 1) prior (mix-SQP)

#log density of N(0, v) at x
lognormpdf(x::Real, v::Real) = -0.5 * (log(2π) + log(v) + x^2 / v)

#normalising constant of the Student t density with ν degrees of freedom (ν = Inf: the normal)
tconst(ν::Real) = isinf(ν) ? -0.5 * log(2π) : loggamma((ν + 1) / 2) - loggamma(ν / 2) - 0.5 * log(ν * π)

#log density at x of a t with ν degrees of freedom and squared scale v (normal when ν = Inf); c = tconst(ν)
logtpdf(x::Real, v::Real, ν::Real, c::Real) = isinf(ν) ? lognormpdf(x, v) : c - 0.5 * log(v) - (ν + 1) / 2 * log1p(x^2 / (ν * v))

#degrees of freedom of each estimate: a scalar (Inf for the normal likelihood) or one value per gene
gene_df(ν::Real, p) = fill(float(ν), p)
gene_df(ν::AbstractVector, p) = collect(float.(ν))

#grid of prior standard deviations: from min(s) / 10 to 2 * sqrt(max(r² - s²)), multiplicative step mult
function ash_grid(r::AbstractVector{<:Real}, s::AbstractVector{<:Real}; mult::Real = sqrt(2))
    σmin = minimum(s) / 10
    σmax = 2 * sqrt(max(maximum(r .^ 2 .- s .^ 2), 0.0))
    σmax = max(σmax, 8 * σmin)
    K = ceil(Int, log(σmax / σmin) / log(mult)) + 1
    return [0.0; σmin .* mult .^ (0:(K - 1))]
end

#component likelihoods Lⱼₖ = N(rⱼ; 0, λ sⱼ² + σₖ²) scaled by their row maximum (returned in logscale)
#with ν (degrees of freedom of each sⱼ) the normal is replaced by a t with scale √(λ sⱼ² + σₖ²), the usual
#approximation of the convolution of a t likelihood with a normal mixture (exact for the point mass)
function mixture_lik!(Lik, logscale, r, s2, σ2, λ; ν = Inf)
    νj = gene_df(ν, length(r))
    for j in eachindex(r)
        m = -Inf
        c = tconst(νj[j])
        for k in eachindex(σ2)
            Lik[j, k] = logtpdf(r[j], λ * s2[j] + σ2[k], νj[j], c)
            m = max(m, Lik[j, k])
        end
        logscale[j] = m
        for k in eachindex(σ2)
            Lik[j, k] = exp(Lik[j, k] - m)
        end
    end
    return Lik
end

#penalised log-likelihood of the mixture weights (rows in free only)
function mixture_objective(π, Lik, logscale, free, prior)
    den = Lik * π
    obj = sum(k -> (prior[k] - 1) * log(max(π[k], 1.0e-300)), eachindex(π))
    for j in free
        obj += log(max(den[j], 1.0e-300)) + logscale[j]
    end
    return obj
end

#non-negative quadratic program min ½ y'Hy - b'y, y ≥ 0 (H positive definite), by the active-set method of
#Lawson & Hanson (1974): add the variable with the largest negative gradient, solve on the passive set, step back if needed
function nnqp(H, b; maxiter = 10 * length(b))
    #symmetric diagonal scaling (unit diagonal) keeps y ≥ 0 and removes most of the ill-conditioning
    D = 1 ./ sqrt.(diag(H))
    return D .* nnqp_scaled(D .* H .* D', D .* b, maxiter)
end
function nnqp_scaled(H, b, maxiter)
    K = length(b)
    tol = 1.0e-10 * max(maximum(abs, b), 1.0)
    passive = falses(K)
    y = zeros(K)
    w = copy(b) #minus the gradient at y
    for _ in 1:maxiter
        candidates = findall(k -> !passive[k] && w[k] > tol, 1:K)
        isempty(candidates) && break
        passive[candidates[argmax(w[candidates])]] = true
        while true
            P = findall(passive)
            z = zeros(K)
            z[P] .= Symmetric(H[P, P]) \ b[P]
            if all(>(0), z[P])
                y = z
                break
            end
            #move towards z until a passive variable reaches zero, and release it
            neg = filter(k -> z[k] <= 0, P)
            α = minimum(y[k] / (y[k] - z[k]) for k in neg)
            y .+= α .* (z .- y)
            for k in P
                if y[k] <= tol
                    y[k] = 0.0
                    passive[k] = false
                end
            end
        end
        w = b .- H * y
    end
    return y
end

#mix-SQP of Kim Y., Carbonetto P., Stephens M. & Anitescu M., J. Comput. Graph. Stat. (2020): https://doi.org/10.1080/10618600.2019.1689985
#maximise F(x) = Σⱼ log(Lⱼ'x) + Σₖ (priorₖ - 1) log(xₖ) - W Σₖ xₖ over x ≥ 0 (W = number of rows + Σₖ (priorₖ - 1)),
#whose maximiser lies on the simplex; each step solves the quadratic approximation as a non-negative QP, then a line search
function mixture_weights(Lik, logscale, free, prior; π0 = fill(1 / size(Lik, 2), size(Lik, 2)), tol = 1.0e-8, maxiter = 100)
    L = Lik[collect(free), :]
    K = size(L, 2)
    a = prior .- 1.0
    W = size(L, 1) + sum(a)
    F(x) = sum(log, max.(L * x, 1.0e-300)) + sum(k -> a[k] > 0 ? a[k] * log(x[k]) : 0.0, 1:K) - W * sum(x)
    x = 0.99 .* π0 ./ sum(π0) .+ 0.01 / K #strictly positive start
    Fx = F(x)
    for _ in 1:maxiter
        u = max.(L * x, 1.0e-300)
        #gradient of F and minus its Hessian (the Dirichlet terms only where aₖ > 0: other xₖ may reach zero)
        barrier = [a[k] > 0 ? a[k] / x[k] : 0.0 for k in 1:K]
        g = L' * (1 ./ u) .+ barrier .- W
        H = L' * (L ./ u .^ 2)
        H[diagind(H)] .+= barrier ./ max.(x, 1.0e-300) .+ 1.0e-8
        d = nnqp(H, H * x .+ g) .- x
        ascent = dot(g, d) #predicted ascent: zero at a KKT point
        ascent <= tol * W && break
        #step capped so the components with a Dirichlet weight stay positive, then Armijo backtracking
        t = 1.0
        for k in 1:K
            a[k] > 0 && d[k] < 0 && (t = min(t, 0.99 * x[k] / -d[k]))
        end
        Fnew = F(x .+ t .* d)
        while Fnew < Fx + 1.0e-4 * t * ascent && t > 1.0e-12
            t /= 2
            Fnew = F(x .+ t .* d)
        end
        t <= 1.0e-12 && break
        x = x .+ t .* d
        Fx = Fnew
    end
    π = x ./ sum(x)
    return π, mixture_objective(π, Lik, logscale, free, prior)
end

#posterior summaries of βⱼ given rⱼ ~ N(βⱼ, vⱼ) and the prior Σₖ πₖ N(0, σₖ²)
function ash_posterior(r, v, σ2, π; ν = Inf)
    νj = gene_df(ν, length(r))
    p = length(r)
    out = (mean = zeros(p), sd = zeros(p), lfsr = zeros(p), lfdr = zeros(p))
    w = similar(π)
    for j in 1:p
        for k in eachindex(π)
            w[k] = π[k] == 0 ? -Inf : log(π[k]) + logtpdf(r[j], v[j] + σ2[k], νj[j], tconst(νj[j]))
        end
        w .= exp.(w .- maximum(w))
        w ./= sum(w)
        m1 = m2 = pneg = ppos = pzero = 0.0
        for k in eachindex(π)
            if σ2[k] == 0
                pzero += w[k]
            else
                vk = σ2[k] * v[j] / (σ2[k] + v[j])
                mk = σ2[k] / (σ2[k] + v[j]) * r[j]
                m1 += w[k] * mk
                m2 += w[k] * (vk + mk^2)
                pk = cdf(Normal(), -mk / sqrt(vk))
                pneg += w[k] * pk
                ppos += w[k] * (1 - pk)
            end
        end
        out.mean[j] = m1
        out.sd[j] = sqrt(max(m2 - m1^2, 0.0))
        out.lfsr[j] = min(pneg, ppos) + pzero
        out.lfdr[j] = pzero
    end
    return out
end

#q-values (from lfdr) or s-values (from lfsr): running mean of the sorted local rates
function cumulative_rate(x::AbstractVector{<:Real})
    o = sortperm(x)
    out = similar(x, Float64)
    acc = 0.0
    for (i, j) in enumerate(o)
        acc += x[j]
        out[j] = acc / i
    end
    return out
end

struct AshFit
    σ::Vector{Float64} #prior standard deviations (σ[1] = 0 is the point mass)
    π::Vector{Float64} #prior weights
    loglik::Float64 #penalised log-likelihood at the optimum
    mean::Vector{Float64}
    sd::Vector{Float64}
    lfsr::Vector{Float64}
    lfdr::Vector{Float64}
    svalue::Vector{Float64}
    qvalue::Vector{Float64}
end

"""
    ash(r, s; nullweight = 10, grid = ash_grid(r, s))

ash for estimates r with standard errors s, nullweight is the Dirichlet weight of the point mass
"""
function ash(r::AbstractVector{<:Real}, s::AbstractVector{<:Real}; nullweight::Real = 10, grid = ash_grid(r, s), df = Inf)
    @assert length(r) == length(s) && all(>(0), s)
    σ2 = grid .^ 2
    s2 = s .^ 2
    Lik = zeros(length(r), length(σ2))
    logscale = zeros(length(r))
    mixture_lik!(Lik, logscale, r, s2, σ2, 1.0; ν = df)
    prior = [nullweight; ones(length(σ2) - 1)]
    π, obj = mixture_weights(Lik, logscale, eachindex(r), prior)
    post = ash_posterior(r, s2, σ2, π; ν = df)
    return AshFit(collect(grid), π, obj, post.mean, post.sd, post.lfsr, post.lfdr, cumulative_rate(post.lfsr), cumulative_rate(post.lfdr))
end

"""
    evalues(r, v, σ, π)

compound e-values of Ramdas A. & Wang R., Hypothesis testing with e-values (2024), Theorem 9.21 / Ignatiadis N.,
Wang R. & Ramdas A., arXiv:2409.19812 (2024): for rⱼ ~ N(βⱼ, vⱼ) and the prior Σₖ πₖ N(0, σₖ²) (σ₁ = 0 the point
mass), the likelihood ratio of the non-null mixture to the null, eⱼ = qⱼ(rⱼ) / N(rⱼ; 0, vⱼ) with
qⱼ = Σₖ₌₂ πₖ N(0, vⱼ + σₖ²) / (1 - π₁), equivalently eⱼ = π₁ / (1 - π₁) × (1 - lfdrⱼ) / lfdrⱼ;
with an estimated prior they are asymptotically valid compound e-values
"""
function evalues(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, σ::AbstractVector{<:Real}, π::AbstractVector{<:Real}; ν = Inf)
    νj = gene_df(ν, length(r))
    π0 = π[1]
    e = zeros(length(r))
    for j in eachindex(r)
        c = tconst(νj[j])
        l0 = logtpdf(r[j], v[j], νj[j], c)
        num = sum(π[k] * exp(logtpdf(r[j], v[j] + σ[k]^2, νj[j], c) - l0) for k in 2:length(π); init = 0.0)
        e[j] = num / max(1 - π0, 1.0e-300)
    end
    return e
end

"""
    ebh(e; level = 0.05)

e-BH procedure of Wang R. & Ramdas A., J. R. Stat. Soc. B (2022): reject the k* hypotheses with the largest
e-values, k* = max{k : k e₍ₖ₎ / K ≥ 1 / level}; controls the FDR under arbitrary dependence for compound e-values.
Return a Boolean vector of discoveries
"""
function ebh(e::AbstractVector{<:Real}; level::Real = 0.05)
    K = length(e)
    o = sortperm(e; rev = true)
    kstar = 0
    for k in 1:K
        k * e[o[k]] / K >= 1 / level && (kstar = k)
    end
    disc = falses(K)
    disc[o[1:kstar]] .= true
    return disc
end
