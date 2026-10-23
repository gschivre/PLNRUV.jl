using PLNRUV, Aqua, Test, Random, Distributions, Statistics

@testset "PLNRUV" begin
    @testset "read_gct" begin
        file = tempname() * ".gct"
        write(file, "#1.3\n2\t3\t0\t0\nid\tName\tDescription\tS1\tS2\tS3\n1\tG1\tgene one\t0\t5\t7\n2\tG2\tgene two\t1\t2.0\t3\n")
        d = read_gct(file)
        @test d.counts == [0 5 7; 1 2 3]
        @test d.genes == ["G1", "G2"]
        @test d.samples == ["S1", "S2", "S3"]
        write(file, "#1.2\n1\t2\nName\tDescription\tS1\tS2\nG1\tgene one\t4\t0\n")
        @test read_gct(file).counts == [4 0]
        @test filter_genes([0 0; 1 3; 10 20]; min_mean = 2.5) == [3]
    end

    @testset "thinning without confounding" begin
        rng = Xoshiro(1)
        counts = fill(1_000_000, 200, 20) #large counts: thinned / original ≈ thinning probability
        sim = thin_2group(rng, counts; ngenes = 150, nsamp = 16, nctrl = 30, prop_null = 0.5, signal = Normal(0, 2))
        @test size(sim.counts) == (16, 150)
        @test sum(sim.X[:, 2]) == 8
        @test length(sim.controls) == 30 && all(iszero, sim.β[sim.controls])
        @test length(sim.signal) == 60 && all(!iszero, sim.β[sim.signal])
        @test all(iszero, sim.β[sim.null])
        #null genes are untouched
        @test all(==(1_000_000), sim.counts[:, sim.null])
        #log2 ratio group 1 / group 0 recovers β
        g1 = sim.X[:, 2] .== 1
        lfc = [log2(mean(c[g1]) / mean(c[.!g1])) for c in eachcol(sim.counts)]
        @test maximum(abs, lfc - sim.β) < 0.01
        #approximately null controls
        sim2 = thin_2group(rng, counts; nctrl = 50, ctrl_signal = Normal(0, 0.1))
        @test all(!iszero, sim2.β[sim2.controls])
    end

    @testset "confounded group assignment" begin
        rng = Xoshiro(2)
        n, p = 400, 300
        z = randn(rng, n) #one strong latent factor
        counts = [rand(rng, Poisson(exp(4 + 1.5 * z[i] * g))) for g in randn(rng, p), i in 1:n]
        cors = [mean(thin_2group(rng, counts; nsamp = n, target_cor = [0.5]).realized_cor[1]) for _ in 1:20]
        @test abs(mean(abs.(cors)) - 0.5) < 0.1
        sim = thin_2group(rng, counts; nsamp = n, target_cor = [0.0])
        @test abs(sim.realized_cor[1]) < 0.2
        #sequencing depth is removed from the surrogate variables
        depth = exp.(randn(rng, n))
        deep = [rand(rng, Poisson(10 * depth[i])) for _ in 1:p, i in 1:n]
        @test abs(cor(surrogate_factors(deep, 1)[:, 1], log.(vec(sum(deep; dims = 1))))) < 0.5
    end

    @testset "log_ratio" begin
        Z = log_ratio([1 3 7; 0 2 4], [1, 2]; pseudocount = 1)
        @test Z ≈ log.([2 4 8; 1 3 5]) .- log.([sqrt(8), sqrt(3)])
    end

    @testset "ash" begin
        rng = Xoshiro(5)
        p = 2_000
        β = (rand(rng, p) .< 0.2) .* randn(rng, p) .* 2
        s = fill(0.5, p)
        r = β .+ s .* randn(rng, p)
        a = ash(r, s)
        @test sum(a.π) ≈ 1
        @test 0.7 < a.π[1] < 0.9
        @test sqrt(mean(abs2, a.mean - β)) < 0.75 * sqrt(mean(abs2, r - β))
        @test all(0 .<= a.lfsr .<= 1) && all(a.lfdr .<= a.lfsr .+ 1.0e-12)
        @test mean(β[a.svalue .< 0.05] .== 0) < 0.1
        #pure noise: everything shrunk to zero
        z = ash(randn(rng, 500), ones(500))
        @test z.π[1] > 0.95 && maximum(abs, z.mean) < 0.1
        #sparse symmetric effects: the gate keeps the prior unimodal at zero (plain ash)
        @test all(iszero, a.μ)
    end

    @testset "ash with a shifted mode (global change of transcription)" begin
        rng = Xoshiro(11)
        p = 3_000
        β = (rand(rng, p) .< 0.6) .* (1 .+ 0.3 .* randn(rng, p)) #40 % null, 60 % shifted around 1
        s = fill(0.3, p)
        r = β .+ s .* randn(rng, p)
        fdp(f) = mean(β[f.svalue .< 0.05] .== 0)
        auto, zero = ash(r, s), ash(r, s; mode = :zero)
        @test abs(auto.μ[end] - 1) < 0.1 && auto.μ[1] == 0 && auto.σ[1] == 0 #null kept at zero, mode at the shift
        @test abs(auto.π[1] - 0.4) < 0.05 #π₀ recovered
        @test zero.π[1] < auto.π[1] - 0.05 #the unimodal prior takes null mass for the shifted effects
        @test fdp(auto) < 0.06 && fdp(auto) < fdp(zero)
        @test PLNRUV.mode_gain(r, PLNRUV.NoiseModel(s .^ 2), 1.0, 1:p, PLNRUV.ash_grid(r, s), 10) > 0.05
        #MOUTHWASH: the mode is estimated with controls, not without them (it would be the compositional offset)
        c = 300
        βc = [zeros(c); β[1:(p - c)]]
        rc = βc .+ s .* randn(rng, p)
        Γ = 0.1 .* randn(rng, 1, p)
        fc = PLNRUV.mouthwash_fit(rc, s, Γ, 1:c)
        @test abs(fc.g.mode - 1) < 0.1
        @test PLNRUV.mouthwash_fit(rc, s, Γ, Int[]).g.mode == 0
    end

    @testset "compound e-values and e-BH" begin
        rng = Xoshiro(13)
        p = 5_000
        β = (rand(rng, p) .< 0.2) .* randn(rng, p) .* 3
        s = fill(0.5, p)
        r = β .+ s .* randn(rng, p)
        a = ash(r, s)
        e = evalues(r, s .^ 2, a.σ, a.π)
        @test all(>=(0), e)
        #e = π₀ / (1 - π₀) × (1 - lfdr) / lfdr
        @test e ≈ a.π[1] / (1 - a.π[1]) .* (1 .- a.lfdr) ./ a.lfdr rtol = 1.0e-6
        #under the null (a fixed prior) the e-values average to at most one
        null = evalues(0.5 .* randn(rng, 200_000), fill(0.25, 200_000), a.σ, a.π)
        @test mean(null) < 1.05
        disc = ebh(e; level = 0.05)
        @test count(disc) > 100 && mean(β[disc] .== 0) < 0.05
        @test ebh([100.0, 1.0, 1.0, 1.0]; level = 0.1) == [true, false, false, false] #1 × 100 / 4 ≥ 10
        @test !any(ebh(ones(10)))
    end

    @testset "variance shrinkage" begin
        rng = Xoshiro(8)
        p, d = 3_000, 6
        #limma: true variances from a scaled inverse χ² prior, so d₀ and s₀² are recovered
        σ2 = 0.5 .* 10 ./ rand(rng, Chisq(10), p)
        s2 = σ2 .* rand(rng, Chisq(d), p) ./ d
        lim = squeeze_var(s2, d)
        @test 6 < lim.d0 < 16 && 0.35 < lim.s02 < 0.65
        @test mean(abs2, log.(lim.s2) - log.(σ2)) < mean(abs2, log.(s2) - log.(σ2))
        @test isapprox(PLNRUV.trigamma_inverse(PLNRUV.trigamma(3.7)), 3.7; rtol = 1.0e-8)
        #vash: bimodal true variances, where a single inverse-χ² prior is misspecified
        σ2 = [rand(rng) < 0.7 ? 0.05 : 1.0 for _ in 1:p] .* exp.(0.2 .* randn(rng, p))
        s2 = σ2 .* rand(rng, Chisq(d), p) ./ d
        err(x) = mean(abs2, log.(x) - log.(σ2))
        va = vash(s2, d)
        @test sum(va.π) ≈ 1
        @test err(va.s2) < err(squeeze_var(s2, d).s2) < err(s2)
        #trend: variances decreasing with a covariate (mean–dispersion trend), an exchangeable prior biases both ends
        x = 10 .* rand(rng, p)
        σ2 = exp.(-0.5 .* x .+ 0.3 .* randn(rng, p))
        s2 = σ2 .* rand(rng, Chisq(d), p) ./ d
        tr = PLNRUV.shrink_variances(s2, d, :vash; covariate = x)
        ex = PLNRUV.shrink_variances(s2, d, :vash)
        @test err(tr.s2) < err(ex.s2)
        low = x .< 2 #high-variance end: the exchangeable prior pulls them down, the trend does not
        @test median(ex.s2[low] ./ σ2[low]) < median(tr.s2[low] ./ σ2[low])
        @test PLNRUV.shrink_variances(s2, d, :none; covariate = x).s2 == s2
    end

    @testset "MOUTHWASH with confounders and a compositional shift" begin
        rng = Xoshiro(6)
        n, p, c = 40, 600, 60
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        Z = X * [0.0 0.0; 1.2 -0.8] + randn(rng, n, 2)
        Y = X * [randn(rng, 1, p); β'] + Z * randn(rng, 2, p) + 0.3 .* randn(rng, n, p)
        Y .+= X[:, 2] .+ 0.3 .* randn(rng, n)
        nc = (c + 1):p
        for ctrl in (1:c, Int[])
            mw = mouthwash(Y, X, 2; controls = ctrl)
            @test abs(mw.α[1] - 1) < 0.2 #compositional shift recovered
            @test abs(mean(mw.mean[nc] - β[nc])) < 0.05
            @test cor(mw.mean[nc], β[nc]) > 0.95
            disc = nc[mw.svalue[nc] .< 0.05]
            @test length(disc) > 50 && mean(β[disc] .== 0) < 0.1
        end
        @test all(iszero, mouthwash(Y, X, 2; controls = 1:c).mean[1:c])
    end

    @testset "MOUTHWASH: additive variance inflation and mixture of t noise" begin
        @test PLNRUV.inflated_df(10.0, 2.0) == 28.0 #4 + (10 - 4) × 2²
        @test PLNRUV.inflated_df(10.0, 0.8) == 10.0 && PLNRUV.inflated_df(Inf, 2.0) == Inf
        #estimates with t(ν) noise of scale v (from the variance estimate) plus normal noise (λ - 1) v it does not see
        p, c, ν, λ = 4_000, 300, 10.0, 2.0
        function addsim(seed, df)
            rng = Xoshiro(seed)
            v = rand(rng, Gamma(4, 0.05), p)
            β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c)]
            r = β .+ sqrt.(v) .* (rand(rng, TDist(ν), p) .+ sqrt(λ - 1) .* randn(rng, p))
            Γ = 0.1 .* randn(rng, 1, p) #a confounder with no effect (α = 0)
            fit = PLNRUV.mouthwash_fit(r, sqrt.(v), Γ, 1:c; df = df)
            nc = (c + 1):p
            disc = nc[PLNRUV.cumulative_rate(fit.lfsr[nc]) .< 0.05]
            wrong = (β[disc] .== 0) .| (sign.(β[disc]) .!= sign.(fit.mean[disc]))
            return (λ = fit.λ, ndisc = length(disc), fdp = isempty(disc) ? 0.0 : mean(wrong), v = v)
        end
        add = [addsim(s, ν) for s in 1:4]
        nor = [addsim(s, Inf) for s in 1:4]
        @test all(1.6 < a.λ < 2.4 for a in add)
        @test all(a.ndisc > 80 for a in add)
        @test mean(a.fdp for a in add) < min(0.08, mean(a.fdp for a in nor)) #the normal likelihood misses the t tails
        v = add[1].v
        #a one-component mixture is the single t
        N1 = PLNRUV.NoiseModel(v, ν)
        Nm = PLNRUV.NoiseModel([0.5 0.5] .* ones(p), [v v], fill(ν, p, 2))
        @test PLNRUV.lognoise(N1, 1, 0.3, 0.1, 1.5) ≈ PLNRUV.lognoise(Nm, 1, 0.3, 0.1, 1.5)
        #vash posterior mixture: weights sum to one and a mixture likelihood runs through MOUTHWASH
        rng = Xoshiro(9)
        s2 = v .* rand(rng, Chisq(8), p) ./ 8
        va = vash(s2, 8)
        @test all(sum(va.mix.w; dims = 2) .≈ 1) && size(va.mix.w) == size(va.mix.s2) == size(va.mix.df)
        r = sqrt.(v) .* randn(rng, p)
        Γ = 0.1 .* randn(rng, 1, p)
        fm = PLNRUV.mouthwash_fit(r, sqrt.(va.s2), Γ, 1:c; noise = PLNRUV.NoiseModel(va.mix.w, va.mix.s2, va.mix.df))
        @test isfinite(fm.loglik) && 0.5 < fm.λ < 4
    end

    @testset "PLN variational E-step" begin
        for t in (-30.0, -1.0, 0.0, 2.5, 40.0)
            w = PLNRUV.lambertw_exp(t)
            @test isapprox(w + log(w), t; atol = 1.0e-12)
        end
        @test isfinite(PLNRUV.lambertw_exp(800.0)) #exp(800) overflows
        M, S2 = PLNRUV.estep_cell(7, 1.2, -0.4, 0.8, 0.0, 0.5; inner = 50)
        f(m, s2) = 7 * (1.2 + m) - exp(1.2 + m + s2 / 2) - ((m + 0.4)^2 + s2) / (2 * 0.8) + log(s2) / 2
        @test f(M, S2) >= maximum(f(M + δ, S2) for δ in -0.01:0.001:0.01)
        @test f(M, S2) >= maximum(f(M, S2 * e) for e in 0.9:0.01:1.1)
    end

    @testset "PLN with many factors on near-empty genes" begin
        #with few residual degrees of freedom (d = 20 - 2 - 10 - 1 = 7), genes whose cells carry no information
        #used to see σ² multiplied by n / d at every M-step until the fit broke
        rng = Xoshiro(12)
        n, p = 20, 300
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        A = -2.5 .+ 0.5 .* randn(rng, 1, p) .+ 0.4 .* randn(rng, n, p)
        Y = [rand(rng, Poisson(exp(A[i, j]))) for i in 1:n, j in 1:p]
        Y[:, 1:5] .= 0
        Y[1, 1:5] .= 1 #genes with a single read
        fit = pln(Y, X, 10)
        @test all(isfinite, fit.M) && all(isfinite, fit.σ2)
        @test maximum(fit.σ2) < 100
    end

    @testset "PLN + MOUTHWASH on low counts" begin
        rng = Xoshiro(7)
        n, p, c = 30, 800, 100
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        σ = sqrt.(rand(rng, Gamma(4, 0.05), p))
        A = X * [randn(rng, 1, p) .* 1.5; β'] + randn(rng, n, p) .* σ' .+ X[:, 2] #shift, no confounders
        o = log(2.0) .+ 0.4 .* randn(rng, n)
        Y = [rand(rng, Poisson(exp(o[i] + A[i, j]))) for i in 1:n, j in 1:p]
        fit = pln_mouthwash(Y, X, 0; controls = 1:c)
        mw = fit.mouthwash
        nc = (c + 1):p
        @test fit.pln.iterations < 500
        @test 0.7 < mw.λ < 1.5 #the latent-scale standard errors need little inflation
        @test cor(mw.mean[nc], β[nc]) > 0.85
        disc = nc[mw.svalue[nc] .< 0.05]
        @test length(disc) > 50 && mean(β[disc] .== 0) < 0.1
    end

    @testset "number of factors by count splitting" begin
        x, w = PLNRUV.gausshermite(10)
        @test sum(w) ≈ sqrt(π) && sum(w .* x .^ 2) ≈ sqrt(π) / 2 #∫ exp(-x²), ∫ x² exp(-x²)
        rng = Xoshiro(10)
        n, p = 20, 400
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        W = randn(rng, n, 2)
        W .-= X * (X \ W)
        A = X * [randn(rng, 1, p) .+ 3; zeros(1, p)] + W * (0.5 .* randn(rng, 2, p)) + 0.3 .* randn(rng, n, p)
        Y = [rand(rng, Poisson(exp(A[i, j]))) for i in 1:n, j in 1:p]
        #counts around 20: the latent noise dominates the Poisson noise, binomial thinning leaves it shared by the two
        #parts and tends to choose too many factors; the beta-binomial split also splits it
        c = choose_q(Y, X; qs = 0:4, nsplit = 1, rng = Xoshiro(1))
        @test c.qs == 0:4
        @test c.q == 2
        @test last(c.path) == 2 && length(c.path) <= 5
        cb = choose_q(Y, X; qs = 0:4, nsplit = 1, split = :binomial, rng = Xoshiro(1))
        @test cb.q >= 2 && length(cb.path) == 1
        #the beta-binomial split: uncorrelated parts at r = 1 / CV² of the latent rate
        rng = Xoshiro(11)
        σ2, m = 0.5, 20.0
        r = 1 / expm1(σ2)
        y = [rand(rng, Poisson(m * exp(sqrt(σ2) * randn(rng) - σ2 / 2))) for _ in 1:200_000]
        y1, y2 = PLNRUV.split_counts(rng, reshape(y, :, 1), [r], 0.5)
        @test abs(cor(vec(y1), vec(y2))) < 0.01
        y1b, y2b = PLNRUV.split_counts(rng, reshape(y, :, 1), [Inf], 0.5)
        @test cor(vec(y1b), vec(y2b)) > 0.5
        #3 samples per group and genes that are on in some samples and off in others: some fits with factors diverge
        #on the split, those q are rejected instead of stopping the choice
        rng = Xoshiro(1)
        n, p = 6, 1_000
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        W = randn(rng, n, 2)
        W .-= X * (X \ W)
        A = X * [randn(rng, 1, p) .* 0.5 .+ 4; zeros(1, p)] + W * (0.5 .* randn(rng, 2, p)) + 0.3 .* randn(rng, n, p)
        A[:, 1:20] .= 3 .+ 4 .* randn(rng, n, 20)
        Y = [rand(rng, Poisson(exp(A[i, j]))) for i in 1:n, j in 1:p]
        c = choose_q(Y, X; qs = 0:2, nsplit = 1, rng = Xoshiro(1))
        @test c.q in 0:2
    end

    @testset "warm-started vash" begin
        rng = Xoshiro(13)
        p, d = 2_000, 8
        s2 = rand(rng, Gamma(2, 0.1), p) .* rand(rng, Chisq(d), p) ./ d
        v0 = vash(s2, d)
        s2b = s2 .* exp.(0.05 .* randn(rng, p)) #the variances of the next EM iteration
        cold = vash(s2b, d)
        warm = vash(s2b, d; warm = (c = v0.c, π = v0.π))
        @test maximum(abs.(warm.s2 ./ cold.s2 .- 1)) < 1.0e-3
        @test isapprox(warm.c, cold.c; rtol = 1.0e-3)
    end

    @testset "display" begin
        rng = Xoshiro(11)
        n, p, c = 30, 200, 40
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        Y = [rand(rng, Poisson(exp(log(5.0) + X[i, 2] * β[j] + 0.3 * randn(rng)))) for i in 1:n, j in 1:p]
        mw = pln_mouthwash(Y, X, 1; controls = 1:c).mouthwash
        txt = sprint(show, MIME"text/plain"(), mw)
        @test occursin("Confounders", txt) && occursin("Discoveries", txt)
    end

    @testset "Code quality (Aqua)" begin
        Aqua.test_all(PLNRUV)
    end
end
