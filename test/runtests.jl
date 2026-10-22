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

    @testset "RUVB with confounders and a compositional shift" begin
        rng = Xoshiro(3)
        n, p, c = 40, 200, 50
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        Z = X * [0.0 0.0; 1.2 -0.8] + randn(rng, n, 2)
        Y = X * [randn(rng, 1, p); β'] + Z * randn(rng, 2, p) + 0.3 .* randn(rng, n, p)
        Y .+= X[:, 2] .+ 0.3 .* randn(rng, n) #per-sample offset correlated with the group
        fit = ruvb(Y, X, 1:c, 2; burnin = 500, nsample = 100, thin = 5, nchain = 2, rng = Xoshiro(4))
        s = summarize(fit)
        nc = (c + 1):p
        @test size(fit.β) == (1, p, 100, 2)
        @test all(isnan, s.mean[1, 1:c])
        @test cor(s.mean[1, nc], β[nc]) > 0.9
        @test abs(mean(s.mean[1, nc] - β[nc])) < 0.1
        @test mean(s.lower[1, nc] .<= β[nc] .<= s.upper[1, nc]) > 0.85
        @test all(0 .<= s.lfsr[1, nc] .<= 0.5)
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

    @testset "PLN-RUVB on low counts" begin
        rng = Xoshiro(9)
        n, p, c = 30, 300, 60
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        σ = sqrt.(rand(rng, Gamma(4, 0.05), p))
        A = X * [randn(rng, 1, p) .* 1.5; β'] + randn(rng, n, p) .* σ' .+ X[:, 2]
        o = log(2.0) .+ 0.4 .* randn(rng, n)
        Y = [rand(rng, Poisson(exp(o[i] + A[i, j]))) for i in 1:n, j in 1:p]
        out = pln_ruvb(Y, X, 1:c, 1; burnin = 300, nsample = 100, thin = 2, nchain = 2, rng = Xoshiro(10))
        s = summarize(out)
        nc = (c + 1):p
        @test all(0.8 .< out.acceptance .< 1)
        @test all(isnan, s.mean[1, 1:c])
        @test mean(s.lower[1, nc] .<= β[nc] .<= s.upper[1, nc]) > 0.85
        @test cor(s.ash_mean[1, nc], β[nc]) > 0.85
        disc = nc[s.ash_svalue[1, nc] .< 0.05]
        @test length(disc) > 30 && mean(β[disc] .== 0) < 0.1
    end

    @testset "convergence diagnostics and display" begin
        rng = Xoshiro(11)
        iid = randn(rng, 1_000, 4)
        @test abs(rhat(iid) - 1) < 0.01
        @test 3_000 < ess(iid) < 5_000
        stuck = randn(rng, 1_000, 4) .+ [0 0 3 3]
        @test rhat(stuck) > 1.5
        ar = zeros(1_000, 4) #AR(1) with ρ = 0.9: ESS ≈ N (1 - ρ) / (1 + ρ)
        for c in 1:4, t in 2:1_000
            ar[t, c] = 0.9 * ar[t - 1, c] + randn(rng)
        end
        @test 100 < ess(ar) < 400
        #display of both branches
        n, p, c = 30, 200, 40
        X = [ones(n) repeat([0.0, 1.0], n ÷ 2)]
        β = [zeros(c); (rand(rng, p - c) .< 0.3) .* randn(rng, p - c) .* 1.5]
        Y = [rand(rng, Poisson(exp(log(5.0) + X[i, 2] * β[j] + 0.3 * randn(rng)))) for i in 1:n, j in 1:p]
        fit = pln_ruvb(Y, X, 1:c, 1; burnin = 100, nsample = 50, thin = 1, nchain = 2, rng = Xoshiro(12))
        txt = sprint(show, MIME"text/plain"(), fit)
        @test occursin("Convergence", txt) && occursin("acceptance", txt) && occursin("Top 10 genes", txt)
        @test occursin("RUVBFit(200 genes", sprint(show, fit))
        mw = pln_mouthwash(Y, X, 1; controls = 1:c).mouthwash
        txt = sprint(show, MIME"text/plain"(), mw)
        @test occursin("Confounders", txt) && occursin("Discoveries", txt)
    end

    @testset "Code quality (Aqua)" begin
        Aqua.test_all(PLNRUV)
    end
end
