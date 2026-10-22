"""
    PLNRUV

Poisson-lognormal differential expression with confounder adjustment, combining CATE of Wang J. et al.,
Ann. Stat. (2017): https://doi.org/10.1214/16-AOS1511 and the RUVB / MOUTHWASH generalisations of
Gerard D. & Stephens M. The compositional offset of sequencing data is a latent factor with loadings fixed to one.
"""
module PLNRUV

using Distributions, LinearAlgebra, Random, SpecialFunctions, Statistics

export read_gct, filter_genes
export surrogate_factors, assign_groups, thin_counts, thin_2group
export log_ratio
export ruvb, summarize, lfsr
export ash, mouthwash, factor_analysis
export pln, pln_mouthwash, pln_ruvb
export squeeze_var, vash

include("io.jl")
include("simulation.jl")
include("transforms.jl")
include("rotation.jl")
include("ruvb.jl")
include("ash.jl")
include("varshrink.jl")
include("mouthwash.jl")
include("pln.jl")
include("plnruvb.jl")

end
