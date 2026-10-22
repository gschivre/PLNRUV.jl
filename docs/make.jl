# Build the manual:
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'   (once; PLNRUV is a [sources] path dependency)
#   julia --project=docs docs/make.jl
# The pages are in docs/src; the manual is written to docs/build/index.html.

using Documenter, PLNRUV

DocMeta.setdocmeta!(PLNRUV, :DocTestSetup, :(using PLNRUV); recursive = true)

makedocs(;
    sitename = "PLNRUV.jl",
    modules = [PLNRUV],
    authors = "Geoffrey Schivre",
    remotes = nothing, #not a git repository yet: no source links
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true", mathengine = Documenter.MathJax3(),
        size_threshold_ignore = ["api.md"]
    ),
    pages = [
        "Home" => "index.md",
        "Method" => [
            "The problem" => "problem.md",
            "The model" => "model.md",
            "Separating the confounders" => "confounders.md",
            "Empirical Bayes" => "empirical_bayes.md",
            "Fitting the PLN layer" => "pln.md",
            "The fully Bayesian branch" => "bayes.md",
        ],
        "Simulating by binomial thinning" => "simulation.md",
        "Validation" => "validation.md",
        "Assumptions and limits" => "limits.md",
        "Notation and references" => "notation.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)
