# Build the manual:
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'   (once; PLNRUV is a [sources] path dependency)
#   julia --project=docs docs/make.jl
# The pages are in docs/src; the manual is written to docs/build/index.html.
# Deployment to GitHub Pages runs only when DOCS_DEPLOY=true (a repository variable, for when Pages is available).

using Documenter, DocumenterCitations, PLNRUV

DocMeta.setdocmeta!(PLNRUV, :DocTestSetup, :(using PLNRUV); recursive = true)

# References live in docs/src/refs.bib; cite them in the pages with [key](@cite) or [key](@citet).
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

makedocs(;
    plugins = [bib],
    sitename = "PLNRUV.jl",
    modules = [PLNRUV],
    authors = "Geoffrey Schivre",
    repo = Remotes.GitHub("gschivre", "PLNRUV.jl"),
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

if get(ENV, "DOCS_DEPLOY", "false") == "true"
    deploydocs(; repo = "github.com/gschivre/PLNRUV.jl.git", devbranch = "main")
end
