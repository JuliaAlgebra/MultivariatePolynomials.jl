using Documenter, MultivariatePolynomials

makedocs(
    format = :html,
    sitename = "MultivariatePolynomials",
    pages = [
        "Introduction" => "index.md",
        "Reference" => "apireference.md",
    ]
)

deploydocs(
    repo   = "github.com/blegat/MultivariatePolynomials.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
