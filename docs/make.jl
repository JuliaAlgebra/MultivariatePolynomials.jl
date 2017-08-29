using Documenter, MultivariatePolynomials

makedocs(
    format = :html,
    sitename = "MultivariatePolynomials",
    pages = [
        "Introduction" => "index.md",
        "Types" => "types.md",
        "Substitution" => "substitution.md",
        "Differentiation" => "differentiation.md",
        "Division" => "division.md",
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
