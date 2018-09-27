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
    repo   = "github.com/JuliaAlgebra/MultivariatePolynomials.jl.git",
    target = "build",
    osname = "linux",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing,
)
