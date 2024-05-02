using Documenter, MultivariatePolynomials

makedocs(
    sitename = "MultivariatePolynomials",

    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,

    pages = [
        "Introduction" => "index.md",
        "Types" => "types.md",
        "Substitution" => "substitution.md",
        "Differentiation" => "differentiation.md",
        "Division" => "division.md",
        "Internal" => "internal.md",
    ],

    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [MultivariatePolynomials],
)

deploydocs(
    repo   = "github.com/JuliaAlgebra/MultivariatePolynomials.jl.git",
)
