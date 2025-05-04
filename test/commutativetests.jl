for file in readdir(joinpath(@__DIR__, "commutative"))
    if file == "independent.jl"
        # It's included by `gcd.jl`
        continue
    end
    if file == "complex.jl" && !isdefined(Mod, Symbol("@complex_polyvar"))
        continue
    end
    include(joinpath(@__DIR__, "commutative", file))
end
