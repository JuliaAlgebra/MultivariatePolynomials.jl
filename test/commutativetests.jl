for file in readdir(joinpath(@__DIR__, "commutative"))
    if file == "complex.jl" && !isdefined(Mod, Symbol("@complex_polyvar"))
        continue
    end
    if file == "mutable_arithmetics.jl"
        continue
    end
    include(joinpath(@__DIR__, "commutative", file))
end
