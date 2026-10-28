for file in
    filter(endswith(".jl"), readdir(joinpath(@__DIR__, "noncommutative")))
    include(joinpath(@__DIR__, "noncommutative", file))
end
