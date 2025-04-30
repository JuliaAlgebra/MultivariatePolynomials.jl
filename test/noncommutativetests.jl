for file in readdir(joinpath(@__DIR__, "noncommutative"))
    include(joinpath(@__DIR__, "noncommutative", file))
end
