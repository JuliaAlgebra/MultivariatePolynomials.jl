facts("Algebraic set") do
    @polyvar x y
    V = AlgebraicSet{Int}()
    @fact_throws ArgumentError addinequality!(V, x*y)
    S = BasicSemialgebraicSet{Int}()
    addequality!(S, x)
    addinequality!(S, x^2*y)
    addequality!(S, 2*x^2*y)
    @fact typeof(Int32(2)*x^2*y) --> Term{Int32}
    addequality!(S, Int32(2)*x^2*y)
    addinequality!(S, 1.0x^2*y)
    addequality!(S, (6//3)*x^2*y+y)
    @fact_throws InexactError addinequality!(S, 1.5x^2*y)
    S = BasicSemialgebraicSet{Float64}(S)
    @fact typeof(S) --> BasicSemialgebraicSet{Float64}
    addinequality!(S, 1.5x+y)
end
