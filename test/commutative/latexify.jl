using Latexify

@testset "Latexify Extension" begin
    
    Mod.@polyvar y[1:10]

    #test no error is thrown
    @test latexify(y[1]) isa Latexify.LaTeXStrings.LaTeXString #variable
    @test latexify(prod(y)) isa Latexify.LaTeXStrings.LaTeXString #monomial
    @test latexify(sqrt(3) * prod(y)) isa Latexify.LaTeXStrings.LaTeXString #term
    @test latexify(sum(sqrt(i) * y[i]^i for i=1:10)) isa Latexify.LaTeXStrings.LaTeXString #polynomial
end

