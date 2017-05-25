@testset "Example 2" begin
    @polyvar x[1:3]
    p = sum(x .* x)
    @test p == x[1]^2 + x[2]^2 + x[3]^2
    @test subs(p, x[1]=>2, x[3]=>3) == x[2]^2 + 13
    A = [1.0 2 3; 4 5 6; 7 8 9]
    @test p(x=>A*vec(x)) == (x[1] + 2x[2] + 3x[3])^2 + (4x[1] + 5x[2] + 6x[3])^2 + (7x[1] + 8x[2] + 9x[3])^2
end
