@testset "Promotion" begin
    Mod.@polyvar x y
    @inferred x*y+x
    @test [x, x*y+x, x] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([1, x/y, x]) <: RationalPoly{<:AbstractTerm{Int}, <:AbstractTerm{Int}}
    @test [(x^2-2x+2) x; x x^2] isa Matrix{<:AbstractPolynomial{Int}}
    @test [2.0x, 3x] isa Vector{<:AbstractTerm{Float64}}
    @inferred Any[x*y, x+y]
    @test Any[x*y, x+y] isa Vector{Any}
    @test [x*y, x+y] isa Vector{<:AbstractPolynomial{Int}}
    @test [2x*y, x+y] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([2.0x, x/y]) <: RationalPoly{<:AbstractTerm{Float64}, <:AbstractTerm{Int}}
    @test eltype([2.0x, x/y, 1y]) <: RationalPoly{<:AbstractTerm{Float64}, <:AbstractTerm{Int}}
    @test eltype([2x+y, x/2.0y, x+1y]) <: RationalPoly{<:AbstractPolynomial{Int}, <:AbstractTerm{Float64}}
    @test eltype([-1/x^2, 1]) <: RationalPoly{<:AbstractTerm{Int}, <:AbstractTerm{Int}}

    X = [x, y]
    Y = [1 2; 3 4] * X
    @test Y[1] == x + 2y
    @test Y[2] == 3x + 4y
    Y = X' * [1 2; 3 4]
    @test Y[1] == x + 3y
    @test Y[2] == 2x + 4y
    @test dot(X, [1 2; 3 4] * X) == x^2 + 5x*y + 4y^2
end

struct A end
Base.zero(::Type{A}) = A()
Base.:*(::A, ::A) = B()
struct B end
Base.zero(::Type{B}) = B()
Base.:+(::B, ::B) = C()
struct C end
Base.:+(::B, ::C) = C()
Base.:+(::C, ::B) = C()
Base.:+(::C, ::C) = C()
@testset "promote_operation with tricky types" begin
    Mod.@polyvar x
    @test MA.promote_operation(*, polynomialtype(x, A), polynomialtype(x, A)) == polynomialtype(x, C)
end
