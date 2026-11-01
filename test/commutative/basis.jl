if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

import StarAlgebras as SA

# Basis indexed by monomials other than the monomial basis, to test the methods
# that are shared by the elements of all bases.
struct TestBasis <: MP.AbstractMonomialIndexed end

@testset "Basis elements" begin
    Mod.@polyvar x y z
    mono = x^2 * y
    @test (@inferred monomial_type(mono)) == typeof(mono)
    @test (@inferred monomial_type(typeof(mono))) == typeof(mono)
    p = MP.BasisElement{TestBasis}(mono)
    @test monomial(p) == mono
    @test monomial_type(p) == typeof(mono)
    @test monomial_type(typeof(p)) == typeof(mono)
    for f in (SA.Term(2, p), polynomial(p))
        @test (@inferred monomial_type(typeof(f))) == typeof(mono)
        @test (@inferred monomial_type(f)) == typeof(mono)
        @test (@inferred constant_monomial(typeof(f))) ==
              constant_monomial(mono)
        @test (@inferred variables([f])) == variables(mono)
        @test (@inferred variables(typeof(f)[])) ==
              variables(constant_monomial(typeof(mono)))
    end
    @test hash(p) == hash(MP.BasisElement{TestBasis}(x^2 * y))
    @test hash(p) != hash(mono)
    q = MP.BasisElement{TestBasis}(mono * z^0)
    @test q == p
    @test hash(q) == hash(p)
    @test !isconstant(p)
    @test isconstant(MP.BasisElement{TestBasis}(constant_monomial(mono)))
    @test sprint(print, p) == "TestBasis(x^2*y)"
    @test adjoint(p) == p
    @test MP.left_constant_mult(2, p) == SA.Term(2, p)
    @test MP.right_constant_mult(p, 2) == SA.Term(2, p)
    for z in [zero(p), zero(typeof(p))]
        @test iszero(z)
        @test SA.basis(z) isa MP.FullBasis{TestBasis}
    end
end
