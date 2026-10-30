import StarAlgebras as SA

# Basis indexed by monomials other than the monomial basis, to test the methods
# that are shared by the elements of all bases.
struct TestBasis <: MP.AbstractMonomialIndexed end

@testset "Basis elements" begin
    Mod.@polyvar x y
    mono = x^2 * y
    p = MP.Polynomial{TestBasis}(mono)
    @test monomial(p) == mono
    @test monomial_type(p) == typeof(mono)
    @test monomial_type(typeof(p)) == typeof(mono)
    @test hash(p) == hash(MP.Polynomial{TestBasis}(x^2 * y))
    @test hash(p) != hash(mono)
    @test !isconstant(p)
    @test isconstant(MP.Polynomial{TestBasis}(constant_monomial(mono)))
    @test sprint(print, p) == "TestBasis(x^2*y)"
    @test adjoint(p) == p
    @test MP.left_constant_mult(2, p) == SA.Term(2, p)
    @test MP.right_constant_mult(p, 2) == SA.Term(2, p)
    for z in [zero(p), zero(typeof(p))]
        @test iszero(z)
        @test SA.basis(z) isa MP.FullBasis{TestBasis}
    end
end
