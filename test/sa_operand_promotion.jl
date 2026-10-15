module TestSAOperandPromotion

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import StarAlgebras as SA
import LinearAlgebra

function same_basis_promotion_allocations(output, x)
    SA._promote_operand!(+, output, x)
    return @allocated SA._promote_operand!(+, output, x)
end

@testset "Scalar addition prepares an implicit basis" begin
    DP.@polyvar x
    p = 2x + 1
    @test (@inferred MP.implicit(p)) === p
    t = 2x
    @test (@inferred MP.implicit(t)) === t
    basis = SA.SubBasis(MP.FullBasis{MP.Monomial}([x]), [[1]])
    q = SA.AlgebraElement([2], MP.algebra(basis))
    for a in (q, SA.Term(parent(q), 1, 2))
        @test (@inferred a + 3) == 2x + 3
        @test (@inferred 3 + a) == 3 + 2x
        @test (@inferred a - 3) == 2x - 3
        @test (@inferred 3 - a) == 3 - 2x
    end
    @test SA.coeffs(q) == [2]
end

@testset "Mutable destinations promote bases" begin
    DP.@polyvar x y z
    for T in (Int, BigInt)
        a = convert(DP.Polynomial{T}, 2x + 2)
        b = convert(DP.Polynomial{T}, y + 1)
        prototype = convert(DP.Polynomial{T}, x + y + z)
        alg = parent(prototype)
        originals = MP.MA.mutable_copy.((a, b))
        same_basis = first(SA.promote_bases(a, alg))
        @test (@inferred SA._promote_operand!(+, prototype, same_basis)) ===
              same_basis
        @test same_basis_promotion_allocations(prototype, same_basis) == 0
        for op in (+, -, *)
            output = zero(prototype)
            @test (@inferred MP.MA.operate_to!(output, op, a, b)) === output
            @test output == op(a, b)
            @test parent(output) === alg
        end
        output = copy(same_basis)
        @test (@inferred MP.MA.mul!!(output, b)) === output
        @test output == a * b
        @test parent(output) === parent(same_basis)
        output = copy(same_basis)
        right = SA.AlgebraElement(SA.coeffs(output), alg)
        @test right !== output
        @test SA.coeffs(right) === SA.coeffs(output)
        @test (@inferred MP.MA.mul!!(output, right)) === output
        @test output == a * a
        @test parent(output) === parent(same_basis)
        output = copy(a)
        @test (@inferred MP.MA.operate!(*, output, b)) === output
        @test output == a * b
        @test parent(output) == parent(a * b)
        @test parent(output) !== parent(a)
        output = zero(prototype)
        @test (@inferred MP.MA.operate_to!(output, -, a)) === output
        @test output == -a
        @test (@inferred MP.MA.operate_to!(output, *, 2, a)) === output
        @test output == 2a
        for op in (*, /, //, div)
            @test (@inferred MP.MA.operate_to!(output, op, a, 2)) === output
            @test output == op(a, 2)
        end
        output = zero(prototype)
        @test (@inferred MP.MA.add!!(output, 2x)) === output
        @test output == 2x
        for op in (MP.MA.add_mul, MP.MA.sub_mul),
            (left, right) in ((a, b), (2x, b), (a, 3y), (2, b), (b, 2))

            output = zero(prototype)
            @test (@inferred MP.MA.operate_to!(output, op, a, left, right)) ===
                  output
            @test output == op(a, left, right)
            @test parent(output) === alg
            expected = op(output, left, right)
            @test (@inferred MP.MA.operate!!(op, output, left, right)) ===
                  output
            @test output == expected
        end
        for operand in (b, 3y)
            output = MP.MA.mutable_copy(a)
            @test (@inferred MP.MA.add!!(output, operand)) === output
            @test output == a + operand
            @test parent(output) == parent(a + operand)
            @test parent(output) !== parent(a)
        end
        for (left, right) in ((a, b), (2, b), (b, 2))
            output = MP.MA.mutable_copy(a)
            @test (@inferred MP.MA.operate!(
                MP.MA.add_mul,
                output,
                left,
                right,
            )) === output
            @test output == a + left * right
        end
        @test (a, b) == originals
        @test parent(a) != alg && parent(b) != alg
    end
end

@testset "Prepared zeros reuse MA reductions" begin
    DP.@polyvar x y z
    for T in (Int, BigInt)
        terms = [T(2) * x, T(3) * y]
        polys = [T(1) * y + 1, T(1) * z + 2]
        originals = MP.MA.mutable_copy.((terms, polys))
        P = DP.Polynomial{T}
        left, right = @inferred SA.promote_bases(terms[1], terms[2])
        @test parent(left) == parent(right)
        @test (left, right) == (terms[1], terms[2])
        acc = @inferred MP.MA.zero_for(P, terms, polys)
        other = @inferred MP.MA.zero_for(P, terms, polys)
        @test parent(acc) == parent(x + y + z)
        @test iszero(acc) && iszero(other)
        @test MP.MA.add!!(acc, terms[1]) === acc
        @test iszero(other)
        for a in (terms, polys)
            expected = a[1] + a[2]
            @test (@inferred sum(a)) == expected
            @test (@inferred MP.MA.operate(sum, a)) == expected
            @test parent(sum(a)) == parent(expected)
            @test sum(reshape(a, 1, 2); dims = 2) == fill(expected, 1, 1)
            seed = T(1) * (x + y + z)
            original_seed = MP.MA.mutable_copy(seed)
            @test sum(a; init = seed) == expected + seed
            @test seed == original_seed
        end
        for op in (MP.MA.add_mul, MP.MA.sub_mul)
            expected = op(zero(T(1) * (x + y + z)), terms[1], polys[1])
            expected = op(expected, terms[2], polys[2])
            @test (@inferred MP.MA.fused_map_reduce(op, terms, polys)) ==
                  expected
        end
        @test (@inferred LinearAlgebra.dot(terms, polys)) ==
              terms[1] * polys[1] + terms[2] * polys[2]
        init = MP.MA.zero_for(P, terms, polys)
        @test MP.MA.fused_map_reduce(MP.MA.add_mul, terms, polys; init) === init
        @test_throws DimensionMismatch MP.MA.fused_map_reduce(
            MP.MA.add_mul,
            terms,
            polys[1:1],
        )
        @test (terms, polys) == originals
        @test iszero(@inferred sum(P[]))
        @test iszero(@inferred MP.MA.fused_map_reduce(MP.MA.add_mul, P[], P[]))
        # Even a zero operand contributes its parent to output preparation.
        result = sum([polys[1], zero(T(1) * x + T(0))])
        @test result == polys[1]
        @test parent(result) == parent(x + y)
    end
    @test sum([x, y]) == x + y
    a, b = [x, y], [y + 1, z + 2]
    @test (@inferred LinearAlgebra.dot(a, b)) == x * (y + 1) + y * (z + 2)
    @test (@inferred LinearAlgebra.dot(b, a)) == (y + 1) * x + (z + 2) * y
    @test LinearAlgebra.dot([0.5, 2.5], a) == 0.5x + 2.5y
    @test LinearAlgebra.dot(a, [0.5, 2.5]) == 0.5x + 2.5y
    DP.@complex_polyvar w
    terms = DP.Term{Complex{Int}}[(1+im)*w, 2x]
    polys = [y + 1, z + 2]
    expected = (1 - im) * conj(w) * (y + 1) + 2x * (z + 2)
    @test (@inferred LinearAlgebra.dot(terms, polys)) == expected
    @test (@inferred MP.MA.fused_map_reduce(MP.MA.add_dot, terms, polys)) ==
          expected
    @test terms' * polys == expected
end

end
