module TestSAPolynomial

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import StarAlgebras as SA
import LinearAlgebra

@testset "Shared polynomial representation" begin
    DP.@polyvar x y
    p = 2x + y
    @test x isa MP.AbstractPolynomialLike
    @test x^2 isa MP.AbstractTermLike
    @test 2x isa MP.AbstractTerm
    @test p isa MP.AbstractPolynomial
    @test typeof(p) === DP.Polynomial{Int}
    @test MP.term_type(x) === DP.Term{Int}
    @test MP.polynomial_type(x) === typeof(p)
    @test MP.polynomial_type(p, Float64) === DP.Polynomial{Float64}
    @test MP.coefficient_type(x^2) === Int
    @test MP.coefficients(p) == [1, 2]
    @test MP.coefficient.(MP.terms(p)) == [1, 2]
    @test length(MP.monomials(p)) == 2
    @test MP.polynomial(p) === p
    @test MP.polynomial!(p) === p
    @test MP.polynomial!(2x, Float64) isa DP.Polynomial{Float64}
    @test sprint(show, typeof(p)) == "DynamicPolynomials.Polynomial{$Int}"
    @test sprint(show, typeof(2x)) == "DynamicPolynomials.Term{$Int}"
    @test promote_type(typeof(x), String) === Any
    @test_throws MethodError p * "x"
    @test_throws InexactError p * 0.5
    @test MP.coefficients(p / 2) == [0.5, 1.0]
    @test MP.coefficients(p // 2) == [1 // 2, 1 // 1]
    @test div(p, 2) == MP.polynomial(x)
    @test_throws InexactError p / 2.5
    @test iszero(p - (x + y) - x)
    @test iszero((x + y)^2 - (x^2 + 2x*y + y^2))
    @test iszero(sum([x, y]) - (x + y))
    @test iszero(sum([p, p]) - 2p)
    @test MP.coefficients(p) == [1, 2]
    @test iszero(sum(DP.Polynomial{Int}[]))
    @test !isone(zero(p))
    @test SA.coeffs(zero(p)) != SA.coeffs(p)
    @test MP.coefficients(MP.differentiate(p, x)) == [2]

    @test MP.coefficient(p, x) == 2
    @test MP.coefficient(x, y) == 0
    @test MP.coefficients(x, [x, y]) == [1, 0]
    @test MP.MA.scaling(one(x)) == 1
    @test p - y == MP.polynomial(2x)
    @test MP.subs(p, x => 2, y => 3) == 7
    @test iszero(MP.subs(p, x => 2) - (y + 4))
    @test MP.monomial([x, y], [1, 0]) == MP.monomial(x)
    @test hash(MP.monomial([x, y], [1, 0])) == hash(MP.monomial(x))

    @test p(x => 2, y => 3) == 7
    @test (2x)(x => 3) == 6
    @test MP.coefficient_type([x, y]) === Int

    DP.@polyvar u v monomial_order=MP.LexOrder
    q = u + v
    @test !(q isa DP.Polynomial)
    @test MP.polynomial_type(u) === typeof(q)
    @test iszero(sum([u, v]) - q)
end

@testset "Ordered noncommutative monomials" begin
    DP.@ncpolyvar x y z
    xy, yx = @inferred(x * y), @inferred(y * x)
    @test xy != yx
    @test MP.variables(xy) == [x, y]
    @test MP.variables(yx) == [y, x]
    @test MP.exponents(xy) == MP.exponents(yx) == [1, 1]
    @test x * x == x^2
    for (word, vars, exps) in (
        (x * (y * x), [x, y, x], [1, 1, 1]),
        ((x * y^2) * (y^3 * z), [x, y, z], [1, 5, 1]),
        ((x * x * y) * (y * x), [x, y, x], [2, 2, 1]),
        ((x * y)^2, [x, y, x, y], [1, 1, 1, 1]),
        ((x * y)^3, [x, y, x, y, x, y], ones(Int, 6)),
        (x^1000000 * y, [x, y], [1000000, 1]),
    )
        @test MP.variables(word) == vars
        @test MP.exponents(word) == exps
    end
    @test x^0 * y == y
    @test y * x^0 == y
    @test isone((x * y)^0)
    @test (x * y)^1 == xy
    copied = xy^1
    @test MP.variables(copied) !== MP.variables(xy)
    @test MP.exponents(copied) !== MP.exponents(xy)
    MP.variables(copied)[1] = z
    MP.exponents(copied)[2] = 3
    @test MP.variables(xy) == [x, y]
    @test MP.exponents(xy) == [1, 1]
    @test (@inferred MP.degree(x * y^4 * x^2, x)) == 3

    a, b = @inferred SA.promote_bases(xy, yx)
    @test MP.variables(a) == MP.variables(b) == [x, y, x]
    @test MP.exponents(a) == [1, 1, 0]
    @test MP.exponents(b) == [0, 1, 1]
    @test a == xy && b == yx
    @test isequal(a, xy) && hash(a) == hash(xy)
    @test MP.exponents(xy) == MP.exponents(yx) == [1, 1]
    @test (@inferred MP.exponents(x^2 * y * x, [x, y, x, z])) == [2, 1, 1, 0]
    @test MP.exponents(a, MP.variables(a)) == MP.exponents(a)
    @test MP.exponents(a, MP.variables(a)) !== MP.exponents(a)
    @test_throws ArgumentError MP.exponents(yx, [x, y])

    p = @inferred x + y
    @test typeof(2x) === MP.term_type(x)
    @test typeof(p) === MP.polynomial_type(x)
    @test p == y + x
    @test (@inferred sum([x, y])) == p
    @test (@inferred LinearAlgebra.dot([x, y], [x, y])) == x^2 + y^2
    @test (x + 1) * (x + 2) == x^2 + 3x + 2
    p = MP.term(2, xy) + MP.term(3, yx)
    @test MP.coefficient(p, xy) == 2
    @test MP.coefficient(p, yx) == 3

    DP.@polyvar u v w
    @test (@inferred MP.exponents(u^2 * v, [u, v, w])) == [2, 1, 0]
    @test_throws ArgumentError MP.exponents(u * v, [u])
end

@testset "Parent-aware sums" begin
    DP.@polyvar x y z
    left, right = @inferred SA.promote_bases(2x, 3y)
    @test parent(left) == parent(right)
    @test left == 2x
    @test right == 3y
    for a in ([x, y], [2x, 3y], [x + 1, y + 2])
        originals = deepcopy(a)
        expected = a[1] + a[2]
        @test (@inferred sum(a)) == expected
        @test (@inferred MP.MA.operate(sum, a)) == expected
        @test sum(a; init = z + 1) == expected + z + 1
        @test sum(a; init = 0.5) ==
              convert(DP.Polynomial{Float64}, expected) + 0.5
        @test a == originals
        @test parent(sum(a)) == parent(expected)
    end

    p = x + y
    a = [p, p]
    original = deepcopy(p)
    @test sum(a; init = p) == 3p
    @test p == original
    @test iszero(sum(DP.Polynomial{Int}[]))
    @test iszero(sum(DP.Term{Int}[]))
    @test iszero(sum(typeof(x)[]))
    @test sum(DP.Polynomial{Int}[]; init = p) == p
    @test parent(sum(DP.Polynomial{Int}[]; init = p)) === parent(p)

    a = [x y; z x]
    @test sum(a; dims = 1) == [x + z y + x]
    @test sum(a; dims = 2) == reshape([x + y, z + x], 2, 1)
    @test sum(a; dims = (1, 2)) == fill(2x + y + z, 1, 1)

    @test MP.polynomial([2x, 3y, -2x]) == 3y
    @test MP.differentiate(x^2 + 2x*y + y, x) == 2x + 2y
    DP.@complex_polyvar w
    @test adjoint((1 + im) * w + x) == (1 - im) * conj(w) + x
end

@testset "Parent-aware dot products" begin
    DP.@polyvar x y z
    for (a, b) in (
        ([x, y], [y, z]),
        ([2x, 3y], [4y, 5z]),
        ([x + 1, y + 2], [y - 1, z + 3]),
        ([2x, 3y], [y - 1, z + 3]),
        ([x + 1, y + 2], [4y, 5z]),
        ([x + 1, y + 2], [big(2) * y, big(3) * z]),
        (
            [big(2) * x + 1, big(3) * y + 2],
            [big(4) * y + 1, big(5) * z + 2],
        ),
    )
        originals = MP.MA.mutable_copy.((a, b))
        expected =
            LinearAlgebra.dot(a[1], b[1]) + LinearAlgebra.dot(a[2], b[2])
        @test (@inferred LinearAlgebra.dot(a, b)) == expected
        @test (@inferred MP.MA.operate(LinearAlgebra.dot, a, b)) == expected
        @test a' * b == expected
        @test LinearAlgebra.dot(reshape(a, 1, 2), reshape(b, 2, 1)) == expected
        @test (a, b) == originals
        @test_throws DimensionMismatch LinearAlgebra.dot(a, b[1:1])
    end
    for T in (Int, BigInt, Complex{Int}), P in (DP.Polynomial{T}, DP.Term{T})
        @test (@inferred LinearAlgebra.dot(P[], P[])) isa DP.Polynomial{T}
        result = LinearAlgebra.dot(P[], P[])
        @test iszero(result)
        @test MP.coefficient_type(result) === T
    end
    @test iszero(LinearAlgebra.dot(typeof(x)[], typeof(x)[]))

    DP.@complex_polyvar w
    a = DP.Polynomial{Complex{Int}}[(1 + 2im) * w + x, w + y]
    b = DP.Polynomial{Complex{Int}}[z + w, (2 - im) * w + x]
    originals = MP.MA.mutable_copy.((a, b))
    @test (@inferred LinearAlgebra.dot(a, b)) ==
          ((1 - 2im) * conj(w) + x) * b[1] + (conj(w) + y) * b[2]
    @test a' * b == LinearAlgebra.dot(a, b)
    @test transpose(a) * b == a[1] * b[1] + a[2] * b[2]
    @test (a, b) == originals
    @test LinearAlgebra.dot([(1 + im) * w, 2x], [3y, 1z]) ==
          3 * (1 - im) * conj(w) * y + 2x * z
end

@testset "Mixed numeric and polynomial dot products" begin
    DP.@polyvar x y
    for polynomials in (
        [x, y],
        [x^2, y^2],
        [2x, 3y],
        [x + 1, y + 2],
        [big(2) * x, big(3) * y],
        [big(2) * x + 1, big(3) * y + 2],
    ), numbers in ([2, 3], [2.0, 3.0]), left in (false, true)

        a, b = left ? (numbers, polynomials) : (polynomials, numbers)
        originals = MP.MA.mutable_copy.((a, b))
        expected =
            LinearAlgebra.dot(a[1], b[1]) + LinearAlgebra.dot(a[2], b[2])
        @test (@inferred LinearAlgebra.dot(a, b)) == expected
        @test (@inferred MP.MA.operate(LinearAlgebra.dot, a, b)) == expected
        @test a' * b == expected
        @test transpose(a) * b == a[1] * b[1] + a[2] * b[2]
        @test LinearAlgebra.dot(reshape(a, 1, 2), reshape(b, 2, 1)) == expected
        @test (a, b) == originals
        @test_throws DimensionMismatch LinearAlgebra.dot(a, b[1:1])
        @test_throws DimensionMismatch LinearAlgebra.dot(a[1:0], b)
    end

    for monomials in ([x, y], [x^2, y^2])
        expected = 0.5 * monomials[1] + 1.5 * monomials[2]
        @test (@inferred LinearAlgebra.dot([0.5, 1.5], monomials)) == expected
        @test (@inferred LinearAlgebra.dot(monomials, [0.5, 1.5])) == expected
    end
    for polynomials in ([2x, 3y], [x + 1, y + 2])
        @test_throws InexactError LinearAlgebra.dot([0.5, 1.5], polynomials)
        @test_throws InexactError LinearAlgebra.dot(polynomials, [0.5, 1.5])
        @test MP.coefficient_type(LinearAlgebra.dot([2.0, 3.0], polynomials)) ===
              Int
    end
    numbers = BigInt[2, 3]
    result = @inferred LinearAlgebra.dot(numbers, [x, y])
    MP.MA.operate!(+, first(MP.coefficients(result)), big(1))
    @test numbers == [2, 3]

    for T in (Int, BigInt, Complex{Int}), P in (DP.Polynomial{T}, DP.Term{T})
        for (a, b) in ((Float64[], P[]), (P[], Float64[]))
            result = @inferred LinearAlgebra.dot(a, b)
            @test result isa DP.Polynomial{T}
            @test iszero(result)
        end
    end
    for T in (Float64, BigInt), P in (typeof(x), typeof(x^2))
        for (a, b) in ((T[], P[]), (P[], T[]))
            result = @inferred LinearAlgebra.dot(a, b)
            @test result isa DP.Polynomial{T}
            @test iszero(result)
        end
    end

    DP.@complex_polyvar w
    numbers = [1 + 2im, 3 - im]
    for polynomials in (
        [w, x],
        [w^2, x^2],
        [(2 + im) * w, (3 - im) * x],
        [(2 + im) * w + x, (3 - im) * x + 1],
    ), left in (false, true)

        a, b = left ? (numbers, polynomials) : (polynomials, numbers)
        originals = MP.MA.mutable_copy.((a, b))
        expected =
            LinearAlgebra.dot(a[1], b[1]) + LinearAlgebra.dot(a[2], b[2])
        @test (@inferred LinearAlgebra.dot(a, b)) == expected
        @test a' * b == expected
        @test transpose(a) * b == a[1] * b[1] + a[2] * b[2]
        @test (a, b) == originals
    end
end

@testset "Matrix products in a common basis" begin
    DP.@polyvar x y z
    for (A, b) in (
        ([1 2; 3 4], [x, y]),
        ([0.5 1.5; 2.5 3.5], [x^2, y^2]),
        ([1 2; 3 4], [2x, 3y]),
        ([2.0 3.0; 4.0 5.0], [x + 1, y + 2]),
        ([1 2; 3 4], [big(2) * x + 1, big(3) * y + 2]),
        ([x y; z x], [0.5, 1.5]),
        ([x y; z x], [y, z]),
        ([x^2 y^2; z^2 x*y], [y, z]),
        ([2x 3y; 4z 1x], [x, z]),
        ([x y; z x], [2y, 3z]),
        ([x + 1 y + 2; z + 3 x + y], [x, y]),
        ([x y; z x], [y + 1, z + 2]),
        ([x + 1 y + 2; z + 3 x + y], [2, 3]),
        ([x + 1 y + 2; z + 3 x + y], [y + 1, z + 2]),
    )
        for B in (b, hcat(b, reverse(b), b))
            originals = MP.MA.mutable_copy.((A, B))
            expected = if B isa AbstractVector
                [A[i, 1] * B[1] + A[i, 2] * B[2] for i in 1:2]
            else
                [A[i, 1] * B[1, j] + A[i, 2] * B[2, j] for i in 1:2, j in 1:3]
            end
            result = @inferred A * B
            @test result == expected
            @test typeof(result) === typeof(expected)
            @test (@inferred MP.MA.operate(*, A, B)) == expected
            @test all(parent(r) == parent(first(result)) for r in result)
            output = similar(result)
            @test (@inferred LinearAlgebra.mul!(output, A, B)) === output
            @test output == expected
            @test (A, B) == originals
            bad = B isa AbstractVector ? B[1:1] : B[1:1, :]
            @test_throws DimensionMismatch A * bad
        end
    end

    DP.@complex_polyvar w
    A = DP.Polynomial{Complex{Int}}[
        (1 + im) * w + x w + y
        w + z (2 - im) * w + x
    ]
    b = [1 + im, 2 - im]
    for matrix in (A, transpose(A), adjoint(A), view(A, :, :))
        expected = [matrix[i, 1] * b[1] + matrix[i, 2] * b[2] for i in 1:2]
        @test (@inferred matrix * b) == expected
        @test (@inferred MP.MA.operate(*, matrix, b)) == expected
    end
    N = [1 + im 2 - im; 2 + im 1 - im]
    for (left, right) in (
        (transpose(A), N),
        (adjoint(A), N),
        (view(A, :, :), N),
        (N, transpose(A)),
        (N, adjoint(A)),
        (N, view(A, :, :)),
        (transpose(A), transpose(A)),
        (adjoint(A), adjoint(A)),
        (transpose(A), adjoint(A)),
        (adjoint(A), transpose(A)),
    )
        expected = [
            left[i, 1] * right[1, j] + left[i, 2] * right[2, j] for
            i in 1:2, j in 1:2
        ]
        @test (@inferred left * right) == expected
        @test (@inferred MP.MA.operate(*, left, right)) == expected
        output = similar(expected)
        @test (@inferred LinearAlgebra.mul!(output, left, right)) === output
        @test output == expected
    end

    for P in (typeof(x), DP.Term{Int}, DP.Polynomial{Int})
        result = @inferred zeros(Float64, 2, 0) * P[]
        @test length(result) == 2
        @test all(iszero, result)
        @test eltype(result) === MP.polynomial_type(
            P,
            P <: MP.AbstractMonomialLike ? Float64 : Int,
        )
        result = @inferred zeros(Float64, 2, 0) * Matrix{P}(undef, 0, 3)
        @test size(result) == (2, 3)
        @test all(iszero, result)
        @test eltype(result) === MP.polynomial_type(
            P,
            P <: MP.AbstractMonomialLike ? Float64 : Int,
        )
        @test isempty(@inferred zeros(Int, 0, 2) * Matrix{P}(undef, 2, 0))
    end
    @test isempty(@inferred zeros(Int, 0, 2) * [x, y])
    @test_throws InexactError [0.5 1.5] * [2x, 3y]
    @test_throws InexactError [0.5 1.5] * [x + 1, y + 2]
    @test_throws InexactError [0.5 1.5] * [2x 3y; 3y 2x]
    @test_throws InexactError [x + 1 y + 2] * [0.5 1.5; 0.5 1.5]
    @test MP.polynomial([1 2; 3 4], [x, y]) == x^2 + 5x*y + 4y^2
    @test (@inferred MP.polynomial([0.5 1.5; 0.5 3.0], [x, y])) ==
          0.5x^2 + 2x*y + 3y^2

    A = [x + 1 y + 2; x + y z + 3]
    originals = MP.MA.mutable_copy(A)
    for B in (A, view(A, :, :))
        @test_throws ArgumentError LinearAlgebra.mul!(A, B, [1 2; 3 4])
        @test_throws ArgumentError LinearAlgebra.mul!(A, [1 2; 3 4], B)
        @test A == originals
    end
    output = [x + 1 y + 2]
    @test_throws DimensionMismatch LinearAlgebra.mul!(output, A, A)
    @test output == [x + 1 y + 2]
    @test_throws DimensionMismatch LinearAlgebra.mul!(output, [1 2], A[1:1, :])
    @test output == [x + 1 y + 2]

    result =
        @inferred zeros(Int, 2, 0) * Matrix{DP.Polynomial{BigInt}}(undef, 0, 2)
    MP.MA.operate!(+, result[1], one(result[1]))
    @test result == [1 0; 0 0]
end

@testset "Equality and adjoints" begin
    DP.@polyvar x y
    p = x + y
    for comp in (==, isequal)
        @test comp(zero(p), 0) && comp(0, zero(p))
        @test comp(one(p), 1) && comp(1, one(p))
        @test comp(MP.polynomial(x), x) && comp(x, MP.polynomial(x))
        @test comp(MP.term(2, one(x)), 2) && comp(2, MP.term(2, one(x)))
        @test comp(one(x), 1) && comp(1, one(x))
        @test !comp(p, 0) && !comp(0, p)
    end
    @test adjoint(p) == transpose(p) == p
    constant = MP.term(2, one(x))
    @test hash(constant) == hash(MP.polynomial(constant)) == hash(2)
    @test haskey(Dict{Any,Bool}(constant => true), 2)
    t = (1 + 2im) * x
    @test adjoint(t) == (1 - 2im) * x
    @test transpose(t) == t
    c = t + (3 - im) * y
    @test adjoint(c) == (1 - 2im) * x + (3 + im) * y
    @test transpose(c) == c
    @test LinearAlgebra.dot(p, p) == p * p
    @test LinearAlgebra.dot(t, t) == 5 * x^2
    @test LinearAlgebra.dot([p, p], [p, p]) == 2 * p * p
    @test LinearAlgebra.dot([t, t], [t, t]) == 10 * x^2
end

@testset "Fused polynomial products" begin
    DP.@polyvar x y
    for T in (Int, BigInt), op in (MP.MA.add_mul, MP.MA.sub_mul)
        a, b = SA.promote_bases(
            convert(DP.Polynomial{T}, x + 1),
            convert(DP.Polynomial{T}, y + 2),
        )
        for (left, right) in ((a, b), (last(MP.terms(a)), last(MP.terms(b))))
            f = convert(DP.Polynomial{T}, x + y)
            originals = MP.MA.copy_if_mutable.((f, left, right))
            expected = op(f, left, right)
            out = zero(f)
            @test (@inferred MP.MA.operate_to!!(out, op, f, left, right)) === out
            @test out == expected
            @test (f, left, right) == originals
            @test (@inferred MP.MA.operate!!(op, f, left, right)) === f
            @test f == expected
            @test (left, right) == originals[2:3]
        end
        f = MP.MA.mutable_copy(a)
        expected = op(f, f, b)
        @test (@inferred MP.MA.operate!!(op, f, f, b)) === f
        @test f == expected
    end
end

@testset "Fused scalar products" begin
    DP.@polyvar x y
    for T in (Int, BigInt), op in (MP.MA.add_mul, MP.MA.sub_mul)
        g = convert(DP.Polynomial{T}, x + 2y)
        for p in (g, last(MP.terms(g))), (a, b) in ((2.0, p), (p, 2.0))
            f = one(g)
            expected = op(f, a, b)
            original = MP.MA.mutable_copy(g)
            out = zero(f)
            @test (@inferred MP.MA.operate_to!!(out, op, f, a, b)) === out
            @test out == expected
            @test (@inferred MP.MA.operate!!(op, f, a, b)) === f
            @test f == expected
            @test g == original
            # Mutating an inserted exponent vector must not modify the source.
            fill!(last(keys(SA.coeffs(f))), 0)
            @test g == original
        end
    end
end

@testset "Matching coefficient type" begin
    DP.@polyvar x
    a = [1 2; 3 4]
    t = MP.term(a, x)
    c = one(x)
    for p in (t + a, a + t)
        @test MP.coefficient(p, x) == a
        @test MP.coefficient(p, c) == a
    end
    @test MP.coefficient(t - a, c) == -a
    @test MP.coefficient(a - t, x) == -a
    @test MP.term(a, c) == a
    @test MP.polynomial(MP.term(a, c)) == a
end

@testset "Coefficient types in division and GCD" begin
    DP.@polyvar x y
    algo = MP.SubresultantAlgorithm()
    immutable = MP.MA.IsNotMutable()
    for (p, T, field) in (
        (x, Int, Rational{Int}),
        (x^2, Int, Rational{Int}),
        (2x, Int, Rational{Int}),
        (2x + 1, Int, Rational{Int}),
        ((2 // 3) * x, Rational{Int}, Rational{Int}),
        ((2 // 3) * x + 1, Rational{Int}, Rational{Int}),
        (2.0x, Float64, Float64),
        (2.0x + 1, Float64, Float64),
    )
        P = typeof(p)
        for op in (div, rem)
            @test MP.MA.promote_operation(op, P, P) === DP.Polynomial{field}
        end
        for op in (MP.pseudo_rem, MP.rem_or_pseudo_rem)
            @test MP.MA.promote_operation(op, P, P, typeof(algo)) ===
                  DP.Polynomial{T}
        end
    end
    for p in (x, x^2)
        @test (@inferred MP.content(p, algo, immutable)) === 1
        @test MP.primitive_part(p, algo, immutable) == p
    end
    @test MP.content(6x, algo, immutable) == 6
    @test MP.content(6x + 9, algo, immutable) == 3
    @test MP.content(zero(x + y), algo, immutable) == 0
    @test MP.primitive_part(6x + 9, algo, immutable) == 2x + 3
    for p in (2.5x, 2.5x + 1)
        @test MP.content(p, algo, immutable) === 1.0
        @test MP.primitive_part(p, algo, immutable) === p
    end
    complex_term = (2.0 + 3.0im) * x
    @test MP.primitive_part(complex_term, algo, immutable) === complex_term
    @test gcd(x^2, x^3) == x^2
    @test gcd(6x, 9x) == 3x
    @test isempty(Test.detect_unbound_args(MP))
end

@testset "Monomial divisibility" begin
    DP.@polyvar x y z
    for (a, b, expected) in (
        (x, x, true),
        (x, y, false),
        (x, x^2, true),
        (x^2, x, false),
        (x^2, x^2, true),
        (x^2, y^2, false),
        (x * y, x^2 * y^3, true),
        (x * y^3, x^2 * y, false),
        (y, x * y * z, true),
        (x * z, y^2, false),
        (x * y, x^2, false),
        (y * z, y^2, false),
        (one(x), x^2, true),
        (x^2, one(x), false),
        (one(x), one(y), true),
        (MP.monomial([x, y, z], [0, 2, 0]), y^3, true),
        (y, MP.monomial([x, y, z], [0, 2, 0]), true),
        (MP.term(3, x * y), MP.term(2, x^2 * y), true),
        (MP.term(2, x^2 * y), MP.term(3, x * y), false),
    )
        originals = deepcopy((a, b))
        @test MP.divides(a, b) === expected
        @test (a, b) == originals
    end
    @test (@inferred MP.divides(x, x^2))
    DP.@ncpolyvar u v
    @test MP.divides(u, u)
    @test !MP.divides(u, v)
    @test_throws ErrorException MP.divides(u, u^2)
    @test_throws ErrorException MP.divides(u^2, v^2)
end

@testset "Polynomial conversion" begin
    DP.@polyvar x
    R = DP.Polynomial{Rational{Int}}
    for p in (x, x^2, 2x, 2x + 1, 3, 0)
        converted = convert(R, p)
        @test converted isa R
        @test converted == p
        @test convert(R, converted) === converted
    end
    p = 2x + 1
    @test parent(convert(R, p)) === parent(p)
    @test_throws InexactError convert(DP.Polynomial{Int}, (1 // 2) * x + 1)
    @test_throws InexactError convert(DP.Polynomial{Int}, (1 // 2) * x)
    @test_throws InexactError convert(DP.Polynomial{Int}, 1 // 2)
end

@testset "Polynomial conversion in division and GCD" begin
    DP.@polyvar x y
    for p in (x, x^2, 2x, 2x + 1)
        original = deepcopy(p)
        q, r = divrem(p, typeof(p)[])
        @test isempty(q)
        @test r == p
        @test r isa DP.Polynomial{Rational{Int}}
        MP.MA.operate!(zero, r)
        @test p == original
    end
    q, r = divrem(zero(2x), 2x)
    @test q isa DP.Polynomial{Rational{Int}}
    @test iszero(q) && iszero(r)

    algo = MP.GeneralizedEuclideanAlgorithm()
    g = MP.primitive_univariate_gcd!(2x, zero(2x), algo)
    @test g isa DP.Polynomial{Int}
    @test g == 2x
    @test MP.inflate(3, one(x), one(x)) == 3
    nested = MP.term(2x, MP.monomial(y))
    c = MP.content(nested, algo, MP.MA.IsNotMutable())
    @test c isa DP.Term{Int}
    @test c == 2x
end

@testset "Leading-term removal and restoration" begin
    DP.@polyvar x
    p = big(2) * x + 3
    owned = MP.MA.copy_if_mutable(p)
    @test parent(owned) === parent(p)
    @test owned !== p
    @test last(keys(SA.coeffs(owned))) !== last(keys(SA.coeffs(p)))
    MP.MA.operate!(+, MP.leading_coefficient(owned), 1)
    last(keys(SA.coeffs(owned)))[1] = 2
    @test p == big(2) * x + 3
    @test MP.MA.operate!!(SA.remove_leading_term, owned) === owned
    @test owned == 3
    for p in (x, x^2, 2x)
        remainder = MP.MA.operate(SA.remove_leading_term, p)
        @test typeof(remainder) === typeof(zero(p))
        @test iszero(remainder)
    end
    for (p, expected) in (
        (zero(x + 1), 0),
        (MP.polynomial(2x), 0),
        (x^2 + x + 1, x + 1),
        (big(2) * x + 3, 3),
    )
        original = deepcopy(p)
        leading = MP.leading_term(p)
        remainder = SA.remove_leading_term(p)
        @test remainder == expected
        @test typeof(remainder) ===
              MP.MA.promote_operation(SA.remove_leading_term, typeof(p)) ===
              typeof(p)
        @test parent(remainder) === parent(p)
        @test p == original
        @test MP.MA.operate!(SA.remove_leading_term, p) === p
        @test p == expected
        @test all(c -> !iszero(c), values(SA.coeffs(p)))
        @test MP.MA.operate!(MP.unsafe_restore_leading_term, p, leading) === p
        @test p == original
        if !iszero(leading)
            @test MP.leading_coefficient(p) === MP.coefficient(leading)
        end
    end

    basis = SA.SubBasis(MP.FullBasis{MP.Monomial}([x]), [[0], [1], [2]])
    for c in ([1, 0, 3], SA.SparseArrays.sparsevec([1, 3], [1, 3], 3))
        p = SA.AlgebraElement(c, MP.algebra(basis))
        leading = MP.leading_term(p)
        @test SA.remove_leading_term(p) == 1
        @test MP.MA.operate!(SA.remove_leading_term, p) === p
        @test p == 1
        @test MP.MA.operate!(MP.unsafe_restore_leading_term, p, leading) === p
        @test p == 3x^2 + 1
    end

    coefficient = [1 2; 3 4]
    p = MP.polynomial(MP.term(coefficient, x))
    leading = MP.leading_term(p)
    MP.MA.operate!(SA.remove_leading_term, p)
    @test iszero(p)
    @test MP.coefficient(leading) == [1 2; 3 4]
    MP.MA.operate!(MP.unsafe_restore_leading_term, p, leading)
    @test MP.leading_coefficient(p) === coefficient
end

@testset "Shared coefficient mapping" begin
    DP.@polyvar x
    for T in (Int, BigInt)
        p = convert(DP.Polynomial{T}, 2x^2 + 3x + 2)
        original = deepcopy(p)
        mapped = @inferred MP.map_coefficients(c -> c // 2, p)
        @test mapped == x^2 + (3 // 2) * x + 1
        @test MP.coefficient_type(mapped) === Rational{T}
        @test parent(mapped) === parent(p)
        @test SA.coeffs(mapped) !== SA.coeffs(p)
        @test MP.map_coefficients(c -> c - 2, p) == x
        @test MP.map_coefficients(c -> 2c, p; nonzero = true) == 2p
        @test p == original

        z = zero(p)
        mapped = @inferred MP.map_coefficients(c -> c // 2, z)
        @test iszero(mapped)
        @test MP.coefficient_type(mapped) === Rational{T}
        @test parent(mapped) === parent(z)

        output = zero(p)
        @test SA.map_coefficients_to!(output, c -> c - 2, p) === output
        @test output == x
        @test p == original
        @test MP.MA.operate!(MP.right_constant_mult, p, T(2)) === p
        @test p == 2original
        @test MP.right_constant_div_multiple(p, T(2), MP.MA.IsMutable()) === p
        @test p == original
        @test SA.map_coefficients!(zero, p) === p
        @test iszero(p)

        f = convert(DP.Polynomial{T}, x^3 + 2x + 1)
        g = convert(DP.Polynomial{T}, 2x^2 + 1)
        g_original = deepcopy(g)
        @test MP.MA.operate!(
            MP.pseudo_rem,
            f,
            g,
            MP.GeneralizedEuclideanAlgorithm(),
        ) === f
        @test f == 3x + 2
        @test g == g_original
    end

    basis = SA.SubBasis(MP.FullBasis{MP.Monomial}([x]), [[0], [1], [2]])
    for c in ([2, 0, 4], SA.SparseArrays.sparsevec([1, 3], [2, 4], 3))
        p = SA.AlgebraElement(c, MP.algebra(basis))
        mapped = @inferred MP.map_coefficients(c -> c / 2 + 1, p)
        @test parent(mapped) === parent(p)
        @test SA.coeffs(mapped) == [2.0, 0.0, 3.0]
        @test SA.coeffs(p) == [2, 0, 4]
    end

    p = x + 1
    c = SA.SparseCoefficients(([0], [1]), (2, 3), SA.coeffs(p).isless)
    p = SA.AlgebraElement(c, parent(p))
    mapped = @inferred MP.map_coefficients(c -> c / 2, p)
    @test mapped == 1 + 1.5x
    @test mapped isa DP.Polynomial{Float64}
    @test parent(mapped) === parent(p)
    @test values(c) == (2, 3)

    a = [1 2; 3 4]
    p = MP.polynomial(MP.term(a, x))
    q = @inferred transpose(p)
    @test MP.coefficient(q, x) == transpose(a)
    @test parent(q) === parent(p)
    @test MP.coefficient_type(q) === typeof(transpose(a))
    @test MP.coefficient(p, x) == a
end

@testset "Coefficient mapping deprecations" begin
    DP.@polyvar x
    p = x + 1
    @test_deprecated r"SA\.map_coefficients!" MP.map_coefficients!(
        identity,
        p;
        nonzero = true,
    )
    @test_deprecated r"SA\.map_coefficients_to!" MP.map_coefficients_to!(
        p,
        identity,
        p;
        nonzero = true,
    )
end

@testset "Remainder operation routing" begin
    DP.@polyvar x y
    algo = MP.GeneralizedEuclideanAlgorithm()
    for T in (Int, BigInt, Rational{Int}, Float64)
        f = convert(DP.Polynomial{T}, x^3 + 2x + 1)
        g = convert(DP.Polynomial{T}, 2x^2 + 1)
        originals = deepcopy((f, g))
        for op in (MP.pseudo_rem, MP.rem_or_pseudo_rem)
            expected =
                op === MP.pseudo_rem || T <: Integer ? 3x + 2 : (3 // 2) * x + 1
            r = @inferred op(f, g, algo)
            @test r == expected
            @test typeof(r) === MP.MA.promote_operation(
                op,
                typeof(f),
                typeof(g),
                typeof(algo),
            )
            @test (f, g) == originals

            owned = MP.MA.mutable_copy(f)
            @test (@inferred MP.MA.operate!!(op, owned, g, algo)) === owned
            @test owned == expected
            @test g == originals[2]
        end
    end

    for op in (MP.pseudo_rem, MP.rem_or_pseudo_rem)
        for (f, g, expected) in (
            (x, 2x + 1, -1),
            (x^2, 2x, 0),
            (2x, x^2 + 1, 2x),
            (x^3 + 2x + 1, convert(DP.Polynomial{BigInt}, 2x^2 + 1), 3x + 2),
        )
            originals = deepcopy((f, g))
            r = @inferred MP.MA.operate!!(op, f, g, algo)
            @test r == expected
            @test (f, g) == originals
        end

        f = x^3 + 2x + 1
        for g in (f, SA.AlgebraElement(SA.coeffs(f), parent(f)))
            original = deepcopy(f)
            @test iszero(MP.MA.operate!!(op, f, g, algo))
            @test f == original
        end

        f = x^3 + 2x + 1
        g = 2x^2 + 1
        buffer = MP.MA.buffer_for(op, typeof(f), typeof(g), typeof(algo))
        @test MP.MA.buffered_operate!!(buffer, op, f, g, algo) === f
        @test f == 3x + 2
        @test g == 2x^2 + 1

        f = x^3 + 2x + 1
        originals = deepcopy((f, g))
        # A missing operate_to! method must not silently allocate.
        @test_throws ErrorException MP.MA.operate_to!!(zero(f), op, f, g, algo)
        @test (f, g) == originals

        g = 2x^2 + 1 + zero(x + y)
        originals = deepcopy((f, g))
        @test_throws ArgumentError MP.MA.operate!!(op, f, g, algo)
        @test (f, g) == originals
        @test op(f, g, algo) == 3x + 2
        @test (f, g) == originals
    end

    f = x^3 + 2x + 1
    g = 2x^2 + 1
    originals = deepcopy((f, g))
    r = @inferred MP.MA.operate!!(rem, f, g, algo)
    @test r == (3 // 2) * x + 1
    @test (f, g) == originals
    @test MP.MA.operate!!(rem, r, g, algo) === r
end

@testset "In-place quotient term accumulation" begin
    DP.@polyvar x y
    for T in (Int, BigInt)
        q = convert(DP.Polynomial{T}, x^2 + 1)
        t = MP.term(T(2), x)
        original = deepcopy(t)
        @test (@inferred MP.MA.add!!(q, t)) === q
        @test q == x^2 + 2x + 1
        @test MP.MA.add!!(q, -t) === q
        @test q == x^2 + 1
        @test t == original

        q = zero(q)
        @test MP.MA.add!!(q, t) === q
        stored = only(keys(SA.coeffs(q)))
        @test stored !== t.index
        @test stored == t.index
        stored[1] += 1
        @test t == original
        if T === BigInt
            MP.MA.operate!(+, MP.leading_coefficient(q), 1)
            @test t == original
        end
    end

    q = x + 1
    t = 0.5x
    original = deepcopy((q, t))
    result = @inferred MP.MA.add!!(q, t)
    @test result == 1.5x + 1
    @test result isa DP.Polynomial{Float64}
    @test (q, t) == original

    q = x + 1
    t = 2y
    original = deepcopy((q, t))
    @test_throws ArgumentError MP.MA.add!!(q, t)
    @test (q, t) == original
    @test_throws ArgumentError MP.MA.operate!(+, q, t)
    @test (q, t) == original
end

@testset "In-place polynomial addition" begin
    DP.@polyvar x y
    for T in (Int, BigInt)
        f = convert(DP.Polynomial{T}, x^2 + 1)
        g = convert(DP.Polynomial{T}, 2x - 1)
        original = deepcopy(g)
        @test (@inferred MP.MA.add!!(f, g)) === f
        @test f == x^2 + 2x
        @test g == original
        alias = SA.AlgebraElement(SA.coeffs(f), parent(f))
        @test (@inferred MP.MA.add!!(f, alias)) === f
        @test f == 2x^2 + 4x
        @test g == original
    end

    f, g = x + 1, 0.5x - 1
    originals = deepcopy((f, g))
    result = @inferred MP.MA.add!!(f, g)
    @test result == 1.5x
    @test result isa DP.Polynomial{Float64}
    @test (f, g) == originals

    g = y + 1
    originals = deepcopy((f, g))
    @test_throws ArgumentError MP.MA.add!!(f, g)
    @test (f, g) == originals
end

@testset "Polynomial division" begin
    DP.@polyvar x y
    for T in (Int, Rational{Int}, Float64, BigInt)
        for (f, g, q, r) in
            ((x^2 - 1, x - 1, x + 1, 0), (x^3 + 2x + 1, x^2 + 1, x, x + 1))
            f = convert(DP.Polynomial{T}, f)
            g = convert(DP.Polynomial{T}, g)
            originals = deepcopy((f, g))
            quotient, remainder = divrem(f, g)
            @test quotient == q
            @test remainder == r
            @test quotient * g + remainder == f
            @test (f, g) == originals
        end
    end
    f = x^2 * y + x + 1
    divisors = [MP.polynomial(x * y)]
    originals = deepcopy((f, divisors))
    quotients, remainder = divrem(f, divisors)
    @test quotients == [x]
    @test remainder == x + 1
    @test (f, divisors) == originals

    g = convert(DP.Polynomial{Rational{Int}}, 3x^2 + 1)
    f = 2 * one(g)
    original = deepcopy(g)
    @test MP.MA.operate!(rem, f, g, MP.GeneralizedEuclideanAlgorithm()) === f
    @test f == 2
    @test g == original
end

@testset "Ordered term products in division" begin
    DP.@polyvar x y z
    f, g, t = x^2 + x + 1, x + 1, 2x
    @test SA.term_product_style(SA.mstructure(f), SA.coeffs(f).isless) isa
          SA.OrderedTermProduct
    buffer = MP.MA.buffer_for(MP.MA.sub_mul, typeof(f), typeof(t), typeof(g))
    @test buffer === nothing
    @test MP.MA.buffered_operate!!(buffer, MP.MA.sub_mul, f, t, g) === f
    @test f == 1 - x - x^2
    @test g == x + 1

    # Coefficient promotion still uses allocating arithmetic.
    f, g, t = x + 1, x + 1, (1 // 2) * x
    original = deepcopy(f)
    result = MP.MA.operate!!(MP.MA.sub_mul, f, t, g)
    @test result == f - t * g
    @test result isa DP.Polynomial{Rational{Int}}
    @test f == original
    @test_throws ArgumentError MP.MA.operate!!(MP.MA.sub_mul, f, 2y, y + 1)
    @test f == original
    @test_throws ArgumentError MP.MA.operate_to!!(
        f,
        MP.MA.sub_mul,
        y + 1,
        2y,
        y + 1,
    )
    @test f == original

    f, g = x^2 + y, x + 1
    originals = deepcopy((f, g))
    @test divrem(f, g) == (x - 1, y + 1)
    @test (f, g) == originals
    f, divisors = x^2 + y^2 + z, [x + 1, y + 1]
    originals = deepcopy((f, divisors))
    q, r = divrem(f, divisors)
    @test q == [x - 1, y - 1]
    @test r == z + 2
    @test (f, divisors) == originals
    @test MP.div_multiple((x + y) * (x + 1), x + 1) == x + y

    f = convert(DP.Polynomial{Rational{Int}}, x^3 + 2x + 1)
    g = convert(DP.Polynomial{Rational{Int}}, x^2 + 1)
    original = deepcopy(g)
    @test MP.MA.operate!(rem, f, g, MP.GeneralizedEuclideanAlgorithm()) === f
    @test f == x + 1
    @test g == original
end

@testset "GCD result type promotion" begin
    DP.@polyvar x y
    for (a, b) in (
        (x, y),
        (x^2, x),
        (6x^2, 9x),
        (2x, y),
        (2.0x, 3y),
        ((2 // 3) * x, 3y),
        ((2 + 3im) * x, 3y),
        ((2 + 3im) * x, 3.0y),
        (MP.term(2x, MP.monomial(y)), MP.term(3x, MP.monomial(y))),
    )
        for (p, q) in ((a, b), (b, a))
            @test MP.MA.promote_operation(gcd, typeof(p), typeof(q)) ===
                  typeof(gcd(p, q))
            for algo in
                (MP.GeneralizedEuclideanAlgorithm(), MP.SubresultantAlgorithm())
                @test MP.MA.promote_operation(
                    gcd,
                    typeof(p),
                    typeof(q),
                    typeof(algo),
                ) === typeof(gcd(p, q, algo))
            end
        end
    end
    @test MP.MA.promote_operation(gcd, typeof(x + 1), typeof(x + 1)) ===
          DP.Polynomial{Int}
    nested = MP.term(x, MP.monomial(y))
    c = MP.content(nested, MP.SubresultantAlgorithm(), MP.MA.IsNotMutable())
    @test c isa MP.AbstractMonomial
    @test c == x
end

end # module
