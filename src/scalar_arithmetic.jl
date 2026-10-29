"""
    MultivariatePolynomials.@scalar_arithmetic Scalar

Allow values of `Scalar` to act as coefficients of bare variables and monomials.
Multiplication returns a term with the scalar as its coefficient. Addition and
subtraction delegate to StarAlgebras' term arithmetic. This also defines scalar
dot products, scalar promotion and conversion, and the MutableArithmetics
promotion, fused accumulation, and dot product methods.
Equality (`==` and `isequal`) with monomials, terms, polynomials, and rational
polynomials uses the corresponding coefficient comparison.
Each operation requires the corresponding coefficient operations and their
MutableArithmetics promotion methods.

Use this unexported macro at module scope, typically in a package extension.
Register disjoint scalar families. Register the same family separately with
`StarAlgebras.@scalar_arithmetic` to support arithmetic with terms and
polynomials. The JuMP and MathOptInterface extensions of both packages do this
automatically.
"""
macro scalar_arithmetic(Scalar)
    S = esc(Scalar)
    code = quote
        function Base.promote_rule(
            ::Type{P},
            ::Type{T},
        ) where {P<:AbstractMonomialLike,T<:$S}
            C = MA.promote_operation(+, coefficient_type(P), T)
            return similar_type(P, C)
        end

        function Base.promote_rule(
            ::Type{RationalPoly{NT,DT}},
            ::Type{T},
        ) where {NT,DT,T<:$S}
            return RationalPoly{promote_type(NT, T),promote_type(DT, term_type(DT))}
        end

        function Base.convert(
            ::Type{P},
            c::$S,
        ) where {P<:Union{AbstractTerm,AbstractPolynomial}}
            return convert(P, constant_term(c, P))
        end

        function Base.convert(::Type{M}, c::$S) where {M<:AbstractMonomial}
            isone(c) || throw(InexactError(:convert, M, c))
            return one(M)
        end

        function Base.convert(::Type{RationalPoly{NT,DT}}, c::$S) where {NT,DT}
            return convert(NT, c) / one(DT)
        end

        Base.:*(a::$S, m::_BasisElementLike) = term(a, m)
        Base.:*(m::_BasisElementLike, a::$S) = term(a, m)

        function MA.promote_operation(
            ::typeof(*),
            ::Type{T},
            ::Type{M},
        ) where {T<:$S,M<:_BasisElementLike}
            return term_type(M, T)
        end

        function MA.promote_operation(
            ::typeof(*),
            ::Type{M},
            ::Type{T},
        ) where {M<:_BasisElementLike,T<:$S}
            return term_type(M, T)
        end

        LinearAlgebra.dot(a::$S, p::AbstractMonomialLike) = a' * p
        LinearAlgebra.dot(p::AbstractMonomialLike, a::$S) = p' * a

        function LinearAlgebra._dot_nonrecursive(
            a::AbstractArray{<:$S},
            b::AbstractArray{<:AbstractMonomialLike},
        )
            return _sum_monomial_products(a, b)
        end

        function LinearAlgebra._dot_nonrecursive(
            a::AbstractArray{<:AbstractMonomialLike},
            b::AbstractArray{<:$S},
        )
            return _sum_monomial_products(b, a)
        end
    end
    for comp in (:(==), :isequal)
        push!(
            code.args,
            quote
                function Base.$comp(
                    p::Union{AbstractPolynomial,AbstractTerm},
                    a::$S,
                )
                    return right_constant_eq(p, a; comp = $comp)
                end

                function Base.$comp(
                    a::$S,
                    p::Union{AbstractPolynomial,AbstractTerm},
                )
                    return left_constant_eq(a, p; comp = $comp)
                end

                Base.$comp(m::AbstractMonomialLike, a::$S) =
                    right_constant_eq(m, a; comp = $comp)
                Base.$comp(a::$S, m::AbstractMonomialLike) =
                    left_constant_eq(a, m; comp = $comp)
                Base.$comp(a::$S, q::RationalPoly) = $comp(a * q.den, q.num)
                Base.$comp(q::RationalPoly, a::$S) = $comp(a, q)
            end,
        )
    end
    for left in (true, false)
        L, R = left ? (S, :AbstractMonomialLike) : (:AbstractMonomialLike, S)
        A, B = left ? (:A, :(term_type(B))) : (:(term_type(A)), :B)
        # If coefficients widen, operate!! falls back to operate. Both paths
        # multiply the scalar by a term before adding or subtracting it.
        for operate in (:operate, :operate!)
            push!(
                code.args,
                quote
                    function MA.$operate(
                        op::MA.AddSubMul,
                        p::AbstractPolynomial,
                        a::$L,
                        b::$R,
                    )
                        return MA.$operate(
                            op,
                            p,
                            _sum_product_operand(a),
                            _sum_product_operand(b),
                        )
                    end
                end,
            )
        end
        push!(
            code.args,
            quote
                function MA.operate(
                    ::typeof(LinearAlgebra.dot),
                    a::AbstractArray{<:$L},
                    b::AbstractArray{<:$R},
                )
                    return LinearAlgebra._dot_nonrecursive(adjoint.(a), b)
                end

                function MA.promote_operation(
                    op::Union{typeof(+),typeof(-)},
                    ::Type{A},
                    ::Type{B},
                ) where {A<:$L,B<:$R}
                    return MA.promote_operation(op, $A, $B)
                end
            end,
        )
        for op in (:+, :-)
            push!(
                code.args,
                quote
                    function Base.$op(a::$L, b::$R)
                        return $op(
                            _sum_product_operand(a),
                            _sum_product_operand(b),
                        )
                    end
                end,
            )
        end
    end
    return code
end

@scalar_arithmetic Number
