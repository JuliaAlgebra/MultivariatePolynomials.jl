function Base.round(p::_APL; args...)
    # round(0.1) is zero so we cannot use `nonzero=true`
    return map_coefficients(p) do term
        return round(term; args...)
    end
end

function Base.div(p::_APL, α::Number, args...)
    return map_coefficients(p) do term
        return div(term, α, args...)
    end
end

"""
    divides(t1::AbstractTermLike, t2::AbstractTermLike)

Returns whether the monomial of t1 divides the monomial of t2.

### Examples

Calling `divides(2x^2y, 3xy)` should return false because `x^2y` does not divide `xy` since `x` has a degree 2 in `x^2y` which is greater than the degree of `x` on `xy`.
However, calling `divides(3xy, 2x^2y)` should return true.
"""
function divides(t1::AbstractTermLike, t2::AbstractTermLike)
    return divides(monomial(t1), monomial(t2))
end
divides(t1::AbstractVariable, t2::AbstractVariable) = t1 == t2

function Base.gcd(m1::AbstractMonomialLike, m2::AbstractMonomialLike)
    return map_exponents(min, m1, m2)
end
function Base.lcm(m1::AbstractMonomialLike, m2::AbstractMonomialLike)
    return map_exponents(max, m1, m2)
end

struct Field end
struct UniqueFactorizationDomain end
const UFD = UniqueFactorizationDomain

algebraic_structure(::Type{<:Integer}) = UFD()
algebraic_structure(::Type{<:AbstractPolynomialLike}) = UFD()
# `Rational`, `AbstractFloat`, JuMP expressions, etc... are fields
algebraic_structure(::Type) = Field()
_field_absorb(::UFD, ::UFD) = UFD()
_field_absorb(::UFD, ::Field) = Field()
_field_absorb(::Field, ::UFD) = Field()
_field_absorb(::Field, ::Field) = Field()

"""
    div_multiple(a, b, ma::MA.MutableTrait)

Return the division of `a` by `b` assuming that `a` is a multiple of `b`.
If `a` is not a multiple of `b` then this function may return anything.
"""
div_multiple(::Field, a, b, ma::MA.MutableTrait) = a / b
div_multiple(::UFD, a, b, ma::MA.IsMutable) = MA.operate!!(div, a, b)
div_multiple(::UFD, a, b, ma::MA.IsNotMutable) = div(a, b)
function div_multiple(a, b, ma::MA.MutableTrait = MA.IsNotMutable())
    return div_multiple(
        algebraic_structure(promote_type(typeof(a), typeof(b))),
        a,
        b,
        ma,
    )
end
function div_multiple(
    m1::AbstractMonomialLike,
    m2::AbstractMonomialLike,
    ::MA.MutableTrait = MA.IsNotMutable(),
)
    return map_exponents(-, m1, m2)
end
function div_multiple(
    t::AbstractTerm,
    m::AbstractMonomial,
    mt::MA.MutableTrait = MA.IsNotMutable(),
)
    return term(_copy(coefficient(t), mt), div_multiple(monomial(t), m))
end
function div_multiple(
    t1::AbstractTermLike,
    t2::AbstractTermLike,
    m1::MA.MutableTrait = MA.IsNotMutable(),
)
    return term(
        div_multiple(coefficient(t1), coefficient(t2), m1),
        div_multiple(monomial(t1), monomial(t2)),
    )
end
function right_constant_div_multiple(
    f::_APL,
    g,
    mf::MA.MutableTrait = MA.IsNotMutable(),
)
    if isone(g)
        return _copy(f, mf)
    end
    return map_coefficients(
        coef -> div_multiple(coef, g, mf),
        f,
        mf;
        nonzero = true,
    )
end
function div_multiple(
    f::_APL,
    g::AbstractMonomialLike,
    mf::MA.MutableTrait = MA.IsNotMutable(),
)
    if isconstant(g)
        return _copy(f, mf)
    end
    return map_exponents(-, f, g, mf)
end
function div_multiple(
    f::_APL,
    g::AbstractTermLike,
    mf::MA.MutableTrait = MA.IsNotMutable(),
)
    f = right_constant_div_multiple(f, coefficient(g), mf)
    return div_multiple(f, monomial(g), MA.IsMutable())
end
function div_multiple(f::_APL, g::_APL, mf::MA.MutableTrait = MA.IsNotMutable())
    lt = leading_term(g)
    if nterms(g) == 1
        return div_multiple(f, lt, mf)
    end
    rf = _copy(f, mf)
    rg = remove_leading_term(g)
    q = zero(rf)
    while !iszero(rf)
        ltf = leading_term(rf)
        if !divides(lt, ltf)
            # In floating point arithmetics, it may happen
            # that `rf` is not zero even if it cannot be reduced further.
            # As `div_multiple` assumes that `g` divides `f`, we know that
            # `rf` is approximately zero anyway.
            break
        end
        qt = div_multiple(ltf, lt)
        q = MA.add!!(q, qt)
        rf = MA.operate!!(remove_leading_term, rf)
        rf = MA.operate!!(MA.sub_mul, rf, qt, rg)
    end
    return q
end

function Base.div(f::_APL, g::Union{_APL,AbstractVector{<:_APL}}; kwargs...)
    return divrem(f, g; kwargs...)[1]
end
function Base.rem(f::_APL, g::Union{_APL,AbstractVector{<:_APL}}; kwargs...)
    return divrem(f, g; kwargs...)[2]
end

function pseudo_divrem(f::_APL{S}, g::_APL{T}, algo) where {S,T}
    return _pseudo_divrem(
        algebraic_structure(MA.promote_operation(-, S, T)),
        f,
        g,
        algo,
    )
end

function _pseudo_divrem(::Field, f::_APL, g::_APL, algo)
    q, r = divrem(f, g)
    return one(q), q, r
end

function _pseudo_divrem(::UFD, f::_APL, g::_APL, algo)
    ltg = leading_term(g)
    rg = remove_leading_term(g)
    ltf = leading_term(f)
    if iszero(f) || !divides(monomial(ltg), ltf)
        return one(f), zero(f), zero(f)
    else
        st = constant_term(coefficient(ltg), f)
        new_f = st * remove_leading_term(f)
        qt = term(coefficient(ltf), div_multiple(monomial(ltf), monomial(ltg)))
        new_g = qt * rg
        # Check with `::` that we don't have any type unstability on this variable.
        return convert(typeof(f), st),
        convert(typeof(f), qt),
        (new_f - new_g)::typeof(f)
    end
end

"""
    pseudo_rem(f::_APL, g::_APL, algo)

Return the pseudo remainder of `f` modulo `g` as defined in [Knu14, Algorithm R, p. 425].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function pseudo_rem(f::_APL, g::_APL, algo)
    return MA.operate!!(pseudo_rem, MA.mutable_copy(f), g, algo)
end

function MA.promote_operation(
    ::typeof(pseudo_rem),
    ::Type{P},
    ::Type{Q},
    ::Type{A},
) where {T,S,P<:_APL{T},Q<:_APL{S},A}
    U1 = MA.promote_operation(*, S, T)
    U2 = MA.promote_operation(*, T, S)
    # `promote_type(P, Q)` is needed for TypedPolynomials in case they use different variables
    return polynomial_type(promote_type(P, Q), MA.promote_operation(-, U1, U2))
end

function MA.buffer_for(::typeof(pseudo_rem), F::Type, G::Type, ::Type)
    return MA.buffer_for(MA.sub_mul, F, term_type(F), G)
end

function _prepare_s_poly!(::typeof(pseudo_rem), f, ltf, ltg)
    MA.operate!(right_constant_mult, f, coefficient(ltg))
    return term(coefficient(ltf), div_multiple(monomial(ltf), monomial(ltg)))
end

function _prepare_s_poly!(::typeof(rem), ::_APL, ltf, ltg)
    return div_multiple(ltf, ltg)
end

function MA.operate!(
    op::Union{typeof(rem),typeof(pseudo_rem)},
    f::_APL,
    g::_APL,
    algo,
)
    return MA.buffered_operate!(nothing, op, f, g, algo)::typeof(f)
end

# TODO As suggested in [Knu14, Algorithm R, p. 426] (univariate case only), if
#      `deg(f) = n` and `deg(g) = m`, during the first `n - m - t` steps, the
#      coefficient of the `t`th power will simply be multiplied by the leading
#      coefficient of `g` so we could speed up by multiplying this coefficient
#      to the power `n - m - t` using `Base.power_by_squaring`.
#      However, since there might be missing terms, so we don't know in advance
#      the needed power but we could keep track of it.
function MA.buffered_operate!(
    buffer,
    op::Union{typeof(rem),typeof(pseudo_rem)},
    f::_APL,
    g::_APL,
    algo,
)
    ltg = leading_term(g)
    ltf = leading_term(f)
    # This only makes sense in the univariate case but it's only used for univariate gcd anyway
    skipped_divisions = maxdegree(f) - maxdegree(g) + 1
    MA.operate!(remove_leading_term, g)
    while !iszero(f)
        if isapproxzero(ltf) # TODO `, kwargs...)`
            MA.operate!(remove_leading_term, f)
        elseif !divides(monomial(ltg), ltf)
            # Since the monomials are sorted in decreasing order,
            # lm is larger than all of them hence it cannot divide any of them
            # This is always the case for univariate.
            # TODO We could also do early termination for Lex order even if `>` returns `false` here
            if monomial(ltg) > monomial(ltf)
                break
            end
            MA.operate!(remove_leading_term, f)
        else
            MA.operate!(remove_leading_term, f)
            t = _prepare_s_poly!(op, f, ltf, ltg)
            skipped_divisions -= 1
            MA.buffered_operate!(buffer, MA.sub_mul, f, t, g)
        end
        if op === pseudo_rem && _primitive_rem(algo)
            f = primitive_part(f, algo, MA.IsMutable())::typeof(f)
        end
        if _skip_last(algo) && maxdegree(f) == maxdegree(g)
            break
        end
        ltf = leading_term(f)
    end
    # Add it back as we cannot modify `g`
    MA.operate!(unsafe_restore_leading_term, g, ltg)
    _set_skipped_divisions!(algo, skipped_divisions)
    return f
end

"""
    rem_or_pseudo_rem(f::_APL, g::_APL, algo)

If the coefficient type is a field, return `rem`, otherwise, return [`pseudo_rem`](ref).
"""
function rem_or_pseudo_rem(f::_APL, g::_APL, algo)
    return MA.operate!!(rem_or_pseudo_rem, MA.mutable_copy(f), g, algo)
end

_op(::Field) = rem
_op(::UFD) = pseudo_rem

function MA.operate!(
    ::typeof(rem_or_pseudo_rem),
    f::_APL{S},
    g::_APL{T},
    algo,
) where {S,T}
    return MA.operate!(
        _op(algebraic_structure(MA.promote_operation(-, S, T))),
        f,
        g,
        algo,
    )::typeof(f)
end

function MA.buffered_operate!(
    buffer,
    ::typeof(rem_or_pseudo_rem),
    f::_APL{S},
    g::_APL{T},
    algo,
) where {S,T}
    return MA.buffered_operate!(
        buffer,
        _op(algebraic_structure(MA.promote_operation(-, S, T))),
        f,
        g,
        algo,
    )::typeof(f)
end

function MA.buffer_for(
    ::typeof(rem_or_pseudo_rem),
    F::Type{<:_APL{S}},
    G::Type{<:_APL{T}},
    A::Type,
) where {S,T}
    return MA.buffer_for(
        _op(algebraic_structure(MA.promote_operation(-, S, T))),
        F,
        G,
        A,
    )
end

function MA.promote_operation(
    ::typeof(rem_or_pseudo_rem),
    ::Type{P},
    ::Type{Q},
    ::Type{A},
) where {T,S,P<:_APL{T},Q<:_APL{S},A}
    return _promote_operation_rem_or_pseudo_rem(
        algebraic_structure(MA.promote_operation(-, S, T)),
        P,
        Q,
        A,
    )
end
function _promote_operation_rem_or_pseudo_rem(
    ::Field,
    ::Type{P},
    ::Type{Q},
    ::Type{A},
) where {P<:_APL,Q<:_APL,A}
    return MA.promote_operation(rem, P, Q)
end
function _promote_operation_rem_or_pseudo_rem(
    ::UFD,
    ::Type{P},
    ::Type{Q},
    ::Type{A},
) where {P<:_APL,Q<:_APL,A}
    return MA.promote_operation(pseudo_rem, P, Q, A)
end

function MA.promote_operation(
    ::typeof(rem),
    ::Type{P},
    ::Type{Q},
    ::Type{A},
) where {P<:_APL,Q<:_APL,A}
    return MA.promote_operation(rem, P, Q)
end
function MA.promote_operation(
    ::Union{typeof(div),typeof(rem)},
    ::Type{P},
    ::Type{Q},
) where {T,S,P<:_APL{T},Q<:_APL{S}}
    U = MA.promote_operation(/, T, S)
    # `promote_type(P, Q)` is needed for TypedPolynomials in case they use different variables
    return polynomial_type(promote_type(P, Q), MA.promote_operation(-, U, U))
end
function Base.divrem(f::_APL{T}, g::_APL{S}; kwargs...) where {T,S}
    rf = convert(
        MA.promote_operation(div, typeof(f), typeof(g)),
        MA.copy_if_mutable(f),
    )
    q = zero(rf)
    r = zero(rf)
    lt = leading_term(g)
    rg = remove_leading_term(g)
    lm = monomial(lt)
    while !iszero(rf)
        ltf = leading_term(rf)
        if isapproxzero(ltf; kwargs...)
            rf = MA.operate!!(remove_leading_term, rf)
        elseif divides(lm, ltf)
            qt = div_multiple(ltf, lt)
            q = MA.add!!(q, qt)
            rf = MA.operate!!(remove_leading_term, rf)
            rf = MA.operate!!(MA.sub_mul, rf, qt, rg)
        elseif lm > monomial(ltf)
            # Since the monomials are sorted in decreasing order,
            # lm is larger than all of them hence it cannot divide any of them
            r = MA.add!!(r, rf)
            break
        else
            r = MA.add!!(r, ltf)
            rf = MA.operate!!(remove_leading_term, rf)
        end
    end
    return q, r
end
function Base.divrem(
    f::_APL{T},
    g::AbstractVector{<:_APL{S}};
    kwargs...,
) where {T,S}
    rf = convert(
        MA.promote_operation(div, typeof(f), eltype(g)),
        MA.copy_if_mutable(f),
    )
    r = zero(rf)
    q = similar(g, typeof(rf))
    for i in eachindex(q)
        q[i] = zero(rf)
    end
    lt = leading_term.(g)
    rg = remove_leading_term.(g)
    lm = monomial.(lt)
    useful = BitSet(eachindex(g))
    while !iszero(rf)
        ltf = leading_term(rf)
        if isapproxzero(ltf; kwargs...)
            rf = MA.operate!!(remove_leading_term, rf)
            continue
        end
        divisionoccured = false
        for i in useful
            if divides(lm[i], ltf)
                qt = div_multiple(ltf, lt[i])
                q[i] = MA.add!!(q[i], qt)
                rf = MA.operate!!(remove_leading_term, rf)
                rf = MA.operate!!(MA.sub_mul, rf, qt, rg[i])
                divisionoccured = true
                break
            elseif lm[i] > monomial(ltf)
                # Since the monomials are sorted in decreasing order,
                # lm is larger than all of them hence it cannot divide any of them
                delete!(useful, i)
            end
        end
        if !divisionoccured
            if isempty(useful)
                r = MA.add!!(r, rf)
                break
            else
                r = MA.add!!(r, ltf)
                rf = MA.operate!!(remove_leading_term, rf)
            end
        end
    end
    return q, r
end
