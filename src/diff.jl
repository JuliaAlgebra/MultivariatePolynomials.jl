# I do not use it but I import the function to add a method
export differentiate

# Fallback for everything else
_diff_promote_op(::Type{T}, ::Type{<:AbstractVariable}) where T = T
differentiate(α::T, v::AbstractVariable) where T = zero(T)

_diff_promote_op(::Type{<:AbstractVariable}, ::Type{<:AbstractVariable}) = Int
differentiate(v1::AbstractVariable, v2::AbstractVariable) = v1 == v2 ? 1 : 0

_diff_promote_op(::Type{TT}, ::Type{<:AbstractVariable}) where {T, TT<:AbstractTermLike{T}} = changecoefficienttype(TT, Base.promote_op(*, T, Int))
differentiate(t::AbstractTermLike, v::AbstractVariable) = coefficient(t) * differentiate(monomial(t), v)

_diff_promote_op(::Type{PT}, ::Type{<:AbstractVariable}) where {T, PT<:APL{T}} = polynomialtype(PT, Base.promote_op(*, T, Int))
# The polynomial function will take care of removing the zeros
differentiate(p::APL, v::AbstractVariable) = polynomial(differentiate.(terms(p), v), SortedState())

differentiate(p::RationalPoly, v::AbstractVariable) = (differentiate(p.num, v) * p.den - p.num * differentiate(p.den, v)) / p.den^2

const ARPL = Union{APL, RationalPoly}

_vec_diff_promote_op(::Type{PT}, ::AbstractVector{VT}) where {PT, VT} = _diff_promote_op(PT, VT)
_vec_diff_promote_op(::Type{PT}, ::NTuple{N, VT}) where {PT, N, VT}   = _diff_promote_op(PT, VT)
_vec_diff_promote_op(::Type{PT}, ::VT, xs...) where {PT, VT}          = _diff_promote_op(PT, VT)
_vec_diff_promote_op(::Type{PT}, xs::Tuple) where PT = _vec_diff_promote_op(PT, xs...)

# even if I annotate with ::Array{_diff_promote_op(T, PolyVar{C}), N+1}, it cannot detect the type since it seems to be unable to determine the dimension N+1 :(
function differentiate(ps::AbstractArray{PT, N}, xs) where {N, PT<:ARPL}
    qs = Array{_vec_diff_promote_op(PT, xs), N+1}(length(xs), size(ps)...)
    for (i, x) in enumerate(xs)
        for j in linearindices(ps)
            J = ind2sub(ps, j)
            qs[i, J...] = differentiate(ps[J...], x)
        end
    end
    qs
end
differentiate(ps::AbstractArray{<:ARPL}, v::AbstractVariable) = differentiate.(ps, v)

differentiate(p::ARPL, xs) = [differentiate(p, x) for x in xs]

# differentiate(p, [x, y]) with TypedPolynomials promote x to a Monomial
differentiate(p::ARPL, m::AbstractMonomial) = differentiate(p, variable(m))

# In Julia v0.5, Base.promote_op returns Any for PolyVar, Monomial and MatPolynomial
# Even on Julia v0.6 and Polynomial, Base.promote_op returns Any...
_diff_promote_op(::Type{PT}, ::Type{VT}) where {PT, VT} = Base.promote_op(differentiate, PT, VT)
_diff_promote_op(::Type{MT}, ::Type{<:AbstractVariable}) where {MT<:AbstractMonomialLike} = termtype(MT, Int)

function differentiate(p, x, deg::Int)
    if deg < 0
        throw(DomainError())
    elseif deg == 0
        # Need the conversion with promote_op to be type stable for PolyVar, Monomial and MatPolynomial
        return convert(_diff_promote_op(typeof(p), typeof(x)), p)
    else
        return differentiate(differentiate(p, x), x, deg-1)
    end
end
