# I do not use it but I import the function to add a method
export differentiate

differentiate(p::MatPolynomial, v) = differentiate(polynomial(p), v)

differentiate(p::RationalPoly, v::AbstractVariable) = (differentiate(p.num, v) * p.den - p.num * differentiate(p.den, v)) / p.den^2

const ARPL = Union{APL, RationalPoly}

_vec_diff_promote_op{PT, VT}(::Type{PT}, ::AbstractVector{VT}) = _diff_promote_op(PT, VT)
_vec_diff_promote_op{PT, N, VT}(::Type{PT}, ::NTuple{N, VT})   = _diff_promote_op(PT, VT)
_vec_diff_promote_op{PT, }(::Type{PT}, xs::Tuple)         = reduce(_diff_promote_op, PT, xs)

# even if I annotate with ::Array{_diff_promote_op(T, PolyVar{C}), N+1}, it cannot detect the type since it seems to be unable to determine the dimension N+1 :(
function differentiate{N, PT<:ARPL}(ps::AbstractArray{PT, N}, xs)
    qs = Array{_vec_diff_promote_op(PT, xs), N+1}(length(xs), size(ps)...)
    for (i, x) in enumerate(xs)
        for j in linearindices(ps)
            J = ind2sub(ps, j)
            qs[i, J...] = differentiate(ps[J...], x)
        end
    end
    qs
end
differentiate(p::ARPL, xs) = [differentiate(p, x) for x in xs]

# In Julia v0.5, Base.promote_op returns Any for PolyVar, Monomial and MatPolynomial
# Even on Julia v0.6 and Polynomial, Base.promote_op returns Any...
_diff_promote_op{PT, VT}(::Type{PT}, ::Type{VT}) = Base.promote_op(differentiate, PT, VT)
_diff_promote_op{MT<:AbstractMonomialLike}(::Type{MT}, ::Type{<:AbstractVariable}) = termtype(MT, Int)
_diff_promote_op{T, TT<:AbstractTerm{T}}(::Type{TT}, ::Type{<:AbstractVariable}) = changecoefficienttype(TT, Base.promote_op(*, T, Int))
_diff_promote_op{T, PT<:APL{T}}(::Type{PT}, ::Type{<:AbstractVariable}) = polynomialtype(PT, Base.promote_op(*, T, Int))

function differentiate{T<:ARPL}(p::Union{T, AbstractArray{T}}, x, deg::Int)
    if deg < 0
        throw(DomainError())
    elseif deg == 0
        # Need the conversion with promote_op to be type stable for PolyVar, Monomial and MatPolynomial
        return _diff_promote_op(typeof(p), typeof(x))(p)
    else
        return differentiate(differentiate(p, x), x, deg-1)
    end
end
