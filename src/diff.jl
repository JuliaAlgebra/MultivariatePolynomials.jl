# I do not use it but I import the function to add a method
export differentiate

differentiate(p::PolyVar, x)  = differentiate(Term(p), x)
differentiate(p::Monomial, x) = differentiate(Term(p), x)

function differentiate{C, T}(t::Term{C, T}, x::PolyVar{C})
    i = findfirst(vars(t), x)
    if i == 0 || t.x.z[i] == 0
        S = Base.promote_op(*, T, Int)
        zero(Term{C, S})
    else
        z = copy(t.x.z)
        z[i] -= 1
        Term(t.α * t.x.z[i], Monomial(vars(t), z))
    end
end

function differentiate{C, T}(p::Polynomial{C, T}, x::PolyVar{C})
    # grlex order preserved
    i = findfirst(vars(p), x)
    S = Base.promote_op(*, T, Int)
    if i == 0
        zero(Polynomial{C, S})
    else
        keep = find([z[i] > 0 for z in p.x.Z])
        Z = [copy(p.x.Z[i]) for i in keep]
        a = Vector{S}(length(keep))
        for j in 1:length(Z)
            a[j] = p.a[keep[j]] * Z[j][i]
            Z[j][i] -= 1
        end
        Polynomial(a, MonomialVector(vars(p), Z))
    end
end

differentiate(p::MatPolynomial, x) = differentiate(Polynomial(p), x)

differentiate(p::RationalPoly, x::PolyVar) = (differentiate(p.num, x) * p.den - p.num * differentiate(p.den, x)) / p.den^2

# even if I annotate with ::Array{_diff_promote_op(T, PolyVar{C}), N+1}, it cannot detect the type since it seems to be unable to determine the dimension N+1 :(
function differentiate{N, C, T<:PolyType{C}}(ps::AbstractArray{T, N}, xs::Union{AbstractVector{PolyVar{C}}, Tuple})
    qs = Array{_diff_promote_op(T, PolyVar{C}), N+1}(length(xs), size(ps)...)
    for (i, x) in enumerate(xs)
        for j in linearindices(ps)
            J = ind2sub(ps, j)
            qs[i, J...] = differentiate(ps[J...], x)
        end
    end
    qs
end
function differentiate{C}(p::PolyType{C}, xs::Union{AbstractVector{PolyVar{C}}, Tuple})
    [differentiate(p, x) for x in xs]
end

# In Julia v0.5, Base.promote_op returns Any for PolyVar, Monomial and MatPolynomial
# Even on Julia v0.6 and Polynomial, Base.promote_op returns Any...
_diff_promote_op(S, T) = Base.promote_op(differentiate, S, T)
_diff_promote_op{C}(::Union{Type{PolyVar{C}}, Type{Monomial{C}}}, ::Type{PolyVar{C}}) = Term{true, Int}
_diff_promote_op{C, T}(::Union{Type{Polynomial{C, T}}, Type{MatPolynomial{C, T}}}, ::Type{PolyVar{C}}) = Polynomial{true, Base.promote_op(*, T, Int)}

function differentiate{T<:PolyType}(p::Union{T, AbstractArray{T}}, x, deg::Int)
    if deg < 0
        throw(DomainError())
    elseif deg == 0
        # Need the conversion with promote_op to be type stable for PolyVar, Monomial and MatPolynomial
        return _diff_promote_op(typeof(p), typeof(x))(p)
    else
        return differentiate(differentiate(p, x), x, deg-1)
    end
end
