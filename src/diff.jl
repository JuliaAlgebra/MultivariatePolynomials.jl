# I do not use it but I import the function to add a method
export differentiate

function differentiate{C, T}(t::Term{C, T}, x::PolyVar{C})
    i = findfirst(vars(t), x)
    if i == 0 || t.x.z[i] == 0
        zero(t)
    else
        z = copy(t.x.z)
        z[i] -= 1
        Term(t.α * t.x.z[i], Monomial(vars(t), z))
    end
end

function differentiate{C, T}(p::Polynomial{C, T}, x::PolyVar{C})
    # grlex order preserved
    i = findfirst(vars(p), x)
    if i == 0
        zero(p)
    else
        keep = find([z[i] > 0 for z in p.x.Z])
        Z = [copy(p.x.Z[i]) for i in keep]
        S = Base.promote_op(*, T, Int)
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

function differentiate{C}(p::PolyType{C}, xs::Vector{PolyVar{C}})
    [differentiate(p, x) for x in xs]
end

function differentiate(p::Polynomial, x, deg::Int)
    deg < 0 && throw(DomainError("Cannot compute a negative derivative of a polynomial"))
    for i in 1:deg
        p = differentiate(p, x)
    end
    p
end
