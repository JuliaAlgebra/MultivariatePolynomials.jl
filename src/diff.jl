# I do not use it but I import the function to add a method
export differentiate

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

function differentiate{C}(p::Polynomial{C}, xs::Vector{PolyVar{C}})
    [differentiate(p, x) for x in xs]
end

differentiate(p::MatPolynomial, x) = differentiate(Polynomial(p), x)
