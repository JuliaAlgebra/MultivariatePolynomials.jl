function show(io::IO, v::AbstractVariable)
    print(io, name(v))
end

function Base.show(io::IO, t::AbstractTerm)
    cst = isconstant(t)
    if coefficient(t) != 1 || cst
        print(io, coefficient(t))
    end
    if !cst
        print(io, monomial(t))
    end
end

function Base.show(io::IO, p::AbstractPolynomial)
    n = nterms(p)
    for (i, t) in enumerate(terms(p))
        print(io, t)
        if i != n
            print(io, " + ")
        end
    end
end

function Base.show(io::IO, p::SOSDecomposition)
    for (i, q) in enumerate(p)
        print(io, "(")
        print(io, q)
        print(io, ")^2")
        if i != length(p)
            print(io, " + ")
        end
    end
end

function Base.show(io::IO, p::RationalPoly)
    print(io, "(")
    print(io, p.num)
    print(io, ") / (")
    print(io, p.den)
    print(io, ")")
end
