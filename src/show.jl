function show(io::IO, v::AbstractVariable)
    print(io, name(v))
end

isone(x::T) where T = x == one(T)
function show(io::IO, m::AbstractMonomial)
    if isconstant(m)
        print(io, "1")
    else
        for (var, exp) in zip(variables(m), exponents(m))
            if !iszero(exp)
                print(io, var)
                if !isone(exp)
                    print(io, "^", exp)
                end
            end
        end
    end
end

function Base.show(io::IO, t::AbstractTerm)
    if isconstant(t)
        print(io, coefficient(t))
    else
        if !isone(coefficient(t))
            print(io, coefficient(t))
        end
        if !iszero(t)
            print(io, monomial(t))
        end
    end
end

lessthanone(x::Real) = x < 0
lessthanone(x::Complex) = false

printcoeff(io::IO, coeff::Real, monomial) = print(io, abs(coeff) * monomial)
printcoeff(io::IO, coeff::Complex, monomial) = print(io, "(", coeff, ")", monomial)

function Base.show(io::IO, p::AbstractPolynomial{T}) where T
    ts = terms(p)
    if isempty(ts)
        print(io, zero(T))
    else
        print(io, first(ts))
        for t in Iterators.drop(ts, 1)
            if lessthanone(coefficient(t))
                print(io, " - ")
            else
                print(io, " + ")
            end
            printcoeff(io, coefficient(t), monomial(t))
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
