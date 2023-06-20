const NonMutable = Union{Int,Float64}
const MP = MultivariatePolynomials
are_independent(::NonMutable, ::NonMutable) = true
are_independent(::BigInt, ::NonMutable) = true
are_independent(::NonMutable, ::BigInt) = true
function are_independent(a::BigInt, b::BigInt)
    return a !== b
end

function are_independent(a::Integer, b::Rational)
    return are_independent(a, b.num) && are_independent(a, b.den)
end

function are_independent(a::Rational, b::Real)
    return are_independent(a.num, b) && are_independent(a.den, b)
end

function are_independent(a::Number, p::_APL)
    return all(are_independent(a, coef) for coef in MP.coefficients(p))
end

function are_independent(p::_APL, a::Number)
    return are_independent(a, p)
end

function are_independent(p::_APL, q::_APL)
    for t in terms(p)
        for s in terms(q)
            if !are_independent(coefficient(t), coefficient(s))
                return false
            end
        end
    end
    return true
end
