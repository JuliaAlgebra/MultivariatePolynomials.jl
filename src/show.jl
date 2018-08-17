# VARIABLES
function Base.show(io::IO, v::AbstractVariable)
    print_var(io, MIME"text/plain"(), v)
end
function Base.show(io::IO, mime::MIME"text/latex", v::AbstractVariable)
    print_var(io, mime, v)
end

function print_var(io::IO, mime::MIME, var::AbstractVariable)
    print_var(io, mime, name(var))
end
function print_var(io::IO, mime::MIME, var::String)
    m = match(r"([a-zA-Z]+)(?:\[((?:\d,)*\d)\])", var)
    if m === nothing
        print(io, var)
    else
        print(io, m.captures[1])
        print_subscript(io, mime, parse.(Int, split(m.captures[2], ",")))
    end
end
function print_subscript(io::IO, ::MIME"text/latex", index)
    print(io, "_{", join(unicode_subscript.(index), ","), "}")
end
function print_subscript(io::IO, mime, index)
    if length(index) == 1
        print(io, unicode_subscript(index[1]))
    else
        print(io, join(unicode_subscript.(index), "\u208B"))
    end
end

# MONOMIALS

function Base.show(io::IO, mime::MIME"text/latex", m::AbstractMonomial)
    print_monomial(io, mime, m)
end
function Base.show(io::IO, mime::MIME"text/plain", m::AbstractMonomial)
    print_monomial(io, mime, m)
end
function Base.show(io::IO, m::AbstractMonomial)
    print_monomial(io, MIME"text/plain"(), m)
end

function print_monomial(io::IO, mime, m::AbstractMonomial)
    if isconstant(m)
        print(io, '1')
    else
        for (var, exp) in zip(variables(m), exponents(m))
            if !iszero(exp)
                print_var(io, mime, var)
                if !isone(exp)
                    print_exponent(io, mime, exp)
                end
            end
        end
    end
end
#
print_exponent(io::IO, ::MIME"text/latex", exp) = print(io, "^{", exp, "}")
function print_exponent(io::IO, mime, exp)
    print(io, join(unicode_superscript.(reverse(digits(exp)))))
end

# TERM

function Base.show(io::IO, t::AbstractTerm)
    print_term(io, MIME"text/plain"(), t)
end
function Base.show(io::IO, mime::MIME"text/latex", t::AbstractTerm)
    print_term(io, mime, t)
end
function Base.show(io::IO, mime::MIME"text/plain", t::AbstractTerm)
    print_term(io, mime, t)
end

function print_term(io::IO, mime, t::AbstractTerm)
    if isconstant(t)
        print_coefficient(io, coefficient(t))
    else
        if should_print_coefficient(coefficient(t))
            if !should_print_coefficient(-coefficient(t))
                print(io, '-')
            else
                print_coefficient(io, coefficient(t))
            end
        end
        if !iszero(t)
            print(io, monomial(t))
        end
    end
end

should_print_coefficient(x) = true  # By default, just print all coefficients
should_print_coefficient(x::Number) = !isone(x) # For numbers, we omit any "one" coefficients
print_coefficient(io::IO, coeff::Real) = print(io, coeff)
print_coefficient(io::IO, coeff) = print(io, "(", coeff, ")")

# POLYNOMIAL

function Base.show(io::IO, t::AbstractPolynomial)
    print_poly(io, MIME"text/plain"(), t)
end
function Base.show(io::IO, mime::MIME"text/plain", t::AbstractPolynomial)
    print_poly(io, mime, t)
end
function Base.show(io::IO, mime::MIME"text/latex", t::AbstractPolynomial)
    print_poly(io, mime, t)
end

function print_poly(io::IO, mime, p::AbstractPolynomial{T}) where T
    ts = terms(p)
    if isempty(ts)
        print(io, zero(T))
    else
        print_term(io, mime, first(ts))
        for t in Iterators.drop(ts, 1)
            if isnegative(coefficient(t))
                print(io, " - ")
                print_term(io, mime, abs(coefficient(t)) * monomial(t))
            else
                print(io, " + ")
                print_term(io, mime, t)
            end
        end
    end
end

isnegative(x::Real) = x < 0
isnegative(x) = false


function Base.show(io::IO, t::RationalPoly)
    print_ratpoly(io, MIME"text/plain"(), t)
end
function Base.show(io::IO, mime::MIME"text/plain", t::RationalPoly)
    print_ratpoly(io, mime, t)
end
function Base.show(io::IO, mime::MIME"text/latex", t::RationalPoly)
    print_ratpoly(io, mime, t)
end

function print_ratpoly(io::IO, mime, p::RationalPoly)
    print(io, "(")
    show(io, mime, p.num)
    print(io, ") / (")
    show(io, mime, p.den)
    print(io, ")")
end

function unicode_subscript(i)
    if i == 0
        "\u2080"
    elseif i == 1
        "\u2081"
    elseif i == 2
        "\u2082"
    elseif i == 3
        "\u2083"
    elseif i == 4
        "\u2084"
    elseif i == 5
        "\u2085"
    elseif i == 6
        "\u2086"
    elseif i == 7
        "\u2087"
    elseif i == 8
        "\u2088"
    elseif i == 9
        "\u2089"
    end
end

function unicode_superscript(i)
    if i == 0
        "\u2070"
    elseif i == 1
        "\u00B9"
    elseif i == 2
        "\u00B2"
    elseif i == 3
        "\u00B3"
    elseif i == 4
        "\u2074"
    elseif i == 5
        "\u2075"
    elseif i == 6
        "\u2076"
    elseif i == 7
        "\u2077"
    elseif i == 8
        "\u2078"
    elseif i == 9
        "\u2079"
    end
end
