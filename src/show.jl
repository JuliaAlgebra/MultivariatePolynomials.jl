Base.print(io::IO, v::AbstractPolynomialLike) = show(io, MIME"text/print"(), v)

# VARIABLES
Base.show(io::IO, v::AbstractVariable) = print_var(io, MIME"text/plain"(), v)
Base.show(io::IO, mime::MIME"text/plain", v::AbstractVariable) = print_var(io, mime, v)
Base.show(io::IO, mime::MIME"text/latex", v::AbstractVariable) = print_var(io, mime, v)
Base.show(io::IO, mime::MIME"text/print", v::AbstractVariable) = print_var(io, mime, v)

function print_var(io::IO, mime::MIME, var::AbstractVariable)
    base, indices = name_base_indices(var)
    if isempty(indices)
        print(io, base)
    else
        print(io, base)
        print_subscript(io, mime, indices)
    end
end
print_var(io::IO, mime::MIME"text/print", var::AbstractVariable) = print(io, name(var))

function print_subscript(io::IO, ::MIME"text/latex", index)
    print(io, "_{", join(index, ","), "}")
end
function print_subscript(io::IO, mime, indices)
    if length(indices) == 1
        print(io, unicode_subscript(indices[1]))
    else
        print(io, join(unicode_subscript.(indices), "\u208B"))
    end
end

const unicode_subscripts = ("₀","₁","₂","₃","₄","₅","₆","₇","₈","₉")
unicode_subscript(i) = join(unicode_subscripts[d+1] for d in reverse(digits(i)))

# MONOMIALS

Base.show(io::IO, mime::MIME"text/latex", m::AbstractMonomial) = print_monomial(io, mime, m)
Base.show(io::IO, mime::MIME"text/plain", m::AbstractMonomial) = print_monomial(io, mime, m)
Base.show(io::IO, mime::MIME"text/print", m::AbstractMonomial) = print_monomial(io, mime, m)
Base.show(io::IO, m::AbstractMonomial) = print_monomial(io, MIME"text/plain"(), m)

function print_monomial(io::IO, mime, m::AbstractMonomial)
    if isconstant(m)
        print(io, '1')
    else
        printed_var = false
        vars = variables(m)
        n = length(vars)
        for (i, var, exp) in zip(1:n, vars, exponents(m))
            if !iszero(exp)
                if mime isa MIME"text/print" && printed_var && i > 0 &&
                    print(io,"*")
                end
                print_var(io, mime, var)
                printed_var = true
                if !isone(exp)
                    print_exponent(io, mime, exp)
                end
            end
        end
    end
end
#
print_exponent(io::IO, ::MIME"text/latex", exp) = print(io, "^{", exp, "}")
print_exponent(io::IO, ::MIME"text/print", exp) = print(io, "^", exp)
function print_exponent(io::IO, mime, exp)
    print(io, join(unicode_superscript.(reverse(digits(exp)))))
end

const unicode_superscripts = ("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
unicode_superscript(i) = unicode_superscripts[i+1]

# TERM

Base.show(io::IO, t::AbstractTerm) = print_term(io, MIME"text/plain"(), t)
Base.show(io::IO, mime::MIME"text/latex", t::AbstractTerm) = print_term(io, mime, t)
Base.show(io::IO, mime::MIME"text/plain", t::AbstractTerm) = print_term(io, mime, t)
Base.show(io::IO, mime::MIME"text/print", t::AbstractTerm) = print_term(io, mime, t)

function print_term(io::IO, mime, t::AbstractTerm)
    if isconstant(t)
        print_coefficient(io, coefficient(t))
    else
        if should_print_coefficient(coefficient(t))
            if !should_print_coefficient(-coefficient(t))
                print(io, '-')
            else
                print_coefficient(io, coefficient(t))
                if !iszero(t)
                    # Print the a multiplication sign between coefficent and monmomial
                    # depending on the mime type
                    print_maybe_multiplication_sign(io, mime)
                end
            end
        end
        if !iszero(t)
            show(io, mime, monomial(t))
        end
    end
end

"""
    print_maybe_multiplication_sign(io, mime)

Prints a multiplication sign depending on the `mime` type.
"""
print_maybe_multiplication_sign(io::IO, ::MIME"text/print") = print(io, "*")
print_maybe_multiplication_sign(io::IO, mime) = nothing

should_print_coefficient(x) = true  # By default, just print all coefficients
should_print_coefficient(x::Number) = !isone(x) # For numbers, we omit any "one" coefficients
print_coefficient(io::IO, coeff::Real) = print(io, coeff)
print_coefficient(io::IO, coeff) = print(io, "(", coeff, ")")

# POLYNOMIAL

Base.show(io::IO, t::AbstractPolynomial) = print_poly(io, MIME"text/plain"(), t)
Base.show(io::IO, mime::MIME"text/plain", t::AbstractPolynomial) = print_poly(io, mime, t)
Base.show(io::IO, mime::MIME"text/latex", t::AbstractPolynomial) = print_poly(io, mime, t)
Base.show(io::IO, mime::MIME"text/print", t::AbstractPolynomial) = print_poly(io, mime, t)

function print_poly(io::IO, mime, p::AbstractPolynomial{T}) where T
    ts = terms(p)
    if isempty(ts)
        print(io, zero(T))
    else
        print_term(io, mime, first(ts))
        for t in Iterators.drop(ts, 1)
            if isnegative(coefficient(t))
                print(io, " - ")
                show(io, mime, abs(coefficient(t)) * monomial(t))
            else
                print(io, " + ")
                show(io, mime, t)
            end
        end
    end
end

isnegative(x::Real) = x < 0
isnegative(x) = false

# RATIONAL POLY

Base.show(io::IO, t::RationalPoly) = print_ratpoly(io, MIME"text/plain"(), t)
Base.show(io::IO, mime::MIME"text/plain", t::RationalPoly) = print_ratpoly(io, mime, t)
Base.show(io::IO, mime::MIME"text/latex", t::RationalPoly) = print_ratpoly(io, mime, t)
Base.show(io::IO, mime::MIME"text/print", t::RationalPoly) = print_ratpoly(io, mime, t)

function print_ratpoly(io::IO, mime, p::RationalPoly)
    print(io, "(")
    show(io, mime, p.num)
    print(io, ") / (")
    show(io, mime, p.den)
    print(io, ")")
end
