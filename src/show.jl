const TypesWithShow = Union{
    AbstractVariable,
    AbstractMonomial,
    AbstractTerm,
    AbstractPolynomial,
    RationalPoly,
}
function Base.show(io::IO, mime::MIME"text/latex", p::TypesWithShow)
    print(io, "\$\$ ")
    _show(io, mime, p)
    return print(io, " \$\$")
end

# If the MIME is not specified, IJulia thinks that it supports images, ...
# and then use the result of show and tries to interpret it as an svg, ...
# We need the two methods to avoid ambiguity
function Base.show(io::IO, mime::MIME"text/plain", p::TypesWithShow)
    return _show(io, mime, p)
end
function Base.show(io::IO, mime::MIME"text/print", p::TypesWithShow)
    return _show(io, mime, p)
end

Base.print(io::IO, p::TypesWithShow) = show(io, MIME"text/print"(), p)
Base.show(io::IO, p::TypesWithShow) = show(io, MIME"text/plain"(), p)

# VARIABLES
function _show(io::IO, mime::MIME, var::AbstractVariable)
    base, indices = name_base_indices(var)
    if isempty(indices)
        print(io, base)
    else
        print(io, base)
        print_subscript(io, mime, indices)
    end
end
function _show(io::IO, mime::MIME"text/print", var::AbstractVariable)
    return print(io, name(var))
end

function print_subscript(io::IO, ::MIME"text/latex", index)
    return print(io, "_{", join(index, ","), "}")
end
function print_subscript(io::IO, mime, indices)
    if length(indices) == 1
        print(io, unicode_subscript(indices[1]))
    else
        print(io, join(unicode_subscript.(indices), "\u208B"))
    end
end

const unicode_subscripts = ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
unicode_subscript(i) = join(unicode_subscripts[d+1] for d in reverse(digits(i)))

# MONOMIALS
function _show(io::IO, mime, m::AbstractMonomial)
    if isconstant(m)
        print(io, '1')
    else
        printed_var = false
        vars = variables(m)
        n = length(vars)
        for (i, var, exp) in zip(1:n, vars, exponents(m))
            if !iszero(exp)
                if mime isa MIME"text/print" && printed_var && i > 0
                    print(io, "*")
                end
                _show(io, mime, var)
                printed_var = true
                if !isone(exp)
                    print_exponent(io, mime, exp)
                end
            end
        end
    end
end

print_exponent(io::IO, ::MIME"text/latex", exp) = print(io, "^{", exp, "}")
print_exponent(io::IO, ::MIME"text/print", exp) = print(io, "^", exp)
function print_exponent(io::IO, mime, exp)
    return print(io, join(unicode_superscript.(reverse(digits(exp)))))
end

const unicode_superscripts = ("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
unicode_superscript(i) = unicode_superscripts[i+1]

# TERM
function _show(io::IO, mime, t::AbstractTerm)
    if isconstant(t)
        print_coefficient(io, mime, coefficient(t))
    else
        if should_print_coefficient(coefficient(t))
            if !should_print_coefficient(-coefficient(t))
                print(io, '-')
            else
                print_coefficient(io, mime, coefficient(t))
                if !iszero(t)
                    # Print a multiplication sign between coefficent and monmomial
                    # depending on the mime type
                    print_maybe_multiplication_sign(io, mime)
                end
            end
        end
        if !iszero(t)
            _show(io, mime, monomial(t))
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
# `Int`, `Float64` don't support MIME"text/latex".
# We could add a check with `showable` if a `Real` subtype supports it and
# the feature is requested.
print_coefficient(io::IO, mime, coeff::Real) = print(io, coeff)
# Scientific notation does not display well in LaTeX so we rewrite it
function print_coefficient(io::IO, mime::MIME"text/latex", coeff::AbstractFloat)
    s = string(coeff)
    if occursin('e', s)
        s = replace(s, 'e' => " \\cdot 10^{") * '}'
    end
    return print(io, s)
end

function _trim_LaTeX(s::AbstractString)
    i = firstindex(s)
    j = lastindex(s)
    while true
        if i < j && isspace(s[i])
            i = nextind(s, i)
        elseif i < j && isspace(s[j])
            j = prevind(s, j)
        elseif i < j && s[i] == '$' && s[j] == '$'
            i = nextind(s, i)
            j = prevind(s, j)
        elseif i < j && ((s[i:i+1] == "\\(" && s[j-1:j] == "\\)") || (s[i:i+1] == "\\[" && s[j-1:j] == "\\]"))
            i = nextind(s, i, 2)
            j = prevind(s, j, 2)
        else
            return s[i:j]
        end
    end
end

# JuMP expressions supports LaTeX output so `showable` will return `true`
# for them. It is important for anonymous variables to display properly as well:
# https://github.com/jump-dev/SumOfSquares.jl/issues/256
# Since they add `$$` around it, we need to trim it with `_trim_LaTeX`
function print_coefficient(io::IO, mime, coeff)
    print(io, "(")
    if showable(mime, coeff)
        print(io, _trim_LaTeX(sprint(show, mime, coeff)))
    else
        show(io, coeff)
    end
    return print(io, ")")
end

# POLYNOMIAL
function _show(io::IO, mime, p::AbstractPolynomial{T}) where {T}
    ts = terms(p)
    if isempty(ts)
        print(io, zero(T))
    else
        _show(io, mime, first(ts))
        for t in Iterators.drop(ts, 1)
            if isnegative(coefficient(t))
                print(io, " - ")
                _show(io, mime, abs(coefficient(t)) * monomial(t))
            else
                print(io, " + ")
                _show(io, mime, t)
            end
        end
    end
end

isnegative(x::Real) = x < 0
isnegative(x) = false

# RATIONAL POLY
function _show(io::IO, mime, p::RationalPoly)
    print(io, mime isa MIME"text/latex" ? "\\frac{" : "(")
    _show(io, mime, p.num)
    print(io, mime isa MIME"text/latex" ? "}{" : ") / (")
    _show(io, mime, p.den)
    return print(io, mime isa MIME"text/latex" ? "}" : ")")
end
