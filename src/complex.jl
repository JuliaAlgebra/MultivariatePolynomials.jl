export isrealpart,
    isimagpart,
    isconj,
    ordinary_variable,
    degree_complex,
    halfdegree,
    mindegree_complex,
    minhalfdegree,
    maxdegree_complex,
    maxhalfdegree,
    extdegree_complex,
    exthalfdegree

"""
    isreal(x::AbstractVariable)

Return whether a given variable was declared as a real-valued or complex-valued variable (also their conjugates are complex,
but their real and imaginary parts are not).
By default, all variables are real-valued.
"""
Base.isreal(::AbstractVariable) = true

"""
    isrealpart(x::AbstractVariable)

Return whether the given variable is the real part of a complex-valued variable.

See also [`isreal`](@ref), [`isimagpart`](@ref), [`isconj`](@ref).
"""
isrealpart(::AbstractVariable) = false

"""
    isimagpart(x::AbstractVariable)

Return whether the given variable is the imaginary part of a complex-valued variable.

See also [`isreal`](@ref), [`isrealpart`](@ref), [`isconj`](@ref).
"""
isimagpart(::AbstractVariable) = false

"""
    isconj(x::AbstractVariable)

Return whether the given variable is obtained by conjugating a user-defined complex-valued variable.

See also [`isreal`](@ref), [`isrealpart`](@ref), [`isimagpart`](@ref).
"""
isconj(::AbstractVariable) = false

"""
    ordinary_variable(x::Union{AbstractVariable, AbstractVector{<:AbstractVariable}})

Given some (complex-valued) variable that was transformed by conjugation, taking its real part, or taking its
imaginary part, return the original variable as it was defined by the user.

See also [`conj`](@ref), [`real`](@ref), [`imag`](@ref).
"""
ordinary_variable(x::AbstractVariable) = MA.copy_if_mutable(x)

"""
    conj(x::AbstractVariable)

Return the complex conjugate of a given variable if it was declared as a complex variable; else return the
variable unchanged.

See also [`isreal`](@ref), [`isconj`](@ref).
"""
Base.conj(x::AbstractVariable) = MA.copy_if_mutable(x)
@doc """
    conj(x::AbstractPolynomialLike)

Return the complex conjugate of `x` by applying conjugation to all coefficients and variables.
""" conj(::_APL)
@doc """
    conj(x::AbstractVector{<:AbstractMonomial})

Return the complex conjugate of `x` by applying conjugation to monomials.
""" conj(::AbstractVector{<:AbstractMonomial})

"""
    real(x::AbstractVariable)

Return the real part of a given variable if it was declared as a complex variable; else return the variable
unchanged.

See also [`imag`](@ref).
"""
Base.real(x::AbstractVariable) = MA.copy_if_mutable(x)
@doc """
    real(x::AbstractPolynomialLike)

Return the real part of `x` by applying [`real`](@ref) to all coefficients and variables; for this purpose, every
complex-valued variable is decomposed into its real- and imaginary parts.

See also [`imag`](@ref).
""" real(::_APL)
@doc """
    real(x::AbstractVector{<:AbstractMonomial})

Return the real part of `x` by applying [`real`](@ref) to all monomials; for this purpose, every complex-valued variable is
decomposed into its real- and imaginary parts. Note that the result will no longer be a monomial vector.

See also [`imag`](@ref).
""" real(::AbstractVector{<:AbstractMonomial})

"""
    imag(x::AbstractVariable)

Return the imaginary part of a given variable if it was declared as a complex variable; else return zero.

See also [`isreal`](@ref), [`isimagpart`](@ref), [`real`](@ref).
"""
Base.imag(::AbstractVariable) = MA.Zero()
@doc """
    imag(x::AbstractPolynomialLike)

Return the imaginary part of `x` by applying [`imag`](@ref) to all coefficients and variables; for this purpose, every
complex-valued variable is decomposed into its real- and imaginary parts.

See also [`real`](@ref).
""" imag(::_APL)
@doc """
    imag(x::AbstractVector{<:AbstractMonomial})

Return the imaginary part of `x` by applying [`imag`](@ref) to all monomials; for this purpose, every complex-valued variable
is decomposed into its real- and imaginary parts. Note that the result will no longer be a monomial vector.

See also [`real`](@ref).
""" imag(::AbstractVector{<:AbstractMonomial})

# extend to higher-level elements. We make all those type-stable (but we need convert, as the construction method may
# deliver simpler types than the inputs if they were deliberately casted, e.g., term to monomial)
"""
    isreal(p::AbstractPolynomialLike)

Returns `true` if and only if no single variable in `p` was declared as a complex variable (in the sense that [`isreal`](@ref)
applied on them would be `true`) and no coefficient is complex-valued.
"""
function Base.isreal(p::_APL)
    for v in variables(p)
        if !isreal(v) && !iszero(maxdegree(p, v))
            return false
        end
    end
    return all(isreal, coefficients(p))
end
"""
    isreal(p::AbstractVector{<:AbstractMonomial})

Returns `true` if and only if every single monomial in `p` would is real-valued.
"""
Base.isreal(p::AbstractVector{<:AbstractMonomial}) = all(isreal, p)

function Base.conj(x::M) where {M<:AbstractMonomial}
    return isreal(x) ? x :
           convert(
        M,
        reduce(
            *,
            conj(var)^exp for (var, exp) in powers(x);
            init = constant_monomial(x),
        ),
    )
end
function Base.conj(x::V) where {V<:AbstractVector{<:AbstractMonomial}}
    return isreal(x) ? MA.copy_if_mutable(x) : monomial_vector(conj.(x))
end
function Base.conj(x::T) where {T<:AbstractTerm}
    return isreal(x) ? MA.copy_if_mutable(x) :
           convert(T, conj(coefficient(x)) * conj(monomial(x)))
end
function Base.conj(x::P) where {P<:AbstractPolynomial}
    return iszero(x) || isreal(x) ? MA.copy_if_mutable(x) :
           convert(P, polynomial([conj(t) for t in x]))
end

# Real and imaginary parts are harder to realize. The real part of a monomial can easily be a polynomial.
for fun in [:real, :imag]
    eval(
        quote
            function Base.$fun(x::_APL)
                # Note: x may also be a polynomial; we could handle this case separately by performing the replacement in each
                # individual term, then adding them all up. This could potentially lower the overall memory requirement (in case
                # the expansions of the individual terms simplify) at the expense of not being able to exploit optimizations that
                # subst can do by knowing the full polynomial.
                iszero(x) &&
                    return zero(polynomial_type(x, real(coefficient_type(x))))
                # We replace every complex variable by its decomposition into real and imaginary part
                subst_vars = filter(Base.:! ∘ isreal, variables(x))
                # To avoid a stack overflow on promote_type, we'll handle the empty case separately
                full_version =
                    isempty(subst_vars) ? polynomial(x) :
                    subs(
                        x,
                        subst_vars =>
                            [real(var) + 1im * imag(var) for var in subst_vars],
                    )
                # Now everything that is imaginary can be recognized by its coefficient
                return convert(
                    polynomial_type(
                        full_version,
                        real(coefficient_type(full_version)),
                    ),
                    map_coefficients!($fun, full_version),
                )
            end
            function Base.$fun(x::AbstractVector{<:AbstractMonomial})
                return map(Base.$fun, x)
            end
        end,
    )
end

# Also give complex-valued degree definitions. We choose not to overwrite degree, as this will lead to issues in monovecs
# and their sorting. So now there are two ways to calculate degrees: strictly by considering all variables independently,
# and also by looking at their complex structure.
for fn in (:degree_complex, :halfdegree)
    @eval function $fn(t::AbstractTermLike)
        realdeg = 0
        cpdeg = 0
        conjdeg = 0
        for (var, exp) in powers(t)
            if isreal(var)
                realdeg += exp
                (isrealpart(var) || isimagpart(var)) && error(
                    "Cannot calculate complex degrees when real or imaginary parts are present",
                )
            else
                if isconj(var)
                    conjdeg += exp
                else
                    cpdeg += exp
                end
            end
        end
        return $(
            fn === :degree_complex ? :(realdeg) : :(div(realdeg, 2, RoundUp))
        ) + max(cpdeg, conjdeg)
    end
end

"""
    degree_complex(t::AbstractTermLike)

Return the _total complex degree_ of the monomial of the term `t`, i.e., the maximum of the total degree of the declared
variables in `t` and the total degree of the conjugate variables in `t`.
To be well-defined, the monomial must not contain real parts or imaginary parts of variables.
If `x` is a real-valued variable and `z` is complex-valued,
- `degree_complex(x^5) = 5`
- `degree_complex(z^3 * conj(z)^4) = 4` and `degree_complex(z^4 * conj(z)^3) = 4`
- `degree_complex(x^5 * z^3 * conj(z^4)) = 5 + 4 = 9`
"""
degree_complex(t::AbstractTermLike)

"""
    halfdegree(t::AbstractTermLike)

Return the equivalent of `ceil(degree(t)/2)`` for real-valued terms or `degree_complex(t)` for terms with only complex
variables; however, respect any mixing between complex and real-valued variables.
To be well-defined, the monomial must not contain real parts or imaginary parts of variables.
If `x` is a real-valued variable and `z` is complex-valued,
- `halfdegree(x^5) = 3`
- `halfdegree(z^3 * conj(z)^4) = 4` and `halfdegree(z^4 * conj(z)^3) = 4`
- `halfdegree(x^5 * z^3 * conj(z^4)) = 3 + 4 = 7`
"""
halfdegree(t::AbstractTermLike)

"""
    degree_complex(t::AbstractTermLike, v::AbstractVariable)

Returns the exponent of the variable `v` or its conjugate in the monomial of the term `t`, whatever is larger.

See also [`isconj`](@ref).
"""
function degree_complex(t::AbstractTermLike, var::AbstractVariable)
    return degree_complex(monomial(t), var)
end
degree_complex(v::AbstractVariable, var::AbstractVariable) = (v == var ? 1 : 0)
function degree_complex(m::AbstractMonomial, v::AbstractVariable)
    deg = 0
    deg_c = 0
    c_v = conj(v)
    for (var, exp) in powers(m)
        (isrealpart(var) || isimagpart(var)) && error(
            "Cannot calculate complex degrees when real or imaginary parts are present",
        )
        if var == v
            deg += exp
        elseif var == c_v
            deg_c += exp
        end
    end
    return max(deg, deg_c)
end

"""
    mindegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the minimal total complex degree of the monomials of `p`, i.e., `minimum(degree_complex, terms(p))`.

    mindegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Return the minimal complex degree of the monomials of `p` in the variable `v`, i.e., `minimum(degree_complex.(terms(p), v))`.
"""
function mindegree_complex(X::AbstractVector{<:AbstractTermLike}, args...)
    return isempty(X) ? 0 :
           minimum(t -> degree_complex(t, args...), X, init = 0)
end
function mindegree_complex(p::_APL, args...)
    return mindegree_complex(terms(p), args...)
end

"""
    minhalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the minmal half degree of the monomials of `p`, i.e., `minimum(halfdegree, terms(p))`
"""
function minhalfdegree(X::AbstractVector{<:AbstractTermLike}, args...)
    return isempty(X) ? 0 : minimum(halfdegree, X, init = 0)
end
minhalfdegree(p::_APL) = minhalfdegree(terms(p))

"""
    maxdegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the maximal total complex degree of the monomials of `p`, i.e., `maximum(degree_complex, terms(p))`.

    maxdegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Return the maximal complex degree of the monomials of `p` in the variable `v`, i.e., `maximum(degree_complex.(terms(p), v))`.
"""
function maxdegree_complex(X::AbstractVector{<:AbstractTermLike}, args...)
    return mapreduce(t -> degree_complex(t, args...), max, X, init = 0)
end
function maxdegree_complex(p::_APL, args...)
    return maxdegree_complex(terms(p), args...)
end

"""
    maxhalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the maximal half degree of the monomials of `p`, i.e., `maximum(halfdegree, terms(p))`
"""
function maxhalfdegree(X::AbstractVector{<:AbstractTermLike}, args...)
    return isempty(X) ? 0 : maximum(halfdegree, X, init = 0)
end
maxhalfdegree(p::_APL) = maxhalfdegree(terms(p))

"""
    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Returns the extremal total complex degrees of the monomials of `p`, i.e., `(mindegree_complex(p), maxdegree_complex(p))`.

    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Returns the extremal complex degrees of the monomials of `p` in the variable `v`, i.e.,
`(mindegree_complex(p, v), maxdegree_complex(p, v))`.
"""
function extdegree_complex(
    p::Union{_APL,AbstractVector{<:AbstractTermLike}},
    args...,
)
    return (mindegree_complex(p, args...), maxdegree_complex(p, args...))
end

"""
    exthalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the extremal half degree of the monomials of `p`, i.e., `(minhalfdegree(p), maxhalfdegree(p))`
"""
function exthalfdegree(p::Union{_APL,AbstractVector{<:AbstractTermLike}})
    return (minhalfdegree(p), maxhalfdegree(p))
end

function ordinary_variable(x::AbstractVector{<:AbstractVariable})
    # let's assume the number of elements in x is small, else a conversion to a dict (probably better OrderedDict) would
    # be better
    results = similar(x, 0)
    sizehint!(results, length(x))
    j = 0
    for el in x
        ov = ordinary_variable(el)
        found = false
        for i in 1:j
            if results[i] == ov
                found = true
                break
            end
        end
        if !found
            push!(results, ov)
            j += 1
        end
    end
    return results
end
