export iscomplex, isrealpart, isimagpart, isconj, ordinary_variable, degree_complex, halfdegree, mindegree_complex,
    minhalfdegree, maxdegree_complex, maxhalfdegree, extdegree_complex, exthalfdegree

"""
    iscomplex(x::AbstractVariable)

Return whether a given variable was declared as a complex-valued variable (also their conjugates are complex, but their real
and imaginary parts are not).
By default, all variables are real-valued.
"""
iscomplex(::AbstractVariable) = false

Base.isreal(v::Union{<:AbstractPolynomialLike,<:AbstractVector{<:AbstractMonomial}}) = !iscomplex(v)

"""
    isrealpart(x::AbstractVariable)

Return whether the given variable is the real part of a complex-valued variable.

See also [`iscomplex`](@ref iscomplex), [`isimagpart`](@ref isimagpart), [`isconj`](@ref isconj).
"""
isrealpart(::AbstractVariable) = false

"""
    isimagpart(x::AbstractVariable)

Return whether the given variable is the imaginary part of a complex-valued variable.

See also [`iscomplex`](@ref iscomplex), [`isrealpart`](@ref isrealpart), [`isconj`](@ref isconj).
"""
isimagpart(::AbstractVariable) = false

"""
    isconj(x::AbstractVariable)

Return whether the given variable is obtained by conjugating a user-defined complex-valued variable.

See also [`iscomplex`](@ref iscomplex), [`isrealpart`](@ref isrealpart), [`isimagpart`](@ref isimagpart).
"""
isconj(::AbstractVariable) = false

"""
    ordinary_variable(x::Union{AbstractVariable, AbstractVector{<:AbstractVariable}})

Given some (complex-valued) variable that was transformed by conjugation, taking its real part, or taking its
imaginary part, return the original variable as it was defined by the user.

See also [`conj`](@ref conj), [`real`](@ref), [`imag`](@ref).
"""
ordinary_variable(x::AbstractVariable) = x

"""
    conj(x::AbstractVariable)

Return the complex conjugate of a given variable if it was declared as a complex variable; else return the
variable unchanged.

    conj(x::AbstractMonomial)
    conj(x::AbstractVector{<:AbstractMonomial})
    conj(x::AbstractTerm)
    conj(x::AbstractPolynomial)

Return the complex conjugate of `x` by applying conjugation to all coefficients and variables.

See also [`iscomplex`](@ref iscomplex), [`isconj`](@ref isconj).
"""
Base.conj(x::AbstractVariable) = x

"""
    real(x::AbstractVariable)

Return the real part of a given variable if it was declared as a complex variable; else return the variable
unchanged.

    real(x::AbstractMonomial)
    real(x::AbstractVector{<:AbstractMonomial})
    real(x::AbstractTerm)
    real(x::AbstractPolynomial)

Return the real part of `x` by applying `real` to all coefficients and variables; for this purpose, every complex-valued
variable is decomposed into its real- and imaginary parts.

See also [`iscomplex`](@ref iscomplex), [`isrealpart`](@ref isrealpart), [`imag`](@ref).
"""
Base.real(x::AbstractVariable) = x

"""
    imag(x::AbstractVariable)

Return the imaginary part of a given variable if it was declared as a complex variable; else return zero.

    imag(x::AbstractMonomial)
    imag(x::AbstractVector{<:AbstractMonomial})
    imag(x::AbstractTerm)
    imag(x::AbstractPolynomial)

Return the imaginary part of `x` by applying `imag` to all coefficients and variables; for this purpose, every complex-valued
variable is decomposed into its real- and imaginary parts.

See also [`iscomplex`](@ref iscomplex), [`isimagpart`](@ref isimagpart), [`real`](@ref).
"""
Base.imag(::AbstractVariable) = MA.Zero()

# extend to higher-level elements. We make all those type-stable (but we need convert, as the construction method may
# deliver simpler types than the inputs if they were deliberately casted, e.g., term to monomial)
function iscomplex(p::AbstractPolynomialLike)
    for v in variables(p)
        if !iszero(maxdegree(p, v)) && iscomplex(v)
            return true
        end
    end
    return false
end
MP.iscomplex(p::AbstractVector{<:AbstractMonomial}) = any(iscomplex, p)

Base.conj(x::M) where {M<:AbstractMonomial} = isreal(x) ? x :
    convert(M, reduce(*, conj(var)^exp for (var, exp) in zip(variables(x), exponents(x)); init=constantmonomial(x)))
Base.conj(x::V) where {V<:AbstractVector{<:AbstractMonomial}} = isreal(x) ? x : monovec(conj.(x))
Base.conj(x::T) where {T<:AbstractTerm} = isreal(x) ? x : convert(T, conj(coefficient(x)) * conj(monomial(x)))
Base.conj(x::P) where {P<:AbstractPolynomial} = iszero(nterms(x)) || isreal(x) ? x : convert(P, sum(conj(t) for t in x))

# Real and imaginary parts are harder to realize. The real part of a monomial can easily be a polynomial.
for fun in [:real, :imag]
    eval(quote
        function Base.$fun(x::APL)
            iszero(nterms(x)) && return zero(polynomialtype(x, real(coefficienttype(x))))
            # We replace every complex variable by its decomposition into real and imaginary part
            substs = [var => real(var) + 1im * imag(var) for var in variables(x) if iscomplex(var)]
            # subs throws an error if it doesn't receive at least one substitution
            full_version = isempty(substs) ? polynomial(x) : subs(x, substs...)
            # Now everything that is imaginary can be recognized by its coefficient
            return convert(
                polynomialtype(full_version, real(coefficienttype(full_version))),
                mapcoefficients!($fun, full_version)
            )
        end
        Base.$fun(x::AbstractVector{<:AbstractMonomial}) = map(Base.$fun, x)
    end)
end

# Also give complex-valued degree definitions. We choose not to overwrite degree, as this will lead to issues in monovecs
# and their sorting. So now there are two ways to calculate degrees: strictly by considering all variables independently,
# and also by looking at their complex structure.
"""
    degree_complex(t::AbstractTermLike)

Return the _total complex degree_ of the monomial of the term `t`, i.e., the maximum of the total degree of the declared
variables in `t` and the total degree of the conjugate variables in `t`.
To be well-defined, the monomial must not contain real parts or imaginary parts of variables.

    degree_complex(t::AbstractTermLike, v::AbstractVariable)

Returns the exponent of the variable `v` or its conjugate in the monomial of the term `t`, whatever is larger.

See also [`isconj`](@ref isconj).
"""
function degree_complex(t::AbstractTermLike)
    vars = variables(t)
    @assert(!any(isrealpart, vars) && !any(isimagpart, vars))
    grouping = isconj.(vars)
    exps = exponents(t)
    return max(sum(exps[grouping]), sum(exps[map(!, grouping)]))
end
degree_complex(t::AbstractTermLike, var::AbstractVariable) = degree_complex(monomial(t), var)
degree_complex(v::AbstractVariable, var::AbstractVariable) = (v == var ? 1 : 0)
function degree_complex(m::AbstractMonomial, v::AbstractVariable)
    deg = 0
    deg_c = 0
    c_v = conj(v)
    for (var, exp) in powers(m)
        @assert(!isrealpart(var) && !isimagpart(var))
        if var == v
            deg += exp
        elseif var == c_v
            deg_c += exp
        end
    end
    return max(deg, deg_c)
end

"""
    halfdegree(t::AbstractTermLike)

Return the equivalent of `ceil(degree(t)/2)`` for real-valued terms or `degree_complex(t)` for terms with only complex
variables; however, respect any mixing between complex and real-valued variables.
"""
function halfdegree(t::AbstractTermLike)
    realdeg = 0
    cpdeg = 0
    conjdeg = 0
    for (var, exp) in powers(t)
        if iscomplex(var)
            if isconj(var)
                conjdeg += exp
            else
                @assert(!isrealpart(var) && !isimagpart(var))
                cpdeg += exp
            end
        else
            realdeg += exp
        end
    end
    return ((realdeg + 1) >> 1) + max(cpdeg, conjdeg)
end

"""
    mindegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the minimal total complex degree of the monomials of `p`, i.e., `minimum(degree_complex, terms(p))`.

    mindegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Return the minimal complex degree of the monomials of `p` in the variable `v`, i.e., `minimum(degree_complex.(terms(p), v))`.
"""
mindegree_complex(X::AbstractVector{<:AbstractTermLike}, args...) =
    isempty(X) ? 0 : minimum(t -> degree_complex(t, args...), X, init=0)
mindegree_complex(p::AbstractPolynomialLike, args...) = mindegree_complex(terms(p), args...)

"""
    minhalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the minmal half degree of the monomials of `p`, i.e., `minimum(halfdegree, terms(p))`
"""
minhalfdegree(X::AbstractVector{<:AbstractTermLike}, args...) = isempty(X) ? 0 : minimum(halfdegree, X, init=0)
minhalfdegree(p::AbstractPolynomialLike) = minhalfdegree(terms(p))

"""
    maxdegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the maximal total complex degree of the monomials of `p`, i.e., `maximum(degree_complex, terms(p))`.

    maxdegree_complex(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Return the maximal complex degree of the monomials of `p` in the variable `v`, i.e., `maximum(degree_complex.(terms(p), v))`.
"""
maxdegree_complex(X::AbstractVector{<:AbstractTermLike}, args...) =
    mapreduce(t -> degree_complex(t, args...), max, X, init=0)
maxdegree_complex(p::AbstractPolynomialLike, args...) = maxdegree_complex(terms(p), args...)

"""
    maxhalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the maximal half degree of the monomials of `p`, i.e., `maximum(halfdegree, terms(p))`
"""
maxhalfdegree(X::AbstractVector{<:AbstractTermLike}, args...) = isempty(X) ? 0 : maximum(halfdegree, X, init=0)
maxhalfdegree(p::AbstractPolynomialLike) = maxhalfdegree(terms(p))

"""
    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Returns the extremal total complex degrees of the monomials of `p`, i.e., `(mindegree_complex(p), maxdegree_complex(p))`.

    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Returns the extremal complex degrees of the monomials of `p` in the variable `v`, i.e.,
`(mindegree_complex(p, v), maxdegree_complex(p, v))`.
"""
extdegree_complex(p::Union{AbstractPolynomialLike,AbstractVector{<:AbstractTermLike}}, args...) =
    (mindegree_complex(p, args...), maxdegree_complex(p, args...))

"""
    exthalfdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Return the extremal half degree of the monomials of `p`, i.e., `(minhalfdegree(p), maxhalfdegree(p))`
"""
exthalfdegree(p::Union{AbstractPolynomialLike,AbstractVector{<:AbstractTermLike}}) = (minhalfdegree(p), maxhalfdegree(p))

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