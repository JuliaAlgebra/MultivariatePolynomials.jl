@deprecate constantmonomial constant_monomial
@deprecate convertconstant convert_constant
@deprecate termtype term_type
@deprecate polynomialtype polynomial_type
@deprecate jointerms join_terms
@deprecate jointerms! join_terms!
@deprecate mapcoefficients map_coefficients
@deprecate mapcoefficients! map_coefficients!
@deprecate mapexponents map_exponents
@deprecate mapexponents! map_exponents!
@deprecate monomialtype monomial_type
@deprecate monovectype monomial_vector_type
@deprecate emptymonovec empty_monomial_vector
@deprecate sortmonovec sort_monomial_vector
@deprecate mergemonovec merge_monomial_vectors
@deprecate monovec monomial_vector
@deprecate leading_term leadingterm
@deprecate leading_monomial leading_monomial
@deprecate leading_coefficient leading_coefficient
@deprecate removeleadingterm remove_leading_term
@deprecate removemonomials remove_monomials
@deprecate mergesorted merge_sorted
@deprecate mergesorted! merge_sorted!
@deprecate coefficienttype coefficient_type
@deprecate constantterm constant_term
@deprecate zeroterm zero_term
@deprecate similarvariable similar_variable
@deprecate pairzip pair_zip
@deprecate tuplezip tuple_zip

function changecoefficienttype(::Type{P}, ::Type{T}) where {P,T}
    Base.depwarn("`changecoefficienttype(P::Type, T::Type)` is deprecated, use `similar_type(P, T)` instead", :changecoefficienttype)
    return similar_type(P, T)
end

function changecoefficienttype(poly::APL, ::Type{T}) where {T}
    Base.depwarn("`changecoefficienttype(poly, T::Type)` is deprecated, use `similar(poly, T)` instead", :changecoefficienttype)
    return similar(poly, T)
end

function mapcoefficientsnz(f::F, p::APL) where {F<:Function}
    Base.depwarn("`mapcoefficientsnz(f, p)` is deprecated, use `map_coefficients(f, p, nonzero = true)` instead", :map_coefficientsnz)
    return map_coefficients(f, p, nonzero = true)
end

function mapcoefficientsnz_to!(output::APL, f::F, p::APL) where {F<:Function}
    Base.depwarn("`mapcoefficientsnz_to!(output, f, p)` is deprecated, use `map_coefficients_to!(output, f, p, nonzero = true)` instead", :map_coefficientsnz_to!)
    return map_coefficients_to!(output, f, p, nonzero = true)
end
