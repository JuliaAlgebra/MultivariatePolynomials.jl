# From MultivariateBases/src/MultivariateBases.jl — algebra() function

function algebra(basis::Union{FullBasis{B,V},SubBasis{B,V}}) where {B,V}
    if is_commutative(V)
        return SA.StarAlgebra(Variables{B}(variables(basis)), MStruct(basis))
    end
    if B != Monomial
        error(
            "Only the `Monomial` basis is supported with noncommutative variables",
        )
    end
    # Use the parent basis as object, like `algebra(::_NoncommutativeBasis)`,
    # so that this algebra can be promoted with the one of the parent basis.
    return SA.StarAlgebra(implicit_basis(basis), SA.DiracMStructure(basis, *))
end

function MA.promote_operation(
    ::typeof(algebra),
    BT::Type{<:Union{FullBasis{B,V,E},SubBasis{B,V,E}}},
) where {B,V,E}
    P = BasisElement{B,V,E}
    if is_commutative(V)
        MS = MStruct{B,V,E,SA.key_type(BT),BT}
        return SA.StarAlgebra{Variables{B,V},P,MS}
    end
    O = MA.promote_operation(implicit_basis, BT)
    MS = SA.DiracMStructure{P,SA.key_type(BT),BT,typeof(*)}
    return SA.StarAlgebra{O,P,MS}
end
