# From MultivariateBases/src/MultivariateBases.jl — algebra() function

function algebra(
    basis::Union{FullBasis{B,V},SubBasis{B,V}},
) where {B,V}
    mstruct = if is_commutative(V)
        MStruct(basis)
    else
        if B != Monomial
            error(
                "Only the `Monomial` basis is supported with noncommutative variables",
            )
        end
        SA.DiracMStructure(basis, *)
    end
    return SA.StarAlgebra(Variables{B}(variables(basis)), mstruct)
end

function MA.promote_operation(
    ::typeof(algebra),
    BT::Type{
        <:Union{
            FullBasis{B,V,E},
            SubBasis{B,V,E},
        },
    },
) where {B,V,E}
    P = Polynomial{B,V,E}
    MS = if is_commutative(V)
        MStruct{B,V,E,SA.key_type(BT),BT}
    else
        SA.DiracMStructure{P,SA.key_type(BT),BT,typeof(*)}
    end
    return SA.StarAlgebra{Variables{B,V},P,MS}
end
