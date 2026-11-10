# From MultivariateBases/src/MultivariateBases.jl — algebra() function

function algebra(basis::Union{FullBasis{B,V},SubBasis{B,V}}) where {B,V}
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

const _FullBasisAlgebra{B,V} = SA.StarAlgebra{
    Variables{B,V},
    <:Any,
    <:MStruct{B,V,<:Any,<:Any,<:FullBasis{B,V}},
}

# Merge the variables of all the operands at once instead of promoting their
# parents one by one which, e.g., merges `O(n)` variables `n` times for the
# sum of `n` terms with different variables.
function _common_parent(operands, ::Type{<:_FullBasisAlgebra{B}}) where {B}
    if isempty(operands)
        return nothing
    end
    vars = reduce(vcat, [variables(SA.object(parent(x))) for x in operands])
    return algebra(FullBasis{B}(unique!(sort!(vars; rev = true))))
end

function SA.common_parent(
    operands::AbstractArray{<:SA.Term{<:Any,A}},
) where {A<:_FullBasisAlgebra{<:Any,<:AbstractVector}}
    return _common_parent(operands, A)
end

function SA.common_parent(
    operands::AbstractArray{<:SA.AlgebraElement{<:Any,A}},
) where {A<:_FullBasisAlgebra{<:Any,<:AbstractVector}}
    return _common_parent(operands, A)
end

function MA.promote_operation(
    ::typeof(algebra),
    BT::Type{<:Union{FullBasis{B,V,E},SubBasis{B,V,E}}},
) where {B,V,E}
    P = BasisElement{B,V,E}
    MS = if is_commutative(V)
        MStruct{B,V,E,SA.key_type(BT),BT}
    else
        SA.DiracMStructure{P,SA.key_type(BT),BT,typeof(*)}
    end
    return SA.StarAlgebra{Variables{B,V},P,MS}
end
