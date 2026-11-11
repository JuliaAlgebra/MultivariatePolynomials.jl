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
    if is_commutative(V)
        MS = MStruct{B,V,E,SA.key_type(BT),BT}
        return SA.StarAlgebra{Variables{B,V},P,MS}
    end
    O = MA.promote_operation(implicit_basis, BT)
    MS = SA.DiracMStructure{P,SA.key_type(BT),BT,typeof(*)}
    return SA.StarAlgebra{O,P,MS}
end
