facts("Vector to Monomial Vector") do
    @polyvar x
    @fact MonomialVector([1]) --> MonomialVector(PolyVar[], [Int[]])
    @fact MonomialVector([x]) --> MonomialVector([x], [[1]])
end
