facts("Monomial") do
    @polyvar x
    @fact_throws ArgumentError Monomial([x], [1,0])
end
facts("Vector to Monomial Vector") do
    @polyvar x
    @fact MonomialVector([1]) --> MonomialVector(PolyVar[], [Int[]])
    @fact MonomialVector([x]) --> MonomialVector([x], [[1]])
end
