facts("Differentiation") do
  @polyvar x y
  @fact differentiate(MatPolynomial{Int}((i,j)->1, [x]), y) --> 0
end
