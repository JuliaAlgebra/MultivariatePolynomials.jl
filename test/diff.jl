facts("Differentiation") do
    @polyvar x y
    @fact differentiate(MatPolynomial{Int}((i,j)->1, [x]), y) --> 0
    @fact differentiate(x*y + 3y^2 , [x, y]) --> [y, x+6y]
end
