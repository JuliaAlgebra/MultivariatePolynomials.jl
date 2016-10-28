facts("Comparison") do
  @polyvar x y
  @fact x*y == x --> false
  @fact x == Monomial(x) --> true
  @fact Monomial([x,y], [1,0]) == x --> true
  @fact x == Monomial([x,y], [0,1]) --> false
  @fact MonomialVector([x,y], [[0,0],[1,0]]) == MonomialVector([x], [[0],[1]]) --> true
end
