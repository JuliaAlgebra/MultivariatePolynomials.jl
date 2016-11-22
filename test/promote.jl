facts("Promotion") do
  @polyvar x y
  @fact typeof([x, x*y+x, x]) --> Vector{VecPolynomial{Int}}
  @fact typeof([1, x/y, x]) --> Vector{RationalPoly{Int, Int}}
  @fact typeof([(x^2-2x+2) x; x x^2]) --> Matrix{VecPolynomial{Int}}
end
