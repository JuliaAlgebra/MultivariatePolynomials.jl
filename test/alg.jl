facts("Algebra") do
  @polyvar x
  @fact 2 .- ((1.+(-x)) .* 4) ./ 2 == x.^2 .* (1 ./ x) .* 2 --> true
  @fact dot(0, x^2 - 2*x^2) --> dot((x^2 - x)', x^2 - x^2)
end
