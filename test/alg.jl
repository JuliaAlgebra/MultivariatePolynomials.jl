facts("Algebra") do
  @polyvar x
  @fact 2 .- ((1.+(-x)) .* 4) ./ 2 == x.^2 .* (1 ./ x) .* 2 --> true
end
