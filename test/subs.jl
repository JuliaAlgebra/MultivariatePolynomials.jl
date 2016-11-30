facts("Substitution") do
    @polyvar x[1:3]

    a = (x[1])([x[2]], [x[1]])
    b = x[2]
    @fact (x[1])([x[2]], [x[1]]) == x[2] --> true

    p = x[1]*x[2]*x[3]
    @fact Int(p([1, 2, 3], x)) --> 6

    p = x[1]^2 + x[1]*x[3] - 3
    @fact Int(p([5, x[1]*x[2], 4], x)) --> 42

    p = x[1]^2 + x[2]^2
    q = p([1 -1; 1 1] * x[1:2], x[1:2])
    @fact q == 2p --> true

    q = (x[1] + 1) / (x[1] + 2)
    @fact isapproxzero(q([-1], [x[1]])) --> true
    @fact isapproxzero(q([1], [x[1]])) --> false
    @fact isapprox(q([1], [x[1]]), 2/3) --> true

    P = [1 2 3; 2 4 5; 3 5 6]
    p = MatPolynomial(P, x)
    @fact p(ones(3), x) --> 31
end
