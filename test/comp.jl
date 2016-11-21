facts("Monomial equality") do
    @polyvar x y
    @fact x*y == x --> false
    @fact x == Monomial(x) --> true
    @fact Monomial([x,y], [1,0]) == x --> true
    @fact x == Monomial([x,y], [0,1]) --> false
    @fact MonomialVector([x,y], [[0,0],[1,0]]) == MonomialVector([x], [[0],[1]]) --> true
end
facts("Graded Lex Order") do
    @polyvar x y z
    @fact x > y > z --> true
    @fact x^2*y > y^3 > z --> true
    # Examples from p. 58, 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
    @fact x^1*y^2*z^3 > x^3*y^2 --> true
    @fact x^1*y^2*z^3 < x^3*y^2 --> false
    @fact x^1*y^2*z^4 > x^1*y^1*z^5 --> true
    @fact x^1*y^2*z^4 < x^1*y^1*z^5 --> false
    @fact x^4*y^7*z > x^4*y^2*z^3 --> true
    @fact x^4*y^7*z < x^4*y^2*z^3 --> false
    @fact x*y^5*z^2 < x^4*y*z^3 --> true
    @fact x*y^5*z^2 > x^4*y*z^3 --> false
    @fact x^5*y*z > x^4*y*z^2 --> true
    @fact x^5*y*z < x^4*y*z^2 --> false
end
facts("Polynomial equality") do
    @polyvar x y
    @fact 2*x*y + 3*y^2 --> 3*y^2 + 2*y*x
    @fact 3*x*y + 2*y^2 != 3*y^2 + 2*y*x --> true
    @fact isapprox((2-1e-3)*x*y + (3+1e-3)*y^2, 3*y^2 + 2*y*x, rtol=1e-2) --> true
    @fact isapprox((2-1e-3)*x*y + (3+1e-1)*y^2, 3*y^2 + 2*y*x, rtol=1e-2) --> false
    @fact isapprox(1e-3*x*y + 3*y^2 + x^2, x^2 + 3*y^2, rtol=1e-2, ztol=1e-2) --> true
    @fact isapprox(3*y^2 + x^2, x^2 + 1e-3*x*y + 3*y^2, rtol=1e-2, ztol=1e-2) --> true
end
