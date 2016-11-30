facts("Monomial") do
    @polyvar x
    @fact_throws ArgumentError Monomial([x], [1,0])
end
facts("Degree to Monomial Vector") do
    @polyvar x y z
    @fact monomials([x,y], 2, z->z[1] <= 1) --> [x*y,y^2]
    @fact MonomialVector([x,y], 1:3) --> MonomialVector([x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y])
    X = MonomialVector([x,y,z], 5)
    @fact issorted(X, rev=true) --> true
    @fact X --> MonomialVector([x^5, x^4*y, x^4*z, x^3*y^2, x^3*y*z, x^3*z^2, x^2*y^3, x^2*y^2*z, x^2*y*z^2, x^2*z^3, x*y^4, x*y^3*z, x*y^2*z^2, x*y*z^3, x*z^4, y^5, y^4*z, y^3*z^2, y^2*z^3, y*z^4, z^5])
end
facts("Vector to Monomial Vector") do
    @polyvar x
    @fact MonomialVector([1]) --> MonomialVector(PolyVar[], [Int[]])
    @fact MonomialVector([x]) --> MonomialVector([x], [[1]])
end
module newmodule
    using FactCheck
    import MultivariatePolynomials
    facts("Polyvar macro hygiene") do
        # Verify that the @polyvar macro works when the package has been activated
        # with `import` instead of `using`.
        MultivariatePolynomials.@polyvar x y
        @fact isa(x, MultivariatePolynomials.PolyVar) --> true
        @fact isa(y, MultivariatePolynomials.PolyVar) --> true
    end
end
