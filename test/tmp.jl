function mapexponents_test()
    @polyvar x y
    a = x^2
    b = x * y
    MultivariatePolynomials.mapexponents!(+, a, b)
    dump(a)
    @test variables(a) == variables(b)
end
mapexponents_test()
