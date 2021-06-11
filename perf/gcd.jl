using DynamicPolynomials

function bench()
    o = one(Rational{BigInt})
    @polyvar x y z
    a = (o * x + o * y^2) * (o * z^3 + o * y^2 + o * x)
    @show a
    b = (o * x + o * y + o * z) * (o * x^2 + o * y)
    @show b
    @time gcd(a, b)
end
