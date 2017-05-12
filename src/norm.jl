Base.norm(p::Polynomial, r::Int) = norm(p.a, r)

Base.norm(p::PolyVar, r)       = norm(Polynomial(p), r)
Base.norm(p::Term, r)          = norm(Polynomial(p), r)
Base.norm(p::MatPolynomial, r) = norm(Polynomial(p), r)

Base.norm(p::Polynomial)    = norm(p, 2)
Base.norm(p::PolyVar)       = norm(p, 2)
Base.norm(p::Term)          = norm(p, 2)
Base.norm(p::MatPolynomial) = norm(p, 2)
