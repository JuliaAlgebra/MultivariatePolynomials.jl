module MultivariatePolynomialsJuMPExt

import JuMP
import MultivariatePolynomials as MP
import StarAlgebras as SA

MP.@scalar_arithmetic JuMP.AbstractJuMPScalar

function JuMP.function_string(
    mime::MIME"text/plain",
    p::MP.AbstractMonomialLike,
)
    return sprint(show, mime, p)
end

function JuMP.function_string(
    mime::MIME"text/latex",
    p::MP.AbstractMonomialLike,
)
    return SA.trim_LaTeX(mime, sprint(show, mime, p))
end

end
