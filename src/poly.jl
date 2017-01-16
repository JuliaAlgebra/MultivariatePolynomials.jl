export Term, VecPolynomial, MatPolynomial, SOSDecomposition, TermType
export NCTerm, NCVecPolynomial, NCMatPolynomial, NCSOSDecomposition, NCTermType

for (TT, TC, T, P, MP, SOSD, PV, M, MV) in [(:TermType, :TermContainer, :Term, :VecPolynomial, :MatPolynomial, :SOSDecomposition, :PolyVar, :Monomial, :MonomialVector), (:NCTermType, :NCTermContainer, :NCTerm, :NCVecPolynomial, :NCMatPolynomial, :NCSOSDecomposition, :NCPolyVar, :NCMonomial, :NCMonomialVector)]
    @eval begin
        abstract $TT{T} <: PolyType
        zero{T}(t::$TT{T}) = $P(T[], $MV(vars(t), Vector{Vector{Int}}()))
        zero{T<:$TT}(::Type{T}) = $P(eltype(T)[], $MV())
        one{T}(t::$TT{T}) = $P([one(T)], $MV(vars(t), [zeros(Int, length(vars(t)))]))
        one{T<:$TT}(::Type{T}) = $P([one(eltype(T))], $MV($PV[], [Int[]]))

        abstract $TC{T} <: $TT{T}
        eltype{T}(::Type{$TC{T}}) = T
        eltype{T}(::Type{$TT{T}}) = T

        type $T{T} <: $TC{T}
            α::T
            x::$M
        end
        $T(t::$T) = t
        $T(x::$M) = $T{Int}(x)
        $T(x::$PV) = $T($M(x))
        $TC{T<:Union{$M,$PV}}(x::T) = $T(x)
        $TC{T}(t::$TC{T}) = t
        $T{T}(α::T) = $T{T}(α, $M())
        Base.convert{T}(::Type{$T{T}}, t::$T{T}) = t
        Base.convert{T}(::Type{$T{T}}, t::$T) = $T{T}(T(t.α), t.x)
        Base.convert{T}(::Type{$T{T}}, x::$M) = $T{T}(one(T), x)
        Base.convert{T}(::Type{$T{T}}, x::$PV) = $T{T}($M(x))
        Base.convert{T}(::Type{$T{T}}, α) = $T(T(α))

        Base.convert{T}(::Type{$TC{T}}, x::Union{$M,$PV}) = $T{T}(x)
        Base.convert{T}(::Type{$TC{T}}, α::T) = $T{T}(α, $M())
        Base.convert{S,T}(::Type{$TC{T}}, α::S) = $TC{T}(T(α))
        $TC{T}(α::T) = $TC{T}(α)

        Base.convert{T}(::Type{$TC{T}}, t::$T) = $T{T}(t)
        Base.convert{T<:$TC}(::Type{T}, t::$T) = $T{eltype(T)}(t)

        zero{T}(t::$T{T}) = $T(zero(T), t.x)
        zero{T}(::Type{$T{T}}) = $T(zero(T), $M())
        one{T}(t::$T{T}) = $T(one(T), $M(t.x.vars, zeros(Int, length(t.x.vars))))
        one{T}(::Type{$T{T}}) = $T(one(T), $M())

        # Invariant:
        # a and x might be empty: meaning it is the zero polynomial
        # a does not contain any zeros
        # x is increasing in the monomial order (i.e. grlex)
        type $P{T} <: $TC{T}
            a::Vector{T}
            x::$MV

            function $P(a::Vector{T}, x::$MV)
                if length(a) != length(x)
                    error("There should be as many coefficient than monomials")
                end
                zeroidx = Int[]
                for (i,α) in enumerate(a)
                    if iszero(α)
                        push!(zeroidx, i)
                    end
                end
                if !isempty(zeroidx)
                    isnz = ones(Bool, length(a))
                    isnz[zeroidx] = false
                    nzidx = find(isnz)
                    a = a[nzidx]
                    x = x[nzidx]
                end
                new(a, x)
            end
        end
        (::Type{$P{T}}){T}(a::Vector{T}, x::Vector) = $P{T}(a, $MV(x))
        (::Type{$P{T}}){S,T}(a::Vector{S}, x::Vector) = $P{T}(Vector{T}(a), $MV(x))

        Base.copy{T}(p::$P{T}) = $P{T}(copy(p.a), copy(p.x))

        $P{T}(a::Vector{T}, x::$MV) = $P{T}(a, x)
        $P(a::Vector, x::Vector) = $P(a, $MV(x))

        $P(x) = $P($T(x))
        $P{T}(p::$P{T}) = p
        $P{T}(t::$T{T}) = $P{T}([t.α], [t.x])
        Base.convert{T}(::Type{$P{T}}, x) = $P($T{T}(x))
        Base.convert{T}(::Type{$P{T}}, t::$T) = $P{T}([T(t.α)], [t.x])
        Base.convert{T}(::Type{$P{T}}, p::$P{T}) = p
        Base.convert{S,T}(::Type{$P{T}}, p::$P{S}) = $P(Vector{T}(p.a), p.x)

        Base.convert{T}(::Type{$TC{T}}, p::$P) = $P{T}(p)

        function (::Type{$P{T}}){T}(f::Function, x::$MV)
            a = T[f(i) for i in 1:length(x)]
            $P{T}(a, x)
        end
        (::Type{$P{T}}){T}(f::Function, x::Vector) = $P{T}(f, $MV(x))

        function vecpolynomialclean{T}(vars::Vector{$PV}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
            a, Z = removedups(adup, Zdup)
            $P(a, $MV(vars, Z))
        end

        eltype{T}(::Type{$P{T}}) = T
        getindex(p::$P, I::Int) = $T(p.a[I[1]], p.x[I[1]])
        removemonomials(p::$P, x::Vector) = removemonomials(p, $MV(x))

        type $MP{T} <: $TT{T}
            Q::Vector{T}
            x::$MV
        end

        function (::Type{$MP{T}}){T}(f::Function, x::$MV)
            $MP{T}(trimat(T, f, length(x)), x)
        end
        (::Type{$MP{T}}){T}(f::Function, x::Vector) = $MP{T}(f, $MV(x))

        function $MP{T}(Q::Matrix{T}, x::$MV)
            $MP{T}((i,j) -> Q[i,j], x)
        end
        $MP(Q::Matrix, x::Vector) = $MP(Q, $MV(x))

        function $P{T}(p::$MP{T})
            if isempty(p.Q)
                zero($P{T})
            else
                n = length(p.x)
                N = trimap(n, n, n)
                Z = Vector{Vector{Int}}(N)
                U = typeof(2*p.Q[1] + p.Q[1])
                a = Vector{U}(N)
                for i in 1:n
                    for j in i:n
                        k = trimap(i, j, n)
                        Z[k] = p.x.Z[i] + p.x.Z[j]
                        if i == j
                            a[k] = p.Q[k]
                        else
                            a[k] = 2*p.Q[k]
                        end
                    end
                end
                vecpolynomialclean(p.x.vars, a, Z)
            end
        end
        $TC(p::$MP) = $P(p)

        type $SOSD{T} <: $TT{T}
            ps::Vector{$P{T}}
            function $SOSD(ps::Vector{$P{T}})
                new(ps)
            end
        end
        function (::Type{$SOSD{T}}){T}(ps::Vector)
            $SOSD(Vector{$P{T}}(ps))
        end
        function $SOSD(ps::Vector)
            T = reduce(promote_type, Int, map(eltype, ps))
            $SOSD{T}(ps)
        end

        function $SOSD{T}(p::$MP{T})
            n = length(p.x)
            # TODO LDL^T factorization for SDP is missing in Julia
            # it would be nice to have though
            A = getmat(p)
            Q = chol(A)
            ps = [$P(Q[i,:], p.x) for i in 1:n]
            $SOSD(ps)
        end
    end
end

zero(::Type{PolyType}) = zero(VecPolynomial{Int})
one(::Type{PolyType}) = one(VecPolynomial{Int})
zero(p::PolyType) = zero(PolyType)
one(p::PolyType) = one(PolyType)

zero(::Type{NCPolyVar}) = zero(NCVecPolynomial{Int})
zero(::Type{NCMonomial}) = zero(NCVecPolynomial{Int})
one(::Type{NCPolyVar}) = one(NCVecPolynomial{Int})
one(::Type{NCMonomial}) = one(NCVecPolynomial{Int})
zero{PVM <: Union{NCPolyVar, NCMonomial}}(::PVM) = zero(PVM)
one{PVM <: Union{NCPolyVar, NCMonomial}}(::PVM) = one(PVM)

export getmat, monomials, removemonomials
typealias AbstractTermContainer{T} Union{TermContainer{T}, NCTermContainer{T}}
typealias AbstractTermType{T} Union{TermType{T}, NCTermType{T}}
typealias AbstractTerm{T} Union{Term{T}, NCTerm{T}}

eltype{T}(p::AbstractTermType{T}) = T

Base.convert(::Type{Any}, t::AbstractTerm) = t
Base.copy{T<:AbstractTerm}(t::T) = T(copy(t.α), copy(t.x))

vars(t::AbstractTerm) = vars(t.x)

eltype{T}(::Type{AbstractTerm{T}}) = T
length(::AbstractTerm) = 1
isempty(::AbstractTerm) = false
start(::AbstractTerm) = false
done(::AbstractTerm, state) = state
next(t::AbstractTerm, state) = (t, true)
getindex(t::AbstractTerm, I::Int) = t

typealias AbstractVecPolynomial{T} Union{VecPolynomial{T}, NCVecPolynomial{T}}
Base.convert(::Type{Any}, p::AbstractVecPolynomial) = p

Base.convert(::Type{PolyType}, p::AbstractTermContainer) = p
function Base.convert{S}(::Type{S}, p::AbstractTermContainer)
    s = zero(S)
    for t in p
        if sum(abs(t.x.z)) > 0
            # The polynomial is not constant
            throw(InexactError())
        end
        s += S(t.α)
    end
    s
end


function removedups{T}(adup::Vector{T}, Zdup::Vector{Vector{Int}})
    σ = sortperm(Zdup, rev=true)
    Z = Vector{Vector{Int}}()
    a = Vector{T}()
    i = 0
    j = 1
    while j <= length(adup)
        k = σ[j]
        if j == 1 || Zdup[k] != Zdup[σ[j-1]]
            push!(Z, Zdup[k])
            push!(a, adup[k])
            i += 1
        else
            a[i] += adup[k]
        end
        j += 1
    end
    a, Z
end

vars(p::AbstractVecPolynomial) = vars(p.x)

length(p::AbstractVecPolynomial) = length(p.a)
isempty(p::AbstractVecPolynomial) = length(p) > 0
start(::AbstractVecPolynomial) = 1
done(p::AbstractVecPolynomial, state) = length(p) < state
next(p::AbstractVecPolynomial, state) = (p[state], state+1)

extdeg(p::AbstractVecPolynomial) = extdeg(p.x)
mindeg(p::AbstractVecPolynomial) = mindeg(p.x)
maxdeg(p::AbstractVecPolynomial) = maxdeg(p.x)

monomials(p::AbstractVecPolynomial) = p.x
function removemonomials(p::AbstractVecPolynomial, x::AbstractMonomialVector)
    # use the fact that monomials are sorted to do this O(n) instead of O(n^2)
    j = 1
    I = Int[]
    for (i,t) in enumerate(p)
        while j <= length(x) && x[j] < t.x
            j += 1
        end
        if j > length(x) || x[j] != t.x
            push!(I, i)
        end
    end
    P(p.a[I], p.x[I])
end

typealias AbstractMatPolynomial{T} Union{MatPolynomial{T}, NCMatPolynomial{T}}

function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function trimat{T}(::Type{T}, f, n)
    Q = Vector{T}(trimap(n, n, n))
    for i in 1:n
        for j in i:n
            Q[trimap(i, j, n)] = f(i, j)
        end
    end
    Q
end

function getindex(p::AbstractMatPolynomial, I::NTuple{2,Int})
    i, j = I
    if i < j
        i, j = (j, i)
    end
    n = length(p.x)
    p.Q[trimap(i,j,n)]
end

function getmat{T}(p::AbstractMatPolynomial{T})
    n = length(p.x)
    A = Matrix{T}(n, n)
    for i in 1:n, j in i:n
        A[j,i] = A[i,j] = p.Q[trimap(i,j,n)]
    end
    A
end

typealias AbstractSOSDecomposition{T} Union{SOSDecomposition{T}, NCSOSDecomposition{T}}

length(p::AbstractSOSDecomposition) = length(p.ps)
isempty(p::AbstractSOSDecomposition) = isempty(p.ps)
start(p::AbstractSOSDecomposition) = start(p.ps)
done(p::AbstractSOSDecomposition, state) = done(p.ps, state)
next(p::AbstractSOSDecomposition, state) = next(p.ps, state)
