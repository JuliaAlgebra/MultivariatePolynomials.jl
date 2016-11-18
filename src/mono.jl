export PolyVar, Monomial, MonomialVector, @polyvar, monomials, polyvecvar, vars

import Base.copy

abstract PolyType
abstract MonomialContainer <: PolyType

function polyvecvar(prefix, idxset)
  [PolyVar("$(prefix * string(i))") for i in idxset]
end

function buildpolyvar(var)
  if isa(var, Symbol)
    esc(:($var = PolyVar($"$var")))
  else
    isa(var, Expr) || error("Expected $var to be a variable name")
    Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
    length(var.args) == 2 || error("Expected $var to have one index set")
    #tmp = gensym()
    varname = var.args[1]
    prefix = string(var.args[1])
    idxset = var.args[2]
    esc(:($varname = polyvecvar($prefix, $idxset)))
  end
end

# Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
macro polyvar(args...)
  reduce((x,y) -> :($x; $y), :(), [buildpolyvar(arg) for arg in args])
end

# TODO equality should be between name ?
immutable PolyVar <: MonomialContainer
  name::AbstractString
end
copy(x::PolyVar) = x

vars(x::PolyVar) = [x]
nvars(::PolyVar) = 1

length(::PolyVar) = 1
isempty(::PolyVar) = false
start(::PolyVar) = false
done(::PolyVar, state) = state
next(x::PolyVar, state) = (Monomial(x), true)

# Invariant:
# vars is increasing
# z may contain 0's (otherwise, getindex of MonomialVector would be inefficient)
type Monomial <: MonomialContainer
  vars::Vector{PolyVar}
  z::Vector{Int}

  function Monomial(vars::Vector{PolyVar}, z::Vector{Int})
    if length(vars) != length(z)
      error("There should be as many vars than exponents")
    end
    new(vars, z)
  end
end
Monomial() = Monomial(PolyVar[], Int[])
deg(x::Monomial) = sum(x.z)

Base.convert(::Type{Monomial}, x::PolyVar) = Monomial([x], [1])

copy(m::Monomial) = Monomial(copy(m.vars), copy(m.z))

length(::Monomial) = 1
isempty(::Monomial) = false
start(::Monomial) = false
done(::Monomial, state) = state
next(x::Monomial, state) = (x, true)

# Invariant: Always sorted and no zero vector
type MonomialVector <: MonomialContainer
  vars::Vector{PolyVar}
  Z::Vector{Vector{Int}}

  function MonomialVector(vars::Vector{PolyVar}, Z::Vector{Vector{Int}})
    for z in Z
      if length(vars) != length(z)
        error("There should be as many vars than exponents")
      end
    end
    new(vars, Z)
  end
end

copy(m::MonomialVector) = MonomialVector(copy(m.vars), copy(m.Z))

function getindex(x::MonomialVector, I)
  MonomialVector(x.vars, x.Z[I])
end
function getindex(x::MonomialVector, i::Integer)
  Monomial(x.vars, x.Z[i])
end
length(x::MonomialVector) = length(x.Z)
isempty(x::MonomialVector) = length(x) == 0
start(::MonomialVector) = 1
done(x::MonomialVector, state) = length(x) < state
next(x::MonomialVector, state) = (Monomial(x.vars, x.Z[state]), state+1)

vars{T<:Union{Monomial,MonomialVector}}(x::T) = x.vars

nvars(x::MonomialContainer) = length(vars(x))

function getZfordegs(n, degs, filter::Function)
  Z = Vector{Vector{Int}}()
  for deg in degs
    z = zeros(Int, n)
    z[1] = deg
    while true
      if filter(z)
        push!(Z, z)
        z = copy(z)
      end
      if z[end] == deg
        break
      end
      sum = 1
      for j in (n-1):-1:1
        if z[j] != 0
          z[j] -= 1
          z[j+1] += sum
          break
        else
          sum += z[j+1]
          z[j+1] = 0
        end
      end
    end
  end
  Z
end

function MonomialVector(vars::Vector{PolyVar}, degs, filter::Function = x->true)
  MonomialVector(vars, getZfordegs(length(vars), degs, filter))
end
MonomialVector(vars::Vector{PolyVar}, degs::Int, filter::Function = x->true) = MonomialVector(vars, [degs], filter)
function monomials(vars::Vector{PolyVar}, degs, filter::Function = x->true)
  Z = getZfordegs(length(vars), degs, filter)
  [Monomial(vars, z) for z in Z]
end
monomials(vars::Vector{PolyVar}, degs::Int, filter::Function = x->true) = monomials(vars, [degs], filter)


function myunion(varsvec::Vector{Vector{PolyVar}})
  n = length(varsvec)
  is = ones(Int, n)
  maps = [ zeros(Int, length(vars)) for vars in varsvec ]
  nonempty = IntSet(find([!isempty(vars) for vars in varsvec]))
  vars = Vector{PolyVar}()
  while !isempty(nonempty)
    imin = 0
    for i in nonempty
      if imin == 0 || varsvec[i][is[i]] < varsvec[imin][is[imin]]
        imin = i
      end
    end
    var = varsvec[imin][is[imin]]
    push!(vars, var)
    for i in nonempty
      if var == varsvec[i][is[i]]
        maps[i][is[i]] = length(vars)
        if is[i] == length(varsvec[i])
          pop!(nonempty, i)
        else
          is[i] += 1
        end
      end
    end
  end
  vars, maps
end

function MonomialVector(X::Vector)
  varsvec = Vector{PolyVar}[ vars(x) for x in X ]
  allvars, maps = myunion(varsvec)
  nvars = length(allvars)
  n = sum([length(x) for x in X])
  Z = [zeros(Int, nvars) for i in 1:n]
  offset = 0
  for (i, x) in enumerate(X)
    for (j, m) in enumerate(x)
      Z[offset+j][maps[i]] = m.z
    end
    offset += length(x)
  end
  MonomialVector(allvars, Z)
end
