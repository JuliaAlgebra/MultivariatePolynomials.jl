export Measure, NCMeasure, zeta, ζ

abstract AbstractMeasure

for (MS, MV, PV) in [(:Measure, :MonomialVector, :PolyVar), (:NCMeasure, :NCMonomialVector, :NCPolyVar)]
    @eval begin
        type $MS{T} <: AbstractMeasure
          a::Vector{T}
          x::$MV

          function $MS(a::Vector{T}, x::$MV)
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

        $MS{T}(a::Vector{T}, x::$MV) = $MS{T}(a, x)
        function $MS(a::Vector, x::Vector)
            if length(a) != length(x)
                error("There should be as many coefficient than monomials")
            end
            perm, X = sortmonovec($PV, x)
            $MS(a[perm], X)
        end

        function ζ{T}(v::Vector{T}, x::$MV, varorder::Vector{$PV})
            $MS(T[m(v, varorder) for m in x], x)
        end
    end
end
