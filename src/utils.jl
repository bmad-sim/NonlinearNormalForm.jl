# =================================================================================== #
# Helper functions
getdef(m::Union{TaylorMap,VectorField}) = TI.getdef(first(m.x))
nvars(m::Union{TaylorMap,VectorField}) = length(m.x0)
nparams(m::Union{TaylorMap,VectorField}) = length(m.x) - length(m.x0)
ndiffs(m::Union{TaylorMap,VectorField}) = length(m.x)
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))
# =================================================================================== #
# Similar eltype
# Returns the equivalent of container A to instead have eltype G and dims M
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:Array{T,N} where{T},G}
  arrtype = Base.typename(A).wrapper
  return arrtype{G, M}
end

# Default fallback:
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{ndims(A)}) where {M,A<:AbstractArray,G}
  return similar_eltype(Array{eltype(A),M}, G, Val{M})
end

# StaticArrays definition
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:StaticArray{N,T} where{T},G}
  return similar_type(A, G, Size(M))
end

# =================================================================================== #
# Coast check

function coastidx(m)
  return -1
  nv = numvars(m)
  for i in nv-1:nv # check only the last two planes
    if abs(m.x[i][0]) < NonlinearNormalForm.coast_threshold
      cycleidx = GTPSA.cycle!(m.x[i], 0, 0, C_NULL, C_NULL)
      if cycleidx == i && abs(m.x[i][i] - 1) < NonlinearNormalForm.coast_threshold
        return i
      end
    end
  end

  return -1
end