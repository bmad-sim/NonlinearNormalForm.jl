# =================================================================================== #
# Helper functions
getdef(m::Union{TaylorMap,VectorField}) = TI.getdef(first(m.x))
nvars(m::Union{TaylorMap,VectorField}) = length(m.x0)
nparams(m::Union{TaylorMap,VectorField}) = length(m.x) - length(m.x0)
ndiffs(m::Union{TaylorMap,VectorField}) = length(m.x)
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))
# =================================================================================== #
# Array type promotion

function _promote_array_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:Array{T,N} where{T}, G<:Number}
  arrtype = Base.typename(A).wrapper
  neweltype = promote_type(eltype(A), G)
  return arrtype{neweltype, M}
end
# Default fallback:
function _promote_array_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{ndims(A)}) where {M,A,G<:Number}
  return _promote_array_eltype(Array{eltype(A),M}, G, Val{M})
end

function _real_array_eltype(::Type{A}, ::Type{Val{M}}=Val{N}) where {M,N, A<:Array{T,N} where{T}}
  arrtype = Base.typename(A).wrapper
  neweltype = real(eltype(A))
  return arrtype{neweltype, M}
end
# Default fallback:
function _real_array_eltype(::Type{A}, ::Type{Val{M}}=Val{ndims(A)})  where {M,A}
  return _real_array_eltype(Array{eltype(A),M}, Val{M})
end

#_vec_type(::Type{A}) where {T,N,A<:AbstractArray{T,N}} = (Base.typename(A).wrapper){T, 1}
#_mat_type(::Type{A}) where {T,N,A<:AbstractArray{T,N}} = (Base.typename(A).wrapper){T, 2}

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