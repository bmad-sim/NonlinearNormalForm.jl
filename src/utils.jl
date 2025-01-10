# =================================================================================== #
# GTPSA helper functions

maxord(m::Union{TaylorMap,VectorField}) = unsafe_load(getdesc(m).desc).mo
prmord(m::Union{TaylorMap,VectorField}) = unsafe_load(getdesc(m).desc).po
vpords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numnn(m))
vords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numvars(m))
pords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numparams(m))

# GTPSA provides these functions for only pure TPS/ComplexTPSs and Descriptor
getdesc(m::Union{TaylorMap{S,T,U,V},VectorField{T,U}}) where {S,T,U,V} = getdesc(first(m.x))
numvars(m::Union{TaylorMap{S,T,U,V},VectorField{T,U}}) where {S,T,U,V} = numvars(first(m.x))
numparams(m::Union{TaylorMap{S,T,U,V},VectorField{T,U}}) where {S,T,U,V} = numparams(first(m.x))
numnn(m::Union{TaylorMap{S,T,U,V},VectorField{T,U}}) where {S,T,U,V} = numnn(first(m.x))

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