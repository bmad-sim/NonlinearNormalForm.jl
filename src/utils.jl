#=

Contains utility routines for getting TPSA info from a map/vector field,
getting the jacobian or transpose-jacobian (jacobiant) from a map/vector 
field array type promotion, checking if the last plane of a map is coasting.

=#
# =================================================================================== #
# Helper functions
getdef(m::Union{TaylorMap,VectorField}) = TI.getdef(first(m.x))
nvars(m::Union{TaylorMap,VectorField}) = length(m.x0)
nparams(m::Union{TaylorMap,VectorField}) = length(m.x) - length(m.x0)
ndiffs(m::Union{TaylorMap,VectorField}) = length(m.x)
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))
# =================================================================================== #
# Jacobian/jacobiant

function jacobian(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  nv = nvars(m)
  M = similar(X0, (nv,nv))
  for i in 1:nv^2
    M[i] = TI.geti(m.x[mod1(i,nv)], floor(Int, (i-1)/nv))
  end
  return M
end

function jacobian(m::TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S}
  nv = nvars(m)
  return StaticArrays.sacollect(SMatrix{nv,nv,eltype(X0)}, TI.geti(m.x[mod1(i,nv)], floor(Int, (i-1)/nv)) for i in 1:nv^2)
end

# =================================================================================== #
# Get/set scalar part of orbital ray
function getscalar(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  nv = nvars(m)
  return map(t->TI.geti(t, 0), view(m.x, 1:nv))
end

function getscalar(m::TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S}
  nv = nvars(m)
  return StaticArrays.sacollect(SVector{nv,eltype(X0)}, TI.geti(m.x[i], 0) for i in 1:nv)
end

function setscalar!(
  m::TaylorMap{X0,X,Q,S}, 
  xs::Number; 
  scl0::Union{Nothing,Number}=nothing,
  scl1::Number=1
) where {X0,X,Q,S}
  nv = nvars(m)
  if isnothing(scl0)
    for i in 1:nv
      TI.seti!(m.x[i], scl1*xs, 0)
    end
  else
    for i in 1:nv
      TI.seti!(m.x[i], TI.geti(m.x[i], 0)*scl0 + scl1*xs, 0)
    end
  end
end

function setscalar!(
  m::TaylorMap{X0,X,Q,S}, 
  xs::AbstractArray; 
  scl0::Union{Nothing,Number}=nothing,
  scl1::Number=1
) where {X0,X,Q,S}

  nv = nvars(m)
  if isnothing(scl0)
    for i in 1:nv
      TI.seti!(m.x[i], xs[i], 0)
    end
  else
    for i in 1:nv
      TI.seti!(m.x[i], TI.geti(m.x[i], 0)*scl0 + scl1*xs[i], 0)
    end
  end
end


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

# =================================================================================== #
