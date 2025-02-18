#=

Contains utility routines for getting TPSA info from a map/vector field,
getting the jacobian or transpose-jacobian (jacobiant) from a map/vector 
field array type promotion, checking if the last plane of a map is coasting.

=#
# =================================================================================== #
# Helper functions
getinit(m::Union{TaylorMap,VectorField}) = TI.getinit(first(m.x))
ndiffs(m::Union{TaylorMap,VectorField}) = length(m.x)
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))

# NNF-specific helpers:
nvars(m::Union{TaylorMap,VectorField}) = length(m.x0)
nparams(m::Union{TaylorMap,VectorField}) = length(m.x) - length(m.x0)
# =================================================================================== #
# Jacobian/jacobiant

function jacobian(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  nv = nvars(m)
  M = similar(X0, (nv,nv))
  for i in 1:nv^2
    M[i] = TI.geti(m.x[mod1(i,nv)], floor(Int, (i-1)/nv)+1)
  end
  return M
end

function jacobian(m::TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S}
  nv = nvars(m)
  return StaticArrays.sacollect(SMatrix{nv,nv,eltype(X0)}, TI.geti(m.x[mod1(i,nv)], floor(Int, (i-1)/nv)+1) for i in 1:nv^2)
end

function jacobiant(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  nv = nvars(m)
  M = similar(X0, (nv,nv))
  for i in 1:nv^2
    M[i] = TI.geti(m.x[floor(Int, (i-1)/nv)+1], mod1(i,nv))
  end
  return M
end

function jacobiant(m::TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S}
  nv = nvars(m)
  return StaticArrays.sacollect(SMatrix{nv,nv,eltype(X0)}, TI.geti(m.x[floor(Int, (i-1)/nv)+1], mod1(i,nv)) for i in 1:nv^2)
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
  return Array{G, M}
end

# Default fallback:
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{ndims(A)}) where {M,A<:AbstractArray,G}
  return similar_eltype(Array{eltype(A),M}, G, Val{M})
end

# StaticArrays definition
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:StaticArray{N,T} where{T},G}
  return similar_type(A, G, Size(M))
end

# Quaternion
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{1}) where {M,A<:Quaternion,G}
  M == 1 || error("Quaternion type can only have dimension 1")
  return Quaternion{promote_type(eltype(A),G)}
end

# =================================================================================== #
# Coast check

function coastidx(m)
  nv = nvars(m)
  for i in nv-1:nv # check only the last two planes
    if abs(TI.geti(m.x[i], 0)) < COAST # if scalar part is 0
      cycleidx = TI.cycle!(m.x[i], 0)
      if cycleidx == i && abs(TI.geti(m.x[i], i) - 1) < COAST
        return i
      end
    end
  end

  return -1
end

# =================================================================================== #
# Context-dependent skew symmetric matrix S
"""
Generic symplectic skew symmetric S matrix (size inferred from 
other matrix in operations) using SkewLinearAlgebra's `JMatrix`
"""
struct SymplecticS end

"""
Generic symplectic skew symmetric S matrix (size inferred from 
other matrix in operations) using SkewLinearAlgebra's `JMatrix`
"""
const S = SymplecticS()

(S::SymplecticS)(n::Integer) = JMatrix{Int8,+1}(n)

for op = (:+, :-, :*, :/)
@eval begin
Base.$op(S::SymplecticS,M) = Base.$op(JMatrix{Int8,+1}(size(M,1)), M)
Base.$op(M,S::SymplecticS) = Base.$op(M, JMatrix{Int8,+1}(size(M,2)))
end
end

"""
    checksymp(M)

Returns `tranpose(M)*S*M - S`, where `S` is the skew-symmetric matrix 
`S = blkdiag([0 1; -1 0], ...)`. If `M` is symplectic, then the result should be a matrix 
containing all zeros. The non-symplectic parts of the matrix can be identified 
by those nonzero elements in the result.
"""
function checksymp(M)
  res = transpose(M)*S*M-S
  return res
end
# =================================================================================== #