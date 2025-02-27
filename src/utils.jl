#=

Contains utility routines for getting TPSA info from a map/vector field,
getting the jacobian or transpose-jacobian (jacobiant) from a map/vector 
field array type promotion, checking if the last plane of a map is coasting.

=#
# =================================================================================== #
# Helper functions
getinit(m::Union{TaylorMap,VectorField}) = TI.getinit(first(m.x))
ndiffs(m::TaylorMap) = length(m.x)
ndiffs(F::VectorField) = TI.ndiffs(first(F.x))
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))

# NNF-specific helpers:
nvars(m::TaylorMap) = length(m.x0)
nvars(F::VectorField) = length(F.x)
nparams(m::Union{TaylorMap,VectorField}) = ndiffs(m) - nvars(m)
nhvars(m::Union{TaylorMap,VectorField}) = iseven(nvars(m)) ? nvars(m) : nvars(m)-1 # number of "harmonic" variables

iscoasting(m::Union{TaylorMap,VectorField}) = !iseven(nvars(m))
# =================================================================================== #
# Jacobian/jacobiant
# Useful options should be 1) harmonic variables, 2) variables, and 3) variables + parameters
abstract type OptionType end
struct HarmonicVariables <: OptionType end  # e.g. 4x4 matrix
struct Variables <: OptionType end # e.g. 4x4 or 5x5 (coasting beam) matrix
struct HarmonicVariablesAndParameters <: OptionType end # e.g. 4 x np matrix
struct VariablesAndParameters <: OptionType end # eg. 5 x np (coasting beam) matrix
struct Parameters <: OptionType end # e.g. 5 x np matrix
struct All <: OptionType end

const HVARS = HarmonicVariables()
const VARS = Variables()
const HVARS_PARAMS = HarmonicVariablesAndParameters()
const VARS_PARAMS = VariablesAndParameters()
const PARAMS = Parameters()
const ALL = All()

@inline function getjacsize(m::TaylorMap, ::T) where {T<:OptionType}
  if T == HarmonicVariables
    nrows = nhvars(m)
    ncols = nhvars(m)
  elseif T == Variables
    nrows = nvars(m)
    ncols = nvars(m)
  elseif T == HarmonicVariablesAndParameters
    nrows = nhvars(m)
    ncols = ndiffs(m)
  elseif T == VariablesAndParameters
    nrows = nvars(m)
    ncols = ndiffs(m)
  elseif T == Parameters
    nrows = nvars(m)
    ncols = nparams(m)
  elseif T == All
    nrows = ndiffs(m)
    ncols = ndiffs(m)
  else
    error("Invalid option type")
  end
  return nrows, ncols
end

function jacobian(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X,Q,S,T<:OptionType}
  nrows, ncols = getjacsize(m, which)
  col_start = T == Parameters ? nvars(m)+1 : 1
  M = similar(X0, (nrows,ncols))
  for col in col_start:col_start+ncols-1
    for row in 1:nrows
      M[row,col-col_start+1] = TI.geti(m.x[row], col)
    end
  end
  return M
end

function jacobian(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X<:StaticArray,Q,S,T<:OptionType}
  nrows, ncols = getjacsize(m, which)
  col_start = T == Parameters ? nvars(m)+1 : 1
  return StaticArrays.sacollect(SMatrix{nrows,ncols,eltype(X0)}, 
  TI.geti(m.x[row], col) for col in col_start:col_start+ncols-1 for row in 1:nrows)
end

function jacobiant(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X,Q,S,T<:OptionType}
  nrowst, ncolst = getjacsize(m, which)
  colt_start = T == Parameters ? nvars(m)+1 : 1
  M = similar(X0, (ncolst,nrowst))
  for colt in colt_start:colt_start+ncolst-1
    for rowt in 1:nrowst
      M[colt-colt_start+1,rowt] = TI.geti(m.x[rowt], colt)
    end
  end
  return M
end

function jacobiant(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X<:StaticArray,Q,S,T<:OptionType}
  return transpose(jacobian(m, which))
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