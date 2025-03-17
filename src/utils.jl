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
struct CoastVariables <: OptionType end # e.g. 1x5 matrix
struct Variables <: OptionType end # e.g. 4x4 or 5x5 (coasting beam) matrix
struct HarmonicVariablesAndParameters <: OptionType end # e.g. 4 x (nv+np) matrix
struct CoastVariablesAndParameters <: OptionType end # e.g. 1 x (nv+np) matrix
struct VariablesAndParameters <: OptionType end # eg. 5 x np (coasting beam) matrix
struct Parameters <: OptionType end # e.g. 5 x np matrix
struct HarmonicParameters <: OptionType end # e.g. 4 x np matrix
struct CoastParameters <: OptionType end # e.g. 1 x np matrix
struct All <: OptionType end

const HVARS = HarmonicVariables()
const CVARS = CoastVariables()
const VARS = Variables()
const HVARS_PARAMS = HarmonicVariablesAndParameters()
const CVARS_PARAMS = CoastVariablesAndParameters()
const VARS_PARAMS = VariablesAndParameters()
const PARAMS = Parameters()
const HPARAMS = HarmonicParameters()
const CPARAMS = CoastParameters()
const ALL = All()

@inline function getjacinfo(m::TaylorMap, ::T) where {T<:OptionType}
  if T == HarmonicVariables
    nrows = nhvars(m)
    ncols = nhvars(m)
    row_start = 1
    col_start = 1
  elseif T == CoastVariables
    nrows = 1
    ncols = nhvars(m)
    row_start = nvars(m)
    col_start = 1
  elseif T == Variables
    nrows = nvars(m)
    ncols = nvars(m)
    row_start = 1
    col_start = 1
  elseif T == HarmonicVariablesAndParameters
    nrows = nhvars(m)
    ncols = ndiffs(m)
    row_start = 1
    col_start = 1
  elseif T == CoastVariablesAndParameters
    nrows = 1
    ncols = ndiffs(m)
    row_start = nvars(m)
    col_start = 1
  elseif T == VariablesAndParameters
    nrows = nvars(m)
    ncols = ndiffs(m)
    row_start = 1
    col_start = 1
  elseif T == Parameters
    nrows = nvars(m)
    ncols = nparams(m)
    row_start = 1
    col_start = nvars(m)+1
  elseif T == HarmonicParameters
    nrows = nhvars(m)
    ncols = nparams(m)
    row_start = 1
    col_start = nvars(m)+1
  elseif T == CoastParameters
    nrows = 1
    ncols = nparams(m)
    row_start = nvars(m)
    col_start = nvars(m)+1
  elseif T == All
    nrows = ndiffs(m)
    ncols = ndiffs(m)
    row_start = 1
    col_start = 1
  else
    error("Invalid option type")
  end
  return nrows, ncols, row_start, col_start
end

function jacobian(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X,Q,S,T<:OptionType}
  nrows, ncols, row_start, col_start = getjacinfo(m, which)
  M = similar(m.x0, nrows, ncols)
  for col in col_start:col_start+ncols-1
    for row in row_start:row_start+nrows-1
      M[row-row_start+1,col-col_start+1] = TI.geti(m.x[row], col)
    end
  end
  return M
end

function jacobian(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X<:StaticArray,Q,S,T<:OptionType}
  nrows, ncols, row_start, col_start = getjacinfo(m, which)
  return StaticArrays.sacollect(SMatrix{nrows,ncols,eltype(X0)}, 
  TI.geti(m.x[row], col) for col in col_start:col_start+ncols-1 for row in row_start:row_start+nrows-1)
end

function jacobiant(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X,Q,S,T<:OptionType}
  nrowst, ncolst, rowt_start, colt_start = getjacinfo(m, which)
  M = similar(m.x0, ncolst,nrowst)
  for colt in colt_start:colt_start+ncolst-1
    for rowt in rowt_start:rowt_start+nrowst-1
      M[colt-colt_start+1,rowt-rowt_start+1] = TI.geti(m.x[rowt], colt)
    end
  end
  return M
end

function jacobiant(m::TaylorMap{X0,X,Q,S}, which::T=VARS) where {X0,X<:StaticArray,Q,S,T<:OptionType}
  return transpose(jacobian(m, which))
end

# =================================================================================== #
# Phasor basis transformations
# inv(c) is from_phasor
function ci_jacobian(m::TaylorMap{X0,X,Q,S}, ::T=VARS) where {X0,X,Q,S,T<:Union{HarmonicVariables,Variables}}
  nv = nvars(m)
  nhv = nhvars(m)
  n = T == HarmonicVariables ? nhv : nv
  CI = similar(m.x0, complex(eltype(m.x0)), n, n)
  CI .= 0
  for i in 1:Int(nhv/2)
    # x_new = 1/sqrt(2)*(x+im*p)
    CI[2*i-1,2*i-1] = 1/sqrt(2)
    CI[2*i-1,2*i] = complex(0,1/sqrt(2))

    # p_new = 1/sqrt(2)*(x-im*p)
    CI[2*i,2*i-1] = 1/sqrt(2)
    CI[2*i,2*i] = complex(0,-1/sqrt(2))
  end

  if T == Variables && iscoasting(m)
    CI[nv,nv] = 1
  end
  return CI
end

function ci_jacobian(m::TaylorMap{X0,X,Q,S}, ::T=VARS) where {X0<:StaticArray,X,Q,S,T<:Union{HarmonicVariables,Variables}}
  if T == Variables && iscoasting(m)
    nv = nvars(m)
    return StaticArrays.sacollect(SMatrix{nv,nv,eltype(X0)}, 
    begin
      if fld1(col,2) != fld1(row,2) 
        0
      elseif row == nv && col == nv
        1
      else # then we are in the block
        if mod1(col,2) == 1 # First column of ci is just 1/sqrt(2)
            1/sqrt(2)
        else # second column of ci is either im/sqrt(2) or -im/sqrt(2)
          if mod1(row,2) == 1
            complex(0,1/sqrt(2))
          else
            complex(0,-1/sqrt(2))
          end
        end
      end
    end for col in 1:nv for row in 1:nv)
  else
    nhv = nhvars(m)
    return StaticArrays.sacollect(SMatrix{nhv,nhv,eltype(X0)}, 
    begin
      if fld1(col,2) != fld1(row,2) 
        0
      else # then we are in the block
        if mod1(col,2) == 1 # First column of ci is just 1/sqrt(2)
            1/sqrt(2)
        else # second column of ci is either im/sqrt(2) or -im/sqrt(2)
          if mod1(row,2) == 1
            complex(0,1/sqrt(2))
          else
            complex(0,-1/sqrt(2))
          end
        end
      end
    end for col in 1:nhv for row in 1:nhv)
  end
end

# c is to_phasor
function c_jacobian(m::TaylorMap{X0,X,Q,S}, ::T=VARS) where {X0,X,Q,S,T<:Union{HarmonicVariables,Variables}}
  nv = nvars(m)
  nhv = nhvars(m)
  n = T == HarmonicVariables ? nhv : nv
  C = similar(m.x0, complex(eltype(m.x0)), n, n)
  C .= 0
  for i in 1:Int(nhv/2)
    C[2*i-1,2*i-1] = 1/sqrt(2)
    C[2*i-1,2*i] = 1/sqrt(2)

    C[2*i,2*i-1] = complex(0,-1/sqrt(2))
    C[2*i,2*i] = complex(0,1/sqrt(2))
  end

  if iscoasting(m) && T == Variables
    C[nv,nv] = 1
  end
  return C
end

function c_jacobian(m::TaylorMap{X0,X,Q,S}, ::T=VARS) where {X0<:StaticArray,X,Q,S,T<:Union{HarmonicVariables,Variables}}
  if T == Variables && iscoasting(m)
    nv = nvars(m)
    return StaticArrays.sacollect(SMatrix{nv,nv,eltype(X0)}, 
    begin
      if fld1(col,2) != fld1(row,2) 
        0
      elseif row == nv && col == nv
        1
      else # then we are in the block
        if mod1(row,2) == 1 # First row of c is just 1/sqrt(2)
            1/sqrt(2)
        else # second row of c is either -im/sqrt(2) or im/sqrt(2)
          if mod1(col,2) == 1
            complex(0,-1/sqrt(2))
          else
            complex(0,1/sqrt(2))
          end
        end
      end
    end for col in 1:nv for row in 1:nv)
  else
    nhv = nhvars(m)
    return StaticArrays.sacollect(SMatrix{nhv,nhv,eltype(X0)}, 
    begin
      if fld1(col,2) != fld1(row,2) 
        0
      else # then we are in the block
        if mod1(row,2) == 1 # First row of c is just 1/sqrt(2)
          1/sqrt(2)
        else # second row of c is either -im/sqrt(2) or im/sqrt(2)
          if mod1(col,2) == 1
            complex(0,-1/sqrt(2))
          else
            complex(0,1/sqrt(2))
          end
        end
      end
    end for col in 1:nhv for row in 1:nhv)
  end
end

function c_map(m::DAMap)
  c=zero(complex(typeof(m)), m)
  setray!(c.x, x_matrix=c_jacobian(m))
  if !isnothing(c.q)
    setquat!(c.q, q=I)
  end
  return c
end

function ci_map(m::DAMap)
  cinv=zero(complex(typeof(m)), m)
  setray!(cinv.x, x_matrix=ci_jacobian(m))
  if !isnothing(cinv.q)
    setquat!(cinv.q, q=I)
  end
  return cinv
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