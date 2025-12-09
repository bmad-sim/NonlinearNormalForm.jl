#=

Contains utility routines for getting TPSA info from a map/vector field,
getting the jacobian or transpose-jacobian (jacobiant) from a map/vector 
field array type promotion, checking if the last plane of a map is coasting.

=#
# =================================================================================== #
# Helper functions
getinit(m::Union{TaylorMap,VectorField}) = TI.getinit(first(m.v))
ndiffs(m::TaylorMap) = length(m.v)
ndiffs(m::TaylorMap{<:Any,V}) where {V} = length(V)
ndiffs(F::VectorField) = TI.ndiffs(first(F.v))
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.v))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.v))

# NNF-specific helpers:
nvars(m::TaylorMap) = length(m.v0)
nvars(m::TaylorMap{V0}) where {V0<:StaticArray} = length(V0)
nvars(F::VectorField) = length(F.v)
nvars(F::VectorField{V0}) where {V0<:StaticArray} = length(V0)
nparams(m::Union{TaylorMap,VectorField}) = ndiffs(m) - nvars(m)
nhvars(m::Union{TaylorMap,VectorField}) = iseven(nvars(m)) ? nvars(m) : nvars(m)-1 # number of "harmonic" variables

iscoasting(m::Union{TaylorMap,VectorField}) = !iseven(nvars(m))
# =================================================================================== #
# Tool to "factor" out a first order dependence on a particular variable
# e.g. for a TPSA with monomial [1,2,3,4], for 1 will return same TPSA but 
# that monomial is now [0,2,3,4]
# This is useful for lattice functions (first order) with nonlinear parameter dependence

factor_out(t, var::Int) = (out = zero(t); factor_out!(out, t, var))

function factor_out!(out, t, var::Int)
  TI.is_tps_type(typeof(t)) isa TI.IsTPSType || error("Function only accepts TPS types")
  TI.is_tps_type(typeof(out)) isa TI.IsTPSType || error("Function only accepts TPS types")
  nn = ndiffs(t)
  v = Ref{TI.numtype(t)}()
  tmpmono = zeros(UInt8, nn) 

  idx = TI.cycle!(t, 0, mono=tmpmono, val=v)
  while idx > 0
    if tmpmono[var] != 0
      tmpmono[var] -= 1
      TI.setm!(out, v[], tmpmono)
    end
    idx = TI.cycle!(t, idx, mono=tmpmono, val=v)
  end
  return out
end

factor_in(t, var::Int, n::Int=1) = (out = zero(t); factor_in!(out, t, var, n))

function factor_in!(out, t, var::Int, n::Int=1)
  TI.is_tps_type(typeof(t)) isa TI.IsTPSType || error("Function only accepts TPS types")
  TI.is_tps_type(typeof(out)) isa TI.IsTPSType || error("Function only accepts TPS types")
  nn = ndiffs(t)
  v = Ref{TI.numtype(t)}()
  tmpmono = zeros(UInt8, nn) 

  # First do scalar part
  tmpmono[var] += n
  TI.setm!(out, TI.geti(t, 0), tmpmono)

  idx = TI.cycle!(t, 0, mono=tmpmono, val=v)
  while idx > 0
    tmpmono[var] += n
    TI.setm!(out, v[], tmpmono)
    idx = TI.cycle!(t, idx, mono=tmpmono, val=v)
  end
  return out
end

fast_var_par(t, var::Int, nv::Int) = (out = zero(t); fast_var_par!(out, t, var, nv); return out)

# This will par out the polynomial with only 
# 1st order dependence on variable, nonlinear 
# parameter dependence
function fast_var_par!(out, t, var::Int,  nv::Int)
  TI.is_tps_type(typeof(t)) isa TI.IsTPSType || error("Function only accepts TPS types")
  TI.is_tps_type(typeof(out)) isa TI.IsTPSType || error("Function only accepts TPS types")
  nn = ndiffs(t)
  np = nn-nv
  v = Ref{TI.numtype(t)}()
  tmpmono = zeros(UInt8, nn) 

  idx = TI.cycle!(t, 0, mono=tmpmono, val=v)
  while idx > 0
    if tmpmono[var] == 1 && all(t->t==0, view(tmpmono, 1:(var-1))) && all(t->t==0, view(tmpmono, (var+1):nv))
      tmpmono[var] -= 1
      TI.setm!(out, v[], tmpmono)
    end
    idx = TI.cycle!(t, idx, mono=tmpmono, val=v)
  end
  return out
end

# =================================================================================== #
# Jacobian/jacobiant
# Useful options should be 1) harmonic variables, 2) variables, and 3) variables + parameters
abstract type OptionType end
struct HarmonicVariables <: OptionType end  # e.g. 4x4 matrix
struct CoastVariables <: OptionType end # e.g. 1x5 matrix
struct Variables <: OptionType end # e.g. 4x4 or 5x5 (coasting beam) matrix
struct HarmonicVariablesAndParameters <: OptionType end # e.g. 4 x (nv+np) matrix
struct VariablesAndCoastParameter <: OptionType end # e.g. 5 x 6 matrix 
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
const VARS_CPARAM  = VariablesAndCoastParameter()
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
  elseif T == VariablesAndCoastParameter
    nrows = nvars(m) + (isodd(nvars(m)) ? 1 : 0)
    ncols = nrows
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

function jacobian(m::TaylorMap{V0,V,Q,S}, which::T=VARS) where {V0,V,Q,S,T<:OptionType}
  nrows, ncols, row_start, col_start = getjacinfo(m, which)
  M = similar(m.v0, nrows, ncols)
  for col in col_start:col_start+ncols-1
    for row in row_start:row_start+nrows-1
      M[row-row_start+1,col-col_start+1] = TI.geti(m.v[row], col)
    end
  end
  return M
end

function jacobian(m::TaylorMap{V0,V,Q,S}, which::T=VARS) where {V0,V<:StaticArray,Q,S,T<:OptionType}
  nrows, ncols, row_start, col_start = getjacinfo(m, which)
  return StaticArrays.sacollect(SMatrix{nrows,ncols,eltype(V0)}, 
  TI.geti(m.v[row], col) for col in col_start:col_start+ncols-1 for row in row_start:row_start+nrows-1)
end

function jacobiant(m::TaylorMap{V0,V,Q,S}, which::T=VARS) where {V0,V,Q,S,T<:OptionType}
  nrowst, ncolst, rowt_start, colt_start = getjacinfo(m, which)
  M = similar(m.v0, ncolst,nrowst)
  for colt in colt_start:colt_start+ncolst-1
    for rowt in rowt_start:rowt_start+nrowst-1
      M[colt-colt_start+1,rowt-rowt_start+1] = TI.geti(m.v[rowt], colt)
    end
  end
  return M
end

function jacobiant(m::TaylorMap{V0,V,Q,S}, which::T=VARS) where {V0,V<:StaticArray,Q,S,T<:OptionType}
  return transpose(jacobian(m, which))
end

# =================================================================================== #
# Phasor basis transformations
# inv(c) is from_phasor
function ci_jacobian(m::TaylorMap{V0,V,Q,S}, ::T=VARS) where {V0,V,Q,S,T<:Union{HarmonicVariables,Variables}}
  nv = nvars(m)
  nhv = nhvars(m)
  n = T == HarmonicVariables ? nhv : nv
  CI = similar(m.v0, complex(eltype(m.v0)), n, n)
  CI .= 0
  for i in 1:Int(nhv/2)
    # x_new = 1/sqrt(2)*(v+im*p)
    CI[2*i-1,2*i-1] = 1/sqrt(2)
    CI[2*i-1,2*i] = complex(0,1/sqrt(2))

    # p_new = 1/sqrt(2)*(v-im*p)
    CI[2*i,2*i-1] = 1/sqrt(2)
    CI[2*i,2*i] = complex(0,-1/sqrt(2))
  end

  if T == Variables && iscoasting(m)
    CI[nv,nv] = 1
  end
  return CI
end

function ci_jacobian(m::TaylorMap{V0,V,Q,S}, ::T=VARS) where {V0<:StaticArray,V,Q,S,T<:Union{HarmonicVariables,Variables}}
  if T == Variables && iscoasting(m)
    nv = nvars(m)
    return StaticArrays.sacollect(SMatrix{nv,nv,complex(eltype(V0))}, 
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
    return StaticArrays.sacollect(SMatrix{nhv,nhv,complex(eltype(V0))}, 
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
function c_jacobian(m::TaylorMap{V0,V,Q,S}, ::T=VARS) where {V0,V,Q,S,T<:Union{HarmonicVariables,Variables}}
  nv = nvars(m)
  nhv = nhvars(m)
  n = T == HarmonicVariables ? nhv : nv
  C = similar(m.v0, complex(eltype(m.v0)), n, n)
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

function c_jacobian(m::TaylorMap{V0,V,Q,S}, ::T=VARS) where {V0<:StaticArray,V,Q,S,T<:Union{HarmonicVariables,Variables}}
  if T == Variables && iscoasting(m)
    nv = nvars(m)
    return StaticArrays.sacollect(SMatrix{nv,nv,complex(eltype(V0))}, 
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
    return StaticArrays.sacollect(SMatrix{nhv,nhv,complex(eltype(V0))}, 
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
  setray!(c.v, v_matrix=c_jacobian(m))
  if !isnothing(c.q)
    setquat!(c.q, q=I)
  end
  return c
end

function ci_map(m::DAMap)
  cinv=zero(complex(typeof(m)), m)
  setray!(cinv.v, v_matrix=ci_jacobian(m))
  if !isnothing(cinv.q)
    setquat!(cinv.q, q=I)
  end
  return cinv
end

# =================================================================================== #
# Get/set scalar part of orbital ray
function getscalar(m::TaylorMap{V0,V,Q,S}) where {V0,V,Q,S}
  nv = nvars(m)
  return map(t->TI.geti(t, 0), view(m.v, 1:nv))
end

function getscalar(m::TaylorMap{V0,V,Q,S}) where {V0,V<:StaticArray,Q,S}
  nv = nvars(m)
  return StaticArrays.sacollect(SVector{nv,eltype(V0)}, TI.geti(m.v[i], 0) for i in 1:nv)
end

function setscalar!(
  m::TaylorMap{V0,V,Q,S}, 
  xs::Number; 
  scl0::Union{Nothing,Number}=nothing,
  scl1::Number=1
) where {V0,V,Q,S}
  nv = nvars(m)
  if isnothing(scl0)
    for i in 1:nv
      TI.seti!(m.v[i], scl1*xs, 0)
    end
  else
    for i in 1:nv
      TI.seti!(m.v[i], TI.geti(m.v[i], 0)*scl0 + scl1*xs, 0)
    end
  end
end

function setscalar!(
  m::TaylorMap{V0,V,Q,S}, 
  xs::AbstractArray; 
  scl0::Union{Nothing,Number}=nothing,
  scl1::Number=1
) where {V0,V,Q,S}

  nv = nvars(m)
  if isnothing(scl0)
    for i in 1:nv
      TI.seti!(m.v[i], xs[i], 0)
    end
  else
    for i in 1:nv
      TI.seti!(m.v[i], TI.geti(m.v[i], 0)*scl0 + scl1*xs[i], 0)
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
  return Quaternion{G}
end

# =================================================================================== #
# Coast check

function check_coast(v)
  if abs(TI.geti(last(v), 0)) ≈ 0 # if scalar part is 0
    cycleidx = TI.cycle!(last(v), 0)
    if cycleidx == length(v) && abs(TI.geti(last(v), length(v)) - 1) ≈ 0
      return true
    end
  end
  return false
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

function checksymp(m::TaylorMap)
  nv = nvars(m)
  if iscoasting(m)
    nv += 1
  end
  Sij = S(nv)
  mo = maxord(m)
  out = zeros(eltype(m.v), nv, nv)
  for i in 1:nv
    for j in 1:nv
      grad1 = [TI.deriv(m.v[j], ell) for ell in 1:nv]
      grad2 = [TI.deriv(m.v[i], k) for k in 1:nv]
      # per thesis, for given k (index) we sum
      s1 = [sum([TI.cutord(Sij[k,ell]*grad1[ell],mo) for ell in 1:nv]) for k in 1:nv]
      out[i,j] = sum([TI.cutord(grad2[k]*s1[k],mo) for k in 1:nv])
    end
  end
  return out-S
end
# =================================================================================== #