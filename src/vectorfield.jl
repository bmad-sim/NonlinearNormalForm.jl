#=

Defines the VectorField type, promotion rules, and constructors.

=#

# =================================================================================== #
# Type

"""
    VectorField{X,Q}

Lie operator which acts on `DAMap`s e.g. dM/dt = FM where F is the `VectorField` and 
`M` is the `DAMap`. `F` can be promoted to a `DAMap` using `exp(F)`.

# Fields
- `x::X` -- Orbital ray as truncated power series, expansion with scalar part equal to EXIT coordinates of map
- `q::Q` -- `Quaternion` as truncated power series if spin is included, else `nothing`

# Type Requirements
- `X <: AbstractVector{<:TPS}`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(X)`
"""
struct VectorField{X<:AbstractVector,Q<:Union{Quaternion,Nothing},}
  x::X
  q::Q     
  function VectorField(x, q)
    vf = new{typeof(x),typeof(q)}(x, q)
    checkvfsanity(vf)
    return vf
  end
end

@inline function checkvfsanity(::VectorField{X,Q}) where {X,Q}
  # Static checks:
  TI.is_tps_type(eltype(X)) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")
  Q == Nothing || eltype(Q) == eltype(X) || error("Quaternion number type $(eltype(Q)) must be $(eltype(X)) (equal to orbital ray)")

  # Runtime checks:
  nvars(first(m.x)) == length(m.x) || error("VectorField orbital ray length disagrees with number of variables in TPSA")
  Q == Nothing || getdef(first(m.q)) == getdef(first(m.x)) || error("Quaternion TPSA definition disagrees with orbital ray TPSA definition")
end

# =================================================================================== #
# Field initialization functions.

# These may be overrided by external array packages.
function init_vf_x(::Type{X}, def::AbstractTPSADef) where {X<:AbstractVector}
  nv = nvars(def)
  nn = ndiffs(def)
  x = similar(X, nv)
  for i in 1:nv
    x[i] = TI.init_tps(TI.numtype(eltype(X)), def) 
  end
  return x
end

init_vf_q(::Type{Q}, def::AbstractTPSADef) where {Q<:Union{Nothing,Quaternion}} = init_map_Q(Q, def)

# =================================================================================== #
# StaticArray field initialization functions.

# Consistency checks are made by the `checkmapsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA
function init_vf_x(::Type{X}, def::AbstractTPSADef) where {X<:StaticVector}
  nv = nvars(def)
  x = StaticArrays.sacollect(X, TI.init_tps(TI.numtype(eltype(X)), def) for i in 1:length(X))
  return x
end

# =================================================================================== #
# Promotion rules

# VectorField:
function promote_rule(::Type{VectorField{X,Q}}, ::Type{G}) where {X,Q,G<:Union{Number,Complex}}
  out_X = similar_eltype(X, promote_type(eltype(X), G))
  out_Q = similar_eltype(Q, promote_type(eltype(Q), G))
  return VectorField{out_X,out_Q}
end

# --- complex type ---
function complex(type::Type{VectorField{X,Q}}) where {X,Q}
  return promote_type(type, complex(TI.numtype(eltype(X))))
end

# --- real type ---
function real(::Type{VectorField{X,Q}}) where {X,Q}
  out_X = similar_eltype(X, real(eltype(X)))
  out_Q = isnothing(Q) ? Nothing : similar_eltype(Q, real(eltype(Q)))
  return VectorField{out_X,out_Q}
end

# =================================================================================== #
# Constructors


# =================================================================================== #









