#=

Defines the types used throughout the package. Specifically, the abstract 
type TaylorMap, and concrete types DAMap and TPSAMap (which differ only in 
concatenation and inversion rules). 

=#

"""
    TaylorMap{X0,X,Q,S}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. For a periodic system, 
the `DAMap` is expanded around the closed orbit, while the `TPSAMap` may be expanded around 
any arbitrary trajectory. For normal form analysis of periodic maps, using a `DAMap` ensures 
no truncation error up to the chosen truncation order.

# Fields
- `x0::X0` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `x::X`   -- Orbital ray as truncated power series, expansion around `x0`, with scalar part equal to EXIT coordinates of map
- `q::Q`   -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `s::S`   -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

# Type Requirements
- `X0 <: AbstractVector{<:Number}` where `ismutabletype(X0) == true` 
- `X <: AbstractVector{<:TPS}` where `GTPSA.numtype(X) == eltype(X0)`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(X)`
- `S <: Union{AbstractMatrix{<:Number},Nothing}` where `S != Nothing` then `eltype(S) == GTPSA.numtype(eltype(X))` AND `ismutabletype(S) == true`

Because the `TPS` type is `mutable` and `GTPSA.jl` provides in-place functions for modifying 
`TPS`s, at the lowest level, all operations on `TaylorMap`s are in-place for performance. Therefore, 
the `x` and `q` arrays which contain `TPS`s may be `immutable`, e.g. the orbital ray `x` may be an 
`SVector` from the `StaticArrays.jl` package, and the `Quaternion` type which is taken from 
`ReferenceFrameRotations.jl` is already `immutable`. The default for the orbital ray is `SVector`.
The `x0` and `s` arrays contain `immutable` number types, and so these arrays MUST be `mutable`.
"""
abstract type TaylorMap{X0<:AbstractVector{<:Number},X<:AbstractVector{<:TPS},Q<:Union{Quaternion{<:TPS},Nothing},S<:Union{AbstractMatrix{<:Number},Nothing}} end 

@inline function checkmapsanity(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  ismutabletype(X0) || error("Reference orbit array must be mutable: $X0 is immutable.")
  eltype(X0) == GTPSA.numtype(eltype(X)) || error("Reference orbit number type $(eltype(X0)) must be $(GTPSA.numtype(eltype(X))) (equal to the type of the orbital ray scalar part)")
  Q == Nothing || eltype(Q) == eltype(X) || error("Quaternion number type $(eltype(Q)) must be $(eltype(X)) (equal to orbital ray)")
  S == Nothing || ismutabletype(S) || error("Stochastic envelope matrix must be mutable: $S is immutable.")
  S == Nothing || eltype(S) == GTPSA.numtype(eltype(X)) || error("Stochastic envelope matrix number type $(eltype(S)) must be $(GTPSA.numtype(eltype(X))) (equal to the type of the orbital ray scalar part)")
end

struct DAMap{X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  x0::X0     # Entrance value of map
  x::X     # Expansion around x0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function DAMap(x0, x, q, s)
    m = new{typeof(x0),typeof(x),typeof(q),typeof(s)}(x0, x, q, s)
    checkmapsanity(m)
    return m
  end
end

struct TPSAMap{X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  x0::X0     # Entrance value of map
  x::X      # Expansion around x0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function TPSAMap(x0, x, q, s)
    m = new{typeof(x0),typeof(x),typeof(q),typeof(s)}(x0, x, q, s)
    checkmapsanity(m)
    return m
  end
end

"""
    VectorField{X,Q}

Lie operator which acts on `DAMap`s e.g. dM/dt = FM where F is the `VectorField` and 
`M` is the `DAMap`. `F` can be promoted to a `DAMap` using `exp(F)`.

# Fields
- `x::X` -- Orbital ray as truncated power series, expansion around `x0`, with scalar part equal to EXIT coordinates of map
- `q::Q` -- `Quaternion` as truncated power series if spin is included, else `nothing`

# Type Requirements
- `X <: AbstractVector{<:TPS}`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(X)`
"""
struct VectorField{X<:AbstractVector{<:TPS},Q<:Union{Quaternion{<:TPS},Nothing}}
  x::X
  q::Q     
  function VectorField(x, q)
    vf = new{typeof(x),typeof(q)}(x, q)
    checkvfsanity(vf)
    return vf
  end
end

@inline function checkmapsanity(m::VectorField{X,Q}) where {X,Q}
  Q == Nothing || eltype(Q) == eltype(X) || error("Quaternion number type $(eltype(Q)) must be $(eltype(X)) (equal to orbital ray)")
end
 

const UseType = Union{Descriptor, TPS, DAMap, TPSAMap, VectorField}

# =================================================================================== #
# Field promotion rules

# Promotion utility functions for easy override by external array packages
# Should promote the element type and return it as a vector (x,x0) or matrix (s)
# with dimensionality specified by Val.
promote_x0_type(::Type{X0}, ::Type{G}) where {X0,G<:Union{Number,Complex}} = _promote_array_eltype(X0, G, Val{1})
promote_x_type(::Type{X}, ::Type{G}) where {X,G<:Union{Number,Complex}} = _promote_array_eltype(X, G, Val{1})
promote_q_type(::Type{Q}, ::Type{G}) where {Q,G<:Union{Number,Complex}} = Quaternion{promote_type(eltype(Q), G)}
promote_s_type(::Type{S}, ::Type{G}) where {S,G<:Union{Number,Complex}} = _promote_array_eltype(S, G, Val{2})
real_x0_type(::Type{X0}) where {X0} = _real_array_eltype(X0, Val{1})
real_x_type(::Type{X}) where {X} = _real_array_eltype(X, Val{1})
real_q_type(::Type{Q}) where {Q} = Quaternion{real(eltype(Q))} 
real_s_type(::Type{S}) where {S} = _real_array_eltype(S, Val{2})

# Fallback for Nothing type of quaternion and stochastic envelope matrix
promote_q_type(::Type{Nothing}, ::Type{G}) where {G<:Union{Number,Complex}} = Nothing
promote_s_type(::Type{Nothing}, ::Type{G}) where {G<:Union{Number,Complex}} = Nothing
real_q_type(::Type{Nothing}) = Nothing
real_s_type(::Type{Nothing}) = Nothing

for t = (:DAMap, :TPSAMap)
@eval begin    

function promote_rule(::Type{$t{X0,X,Q,S}}, ::Type{G}) where {X0,X,Q,S,G<:Union{Number,Complex}}
  outS = promote_x0_type(X0, G)
  outT = promote_x_type(X, G)
  outU = promote_q_type(Q, G)
  outV = promote_s_type(S, G)
  return $t{outS,outT,outU,outV}
end

# --- complex type ---
function complex(type::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  return promote_type(type, complex(GTPSA.numtype(eltype(X))))
end

# --- real type ---
function real(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  outS = real_x0_type(X0)
  outT = real_x_type(X)
  outU = real_q_type(Q)
  outV = real_s_type(S)
  return $t{outS,outT,outU,outV}
end

end
end

function promote_rule(::Type{VectorField{X,Q}}, ::Type{G}) where {X,Q,G<:Union{Number,Complex}}
  outT = promote_x_type(X, G)
  outU = promote_q_type(Q, G)
  return VectorField{outT,outU}
end

# --- complex type ---
function complex(type::Type{VectorField{X,Q}}) where {X,Q}
  return promote_type(type, complex(GTPSA.numtype(eltype(X))))
end

# --- real type ---
function real(::Type{VectorField{X,Q}}) where {X,Q}
  outT = real_x_type(X)
  outU = real_q_type(Q)
  return VectorField{outT,outU}
end


# =================================================================================== #
# Field initialization functions.

# These may be overrided by external array packages.
function init_x0(::Type{X0}, use::UseType) where {X0}
  nv = numvars(use)
  x0 = similar(X0, nv)
  x0 .= 0
  return x0
end

function init_x(::Type{X}, use::UseType) where {X}
  desc = getdesc(use)
  nv = numvars(desc)
  nn = numnn(desc)
  x = similar(X, nn)
  for i in 1:nv
    x[i] = eltype(x)(use=desc)
  end
  # use same parameters if use isa TaylorMap and eltype(x) == eltype(use.x)
  if use isa TaylorMap && eltype(x) == eltype(use.x)
    x[nv+1:nn] .= view(use.x, nv+1:nn)
  else # allocate
    if eltype(x) == TPS{Float64}
      x[nv+1:nn] .= params(desc)
    else
      x[nv+1:nn] .= complexparams(desc)
    end
  end
  return x
end

function init_q(::Type{Q}, use::UseType) where {Q}
  if Q != Nothing
    desc = getdesc(use)
    q0 = eltype(Q)(use=desc)
    q1 = eltype(Q)(use=desc)
    q2 = eltype(Q)(use=desc)
    q3 = eltype(Q)(use=desc)
    q = Quaternion(q0,q1,q2,q3)
  else
    q = nothing
  end
  return q
end

function init_s(::Type{S}, use::UseType) where {S}
  if S != Nothing
    nv = numvars(use)
    s = similar(S, nv, nv)
    s .= 0
  else
    s = nothing
  end
  return s
end






