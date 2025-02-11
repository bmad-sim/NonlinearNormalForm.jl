#=

Defines the types abstract type TaylorMap, and concrete types DAMap 
and TPSAMap (which differ only in concatenation and inversion rules), 
promotion rules, and constructors.

=#

# =================================================================================== #
# Types

"""
    TaylorMap{X0,X,Q,S}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. `DAMap`s have a coordinate 
system chosen so that the expansion point is always around zero, e.g. Δx = x, while `TPSAMap`s 
can have any coordinate system/expansion point, but therefore will accrue truncation error when 
trying to compose two `TPSAMap`s with differing expansion points. For normal form analysis of 
periodic maps, using a `DAMap` ensures no truncation error up to the chosen truncation order.

Any truncated power series (TPS) type supported by `TPSAInterface.jl` is allowed for use in 
a `TaylorMap`. Henceforth we will generically refer to this type as a `TPS`

# Fields
- `x0::X0` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `x::X`   -- Orbital ray as truncated power series, expansion around `x0`, with scalar part equal to EXIT coordinates of map
- `q::Q`   -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `s::S`   -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

# Type Requirements
- `X0 <: AbstractVector{<:Number}` where `ismutabletype(X0) == true` 
- `X <: AbstractVector{<:TPS}` where `TPSAInterface.numtype(X) == eltype(X0)`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(X)`
- `S <: Union{AbstractMatrix{<:Number},Nothing}` where `S != Nothing` then `eltype(S) == TPSAInterface.numtype(eltype(X))` AND `ismutabletype(S) == true`

Because the TPS type is `mutable` and `TPSAInterface.jl` provides in-place functions for modifying 
TPSs, at the lowest level, all operations on `TaylorMap`s are in-place for performance. Therefore, 
the `x` and `q` arrays which contain TPSs may be `immutable`, e.g. the orbital ray `x` may be an 
`SVector` from the `StaticArrays.jl` package, and the `Quaternion` type which is taken from 
`ReferenceFrameRotations.jl` is already `immutable`. The default for the orbital ray is `SVector`.
The `x0` and `s` arrays contain `immutable` number types, and so these arrays MUST be `mutable`.
"""
abstract type TaylorMap{X0<:AbstractVector,X<:AbstractVector,Q<:Union{Quaternion,Nothing},S<:Union{AbstractMatrix,Nothing}} end 

@inline function checkmapsanity(m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S}
  # Static checks:
  ismutabletype(X0) || error("Reference orbit array must be mutable: $X0 is immutable.")
  TI.is_tps_type(eltype(X)) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")
  eltype(X0) == TI.numtype(eltype(X)) || error("Reference orbit number type $(eltype(X0)) must be $(TI.numtype(eltype(X))) (equal to the type of the orbital ray scalar part)")
  Q == Nothing || eltype(Q) == eltype(X) || error("Quaternion number type $(eltype(Q)) must be $(eltype(X)) (equal to orbital ray)")
  S == Nothing || ismutabletype(S) || error("Stochastic envelope matrix must be mutable: $S is immutable.")
  S == Nothing || eltype(S) == TI.numtype(eltype(X)) || error("Stochastic envelope matrix number type $(eltype(S)) must be $(TI.numtype(eltype(X))) (equal to the type of the orbital ray scalar part)")
  
  # Runtime checks:
  nvars(first(m.x)) == length(m.x0) || error("Reference orbit array length disagrees with number of variables in TPSA")
  ndiffs(first(m.x)) == length(m.x) || error("Orbital ray length disagrees with number of differentials (variables + parameters) in TPSA")
  Q == Nothing || getdef(first(m.q)) == getdef(first(m.x)) || error("Quaternion TPSA definition disagrees with orbital ray TPSA definition")
  S == Nothing || size(S) == (length(m.x0), length(m.x0)) || error("Size of stochastic matrix disagrees with number of variables in TPSA")
end

struct DAMap{X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  x0::X0    # Entrance value of map
  x::X      # Expansion around x0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function DAMap(x0, x, q, s)
    m = new{typeof(x0),typeof(x),typeof(q),typeof(s)}(x0, x, q, s)
    checkmapsanity(m)
    return m
  end
end

struct TPSAMap{X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  x0::X0    # Entrance value of map
  x::X      # Expansion around x0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function TPSAMap(x0, x, q, s)
    m = new{typeof(x0),typeof(x),typeof(q),typeof(s)}(x0, x, q, s)
    checkmapsanity(m)
    return m
  end
end

# =================================================================================== #
# Field initialization functions.
# Note that even when a TPSA definition is inferrable from the types (X0, X, ...),
# we still explicitly pass the def in case it is not. This is redundant, however ensures 
# flexibility for different TPSA packages/modes.

# These may be overrided by external array packages.
function init_map_x0(::Type{X0}, def::AbstractTPSADef) where {X0<:AbstractVector}
  nv = nvars(def)
  x0 = similar(X0, nv)
  x0 .= 0
  return x0
end

function init_map_x(::Type{X}, def::AbstractTPSADef, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:AbstractVector}
  nv = nvars(def)
  nn = ndiffs(def)
  x = similar(X, nn)
  for i in 1:nv
    x[i] = TI.init_tps(TI.numtype(eltype(X)), def) 
  end
  # use same parameters if reuse isa TaylorMap and eltype(x) == eltype(reuse.x)
  if reuse isa TaylorMap && eltype(x) == eltype(reuse.x) && def == getdef(reuse)
    x[nv+1:nn] .= view(reuse.x, nv+1:nn)
  else # allocate
    for i in nv+1:nn
      x[i] = TI.init_tps(TI.numtype(eltype(X)), def) 
      TI.seti!(x[i], 1, i)
    end
  end
  return x
end

function init_map_q(::Type{Q}, def::AbstractTPSADef) where {Q<:Union{Nothing,Quaternion}}
  if Q != Nothing
    q0 = TI.init_tps(TI.numtype(eltype(Q)), def) 
    q1 = TI.init_tps(TI.numtype(eltype(Q)), def) 
    q2 = TI.init_tps(TI.numtype(eltype(Q)), def) 
    q3 = TI.init_tps(TI.numtype(eltype(Q)), def) 
    q = Quaternion(q0,q1,q2,q3)
  else
    q = nothing
  end
  return q
end

function init_map_s(::Type{S}, def::AbstractTPSADef) where {S<:Union{Nothing,AbstractMatrix}}
  if S != Nothing
    nv = nvars(def)
    s = similar(S, nv, nv)
    s .= 0
  else
    s = nothing
  end
  return s
end

# =================================================================================== #
# StaticArray field initialization functions.

# Consistency checks are made by the `checkmapsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA
function init_map_x0(a::Type{X0}, ::AbstractTPSADef) where {X0<:StaticVector}
  x0 = StaticArrays.sacollect(X0, 0 for i in 1:length(a))
  return x0
end

function init_map_x(::Type{X}, def::AbstractTPSADef, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:StaticVector}
  nv = nvars(def)
  # reuse parameters if applicable
  if reuse isa TaylorMap && eltype(X) == eltype(reuse.x) && def == getdef(reuse)
    x = StaticArrays.sacollect(X, (i <= nv ? TI.init_tps(TI.numtype(eltype(X)), def) :  reuse.x[i]) for i in 1:length(X))
  else # allocate
    x = StaticArrays.sacollect(X, (i <= nv ? 
                                    TI.init_tps(TI.numtype(eltype(X)), def) : 
                                    (t = TI.init_tps(TI.numtype(eltype(X)), def); TI.seti!(t, 1, i); t)) for i in 1:length(X))
  end
  return x
end

function init_map_s(::Type{S}, ::AbstractTPSADef) where {S<:StaticMatrix}
  s = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  return s
end

# =================================================================================== #
# Promotion rules

for t = (:DAMap, :TPSAMap)
@eval begin    

function promote_rule(::Type{$t{X0,X,Q,S}}, ::Type{G}) where {X0,X,Q,S,G<:Union{Number,Complex}}
  out_X0 = similar_eltype(X0, promote_type(eltype(X0), G))
  out_X = similar_eltype(X, promote_type(eltype(X), G))
  out_Q = isnothing(Q) ? Nothing : similar_eltype(Q, promote_type(eltype(Q), G))
  out_S = isnothing(S) ? Nothing : similar_eltype(S, promote_type(eltype(S), G))
  return $t{out_X0,out_X,out_Q,out_S}
end

# --- complex type ---
function complex(type::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  return promote_type(type, complex(TI.numtype(eltype(X))))
end

# --- real type ---
function real(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  out_X0 = similar_eltype(X0, real(eltype(X0)))
  out_X = similar_eltype(X, real(eltype(X)))
  out_Q = isnothing(Q) ? Nothing : similar_eltype(Q, real(eltype(Q)))
  out_S = isnothing(S) ? Nothing : similar_eltype(S, real(eltype(S)))
  return $t{out_X0,out_X,out_Q,out_S}
end

end
end


# =================================================================================== #
# Constructors
for t = (:DAMap, :TPSAMap)
@eval begin

# Lowest-level, internal
function _zero(::Type{$t{X0,X,Q,S}}, def::AbstractTPSADef, reuse::Union{TaylorMap,Nothing}) where {X0,X,Q,S}
  out_x0 = init_map_x0(X0, def)
  out_x = init_map_x(X, def, reuse)
  out_q = init_map_q(Q, def)
  out_s = init_map_s(S, def)
  out_m = $t(out_x0, out_x, out_q, out_s)
  return out_m
end

function zero(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  return _zero($t{X0,X,Q,S}, getdef(eltype(X)), nothing)
end

# Copy ctor including optional TPSA def change
function $t(m::TaylorMap; def::AbstractTPSADef=getdef(m))
  checktpsas(m, def)
  X0 = typeof(m.x0)
  X = similar_eltype(typeof(m.x), TI.init_tps_type(eltype(X0), def))
  Q = isnothing(m.q) ? Nothing : Quaternion{TI.init_tps_type(eltype(X0), def)}
  S = typeof(m.s)
  out_m = _zero($t{X0,X,Q,S}, def, m)
  copy!(out_m, m)
  return out_m
end

# Kwarg ctor:
function $t(;
  def::Union{AbstractTPSADef,Nothing}=nothing,
  x0::Union{AbstractVector,Nothing}=nothing,
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
  s::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
  stochastic::Union{Bool,Nothing}=nothing,
) 
  if isnothing(def)
    if !isnothing(x) && TI.is_tps_type(eltype(x)) isa TI.IsTPSType
      def = getdef(first(x))
    elseif !isnothing(q) && TI.is_tps_type(eltype(q)) isa TI.IsTPSType
      def = getdef(first(q))
    else
      error("No TPSA definition has been provided, nor is one inferrable from the input arguments!")
    end
  end


  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? TI.numtype(eltype(t)) : Float64), (x0, x, q, s))...)
  TW = TI.init_tps_type(W, def)
  X0 = isnothing(x0) ? @_DEFAULT_X0(nvars(def)){W} : similar_eltype(typeof(x0), W)
  X = isnothing(x) ? @_DEFAULT_X(ndiffs(def)){TW} : similar_eltype(typeof(x), TW)
  
  if isnothing(spin)
    Q = isnothing(q) && isnothing(q_map) ? Nothing : Quaternion{TW}
  elseif spin
    Q = Quaternion{TW}
  else
    Q = Nothing
  end

  if isnothing(stochastic)
    S = isnothing(s) ? Nothing : similar_eltype(typeof(s), W)
  elseif stochastic
    S = isnothing(s) ? @_DEFAULT_S(nvars(def)){W} : similar_eltype(typeof(s), W)
  else
    S = Nothing
  end

  # Construct map:
  out_x0 = init_map_x0(X0, def)
  out_x = init_map_x(X, def)
  out_q = init_map_q(Q, def)
  out_s = init_map_s(S, def)

  out_m = $t(out_x0, out_x, out_q, out_s)


  if !isnothing(x0)
    out_m.x0 .= x0
  end

  setray!(out_m.x, x=x, x_matrix=x_matrix)
  if !isnothing(out_m.q) && !isnothing(q)
    setquat!(out_m.q, q=q, q_map=q_map)
  end

  if !isnothing(out_m.s) && !isnothing(s)
    out_m.s .= s
  end

  return out_m
end

end
end


"""
    zero(m::TaylorMap)

Creates a zero `m` with the same properties as `m` including TPSA definiton,
spin, and stochasticity.
"""
function zero(m::TaylorMap) 
  return _zero(typeof(m), getdef(m), m)
end

"""
    one(m::TaylorMap)
  
Creates an identity `m` with the same properties as `m`, including GTPSA
`Descriptor`, spin, and stochasticity.
"""
function one(m::TaylorMap)
  out_m = zero(m)
  nv = nvars(m)

  for i in 1:nv
    TI.seti!(out_m.x[i], 1, i)
  end

  if !isnothing(m.q)
    TI.seti!(out_m.q.q0, 1, 0)
  end

  return out_m
end

# =================================================================================== #
# composition
for t = (:DAMap, :TPSAMap)
@eval begin
  
function _compose!(m::$t, m2::$t, m1::$t; dospin::Bool=true, dostochastic::Bool=true, work_Q::Union{Nothing,Quaternion}=prep_work_Q(m))
  @assert !(m === m1) "Cannot compose_it!(m, m2, m1) with m === m1"
  @assert !(m === m2) "Cannot compose_it!(m, m2, m1) with m === m2"

  m.x0 .= m1.x0

  TI.compose!(view(m.x, 1:nv), view(m2.x, 1:nv), view(m1.x, 1:nn))

  # Spin:
  if !isnothing(m.q) && dospin
    TI.compose!(m.q.q0, m2.q.q0, m1.x)
    TI.compose!(m.q.q1, m2.q.q1, m1.x)
    TI.compose!(m.q.q2, m2.q.q2, m1.x)
    TI.compose!(m.q.q3, m2.q.q3, m1.x)
    mul!(m.q, m.q, m1.q) # This will be made faster eventually
  end 

  # Stochastic
  # MAKE THIS FASTER!
  if !isnothing(m.s) && dostochastic
    M2 = jacobian(m2)   
    m.s .= M2*m1.s*transpose(M2) + m2.s
  end

  return m
end

end
end
# =================================================================================== #