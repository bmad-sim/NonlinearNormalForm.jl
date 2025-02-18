#=

Defines the types abstract type TaylorMap, and concrete types DAMap 
and TPSAMap (which differ only in concatenation and inversion rules), 
promotion rules, constructors, and map-specific operators.

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
  ndiffs(first(m.x)) == length(m.x) || error("Orbital ray length disagrees with number of differentials in TPSA")
  Q == Nothing || getinit(first(m.q)) == getinit(first(m.x)) || error("Quaternion TPSA definition disagrees with orbital ray TPSA definition")
  S == Nothing || size(S) == (length(m.x0), length(m.x0)) || error("Size of stochastic matrix disagrees with number of variables in map")
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
# we still explicitly pass the init in case it is not. This is redundant, however ensures 
# flexibility for different TPSA packages/modes.

# These may be overrided by external array packages.
function init_map_x0(::Type{X0}, init::AbstractTPSAInit, nv::Integer) where {X0<:AbstractVector}
  x0 = similar(X0, nv)
  x0 .= 0
  return x0
end

function init_map_x(::Type{X}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:AbstractVector}
  nn = ndiffs(init)
  x = similar(X, nn)
  for i in 1:nv
    x[i] = TI.init_tps(TI.numtype(eltype(X)), init) 
  end
  # use same parameters if reuse isa TaylorMap and eltype(x) == eltype(reuse.x)
  if reuse isa TaylorMap && eltype(x) == eltype(reuse.x) && init == getinit(reuse)
    x[nv+1:nn] .= view(reuse.x, nv+1:nn)
  else # allocate
    for i in nv+1:nn
      x[i] = TI.init_tps(TI.numtype(eltype(X)), init) 
      TI.seti!(x[i], 1, i)
    end
  end
  return x
end

function init_map_q(::Type{Q}, init::AbstractTPSAInit) where {Q<:Union{Nothing,Quaternion}}
  if Q != Nothing
    q0 = TI.init_tps(TI.numtype(eltype(Q)), init) 
    q1 = TI.init_tps(TI.numtype(eltype(Q)), init) 
    q2 = TI.init_tps(TI.numtype(eltype(Q)), init) 
    q3 = TI.init_tps(TI.numtype(eltype(Q)), init) 
    q = Quaternion(q0,q1,q2,q3)
  else
    q = nothing
  end
  return q
end

function init_map_s(::Type{S}, init::AbstractTPSAInit, nv::Integer) where {S<:Union{Nothing,AbstractMatrix}}
  if S != Nothing
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
function init_map_x0(a::Type{X0}, ::AbstractTPSAInit, nv::Integer) where {X0<:StaticVector}
  x0 = StaticArrays.sacollect(X0, 0 for i in 1:length(a))
  return x0
end

function init_map_x(::Type{X}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:StaticVector}
  # reuse parameters if applicable
  if reuse isa TaylorMap && eltype(X) == eltype(reuse.x) && init == getinit(reuse)
    x = StaticArrays.sacollect(X, (i <= nv ? TI.init_tps(TI.numtype(eltype(X)), init) :  reuse.x[i]) for i in 1:length(X))
  else # allocate
    x = StaticArrays.sacollect(X, (i <= nv ? 
                                    TI.init_tps(TI.numtype(eltype(X)), init) : 
                                    (t = TI.init_tps(TI.numtype(eltype(X)), init); TI.seti!(t, 1, i); t)) for i in 1:length(X))
  end
  return x
end

function init_map_s(::Type{S}, ::AbstractTPSAInit, nv::Integer) where {S<:StaticMatrix}
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
  out_Q = Q == Nothing ? Nothing : similar_eltype(Q, promote_type(eltype(Q), G))
  out_S = S == Nothing ? Nothing : similar_eltype(S, promote_type(eltype(S), G))
  return $t{out_X0,out_X,out_Q,out_S}
end

function promote_rule(::Type{$t{X01,X1,Q1,S1}}, ::Type{$t{X02,X2,Q2,S2}}) where {X01,X02,X1,X2,Q1,Q2,S1,S2}
  out_X0 = similar_eltype(X01, promote_type(eltype(X01), eltype(X02)))
  out_X = similar_eltype(X1, promote_type(eltype(X1), eltype(X2)))
  !xor(Q1==Nothing, Q2==Nothing) || error("Cannot promote $(t)s: one includes spin while the other does not")
  out_Q = Q1 == Nothing ? Nothing : similar_eltype(Q1, promote_type(eltype(Q1),eltype(Q2)))
  !xor(S1==Nothing, S2==Nothing) || error("Cannot promote $(t)s: one includes stochastic kicks while the other does not")
  out_S = S1 == Nothing ? Nothing : similar_eltype(S1, promote_type(eltype(S1),eltype(S2)))
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
  out_Q = Q == Nothing ? Nothing : similar_eltype(Q, real(eltype(Q)))
  out_S = S == Nothing ? Nothing : similar_eltype(S, real(eltype(S)))
  return $t{out_X0,out_X,out_Q,out_S}
end

end
end


# =================================================================================== #
# Constructors
for t = (:DAMap, :TPSAMap)
@eval begin

# Lowest-level, internal
function _zero(::Type{$t{X0,X,Q,S}}, init::AbstractTPSAInit, nv::Integer, reuse::Union{TaylorMap,Nothing}=nothing) where {X0,X,Q,S}
  out_x0 = init_map_x0(X0, init, nv)
  out_x = init_map_x(X, init, nv, reuse)
  out_q = init_map_q(Q, init)
  out_s = init_map_s(S, init, nv)
  out_m = $t(out_x0, out_x, out_q, out_s)
  return out_m
end

#=
function zero(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  return _zero($t{X0,X,Q,S}, getinit(eltype(X)), nothing)
end

function one(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  out_m = _zero($t{X0,X,Q,S}, getinit(eltype(X)), nothing)
  nv = nvars(out_m)

  for i in 1:nv
    TI.seti!(out_m.x[i], 1, i)
  end

  if !isnothing(m.q)
    TI.seti!(out_m.q.q0, 1, 0)
  end
  return out_m
end
=#
# Explicit type specification
# Init change would be static (in type)
function $t{X0,X,Q,S}(m::TaylorMap) where {X0,X,Q,S}
  out_m = zero($t{X0,X,Q,S}, m) 
  copy!(out_m, m)
  return out_m
end

# zero but new type potentially
function zero(::Type{$t{X0,X,Q,S}}, m::TaylorMap) where {X0,X,Q,S}
  return _zero($t{X0,X,Q,S}, getinit(eltype(X)), nvars(m), m)
end

zero(::Type{$t}, m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S} = zero($t{X0,X,Q,S}, m)

# Copy ctor including optional TPSA init change
function $t(m::TaylorMap; init::AbstractTPSAInit=getinit(m))
  X0 = typeof(m.x0)
  X = similar_eltype(typeof(m.x), TI.init_tps_type(eltype(X0), init))
  Q = isnothing(m.q) ? Nothing : Quaternion{TI.init_tps_type(eltype(X0), init)}
  S = typeof(m.s)
  out_m = _zero($t{X0,X,Q,S}, init, nvars(m), m)
  copy!(out_m, m)
  return out_m
end

# Kwarg ctor:
function $t(;
  init::Union{AbstractTPSAInit,Nothing}=nothing,
  nv::Union{Integer,Nothing}=nothing,
  np::Union{Integer,Nothing}=nothing,
  x0::Union{AbstractVector,Nothing}=nothing,
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
  s::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
  stochastic::Union{Bool,Nothing}=nothing,
) 
  if isnothing(init)
    if !isnothing(x) && TI.is_tps_type(eltype(x)) isa TI.IsTPSType
      init = getinit(first(x))
    elseif !isnothing(q) && TI.is_tps_type(eltype(q)) isa TI.IsTPSType
      init = getinit(first(q))
    else
      error("No TPSA definition has been provided, nor is one inferrable from the input arguments")
    end
  end

  # Try to infer nv if not provided
  if isnothing(nv)
    if !isnothing(x0)
      nv = length(x0)
    elseif !isnothing(x)
      nv = length(x)
    elseif x_matrix isa AbstractMatrix
      nv = size(x_matrix, 1)
    elseif !isnothing(s)
      nv = size(s, 1)
    else
      nv = DEFAULT_NVARS
    end
  end

  if isnothing(np)
    np = 0
  end
  
  nn = nv+np
  # Check if nv+np agrees with ndiffs in init
  nn == ndiffs(init) || error("Number of variables + parameters does not agree with the number of differentials in the TPSA")

  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? TI.numtype(eltype(t)) : Float64), (x0, x, q, s))...)
  TW = TI.init_tps_type(W, init)
  X0 = isnothing(x0) ? @_DEFAULT_X0(nv){W} : similar_eltype(typeof(x0), W)
  X = isnothing(x) ? @_DEFAULT_X(nn){TW} : similar_eltype(typeof(x), TW)
  
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
    S = isnothing(s) ? @_DEFAULT_S(nv){W} : similar_eltype(typeof(s), W)
  else
    S = Nothing
  end

  # Construct map:
  out_m = _zero($t{X0,X,Q,S}, init, nv)

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
  return _zero(typeof(m), getinit(m), nvars(m), m)
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

# --- internal composer used by both TPSAMap and DAMap ---
function _compose!(
  m::TaylorMap, 
  m2::TaylorMap, 
  m1::TaylorMap,
  work_q::Union{Nothing,Quaternion}, 
  do_spin::Bool, 
  do_stochastic::Bool
)
  m.x0 .= m1.x0
  nv = nvars(m)

  TI.compose!(view(m.x, 1:nv), view(m2.x, 1:nv), m1.x)

  # Spin:
  if !isnothing(m.q)
    if do_spin
      # TO-DO: use MQuaternion (mutable quaternion) so only 1 vectorized compose! call
      TI.compose!(work_q.q0, m2.q.q0, m1.x)
      TI.compose!(work_q.q1, m2.q.q1, m1.x)
      TI.compose!(work_q.q2, m2.q.q2, m1.x)
      TI.compose!(work_q.q3, m2.q.q3, m1.x)
      mul!(m.q, work_q, m1.q) 
    else
      TI.clear!(m.q.q0)
      TI.clear!(m.q.q1)
      TI.clear!(m.q.q2)
      TI.clear!(m.q.q3)
    end
  end 

  # Stochastic (fast with StaticArrays)
  if !isnothing(m.s)
    if do_stochastic
      M2 = jacobian(m2)   
      m.s .= M2*m1.s*transpose(M2) + m2.s
    else
      m.s .= 0
    end
  end

  return m
end

# --- DAMap ---
function compose!(
  m::DAMap, 
  m2::DAMap, 
  m1::DAMap; 
  work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q),
  do_spin::Bool=true, 
  do_stochastic::Bool=true,
  keep_scalar::Bool=true, 
)
  checkinplace(m, m2, m1)
  !(m === m1) || error("Cannot compose!(m, m2, m1) with m === m1")
  !(m === m2) || error("Cannot compose!(m, m2, m1) with m === m2")

  if keep_scalar
    m1_out = getscalar(m1)
  end

  setscalar!(m1, 0)

  _compose!(m, m2, m1, work_q, do_spin, do_stochastic)

  # Put back the m1_out and if m1 === m2, also add to outx
  if keep_scalar
    setscalar!(m1, m1_out)
    if m1 === m2
      setscalar!(m, m1_out, scl0=1)
    end
  end

  return m
end

# ---TPSAMap ---
function compose!(
  m::TPSAMap, 
  m2::TPSAMap, 
  m1::TPSAMap; 
  work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q),
  do_spin::Bool=true, 
  do_stochastic::Bool=true,
  keep_scalar=nothing # For TPSAMap, this is just a placeholder, scalar part is never dropped
)
  isnothing(keep_scalar) || error("The keep_scalar kwarg is invalid for TPSAMap, because the scalar part is never dropped")
  checkinplace(m, m2, m1)
  !(m === m1) || error("Cannot compose!(m, m2, m1) with m === m1")
  !(m === m2) || error("Cannot compose!(m, m2, m1) with m === m2")

  # TPSAMap setup:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  setscalar!(m1, m2.x0, scl0=1, scl1=-1)

  _compose!(m, m2, m1, work_q, do_spin, do_stochastic)

  # Now fix m1 and if m2 === m1, add to output too:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)

  setscalar!(m1, m2.x0, scl0=1)
  if m1 === m2
    setscalar!(m, m2.x0, scl0=1)
  end

  return m
end

# =================================================================================== #
# Inversion
function inv!(
  m::TaylorMap,
  m1::TaylorMap; 
  do_spin::Bool=true, 
  work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q)
)

  checkinplace(m, m1)
  !(m === m1) || error("Cannot inv!(m, m1) with m === m1")

  TI.inv!(m.x, m1.x)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  if !isnothing(m.q)
    if do_spin
      inv!(work_q, m1.q)
      compose!(m.q, work_q, m.x)
      # TO-DO: use MQuaternion (mutable quaternion) so only 1 vectorized compose! call
      TI.compose!(m.q.q0, work_q.q.q0, m.x)
      TI.compose!(m.q.q1, work_q.q.q1, m.x)
      TI.compose!(m.q.q2, work_q.q.q2, m.x)
      TI.compose!(m.q.q3, work_q.q.q3, m.x)
    else
      TI.clear!(m.q.q0)
      TI.clear!(m.q.q1)
      TI.clear!(m.q.q2)
      TI.clear!(m.q.q3)
    end
  end

  m.x0 .= getscalar(m1)

  setscalar!(m, m1.x0)

  return m
end

# =================================================================================== #
# Map powers/pow!

# --- internal power used by both TPSAMap and DAMap ---
function _pow!(
  m::TaylorMap, 
  m1::TaylorMap, 
  n::Integer,
  work_m::Union{Nothing,TaylorMap},
  work_q::Union{Nothing,Quaternion},
  keep_scalar::Union{Bool,Nothing}
)
  checkinplace(m, m1)
  !(m === m1) || error("Cannot pow!(m, m1) with m === m1")
  if n == 0
    clear!(m)
    m.x0 .= m1.x0
    setray!(m, x_matrix=I)
    if !isnothing(m.q)
      setquat!(m, q=I)
    end
    return m
  end

  pown = abs(n)-1
  tmp1 = m
  tmp2 = work_m
  if n > 0
    if isodd(pown) # e.g. n = 2 corresponds to pown = 1
      copy!(tmp2, m1)
      for _ in 1:pown
        compose!(tmp1, tmp2, m1, keep_scalar=keep_scalar, work_q=work_q)
        tmp2, tmp1 = tmp1, tmp2
      end
    else  # e.g. n = 3 corresponds to pown = 2
      copy!(tmp1, m1)
      for _ in 1:pown
        compose!(tmp2, tmp1, m1, keep_scalar=keep_scalar, work_q=work_q)
        tmp2, tmp1 = tmp1, tmp2
      end
    end
  else
    # If n < 0 then we need one more step to invert the map, so opposite of above
    if isodd(pown) 
      copy!(tmp1, m1)
      for _ in 1:pown
        compose!(tmp2, tmp1, m1, keep_scalar=keep_scalar, work_q=work_q)
        tmp2, tmp1 = tmp1, tmp2
      end
    else
      copy!(tmp2, m1)
      for _ in 1:pown
        compose!(tmp1, tmp2, m1, keep_scalar=keep_scalar, work_q=work_q)
        tmp2, tmp1 = tmp1, tmp2
      end
    end
    inv!(m, work_m, work_q=work_q)
  end
  return m
end


function pow!(
  m::DAMap, 
  m1::DAMap, 
  n::Integer; 
  work_m::Union{DAMap,Nothing}=zero(m), 
  work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q)
)
  checkinplace(m, m1, work_m)
  !(m === m1) || error("Cannot pow!(m, m1) with m === m1")
  !(m === work_m) || error("Cannot pow!(m, m1; work_m=work_m) with m === work_m")
  m1_out = getscalar(m1)
  _pow!(m, m1, n, work_m, work_q, false)
  setscalar!(m1, m1_out)
  return m
end

function pow!(  
  m::TPSAMap, 
  m1::TPSAMap, 
  n::Integer; 
  work_m::Union{TPSAMap,Nothing}=zero(m), 
  work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q)
)
  checkinplace(m, m1, work_m)
  !(m === m1) || error("Cannot pow!(m, m1) with m === m1")
  !(m === work_m) || error("Cannot pow!(m, m1; work_m=work_m) with m === work_m")
  _pow!(m, m1, n, work_m, work_q, nothing)
  return m
end

# =================================================================================== #
# Composition and inversion out-of-place operators
for t = (:DAMap, :TPSAMap)
@eval begin

function ∘(m2::$t, m1::$t)
  m2prom, m1prom = promote(m2, m1)
  m = zero(m2prom)
  compose!(m, m2prom, m1prom)
  return m
end

# When composing a TPS scalar/vector function w a map, use orbital part of map:
function ∘(m2, m1::$t)
  TI.is_tps_type(eltype(m2)) isa TI.IsTPSType || error("Cannot compose: $(eltype(m2)) is not a TPS type supported by TPSAInterface.jl")
  T = promote_type(eltype(m1.x), eltype(m2))
  T == eltype(m1.x) ? m1xprom = m1.x : m1xprom = T.(m1.x)
  T == eltype(m2) ? m2prom = m2 : m2prom = T.(m2)
  m = zero(m2prom)
  TI.compose!(m, m2prom, m1xprom)
  return m
end

# necessary because inv(m)^n != inv(m^n) but Julia defaults to the first
literal_pow(::typeof(^), m::$t{X0,X,Q,S}, vn::Val{n}) where {X0,X,Q,S,n} = ^(m, n) 

literal_pow(::typeof(^), m::$t{X0,X,Q,S}, vn::Val{1}) where {X0,X,Q,S} = copy(m)
literal_pow(::typeof(^), m::$t{X0,X,Q,S}, vn::Val{0}) where {X0,X,Q,S} = one(m)
literal_pow(::typeof(^), m::$t{X0,X,Q,S}, vn::Val{-1}) where {X0,X,Q,S} = inv(m)
inv(m::$t; do_spin::Bool=true) = (out_m = zero(m); inv!(out_m, m, do_spin=do_spin); return out_m)
^(m::$t, n::Integer) = (out_m = zero(m); pow!(out_m, m, n); return out_m)

# Also allow * for simplicity and \ and / 
*(m2, m1::$t) = ∘(m2, m1)
/(m2::$t, m1::$t) = m2 ∘ inv(m1) 
\(m2::$t, m1::$t) = inv(m2) ∘ m1

# Uniform scaling for * (∘) and /, \
∘(m::$t, J::UniformScaling) = $t(m)
∘(J::UniformScaling, m::$t) = $t(m)
*(m::$t, J::UniformScaling) = $t(m)
*(J::UniformScaling, m::$t) = $t(m)
/(m::$t, J::UniformScaling) = $t(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = $t(m)

#compose!(m::$t, m1::$t, J::UniformScaling; kwargs...) = copy!(m, m1)
#compose!(m::$t, J::UniformScaling, m1::$t; kwargs...) = copy!(m, m1)

end
end