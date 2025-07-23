# Field initialization functions.
# Note that even when a TPSA definition is inferrable from the types (V0, V, ...),
# we still explicitly pass the init in case it is not. This is redundant, however ensures 
# flexibility for different TPSA packages/modes.

# These may be overrided by external array packages.
#---
# StaticArray field initialization functions.

# Consistency checks are made by the `checkmapsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA



#=

This file defines the types abstract type TaylorMap, and concrete types DAMap 
and TPSAMap (which differ only in concatenation and inversion rules), 
promotion rules, constructors, and map-specific operators.

=#

# =================================================================================== #
# Types

"""
    abstract type TaylorMap{V0,V,Q,S}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. `DAMap`s have a coordinate 
system chosen so that the expansion point is always around zero, e.g. Δv = v, while `TPSAMap`s 
can have any coordinate system/expansion point, but therefore will accrue truncation error when 
trying to compose two `TPSAMap`s with differing expansion points. For normal form analysis of 
periodic maps, using a `DAMap` ensures no truncation error up to the chosen truncation order.

Any truncated power series (TPS) type supported by `TPSAInterface.jl` is allowed for use in 
a `TaylorMap`. Henceforth we will generically refer to this type as a `TPS`

## Fields
- `v0::V0` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `v::V`   -- Orbital ray as truncated power series, expansion around `v0`, with scalar part equal to EXIT coordinates of map
- `q::Q`   -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `s::S`   -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

## Type Requirements
- `V0 <: AbstractVector{<:Number}` where `ismutabletype(V0) == true` 
- `V <: AbstractVector{<:TPS}` where `TPSAInterface.numtype(V) == eltype(V0)`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(V)`
- `S <: Union{AbstractMatrix{<:Number},Nothing}` where `S != Nothing` then `eltype(S) == TPSAInterface.numtype(eltype(V))` AND `ismutabletype(S) == true`

Because the TPS type is `mutable` and `TPSAInterface.jl` provides in-place functions for modifying 
TPSs, at the lowest level, all operations on `TaylorMap`s are in-place for performance. Therefore, 
the `v` and `q` arrays which contain TPSs may be `immutable`, e.g. the orbital ray `v` may be an 
`SVector` from the `StaticArrays.jl` package, and the `Quaternion` type which is taken from 
`ReferenceFrameRotations.jl` is already `immutable`. The default for the orbital ray is `SVector`.
The `v0` and `s` arrays contain `immutable` number types, and so these arrays MUST be `mutable`.
"""
abstract type TaylorMap{V0<:AbstractVector,V<:AbstractVector,Q<:Union{Quaternion,Nothing},S<:Union{AbstractMatrix,Nothing}} end 

@inline function checkmapsanity(m::TaylorMap{V0,V,Q,S}) where {V0,V,Q,S}
  # Static checks:
  ismutabletype(V0) || error("Reference orbit array must be mutable: $V0 is immutable.")
  TI.is_tps_type(eltype(V)) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")
  eltype(V0) == TI.numtype(eltype(V)) || error("Reference orbit number type $(eltype(V0)) must be $(TI.numtype(eltype(V))) (equal to the type of the orbital ray scalar part)")
  Q == Nothing || eltype(Q) == eltype(V) || error("Quaternion number type $(eltype(Q)) must be $(eltype(V)) (equal to orbital ray)")
  S == Nothing || ismutabletype(S) || error("Stochastic envelope matrix must be mutable: $S is immutable.")
  S == Nothing || eltype(S) == TI.numtype(eltype(V)) || error("Stochastic envelope matrix number type $(eltype(S)) must be $(TI.numtype(eltype(V))) (equal to the type of the orbital ray scalar part)")
  
  # Runtime checks:
  ndiffs(first(m.v)) == length(m.v) || error("Orbital ray length disagrees with number of differentials in TPSA")
  iseven(length(m.v0)) || ndiffs(first(m.v)) > length(m.v0) || error("Coasting plane requires at least one differential (as a parameter) which corresponds to the energy-like canonical variable")
  Q == Nothing || getinit(first(m.q)) == getinit(first(m.v)) || error("Quaternion TPSA definition disagrees with orbital ray TPSA definition")
  S == Nothing || size(m.s) == (length(m.v0), length(m.v0)) || error("Size of stochastic matrix disagrees with number of variables in map")
end

# =================================================================================== #
# DAMap

"""
    struct DAMap{V0,V,Q,S} <: TaylorMap{V0,V,Q,S}

Map such that the expansion point is around zero. 
See the `TaylorMap` documentation for details.

## Fields
- `v0::V0` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `v::V`   -- Orbital ray as truncated power series, expansion around `v0`, with scalar part equal to EXIT coordinates of map
- `q::Q`   -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `s::S`   -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

"""
struct DAMap{V0,V,Q,S} <: TaylorMap{V0,V,Q,S}
  v0::V0    # Entrance value of map
  v::V      # Expansion around v0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function DAMap(v0, v, q, s)
    m = new{typeof(v0),typeof(v),typeof(q),typeof(s)}(v0, v, q, s)
    checkmapsanity(m)
    return m
  end
end

# =================================================================================== #
# TPSAMap

"""
    struct TPSAMap{V0,V,Q,S} <: TaylorMap{V0,V,Q,S}

Map where the expansion point does not have to be around zero. 
Includes feed down error if composing maps with different expansion points.
See the `TaylorMap` documentation for details.

## Fields
- `v0::V0` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `v::V`   -- Orbital ray as truncated power series, expansion around `v0`, with scalar part equal to EXIT coordinates of map
- `q::Q`   -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `s::S`   -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

"""
struct TPSAMap{V0,V,Q,S} <: TaylorMap{V0,V,Q,S}
  v0::V0    # Entrance value of map
  v::V      # Expansion around v0, with scalar part equal to EXIT value of map
  q::Q      # Quaternion for spin
  s::S      # Envelope for stochasticity

  function TPSAMap(v0, v, q, s)
    m = new{typeof(v0),typeof(v),typeof(q),typeof(s)}(v0, v, q, s)
    checkmapsanity(m)
    return m
  end
end

# =================================================================================== #
# init_map_x0

"""
   init_map_x0(::Type{V0}, nv::Integer) where {V0<:AbstractVector}
   init_map_x0(::Type{V0}, nv::Integer) where {V0<:StaticVector}

Returns an array of type `V0` where each element of the array is initialized to zero. 

For the static case, the `nv` argument is ignored and the length of the returned array is the same as `V0`. 
For the general case (since the length of the vector is unknown), the length is `nv`.
""" init_map_x0

function init_map_x0(::Type{V0}, nv::Integer) where {V0<:AbstractVector}
  v0 = similar(V0, nv)
  v0 .= 0
  return v0
end

function init_map_x0(::Type{V0}, nv::Integer) where {V0<:StaticVector}
  v0 = StaticArrays.sacollect(V0, 0 for i in 1:length(V0))
  return v0
end

# =================================================================================== #
# init_map_x

"""
    init_map_x(::Type{V}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {V<:AbstractVector}
    init_map_x(::Type{V}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {V<:StaticVector}

Returns an array of of type `V` with each element of the array initialized by `init`.

For the static case, the `nv` argument is ignored and the length of the returned array is the same as `V0`. 
For the general case (since the length of the vector is unknown), the length is `nv`.
""" init_map_x

function init_map_x(::Type{V}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {V<:AbstractVector}
  nn = ndiffs(init)     # Number of variables+parameters in TPSA
  v = similar(V, nn)
  for i in 1:nv
    v[i] = TI.init_tps(TI.numtype(eltype(V)), init) 
  end
  # use same parameters if reuse isa TaylorMap and eltype(v) == eltype(reuse.v)
  if reuse isa TaylorMap && eltype(v) == eltype(reuse.v) && init == getinit(reuse)
    v[nv+1:nn] .= view(reuse.v, nv+1:nn)
  else # allocate
    for i in nv+1:nn
      v[i] = TI.init_tps(TI.numtype(eltype(V)), init) 
      TI.seti!(v[i], 1, i)
    end
  end
  return v
end

function init_map_x(::Type{V}, init::AbstractTPSAInit, nv::Integer, reuse::Union{Nothing,TaylorMap}=nothing) where {V<:StaticVector}
  # reuse parameters if applicable
  if reuse isa TaylorMap && eltype(V) == eltype(reuse.v) && init == getinit(reuse)
    v = StaticArrays.sacollect(V, (i <= nv ? TI.init_tps(TI.numtype(eltype(V)), init) :  reuse.v[i]) for i in 1:length(V))
  else # allocate
    v = StaticArrays.sacollect(V, (i <= nv ? 
                                    TI.init_tps(TI.numtype(eltype(V)), init) : 
                                    (t = TI.init_tps(TI.numtype(eltype(V)), init); TI.seti!(t, 1, i); t)) for i in 1:length(V))
  end
  return v
end

# =================================================================================== #
# init_map_q

"""
    init_map_q(::Type{Q}, init::AbstractTPSAInit) where {Q<:Union{Nothing,Quaternion}}
"""
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

# =================================================================================== #
# init_map_s

"""
    init_map_s(::Type{S}, nv::Integer) where {S<:Union{Nothing,AbstractMatrix}}
    init_map_s(::Type{S}, nv::Integer) where {S<:StaticMatrix}
"""
function init_map_s(::Type{S}, nv::Integer) where {S<:Union{Nothing,AbstractMatrix}}
  if S != Nothing
    s = similar(S, nv, nv)
    s .= 0
  else
    s = nothing
  end
  return s
end

function init_map_s(::Type{S}, nv::Integer) where {S<:StaticMatrix}
  s = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  return s
end

# =================================================================================== #
# Promotion rules

for t = (:DAMap, :TPSAMap)
  @eval begin    

    function promote_rule(::Type{$t{V0,V,Q,S}}, ::Type{G}) where {V0,V,Q,S,G<:Union{Number,Complex}}
      out_X0 = similar_eltype(V0, promote_type(eltype(V0), G))
      out_X = similar_eltype(V, promote_type(eltype(V), G))
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
    function complex(type::Type{$t{V0,V,Q,S}}) where {V0,V,Q,S}
      return promote_type(type, complex(TI.numtype(eltype(V))))
    end

    # --- real type ---
    function real(::Type{$t{V0,V,Q,S}}) where {V0,V,Q,S}
      out_X0 = similar_eltype(V0, real(eltype(V0)))
      out_X = similar_eltype(V, real(eltype(V)))
      out_Q = Q == Nothing ? Nothing : similar_eltype(Q, real(eltype(Q)))
      out_S = S == Nothing ? Nothing : similar_eltype(S, real(eltype(S)))
      return $t{out_X0,out_X,out_Q,out_S}
    end

  end  # @eval
end

# =================================================================================== #
# Constructors

for t in (:DAMap, :TPSAMap)
  @eval begin

    # Lowest-level, internal
    function _zero(::Type{$t{V0,V,Q,S}}, init::AbstractTPSAInit, nv::Integer, reuse::Union{TaylorMap,Nothing}=nothing) where {V0,V,Q,S}
      out_v0 = init_map_x0(V0, nv)
      out_v = init_map_x(V, init, nv, reuse)
      out_q = init_map_q(Q, init)
      out_s = init_map_s(S, nv)
      out_m = $t(out_v0, out_v, out_q, out_s)
      return out_m
    end

    # Explicit type specification
    # Init change would be static (in type)
    function $t{V0,V,Q,S}(m::TaylorMap) where {V0,V,Q,S}
      out_m = zero($t{V0,V,Q,S}, m) 
      copy!(out_m, m)
      return out_m
    end

    # zero but new type potentially
    function zero(::Type{$t{V0,V,Q,S}}, m::TaylorMap) where {V0,V,Q,S}
      return _zero($t{V0,V,Q,S}, getinit(eltype(V)), nvars(m), m)
    end

    zero(::Type{$t}, m::TaylorMap{V0,V,Q,S}) where {V0,V,Q,S} = zero($t{V0,V,Q,S}, m)

    # Copy ctor including optional TPSA init change
    function $t(m::TaylorMap; init::AbstractTPSAInit=getinit(m))
      V0 = typeof(m.v0)
      V = similar_eltype(typeof(m.v), TI.init_tps_type(eltype(V0), init))
      Q = isnothing(m.q) ? Nothing : Quaternion{TI.init_tps_type(eltype(V0), init)}
      S = typeof(m.s)
      out_m = _zero($t{V0,V,Q,S}, init, nvars(m), m)
      copy!(out_m, m)
      return out_m
    end

    # Kwarg constructor:
    function $t(;
          init::Union{AbstractTPSAInit,Nothing}=nothing,
          nv::Union{Integer,Nothing}=nothing,
          np::Union{Integer,Nothing}=nothing,
          v0::Union{AbstractVector,Nothing}=nothing,
          v::Union{AbstractVector,Nothing}=nothing,
          v_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
          q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
          q_map::Union{AbstractMatrix,Nothing}=nothing,
          s::Union{AbstractMatrix,Nothing}=nothing,
          spin::Union{Bool,Nothing}=nothing,
          stochastic::Union{Bool,Nothing}=nothing)

      if isnothing(init)
        if !isnothing(v) && TI.is_tps_type(eltype(v)) isa TI.IsTPSType
          init = getinit(first(v))
        elseif !isnothing(q) && TI.is_tps_type(eltype(q)) isa TI.IsTPSType
          init = getinit(first(q))
        else
          error("No TPSA definition has been provided, nor is one inferrable from the input arguments")
        end
      end

      # Try to infer nv if not provided
      if isnothing(nv)
        coast = false
        if !isnothing(v0)
          nv = length(v0)
        elseif !isnothing(v)
          # Coast check
          coast = check_coast(v)
        elseif v_matrix isa AbstractMatrix
          nv = size(v_matrix, 1)
          if all(t->t ≈ 0, view(v_matrix, nv, 1:nv-1)) && v_matrix[nv,nv] ≈ 1
            coast = true
          end
        elseif !isnothing(s)
          nv = size(s, 1)
        else
          nv = DEFAULT_NVARS
        end

        if coast
          nv = length(v) - 1
        else
          nv = length(v)
        end
      else
        coast = isodd(nv)
      end

      if isnothing(np)
        np = ndiffs(init) - nv
      end

      nn = nv+np
      # Check if nv+np agrees with ndiffs in init
      nn == ndiffs(init) || error("Number of variables + parameters does not agree with the number of differentials in the TPSA")

      # Assemble types:
      W = promote_type(map(t->(!isnothing(t) ? TI.numtype(eltype(t)) : Float64), (v0, v, q, s))...)
      TW = TI.init_tps_type(W, init)
      V0 =  @_DEFAULT_X0(nv){W} #isnothing(v0) ? @_DEFAULT_X0(nv){W} : similar_eltype(typeof(v0), W)
      V = @_DEFAULT_X(nn){TW} #isnothing(v) ? @_DEFAULT_X(nn){TW} : similar_eltype(typeof(v), TW)

      if isnothing(spin)
        Q = isnothing(q) && isnothing(q_map) ? Nothing : Quaternion{TW}
      elseif spin
        Q = Quaternion{TW}
      else
        Q = Nothing
      end

      if isnothing(stochastic)
        S = isnothing(s) ? Nothing : @_DEFAULT_S(nv){W} #similar_eltype(typeof(s), W)
      elseif stochastic
        S = @_DEFAULT_S(nv){W} #isnothing(s) ? @_DEFAULT_S(nv){W} : similar_eltype(typeof(s), W)
       else
        S = Nothing
      end

      #return $t{V0,V,Q,S}, init, nv

      # Construct map:
      out_m = _zero($t{V0,V,Q,S}, init, nv)

      if !isnothing(v0)
        out_m.v0 .= v0
      end

      setray!(out_m.v, v=v, v_matrix=v_matrix)
      if !isnothing(out_m.q) && !isnothing(q)
        setquat!(out_m.q, q=q, q_map=q_map)
      end

      if !isnothing(out_m.s) && !isnothing(s)
        out_m.s .= s
      end

      return out_m
    end

  end  # @eval
end

# =================================================================================== #

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
    TI.seti!(out_m.v[i], 1, i)
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
  m.v0 .= m1.v0
  nv = nvars(m)

  TI.compose!(view(m.v, 1:nv), view(m2.v, 1:nv), m1.v)

  # Spin:
  if !isnothing(m.q)
    if do_spin
      # m.q = m2.q(m1.v)*m1.q
      TI.compose!(work_q, m2.q, m1.v)
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
  # For TPSA Map concatenation, we need to subtract w_0 (m2 v0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 v0)
  setscalar!(m1, m2.v0, scl0=1, scl1=-1)

  _compose!(m, m2, m1, work_q, do_spin, do_stochastic)

  # Now fix m1 and if m2 === m1, add to output too:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 v0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 v0)

  setscalar!(m1, m2.v0, scl0=1)
  if m1 === m2
    setscalar!(m, m2.v0, scl0=1)
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

  TI.inv!(m.v, m1.v)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  if !isnothing(m.q)
    if do_spin
      inv!(work_q, m1.q)
      TI.compose!(m.q, work_q, m.v)
    else
      TI.clear!(m.q.q0)
      TI.clear!(m.q.q1)
      TI.clear!(m.q.q2)
      TI.clear!(m.q.q3)
    end
  end

  m.v0 .= getscalar(m1)

  setscalar!(m, m1.v0)

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
  !(m === m1) || error("Cannot pow!(m, m1) with m === m1")
  if n == 0
    clear!(m)
    m.v0 .= m1.v0
    setray!(m, v_matrix=I)
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
       work_q::Union{Nothing,Quaternion}=isnothing(m.q) ? nothing : zero(m.q))
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
    T = promote_type(eltype(m1.v), eltype(m2))
    T == eltype(m1.v) ? m1xprom = m1.v : m1xprom = T.(m1.v)
    T == eltype(m2) ? m2prom = m2 : m2prom = T.(m2)
    m = zero(m2prom)
    TI.compose!(m, m2prom, m1xprom)
    return m
  end

  # necessary because inv(m)^n != inv(m^n) but Julia defaults to the first
  literal_pow(::typeof(^), m::$t{V0,V,Q,S}, vn::Val{n}) where {V0,V,Q,S,n} = ^(m, n) 

  literal_pow(::typeof(^), m::$t{V0,V,Q,S}, vn::Val{1}) where {V0,V,Q,S} = copy(m)
  literal_pow(::typeof(^), m::$t{V0,V,Q,S}, vn::Val{0}) where {V0,V,Q,S} = one(m)
  literal_pow(::typeof(^), m::$t{V0,V,Q,S}, vn::Val{-1}) where {V0,V,Q,S} = inv(m)
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