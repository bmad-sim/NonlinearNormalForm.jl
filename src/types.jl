#=

Defines the types used throughout the package. Specifically, the abstract 
type TaylorMap, and concrete types DAMap and TPSAMap (which differ only in 
concatenation and inversion rules). 

=#

"""
    TaylorMap{S,T,U,V}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

# Fields
- `x0::S` -- Reference orbit. The entrance coordinates of the map as scalars, or equivalently the Taylor map expansion point.
- `x::T`  -- Orbital ray as truncated power series, expansion around `x0`, with scalar part equal to EXIT coordinates of map
- `Q::U`  -- `Quaternion` as truncated power series if spin is included, else `nothing`
- `E::V`  -- Matrix of the envelope for stochastic kicks as scalars if included, else `nothing`

# Type Requirements
- `S <: AbstractVector{<:Number}` where `ismutabletype(S) == true` 
- `T <: AbstractVector{<:TPS}` where `GTPSA.numtype(T) == eltype(S)`
- `U <: Union{Quaternion{<:TPS},Nothing}` and if `U != Nothing` then `eltype(U) == eltype(T)`
- `V <: Union{AbstractMatrix{<:Number},Nothing}` where `V != Nothing` then `eltype(V) == GTPSA.numtype(eltype(T))` AND `ismutabletype(V) == true`

Because the `TPS` type is `mutable` and `GTPSA.jl` provides in-place functions for modifying 
`TPS`s, at the lowest level, all operations on `TaylorMap`s are in-place for performance. Therefore, 
the `x` and `Q` arrays which contain `TPS`s may be `immutable`, e.g. the orbital ray `x` may be an 
`SVector` from the `StaticArrays.jl` package, and the `Quaternion` type which is taken from 
`ReferenceFrameRotations.jl` is already `immutable`. The default for the orbital ray is `SVector`.
The `x0` and `E` arrays contain `immutable` number types, and so these arrays MUST be `mutable`.
"""
abstract type TaylorMap{S<:AbstractVector{<:Number},T<:AbstractVector{<:TPS},U<:Union{Quaternion{<:TPS},Nothing},V<:Union{AbstractMatrix{<:Number},Nothing}} end 

@inline function checkmapsanity(m::TaylorMap{S,T,U,V}) where {S,T,U,V}
  ismutabletype(S) || error("Reference orbit array must be mutable: $S is immutable.")
  eltype(S) == GTPSA.numtype(eltype(T)) || error("Reference orbit number type $(eltype(S)) must be $(GTPSA.numtype(eltype(T))) (equal to the type of the orbital ray scalar part)")
  U == Nothing || eltype(U) == eltype(T) || error("Quaternion number type $(eltype(U)) must be $(eltype(T)) (equal to orbital ray)")
  V == Nothing || ismutabletype(V) || error("Stochastic envelope matrix must be mutable: $V is immutable.")
  V == Nothing || eltype(V) == GTPSA.numtype(eltype(T)) || error("Stochastic envelope matrix number type $(eltype(V)) must be $(GTPSA.numtype(eltype(T))) (equal to the type of the orbital ray scalar part)")
end

"""
    DAMap{S,T,U,V} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverts as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T,U,V} <: TaylorMap{S,T,U,V}
  x0::S     # Entrance value of map
  x::T      # Expansion around x0, with scalar part equal to EXIT value of map
  Q::U      # Quaternion for spin
  E::V      # Envelope for stochasticity

  function DAMap(x0, x, Q, E)
    m = new{typeof(x0),typeof(x),typeof(Q),typeof(E)}(x0, x, Q, E)
    checkmapsanity(m)
    return m
  end
end

"""
    TPSAMap{S,T,U,V} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverts as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T,U,V} <: TaylorMap{S,T,U,V}
  x0::S     # Entrance value of map
  x::T      # Expansion around x0, with scalar part equal to EXIT value of map
  Q::U      # Quaternion for spin
  E::V      # Envelope for stochasticity

  function TPSAMap(x0, x, Q, E)
    m = new{typeof(x0),typeof(x),typeof(Q),typeof(E)}(x0, x, Q, E)
    checkmapsanity(m)
    return m
  end
end

# Promotion utility functions for easy override by external array packages
# Should promote the element type and return it as a vector (x,x0) or matrix (E)
# with dimensionality specified by Val.
promote_x0_type(::Type{S}, ::Type{G}) where {S,G<:Union{Number,Complex}} = _promote_array_eltype(S, G, Val{1})
promote_x_type(::Type{T}, ::Type{G}) where {T,G<:Union{Number,Complex}} = _promote_array_eltype(T, G, Val{1})
promote_Q_type(::Type{U}, ::Type{G}) where {U,G<:Union{Number,Complex}} = U != Nothing ? Quaternion{promote_type(eltype(U), G)} : Nothing
promote_E_type(::Type{V}, ::Type{G}) where {V,G<:Union{Number,Complex}} = V != Nothing ? _promote_array_eltype(V, G, Val{2}) : Nothing
real_x0_type(::Type{S}) where {S} = _real_array_eltype(S, Val{1})
real_x_type(::Type{T}) where {T} = _real_array_eltype(T, Val{1})
real_Q_type(::Type{U}) where {U} = U != Nothing ? Quaternion{real(eltype(U))} : Nothing 
real_E_type(::Type{V}) where {V} = V != Nothing ? _real_array_eltype(V, Val{2}) : Nothing

for t = (:DAMap, :TPSAMap)
@eval begin    

function promote_rule(::Type{$t{S,T,U,V}}, ::Type{G}) where {S,T,U,V,G<:Union{Number,Complex}}
  outS = promote_x0_type(S, G)
  outT = promote_x_type(T, G)
  outU = promote_Q_type(U, G)
  outV = promote_E_type(V, G)
  return $t{outS,outT,outU,outV}
end

# --- complex type ---
function complex(type::Type{$t{S,T,U,V}}) where {S,T,U,V}
  return promote_type(type, complex(GTPSA.numtype(eltype(T))))
end

# --- real type ---
function real(::Type{$t{S,T,U,V}}) where {S,T,U,V}
  outS = real_x0_type(S)
  outT = real_x_type(T)
  outU = real_Q_type(U)
  outV = real_E_type(V)
  return $t{outS,outT,outU,outV}
end

end
end


"""
    VectorField{T,U}

Lie operator which acts on `DAMap`s e.g. dM/dt = FM where F is the `VectorField` and 
`M` is the `DAMap`. `F` can be promoted to a `DAMap` using `exp(F)`.

# Fields
- `x::T`  -- Orbital ray as truncated power series, expansion around `x0`, with scalar part equal to EXIT coordinates of map
- `Q::U`  -- `Quaternion` as truncated power series if spin is included, else `nothing`

# Type Requirements
- `T <: AbstractVector{<:TPS}`
- `U <: Union{Quaternion{<:TPS},Nothing}` and if `U != Nothing` then `eltype(U) == eltype(T)`
"""
struct VectorField{T<:AbstractVector{<:TPS},U<:Union{Quaternion{<:TPS},Nothing}}
  x::T
  Q::U     
  function VectorField(x, Q)
    m = new{typeof(x),typeof(Q)}(x, Q)
    typeof(Q) == Nothing || eltype(Q) == eltype(x) || error("Quaternion number type $(eltype(Q)) must be $(eltype(x)) (equal to orbital ray)")
    return m
  end
end


function promote_rule(::Type{VectorField{T,U}}, ::Type{G}) where {T,U,G<:Union{Number,Complex}}
  outT = promote_x_type(T, G)
  outU = promote_Q_type(U, G)
  return VectorField{outT,outU}
end

# --- complex type ---
function complex(type::Type{VectorField{T,U}}) where {T,U}
  return promote_type(type, complex(GTPSA.numtype(eltype(T))))
end

# --- real type ---
function real(::Type{VectorField{T,U}}) where {T,U}
  outT = real_x_type(T)
  outU = real_Q_type(U)
  return VectorField{outT,outU}
end

const UseType = Union{Descriptor, TPS, DAMap, TPSAMap, VectorField}







