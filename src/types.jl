
"""
    TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

All `TaylorMap`s contain `x0` and `x` as the entrance coordinates and transfer map 
as a truncated power series respectively. If spin is included, a field `Q` containing 
a `Quaternion` as a truncated power series is included, else `Q` is `nothing`. If 
radiation is included, a field `E` contains a matrix of the envelope for stochastic 
radiation, else `E` is nothing.

### Fields
- `x0` -- Entrance coordinates of the map, Taylor expansion point
- `x`  -- Orbital ray as a truncated power series, expansion around `x0` + scalar part equal to EXIT coordinates of map
- `Q`  -- `Quaternion` as a truncated power series if spin is included, else `nothing`
- `E`  -- Matrix of the envelope for stochastic radiation if included, else `nothing`
"""
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} end 

"""
    DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverses as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

"""
    TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverses as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}  <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

"""
    Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}

Parametric type used for tracking. The orbital/spin part can contain 
either scalars or `TPS`s. If spin is included, the field `Q` is a `Quaternion` 
with the same type as `x`, else `Q` is `nothing`. If radiation is included, 
the field `E` contains the stochastic matrix with `eltype(E)==S`, else `E` is `nothing`.
"""
struct Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}
  x0::Vector{S}   # Entrance coordinates
  x::Vector{T}    # Out coordinates
  Q::U            # Quaternion
  E::V            # Stochastic matrix
end

"""


Lie operator to act on maps. Can be turned into a map with exp(:F:)
"""
struct VectorField{T<:Union{TPS,ComplexTPS}, U<:Union{Quaternion{T},Nothing}}
  x::Vector{T}  
  Q::U           
end

UseType = Union{Descriptor, TPS, ComplexTPS, TaylorMap, Probe{<:Any,Union{TPS,ComplexTPS},<:Any,<:Any}, VectorField, Nothing}