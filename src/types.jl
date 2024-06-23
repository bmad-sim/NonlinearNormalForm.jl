
"""
    TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

All `TaylorMap`s contain `x0` and `x` as the entrance coordinates and transfer map 
as a truncated power series respectively. If spin is included, a field `Q` containing 
a `Quaternion` as a truncated power series is included, else `Q` is `nothing`. If 
stochasticity is included, a field `E` contains a matrix of the envelope for the stochastic
kicks, else `E` is nothing.

If all planes are exhibiting pseudo-harmonic oscillations, then `idpt` is `nothing`. 
If the last plane is coasting, then `idpt` specifies which variable in the last plane 
is constant (energy-like): `idpt=false` if the variable with index `NV-1` is constant, or 
`idpt=true` is the variable with index `NV` is constant.

### Fields
- `x0`   -- Entrance coordinates of the map, Taylor expansion point
- `x`    -- Orbital ray as a truncated power series, expansion around `x0` + scalar part equal to EXIT coordinates of map
- `Q`    -- `Quaternion` as a truncated power series if spin is included, else `nothing`
- `E`    -- Matrix of the envelope for stochastic kicks if included, else `nothing`
- `idpt` -- If the last plane is coasting, then `idpt=false` if the first variable in last plane is constant (energy-like) or `true` if the second. If no coasting, then `idpt=nothing`
"""
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}} end 

"""
    DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}

`TaylorMap` that composes and inverses as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochasticity
  idpt::W          # Specifies index of constant (energy-like) variable
end

"""
    TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}

`TaylorMap` that composes and inverses as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochasticity
  idpt::W          # Specifies index of constant (energy-like) variable
end

"""
    Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}}

Parametric type used for tracking. The orbital/spin part can contain 
either scalars or `TPS`s. If spin is included, the field `Q` is a `Quaternion` 
with the same type as `x`, else `Q` is `nothing`. If radiation is included, 
the field `E` contains the stochastic matrix with `eltype(E)==S`, else `E` is `nothing`.

If all planes are exhibiting pseudo-harmonic oscillations, then `idpt` is `nothing`. 
If the last plane is coasting, then `idpt` specifies which variable in the last plane 
is constant (energy-like): `idpt=false` if the variable with index `NV-1` is constant, or 
`idpt=true` is the variable with index `NV` is constant.
"""
struct Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing},W<:Union{Nothing,Bool}}
  x0::Vector{S}   # Entrance coordinates
  x::Vector{T}    # Out coordinates
  Q::U            # Quaternion
  E::V            # Stochastic matrix
  idpt::W         # Specifies index of constant (energy-like) variable
end

"""


Lie operator to act on maps. Can be turned into a map with exp(:F:)
"""
struct VectorField{T<:Union{TPS,ComplexTPS}, U<:Union{Quaternion{T},Nothing}}
  x::Vector{T}  
  Q::U           
end



UseType = Union{Descriptor, TPS, ComplexTPS, TaylorMap, Probe{<:Any,Union{TPS,ComplexTPS},<:Any,<:Any}, VectorField, Nothing}
#=
function promote_rule(m1::Type{TaylorMap{S1,T1,U,V1,W}}, m2::Type{TaylorMap{S2,T2,U,V2,W}}) where {S1,S2,T1,T2,U,V1,V2,W}
  S = promote_rule(eltype(m1.x0), eltype(m2.x0))
  T = promote_rule(eltype(m1.x), eltype(m2.x))
  if !isnothing
end =#

