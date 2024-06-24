
"""
    TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

All `TaylorMap`s contain `x0` and `x` as the entrance coordinates and transfer map 
as a truncated power series respectively. If spin is included, a field `Q` containing 
a `Quaternion` as a truncated power series is included, else `Q` is `nothing`. If 
stochasticity is included, a field `E` contains a matrix of the envelope for the FD
kicks, else `E` is nothing.

If all planes are exhibiting pseudo-harmonic oscillations, then `idpt` is `nothing`. 
If the last plane is coasting, then `idpt` specifies which variable in the last plane 
is constant (energy-like): `idpt=false` if the variable with index `NV-1` is constant, or 
`idpt=true` is the variable with index `NV` is constant.

### Fields
- `x0`   -- Entrance coordinates of the map, Taylor expansion point
- `x`    -- Orbital ray as a truncated power series, expansion around `x0` + scalar part equal to EXIT coordinates of map
- `Q`    -- `Quaternion` as a truncated power series if spin is included, else `nothing`
- `E`    -- Matrix of the envelope for FD kicks if included, else `nothing`
- `idpt` -- If the last plane is coasting, then `idpt=false` if the first variable in last plane is constant (energy-like) or `true` if the second. If no coasting, then `idpt=nothing`
"""
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}} end 

"""
    DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}

`TaylorMap` that composes and inverses as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochasticity
  idpt::W          # Specifies index of constant (energy-like) variable
end

"""
    TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}

`TaylorMap` that composes and inverses as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}} <: TaylorMap{S,T,U,V,W}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochasticity
  idpt::W          # Specifies index of constant (energy-like) variable
end

"""
    Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

Parametric type used for tracking. The orbital/spin part can contain 
either scalars or `TPS`s. If spin is included, the field `Q` is a `Quaternion` 
with the same type as `x`, else `Q` is `nothing`. If radiation is included, 
the field `E` contains the FD matrix with `eltype(E)==S`, else `E` is `nothing`.

If all planes are exhibiting pseudo-harmonic oscillations, then `idpt` is `nothing`. 
If the last plane is coasting, then `idpt` specifies which variable in the last plane 
is constant (energy-like): `idpt=false` if the variable with index `NV-1` is constant, or 
`idpt=true` is the variable with index `NV` is constant.
"""
struct Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}
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



#UseType = Union{Descriptor, TPS, ComplexTPS, DAMap, TPSAMap, Probe{<:Any,Union{TPS,ComplexTPS},<:Any,<:Any}, VectorField, Nothing}

for t = (:DAMap, :TPSAMap)
@eval begin    

function promote_rule(::Type{$t{S,T,U,V,W}}, ::Type{G}) where {S,T,U,V,W,G<:Union{Number,Complex}}
  outS = promote_type(S,numtype(T),G)
  outT = promote_type(T,G)
  U != Nothing ? outU = Quaternion{promote_type(eltype(U), G)} : outU = Nothing
  V != Nothing ? outV = promote_type(Matrix{G},V) : outV = Nothing
  return $t{outS,outT,outU,outV,W}
end

function promote_rule(::Type{$t{S1,T1,U1,V1,W}}, ::Type{$t{S2,T2,U2,V2,W}}) where {S1,S2,T1,T2,U1,U2,V1,V2,W} 
  outS = promote_type(S1, S2, numtype(T1), numtype(T2))
  outT = promote_type(T1, T2)
  U1 != Nothing ? outU = Quaternion{promote_type(eltype(U2),eltype(U2))} : outU = Nothing
  V1 != Nothing ? outV = promote_type(V1,V2) : outV = Nothing
  return $t{outS,outT,outU,outV,W}
end 

end
end



