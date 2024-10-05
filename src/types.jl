
"""
    TaylorMap{S,T<:TPS,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

All `TaylorMap`s contain `x0` and `x` as the entrance coordinates and transfer map 
as a truncated power series respectively. If spin is included, a field `Q` containing 
a `Quaternion` as a truncated power series is included, else `Q` is `nothing`. If 
stochasticity is included, a field `E` contains a matrix of the envelope for the 
fluctuations, else `E` is nothing.

### Fields
- `x0`   -- Entrance coordinates of the map, Taylor expansion point
- `x`    -- Orbital ray as a truncated power series, expansion around `x0` + scalar part equal to EXIT coordinates of map
- `Q`    -- `Quaternion` as a truncated power series if spin is included, else `nothing`
- `E`    -- Matrix of the envelope for stochastic kicks if included, else `nothing`

"""
abstract type TaylorMap{S<:Vector,T<:Vector,U<:Union{Quaternion,Nothing},V<:Union{Matrix,Nothing}} end 

"""
    DAMap{S,T,U,V} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverses as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T,U,V} <: TaylorMap{S,T,U,V}
  x0::S     # Entrance value of map
  x::T      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
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

`TaylorMap` that composes and inverses as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T,U,V} <: TaylorMap{S,T,U,V}
  x0::S     # Entrance value of map
  x::T      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U      # Quaternion for spin
  E::V      # Envelope for stochasticity

  function TPSAMap(x0, x, Q, E)
    m = new{typeof(x0),typeof(x),typeof(Q),typeof(E)}(x0, x, Q, E)
    checkmapsanity(m)
    return m
  end
end

@inline function checkmapsanity(m::TaylorMap)
  eltype(m.x0) == eltype(eltype(m.x)) || error("Reference orbit type $(eltype(m.x0)) must be $(eltype(eltype(m.x))) (equal to scalar of orbital)")
  isnothing(m.Q) || eltype(m.Q) == eltype(m.x) || error("Quaternion type $(eltype(m.Q)) must be $(eltype(m.x)) (equal to orbital)")
  isnothing(m.E)|| eltype(m.E) == eltype(eltype(m.x)) || error("Stochastic matrix type $(eltype(m.E)) must be $(eltype(eltype(m.x))) (equal to scalar of orbital)")
end

"""
    Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

Parametric type used for tracking. The orbital/spin part can contain 
either scalars or `TPS`s. If spin is included, the field `Q` is a `Quaternion` 
with the same type as `x`, else `Q` is `nothing`. If radiation is included, 
the field `E` contains the FD matrix with `eltype(E)==S`, else `E` is `nothing`.
"""
struct Probe{S,T,U<:Union{Quaternion,Nothing},V<:Union{Matrix,Nothing}}
  x0::Vector{S}   # Entrance coordinates
  x::Vector{T}    # Out coordinates
  Q::U            # Quaternion
  E::V            # Stochastic matrix
end

"""


Lie operator to act on maps. Can be turned into a map with exp(:F:)
"""
struct VectorField{T<:Vector, U<:Union{Quaternion,Nothing}}
  x::T
  Q::U           
end

const UseType = Union{Descriptor, TPS, DAMap, TPSAMap, Probe{<:Any,TPS,<:Any,<:Any}, VectorField, Nothing}

for t = (:DAMap, :TPSAMap)
@eval begin    

function promote_rule(::Type{$t{S,T,U,V}}, ::Type{G}) where {S,T,U,V,G<:Union{Number,Complex}}
  outS = Vector{promote_type(eltype(S),eltype(eltype(T)),G)}
  outT = Vector{promote_type(eltype(T),G)}
  U != Nothing ? outU = Quaternion{promote_type(eltype(U), G)} : outU = Nothing
  V != Nothing ? outV = Matrix{promote_type(G,eltype(V))} : outV = Nothing
  return $t{outS,outT,outU,outV}
end

# Currently promote_type in promotion.jl gives
# promote_type(::Type{T}, ::Type{T}) where {T} = T
# and does not even call promote_rule, therefore this is never reached
# Therefore I will required the reference orbit to have the same eltype as the 
# TPS at construction.
function promote_rule(::Type{$t{S1,T1,U1,V1}}, ::Type{$t{S2,T2,U2,V2}}) where {S1,S2,T1,T2,U1,U2,V1,V2} 
  outS = Vector{promote_type(eltype(eltype(T1)), eltype(eltype(T2)))}
  outT = Vector{promote_type(eltype(T1), eltype(T2))}
  U1 != Nothing ? outU = Quaternion{promote_type(eltype(U2),eltype(U2))} : outU = Nothing
  V1 != Nothing ? outV = promote_type(V1,V2) : outV = Nothing
  return $t{outS,outT,outU,outV}
end 

end
end



