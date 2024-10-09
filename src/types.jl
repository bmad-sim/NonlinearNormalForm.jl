#=

Defines the types used throughout the package. Specifically, the abstract 
type TaylorMap, and concrete types DAMap and TPSAMap (which differ only in 
concatenation and inversion rules). 

=#
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

# --- complex type ---
function complex(::Type{$t{S,T,U,V}}) where {S,T,U,V}
  return $t{Vector{ComplexF64},Vector{ComplexTPS64},U == Nothing ? Nothing : Quaternion{ComplexTPS64}, V == Nothing ? Nothing : Matrix{ComplexF64}}
end

# --- real type ---
function real(::Type{$t{S,T,U,V}}) where {S,T,U,V}
  return $t{Vector{Float64},Vector{TPS64},U == Nothing ? Nothing : Quaternion{TPS64}, V == Nothing ? Nothing : Matrix{Float64}}
end

end
end


"""
    VectorField{T<:Vector, U<:Union{Quaternion,Nothing}}

Lie operator to act on maps. Can be turned into a map with exp(:F:)
"""
struct VectorField{T<:Vector, U<:Union{Quaternion,Nothing}}
  x::T
  Q::U           
end

# --- complex type ---
function complex(::Type{VectorField{T,U}}) where {T,U}
  return VectorField{ComplexTPS64, U == Nothing ? Nothing : Quaternion{ComplexTPS64}}
end

# --- real type ---
function real(::Type{VectorField{T,U}}) where {T,U}
  return VectorField{Vector{TPS64},U == Nothing ? Nothing : Quaternion{TPS64}}
end

function promote_rule(::Type{VectorField{T,U}}, ::Type{G}) where {T,U,G<:Union{Number,Complex}}
  outT = Vector{promote_type(eltype(T),G)}
  U != Nothing ? outU = Quaternion{promote_type(eltype(U), G)} : outU = Nothing
  return VectorField{outT,outU}
end


const UseType = Union{Descriptor, TPS, DAMap, TPSAMap, VectorField, Nothing}







