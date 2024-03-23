#=
Probe is used for tracking. As a parametric type, 
the orbital/spin part can contain either a Float64 or 
a real Truncated Power Series. 

Because of Julia's JIT + multiple dispatch, we actually
could use only one Probe both for particle tracking and 
TPS tracking (no slowdowns expected).

Parametric type as generic as possible to eventually allow 
differentiating through.
=#
struct Probe{S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}
  x0::Vector{S}   # Entrance coordinates
  x::Vector{T}    # Out coordinates
  Q::U            # Quaternion
  E::V            # Stochastic matrix
end

# Type stable constructor
# spin=false or radiation=false is an error in order to keep type stable 
# basically if spin isa Bool then spin track 
function Probe(x::Vector{T}; x0::Vector=zeros(length(x)), Q::U=nothing, E::V=nothing, spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing) where {T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}
  length(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")

  if isnothing(spin)
    Q1 = Q
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(first(x))
    else
      Q1 = Q1
    end
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing")
  end

  if isnothing(radiation)
    E1 = E
  elseif radiation
    if isnothing(E)
      E1 = zeros(eltype(x0), length(x), length(x)) 
    else
      (length(x),length(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      E1 = E
    end
  else
    error("For no radiation, please omit the radiation kwarg or set radiation=nothing")
  end

  return Probe{eltype(x0),T,typeof(Q1),typeof(E)}(x0, x, Q1, E)
end


# Copy ctor:
function Probe(p::Probe{S,T,U,V}) where {S,T,U,V}
  return Probe{S,T,U,V}(copy(p.x0), deepcopy(p.x), deepcopy(p.Q), copy(p.E))
end

==(p1::Probe, p2::Probe) = (p1.x0 == p2.x0 && p1.x == p2.x && p1.q == p2.q && p1.E == p2.E)
