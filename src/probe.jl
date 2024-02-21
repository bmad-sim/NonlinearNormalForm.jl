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
struct Probe{S<:Real, T<:Real, U<:Real, V<:Real}
  x0::Vector{S}          # Entrance coordinates
  x::Vector{T}           # Out coordinates
  q::Quaternion{U}       # Quaternion
  E::Matrix{V}           # Stochastic matrix

  function Probe{S,T,U,V}(x0, x, q, E) where {S<:Real, T<:Real, U<:Real, V<:Real}
    length(x) == length(x0) || error("Length of orbital ray != length of coordinate system origin vector!")
    (length(x),length(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
    return new(x0, x, q, E)
  end
end

# ctors require phase space ray and optionally x0 q and E as keyword argument
function Probe(x::Vector{T}; x0::Vector{S}=zeros(length(x)), q::Quaternion{U}=Quaternion(1.,0.,0.,0.), E::Matrix{V}=zeros(length(x), length(x))) where {S<:Real, T<:Real, U<:Real, V<:Real}
  return Probe{S,T,U,V}(x0, x, q, E)
end

# For now orbit TPS always -> spin TPS
function Probe(x::Vector{T}; x0::Vector{S}=zeros(numvars(x)), q::Quaternion{U}=unit_quat(first(x)), E::Matrix{V}=zeros(numvars(x), numvars(x))) where {S<:Real, T<:TPS, U<:TPS, V<:Real}
  numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
  getdesc(x) == getdesc(q) || error("Orbital ray descriptor inconsistent with quaternion descriptor!")
  return Probe{S,T,U,V}(x0, x, q, E)
end

# Copy ctor:
function Probe(p::Probe)
  return deepcopy(p)
end

==(p1::Probe, p2::Probe) = (p1.x0 == p2.x0 && p1.x == p2.x && p1.q == p2.q && p1.E == p2.E)
