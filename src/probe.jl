#=
Probe is used for tracking. As a parametric type, 
the orbital/spin part can contain either a Float64 or 
a real Truncated Power Series. 

Because of Julia's JIT + multiple dispatch, we actually
could use only one Probe both for particle tracking and 
TPS tracking (no slowdowns expected).
=#
struct Probe{T <: Union{Float64, TPS}}
  x::Vector{T}           # Out coordinates
  x0::Vector{Float64}    # Entrance coordinates
  q::Quaternion{T}       # Quaternion
  E::Matrix{Float64}     # Stochastic matrix
end

# Quasi-keyword argument constructor for easy expansion/modification
function Probe(; x::Vector{T}, x0::Vector{Float64}, q::Quaternion{T}, E::Matrix{Float64}) where T <: Union{Float64, TPS}
  length(x) == length(x0) || error("Length of coordinates vector != length of coordinate system origin vector!")
  (length(x),length(x)) == size(E) || error("Length of coordinates vector inconsistent with size of stochastic matrix!")
  Probe(x, x0, q, E)
end

# ctors require phase space ray and optionally x0 q and E as keyword argument
function Probe(x::Vector{<:Real}; x0::Vector{<:Real}=zeros(length(x)), q::Quaternion{<:Real}=Quaternion((0.,0.,0.,0.)), E::Matrix{<:Real}=zeros(length(x), length(x)))
  return Probe(x=convert(Vector{Float64}, x), x0=convert(Vector{Float64}, x0), q=Quaternion{Float64}(q.q), E=convert(Matrix{Float64}, E))
end

function Probe(x::Vector{TPS}; x0::Vector{<:Real}=zeros(numvars(first(x))), q::Quaternion{TPS}=Quaternion((zero(first(x)), zero(first(x)), zero(first(x)), zero(first(x)))), E::Matrix{<:Real}=zeros(numvars(first(x)), numvars(first(x))))
  numvars(first(x)) == length(x) || error("Length of coordinates vector inconsistent with number of variables in GTPSA!")
  Probe(x=x, x0=convert(Vector{Float64}, x0), q=Quaternion{TPS}(q.q), E=convert(Matrix{Float64}, E))
end

# Copy ctor:
function Probe(p::Probe)
  return deepcopy(p)
end



