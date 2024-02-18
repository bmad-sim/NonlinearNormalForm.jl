#=
Probe is used for tracking. As a parametric type, 
the orbital/spin part can contain either a Float64 or 
a real Truncated Power Series. 
=#
struct Probe{T <: Union{Float64, TPS}}
  v::Vector{T} 
  q::Quaternion{T}
  x0::Vector{Float64}
  E::Matrix{Float64}
end

function Probe(v::Vector{<:Real}; q::Quaternion{<:Real}=Quaternion((0.,0.,0.,0.)), x0::Vector{<:Real}=zeros(length(v)), E::Matrix{<:Real}=zeros(length(v), length(v)))
  Probe(convert(Vector{Float64}, v), Quaternion(convert(NTuple{4,Float64}, q.q)), convert(Vector{Float64},x0), convert(Matrix{Float64},E))
end

function Probe(v::Vector{TPS}; q::Quaternion{TPS}=Quaternion((zero(first(v)), zero(first(v)), zero(first(v)), zero(first(v)))), x0::Vector{<:Real}=zeros(length(v)), E::Matrix{<:Real}=zeros(length(v), length(v)))
  numvars(first(v)) == length(v) || error("Orbital GTPSA ray length != number of variables in GTPSA!")
  length(v) == length(x0) || error("Length of initial coordinates vector != length of orbital ray!")
  Probe(v, q, convert(Vector{Float64},x0), convert(Matrix{Float64},E))
end

# Copy ctor:
function Probe(p::Probe)
  return deepcopy(p)
end



