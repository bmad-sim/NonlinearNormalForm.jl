#=
Probe is used for tracking. As a parametric type, 
the orbital/spin part can contain either a Float64 or 
a real Truncated Power Series. 

Because of Julia's JIT + multiple dispatch, we actually
could use only one Probe both for particle tracking and 
TPS tracking (no slowdowns expected).
=#
struct Probe{T <: Union{Float64, TPS}}
  x0::Vector{Float64}   
  v::Vector{T} 
  q::Quaternion{T}
  E::Matrix{Float64}
end

# ctors accept only phase space ray for now
function Probe(vec::Vector{<:Real}; q::Quaternion{<:Real}=Quaternion((0.,0.,0.,0.)), E::Matrix{<:Real}=zeros(length(vec), length(vec)))
  return Probe(convert(Vector{Float64}, vec), zeros(length(vec)), Quaternion{Float64}(q.q), convert(Matrix{Float64}, E))
end

function Probe(vec::Vector{TPS}; q::Quaternion{TPS}=Quaternion((zero(first(vec)), zero(first(vec)), zero(first(vec)), zero(first(vec)))), E::Matrix{<:Real}=zeros(numvars(first(vec)), numvars(first(vec))))
  x0 = map(x->x[0], vec)
  for t in vec t[0] = 0 end
  Probe(x0, vec, Quaternion{TPS}(q.q), convert(Matrix{Float64}, E))
end

# Copy ctor:
function Probe(p::Probe)
  return deepcopy(p)
end



