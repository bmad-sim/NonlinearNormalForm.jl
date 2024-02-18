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

# ctors accept phase space ray 
function Probe(vec::Vector{<:Real}; q::Quaternion{<:Real}=Quaternion((1.,0.,0.,0.)), E::Matrix{<:Real}=zeros(length(vec), length(vec)))
  return Probe(convert(Vector{Float64}, vec), zeros(length(vec)), Quaternion{Float64}(q.q), convert(Matrix{Float64}, E))
end

function Probe(vec::Vector{TPS}; q::Quaternion{TPS}=Quaternion((one(first(vec)), zero(first(vec)), zero(first(vec)), zero(first(vec)))), E::Matrix{<:Real}=zeros(numvars(first(vec)), numvars(first(vec))))
  x0 = map(x->x[0], vec)
  for t in vec t[0] = 0 end
  Probe(x0, vec, Quaternion{TPS}(q.q), convert(Matrix{Float64}, E))
end

#=
# kwarg constructors with correct default values:
# This ctor used for scalar Probe
function Probe(;x0::Vector{<:Real}, v::Vector{<:Real}=zeros(length(v)), q::Quaternion{<:Real}=Quaternion((1.,0.,0.,0.)), E::Matrix{<:Real}=zeros(length(x0), length(x0)))
  length(x0) == length(v) || error("Length of initial ray value != length orbital ray!")
  (length(x0), length(x0)) == size(E) || error("Number of variables inconsistent with stochastic matrix size!")
  return Probe(convert(Vector{Float64}, x0), convert(Vector{Float64}, v), Quaternion{Float64}(q.q), convert(Matrix{Float64}, E))
end

# This ctor is used for TPS Probe
function Probe(;v::Vector{TPS}, x0::Vector{<:Real}=zeros(numvars(first(v))), q::Quaternion{<:Real}=Quaternion(one(first(v)), zero(first(v)), zero(first(v)), zero(first(v))), E::Matrix{<:Real}=zeros(numvars(first(v)), numvars(first(v))))
  length(v) == numvars(first(v)) || error("Length of orbital ray != number of variables in GTPSA!")
  numvars(first(v)) == length(x0) || error("Length of initial ray value != number of variables in GTPSA!")
  (numvars(first(v)), numvars(first(v))) == size(E) || error("Number of variables inconsistent with stochastic matrix size!")
  scalar(first(v)) != 0 
  return Probe(convert(Vector{Float64}, x0), v, Quaternion{TPS}(q.q), convert(Matrix{Float64}, E))
end
=#

# Copy ctor:
function Probe(p::Probe)
  return deepcopy(p)
end



