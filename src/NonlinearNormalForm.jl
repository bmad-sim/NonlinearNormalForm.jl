module NonlinearNormalForm

using Reexport
using GTPSA
using LinearAlgebra

@reexport using GTPSA

import Base: getindex,
             setindex!,
             convert,
             ∘,
             promote_rule,
             show,
             complex

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, checksympm


abstract type TaylorMap end 

# helper functions
getdesc(m::TaylorMap) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d))
numvars(m::TaylorMap) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d)).nv
numparams(m::TaylorMap) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d)).np

getdesc(t::Union{TPS,ComplexTPS}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d))
numvars(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).nv
numparams(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).np

getdesc(d::Descriptor) = d
numvars(d::Descriptor) = unsafe_load(d.desc).nv
numparams(d::Descriptor) = unsafe_load(d.desc).nv


#=
zero(m::TaylorMap) = (typeof(m))(use=m)

# Identity map:
function one(m::Union{TaylorMap,Nothing}=nothing)
  return low_one(m, use)
end

function low_one(m::TaylorMap, use::Nothing)
  nv = numvars(m)
  x0 = zeros(typeof(first(m.x0)), nv)
  v = Vector{typeof(first(m.v))}(undef, nv)
  for i=1:nv
    t = TPS(use=first(m.v))
    GTPSA.mad_tpsa_seti!(t.tpsa, Cint(i), 0.0, 1.0)
    @inbounds v[i] = t
  end
  return (typeof(m))(x0, v)
end

function low_one(m::TaylorMap, use::Nothing)

end
=#
#=
# Use latest Descriptor
function one(type::Union{Type{DAMap},Type{GTPSAMap},Type{DAMap{S,T}}, Type{GTPSAMap{S,T}}}) where {S <: Union{Float64,ComplexF64},T <: Union{TPS,ComplexTPS}}
  desc = unsafe_load(GTPSA.desc_current.desc)
  nv = desc.nv
  x0 = zeros(nv) 
  v = Vector{TPS}(undef, nv)
  for i=1:nv
    t = TPS(use=GTPSA.desc_current)
    GTPSA.mad_tpsa_seti!(t.tpsa, Cint(i), 0.0, 1.0)
    @inbounds v[i] = t
  end
  return (type)(x0, v)
end
=#

function checksympm(M::Matrix)
  nv = size(M)[1]
  J = zeros(nv, nv)
  for i=1:2:nv
    J[i:i+1,i:i+1] = [0 1; -1 0];
  end
  res = transpose(M)*J*M
  return sum(abs.(res - J))
end

function checksymp(m::TaylorMap)
  nv = numvars(m)
  J = zeros(nv, nv)
  for i=1:2:nv
    J[i:i+1,i:i+1] = [0 1; -1 0];
  end
  M = jacobian(m.v)
  res = transpose(M)*J*M
  return sum(abs.(res - J))
end

include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("show.jl")

end
