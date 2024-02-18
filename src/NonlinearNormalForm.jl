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

export DAMap, GTPSAMap, Probe, TaylorMap, Descriptor,

       checksymp


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

# Only difference is in map composition
function ∘(m2::DAMap,m1::DAMap)
  ref = Vector{ComplexF64}(undef, length(m1.v))
  for i=1:length(m1.v)
    @inbounds ref[i] = m1.v[i][0]
    @inbounds m1.v[i][0] = 0
  end
  outv = ∘(m2.v, vcat(m1.v, complexparams(use=first(m1.v))...))
  for i=1:length(m1.v)
    @inbounds m1.v[i][0] = ref[i]
  end
  return DAMap(zeros(ComplexF64, length(m1.v)), outv)
end

function ∘(m2::GTPSAMap,m1::GTPSAMap)
  return GTPSAMap(zeros(ComplexF64, length(m1.v)),∘(m2.v, vcat(m1.v, complexparams(use=first(m1.v))...)))
end

function checksymp(m::TaylorMap, tol)
  nv = numvars(m)
  J = zeros(nv, nv)
  for i=1:2:nv
    J[i:i+1,i:i+1] = [0 1; -1 0];
  end
  M = jacobian(m.v)
  res = tranpose(M)*J*M
  return sum(abs.(res - J))
end


include("show.jl")

end
