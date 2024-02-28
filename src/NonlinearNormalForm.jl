module NonlinearNormalForm

import Base: ∘,
             ^,
             show,
             convert,
             inv,
             ==

using LinearAlgebra,
      Printf,
      Reexport

@reexport using GTPSA

import LinearAlgebra: norm, dot

import GTPSA: numtype, 
              Desc, 
              RTPSA, 
              CTPSA,
              compose!,
              minv!,
              jacobian

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, qmul!,
        normalize!, dot

include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("show.jl")

# helper functions
getdesc(m::Union{Probe{<:Real,TPS,TPS,<:Real},<:TaylorMap}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d))
numvars(m::Union{Probe{<:Real,TPS,TPS,<:Real},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nv
numparams(m::Union{Probe{<:Real,TPS,TPS,<:Real},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).np

getdesc(t::Union{TPS,ComplexTPS}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d))
numvars(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).nv
numparams(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).np

getdesc(t::Vector{<:Union{TPS,ComplexTPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d))
numvars(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nv
numparams(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).np

getdesc(m::Quaternion{<:Union{ComplexTPS,TPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d))
numvars(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nv
numparams(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).np

getdesc(d::Descriptor) = d
numvars(d::Descriptor) = unsafe_load(d.desc).nv
numparams(d::Descriptor) = unsafe_load(d.desc).np

getdesc(d::Nothing) = nothing

function checksymp(M::Matrix{T}) where T<:Number
  s = size(M)
  nv = first(s)
  nv == last(s) || error("Non-square matrix!")
  iseven(nv) || error("Matrix contains odd number of rows/columns!")
  J = zeros(nv,nv)
  for i=1:2:nv
    J[i:i+1,i:i+1] = [0 1; -1 0];
  end
  res = transpose(M)*J*M
  return sum(abs.(res - J))
end

#hessian(m::TaylorMap,include_params=false) = hessian(m.x[1:numvars(m)],include_params=include_params)
jacobian(m::TaylorMap,include_params=false) = jacobian(m.x[1:numvars(m)],include_params=include_params)
checksymp(m::TaylorMap) = checksymp(jacobian(m))


end
