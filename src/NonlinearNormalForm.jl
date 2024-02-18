module NonlinearNormalForm

import Base: ∘,
             show,
             convert

using GTPSA,
      LinearAlgebra,
      Printf,
      Reexport

@reexport using GTPSA

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
