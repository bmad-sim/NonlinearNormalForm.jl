module NonlinearNormalForm

using Reexport
using GTPSA

@reexport using GTPSA

import Base: getindex,
             setindex!,
             convert,
             ∘,
             one,
             promote_rule,
             show

export ComplexDAMap, DAMap, GTPSAMap, ComplexGTPSAMap, TaylorMap, Descriptor

abstract type TaylorMap{S<:Union{Float64,ComplexF64}, T<:Union{TPS,ComplexTPS}} end

struct DAMap{S,T} <: TaylorMap{S,T}
  x0::Vector{S}   # ENTRANCE VALUE OF THE MAP!
  v::Vector{T}    # Expansion around x0, with scalar part equal to EXIT value of map
end

struct GTPSAMap{S,T} <: TaylorMap{S,T}
  x0::Vector{S} 
  v::Vector{T}
end

#const ComplexGTPSAMap = GTPSAMap{ComplexF64, ComplexTPS}
#const ComplexDAMap = DAMap{ComplexF64, ComplexTPS}

#promote_rule(t1::Type{DAMap{<:Real,TPS}}, t2::Type{DAMap{<:Real,TPS}}) = t1{Float64,TPS}
#promote_rule(t1::Type{DAMap{<:Number,Union{TPS,ComplexTPS}}},t2::Type{DAMap{<:Number,Union{TPS,ComplexTPS}}}) = t1{ComplexF64,ComplexTPS}

#promote_rule(t1::Type{GTPSAMap{<:Real,TPS}}, t2::Type{GTPSAMap{<:Real,TPS}}) = t1{Float64,TPS}
#promote_rule(t1::Type{GTPSAMap{<:Number,Union{TPS,ComplexTPS}}},t2::Type{GTPSAMap{<:Number,Union{TPS,ComplexTPS}}}) = t1{ComplexF64,ComplexTPS}

# Disallow just adding scalars for now, maybe forever:
# promote_rule(::Type{DAMap{TPS}}, ::Type{Vector{<:Union{Type{AbstractFloat}, Type{Integer}, Type{Rational}, Type{AbstractIrrational}}}}) = DAMap{TPS}
# promote_rule(::Type{DAMap{ComplexTPS}}, ::Type{Vector{<:Union{Type{Complex},Type{AbstractFloat}, Type{Integer}, Type{Rational}, Type{AbstractIrrational}}}}) = DAMap{ComplexTPS}

function DAMap(m::Union{TaylorMap,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(DAMap, m, use)
end

function GTPSAMap(m::Union{TaylorMap,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(GTPSAMap, m, use)
end

# Create identity map using specified Descriptor
function low_map(type::Type, m::Nothing, use::Descriptor)
  desc = unsafe_load(use.desc)
  nv = desc.nv
  np = desc.np
  x0 = zeros(nv) 
  v = Vector{TPS}(undef, nv+np)
  for i=1:nv+np
    t = TPS(use=use)
    GTPSA.mad_tpsa_seti!(t.tpsa, Cint(i), 0.0, 1.0)
    v[i] = t
  end
  return (type)(x0, v)
end

# Use latest descriptor
function low_map(type::Type, m::Nothing, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# Copy ctor: use DAMap/GTPSA map Descriptor
function low_map(type::Type, m::TaylorMap, use::Nothing)
  return (type)(deepcopy(m.x0), deepcopy(m.v))
end

function low_map(type::Type, m::TaylorMap, use::Union{Descriptor,TaylorMap})
  error("Changing descriptor of maps not yet implemented.")
end

# Use descriptor in Vector of TPSA. This probably shouldn't be done.
# Do we need a probe?
#= NOT ALLOWED!
function low_map(type::Type, m::Vector{<:Union{TPS,ComplexTPS}}, use::Nothing)
  nv = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m).tpsa).d)).nv
  length(m) == nv || error("Length of $(typeof(m)) != Number of variables")
  return (type)(zeros(GTPSA.numtype(first(m)), nv), deepcopy(m))
end
=#

getdesc(m::TaylorMap) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d))
numvars(m::TaylorMap) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d)).nv
numparams(m::TaylorMap) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.v).tpsa).d)).np



# Only difference is in map composition
function ∘(d1::DAMap,d2::DAMap)
  ref = Vector{eltype(d2.v)}(undef, length(d2.v))
  for i=1:length(d2.v)
    ref[i] = d2.v[i][0]
    d2.v[i][0] = 0
  end
  out = ∘(d1.v, d2.v)
  for i=1:length(d2.v)
    d2.v[i][0] = ref[i]
  end
  return out
end

include("show.jl")



end
