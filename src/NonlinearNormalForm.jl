module NonlinearNormalForm

import Base: getindex,
             setindex!,
             convert,
             ∘,
             promote_rule

using GTPSA

export ComplexDAMap, DAMap

struct DAMap{T<:Union{TPS,ComplexTPS}}
  v::Vector{T}
end

const ComplexDAMap = DAMap{ComplexTPS}

promote_rule(::Type{DAMap{TPS}}, ::Type{DAMap{TPS}}) = DAMap{TPS}
# promote_rule(::Type{ComplexTPS}, ::Union{Type{Complex},Type{AbstractFloat}, Type{Integer}, Type{Rational}, Type{AbstractIrrational}}) = ComplexTPS
# promote_rule(::Type{TPS}, ::Union{Type{ComplexTPS}, Type{Complex}}) = ComplexTPS

#= Simple helper functions
get_cdesc(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d))
get_cdesc(m::Vector{<:Union{TPS,ComplexTPS}}) = get_cdesc(first(m))
get_cdesc(da::DAMap) = unsafe_load(first(da.v).desc)
=#

function DAMap(da::Union{DAMap,Vector{<:Union{TPS,ComplexTPS,Number}},Nothing}; use::Union{Descriptor,TPS,ComplexTPS,Nothing}=nothing)
  low_DAMap(da, use)
end

# Use latest descriptor
function low_DAMap(da::Nothing, use::Nothing)
  return DAMap(vars(GTPSA.desc_current))
end

# Copy ctor: use DAMap Descriptor
function low_DAMap(da::DAMap, use::Nothing)
  return DAMap(deepcopy(da.v))
end

# Use descriptor in TPSA map
function low_DAMap(m::Union{Vector{<:Union{TPS,ComplexTPS}}}, use::Nothing)
  nv = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).nv
  length(m) == nv || error("Length of $(typeof(m)) != Number of variables")
  return DAMap(deepcopy(m))
end

# Indexing should be like TPSA map, note will return TPSA
getindex(m::Union{DAMap,ComplexDAMap}, i::Integer) = m.v[i]
setindex!(m::Union{DAMap,ComplexDAMap}, t::Number, i::Integer) = setindex!(m.v, t, i)


# Only difference is in map composition
function ∘(d1::Union{DAMap,ComplexDAMap},d2::Union{DAMap,ComplexDAMap})
  ref = Vector{eltype(d2.v)}(undef, length(d2.v))
  for i=1:length(d2.v)
    ref[i] = d2.v[i][0]
    d2.v[i][0] = 0
  end
  tmp = ∘(d1.v, d2.v)
  
  
end



#=
# Operators:
for op = (:+, :-, :*, :/, :^)   #, :± :∓ , :⨰ , :⨱, :⤊ 
  eval(quote
    Base.$op(d1::Union{DAMap,ComplexDAMap},d2::Union{DAMap,ComplexDAMap}) = GTPSA.$op(d1.v, d2.v)
    Base.$op(d1::Union{DAMap,ComplexDAMap}, a::Number) = GTPSA.$op(d1.v, a)
    Base.$op(a::Number, d1::Union{DAMap,ComplexDAMap}) = GTPSA.$op(a, d1.v)
  end)
end
Base.∘(d1::Union{DAMap,ComplexDAMap},d2::Union{DAMap,ComplexDAMap}) = ∘(d1.v, d2.v)
for un = (:abs, :sqrt, :exp, :log, :sin, :cos, :tan, :csc, :sec, :cot, :sinc, :sinh, :cosh, 
          :tanh, :csch, :sech, :coth, :asin, :acos, :atan, :acsc, :asec, :acot, :asinh, :acosh, 
          :atanh, :acsch, :asech, :acoth)

+(d1::DAMap, d2::ComplexDAMap) = +(d1.v, d2.v)
=#  



end
