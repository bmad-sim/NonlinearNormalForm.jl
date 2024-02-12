module NonlinearNormalForm

import Base: getindex,
             setindex!,
             convert

using GTPSA

export ComplexDAMap, DAMap

struct DAMap
  ref::Vector{Float64}
  v::Vector{TPS}
end

struct ComplexDAMap
  ref::Vector{ComplexF64}
  v::Vector{ComplexTPS}
end

convert(DAMap, m::Vector{TPS}) = DAMap(m)
convert(ComplexDAMap, m::Vector{ComplexTPS}) = ComplexDAMap(m)

function DAMap(d::Descriptor=GTPSA.desc_current)
  v = vars(d)
  nv = length(v)
  ref = zeros(Float64, nv)
  return DAMap(ref, v)
end

function DAMap(m::Vector{TPS})
  nv = GTPSA.get_cdesc(first(m)).nv
  length(m) == nv || error("Length of $(typeof(m)) != Number of variables")
  ref = Vector{Float64}(undef, nv)
  da = DAMap(ref, m)
  rezero!(da)
  return da
end

function ComplexDAMap(d::Descriptor=GTPSA.desc_current)
  v = complexvars(d)
  nv = length(v)
  ref = zeros(ComplexF64, nv)
  return ComplexDAMap(ref, v)
end

function ComplexDAMap(m::Vector{ComplexTPS})
  nv = GTPSA.get_cdesc(first(m)).nv
  length(m) == nv || error("Length of $(typeof(m)) != Number of variables")
  ref = Vector{ComplexF64}(undef, nv)
  da = ComplexDAMap(ref, m)
  rezero!(da)
  return da
end

# Indexing should be like TPSA
getindex(m::Union{DAMap,ComplexDAMap}, i::Integer) = m.v[i]
setindex!(m::Union{DAMap,ComplexDAMap}, t::Number, i::Integer) = setindex!(m.v, t, i)

# For all operations, we need something to get the zeroth order part after and put into ref
function rezero!(da::Union{DAMap,ComplexDAMap})
  for i in LinearIndices(da.v)
    @inbounds da.ref[i] = da.v[i][0]
    @inbounds da.v[i][0] = 0
  end
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
