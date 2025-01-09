module NonlinearNormalFormStaticArraysExt
import NonlinearNormalForm as NNF
using StaticArrays

NNF.numvars(::NNF.TaylorMap{S,T,U,V}) where {S<:StaticArray,T,U,V} = length(S)
NNF.numnn(::NNF.TaylorMap{S,T,U,V}) where {S,T<:StaticArray,U,V} = length(T)
NNF.numparams(::NNF.TaylorMap{S,T,U,V}) where {S<:StaticArray,T<:StaticArray,U,V} = length(S) - length(T)

NNF.promote_x0_type(::Type{S}, ::Type{G}) where {S<:StaticArray,G<:Union{Number,Complex}} = similar_type(S, promote_type(eltype(S), G))
NNF.promote_x_type(::Type{T}, ::Type{G}) where {T<:StaticArray,G<:Union{Number,Complex}} = similar_type(T, promote_type(eltype(T), G), Size(Size(T)[1]))
NNF.promote_E_type(::Type{V}, ::Type{G}) where {V<:StaticArray,G<:Union{Number,Complex}} = V != Nothing ? similar_type(V, promote_type(eltype(V), G)) : Nothing

function NNF.init_x0(::Type{S}, use) where {S<:StaticArray}
  nv = NNF.numvars(use)
  length(S) == nv || error("Incorrect length for reference orbit: received $(length(S)), require $nv.")
  x0 = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  return x0
end

function NNF.init_x(::Type{T}, use) where {T<:StaticArray}
  desc = NNF.getdesc(use)
  nv = NNF.numvars(desc)
  nn = NNF.numnn(desc)
  length(T) == nn || error("Incorrect length for orbital ray: received $(length(T)), require $nn.")
  if use isa NNF.TaylorMap && eltype(T) == eltype(use.x)
    x = StaticArrays.sacollect(T, (i <= nv ? eltype(T)(use=desc) :  use.x[i]) for i in 1:length(T))
  else
    x = StaticArrays.sacollect(T, eltype(T)(use=desc) for i in 1:length(T))
    for i in nv+1:nn
      x[i][i] = 1
    end
  end
  return x
end

function NNF.init_E(::Type{V}, use) where {V<:StaticArray}
  if V != Nothing
    nv = NNF.numvars(use)
    size(V) == (nv,nv) || error("Incorrect size for stochastic matrix: received $(size(V)), require $((nv,nv)).")
    E = StaticArrays.sacollect(V, 0 for i in 1:length(V))
  else
    E = nothing
  end
  return E
end

end