module NonlinearNormalFormStaticArraysExt
import NonlinearNormalForm as NNF
using StaticArrays

# Static getters
NNF.numvars(::NNF.TaylorMap{X0,X,Q,S}) where {X0<:StaticArray,X,Q,S} = length(X0)
NNF.numvars(::NNF.TaylorMap{X0,X,Q,S}) where {X0,X,Q,S<:StaticArray} = size(S,1)
NNF.numvars(::NNF.TaylorMap{X0,X,Q,S}) where {X0<:StaticArray,X,Q,S<:StaticArray} = length(X0)
NNF.numnn(::NNF.TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S} = length(X)
NNF.numparams(::NNF.TaylorMap{X0,X,Q,S}) where {X0<:StaticArray,X<:StaticArray,Q,S} = length(X) - length(X0)
NNF.numparams(::NNF.TaylorMap{X0,X,Q,S}) where {X0,X<:StaticArray,Q,S<:StaticArray} = length(X) - length(S)
NNF.numparams(::NNF.TaylorMap{X0,X,Q,S}) where {X0<:StaticArray,X<:StaticArray,Q,S<:StaticArray} = length(X) - length(X0)

# Field promotion rules
NNF.promote_x0_type(::Type{X0}, ::Type{G}) where {X0<:StaticArray,G<:Union{Number,Complex}} = MVector{Size(X0)[1], promote_type(eltype(X0), G)} 
NNF.promote_x_type(::Type{X}, ::Type{G}) where {X<:StaticArray,G<:Union{Number,Complex}} = similar_type(X, promote_type(eltype(X), G), Size(Size(X)[1]))
NNF.promote_s_type(::Type{S}, ::Type{G}) where {S<:StaticArray,G<:Union{Number,Complex}} = MMatrix{Size(S)[1], Size(S)[1], promote_type(eltype(S), G)} 
NNF.real_x0_type(::Type{X0}) where {X0<:StaticArray} = real(X0)
NNF.real_x_type(::Type{X}) where {X<:StaticArray} = real(X)
NNF.real_s_type(::Type{S}) where {S<:StaticArray} = real(S)


# Field initialization functions
function NNF.init_x0(::Type{X0}, use::NNF.UseType) where {X0<:StaticArray}
  nv = NNF.numvars(use)
  length(X0) == nv || error("Incorrect length for reference orbit: received $(length(X0)), require $nv.")
  x0 = StaticArrays.sacollect(X0, 0 for i in 1:length(X0))
  return x0
end

function NNF.init_x(::Type{X}, use::NNF.UseType) where {X<:StaticArray}
  desc = NNF.getdesc(use)
  nv = NNF.numvars(desc)
  nn = NNF.numnn(desc)
  length(X) == nn || error("Incorrect length for orbital ray: received $(length(X)), require $nn.")
  if use isa NNF.TaylorMap && eltype(X) == eltype(use.x)
    x = StaticArrays.sacollect(X, (i <= nv ? eltype(X)(use=desc) :  use.x[i]) for i in 1:length(X))
  else
    x = StaticArrays.sacollect(X, eltype(X)(use=desc) for i in 1:length(X))
    for i in nv+1:nn
      x[i][i] = 1
    end
  end
  return x
end

function NNF.init_s(::Type{S}, use::NNF.UseType) where {S<:StaticArray}
  if S != Nothing
    nv = NNF.numvars(use)
    size(S) == (nv,nv) || error("Incorrect size for stochastic matrix: received $(size(S)), require $((nv,nv)).")
    s = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  else
    s = nothing
  end
  return s
end

end