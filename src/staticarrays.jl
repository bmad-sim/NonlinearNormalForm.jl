# Field promotion rules
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:StaticArray{N,T} where{T},G}
  return similar_type(A, G, Size(M))
end

# Field initialization functions
# Consistency checks are made by the `checkmapsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA
function init_x0(a::Type{X0}, ::AbstractTPSADef) where {X0<:StaticVector}
  x0 = StaticArrays.sacollect(X0, 0 for i in 1:length(a))
  return x0
end

function init_x(::Type{X}, def::AbstractTPSADef, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:StaticVector}
  nv = nvars(def)
  # reuse parameters if applicable
  if reuse isa TaylorMap && eltype(X) == eltype(reuse.x) && def == getdef(reuse)
    x = StaticArrays.sacollect(X, (i <= nv ? TI.init_tps(TI.numtype(eltype(X)), def) :  reuse.x[i]) for i in 1:length(X))
  else # allocate
    x = StaticArrays.sacollect(X, (i <= nv ? 
                                    TI.init_tps(TI.numtype(eltype(X)), def) : 
                                    (t = TI.init_tps(TI.numtype(eltype(X)), def); TI.seti!(t, 1, i); t)) for i in 1:length(X))
  end
  return x
end

function init_s(::Type{S}, ::AbstractTPSADef) where {S<:StaticMatrix}
  s = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  return s
end
