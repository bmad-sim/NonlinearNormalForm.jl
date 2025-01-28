# Field promotion rules
promote_x0_type(::Type{X0}, ::Type{G}) where {X0<:StaticArray,G<:Union{Number,Complex}} = MVector{Size(X0)[1], promote_type(eltype(X0), G)} 
promote_x_type(::Type{X}, ::Type{G}) where {X<:StaticArray,G<:Union{Number,Complex}} = similar_type(X, promote_type(eltype(X), G), Size(Size(X)[1]))
promote_s_type(::Type{S}, ::Type{G}) where {S<:StaticArray,G<:Union{Number,Complex}} = MMatrix{Size(S)[1], Size(S)[1], promote_type(eltype(S), G)} 
real_x0_type(::Type{X0}) where {X0<:StaticArray} = real(X0)
real_x_type(::Type{X}) where {X<:StaticArray} = real(X)
real_s_type(::Type{S}) where {S<:StaticArray} = real(S)

# Field initialization functions
# Consistency checks are made by the `checkmapsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA
function init_x0(a::Type{X0}, def::AbstractTPSADef) where {X0<:StaticArray}
  x0 = StaticArrays.sacollect(X0, 0 for i in 1:length(a))
  return x0
end

function init_x(::Type{X}, def::AbstractTPSADef, reuse::Union{Nothing,TaylorMap}=nothing) where {X<:StaticArray}
  nv = nvars(def)
  nn = length(X)
  # reuse parameters if applicable
  if reuse isa TaylorMap && eltype(X) == eltype(reuse.x)
    x = StaticArrays.sacollect(X, (i <= nv ? TI.init_tps(TI.numtype(eltype(X)), def) :  reuse.x[i]) for i in 1:length(X))
  else # allocate
    x = StaticArrays.sacollect(X, TI.init_tps(TI.numtype(eltype(X)), def) for i in 1:length(X))
    for i in nv+1:nn
      TI.seti!(x[i], 1, i)
    end
  end
  return x
end

function init_s(::Type{S}, def::AbstractTPSADef) where {S<:StaticArray}
  if S != Nothing
    s = StaticArrays.sacollect(S, 0 for i in 1:length(S))
  else
    s = nothing
  end
  return s
end
