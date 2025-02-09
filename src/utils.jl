# =================================================================================== #
# Helper functions
getdef(m::Union{TaylorMap,VectorField}) = TI.getdef(first(m.x))
nvars(m::Union{TaylorMap,VectorField}) = length(m.x0)
nparams(m::Union{TaylorMap,VectorField}) = length(m.x) - length(m.x0)
ndiffs(m::Union{TaylorMap,VectorField}) = length(m.x)
nmonos(m::Union{TaylorMap,VectorField}) = TI.nmonos(first(m.x))
maxord(m::Union{TaylorMap,VectorField}) = TI.maxord(first(m.x))
# =================================================================================== #
# Similar eltype
# Returns the equivalent of container A to instead have eltype G and dims M
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{N}) where {M,N,A<:Array{T,N} where{T},G}
  arrtype = Base.typename(A).wrapper
  return arrtype{G, M}
end

# Default fallback:
function similar_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{ndims(A)}) where {M,A<:AbstractArray,G}
  return similar_eltype(Array{eltype(A),M}, G, Val{M})
end
#=
# =================================================================================== #
# Array eltype modification
# Returns the equivalent of container A to instead have eltype f(eltype(A),args...) and dims M
function similar_f_eltype(::Type{A}, ::F, ::Type{Val{M}}=Val{N}, args...) where {M,N,A<:Array{T,N} where {T},F}
  similar_eltype()
end


function _promote_array_eltype(::Type{A}, ::Type{G}, ::Type{Val{M}}=Val{ndims(A)}) where {M,A<:AbstractArray,G}
  return similar_eltype(A, promote_type(eltype(A), G), Val{M})
end

function _f_array_eltype(::Type{A}, ::Type{Val{M}}=Val{N}) where {M,N,A<:Array{T,N} where{T}}
  arrtype = Base.typename(A).wrapper
  neweltype = real(eltype(A))
  return arrtype{neweltype, M}
end
# Default fallback:
function _real_array_eltype(::Type{A}, ::Type{Val{M}}=Val{ndims(A)})  where {M,A}
  return _real_array_eltype(Array{eltype(A),M}, Val{M})
end
# =================================================================================== #
# Array TPSA def change
function _change_array_def(::Type{A}, def::AbstractTPSADef, ::Type{Val{M}}=Val{N}) where {M,N,A<:Array{T,N} where{T}}
  arrtype = Base.typename(A).wrapper
  neweltype = TI.init_tps_type(TI.numtype(eltype(A)), def)
  return arrtype{neweltype, M}
end

# Default fallback:
function _change_array_def(::Type{A}, def::AbstractTPSADef, ::Type{Val{M}}=Val{ndims(A)}) where {M,A}
  return _change_array_def(Array{eltype(A),M}, G, Val{M})
end
#_vec_type(::Type{A}) where {T,N,A<:AbstractArray{T,N}} = (Base.typename(A).wrapper){T, 1}
#_mat_type(::Type{A}) where {T,N,A<:AbstractArray{T,N}} = (Base.typename(A).wrapper){T, 2}
=#
# =================================================================================== #
# Coast check

function coastidx(m)
  return -1
  nv = numvars(m)
  for i in nv-1:nv # check only the last two planes
    if abs(m.x[i][0]) < NonlinearNormalForm.coast_threshold
      cycleidx = GTPSA.cycle!(m.x[i], 0, 0, C_NULL, C_NULL)
      if cycleidx == i && abs(m.x[i][i] - 1) < NonlinearNormalForm.coast_threshold
        return i
      end
    end
  end

  return -1
end