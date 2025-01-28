# =================================================================================== #
# General convenience setters

function setray!(
  r::AbstractArray{<:T}; 
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
) where {T}
  TI.is_tps_type(T) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")

  NV = nvars(first(r))
  if !isnothing(x)
    length(x) <= NV || error("Length of input vector `x` cannot be greater than the number of variables in the TPSA!")
    foreach((out_xi, xi)->copy!(out_xi, xi), view(r, 1:length(x)), x)
  end

  if !isnothing(x_matrix)
    if x_matrix isa AbstractMatrix # Map as a matrix:
      Base.require_one_based_indexing(x_matrix)
      size(x_matrix,1) <= NV || error("Number of rows of `x_matrix` cannot be greater than the number of variables in the TPSA!")
      size(x_matrix,2) <= nmonos(first(r))-1 || error("Number of columns of `x_matrix` cannot be greater than the number of monomial coefficients in the TPSA!")
      for varidx in 1:size(x_matrix,1)
        for monoidx in 1:size(x_matrix,2)
          TI.seti!(r[varidx], x_matrix[varidx,monoidx], monoidx)
        end
      end
    else # Uniform scaling: Making identity map
      for varidx in 1:NV
        TI.seti!(r[varidx], 1, varidx)
      end
    end
  end

  return r
end

function setquat!(
  quat::Quaternion{T}; 
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
) where {T}
  TI.is_tps_type(T) isa TI.IsTPSType || error("Quaternion element type must be a truncated power series type supported by `TPSAInterface.jl`")

  if !isnothing(q)
    if q isa AbstractVector || q isa Quaternion #  TPSA map or scalar part provided:
      length(q) <= 4 || error("Length of input vector `q` cannot be greater than 4!")
      foreach((out_qi, qi)->copy!(out_qi, qi), quat, q)
    else # Uniform scaling: Making identity quaternion
      TI.copy!(quat.q0, 1)
      TI.copy!(quat.q1, 0)
      TI.copy!(quat.q2, 0)
      TI.copy!(quat.q3, 0)
    end
  end

  if !isnothing(q_map)
    Base.require_one_based_indexing(q_map)
    size(q_map,1) <= 4 || error("Number of rows of `q_map` cannot be greater than 4!")
    size(q_map,2) <= nmonos(first(quat))-1 || error("Number of columns of `q_map` cannot be greater than the number of monomial coefficients in the GTPSA!")
    for qidx in 1:size(q_map,1)
      for monoidx in 1:size(q_map,2)
        TI.seti!(quat[qidx], q_map[qidx,monoidx], monoidx)
      end
    end
  end

  return quat
end