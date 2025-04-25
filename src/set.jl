# =================================================================================== #
# General convenience setters

function setray!(
  r::AbstractArray{<:T}; 
  v::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  x_matrix_offset::Integer=0,
) where {T}
  TI.is_tps_type(T) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")

  length(r) <= ndiffs(first(r)) || error("Length of output orbital ray `r` cannot be greater than the number of differentials in the TPSA!")

  if !isnothing(v)
    length(v) <= length(r) || error("Length of input vector `v` cannot be greater than the length of output vector `r`!")
    foreach((out_xi, xi)->copy!(out_xi, xi), view(r, 1:length(v)), v)
  end

  if !isnothing(x_matrix)
    if x_matrix isa AbstractMatrix # Map as a matrix:
      Base.require_one_based_indexing(x_matrix)
      size(x_matrix,1) <= length(r) || error("Number of rows of `x_matrix` cannot be greater than the length of output vector `r`!")
      size(x_matrix,2) <= nmonos(first(r))-1-x_matrix_offset || error("Number of columns of `x_matrix` cannot be greater than the number of monomial coefficients in the TPSA - x_matrix_offset")
      for varidx in 1:size(x_matrix,1)
        for monoidx in 1:size(x_matrix,2)
          TI.seti!(r[varidx], x_matrix[varidx,monoidx], monoidx+x_matrix_offset)
        end
      end
    else # Uniform scaling: Making identity map
      x_matrix_offset == 0 || error("`x_matrix_offset` must be zero for `UniformScaling` `x_matrix`")
      for varidx in 1:length(r)
        for monoidx in 1:length(r)
          if varidx == monoidx
            TI.seti!(r[varidx], 1, varidx)
          else
            TI.seti!(r[varidx], 0, monoidx)
          end
        end
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