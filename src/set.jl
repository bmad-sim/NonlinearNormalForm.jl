# =================================================================================== #
# General convenience setters

function setray!(
  a::AbstractArray{TPS{T}}; 
  x::Union{AbstractArray,Nothing}=nothing,
  x_matrix::Union{AbstractArray,Nothing}=nothing,
) where {T<:TPS}

  if !isnothing(x)
    length(x) <= length(a) || error("Length of input vector `x` cannot be greater than length of destination `a`!")
    foreach((out_xi, xi)->setTPS!(out_xi, xi, change=true), view(a, 1:length(x)), x)
  end

  if !isnothing(x_matrix)
    if x_matrix isa AbstractMatrix # Map as a matrix:
      size(x_matrix,1) <= length(a) || error("Number of rows of `x_matrix` cannot be greater than length of destination `a`!")
      size(x_matrix,2) <= GTPSA.numcoefs(first(a))-1 || error("Number of columns of `x_matrix` cannot be greater than the number of coefficients in `a` GTPSA!")
      for varidx in 1:size(x_matrix,1)
        a[varidx][1:size(x_matrix,2)] = view(x_matrix, varidx, :)
      end
    else # Uniform scaling: Making identity map
      for varidx in 1:length(a)
        a[varidx][varidx] = 1
      end
    end
  end
  return a
end

function setQ!(
  a::AbstractArray{TPS{T}}; 
  x::Union{AbstractArray,Nothing}=nothing,
  x_matrix::Union{AbstractArray,Nothing}=nothing,
) where {T<:TPS}

  if !isnothing(m.Q)
    if !isnothing(Q)
      if Q isa AbstractVector || Q isa Quaternion #  TPSA map or scalar part provided:
        length(Q) <= 4 || error("Length of input vector `Q` cannot be greater than 4!")
        foreach((out_Qi, Qi)->setTPS!(out_Qi, Qi, change=true), m.Q, Q)
      else # Uniform scaling: Making identity quaternion:
        m.Q.q0[0] = 1
      end
    end

    if !isnothing(Q_map)
      size(Q_map,1) <= 4 || error("Number of rows of `Q_map` cannot be greater than 4!")
      size(Q_map,2) <= GTPSA.numcoefs(first(m.Q))-1 || error("Number of columns of `Q_map` cannot be greater than the number of coefficients in `use` GTPSA!")
      for qidx in 1:size(Q_map,1)
        m.Q[qidx][1:size(Q_map,2)] = view(Q_map, qidx, :)
      end
    end
  end
  return a
end