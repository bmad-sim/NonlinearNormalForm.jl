function prep_comp_work_low(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  #prepare work
  nn = numnn(m1)
  nv = numvars(m1)
  outx_low = Vector{lowtype(T)}(undef, nv)
  m2x_low = Vector{lowtype(T)}(undef, nv)
  m1x_low = Vector{lowtype(T)}(undef, nn)
  if !isnothing(m1.Q)
    if nv >= 4 # Reuse container
      outQ_low = outx_low
      m2Q_low = m2x_low
    else
      outQ_low = Vector{lowtype(T)}(undef, 4)
      m2Q_low = Vector{lowtype(T)}(undef, 4)
    end
    work_low = (outx_low, m2x_low, m1x_low, outQ_low, m2Q_low)
  else
    work_low = (outx_low, m2x_low, m1x_low)
  end

  return work_low
end

function prep_comp_inv_work_low(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  #prepare work
  nn = numnn(m1)
  nv = numvars(m1)
  outx_low = Vector{lowtype(T)}(undef, nn)
  m2x_low = Vector{lowtype(T)}(undef, nv)
  m1x_low = Vector{lowtype(T)}(undef, nn)
  if !isnothing(m1.Q)
    if nv >= 4 # Reuse container
      outQ_low = outx_low
      m2Q_low = m2x_low
    else
      outQ_low = Vector{lowtype(T)}(undef, 4)
      m2Q_low = Vector{lowtype(T)}(undef, 4)
    end
    comp_work_low = (outx_low, m2x_low, m1x_low, outQ_low, m2Q_low)
    inv_work_low = (outx_low,m1x_low,m2Q_low)
  else
    comp_work_low = (outx_low, m2x_low, m1x_low)
    inv_work_low = (outx_low,m1x_low)
  end
  return comp_work_low, inv_work_low
end

function prep_inv_work_low(m::TaylorMap{S,T,U,V}) where {S,T,U,V}

end