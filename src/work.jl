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

function prep_comp_work_prom(m::TaylorMap, m2::TaylorMap, m1::TaylorMap)
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  if eltype(m.x) != eltype(m1.x)
    m1x_prom = Vector{ComplexTPS}(undef, nn)
    for i=1:nn  # Allocate
      @inbounds m1x_prom[i] = ComplexTPS(use=desc)
    end
    return (m1x_prom,)
  elseif eltype(m.x) != eltype(m2.x)
    m2x_prom = Vector{ComplexTPS}(undef, nv)
    for i=1:nv
      @inbounds m2x_prom[i] = ComplexTPS(use=desc)
    end
    if !isnothing(m.Q)
      m2Q_prom = Vector{ComplexTPS}(undef, 4)
      @inbounds m2Q_prom[1] = ComplexTPS(use=desc)
      @inbounds m2Q_prom[2] = ComplexTPS(use=desc)
      @inbounds m2Q_prom[3] = ComplexTPS(use=desc)
      @inbounds m2Q_prom[4] = ComplexTPS(use=desc)
      return (m2x_prom, m2Q_prom)
    else
      return (m2x_prom,)
    end
  else
    return nothing
  end
end

function prep_comp_inv_work_low(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  #prepare work
  nn = numnn(m1)
  nv = numvars(m1)
  outx_low = Vector{lowtype(T)}(undef, nn)   # change
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

function prep_inv_work_low(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  nn = numnn(m1)
  outx_low = Vector{lowtype(T)}(undef, nn)    # SHOULD ONLY NEED TO BE NV  BUT GTPSA BUG
  m1x_low = Vector{lowtype(T)}(undef, nn)
  if !isnothing(m1.Q)
    if nn >= 4   # reuse
      outQ_low = m1x_low
    else
      outQ_low = Vector{lowtype(T)}(undef, 4)
    end
    return (outx_low, m1x_low, outQ_low)
  else
    return (outx_low,m1x_low)
  end
end

function prep_work_ref(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  return Vector{numtype(eltype(m1.x))}(undef, numvars(m1))
end

function prep_lb_work_low(F::VectorField{T,U}) where {T,U}
  nv = numvars(F)
  Fx_low = Vector{lowtype(T)}(undef, nv)
  Hx_low = Vector{lowtype(T)}(undef, nv)
  Gx_low = Vector{lowtype(T)}(undef, nv)
  return (Fx_low, Hx_low, Gx_low)
end

function prep_vf_work_Q(F::VectorField{T,U}) where {T,U}
  if !isnothing(F.Q)
    return Quaternion(first(F.Q))
  else
    return nothing
  end
end
