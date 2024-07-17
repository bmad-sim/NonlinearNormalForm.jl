function prep_comp_work_prom(m::TaylorMap, m2::TaylorMap, m1::TaylorMap)
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  if eltype(m.x) != eltype(m1.x)
    m1x_prom = Vector{ComplexTPS64}(undef, nn)
    for i=1:nn  # Allocate
      @inbounds m1x_prom[i] = ComplexTPS64(use=desc)
    end
    return (m1x_prom,)
  elseif eltype(m.x) != eltype(m2.x)
    m2x_prom = Vector{ComplexTPS64}(undef, nv)
    for i=1:nv
      @inbounds m2x_prom[i] = ComplexTPS64(use=desc)
    end
    if !isnothing(m.Q)
      m2Q_prom = Vector{ComplexTPS64}(undef, 4)
      @inbounds m2Q_prom[1] = ComplexTPS64(use=desc)
      @inbounds m2Q_prom[2] = ComplexTPS64(use=desc)
      @inbounds m2Q_prom[3] = ComplexTPS64(use=desc)
      @inbounds m2Q_prom[4] = ComplexTPS64(use=desc)
      return (m2x_prom, m2Q_prom)
    else
      return (m2x_prom,)
    end
  else
    return nothing
  end
end

function prep_work_ref(m1::TaylorMap)
  return Vector{eltype(eltype(m1.x))}(undef, numvars(m1))
end

function prep_vf_work_Q(F::VectorField)
  if !isnothing(F.Q)
    T = eltype(F.Q)
    use = getdesc(F)
    return Quaternion(T(use=use), T(use=use), T(use=use), T(use=use))
  else
    return nothing
  end
end

function prep_log_work(m::DAMap{S,T,U,V,W}) where {S,T,U,V,W}
  return (zero(m), zero(m), zero(m), zero(VectorField{T,U},use=m), zero(VectorField{T,U},use=m))
end