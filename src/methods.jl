#=

Non-arithmetic functions acting on both maps and vector fields.

=#

# --- clear! ---
function clear!(m::Union{TaylorMap,VectorField})
  nv = nvars(m)
  
  for i in 1:nv
    TI.clear!(m.v[i])
  end
  if !isnothing(m.q)
    TI.clear!(m.q.q0)
    TI.clear!(m.q.q1)
    TI.clear!(m.q.q2)
    TI.clear!(m.q.q3)
  end

  if m isa TaylorMap
    m.v0 .= 0
    if !isnothing(m.s)
      m.s .= 0
    end
  end
  return
end

# --- copy!/copy ---
function copy!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1)
  if m1 isa TaylorMap && m isa TaylorMap
    m.v0 .= m1.v0
  end
  nv = nvars(m)
  foreach((xi, x1i)->TI.copy!(xi, x1i), view(m.v, 1:nv), m1.v)
  if m isa TaylorMap && !isnothing(m.q)
    foreach((qi, q1i)->TI.copy!(qi, q1i), m.q, m1.q)
  end
  if m isa TaylorMap && !isnothing(m.s) && m1 isa TaylorMap
    m.s .= m1.s
  end
  return m
end

copy(m::Union{TaylorMap,VectorField}) = (out_m = zero(m); copy!(out_m, m); return out_m)

# --- convert ---
convert(::Type{T}, m::Union{TaylorMap,VectorField}) where {T<:Union{TaylorMap,VectorField}} = T(m) 
convert(::Type{T}, m::T) where {T<:Union{TaylorMap,VectorField}} = m

# --- norm ---
function norm(m::Union{TaylorMap,VectorField})
  nrm = 0.

  nv = nvars(m)
  for i in 1:nv
    nrm += TI.norm_tps(m.v[i])
  end
  
  if !isnothing(m.q)
    nrm += TI.norm_tps(m.q.q0)
    nrm += TI.norm_tps(m.q.q1)
    nrm += TI.norm_tps(m.q.q2)
    nrm += TI.norm_tps(m.q.q3)
  end

  return nrm
end

# --- real/imag ---
function real!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1) # No checkinplace because m can be real
  if m1 isa TaylorMap && m isa TaylorMap
    @. m.v0 = real(m1.v0)
  end
  nv = nvars(m)
  foreach((xi, x1i)->TI.real!(xi, x1i), view(m.v, 1:nv), m1.v)
  if m isa TaylorMap && !isnothing(m.q)
    foreach((qi, q1i)->TI.real!(qi, q1i), m.q, m1.q)
  end
  if m isa TaylorMap && !isnothing(m.s) && m1 isa TaylorMap
    @. m.s = real(m1.s)
  end
  return m
end

function imag!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1) # No checkinplace because m can be real
  if m1 isa TaylorMap && m isa TaylorMap
    @. m.v0 = imag(m1.v0)
  end
  nv = nvars(m)
  foreach((xi, x1i)->TI.imag!(xi, x1i), view(m.v, 1:nv), m1.v)
  if m isa TaylorMap && !isnothing(m.q)
    foreach((qi, q1i)->TI.imag!(qi, q1i), m.q, m1.q)
  end
  if m isa TaylorMap && !isnothing(m.s) && m1 isa TaylorMap
    @. m.s = imag(m1.s)
  end
  return m
end

function real(m::TaylorMap)
  out_m = _zero(real(typeof(m)), getinit(m), nvars(m), m)
  real!(out_m, m)
  return out_m
end

function imag(m::TaylorMap)
  out_m = _zero(real(typeof(m)), getinit(m), nvars(m), m)
  imag!(out_m, m)
  return out_m
end

function complex(m::TaylorMap)
  out_m = _zero(complex(typeof(m)), getinit(m), nvars(m), m)
  copy!(out_m, m)
  return out_m
end


# --- cutord ---
function cutord!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; do_spin::Bool=true)
  checkinplace(m, m1)

  nv = nvars(m)
  
  for i in 1:nv
    TI.cutord!(m.v[i], m1.v[i], order)
  end

  if !isnothing(m1.q) && do_spin
    TI.cutord!(m.q.q0, m1.q.q0, spin_order)
    TI.cutord!(m.q.q1, m1.q.q1, spin_order)
    TI.cutord!(m.q.q2, m1.q.q2, spin_order)
    TI.cutord!(m.q.q3, m1.q.q3, spin_order)
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.v0 .= m1.v0
      if !isnothing(m1.s)
        m.s .= m1.s
      end
    else
      m.v0 .= 0
      if !isnothing(m1.s)
        m.s .= 0
      end
    end
  end
  
  return
end

function cutord(m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; do_spin::Bool=true)
  m = zero(m1)
  cutord!(m, m1, order, spin_order, do_spin=do_spin)
  return m
end

# --- getord ---
function getord!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; do_spin::Bool=true)
  checkinplace(m, m1)
  
  nv = nvars(m)

  for i in 1:nv
    TI.getord!(m.v[i], m1.v[i], order)
  end

  if !isnothing(m1.q) && do_spin
    TI.getord!(m.q.q0, m1.q.q0, spin_order)
    TI.getord!(m.q.q1, m1.q.q1, spin_order)
    TI.getord!(m.q.q2, m1.q.q2, spin_order)
    TI.getord!(m.q.q3, m1.q.q3, spin_order)
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.v0 .= m1.v0
      if !isnothing(m1.s)
        m.s .= m1.s
      end
    else
      m.v0 .= 0
      if !isnothing(m1.s)
        m.s .= 0
      end
    end
  end

  return
end

function getord(m1::TaylorMap, order::Integer, spin_order::Integer=order; do_spin::Bool=true)
  m = zero(m1)
  getord!(m, m1, order, spin_order, do_spin=do_spin)
  return m
end