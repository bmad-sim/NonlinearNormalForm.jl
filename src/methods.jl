#=

Non-arithmetic functions acting on both maps and vector fields.

=#


# --- norm ---
function norm(m::Union{TaylorMap,VectorField})
  nrm = 0.

  nv = numvars(m)
  for i=1:nv
    nrm += normTPS(m.x[i])
  end
  
  if !isnothing(m.Q)
    nrm += normTPS(m.Q.q0)
    nrm += normTPS(m.Q.q1)
    nrm += normTPS(m.Q.q2)
    nrm += normTPS(m.Q.q3)
  end

  return nrm
end

# --- complex ---
function complex(m::Union{TaylorMap,VectorField})
  outm = zero_op(m,im)
  copy!(outm, m)
  return outm
end


# --- copy --- 
function copy!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1)
  nv = numvars(m)

  for i=1:nv
    copy!(m.x[i], m1.x[i])
  end

  if !isnothing(m1.Q)
    copy!(m.Q.q0, m1.Q.q0)
    copy!(m.Q.q1, m1.Q.q1)
    copy!(m.Q.q2, m1.Q.q2)
    copy!(m.Q.q3, m1.Q.q3)
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.x0 .= m1.x0
      if !isnothing(m1.E)
        m.E .= m1.E
      end
    else
      m.x0 .= 0
      if !isnothing(m1.E)
        m.E .= 0
      end
    end
  end

  return m
end

copy(m1::Union{TaylorMap,VectorField}) = (m = zero(m1); copy!(m, m1); return m)

# --- real/imag ---
function real!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1) # no checkinplace because m can be real
  nv = numvars(m)

  for i=1:nv
    @FastGTPSA! m.x[i] = real(m1.x[i])
  end

  if !isnothing(m1.Q)
    @FastGTPSA! begin
      m.Q.q0 = real(m1.Q.q0)
      m.Q.q1 = real(m1.Q.q1)
      m.Q.q2 = real(m1.Q.q2)
      m.Q.q3 = real(m1.Q.q3)
    end
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.x0 .= real.(m1.x0)
      if !isnothing(m1.E)
        m.E .= real.(m1.E)
      end
    else
      m.x0 .= 0
      if !isnothing(m1.E)
        m.E .= 0
      end
    end
  end

  return m
end

function imag!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1) # no checkinplace because m can be real

  m.x0 .= imag.(m1.x0)
  nv = numvars(m)

  for i=1:nv
    @FastGTPSA! m.x[i] = imag(m1.x[i])
  end

  if !isnothing(m1.Q)
    @FastGTPSA! begin
      m.Q.q0 = imag(m1.Q.q0)
      m.Q.q1 = imag(m1.Q.q1)
      m.Q.q2 = imag(m1.Q.q2)
      m.Q.q3 = imag(m1.Q.q3)
    end
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.x0 .= imag.(m1.x0)
      if !isnothing(m1.E)
        m.E .= imag.(m1.E)
      end
    else
      m.x0 .= 0
      if !isnothing(m1.E)
        m.E .= 0
      end
    end
  end

  return m
end


real(m1::Union{TaylorMap,VectorField}) = (m = zero(real(typeof(m1)), use=m1); real!(m, m1); return m)
imag(m1::Union{TaylorMap,VectorField}) = (m = zero(real(typeof(m1)), use=m1); imag!(m, m1); return m)

# --- clear! ---
function clear!(m::Union{TaylorMap,VectorField})
  nv = numvars(m)
  
  for i=1:nv
    clear!(m.x[i])
  end
  if !isnothing(m.Q)
    clear!(m.Q.q0)
    clear!(m.Q.q1)
    clear!(m.Q.q2)
    clear!(m.Q.q3)
  end

  if m isa TaylorMap
    m.x0 .= 0
    if !isnothing(m.E)
      m.E .= 0
    end
  end
  return
end


# --- cutord ---
function cutord!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  checkinplace(m, m1)

  nv = numvars(m)
  
  for i=1:nv
    cutord!(m.x[i], m1.x[i], order)
  end

  if !isnothing(m1.Q) && dospin
    cutord!(m.Q.q0, m1.Q.q0, spin_order)
    cutord!(m.Q.q1, m1.Q.q1, spin_order)
    cutord!(m.Q.q2, m1.Q.q2, spin_order)
    cutord!(m.Q.q3, m1.Q.q3, spin_order)
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.x0 .= m1.x0
      if !isnothing(m1.E)
        m.E .= m1.E
      end
    else
      m.x0 .= 0
      if !isnothing(m1.E)
        m.E .= 0
      end
    end
  end
  
  return
end

function cutord(m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  m = zero(m1)
  cutord!(m, m1, order, spin_order, dospin=dospin)
  return m
end

# --- getord ---
function getord!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  checkinplace(m, m1)
  
  nv = numvars(m)

  for i=1:nv
    getord!(m.x[i], m1.x[i], order)
  end

  if !isnothing(m1.Q) && dospin
    getord!(m.Q.q0, m1.Q.q0, spin_order)
    getord!(m.Q.q1, m1.Q.q1, spin_order)
    getord!(m.Q.q2, m1.Q.q2, spin_order)
    getord!(m.Q.q3, m1.Q.q3, spin_order)
  end

  if m isa TaylorMap 
    if m1 isa TaylorMap
      m.x0 .= m1.x0
      if !isnothing(m1.E)
        m.E .= m1.E
      end
    else
      m.x0 .= 0
      if !isnothing(m1.E)
        m.E .= 0
      end
    end
  end

  return
end

function getord(m1::TaylorMap, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  m = zero(m1)
  getord!(m, m1, order, spin_order, dospin=dospin)
  return m
end