# --- norm ---
function norm(m::Union{TaylorMap,VectorField})
  nrm = 0.

  nv = numvars(m)
  for i=1:nv
    @inbounds nrm += normTPS(m.x[i])
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
for t = (:DAMap, :TPSAMap)
@eval begin    
function complex(m::$t)
  outm = zero_op(m,im)
  copy!(outm, m)
  return outm
end

# --- complex type ---
function complex(::Type{$t{S,T,U,V,W}}) where {S,T,U,V,W}
  return $t{Vector{ComplexF64},Vector{ComplexTPS64},U == Nothing ? Nothing : Quaternion{ComplexTPS64}, V == Nothing ? Nothing : Matrix{ComplexF64}, W}
end

end
end

# --- copy! ---
function copy!(m::TaylorMap, m1::TaylorMap)
  m.x0 .= m1.x0
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  for i=1:nv
    copy!(m.x[i], m1.x[i])
  end

  #@inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q)
    copy!(m.Q.q0, m1.Q.q0)
    copy!(m.Q.q1, m1.Q.q1)
    copy!(m.Q.q2, m1.Q.q2)
    copy!(m.Q.q3, m1.Q.q3)
  end

  if !isnothing(m1.E)
    m.E .= m1.E
  end

  return m
end

# --- clear ---
function clear!(m::TaylorMap)
  nv = numvars(m)
  m.x0 .= 0
  for i=1:nv
    @inbounds clear!(m.x[i])
  end
  if !isnothing(m.Q)
    clear!(m.Q.q0)
    clear!(m.Q.q1)
    clear!(m.Q.q2)
    clear!(m.Q.q3)
  end
  if !isnothing(m.E)
    m.E .= 0
  end
  return
end

# --- jacobian/jacobiant --- 
jacobian(m::TaylorMap;include_params=false) = GTPSA.jacobian(view(m.x, 1:numvars(m)),include_params=include_params)
jacobiant(m::TaylorMap;include_params=false) = GTPSA.jacobiant(view(m.x, 1:numvars(m)), include_params=include_params)

# --- checksymp ---
checksymp(m::TaylorMap) = checksymp(GTPSA.jacobian(m))

# --- cutord ---
function cutord(m1::TaylorMap, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  m = zero(m1)
  cutord!(m, m1, order, spin_order, dospin=dospin)
  return m
end

function cutord!(m::TaylorMap, m1::TaylorMap, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  checkinplace(m, m1)

  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  m.x0 .= m1.x0
  for i=1:nv
    @inbounds cutord!(m.x[i], m1.x[i], order)
  end
  # add immutable parameters to outx
  #@inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q) && dospin
    cutord!(m.Q.q0, m1.Q.q0, spin_order)
    cutord!(m.Q.q1, m1.Q.q1, spin_order)
    cutord!(m.Q.q2, m1.Q.q2, spin_order)
    cutord!(m.Q.q3, m1.Q.q3, spin_order)
  end
  if !isnothing(m1.E)
    m.E .= m.E
  end
  return
end

# --- getord ---
function getord(m1::TaylorMap, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  m = zero(m1)
  getord!(m, m1, order, spin_order, dospin=dospin)
  return m
end

function getord!(m::TaylorMap, m1::TaylorMap, order::Integer, spin_order::Integer=order; dospin::Bool=true)
  checkinplace(m, m1)
  
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  m.x0 .= m1.x0
  for i=1:nv
    @inbounds getord!(m.x[i], m1.x[i], order)
  end
  # add immutable parameters to outx
  #@inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q) && dospin
    getord!(m.Q.q0, m1.Q.q0, spin_order)
    getord!(m.Q.q1, m1.Q.q1, spin_order)
    getord!(m.Q.q2, m1.Q.q2, spin_order)
    getord!(m.Q.q3, m1.Q.q3, spin_order)
  end
  if !isnothing(m1.E)
    m.E .= m.E
  end
  return
end

# --- setmatrix! ---

function setmatrix!(m::DAMap, M::AbstractMatrix)
  Base.require_one_based_indexing(M)
  nv = numvars(m)
  nn = numnn(m)

  nv >= size(M,1) || error("Number of rows in matrix > number of variables in GTPSA!")

  for i=1:size(M,1)
    for j=1:size(M,2)
      @inbounds m.x[i][j] = M[i,j]
    end
  end
end