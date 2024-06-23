# --- norm ---
function norm(m::Union{TaylorMap,VectorField})
  nrm = zero(numtype(m))

  nv = numvars(m)
  for i=1:nv
    @inbounds nrm += norm(m.x[i])
  end
  
  if !isnothing(m.Q)
    for i=1:4
      @inbounds nrm += norm(m.Q.q[i])
    end
  end

  return nrm
end

# --- complex ---
for t = (:DAMap, :TPSAMap)
@eval begin    
function complex(m::$t)
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  x0 = map(t->complex(t), m.x0)
  
  x = Vector{ComplexTPS}(undef, nn)
  for i=1:nv
    @inbounds x[i] = ComplexTPS(m.x[i],use=desc)
  end

  # use same parameters if complex already
  T = eltype(m.x)
  if T == ComplexTPS
    @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)
  else
    @inbounds x[nv+1:nn] .= complexparams(getdesc(first(x)))
  end

  if !isnothing(m.Q)
    q = Vector{ComplexTPS}(undef, 4)
    for i=1:4
      @inbounds q[i] = ComplexTPS(m.Q.q[i],use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = map(t->complex(t), m.E)
  else
    E = nothing
  end
  return $t(x0, x, Q, E, m.idpt)
end

# --- complex type ---
function complex(::Type{$t{S,T,U,V,W}}) where {S,T,U,V,W}
  return $t{ComplexF64,ComplexTPS,U == Nothing ? Nothing : Quaternion{ComplexTPS}, V == Nothing ? Nothing : Matrix{ComplexF64}, W}
end

end
end

# --- copy! ---
function copy!(m::TaylorMap, m1::TaylorMap)
  checkinplace(m, m1)

  m.x0 .= m1.x0
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  for i=1:nv
    copy!(m.x[i], m1.x[i])
  end

  @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q)
    for i=1:4
      @inbounds copy!(m.Q.q[i], m1.Q.q[i])
    end
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
    for i=1:4
      @inbounds clear!(m.Q.q[i])
    end
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
  @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q) && dospin
    for i=1:4
      @inbounds cutord!(m.Q.q[i], m1.Q.q[i], spin_order)
    end
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
  @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q) && dospin
    for i=1:4
      @inbounds getord!(m.Q.q[i], m1.Q.q[i], spin_order)
    end
  end
  if !isnothing(m1.E)
    m.E .= m.E
  end
  return
end