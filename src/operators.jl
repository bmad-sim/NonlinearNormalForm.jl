# --- powers ---

function ^(m1::DAMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  nv = numvars(m1)
  # Do it
  if n>0
    #prepare work
    work_low = prep_comp_work_low(m1)
    ref = prep_work_ref(m1)

    # do it:
    m = DAMap(m1)
    map!(t->t[0], ref, view(m.x,1:nv))
    for i=1:(n-1)
      compose!(m,m,m1, keep_scalar=false,work_low=work_low)
    end
    for i=1:nv
        @inbounds m1.x[i][0] = ref[i]
        #@inbounds m.x[i][0] += ref[i]
    end
    return m
  elseif n<0
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m1)
    ref = prep_work_ref(m1)
    #println(typeof(ref))

    m = DAMap(m1)
    map!(t->t[0], ref, view(m.x,1:nv))
    for i=1:(-n-1)
      compose!(m,m,m1, keep_scalar=false,work_low=comp_work_low)
    end
    for i=1:nv
        @inbounds m1.x[i][0] = ref[i]
        #@inbounds m.x[i][0] += ref[i]
    end
    inv!(m,m,work_ref=ref,work_low=inv_work_low)
    return m
  else
    return DAMap(m1)
  end
end

function ^(m1::TPSAMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  # Do it
  if n>0
    work_low = prep_comp_work_low(m1)
    m = TPSAMap(m1)
    for i=1:(n-1)
      compose!(m,m,m1,work_low=work_low)
    end
    return m
  elseif n<0
    #prepare work
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m1)
    m = TPSAMap(m1)
    for i=1:(-n-1)
      compose!(m,m,m1,work_low=comp_work_low)
    end
    inv!(m,m,work_low=inv_work_low)
    return m
  else
    return TPSAMap(m1)
  end
end

# --- others ---

for t = (:DAMap, :TPSAMap)
@eval begin

function ∘(m2::$t,m1::$t)
  @assert !isnothing(m1.Q) && !isnothing(m2.Q) || m1.Q == m2.Q "Cannot compose: one map includes spin, other does not"
  @assert !isnothing(m1.E) && !isnothing(m2.E) || m1.E == m2.E "Cannot compose: one map includes radiation, other does not"

  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)

  outT = promote_type(eltype(m2.x),eltype(m1.x))
  
  # set up outx0
  outx0 = Vector{numtype(outT)}(undef, nv)

  # Set up outx:
  outx = Vector{outT}(undef, nn)
  for i=1:nv  # no need to allocate immutable parameters taken care of inside compose_it!
      @inbounds outx[i] = outT(use=desc)
  end

  # set up quaternion out:
  if !isnothing(m1.Q)
    outq = Vector{outT}(undef, 4)
    for i=1:4
      @inbounds outq[i] = outT(use=desc)
    end
    outQ = Quaternion(outq)
  else
    outQ = nothing
  end

  # set up radiation out
  if isnothing(m1.E)
    outE = nothing
  else
    outE = Matrix{numtype(outT)}(undef, nv, nv)
  end

  m = $t(outx0, outx, outQ, outE)
  compose!(m, m2, m1)
  
  return m
end

# Basic operators
function +(m2::$t{S,T,U,V},m1::$t{S,T,U,V}) where {S,T,U,V}
  if xor(isnothing(m2.Q), isnothing(m1.Q))
    error("Cannot +: one map includes spin, the other does not")
  end
  if xor(isnothing(m2.E), isnothing(m1.E))
    error("Cannot +: one map includes radiation, the other does not")
  end

  if !isnothing(m2.Q)
    Q1 = Quaternion(m2.Q.q+m1.Q.q)
  else
    Q1 = nothing
  end

  if !isnothing(m2.E)
    E1 =  m2.E+m1.E
  else
    E1 = nothing
  end

  return $t(m2.x0+m1.x0, m2.x+m1.x, Q1, E1)
end

function -(m2::$t,m1::$t)
  if xor(isnothing(m2.Q), isnothing(m1.Q))
    error("Cannot -: one map includes spin, the other does not")
  end
  if xor(isnothing(m2.E), isnothing(m1.E))
    error("Cannot -: one map includes radiation, the other does not")
  end

  if !isnothing(m2.Q)
    Q1 = Quaternion(m2.Q.q-m1.Q.q)
  else
    Q1 = nothing
  end

  if !isnothing(m2.E)
    E1 =  m2.E-m1.E
  else
    E1 = nothing
  end

  return $t(m2.x0-m1.x0, m2.x-m1.x, Q1, E1)
end

literal_pow(::typeof(^), m::$t{S,T,U,V}, vn::Val{n}) where {S,T,U,V,n} = ^(m,n)

# Also allow * for simpliticty and \ and / because why not
*(m2::$t,m1::$t) = ∘(m2,m1)
/(m2::$t,m1::$t) = m2∘inv(m1) 
\(m2::$t,m1::$t) = inv(m2)∘m1

# Uniform scaling (identity map)
+(m::$t, J::UniformScaling) = m + one(m)
+(J::UniformScaling, m::$t) = +(m,J)
-(m::$t, J::UniformScaling) = m - one(m)
-(J::UniformScaling, m::$t) = one(m) - m
∘(m::$t, J::UniformScaling) = DAMap(m)
∘(J::UniformScaling, m::$t) = ∘(m,J)
*(m::$t, J::UniformScaling) = ∘(m,J)
*(J::UniformScaling, m::$t) = ∘(m,J)
/(m::$t, J::UniformScaling) = DAMap(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = DAMap(m)
end
end