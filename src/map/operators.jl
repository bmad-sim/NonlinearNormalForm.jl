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


for t = (:DAMap, :TPSAMap)
@eval begin


# --- compose ---
"""
    ∘(m2::$($t),m1::$($t)) -> $($t)

$($t) composition, $( $t == DAMap ? "ignoring the scalar part of `m1`" : "including the scalar part of `m1`")
"""
∘(m2::$t, m1::$t) = compose(m2, m1)

# --- add ---
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

# --- subtract ---
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
∘(m::$t, J::UniformScaling) = $t(m)
∘(J::UniformScaling, m::$t) = ∘(m,J)
*(m::$t, J::UniformScaling) = ∘(m,J)
*(J::UniformScaling, m::$t) = ∘(m,J)
/(m::$t, J::UniformScaling) = $t(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = $t(m)
end
end