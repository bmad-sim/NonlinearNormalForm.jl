# --- map powers ---
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

literal_pow(::typeof(^), m::$t{S,T,U,V}, vn::Val{n}) where {S,T,U,V,n} = ^(m,n)

# Also allow * for simpliticty and \ and / because why not
*(m2::$t,m1::$t) = ∘(m2,m1)
/(m2::$t,m1::$t) = m2∘inv(m1) 
\(m2::$t,m1::$t) = inv(m2)∘m1

# Uniform scaling for * (∘) and /, \
compose!(m::$t, m1::$t, J::UniformScaling) = copy!(m, m1)
compose!(m::$t, J::UniformScaling, m1::$t) = copy!(m, m1)
∘(m::$t, J::UniformScaling) = $t(m)
∘(J::UniformScaling, m::$t) = $t(m)
*(m::$t, J::UniformScaling) = $t(m)
*(J::UniformScaling, m::$t) = $t(m)
/(m::$t, J::UniformScaling) = $t(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = $t(m)
end
end

# --- add and subtract for maps/UniformScaling (other operators defined already) ---
# For our in-place operators, we assume that `m` is an allocated map, so the 
# immutable parameters are already in there and we do not need to worry about them
for ops = (("add!", :+), ("sub!",:-))
@eval begin

function $(Meta.parse(ops[1]))(m::TaylorMap, m1::TaylorMap, m2::TaylorMap; dospin::Bool=true)
  nv = numvars(m)

  m.x0 .= m1.x0

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(m.x[i], m1.x[i], m2.x[i])
  end

  if !isnothing(m.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(m.Q.q[i], m1.Q.q[i], m2.Q.q[i])
    end
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end
  return
end

function $(Meta.parse(ops[1]))(m::TaylorMap, J::UniformScaling, m1::TaylorMap; dospin::Bool=true)
  nv = numvars(m)

  m.x0 .= m1.x0

  for i=1:nv
    @inbounds copy!(m.x[i], m1.x[i])
    @inbounds m.x[i][i] = $(ops[2])(1, m.x[i][i])
  end

  if !isnothing(m.Q) && dospin
    for i=1:4
      @inbounds copy!(m.Q.q[i], m1.Q.q[i])
    end
    m.Q.q[1][0] = $(ops[2])(1, m.Q.q[1][0])
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end
  return
end

function $(Meta.parse(ops[1]))(m::TaylorMap, m1::TaylorMap, J::UniformScaling; dospin::Bool=true)
  nv = numvars(m)

  m.x0 .= m1.x0

  for i=1:nv
    @inbounds copy!(m.x[i], m1.x[i])
    @inbounds m.x[i][i] = $(ops[2])(m.x[i][i], 1)
  end

  if !isnothing(m.Q) && dospin
    for i=1:4
      @inbounds copy!(m.Q.q[i], m1.Q.q[i])
    end
    m.Q.q[1][0] = $(ops[2])(m.Q.q[1][0], 1)
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end
  return
end

function $(ops[2])(m1::TaylorMap, m2::TaylorMap)
  # Promote if necessary:
  if eltype(m1.x) == ComplexTPS
    m = zero(m1)
  else
    m = zero(m2)
  end
  $(Meta.parse(ops[1]))(m, m1, m2)
  return m
end

function $(ops[2])(m1::TaylorMap, J::UniformScaling)
  m = zero(m1)
  $(Meta.parse(ops[1]))(m, m1, J)
  return m
end


function $(ops[2])(J::UniformScaling, m1::TaylorMap)
  m = zero(m1)
  $(Meta.parse(ops[1]))(m, J, m1)
  return m
end

end
end


# --- add, subtract, multiply, divide for scalars ---
for ops = (("add!", :+), ("sub!",:-), ("mul!",:*), ("div!",:/))
@eval begin
function $(Meta.parse(ops[1]))(m::TaylorMap, a::Number, m1::TaylorMap; dospin::Bool=true)

  m.x0 .= m1.x0

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(m.x[i], a, m1.x[i])
  end

  if !isnothing(m.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(m.Q.q[i], a, m1.Q.q[i])
    end
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end
  return
end

function $(Meta.parse(ops[1]))(m::TaylorMap, m1::TaylorMap, a::Number; dospin::Bool=true)
  nv = numvars(m)

  m.x0 .= m1.x0

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(m.x[i], m1.x[i], a)
  end

  if !isnothing(m.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(m.Q.q[i], m1.Q.q[i], a)
    end
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end
  return
end

function $(ops[2])(a::Number, m1::TaylorMap)
  if a isa Complex
    m = zero(complex(typeof(m1)), use=m1)
  else
    m = zero(m1)
  end
  $(Meta.parse(ops[1]))(m, a, m1)
  return m
end

function $(ops[2])(m1::TaylorMap, a::Number)
  if a isa Complex
    m = zero(complex(typeof(m1)), use=m1)
  else
    m = zero(m1)
  end
  $(Meta.parse(ops[1]))(m, m1, a)
  return m
end

end
end