#=

Arithmetic operators for maps/vector fields with maps/vector fields 
or UniformScaling/scalars. The destination map/vf is assumed to be 
allocated, and an error will be thrown by GTPSA if the Descriptors 
do not match. If the Descriptors do match, then because of the defined 
map ctors we can guarantee the immutable parameters are there if map.

Contains add!, sub!, mul!, div!,  +, -, *, /. Note that mul! and div! 
for maps with non-scalar types are "special" (e.g. compose! and inv!) 
and so are not defined.

??? HOW TO DEAL WITH STOCHASTIC MATRIX FOR ADD SUB WITH MAPS/VFS AND 
ADD SUB MUL DIV WITHS SCALARS

=#

# --- unary ---
function +(m1::Union{TaylorMap,VectorField})
  F = zero(m1)
  copy!(m,m1)
  return m
end

function -(m1::Union{TaylorMap,VectorField})
  m = -1*m1
  return m
end


# --- add and subtract for VectorFields/Maps/UniformScaling ---
for ops = (("add!", :+), ("sub!",:-))
@eval begin

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, m2::Union{TaylorMap,VectorField}; do_spin::Bool=true)
  checkinplace(m, m1, m2, internal_promotion=true)
  
  nv = nvars(m)

  for i in 1:nv
    TI.$(Meta.parse(ops[1]))(m.v[i], m1.v[i], m2.v[i])
  end

  if !isnothing(m.q) && do_spin
    TI.$(Meta.parse(ops[1]))(m.q.q0, m1.q.q0, m2.q.q0)
    TI.$(Meta.parse(ops[1]))(m.q.q1, m1.q.q1, m2.q.q1)
    TI.$(Meta.parse(ops[1]))(m.q.q2, m1.q.q2, m2.q.q2)
    TI.$(Meta.parse(ops[1]))(m.q.q3, m1.q.q3, m2.q.q3)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.v0 .= m1.v0
    elseif m2 isa TaylorMap
      m.v0 .= m2.v0
    else
      m.v0 .= 0
    end

    #=
    if !isnothing(m.E)
      if m1 isa TaylorMap
        if m2 isa TaylorMap
          @. m.E = $(ops[2])(m1.E, m2.E)
        else
          @. m.E = $(ops[2])(m1.E, 0)
        end
      elseif m2 isa TaylorMap
        @. m.E = $(ops[2])(0, m2.E)
      else
        m.E .= 0
      end
    end
    =#
  end

  return m
end

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, J::UniformScaling, m1::Union{TaylorMap,VectorField}; do_spin::Bool=true)
  checkinplace(m, m1, internal_promotion=true)
  
  nv = nvars(m)

  for i in 1:nv
    copy!(m.v[i], m1.v[i])
    TI.seti!(m.v[i], $(ops[2])(1, TI.geti(m.v[i], i)), i)
  end

  if !isnothing(m.q) && do_spin
    copy!(m.q.q0, m1.q.q0)
    copy!(m.q.q1, m1.q.q1)
    copy!(m.q.q2, m1.q.q2)
    copy!(m.q.q3, m1.q.q3)
    TI.seti!(m.q.q0, $(ops[2])(1, TI.geti(m.q.q0, 0)), 0)
  end


  if m isa TaylorMap
    if m1 isa TaylorMap
      m.v0 .= m1.v0
    else
      m.v0 .= 0
    end
    #=
    if !isnothing(m.E)
      if m1 isa TaylorMap
        @. m.E = $(ops[2])(0, m1.E)
      else
        m.E .= 0
      end
    end
    =#
  end

  return m
end

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, J::UniformScaling; do_spin::Bool=true)
  checkinplace(m, m1, internal_promotion=true)

  nv = nvars(m)

  for i in 1:nv
    copy!(m.v[i], m1.v[i])
    TI.seti!(m.v[i], $(ops[2])(TI.geti(m.v[i], i), 1), i)
  end

  if !isnothing(m.q) && do_spin
    copy!(m.q.q0, m1.q.q0)
    copy!(m.q.q1, m1.q.q1)
    copy!(m.q.q2, m1.q.q2)
    copy!(m.q.q3, m1.q.q3)
    TI.seti!(m.q.q0, $(ops[2])(TI.geti(m.q.q0, 0), 1), 0)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.v0 .= m1.v0
    else
      m.v0 .= 0
    end
    #=
    if !isnothing(m.E)
      if m1 isa TaylorMap
        @. m.E = $(ops[2])(m1.E, 0)
      else
        m.E .= 0
      end
    end
    =#
  end

  return m
end

function $(ops[2])(m1::T, m2::T) where {T<:Union{TaylorMap,VectorField}}
  m = zero(promote_type(typeof(m1), typeof(m2)), m1)
  $(Meta.parse(ops[1]))(m, m1, m2)
  return m
end

function $(ops[2])(m1::T, J::UniformScaling) where {T<:Union{TaylorMap,VectorField}}
  m = zero(m1)
  $(Meta.parse(ops[1]))(m, m1, J)
  return m
end


function $(ops[2])(J::UniformScaling, m1::T) where {T<:Union{TaylorMap,VectorField}}
  m = zero(m1)
  $(Meta.parse(ops[1]))(m, J, m1)
  return m
end

end
end


# --- add, subtract, multiply, divide for scalars ---
for ops = (("add!", :+), ("sub!",:-), ("mul!",:*), ("div!",:/))
@eval begin
function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, a::Number, m1::Union{TaylorMap,VectorField}; do_spin::Bool=true)
  checkinplace(m, a, m1, internal_promotion=true)
  
  nv = nvars(m)

  for i in 1:nv
    TI.$(Meta.parse(ops[1]))(m.v[i], a, m1.v[i])
  end

  if !isnothing(m.q) && do_spin
    TI.$(Meta.parse(ops[1]))(m.q.q0, a, m1.q.q0)
    TI.$(Meta.parse(ops[1]))(m.q.q1, a, m1.q.q1)
    TI.$(Meta.parse(ops[1]))(m.q.q2, a, m1.q.q2)
    TI.$(Meta.parse(ops[1]))(m.q.q3, a, m1.q.q3)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.v0 .= m1.v0
    else
      m.v0 .= 0
    end
    #=
    if !isnothing(m.E)
      if m1 isa TaylorMap
        @. m.E = $(ops[2])(0, m1.E)
      else
        m.E .= 0
      end
    end
    =#
  end

  return m
end

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, a::Number; do_spin::Bool=true)
  checkinplace(m, a, m1, internal_promotion=true)
  
  nv = nvars(m)

  for i in 1:nv
    TI.$(Meta.parse(ops[1]))(m.v[i], m1.v[i], a)
  end

  if !isnothing(m.q) && do_spin
    TI.$(Meta.parse(ops[1]))(m.q.q0, m1.q.q0, a)
    TI.$(Meta.parse(ops[1]))(m.q.q1, m1.q.q1, a)
    TI.$(Meta.parse(ops[1]))(m.q.q2, m1.q.q2, a)
    TI.$(Meta.parse(ops[1]))(m.q.q3, m1.q.q3, a)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.v0 .= m1.v0
    else
      m.v0 .= 0
    end
    #=
    if !isnothing(m.E)
      if m1 isa TaylorMap
        @. m.E = $(ops[2])(m1.E, 0)
      else
        m.E .= 0
      end
    end
    =#
  end
  
  return m
end

function $(ops[2])(a::Number, m1::Union{TaylorMap,VectorField})
  m = zero(promote_type(typeof(m1), typeof(a)), m1)
  $(Meta.parse(ops[1]))(m, a, m1)
  return m
end

function $(ops[2])(m1::Union{TaylorMap,VectorField}, a::Number)
  m = zero(promote_type(typeof(m1), typeof(a)), m1)
  $(Meta.parse(ops[1]))(m, m1, a)
  return m
end

end
end