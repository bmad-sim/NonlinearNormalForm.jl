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

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, m2::Union{TaylorMap,VectorField}; dospin::Bool=true)
  checkinplace(m, m1, m2)
  
  nv = numvars(m)

  for i=1:nv
    $(Meta.parse(ops[1]))(m.x[i], m1.x[i], m2.x[i])
  end

  if !isnothing(m.Q) && dospin
    $(Meta.parse(ops[1]))(m.Q.q0, m1.Q.q0, m2.Q.q0)
    $(Meta.parse(ops[1]))(m.Q.q1, m1.Q.q1, m2.Q.q1)
    $(Meta.parse(ops[1]))(m.Q.q2, m1.Q.q2, m2.Q.q2)
    $(Meta.parse(ops[1]))(m.Q.q3, m1.Q.q3, m2.Q.q3)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.x0 .= m1.x0
    elseif m2 isa TaylorMap
      m.x0 .= m2.x0
    else
      m.x0 .= 0
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

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, J::UniformScaling, m1::Union{TaylorMap,VectorField}; dospin::Bool=true)
  checkinplace(m, m1)
  
  nv = numvars(m)

  for i=1:nv
    copy!(m.x[i], m1.x[i])
    m.x[i][i] = $(ops[2])(1, m.x[i][i])
  end

  if !isnothing(m.Q) && dospin
    copy!(m.Q.q0, m1.Q.q0)
    copy!(m.Q.q1, m1.Q.q1)
    copy!(m.Q.q2, m1.Q.q2)
    copy!(m.Q.q3, m1.Q.q3)
    m.Q.q0[0] = $(ops[2])(1, m.Q.q0[0])
  end


  if m isa TaylorMap
    if m1 isa TaylorMap
      m.x0 .= m1.x0
    else
      m.x0 .= 0
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

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, J::UniformScaling; dospin::Bool=true)
  checkinplace(m, m1)

  nv = numvars(m)

  for i=1:nv
    copy!(m.x[i], m1.x[i])
    m.x[i][i] = $(ops[2])(m.x[i][i], 1)
  end

  if !isnothing(m.Q) && dospin
    copy!(m.Q.q0, m1.Q.q0)
    copy!(m.Q.q1, m1.Q.q1)
    copy!(m.Q.q2, m1.Q.q2)
    copy!(m.Q.q3, m1.Q.q3)
    m.Q.q0[0] = $(ops[2])(m.Q.q0[0], 1)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.x0 .= m1.x0
    else
      m.x0 .= 0
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
  m = zero_op(m1, m2)
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
function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, a::Number, m1::Union{TaylorMap,VectorField}; dospin::Bool=true)
  checkinplace(m, a, m1)
  
  nv = numvars(m)

  for i=1:nv
    $(Meta.parse(ops[1]))(m.x[i], a, m1.x[i])
  end

  if !isnothing(m.Q) && dospin
    $(Meta.parse(ops[1]))(m.Q.q0, a, m1.Q.q0)
    $(Meta.parse(ops[1]))(m.Q.q1, a, m1.Q.q1)
    $(Meta.parse(ops[1]))(m.Q.q2, a, m1.Q.q2)
    $(Meta.parse(ops[1]))(m.Q.q3, a, m1.Q.q3)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.x0 .= m1.x0
    else
      m.x0 .= 0
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

function $(Meta.parse(ops[1]))(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField}, a::Number; dospin::Bool=true)
  checkinplace(m, a, m1)
  
  nv = numvars(m)

  for i=1:nv
    $(Meta.parse(ops[1]))(m.x[i], m1.x[i], a)
  end

  if !isnothing(m.Q) && dospin
    $(Meta.parse(ops[1]))(m.Q.q0, m1.Q.q0, a)
    $(Meta.parse(ops[1]))(m.Q.q1, m1.Q.q1, a)
    $(Meta.parse(ops[1]))(m.Q.q2, m1.Q.q2, a)
    $(Meta.parse(ops[1]))(m.Q.q3, m1.Q.q3, a)
  end

  if m isa TaylorMap
    if m1 isa TaylorMap
      m.x0 .= m1.x0
    else
      m.x0 .= 0
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
  m = zero_op(m1, a)
  $(Meta.parse(ops[1]))(m, a, m1)
  return m
end

function $(ops[2])(m1::Union{TaylorMap,VectorField}, a::Number)
  m = zero_op(m1, a)
  $(Meta.parse(ops[1]))(m, m1, a)
  return m
end

end
end