# --- unary ---
function +(F1::VectorField)
  F = zero(F1)
  copy!(F,F1)
  return F
end

function -(F1::VectorField)
  F = -1*F1
  return F
end


# --- add and subtract for VectorFields/Maps/UniformScaling (mul defined already for VectorField*Map) ---
for ops = (("add!", :+), ("sub!",:-))
@eval begin
  
function $(Meta.parse(ops[1]))(F::VectorField, F1::Union{VectorField,TaylorMap}, F2::Union{VectorField,TaylorMap}; dospin::Bool=true)
  checkinplace(F, F1, F2)
  
  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], F1.x[i], F2.x[i])
  end

  if !isnothing(F.Q) && dospin
    $(Meta.parse(ops[1]))(F.Q.q0, F1.Q.q0, F2.Q.q0)
    $(Meta.parse(ops[1]))(F.Q.q1, F1.Q.q1, F2.Q.q1)
    $(Meta.parse(ops[1]))(F.Q.q2, F1.Q.q2, F2.Q.q2)
    $(Meta.parse(ops[1]))(F.Q.q3, F1.Q.q3, F2.Q.q3)
  end

  return
end

function $(Meta.parse(ops[1]))(F::VectorField, J::UniformScaling, F1::Union{VectorField,TaylorMap}; dospin::Bool=true)
  checkinplace(F, F1)
  
  nv = numvars(F)

  for i=1:nv
    @inbounds copy!(F.x[i], F1.x[i])
    @inbounds F.x[i][i] = $(ops[2])(1, F.x[i][i])
  end

  if !isnothing(F.Q) && dospin
    copy!(F.Q.q0, F1.Q.q0)
    copy!(F.Q.q1, F1.Q.q1)
    copy!(F.Q.q2, F1.Q.q2)
    copy!(F.Q.q3, F1.Q.q3)
    F.Q.q0[0] = $(ops[2])(1, F.Q.q0[0])
  end

  return
end

function $(Meta.parse(ops[1]))(F::VectorField, F1::Union{VectorField,TaylorMap}, J::UniformScaling; dospin::Bool=true)
  checkinplace(F, F1)
  
  nv = numvars(F)

  for i=1:nv
    @inbounds copy!(F.x[i], F1.x[i])
    @inbounds F.x[i][i] = $(ops[2])(F.x[i][i], 1)
  end

  if !isnothing(F.Q) && dospin
    copy!(F.Q.q0, F1.Q.q0)
    copy!(F.Q.q1, F1.Q.q1)
    copy!(F.Q.q2, F1.Q.q2)
    copy!(F.Q.q3, F1.Q.q3)
    F.Q.q0[0] = $(ops[2])(F.Q.q0[0], 1)
  end

  return
end

function $(ops[2])(F1::VectorField, F2::VectorField)
  checkop(F1, F2)

  # Promote if necessary:
  if eltype(F1.x) == ComplexTPS64
    F = zero(F1)
  else
    F = zero(F2)
  end
  $(Meta.parse(ops[1]))(F, F1, F2)
  return F
end

function $(ops[2])(F1::VectorField, J::UniformScaling)
  F = zero(F1)
  $(Meta.parse(ops[1]))(F, F1, J)
  return F
end


function $(ops[2])(J::UniformScaling, F1::VectorField)
  F = zero(F1)
  $(Meta.parse(ops[1]))(F, J, F1)
  return F
end
  
end
end



# --- add, subtract, multiply, divide for scalars ---
for ops = (("add!", :+), ("sub!",:-), ("mul!",:*), ("div!",:/))
@eval begin
function $(Meta.parse(ops[1]))(F::VectorField, a::Number, F1::VectorField; dospin::Bool=true)
  checkinplace(F, F1, a)
  
  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], a, F1.x[i])
  end

  if !isnothing(F.Q) && dospin
    $(Meta.parse(ops[1]))(F.Q.q0, a, F1.Q.q0)
    $(Meta.parse(ops[1]))(F.Q.q1, a, F1.Q.q1)
    $(Meta.parse(ops[1]))(F.Q.q2, a, F1.Q.q2)
    $(Meta.parse(ops[1]))(F.Q.q3, a, F1.Q.q3)
  end
  return
end

function $(Meta.parse(ops[1]))(F::VectorField, F1::VectorField, a::Number; dospin::Bool=true)
  checkinplace(F, F1, a)

  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], F1.x[i], a)
  end

  if !isnothing(F.Q) && dospin
    $(Meta.parse(ops[1]))(F.Q.q0, F1.Q.q0, a)
    $(Meta.parse(ops[1]))(F.Q.q1, F1.Q.q1, a)
    $(Meta.parse(ops[1]))(F.Q.q2, F1.Q.q2, a)
    $(Meta.parse(ops[1]))(F.Q.q3, F1.Q.q3, a)
  end

  return
end

function $(ops[2])(a::Number, F1::VectorField)
  if a isa Complex
    F = zero(complex(typeof(F1)), use=F1)
  else
    F = zero(F1)
  end
  $(Meta.parse(ops[1]))(F, a, F1)
  return F
end

function $(ops[2])(F1::VectorField, a::Number)
  if a isa Complex
    F = zero(complex(typeof(F1)), use=F1)
  else
    F = zero(F1)
  end
  $(Meta.parse(ops[1]))(F, F1, a)
  return F
end

end
end