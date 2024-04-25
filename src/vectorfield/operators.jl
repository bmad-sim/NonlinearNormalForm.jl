# --- add and subtract for VectorFields/Maps/UniformScaling (mul defined already for VectorField*Map) ---
for ops = (("add!", :+), ("sub!",:-))
@eval begin
  
function $(Meta.parse(ops[1]))(F::VectorField, F1::Union{VectorField,TaylorMap}, F2::Union{VectorField,TaylorMap}; dospin::Bool=true)
  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], F1.x[i], F2.x[i])
  end

  if !isnothing(F.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(F.Q.q[i], F1.Q.q[i], F2.Q.q[i])
    end
  end

  return
end

function $(Meta.parse(ops[1]))(F::VectorField, J::UniformScaling, F1::Union{VectorField,TaylorMap}; dospin::Bool=true)
  nv = numvars(F)

  for i=1:nv
    @inbounds copy!(F.x[i], F1.x[i])
    @inbounds F.x[i][i] = $(ops[2])(1, F.x[i][i])
  end

  if !isnothing(F.Q) && dospin
    for i=1:4
      @inbounds copy!(F.Q.q[i], F1.Q.q[i])
    end
    F.Q.q[1][0] = $(ops[2])(1, F.Q.q[1][0])
  end

  return
end

function $(Meta.parse(ops[1]))(F::VectorField, F1::Union{VectorField,TaylorMap}, J::UniformScaling; dospin::Bool=true)
  nv = numvars(F)

  for i=1:nv
    @inbounds copy!(F.x[i], F1.x[i])
    @inbounds F.x[i][i] = $(ops[2])(F.x[i][i], 1)
  end

  if !isnothing(F.Q) && dospin
    for i=1:4
      @inbounds copy!(F.Q.q[i], F1.Q.q[i])
    end
    F.Q.q[1][0] = $(ops[2])(F.Q.q[1][0], 1)
  end

  return
end

function $(ops[2])(F1::VectorField, F2::VectorField)
  # Promote if necessary:
  if eltype(F1.x) == ComplexTPS
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
  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], a, F1.x[i])
  end

  if !isnothing(F.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(F.Q.q[i], a, F1.Q.q[i])
    end
  end
  return
end

function $(Meta.parse(ops[1]))(F::VectorField, F1::VectorField, a::Number; dospin::Bool=true)
  nv = numvars(F)

  for i=1:nv
    @inbounds $(Meta.parse(ops[1]))(F.x[i], F1.x[i], a)
  end

  if !isnothing(F.Q) && dospin
    for i=1:4
      @inbounds $(Meta.parse(ops[1]))(F.Q.q[i], F1.Q.q[i], a)
    end
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