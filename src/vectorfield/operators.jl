# Operators with scalars:
for ops = (("add!", :+), ("sub!",:-), ("mul!",:*), ("div!",:/))
@eval begin
function $(Meta.parse(ops[1]))(F::VectorField, a, F1::VectorField)
  nv = numvars(F)

  for i=1:nv
    $(Meta.parse(ops[1]))(F.x[i], a, F1.x[i])
  end

  if !isnothing(F.Q)
    for i=1:4
      $(Meta.parse(ops[1]))(F.Q.q[i], a, F1.Q.q[i])
    end
  end
  return
end

function $(Meta.parse(ops[1]))(F::VectorField, F1::VectorField, a)
  nv = numvars(F)

  for i=1:nv
    $(Meta.parse(ops[1]))(F.x[i], F1.x[i], a)
  end

  if !isnothing(F.Q)
    for i=1:4
      $(Meta.parse(ops[1]))(F.Q.q[i], F1.Q.q[i], a)
    end
  end

  return
end

function $(ops[2])(a::Number, F1::VectorField)
  F = zero(F1)
  $(Meta.parse(ops[1]))(F, a, F1)
  return F
end

function $(ops[2])(F1::VectorField, a::Number)
  F = zero(F1)
  $(Meta.parse(ops[1]))(F, F1, a)
  return F
end

end
end