# FIX!!!! Define operators:
for ops = (("add!", :+), ("sub!",:-), ("mul!",:*), ("div!",:/))
  @eval begin
  function $(Meta.parse(ops[1]))(m::TaylorMap, a, m1::TaylorMap)
    nv = numvars(m)
    nn = numnn(m)
  
    m.x0 .= m1.x0
  
    for i=1:nv
      $(Meta.parse(ops[1]))(m.x[i], a, m1.x[i])
    end
    m.x[nv+1:nn] .= view(m.x, nv+1:nn)
  
    if !isnothing(m.Q)
      for i=1:4
        $(Meta.parse(ops[1]))(m.Q.q[i], a, m1.Q.q[i])
      end
    end
  
    if !isnothing(m.E)
      m.E .= m1.E
    end
    return
  end
  
  function $(Meta.parse(ops[1]))(m::TaylorMap, m1::TaylorMap, a)
    nv = numvars(m)
    nn = numnn(m)
  
    m.x0 .= m1.x0
  
    for i=1:nv
      $(Meta.parse(ops[1]))(m.x[i], m1.x[i], a)
    end
    m.x[nv+1:nn] .= view(m.x, nv+1:nn)
  
    if !isnothing(m.Q)
      for i=1:4
        $(Meta.parse(ops[1]))(m.Q.q[i], m1.Q.q[i], a)
      end
    end
  
    if !isnothing(m.E)
      m.E .= m1.E
    end
    return
  end