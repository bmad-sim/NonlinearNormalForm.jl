function read_fpp_map(file)
  data = readdlm(file, skipblanks=true)
  nv = data[findfirst(t->t=="Dimensional", data)- CartesianIndex(0,1)]
  no = data[findfirst(t->t=="NO", data) + CartesianIndex(0,2)]
  np = data[findfirst(t->t=="NV", data) + CartesianIndex(0,2)] - nv

  # Make the TPSA
  d = Descriptor(nv, no, np, no)
  m = DAMap(repeat([ComplexTPS(use=d)], nv), Q=Quaternion([ComplexTPS(1,use=d), repeat([ComplexTPS(use=d)], 3)...]))

  idx=3
  data=data[3:end,:]
  # Now fill
  for i=1:nv
    idx = findfirst(x->(x isa Integer), data[:,1])
    println(idx)
    println(data[idx,1])
    count = 0
    while data[idx,1] >= 0
      a = data[idx,2]
      b = data[idx,3]
      ords = data[idx,4:end]
      m.x[i][ords...] = a + im*b
      idx += 1
      count += 1
    end
    if count != 0 && -data[idx,1] != count
      println(m)
      println(data[idx,1])
      println(count)
      error("This should not have been reached! Incorrect number of monomials read for variable $(i)")
    end
    idx += 1
    data=data[idx:end,:]
  end

  # spin?
  idx = findfirst(t->t=="c_quaternion", data[:,1])
  if idx > 0
    for i=1:4
      idx = findfirst(x->(x isa Integer), data[:,1])
      println(idx)
      println(data[idx,1])
      count = 0
      while data[idx,1] >= 0
        a = data[idx,2]
        b = data[idx,3]
        ords = data[idx,4:end]
        m.Q.q[i][ords...] = a + im*b
        idx += 1
        count += 1
      end
      if count!=0 && -data[idx,1] != count
        println(m)
        println(data[idx,1])
        println(count)
        error("This should not have been reached! Incorrect number of monomials read for variable $(i)")
      end
      idx += 1
      data=data[idx:end,:]
    end
  end 
  return m
end