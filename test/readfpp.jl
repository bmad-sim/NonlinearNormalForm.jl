
using DelimitedFiles
function read_fpp_map(file; stochastic::Union{Nothing,Bool}=nothing,spin::Union{Nothing,Bool}=nothing,coast::Bool=false) 
  data = readdlm(file, skipblanks=true)
  nv = data[findfirst(t->t=="Dimensional", data)- CartesianIndex(0,1)]
  no = data[findfirst(t->t=="NO", data) + CartesianIndex(0,2)]
  np = data[findfirst(t->t=="NV", data) + CartesianIndex(0,2)] - nv
  nn = nv+np

  # Check if map has stochasticity
  stoch_idx = findfirst(t->t=="Stochastic", data)
  if isnothing(stochastic) # Default
    if !isnothing(stoch_idx) && stoch_idx[2]== 1 # stochastic in first column meaning no "No"
      stochastic = true
    else
      stochastic = nothing
    end
  elseif stochastic == true
    if !isnothing(stoch_idx) && stoch_idx[2]== 1 # stochastic in first column meaning no "No"
      stochastic = true
    else
      error("Cannot include stochastic: no stochastic matrix detected")
    end
  else
    stochastic = nothing
  end

  # Check if map has spin
  spin_idx = findfirst(t->t=="c_quaternion", data)
  if isnothing(spin) # Default
    if !isnothing(spin_idx) && spin_idx[2]== 1 # c_quaternion in first column meaning no "No"
      spin = true
    else
      spin = nothing
    end
  elseif spin == true
    if !isnothing(spin_idx) && spin_idx[2]== 1 # c_quaternion in first column meaning no "No"
      spin = true
    else
      error("Cannot include spin: no quaternion detected")
    end
  else
    spin = nothing
  end


  # Make the TPSA
  
  if coast
    nvt = nv-1
    npt = np+1
  else
    nvt = nv
    npt = np
  end
  d = InitGTPSA{Descriptor(nvt, no, npt, no),Nothing}()
  # for old GTPSA compatibility with TI:
  #d = InitGTPSA{Nothing,Descriptor}(dynamic_descriptor=Descriptor(nvt, no, npt, no))
  m = complex(DAMap(init=d,nv=nvt,np=npt,stochastic=stochastic,spin=spin))

  idx=3
  data=data[3:end,:]
  # Now fill
  for i=1:nv
    idx = findfirst(x->(x isa Integer), data[:,1])
    count = 0
    while data[idx,1] >= 0
      a = data[idx,2]
      b = data[idx,3]
      ords = data[idx,4:4+nn-1]
      #println(ords)
      TI.setm!(m.x[i], a + im*b, collect(Int, ords))
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
  # dont forget params
  #m.x[nv+1:nn] .= complexparams(d)


  if !isnothing(spin)
    idx = findfirst(t->t=="c_quaternion", data[:,1])
    if data[idx,3] == "identity"
      TI.copy!(m.q.q0, 1)
      TI.copy!(m.q.q1, 0)
      TI.copy!(m.q.q2, 0)
      TI.copy!(m.q.q3, 0)
    else
      for qi in m.q
        idx = findfirst(x->(x isa Integer), data[:,1])
        count = 0
        while data[idx,1] >= 0
          a = data[idx,2]
          b = data[idx,3]
          ords = data[idx,4:4+nn-1]
          #println(ords)
          TI.setm!(qi, a + im*b, collect(Int,ords))
          idx += 1
          count += 1
        end
        if count!=0 && -data[idx,1] != count
          println(m)
          println(data[idx,1])
          println(count)
          error("This should not have been reached! Incorrect number of monomials read for variable $(qi)")
        end
        idx += 1
        data=data[idx:end,:]
      end
    end 
  end

  # stochastic
  if !isnothing(stochastic)
    idx = findfirst(t->t=="Stochastic", data)[1]
    data = data[idx+1:end,:]
    for i=1:size(data, 1)
      row = data[i,1]
      col = data[i,2]
      numstr = data[i,3]
      comidx = findfirst(",", numstr)[1]
      num = complex(parse(Float64, numstr[2:comidx-1]), parse(Float64,numstr[comidx+1:end-1]))
      m.s[row,col] = num
    end
  end

  return m
end


