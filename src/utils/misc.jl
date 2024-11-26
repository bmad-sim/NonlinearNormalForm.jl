#=

coastidx is also defined to quickly check for the coasting plane case, returning 
the index of which variable is constant (energy-like). Note that the last plane 
must be the coasting plane.

=#

# --- helper functions to get the numvars/numparams/descriptor in any context ---

# GTPSA provides these functions for only pure TPS/ComplexTPSs
getdesc(m::Union{TaylorMap,VectorField}) = getdesc(first(m.x))
numvars(m::Union{TaylorMap,VectorField}) = numvars(first(m.x))
numparams(m::Union{TaylorMap,VectorField}) = numparams(first(m.x))
numnn(m::Union{TaylorMap,VectorField}) = numnn(first(m.x))
#=
getdesc(m::Quaternion{<:Union{ComplexTPS,TPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d))
numvars(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nv
numparams(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).np
numnn(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nn
=#
eltype(m::Union{TaylorMap{<:Number,T,<:Any,<:Any},VectorField{T,<:Any}}) where {T<:TPS} = eltype(T)

maxord(m::Union{TaylorMap,VectorField}) = unsafe_load(getdesc(m).desc).mo
prmord(m::Union{TaylorMap,VectorField}) = unsafe_load(getdesc(m).desc).po
vpords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numnn(m))
vords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numvars(m))
pords(m::Union{TaylorMap,VectorField}) = unsafe_wrap(Vector{UInt8}, unsafe_load(getdesc(m).desc).no, numparams(m))

function coastidx(m)
  nv = numvars(m)
  for i in nv-1:nv # check only the last two planes
    if abs(m.x[i][0]) < NonlinearNormalForm.coast_threshold
      cycleidx = GTPSA.cycle!(m.x[i], 0, 0, C_NULL, C_NULL)
      if cycleidx == i && abs(m.x[i][i] - 1) < NonlinearNormalForm.coast_threshold
        return i
      end
    end
  end

  return -1
end


#=
# --- random symplectic map ---
function rand(t::Union{Type{DAMap},Type{TPSAMap}}; spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TPS,ComplexTPS}=GTPSA.desc_current, ndpt::Union{Nothing,Integer}=nothing)
  if isnothing(spin)
    U = Nothing
  else
    if spin
      U = Quaternion{TPS}
    else
      U = Nothing
    end
  end

  if isnothing(FD)
    V = Nothing
  else
    if FD
      V = Matrix{Float64}
    else
      V = Nothing
    end
  end

  return rand(t{Float64,TPS,U,V},use=use)
end


"""
Generate map symplectic up to order in Descriptor
"""
function rand(t::Union{Type{DAMap},Type{TPSAMap}}; require_stable::Bool=true, use::Union{Descriptor,TPS,ComplexTPS}=GTPSA.desc_current, spin::Bool=false, FD::Bool=false) where {S,T,U,V}
  desc = getdesc(use)
  desc.desc != C_NULL || error("No Descriptor defined!")

  # Must create a descriptor that is one order higher than use to create VectorField
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  no = unsafe_wrap(Vector{UInt8}, unsafe_load(desc.desc).no, nn)
  mo = unsafe_load(desc.desc).mo

  if np > 0
    po = unsafe_load(desc.desc).po
    @views dtmp = Descriptor(no[1:nv].+1, mo+1, no[nv+1:nn].+1, po+1)
  else
    @views dtmp = Descriptor(no[1:nv].+1, mo+1)
  end
  
  if require_stable # make hamiltonian in phasors basis then reverse
    h = ComplexTPS(use=dtmp)
    len = length(h)
    #  random coefficients for hamiltonian except 0th and 1st order terms
    for i=nv+1:len-1
      h[i] = rand(eltype(T))
    end
    c = from_phasor(DAMap(spin=spin,FD=FD))
    hc = h * c
  else
    h = TPS(use=dtmp)
    len = length(h)
  
    #  random coefficients for hamiltonian except 0th and 1st order terms
    for i=nv+1:len-1
      h[i] = rand(eltype(T))
    end
    hc = h
  end
  
  F = VectorField(h)
  mtmp = exp(F)

  # Make a copy and change the descriptor
  if t == DAMap{S,T,U,V}
    m = DAMap(mtmp,use=desc)
  else
    m = TPSAMap(mtmp,use=desc)
  end

  # Reset GTPSA descriptor
  GTPSA.desc_current=desc

  if U != Nothing
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = rand(T,use=desc)
    end
    q = q/sqrt(dot(q,q))
    m.Q.q .= q
  end

  return m
end

=#




function read_fpp_map(file; stochastic::Union{Nothing,Bool}=nothing,spin::Union{Nothing,Bool}=nothing) 
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
  d = Descriptor(nv, no, np, no)
  m = complex(DAMap(use=d,stochastic=stochastic,spin=spin)) #repeat([ComplexTPS(use=d)], nv), Q=Quaternion([ComplexTPS(1,use=d), repeat([ComplexTPS(use=d)], 3)...]))

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
      m.x[i][collect(Int, ords)] = a + im*b
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
  m.x[nv+1:nn] .= complexparams(d)


  if !isnothing(spin)
    idx = findfirst(t->t=="c_quaternion", data[:,1])
    if data[idx,3] == "identity"
      setTPS!(m.Q.q0, 1)
      setTPS!(m.Q.q1, 0)
      setTPS!(m.Q.q2, 0)
      setTPS!(m.Q.q3, 0)
    else
      for qi in m.Q
        idx = findfirst(x->(x isa Integer), data[:,1])
        count = 0
        while data[idx,1] >= 0
          a = data[idx,2]
          b = data[idx,3]
          ords = data[idx,4:4+nn-1]
          #println(ords)
          qi[collect(Int,ords)] = a + im*b
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
      m.E[row,col] = num
    end
  end
  #return m


  return m
end


