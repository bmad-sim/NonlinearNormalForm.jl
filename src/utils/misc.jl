# --- helper functions to get the numvars/numparams/descriptor in any context ---

# GTPSA provides these functions for only pure TPS/ComplexTPSs
getdesc(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap,VectorField}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d))
numvars(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap,VectorField}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nv
numparams(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap,VectorField}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).np
numnn(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap,VectorField}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nn

getdesc(t::Vector{<:Union{TPS,ComplexTPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d))
numvars(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nv
numparams(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).np
numnn(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nn

getdesc(m::Quaternion{<:Union{ComplexTPS,TPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d))
numvars(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nv
numparams(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).np
numnn(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nn

numtype(m::Union{Probe{<:Real,T,<:Any,<:Any},<:TaylorMap{<:Number,T,<:Any,<:Any},VectorField{T,<:Any}}) where {T<:Union{TPS,ComplexTPS}} = numtype(T)


# --- random symplectic map ---
function rand(t::Union{Type{DAMap},Type{TPSAMap}}; spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TPS,ComplexTPS}=GTPSA.desc_current)
  if isnothing(spin)
    U = Nothing
  else
    if spin
      U = Quaternion{TPS}
    else
      U = Nothing
    end
  end

  if isnothing(radiation)
    V = Nothing
  else
    if spin
      V = Float64
    else
      V = Nothing
    end
  end

  return rand(t{Float64,TPS,U,V},use=use)
end


"""
Generate map symplectic up to order in Descriptor
"""
function rand(t::Union{Type{DAMap{S,T,U,V}},Type{TPSAMap{S,T,U,V}}}; require_stable::Bool=false, use::Union{Descriptor,TPS,ComplexTPS}=GTPSA.desc_current) where {S,T,U,V}
  if require_stable # make hamiltonian in phasors basis then reverse
    
    
  end
  
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
  
  # If require_stable, we need to go into phasors basis



  h = T(use=dtmp)
  len = length(h)

  #  random coefficients for hamiltonian except 0th and 1st order terms
  for i=nv+1:len-1
    h[i] = rand(numtype(T))
  end

  F = VectorField{T,U}(h)
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






function read_fpp_map(file)
  data = readdlm(file, skipblanks=true)
  nv = data[findfirst(t->t=="Dimensional", data)- CartesianIndex(0,1)]
  no = data[findfirst(t->t=="NO", data) + CartesianIndex(0,2)]
  np = data[findfirst(t->t=="NV", data) + CartesianIndex(0,2)] - nv
  nn = nv+np
  # Make the TPSA
  d = Descriptor(nv, no, np, no)
  m = complex(DAMap()) #repeat([ComplexTPS(use=d)], nv), Q=Quaternion([ComplexTPS(1,use=d), repeat([ComplexTPS(use=d)], 3)...]))

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
  return m


  # spin?
  idx = findfirst(t->t=="c_quaternion", data[:,1])
  if idx > 0
    if data[idx,3] == "identity"
      m.Q.q .= [1.0, 0., 0., 0.]
    else
      for i=1:4
        idx = findfirst(x->(x isa Integer), data[:,1])
        count = 0
        while data[idx,1] >= 0
          a = data[idx,2]
          b = data[idx,3]
          ords = data[idx,4:end]
          m.Q.q[i][collect(Int,ords)] = a + im*b
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
  end
  return m
end


