# helper functions
getdesc(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d))
numvars(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nv
numparams(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).np
numnn(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nn

getdesc(t::Union{TPS,ComplexTPS}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d))
numvars(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).nv
numparams(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).np
numnn(t::Union{TPS,ComplexTPS}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(t.tpsa).d)).nn

getdesc(t::Vector{<:Union{TPS,ComplexTPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d))
numvars(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nv
numparams(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).np
numnn(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nn

getdesc(m::Quaternion{<:Union{ComplexTPS,TPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d))
numvars(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nv
numparams(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).np
numnn(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nn

getdesc(d::Descriptor) = d
numvars(d::Descriptor) = unsafe_load(d.desc).nv
numparams(d::Descriptor) = unsafe_load(d.desc).np
numnn(d::Descriptor) = unsafe_load(d.desc).nn

getdesc(d::Nothing) = nothing

numtype(t::TPS) = Float64
numtype(ct::ComplexTPS) = ComplexF64
numtype(::Type{TPS}) = Float64
numtype(::Type{ComplexTPS}) = ComplexF64

lowtype(t::TPS) = Ptr{RTPSA}
lowtype(ct::ComplexTPS) = Ptr{CTPSA}
lowtype(::Type{TPS}) = Ptr{RTPSA}
lowtype(::Type{ComplexTPS}) = Ptr{CTPSA}


function checksymp(M::Matrix{T}) where T<:Number
  s = size(M)
  nv = first(s)
  nv == last(s) || error("Non-square matrix!")
  iseven(nv) || error("Matrix contains odd number of rows/columns!")
  J = zeros(nv,nv)
  for i=1:2:nv
    J[i:i+1,i:i+1] = [0 1; -1 0];
  end
  res = transpose(M)*J*M-J
  return res
end

jacobian(m::TaylorMap;include_params=false) = jacobian(m.x[1:numvars(m)],include_params=include_params)
checksymp(m::TaylorMap) = checksymp(jacobian(m))

function read_fpp_map(file)
  data = readdlm(file, skipblanks=true)
  nv = data[findfirst(t->t=="Dimensional", data)- CartesianIndex(0,1)]
  no = data[findfirst(t->t=="NO", data) + CartesianIndex(0,2)]
  np = data[findfirst(t->t=="NV", data) + CartesianIndex(0,2)] - nv
  nn = nv+np
  # Make the TPSA
  d = Descriptor(nv, no, np, no)
  m = DAMap(repeat([ComplexTPS(use=d)], nv), Q=Quaternion([ComplexTPS(1,use=d), repeat([ComplexTPS(use=d)], 3)...]))

  idx=3
  data=data[3:end,:]
  # Now fill
  for i=1:nv
    idx = findfirst(x->(x isa Integer), data[:,1])
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
  # dont forget params
  m.x[nv+1:nn] = complexparams(d)


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
  end
  return m
end