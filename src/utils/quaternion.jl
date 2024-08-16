"""
    Quaternion{T <: Number}

Lightweight quaternion implementation for simulations.
Mutable in order to compose in GTPSA.
"""
mutable struct Quaternion{T <: Number} <: FieldVector{4, T}
  q0::T
  q1::T
  q2::T
  q3::T
end

# This file provides special fast functions for Quaternion{<:TPS}

"""
    mul!(Q::Quaternion, Q1::Quaternion, Q2::Quaternion)

Sets `Q` equal to `Q1*Q2`, with zero allocations. `Q` must not be 
aliased with either `Q1` or `Q2`.
"""
function mul!(Q::Quaternion{<:TPS}, Q1::Quaternion{<:TPS}, Q2::Quaternion{<:TPS})
  @assert !(Q === Q1) && !(Q === Q2) "Aliasing Q with Q1 or Q2 not allowed!"
  @FastGTPSA! begin
  Q.q0 = Q1.q0 * Q2.q0 - Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3
  Q.q1 = Q1.q0 * Q2.q1 + Q1.q1 * Q2.q0 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2
  Q.q2 = Q1.q0 * Q2.q2 - Q1.q1 * Q2.q3 + Q1.q2 * Q2.q0 + Q1.q3 * Q2.q1
  Q.q3 = Q1.q0 * Q2.q3 + Q1.q1 * Q2.q2 - Q1.q2 * Q2.q1 + Q1.q3 * Q2.q0
  end
  return Q
end

function *(Q1::Quaternion, Q2::Quaternion)
  @FastGTPSA begin
  q0 = Q1.q0 * Q2.q0 - Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3
  q1 = Q1.q0 * Q2.q1 + Q1.q1 * Q2.q0 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2
  q2 = Q1.q0 * Q2.q2 - Q1.q1 * Q2.q3 + Q1.q2 * Q2.q0 + Q1.q3 * Q2.q1
  q3 = Q1.q0 * Q2.q3 + Q1.q1 * Q2.q2 - Q1.q2 * Q2.q1 + Q1.q3 * Q2.q0
  end
  return Quaternion(q0,q1,q2,q3)
end

norm(Q1::Quaternion{<:TPS}) = @FastGTPSA sqrt(Q1.q0^2 + Q1.q1^2 + Q1.q2^2 + Q1.q3^2)

function inv!(Q::Quaternion{<:TPS}, Q1::Quaternion{<:TPS})
  @assert !(Q === Q1) "Aliasing Q with Q1 not allowed!"
  @FastGTPSA! begin
  Q.q0 = Q1.q0 * Q1.q0 + Q1.q1 * Q1.q1 + Q1.q2 * Q1.q2 + Q1.q3 * Q1.q3
  Q.q1 = -Q1.q1/Q.q0
  Q.q2 = -Q1.q2/Q.q0
  Q.q3 = -Q1.q3/Q.q0
  Q.q0 = Q1.q0/Q.q0
  end
  return
end

function inv(Q1::Quaternion)
  @FastGTPSA begin
  q0 = Q1.q0 * Q1.q0 + Q1.q1 * Q1.q1 + Q1.q2 * Q1.q2 + Q1.q3 * Q1.q3
  q1 = -Q1.q1/q0
  q2 = -Q1.q2/q0
  q3 = -Q1.q3/q0
  end
  @FastGTPSA! q0 = Q1.q0/q0
  return Quaternion(q0,q1,q2,q3)
end

# This needs to be optimized
function exp(Q1::Quaternion{T}) where {T<:TPS}
  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*4
  nrm_ =Inf
  conv = false
  slow = false

  Qout = Quaternion{T}(1,0,0,0)
  Q = Quaternion{T}(1,0,0,0)
  for j=1:nmax
    Q = Q*Q1/j
    Qout .+= Q

    nrm = norm(normTPS(Q.q0) + normTPS(Q.q1) + normTPS(Q.q2) + normTPS(Q.q2))
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      if slow
        @warn "exp! slow convergence: required n = $(j) iterations"
      end
      return Qout
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  @warn "exp convergence not reached for $nmax iterations"
  return 
end

function show(io::IO, Q::Quaternion{<:TPS})
  GTPSA.show_map!(io, collect(Q), Ref{Int}(0), false, [" q0:"," q1:"," q2:"," q3:"])
end

show(io::IO, ::MIME"text/plain", Q::Quaternion{<:TPS}) = GTPSA.show_map!(io, collect(Q), Ref{Int}(0), false, [" q0:"," q1:"," q2:"," q3:"])



#convert(::Type{Quaternion{S}}, Q::Quaternion{T}) where {S,T} = Quaternion(map(x->(S)(x), Q.q))
#Quaternion(Q::Quaternion{T}) where T <: Number = Quaternion(map(x->T(x), Q.q))
#Quaternion(q0, q1, q2, q3) = Quaternion(MVector{4}(promote(q0, q1, q2, q3)))
#Quaternion(q::AbstractVector{<:Number}) = Quaternion(MVector{4}(q))

#Quaternion(t::T) where T <: Number = Quaternion([one(t), zero(t), zero(t), zero(t)])
#Quaternion(::Nothing) = Quaternion{Nothing}(Nothing[])

#eltype(::Type{Quaternion{T}}) where T = T
#eltype(q::Quaternion{T}) where T = T


#=

function mul!(Q::Quaternion, Q1::Quaternion, Q2::Quaternion)
  q0 = Q1.q0 * Q2.q0 - Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3
  q1 = Q1.q0 * Q2.q1 + Q1.q1 * Q2.q0 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2
  q2 = Q1.q0 * Q2.q2 - Q1.q1 * Q2.q3 + Q1.q2 * Q2.q0 + Q1.q3 * Q2.q1
  q3 = Q1.q0 * Q2.q3 + Q1.q1 * Q2.q2 - Q1.q2 * Q2.q1 + Q1.q3 * Q2.q0

  copy!(Q.q0, q0)
  copy!(Q.q1, q1)
  copy!(Q.q2, q2)
  copy!(Q.q3, q3)
  return Q
end
function *(Q1::Quaternion, Q2::Quaternion)
  q0 = Q1.q0 * Q2.q0 - Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3
  q1 = Q1.q0 * Q2.q1 + Q1.q1 * Q2.q0 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2
  q2 = Q1.q0 * Q2.q2 - Q1.q1 * Q2.q3 + Q1.q2 * Q2.q0 + Q1.q3 * Q2.q1
  q3 = Q1.q0 * Q2.q3 + Q1.q1 * Q2.q2 - Q1.q2 * Q2.q1 + Q1.q3 * Q2.q0
  return Quaternion(q0,q1,q2,q3)
end

function inv(Q1::Quaternion)
  nrm = Q1.q0 * Q1.q0 + Q1.q1 * Q1.q1 + Q1.q2 * Q1.q2 + Q1.q3 * Q1.q3
  q0 = Q1.q0/nrm
  q1 = -Q1.q1/nrm
  q2 = -Q1.q2/nrm
  q3 = -Q1.q3/nrm
  return Quaternion(q0,q1,q2,q3)
end

function inv!(Q::Quaternion, Q1::Quaternion)
  nrm = Q1.q0 * Q1.q0 + Q1.q1 * Q1.q1 + Q1.q2 * Q1.q2 + Q1.q3 * Q1.q3

  div!(Q.q0, Q1.q0, nrm)
  div!(Q.q1, Q1.q1, nrm)
  div!(Q.q2, Q1.q2, nrm)
  div!(Q.q3, Q1.q3, nrm)

  mul!(Q.q1, Q.q1, -1)
  mul!(Q.q2, Q.q2, -1)
  mul!(Q.q3, Q.q3, -1)
  return
end

=#


#mad_compose!(na, ma::Quaternion{TPS{Float64}},    nb, mb::AbstractVector{TPS{Float64}},    mc::Quaternion{TPS{Float64}}) = GTPSA.mad_tpsa_compose!(Cint(na), ma, Cint(nb), mb, mc)
#mad_compose!(na, ma::Quaternion{TPS{ComplexF64}}, nb, mb::AbstractVector{TPS{ComplexF64}}, mc::Quaternion{TPS{ComplexF64}}) = GTPSA.mad_ctpsa_compose!(Cint(na), ma, Cint(nb), mb, mc)
#Base.unsafe_convert(::Type{Ptr{TPS{T}}}, Q::Quaternion{<:TPS{<:T}}) where {T} = unsafe_convert(Ptr{TPS{T}}, Ref(Q))


#=
## THIS IS DEPRECATED! WE NOW USE 
# ReferenceFrameRotations.jl Quaternion implementation



"""
    mul!(Q3::Quaternion, Q1::Quaternion{S}, Q2::Quaternion{T}) where {S,T}

Multiplies the quaternions `Q1*Q2` and stores the result in-place in `Q3`. 
Aliasing of all three arguments is allowed.
"""
function mul!(Q3::Quaternion, Q1::Quaternion{S}, Q2::Quaternion{T}) where {S,T}
  q1 = Q1.q
  q2 = Q2.q
  A =  (q1[1]+q1[2])*(q2[1]+q2[2])
  B =  (q1[4]-q1[3])*(q2[3]-q2[4])
  C =  (q1[1]-q1[2])*(q2[3]+q2[4]) 
  D =  (q1[3]+q1[4])*(q2[1]-q2[2])
  E =  (q1[2]+q1[4])*(q2[2]+q2[3])
  F =  (q1[2]-q1[4])*(q2[2]-q2[3])
  G =  (q1[1]+q1[3])*(q2[1]-q2[4])
  H =  (q1[1]-q1[3])*(q2[1]+q2[4])
  
  Q3.q0 =  B+(-E-F+G+H)/2
  Q3.q1 =  A-(E+F+G+H)/2
  Q3.q2 =  C+(E-F+G-H)/2
  Q3.q3 =  D+(E-F-G+H)/2
end

+(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.q .+ Q2.q)
-(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.q .- Q2.q)



function *(Q1::Quaternion, Q2::Quaternion)
  q1 = Q1.q
  q2 = Q2.q
  E = (q1[2]+q1[4])*(q2[2]+q2[3])
  F = (q1[2]-q1[4])*(q2[2]-q2[3])
  G = (q1[1]+q1[3])*(q2[1]-q2[4])
  H = (q1[1]-q1[3])*(q2[1]+q2[4])
  
  out1 = (q1[4]-q1[3])*(q2[3]-q2[4])+(-E-F+G+H)/2
  out2 = (q1[1]+q1[2])*(q2[1]+q2[2])-(E+F+G+H)/2
  out3 = (q1[1]-q1[2])*(q2[3]+q2[4])+(E-F+G-H)/2
  out4 = (q1[3]+q1[4])*(q2[1]-q2[2])+(E-F-G+H)/2

  return Quaternion([out1,out2,out3,out4])
end


function to_SO3(Q1::Quaternion)
  sq1 = Q1.q0 * Q1.q0
  sqx = Q1.q1 * Q1.q1
  sqy = Q1.q2 * Q1.q2
  sqz = Q1.q3 * Q1.q3
  rmat = Matrix{eltype(Q1.q)}(undef,3,3)

  # invs (inverse square length) is only required if quaternion is not already normalised                                       

  invs = 1 / (sqx + sqy + sqz + sq1)
  rmat[1,1] = ( sqx - sqy - sqz + sq1) * invs   # since sq1 + sqx + sqy + sqz =1/invs * invs                                     
  rmat[2,2] = (-sqx + sqy - sqz + sq1) * invs
  rmat[3,3] = (-sqx - sqy + sqz + sq1) * invs

  tmp1 = Q1.q1 * Q1.q2
  tmp2 = Q1.q3 * Q1.q0
  rmat[2,1] = 2 * (tmp1 + tmp2) * invs
  rmat[1,2] = 2 * (tmp1 - tmp2) * invs

  tmp1 = Q1.q1 * Q1.q3
  tmp2 = Q1.q2 * Q1.q0
  rmat[3,1] = 2 * (tmp1 - tmp2) * invs
  rmat[1,3] = 2 * (tmp1 + tmp2) * invs
  tmp1 = Q1.q2 * Q1.q3
  tmp2 = Q1.q1 * Q1.q0
  rmat[3,2] = 2 * (tmp1 + tmp2) * invs
  rmat[2,3] = 2 * (tmp1 - tmp2) * invs
  return rmat
end
=#