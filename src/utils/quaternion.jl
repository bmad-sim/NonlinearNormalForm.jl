"""
    Quaternion{T <: Number}

Lightweight quaternion implementation for simulations.
"""
struct Quaternion{T <: Number}
  q::Vector{T}
end

convert(::Type{Quaternion{S}}, Q::Quaternion{T}) where {S,T} = Quaternion(map(x->(S)(x), Q.q))
Quaternion(Q::Quaternion{T}) where T <: Number = Quaternion(map(x->T(x), Q.q))

Quaternion(t::T) where T <: Number = Quaternion([one(t), zero(t), zero(t), zero(t)])
Quaternion(::Nothing) = Quaternion{Nothing}(Nothing[])

==(Q1::Quaternion, Q2::Quaternion) = Q1.q == Q2.q

"""
    mul!(Q3::Quaternion, Q1::Quaternion{S}, Q2::Quaternion{T}) where {S,T}

Multiplies the quaternions `Q1*Q2` and stores the result in-place in `Q3`. 
Aliasing of all three arguments is allowed.
"""
function mul!(Q3::Quaternion, Q1::Quaternion{S}, Q2::Quaternion{T}) where {S,T}
  q1 = Q1.q
  q2 = Q2.q
  A = @FastGTPSA  (q1[1]+q1[2])*(q2[1]+q2[2])
  B = @FastGTPSA  (q1[4]-q1[3])*(q2[3]-q2[4])
  C = @FastGTPSA  (q1[1]-q1[2])*(q2[3]+q2[4]) 
  D = @FastGTPSA  (q1[3]+q1[4])*(q2[1]-q2[2])
  E = @FastGTPSA  (q1[2]+q1[4])*(q2[2]+q2[3])
  F = @FastGTPSA  (q1[2]-q1[4])*(q2[2]-q2[3])
  G = @FastGTPSA  (q1[1]+q1[3])*(q2[1]-q2[4])
  H = @FastGTPSA  (q1[1]-q1[3])*(q2[1]+q2[4])
  
  Q3.q[1] = @FastGTPSA  B+(-E-F+G+H)/2
  Q3.q[2] = @FastGTPSA  A-(E+F+G+H)/2
  Q3.q[3] = @FastGTPSA  C+(E-F+G-H)/2
  Q3.q[4] = @FastGTPSA  D+(E-F-G+H)/2
end

+(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.q .+ Q2.q)
-(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.q .- Q2.q)

function norm(Q1::Quaternion)
  return sqrt(dot(Q1,Q1))
end

function *(Q1::Quaternion, Q2::Quaternion)
  q1 = Q1.q
  q2 = Q2.q
  E = @FastGTPSA (q1[2]+q1[4])*(q2[2]+q2[3])
  F = @FastGTPSA (q1[2]-q1[4])*(q2[2]-q2[3])
  G = @FastGTPSA (q1[1]+q1[3])*(q2[1]-q2[4])
  H = @FastGTPSA (q1[1]-q1[3])*(q2[1]+q2[4])
  
  out1 = @FastGTPSA (q1[4]-q1[3])*(q2[3]-q2[4])+(-E-F+G+H)/2
  out2 = @FastGTPSA (q1[1]+q1[2])*(q2[1]+q2[2])-(E+F+G+H)/2
  out3 = @FastGTPSA (q1[1]-q1[2])*(q2[3]+q2[4])+(E-F+G-H)/2
  out4 = @FastGTPSA (q1[3]+q1[4])*(q2[1]-q2[2])+(E-F-G+H)/2

  return Quaternion([out1,out2,out3,out4])
end

function dot(Q1::Quaternion, Q2::Quaternion)
  q1 = Q1.q
  q2 = Q2.q
  return  conj(q1[1])*q2[1]+conj(q1[2])*q2[2]+conj(q1[3])*q2[3]+conj(q1[4])*q2[4]
end

function inv(Q1::Quaternion)
  q1 = Q1.q
  out1 = q1./dot(Q1,Q1)
  return Quaternion([out1[1], -out1[2], -out1[3], -out1[4]])
end

function inv!(Q::Quaternion, Q1::Quaternion)
  Q.q .= Q1.q./dot(Q1,Q1)
  @inbounds Q.q[2:4] .*= -1
  return
end

function to_SO3(Q1::Quaternion)
  sq1 = Q1.q[1] * Q1.q[1]
  sqx = Q1.q[2] * Q1.q[2]
  sqy = Q1.q[3] * Q1.q[3]
  sqz = Q1.q[4] * Q1.q[4]
  rmat = Matrix{eltype(Q1.q)}(undef,3,3)

  # invs (inverse square length) is only required if quaternion is not already normalised                                       

  invs = 1 / (sqx + sqy + sqz + sq1)
  rmat[1,1] = ( sqx - sqy - sqz + sq1) * invs   # since sq1 + sqx + sqy + sqz =1/invs * invs                                     
  rmat[2,2] = (-sqx + sqy - sqz + sq1) * invs
  rmat[3,3] = (-sqx - sqy + sqz + sq1) * invs

  tmp1 = Q1.q[2] * Q1.q[3]
  tmp2 = Q1.q[4] * Q1.q[1]
  rmat[2,1] = 2 * (tmp1 + tmp2) * invs
  rmat[1,2] = 2 * (tmp1 - tmp2) * invs

  tmp1 = Q1.q[2] * Q1.q[4]
  tmp2 = Q1.q[3] * Q1.q[1]
  rmat[3,1] = 2 * (tmp1 - tmp2) * invs
  rmat[1,3] = 2 * (tmp1 + tmp2) * invs
  tmp1 = Q1.q[3] * Q1.q[4]
  tmp2 = Q1.q[2] * Q1.q[1]
  rmat[3,2] = 2 * (tmp1 + tmp2) * invs
  rmat[2,3] = 2 * (tmp1 - tmp2) * invs
  return rmat
end
