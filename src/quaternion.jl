#= 
Simple differentiable (bi)quaternion bc Julia 
doesn't seem to have anything...

Or Julia does have stuff but it requires 
Quaternion{T<:Real} + im*Quaternion{T<:Real} 
to fit in with the concrete Complex type
presumably.
=#
struct Quaternion{T <: Number}
  q::Vector{T}
end

convert(::Type{Quaternion{S}}, Q::Quaternion{T}) where {S,T} = Quaternion(map(x->(S)(x), Q1.q))
Quaternion(Q::Quaternion) = Quaternion(deepcopy(Q1.q))

Quaternion(t::T) where T <: Number = Quaternion([one(t), zero(t), zero(t), zero(t)])

function qmul!(Q1::Quaternion{S}, Q2::Quaternion{T}, Q3::Quaternion) where {S,T}
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
  
  Q3.q[1] =  B+(-E-F+G+H)/2
  Q3.q[2] =  A-(E+F+G+H)/2
  Q3.q[3] =  C+(E-F+G-H)/2
  Q3.q[4] =  D+(E-F-G+H)/2
end

function Base.:+(Q1::Quaternion, Q2::Quaternion)
  return Quaternion(Q1.q .+ Q2.q)
end

function LinearAlgebra.norm(Q1::Quaternion)
  return sqrt(dot(Q1,Q1))
end

function Base.:*(Q1::Quaternion, Q2::Quaternion)
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
  
  out1 =  B+(-E-F+G+H)/2
  out2 =  A-(E+F+G+H)/2
  out3 =  C+(E-F+G-H)/2
  out4 =  D+(E-F-G+H)/2

  return Quaternion([out1,out2,out3,out4])
end

function LinearAlgebra.dot(Q1::Quaternion, Q2::Quaternion)
  q1 = Q1.q
  q2 = Q2.q
  return  conj(q1[1])*q2[1]+conj(q1[2])*q2[2]+conj(q1[3])*q2[3]+conj(q1[4])*q2[4]
end

function Base.inv(Q1::Quaternion)
  q1 = Q1.q
  out1 = q1./dot(Q1,Q1)
  return Quaternion([out1[1], -out1[2], -out1[3], -out1[4]])
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