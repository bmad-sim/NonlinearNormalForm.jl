#= 
Simple differentiable (bi)quaternion bc Julia 
doesn't seem to have anything...

Or Julia does have stuff but it requires 
Quaternion{T<:Real} + im*Quaternion{T<:Real} 
to fit in with the stupid concrete Complex type
presumably.
=#
struct Quaternion{T <: Number}
  q::Vector{T}
end

convert(::Type{Quaternion{S}}, Q::Quaternion{T}) where {S,T} = Quaternion(map(x->(S)(x), Q.q))
Quaternion(Q::Quaternion) = Quaternion(deepcopy(Q.q))

Quaternion(t::T) where T <: Number = Quaternion([one(t), zero(t), zero(t), zero(t)])

function qmul!(Q1::Quaternion{S}, Q2::Quaternion{T}, Q3::Quaternion) where {S,T}
  q1 = Q1.q
  q2 = Q2.q
  A = @FastGTPSA (q1[1]+q1[2])*(q2[1]+q2[2])
  B = @FastGTPSA (q1[4]-q1[3])*(q2[3]-q2[4])
  C = @FastGTPSA (q1[1]-q1[2])*(q2[3]+q2[4]) 
  D = @FastGTPSA (q1[3]+q1[4])*(q2[1]-q2[2])
  E = @FastGTPSA (q1[2]+q1[4])*(q2[2]+q2[3])
  F = @FastGTPSA (q1[2]-q1[4])*(q2[2]-q2[3])
  G = @FastGTPSA (q1[1]+q1[3])*(q2[1]-q2[4])
  H = @FastGTPSA (q1[1]-q1[3])*(q2[1]+q2[4])
  
  Q3.q[1] = @FastGTPSA B+(-E-F+G+H)/2
  Q3.q[2] = @FastGTPSA A-(E+F+G+H)/2
  Q3.q[3] = @FastGTPSA C+(E-F+G-H)/2
  Q3.q[4] = @FastGTPSA D+(E-F-G+H)/2
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
  A = @FastGTPSA (q1[1]+q1[2])*(q2[1]+q2[2])
  B = @FastGTPSA (q1[4]-q1[3])*(q2[3]-q2[4])
  C = @FastGTPSA (q1[1]-q1[2])*(q2[3]+q2[4]) 
  D = @FastGTPSA (q1[3]+q1[4])*(q2[1]-q2[2])
  E = @FastGTPSA (q1[2]+q1[4])*(q2[2]+q2[3])
  F = @FastGTPSA (q1[2]-q1[4])*(q2[2]-q2[3])
  G = @FastGTPSA (q1[1]+q1[3])*(q2[1]-q2[4])
  H = @FastGTPSA (q1[1]-q1[3])*(q2[1]+q2[4])
  
  out1 = @FastGTPSA B+(-E-F+G+H)/2
  out2 = @FastGTPSA A-(E+F+G+H)/2
  out3 = @FastGTPSA C+(E-F+G-H)/2
  out4 = @FastGTPSA D+(E-F-G+H)/2

  return Quaternion([out1,out2,out3,out4])
end

function LinearAlgebra.dot(Q1::Quaternion, Q2::Quaternion)
  q1 = Q1.q
  q2 = Q2.q
  return @FastGTPSA conj(q1[1])*q2[1]+conj(q1[2])*q2[2]+conj(q1[3])*q2[3]+conj(q1[4])*q2[4]
end

function Base.inv(Q1::Quaternion)
  q1 = Q1.q
  out1 = q1./dot(Q1,Q1)
  return Quaternion([out1[1], -out1[2], -out1[3], -out1[4]])
end