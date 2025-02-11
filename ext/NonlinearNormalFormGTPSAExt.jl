module NonlinearNormalFormGTPSAExt
import NonlinearNormalForm as NNF
using NonlinearNormalForm: Quaternion
using GTPSA: GTPSA, @FastGTPSA, @FastGTPSA!, TPS

function NNF.mul!(q::Quaternion{<:TPS}, q1::Quaternion{<:TPS}, q2::Quaternion{<:TPS})
  @assert !(q === q1) && !(q === q2) "Aliasing q with q1 or q2 not allowed!"
  @FastGTPSA! begin
  q.q0 = q1.q0 * q2.q0 - q1.q1 * q2.q1 - q1.q2 * q2.q2 - q1.q3 * q2.q3
  q.q1 = q1.q0 * q2.q1 + q1.q1 * q2.q0 + q1.q2 * q2.q3 - q1.q3 * q2.q2
  q.q2 = q1.q0 * q2.q2 - q1.q1 * q2.q3 + q1.q2 * q2.q0 + q1.q3 * q2.q1
  q.q3 = q1.q0 * q2.q3 + q1.q1 * q2.q2 - q1.q2 * q2.q1 + q1.q3 * q2.q0
  end
  return q
end

function NNF.:*(q1::Quaternion{<:TPS}, q2::Quaternion{<:TPS})
  @FastGTPSA begin
  q0 = q1.q0 * q2.q0 - q1.q1 * q2.q1 - q1.q2 * q2.q2 - q1.q3 * q2.q3
  q1 = q1.q0 * q2.q1 + q1.q1 * q2.q0 + q1.q2 * q2.q3 - q1.q3 * q2.q2
  q2 = q1.q0 * q2.q2 - q1.q1 * q2.q3 + q1.q2 * q2.q0 + q1.q3 * q2.q1
  q3 = q1.q0 * q2.q3 + q1.q1 * q2.q2 - q1.q2 * q2.q1 + q1.q3 * q2.q0
  end
  return Quaternion(q0,q1,q2,q3)
end

NNF.norm(q1::Quaternion{<:TPS}) = @FastGTPSA sqrt(q1.q0^2 + q1.q1^2 + q1.q2^2 + q1.q3^2)

function NNF.inv!(q::Quaternion{<:TPS}, q1::Quaternion{<:TPS})
  @assert !(q === q1) "Aliasing q with q1 not allowed!"
  @FastGTPSA! begin
    q.q0 = q1.q0 * q1.q0 + q1.q1 * q1.q1 + q1.q2 * q1.q2 + q1.q3 * q1.q3
    q.q1 = -q1.q1/q.q0
    q.q2 = -q1.q2/q.q0
    q.q3 = -q1.q3/q.q0
    q.q0 = q1.q0/q.q0
  end
  return
end

function NNF.inv(q1::Quaternion{<:TPS})
  @FastGTPSA begin
    q0 = q1.q0 * q1.q0 + q1.q1 * q1.q1 + q1.q2 * q1.q2 + q1.q3 * q1.q3
    q1 = -q1.q1/q0
    q2 = -q1.q2/q0
    q3 = -q1.q3/q0
  end
  @FastGTPSA! q0 = q1.q0/q0
  return Quaternion(q0,q1,q2,q3)
end

function NNF.show(io::IO, q::Quaternion{<:TPS})
  println(io, "$(typeof(q)):")
  GTPSA.show_map!(io, collect(q), Ref{Int}(1), false, [" q0:"," q1:"," q2:"," q3:"])
end

NNF.show(io::IO, ::MIME"text/plain", q::Quaternion{<:TPS}) = (println(io, "$(typeof(q)):"); GTPSA.show_map!(io, collect(q), Ref{Int}(1), false, [" q0:"," q1:"," q2:"," q3:"]))


end