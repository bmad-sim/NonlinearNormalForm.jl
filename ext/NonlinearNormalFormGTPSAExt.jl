module NonlinearNormalFormGTPSAExt
import NonlinearNormalForm as NNF
using NonlinearNormalForm: Quaternion
using GTPSA: GTPSA, @FastGTPSA, @FastGTPSA!, TPS, Descriptor

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
  qf0 = q1.q0 * q2.q0 - q1.q1 * q2.q1 - q1.q2 * q2.q2 - q1.q3 * q2.q3
  qf1 = q1.q0 * q2.q1 + q1.q1 * q2.q0 + q1.q2 * q2.q3 - q1.q3 * q2.q2
  qf2 = q1.q0 * q2.q2 - q1.q1 * q2.q3 + q1.q2 * q2.q0 + q1.q3 * q2.q1
  qf3 = q1.q0 * q2.q3 + q1.q1 * q2.q2 - q1.q2 * q2.q1 + q1.q3 * q2.q0
  end
  return Quaternion(qf0,qf1,qf2,qf3)
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
    qf0 = q1.q0 * q1.q0 + q1.q1 * q1.q1 + q1.q2 * q1.q2 + q1.q3 * q1.q3
    qf1 = -q1.q1/qf0
    qf2 = -q1.q2/qf0
    qf3 = -q1.q3/qf0
  end
  @FastGTPSA! qf0 = q1.q0/qf0
  return Quaternion(qf0,qf1,qf2,qf3)
end

NNF.zero(q1::Quaternion{<:TPS}) = return Quaternion(zero(q1.q0), zero(q1.q1), zero(q1.q2), zero(q1.q3))
NNF.one(q1::Quaternion{<:TPS}) = return Quaternion(one(q1.q0), zero(q1.q1), zero(q1.q2), zero(q1.q3))


function show_quat(io, m::Quaternion{<:TPS{T,D}}) where {T,D}
  lines_used = 0
  println(io, typeof(m), ":") 
  lines_used += 1
  desc = first(m).d
  diffdescs = false
  for i in 1:4
    if !diffdescs && desc != m[i].d
      println(io, "WARNING: Atleast one $(eltype(m)) has a different Descriptor!")
      diffdescs = true
      lines_used += 1
    end
  end
  if !diffdescs && D == GTPSA.Dynamic
    println(io, Descriptor(desc))
    lines_used += 1
  end
  tpsouts = Any[]
  coef_str = "INDEX  COEFFICIENT             ORDER   EXPONENTS"
  println(io, coef_str)
  lines_used += 1
  oversized=false
  #@show lines_used
  for i in 1:4
    lines_used += 1 # For the line ----
    if diffdescs # for the Descriptor
      lines_used += 1
    end
    t = m[i]
    lines = GTPSA.format_tps_contents(t)
    lines = map(lines[2:end-1]) do line
      GTPSA.@sprintf(" q%i:  %s", i-1, line)
    end

    if get(io, :limit, false) && lines_used + length(lines) > displaysize(io)[1]-4
      oversized=true
      #@show lines_used
      #@show lines_used-(displaysize(io)[1]-6)
      #@show length(lines)
     ## @show min(lines_used-(displaysize(io)[1]-6), length(lines))
      lines = lines[1:min(displaysize(io)[1]-4-lines_used, length(lines))] #lines_used-displaysize(io)[1]-4]
      if !isempty(lines) 
        push!(tpsouts, join(lines, '\n')) 
      end
      break
    else
      lines_used += length(lines)
    end
    push!(tpsouts, join(lines, '\n')) 
  end
  
  line_length = findmax(tpsouts) do tpsout
    findmax(length.(split(tpsout, '\n')))[1]
  end
  
  for (tps,tpsout) in zip(m,tpsouts)
    println(io, repeat('-', max(line_length[1], length(coef_str))))
    if diffdescs
      println(io, Descriptor(tps.d))
      #println(io, Descriptor(desc))
    end
    println(io, tpsout)
  end
  if oversized
    print(io, "       ... (Output truncated)")
    return
  end
  return
end
  
NNF.show(io::IO, q::Quaternion{<:TPS}) = show_quat(io, q)
NNF.show(io::IO, ::MIME"text/plain", q::Quaternion{<:TPS}) = show_quat(io, q) # (println(io, "$(typeof(q)):"); GTPSA.show_map!(io, collect(q), Ref{Int}(1), false, [" q0:"," q1:"," q2:"," q3:"]))


end