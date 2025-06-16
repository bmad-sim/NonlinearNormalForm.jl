#=

Special quaternion routines for NonlinearNormalForm

=#

# TO-DO: use MQuaternion (mutable quaternion) when ready
function TI.compose!(q::Quaternion, q1::Quaternion, m1::AbstractArray)
  qt = SA[q.q0, q.q1, q.q2, q.q3]
  q1 = SA[q1.q0, q1.q1, q1.q2, q1.q3]
  TI.compose!(qt, q1, m1)
  return q
end

function mul!(q::Quaternion, q1::Quaternion, q2::Quaternion)
  @assert !(q === q1) && !(q === q2) "Aliasing q with q1 or q2 not allowed!"
  all(qi->TI.IsTPSType(eltype(qi)) isa TI.IsTPSType) || error("mul! only works on TPS types supported by TPSAInterface.jl!")
  TI.copy!(q.q0, q1.q0 * q2.q0 - q1.q1 * q2.q1 - q1.q2 * q2.q2 - q1.q3 * q2.q3)
  TI.copy!(q.q1, q1.q0 * q2.q1 + q1.q1 * q2.q0 + q1.q2 * q2.q3 - q1.q3 * q2.q2)
  TI.copy!(q.q2, q1.q0 * q2.q2 - q1.q1 * q2.q3 + q1.q2 * q2.q0 + q1.q3 * q2.q1)
  TI.copy!(q.q3, q1.q0 * q2.q3 + q1.q1 * q2.q2 - q1.q2 * q2.q1 + q1.q3 * q2.q0)
  return q
end

function inv!(q::Quaternion, q1::Quaternion)
  @assert !(q === q1) "Aliasing q with q1 not allowed!"
  all(qi->TI.IsTPSType(eltype(qi)) isa TI.IsTPSType) || error("inv! only works on TPS types supported by TPSAInterface.jl!")
  TI.copy!(q.q0, q1.q0 * q1.q0 + q1.q1 * q1.q1 + q1.q2 * q1.q2 + q1.q3 * q1.q3)
  TI.copy!(q.q1, -q1.q1/q.q0)
  TI.copy!(q.q2, -q1.q2/q.q0)
  TI.copy!(q.q3, -q1.q3/q.q0)
  TI.copy!(q.q0, q1.q0/q.q0)
  return
end

# TO-DO: Optimize
function exp(q1::Quaternion{T}) where {T}
  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*4
  nrm_ =Inf
  conv = false
  slow = false

  q_out = Quaternion{T}(1,0,0,0)
  q = Quaternion{T}(1,0,0,0)
  for j in 1:nmax
    q = q*q1/j
    q_out += q

    nrm = norm(TI.norm_tps(q.q0) + TI.norm_tps(q.q1) + TI.norm_tps(q.q2) + TI.norm_tps(q.q3))
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      if slow
        @warn "exp! slow convergence: required n = $(j) iterations"
      end
      return q_out
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  @warn "exp convergence not reached for $nmax iterations"
  return 
end
