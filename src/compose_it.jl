for t = (:DAMap, :TPSAMap)
@eval begin
  
"""
    compose_it!(m, m2, m1; work_low=nothing, work_prom=nothing)

Low level composition function, `m = m2 ∘ m1`. Aliasing `m` with `m2` is allowed, but not `m` with `m1`.
Assumes the destination map is properly set up (with correct types promoted if necessary), and 
that `m.x[1:nv]` (and `m.Q.q` if spin) contain allocated TPSs. The parameters part of `m.x` (`m.x[nv+1:nn]`) 
does not need to contain allocated TPSs.

For all compositions, 5 temporary vectors must be generated that contain Ptr{RTPSA} or Ptr{CTPSA}
for each TPS in the map (depending on output type), to pass to the low-level C composition function in GTPSA. 
Three vectors are for the orbital part (`m.x`, `m2.x`, `m1.x` correspondingly referred in the code as 
`outx_low`, `m2x_low`, and `m1x_low`) and two are for the spin part (`m.Q.q` and `m2.Q.q` correspondingly 
referred in the code as `outQ_low` and `m2Q_low`).  These 5 temporaries can be optionally passed as a tuple 
in `work_low`, and must satisfy the following requirements:

work_low[1] = outx_low   # Length >= nv
work_low[2] = m2x_low    # Length >= nv
work_low[3] = m1x_low    # Length >= nv+np
work_low[4] = outQ_low   # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4
work_low[5] = m2Q_low    # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4

Furthermore, for the spin part both of `outx_low` and `m2x_low` could be reused for `outQ_low` and `m2Q_low` 
if `nv >= 4` , however `m1x_low` MUST NOT BE REUSED!

If promotion is occuring, then one of the maps is `ComplexTPS` and the other `TPS`, with output map `ComplexTPS`. 
Note that the spin part is required to agree with the orbital part in terms of type by definition of the `TaylorMap` 
struct. `work_prom` can optionally be passed as a tuple containing the temporary vectors if promotion is occuring:

If `eltype(m.x) != eltype(m1.x)` (then `m1` must be promoted):
work_prom[1] = m1x_prom  # Length >= nv+np, could be = Vector{ComplexTPS}(undef, nv+np)

else if `eltype(m.x) != eltype(m2.x)` (then `m2` must be promoted):
work_prom[1] = m2x_prom  # Length >= nv, could be = Vector{ComplexTPS}(undef, nv)
work_prom[2] = m2Q_prom  # Length >= 4, could be = Vector{ComplexTPS}(undef, 4)

Note that the `ComplexTPS`s in the vector do NOT need to be defined - just the container should be passed.
""" 
function compose_it!(m::$t, m2::$t, m1::$t; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))
  @assert !(m === m1) "Cannot compose_it!(m, m2, m1) with m === m1"
  @assert !isnothing(m1.Q) && !isnothing(m2.Q) || m1.Q == m2.Q "Cannot compose: one map includes spin, other does not"
  @assert !isnothing(m1.E) && !isnothing(m2.E) || m1.E == m2.E "Cannot compose: one map includes radiation, other does not"
  @assert eltype(m.x) == promote_type(eltype(m2.x),eltype(m1.x)) "Cannot compose: output map orbital ray type $(eltype(m.x)) must be $(promote_type(eltype(m2.x),eltype(m1.x)))"
  @assert isnothing(m1.Q) && isnothing(m.Q) || eltype(m.Q.q) == promote_type(eltype(m2.Q.q),eltype(m1.Q.q)) "Cannot compose: output map quaternion type $(eltype(m.x)) must be $(promote_type(eltype(m2.x),eltype(m1.x)))"
  @assert isnothing(m1.E) && isnothing(m.E) || eltype(m.E) == promote_type(eltype(m2.E),eltype(m1.E)) "Cannot compose: output map stochastic matrix type $(eltype(m.E)) must be $(promote_type(eltype(m2.E),eltype(m1.E)))"

  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  outT = eltype(m.x)

  outx_low = work_low[1]
  m2x_low = work_low[2]
  m1x_low = work_low[3]
  @assert length(outx_low) >= nv "Incorrect length for work_low[1] = outx_low. Received $(length(outx_low)), should be >=$nv"
  @assert length(m2x_low) >= nv "Incorrect length for work_low[2] = m2x_low. Received $(length(m2x_low)), should be >=$nv"
  @assert length(m1x_low) >= nn "Incorrect length for work_low[3] = m1x_low. Received $(length(m1x_low)), should be >=$nn"
  if !isnothing(m.Q) && dospin
    outQ_low = work_low[4]
    m2Q_low = work_low[5]
    @assert length(outQ_low) >= 4 "Incorrect length for outQ_low: length(work_low[4]) < 4"
    @assert length(m2Q_low) >= 4 "Incorrect length for m2Q_low: length(work_low[5]) < 4"
    @assert !(outQ_low === m1x_low) "m1x_low === outQ_low !! m1x_low must NOT be reused!"
    @assert !(m2Q_low === m1x_low) "m1x_low === m2Q_low !! m1x_low must NOT be reused!"
  end

  if outT != eltype(m1.x)
    m1x_prom = work_prom[1]
    m2x_prom = nothing
    m2Q_prom = nothing
    @assert length(m1x_prom) >= nn "Incorrect length for work_prom[1] = m1x_prom: Received $(length(m1x_prom)), should be >=$nn"
  elseif outT != eltype(m2.x)
    m1x_prom = nothing
    m2x_prom = work_prom[1]
    if isnothing(m.Q) && dospin
      m2Q_prom = work_prom[2]
      @assert length(m2Q_prom) >= 4 "Incorrect length for work_prom[2] = m2Q_prom: Received $(length(m2Q_prom)), should be >=4"
    end
    @assert length(m2x_prom) >= nv "Incorrect length for work_prom[1] = m2x_prom: Received $(length(m2x_prom)), should be >=$nv"
  else
    m1x_prom = nothing
    m2x_prom = nothing
    m2Q_prom = nothing
  end

  # Deal with x0:
  m.x0 .= m1.x0

  # add immutable parameters to outx
  if eltype(m.x) == eltype(m1.x)
    @inbounds m.x[nv+1:nn] = view(m1.x, nv+1:nn)
  else
    @inbounds m.x[nv+1:nn] = view(m2.x, nv+1:nn)
  end

  # Do the composition, promoting if necessary
  if outT != eltype(m1.x) 
    # Promote to ComplexTPS:
    map!(t->ComplexTPS(t), m1x_prom,  m1.x)   
    map!(t->t.tpsa, m1x_low, m1x_prom)
  else
    map!(t->t.tpsa, m1x_low, m1.x)
  end

  if outT != eltype(m2.x)
    map!(t->ComplexTPS(t), m2x_prom, view(m2.x,1:nv))  
    map!(t->t.tpsa, m2x_low, m2x_prom)
  else
    map!(t->t.tpsa, m2x_low, view(m2.x, 1:nv))
  end

  # go low
  map!(t->t.tpsa, outx_low, view(m.x,1:nv))

  # Orbit:
  GC.@preserve m1x_prom m2x_prom compose!(nv, m2x_low, nv+np, m1x_low, outx_low)

  # Spin:
  if !isnothing(m.Q) && dospin
    if outT != eltype(m2.x)
      map!(t->ComplexTPS(t), m2Q_prom, m2.Q.q) 
      map!(t->t.tpsa, m2Q_low, m2Q_prom)
    else
      map!(t->t.tpsa, m2Q_low, m2.Q.q)
    end
    # Go low
    map!(t->t.tpsa, outQ_low, m.Q.q)
    # Spin (spectator) q(z0)=q2(M(z0))q1(z0)
    # First obtain q2(M(z0))
    GC.@preserve m1x_prom m2Q_prom compose!(Cint(4), m2Q_low, nv+np, m1x_low, outQ_low)
    # Now concatenate
    mul!(m.Q, m1.Q, m.Q)
  end

  return 
end

end
end