for t = (:DAMap, :TPSAMap)
@eval begin
  
"""
    compose_it!(m, m2, m1; dospin::Bool=true, dostochastic::Bool=true, work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))

Low level composition function, `m = m2 ∘ m1`. Aliasing `m` with `m2` is allowed, but not `m` with `m1`.
Assumes the destination map is properly set up (with correct types promoted if necessary), and 
that `m.x[1:nv]` (and `m.Q` if spin) contain allocated TPSs. The parameters part of `m.x` (`m.x[nv+1:nn]`) 
does not need to contain allocated TPSs.

If promotion is occuring, then one of the maps is `ComplexTPS64` and the other `TPS`, with output map `ComplexTPS64`. 
Note that the spin part is required to agree with the orbital part in terms of type by definition of the `TaylorMap` 
struct. `work_prom` can optionally be passed as a tuple containing the temporary `ComplexTPS64`s if promotion is occuring:

If `eltype(m.x) != eltype(m1.x)` (then `m1` must be promoted):
`work_prom[1] = m1x_prom  # Length >= nv+np, Vector{ComplexTPS64}`

else if `eltype(m.x) != eltype(m2.x)` (then `m2` must be promoted):
```
work_prom[1] = m2x_prom  # Length >= nv, Vector{ComplexTPS64}
work_prom[2] = m2Q_prom  # Length >= 4, Vector{ComplexTPS64}
```
Note that the `ComplexTPS64`s in the vector(s) must be allocated and have the same `Descriptor`.

If spin is included, not that the final quaternion concatenation step mul! will create allocations
""" 
function compose_it!(m::$t, m2::$t, m1::$t; dospin::Bool=true, dostochastic::Bool=true, work_Q::Union{Nothing,Quaternion}=prep_work_Q(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))
  checkinplace(m, m2, m1)
  @assert !(m === m1) "Cannot compose_it!(m, m2, m1) with m === m1"

  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  outT = eltype(m.x)

  # Work sanity check:
  if outT != eltype(m1.x)
    x_prom = work_prom[1]
    Q_prom = nothing
    @assert length(x_prom) >= nn "Incorrect length for work_prom[1] = x_prom: Received $(length(x_prom)), should be >=$nn"
  elseif outT != eltype(m2.x)
    x_prom = work_prom[1]
    if !isnothing(m.Q) && dospin
      Q_prom = work_prom[2]
      @assert length(Q_prom) >= 4 "Incorrect length for work_prom[2] = Q_prom: Received $(length(Q_prom)), should be >=4"
    else
      Q_prom = nothing
    end
    @assert length(x_prom) >= nv "Incorrect length for work_prom[1] = x_prom: Received $(length(x_prom)), should be >=$nv"
  else
    x_prom = nothing
    Q_prom = nothing
  end

  # Deal with x0:
  m.x0 .= m1.x0

  # add immutable parameters to outx
  if eltype(m.x) == eltype(m1.x)
    @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)
  else
    @inbounds m.x[nv+1:nn] .= view(m2.x, nv+1:nn)
  end

  # Do the composition, promoting if necessary
  compose!(view(m.x, 1:nv), view(m2.x, 1:nv), view(m1.x, 1:nn), work=x_prom)

  # Spin:
  if !isnothing(m.Q) && dospin
    compose!(work_Q, m2.Q, m1.x, work=Q_prom)
    mul!(m.Q, work_Q, m1.Q) # This will be made faster eventually
  end 

  # Stochastic
  # MAKE THIS FASTER!
  if !isnothing(m.E) && dostochastic
    M2 = jacobian(m2)   
    m.E .= M2*m1.E*transpose(M2) + m2.E
  end

  return 
end

end
end