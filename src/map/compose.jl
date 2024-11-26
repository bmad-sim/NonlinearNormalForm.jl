#=

In-place composition functions compose! for both 
DAMaps and TPSAMaps, as well as the low-level 
compose_it! function which is called by both.

=#

"""
compose!(m::DAMap, m2::DAMap, m1::DAMap; keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, dospin::Bool=true, work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))

In-place `DAMap` composition which calculates `m = m2 ∘ m1`, ignoring the scalar part of `m1`.

Aliasing `m === m2` is allowed, but NOT `m === m1`. `m2 === m1` is fine. The destination map `m` should be 
properly set up (with correct types promoted if necessary), and have `m.x[1:nv]` (and `m.Q` if spin) 
containing allocated TPSs.

If `keep_scalar` is set to `true`, then the scalar part of `m1` is retained. In this case, a temporary 
vector `work_ref` must be used to store the scalar part before calling the low-level `compose_it!`, 
then added back into `m1` after composition. `work_ref` can be optionally provided, or will be created 
internally in this case. Default is `true`.

If `dospin` is `true`, then the quaternion part of the maps will be composed as well. Default is `true`.

See the documentation for `compose_it!` for information on `work_prom`.

### Keyword Arguments
- `keep_scalar` -- Specify whether to keep the scalar part in `m1` or throw it away. If `true`, a temporary vector storing the scalar part must be used. Default is true
- `work_ref` -- If `keep_scalar` is true, the temporary vector can be provided in this keyword argument. If `nothing` is provided, the temporary will be created internally. Default is `nothing`
- `dospin` -- Specify whether or not to include the quaternions in the concatenation. Default is `true`
- `work_prom` -- Temporary vector of allocated `ComplexTPS64`s when there is implicit promotion. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_prom(m, m2, m1)`
"""
function compose!(m::DAMap, m2::DAMap, m1::DAMap; keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, dospin::Bool=true, work_Q::Union{Nothing,Quaternion}=prep_work_Q(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))
  checkinplace(m, m2, m1)
  
  # DAMap setup:
  desc = getdesc(m1)
  nv = numvars(desc)

  if keep_scalar
    if isnothing(work_ref)
      ref = prep_work_ref(m)
    else
      ref=work_ref
    end
    @assert length(ref) >= nv "Incorrect length for work_ref, received $(length(work_ref)) but should be atleast $nv"

    # Take out scalar part and store it
    for i=1:nv
        ref[i] = m1.x[i][0]
        m1.x[i][0] = 0
    end
  else
    for i=1:nv
      m1.x[i][0] = 0
    end
  end

  compose_it!(m, m2, m1, dospin=dospin, work_prom=work_prom)

  # Put back the reference and if m1 === m2, also add to outx
  if keep_scalar
    if m1 === m2
      for i=1:nv
          m1.x[i][0] = ref[i]
          m.x[i][0] += ref[i]
      end
    else
      for i=1:nv
          m1.x[i][0] = ref[i]
      end
    end
  end

  return m
end

"""
    compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; dospin::Bool=true, work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))

In-place `TPSAMap` composition which calculates `m = m2 ∘ m1`, including the scalar part of `m1`.

Aliasing `m === m2` is allowed, but NOT `m === m1`. `m2 === m1` is fine. The destination map `m` should be 
properly set up (with correct types promoted if necessary), and have `m.x[1:nv]` (and `m.Q` if spin) 
containing allocated TPSs.

If `dospin` is `true`, then the quaternion part of the maps will be composed as well. Default is `true`.

See the documentation for `compose_it!` for information on `work_prom`.

### Keyword Arguments
- `dospin` -- Specify whether or not to include the quaternions in the concatenation. Default is `true`
- `work_prom` -- Temporary vector of allocated `ComplexTPS64`s when there is implicit promotion. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_prom(m, m2, m1)`
"""
function compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; dospin::Bool=true, work_Q::Union{Nothing,Quaternion}=prep_work_Q(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS64}}}}=prep_comp_work_prom(m,m2,m1))
  checkinplace(m, m2, m1)
  
  # TPSAMap setup:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  desc = getdesc(m1)
  nv = numvars(desc)
  for i=1:nv
    m1.x[i][0] -= m2.x0[i]
  end

  compose_it!(m, m2, m1, dospin=dospin, work_prom=work_prom)

  # Now fix m1 and if m2 === m1, add to output too:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  if m1 === m2
    for i=1:nv
       m1.x[i][0] += m2.x0[i]
       m.x[i][0] += m2.x0[i]
    end
  else
    for i=1:nv
       m1.x[i][0] += m2.x0[i]
    end
  end

  return m
end


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

"""
    compose(m2::$($t),m1::$($t)) -> $($t)

$($t) composition, which calculates `m2 ∘ m1` $( $t == DAMap ? "ignoring the scalar part of `m1`" : "including the scalar part of `m1`")
"""
function compose(m2::$t,m1::$t)
  m = zero_op(m1, m2)
  compose!(m, m2, m1)
  
  return m
end

end
end


