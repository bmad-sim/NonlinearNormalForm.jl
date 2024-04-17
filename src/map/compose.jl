# compose.jl contains the high-level composition functions for both DAMaps and TPSAMaps

# --- compose! ---

"""
    compose!(m::DAMap, m2::DAMap, m1::DAMap; dospin::Bool=true, keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))

In-place `DAMap` composition which calculates `m = m2 ∘ m1`, ignoring the scalar part of `m1`.

Aliasing `m === m2` is allowed, but NOT `m === m1`. `m2 === m1` is fine. The destination map `m` should be 
properly set up (with correct types promoted if necessary), and have `m.x[1:nv]` (and `m.Q.q` if spin) 
containing allocated TPSs.

If `keep_scalar` is set to `true`, then the scalar part of `m1` is retained. In this case, a temporary 
vector `work_ref` must be used to store the scalar part before calling the low-level `compose_it!`, 
then added back into `m1` after composition. `work_ref` can be optionally provided, or will be created 
internally in this case. Default is `true`.

If `dospin` is `true`, then the quaternion part of the maps will be composed as well. Default is `true`.

See the documentation for `compose_it!` for information on `work_low` and `work_prom`.

### Keyword Arguments
- `keep_scalar` -- Specify whether to keep the scalar part in `m1` or throw it away. If `true`, a temporary vector storing the scalar part must be used. Default is true
- `work_ref` -- If `keep_scalar` is true, the temporary vector can be provided in this keyword argument. If `nothing` is provided, the temporary will be created internally. Default is `nothing`
- `dospin` -- Specify whether or not to include the quaternions in the concatenation. Default is `true`
- `work_low` -- Temporary vector to hold the low-level C pointers. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_low(m)`
- `work_prom` -- Temporary vector of allocated `ComplexTPS`s when there is implicit promotion. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_prom(m, m2, m1)`
"""
function compose!(m::DAMap, m2::DAMap, m1::DAMap; keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))
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
        @inbounds ref[i] = m1.x[i][0]
        @inbounds m1.x[i][0] = 0
    end
  else
    for i=1:nv
      @inbounds m1.x[i][0] = 0
    end
  end

  compose_it!(m, m2, m1, dospin=dospin, work_low=work_low, work_prom=work_prom)

  # Put back the reference and if m1 === m2, also add to outx
  if keep_scalar
    if m1 === m2
      for i=1:nv
          @inbounds m1.x[i][0] = ref[i]
          @inbounds m.x[i][0] += ref[i]
      end
    else
      for i=1:nv
          @inbounds m1.x[i][0] = ref[i]
      end
    end
  end

  return m
end

"""
    compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))

In-place `TPSAMap` composition which calculates `m = m2 ∘ m1`, including the scalar part of `m1`.

Aliasing `m === m2` is allowed, but NOT `m === m1`. `m2 === m1` is fine. The destination map `m` should be 
properly set up (with correct types promoted if necessary), and have `m.x[1:nv]` (and `m.Q.q` if spin) 
containing allocated TPSs.

If `dospin` is `true`, then the quaternion part of the maps will be composed as well. Default is `true`.

See the documentation for `compose_it!` for information on `work_low` and `work_prom`.

### Keyword Arguments
- `dospin` -- Specify whether or not to include the quaternions in the concatenation. Default is `true`
- `work_low` -- Temporary vector to hold the low-level C pointers. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_low(m)`
- `work_prom` -- Temporary vector of allocated `ComplexTPS`s when there is implicit promotion. See the `compose_it!` documentation for more details. Default is output from `prep_comp_work_prom(m, m2, m1)`
"""
function compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))
  # TPSAMap setup:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  desc = getdesc(m1)
  nv = numvars(desc)
  for i=1:nv
    @inbounds m1.x[i][0] -= m2.x0[i]
  end

  compose_it!(m, m2, m1, dospin=dospin, work_low=work_low, work_prom=work_prom)

  # Now fix m1 and if m2 === m1, add to output too:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  if m1 === m2
    for i=1:nv
       @inbounds m1.x[i][0] += m2.x0[i]
       @inbounds m.x[i][0] += m2.x0[i]
    end
  else
    for i=1:nv
       @inbounds m1.x[i][0] += m2.x0[i]
    end
  end

  return m
end

# --- compose ---
for t = (:DAMap, :TPSAMap)
@eval begin
  
"""
    compose(m2::$($t),m1::$($t)) -> $($t)

$($t) composition, which calculates `m2 ∘ m1` $( $t == DAMap ? "ignoring the scalar part of `m1`" : "including the scalar part of `m1`")
"""
function compose(m2::$t,m1::$t)
  @assert !isnothing(m1.Q) && !isnothing(m2.Q) || m1.Q == m2.Q "Cannot compose: one map includes spin, other does not"
  @assert !isnothing(m1.E) && !isnothing(m2.E) || m1.E == m2.E "Cannot compose: one map includes radiation, other does not"

  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)

  outT = promote_type(eltype(m2.x),eltype(m1.x))
  
  # set up outx0
  outx0 = Vector{numtype(outT)}(undef, nv)

  # Set up outx:
  outx = Vector{outT}(undef, nn)
  for i=1:nv  # no need to allocate immutable parameters taken care of inside compose_it!
      @inbounds outx[i] = outT(use=desc)
  end

  # set up quaternion out:
  if !isnothing(m1.Q)
    outq = Vector{outT}(undef, 4)
    for i=1:4
      @inbounds outq[i] = outT(use=desc)
    end
    outQ = Quaternion(outq)
  else
    outQ = nothing
  end

  # set up radiation out
  if isnothing(m1.E)
    outE = nothing
  else
    outE = Matrix{numtype(outT)}(undef, nv, nv)
  end

  m = $t(outx0, outx, outQ, outE)
  compose!(m, m2, m1)
  
  return m
end
end
end