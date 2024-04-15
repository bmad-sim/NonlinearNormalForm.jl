# --- compose ---

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

In-place TPSAMap` composition which calculates `m = m2 ∘ m1`, including the scalar part of `m1`.

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

# --- inverse ---
minv!(na::Cint, ma::Vector{Ptr{RTPSA}}, nb::Cint, mc::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_minv!(na, ma, nb, mc))
minv!(na::Cint, ma::Vector{Ptr{CTPSA}}, nb::Cint, mc::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_minv!(na, ma, nb, mc))

"""
    inv(m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}

Inverts the `TaylorMap`.

### Keyword Arguments
- `dospin` -- Specify whether to invert the quaternion as well, default is `true`
- `work_low` -- Temporary vector to hold the low-level C pointers. Default is output from `prep_inv_work_low`
"""
function inv(m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}
  m = zero(m1)
  inv!(m,m1,dospin=dospin,work_low=work_low)
  return m
end

"""
    inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}

In-place inversion of the `TaylorMap` setting `m = inv(m1)`. Aliasing `m === m1` is allowed, however 
in this case a temporary vector must be used to store the scalar part of `m1` prior to inversion so 
that the entrance/exit coordinates of the map can be properly handled.

### Keyword Arguments
- `dospin` -- Specify whether to invert the quaternion as well, default is `true`
- `work_ref` -- If `m === m1`, then a temporary vector must be used to store the scalar part. If not provided and `m === m1`, this temporary will be created internally. Default is `nothing`
- `work_low` -- Temporary vector to hold the low-level C pointers. Default is output from `prep_inv_work_low`
"""
function inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}
  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)

  outx_low = work_low[1]
  m1x_low = work_low[2]
  @assert length(outx_low) >= nn "Cannot inv!: incorrect length for outx_low = work_low[1]. Received $(length(outx_low)), must be >= $nn"
  @assert length(m1x_low) >= nn "Cannot inv!: incorrect length for m1x_low = work_low[2]. Received $(length(m1x_low)), must be >= $nn"
  if !isnothing(m.Q) && dospin
    outQ_low = work_low[3]
    @assert length(outQ_low) >= 4 "Cannot inv!: incorrect length for outQ_low = work_low[3]. Received $(length(outQ_low)), must be >= 4"
    @assert !(outQ_low === outx_low) "Cannot inv!: outQ_low === outx_low !! outx_low must NOT be reused!"
  end
  
  # if aliasing, must use vector to store x0
  if m1 === m
    if isnothing(work_ref)
      ref = Vector{numtype(T)}(undef, nv)
    else
      ref = work_ref
      @assert length(ref) >= nv "Cannot inv!: incorrect length for ref. Received $(length(ref)), must be >= $nv"
    end
    map!(t->t[0], ref, view(m1.x,1:nv))
  elseif !isnothing(work_ref)
    #@warn "work_ref provided to inv!, but !(m1 === m)"
  end

  # add immutable parameters to outx
  @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  map!(t->t.tpsa, m1x_low, m1.x)
  map!(t->t.tpsa, outx_low, m.x)

  # This C function ignores the scalar part so no need to take it out
  minv!(nn, m1x_low, nv, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  if !isnothing(m.Q) && dospin
    inv!(m.Q, m1.Q)
    map!(t->t.tpsa, outQ_low, m.Q.q)
    compose!(Cint(4), outQ_low, nn, outx_low, outQ_low)
  end

  if m1 === m
    for i=1:nv
      @inbounds m.x[i][0] = m.x0[i]
      @inbounds m.x0[i] = ref[i]
    end
  else
    for i=1:nv
       @inbounds m.x0[i] = m1.x[i][0]
       @inbounds m.x[i][0] = m1.x0[i]
    end
  end
  
  return 
end


for t = (:DAMap, :TPSAMap)
@eval begin  

# --- norm ---
function norm(m::$t)
  n = norm(m.x0) + norm(norm.(m.x))
  if !isnothing(m.Q)
    n += norm(m.Q.q)
  end

  if !isnothing(m.E)
    n += norm(m.E)
  end
  return n
end

# --- complex ---
function complex(m::$t{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  x0 = map(t->complex(t), m.x0)
  
  x = Vector{ComplexTPS}(undef, nn)
  for i=1:nv
    @inbounds x[i] = ComplexTPS(m.x[i],use=desc)
  end

  # use same parameters if complex already
  if T == ComplexTPS
    @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)
  else
    @inbounds x[nv+1:nn] = complexparams(getdesc(first(x)))
  end

  if !isnothing(m.Q)
    q = Vector{ComplexTPS}(undef, 4)
    for i=1:4
      @inbounds q[i] = ComplexTPS(m.Q.q[i],use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = map(t->complex(t), m.E)
  else
    E = nothing
  end
  return $t(x0, x, Q, E)
end

end
end

# --- clear ---
function clear!(m::TaylorMap)
  nv = numvars(m)
  m.x0 .= 0
  for i=1:nv
    @inbounds clear!(m.x[i])
  end
  if !isnothing(m.Q)
    for i=1:4
      @inbounds clear!(m.Q.q[i])
    end
  end
  if !isnothing(m.E)
    m.E .= 0
  end
  return
end

# --- cut ---
function cutord(m1::TaylorMap{S,T,U,V}, order::Integer; dospin::Bool=true) where {S,T,U,V}
  m = zero(m1)
  cutord!(m, m1, order, dospin=dospin)
  return m
end

function cutord!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}, order::Integer; dospin::Bool=true) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  m.x0 .= m1.x0
  for i=1:nv
    @inbounds cutord!(m.x[i], m1.x[i], order)
    #GTPSA.cutord!(m1.x[i].tpsa, m.x[i].tpsa, convert(Cint, ord))
  end
  # add immutable parameters to outx
  @inbounds m.x[nv+1:nn] .= view(m1.x, nv+1:nn)

  if !isnothing(m1.Q) && dospin
    for i=1:4
      @inbounds cutord!(m.Q.q[i], m1.Q.q[i], order)
      #@inbounds GTPSA.cutord!(m1.Q.q[i].tpsa, m.Q.q[i].tpsa, convert(Cint, ord))
    end
  end
  if !isnothing(m1.E)
    m.E .= m.E
  end
  return
end