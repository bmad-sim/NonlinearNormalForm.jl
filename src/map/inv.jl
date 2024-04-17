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
      ref = prep_work_ref(m) 
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