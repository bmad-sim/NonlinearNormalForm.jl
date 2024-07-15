# --- inverse ---
minv!(na, ma::Vector{TPS{Float64}},    nb, mc::Vector{TPS{Float64}})    = GTPSA.mad_tpsa_minv!(Cint(na), ma, Cint(nb), mc)
minv!(na, ma::Vector{TPS{ComplexF64}}, nb, mc::Vector{TPS{ComplexF64}}) = GTPSA.mad_ctpsa_minv!(Cint(na), ma, Cint(nb), mc)

"""
    inv(m1::TaylorMap; dospin::Bool=true)

Inverts the `TaylorMap`.

### Keyword Arguments
- `dospin` -- Specify whether to invert the quaternion as well, default is `true`
"""
function inv(m1::TaylorMap; dospin::Bool=true)
  m = zero(m1)
  inv!(m,m1,dospin=dospin)
  return m
end

#=
function inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low=nothing) where {S,T,U,V}
  checkidpt(m,m1)
  lin = cutord(m1, 2)
  np = numparams(m1)
  nn = numnn(m1)
  nv = numvars(m1)
  M = zeros(eltype(m1), nn, nn)
  M[1:nv,1:nn] = jacobian(lin,include_params=true)
  for i=1:np
    M[nv+i,nv+i] = 1
  end
  lin_inv = DAMap(inv(M)[1:nv,1:nn],use=m1,idpt=m.idpt)
  n1 = lin_inv*m1
  F = log(n1)
  n1 = exp(-F,one(lin_inv))
  copy!(m, n1*lin_inv)
  return
end
=#

"""
    inv!(m::TaylorMap, m1::TaylorMap; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing)

In-place inversion of the `TaylorMap` setting `m = inv(m1)`. Aliasing `m === m1` is allowed, however 
in this case a temporary vector must be used to store the scalar part of `m1` prior to inversion so 
that the entrance/exit coordinates of the map can be properly handled.

### Keyword Arguments
- `dospin` -- Specify whether to invert the quaternion as well, default is `true`
- `work_ref` -- If `m === m1`, then a temporary vector must be used to store the scalar part. If not provided and `m === m1`, this temporary will be created internally. Default is `nothing`
"""
function inv!(m::TaylorMap, m1::TaylorMap; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing)
  checkinplace(m, m1)

  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)
  
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

  # This C function ignores the scalar part so no need to take it out
  minv!(nn, m1.x, nv, m.x)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  if !isnothing(m.Q) && dospin
    inv!(m.Q, m1.Q)
    compose!(m.Q.q, m.Q.q, m.x)
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

