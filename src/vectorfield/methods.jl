# --- copy! ---
function copy!(F::VectorField{T,U}, F1::VectorField{T,U}) where {T,U}
  desc = getdesc(F)
  nv = numvars(desc)

  for i=1:nv
   @inbounds copy!(F.x[i], F1.x[i])
  end

  if !isnothing(F.Q)
    for i=1:4
      @inbounds copy!(F.Q.q[i], F1.Q.q[i])
    end
  end

  return F
end

# --- complex ---
function complex(F::VectorField{T,U}) where {T,U}
  desc = getdesc(F)
  nn = numnn(desc)
  nv = numvars(desc)
  
  x = Vector{ComplexTPS}(undef, nv)
  for i=1:nv
    @inbounds x[i] = ComplexTPS(F.x[i],use=desc)
  end

  if !isnothing(m.Q)
    q = Vector{ComplexTPS}(undef, 4)
    for i=1:4
      @inbounds q[i] = ComplexTPS(F.Q.q[i],use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  return VectorField(x,Q)
end

# --- complex type ---
function complex(::Type{VectorField{T,U}}) where {T,U}
  return VectorField{ComplexTPS, U == Nothing ? Nothing : Quaternion{ComplexTPS}}
end

# --- Lie bracket including spin ---
# GTPSA only provides routines for orbital part:
liebra!(na::Cint, m1::Vector{Ptr{RTPSA}}, m2::Vector{Ptr{RTPSA}}, m3::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_liebra!(na, m1, m2, m3))
liebra!(na::Cint, m1::Vector{Ptr{CTPSA}}, m2::Vector{Ptr{CTPSA}}, m3::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_liebra!(na, m1, m2, m3))

"""
    lb(F::VectorField{T,U}, H::VectorField{T,U}) where {T,U}

Computes the Lie bracket of the vector fields `F` and `H`. Explicitly, and including 
spin (with the lower case letter for the quaternion of the vector field), this is:

`(G,g) = ⟨(F,f) , (H,h)⟩ = (F⋅∇H-H⋅∇F , [h,f]+F⋅∇h-G⋅∇f)`

where `[h,f] = h*f - f*h` is just a quaternion commutator. See Equation 44.52 in the Bmad manual
"""
function lb(F::VectorField{T,U}, H::VectorField{T,U}) where {T,U}
  G = zero(F)
  lb!(G,F,H)
  return G
end


"""
    lb!(G::VectorField{T,U}, F::VectorField{T,U}, H::VectorField{T,U}; work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_vf_work_low(F), work_Q::Union{Nothing,U}=nothing) where {T,U}

Computes the Lie bracket of the vector fields `F` and `H`, and stores the result in `G`. 
See `lb` for more details.

### Keyword Arguments
- `work_low` -- Tuple of 3 vectors of type `lowtype(T)` of length `>=nv`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function lb!(G::VectorField{T,U}, F::VectorField{T,U}, H::VectorField{T,U}; work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_lb_work_low(F), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {T,U}
  nv = numvars(F)
  Fx_low = work_low[1]
  Hx_low = work_low[2]
  Gx_low = work_low[3]
  @assert length(Fx_low) >= nv "Incorrect length for work_low[1] = Fx_low; received $(length(Fx_low)), should be >=$nv"
  @assert length(Hx_low) >= nv "Incorrect length for work_low[2] = Hx_low; received $(length(Hx_low)), should be >=$nv"
  @assert length(Gx_low) >= nv "Incorrect length for work_low[3] = Gx_low; received $(length(Gx_low)), should be >=$nv"
  @assert eltype(Fx_low) == lowtype(T) "Incorrect eltype of work_low[1] = Fx_low. Received $(eltype(Fx_low)), should be $(lowtype(T))"
  @assert eltype(Hx_low) == lowtype(T) "Incorrect eltype of work_low[2] = Hx_low. Received $(eltype(Hx_low)), should be $(lowtype(T))"
  @assert eltype(Gx_low) == lowtype(T) "Incorrect eltype of work_low[3] = Gx_low. Received $(eltype(Gx_low)), should be $(lowtype(T))"
  @assert !(G === F) && !(G === H) "Aliasing any source arguments with the destination in lb! is not allowed"

  # Orbital part (Eq. 44.51 in Bmad manual, handled by GTPSA):
  map!(t->t.tpsa, Fx_low, F.x)
  map!(t->t.tpsa, Hx_low, H.x)
  map!(t->t.tpsa, Gx_low, G.x)
  liebra!(nv, Fx_low, Hx_low, Gx_low)

  # Spin (Eq. 44.51 in Bmad manual, NOT handled by GTPSA):
  if !isnothing(F.Q)
    # first [h,f] (this part involves quaternion multiplication so will be slow)
    mul!(G.Q, H.Q, F.Q)
    mul!(work_Q, F.Q, H.Q)
    for i=1:4
      @inbounds sub!(G.Q.q[i], G.Q.q[i], work_Q[i])
    end

    # then +F⋅∇h 
    tmp = work_Q.Q.q[1]
    for i=1:4
      @inbounds fgrad!(tmp, F.x, H.Q.q[i], work_low=Fx_low)
      @inbounds add!(G.Q.q[i], G.Q.q[i], tmp)
    end
    
    # finally -G⋅∇f
    for i=1:4
      @inbounds fgrad!(tmp, G.x, F.Q.q[i], work_low=Fx_low)
      @inbounds sub!(G.Q.q[i], G.Q.q[i], tmp)
    end
  end
  return
end


# --- Lie operator (VectorField) acting on DAMap ---
# GTPSA provides fgrad with is used, but nothing for spin
"""
    *(F::VectorField{T,U}, m1::DAMap{S,T,U,V}) where {S,T,U,V}

Calculates a Lie operator `VectorField` acting on a `DAMap`. Explicity, if spin is 
included that is `F * m = (F.x, F.Q) * (m.x, m.Q) = (F.x ⋅ ∇ m.x , F.x ⋅ ∇ m.Q + m.Q*F.Q)`
"""
function *(F::VectorField{T,U}, m1::Union{DAMap{S,T,U,V},UniformScaling}) where {S,T,U,V}
  if m1 isa UniformScaling 
    m1 = one(DAMap{numtype(T),T,U,Nothing},use=getdesc(F))
  end
  m = zero(m1)
  mul!(m, F, m1)
  return m
end

"""
    mul!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}

Computes the Lie operator `F` acting on a `DAMap` `m1`, and stores the result in `m`.
Explicity, that is `F * m = (F.x, F.Q) * (m.x, m.Q) = (F.x ⋅ ∇ m.x , F.x ⋅ ∇ m.Q + m.Q*F.Q)`

### Keyword Arguments
- `work_low` -- Vector with eltype `lowtype(T)` and length `>=nv`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function mul!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(F)

  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low. Received $(eltype(work_low)), should be $(lowtype(T))"
  @assert !(m === m1) "Aliasing m === m1 is not allowed for mul!"


  # Orbital part F⋅∇m1 :
  for i=1:nv
    @inbounds fgrad!(m.x[i], F.x, m1.x[i], work_low=work_low)
  end

  # Spin F⋅∇q + qf (quaternion part acts in reverse order):
  if !isnothing(F.Q)
    mul!(work_Q, m1.Q, F.Q)

    for i=1:4
      @inbounds fgrad!(m.Q.q[i], F.x, m1.Q.q[i], work_low=work_low)
    end
    
    for i=1:4
      @inbounds add!(m.Q.q[i], m.Q.q[i], work_Q.q[i])
    end
  end
  return
end


# --- exp(F)*m ---
# While GTPSA provides this for only the orbital part, it separately converges 
# each variable and does not include spin. Here I include spin and do the entire map at once.
function exp(F::VectorField{T,U}, m1::Union{UniformScaling,DAMap{S,T,U,V}}=I) where {S,T,U,V}
  if m1 isa UniformScaling 
    m1 = one(DAMap{numtype(T),T,U,Nothing},use=getdesc(F))
  end
  m = zero(m1)
  exp!(m, F, m1)
  return m
end


"""
    exp!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_maps::Tuple{Vararg{DAMap{S,T,U,V}}}=(zero(m1),zero(m1)), work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}

Computes `exp(F)*m1`, and stores the result in `m`. Explicity, this is
`exp(F)*m1 = m1 + F*m1 + 1/2*F*(F*m1) + 1/6*F*(F*(F*m1)) + ...`, where `*` is
the operation of a `VectorField` on the map. See the documentation for `mul!` for 
more details of this operation.

### Keyword Arguments
- `work_maps` -- Tuple of 2 `DAMap`s of the same type as `m1`
- `work_low` -- Vector with eltype `lowtype(T)` and length `>=nv`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function exp!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_maps::Tuple{Vararg{DAMap{S,T,U,V}}}=(zero(m1),zero(m1)), work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(F)

  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low. Received $(eltype(work_low)), should be $(lowtype(T))"
  tmp = work_maps[1]
  tmp2 = work_maps[2]

  # convergence parameters taken exactly from GTPSA
  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(numtype(T))*nv
  nrm_ =Inf
  conv = false

  copy!(tmp, m1)
  copy!(m, m1)

  for j=1:nmax
    div!(tmp, tmp, j)
    mul!(tmp2, F, tmp, work_low=work_low, work_Q=work_Q)
    add!(m, m, tmp2)
    copy!(tmp, tmp2)
    nrm = 0.
    for i=1:nv
      @inbounds nrm += norm(tmp2.x[i])
    end

    # Check convergence
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      return
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  @warn "exp! convergence not reached for $nmax iterations"
end


# --- log ---

# See Bmad manual Ch. 47 for log definition using Lie operators (no matrices!!)
function log(m1::DAMap{S,T,U,V}) where {S,T,U,V}
  nv = numvars(m1)
  F = zero(VectorField{T,U},use=m1)
  tmp = m1 - I
  F.x = tmp.x[1:nv]
  F.Q  = tmp.Q

  for i=1:1000
    mt = exp(F,m1)
    
  end

  # Start with first approx for vector field F = (t,u) = (M-I, q-1)



  eps2 = -1/2*F*tmp  # eqn 47.6
  T3 = F + eps2

end


"""

"""
function log!(F::VectorField{T,U}, m1::DAMap{S,T,U,V}) where {S,T,U,V}


end



















logpb!(na::Cint, ma::Vector{Ptr{RTPSA}}, mb::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_logpb!(na, ma, mb, mc))
logpb!(na::Cint, ma::Vector{Ptr{CTPSA}}, mb::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_logpb!(na, ma, mb, mc))




"""


All arguments must have correct types. No promotion here yet. I am tired. Maybe soon

### Keyword arguments
- `'work_low` -- Tuple of 3 vectors of length nv with type equal to the lowtype of `F`.

"""
function log!(F::VectorField{T,U}, m2::DAMap{S,T,U,V}, m1::DAMap{S,T,U,V}; work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_log_work_low(F)) where {S,T,U,V}
  nv = numvars(F)

  Fx_low = work_low[1]
  m2x_low = work_low[2]
  m1x_low = work_low[3]

  @assert length(Fx_low) >= nv "Incorrect length for work_low[1] = Fx_low; received $(length(Fx_low)), should be >=$nv"
  @assert length(m2x_low) >= nv "Incorrect length for work_low[1] = m2x_low; received $(length(m2x_low)), should be >=$nv"
  @assert length(m1x_low) >= nv "Incorrect length for work_low[1] = m1x_low; received $(length(m1x_low)), should be >=$nv"
  @assert eltype(outx_low) == lowtype(outT) "Incorrect eltype of work_low[1] = outx_low. Received $(eltype(Fx_low)), should be $(lowtype(T))"
  @assert eltype(m2x_low) == lowtype(outT) "Incorrect eltype of work_low[2] = m2x_low. Received $(eltype(m2x_low)), should be $(lowtype(T))"
  @assert eltype(m1x_low) == lowtype(outT) "Incorrect eltype of work_low[3] = m1x_low. Received $(eltype(m1x_low)), should be $(lowtype(T))"



end


"""
    logpb(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))

Given a final map `mf` and initial map `mi`, this function calculates the vector field `F`
such that `mf=exp(F⋅∇)mi`. If `mi` is not provided, it is assumed to be the identity.

```julia-repl
julia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];

julia> time = 0.01; k = 2; m = 0.01;

julia> h = p^2/(2m) + 1/2*k*x^2;

julia> hf = getvectorfield(h);

julia> map = exppb(-time*hf);

julia> logpb(map) == -time*hf
true
```
"""
function log(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(mf)))
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(mf[1].tpsa).d))
  if length(mf) != desc.nv || length(mi) != desc.nv
    error("Vector length != number of variables in the GTPSA")
  end
  ma1, mb1 = promote(mf, mi)
  m1 = map(t->t.tpsa, ma1) 
  m2 = map(t->t.tpsa, mb1) 
  mc = zero.(ma1)
  m3 = map(t->t.tpsa, mc)
  GC.@preserve ma1 mb1 logpb!(Cint(length(mf)), m1, m2, m3)  
  return mc
end
