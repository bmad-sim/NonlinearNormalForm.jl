# --- copy! ---
function copy!(F::VectorField{T,U}, F1::Union{VectorField{T,U},TaylorMap{S,T,U,V}}) where {S,T,U,V}
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
      @inbounds sub!(G.Q.q[i], G.Q.q[i], work_Q.q[i])
    end

    # then +F⋅∇h 
    tmp = work_Q.q[1]
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

  m.x0 .= m1.x0
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

  if !isnothing(m.E)
    m.E .= m1.E
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
more details of this operation. If the linear part of `F` is zero, then the 
number of iterations is equal to the number of higher orders left.

### Keyword Arguments
- `work_maps` -- Tuple of 2 `DAMap`s of the same type as `m1`
- `work_low` -- Vector with eltype `lowtype(T)` and length `>=nv`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function exp!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_maps::Tuple{Vararg{DAMap{S,T,U,V}}}=(zero(m1),zero(m1)), work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(F)
  #GTPSA.mad_tpsa_exppb!(nv, map(t->t.tpsa, F.x), map(t->t.tpsa, view(m1.x, 1:nv)), map(t->t.tpsa, view(m.x,1:nv)))
  
  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low. Received $(eltype(work_low)), should be $(lowtype(T))"
  tmp = work_maps[1]
  tmp2 = work_maps[2]

  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(numtype(T))*nv
  nrm_ =Inf
  conv = false
  slow = false

  copy!(tmp, m1)
  copy!(m, m1)

  for j=1:nmax
    if j == 25
      slow=true
    end

    div!(tmp, tmp, j)
    mul!(tmp2, F, tmp, work_low=work_low, work_Q=work_Q)
    add!(m, m, tmp2)
    copy!(tmp, tmp2)
    nrm = norm(tmp2)

    # Check convergence
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      if slow
        @warn "exp! slow convergence: required n = $(j) iterations"
      end
      return
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  @warn "exp! convergence not reached for $nmax iterations"
  return 
end


# --- log ---

# See Bmad manual Ch. 47 for calculation of log using Lie operators
# 3 work_maps
# 1 work_Q
# 1 work_low length >= nv
# 2 work_vfs

"""
    log!(F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work::Tuple{DAMap{S,T,U,V},DAMap{S,T,U,V},DAMap{S,T,U,V},VectorField{T,U},VectorField{T,U}}=prep_log_work(m1), work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_lb_work_low(F), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}

Computes the log of the map `m1` - that is, calculates the `VectorField` `F` that 
would represent the map `m1` as a Lie exponent `exp(F)` - and stores the result in `F`.
The map `m1` should be close to the identity for this to converge quickly.

### Keyword Arguments
- `work` -- Tuple of 3 `DAMap`s of the same type as `m1` followed by 2 `VectorField`S
- `work_low` -- Vector with eltype `lowtype(T)` and length `>=nv`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function log!(F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work::Tuple{DAMap{S,T,U,V},DAMap{S,T,U,V},DAMap{S,T,U,V},VectorField{T,U},VectorField{T,U}}=prep_log_work(m1), work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_lb_work_low(F), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(m1)
  nmax = 100
  nrm0 = norm(m1)
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(numtype(T))*nv
  nrm0 = norm(m1)
  epsone = nrm0/1000
  conv = false 
  slow = false
  nrm_ = Inf

  # exp requires 2 work_maps, 1 work_low (>= nv length), 1 work_Q
  # mul (called by exp) requires work_low and work_Q
  # lb requires 3 work_low (>=nv length), and work_Q

  work_maps = (work[1], work[2], work[3])
  work_vfs = (work[4], work[5])

  # Choose the initial guess for the VectorField to be (M,q) - (I,1)
  sub!(F, m1, I)

  # Now we will iterate:
  for j=1:nmax
    if j == 25
      slow=true
    end
    # First, rotate back our guess:
    mul!(F, -1, F)
    exp!(work_maps[3], F, m1, work_maps=work_maps, work_low=work_low[1], work_Q=work_Q)
    mul!(F, -1, F)

    # Now we have (I,1) + (t,u) where (t,u) is the leftover garbage we want to be 0
    # solve for the vector field that gives us our garbage map
    # this will be exp(T)(I,1) = (I,1) + (t,u)
    # to second order in (t,u) is 
    # G = (t,u) + epsilon_2 where 
    # epsilon_2 = (t dot grad t, t dot grad u + u^2)
    # get the garbage as a VectorField and as a map
    sub!(work_maps[3], work_maps[3], I)
    copy!(work_vfs[1], work_maps[3])

    # Let the vector field act on the garbage map
    mul!(work_maps[2], work_vfs[1], work_maps[3], work_low=work_low[1], work_Q=work_Q)

    # Finally get our approximation 
    add!(work_vfs[1], work_vfs[1], work_maps[2])

    # Now we can see that  M = exp(F)exp(G) approximately
    # CBH formula would be exact combination for exp(F)exp(G)=exp(H)
    # but if the leftover stuff is still big then we can just approximately
    # with linear term (adding F+G trivially)

    if norm(work_vfs[1]) < epsone # small! use CBH!
      lb!(work_vfs[2], F, work_vfs[1], work_low=work_low, work_Q=work_Q)
      mul!(work_vfs[2], 0.5, work_vfs[2])
      add!(work_vfs[1], work_vfs[1], work_vfs[2])
      add!(F,F,work_vfs[1])
    else # big! just linear!
      add!(F,F,work_vfs[1])
    end

    nrm = norm(work_vfs[1])/nrm0

    if nrm <= nrm_min2 || conv && nrm >= nrm_
      if slow
        @warn "log! slow convergence: required n = $(j) iterations"
      end
      return
    end

    if nrm <= nrm_min1
      conv = true
    end

    nrm_ = nrm
  end
  @warn "log! convergence not reached for $(nmax) iterations"
  return
end


"""
    log(m1::DAMap{S,T,U,V}) where {S,T,U,V}

Computes the log of the map `m1` - that is, calculates the `VectorField` `F` that 
would represent the map `m1` as a Lie exponent `exp(F)`. The map `m1` should be close 
to the identity for this to converge quickly.
"""
function log(m1::DAMap{S,T,U,V}) where {S,T,U,V}
  F = zero(VectorField{T,U},use=m1)
  log!(F,m1)
  return F
end



