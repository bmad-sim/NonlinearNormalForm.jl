# --- copy! ---
function copy!(F::VectorField, F1::Union{VectorField,TaylorMap})
  checkinplace(F, F1)

  desc = getdesc(F)
  nv = numvars(desc)

  for i=1:nv
   @inbounds copy!(F.x[i], F1.x[i])
  end

  if !isnothing(F.Q)
    copy!(F.Q.q0, F1.Q.q0)
    copy!(F.Q.q1, F1.Q.q1)
    copy!(F.Q.q2, F1.Q.q2)
    copy!(F.Q.q3, F1.Q.q3)
  end

  return F
end

# --- clear! ---

function clear!(F::VectorField)
  nv = numvars(F)
  for i=1:nv
    @inbounds clear!(F.x[i])
  end
  if !isnothing(F.Q)
    clear!(F.Q.q0)
    clear!(F.Q.q1)
    clear!(F.Q.q2)
    clear!(F.Q.q3)
  end
  return
end

# --- complex ---
#=
function complex(F::VectorField)
  desc = getdesc(F)
  nn = numnn(desc)
  nv = numvars(desc)
  
  x = Vector{ComplexTPS64}(undef, nv)
  for i=1:nv
    @inbounds x[i] = ComplexTPS64(F.x[i],use=desc)
  end

  if !isnothing(m.Q)
    q = Vector{ComplexTPS64}(undef, 4)
    for i=1:4
      @inbounds q[i] = ComplexTPS64(F.Q.q[i],use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  return VectorField(x,Q)
end
=#
# --- complex type ---
function complex(::Type{VectorField{T,U}}) where {T,U}
  return VectorField{ComplexTPS64, U == Nothing ? Nothing : Quaternion{ComplexTPS64}}
end

# --- Lie bracket including spin ---
# GTPSA only provides routines for orbital part:
liebra!(na, m1::Vector{TPS{Float64}},    m2::Vector{TPS{Float64}},    m3::Vector{TPS{Float64}})    = GTPSA.mad_tpsa_liebra!(Cint(na), m1, m2, m3)
liebra!(na, m1::Vector{TPS{ComplexF64}}, m2::Vector{TPS{ComplexF64}}, m3::Vector{TPS{ComplexF64}}) = GTPSA.mad_ctpsa_liebra!(Cint(na), m1, m2, m3)

"""
    lb(F::VectorField, H::VectorField)

Computes the Lie bracket of the vector fields `F` and `H`. Explicitly, and including 
spin (with the lower case letter for the quaternion of the vector field), this is:

`(G,g) = ⟨(F,f) , (H,h)⟩ = (F⋅∇H-H⋅∇F , [h,f]+F⋅∇h-G⋅∇f)`

where `[h,f] = h*f - f*h` is just a quaternion commutator. See Equation 44.52 in the Bmad manual
"""
function lb(F::VectorField, H::VectorField)
  checkop(F,H)

  G = zero(F)
  lb!(G,F,H)
  return G
end


"""
    lb!(G::VectorField, F::VectorField, H::VectorField; work_Q::Union{Quaternion,Nothing}=nothing)

Computes the Lie bracket of the vector fields `F` and `H`, and stores the result in `G`. 
See `lb` for more details.

### Keyword Arguments
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function lb!(G::VectorField, F::VectorField, H::VectorField; work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))
  checkinplace(G, F, H)

  T = eltype(G.x)

  nv = numvars(F)
  @assert !(G === F) && !(G === H) "Aliasing any source arguments with the destination in lb! is not allowed"

  # Orbital part (Eq. 44.51 in Bmad manual, handled by GTPSA):
  liebra!(nv, F.x, H.x, G.x)

  # Spin (Eq. 44.51 in Bmad manual, NOT handled by GTPSA):
  if !isnothing(F.Q)
    # first [h,f] (this part involves quaternion multiplication so will be slow)
    mul!(G.Q, H.Q, F.Q)
    mul!(work_Q, F.Q, H.Q)

    sub!(G.Q.q0, G.Q.q0, work_Q.q0)
    sub!(G.Q.q1, G.Q.q1, work_Q.q1)
    sub!(G.Q.q2, G.Q.q2, work_Q.q2)
    sub!(G.Q.q3, G.Q.q3, work_Q.q3)

    # then +F⋅∇h 
    tmp = work_Q.q0
    fgrad!(tmp, F.x, H.Q.q0)
    add!(G.Q.q0, G.Q.q0, tmp)

    fgrad!(tmp, F.x, H.Q.q1)
    add!(G.Q.q1, G.Q.q1, tmp)

    fgrad!(tmp, F.x, H.Q.q2)
    add!(G.Q.q2, G.Q.q2, tmp)

    fgrad!(tmp, F.x, H.Q.q3)
    add!(G.Q.q3, G.Q.q3, tmp)
    
    # finally -G⋅∇f
    fgrad!(tmp, G.x, F.Q.q0)
    sub!(G.Q.q0, G.Q.q0, tmp)

    fgrad!(tmp, G.x, F.Q.q1)
    sub!(G.Q.q1, G.Q.q1, tmp)

    fgrad!(tmp, G.x, F.Q.q2)
    sub!(G.Q.q2, G.Q.q2, tmp)

    fgrad!(tmp, G.x, F.Q.q3)
    sub!(G.Q.q3, G.Q.q3, tmp)
  end
  return
end


# --- Lie operator (VectorField) acting on DAMap ---
# GTPSA provides fgrad with is used, but nothing for spin
"""
    *(F::VectorField, m1::DAMap)

Calculates a Lie operator `VectorField` acting on a `DAMap`. Explicity, if spin is 
included that is `F * m = (F.x, F.Q) * (m.x, m.Q) = (F.x ⋅ ∇ m.x , F.x ⋅ ∇ m.Q + m.Q*F.Q)`
"""
function *(F::VectorField, m1::Union{DAMap,UniformScaling})
  if m1 isa UniformScaling 
    m1 = one(DAMap{Vector{eltype(eltype(F.x))},typeof(F.x),typeof(F.Q),Nothing,Nothing},use=getdesc(F))
  end
  checkop(F, m1)
  m = zero(m1)
  mul!(m, F, m1)
  return m
end

"""
    mul!(m::DAMap, F::VectorField, m1::DAMap; work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))

Computes the Lie operator `F` acting on a `DAMap` `m1`, and stores the result in `m`.
Explicity, that is `F * m = (F.x, F.Q) * (m.x, m.Q) = (F.x ⋅ ∇ m.x , F.x ⋅ ∇ m.Q + m.Q*F.Q)`

### Keyword Arguments\
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function mul!(m::DAMap, F::VectorField, m1::DAMap; work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))
  checkinplace(m, F, m1)

  T = eltype(m.x)
  nv = numvars(F)

  @assert !(m === m1) "Aliasing m === m1 is not allowed for mul!"

  m.x0 .= m1.x0
  # Orbital part F⋅∇m1 :
  for i=1:nv
    @inbounds fgrad!(m.x[i], F.x, m1.x[i])
  end

  # Spin F⋅∇q + qf (quaternion part acts in reverse order):
  if !isnothing(F.Q)
    mul!(work_Q, m1.Q, F.Q)
    fgrad!(m.Q.q0, F.x, m1.Q.q0)
    fgrad!(m.Q.q1, F.x, m1.Q.q1)
    fgrad!(m.Q.q2, F.x, m1.Q.q2)
    fgrad!(m.Q.q3, F.x, m1.Q.q3)

    add!(m.Q.q0, m.Q.q0, work_Q.q0)
    add!(m.Q.q1, m.Q.q1, work_Q.q1)
    add!(m.Q.q2, m.Q.q2, work_Q.q2)
    add!(m.Q.q3, m.Q.q3, work_Q.q3)
  end

  if !isnothing(m.E)
    m.E .= m1.E
  end

  return
end


# --- exp(F)*m ---
# While GTPSA provides this for only the orbital part, it separately converges 
# each variable and does not include spin. Here I include spin and do the entire map at once.
function exp(F::VectorField, m1::Union{UniformScaling,DAMap}=I)
  if m1 isa UniformScaling 
    m1 = one(DAMap{Vector{eltype(eltype(F.x))},typeof(F.x),typeof(F.Q),Nothing,Nothing},use=getdesc(F))
  end
  checkop(F,m1)

  m = zero(m1)
  exp!(m, F, m1)
  return m
end


"""
    exp!(m::DAMap, F::VectorField, m1::DAMap; work_maps::Tuple{Vararg{DAMap}}=(zero(m1),zero(m1)), work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))

Computes `exp(F)*m1`, and stores the result in `m`. Explicity, this is
`exp(F)*m1 = m1 + F*m1 + 1/2*F*(F*m1) + 1/6*F*(F*(F*m1)) + ...`, where `*` is
the operation of a `VectorField` on the map. See the documentation for `mul!` for 
more details of this operation. If the linear part of `F` is zero, then the 
number of iterations is equal to the number of higher orders left.

### Keyword Arguments
- `work_maps` -- Tuple of 2 `DAMap`s of the same type as `m1`
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function exp!(m::DAMap, F::VectorField, m1::DAMap; work_maps::Tuple{Vararg{DAMap}}=(zero(m1),zero(m1)), work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))
  checkinplace(m, F, m1, work_maps...)

  T = eltype(m.x)
  nv = numvars(F)
  #GTPSA.mad_ctpsa_exppb!(nv, map(t->t.tpsa, F.x), map(t->t.tpsa, view(m1.x, 1:nv)), map(t->t.tpsa, view(m.x,1:nv)))
  #return
  tmp = work_maps[1]
  tmp2 = work_maps[2]

  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*nv
  nrm_ =Inf
  conv = false
  slow = false

  copy!(tmp, m1)
  copy!(m, m1)

  for j=1:nmax
    #println("at iteration $j")
    if j == 25
      slow=true
    end

    div!(tmp, tmp, j)
    mul!(tmp2, F, tmp, work_Q=work_Q)
    add!(m, m, tmp2)
    copy!(tmp, tmp2)
    nrm = norm(norm(tmp2))

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
# 2 work_vfs

"""
    log!(F::VectorField, m1::DAMap; work::Tuple{DAMap,DAMap,DAMap,VectorField,VectorField}=prep_log_work(m1), work_Q::Union{Quaternion,Nothing}=prep_work_Q(F)) 

Computes the log of the map `m1` - that is, calculates the `VectorField` `F` that 
would represent the map `m1` as a Lie exponent `exp(F)` - and stores the result in `F`.
The map `m1` should be close to the identity for this to converge quickly.

### Keyword Arguments
- `work` -- Tuple of 3 `DAMap`s of the same type as `m1` followed by 2 `VectorField`S
- `work_Q`   -- `Quaternion{T}` if spin is included in the vector field, else `nothing`
"""
function log!(F::VectorField, m1::DAMap; work::Tuple{DAMap,DAMap,DAMap,VectorField,VectorField}=prep_log_work(m1), work_Q::Union{Quaternion,Nothing}=prep_work_Q(F))
  checkinplace(F, m1, work...)
  
  nv = numvars(m1)
  nmax = 100
  nrm0 = norm(m1)
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*nv
  nrm0 = norm(norm(m1))
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
    exp!(work_maps[3], F, m1, work_maps=work_maps, work_Q=work_Q)
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
    mul!(work_maps[2], work_vfs[1], work_maps[3], work_Q=work_Q)

    # Finally get our approximation 
    add!(work_vfs[1], work_vfs[1], work_maps[2])

    # Now we can see that  M = exp(F)exp(G) approximately
    # CBH formula would be exact combination for exp(F)exp(G)=exp(H)
    # but if the leftover stuff is still big then we can just approximately
    # with linear term (adding F+G trivially)

    if norm(norm(work_vfs[1])) < epsone # small! use CBH!
      lb!(work_vfs[2], F, work_vfs[1], work_Q=work_Q)
      mul!(work_vfs[2], 0.5, work_vfs[2])
      add!(work_vfs[1], work_vfs[1], work_vfs[2])
      add!(F,F,work_vfs[1])
    else # big! just linear!
      add!(F,F,work_vfs[1])
    end

    nrm = norm(norm(work_vfs[1])/nrm0)

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
    log(m1::DAMap)

Computes the log of the map `m1` - that is, calculates the `VectorField` `F` that 
would represent the map `m1` as a Lie exponent `exp(F)`. The map `m1` should be close 
to the identity for this to converge quickly.
"""
function log(m1::DAMap)
  F = zero(VectorField{typeof(m1.x),typeof(m1.Q)},use=m1)
  log!(F,m1)
  return F
end



