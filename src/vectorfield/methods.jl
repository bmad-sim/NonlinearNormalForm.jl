# --- copy! ---
function copy!(F::VectorField{T,U}, F1::VectorField{T,U}) where {T,U}
  desc = getdesc(F)
  nv = numvars(desc)

  for i=1:nv
    copy!(F.x[i], F1.x[i])
  end

  if !isnothing(F.Q)
    for i=1:4
      copy!(F.Q.q[i], F1.Q.q[i])
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

  # Orbital part (Eq. 44.51 in Bmad manual, handled by GTPSA):
  map!(t->t.tpsa, Fx_low, F.x)
  map!(t->t.tpsa, Hx_low, H.x)
  map!(t->t.tpsa, Gx_low, G.x)
  liebra!(nv, Fx_low, Hx_low, Gx_low)

  # Spin (Eq. 44.51 in Bmad manual, NOT handled by GTPSA):
  if !isnothing(F.Q)
    @assert !(G === F) && !(G === H) "Aliasing the solution in lb! with spin is not allowed"
    # first [h,f] (this part involves quaternion multiplication so will be slow)
    mul!(G.Q, H.Q, F.Q)
    mul!(work_Q, F.Q, H.Q)
    for i=1:4
      sub!(G.Q.q[i], G.Q.q[i], work_Q[i])
    end

    # then +F⋅∇h 
    tmp = work_Q.Q.q[1]
    for i=1:4
      fgrad!(tmp, F.x, H.Q.q[i], work_low=Fx_low)
      add!(G.Q.q[i], G.Q.q[i], tmp)
    end
    
    # finally -G⋅∇f
    for i=1:4
      fgrad!(tmp, G.x, F.Q.q[i], work_low=Fx_low)
      sub!(G.Q.q[i], G.Q.q[i], tmp)
    end
  end
  return
end


# --- Lie operator (VectorField) acting on DAMap ---
function *(F::VectorField{T,U}, m1::DAMap{S,T,U,V}) where {S,T,U,V}
  m = zero(m1)
  mul!(m, F, m1)
  return m
end


"""


"""
function mul!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(F)

  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low. Received $(eltype(work_low)), should be $(lowtype(T))"
  @assert !(m === m1) "Aliasing m === m1 is not allowed for mul!"


  # Orbital part F⋅∇m1 :
  for i=1:nv
    fgrad!(m.x[i], F.x, m1.x[i], work_low=work_low)
  end

  # Spin F⋅∇q + qf (quaternion part acts in reverse order):
  if !isnothing(F.Q)
    mul!(work_Q, m1.Q, F.Q)


    for i=1:4
      fgrad!(m.Q.q[i], F.x, m1.Q.q[i], work_low=work_low)
    end

    
    for i=1:4
      add!(m.Q.q[i], m.Q.q[i], work_Q.q[i])
    end
  end
  return
end


# --- exp(F)*m ---
function exp(F::VectorField{T,U}, m::Union{UniformScaling,DAMap{S,T,U,V}}=I) where {S,T,U,V}
  if m isa UniformScaling 
    m = one(DAMap{numtype(T),T,U,Nothing},use=getdesc(F))
  end
  mf = zero(m)
  exp!(mf, F, m)
  return mf
end

function exp!(m::DAMap{S,T,U,V}, F::VectorField{T,U}, m1::DAMap{S,T,U,V}; work_map::DAMap{S,T,U,V}=zero(m1), work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(first(F.x))}(undef, numvars(F)), work_Q::Union{U,Nothing}=prep_vf_work_Q(F)) where {S,T,U,V}
  nv = numvars(F)

  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low. Received $(eltype(work_low)), should be $(lowtype(T))"

 # exppb!(nv, map(t->t.tpsa, F.x), map(t->t.tpsa, m1.x), map(t->t.tpsa,m.x))
#return
  # The convergence checks are taken exactly from GTPSA
  
  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(numtype(T))*nv
  nrm_ =Inf
  conv = false

  tmp = work_map
  tmp2 = zero(m1)

  copy!(tmp, m1)
  copy!(m, m1)

  for j=1:nmax
    div!(tmp, tmp, j)
    mul!(tmp2, F, tmp)
    add!(m, m, tmp2)
    copy!(tmp, tmp2)
    nrm = 0
    for i=1:nv
      nrm += norm(tmp2.x[i])
    end

    # Check convergence
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      #println("converged at j = ", j)
      return
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  #=

  tmp = TPS()
  tmp2 = TPS()

  for i=1:nv
    nrm_ = typemax(Float64)
    conv = false

    copy!(tmp, m1.x[i])
    #println(tmp)
    copy!(m.x[i], m1.x[i])
    nrm = 0
j=1
    while j <= 100
      #println(j)
      #println("F = ", F)
      #println("tmp = ", tmp)
      tmp /= j
      fgrad!(tmp2, F.x, tmp)
      m.x[i] += tmp2
      #println("m.x[$i] = ", m)
      copy!(tmp, tmp2)
      nrm = norm(tmp2)
      #println(tmp2)
      #println(nrm)

      # Check convergence
      if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
       # println("converged at j = ", j)
        #copy!(m, tmp)
        j = 101
      end

      if nrm <= nrm_min1 # convergence reached just refine a bit
        conv = true
      end
      nrm_ = nrm
      j += 1
    end
  end=#
#=



  for i=1:nmax
    work_F = work_F / i
    # Orbital part F⋅∇m1 :
    for i=1:nv
      fgrad!(Cint(nv),map(t->t.tpsa, work_F.x), work_map.x[i].tpsa, tmp2.x[i].tpsa)
      #fgrad!(work_map.x[i], F.x, work_map.x[i], work_low=work_low)
    end

    tmp += tmp2
    copy!(work_map, tmp2)
    #work_map = work_F*work_map
    #tmp += work_map

    #div!(work_F,work_F,i)
    #println("i = ", i)
    #println("F = ", work_F)
    #println("work_map")
    #mul!(tmp,work_F,work_map,work_low=work_low,work_Q=work_Q)
    #add!(m,m,tmp)


    #copy!(work_map,tmp)
    # Check norm of orbital part (also done no quaternion in FPP)
    nrm = 0
    for i=1:nv
      nrm += norm(tmp2.x[i])
    end
    println(nrm)

    # Check convergence
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      copy!(m, tmp)
      return
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
  end

  nrm_ = nrm
  @warn "exp! convergence not reached for $nmax iterations"=#
end


# --- log ---

"""

Instead of using GTPSA's `logpb`, we need a log that will work with `VectorField`s including 
a quaternion. Therefore I am writing my own.

If we formulate the matrix of monomials, the log is easy. But this is gross and slow.

Therefore we use the fast algorithm described in Ch. 47 of Bmad manual


"""
function log!(mf::DAMap{S,T,U,V}, m::DAMap{S,T,U,V}) where {S,T,U,V}

end



















logpb!(na::Cint, ma::Vector{Ptr{RTPSA}}, mb::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_logpb!(na, ma, mb, mc))
logpb!(na::Cint, ma::Vector{Ptr{CTPSA}}, mb::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_logpb!(na, ma, mb, mc))




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
