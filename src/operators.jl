# --- compose ---
function compose!(m::DAMap, m2::DAMap, m1::DAMap; keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=nothing)
  # DAMap setup:
  desc = getdesc(m1)
  nv = numvars(desc)
  if isnothing(work_ref)
    ref = Vector{numtype(eltype(m1.x))}(undef, numvars(m1))
  else
    @assert length(work_ref) >= nv "Incorrect length for work_ref, received $(length(work_ref)) but should be atleast $nv"
    ref = work_ref
  end   

  if keep_scalar
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

  compose_it!(m, m2, m1, work_low=work_low, work_prom=work_prom)

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

function compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=nothing)
  # TPSAMap setup:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  desc = getdesc(m1)
  nv = numvars(desc)
  for i=1:nv
    @inbounds m1.x[i] -= m2.x0[i]
  end

  compose_it!(m, m2, m1, work_low=work_low, work_prom=work_prom)

  # Now fix m1 and if m2 === m1, add to output too:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  if m1 === m2
    for i=1:nv
       @inbounds m1.x[i] += m2.x0[i]
       @inbounds m.x[i] += m2.x0[i]
    end
  else
    for i=1:nv
       @inbounds m1.x[i] += m2.x0[i]
    end
  end

  return m
end

# --- inverse ---
minv!(na::Cint, ma::Vector{Ptr{RTPSA}}, nb::Cint, mc::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_minv!(na, ma, nb, mc))
minv!(na::Cint, ma::Vector{Ptr{CTPSA}}, nb::Cint, mc::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_minv!(na, ma, nb, mc))


function inv(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  m = zero(m1)
  inv!(m,m1)
  return m
end

function inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}; work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing) where {S,T,U,V}
  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)

  # Set up work:
  if isnothing(work_low)
    outx_low = Vector{lowtype(T)}(undef, nn)    # SHOULD ONLY NEED TO BE NV  BUT GTPSA BUG
    m1x_low = Vector{lowtype(T)}(undef, nn)
    if !isnothing(m.Q)
      if nn >= 4   # reuse
        outQ_low = m1x_low
      else
        outQ_low = Vector{lowtype(T)}(undef, 4)
      end
    end
  else
    outx_low = work_low[1]
    m1x_low = work_low[2]
    @assert length(outx_low) >= nn "Cannot inv!: incorrect length for outx_low = work_low[1]. Received $(length(outx_low)), must be >= $nn"
    @assert length(m1x_low) >= nn "Cannot inv!: incorrect length for m1x_low = work_low[2]. Received $(length(m1x_low)), must be >= $nn"
    if !isnothing(m.Q)
      outQ_low = work_low[3]
      @assert length(outQ_low) >= 4 "Cannot inv!: incorrect length for outQ_low = work_low[3]. Received $(length(outQ_low)), must be >= 4"
      @assert !(outQ_low === outx_low) "Cannot inv!: outQ_low === outx_low !! outx_low must NOT be reused!"
    end
  end
  # if aliasing, must pass vector to store x0
  if m1 === m
    if isnothing(work_ref)
      ref = Vector{numtype(T)}(undef, nv)
    else
      ref = work_ref
      @assert length(ref) >= nv "Cannot inv!: incorrect length for ref. Received $(length(ref)), must be >= $nv"
    end
    map!(t->t[0], ref, view(m1.x,1:nv))
  end
  # add immutable parameters to outx
  @inbounds m.x[nv+1:nn] = view(m1.x, nv+1:nn)

  map!(t->t.tpsa, m1x_low, m1.x)
  map!(t->t.tpsa, outx_low, m.x)

  # This C function ignores the scalar part so no need to take it out
  minv!(nn, m1x_low, nv, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  if !isnothing(m.Q)
    inv!(m.Q, m1.Q)
    map!(t->t.tpsa, outQ_low, m.Q.q)
    compose!(Cint(4), outQ_low, nn, outx_low, outQ_low)
  end

  if m1 === m
    for i=1:nv
      @inbounds m.x[i][0] = m1.x0[i]
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


# --- powers ---

function ^(m1::DAMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  nv = numvars(m1)
  # Do it
  if n>0
    #prepare work
    work_low = prep_comp_work_low(m1)
    ref = Vector{numtype(T)}(undef, nv)

    # do it:
    m = DAMap(m1)
    map!(t->t[0], ref, view(m.x,1:nv))
    for i=1:(n-1)
      compose!(m,m,m1, keep_scalar=false,work_low=work_low)
    end
    for i=1:nv
        @inbounds m1.x[i][0] = ref[i]
        #@inbounds m.x[i][0] += ref[i]
    end
    return m
  elseif n<0
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m1)
    ref = Vector{numtype(T)}(undef, nv)

    m = DAMap(m1)
    map!(t->t[0], ref, view(m.x,1:nv))
    for i=1:(-n-1)
      compose!(m,m,m1, keep_scalar=false,work_low=comp_work_low)
    end
    for i=1:nv
        @inbounds m1.x[i][0] = ref[i]
        #@inbounds m.x[i][0] += ref[i]
    end
    inv!(m,m,work_ref=ref,work_low=inv_work_low)
    return m
  else
    return DAMap(m1)
  end
end

function ^(m1::TPSAMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  # Do it
  if n>0
    work_low = prep_comp_work_low(m1)
    m = TPSAMap(m1)
    for i=1:(n-1)
      compose!(m,m,m1,work_low=work_low)
    end
    return m
  elseif n<0
    #prepare work
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m1)
    m = TPSAMap(m1)
    for i=1:(-n-1)
      compose!(m,m,m1,work_low=comp_work_low)
    end
    inv!(m,m,work_low=inv_work_low)
    return m
  else
    return TPSAMap(m1)
  end
end

# --- others ---

for t = (:DAMap, :TPSAMap)
@eval begin

function ∘(m2::$t,m1::$t)
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

# Basic operators
function +(m2::$t{S,T,U,V},m1::$t{S,T,U,V}) where {S,T,U,V}
  if xor(isnothing(m2.Q), isnothing(m1.Q))
    error("Cannot +: one map includes spin, the other does not")
  end
  if xor(isnothing(m2.E), isnothing(m1.E))
    error("Cannot +: one map includes radiation, the other does not")
  end

  if !isnothing(m2.Q)
    Q1 = Quaternion(m2.Q.q+m1.Q.q)
  else
    Q1 = nothing
  end

  if !isnothing(m2.E)
    E1 =  m2.E+m1.E
  else
    E1 = nothing
  end

  return $t(m2.x0+m1.x0, m2.x+m1.x, Q1, E1)
end

function -(m2::$t,m1::$t)
  if xor(isnothing(m2.Q), isnothing(m1.Q))
    error("Cannot -: one map includes spin, the other does not")
  end
  if xor(isnothing(m2.E), isnothing(m1.E))
    error("Cannot -: one map includes radiation, the other does not")
  end

  if !isnothing(m2.Q)
    Q1 = Quaternion(m2.Q.q-m1.Q.q)
  else
    Q1 = nothing
  end

  if !isnothing(m2.E)
    E1 =  m2.E-m1.E
  else
    E1 = nothing
  end

  return $t(m2.x0-m1.x0, m2.x-m1.x, Q1, E1)
end

literal_pow(::typeof(^), m::$t{S,T,U,V}, vn::Val{n}) where {S,T,U,V,n} = ^(m,n)

# Also allow * for simpliticty and \ and / because why not
*(m2::$t,m1::$t) = ∘(m2,m1)
/(m2::$t,m1::$t) = m2∘inv(m1) 
\(m2::$t,m1::$t) = inv(m2)∘m1

# Uniform scaling (identity map)
+(m::$t, J::UniformScaling) = m + one(m)
+(J::UniformScaling, m::$t) = +(m,J)
-(m::$t, J::UniformScaling) = m - one(m)
-(J::UniformScaling, m::$t) = one(m) - m
∘(m::$t, J::UniformScaling) = DAMap(m)
∘(J::UniformScaling, m::$t) = ∘(m,J)
*(m::$t, J::UniformScaling) = ∘(m,J)
*(J::UniformScaling, m::$t) = ∘(m,J)
/(m::$t, J::UniformScaling) = DAMap(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = DAMap(m)

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