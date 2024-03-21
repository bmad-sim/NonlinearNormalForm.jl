#=
TPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.

Once again parametric types, however being a TaylorMap we now require the 
orbit x and spin q to be TPSs.
=#
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V} end 

struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V} <: TaylorMap{S,T,U,V}
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V} <: TaylorMap{S,T,U,V}
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

for t = (:DAMap, :TPSAMap)
@eval begin

# Create a TaylorMap from other TaylorMap
function $t(m::TaylorMap{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S, T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS}, V}
  desc = getdesc(m)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np+nv
  x = Vector{T}(undef, nn)
  q = Vector{U}(undef, 4)
  for i=1:nv
    @inbounds x[i] = T(m.x[i], use=getdesc(use))
  end

  # use same parameters if same descriptor (use=nothing)
  if isnothing(use)
    @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)
  else
    if T == TPS
      @inbounds x[nv+1:nn] = params(getdesc(first(x)))
    else
      @inbounds x[nv+1:nn] = complexparams(getdesc(first(x)))
    end
  end

  for i=1:4
    @inbounds q[i] = U(m.Q.q[i],use=getdesc(use))
  end
  return $t{S,T,U,V}(copy(m.x0), x, Quaternion(q), copy(m.E))
end

# Create TaylorMap from a Probe. Probes do not have parameters tacked on so must allocate new
function $t(p::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S, T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS}, V}
  desc = getdesc(p)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np+nv
  x = Vector{T}(undef, nn)
  q = Vector{U}(undef, 4)
  for i=1:nv
    @inbounds x[i] = T(p.x[i], use=getdesc(use))
  end

  # Must allocate new for parameters
  if T == TPS
    @inbounds x[nv+1:nn] = params(getdesc(first(x)))
  else
    @inbounds x[nv+1:nn] = complexparams(getdesc(first(x)))
  end

  for i=1:4
    @inbounds q[i] = U(p.Q.q[i],use=getdesc(use))
  end
  return $t{S,T,U,V}(copy(p.x0), x, Quaternion(q), copy(p.E))
end

# Create from vector and blank (in this case must check for consistency)
function $t(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::Quaternion{U}=Quaternion(first(x)), E::Matrix{V}=zeros(numtype(first(x)), numvars(x), numvars(x)), use::Union{Descriptor,<:TaylorMap,Probe{S,TPS,TPS,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS},V}
  numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
  numvars(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")
  (numvars(x),numvars(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
  (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
  nv = numvars(x)
  np = numparams(x)
  nn = np+nv
  x1 = Vector{T}(undef, nn)
  x1[1:nv] = map(x->(T)(x, use=getdesc(use)), x)
  if T == TPS
    x1[nv+1:nn] = params(getdesc(first(x)))
  else
    x1[nv+1:nn] = complexparams(getdesc(first(x)))
  end
  Q1 = Quaternion(map(x->(U)(x, use=getdesc(use)), Q.q))
  return $t{eltype(x0),T,U,eltype(E)}(copy(x0), x1, Q1, copy(E))
end

# Basic operators
function +(m2::$t,m1::$t)
  return $t(m2.x0+m1.x0, m2.x+m1.x, Quaternion(m2.Q.q+m1.Q.q), m2.E+m1.E)
end

function -(m2::$t,m1::$t)
  return $t(m2.x0-m1.x0, m2.x-m1.x, Quaternion(m2.Q.q-m1.Q.q), m2.E-m1.E)
end

# zero map (empty but still identity in parameters)
function zero(m::$t)
  desc = getdesc(m)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np+nv
  x = Vector{eltype(m.x)}(undef, nn)
  q = Vector{eltype(m.Q.q)}(undef, 4)
  for i=1:nv
    @inbounds x[i] = zero(m.x[i])
  end

  # use same parameters 
  @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)

  for i=1:4
    @inbounds q[i] = zero(m.Q.q[i])
  end
  return $t(copy(m.x0), x, Quaternion(q), copy(m.E))
end


"""
    compose_it!(m, m2, m1)

Composes the maps  in-place, `m2 ∘ m1`. Aliasing is allowed, however in some cases may 
be slower than not aliasing, especially when `m === m1`. Assumes the destination map is 
properly set up (with correct types promoted if necessary), and that `m.x[1:nv]` and 
`m.Q.q` contain allocated TPSs.
""" 
function compose_it!(m::$t, m2::$t, m1::$t)
  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  
  outT = eltype(m.x)
  outU = eltype(m.Q.q)
  m.x0 .= m1.x0
  m.E .= zeros(nv,nv)

  # add immutable parameters to outx
  if eltype(m.x) == eltype(m1.x)
    @inbounds m.x[nv+1:nn] = view(m1.x, nv+1:nn)
  else
    @inbounds m.x[nv+1:nn] = view(m2.x, nv+1:nn)
  end

  # Do the composition, promoting if necessary
  # --- Orbit ---
  if outT != eltype(m1.x) #T1
    m1x_store = Vector{ComplexTPS}(undef, nn)
    # Promote to ComplexTPS:
    map!(t->ComplexTPS(t), m1x_store,  m1.x) # Must store to prevent GC   
    m1x_low = map(t->t.tpsa, m1x_store)
    if m1 === m
      m1xquat_store = Vector{ComplexTPS}(undef, nn)
      map!(t->ComplexTPS(t), m1xquat_store, m1x_store) 
      m1xquat_low = map(t->t.tpsa, m1xquat_store)
      m1Q = Quaternion(m1.Q)
    else
      m1xquat_store = m1x_store
      m1xquat_low = m1x_low
      m1Q = m1.Q
    end
  else
    m1x_store = nothing
    m1x_low = map(t->t.tpsa, m1.x)
    if m1 === m
      m1xquat_store = Vector{outT}(undef, nn)
      map!(t->(eltype(m1.x))(t), m1xquat_store, m1.x) 
      m1xquat_low = map(t->t.tpsa, m1xquat_store)
      m1Q = Quaternion(m1.Q)
    else
      m1xquat_store = nothing
      m1xquat_low = m1x_low
      m1Q = m1.Q
    end
  end
  if outT != eltype(m2.x)# T2
    m2x_store = Vector{ComplexTPS}(undef, nv)
    # Promote to ComplexTPS:
    map!(t->ComplexTPS(t), m2x_store, view(m2.x,1:nv)) # Must store to prevent GC   
    m2x_low = map(t->t.tpsa, m2x_store)
  else
    m2x_store = nothing
    m2x_low = map(t->t.tpsa, view(m2.x, 1:nv))
  end

  # --- Quaternion ---
  if outU != eltype(m2.Q.q) #U2
    m2Q_store = Vector{ComplexTPS}(undef, 4)
    map!(t->ComplexTPS(t), m2Q_store, m2.Q.q) # Must store to prevent GC
    m2Q_low = map(t->t.tpsa, m2Q_store)
  else
    m2Q_store = nothing
    m2Q_low = map(t->t.tpsa, m2.Q.q)
  end

  # go low
  outx_low = map(t->t.tpsa, view(m.x,1:nv))
  outQ_low = map(t->t.tpsa, m.Q.q)

  # Orbit:
  GC.@preserve m1x_store m2x_store compose!(nv, m2x_low, nv+np, m1x_low, outx_low)

  # Spin (spectator) q(z0)=q2(M(z0))q1(z0)
  # First obtain q2(M(z0))
  GC.@preserve m1xquat_store m2Q_store compose!(Cint(4), m2Q_low, nv+np, m1xquat_low, outQ_low)
  # Now concatenate
  mul!(m.Q, m1Q, m.Q)
  # Make that map
  return 
end

function compose_it(m2::$t, m1::$t)
  desc = getdesc(m1)
  nn = numnn(desc)
  nv = numvars(desc)

  # Set up outx:
  outT = promote_type(eltype(m2.x),eltype(m1.x))
  outx = Vector{outT}(undef, nn)
  for i=1:nv
     @inbounds outx[i] = outT(use=desc)
  end

  # set up quaternion out:
  outU = promote_type(eltype(m2.Q.q), eltype(m1.Q.q))
  outQ = Quaternion{outU}([outU(use=desc), outU(use=desc), outU(use=desc), outU(use=desc)])

  # make a map
  m = $t(Vector{eltype(m1.x0)}(undef, nv), outx, outQ, Matrix{eltype(m1.E)}(undef, nv, nv))

  # compose
  compose_it!(m, m2, m1)

  return m
end

end
end

==(m1::TaylorMap, m2::TaylorMap) = (m1.x0 == m2.x0 && m1.x == m2.x && m1.Q == m2.Q && m1.E == m2.E)

# --- compose ---
function compose!(m::DAMap, m2::DAMap, m1::DAMap)
  # DAMap setup:
  desc = getdesc(m1)
  nv = numvars(desc)
  ref = Vector{numtype(eltype(m1.x))}(undef, nv)
  # Take out scalar part and store it
  for i=1:nv
      @inbounds ref[i] = m1.x[i][0]
      @inbounds m1.x[i][0] = 0
  end

  compose_it!(m, m2, m1)

  # Put back the reference and if m1 === m2, also add to outx
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

  return m
end

function ∘(m2::DAMap,m1::DAMap)
  # DAMap setup:
  desc = getdesc(m1)
  nv = numvars(desc)
  ref = Vector{numtype(eltype(m1.x))}(undef, nv)
  # Take out scalar part and store it
  for i=1:nv
     @inbounds ref[i] = m1.x[i][0]
     @inbounds m1.x[i][0] = 0
  end
  
  m = compose_it(m2,m1)
  
  # Put back the reference and if m1 === m2, also add to outx
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

  return m
end

function ∘(m2::TPSAMap, m1::TPSAMap)
  # TPSAMap setup:
  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  desc = getdesc(m1)
  nv = numvars(desc)
  for i=1:nv
    @inbounds m1.x[i] -= m2.x0[i]
  end

  m = compose_it(m2,m1)
  
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

# Also allow * for simpliticty
*(m2::DAMap,m1::DAMap) = ∘(m2,m1)
*(m2::TPSAMap,m1::TPSAMap) = ∘(m2,m1)

# --- inverse ---
minv!(na::Cint, ma::Vector{Ptr{RTPSA}}, nb::Cint, mc::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_minv!(na, ma, nb, mc))
minv!(na::Cint, ma::Vector{Ptr{CTPSA}}, nb::Cint, mc::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_minv!(na, ma, nb, mc))

function inv(m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np + nv
  outx = Vector{T}(undef, nn)
  for i=1:nv
     @inbounds outx[i] = T(use=desc)
  end
  # add immutable parameters to outx
  @inbounds outx[nv+1:nn] = view(m1.x, nv+1:nn)

  outQ = inv(m1.Q)

  m1x_low = map(t->t.tpsa, m1.x)
  outx_low = map(t->t.tpsa, outx)

  # This C function ignores the scalar part so no need to take it out
  minv!(nv+np, m1x_low, nv, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  outQ_low = map(t->t.tpsa, outQ.q)
  compose!(Cint(4), outQ_low, nv+np, outx_low, outQ_low)

  outx0 = Vector{S}(undef, nv)
  for i=1:nv
     @inbounds outx0[i] = m1.x[i][0]
     @inbounds outx[i] += m1.x0[i]
  end
  
  return (typeof(m1))(outx0, outx, outQ, zeros(nv,nv))
end

function inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np + nv

  # add immutable parameters to outx
  @inbounds m.x[nv+1:nn] = view(m1.x, nv+1:nn)

  inv!(m.Q, m1.Q)

  m1x_low = map(t->t.tpsa, m1.x)
  outx_low = map(t->t.tpsa, m.x)

  # This C function ignores the scalar part so no need to take it out
  minv!(nv+np, m1x_low, nv, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  outQ_low = map(t->t.tpsa, m.Q.q)
  compose!(Cint(4), outQ_low, nv+np, outx_low, outQ_low)
  for i=1:nv
     @inbounds m.x0[i] = m1.x[i][0]
     @inbounds m.x[i] += m1.x0[i]
  end
  
  return 
end

literal_pow(::typeof(^), m::TaylorMap{S,T,U,V}, vn::Val{n}) where {S,T,U,V,n} = ^(m,n)

function ^(m1::TaylorMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  if n>0
    m = m1
    for i=1:(n-1)
      m = m1∘m
    end
    return m
  elseif n<0
    m = m1
    for i=1:(-n-1)
      m = m1∘m
    end
    return inv(m)
  else
    return (typeof(m1)){S,T,U,V}(m1)
  end
end

function cut(m1::TaylorMap{S,T,U,V}, order::Integer) where {S,T,U,V}
  m = zero(m1)
  cut!(m, m1, order)
  return m
end

function cut!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}, order::Integer) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  nn = np + nv
  m.x0 .= m1.x0
  ord = Cint(order)
  for i=1:nv
    GTPSA.cutord!(m1.x[i].tpsa, m.x[i].tpsa, convert(Cint, ord))
  end
  # add immutable parameters to outx
  @inbounds m.x[nv+1:nn] = view(m1.x, nv+1:nn)
  for i=1:4
    @inbounds GTPSA.cutord!(m1.Q.q[i].tpsa, m.Q.q[i].tpsa, convert(Cint, ord))
  end
  m.E .= m.E
  return
end