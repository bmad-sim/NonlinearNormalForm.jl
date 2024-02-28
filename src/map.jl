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
      x0 = deepcopy(m.x0)
      x = map(x->(T)(x, use=getdesc(use)), m.x)
      Q = Quaternion(map(x->(U)(x, use=getdesc(use)), m.Q.q))
      E = deepcopy(m.E)
      return $t{S,T,U,V}(x0, x, Q, E)
    end

    # Create TaylorMap from a Probe (must tag on Parameters)
    function $t(m::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S, T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS}, V}
      x0 = deepcopy(m.x0)
      v = map(x->(T)(x, use=getdesc(use)), m.x)
      if eltype(v) == TPS
        k = params(getdesc(v))
      else
        k = complexparams(getdesc(v))
      end
      x = vcat(v,k)
      Q = Quaternion(map(x->(U)(x, use=getdesc(use)), m.Q.q))
      E = deepcopy(m.E)
      return $t{S,T,U,V}(x0, x, Q, E)
    end
    
    # Create from vector and blank (in this case must check for consistency)
    function $t(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::Quaternion{U}=Quaternion(first(x)), E::Matrix{V}=zeros(numtype(first(x)), numvars(x), numvars(x)), use::Union{Descriptor,<:TaylorMap,Probe{S,TPS,TPS,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS},V}
      numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
      numvars(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")
      (numvars(x),numvars(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      (isnothing(use) && getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      
      x1 = map(x->(T)(x, use=getdesc(use)), x)
      Q1 = Quaternion(map(x->(U)(x, use=getdesc(use)), Q.q))
      
      return $t{eltype(x0),T,U,eltype(E)}(deepcopy(x0), x1, Q1, deepcopy(E))
    end
  end
end


# --- composition ---
function ∘(m2::DAMap{S2,T2,U2,V2},m1::DAMap{S1,T1,U1,V1}) where {S2,T2,U2,V2,S1,T1,U1,V1}
  # all(scalar.(m1.x)+m1.x0 - m2.x0 .< 1e-20) || error("Disconnected DAMaps! Exit coordinates of first map != entrance coordinates of second map!")
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)
  ref = Vector{numtype(T1)}(undef, nv)

  outT = promote_type(T2,T1)

  outx = Vector{outT}(undef,  nv+np)
  # Set up orbit out:
  for i=1:nv
    @inbounds outx[i] = outT(use=desc)
  end

  outU = promote_type(U2,U1)

  outQ = Quaternion{outU}([outU(use=desc), outU(use=desc), outU(use=desc), outU(use=desc)])

  # Take out scalar part and store it
  for i=1:nv
    @inbounds ref[i] = m1.x[i][0]
    @inbounds m1.x[i][0] = 0
  end

  # Do the composition, promoting if necessary
  # --- Orbit ---
  if outT != T1
    # Promote to ComplexTPS:
    m1x_store = map(t->outT(t), view(m1.x, :)) # Must store to prevent GC   
    m1x_low = map(t->t.tpsa, m1x_store)
  else
    m1x_store = nothing
    m1x_low = map(t->t.tpsa, view(m1.x, :))
  end
  if outT != T2
    m2x_store = map(t->outT(t), view(m2.x, 1:nv)) # Must store to prevent GC   
    m2x_low = map(t->t.tpsa, m2x_store)
  else
    m2x_store = nothing
    m2x_low = map(t->t.tpsa, view(m2.x, 1:nv))
  end

  # --- Quaternion ---
  if outU != U2
    m2Q_store = map(t->outU(t), view(m2.Q.q, :)) # Must store to prevent GC
    m2Q_low = map(t->t.tpsa, m2Q_store)
  else
    m2Q_store = nothing
    m2Q_low = map(t->t.tpsa, view(m2.Q.q, :))
  end

  # go low
  outx_low = map(t->t.tpsa, view(outx, 1:nv))
  outQ_low = map(t->t.tpsa, outQ.q)

  # Orbit:
  GC.@preserve m1x_store m2x_store compose!(nv, m2x_low, nv+np, m1x_low, outx_low)

  # Spin (spectator) q(z0)=q2(M(z0))q1(z0)
  # First obtain q2(M(z0))
  GC.@preserve m1x_store m2Q_store compose!(Cint(4), m2Q_low, nv, m1x_low, outQ_low)
  # Now concatenate
  qmul!(outQ, m1.Q, outQ)

  # Add the params to outx:
  if outT == TPS
    k = params(desc)
  else
    k = complexparams(desc)
  end
  outx[nv+1:end] = k
  
  # Put back the reference and if m1 === m2, also add to outx
  if m1 === m2
    for i=1:nv
      @inbounds m1.x[i][0] = ref[i]
      @inbounds outx[i][0] += ref[i]
    end
  else
    for i=1:nv
      @inbounds m1.x[i][0] = ref[i]
    end
  end


  # Make that map!
  return DAMap(deepcopy(m1.x0), outx, outQ, zeros(nv, nv))
end

function ∘(m2::TPSAMap{S2,T2,U2,V2},m1::TPSAMap{S1,T1,U1,V1}) where {S2,T2,U2,V2,S1,T1,U1,V1}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc)

  outT = promote_type(T2,T1)

  outx = Vector{outT}(undef,  nv+np)
  # Set up orbit out:
  for i=1:nv
    @inbounds outx[i] = outT(use=desc)
  end
  
  outU = promote_type(U2,U1)

  outQ = Quaternion{outU}([outU(use=desc), outU(use=desc), outU(use=desc), outU(use=desc)])

  # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  for i=1:nv
    @inbounds m1.x[i] -= m2.x0[i]
  end

  # Do the composition, promoting if necessary
  # --- Orbit ---
  if outT != T1
    # Promote to ComplexTPS:
    m1x_store = map(t->outT(t), view(m1.x, :)) # Must store to prevent GC   
    m1x_low = map(t->t.tpsa, m1x_store)
  else
    m1x_store = nothing
    m1x_low = map(t->t.tpsa, view(m1.x, :))
  end
  if outT != T2
    m2x_store = map(t->outT(t), view(m2.x, 1:nv)) # Must store to prevent GC   
    m2x_low = map(t->t.tpsa, m2x_store)
  else
    m2x_store = nothing
    m2x_low = map(t->t.tpsa, view(m2.x, 1:nv))
  end

  # --- Quaternion ---
  if outU != U2
    m2Q_store = map(t->outU(t), view(m2.Q.q, :)) # Must store to prevent GC
    m2Q_low = map(t->t.tpsa, m2Q_store)
  else
    m2Q_store = nothing
    m2Q_low = map(t->t.tpsa, view(m2.Q.q, :))
  end

  # go low
  outx_low = map(t->t.tpsa, view(outx, 1:nv))
  outQ_low = map(t->t.tpsa, outQ.q)

  # Orbit:
  GC.@preserve m1x_store m2x_store compose!(nv, m2x_low, nv+np, m1x_low, outx_low)

  # Spin (spectator) q(z0)=q2(M(z0))q1(z0)
  # First obtain q2(M(z0))
  GC.@preserve m1x_store m2Q_store compose!(Cint(4), m2Q_low, nv, m1x_low, outQ_low)
  # Now concatenate
  qmul!(outQ, m1.Q, outQ)

  # Add the params to outx:
  if outT == TPS
    k = params(desc)
  else
    k = complexparams(desc)
  end
  outx[nv+1:end] = k

  # Now fix m1 and if m2 === m1, add to output too:
 # For TPSA Map concatenation, we need to subtract w_0 (m2 x0) (Eq. 33)
  # Because we are still expressing in terms of z_0 (m1 x0)
  if m1 === m2
    for i=1:nv
      @inbounds m1.x[i] += m2.x0[i]
      @inbounds outx[i] += m2.x0[i]
    end
  else
    for i=1:nv
      @inbounds m1.x[i] += m2.x0[i]
    end
  end

  # Make that map!
  return TPSAMap(deepcopy(m1.x0), outx, outQ, zeros(nv, nv))
end

# --- inverse ---
function inv(m1::DAMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc) 
  outx = Vector{T}(undef, nv+np)
  for i=1:nv
    @inbounds outx[i] = T(use=desc)
  end
  # 
  outQ = inv(m1.Q)

  m1x_low = map(t->t.tpsa, view(m1.x, 1:nv))
  outx_low = map(t->t.tpsa, view(outx, 1:nv))

  # This C function ignores the scalar part so no need to take it out
  minv!(nv, m1x_low, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  outQ_low = map(t->t.tpsa, outQ.q)
  compose!(Cint(4), outQ_low, nv, outx_low, outQ_low)

  # Add the params to outx:
  if T == TPS
    k = params(desc)
  else
    k = complexparams(desc)
  end

  outx0 = Vector{numtype(T)}(undef, nv)
  for i=1:nv
    @inbounds outx0[i] = m1.x[i][0]
    @inbounds outx[i] += m1.x0[i]
  end
  outx[nv+1:end] = k
  
  return DAMap(outx0, outx, outQ, zeros(nv,nv))
end

function inv(m1::TPSAMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m1)
  nv = numvars(desc)
  np = numparams(desc) 
  outx = Vector{T}(undef, nv+np)
  for i=1:nv
    @inbounds outx[i] = T(use=desc)
  end
  # 
  outQ = inv(m1.Q)

  m1x_low = map(t->t.tpsa, view(m1.x, 1:nv))
  outx_low = map(t->t.tpsa, view(outx, 1:nv))

  # This C function ignores the scalar part so no need to take it out
  minv!(nv, m1x_low, outx_low)

  # Now do quaternion: inverse of q(z0) is q^-1(M^-1(zf))
  outQ_low = map(t->t.tpsa, outQ.q)
  compose!(Cint(4), outQ_low, nv, outx_low, outQ_low)

  # Add the params to outx:
  if T == TPS
    k = params(desc)
  else
    k = complexparams(desc)
  end

  outx0 = Vector{numtype(T)}(undef, nv)
  for i=1:nv
    @inbounds outx0[i] = m1.x[i][0]
    @inbounds outx[i] += m1.x0[i]
  end
  outx[nv+1:end] = k

  return TPSAMap(outx0, outx, outQ, zeros(nv,nv))
end

function ^(m1::TaylorMap{S,T,U,V}, n::Integer) where {S,T,U,V}
  if n>0
    m = m1
    for i=1:n-1
      m = m1∘m
    end
    return m
  elseif n<0
    m = m1
    for i=1:-n+1
      m = m1∘m
    end
    return inv(m)
  else
    return (typeof(m1)){S,T,U,V}(m1)
  end
end

==(m1::TaylorMap, m2::TaylorMap) = (m1.x0 == m2.x0 && m1.x == m2.x && m1.Q == m2.Q && m1.E == m2.E)