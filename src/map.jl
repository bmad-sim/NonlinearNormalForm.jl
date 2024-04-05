#=
TPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.

Once again parametric types, however being a TaylorMap we now require the 
orbit x and spin q to be TPSs. For no spin tracking, U==Nothing, for no radiation,
V == Nothing.
=#
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} end 

struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}  <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

for t = (:DAMap, :TPSAMap)
@eval begin

# Create a TaylorMap from other TaylorMap
function $t(m::TaylorMap{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}
  desc = getdesc(m)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  x = Vector{T}(undef, nn)
  
  for i=1:nv
    @inbounds x[i] = T(m.x[i], use=getdesc(use))
  end

  # use same parameters if same descriptor (use=nothing)
  if isnothing(use) || getdesc(use) == desc
    @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)
  else
    if T == TPS
      @inbounds x[nv+1:nn] = params(getdesc(first(x)))
    else
      @inbounds x[nv+1:nn] = complexparams(getdesc(first(x)))
    end
  end

  if !isnothing(m.Q)
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(m.Q.q[i],use=getdesc(use))
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = copy(m.E)
  else
    E = nothing
  end

  return $t{S,T,U,V}(copy(m.x0), x, Q, E)
end

# Create TaylorMap from a Probe. Probes do not have parameters tacked on so must allocate new
function $t(p::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}
  desc = getdesc(p)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  x = Vector{T}(undef, nn)

  length(p.x) == nv || error("Length of orbital ray ($(length(p.x))) inconsistent with number of variables in GTPSA ($(nv))")
  
  for i=1:nv
    @inbounds x[i] = T(p.x[i], use=getdesc(use))
  end

  if T == TPS
    @inbounds x[nv+1:nn] = params(getdesc(first(x)))
  else
    @inbounds x[nv+1:nn] = complexparams(getdesc(first(x)))
  end


  if !isnothing(p.Q)
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(p.Q.q[i],use=getdesc(use))
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(p.E)
    E = copy(p.E)
  else
    E = nothing
  end

  return $t{S,T,U,V}(copy(p.x0), x, Q, E)
end

# Create an undefined TaylorMap with properties from passed map
# except for parameters which are immutable
function $t(u::UndefInitializer, m::TaylorMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  x0 = Vector{S}(undef, nv)
  x = Vector{T}(undef, nn)

  # use same parameters if same descriptor (use=nothing)
  @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)

  if !isnothing(m.Q)
    q = Vector{T}(undef, 4)
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = Matrix{S}(undef, nv, nv)
  else
    E = nothing
  end

  return $t(x0, x, Q, E)
end

# Create from vector and blank (in this case must check for consistency)
# This constructor is type unstable if the spin or radiation kwargs are set, 
# and should be avoided in performance critical applications
function $t(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,<:TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}
  nv = numvars(x)
  np = numparams(x)
  nn = numnn(x)

  length(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")
  numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")


  x1 = Vector{T}(undef, nn)
  @inbounds x1[1:nv] = map(x->(T)(x, use=getdesc(use)), x)

  if T == TPS
    @inbounds x1[nv+1:nn] = params(getdesc(first(x)))
  else
    @inbounds x1[nv+1:nn] = complexparams(getdesc(first(x)))
  end

  if isnothing(spin)
    (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(Q.q[i],use=getdesc(use))
    end
    Q1 = Quaternion(q)
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(first(x1)) # implicilty uses use descriptor
    else
      (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      q = Vector{T}(undef, 4)
      for i=1:4
        @inbounds q[i] = T(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
  else
    # error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    Q1 = nothing # For type instability
  end

  if isnothing(radiation)
    E1 = E
  elseif radiation
    if isnothing(E)
      E1 = zeros(eltype(x0), nv, nv) 
    else
      (nv,nv) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      E1 = E
    end
  else
    # error("For no radiation, please omit the radiation kwarg or set radiation=nothing") # For type stability
    E1 = nothing # for type instability
  end

  return $t{S,T,typeof(Q1),typeof(E1)}(copy(x0), x1, Q1, E1)
end

# zero map (empty but still identity in parameters)
function zero(m::$t{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  
  x = Vector{T}(undef, nn)
  for i=1:nv
    @inbounds x[i] = T(use=desc)
  end

  # use same parameters 
  @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)

  if !isnothing(m.Q)
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = zeros(eltype(m.E), nv, nv)
  else
    E = nothing
  end

  return $t(zeros(eltype(m.x0), nv), x, Q, E)
end

# identity map
function one(m::$t{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  
  x = Vector{T}(undef, nn)
  if T == ComplexTPS
    for i=1:nv
      @inbounds x[i] = complexmono(i,use=desc)
    end
  else
    for i=1:nv
      @inbounds x[i] = mono(i,use=desc)
    end
  end

  # use same parameters 
  @inbounds x[nv+1:nn] = view(m.x, nv+1:nn)

  if !isnothing(m.Q)
    q = Vector{T}(undef, 4)
    @inbounds q[1] = one(first(x))
    for i=2:4
      @inbounds q[i] = T(use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if !isnothing(m.E)
    E = zeros(eltype(m.E), nv, nv)
  else
    E = nothing
  end

  return $t(zeros(eltype(m.x0), nv), x, Q, E)
end

end
end


==(m1::TaylorMap, m2::TaylorMap) = (m1.x0 == m2.x0 && m1.x == m2.x && m1.Q == m2.Q && m1.E == m2.E)
