for t = (:DAMap, :TPSAMap)
@eval begin

"""
    $($t)(m::Union{TaylorMap{S,T,U,V,W},Probe{S,T,U,V,W}}; use=nothing, idpt::Union{Nothing,Bool}=m.idpt) where {S,T,U,V,W}

Creates a new copy of the passed `TaylorMap` as a `$($t)`. 

If `use` is not specified, then the same GTPSA `Descriptor` as `m` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
copy of `m` as a new $($(t)) will have the same `Descriptor` as in `use.` The total number of variables + 
parameters must agree, however the orders may be different.

If `idpt` is not specified, then the same `idpt` as `m` will be used, else that specified will be used.
"""
function $t(m::Union{TaylorMap{S,T,U,V,W},Probe{S,T,U,V,W}}; use=m, idpt::Union{Nothing,Bool}=m.idpt) where {S,T<:Union{TPS,ComplexTPS},U,V,W}
  numnn(use) == numnn(m) || error("Number of variables + parameters in GTPSAs for `m` and `use` disagree!")
  
  outm = $t{S,T,U,V,typeof(idpt)}(undef, use = use, idpt = idpt) # undefined map but parameters allocated
  nv = numvars(use)

  # allocate variables
  for i=1:nv
    @inbounds outm.x[i] = T(m.x[i], use=getdesc(use))
  end

  # allocate quaternion if present
  if !isnothing(outm.Q)
    for i=1:4
      @inbounds outm.Q.q[i] = T(m.Q.q[i],use=getdesc(use))
    end
  end

  # set the reference orbit properly
  if nv > numvars(m)  # Increasing dimensionality
    outm.x0[1:numvars(m)] .= m.x0
    outm.x0[numvars(m)+1:nv] .= 0
  else    # Reducing or keeping same dimensionality
    outm.x0[1:nv] .= view(m.x0, 1:nv)
  end

  # set the FD matrix properly
  if !isnothing(outm.E)
    if nv > numvars(m)  # Increasing dimensionality
      outm.E[1:numvars(m),1:numvars(m)] .= view(m.E, 1:numvars(m), 1:numvars(m))
      outm.E[numvars(m)+1:nv,:] .= 0
      outm.E[:,numvars(m)+1:nv] .= 0
    else
      outm.E[1:nv,1:nv] .= view(m.E, 1:nv, 1:nv)
    end
  end

  return outm
end

"""
    $($t){S,T,U,V,W}(u::UndefInitializer; use=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}

Creates an undefined `$($t){S,T,U,V,W}` with same `Descriptor` as `use`. The immutable 
parameters will be allocated if `use` is not a `TaylorMap`, else the immutable parameters 
from `use` will be used.
"""
function $t{S,T,U,V,W}(u::UndefInitializer; use=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  desc = getdesc(use)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  x0 = Vector{S}(undef, nv)
  x = Vector{T}(undef, nn)

  

  # use same parameters if use isa TaylorMap and T == eltype(use.x)
  if use isa TaylorMap && T == eltype(use.x)
    @inbounds x[nv+1:nn] .= view(use.x, nv+1:nn)
  else # allocate
    if T == TPS
      @inbounds x[nv+1:nn] .= params(desc)
    else
      @inbounds x[nv+1:nn] .= complexparams(desc)
    end
  end

  if U != Nothing
    q = Vector{eltype(U)}(undef, 4)
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if V != Nothing
    E = Matrix{eltype(V)}(undef, nv, nv)
  else
    E = nothing
  end

  return $t{S,T,U,V,W}(x0, x, Q, E, idpt)
end

"""
    $($t)(;use=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 

Constructs a $($t) with the passed vector of `TPS`/`ComplexTPS` as the orbital ray, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and FD matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `FD` may be set to `true` to construct a $($t) with an identity quaternion/FD 
matrix, or `false` for no spin/FD. Note that setting `spin`/`FD` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`. The `use` kwarg may also be used to change the `Descriptor` of the TPSs, provided the number of variables 
+ parameters agree (orders may be different).
"""
function $t(;use=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 
  if !isnothing(Q)
    if !isnothing(E)
      T = promote_type(TPS,eltype(x0),eltype(x),eltype(Q),eltype(E))
    else
      T = promote_type(TPS,eltype(x0),eltype(x),eltype(Q))
    end
  else
    T = promote_type(TPS,eltype(x0),eltype(x))
  end

  S = numtype(T)
  
  # set up
  if isnothing(spin)
    if isnothing(Q)
      U = Nothing
    else
      U = Quaternion{T}
    end
  elseif spin
    U = Quaternion{T}
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #U = Nothing # For type instability
  end

  if isnothing(FD)
    if isnothing(E)
      V = Nothing
    else
      V = Matrix{S}
    end
  elseif FD
    V = Matrix{S}
  else
    error("For no fluctuation-dissipation, please omit the FD kwarg or set FD=nothing") # For type stability
    #V = Nothing # For type instability
  end

  outm = $t{S,T,U,V,typeof(idpt)}(undef, use = use, idpt = idpt)   

  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  # sanity checks
  length(x0) == nv || error("Number of variables != length of reference orbit vector!")
  length(x) == nv || error("Number of variables in GTPSAs for `x` and `use` disagree!")


  outm.x0 .= x0

  for i=1:nv
    @inbounds outm.x[i] = T(x[i], use=getdesc(use))
  end

  if !isnothing(outm.Q)
    if !isnothing(Q)
      for i=1:4
        @inbounds outm.Q.q[i] = T(Q.q[i], use=getdesc(use))
      end
    else
      for i=1:4
        @inbounds outm.Q.q[i] = T(use=getdesc(use))
      end
      outm.Q.q[1][0] = 1
    end
  end

  if !isnothing(outm.E)
    if !isnothing(E)
      (nv,nv) == size(E) || error("Size of FD matrix inconsistent with number of variables!")
      outm.E .= E
    else
      outm.E .= 0
    end
  end

  return outm
end

"""
    $($t)(M::AbstractMatrix; use=GTPSA.desc_current, x0::Vector{S}=zeros(numtype(T), numvars(use)), Q::U=nothing, E::V=nothing, idpt::W=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) where {S,U<:Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

`M` must represent a matrix with linear indexing.

Constructs a $($t) with the passed matrix of scalars `M` as the linear part of the `TaylorMap`, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and FD matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `FD` may be set to `true` to construct a $($t) with an identity quaternion/FD 
matrix, or `false` for no spin/FD. Note that setting `spin`/`FD` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`.
"""
function $t(M::AbstractMatrix; use=GTPSA.desc_current, x0::Vector=zeros(numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 
  Base.require_one_based_indexing(M)
  nv = numvars(use)
  nn = numnn(use)

  nv >= size(M,1) || error("Number of rows in matrix > number of variables in GTPSA!")

  if eltype(M) <: Complex
    T = ComplexTPS
  else
    T = TPS
  end

  x = Vector{T}(undef, nv)
  for i=1:size(M,1)
    @inbounds x[i] = T(use=getdesc(use))
    for j=1:size(M,2)
      @inbounds x[i][j] = M[i,j]
    end
  end

  return DAMap(use=use,x=x,x0=x0,Q=Q,E=E,idpt=idpt,spin=spin,FD=FD)
end


"""
    zero(m::$($t))

Creates a $($t) with the same GTPSA `Descriptor`, and spin/FD on/off,
as `m` but with all zeros for each quantity (except for the immutable parameters 
in `x[nv+1:nn]`, which will be copied from `m.x`)
"""
function zero(m::$t)
  return zero(typeof(m), use=m, idpt=m.idpt)
end

function zero(::Type{$t{S,T,U,V,W}}; use=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  desc = getdesc(use)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  x0 = zeros(S, nv)

  x = Vector{T}(undef, nn)
  for i=1:nv
    @inbounds x[i] = T(use=desc)
  end

  if use isa Union{TaylorMap,Probe} && eltype(use.x) == T
    # use same parameters 
    @inbounds x[nv+1:nn] .= view(use.x, nv+1:nn)
  else
    # allocate
    if T == TPS
      @inbounds x[nv+1:nn] .= params(desc)
    else
      @inbounds x[nv+1:nn] .= complexparams(desc)
    end
  end

  if U != Nothing
    q = Vector{eltype(U)}(undef, 4)
    for i=1:4
      @inbounds q[i] = eltype(U)(use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if V != Nothing
    E = Matrix{eltype(V)}(undef, nv, nv) #zeros(eltype(V), nv, nv)
  else
    E = nothing
  end

  return $t{S,T,U,V,W}(x0, x, Q, E, idpt)
end

function zero_op(m2::Union{$t,Number}, m1::Union{$t,Number})
  outtype = promote_type(typeof(m1),typeof(m2))

  # If either inputs are ComplexTPS, use those parameters
  if m2 isa $t
    if m1 isa $t
      if eltype(m2.x) == ComplexTPS
        return zero(outtype,use=m2,idpt=m2.idpt)
      else
        return zero(outtype,use=m1,idpt=m1.idpt)
      end
    else
      return zero(outtype,use=m2,idpt=m2.idpt)
    end
  elseif m1 isa $t
    return zero(outtype,use=m1,idpt=m1.idpt)
  else
    error("Cannot create new map based only on scalars")
  end
end


"""
    one(m::$($t))
  
Construct an identity map based on `m`.
"""
function one(m::$t)
  return one(typeof(m), use=m, idpt=m.idpt)
end

function one(t::Type{$t{S,T,U,V,W}}; use=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  m = zero(t, use=use, idpt=idpt)
  nv = numvars(m)
  
  for i=1:nv
    @inbounds m.x[i][i] = 1
  end

  if !isnothing(m.Q)
    @inbounds m.Q.q[1][0] = 1
  end

  return m
end

end
end
#=
function test1(m2,m1)
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

  # set up FD out
  if isnothing(m1.E) && isnothing(m2.E)
    outE = nothing
  else
    outE = Matrix{numtype(outT)}(undef, nv, nv)
  end
  return DAMap(outx0, outx, outQ, outE, m1.idpt)
end=#