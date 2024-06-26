for t = (:DAMap, :TPSAMap)
@eval begin

"""
    $($t)(m::Union{TaylorMap{S,T,U,V,W},Probe{S,T,U,V,W}}; use::UseType=m, idpt::Union{Nothing,Bool}=m.idpt) where {S,T,U,V,W}

Creates a new copy of the passed `TaylorMap` as a `$($t)`. 

If `use` is not specified, then the same GTPSA `Descriptor` as `m` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
copy of `m` as a new $($(t)) will have the same `Descriptor` as in `use.` The total number of variables + 
parameters must agree, however the orders may be different.

If `idpt` is not specified, then the same `idpt` as `m` will be used, else that specified will be used.
"""
function $t(m::Union{TaylorMap{S,T,U,V,W},Probe{S,T,U,V,W}}; use::UseType=m, idpt::Union{Nothing,Bool}=m.idpt) where {S,T,U,V,W}
  numnn(use) == numnn(m) || error("Number of variables + parameters in GTPSAs for `m` and `use` disagree!")
  
  outm = zero($t{S,T,U,V,typeof(idpt)}, use=use, idpt=idpt)
  nv = numvars(use)

  # set variables
  for i=1:nv
    @inbounds GTPSA.change!(outm.x[i], m.x[i])
  end

  # set quaternion
  if !isnothing(outm.Q)
    GTPSA.change!(outm.Q.q0, m.Q.q0)
    GTPSA.change!(outm.Q.q1, m.Q.q1)
    GTPSA.change!(outm.Q.q2, m.Q.q2)
    GTPSA.change!(outm.Q.q3, m.Q.q3)
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
    $($t)(;use::UseType=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 

Constructs a $($t) with the passed vector of `TPS`/`ComplexTPS` as the orbital ray, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and FD matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `FD` may be set to `true` to construct a $($t) with an identity quaternion/FD 
matrix, or `false` for no spin/FD. Note that setting `spin`/`FD` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`. The `use` kwarg may also be used to change the `Descriptor` of the TPSs, provided the number of variables 
+ parameters agree (orders may be different).
"""
function $t(;use::UseType=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 
  Base.require_one_based_indexing(x,x0)

  if !isnothing(Q)
    if !isnothing(E)
      T = Vector{promote_type(TPS,eltype(x0),eltype(x),eltype(Q),eltype(E))}
      Base.require_one_based_indexing(E)
    else
      T = Vector{promote_type(TPS,eltype(x0),eltype(x),eltype(Q))}
    end
  else
    T = Vector{promote_type(TPS,eltype(x0),eltype(x))}
  end

  S = Vector{numtype(eltype(T))}
  
  # set up
  if isnothing(spin)
    if isnothing(Q)
      U = Nothing
    else
      U = Quaternion{eltype(T)}
    end
  elseif spin
    U = Quaternion{eltype(T)}
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #U = Nothing # For type instability
  end

  if isnothing(FD)
    if isnothing(E)
      V = Nothing
    else
      V = Matrix{numtype(eltype(T))}
    end
  elseif FD
    V = Matrix{numtype(eltype(T))}
  else
    error("For no fluctuation-dissipation, please omit the FD kwarg or set FD=nothing") # For type stability
    #V = Nothing # For type instability
  end

  outm = zero($t{S,T,U,V,typeof(idpt)}, use = use, idpt = idpt)   

  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  # sanity checks
  length(x0) <= nv || error("Number of variables $nv != length of reference orbit vector $(length(x0))!")
  length(x) <= nv || error("Number of variables in GTPSAs for `x` and `use` disagree!")


  @views outm.x0 .= x0[1:length(outm.x0)]

  # set variables
  for i=1:nv
    GTPSA.change!(outm.x[i], x[i])
  end

  # set quaternion
  if !isnothing(outm.Q)
    if !isnothing(Q)
      GTPSA.change!(outm.Q.q0, m.Q.q0)
      GTPSA.change!(outm.Q.q1, m.Q.q1)
      GTPSA.change!(outm.Q.q2, m.Q.q2)
      GTPSA.change!(outm.Q.q3, m.Q.q3)
    else
      outm.Q.q0[0] = 1
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
    $($t)(M::AbstractMatrix; use::UseType=GTPSA.desc_current, x0::Vector{S}=zeros(numtype(T), numvars(use)), Q::U=nothing, E::V=nothing, idpt::W=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) where {S,U<:Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

This function could be optimized.

`M` must represent a matrix with linear indexing.

Constructs a $($t) with the passed matrix of scalars `M` as the linear part of the `TaylorMap`, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and FD matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `FD` may be set to `true` to construct a $($t) with an identity quaternion/FD 
matrix, or `false` for no spin/FD. Note that setting `spin`/`FD` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`.
"""
function $t(M::AbstractMatrix; use::UseType=GTPSA.desc_current, x0::Vector=zeros(eltype(M), size(M,1)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, idpt::Union{Bool,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing) 
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

function zero(::Type{$t{S,T,U,V,W}}; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  desc = getdesc(use)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)

  x0 = similar(S, nv) 
  x = similar(T, nn) 
  Base.require_one_based_indexing(x0, x)

  x0 .= 0

  for i=1:nv
    @inbounds x[i] = eltype(x)(use=desc)
  end

  # use same parameters if use isa TaylorMap and eltype(x) == eltype(use.x)
  if use isa TaylorMap && eltype(x) == eltype(use.x)
    @inbounds x[nv+1:nn] .= view(use.x, nv+1:nn)
  else # allocate
    if eltype(x) == TPS
      @inbounds x[nv+1:nn] .= params(desc)
    else
      @inbounds x[nv+1:nn] .= complexparams(desc)
    end
  end

  if U != Nothing
    q0 = eltype(x)(use=desc)
    q1 = eltype(x)(use=desc)
    q2 = eltype(x)(use=desc)
    q3 = eltype(x)(use=desc)
    Q = Quaternion(q0,q1,q2,q3)
  else
    Q = nothing
  end

  if V != Nothing
    E = similar(V, nv, nv)
    E .= 0
  else
    E = nothing
  end
  return $t(x0, x, Q, E, idpt)
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

function one(t::Type{$t{S,T,U,V,W}}; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  m = zero(t, use=use, idpt=idpt)
  nv = numvars(m)
  
  for i=1:nv
    @inbounds m.x[i][i] = 1
  end

  if !isnothing(m.Q)
    @inbounds m.Q.q0[0] = 1
  end

  return m
end

end
end
