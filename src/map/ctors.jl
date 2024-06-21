for t = (:DAMap, :TPSAMap)
@eval begin

"""
    $($t)(m::TaylorMap{S,T,U,V,W}; use::UseType=nothing, idpt::Union{Nothing,Bool}=m.idpt) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

Creates a new copy of the passed `TaylorMap` as a `$($t)`. 

If `use` is not specified, then the same GTPSA `Descriptor` as `m` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
copy of `m` as a new $($(t)) will have the same `Descriptor` as in `use.` The total number of variables + 
parameters must agree, however the orders may be different.

If `idpt` is not specified, then the same `idpt` as `m` will be used, else that specified will be used.
"""
function $t(m::TaylorMap{S,T,U,V,W}; use::UseType=nothing, idpt::Union{Nothing,Bool}=m.idpt) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}
  if isnothing(use)
    use = getdesc(m)
  else
    numnn(use) == numnn(m) || error("Number of variables + parameters in GTPSAs for `m` and `use` disagree!")
  end
  
  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)


  x = Vector{T}(undef, nn)
  for i=1:nv
    x[i] = T(m.x[i], use=getdesc(use))
  end

  # use same parameters if same descriptor (use=nothing)
  if isnothing(use) || getdesc(use) == getdesc(m)
    @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)
  else
    if T == TPS
      @inbounds x[nv+1:nn] .= params(getdesc(first(x)))
    else
      @inbounds x[nv+1:nn] .= complexparams(getdesc(first(x)))
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

  x0 = Vector{S}(undef, nv)
  if nv > numvars(m)  # Increasing dimensionality
    x0[1:numvars(m)] .= m.x0
    x0[numvars(m)+1:nv] .= 0
  else    # Reducing or keeping same dimensionality
    x0[1:nv] .= view(m.x0, 1:nv)
  end

  if !isnothing(m.E)
    E = similar(m.E, nv, nv)
    if nv > numvars(m)  # Increasing dimensionality
      E[1:numvars(m),1:numvars(m)] .= view(m.E, 1:numvars(m), 1:numvars(m))
      E[numvars(m)+1:nv,:] .= 0
      E[:,numvars(m)+1:nv] .= 0
    else
      E[1:nv,1:nv] .= view(m.E, 1:nv, 1:nv)
    end
  else
    E = nothing
  end

  return $t{S,T,U,V,typeof(idpt)}(x0, x, Q, E, idpt)
end

"""
    $($t)(p::Probe{S,T,U,V,W}; use::UseType=nothing, idpt::Union{Nothing,Bool}=p.idpt) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}

Creates a `$($t)` from the `Probe`, which must contain `TPS`s. 

If `use` is not specified, then the same GTPSA `Descriptor` as `p` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
`p` promoted to a $($(t)) will have the same `Descriptor` as in `use.` The total number of variables + 
and parameters must agree, however the orders may be different.

If `idpt` is not specified, then the same `idpt` as `m` will be used, else that specified will be used.
"""
function $t(p::Probe{S,T,U,V,W}; use::UseType=nothing, idpt::Union{Nothing,Bool}=p.idpt) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}
  if isnothing(use)
    use = getdesc(p)
  else
    numnn(use) == numnn(p) || error("Number of variables + parameters in GTPSAs for `p` and `use` disagree!")
  end
  
  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)
  
  x = Vector{T}(undef, nn)

  length(p.x) == numvars(p) || error("Length of orbital ray ($(length(p.x))) inconsistent with number of variables in GTPSA ($(nv))")
  
  for i=1:numvars(p)
    @inbounds x[i] = T(p.x[i], use=getdesc(use))
  end

  for i=numvars(p)+1:nv
    @inbounds x[i] = T(use=getdesc(use))
  end

  if T == TPS
    @inbounds x[nv+1:nn] .= params(getdesc(first(x)))
  else
    @inbounds x[nv+1:nn] .= complexparams(getdesc(first(x)))
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

  x0 = Vector{S}(undef, nv)
  if nv > numvars(p)  # Increasing dimensionality
    x0[1:numvars(p)] .= p.x0
    x0[numvars(p)+1:nv] .= 0
  else    # Reducing or keeping same dimensionality
    x0[1:nv] .= view(p.x0, 1:nv)
  end

  if !isnothing(p.E)
    E = similar(p.E, nv, nv)
    if nv > numvars(p)  # Increasing dimensionality
      E[1:numvars(p),1:numvars(p)] .= view(p.E, 1:numvars(p), 1:numvars(p))
      E[numvars(p)+1:nv,:] .= 0
      E[:,numvars(p)+1:nv] .= 0
    else
      E[1:nv,1:nv] .= view(p.E, 1:nv, 1:nv)
    end
  else
    E = nothing
  end

  return $t{S,T,U,V,typeof(idpt)}(x0, x, Q, E, idpt)
end

"""
    $($t){S,T,U,V,W}(u::UndefInitializer; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}

Creates an undefined `$($t){S,T,U,V}` with same `Descriptor` as `use`. The immutable 
parameters will be allocated if `use` is not a `TaylorMap`, else the immutable parameters 
from `use` will be used.
"""
function $t{S,T,U,V,W}(u::UndefInitializer; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
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
    q = Vector{T}(undef, 4)
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if V != Nothing
    E = Matrix{S}(undef, nv, nv)
  else
    E = nothing
  end

  return $t(x0, x, Q, E, idpt)
end

"""
    $($t)(; x::Union{Vector{<:Union{TPS,ComplexTPS}},Nothing}=nothing, x0::Union{Vector,Nothing}=nothing, Q::Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing}=nothing, E::Union{Matrix,Nothing}=nothing,  spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing, idpt::Union{Nothing,Bool}=nothing, use::UseType=nothing)

Constructs a $($t) with the passed vector of `TPS`/`ComplexTPS` as the orbital ray, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `stochastic` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/stochastic. Note that setting `spin`/`stochastic` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`. The `use` kwarg may also be used to change the `Descriptor` of the TPSs, provided the number of variables 
+ parameters agree (orders may be different).
"""
function $t(; x::Union{Vector{<:Union{TPS,ComplexTPS}},Nothing}=nothing, x0::Union{Vector,Nothing}=nothing, Q::Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing}=nothing, E::Union{Matrix,Nothing}=nothing,  spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing, idpt::Union{Nothing,Bool}=nothing, use::UseType=nothing)
  if isnothing(use)
    use = GTPSA.desc_current
  end

  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  if isnothing(x)
    x1 = Vector{TPS}(undef, nn)
    for i=1:nv
      @inbounds x1[i] = TPS(use=use)
    end
  else
    numvars(use) == numvars(x) || error("Number of variables in GTPSAs for `x` and `use` disagree!")
    numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
    x1 = Vector{eltype(x)}(undef, nn)
    @inbounds x1[1:nv] .= view(x, 1:nv)
  end

  if eltype(x1) == TPS
    @inbounds x1[nv+1:nn] .= params(getdesc(use))
  else
    @inbounds x1[nv+1:nn] .= complexparams(getdesc(use))
  end

  if isnothing(x0)
    x01 = zeros(numtype(first(x1)), nv)
  else
    numvars(x1) == length(x0) || error("Number of variables != length of reference orbit vector!")
    x01 = copy(x0)
  end

  if isnothing(spin)
    if isnothing(Q)
      Q1 = Q
    else
      (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      q = Vector{eltype(x1)}(undef, 4)
      for i=1:4
        @inbounds q[i] = (eltype(x1))(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(first(x1)) # implicilty uses use descriptor
    else
      (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      q = Vector{eltype(x1)}(undef, 4)
      for i=1:4
        @inbounds q[i] = (eltype(x1))(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #Q1 = nothing # For type instability
  end

  if isnothing(stochastic)
    E1 = E
  elseif stochastic
    if isnothing(E)
      E1 = zeros(eltype(x01), nv, nv) 
    else
      (nv,nv) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      E1 = E
    end
  else
    error("For no stochastic, please omit the stochastic kwarg or set stochastic=nothing") # For type stability
    #E1 = nothing # for type instability
  end

  return $t(x01, x1, Q1, E1, idpt)
end

"""
    $($t)(M; use::UseType=GTPSA.desc_current, x0::Vector{S}=zeros(eltype(M), size(M,1)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing, idpt::Union{Nothing,Bool}=nothing) where {S,U<:Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing},V<:Union{Matrix,Nothing}}

`M` must represent a matrix with linear indexing.

Constructs a $($t) with the passed matrix of scalars `M` as the linear part of the `TaylorMap`, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `stochastic` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/stochastic. Note that setting `spin`/`stochastic` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`.
"""
function $t(M; use::UseType=GTPSA.desc_current, x0::Vector{S}=zeros(eltype(M), numvars(use)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing, idpt::Union{Nothing,Bool}=nothing) where {S,U<:Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing},V<:Union{Matrix,Nothing}}
  Base.require_one_based_indexing(M)
  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  nv == length(x0) || error("Number of variables in GTPSA Descriptor != length of reference orbit vector!")
  nv >= size(M,1) || error("Number of rows in transfer matrix > number of variables in GTPSA!")

  if eltype(M) <: Complex
    outT = ComplexTPS
  else
    outT = TPS
  end

  x1 = Vector{outT}(undef, nn)
  for i=1:size(M,1)
    @inbounds x1[i] = (outT)(use=getdesc(use))
    for j=1:size(M,2)
      @inbounds x1[i][j] = M[i,j]
    end
  end

  for i=size(M,1)+1:nv
    @inbounds x1[i] = (outT)(use=getdesc(use))
  end

  if outT == TPS
    @inbounds x1[nv+1:nn] .= params(getdesc(first(x1)))
  else
    @inbounds x1[nv+1:nn] .= complexparams(getdesc(first(x1)))
  end

  if isnothing(spin)
    if isnothing(Q)
      Q1 = Q
    else
      q = Vector{outT}(undef, 4)
      for i=1:4
        @inbounds q[i] = outT(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(first(x1)) # implicilty uses use descriptor
    else
      q = Vector{outT}(undef, 4)
      for i=1:4
        @inbounds q[i] = outT(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #Q1 = nothing # For type instability
  end

  if isnothing(stochastic)
    E1 = E
  elseif stochastic
    if isnothing(E)
      E1 = zeros(eltype(x01), nv, nv) 
    else
      (nv,nv) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      E1 = E
    end
  else
    error("For no stochastic, please omit the stochastic kwarg or set stochastic=nothing") # For type stability
    #E1 = nothing # for type instability
  end
  return $t{eltype(M),outT,typeof(Q1),typeof(E1),typeof(idpt)}(copy(x0), x1, Q1, E1,idpt)
end


"""
    zero(m::$($t))

Creates a $($t) with the same GTPSA `Descriptor`, and spin/stochastic on/off,
as `m` but with all zeros for each quantity (except for the immutable parameters 
in `x[nv+1:nn]`, which will be copied from `m.x`)
"""
function zero(m::$t)
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)
  
  x = Vector{eltype(m.x)}(undef, nn)
  for i=1:nv
    @inbounds x[i] = (eltype(m.x))(use=desc)
  end

  # use same parameters 
  @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)

  if !isnothing(m.Q)
    q = Vector{eltype(m.x)}(undef, 4)
    for i=1:4
      @inbounds q[i] = (eltype(m.x))(use=desc)
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

  return $t(zeros(eltype(m.x0), nv), x, Q, E, m.idpt)
end

function zero(::Type{$t{S,T,U,V,W}}; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
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
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(use=desc)
    end
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if V != Nothing
    E = zeros(S, nv, nv)
  else
    E = nothing
  end

  return $t(x0, x, Q, E, idpt)
end

"""
    one(m::$($t))
  
Construct an identity map based on `m`.
"""
function one(m::$t)
  desc = getdesc(m)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  T = eltype(m.x)
  
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
  @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)

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

  return $t(zeros(eltype(m.x0), nv), x, Q, E, m.idpt)
end

function one(::Type{$t{S,T,U,V,W}}; use::UseType=GTPSA.desc_current, idpt::W=nothing) where {S,T,U,V,W}
  desc = getdesc(use)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  x0 = zeros(S, nv)

  x = Vector{T}(undef, nn)
  for i=1:nv
    @inbounds x[i] = T(use=desc)
    x[i][i] = 1
  end

  if use isa Union{TaylorMap,Probe} && eltype(use.x) == T
    # use same parameters 
    @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)
  else
    # allocate
    if T == TPS
      @inbounds x[nv+1:nn] .= params(desc)
    else
      @inbounds x[nv+1:nn] .= complexparams(desc)
    end
  end

  if U != Nothing
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(use=desc)
    end
    q[1][0] = 1
    Q = Quaternion(q)
  else
    Q = nothing
  end

  if V != Nothing
    E = zeros(S, nv, nv)
  else
    E = nothing
  end

  return $t(x0, x, Q, E, idpt)
end

end
end



