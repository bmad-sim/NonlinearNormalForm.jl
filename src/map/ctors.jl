"""
    TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}

Abstract type for `TPSAMap` and `DAMap` used for normal form analysis. 

All `TaylorMap`s contain `x0` and `x` as the entrance coordinates and transfer map 
as a truncated power series respectively. If spin is included, a field `Q` containing 
a `Quaternion` as a truncated power series is included, else `Q` is `nothing`. If 
radiation is included, a field `E` contains a matrix of the envelope for stochastic 
radiation, else `E` is nothing.

### Fields
- `x0` -- Entrance coordinates of the map, Taylor expansion point
- `x`  -- Orbital ray as a truncated power series, expansion around `x0` + scalar part equal to EXIT coordinates of map
- `Q`  -- `Quaternion` as a truncated power series if spin is included, else `nothing`
- `E`  -- Matrix of the envelope for stochastic radiation if included, else `nothing`
"""
abstract type TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} end 

"""
    DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverses as a `DAMap` (with the scalar part ignored).
See `TaylorMap` for more information.
"""
struct DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

"""
    TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}

`TaylorMap` that composes and inverses as a `TPSAMap` (with the scalar part included).
See `TaylorMap` for more information.
"""
struct TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}  <: TaylorMap{S,T,U,V}
  x0::Vector{S}    # Entrance value of map
  x::Vector{T}     # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  Q::U             # Quaternion for spin
  E::V             # Envelope for stochastic radiation
end

for t = (:DAMap, :TPSAMap)
@eval begin

"""
    $($t)(TaylorMap{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

Creates a new copy of the passed `TaylorMap` as a `$($t)`. 

If `use` is not specified, then the same GTPSA `Descriptor` as `m` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
copy of `m` as a new $($(t)) will have the same `Descriptor` as in `use.` The number of variables 
and parameters must agree, however the orders may be different.
"""
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

  if !isnothing(m.E)
    E = copy(m.E)
  else
    E = nothing
  end

  return $t{S,T,U,V}(copy(m.x0), x, Q, E)
end

"""
    $($t)(p::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

Creates a `$($t)` from the `Probe`, which must contain `TPS`s. 

If `use` is not specified, then the same GTPSA `Descriptor` as `p` will be used. If `use` is 
specified (could be another `Descriptor`, `TaylorMap`, or a `Probe` containing `TPS`s), then the 
`p` promoted to a $($(t)) will have the same `Descriptor` as in `use.` The number of variables 
and parameters must agree, however the orders may be different.
"""
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

  if !isnothing(p.E)
    E = copy(p.E)
  else
    E = nothing
  end

  return $t{S,T,U,V}(copy(p.x0), x, Q, E)
end

"""
    $($t)(u::UndefInitializer, m::TaylorMap{S,T,U,V}) where {S,T,U,V}

Creates an undefined `$($t)` based on `m` (same `Descriptor` and with
radiation/spin on/off.)
"""
function $t(u::UndefInitializer, m::TaylorMap{S,T,U,V}) where {S,T,U,V}
  desc = getdesc(m)
  nv = numvars(desc)
  np = numparams(desc)
  nn = numnn(desc)
  x0 = Vector{S}(undef, nv)
  x = Vector{T}(undef, nn)

  # use same parameters if same descriptor (use=nothing)
  @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)

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

"""
    $($t)(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

Constructs a $($t) with the passed vector of `TPS`/`ComplexTPS` as the orbital ray, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `radiation` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/radiation. Note that setting `spin`/`radiation` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`. The `use` kwarg may also be used to change the `Descriptor` of the TPSs, provided the number of variables 
and parameters agree (orders may be different).
"""
function $t(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}
  nv = numvars(x)
  np = numparams(x)
  nn = numnn(x)

  length(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")
  numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")


  x1 = Vector{T}(undef, nn)
  @inbounds x1[1:nv] = map(x->(T)(x, use=getdesc(use)), x)

  if T == TPS
    @inbounds x1[nv+1:nn] .= params(getdesc(first(x)))
  else
    @inbounds x1[nv+1:nn] .= complexparams(getdesc(first(x)))
  end

  if isnothing(spin)
    if isnothing(Q)
      Q1 = Q
    else
      (!isnothing(use) || getdesc(x) == getdesc(Q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      q = Vector{T}(undef, 4)
      for i=1:4
        @inbounds q[i] = T(Q.q[i],use=getdesc(use))
      end
      Q1 = Quaternion(q)
    end
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

"""
    $($t)(M; x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}

`M` must represent a matrix with linear indexing.

Constructs a $($t) with the passed matrix of scalars `M` as the linear part of the `TaylorMap`, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `radiation` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/radiation. Note that setting `spin`/`radiation` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`.
"""
function $t(M; use::Union{Descriptor,TaylorMap,Probe{S,Union{TPS,ComplexTPS},U,V}}=GTPSA.desc_current, x0::Vector{S}=zeros(eltype(M), size(M,1)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing) where {S,U<:Union{Quaternion{<:Union{TPS,ComplexTPS}},Nothing},V<:Union{Matrix,Nothing}}
  Base.require_one_based_indexing(M)
  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  nv == length(x0) || error("Number of variables in GTPSA Descriptor != length of reference orbit vector!")
  nv == size(M,1) || error("Number of rows in transfer matrix inconsistent with number of variables in GTPSA!")

  if eltype(M) <: Complex
    outT = ComplexTPS
  else
    outT = TPS
  end

  x1 = Vector{outT}(undef, nn)
  for i=1:nv
    @inbounds x1[i] = (outT)(use=getdesc(use))
    for j=1:size(M,2)
      @inbounds x1[i][j] = M[i,j]
    end
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

  return $t{eltype(M),outT,typeof(Q1),typeof(E1)}(copy(x0), x1, Q1, E1)
end


"""
    zero(m::$($t){S,T,U,V}) where {S,T,U,V}

Creates a $($t) with the same GTPSA `Descriptor`, and spin/radiation on/off,
as `m` but with all zeros for each quantity (except for the immutable parameters 
in `x[nv+1:nn]`, which will be copied from `m.x`)
"""
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
  @inbounds x[nv+1:nn] .= view(m.x, nv+1:nn)

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

function zero(::Type{$t{S,T,U,V}}; use::Union{Descriptor,TaylorMap,Probe{<:Any,Union{TPS,ComplexTPS},<:Any,<:Any}}=GTPSA.desc_current) where {S,T,U,V}
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
      x[nv+1:nn] .= params(desc)
    else
      x[nv+1:nn] .= complexparams(desc)
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

  return $t(x0, x, Q, E)
end

"""
    one(m::$($t){S,T,U,V}) where {S,T,U,V}
  
Construct an identity map based on `m`.
"""
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

  return $t(zeros(eltype(m.x0), nv), x, Q, E)
end

function one(::Type{$t{S,T,U,V}}; use::Union{Descriptor,TaylorMap,Probe{S,Union{TPS,ComplexTPS},U,V}}=GTPSA.desc_current) where {S,T,U,V}
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
      x[nv+1:nn] .= params(desc)
    else
      x[nv+1:nn] .= complexparams(desc)
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

  return $t(x0, x, Q, E)
end

end
end



