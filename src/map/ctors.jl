#=

Constructors for the DAMap and TPSAMap. Includes zero, one, 
copy ctor (with possible descriptor change), custom map ctor using 
keyword arguments, map ctor from a matrix, and zero_op which creates 
a zero map of the properly promoted type (useful for out-of-place 
functions).

=#

function init_x0(::Type{S}, use) where {S}
  nv = numvars(use)
  x0 = similar(S, nv)
  x0 .= 0
  return x0
end

function init_x(::Type{T}, use) where {T}
  desc = getdesc(use)
  nv = numvars(desc)
  nn = numnn(desc)
  x = similar(T, nn)
  for i in 1:nv
    x[i] = eltype(x)(use=desc)
  end
  # use same parameters if use isa TaylorMap and eltype(x) == eltype(use.x)
  if use isa TaylorMap && eltype(x) == eltype(use.x)
    x[nv+1:nn] .= view(use.x, nv+1:nn)
  else # allocate
    if eltype(x) == TPS{Float64}
      x[nv+1:nn] .= params(desc)
    else
      x[nv+1:nn] .= complexparams(desc)
    end
  end
  return x
end

function init_Q(::Type{U}, use) where {U}
  if U != Nothing
    desc = getdesc(use)
    q0 = eltype(x)(use=desc)
    q1 = eltype(x)(use=desc)
    q2 = eltype(x)(use=desc)
    q3 = eltype(x)(use=desc)
    Q = Quaternion(q0,q1,q2,q3)
  else
    Q = nothing
  end
  return Q
end

function init_E(::Type{V}, use) where {V}
  if V != Nothing
    nv = numvars(use)
    E = similar(V, nv, nv)
    E .= 0
  else
    E = nothing
  end
  return E
end


for t = (:DAMap, :TPSAMap)
@eval begin

"""
    zero(m::$($t))

Creates a zero $($t) with the same properties as `m`, including GTPSA `Descriptor`,
spin, and stochasticity
"""
function zero(m::$t)
  return zero(typeof(m), use=m)
end

function zero(::Type{$t{S,T,U,V}}; use::UseType=GTPSA.desc_current) where {S,T,U,V}
  getdesc(use).desc != C_NULL || error("GTPSA Descriptor not defined!")
  x0 = init_x0(S, use)
  x = init_x(T, use)
  Q = init_Q(U, use)
  E = init_E(V, use)
  return $t(x0, x, Q, E)
end

"""
    one(m::$($t))
  
Construct an identity map based on `m`.
"""
function one(m::$t)
  return one(typeof(m), use=m)
end

function one(t::Type{$t{S,T,U,V}}; use::UseType=GTPSA.desc_current) where {S,T,U,V}
  m = zero(t, use=use)
  nv = numvars(m)
  
  for i in 1:nv
    m.x[i][i] = 1
  end

  if !isnothing(m.Q)
    m.Q.q0[0] = 1
  end

  return m
end

"""
    $($t)(m::Union{TaylorMap{S,T,U,V}; use::UseType=m) where {S,T,U,V}

Creates a new copy of the passed `TaylorMap` as a `$($t)`. 

If `use` is not specified, then the same GTPSA `Descriptor` as `m` will be used. If `use` is 
specified (could be another `Descriptor` or `TaylorMap`), then the copy of `m` as a new $($(t)) 
will have the same `Descriptor` as in `use.` The total number of variables + parameters must 
agree, however the orders may be different.
"""
function $t(m::Union{TaylorMap{S,T,U,V}}; use::UseType=m) where {S,T,U,V}
  numnn(use) == numnn(m) || error("Number of variables + parameters in GTPSAs for `m` and `use` disagree!")
  
  outm = zero($t{S,T,U,V}, use=use)
  nv = numvars(use)

  # set variables
  for i=1:nv
    setTPS!(outm.x[i], m.x[i], change=true)
  end

  # set quaternion
  if !isnothing(outm.Q)
    setTPS!(outm.Q.q0, m.Q.q0, change=true)
    setTPS!(outm.Q.q1, m.Q.q1, change=true)
    setTPS!(outm.Q.q2, m.Q.q2, change=true)
    setTPS!(outm.Q.q3, m.Q.q3, change=true)
  end
  
  # set the reference orbit properly
  if nv > numvars(m)  # Increasing dimensionality
    outm.x0[1:numvars(m)] .= m.x0
    outm.x0[numvars(m)+1:nv] .= 0
  else    # Reducing or keeping same dimensionality
    outm.x0[1:nv] .= view(m.x0, 1:nv)
  end

  # set the stochastic matrix properly
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
    $($t)(;use::UseType=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(GTPSA.numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing) 

Constructs a $($t) with the passed vector of `TPS`/`ComplexTPS64` as the orbital ray, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `stochastic` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/stochastic. Note that setting `spin`/`stochastic` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`. The `use` kwarg may also be used to change the `Descriptor` of the TPSs, provided the number of variables 
+ parameters agree (orders may be different).
"""
function $t(;use::UseType=GTPSA.desc_current, x::Vector=vars(getdesc(use)), x0::Vector=zeros(GTPSA.numtype(eltype(x)), numvars(use)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing) 
  Base.require_one_based_indexing(x,x0)

  if !isnothing(Q)
    if !isnothing(E)
      T = Vector{promote_type(TPS{Float64},eltype(x0),eltype(x),eltype(Q),eltype(E))}
      Base.require_one_based_indexing(E)
    else
      T = Vector{promote_type(TPS{Float64},eltype(x0),eltype(x),eltype(Q))}
    end
  else
    T = Vector{promote_type(TPS{Float64},eltype(x0),eltype(x))}
  end

  S = Vector{GTPSA.numtype(eltype(T))}
  
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

  if isnothing(stochastic)
    if isnothing(E)
      V = Nothing
    else
      V = Matrix{GTPSA.numtype(eltype(T))}
    end
  elseif stochastic
    V = Matrix{GTPSA.numtype(eltype(T))}
  else
    error("For no fluctuation-dissipation, please omit the stochastic kwarg or set stochastic=nothing") # For type stability
    #V = Nothing # For type instability
  end

  outm = zero($t{S,T,U,V}, use=use)   

  nv = numvars(use)
  np = numparams(use)
  nn = numnn(use)

  # sanity checks
  length(x0) <= nv || error("Number of variables $nv != length of reference orbit vector $(length(x0))!")
  length(x) <= nv || error("Number of variables in GTPSAs for `x` and `use` disagree!")


  @views outm.x0 .= x0[1:length(outm.x0)]

  # set variables
  for i=1:nv
    setTPS!(outm.x[i], x[i], change=true)
  end

  # set quaternion
  if !isnothing(outm.Q)
    if !isnothing(Q)
      setTPS!(outm.Q.q0, Q.q0, change=true)
      setTPS!(outm.Q.q1, Q.q1, change=true)
      setTPS!(outm.Q.q2, Q.q2, change=true)
      setTPS!(outm.Q.q3, Q.q3, change=true)
    else
      outm.Q.q0[0] = 1
    end
  end

  if !isnothing(outm.E) && !isnothing(E)
    (nv,nv) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
    outm.E .= E
  end

  return outm
end

"""
    $($t)(M::AbstractMatrix; use::UseType=GTPSA.desc_current, x0::Vector=zeros(eltype(M), size(M,1)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing) 

`M` must be a matrix with linear indexing.

Constructs a $($t) with the passed matrix of scalars `M` as the linear part of the `TaylorMap`, and optionally the entrance 
coordinates `x0`, `Quaternion` for spin `Q`, and stochastic matrix `E` as keyword arguments. The helper keyword 
arguments `spin` and `stochastic` may be set to `true` to construct a $($t) with an identity quaternion/stochastic 
matrix, or `false` for no spin/stochastic. Note that setting `spin`/`stochastic` to any `Bool` value without `Q` or `E` 
specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA 
`Descriptor`.
"""
function $t(M::AbstractMatrix; use::UseType=GTPSA.desc_current, x0::Vector=zeros(eltype(M), size(M,1)), Q::Union{Quaternion,Nothing}=nothing, E::Union{Matrix,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing, stochastic::Union{Bool,Nothing}=nothing) 
  Base.require_one_based_indexing(M)
  m = DAMap(use=use, x0=x0, Q=Q, E=E, spin=spin, stochastic=stochastic)
  setmatrix!(m, M)
  return m
end


"""

    zero_op(m2::$($t), m1::Union{$($t),Number})

Returns a new zero $($t) with type equal to promoted type of `m1` and `m2`.
`m2` could be a number. This function is necessary to both ensure correct 
usage of the GTPSA descriptor and share the immutable parameters when possible.
"""
function zero_op(m2::$t, m1::Union{$t,Number})
  outtype = promote_type(typeof(m1),typeof(m2))

  # If either inputs are ComplexTPS64, use those parameters
  if m1 isa $t
    if eltype(m2.x) == ComplexTPS64
      return zero(outtype,use=m2)
    else
      return zero(outtype,use=m1)
    end
  else
    return zero(outtype,use=m2)
  end
end

end
end
