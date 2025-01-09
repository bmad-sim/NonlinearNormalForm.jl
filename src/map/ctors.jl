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
    q0 = eltype(U)(use=desc)
    q1 = eltype(U)(use=desc)
    q2 = eltype(U)(use=desc)
    q3 = eltype(U)(use=desc)
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

# Blank lowest-level ctor
function $t{S,T,U,V}(use::UseType=GTPSA.desc_current) where {S,T,U,V}
  getdesc(use).desc != C_NULL || error("GTPSA Descriptor not defined!")
  out_x0 = init_x0(S, use)
  out_x = init_x(T, use)
  out_Q = init_Q(U, use)
  out_E = init_E(V, use)

  return $t(out_x0, out_x, out_Q, out_E)
end

# Copy ctor including optional Descriptor change
function $t(m::TaylorMap{S,T,U,V}; use::UseType=m) where {S,T,U,V}
  numvars(use) == numvars(m) || error("Number of variables in GTPSAs for `m` and `use` disagree!")
  numparams(use) == numparams(m) || error("Number of parameters in GTPSAs for `m` and `use` disagree!") 

  out_m = $t{S,T,U,V}(use)
  out_m.x0 = m.x0
  foreach((out_xi, xi)->setTPS!(out_xi, xi, change=true), view(out_m.x, 1:numvars(use)), m.x)
  if !isnothing(out_m.Q)
    foreach((out_Qi, Qi)->setTPS!(out_Qi, Qi, change=true), out_m.Q, m.Q)
  end
  if !isnothing(out_m.E)
    out_m.E .= m.E
  end
  return out_m
end

function $t(;
  use::UseType=GTPSA.desc_current,
  x0::Union{AbstractVector,Nothing}=nothing,
  x::Union{AbstractVector,AbstractMatrix,UniformScaling,Nothing}=nothing,
  Q::Union{Quaternion,AbstractVector,AbstractMatrix,UniformScaling,Nothing}=nothing,
  E::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
  stochastic::Union{Bool,Nothing}=nothing,
) 
  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? GTPSA.numtype(eltype(t)) : Float64), (x0, x, Q, E))...)
  
  S = isnothing(x0) ? Vector{W} : promote_x0_type(typeof(x0), W)
  T = isnothing(x) ? Vector{TPS{W}} : promote_x_type(typeof(x), TPS{W})

  if isnothing(spin)
    U = isnothing(Q) ? Nothing : Quaternion{TPS{W}}
  elseif spin
    U = Quaternion{TPS{W}}
  else
    U = Nothing
  end

  if isnothing(stochastic)
    V = isnothing(E) ? Nothing : promote_E_type(typeof(E), W)
  elseif stochastic
    V = isnothing(E) ? Matrix{W} : promote_E_type(typeof(E), W)
  else
    V = Nothing
  end

  # Construct map using low level ctor:
  m = $t{S,T,U,V}(use)

  # Set if values provided:
  if !isnothing(x0)
    m.x0 .= x0
  end
  
  if !isnothing(x)
    if x isa AbstractVector #  TPSA map or scalar part provided:
      length(x) <= numvars(use) || error("Length of input vector `x` cannot be greater than the number of variables in `use` GTPSA!")
      foreach((out_xi, xi)->setTPS!(out_xi, xi, change=true), view(m.x, 1:length(x)), x)
    elseif x isa AbstractMatrix # Map as a matrix:
      size(x,1) <= numvars(use) || error("Number of rows of input matrix `x` cannot be greater than the number of variables in `use` GTPSA!")
      size(x,2) <= GTPSA.numcoefs(first(m.x))-1 || error("Number of columns of input matrix `x` cannot be greater than the number of coefficients in `use` GTPSA!")
      for varidx in 1:size(x,1)
        m.x[varidx][1:size(x,2)] = view(x, varidx, :)
      end
    else # Uniform scaling: Making identity map
      for varidx in 1:numvars(use)
        m.x[varidx][varidx] = 1
      end
    end
  end

  if !isnothing(m.Q) && !isnothing(Q)
    if Q isa AbstractVector || Q isa Quaternion #  TPSA map or scalar part provided:
      length(Q) <= 4 || error("Length of input vector `Q` cannot be greater than 4!")
      foreach((out_Qi, Qi)->setTPS!(out_Qi, Qi, change=true), m.Q, Q)
    elseif Q isa AbstractMatrix # Map as a matrix:
      size(Q,1) <= 4 || error("Number of rows of input matrix `Q` cannot be greater than 4!")
      size(Q,2) <= GTPSA.numcoefs(first(m.Q))-1 || error("Number of columns of input matrix `x` cannot be greater than the number of coefficients in `use` GTPSA!")
      for qidx in 1:size(Q,1)
        m.Q[qidx][1:size(Q,2)] = view(Q, qidx, :)
      end
    else # Uniform scaling: Making identity quaternion:
      m.Q.q0[0] = 1
    end
  end

  if !isnothing(m.E) && !isnothing(E)
    m.E .= E
  end

  return m
end


"""
    zero(m::$($t))

Creates a zero $($t) with the same properties as `m`, including GTPSA `Descriptor`,
spin, and stochasticity.
"""
zero(m::$t) = typeof(m)(use=m)

"""
    one(m::$($t))
  
Creates an identity $($t) with the same properties as `m`, including GTPSA
`Descriptor`, spin, and stochasticity.
"""
function one(m::$t)
  out_m = zero(m)

  for i in 1:nv
    out_m.x[i][i] = 1
  end

  if !isnothing(m.Q)
    out_m.Q.q0[0] = 1
  end

  return out_m
end


end
end
