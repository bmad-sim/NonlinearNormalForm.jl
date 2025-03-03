#=

Defines the VectorField type, promotion rules, and constructors.

=#

# =================================================================================== #
# Type

"""
    VectorField{X,Q}

Lie operator which acts on `DAMap`s e.g. dM/dt = FM where F is the `VectorField` and 
`M` is the `DAMap`. `F` can be promoted to a `DAMap` using `exp(F)`.

# Fields
- `x::X` -- Orbital ray as truncated power series, expansion with scalar part equal to EXIT coordinates of map
- `q::Q` -- `Quaternion` as truncated power series if spin is included, else `nothing`

# Type Requirements
- `X <: AbstractVector{<:TPS}`
- `Q <: Union{Quaternion{<:TPS},Nothing}` and if `Q != Nothing` then `eltype(Q) == eltype(X)`
"""
struct VectorField{X<:AbstractVector,Q<:Union{Quaternion,Nothing},}
  x::X
  q::Q     
  function VectorField(x, q)
    F = new{typeof(x),typeof(q)}(x, q)
    checkvfsanity(F)
    return F
  end
end

@inline function checkvfsanity(F::VectorField{X,Q}) where {X,Q}
  # Static checks:
  TI.is_tps_type(eltype(X)) isa TI.IsTPSType || error("Orbital ray element type must be a truncated power series type supported by `TPSAInterface.jl`")
  Q == Nothing || eltype(Q) == eltype(X) || error("Quaternion number type $(eltype(Q)) must be $(eltype(X)) (equal to orbital ray)")

  # Runtime checks:
  length(F.x) <= ndiffs(getinit(first(F.x))) || error("Orbital ray number of variables ($(length(F.x))) exceeds number of differentials in TPSA ($(ndiffs(getinit(F.x))))")
  Q == Nothing || getinit(first(F.q)) == getinit(first(F.x)) || error("Quaternion TPSA definition disagrees with orbital ray TPSA definition")
end

# =================================================================================== #
# Field initialization functions.

# These may be overrided by external array packages.
function init_vf_x(::Type{X}, init::AbstractTPSAInit, nv::Integer) where {X<:AbstractVector}
  x = similar(X, nv)
  for i in 1:nv
    x[i] = TI.init_tps(TI.numtype(eltype(X)), init) 
  end
  return x
end

init_vf_q(::Type{Q}, init::AbstractTPSAInit) where {Q<:Union{Nothing,Quaternion}} = init_map_q(Q, init)

# =================================================================================== #
# StaticArray field initialization functions.

# Consistency checks are made by the `checkvfsanity` run by every map construction,
# so lengths of arrays here are not checked for consistency with the TPSA
function init_vf_x(::Type{X}, init::AbstractTPSAInit, nv::Integer) where {X<:StaticVector}
  x = StaticArrays.sacollect(X, TI.init_tps(TI.numtype(eltype(X)), init) for i in 1:length(X))
  return x
end

# =================================================================================== #
# Promotion rules

function promote_rule(::Type{VectorField{X,Q}}, ::Type{G}) where {X,Q,G<:Union{Number,Complex}}
  out_X = similar_eltype(X, promote_type(eltype(X), G))
  out_Q = Q == Nothing ? Nothing : similar_eltype(Q, promote_type(eltype(Q), G))
  return VectorField{out_X,out_Q}
end

function promote_rule(::Type{VectorField{X1,Q1}}, ::Type{VectorField{X2,Q2}}) where {X1,X2,Q1,Q2}
  out_X = similar_eltype(X1, promote_type(eltype(X1), eltype(X2)))
  !xor(Q1==Nothing, Q2==Nothing) || error("Cannot promote VectorFields: one includes spin while the other does not")
  out_Q = Q1 == Nothing ? Nothing : similar_eltype(Q1, promote_type(eltype(Q1),eltype(Q2)))
  return VectorField{out_X,out_Q}
end

# --- complex type ---
function complex(type::Type{VectorField{X,Q}}) where {X,Q}
  return promote_type(type, complex(TI.numtype(eltype(X))))
end

# --- real type ---
function real(::Type{VectorField{X,Q}}) where {X,Q}
  out_X = similar_eltype(X, real(eltype(X)))
  out_Q = Q == Nothing ? Nothing : similar_eltype(Q, real(eltype(Q)))
  return VectorField{out_X,out_Q}
end

# =================================================================================== #
# Constructors

# Lowest-level, internal
function _zero(::Type{VectorField{X,Q}}, init::AbstractTPSAInit, nv::Integer) where {X,Q}
  x = init_vf_x(X, init, nv)
  q = init_vf_q(Q, init)
  return VectorField(x, q)
end

#=
function zero(::Type{VectorField{X,Q}}) where {X,Q}
  return _zero(VectorField{X,Q}, getinit(eltype(X)))
end
=#
# Explicit type specification
# Def change would be static (in type)
function VectorField{X,Q}(F::VectorField) where {X,Q}
  out_F = zero(VectorField{X,Q}, F)
  copy!(out_F, F)
  return out_F
end

# zero but new type potentially
function zero(::Type{VectorField{X,Q}}, F::Union{VectorField,TaylorMap}) where {X,Q}
  return _zero(VectorField{X,Q}, getinit(eltype(X)), nvars(F))
end

zero(::Type{VectorField}, F::VectorField{X,Q}) where {X,Q} = zero(VectorField{X,Q}, F)

zero(::Type{VectorField}, m::TaylorMap{X0,X,Q,S}) where {X0,X,Q,S} = zero(VectorField{similar_eltype(X0,eltype(X)),Q}, m)

# Copy ctor including optional TPSA init change
function VectorField(F::VectorField; init::AbstractTPSAInit=getinit(F))
  X = similar_eltype(typeof(F.x), TI.init_tps_type(eltype(X0), init))
  Q = isnothing(F.q) ? Nothing : Quaternion{TI.init_tps_type(eltype(X0), init)}
  out_F = _zero(VectorField{X,Q}, init, nvars(F))
  copy!(out_F, F)
  return out_F
end

# Kwarg ctor
function VectorField(;
  init::Union{AbstractTPSAInit,Nothing}=nothing,
  nv::Integer=6,
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
) 
  if isnothing(init)
    if !isnothing(x) && TI.is_tps_type(eltype(x)) isa TI.IsTPSType
      init = getinit(first(x))
    elseif !isnothing(q) && TI.is_tps_type(eltype(q)) isa TI.IsTPSType
      init = getinit(first(q))
    else
      error("No TPSA definition has been provided, nor is one inferrable from the input arguments")
    end
  end

  nv <= ndiffs(init) || error("Number of variables specified for VectorField ($nv) is greater than the number of differentials ($(ndiffs(init))) in the TPSA")

  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? TI.numtype(eltype(t)) : Float64), (x, q))...)
  TW = TI.init_tps_type(W, init)
  X = isnothing(x) ? @_DEFAULT_X(nv){TW} : similar_eltype(typeof(x), TW)
  
  if isnothing(spin)
    Q = isnothing(q) && isnothing(q_map) ? Nothing : Quaternion{TW}
  elseif spin
    Q = Quaternion{TW}
  else
    Q = Nothing
  end

  # Construct vector field:
  out_F = _zero(VectorField{X,Q}, init, nv)

  setray!(out_F.x, x=x, x_matrix=x_matrix)
  if !isnothing(out_F.q) && !isnothing(q)
    setquat!(out_F.q, q=q, q_map=q_map)
  end

  return out_F
end

function zero(F::VectorField)
  return _zero(typeof(F), getinit(F), nvars(F))
end

# =================================================================================== #
# Lie bracket

"""
    lb!(G::VectorField, F::VectorField, H::VectorField; work_q::Union{Nothing,Quaternion}=nothing) -> G

Sets G equal to the Lie bracket of the vector fields `F` and `H`. Explicitly, including 
spin (with the lower case letter for the quaternion of the vector field), this is:

`(G,g) = ⟨(F,f) , (H,h)⟩ = (F⋅∇H-H⋅∇F , [h,f]+F⋅∇h-G⋅∇f)`

where `[h,f] = h*f - f*h` is just a quaternion commutator. See Equation 44.52 in the Bmad manual
"""
function lb!(
  G::VectorField, 
  F::VectorField, 
  H::VectorField; 
  work_q::Union{Nothing,Quaternion}= isnothing(G.q) ? nothing : zero(G.q)
)
  checkinplace(G, F, H)
  @assert !(G === F) && !(G === H) "Aliasing any source arguments with the destination in lb! is not allowed"

  # Orbital part (Eq. 44.51 in Bmad manual, handled by GTPSA):
  TI.liebra!(G.x, F.x, H.x)

  # Spin (Eq. 44.51 in Bmad manual, NOT handled by GTPSA):
  if !isnothing(F.q)
    # first [h,f] (this part involves quaternion multiplication so will be slow)
    mul!(G.q, H.q, F.q)
    mul!(work_q, F.q, H.q)

    TI.sub!(G.q.q0, G.q.q0, work_q.q0)
    TI.sub!(G.q.q1, G.q.q1, work_q.q1)
    TI.sub!(G.q.q2, G.q.q2, work_q.q2)
    TI.sub!(G.q.q3, G.q.q3, work_q.q3)

    # then +F⋅∇h 
    tmp = work_q.q0
    TI.fgrad!(tmp, F.x, H.q.q0)
    TI.add!(G.q.q0, G.q.q0, tmp)

    TI.fgrad!(tmp, F.x, H.q.q1)
    TI.add!(G.q.q1, G.q.q1, tmp)

    TI.fgrad!(tmp, F.x, H.q.q2)
    TI.add!(G.q.q2, G.q.q2, tmp)

    TI.fgrad!(tmp, F.x, H.q.q3)
    TI.add!(G.q.q3, G.q.q3, tmp)
    
    # finally -G⋅∇f
    TI.fgrad!(tmp, G.x, F.q.q0)
    TI.sub!(G.q.q0, G.q.q0, tmp)

    TI.fgrad!(tmp, G.x, F.q.q1)
    TI.sub!(G.q.q1, G.q.q1, tmp)

    TI.fgrad!(tmp, G.x, F.q.q2)
    TI.sub!(G.q.q2, G.q.q2, tmp)

    TI.fgrad!(tmp, G.x, F.q.q3)
    TI.sub!(G.q.q3, G.q.q3, tmp)
  end
  return G
end


# =================================================================================== #
# "Multiplication" (Lie operator acting on a DAMap)
"""
    mul!(m::DAMap, F::VectorField, m1::DAMap; work_Q::Union{Quaternion,Nothing}=isnothing(m.q) ? nothing : zero(m.q)) -> m

Computes the Lie operator `F` acting on a `DAMap` `m1`, and stores the result in `m`.
Explicity, that is `F * m = (F.x, F.q) * (m.x, m.q) = (F.x ⋅ ∇ m.x , F.x ⋅ ∇ m.q + m.q*F.q)`
"""
function mul!(
  m::DAMap, 
  F::VectorField,
  m1::DAMap; 
  work_q::Union{Quaternion,Nothing}= isnothing(m.q) ? nothing : zero(m.q)
)
  checkinplace(m, F, m1)

  nv = nvars(F)

  @assert !(m === m1) "Aliasing m === m1 is not allowed for mul!(m, F, m1)"

  m.x0 .= m1.x0
  # Orbital part F⋅∇m1 :
  for i=1:nv
    TI.fgrad!(m.x[i], F.x, m1.x[i])
  end

  # Spin F⋅∇q + qf (quaternion part acts in reverse order):
  if !isnothing(F.q)
    TI.mul!(work_q, m1.q, F.q)
    TI.fgrad!(m.q.q0, F.x, m1.q.q0)
    TI.fgrad!(m.q.q1, F.x, m1.q.q1)
    TI.fgrad!(m.q.q2, F.x, m1.q.q2)
    TI.fgrad!(m.q.q3, F.x, m1.q.q3)

    TI.add!(m.q.q0, m.q.q0, work_Q.q0)
    TI.add!(m.q.q1, m.q.q1, work_Q.q1)
    TI.add!(m.q.q2, m.q.q2, work_Q.q2)
    TI.add!(m.q.q3, m.q.q3, work_Q.q3)
  end

  if !isnothing(m.s)
    m.s .= m1.s
  end

  return m
end

# =================================================================================== #
# Exponentiation
"""
    exp!(
      m::DAMap, 
      F::VectorField, 
      m1::DAMap; 
      work_maps::NTuple{2, DAMap}=ntuple(t->zero(m1), Val{2}()), 
      work_q::Union{Quaternion,Nothing}=isnothing(m.q) ? nothing : zero(m.q)
    )

Computes `exp(F)*m1`, and stores the result in `m`. Explicity, this is
`exp(F)*m1 = m1 + F*m1 + 1/2*F*(F*m1) + 1/6*F*(F*(F*m1)) + ...`, where `*` is
the operation of a `VectorField` on the map. See the documentation for `mul!` for 
more details of this operation. If the linear part of `F` is zero, then the 
number of iterations is equal to the number of higher orders left.
"""
function exp!(
  m::DAMap, 
  F::VectorField, 
  m1::DAMap; 
  work_maps::NTuple{2, DAMap}=ntuple(t->zero(m1), Val{2}()), 
  work_q::Union{Quaternion,Nothing}=isnothing(m.q) ? nothing : zero(m.q)
)
  checkinplace(m, F, m1, work_maps...)
  nv = nvars(F)
  tmp = work_maps[1]
  tmp2 = work_maps[2]

  nmax = 100
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*nv
  nrm_ =Inf
  conv = false
  slow = false

  copy!(tmp, m1)
  copy!(m, m1)

  for j=1:nmax
    if j == 25
      slow=true
    end

    div!(tmp, tmp, j)
    mul!(tmp2, F, tmp, work_q=work_q)
    add!(m, m, tmp2)
    copy!(tmp, tmp2)
    nrm = norm(tmp2)

    # Check convergence
    if nrm <= nrm_min2 || conv && nrm >= nrm_ # done
      if slow
        @warn "exp! slow convergence: required n = $(j) iterations"
      end
      return m
    end

    if nrm <= nrm_min1 # convergence reached just refine a bit
      conv = true
    end
    nrm_ = nrm
  end
  @warn "exp! convergence not reached for $nmax iterations"
  return m
end

# =================================================================================== #
# Logarithm

"""
    log!(
      F::VectorField, 
      m1::DAMap; 
      work_maps::NTuple{3, DAMap}=ntuple(t->zero(m1), Val{3}()),
      work_vfs::NTuple{2, VectorField}=ntuple(t->zero(F), Val{2}()),
      work_q::Union{Quaternion,Nothing}=isnothing(m1.q) ? nothing : zero(m1.q)
    )

Computes the log of the map `m1` - that is, calculates the `VectorField` `F` that 
would represent the map `m1` as a Lie exponent `exp(F)` - and stores the result in `F`.
The map `m1` should be close to the identity for this to converge quickly.
"""
function log!(
  F::VectorField, 
  m1::DAMap; 
  work_maps::NTuple{3, DAMap}=ntuple(t->zero(m1), Val{3}()),
  work_vfs::NTuple{2, VectorField}=ntuple(t->zero(F), Val{2}()),
  work_q::Union{Quaternion,Nothing}=isnothing(m1.q) ? nothing : zero(m1.q)
)
  checkinplace(F, m1, work_maps..., work_vfs...)
  
  nv = nvars(m1)
  nmax = 100
  nrm0 = norm(m1)
  nrm_min1 = 1e-9
  nrm_min2 = 100*eps(Float64)*nv
  nrm0 = norm(m1)
  epsone = nrm0/1000
  conv = false 
  slow = false
  nrm_ = Inf

  # Choose the initial guess for the VectorField to be (M,q) - (I,1)
  sub!(F, m1, I)

  # Now we will iterate:
  for j=1:nmax
    if j == 25
      slow=true
    end
    # First, rotate back our guess:
    mul!(F, -1, F)
    exp!(work_maps[3], F, m1, work_maps=work_maps[1:2], work_q=work_q)
    mul!(F, -1, F)

    # Now we have (I,1) + (t,u) where (t,u) is the leftover garbage we want to be 0
    # solve for the vector field that gives us our garbage map
    # this will be exp(T)(I,1) = (I,1) + (t,u)
    # to second order in (t,u) is 
    # G = (t,u) + epsilon_2 where 
    # epsilon_2 = (t dot grad t, t dot grad u + u^2)
    # get the garbage as a VectorField and as a map
    sub!(work_maps[3], work_maps[3], I)
    copy!(work_vfs[1], work_maps[3])

    # Let the vector field act on the garbage map
    mul!(work_maps[2], work_vfs[1], work_maps[3], work_q=work_q)

    # Finally get our approximation 
    add!(work_vfs[1], work_vfs[1], work_maps[2])

    # Now we can see that  M = exp(F)exp(G) approximately
    # CBH formula would be exact combination for exp(F)exp(G)=exp(H)
    # but if the leftover stuff is still big then we can just approximately
    # with linear term (adding F+G trivially)

    if norm(work_vfs[1]) < epsone # small! use CBH!
      lb!(work_vfs[2], F, work_vfs[1], work_q=work_q)
      mul!(work_vfs[2], 0.5, work_vfs[2])
      add!(work_vfs[1], work_vfs[1], work_vfs[2])
      add!(F,F,work_vfs[1])
    else # big! just linear!
      add!(F,F,work_vfs[1])
    end

    nrm = norm(work_vfs[1])/nrm0

    if nrm <= nrm_min2 || conv && nrm >= nrm_
      if slow
        @warn "log! slow convergence: required n = $(j) iterations"
      end
      return
    end

    if nrm <= nrm_min1
      conv = true
    end

    nrm_ = nrm
  end
  @warn "log! convergence not reached for $(nmax) iterations"
  return
end

# =================================================================================== #
# Out-of-place operators

function lb(F::VectorField, H::VectorField)
  Fprom, Gprom = promote(F,G)
  H = zero(Fprom)
  lb!(H, Fprom, Gprom)
  return H
end


function *(F::VectorField, m1::Union{DAMap,UniformScaling})
  zer = zero(TI.numtype(eltype(m1.x)))
  Fprom = promote(F,zer)
  m1prom = promote(m1,zer)
  m = zero(m1prom)
  mul!(m, Fprom, m1prom)
  return m
end

function exp(F::VectorField, m1::Union{UniformScaling,DAMap}=I)
  if m1 isa UniformScaling
    error("Not currently supported please provide a map")
    # TO-DO: How/should I store nparams in VectorField?
    # It is just equal to nvars - ndiffs, though this will be slow...
    m1prom = one(DAMap{typeof(F.x),typeof(F.q)})
    Fprom = F
  else
    if TI.numtype(eltype(m1.x)) != TI.numtype(eltype(F.x))
      if promote_type(typeof(F), TI.numtype(eltype(m1.x))) != typeof(F)
        Fprom = zero(promote_type(typeof(F), TI.numtype(eltype(m1.x))), F)
        copy!(Fprom, F)
      else
        m1prom = zero(promote_type(typeof(m1), TI.numtype(eltype(m1.x))), m1)
        copy!(m1prom, m1)
      end
    else
      Fprom = F
      m1prom = m1
    end
  end
  m = zero(m1prom)
  exp!(m, Fprom, m1prom)
  return m
end

function log(m1::DAMap)
  F = zero(VectorField, m1)
  log!(F, m1)
  return F
end