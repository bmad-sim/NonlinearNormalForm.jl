# --- getvectorfield ---
vec2fld!(na::Cint, tpsa::TPS{Float64}, m::Vector{TPS{Float64}}) = (@inline; GTPSA.mad_tpsa_vec2fld!(na, tpsa, m))
vec2fld!(na::Cint, ctpsa::TPS{ComplexF64}, m::Vector{TPS{ComplexF64}}) = (@inline; GTPSA.mad_ctpsa_vec2fld!(na, ctpsa, m))

"""
VectorField(h::TPS; Q::Union{Quaternion,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing)

Constructs a `VectorField` from the passed Hamiltonian `h`. Explicity, 
for `h`, constructs a vector field `F` such that
  
`F.x = [-∂h/∂p₁, ∂h/∂q₁, ...]`
"""
function VectorField(h::TPS; Q::Union{Quaternion,Nothing}=nothing, spin::Union{Bool,Nothing}=nothing)
  if !isnothing(Q)
    T = Vector{promote_type(typeof(h),eltype(Q))}
  else
    T = Vector{typeof(h)}
  end

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

  outF = zero(VectorField{T,U},use=h)
  return outF

  nv = numvars(h)
  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), require >=$nv"

  map!(t->t.tpsa, work_low, outF.x)
  vec2fld!(nv, h.tpsa, work_low)

  return outF
end

"""
    VectorField(m::DAMap)

Creates a `VectorField` from the map
"""
function VectorField(m::DAMap)
  F =  zero(VectorField{typeof(m.x),typeof(m.Q)}, use=m) 
  for i=1:numvars(F)
    F.x[i] = copy(m.x[i])
  end

  if !isnothing(m.Q)
    F.Q.q0 = copy(m.Q.q0)
    F.Q.q1 = copy(m.Q.q1)
    F.Q.q2 = copy(m.Q.q2)
    F.Q.q3 = copy(m.Q.q3)
  end
  return F
end


# --- zero ---
function zero(F::VectorField{T,U}) where {T,U}
  return zero(typeof(F), use=F)
end

function zero(::Type{VectorField{T,U}}; use::UseType=GTPSA.desc_current) where {T,U}
  desc = getdesc(use)
  nv = numvars(desc)
  np = numparams(desc)

  x = similar(T, nv) 
  Base.require_one_based_indexing(x)

  for i=1:nv
    @inbounds x[i] = eltype(x)(use=desc)
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

  return VectorField(x, Q)
end

#= --- one ---
"""
    one(F::VectorField)
  
Construct a 
"""
function one(F::VectorField)
  return one(typeof(F), use=F, idpt=F.idpt)
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
=#