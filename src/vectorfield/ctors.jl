#=

Constructors for the VectorField. Includes zero, copy ctor (with possible 
descriptor change), a ctor from a hamiltonian, and zero_op.

=#

# --- zero ---
function zero(F::VectorField{T,U}) where {T,U}
  return zero(typeof(F), use=F)
end

function zero(::Type{VectorField{T,U}}; use::UseType=GTPSA.desc_current) where {T,U}
  desc = getdesc(use)
  nv = numvars(desc)
  x = similar(T, nv) 
  Base.require_one_based_indexing(x)

  for i=1:nv
    x[i] = eltype(x)(use=desc)
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


# --- copy ctor (map or vector field) --- 
"""
    VectorField(m::Union{TaylorMap,VectorField}, use::UseType=m)

Creates a copy `VectorField` from the map or vector field with 
possible GTPSA descriptor change if specified.
"""
function VectorField(m::Union{TaylorMap,VectorField}; use::UseType=m)
  F = zero(VectorField{typeof(m.x),typeof(m.Q)}, use=use) 
  nv = numvars(F)
  for i=1:nv
    setTPS!(F.x[i], m.x[i], change=true)
  end

  # set quaternion
  if !isnothing(F.Q)
    setTPS!(F.Q.q0, m.Q.q0, change=true)
    setTPS!(F.Q.q1, m.Q.q1, change=true)
    setTPS!(F.Q.q2, m.Q.q2, change=true)
    setTPS!(F.Q.q3, m.Q.q3, change=true)
  end
  return F
end



# --- from hamiltonian (getvectorfield) ---
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

  nv = numvars(h)
  GTPSA.vec2fld!(nv, h.tpsa, outF.x)

  if !isnothing(outF.Q) && !isnothing(Q)
    outF.Q .= Q
  end

  return outF
end


"""

    zero_op(F2::VectorField, F1::Union{VectorField,Number})

Returns a new zero VectorField with type equal to promoted type of `F1` and `F2`.
`F2` could be a number. This function is easy for VectorField because VectorFields 
do not carry the immutable parameters like maps do.
"""
function zero_op(F2::VectorField, F1::Union{VectorField,Number})
  outtype = promote_type(typeof(F1),typeof(F2))
  return zero(outtype, use=F2)
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