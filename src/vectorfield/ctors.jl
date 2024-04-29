# --- getvectorfield ---
vec2fld!(na::Cint, tpsa::Ptr{RTPSA}, m::Vector{Ptr{RTPSA}}) = (@inline; GTPSA.mad_tpsa_vec2fld!(na, tpsa, m))
vec2fld!(na::Cint, ctpsa::Ptr{CTPSA}, m::Vector{Ptr{CTPSA}}) = (@inline; GTPSA.mad_ctpsa_vec2fld!(na, ctpsa, m))

"""
    VectorField(h::T; Q::U=nothing, spin::Union{Bool,Nothing}=nothing, work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(h)}(undef,numvars(h))) where {T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing}}

Constructs a `VectorField` from the passed Hamiltonian `h`. Explicity, 
for `h`, constructs a vector field `F` such that
  
`F.x = [∂h/∂p₁, -∂h/∂q₁, ...]`
"""
function VectorField(h::T; Q::U=nothing, spin::Union{Bool,Nothing}=nothing, work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(h)}(undef,numvars(h))) where {T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing}}
  nv = numvars(h)
  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), require >=$nv"

  x = Vector{T}(undef, nv)
  for i=1:nv
    @inbounds x[i] = T(use=h)
  end

  map!(t->t.tpsa, work_low, x)
  vec2fld!(nv, h.tpsa, work_low)

  if isnothing(spin)
    Q1 = Q
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(h)
    else
      Q1 = Q1
    end
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #Q1 = nothing # For type instability
  end
  return VectorField(x, Q1)
end

function VectorField{T,U}(h::T; Q::Union{U,Nothing}=nothing, work_low::Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}=Vector{lowtype(h)}(undef,numvars(h))) where {T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing}}
  nv = numvars(h)
  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), require >=$nv"

  x = Vector{T}(undef, nv)
  for i=1:nv
    @inbounds x[i] = T(use=h)
  end

  map!(t->t.tpsa, work_low, x)
  vec2fld!(nv, h.tpsa, work_low)

  if U != Nothing
    if isnothing(Q)
      Q1 = Quaternion(h)
    else
      Q1 = Q
    end
  else
    if !isnothing(Q)
      error("Quaternion cannot be specified in constructor for VectorField{$T,$U}")
    end
    Q1 = nothing
  end

  return VectorField(x, Q1)
end

"""
    VectorField{T,U}(u::UndefInitializer; use::UseType=GTPSA.desc_current) where {T,U}

Creates an undefined `VectorField{T,U}` with same number of variables as `use`.
"""
function VectorField{T,U}(u::UndefInitializer; use::UseType=GTPSA.desc_current) where {T,U}
  desc = getdesc(use)
  nv = numvars(desc)
  x = Vector{T}(undef, nv)

  if U != Nothing
    q = Vector{T}(undef, 4)
    Q = Quaternion(q)
  else
    Q = nothing
  end

  return VectorField(x,Q)
end

"""
    VectorField{T,U}(u::UndefInitializer, nv::Integer) where {T,U}

Creates an undefined `VectorField{T,U}` with specified number of variables `nv`.
"""
function VectorField{T,U}(u::UndefInitializer, nv::Integer) where {T,U}
  x = Vector{T}(undef, nv)

  if U != Nothing
    q = Vector{T}(undef, 4)
    Q = Quaternion(q)
  else
    Q = nothing
  end

  return VectorField(x,Q)
end


#=
"""
    VectorField(m::DAMap{S,T,U,V}) where {S,T,U,V}

Create a `VectorField` from the orbital (and quaternion) part
of the `DAMap`.
"""
function VectorField(m::DAMap{S,T,U,V}) where {S,T,U,V}


end=#

# --- zero ---
function zero(F::VectorField{T,U}) where {T,U}
  desc = getdesc(F)
  nv = numvars(desc)
  x = Vector{T}(undef, nv)
  for i=1:nv
    @inbounds x[i] = T(use=desc)
  end

  if isnothing(F.Q)
    Q = nothing
  else
    q = Vector{T}(undef, 4)
    for i=1:4
      @inbounds q[i] = T(use=desc)
    end
    Q =  Quaternion(q)
  end
  return VectorField(x,Q)
end

function zero(::Type{VectorField{T,U}}; use::Union{Descriptor,TPS,ComplexTPS,TaylorMap,Probe{<:Any,Union{TPS,ComplexTPS},<:Any,<:Any}}=GTPSA.desc_current) where {T,U}
  desc = getdesc(use)
  nn = numnn(desc)
  nv = numvars(desc)
  np = numparams(desc)

  x = Vector{T}(undef, nn)
  for i=1:nv
    @inbounds x[i] = T(use=desc)
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

  return VectorField{T,U}(x, Q)
end