function Probe(x::Vector{T}; x0::Vector{S}=zeros(length(x)), Q::U=nothing, E::V=nothing, spin::Union{Bool,Nothing}=nothing, FD::Union{Bool,Nothing}=nothing, idpt::W=nothing) where {S,T,U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing},W<:Union{Nothing,Bool}}
  length(x) == length(x0) || error("Length of orbital ray != length of reference orbit vector!")

  if isnothing(spin)
    Q1 = Q
  elseif spin
    if isnothing(Q)
      Q1 = Quaternion(first(x))
    else
      Q1 = Q1
    end
  else
    error("For no spin tracking, please omit the spin kwarg or set spin=nothing") # For type stability
    #Q1 = nothing # For type instability
  end

  if isnothing(FD)
    E1 = E
  elseif FD
    if isnothing(E)
      E1 = zeros(eltype(x0), length(x), length(x)) 
    else
      (length(x),length(x)) == size(E) || error("Size of FD matrix inconsistent with number of variables!")
      E1 = E
    end
  else
    error("For no stochasticity, please omit the FD kwarg or set FD=nothing") # For type stability
    #E1 = nothing # for type instability
  end

  return Probe{S,T,typeof(Q1),typeof(E1),W}(x0, x, Q1, E, idpt)
end


# Copy ctor:
function Probe(p::Probe{S,T,U,V,W}) where {S,T,U,V,W}

  x = Vector{T}(undef,length(p.x))

  for i=1:length(p.x)
    @inbounds x[i] = T(p.x[i])
  end

  if isnothing(p.Q)
    Q = nothing
  else
    q = Vector{T}(undef,4)
    for i=1:4
      @inbounds q[i] = T(p.Q.q[i])
    end
    Q = Quaternion(q)
  end

  if isnothing(p.E)
    E = nothing
  else
    E = copy(m.E)
  end

  return Probe{S,T,U,V,W}(copy(p.x0), x, Q, E, p.idpt)
end

==(p1::Probe, p2::Probe) = (p1.x0 == p2.x0 && p1.x == p2.x && p1.q == p2.q && p1.E == p2.E && p1.idpt == p2.idpt)
