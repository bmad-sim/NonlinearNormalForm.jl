#=
TPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.

Once again parametric types, however being a TaylorMap we now require the 
orbit x and spin q to be TPSs.
=#
abstract type TaylorMap end 

struct DAMap{S<:Number,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V<:Number} <: TaylorMap
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

struct TPSAMap{S<:Number,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V<:Number} <: TaylorMap
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

# List of acceptable Probes to create TaylorMaps from (e.g. including/excluding orbit/spin):
# For now I will just assume both orbit and spin are TPSs but can change later:
const TaylorProbe{S<:Real,V<:Real,W<:Real} = Probe{S,TPS,TPS,V} #Union{Probe{S,TPS,TPS,V}, Probe{S,W,TPS,V}, Probe{S,TPS,W,V}}

for t = (:DAMap, :TPSAMap)
  @eval begin
    # Function to create a TaylorMap from a TaylorProbe (assuming consistent definition)
    function $t(p::TaylorProbe; use::Union{Descriptor,<:TaylorMap,TaylorProbe,Nothing}=nothing)
      return $t(deepcopy(p.x0), (eltype(p.x)).(p.x,use=getdesc(use)), Quaternion((eltype(p.q.q)).(p.q.q, use=getdesc(use))), deepcopy(p.E))
    end

    # Copy ctor/change ctor (again assumes consistent definition)
    function $t(m::TaylorMap; use::Union{Descriptor,<:TaylorMap,TaylorProbe,Nothing}=nothing)
      return $t(deepcopy(m.x0),(eltype(p.x)).(m.x,use=getdesc(use)), Quaternion((eltype(p.x)).(m.q.q, use=getdesc(use))), deepcopy(m.E))
    end

    # For some reason type unstable if TPS type for quaternion? huh?
    # Create from vector and blank (in this case must check for consistency)
    function $t(x::Vector{<:Union{TPS,ComplexTPS}}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{<:Number}=zeros(numtype(first(x)), numvars(x)), q::Quaternion{<:Union{TPS,ComplexTPS}}, E::Matrix{<:Number}=zeros(numtype(first(x)), numvars(x), numvars(x)), use::Union{Descriptor,<:TaylorMap,TaylorProbe,Nothing}=nothing)
      numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
      numvars(x) == length(x0) || error("Length of orbital ray != length of coordinate system origin vector!")
      (numvars(x),numvars(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      x1 = (eltype(x)).(x, use=getdesc(use))
      if isnothing(q)
        q1 = unit_quat(first(x1))
      else
        q1 = Quaternion((eltype(x)).(q.q, use=getdesc(use)))
      end
      
      return $t(x0, x, q1, E)
    end
  end
end

# --- Operators --- 

# Composition
function ∘(m2::DAMap,m1::DAMap)
  scalar.(m1) == m2.x0 || error("Disconnected DAMaps! Scalar part of first map != entrance coordinates of second map!")
  
  nv = numvars(m1)
  ref = Vector{ComplexF64}(undef, nv)
  for i=1:nv
    @inbounds ref[i] = m1.x[i][0]
    @inbounds m1.x[i][0] = 0
  end
  outx = ∘(m2.x, vcat(m1.x, complexparams(use=getdesc(m1))...))
  for i=1:nv
    @inbounds m1.x[i][0] = ref[i]
  end
  return DAMap(outx, deepcopy(m1.x0), Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

function ∘(m2::TPSAMap,m1::TPSAMap)
  nv = numvars(m1)
  outx = ∘(m2.v, vcat((m1.x0+m1.v), complexparams(use=getdesc(m1))...))
  return TPSAMap(x0=m1.x0, x=outx, q=Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), E=zeros(nv, nv))
end

# Inverse
function inv(m1::DAMap)
  nv = numvars(m1)
  outx = map(x->x[0], m1.v)
  for t in m1.v t[0] = 0 end
  outv = inv(m1.v)
  return DAMap(outx, outv, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

