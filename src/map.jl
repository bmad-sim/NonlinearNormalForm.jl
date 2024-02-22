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
abstract type TaylorMap{S<:Number,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V<:Number} end 

struct DAMap{S<:Number,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V<:Number} <: TaylorMap{S,T,U,V}
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

struct TPSAMap{S<:Number,T<:Union{TPS,ComplexTPS},U<:Union{TPS,ComplexTPS},V<:Number} <: TaylorMap{S,T,U,V}
  x0::Vector{S}     # Entrance value of map
  x::Vector{T}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  q::Quaternion{U}  # Quaternion for spin
  E::Matrix{V}      # Envelope for stochastic radiation
end

# List of acceptable Probes to create TaylorMaps from (e.g. including/excluding orbit/spin):
# For now I will just assume both orbit and spin are TPSs but can change later:
#const TaylorProbe{S<:Real,V<:Real} = Probe{S,TPS,TPS,V} #Union{Probe{S,TPS,TPS,V}, Probe{S,W,TPS,V}, Probe{S,TPS,W,V}}

for t = (:DAMap, :TPSAMap)
  @eval begin
    # Function to create a TaylorMap from a valid Probe or other TaylorMap
    function $t(m::Union{TaylorMap{S,T,U,V},Probe{S,T,U,V}}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S<:Real, T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS}, V<:Real}
      x0 = deepcopy(m.x0)
      x = map(x->(T)(x, use=getdesc(use)), m.x)
      q = Quaternion(map(x->(U)(x, use=getdesc(use)), m.q.q))
      E = deepcopy(m.E)
      return $t{S,T,U,V}(x0, x, q, E)
    end
    
    # Create from vector and blank (in this case must check for consistency)
    function $t(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), q::Quaternion{U}=unit_quat(first(x)), E::Matrix{V}=zeros(numtype(first(x)), numvars(x), numvars(x)), use::Union{Descriptor,<:TaylorMap,Probe{<:Real,TPS,TPS,<:Real},Nothing}=nothing) where {S<:Number,T<:Union{TPS,ComplexTPS}, U<:Union{TPS,ComplexTPS},V<:Number}
      numvars(x) == length(x) || error("Length of orbital ray inconsistent with number of variables in GTPSA!")
      numvars(x) == length(x0) || error("Length of orbital ray != length of coordinate system origin vector!")
      (numvars(x),numvars(x)) == size(E) || error("Size of stochastic matrix inconsistent with number of variables!")
      (isnothing(use) && getdesc(x) == getdesc(q)) || error("Orbital ray Descriptor different from quaternion Descriptor!")
      
      x1 = map(x->(T)(x, use=getdesc(use)), x)
      #(eltype(x)).(x, use=getdesc(use))
      q1 = Quaternion(map(x->(U)(x, use=getdesc(use)), q.q))
      
      return $t{eltype(x0),T,U,eltype(E)}(deepcopy(x0), x1, q1, deepcopy(E))
    end
  end
end


# --- Operators --- 

# --- Compose ---
function ∘(m2::DAMap,m1::DAMap{S,T,U,V}) where {S,T,U,V}
  scalar.(m1.x) == m2.x0 || error("Disconnected DAMaps! Scalar part of first map != entrance coordinates of second map!")
  nv = numvars(m1)
  ref = Vector{numtype(T)}(undef, nv)
  for i=1:nv
    @inbounds ref[i] = m1.x[i][0]
    @inbounds m1.x[i][0] = 0
  end
  # What kind to get? Complex or Real?
  # there is probably a better way to do this 
  # but this is type stable
  # maybe we should make params also just take a TPS 
  # and return that TPSs params with the correct type
  if eltype(m1.x) == TPS
    k = params(use=getdesc(m1.x))
  else
    k = complexparams(use=getdesc(m1.x))
  end
  tmp = vcat(m1.x, k...)
  outx = ∘(m2.x, tmp)
  for i=1:nv
    @inbounds m1.x[i][0] = ref[i]
  end
  return DAMap(deepcopy(m1.x0), outx, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

function ∘(m2::TPSAMap,m1::TPSAMap)
  nv = numvars(m1)
  if eltype(m1.x) == TPS
    k = params(use=getdesc(m1.x))
  else
    k = complexparams(use=getdesc(m1.x))
  end
  outx = ∘(m2.x, vcat((m1.x0+m1.x), k...))
  return TPSAMap(deepcopy(m1.x0), outx, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

# --- Inverse ---
function inv(m1::DAMap)
  nv = numvars(m1)
  outx = map(x->x[0], m1.v)
  for t in m1.v t[0] = 0 end
  outv = inv(m1.v)
  return DAMap(outx, outv, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

