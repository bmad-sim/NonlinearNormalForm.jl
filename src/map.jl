#=
TPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.
=#
struct DAMap <: TaylorMap
  x0::Vector{ComplexF64}     # Entrance value of map
  v::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map
  q::Quaternion{ComplexTPS}  # Quaternion for spin
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
end

struct TPSAMap <: TaylorMap
  x0::Vector{ComplexF64}     # Entrance value of map
  v::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map
  q::Quaternion{ComplexTPS}  # Quaternion for spin
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
end

function DAMap(m::Union{TaylorMap,Probe,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(DAMap, m, use)
end

function TPSAMap(m::Union{TaylorMap,Probe,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(TPSAMap, m, use)
end

# --- Probe ---
# If Probe contains TPS, then we will use that Descriptor:
function low_map(type::Type, m::Probe{TPS}, use::Nothing)
  return (type)(deepcopy(m.x0), deepcopy(m.v), deepcopy(m.q), deepcopy(m.E))
end

# If Probe contains TPS and Descriptor is specified, change Descriptor
function low_map(type::Type, m::Probe{TPS}, use::Union{Descriptor,TaylorMap})
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use)), m.v), Quaternion(map(x->ComplexTPS(x,use=getdesc(use)), m.q.q)), map(x->ComplexTPS(x,use=getdesc(use)), deepcopy(m.E))) # Error checking done in GTPSA
end

# If Probe does not contain TPS, use specified Descriptor:
function low_map(type::Type, m::Probe{Float64}, use::Union{Descriptor,TaylorMap})
  length(m.x0) == numvars(use) || error("Number of variables in GTPSA != number of variables in Probe!")
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use)), m.v), Quaternion(map(x->ComplexTPS(x,use=getdesc(use)), m.q.q)), map(x->ComplexTPS(x,use=getdesc(use)), deepcopy(m.E))) # Error checking done in GTPSA
end

# If probe does not contain TPS and no Descriptor specified, use GTPSA.desc_current
function low_map(type::Type, m::Probe{Float64}, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# --- TaylorMap ---

# Use Descriptor as in Map, copy
function low_map(type::Type, m::TaylorMap, use::Nothing)
  return (type)(deepcopy(m.x0), deepcopy(m.v), deepcopy(m.q), deepcopy(m.E))
end

# Change Descriptor
function low_map(type::Type, m::TaylorMap, use::Union{Descriptor,TaylorMap})
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use)), m.v), Quaternion(map(x->ComplexTPS(x,use=getdesc(use)), m.q.q)), map(x->ComplexTPS(x,use=getdesc(use)), deepcopy(m.E))) # Error checking done in GTPSA
end

# --- zero ---
function low_map(type::Type, m::Nothing, use::Union{Descriptor,TaylorMap})
  nv = numvars(use)
  x0 = zeros(ComplexF64, nv) 
  v = Vector{ComplexTPS}(undef, nv)
  for i in eachindex(v)
    v[i] = ComplexTPS(use=getdesc(use))
  end
  return (type)(x0, v, Quaternion((ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)))), zeros(nv, nv))
end

function low_map(type::Type, m::Nothing, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# --- Operators --- 

# Only difference is in map composition
function ∘(m2::DAMap,m1::DAMap)
  nv = numvars(m1)
  ref = Vector{ComplexF64}(undef, nv)
  for i=1:nv
    @inbounds ref[i] = m1.v[i][0]
    @inbounds m1.v[i][0] = 0
  end
  outv = ∘(m2.v, vcat(m1.v, complexparams(use=getdesc(m1))...))
  for i=1:nv
    @inbounds m1.v[i][0] = ref[i]
  end
  return DAMap(deepcopy(m1.x0), outv, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end

function ∘(m2::TPSAMap,m1::TPSAMap)
  nv = numvars(m1)
  comp = ∘(m2.x0+m2.v, vcat((m1.x0+m1.v), complexparams(use=getdesc(m1))...))
  x0 = map(x->x[0], comp)
  for t in comp t[0] = 0 end
  return TPSAMap(x0, comp, Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), zeros(nv, nv))
end