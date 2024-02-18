#=
GTPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.
=#
struct DAMap <: TaylorMap
  x0::Vector{ComplexF64}     # Entrance value of map
  v::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
  q::Quaternion{ComplexTPS}  # Quaternion for spin
end

struct GTPSAMap <: TaylorMap
  x0::Vector{ComplexF64}     # Entrance value of map
  v::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
  q::Quaternion{ComplexTPS}  # Quaternion for spin
end

function DAMap(m::Union{TaylorMap,Probe,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(DAMap, m, use)
end

function GTPSAMap(m::Union{TaylorMap,Probe,Nothing}=nothing; use::Union{Descriptor,TaylorMap,Nothing}=nothing)
  return low_map(GTPSAMap, m, use)
end

# --- Probe ---
# If Probe contains TPS, then we will use that Descriptor:
function low_map(type::Type, m::Probe{TPS}, use::Nothing)
  return deepcopy(m)
end

# If Probe contains TPS and Descriptor is specified, change Descriptor
function low_map(type::Type, m::Probe{TPS}, use::Union{Descriptor,TaylorMap})
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use), m.v)) # Error checking done in GTPSA
end

# If Probe does not contain TPS, use specified Descriptor:
function low_map(type::Type, m::Probe{Float64}, use::Union{Descriptor,TaylorMap})
  length(m.x0) == numvars(use) || error("Number of variables in GTPSA != number of variables in Probe!")
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use)), m.v))
end

# If probe does not contain TPS and no Descriptor specified, use GTPSA.desc_current
function low_map(type::Type, m::Probe{Float64}, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# --- TaylorMap ---

# Use Descriptor as in Map, copy
function low_map(type::Type, m::TaylorMap, use::Nothing)
  return (type)(deepcopy(m.x0), deepcopy(m.v))
end

# Change Descriptor
function low_map(type::Type, m::TaylorMap, use::Union{Descriptor,TaylorMap})
  return (type)(deepcopy(m.x0), map(x->ComplexTPS(x,use=getdesc(use)), m.v)) # error checking done in GTPSA
end

# --- zero ---
function low_map(type::Type, m::Nothing, use::Union{Descriptor,TaylorMap})
  nv = numvars(use)
  x0 = zeros(ComplexF64, nv) 
  v = Vector{ComplexTPS}(undef, nv)
  for i=1:nv
    @inbounds v[i] = ComplexTPS(use=getdesc(use))
  end
  return (type)(x0, v)
end

function low_map(type::Type, m::Nothing, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# --- Operators --- 