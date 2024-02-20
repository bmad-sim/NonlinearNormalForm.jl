#=
TPSAMap and DAMap used for normal form analysis. Maps can be 
constructed from Probes, GTPSA Descriptors, or other Maps. If 
no argument is provided, then GTPSA.desc_current is used to generate
a zero map (all coefficients zero with zero entrance value).

The GTPSA Descriptor can be changed when constructing a Map from a
Probe. The numbers of variables and parameters in the GTPSAs must agree.
=#
struct DAMap <: TaylorMap
  x::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  x0::Vector{ComplexF64}     # Entrance value of map
  q::Quaternion{ComplexTPS}  # Quaternion for spin
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
end

struct TPSAMap <: TaylorMap
  x::Vector{ComplexTPS}      # Expansion around x0, with scalar part equal to EXIT value of map wrt initial coordinates x0
  x0::Vector{ComplexF64}     # Entrance value of map
  q::Quaternion{ComplexTPS}  # Quaternion for spin
  E::Matrix{ComplexF64}      # Envelope for stochastic radiation
end

# maps can default to empty using GTPSA.desc_current, unlike Probes which are also used in tracking regular

# Quasi-keyword argument constructor for easy expansion/modification
function DAMap(; x::Vector{ComplexTPS}, x0::Vector{ComplexF64}, q::Quaternion{ComplexTPS}, E::Matrix{ComplexF64})
  return DAMap(x, x0, q, E)
end

# Quasi-keyword argument constructor for easy expansion/modification
function TPSAMap(; x::Vector{ComplexTPS}, x0::Vector{ComplexF64}, q::Quaternion{ComplexTPS}, E::Matrix{ComplexF64})
  return TPSAMap(x, x0, q, E)
end

# --- Definitions ---
function DAMap(m::Union{<:TaylorMap,Probe,Vector{<:Union{TPS,ComplexTPS}}}; use::Union{Descriptor,<:TaylorMap,Nothing}=nothing)
  return low_map(DAMap, m, use)
end

function TPSAMap(m::Union{<:TaylorMap,Probe}; use::Union{Descriptor,<:TaylorMap,Nothing}=nothing)
  return low_map(TPSAMap, m, use)
end

# --- Low-level ctors ---
function low_map(type::Type, m::Union{<:TaylorMap,Probe{TPS}}, use::Nothing)
  return (type)(x=ComplexTPS.(m.x),x0=convert(Vector{ComplexF64}, m.x0), q=Quaternion(ComplexTPS.(m.q.q)), E=convert(Matrix{ComplexF64}, m.E))
end

function low_map(type::Type, m::Probe{Float64}, use::Union{Descriptor,<:TaylorMap})
  length(m.x0) == numvars(use) || error("Number of variables in GTPSA != number of variables in Probe!")
  return (type)(x=ComplexTPS.(m.x,use=getdesc(use)),x0=convert(Vector{ComplexF64}, m.x0), q=Quaternion(ComplexTPS.(m.q.q,use=getdesc(use))), E=convert(Matrix{ComplexF64}, m.E)) # Error checking done in GTPSA
end

function low_map(type::Type, m::Probe{Float64}, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

function low_map(type::Type, m::Nothing, use::Union{Descriptor,TaylorMap})
  nv = numvars(use)
  x0 = zeros(ComplexF64, nv) 
  x = Vector{ComplexTPS}(undef, nv)
  for i in eachindex(x)
    x[i] = ComplexTPS(use=getdesc(use))
  end
  return (type)(x=x,x0=x0,q=Quaternion((ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)), ComplexTPS(use=getdesc(use)))), E=zeros(nv, nv))
end

function low_map(type::Type, m::Nothing, use::Nothing)
  return low_map(type, m, GTPSA.desc_current)
end

# Change Descriptor:
function low_map(type::Type, m::Union{<:TaylorMap,Probe{TPS}}, use::Union{Descriptor,<:TaylorMap})
  return (type)(x=ComplexTPS.(m.x, use=getdesc(use)),x0=convert(Vector{ComplexF64}, m.x0), q=Quaternion(ComplexTPS.(m.q.q, use=getdesc(use))), E=convert(Matrix{ComplexF64}, m.E))
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
  return DAMap(x=outx, x0=deepcopy(m1.x0), q=Quaternion((ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)), ComplexTPS(use=getdesc(m1)))), E=zeros(nv, nv))
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

