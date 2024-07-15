# From GTPSA:
#=
# --- Poisson bracket ---
# Low-level calls
poisbra!(tpsa1::TPS{Float64},     tpsa2::TPS{Float64},     tpsa::TPS{Float64},     nv::Cint) = (@inline; GTPSA.mad_tpsa_poisbra!(tpsa1, tpsa2, tpsa, nv))
poisbra!(tpsa1::TPS{Float64},     ctpsa1::TPS{ComplexF64}, ctpsa::TPS{ComplexF64}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_poisbrat!(ctpsa1,tpsa1,ctpsa, nv))
poisbra!(ctpsa1::TPS{ComplexF64}, tpsa1::TPS{Float64},     ctpsa::TPS{ComplexF64}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_tpoisbra!(tpsa1, ctpsa1, ctpsa, nv))
poisbra!(ctpsa1::TPS{ComplexF64}, ctpsa2::TPS{ComplexF64}, ctpsa::TPS{ComplexF64}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_poisbra!(ctpsa1,ctpsa2, ctpsa, nv))

"""
    pb!(h::Union{TPS,ComplexTPS}, f::Union{TPS,ComplexTPS}, g::Union{TPS,ComplexTPS})

In-place Poisson bracket, see the documentation for `pb` for more details.
"""
function pb!(h::Union{TPS,ComplexTPS}, f::Union{TPS,ComplexTPS}, g::Union{TPS,ComplexTPS})
  poisbra!(f.tpsa,g.tpsa,h.tpsa, numvars(f))
  return
end

"""
    pb(f::Union{TPS, ComplexTPS}, g::Union{TPS, ComplexTPS})

Assuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically-
conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), computes the Poisson bracket 
of the scalar functions `f` and `g`. The Poisson bracket of two functions `{f, g}` is defined as 
`Σᵢ (∂f/∂qᵢ)(∂g/∂pᵢ) - (∂g/∂qᵢ)(∂f/∂pᵢ)`.

# Examples
```julia-repl
julia> d = Descriptor(4,10);

julia> x = vars(d);

julia> f = (x[1]^2 + x[2]^2)/2 + (x[3]^2 + x[4]^2)/2;

julia> pb(f,x[1])
TPS:
  Coefficient              Order     Exponent
  -1.0000000000000000e+00    1        0    1    0    0


julia> pb(f,x[2])
TPS:
  Coefficient              Order     Exponent
   1.0000000000000000e+00    1        1    0    0    0


julia> pb(f,x[3])
TPS:
  Coefficient              Order     Exponent
  -1.0000000000000000e+00    1        0    0    0    1


julia> pb(f,x[4])
TPS:
  Coefficient              Order     Exponent
   1.0000000000000000e+00    1        0    0    1    0
```
"""
function pb(f::Union{TPS, ComplexTPS}, g::Union{TPS, ComplexTPS})
  h = promote_type(typeof(f),typeof(g))(use=f)
  pb!(h, f, g)
  return h
end
=#

# --- F . grad ---
function fgrad!(g::T, F::AbstractVector{<:T}, h::T) where {T<:TPS{<:Union{Float64,ComplexF64}}}
  Base.require_one_based_indexing(F)
  nv = numvars(h)
  @assert length(F) == nv "Incorrect length of F; received $(length(F)), should be $nv"
  @assert !(g === h) "Aliasing g === h not allowed for fgrad!"
  if T != ComplexTPS
    GTPSA.mad_tpsa_fgrad!(Cint(length(F)), F, h, g)
  else
    GTPSA.mad_ctpsa_fgrad!(Cint(length(F)), F, h, g)
  end
  return g
end

"""
    fgrad(F::AbstractVector{<:T}, h::T) where {T<:TPS{<:Union{Float64,ComplexF64}}}    

Calculates `F⋅∇h`.
"""
function fgrad(F::AbstractVector{<:T}, h::T) where {T<:TPS{<:Union{Float64,ComplexF64}}}
  g = zero(h)
  fgrad!(g, F, h)
  return g
end



#=





# --- gethamiltonian ---
fld2vec!(na::Cint, ma::Vector{TPS{Float64}}, tpsa::TPS{Float64}) = (@inline; mad_tpsa_fld2vec!(na, ma, tpsa))
fld2vec!(na::Cint, ma::Vector{TPS{ComplexF64}},  ctpsa::TPS{ComplexF64}) = (@inline; mad_ctpsa_fld2vec!(na, ma, ctpsa))

"""
    gethamiltonian(F::Vector{<:Union{TPS,ComplexTPS}})

Assuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically-
conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), this function calculates the Hamiltonian 
from a vector field `F` that can be obtained from a Hamiltonian (e.g. by `getvectorfield`). Explicitly, 
`∫ F₁ dp₁ - ∫ F₂ dq₁ + ... + ∫ F₂ₙ₋₁ dpₙ - ∫ F₂ₙ dqₙ `

# Example
```julia-repl
julia> d = Descriptor(2,10); x = vars();

julia> h = (x[1]^2 + x[2]^2)/2
TPS:
 Coefficient                Order   Exponent
  5.0000000000000000e-01      2      2   0
  5.0000000000000000e-01      2      0   2


julia> F = getvectorfield(h)
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:  -1.0000000000000000e+00      1      0   1
-------------------------------------------------
   2:   1.0000000000000000e+00      1      1   0


julia> gethamiltonian(F)
TPS:
 Coefficient                Order   Exponent
  5.0000000000000000e-01      2      2   0
  5.0000000000000000e-01      2      0   2
```
"""
function gethamiltonian(F::Vector{<:Union{TPS,ComplexTPS}})
  descF = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(F[1].tpsa).d))
  if length(F) != descF.nv
    error("Vector length != number of variables in the GTPSA")
  end
  h = zero(F[1])
  m1 = map(t->t.tpsa, F)
  fld2vec!(Cint(length(m)), m1, h.tpsa)
  return h
end










# --- partial inversion ---
pminv!(na::Cint, ma::Vector{TPS{Float64}}, mc::Vector{TPS{Float64}}, select::Vector{Cint}) = (@inline; mad_tpsa_pminv!(na, ma, mc, select))
pminv!(na::Cint, ma::Vector{TPS{ComplexF64}}, mc::Vector{TPS{ComplexF64}}, select::Vector{Cint}) = (@inline; mad_ctpsa_pminv!(na, ma, mc, select))

"""
    ptinv(ma::Vector{<:Union{TPS,ComplexTPS}}, vars::Vector{<:Integer})

Partially-inverts the map `ma`, inverting only the variables specified by index
in `vars`.
"""
function ptinv(ma::Vector{<:Union{TPS,ComplexTPS}}, vars::Vector{<:Integer})
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(ma[1].tpsa).d))
  if length(ma) != desc.nv
    error("Map length != number of variables in the GTPSA")
  end
  mc = zero.(ma)
  ma1 = map(x->x.tpsa, ma)
  mc1 = map(x->x.tpsa, mc)
  na = Cint(length(ma))
  select = zeros(Cint, na)
  select[vars] .= Cint(1)
  pminv!(na, ma1, mc1, select)
  return mc
end
=#