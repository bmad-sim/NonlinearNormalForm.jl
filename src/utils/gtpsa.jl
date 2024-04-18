# From GTPSA:
# --- Poisson bracket ---
# Low-level calls
poisbra!(tpsa1::Ptr{RTPSA}, tpsa2::Ptr{RTPSA}, tpsa::Ptr{RTPSA}, nv::Cint) = (@inline; GTPSA.mad_tpsa_poisbra!(tpsa1, tpsa2, tpsa, nv))
poisbra!(tpsa1::Ptr{RTPSA}, ctpsa1::Ptr{CTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_poisbrat!(ctpsa1,tpsa1,ctpsa, nv))
poisbra!(ctpsa1::Ptr{CTPSA}, tpsa1::Ptr{RTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_tpoisbra!(tpsa1, ctpsa1, ctpsa, nv))
poisbra!(ctpsa1::Ptr{CTPSA}, ctpsa2::Ptr{CTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; GTPSA.mad_ctpsa_poisbra!(ctpsa1,ctpsa2, ctpsa, nv))

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


# --- F . grad ---
fgrad!(na::Cint, ma::Vector{Ptr{RTPSA}}, b::Ptr{RTPSA}, c::Ptr{RTPSA}) = (@inline; mad_tpsa_fgrad!(na, ma, b, c))
fgrad!(na::Cint, ma::Vector{Ptr{CTPSA}}, b::Ptr{CTPSA}, c::Ptr{CTPSA}) = (@inline; mad_ctpsa_fgrad!(na, ma, b, c))


function fgrad!(g::T, F::Vector{<:T}, h::T; work_low::Vector{<:Union{TPS,ComplexTPS}}=Vector{lowtype(h)}(undef, numvars(h))) where {T<:Union{TPS,ComplexTPS}}
  nv = numvars(h)
  @assert length(F) == nv "Incorrect length of F; received $(length(F)), should be $nv"
  @assert length(work_low) >= nv "Incorrect length for work_low; received $(length(work_low)), should be >=$nv"
  @assert eltype(work_low) == lowtype(T) "Incorrect eltype of work_low; received $(eltype(work_low)), should be $(lowtype(T))"
  map!(t->t.tpsa, work_low, F)
  fgrad!(nv, work_low, h, g)
  return
end

"""
    fgrad(F::Vector{<:T}, h::T) where {T<:Union{TPS,ComplexTPS}}    

Calculates `F⋅∇g`.
"""
function fgrad(F::Vector{<:T}, h::T) where {T<:Union{TPS,ComplexTPS}}
  g = zero(h)
  fgrad!(g, F, h)
  return g
end














# --- getvectorfield ---
vec2fld!(na::Cint, tpsa::Ptr{RTPSA}, m::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_vec2fld!(na, tpsa, m))
vec2fld!(na::Cint, ctpsa::Ptr{CTPSA}, m::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_vec2fld!(na, ctpsa, m))

"""
    getvectorfield(h::Union{TPS,ComplexTPS})::Vector{<:typeof(h)}

Assuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically-
conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), calculates the vector field (Hamilton's 
equations) from the passed Hamiltonian, defined as `[∂h/∂p₁, -∂h/∂q₁, ...]`

# Example
```julia-repl
julia> d = Descriptor(2,10); x = vars();

julia> h = (x[1]^2 + x[2]^2)/2
TPS:
 Coefficient                Order   Exponent
  5.0000000000000000e-01      2      2   0
  5.0000000000000000e-01      2      0   2


julia> getvectorfield(h)
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:  -1.0000000000000000e+00      1      0   1
-------------------------------------------------
   2:   1.0000000000000000e+00      1      1   0
```
"""
function getvectorfield(h::Union{TPS,ComplexTPS})::Vector{<:typeof(h)}
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(h.tpsa).d))
  na = desc.nv
  mc = Vector{typeof(h)}(undef, na)
  for i in eachindex(mc)
    mc[i] = zero(h)
  end
  m = map(t->t.tpsa, mc)
  vec2fld!(na, h.tpsa, m)
  return mc
end

# --- gethamiltonian ---
fld2vec!(na::Cint, ma::Vector{Ptr{RTPSA}}, tpsa::Ptr{RTPSA}) = (@inline; mad_tpsa_fld2vec!(na, ma, tpsa))
fld2vec!(na::Cint, ma::Vector{Ptr{CTPSA}},  ctpsa::Ptr{CTPSA}) = (@inline; mad_ctpsa_fld2vec!(na, ma, ctpsa))

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


# --- exp(F . grad) m ---
exppb!(na::Cint, ma::Vector{Ptr{RTPSA}}, mb::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_exppb!(na, ma, mb, mc))
exppb!(na::Cint, ma::Vector{Ptr{CTPSA}}, mb::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_exppb!(na, ma, mb, mc))

"""
    exppb(F::Vector{<:Union{TPS,ComplexTPS}}, m::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))

Calculates `exp(F⋅∇)m = m + F⋅∇m + (F⋅∇)²m/2! + ...`. If `m` is not provided, it is assumed 
to be the identity. 

# Example

```julia-repl
julia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];

julia> time = 0.01; k = 2; m = 0.01;

julia> h = p^2/(2m) + 1/2*k*x^2;

julia> hf = getvectorfield(h);

julia> map = exppb(-time*hf, [x, p])
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:   9.9001665555952290e-01      1      1   0
   1:   9.9666999841313930e-01      1      0   1
-------------------------------------------------
   2:  -1.9933399968262787e-02      1      1   0
   2:   9.9001665555952378e-01      1      0   1
```
"""
function exppb(F::Vector{<:Union{TPS,ComplexTPS}}, m::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))
  descF = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(F[1].tpsa).d))
  if length(F) != descF.nv
    error("Vector length != number of variables in the GTPSA")
  end
  ma1, mb1 = promote(F, m)
  m1 = map(t->t.tpsa, ma1) 
  m2 = map(t->t.tpsa, mb1) 
  mc = zero.(ma1)
  m3 = map(t->t.tpsa, mc)
  GC.@preserve ma1 mb1 exppb!(Cint(length(F)), m1, m2, m3)  
  return mc
end




# --- mnrm ---
mnrm(na::Cint, ma::Vector{Ptr{RTPSA}})::Float64 = mad_tpsa_mnrm(na, ma)
mnrm(na::Cint, ma::Vector{Ptr{CTPSA}})::ComplexF64 = mad_ctpsa_mnrm(na, ma)

#=
"""
    norm(ma::Vector{<:Union{TPS,ComplexTPS}})

Calculates the norm of the map `ma`, defined as `sum(norm.(ma))` or the 
sum of the absolute value of all coefficients in each TPS.
"""
function norm(ma::Vector{<:Union{TPS,ComplexTPS}})
  return mnrm(Cint(length(ma)), map(x->x.tpsa, ma))
end
=#


# --- partial inversion ---
pminv!(na::Cint, ma::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}, select::Vector{Cint}) = (@inline; mad_tpsa_pminv!(na, ma, mc, select))
pminv!(na::Cint, ma::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}, select::Vector{Cint}) = (@inline; mad_ctpsa_pminv!(na, ma, mc, select))

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
