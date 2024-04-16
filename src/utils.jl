# helper functions
getdesc(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d))
numvars(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nv
numparams(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).np
numnn(m::Union{Probe{<:Real,<:Union{TPS,ComplexTPS},<:Any,<:Any},<:TaylorMap}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.x).tpsa).d)).nn

getdesc(t::Vector{<:Union{TPS,ComplexTPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d))
numvars(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nv
numparams(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).np
numnn(t::Vector{<:Union{TPS,ComplexTPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(t).tpsa).d)).nn

getdesc(m::Quaternion{<:Union{ComplexTPS,TPS}}) = Descriptor(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d))
numvars(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nv
numparams(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).np
numnn(m::Quaternion{<:Union{ComplexTPS,TPS}}) = unsafe_load(Base.unsafe_convert(Ptr{GTPSA.Desc}, unsafe_load(first(m.q).tpsa).d)).nn

numtype(m::Union{Probe{<:Real,T,<:Any,<:Any},<:TaylorMap{<:Number,T,<:Any,<:Any}}) where {T<:Union{TPS,ComplexTPS}} = numtype(T)

jacobian(m::TaylorMap;include_params=false) = jacobian(view(m.x, 1:numvars(m)),include_params=include_params)
jacobiant(m::TaylorMap;include_params=false) = jacobiant(view(m.x, 1:numvars(m)), include_params=include_params)
checksymp(m::TaylorMap) = checksymp(jacobian(m))


function read_fpp_map(file)
  data = readdlm(file, skipblanks=true)
  nv = data[findfirst(t->t=="Dimensional", data)- CartesianIndex(0,1)]
  no = data[findfirst(t->t=="NO", data) + CartesianIndex(0,2)]
  np = data[findfirst(t->t=="NV", data) + CartesianIndex(0,2)] - nv
  nn = nv+np
  # Make the TPSA
  d = Descriptor(nv, no, np, no)
  m = complex(DAMap(spin=true)) #repeat([ComplexTPS(use=d)], nv), Q=Quaternion([ComplexTPS(1,use=d), repeat([ComplexTPS(use=d)], 3)...]))

  idx=3
  data=data[3:end,:]
  # Now fill
  for i=1:nv
    idx = findfirst(x->(x isa Integer), data[:,1])
    count = 0
    while data[idx,1] >= 0
      a = data[idx,2]
      b = data[idx,3]
      ords = data[idx,4:end]
      m.x[i][collect(Int, ords)] = a + im*b
      idx += 1
      count += 1
    end
    if count != 0 && -data[idx,1] != count
      println(m)
      println(data[idx,1])
      println(count)
      error("This should not have been reached! Incorrect number of monomials read for variable $(i)")
    end
    idx += 1
    data=data[idx:end,:]
  end
  # dont forget params
  m.x[nv+1:nn] .= complexparams(d)


  # spin?
  idx = findfirst(t->t=="c_quaternion", data[:,1])
  if idx > 0
    if data[idx,3] == "identity"
      m.Q.q .= [1.0, 0., 0., 0.]
    else
      for i=1:4
        idx = findfirst(x->(x isa Integer), data[:,1])
        count = 0
        while data[idx,1] >= 0
          a = data[idx,2]
          b = data[idx,3]
          ords = data[idx,4:end]
          m.Q.q[i][collect(Int,ords)] = a + im*b
          idx += 1
          count += 1
        end
        if count!=0 && -data[idx,1] != count
          println(m)
          println(data[idx,1])
          println(count)
          error("This should not have been reached! Incorrect number of monomials read for variable $(i)")
        end
        idx += 1
        data=data[idx:end,:]
      end
    end 
  end
  return m
end


# From GTPSA:
# --- Poisson bracket ---
# Low-level calls
poisbra!(tpsa1::Ptr{RTPSA}, tpsa2::Ptr{RTPSA}, tpsa::Ptr{RTPSA}, nv::Cint) = (@inline; mad_tpsa_poisbra!(tpsa1, tpsa2, tpsa, nv))
poisbra!(tpsa1::Ptr{RTPSA}, ctpsa1::Ptr{CTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; mad_ctpsa_poisbrat!(ctpsa1,tpsa1,ctpsa, nv))
poisbra!(ctpsa1::Ptr{CTPSA}, tpsa1::Ptr{RTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; mad_ctpsa_tpoisbra!(tpsa1, ctpsa1, ctpsa, nv))
poisbra!(ctpsa1::Ptr{CTPSA}, ctpsa2::Ptr{CTPSA}, ctpsa::Ptr{CTPSA}, nv::Cint) = (@inline; mad_ctpsa_poisbra!(ctpsa1,ctpsa2, ctpsa, nv))

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
  t = promote_type(typeof(f),typeof(g))(use=f)
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(f.tpsa).d))
  poisbra!(f.tpsa,g.tpsa,t.tpsa, desc.nv)
  return t
end

# --- Lie bracket ---
liebra!(na::Cint, m1::Vector{Ptr{RTPSA}}, m2::Vector{Ptr{RTPSA}}, m3::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_liebra!(na, m1, m2, m3))
liebra!(na::Cint, m1::Vector{Ptr{CTPSA}}, m2::Vector{Ptr{CTPSA}}, m3::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_liebra!(na, m1, m2, m3))

"""
    lb(A::Vector{<:Union{TPS,ComplexTPS}}, F::Vector{<:Union{TPS,ComplexTPS}})

Computes the Lie bracket of the vector functions `A` and `F`, defined over N variables as 
`Σᵢᴺ Aᵢ (∂F/∂xᵢ) - Fᵢ (∂A/∂xᵢ)`

# Example
```julia-repl
julia> d = Descriptor(2,10); x = vars();

julia> A = [-x[2], x[1]]
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:  -1.0000000000000000e+00      1      0   1
-------------------------------------------------
   2:   1.0000000000000000e+00      1      1   0


julia> F = [-x[1]^2, 2*x[1]*x[2]]
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:  -1.0000000000000000e+00      2      2   0
-------------------------------------------------
   2:   2.0000000000000000e+00      2      1   1


julia> lb(A,F)
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:   4.0000000000000000e+00      2      1   1
-------------------------------------------------
   2:   3.0000000000000000e+00      2      2   0
   2:  -2.0000000000000000e+00      2      0   2
```
"""
function lb(A::Vector{<:Union{TPS,ComplexTPS}}, F::Vector{<:Union{TPS,ComplexTPS}})
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(A[1].tpsa).d))
  if length(A) != desc.nv || length(F) != desc.nv
    error("Vector length != number of variables in the GTPSA")
  end
  A1, F1 = promote(A, F)
  m1 = map(t->t.tpsa, A1)  
  m2 = map(t->t.tpsa, F1) 
  mc = zero.(A1)
  m3 = map(t->t.tpsa, mc)
  GC.@preserve A1 F1 liebra!(Cint(length(A)), m1, m2, m3)     
  return mc
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

# --- logpb ---
logpb!(na::Cint, ma::Vector{Ptr{RTPSA}}, mb::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_logpb!(na, ma, mb, mc))
logpb!(na::Cint, ma::Vector{Ptr{CTPSA}}, mb::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_logpb!(na, ma, mb, mc))

"""
    logpb(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))

Given a final map `mf` and initial map `mi`, this function calculates the vector field `F`
such that `mf=exp(F⋅∇)mi`. If `mi` is not provided, it is assumed to be the identity.

```julia-repl
julia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];

julia> time = 0.01; k = 2; m = 0.01;

julia> h = p^2/(2m) + 1/2*k*x^2;

julia> hf = getvectorfield(h);

julia> map = exppb(-time*hf);

julia> logpb(map) == -time*hf
true
```
"""
function logpb(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(mf)))
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(mf[1].tpsa).d))
  if length(mf) != desc.nv || length(mi) != desc.nv
    error("Vector length != number of variables in the GTPSA")
  end
  ma1, mb1 = promote(mf, mi)
  m1 = map(t->t.tpsa, ma1) 
  m2 = map(t->t.tpsa, mb1) 
  mc = zero.(ma1)
  m3 = map(t->t.tpsa, mc)
  GC.@preserve ma1 mb1 logpb!(Cint(length(mf)), m1, m2, m3)  
  return mc
end

# --- F . grad ---
fgrad!(na::Cint, ma::Vector{Ptr{RTPSA}}, b::Ptr{RTPSA}, c::Ptr{RTPSA}) = (@inline; mad_tpsa_fgrad!(na, ma, b, c))
fgrad!(na::Cint, ma::Vector{Ptr{CTPSA}}, b::Ptr{CTPSA}, c::Ptr{CTPSA}) = (@inline; mad_ctpsa_fgrad!(na, ma, b, c))

"""
    fgrad(F::Vector{<:Union{TPS,ComplexTPS}}, g::Union{TPS,ComplexTPS})

Calculates `F⋅∇g`.
"""
function fgrad(F::Vector{<:Union{TPS,ComplexTPS}}, g::Union{TPS,ComplexTPS})
  descF = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(F[1].tpsa).d))
  if length(F) != descF.nv
    error("Vector length != number of variables in the GTPSA")
  end
  type = promote_type(typeof(F[1]), typeof(g))
  ma1 = convert(Vector{type},  F)
  b1 = convert(type, g)
  m1 = map(t->t.tpsa, ma1) 
  c = zero(b1)
  GC.@preserve ma1 fgrad!(Cint(length(ma)), m1, b1.tpsa, c.tpsa)  
  return c
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
# --- map inversion ---
minv!(na::Cint, ma::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_minv!(na, ma, mc))
minv!(na::Cint, ma::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_minv!(na, ma, mc))

"""
    inv(ma::Vector{<:Union{TPS,ComplexTPS}})

Inverts the map `ma` such that `ma ∘ inv(ma) = 1` in the variables.

# Example

```julia-repl

julia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];

julia> time = 0.01; k = 2; m = 0.01;

julia> h = p^2/(2m) + 1/2*k*x^2;

julia> hf = getvectorfield(h);

julia> map = exppb(-time*hf, [x, p]);

julia> map ∘ inv(map)
2-element Vector{TPS}:
  Out  Coefficient                Order   Exponent
-------------------------------------------------
   1:   1.0000000000000000e+00      1      1   0
-------------------------------------------------
   2:   1.0000000000000002e+00      1      0   1
```
"""
function inv(ma::Vector{<:Union{TPS,ComplexTPS}})
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(ma[1].tpsa).d))
  if length(ma) != desc.nv
    error("Map length != number of variables in the GTPSA")
  end
  mc = zero.(ma)
  ma1 = map(x->x.tpsa, ma)
  mc1 = map(x->x.tpsa, mc)
  minv!(Cint(length(ma)), ma1, mc1)
  return mc
end

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
