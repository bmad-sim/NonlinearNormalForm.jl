function copy!(m::Union{TaylorMap,VectorField}, m1::Union{TaylorMap,VectorField})
  checkstates(m, m1)
  if m1 isa TaylorMap && m isa TaylorMap
    m.x0 .= m1.x0
  end
  NV = nvars(m)
  foreach((xi, x1i)->TI.copy!(xi, x1i), view(m.x, 1:NV), m1.x)
  if !isnothing(m.q)
    foreach((qi, q1i)->TI.copy!(qi, q1i), m.q, m1.q)
  end
  if !isnothing(m.s) && m1 isa TaylorMap
    m.s .= m1.s
  end
  return m
end

copy(m::Union{TaylorMap,VectorField}) = (out_m = zero(m); copy!(out_m, m); return out_m)

"""
    zero(m::TaylorMap)

Creates a zero `m` with the same properties as `m` including TPSA definiton,
spin, and stochasticity.
"""
function zero(m::TaylorMap) 
  return _zero_map(typeof(m), getdef(m), m)
end

for t = (:DAMap, :TPSAMap)
@eval begin

# Lowest-level, internal
function _zero_map(::Type{$t{X0,X,Q,S}}, def::AbstractTPSADef, reuse::Union{TaylorMap,Nothing}) where {X0,X,Q,S}
  out_x0 = init_x0(X0, def)
  out_x = init_x(X, def, reuse)
  out_q = init_q(Q, def)
  out_s = init_s(S, def)
  out_m = $t(out_x0, out_x, out_q, out_s)
  return out_m
end

function zero(::Type{$t{X0,X,Q,S}}) where {X0,X,Q,S}
  return _zero_map($t{X0,X,Q,S}, getdef(eltype(X)), nothing)
end
#=
# Explicit type specification
# Def change would be static (in type)
# The current problem with this is that dynamic gtpsa resolution would have a problem
# because GPSA.desc_current will be used
function $t{X0,X,Q,S}(m::Union{TaylorMap,Nothing}=nothing) where {X0,X,Q,S}
  out_m = _zero_map($t{X0,X,Q,S}, getdef(eltype(X)), m)
  if !isnothing(m)
    copy!(out_m, m)
  end
  return out_m
end
=#
# Copy ctor including optional TPSA def change
function $t(m::TaylorMap; def::AbstractTPSADef=getdef(m))
  checktpsas(m, def)
  X0 = typeof(m.x0)
  X = similar_eltype(typeof(m.x), TI.init_tps_type(eltype(X0), def))
  Q = isnothing(m.q) ? Nothing : Quaternion{TI.init_tps_type(eltype(X0), def)}
  S = typeof(m.s)
  out_m = _zero_map($t{X0,X,Q,S}, def, m)
  copy!(out_m, m)
  return out_m
end

# Kwarg-esque ctor:
function $t(;
  def::Union{AbstractTPSADef,Nothing}=nothing,
  x0::Union{AbstractVector,Nothing}=nothing,
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
  s::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
  stochastic::Union{Bool,Nothing}=nothing,
) 
  if isnothing(def)
    if !isnothing(x) && TI.is_tps_type(eltype(x)) isa TI.IsTPSType
      def = getdef(first(x))
    elseif !isnothing(q) && TI.is_tps_type(eltype(q)) isa TI.IsTPSType
      def = getdef(first(q))
    else
      error("No TPSA definition has been provided, nor is one inferrable from the input arguments!")
    end
  end


  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? TI.numtype(eltype(t)) : Float64), (x0, x, q, s))...)
  TW = TI.init_tps_type(W, def)
  X0 = isnothing(x0) ? @_DEFAULT_X0(nvars(def)){W} : similar_eltype(typeof(x0), W)
  X = isnothing(x) ? @_DEFAULT_X(ndiffs(def)){TW} : similar_eltype(typeof(x), TW)
  
  if isnothing(spin)
    Q = isnothing(q) && isnothing(q_map) ? Nothing : Quaternion{TW}
  elseif spin
    Q = Quaternion{TW}
  else
    Q = Nothing
  end

  if isnothing(stochastic)
    S = isnothing(s) ? Nothing : similar_eltype(typeof(s), W)
  elseif stochastic
    S = isnothing(s) ? @_DEFAULT_S(nvars(def)){W} : similar_eltype(typeof(s), W)
  else
    S = Nothing
  end

  # Construct map:
  out_x0 = init_x0(X0, def)
  out_x = init_x(X, def)
  out_q = init_q(Q, def)
  out_s = init_s(S, def)

  out_m = $t(out_x0, out_x, out_q, out_s)


  if !isnothing(x0)
    out_m.x0 .= x0
  end

  setray!(out_m.x, x=x, x_matrix=x_matrix)
  if !isnothing(out_m.q) && !isnothing(q)
    setquat!(out_m.q, q=q, q_map=q_map)
  end

  if !isnothing(out_m.s) && !isnothing(s)
    out_m.s .= s
  end

  return out_m
end

end
end
#=

"""
    zero(m::TaylorMap)

Creates a zero `m` with the same properties as `m` including GTPSA `Descriptor`,
spin, and stochasticity.
"""
zero(m::TaylorMap) = typeof(m)(m)
=#
"""
    one(m::TaylorMap)
  
Creates an identity `m` with the same properties as `m`, including GTPSA
`Descriptor`, spin, and stochasticity.
"""
function one(m::TaylorMap)
  out_m = zero(m)
  nv = nvars(m)

  for i in 1:nv
    TI.seti!(out_m.x[i], 1, i)
  end

  if !isnothing(m.q)
    TI.seti!(out_m.q.q0, 1, 0)
  end

  return out_m
end



#=

"""
    $($t){X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  
A `$($t)` is a `TaylorMap` which composes and inverts $( $t == DAMap ? "with the scalar part IGNORED" : "with the scalar part INCLUDED").
Both `DAMap` and `TPSAMap` have the exact same internal structures. See the documentation for 
`TaylorMap` for more details of the differences.

# Constructors

1. Explicit type specification and optional copy:

```julia
$($t){X0,X,Q,S}(m::Union{TaylorMap,Nothing}=nothing) where {X0,X,Q,S}
```

where the types `{X0,X,Q,S}` of each field are specified (see `TaylorMap`). The TPSA definition 
is inferred from `eltype(X)`, and must have the same number of variables + parameters as `m` if 
provided. The result is a `$($t){X0,X,Q,S}` with values equal to `m`, up to any truncations, if 
provided.

------

2. A copy constructor with possible change of `Descriptor` via

```julia
$($t)(m::Union{TaylorMap,Nothing}=nothing; def::Union{AbstractTPSADef,Nothing}=nothing)
```

The number variables and number of parameters must agree for a `Descriptor` change.

------

3. The general purpose constructor

```julia
$($t)(;
  use::UseType=GTPSA.desc_current,
  x0::Union{AbstractVector,Nothing}=nothing,
  x::Union{AbstractVector,Nothing}=nothing,
  x_matrix::Union{AbstractMatrix,UniformScaling,Nothing}=nothing,
  q::Union{Quaternion,AbstractVector,UniformScaling,Nothing}=nothing,
  q_map::Union{AbstractMatrix,Nothing}=nothing,
  s::Union{AbstractMatrix,Nothing}=nothing,
  spin::Union{Bool,Nothing}=nothing,
  stochastic::Union{Bool,Nothing}=nothing,
) 
```
If provided, the keyword arguments have the following behaviors:
### Keyword Arguments
- `use`        -- Specifies which GTPSA `Descriptor` to use. Must be a `Descriptor`, `TPS`, `TaylorMap`, or `VectorField`
- `x0`         -- Sets the values of the reference orbit equal to `x0`. Must be an `AbstractVector`
- `x`          -- Sets the orbital ray equal to `x`. Must be an `AbstractVector`
- `x_matrix`   -- Sets the linear (and optionally nonlinear) part of the orbital ray equal to the transfer matrix `x_matrix`. This will be done after setting `x`. Must be an `AbstractMatrix` or `UniformScaling`
- `q`          -- Sets the quaternion equal to `q`. Must be a `Quaternion`, `AbstractVector`, or `UniformScaling`
- `q_map`      -- Sets the linear (and optionally nonlinear) part of the quaternion map equal to the transfer matrix `q_map`. This will be done after setting `q`. Must be an `AbstractMatrix`
- `s`          -- Sets the stochastic kicks envelope matrix equal to `s`
- `spin`       -- Boolean which specifies if spin tracking should be included. This is type-unstable and should not be used in performance-critical places
- `stochastic` -- Boolean which specifies if the envelope matrix of stochastic kicks should be included. This is type-unstable and should not be used in performance-critical places

All of the setter keyword arguments do NOT use the directly passed inputs in the map, rather set the 
values equal to those after construction of the map. So for example:

```julia
julia> v = collect(Float64, 1:6);

julia> m = DAMap(x0=v);

julia> m.x0 === v
false
```
"""
$t

=#