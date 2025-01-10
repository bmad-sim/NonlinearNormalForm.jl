for t = (:DAMap, :TPSAMap)
@eval begin

# Blank lowest-level ctor
function $t{X0,X,Q,S}(use::UseType=GTPSA.desc_current) where {X0,X,Q,S}
  getdesc(use).desc != C_NULL || error("GTPSA Descriptor not defined!")
  out_x0 = init_x0(X0, use)
  out_x = init_x(X, use)
  out_q = init_q(Q, use)
  out_s = init_s(S, use)

  return $t(out_x0, out_x, out_q, out_s)
end

# Copy ctor including optional Descriptor change
function $t(m::TaylorMap{X0,X,Q,S}; use::UseType=m) where {X0,X,Q,S}
  numvars(use) == numvars(m) || error("Number of variables in GTPSAs for `m` and `use` disagree!")
  numparams(use) == numparams(m) || error("Number of parameters in GTPSAs for `m` and `use` disagree!") 

  out_m = $t{X0,X,Q,S}(use)
  out_m.x0 = m.x0
  foreach((out_xi, xi)->setTPS!(out_xi, xi, change=true), view(out_m.x, 1:numvars(use)), m.x)
  if !isnothing(out_m.q)
    foreach((out_qi, Qi)->setTPS!(out_qi, Qi, change=true), out_m.q, m.q)
  end
  if !isnothing(out_m.s)
    out_m.s .= m.s
  end
  return out_m
end


"""
    $($t){X0,X,Q,S} <: TaylorMap{X0,X,Q,S}
  
A `$($t)` is a `TaylorMap` which composes and inverts $( $t == DAMap ? "with the scalar part IGNORED" : "with the scalar part INCLUDED").
Both `DAMap` and `TPSAMap` have the exact same internal structures. See the documentation for 
`TaylorMap` for more details of the differences.

# Constructors

1. At the lowest level, one may call

```julia
$($t){X0,X,Q,S}(use::UseType=GTPSA.desc_current) where {X0,X,Q,S}
```

where the types `{X0,X,Q,S}` of each field are specified (see `TaylorMap`), and the GTPSA 
`Descriptor` defined in `use` is used (see `UseType` for allowed types for `use`). A `TaylorMap` 
type is generally recommended for performance reasons as the already allocated `TPS`s for the 
parameters can be reused.

------

2. A copy constructor with possible change of `Descriptor` via

```julia
$($t)(m::TaylorMap{X0,X,Q,S}; use::UseType=m) where {X0,X,Q,S}
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

function $t(;
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
  # Assemble types:
  W = promote_type(map(t->(!isnothing(t) ? GTPSA.numtype(eltype(t)) : Float64), (x0, x, q, s))...)
  
  X0 = isnothing(x0) ? Vector{W} : promote_x0_type(typeof(x0), W)
  X = isnothing(x) ? Vector{TPS{W}} : promote_x_type(typeof(x), TPS{W})

  if isnothing(spin)
    Q = isnothing(q) && isnothing(q_map) ? Nothing : Quaternion{TPS{W}}
  elseif spin
    Q = Quaternion{TPS{W}}
  else
    Q = Nothing
  end

  if isnothing(stochastic)
    S = isnothing(s) ? Nothing : promote_s_type(typeof(s), W)
  elseif stochastic
    S = isnothing(s) ? Matrix{W} : promote_s_type(typeof(s), W)
  else
    S = Nothing
  end

  # Construct map using low level ctor:
  m = $t{X0,X,Q,S}(use)

  # Set if values provided:
  if !isnothing(x0)
    m.x0 .= x0
  end
  
  if !isnothing(x)
    length(x) <= numvars(use) || error("Length of input vector `x` cannot be greater than the number of variables in `use` GTPSA!")
    foreach((out_xi, xi)->setTPS!(out_xi, xi, change=true), view(m.x, 1:length(x)), x)
  end

  if !isnothing(x_matrix)
    if x isa AbstractMatrix # Map as a matrix:
      size(x_matrix,1) <= numvars(use) || error("Number of rows of `x_matrix` cannot be greater than the number of variables in `use` GTPSA!")
      size(x_matrix,2) <= GTPSA.numcoefs(first(m.x))-1 || error("Number of columns of `x_matrix` cannot be greater than the number of coefficients in `use` GTPSA!")
      for varidx in 1:size(x_matrix,1)
        m.x[varidx][1:size(x_matrix,2)] = view(x_matrix, varidx, :)
      end
    else # Uniform scaling: Making identity map
      for varidx in 1:numvars(use)
        m.x[varidx][varidx] = 1
      end
    end
  end

  if !isnothing(m.s) && !isnothing(s)
    m.s .= s
  end

  return m
end

end
end


"""
    zero(m::TaylorMap)

Creates a zero `m` with the same properties as `m` including GTPSA `Descriptor`,
spin, and stochasticity.
"""
zero(m::TaylorMap) = typeof(m)(m)

"""
    one(m::TaylorMap)
  
Creates an identity `m` with the same properties as `m`, including GTPSA
`Descriptor`, spin, and stochasticity.
"""
function one(m::TaylorMap)
  out_m = zero(m)
  nv = numvars(m)

  for i in 1:nv
    out_m.x[i][i] = 1
  end

  if !isnothing(m.q)
    out_m.q.q0[0] = 1
  end

  return out_m
end
