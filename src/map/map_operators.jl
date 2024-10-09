#=

Compose and inv using arithmetic operators and 
with uniform scaling, and map powers.

=#


# --- compose ---
for t = (:DAMap, :TPSAMap)
@eval begin
"""
    ∘(m2::$($t), m1::$($t)) -> $($t)

$($t) composition, $( $t == DAMap ? "ignoring" : "including") the scalar part of `m1`
"""
∘(m2::$t, m1::$t) = compose(m2, m1)
# When composing a TPS scalar/vector function w a map, use orbital part of map:
∘(m2::Union{TPS, AbstractVector{<:TPS}}, m1::$t) = GTPSA.compose(m2, m1.x)

literal_pow(::typeof(^), m::$t{S,T,U,V}, vn::Val{n}) where {S,T,U,V,n} = ^(m,n)

# Also allow * for simpliticty and \ and / because why not
*(m2::$t, m1::$t) = ∘(m2, m1)
*(m2::Union{TPS,AbstractVector{<:TPS{<:Union{Float64,ComplexF64}}}}, m1::$t) = GTPSA.compose(m2, m1.x)
/(m2::$t, m1::$t) = m2 ∘ inv(m1) 
\(m2::$t, m1::$t) = inv(m2) ∘ m1

# Uniform scaling for * (∘) and /, \
∘(m::$t, J::UniformScaling) = $t(m)
∘(J::UniformScaling, m::$t) = $t(m)
*(m::$t, J::UniformScaling) = $t(m)
*(J::UniformScaling, m::$t) = $t(m)
/(m::$t, J::UniformScaling) = $t(m)
/(J::UniformScaling, m::$t) = inv(m)
\(m::$t, J::UniformScaling) = inv(m)
\(J::UniformScaling, m::$t) = $t(m)

compose!(m::$t, m1::$t, J::UniformScaling) = copy!(m, m1)
compose!(m::$t, J::UniformScaling, m1::$t) = copy!(m, m1)

end
end



# --- map powers/pow! ---
function pow!(m::DAMap, m1::DAMap, n::Integer)
  checkinplace(m, m1)
  nv = numvars(m1)
  #prepare work
  ref = prep_work_ref(m1)
  
  copy!(m, m1)

  # Store the reference
  map!(t->t[0], ref, view(m.x,1:nv))
  for i=1:(abs(n)-1)
    compose!(m,m,m1, keep_scalar=false)
  end
  
  # Reset the reference (thrown away in compose!)
  for i=1:nv
      m1.x[i][0] = ref[i]
  end

  # Invert and use the container
  if n < 0
    inv!(m,m,work_ref=ref)
  end

  return m
end

function pow!(m::TPSAMap, m1::TPSAMap, n::Integer)
  checkinplace(m, m1)
  copy!(m, m1)
  for i=1:(abs(n)-1)
    compose!(m,m,m1)
  end
  if n < 0
    inv!(m, m)
  end
  return m
end

function ^(m1::TaylorMap, n::Integer)
  m = zero(m1)
  pow!(m, m1, n)
  return m
end
