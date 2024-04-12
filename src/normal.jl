function normal(m::DAMap)
  tmp1 = zero(m)
  tmp2 = zero(m)
  comp_work_low, inv_work_low = prep_comp_inv_work_low(m)

  a0 = tmp1
  work_map = tmp2

  # 1: Go to the parameter dependent fixed point
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  # Similarity transformation to make parameter part of jacobian = 0 (inv(a0)∘m∘a0)
  compose!(work_map, m, a0, work_low=comp_work_low,keep_scalar=false)
  inv!(a0, a0, work_low=inv_work_low)
  compose!(a0, a0, work_map, work_low=comp_work_low,keep_scalar=false)

  m0 = a0
  # 2: do the linear normal form exactly 
  return linear_a(m0)

end

function gofix!(a0::DAMap, m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  desc = getdesc(m)
  nv = numvars(desc)
  clear!(a0)

  if isnothing(comp_work_low) && isnothing(inv_work_low)
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m)
  elseif isnothing(comp_work_low)
    comp_work_low = prep_comp_work_low(m)
  elseif isnothing(inv_work_low)
    inv_work_low = prep_inv_work(m)
  end

  # 1: v = map-identity in harmonic planes, identity in spin
  for i=1:nv
    @inbounds copy!(a0.x[i], m.x[i])
    @inbounds a0.x[i][i] -= 1
  end

  # 2: map is cut to order 2 or above
  cutord!(work_map,a0,order+1,dospin=false)

  # 3: map is inverted at least to order 1:
  inv!(a0,work_map,work_low=inv_work_low,dospin=false)

  # 4: a map x is created with dimension nv
  clear!(work_map)
  # x is zero except for the parameters and delta if coasting
  compose!(a0,a0,work_map,work_low=comp_work_low,dospin=false,keep_scalar=false)

  # 5: add back in identity
  for i=1:nv
    @inbounds a0.x0[i] = m.x0[i]
    @inbounds a0.x[i][i] += 1
  end
  
  if !isnothing(m.Q)
    a0.Q.q[1][0] = 1
  end

  return a0
end

function gofix(m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  a0 = zero(m)
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  return a0
end

#=
# -1 i inv(p)*m*p
function simil!(p::DAMap, m::DAMap,; inverse::Bool=false, work_map::DAMap=zero(m))
  compose!(work_map, m, p, work_low=comp_work_low,keep_scalar=false)
  inv!(p, p, work_low=inv_work_low)
  compose!(p, p, work_map, work_low=comp_work_low,keep_scalar=false)
end
=#
function linear_a(m0::DAMap)
  # We now get the eigenvectors of the compositional map ℳ (which acts on functions of phase space instead 
  # of phase space) in the linear regime. basically f∘ζ = ℳf where ζ is the linear map (vector). f=f(x,p) could 
  # be a vector or scalar function. Assuming f(x,p) = v₁x + v₂p we see that if (written as a transfer matrix) 
  # ζ = [a b; c d] then ℳf = v₁(ax + bp) + v₂(cx + dp) =  (v₁a+v₂c)x + (v₁b+v₂d)p = v̄₁x + v̄₂p . ℳ acts on 
  # functions so essentially we have [v̄₁, v̄₂] = [a c; b d]*[v₁, v₂] =  Mᵀ * [v₁, v₂] . See Etienne's yellow 
  # book Eq 2.39.

  nhv = numvars(m0)  # Number harmonic variables
  nhpl = Int(nhv/2) # Number harmonic variables / 2 = number harmonic planes
  Mt = zeros(numtype(m0), nhv, nhv) #Matrix{numtype(m0)}(undef, nhv, nhv)
  planes = Vector{Int}(undef, nhpl)
  fm = zeros(numtype(m0), nhv, numnn(m0)) #Matrix{complex(numtype(m0))}(undef, nhv, numnn(m0))

  @views jacobiant!(Mt, m0.x[1:nhv])  # parameters never included
  F = mat_eigen!(Mt)

  return F
end

