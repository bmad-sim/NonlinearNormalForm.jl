function normal(m::DAMap)
  
  work_map = zero(m)
  comp_work_low, inv_work_low = prep_comp_inv_work_low(m)

  a0 = zero(m)

  # 1: Go to the parameter dependent fixed point
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)

  # Similarity transformation to make parameter part of jacobian = 0 (inv(a0)∘m∘a0)
  #m1 = zero(m)
  #compose!(work_map, m, a0, work_low=comp_work_low)
  #inv!(a0, a0, work_low=inv_work_low)
  #compose!(m1, a0, work_map, work_low=comp_work_low)
  return a0
  #a0 = gofix(m)
  #m = inv(a0)∘m∘a0 # similarity transformation to make parameter part of jacobian 0


  # 2: do the linear normal form exactly 

end

function gofix(m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  a0 = zero(m)
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  return a0
end


function gofix!(a0::DAMap, m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  desc = getdesc(m)
  nv = numvars(desc)
  #clear!(a0)

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
  cut!(a0,a0,order+1,dospin=false)

  # 3: map is inverted at least to order 1:
  inv!(a0,a0,work_low=inv_work_low,dospin=false)

  # 4: a map x is created with dimension nv
  # x is zero except for the parameters and delta if coasting
  compose!(a0,a0,work_map,work_low=comp_work_low,dospin=false)

  # 5: add back in identity
  for i=1:nv
    @inbounds a0.x0[i] = m.x0[i]
    @inbounds a0.x[i][i] += 1
  end

  a0.Q.q[1][0] = 1

  return a0
end

function linear_a(xy::DAMap)
  fm0 = transpose(jacobian(xy)) # OPTIMIZE IN FUTURE
  eigen(fm0)

end