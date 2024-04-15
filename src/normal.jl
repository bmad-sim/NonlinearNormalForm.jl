function testallocs!(m, tmp1, tmp2, comp_work_low, inv_work_low, work_ref)
  a0 = tmp1
  work_map = tmp2

  # 1: Go to the parameter dependent fixed point
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  # Similarity transformation to make parameter part of jacobian = 0 (inv(a0)∘m∘a0)
  compose!(work_map, m, a0, work_low=comp_work_low,keep_scalar=false)
  inv!(a0, a0, work_ref=work_ref, work_low=inv_work_low,dospin=false)
  compose!(a0, a0, work_map, work_low=comp_work_low,keep_scalar=false)
  return
end

function normal(m::DAMap)
  tmp1 = zero(m)
  tmp2 = zero(m)
  tmp3 = zero(m)
  comp_work_low, inv_work_low = prep_comp_inv_work_low(m)
  #work_ref = prep_work_ref(m)

  # 1: Go to the parameter dependent fixed point
  gofix!(tmp1, m, 1, work_map=tmp2, comp_work_low=comp_work_low, inv_work_low=inv_work_low)

  # Similarity transformation to make parameter part of jacobian = 0 (inv(a0)∘m∘a0)
  compose!(tmp2, m, tmp1, work_low=comp_work_low,keep_scalar=false)
  inv!(tmp3, tmp1, work_low=inv_work_low,dospin=false)
  compose!(tmp1, tmp3, tmp2, work_low=comp_work_low,keep_scalar=false)

  # tmp1 = m0 now 

  # 2: do the linear normal form exactly 
  linear_a!(tmp2, tmp1, inverse=true)

  # tmp1 is still m0
  # tmp2 is inv(a1)

  # Now normalize linear map inv(a1)*m0*a1
  compose!(tmp3, tmp2, tmp1, work_low=comp_work_low,keep_scalar=false)
  inv!(tmp1, tmp2, work_low=inv_work_low)
  compose!(tmp2, tmp3, tmp1, work_low=comp_work_low,keep_scalar=false)

  # Now tmp2 is m1
  return tmp2
end

function gofix!(a0::DAMap, m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  @assert !(a0 === m) "Aliasing `a0 === m` is not allowed."
  
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

function linear_a!(a1::DAMap, m0::DAMap; inverse=false)
  @assert !(a1 === m0) "Aliasing `a1 === m0` is not allowed."
  
  # We now get the eigenvectors of the compositional map ℳ (which acts on functions of phase space instead 
  # of phase space) in the linear regime. basically f∘ζ = ℳf where ζ is the linear map (vector). f=f(x,p) could 
  # be a vector or scalar function. Assuming f(x,p) = v₁x + v₂p we see that if (written as a transfer matrix) 
  # ζ = [a b; c d] then ℳf = v₁(ax + bp) + v₂(cx + dp) =  (v₁a+v₂c)x + (v₁b+v₂d)p = v̄₁x + v̄₂p . ℳ acts on 
  # functions so essentially we have [v̄₁, v̄₂] = [a c; b d]*[v₁, v₂] =  Mᵀ * [v₁, v₂] . See Etienne's yellow 
  # book Eq 2.39.

  nhv = numvars(m0) # Number harmonic variables
  nhpl = Int(nhv/2) # Number harmonic variables / 2 = number harmonic planes

  work_matrix=zeros(numtype(m0), numvars(m0), numvars(m0))

  #@assert size(work_matrix) == (nhv, nhv) "Incorrect size for work_matrix: received $(size(work_matrix)), require $(tuple(nhv, nhv))"
  #work_matrix .= 0

  @views jacobiant!(work_matrix, m0.x[1:nhv])  # parameters never included
  F = mat_eigen!(work_matrix, phase_modes=false) # no need to phase modes. just a rotation
  
  for i=1:nhv
    for j=1:nhpl
      work_matrix[2*j-1,i] = sqrt(2)*real(F.vectors[i,2*j-1])
      work_matrix[2*j,i] = sqrt(2)*imag(F.vectors[i,2*j-1])
    end
  end

  #work_matrix = [0 -1; 1 0]

  if !inverse # yes this is intentional... the inverse normalizing linear map requires NO inverse step here, where the non-inverse DOES
    work_matrix = inv(work_matrix) # no in-place inverter in Julia, and using this vs. GTPSA is easier + same speed
  end

  clear!(a1)
  for i=1:nhv
    for j=1:nhv
      @inbounds a1.x[i][j] = work_matrix[i,j]
    end
  end
  return a1
end



#=
# Forward is inv(p)*m*p
# Reverse is p*m*inv(p)
# Default is forward
function simil!(out::DAMap, p::DAMap, m::DAMap; reverse::Bool=false, work_map::DAMap=zero(m), comp_work_low=nothing, inv_work_low=nothing, keep_scalar=true)
  compose!(work_map, m, p, work_low=comp_work_low,keep_scalar=keep_scalar)
  inv!(p, p, work_low=inv_work_low)
  compose!(p, p, work_map, work_low=comp_work_low,keep_scalar=keep_scalar)
  
  if reverse
    compose!(work_map, p, m, work_low=comp_work_low,keep_scalar=false)
    #inv!()
  else

  end
end
=#