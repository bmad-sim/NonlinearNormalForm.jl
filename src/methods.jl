function normal(m::DAMap)
  # 1: Go to the parameter dependent fixed point
  a0 = gofix(m)
  m = inv(a0)∘m∘a0 # similarity transformation to make parameter part of jacobian 0

  # 2: do the linear normal form exactly



end

"""

Finds the parameter dependent fixed point to the specified order.
Does not currently allow coasting beam.
"""
function gofix(xy::DAMap, order=1)
  desc = getdesc(xy)
  nv = numvars(desc)

  # 1: v = map-identity in harmonic planes, identity in spin
  v = zero(xy)
  for i=1:nv
    @inbounds v.x[i] =  xy.x[i] - mono(i,use=desc)
  end
  v.Q.q[1] = 1

  # 2: map is cut to order 2 or above
  cut!(v,v,order+1)

  # 3: map is inverted at least to order 1:
  inv!(v,v)

  # 4: a map x is created with dimension nv
  # x is zero except for the parameters and delta if coasting
  x = zero(v) # default identity in parameters
  x.Q.q[1] = 1 # identity in spin
  compose!(v,v,x)

  # 5: add back in identity
  for i=1:nv
    @inbounds v.x[i] += mono(i,use=desc)
  end

  return v
end


function linear_a(xy::DAMap)
  fm0 = transpose(jacobian(xy)) # OPTIMIZE IN FUTURE
  eigen(fm0)

end