"""

Finds the parameter dependent fixed point to the specified order.
Does not currently allow coasting beam.

(a1-I) = (M-I)^-1*(k)

(M-I)*(a1-I) = k
M*a1 - M - a1 - I = k

M*a1 = k +a1 + M + I

M*a1 = a1+k

"""
function gofix(xy::DAMap, order=1)
  desc = getdesc(xy)
  nv = numvars(desc)

  # 1: v = map-identity in harmonic planes
  v = DAMap(xy)
  for i=1:nv
    @inbounds v.x[i] -= mono(i,use=desc)
  end

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