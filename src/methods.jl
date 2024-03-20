"""

Finds the parameter dependent fixed point to the specified order.
"""
function gofix(xy::DAMap, order=1)
  desc = getdesc(xy)
  nv = numvars(desc)

  # 1: v = map-identity in harmonic planes
  v = DAMap(xy)
  v.x[1:nv] .-= vars(desc)
  #println(read_fpp_map("gofix1.map").x - v.x)

  # 2: map is cut to order 2 or above
  cut!(v,v,order+1)
  print(read_fpp_map("gofix2.map").x - v.x)

  # 3: map is inverted at least to order 1:
  w = inv(v)

  # 4: a map x is created with dimension nv
  # x is zero except for the parameters and delta if coasting
  x = zero(w) # default identity in parameters
  x.Q.q[1] = 1 # identity in spin
  a1 = w∘x

  # 5: add back in idenrttiy
  a1.x[1:nv] .+= vars(desc)

  return a1
end