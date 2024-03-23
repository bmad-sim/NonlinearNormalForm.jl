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
  m = cut(xy, order+1)
  return (m-I)^-1∘zero(m)+I
end


function linear_a(xy::DAMap)
  fm0 = transpose(jacobian(xy)) # OPTIMIZE IN FUTURE
  eigen(fm0)

end