struct LinearDAMap{V0,V,Q,S} <: TaylorMap{V0,V,Q,S}
  v0::V0  # vector with length nvars
  v::V    # ndiffs x ndiffs matrix
  q::Q    # optional quaternion matrix
  s::S    # matrix of second order moments
end

ndiffs(m::LinearDAMap) = size(m.v, 2)
nvars(m::LinearDAMap) = length(m.v0)
nparams(m::LinearDAMap) = ndiffs(m) - nvars(m)
nhvars(m::LinearDAMap) = iseven(nvars(m)) ? nvars(m) : nvars(m)-1 # number of "harmonic" variables
iscoasting(m::LinearDAMap) = !iseven(nvars(m))

function ∘(m2::LinearDAMap, m1::LinearDAMap)
  m.v0 .= m1.v0
  nv = nvars(m)

  # Spin:
  if !isnothing(m.q)
    error("Spin for LinearDAMap not implemented yet; use DAMap")
    if do_spin
      # m.q = m2.q(m1.v)*m1.q
    else

    end
  end 
  
  # Stochastic (fast with StaticArrays)
  if !isnothing(m2.s)
    if do_stochastic
      M2 = jacobian(m2)
      m.s .= M2*m1.s*transpose(M2) + m2.s
    else
      m.s .= 0
    end
  end

  return LinearDAMap(m1.v0, m2*m1, )
end