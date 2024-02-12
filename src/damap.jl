using GTPSA
using ForwardDiff
using BenchmarkTools: @btime, @benchmark

function track_qf(z0, k1, hkick)::typeof(z0)
  L = 0.5
  z1 =  cos(sqrt(k1)*L)*z0[1]           + 1. /sqrt(k1)*sin(sqrt(k1)*L)*z0[2]
  z2 =  -sqrt(k1)*sin(sqrt(k1)*L)*z0[1] + cos(sqrt(k1)*L)*z0[2] + hkick
  z3 =  cosh(sqrt(k1)*L)*z0[3]          + 1. /sqrt(k1)*sinh(sqrt(k1)*L)*z0[4]
  z4 =  sqrt(k1)*sinh(sqrt(k1)*L)*z0[3] + cosh(sqrt(k1)*L)*z0[4]
  return [z1,z2,z3,z4]
end

function track_qd(z0, k1, vkick)::typeof(z0)
  L = 0.5
  z1 =  cosh(sqrt(k1)*L)*z0[1]          + 1. /sqrt(k1)*sinh(sqrt(k1)*L)*z0[2]
  z2 =  sqrt(k1)*sinh(sqrt(k1)*L)*z0[1] + cosh(sqrt(k1)*L)*z0[2]
  z3 =  cos(sqrt(k1)*L)*z0[3]           + 1. /sqrt(k1)*sin(sqrt(k1)*L)*z0[4]
  z4 =  -sqrt(k1)*sin(sqrt(k1)*L)*z0[3] + cos(sqrt(k1)*L)*z0[4] + vkick
  return [z1,z2,z3,z4] 
end

function track_drift(z0)::typeof(z0)
  L = 0.75
  z1 =  z0[1]+z0[2]*L
  z3 =  z0[3]+z0[4]*L
  return [z1,z0[2],z3 , z0[4]]
end

function track_sextupole(z0, k2l)::typeof(z0)
  z2 =  z0[2]-k2l*(z0[1]^2 - z0[3]^2)
  z4 =  z0[4]+k2l/2.0*z0[1]*z0[3]
  return  [z0[1], z2, z0[3], z4]
end

function track_fodo(z0, k1, k2l, kick)::typeof(z0)
  z1 = track_qf(z0, k1, kick)
  z2 = track_sextupole(z1, k2l)
  z3 = track_drift(z2)
  z4 = track_qd(z3, k1, 0)
  z5 = track_sextupole(z4, -k2l)
  z6 = track_drift(z5)
  return z6
end

function track_ring(z0, k1=0.36, k2l=1.2, kick=zeros(50))::typeof(z0)
  for i=1:50
    z0 = track_fodo(z0, k1, k2l, kick[i])
  end
  return z0
end

function benchmark_GTPSA()
  d = Descriptor(4,2,52,2)
  z = vars(d)
  k = params(d)
  map = track_ring([z[1], z[2], z[3], z[4]], 0.36+k[1], k[2], k[3:end])
  return map
end

function benchmark_ForwardDiff()
  m(z) = track_ring([z[1], z[2], z[3], z[4]], 0.36+z[5], z[6], z[7:end])
  j = Array{Float64}(undef,4,56)
  h = Array{Float64}(undef,224,56)
  #c = Array{Float64}(undef,12544,56)
  ForwardDiff.jacobian!(j, m, zeros(56))
  ForwardDiff.jacobian!(h, z->ForwardDiff.jacobian(z->m(z), z), zeros(56))
  #ForwardDiff.jacobian!(c, z->ForwardDiff.jacobian(z->ForwardDiff.jacobian(z->m(z), z), z), zeros(56))
  return j, h #, c
end
