using Revise
using BenchmarkTools: @btime, @benchmark
using NonlinearNormalForm

undef_init(vec::Vector) = Vector{eltype(vec)}(undef, length(vec))

function track_qf(p::Probe, k1, hkick)
  z0 = p.x0 + p.v
  z = undef_init(z0)
  
  lbend=0.1
  L  = @FastGTPSA 0.5/(1.0+z0[6])
  h  = @FastGTPSA -L*(z0[2]^2+k1*z0[1]^2+ z0[4]^2-k1*z0[3]^2)/(1.0+z0[6])/2.0
  z[1] = @FastGTPSA cos(sqrt(k1)*L)*z0[1]+1/sqrt(k1)*sin(sqrt(k1)*L)*z0[2]
  z[2] = @FastGTPSA -sqrt(k1)*sin(sqrt(k1)*L)*z0[1]+cos(sqrt(k1)*L)*z0[2]+hkick+lbend*z0[6]
  z[3] = @FastGTPSA cosh(sqrt(k1)*L)*z0[3]+1/sqrt(k1)*sinh(sqrt(k1)*L)*z0[4]
  z[4] = @FastGTPSA sqrt(k1)*sinh(sqrt(k1)*L)*z0[3]+cosh(sqrt(k1)*L)*z0[4]
  z[5] = @FastGTPSA z0[5]+h-lbend*z[1]
  z[6] = +z0[6]
  return Probe(z)
end

function track_qd(p::Probe, k1, vkick)
  z0 = p.x0 + p.v
  z = undef_init(z0)

  L  = @FastGTPSA 0.5/(1.0+z0[6])
  h  = @FastGTPSA -L*(z0[2]^2-k1*z0[1]^2+z0[4]^2+k1*z0[3]^2)/(1.0+z0[6])/2.0
  z[1] = @FastGTPSA cosh(sqrt(k1)*L)*z0[1]+1/sqrt(k1)*sinh(sqrt(k1)*L)*z0[2]
  z[2] = @FastGTPSA sqrt(k1)*sinh(sqrt(k1)*L)*z0[1]+cosh(sqrt(k1)*L)*z0[2]
  z[3] = @FastGTPSA cos(sqrt(k1)*L)*z0[3]+1/sqrt(k1)*sin(sqrt(k1)*L)*z0[4]
  z[4] = @FastGTPSA -sqrt(k1)*sin(sqrt(k1)*L)*z0[3]+cos(sqrt(k1)*L)*z0[4]+vkick 
  z[5] = @FastGTPSA z0[5]+h
  z[6] = +z0[6]
  return Probe(z)
end

function track_drift(p::Probe)
  z0 = p.x0 + p.v
  z = undef_init(z0)

  L = 0.75
  z[1] = @FastGTPSA z0[1]+z0[2]*L/(1.0+z0[6])
  z[2] = +z0[2]
  z[3] = @FastGTPSA z0[3]+z0[4]*L/(1.0+z0[6])
  z[4] = +z0[4]
  z[5] = @FastGTPSA z0[5]-L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  z[6] = +z0[6] 
  return Probe(z)
end


function track_cav(p::Probe)
  z0 = p.x0 + p.v
  z = undef_init(z0)

  z[1] = +z0[1]
  z[2] = +z0[2]
  z[3] = +z0[3]
  z[4] = +z0[4]
  z[5] = +z0[5]
  z[6] = @FastGTPSA z0[6]+0.0001*z0[5]
  return Probe(z)
end

function track_sextupole(p::Probe, k2l)
  z0 = p.x0 + p.v
  z = undef_init(z0)
  
  z[1] = +z0[1]
  z[2] = @FastGTPSA z0[2]-k2l*(z0[1]^2-z0[3]^2)
  z[3] = +z0[3]
  z[4] = @FastGTPSA z0[4]+k2l*2.0*z0[1]*z0[3]
  z[5] = +z0[5]
  z[6] = +z0[6]
  return Probe(z)
end

function track_fodo(z0, k1, k2l, kick)
  z1 = track_qf(z0, k1, kick)
  z2 = track_sextupole(z1, k2l)
  z3 = track_drift(z2)
  z4 = track_qd(z3, k1, 0)
  z5 = track_sextupole(z4, -k2l)
  z6 = track_drift(z5)
  return z6
end

function track_ring(z0, k1=0.36, k2l=1.2, kick=zeros(50))
  for i=1:2
    z0 = track_fodo(z0, k1, k2l, kick[i])
  end
# cavity on or off
  z0 = track_cav(z0)
  return z0
end


closed_orbit = 1e-3

d = Descriptor(6,2,2,2);
x = vars();
k = params()

x0 = repeat([closed_orbit], 6)

xs = Probe(x+x0)
k1 = 0.36;
k2l = 1.2;
xs = track_ring(xs, k1, k2l, [k[1], k[2], zeros(TPS,48)...])






#[k[1], k[2], zeros(TPS, 50)...])


#=
m = DAMap(x0, vars())

m = DAMap()

xs = vars();
k = params();

xs .+= closed_orbit;




m = DAMap(xs)

print(m)
=#