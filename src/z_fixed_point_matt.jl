using Revise
using NonlinearNormalForm

function track_qf(z0, k1, hkick)

 lbend=0.1



  L = 0.5/(1.0+z0[6])
    h=-L*(z0[2]^2+k1*z0[1]^2+ z0[4]^2-k1*z0[3]^2)/(1.0+z0[6])/2.0
  z1 = @FastGTPSA cos(sqrt(k1)*L)*z0[1]           + 1. /sqrt(k1)*sin(sqrt(k1)*L)*z0[2]
  z2 = @FastGTPSA -sqrt(k1)*sin(sqrt(k1)*L)*z0[1] + cos(sqrt(k1)*L)*z0[2] + hkick 
  z3 = @FastGTPSA cosh(sqrt(k1)*L)*z0[3]          + 1. /sqrt(k1)*sinh(sqrt(k1)*L)*z0[4]
  z4 = @FastGTPSA sqrt(k1)*sinh(sqrt(k1)*L)*z0[3] + cosh(sqrt(k1)*L)*z0[4]
  z5 = @FastGTPSA z0[5]+h
  z6 = @FastGTPSA z0[6]
  z2 = z2 +lbend*z6
  z5 = z5 -lbend*z1

  return [z1,z2,z3,z4,z5,z6]
end

function track_qd(z0, k1, vkick)
  L = 0.5/(1.0+z0[6])
    h=-L*(z0[2]^2-k1*z0[1]^2+ z0[4]^2+k1*z0[3]^2)/(1.0+z0[6])/2.0
  z1 = @FastGTPSA cosh(sqrt(k1)*L)*z0[1]          + 1. /sqrt(k1)*sinh(sqrt(k1)*L)*z0[2]
  z2 = @FastGTPSA sqrt(k1)*sinh(sqrt(k1)*L)*z0[1] + cosh(sqrt(k1)*L)*z0[2]
  z3 = @FastGTPSA cos(sqrt(k1)*L)*z0[3]           + 1. /sqrt(k1)*sin(sqrt(k1)*L)*z0[4]
  z4 = @FastGTPSA -sqrt(k1)*sin(sqrt(k1)*L)*z0[3] + cos(sqrt(k1)*L)*z0[4] + vkick 
  z5 = @FastGTPSA z0[5]+h
  z6 = @FastGTPSA z0[6]
  return [z1,z2,z3,z4,z5,z6] 
end

function track_drift(z0)
  L = 0.75
  z1 = @FastGTPSA z0[1]+z0[2]*L/(1.0+z0[6])
  z3 = @FastGTPSA z0[3]+z0[4]*L/(1.0+z0[6])
  z5=  z0[5] -L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  z6 =  z0[6] 
  return [z1,z0[2],z3 , z0[4],z5,z6]
end


function track_cav(z0)
  z6 =  z0[6] + 0.0001*z0[5]
  return [z0[1],z0[2],z0[3] , z0[4],z0[5],z6]
end

function track_sextupole(z0, k2l)
  z2 = @FastGTPSA z0[2]-k2l*(z0[1]^2 - z0[3]^2)
  z4 = @FastGTPSA z0[4]+k2l/2.0*z0[1]*z0[3]
  return  [z0[1], z2, z0[3], z4, z0[5], z0[6]]
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
    z0=track_cav(z0)
  return z0
end


closed_orbit = 1e-3

d = Descriptor(6,2,2,2);
xs = Probe(repeat([closed_orbit], 6))

x0 = repeat([closed_orbit], 6);
m = DAMap(x0, vars())

m = DAMap()

xs = vars();
k = params();

xs .+= closed_orbit;
k1 = 0.36;
k2l = 1.2;

xs = track_ring(xs, k1, k2l, [k[1], k[2], zeros(TPS, 50)...])

m = DAMap(xs)

print(m)