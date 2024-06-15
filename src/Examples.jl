module Examples
using NonlinearNormalForm
using Optim
export track_drift,
       track_fodo,
       track_qd,
       track_qf,
       track_ring,
       track_ring0,
       track_sextupole

const a = 0.00115965218128 
const gamma_0 = 40.5/a
rf_on::Bool = true
#radiation_on::Bool = false
N_fodo::Int = 50

function track_qf(p::Probe, k1, hkick)
  z0 = p.x
  z = Vector{promote_type(eltype(z0),typeof(k1),typeof(hkick))}(undef, length(z0))
  
  lbend=0.1#2*pi/Examples.N_fodo
  L  =  0.5/(1.0+z0[6])
  h  =  -L*(z0[2]^2+k1*z0[1]^2+ z0[4]^2-k1*z0[3]^2)/(1.0+z0[6])/2.0
  z[1] =  cos(sqrt(k1)*L)*z0[1]+1/sqrt(k1)*sin(sqrt(k1)*L)*z0[2]
  z[2] =  -sqrt(k1)*sin(sqrt(k1)*L)*z0[1]+cos(sqrt(k1)*L)*z0[2]+hkick+lbend*z0[6]
  z[3] =  cosh(sqrt(k1)*L)*z0[3]+1/sqrt(k1)*sinh(sqrt(k1)*L)*z0[4]
  z[4] =  sqrt(k1)*sinh(sqrt(k1)*L)*z0[3]+cosh(sqrt(k1)*L)*z0[4]
  z[5] =  z0[5]+h-lbend*z[1]
  z[6] = +z0[6]

  # Spin:
  if !isnothing(p.Q)
    gx =  lbend/L
    sf =  sin(a*lbend*gamma_0/2)
    cf =  cos(a*lbend*gamma_0/2)
    chi = 1+a*gamma_0
    psi = gamma_0^2-1
    zeta = gamma_0-1
    alf =  2*(a^2*gamma_0^2*gx^2+k1)
    bet =  a*gx*k1*(gamma_0*chi-zeta)
    kx =  k1+gx^2
    w_x =  sqrt(kx)
    w_y = sqrt(k1)
    sx =  sin(L*w_x)
    cx =  cos(L*w_x)
    sy =  sinh(L*w_y)
    cy =  cosh(L*w_y)
    sig =  w_y*(k1 + a*k1*gamma_0 + a^2*gx^2*zeta*gamma_0)
    xi =  w_y*(k1*chi + a^2*gx^2*zeta*gamma_0)
    
    q0 =  cf - z0[1]*kx*chi/(2*w_x)*sx*sf - z0[2]*kx*chi/(2*w_x^2)*(1-cx)*sf + z0[6]*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*sf
    q1 =  z0[3]*-1/alf*(bet*(1+cy)*sf + sig*sy*cf) + z0[4]*1/(w_y*alf)*(xi*(1-cy)*cf-bet*sy*sf)
    q2 =  -sf - z0[1]*kx*chi/(2*w_x)*sx*cf - z0[2]*kx*chi/(2*w_x^2)*(1-cx)*cf + z0[6]*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*cf
    q3 =  z0[3]/alf*(bet*(1-cy)*cf + sig*sy*sf) + z0[4]/(w_y*alf)*(xi*(1+cy)*sf-bet*sy*cf)

    q = [q0,q1,q2,q3]
    q = q/sqrt(dot(q,q))
    Q = Quaternion{promote_type(eltype(z0),typeof(k1),typeof(hkick))}(q)
    mul!(Q,Q,p.Q)
    Qkick = Quaternion{promote_type(eltype(z0),typeof(k1),typeof(hkick))}([cos(hkick*(1+a*gamma_0)/2),0, sin(hkick*(1+a*gamma_0)/2), 0])
    mul!(Q, Qkick, Q)
  else 
    Q = nothing
  end

  if !isnothing(p.E)
    lrad1=0.2
    lrad2=0.2
    lrad3=0.2
    z[1] +=  exp(lrad1*(1.0+z[1]^2))*z[1] 
    z[3] +=  exp(lrad2*(1.0+z[3]^2))*z[3] 
    z[5] +=  exp(lrad3*(1.0+z[5]^2))*z[5] 
  end

  return Probe(z,x0=p.x0,Q=Q)
end

function track_qd(p::Probe, k1, vkick)
  z0 = p.x
  z = Vector{promote_type(eltype(z0),typeof(k1),typeof(vkick))}(undef, length(z0))

  L  =  0.5/(1.0+z0[6])
  h  =  -L*(z0[2]^2-k1*z0[1]^2+z0[4]^2+k1*z0[3]^2)/(1.0+z0[6])/2.0
  z[1] =  cosh(sqrt(k1)*L)*z0[1]+1/sqrt(k1)*sinh(sqrt(k1)*L)*z0[2]
  z[2] =  sqrt(k1)*sinh(sqrt(k1)*L)*z0[1]+cosh(sqrt(k1)*L)*z0[2]
  z[3] =  cos(sqrt(k1)*L)*z0[3]+1/sqrt(k1)*sin(sqrt(k1)*L)*z0[4]
  z[4] =  -sqrt(k1)*sin(sqrt(k1)*L)*z0[3]+cos(sqrt(k1)*L)*z0[4]+vkick 
  z[5] =  z0[5]+h
  z[6] = +z0[6]

  # Spin:
  if !isnothing(p.Q)
    gx = 0.
    sf = 0.
    cf = 1.0
    chi = 1+a*gamma_0
    psi = gamma_0^2-1
    zeta = gamma_0-1
    alf = 2*k1
    bet = 0
    kx = k1
    w_x =  sqrt(kx)
    w_y = w_x
    sx =  sinh(L*w_x)
    cx =  cosh(L*w_x)
    sy =  sin(L*w_y)
    cy =  cos(L*w_y)
    sig =  w_y*k1*chi
    xi = sig
    
    q0 =  cf - z0[1]*kx*chi/(2*w_x)*sx*sf + z0[2]*kx*chi/(2*w_x^2)*(1-cx)*sf + z0[6]*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*sf
    q1 =  z0[3]*-1/alf*(bet*(1+cy)*sf - sig*sy*cf) + z0[4]*1/(w_y*alf)*(xi*(1-cy)*cf-bet*sy*sf)
    q2 =  -sf - z0[1]*kx*chi/(2*w_x)*sx*cf + z0[2]*kx*chi/(2*w_x^2)*(1-cx)*cf + z0[6]*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*cf
    q3 =  z0[3]/alf*(bet*(1-cy)*cf - sig*sy*sf) + z0[4]/(w_y*alf)*(xi*(1+cy)*sf-bet*sy*cf)

    q = [q0,q1,q2,q3]
    q = q/sqrt(dot(q,q))
    Q = Quaternion{promote_type(eltype(z0),typeof(k1),typeof(vkick))}(q)
    mul!(Q, Q, p.Q)
    Qkick = Quaternion{promote_type(eltype(z0),typeof(k1),typeof(vkick))}([cos(vkick*(1+a*gamma_0)/2), sin(vkick*(1+a*gamma_0)/2), 0, 0])
    mul!(Q, Qkick, Q)
  else 
    Q = nothing
  end

  return Probe(z,x0=p.x0,Q=Q)
end

function track_drift(p::Probe)
  z0 = p.x
  z = Vector{eltype(z0)}(undef, length(z0))

  L = 0.75
  z[1] =  z0[1]+z0[2]*L/(1.0+z0[6])
  z[2] = +z0[2]
  z[3] =  z0[3]+z0[4]*L/(1.0+z0[6])
  z[4] = +z0[4]
  z[5] =  z0[5]-L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  z[6] = +z0[6] 
  
  return Probe(z,x0=p.x0,Q=p.Q)
end


function track_cav(p::Probe)
  z0 = p.x
  z = Vector{eltype(z0)}(undef, length(z0))

  z[1] = +z0[1]
  z[2] = +z0[2]
  z[3] = +z0[3]
  z[4] = +z0[4]
  z[5] = +z0[5]
  z[6] =  z0[6]+0.0001*z0[5]

  return Probe(z,x0=p.x0,Q=p.Q)
end

function track_sextupole(p::Probe, k2l)
  z0 = p.x
  z = Vector{promote_type(eltype(z0),typeof(k2l))}(undef, length(z0))
  
  z[1] = +z0[1]
  z[2] =  z0[2]-k2l*(z0[1]^2-z0[3]^2)
  z[3] = +z0[3]
  z[4] =  z0[4]+k2l*2.0*z0[1]*z0[3]
  z[5] = +z0[5]
  z[6] = +z0[6]  
  # Spin:
  if !isnothing(p.Q)
    q = [1,-k2l*2.0*z0[1]*z0[3]*(1+a*gamma_0)/2,-k2l*(z0[1]^2-z0[3]^2)*(1+a*gamma_0)/2,0]
    q = q/sqrt(dot(q,q))
    Q = Quaternion{promote_type(eltype(z0),typeof(k2l))}(q)
    mul!(Q, Q, p.Q)
  else 
    Q = nothing
  end

  return Probe(z,x0=p.x0,Q=Q)
end

function track_fodo(p0, k1, k2l, hkick, vkick)
  p1 = track_qf(p0, k1, hkick)
  p2 = track_sextupole(p1, k2l)
  p3 = track_drift(p2)
  p4 = track_qd(p3, k1, vkick)
  p5 = track_sextupole(p4, -k2l)
  p6 = track_drift(p5)
  return p6
end

function track_ring(p0; N_fodo=50, k1=0.36, k2l=1.2, hkicks=zeros(N_fodo), vkicks=zeros(N_fodo))
  for i=1:N_fodo
    p0 = track_fodo(p0, k1, k2l, hkicks[i], vkicks[i])
  end
  if Examples.rf_on
    p0 = track_cav(p0)
  end
  return p0
end

function track_ring0()
  vkicks = [0, zeros(49)...]  # first vertical coil has strength 1e-4 

  # Calculate closed orbit:
  orbit(z) = norm(track_ring(Probe(z), vkicks=vkicks).x - z)
  x0 = optimize(orbit, zeros(6),g_tol=1e-25).minimizer 

  # Now expand along closed orbit
  d = Descriptor(6,3,2,3)
  x = vars()
  k = params()

  p = Probe(x+x0, x0=x0)
  p = track_ring(p,k1=0.36+k[2], vkicks=[vkicks[1]+k[1], 0, zeros(TPS,48)...]) # first and second coil are knobs

  # Make DAMap
  m1 = DAMap(p)
  #m = m1^-1∘m1 
  #println(m)
  return m1

  # println(m1)

  # DAMap concatenation and inversion:
  # m = m1^-1∘m1  # or inv(m1)∘m1
  # println(m)

  # TPSAMap concatenation and inversion:
  # mt1 = TPSAMap(m1)
  # mt = mt1^-1∘mt1  # or inv(mt1)∘mt1
end

end