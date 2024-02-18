function z_deriv()
  d = Descriptor(2,4) # nv = 2, no = 4
  x = vars() # Gets vector of variables

  # Creates 2 x_1 x_ 2 ^2+3 x_1x_2 +4
  f = 2*x[1]*x[2]^2 + 3*x[1]*x[2] + 4
  # or f = 2*mono([1, 2]) + 3*mono([1, 1]) + 4               # Indexing by order
  # or f = 2*mono([1=>1, 2=>2]) + 3*mono([1=>1, 2=>1]) + 4   # Indexing by "sparse monomial"

  df = ∂(f,2) #df/dx_2    Creates 4 x_1x_2 + 3 x_1
  # or df = ∂(f, [0, 1])    # Order
  # or df = ∂(f, [2=>1])    # Sparse Monomial

  idf = ∫(df, 2) # int df/dx_2    recreates 2 x_1 x_ 2 ^2+3 x_1x_2 WITHOUT 4

  print(f)
  print(df)
  print(idf)
end

function z_sub_i()
  d = Descriptor(2,4) # nv = 2, no = 4
  x = vars() # Gets vector of variables

  # Creates 2 x_1 x_ 2 ^2 + 3 x_1^3x_2 + 4 x_2^4 + 5
  f = 2*x[1]*x[2]^2 + 3*x[1]^3*x[2] + 4*x[2]^4 + 5
  # or f = 2*mono([1, 2]) + 3*mono([3, 1]) + 4*mono([0, 4]) + 5
  # or f = 2*mono([1=>1, 2=>2]) + 3*mono([1=>3, 2=>1]) + 4*mono([2=>4]) + 5

  g = getord(f, 4) #  Creates 3 x_1^3x_2 + 4 x_2^4
  print(f)
  print(g)
end

function z_cut()
  d = Descriptor(2,4) # nv = 2, no = 4
  x = vars() # Gets vector of variables
  
  # Creates 2.d0 x_1 x_ 2 + 3.d0x_1^2x_2 + 4.d0x_2^4 + 5.d0
  f = 2*x[1]*x[2] + 3*x[1]^2*x[2] + 4*x[2]^4 + 5
  # or f = 2*mono([1, 1]) + 3*mono([2, 1]) + 4*mono([0, 4]) + 5 
  # or f = 2*mono([1=>1, 2=>1) + 3*mono([1=>2, 2=>1]) + 4*mono([2=>4]) + 5
  
  g = cutord(f, 3) # Creates 2.d0 x_1 x_2 + 5.d0
  print(f)
  print(g)
  
  g = cutord(f, 1) # Creates 5.d0
  print(f)
  print(g)
end

function z_sub_j()
  d = Descriptor(2,4) # nv = 2, no = 4
  x = vars()

  f = 2*x[1]*x[2]^2 + 5
  
  # Peeks coefficient 2.0:
  r = f[1,2]        # Indexing by order
  r = f[1=>1, 2=>2] # Indexing by sparse monomial
  println(f)
  println(r)
  
  f = 2*x[1]*x[2]^2 + 4*x[1]^3*x[2] + 5
  
  # Peeks coefficient 4.0:
  r = f[3,1]         # Indexing by order
  r = f[1=>3, 2=>1]  # Indexing by sparse monomial
  println(f)
  println(r)
  
  f = 2*x[1]*x[2]^2 + 4*x[1]*x[2]^3 + 3*x[2] + 5
  
  # Peeks coefficient 4.0
  r = f[1,3]         # Indexing by order
  r = f[1=>1, 2=>3]  # Indexing by sparse monomial
  println(f)
  println(r)
end

function z_var()
  d = Descriptor(2,4)
  x = vars(use=d)
  
  f = 3 + x[1]
  print(f)
  
  r = [5, 6]
  f = r[1] + r[2]*x[2]
  print(f)
  
  # Could also have
  r = [5im, 6im]          # im is imaginary i in julia
  f = r[1] + r[2]*x[2]    # Promotion to ComplexTPS automatically done
  print(f)
end

function z_mono()
  d = Descriptor(2, 4)
  j = [1,2]
  k = [3,1]

  f = 2*mono(j) + 4*mono(k)
  print(f)

  f = 2*mono([1, 2]) + 3*mono([1, 1])
  print(f)
  
  i = 2
  x = vars(use=d)
  f = 2*x[i] + 3*x[i-1]
  # or f = 2*mono(i) + 3*mono(i-1)
  print(f)
end

function z_par()
  d = Descriptor(4, 6)
  j = [1,2,2,1]
  jj = [1,1,0,3]
  n = 2
  k = j[1:n]

  f = 2*mono(j) + 3*mono(jj) + 5
  g =  f[k...,:] # or par(f, k)
  print(f)
  print(g)

  f = 2*mono([1,1,1,1]) + 4*mono([2,1,1]) + 1*mono([1,1,1,2])
  g = par(f, [1,1])

  print(f)
  print(g)

  f = 2*mono([2,0,1,1]) + 3*mono([2,0,1,2]) + 6*mono([0,0,1]) + 5
  g = par(f, [1=>2, 3=>1])
  print(f)
  print(g)
end

function z_parT()
  error("parT not implemented yet.")
end

function z_shift()
  println("<= not implemented yet.")
end

function z_pseudo_deriv()
  println("k not implemented yet.")
end

function z_poisson()
  d = Descriptor(2,4) # nv = 2, no = 4
  x = vars()[1] # or mono(1)
  p = vars()[2] # or mono(2)
  
  time = 0.01 # seconds
  k = 2
  m = 0.01

  # Hamiltonian Lie Operator for a spring
  h = -time*(p^2/m + k*x^2)/2
  hf = getvectorfield(h)
  display(hf)

  map = exppb(hf, [x, p])    # vars is identity map
  
  println("[x_0, p_0] = 1")
  display([x, p])
  display(map)
  pb = pb(x, p)
  display(pb)

  # Preservation of Poisson Bracket after a time of 0.01
  pb = pb(map[1], map[2])
  println("[x_t, p_t] = 1")
  display(pb)
end


function z_cut_sub()
  error("Random taylor map generator not implemented yet")
end

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
  z5=  z0[5] -L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])/2.0
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
