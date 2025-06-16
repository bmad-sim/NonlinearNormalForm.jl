#=
struct NormalForm{A,V,S}
  a::A    # Normalizing transformation
  eg::V   # 
  Σ::S
end
=#
function normal(m::DAMap, order::Integer=maxord(m); res=nothing, spin_res=nothing)
  nn = ndiffs(m)
  nhv = nhvars(m)
  nv = nvars(m)
  np = nparams(m)
  mo = order
  coast = iscoasting(m)

  # 1: Go to parameter-dependent fixed point to first order ONLY!
  # Higher orders will be taken care of in nonlinear part
  #= Idea is following:

  This may be useful: https://www.statlect.com/matrix-algebra/determinant-of-block-matrix

  We have a linear map M = [Mz Mk; 0 I] where Mz corresponds to orbital part and Mk to parameters
  Note that parameters always identity and do not depend on orbital (hence 0 and I)

  Now we seek a transformation A = [A11 A12; 0 I] such that A^-1*M*A = [Mz 0; 0 I] (to make M no longer dependent on parameters)
  We can rewrite this as M*A = A*[Mz 0; 0 I]  as

  [Mz*A11, Mz*A12 + Mk;]   =    [A11*Mz,  A12;]
  [0,      I          ;]        [0,       I  ;]      

  The transformation must be symplectic with determinant 1 and so A11 = I
  We are left with 

  Mz*A12 + Mk = A12
  ->
  A12 = -(Mz-I)^-1*Mk

  For coasting beam, we have

  M = [Mhz   0    Mhk ]
      [Vthz' 1    Vtk']
      [0     0    I   ]

  We seek a transformation

  A = [I     0   Ahk]
      [Ath   1   Atk]
      [0     0   I  ]

  to put the map into the form

  A^-1*M*A = [Mhz  0   0   ]
             [0    1   Utk']
             [0    0   I   ]

  The reason for this is because there is no "time-like" fixed point with coasting beam. However 
  we do transform to a space where the time-like coordinate does not depend on the variables at all.

  =#
  if np > 0
    M_hvars = jacobian(m, HVARS)  
    M_hparams = jacobian(m, HPARAMS)
    A0_12 = -inv(M_hvars-I)*M_hparams
    a0 = one(m)
    setray!(a0.v, v_matrix=A0_12, v_matrix_offset=nv) # offset to parameters part

    # if we are coasting, then we need to include the effect of the variables on the time-like coordinate
    # to ensure Poisson bracket does not change
    if coast
      nt = nv
      ndpt = nv+1

      # Parameters part:
      M_cvars = jacobian(m, CVARS)
      M_cparams = jacobian(m, CPARAMS)
      coastrow = M_cvars*A0_12 + M_cparams
      setray!(view(a0.v, nt:nt), v_matrix=coastrow, v_matrix_offset=nv)

      # Variables part:
      # Note that here, because we are only doing first-order, we have exp(F) = I + F
      # And so we can capture this just by adding the Poisson bracket. 
      # But in the nonlinear part (where we "factorize" we have to actually exponentiate
      # the Poisson bracket captured in F to handle this properly.)
      for i in 1:Int(nhv/2)
        TI.seti!(a0.v[nt], TI.geti(a0.v[2*i-1], ndpt), 2*i)
        TI.seti!(a0.v[nt], -TI.geti(a0.v[2*i], ndpt), 2*i-1)
      end
    end
  else
    a0 = I
  end

  if mo == 0
    return a0
  end

  # 2: Do the linear normal form exactly
  # We calculate the matrix A1 such that inv(A1)*M0*A1 is a rotation matrix
  # of each oscillating plane
  Mt = jacobiant(m, HVARS) # parameters (+ coast) never included
  F = mat_eigen(Mt, phase_modes=false) # Returns eigenvectors with vⱼ'*S*vⱼ = +im for odd j, and -im for even j
  a1_inv_matrix = real(c_jacobian(m, HVARS)*transpose(F.vectors))

  a1 = one(m)
  setray!(a1.v, v_matrix=inv(a1_inv_matrix))
  
  if mo == 1 
    return a0 ∘ a1
  end


  # Continue:
  a1i = one(m)
  setray!(a1i.v, v_matrix=a1_inv_matrix)

  c = c_map(m)
  ci = ci_map(m)
  m0 = inv(a0) ∘ m ∘ a0
  m1 = a1i ∘ m0 ∘ a1

  # 3: Go into phasor's basis
  m1 = ci ∘ m1 ∘ c

  # ---- Nonlinear -----
  # 4: Nonlinear algorithm
  # order by order
  # check if tune shift and kill
  R_inv = inv(jacobian(m1, HVARS)) # R is diagonal matrix
  ri = one(m1)
  setray!(ri.v, v_matrix=R_inv)

  # eg = tunes POTENTIALLY INCLUDING DAMPING!!
  # Cannot just use eigenvalues here, need to get from linear normalized map
  eg = diag(jacobian(ri, HVARS))

  an = one(m1)
  ker = one(m1) # current kernel (containing only tune shifts and maybe resonance)
  F =  zero(VectorField, m1)  # VectorField of monomials to remove for this particular order 
  Fker =  zero(VectorField, m1)  # VectorField of tune shifts to keep for this particular order

  ords = similar(m1.v, Int) # monomial orders
  tmpmono = zeros(UInt8, nn) # same as ords but GTPSA v1.4.0 only compatible with this (UInt8 and Vector type) right now

  id = one(m1)

  v = Ref{TI.numtype(eltype(m1.v))}() # monomial value 
  for i in 2:mo
    clear!(F)
    clear!(Fker)

    # Here, we have 
    # m1 = ℛexp(K)(ℐ+ϵ²𝒞₂) 
    # get rid of ℛ: below does:
    nonl = m1 ∘ ri
    nonl = nonl ∘ inv(ker) # get rid of kernel
    # compose!(tmps[1], m1, ri, do_stochastic=false)
    # inv!(tmps[2], ker)
    # compose!(tmps[3], tmps[1], tmps[2], do_stochastic=false) # get rid of kernel

    nonl = getord(nonl, i)
    #getord!(tmps[3], tmps[3], i) # Get only the leading order to stay in symplectic group 
    # nonl = tmps[3]
    # now nonl = ϵ²𝒞₂

    # For each variable in the nonlinear map
    for j in 1:nv
      idx = TI.cycle!(nonl.v[j], 0, mono=tmpmono, val=v)
      while idx > 0
        ords .= tmpmono
        # Tune shifts should be left in the map (kernel) because we cannot remove them
        # if there is damping, we technically could remove them
        if !is_tune_shift(j, ords, nhv) && !is_orbital_resonance(j, ords, nhv, res, spin_res) # then remove it
          ords[j] -= 1
          lam = 1
          for k in 1:nhv # ignore coasting plane
            lam *= eg[k]^ords[k]
          end
          TI.setm!(F.v[j], v[]/(1-lam), tmpmono)
        else # cannot remove it - add to kernel (if close to res and res is specified it will nearly blow up)
          TI.setm!(Fker.v[j], v[], tmpmono)
        end
        idx = TI.cycle!(nonl.v[j], idx, mono=tmpmono, val=v)
      end
    end
    
    kert = exp(Fker,id) # I + Fker 
    ker = kert ∘ ker


    ant = exp(F,id) # I + F
    an = an ∘ ant
    m1 = inv(ant) ∘ m1 ∘ ant
    
#=
    exp!(tmps[3], Fker, id, work_maps=(tmps[1],tmps[2]))
    compose!(tmps[2], tmps[3], ker, do_stochastic=false)
    copy!(ker, tmps[2])

    exp!(tmps[3], F, id, work_maps=(tmps[1],tmps[2]))
    compose!(tmps[2], an, tmps[3], do_stochastic=false)
    copy!(an, tmps[2])
    inv!(tmps[2], tmps[3])
    compose!(tmps[1], tmps[2], m1, do_stochastic=false)
    compose!(m1, tmps[1], tmps[3], do_stochastic=false)
  =#
  end

  a = a0 ∘ a1 ∘ c ∘ an ∘ ci

  # Fully nonlinear map in regular basis

  if !isnothing(m.q)
    
    # First get the 0th order quaternion to make 
    # orbital quaternion be solely a rotation around 
    # vertical to 0th order:

    # The choice here is 
    # We do a cross product of n0 and yhat:
    n0 = TI.scalar.(vect(m1.q))
    n0 = n0/norm(n0) # normalize to 1
    alpha = acos(n0[2])

    nrm = sqrt(n0[1]^2+n0[3]^2)
    qr = Quaternion(cos(alpha/2), sin(alpha/2)*n0[3]/nrm, 0, -sin(alpha/2)*n0[1]/nrm)
    # now concatenate m1:
    aspin = one(m)
    setquat!(aspin.q, q=qr)
    m1 = inv(aspin) ∘ m1 ∘ aspin # == Quaternion(scalar.(Qr*m.q*inv(Qr)))
    nu0 = 2*acos(TI.scalar(m1.q.q0))  # closed orbit spin tune ( i guess we could have gotten this earlier?)
    # it is equal to  2*acos(Quaternion(scalar.(m.q)).q0) (i.e. before transforming)
    # Now we start killing the spin. The first step is to start with a map 
    # (identity in orbital because we are done with orbital) that does this zero order rotation
    qr_inv = one(m)
    setquat!(qr_inv.q, q=inv(Quaternion(TI.scalar.(m1.q))))
    # Now store analogous to eg -> egspin
    egspin = SA[cos(nu0)+im*sin(nu0), 1, cos(nu0)-im*sin(nu0)]

    for i in 1:mo
      # get rid of ℛ:
      # linandnonl = m1.q(I)*qr_inv.q
      linandnonl = m1 ∘ qr_inv
      
      # leaves 1st, 2nd, 3rd, etc terms
      
      # linandnonl contains identity quaternion + delta
      # the quaternion can be written as exp(delta) = 1+delta+delta^2 etc
      # We see at the first iteration, linandnonl.q.q0 contains only second order stuff
      # (not first order). i.e. to leading order the scalar part is cleaned up but the vector
      # part is not (first order stuff here). At the next iteration the q0 part will be only third 
      # order and vector part second order, etc. q0 cleans up itself
      
      # vector part:
      # _s = _something bc I'm only partially understanding this 
      n_s = vect(linandnonl.q)

      # Now we go into eigen-operators of spin 
      # so we can identify the terms to kill much more easily
      nr0_s = MVector(n_s[1]-im*n_s[3], n_s[2], n_s[1]+im*n_s[3])
      nr_s = zero.(nr0_s)
      TI.compose!(nr_s, nr0_s, ri.v)
      na = nr0_s
      TI.clear!.(na)

      # now kill the terms
      for j in 1:3
        idx = TI.cycle!(nr_s[j], 0, mono=tmpmono, val=v)
        while idx > 0
          ords .= tmpmono
          # We remove every term in v and z, tune shifts will only be left in y component
          # because of how we defined everything
          # NOTE SPIN RESONANCES WILL HAPPEN WHEN j != 2 SO WE CHECK THAT FIRST!!
          if !is_spin_resonance(j, ords, nhv, res, spin_res) && (j != 2 || !is_tune_shift(j, ords, nhv, true)) # then remove it, note spin components are like hamiltonian
            lam = egspin[j]
            for k in 1:nhv # ignore coasting plane
              lam *= eg[k]^ords[k]
            end
            TI.setm!(na[j], v[]/(1-lam), tmpmono)
          end
          idx = TI.cycle!(nr_s[j], idx, mono=tmpmono, val=v)
        end
      end
      # Now exit the basis and exponentiate
      na = SA[(na[1]+na[3])/2, na[2], im*(na[1]-na[3])/2]
      qnr = one(m)
      setquat!(qnr.q, q=exp(Quaternion(0,na)))

      aspin = aspin ∘ qnr # put in normalizing map
      m1 = inv(qnr) ∘ m1 ∘ qnr # kill the terms in m1
    end
    a = a ∘ c ∘ aspin ∘ ci
  end

  return  a #real(a)
end


# to get dbeta/ddelta, first go to fully nonlinear parameter dependent fixed point
# then calculate lattice functions. Lattice functions will be TPSA and then you 
# can extract dbeta/ddelta
function factorize(a)
  nv = nvars(a)
  nhv = nhvars(a)
  nn = ndiffs(a)
  mo = maxord(a)
  coast = iscoasting(a)

  if !isnothing(a.q)
    as = one(a)
    tmp = inv(a)
    setquat!(tmp.q, q=I)
    TI.compose!(as.q, a.q, tmp.v)
  end

  # Simply get parameter dependent part
  a0 = a ∘ zero(a)

  if coast
    # At the moment, modulations are NOT implemented  in this package.
    # If modulations AND coasting plane is used, then this computation may disagree with FPP
    # I am not sure which is correct nor if it has been tested in FPP... so if someone wants 
    # to implement modulations make sure there is no coasting
    nt = nv
    ndpt = nv+1
    id = one(a)
    vf = VectorField(v=view(a0.v, 1:nv), q=a0.q)
    TI.clear!(vf.v[nt])
    t1 = zero(vf.v[nt])
    t2 = zero(t1)
    tm = zero(t1)
    # set the timelike variable so poisson bracket does not change
    for i in 1:Int(nhv/2)
      TI.clear!(tm)
      TI.deriv!(t1, vf.v[2*i-1], ndpt)
      TI.seti!(tm, 1, 2*i)
      TI.mul!(t2, t1, tm)
      TI.add!(vf.v[nt], vf.v[nt], t2)

      TI.deriv!(t1, vf.v[2*i], ndpt)
      TI.clear!(tm)
      TI.seti!(tm, -1, 2*i-1)
      TI.mul!(t2, t1, tm)
      TI.add!(vf.v[nt], vf.v[nt], t2)
    end
    a0 = exp(vf,exp(-vf,id) ∘ a0 + I)
  else
    add!(a0, a0, I)
  end
  a1t = inv(a0)*a 
  a1 = zero(a1t)
  ords = zeros(UInt8, nn)
  tmp = zero(a1.v[1])

  # Gets the linear part but retains nonlinear parameter dependence
  v = Ref{TI.numtype(eltype(a1t.v))}() # monomial value 
  for i in 1:nhv
    TI.clear!(tmp)
    ords .= 0
    idx = TI.cycle!(a1t.v[i], 0, mono=ords, val=v)
    while idx > -1
      if sum(view(ords, 1:nhv)) == 1
        TI.setm!(tmp, v[], ords)
      end
      idx = TI.cycle!(a1t.v[i], idx, mono=ords, val=v)
    end
    TI.copy!(a1.v[i], tmp)
  end

  # for the coasting part need to remove quadratic orbital part:
  if coast
    TI.clear!(tmp)
    ords .= 0
    nt = nv
    idx = TI.cycle!(a1t.v[nt], 0, mono=ords, val=v)
    while idx > -1
      if sum(view(ords, 1:nhv)) <= 2
        TI.setm!(tmp, v[], ords)
      end
      idx = TI.cycle!(a1t.v[nt], idx, mono=ords, val=v)
    end
    TI.seti!(tmp, 1, nt)
    TI.copy!(a1.v[nt], tmp)
  end

  if !isnothing(a.q)
    setquat!(a1.q, q=I)
  end

  a2 = inv(a1) ∘ a1t

  if !isnothing(a.q)
    setquat!(a2.q, q=I)
  end

  if !isnothing(a.q)
    return (; as=as, a0=a0, a1=a1, a2=a2)
  else
    return (; a0=a0, a1=a1, a2=a2)
  end
end


"""
    is_tune_shift(varidx, ords, nhv)

Checks if the monomial corresponds to a tune shift.

### Input
- `varidx`      -- Current variable index (e.g. 1 is v, 2 is px, etc)
- `ords`        -- Array of monomial index as orders
- `nhv`         -- Number harmonic variables
- `hamiltonian` -- Default is false, if the monomial is in a vector field and not a hamitlonian then this should be false.
"""
function is_tune_shift(varidx, ords, nhv, hamiltonian=false)
  if !hamiltonian
    ords[varidx] -= 1 # have to subtract because using vectorfield and not hamiltonian
  end
  t = 0
 
  # remove it if there are equal powers of hhbar = J

  for k in 1:2:nhv # ignore coasting plane
    t += abs(ords[k]-ords[k+1]) 
  end

  if !hamiltonian
    ords[varidx] += 1 
  end

  if t != 0
    return false # remove it
  else
    return true # keep it in
  end
end

"""
    is_orbital_resonance(varidx, ords, nhv, res)

Checks if the monomial corresponds to a particular resonance 
(and resonance in the same family)

Note that for 3*Q_x = integer, we also have 6*Q_x though where 

lambda*(3*Q_x) = integer where lambda*3 <= MAXORD + 1 according to notes
are in the same family and so these must be left in too. These are the same family

each column of resonance corresponds to one res in the family

Each row is a multiple of previous

So For example for the 1*Q_1 + 2*Q_2 + 3*Q_3 = integer resonance:
m = 
[ 1  2  3  ...;
  2  4  6  ...;
  3  6  9  ...];

"""
function is_orbital_resonance(varidx, ords, nhv, res, spin_res)
  if isnothing(res)
    return false
  end

  ords[varidx] -= 1

  for curresidx in 1:size(res, 2) # for each res in the family
    if !isnothing(spin_res) && spin_res[curresidx] != 0
      return false # spin res not orbital res
    end
    t1 = 0
    t2 = 0
    
    for k in 1:2:nhv # ignore coasting plane
      t1 += abs(ords[k]-ords[k+1]+res[Int((k+1)/2),curresidx]) 
      t2 += abs(ords[k]-ords[k+1]-res[Int((k+1)/2),curresidx]) 
    end
    if t1 == 0 || t2 == 0
      ords[varidx] += 1
      return true # keep it in!
    end
  end
  ords[varidx] += 1
  return false
end

"""

nu_s + j dot Q = n

0*Qx + 1*Qy + 2*Qs = n

IMPORTANT::::
For spin resonances, there is only one resonance
in each resonance family:

m = [0 ;
     1  ]

ms = [1];

In the code, check - of m + ms
"""
function is_spin_resonance(spinidx, ords, nhv, res, spin_res)
  if isnothing(res) && isnothing(spin_res)
    return false
  end

  size(res, 2) == length(spin_res) || error("Number of resonances in spin_res != number of resonances in res")

  for curresidx in 1:size(res, 2) # for each res in the family
    t1 = 0
    t2 = 0
    
    for k in  1:2:nhv # ignore coasting plane
      t1 += abs(ords[k]-ords[k+1]+res[Int((k+1)/2),curresidx]) 
      t2 += abs(ords[k]-ords[k+1]-res[Int((k+1)/2),curresidx]) 
    end

    m = spin_res[curresidx]
    if spinidx == 1      # i
      if m > 0 && t1 == 0
        return true
      elseif m < 0 && t2 == 0
        return true
      end
    elseif spinidx == 3  # k
      if m > 0 && t2 == 0
        return true
      elseif m < 0 && t1 == 0
        return true
      end
    else                 # j vertical 
      if t1 + abs(m) == 0 || t2 + abs(m) == 0
        return true
      end
    end
  end
  return false
end


function equilibrium_moments(m::DAMap, a::DAMap=normal(m,1))
  !isnothing(m.s) || error("Map does not have stochasticity")
  !all(m.s .== 0) || error("No FD fluctuations in map (m.s .== 0)")

  # Moments Σ transform with map like MΣMᵀ + B 
  # e.g. for m2 ∘ m1:  m.s .= M2*m1.s*transpose(M2) + m2.s
  # This is only linear because this is very complex to higher order not necessary
  # For now we do not include parameters but I'd like to add this later
  # To include parameters later, we would:
  # First must go to parameter-dependent fixed point to all orders (obtained from factorization after 
  # normal form) and then can get the a1 matrix around there
  # tracking code would have to give FD matrix as a function of the parameters

  # Let B = m.s (FD part)
  # We want to find Σ such that Σ = MΣMᵀ + B  (fixed point)
  # very easy to do in phasors basis
  # fixed point transformation does nothing (note a0.s = Bₐ₀ = 0 of course and a0.v is identity in variable but not in parameters)
  # When including parameters, fixed point transformation would have to be fully 
  # nonlinear, probably obtained from factorized `a`
  
  # For now because excluding parameters I do not need a fixed point transformation

  # MΣMᵀ + B  = Σ
  # A₁MA₁⁻¹ A₁Σ A₁ᵀ A₁⁻¹ᵀMᵀA₁ᵀ  + A₁BA₁ᵀ = A₁ΣA₁ᵀ
  # We have R =  A₁MA₁⁻¹ , β = A₁BA₁ᵀ , σ = A₁ΣA₁ᵀ so
  # RσRᵀ + β = σ
  # Now go to phasors and CRC⁻¹=Λ=Λᵀ is diagonal - simple equation
  # Let s = CσC⁻¹ and b = CβCᵀ:
  # ΛsΛ + b = s
  # solve simply

  # for m2 ∘ m1:  m.s .= M2*m1.s*transpose(M2) + m2.s
  C = c_jacobian(m, HVARS)
  CI = ci_jacobian(m, HVARS)
  A = jacobian(a, HVARS)
  AI = inv(A)
  Λ = CI * AI * jacobian(m, HVARS) * A * C
  b = (CI*AI)*m.s*transpose(CI*AI)

  Σc = (1 ./(1 .- diag(Λ)*transpose(diag(Λ)))) .* b

  # Now take it out of phasor basis
  return real((A*C)*Σc*transpose(A*C))
end

# making the 12, 34, 56 elements 0 in the normalizing map
# and returns the phase added to do so
# only canonizes linear part, see c_full_canonise for nonlinear
# modify this to return canonizing map

# Returns the rotation map to put a in Courant Snyder form 
# The phase advance can be acquired from this map by atan(r11,r22), etc etc
function fast_canonize(a::DAMap, damping::Bool=!isnothing(a.s); phase=nothing)
  nv = nvars(a)
  nhv = nhvars(a)
  coast = iscoasting(a)

  canonizer = zero(a)

  a_matrix = real.(jacobian(a, VARS))
  

  #phase = zeros(nvars(a))

  for i in 1:Int(nhv/2) # for each harmonic oscillator
    t = sqrt(a_matrix[2*i-1,2*i-1]^2 + a_matrix[2*i-1,2*i]^2)
    cphi = a_matrix[2*i-1,2*i-1]/t
    sphi = a_matrix[2*i-1,2*i]/t
    if sphi*a_matrix[2*i-1,2*i] + cphi*a_matrix[2*i-1,2*i-1] < 0
      cphi = -cphi
      sphi = -sphi
    end

    TI.seti!(canonizer.v[2*i-1],  cphi, 2*i-1)
    TI.seti!(canonizer.v[2*i],    cphi, 2*i)
    TI.seti!(canonizer.v[2*i-1], -sphi, 2*i)
    TI.seti!(canonizer.v[2*i],    sphi, 2*i-1)

    if !isnothing(phase)
      phase[i] += atan(sphi,cphi)/(2*pi)
    end
  end

  if coast
    nt = nv
    ndpt = nv + 1
    TI.seti!(canonizer.v[nt], 1, nt)
    TI.seti!(canonizer.v[ndpt], 1, ndpt)
    TI.seti!(canonizer.v[nt], -TI.geti(a.v[nt], ndpt), ndpt)
    if !isnothing(phase)
      phase[end] += -TI.geti(a.v[nt], ndpt)
    end
  end
  #a_rot = a_matrix*ri
  #return a_rot

  # Now we have rotated a so that a_12, a_34, a_56, etc are 0 (Courant Snyder)
  # But if we have damping, we also have
  # A*S*transpose(A) != S
  # We can multiply the normalizing map A by some dilation to make it so that, 
  # even though we don't have exactly A*S*transpose(A) == S, that 
  # we atleast have (A*S*transpose(A))[1,2] == 1, (A*S*transpose(A))[2,1] == -1, etc 

  # note that with damping we have M as
  # A*Λ*R*A^-1  where R is the amplitude dependent rotation (diagonal matrix with 
  # complex values on unit circle) and Λ is a diagonal matrix with real values 
  # which correspond to the damping (same in each plane, Diagonal(lambda1, lambda1, lambda2, lambda2, etc)

  if damping
    error("need to finish")
    damp = zeros(Int(nvars(a)/2))
    tmp = zeros(Int(nvars(a)/2), Int(nvars(a)/2))
    for i in 1:Int(nvars(a)/2)
      tmp[i,i] = a_rot[2*i-1,2*i-1]*a_rot[2*i,2*i]-a_rot[2*i-1,2*i]*a_rot[2*i,2*i-1]
      for j in 1:Int(nvars(a)/2)
        if i != j
          tmp[i,j] = a_rot[2*i-1,2*j-1]*a_rot[2*i,2*j]-a_rot[2*i-1,2*j]*a_rot[2*i,2*j-1]
        end
      end
    end
    tmp = inv(tmp)
    damp .= sqrt.(sum(tmp, dims=2))
    a_rot = a_rot*Diagonal(repeat(damp,inner=2))
    damp .= log.(damp)
  end


  return canonizer #a_rot
end

# This can help give you the fixed point
function calc_Hr(m, n, res)
  a = n.a
  c = c_map(m)
  R = inv(c)*inv(a)*m*a*c
  
  # This could be old ----------------- 
  # Now we want to refactorize the map as
  # R = exp(-:p*2*pi/|mr|^2 * mr dot J :)*Rc
  # Eq 5.16 in EYB
  # Then log(Rc) is :Hr: which you can use to calculate the fixed point
  # --------------------------------
  
  # First we have to rotate back 
  # We can use any res in the res family, so just use first
  # need to get tunes around -0.5, 0.5
  mu = imag.(log.(n.evals[1:2:end]))

  s = 0
  mr = res[:,1]
  for i in eachindex(mu)
    if abs(mu[i]) > pi
      sgn = sign(mu[i])
      mu[i] = sgn*-2*pi+mu[i]
    end
    s += mr[i]*mu[i]/(2*pi)
  end
  p = round(s)

  # This is more accurate description -----
  # We first have
  # R = exp(-mu dot J)*Nonlinear
  # We want to split this into components parallel and perpendicular 
  # to plane of res, e.g.

  # -mu dot J = -(mu dot mr)/(mr^2)*mr dot J - (mu dot a)/(a^2)*a dot J  
  # where the first term is the parallel to res and second is
  # perpendicular to res (a's are vectors perpendicular to mr)

  # First need to compute -mu dot J as a VectorField
  # Note 3.51 and 3.51 in EBB are incorrect
  # going from vectorfield to hamiltonian etc are NEGATIVE hamilton's equations
  # but some extra stuff because phasors
  # Usually we have dx/dt = -[H,v] = -dH/dp   and dp/dt = -[H,p] = dH/dx 
  # Now with phasors we have dh/dt = -*H,h* = -i*dH/dhbar  and dhbar/dt = -*H,hbar* = i*dH/dh
  # ( eq 44.65 in special Bmad manual)
  # So per eq 44.67 in special Bmad manual we see the vector field needs a factor of i
  # as opposed to regular

  coef = dot(mr, mu)/norm(mr)^2

  F =  zero(VectorField{typeof(m.v),typeof(m.q)}, use=m) 
  for i in eachindex(mu)
    F.v[2*i-1][2*i-1] =  im*mu[i] - (im*coef*mr[i]) 
    F.v[2*i][2*i] =  -im*mu[i] - (-im*coef*mr[i]) 
  end

  # exp(:F:) is now NEGATIVE the perpendicular part
  # this is negative the first term in Eq 5.18 in EYB
  #eturn F
  N_c = exp(F,R)
  #return N_c

  # now remove co-moving map - this is the epsilon thing
  # negative the second term in Eq. 5.18 in EYB
  clear!(F)
  for i in eachindex(mr)
    F.v[2*i-1][2*i-1] = im*p*mr[i]*2*pi/norm(mr)^2
    F.v[2*i][2*i] = -im*p*mr[i]*2*pi/norm(mr)^2
  end
  #return exp(F,N_c)

  # The left over stuff will now give exp(:Hr:)
  #return exp(F,N_c)
  Hr_vec = log(exp(F,N_c))
  return Hr_vec

end


