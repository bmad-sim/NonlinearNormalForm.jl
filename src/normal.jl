struct NormalForm
  a
  evals
end

function normal(m::DAMap, resonance=nothing)
  nn = numnn(m)
  
  # 1: Go to parameter-dependent fixed point to first order ONLY!
  # Higher orders will be taken care of in nonlinear part
 # zero(m) is zero in variables but identity in parameters
  if !isnothing(m.idpt) # if coasting, set number of variables executing pseudo-harmonic oscillations
    nhv = numvars(m)-2
    eye = zero(m)
    setmatrix!(eye, I(nhv))
    #eye = DAMap(I(nhv),use=m,idpt=m.idpt)
    ndpt = numvars(m)-1+m.idpt # energy like variable index
    sgn = 1-2*m.idpt
    nt = ndpt+sgn # timelike variable index
    zer = zero(m); zer.x[nt][nt]=1; zer.x[ndpt][ndpt] = 1
    a0 = (cutord(m,2)-eye)^-1 * zer + eye
    # ensure poisson bracket does not change
    for i=1:Int(nhv/2)
      a0.x[nt] += sgn*a0.x[2*i][ndpt]*mono(2*i-1,use=getdesc(m)) - sgn*a0.x[2*i-1][ndpt]*mono(2*i,use=getdesc(m))
    end
  else
    nhv = numvars(m)
    a0 = inv(cutord(m,2)-I, dospin=false) ∘ zero(m) + I 
  end

  m0 = a0^-1 ∘ m ∘ a0

  # 2: Do the linear normal form exactly
  Mt = GTPSA.jacobiant(m0.x[1:nhv])[1:nhv,:]  # parameters (+ coast) never included
  F = mat_eigen(Mt, phase_modes=false) # Returns eigenvectors with vⱼ'*S*vⱼ = +im for odd j, and -im for even j
  a1_inv_matrix = zeros(numvars(m),numvars(m))
  for i=1:nhv
    for j=1:Int(nhv/2)  # See Eqs. 2.48, 2.49 in EYB
      a1_inv_matrix[2*j-1,i] = sqrt(2)*real(F.vectors[i,2*j-1])  # See Eq. 3.74 in EBB for factor of 2
      a1_inv_matrix[2*j,i] = sqrt(2)*imag(F.vectors[i,2*j-1])  
    end
  end
  if !isnothing(m.idpt)
    a1_inv_matrix[nhv+1,nhv+1] = 1; a1_inv_matrix[nhv+2,nhv+2] = 1;
  end

  a1 = zero(promote_type(eltype(a1_inv_matrix),typeof(m0)),use=m0,idpt=m.idpt)
  setmatrix!(a1, inv(a1_inv_matrix))

  a1_mat = fast_canonize(a1)
  clear!(a1)
  setmatrix!(a1, a1_mat)
  if !isnothing(m.Q)
    a1.Q.q0[0] = 1
  end



  m1 = a1^-1 ∘ m0 ∘ a1

  # 3: Go into phasor's basis
  c = to_phasor(m1)
  m1 = c^-1 ∘ m1 ∘ c


  # ---- Nonlinear -----
  # 4: Nonlinear algorithm
  # order by order
  # check if tune shift and kill
  R_inv = inv(getord(m1, 1, 0, dospin=false), dospin=false) # R is diagonal matrix
  if !isnothing(R_inv.Q)
    R_inv.Q.q0[0] = 1
  end

  # Store the tunes
  eg = Vector{numtype(eltype(R_inv.x))}(undef, numvars(m))
  for i=1:nhv
    eg[i] = R_inv.x[i][i]
  end

  mo = maxord(m)

  an = one(m1)
  ker = one(m1)
  # Kernel - these are tune shifts we leave in the map
  for i = 2:mo
    # Here, we have 
    # m1 = ℛexp(K)(ℐ+ϵ²𝒞₂) 
    # get rid of ℛ:
    nonl = m1 ∘ R_inv
    nonl = nonl ∘ inv(ker) # get rid of kernel
    nonl = getord(nonl, i)  # Get only the leading order to stay in symplectic group
    # now nonl = ϵ²𝒞₂

    F =  zero(VectorField{typeof(m1.x),typeof(m1.Q)}, use=m1)  # temporary to later exponentiate
    Fker =  zero(VectorField{typeof(m1.x),typeof(m1.Q)}, use=m1)  # temporary to later exponentiate

    # For each variable in the nonlinear map
    for j=1:numvars(m)
      v = Ref{ComplexF64}()
      ords = Vector{UInt8}(undef, nn)
      idx = GTPSA.cycle!(nonl.x[j].tpsa, Cint(j), nn, ords, v)
      while idx > 0
        # Tune shifts should be left in the map (kernel) because we cannot remove them
        # if there is damping, we technically could remove them
        if !is_tune_shift(j, ords, nhv) && !is_resonance(j, ords, nhv, resonance) # then remove it
          je = convert(Vector{Int}, ords)
          je[j] -= 1
          lam = 1
          for k = 1:nhv # ignore coasting plane
            lam *= eg[k]^je[k]
          end
          F.x[j] += mono(ords,use=getdesc(m1))*v[]/(1-lam)
        else # cannot remove it - add to kernel (if close to resonance and resonance is specified it will nearly blow up)
          Fker.x[j] += mono(ords,use=getdesc(m1))*v[]
        end

        idx = GTPSA.cycle!(nonl.x[j].tpsa, Cint(idx), nn, ords, v)
      end
    end
    kert = exp(Fker,one(m)) #I + Fker
    ant = exp(F,one(m)) #I + F

    ker = kert ∘ ker
    an = an ∘ ant
    m1 = inv(ant) ∘ m1 ∘ ant
  end

  an = c∘an∘c^-1
  
  return NormalForm(a0 ∘ a1 ∘ an, eg)
end




"""
    is_tune_shift(varidx, ords, nhv)

Checks if the monomial corresponds to a tune shift.

### Input
- `varidx` -- Current variable index (e.g. 1 is x, 2 is px, etc)
- `ords`   -- Array of monomial index as orders
- `nhv`    -- Number harmonic variables
"""
function is_tune_shift(varidx, ords, nhv)
  je = convert(Vector{Int}, ords)
  je[varidx] -= 1
  t = 0
  
  for k = 1:2:nhv # ignore coasting plane
    t += abs(je[k]-je[k+1]) 
  end

  if t != 0
    return false # remove it
  else
    return true # keep it in
  end
end

"""
    is_resonance(varidx, ords, nhv, resonance)

Checks if the monomial corresponds to a particular resonance 
(and resonance in the same family)

Note that for 3*Q_x = integer, we also have 6*Q_x though where 

lambda*(3*Q_x) = integer where lambda*3 <= MAXORD + 1 according to notes
are in the same family and so these must be left in too. These are the same family

each column of resonance corresponds to one resonance in the family

Each row is a multiple of previous

So For example for the 1*Q_1 + 2*Q_2 + 3*Q_3 = integer resonance:
m = 
[ 1  2  3  ...;
  2  4  6  ...;
  3  6  9  ...];

"""
function is_resonance(varidx, ords, nhv, resonance)
  if isnothing(resonance)
    return false
  end
  je = convert(Vector{Int}, ords)
  je[varidx] -= 1

  for curresidx=1:size(resonance, 2) # for each resonance in the family
    t1 = 0
    t2 = 0
    
    for k = 1:2:nhv # ignore coasting plane
      t1 += abs(je[k]-je[k+1]+resonance[Int((k+1)/2),curresidx]) 
      t2 += abs(je[k]-je[k+1]-resonance[Int((k+1)/2),curresidx]) 
    end
    if t1 == 0 || t2 == 0
      return true # keep it in!
    end
  end

  return false
end

function equilibrium_moments(m::DAMap, a::DAMap)
  !isnothing(m.E) || error("Map does not have stochasticity")
  !all(m.E .== 0) || error("No FD fluctuations in map (m.E .== 0)")

  # Moments Σ transform with map like MΣMᵀ + B 
  # This is only linear because this is very complex to higher order not necessary
  # For now we do not include parameters but I'd like to add this later
  # To include parameters later, we would:
  # First must go to parameter-dependent fixed point to all orders (obtained from factorization after 
  # normal form) and then can get the a1 matrix around there
  # tracking code would have to give FD matrix as a function of the parameters

  # Let B = m.E (FD part)
  # We want to find Σ such that Σ = MΣMᵀ + B  (fixed point)
  # very easy to do in phasors basis
  # fixed point transformation does nothing (note a0.E = Bₐ₀ = 0 of course and a0.x is identity in variable but not in parameters)
  # When including parameters, fixed point transformation would have to be fully 
  # nonlinear, probably obtained from factorized `a`
  
  # For now because excluding parameters I do not need a fixed point transformatio

  # MΣMᵀ + B  = Σ
  # A₁MA₁⁻¹ A₁Σ A₁ᵀ A₁⁻¹ᵀMᵀA₁ᵀ  + A₁BA₁ᵀ = A₁ΣA₁ᵀ
  # We have R =  A₁MA₁⁻¹ , β = A₁BA₁ᵀ , σ = A₁ΣA₁ᵀ so
  # RσRᵀ + β = σ
  # Now go to phasors and CRC⁻¹=Λ=Λᵀ is diagonal - simple equation
  # Let s = CσC⁻¹ and b = CβCᵀ:
  # ΛsΛ + b = s
  # solve simply

  c = to_phasor(m)

  R = inv(c)*inv(a)*m*a*c
  b = R.E
  Λ = GTPSA.jacobian(R)

  #return b

  s = zeros(ComplexF64, numvars(m),numvars(m)) # beam envelope in phasors basis
  emits = zeros(Float64, Int(numvars(m)/2)) # emittances

  for i=1:numvars(m)
    for j=1:numvars(m)
      s[i,j] = 1/(1-Λ[i,i]*Λ[j,j])*b[i,j]
    end
  end

  R.E .= s
  return (a*c*R*inv(c)*inv(a)).E

end

# making the 12, 34, 56 elements 0 in the normalizing map
# and returns the phase added to do so
function fast_canonize(a::DAMap, damping::Bool=!isnothing(a.E))
  # Basically rotates a so that we are in Courant-Snyder form of a

  a_matrix = real.(GTPSA.jacobian(a))
  ri = zero(a_matrix)
  #ri .= 0

  if !isnothing(a.idpt) # if coasting, set number of variables executing pseudo-harmonic oscillations
    nhv = numvars(a)-2
  else
    nhv = numvars(a)
  end

  phase = zeros(numvars(a))

  for i=1:Int(nhv/2) # for each harmonic oscillator
    t = sqrt(a_matrix[2*i-1,2*i-1]^2 + a_matrix[2*i-1,2*i]^2)
    cphi = a_matrix[2*i-1,2*i-1]/t
    sphi = a_matrix[2*i-1,2*i]/t
    if sphi*a_matrix[2*i-1,2*i] + cphi*a_matrix[2*i-1,2*i-1] < 0
      cphi = -cphi
      sphi = -sphi
    end

    ri[2*i-1,2*i-1] =  cphi 
    ri[2*i,2*i]     =  cphi 
    ri[2*i-1,2*i]   = -sphi  
    ri[2*i,2*i-1]   =  sphi  

    phase[i] += atan(sphi,cphi)/(2*pi)
  end

  if !isnothing(a.idpt)
    ndpt = numvars(a)-1+a.idpt
    sgn = 1-2*a.idpt
    nt = ndpt+sgn
    ri[nt,nt] = 1
    ri[ndpt,ndpt] = 1
    ri[nt,ndpt] = -a_matrix[nt,ndpt]
    phase[end] += a_matrix[nt,ndpt]
  end
  #return ri
  a_rot = a_matrix*ri
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
    damp = zeros(Int(numvars(a)/2))
    tmp = zeros(Int(numvars(a)/2), Int(numvars(a)/2))
    for i=1:Int(numvars(a)/2)
      tmp[i,i] = a_rot[2*i-1,2*i-1]*a_rot[2*i,2*i]-a_rot[2*i-1,2*i]*a_rot[2*i,2*i-1]
      for j=1:Int(numvars(a)/2)
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


  return a_rot
end

# This can help give you the fixed point
function calc_Hr(m, n, resonance)
  a = n.a
  c = to_phasor(m)
  R = inv(c)*inv(a)*m*a*c
  
  # This could be old ----------------- 
  # Now we want to refactorize the map as
  # R = exp(-:p*2*pi/|mr|^2 * mr dot J :)*Rc
  # Eq 5.16 in EYB
  # Then log(Rc) is :Hr: which you can use to calculate the fixed point
  # --------------------------------
  
  # First we have to rotate back 
  # We can use any resonance in the resonance family, so just use first
  # need to get tunes around -0.5, 0.5
  mu = imag.(log.(n.evals[1:2:end]))

  s = 0
  mr = resonance[:,1]
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
  # to plane of resonance, e.g.

  # -mu dot J = -(mu dot mr)/(mr^2)*mr dot J - (mu dot a)/(a^2)*a dot J  
  # where the first term is the parallel to resonance and second is
  # perpendicular to resonance (a's are vectors perpendicular to mr)

  # First need to compute -mu dot J as a VectorField
  # Note 3.51 and 3.51 in EBB are incorrect
  # going from vectorfield to hamiltonian etc are NEGATIVE hamilton's equations
  # but some extra stuff because phasors
  # Usually we have dx/dt = -[H,x] = -dH/dp   and dp/dt = -[H,p] = dH/dx 
  # Now with phasors we have dh/dt = -*H,h* = -i*dH/dhbar  and dhbar/dt = -*H,hbar* = i*dH/dh
  # ( eq 44.65 in special Bmad manual)
  # So per eq 44.67 in special Bmad manual we see the vector field needs a factor of i
  # as opposed to regular

  coef = dot(mr, mu)/norm(mr)^2

  F =  zero(VectorField{typeof(m.x),typeof(m.Q)}, use=m) 
  for i in eachindex(mu)
    F.x[2*i-1][2*i-1] =  im*mu[i] - (im*coef*mr[i]) 
    F.x[2*i][2*i] =  -im*mu[i] - (-im*coef*mr[i]) 
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
    F.x[2*i-1][2*i-1] = im*p*mr[i]*2*pi/norm(mr)^2
    F.x[2*i][2*i] = -im*p*mr[i]*2*pi/norm(mr)^2
  end
  #return exp(F,N_c)

  # The left over stuff will now give exp(:Hr:)
  #return exp(F,N_c)
  Hr_vec = log(exp(F,N_c))
  return Hr_vec

end
















# 2x faster but not as pretty
function normal_fast(m::DAMap{S,T,U,V}) where {S,T,U,V}
  tmp1 = zero(m)
  tmp2 = zero(m)
  tmp3 = zero(m)
  comp_work_low, inv_work_low = prep_comp_inv_work_low(m)

  if numparams(m) > 0
    # 1: Go to the parameter dependent fixed point ------------------------------------------
    gofix!(tmp1, m, 1, work_map=tmp2, comp_work_low=comp_work_low, inv_work_low=inv_work_low)

    # tmp1 is a0 (linear transformation to parameter-dependent fixed point)

    # Similarity transformation to make parameter part of jacobian = 0 (m0 = inv(a0)∘m∘a0)
    compose!(tmp2, m, tmp1, work_low=comp_work_low,keep_scalar=false)
    inv!(tmp3, tmp1, work_low=inv_work_low,dospin=false)
    compose!(tmp1, tmp3, tmp2, work_low=comp_work_low,keep_scalar=false)

    # tmp1 is m0 (at parameter-dependent fixed point)
  else
    copy!(tmp1, m)
  end
  
  # 2: do the linear normal form exactly --------------------------------------------------
  linear_a!(tmp2, tmp1, inverse=true)

  # tmp1 is still m0
  # tmp2 is inv(a1)
  # Now normalize linear map inv(a1)*m0*a1
  compose!(tmp3, tmp2, tmp1, work_low=comp_work_low,keep_scalar=false)
  inv!(tmp1, tmp2, work_low=inv_work_low,dospin=false)
  compose!(tmp2, tmp3, tmp1, work_low=comp_work_low,keep_scalar=false)

  # tmp2 is m1 , if m is complex then w
  m1 = tmp2
  
  # 3: Go into phasors' basis = c*m1*inv(c) ------------------------------------------------
  # The nonlinear part of the normal form should be in the phasor's basis (see Eq. 3.88 in 
  # Etienne's blue book)
  # We now need complex stuff if we don't already have it
  if eltype(m.x) != ComplexTPS
    ctmp1 = zero(complex(typeof(m)),use=m)
    ctmp2 = zero(ctmp1)
    comp_work_low, inv_work_low = prep_comp_inv_work_low(ctmp1)
  else
    ctmp1 = tmp1
    ctmp2 = tmp3
  end
  
  work_prom = prep_comp_work_prom(ctmp2, ctmp1, m1) # will be nothing if there is no promotion

  to_phasor!(ctmp1,m1) # ctmp1 = c
  compose!(ctmp2, ctmp1, m1, work_low=comp_work_low,keep_scalar=false,work_prom=work_prom)
  from_phasor!(ctmp1, m1) # ctmp1 = inv(c)
  compose!(ctmp2, ctmp2, ctmp1, work_low=comp_work_low,keep_scalar=false)

  # ctmp2 is map in phasors basis
  return ctmp2
  # 4: Nonlinear algorithm --------------------------------------------------------------------
  # Now we have m1, in the phasor's basis. Therefore the nonlinear part of the canonical 
  # transformation a will have a factored Lie representation:
  #
  #         a₂∘a₃∘... = ...exp(F₃)exp(F₂)I
  # 
  # Next step is to get only the nonlinear part of ctmp1
  # this can be done by inverting the linear part of ctmp1 and composing
  mo = getdesc(m).mo
  nv = numvars(m)
  n = Vector{VectorField{T,U}}(undef, mo)  # For now, analog to c_factored_lie, likely will make separate type
  for i=2:mo
    # Get current order + identity
    getord!(ctmp1, ctmp2, i)
    add!(ctmp1,ctmp1,I)

    # get vector field
    G = log(ctmp1)
    n[i] = zero(G)

    # Now kill the terms if not a tune shift!!
    for j=1:nv


    end

  end

  return ctmp2
end

function gofix!(a0::DAMap, m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  checkidpt(a0,m,work_map)
  @assert !(a0 === m) "Aliasing `a0 === m` is not allowed."
  
  desc = getdesc(m)
  nv = numvars(desc)
  clear!(a0)

  if isnothing(comp_work_low) && isnothing(inv_work_low)
    comp_work_low, inv_work_low = prep_comp_inv_work_low(m)
  elseif isnothing(comp_work_low)
    comp_work_low = prep_comp_work_low(m)
  elseif isnothing(inv_work_low)
    inv_work_low = prep_inv_work(m)
  end

  # 1: v = map-identity in harmonic planes, identity in spin
  sub!(a0, m, I, dospin=false)
  if !isnothing(a0.Q)
    a0.Q.q[1][0] = 0
  end

  # 2: map is cut to order 2 or above
  cutord!(work_map,a0,order+1,dospin=false)

  # 3: map is inverted at least to order 1:
  inv!(a0,work_map,work_low=inv_work_low,dospin=false)

  # 4: a map x is created with dimension nv
  clear!(work_map)
  # x is zero except for the parameters and delta if coasting
  compose!(a0,a0,work_map,work_low=comp_work_low,dospin=false,keep_scalar=false)

  # 5: add back in identity
  add!(a0, a0, I, dospin=false)

  a0.x0 .= m.x0

  if !isnothing(m.Q)
    a0.Q.q[1][0] = 1
  end

  return a0
end

function testallocs!(m, tmp1, tmp2, comp_work_low, inv_work_low, work_ref)
  a0 = tmp1
  work_map = tmp2

  # 1: Go to the parameter dependent fixed point
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  # Similarity transformation to make parameter part of jacobian = 0 (inv(a0)∘m∘a0)
  compose!(work_map, m, a0, work_low=comp_work_low,keep_scalar=false)
  inv!(a0, a0, work_ref=work_ref, work_low=inv_work_low,dospin=false)
  compose!(a0, a0, work_map, work_low=comp_work_low,keep_scalar=false)
  return
end

function gofix(m::DAMap, order=1; work_map::DAMap=zero(m), comp_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing, inv_work_low::Union{Nothing,Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}}=nothing)
  a0 = zero(m)
  gofix!(a0, m, 1, work_map=work_map, comp_work_low=comp_work_low, inv_work_low=inv_work_low)
  return a0
end

function linear_a!(a1::DAMap, m0::DAMap; inverse=false)
  checkidpt(a1,m0)
  @assert !(a1 === m0) "Aliasing `a1 === m0` is not allowed."
  
  # We now get the eigenvectors of the compositional map ℳ (which acts on functions of phase space instead 
  # of phase space) in the linear regime. basically f∘ζ = ℳf where ζ is the linear map (vector). f=f(x,p) could 
  # be a vector or scalar function. Assuming f(x,p) = v₁x + v₂p we see that if (written as a transfer matrix) 
  # ζ = [a b; c d] then ℳf = v₁(ax + bp) + v₂(cx + dp) =  (v₁a+v₂c)x + (v₁b+v₂d)p = v̄₁x + v̄₂p . ℳ acts on 
  # functions so essentially we have [v̄₁, v̄₂] = [a c; b d]*[v₁, v₂] =  Mᵀ * [v₁, v₂] . See Etienne's yellow 
  # book Eq 2.39.

  nhv = numvars(m0) # Number harmonic variables
  nhpl = Int(nhv/2) # Number harmonic variables / 2 = number harmonic planes

  work_matrix=zeros(numtype(m0), numvars(m0), numvars(m0))

  @views GTPSA.jacobiant!(work_matrix, m0.x[1:nhv])  # parameters never included
  F = mat_eigen!(work_matrix, phase_modes=false) # no need to phase modes. just a rotation
  
  for i=1:nhv
    for j=1:nhpl
      work_matrix[2*j-1,i] = sqrt(2)*real(F.vectors[i,2*j-1])  # See Eq. 3.74 in EBB for factor of 2
      work_matrix[2*j,i] = sqrt(2)*imag(F.vectors[i,2*j-1])
    end
  end

  if !inverse # yes this is intentional... the inverse normalizing linear map requires NO inverse step here, where the non-inverse DOES
    work_matrix = inv(work_matrix) # no in-place inverter in Julia, and using this vs. GTPSA is easier + same speed
  end

  clear!(a1)
  for i=1:nhv
    for j=1:nhv
      @inbounds a1.x[i][j] = work_matrix[i,j]
    end
  end

  # Make spin identity
  if !isnothing(m0.Q)
    a1.Q.q[1][0] = 1
  end

  return a1
end

function linear_a(m0::DAMap; inverse=false)
  a1 = zero(m0)
  linear_a!(a1, m0, inverse=inverse)
  return a1
end

# computes c^-1
function from_phasor!(cinv::DAMap, m::DAMap)
  checkidpt(cinv,m)
  clear!(cinv)

  nhv = numvars(m)
  if !isnothing(m.idpt)
    nhv -= 2
    cinv.x[nhv+1][nhv+1] = 1
    cinv.x[nhv+2][nhv+2] = 1
  end

  for i=1:Int(nhv/2)
    # x_new = 1/sqrt(2)*(x+im*p)
    cinv.x[2*i-1][2*i-1] = 1/sqrt(2)
    cinv.x[2*i-1][2*i]   = complex(0,1/sqrt(2))

    # p_new = 1/sqrt(2)*(x-im*p)
    cinv.x[2*i][2*i-1] = 1/sqrt(2)
    cinv.x[2*i][2*i]   = complex(0,-1/sqrt(2))
  end

  if !isnothing(cinv.Q)
    cinv.Q.q0[0] = 1
  end

  return
end

function from_phasor(m::DAMap)
  cinv=zero(complex(typeof(m)),use=m,idpt=m.idpt);
  from_phasor!(cinv,m);
  return cinv
end

# computes c
function to_phasor!(c::DAMap, m::DAMap)
  checkidpt(c,m)
  clear!(c)

  nhv = numvars(m)
  if !isnothing(m.idpt)
    nhv -= 2
    c.x[nhv+1][nhv+1] = 1
    c.x[nhv+2][nhv+2] = 1
  end


  for i=1:Int(nhv/2)
    c.x[2*i-1][2*i-1] = 1/sqrt(2)
    c.x[2*i-1][2*i]   = 1/sqrt(2)

    c.x[2*i][2*i-1] = complex(0,-1/sqrt(2))
    c.x[2*i][2*i]   = complex(0,1/sqrt(2))
  end

  if !isnothing(c.Q)
    c.Q.q0[0] = 1
  end

  return
end

function to_phasor(m::DAMap)
  c=zero(complex(typeof(m)),use=m,idpt=m.idpt);
  to_phasor!(c,m);
  return c
end


#=
# Forward is inv(p)*m*p
# Reverse is p*m*inv(p)
# Default is forward
function simil!(out::DAMap, p::DAMap, m::DAMap; reverse::Bool=false, work_map::DAMap=zero(m), comp_work_low=nothing, inv_work_low=nothing, keep_scalar=true)
  compose!(work_map, m, p, work_low=comp_work_low,keep_scalar=keep_scalar)
  inv!(p, p, work_low=inv_work_low)
  compose!(p, p, work_map, work_low=comp_work_low,keep_scalar=keep_scalar)
  
  if reverse
    compose!(work_map, p, m, work_low=comp_work_low,keep_scalar=false)
    #inv!()
  else

  end
end
=#