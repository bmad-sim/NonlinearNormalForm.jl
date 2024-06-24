function normal(m::DAMap)
  nn = numnn(m)
  
  # 1: Go to parameter-dependent fixed point to first order ONLY!
  # Higher orders will be taken care of in nonlinear part
 # zero(m) is zero in variables but identity in parameters
  if !isnothing(m.idpt) # if coasting, set number of variables executing pseudo-harmonic oscillations
    nhv = numvars(m)-2
    eye = DAMap(I(nhv),use=m,idpt=m.idpt)
    ndpt = numvars(m)-1+m.idpt
    sgn = 1-2*m.idpt
    nt = ndpt+sgn
    zer = zero(m); zer.x[nt][nt]=1; zer.x[ndpt][ndpt] = 1
    a0 = (cutord(m,2)-eye)^-1 * zer + eye
    println(a0.x[5])
    for i=1:Int(nhv/2)
      a0.x[nt] += sgn*a0.x[2*i][ndpt]*mono(2*i-1,use=getdesc(m)) - sgn*a0.x[2*i-1][ndpt]*mono(2*i,use=getdesc(m))
    end
  else
    nhv = numvars(m)
    a0 = (cutord(m,2)-I)^-1 ∘ zero(m) + I 
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
    R_inv.Q.q[1][0] = 1
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
    nonl = nonl ∘ inv(ker) # get rid of kernal
    nonl = getord(nonl, i)  # Get only the leading order to stay in symplectic group
    # now nonl = ϵ²𝒞₂

    F =  zero(VectorField{eltype(m1.x),typeof(m1.Q)}, use=m1)  # temporary to later exponentiate
    Fker =  zero(VectorField{eltype(m1.x),typeof(m1.Q)}, use=m1)  # temporary to later exponentiate

    # For each variable in the nonlinear map
    for j=1:numvars(m)
      v = Ref{ComplexF64}()
      ords = Vector{UInt8}(undef, nn)
      idx = GTPSA.cycle!(nonl.x[j].tpsa, Cint(j), nn, ords, v)
      while idx > 0
        # Tune shifts should be left in the map (kernel) because we cannot remove them
        # if there is damping, we technically could remove them
        je = convert(Vector{Int}, ords)
        je[j] -= 1
        t = 0
        for k = 1:2:nhv
          t += abs(je[k]-je[k+1])
        end

        if t != 0  # then remove it
          lam = 1
          for k = 1:nhv # ignore coasting plane
            lam *= eg[k]^je[k]
          end
          F.x[j] += mono(ords,use=getdesc(m1))*v[]/(1-lam)
        else # cannot remove it - add to kernel
          je[j] += 1
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
  
  return a0 ∘ a1 ∘ an
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
    cinv.Q.q[1][0] = 1
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
    c.Q.q[1][0] = 1
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