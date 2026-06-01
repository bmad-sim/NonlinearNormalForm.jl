function compute_bengtsson(a0, a1, m)
  # The return type will be a Dict with key of type NTuple{nhv,Int}
  # and value of type complex scalar if np == 0 or complex TPS if np > 0.
  # This is type stable bc this is stored in the type.
  nhv = nhvars(m)
  nv = nvars(m)
  np = nparams(m)
  nn = ndiffs(m)

  default = np == 0 ? complex(zero(TI.numtype(first(m.v)))) : complex(zero(first(m.v)))
  valtype = typeof(default)
  h = DefaultOrderedDict{NTuple{nhv,Int},valtype}(default)
  a1t = one(a1)
  setray!(a1t.v, v_matrix=jacobian(a1, HVARS))
  #@show jacobian(a1t, ALL)
  #@show jacobian(a0, HVARS)

  a_cs = a0 ∘ a1t # cutord(a, 2) #cutord(a0 ∘ a1 ∘ a2, 2)
  #cutord(a0, 2) ∘ a1t #a0 ∘ a1 #t
#  cutord(a1, 2)
#  cutord(a0 ∘ a1t, 2)#t #cutord(a1, 2) #t #cutord(a1, 2) #cutord(a0 ∘ a1, 2)
#@show jacobian(a_cs, ALL)
#error()
  #@show jacobian(cutord(a1, 2), ALL)
 # error()
  c = c_map(a0)
  ci = ci_map(a0)

  r =  ci ∘ inv(a_cs) ∘ m ∘ a_cs ∘ c
#=
  # Gets the linear part but retains nonlinear parameter dependence
  v = Ref{TI.numtype(eltype(r.v))}() # monomial value 
  r_lin = zero(r)
  tmp = zero(r.v[1])
  ords = zeros(UInt8, nn)

  for i in 1:nhv
    TI.clear!(tmp)
    ords .= 0
    idx = TI.cycle!(r.v[i], 0, mono=ords, val=v)
    while idx > -1
      if sum(view(ords, 1:nhv)) <= 1
        TI.setm!(tmp, v[], ords)
      end
      idx = TI.cycle!(r.v[i], idx, mono=ords, val=v)
    end
    TI.copy_tps!(r_lin.v[i], tmp)
  end
  =#
#=
  # for the coasting part need to remove quadratic orbital part:
  if nv==5
    TI.clear!(tmp)
    ords .= 0
    nt = nv
    idx = TI.cycle!(r.v[nt], 0, mono=ords, val=v)
    while idx > -1
      if sum(view(ords, 1:nhv)) <= 2
        TI.setm!(tmp, v[], ords)
      end
      idx = TI.cycle!(r.v[nt], idx, mono=ords, val=v)
    end
    TI.seti!(tmp, 1, nt)
    TI.copy_tps!(r_lin.v[nt], tmp)
  end
  =#
  #@show (checksymp(r_lin))
  #@show jacobian(r, ALL)
  #@show jacobian(r_lin, ALL)
  #error()
#@show norm(jacobian(r, ALL) - jacobian(r_lin, ALL))
  #=
  rlin = zero(r)
  mono = zeros(UInt8, nv)
  for i in 1:nv
    for j in 1:nv
      mono[j] = 1
      TI.add!(rlin.v[i], rlin.v[i], r.v[i][[mono...,0]])
      TI.add!(rlin.v[i], rlin.v[i], r.v[i][[mono...,1]])
      mono[j] = 0
      #ast_var_slice!(rlin.v[i], r.v[i], j, nv)
    end
  end
  =#
  #r0, r1, r2 = factorize(r)
  #r = ci ∘ inv(a1) ∘ inv(a0) ∘ m ∘ a0 ∘ a1 ∘ c
  #@show jacobian(r, ALL)
  R_inv = inv(jacobian(r, HVARS)) # R is diagonal matrix
  #@show R_inv
  #@show jacobian(inv(rlin), ALL)
  ri = one(r)
  setray!(ri.v, v_matrix=R_inv)
  # Remove the linear part:
  r_nonl = r ∘ ri

  # Log:
  h_vf = log(r_nonl)

  # now extract RDTs and fill the Dict without elevating actual TPSA order
  tmpmono = zeros(UInt8, nn) # same as ords but GTPSA v1.4.0 only compatible with this (UInt8 and Vector type) right now
  v = Ref{TI.numtype(eltype(h_vf.v))}() # monomial value 

  # in phasors ; h ; = im*S*grad(h) \dot \nabla
  
  # We have h_vf = im*S*grad(h)

  for i in 1:nhv
    # i odd  -> qbar component (a = 2*i-1, bumps slot i+1 (pbar),  sign +im
    # i even -> pbar component (a = 2*i),  bumps slot i-1 (qbar),  sign -im
    if isodd(i)
      bump = i + 1
      sgn  = +im
    else
      bump = i - 1
      sgn  = -im
    end

    idx = TI.cycle!(h_vf.v[i], 0, mono=tmpmono, val=v)
    while idx > 0
      fac = 1
      for j in 1:nhv
        fac += Int(tmpmono[j])
      end
      tmpmono[bump] += 1
      contribution = (sgn * v[]) / fac # (fac * sqrt(2)^fac)
      h_idx = ntuple(k -> Int(tmpmono[k]), Val(nhv))
      if TI.is_tps_type(valtype) isa TI.IsTPSType # parameters included
        if haskey(h, h_idx)
          cur = h[h_idx]
        else
          cur = zero(h_vf.v[1])
          h[h_idx] = cur
        end
        tmpmono[1:nhv] .= 0
        cur_v = TI.getm(cur, tmpmono)
        TI.setm!(cur, cur_v + contribution, tmpmono)
      else # All scalar, just add
        h[h_idx] = get(h, h_idx, zero(contribution)) + contribution
      end
      idx = TI.cycle!(h_vf.v[i], idx, mono=tmpmono, val=v)
    end
  end

  return h
end
