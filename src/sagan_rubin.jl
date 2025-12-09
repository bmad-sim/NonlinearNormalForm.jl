function compute_sagan_rubin(a1::DAMap{V0}, ::Val{linear}=Val{false}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix

  # If we are NOT coasting, then we will return the H[3][:,5] and H[3][:,6]
  # column vectors (eta and zeta in Etienne's blue book). These are all that 
  # is needed for expressions correct to leading order in the coupling. If 
  # that is not enough, then the de Moivre lattice functions should be used.
  # Note that in such an approximation the longitudinal tune will be the 
  # phase slip

  # If we are coasting, then we will not return any dispersion stuff, since 
  # that is all well-defined in the parameter-dependent fixed point.
  nv = nvars(a1)
  nhv = nhvars(a1)
  coast = iscoasting(a1)
  if linear
    if !coast
      let a1_mat = jacobian(a1, HVARS), a1i_mat = inv(a1_mat), nvm = nv-2
        H = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nv/2))
        gamma_c = sqrt(H[1][1,1])
        C = nv > 2 ? StaticArrays.sacollect(SMatrix{2,2,eltype(H[1])}, H[1][i,j] for i in 3:4 for j in 3:4) : 0
        Ct = nv > 2 ? SA[C[2,2] -C[1,2]; -C[2,1] C[1,1]] : 0
        Vi = gamma_c*I + vcat(hcat(zero(C), -C), hcat(Ct, zero(Ct)))
        N = Vi*StaticArrays.sacollect(SMatrix{nvm,nvm,eltype(a1_mat)}, a1_mat[i,j] for i in 1:nvm for j in 1:nvm)
        beta = StaticArrays.sacollect(SVector{Int(nvm/2),eltype(N)}, N[i,i]^2 for i in 1:2:Int(nvm))
        alpha = StaticArrays.sacollect(SVector{Int(nvm/2),eltype(N)}, -N[i,i+1]*N[i,i] for i in 1:2:Int(nvm))
        # Note that eta and zeta here are in x,px,y,py, not in a,pa,b,pb
        eta = StaticArrays.sacollect(SVector{nvm,eltype(H[1])}, H[end][i,nv] for i in 1:(nvm))
        zeta = StaticArrays.sacollect(SVector{nvm,eltype(H[1])}, H[end][i,nv-1] for i in 1:(nvm))
        # I don't really understand the utility of this, but for consistency with Sagan-Rubin
        # multiply by Vi to put in "normalized" coordinates
        return (; beta=beta, alpha=alpha, eta=Vi*eta, zeta=Vi*zeta, gamma_c=gamma_c, C=C)
      end
    else
      let a1_mat = jacobian(a1, HVARS), a1i_mat = inv(a1_mat)
        H = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2))
        gamma_c = sqrt(H[1][1,1])
        C = nv > 2 ? StaticArrays.sacollect(SMatrix{2,2,eltype(H[1])}, H[1][i,j] for i in 3:4 for j in 3:4) : 0
        Ct = nv > 2 ? SA[C[2,2] -C[1,2]; -C[2,1] C[1,1]] : 0
        Vi = gamma_c*I + vcat(hcat(zero(C), -C), hcat(Ct, zero(Ct)))
        N = Vi*StaticArrays.sacollect(SMatrix{nhv,nhv,eltype(a1_mat)}, a1_mat[i,j] for i in 1:nhv for j in 1:nhv)
        beta = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, N[i,i]^2 for i in 1:2:nhv)
        alpha = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, -N[i,i+1]*N[i,i] for i in 1:2:nhv)
        return (; beta=beta, alpha=alpha, gamma_c=gamma_c, C=C)
      end
    end
  else
    if !coast
      let nvm = nv-2
        tmp = zero(a1)
        a1i = inv(a1)
        H = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
          setray!(tmp.v, v_matrix=ip_mat(a1, i))
          a1∘tmp∘a1i
        end for i in 1:Int(nhv/2)
        )
        mo = maxord(a1)
        gamma_c = TI.cutord(sqrt(factor_out(H[1].v[1], 1)), mo)
        C = nv > 2 ? StaticArrays.sacollect(SMatrix{2,2,eltype(H[1].v)}, factor_out(H[1].v[i], j) for i in 3:4 for j in 3:4) : 0
        Ct = nv > 2 ? SA[C[2,2] -C[1,2]; -C[2,1] C[1,1]] : 0
        Vi = gamma_c*I + vcat(hcat(zero(C), -C), hcat(Ct, zero(Ct)))
        N = TI.cutord.(Vi*StaticArrays.sacollect(SMatrix{nvm,nvm,eltype(a1.v)}, factor_out(a1.v[i], j) for i in 1:nvm for j in 1:nvm), mo)
        beta = StaticArrays.sacollect(SVector{Int(nvm/2),eltype(N)}, TI.cutord(N[i,i]^2, mo) for i in 1:2:Int(nvm))
        alpha = StaticArrays.sacollect(SVector{Int(nvm/2),eltype(N)}, TI.cutord(-N[i,i+1]*N[i,i], mo) for i in 1:2:Int(nvm))
        # Note that eta and zeta here are in x,px,y,py, not in a,pa,b,pb
        eta = StaticArrays.sacollect(SVector{nvm,eltype(H[1].v)}, factor_out(H[end].v[i], nv) for i in 1:(nvm))
        zeta = StaticArrays.sacollect(SVector{nvm,eltype(H[1].v)}, factor_out(H[end].v[i], nv-1) for i in 1:(nvm))
        # I don't really understand the utility of this, but for consistency with Sagan-Rubin
        # multiply by Vi to put in "normalized" coordinates
        return (; beta=beta, alpha=alpha, eta=Vi*eta, zeta=Vi*zeta, gamma_c=gamma_c, C=C)
      end
    else
      let
        tmp = zero(a1)
        a1i = inv(a1)
        H = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
          setray!(tmp.v, v_matrix=ip_mat(a1, i))
          a1∘tmp∘a1i
        end for i in 1:Int(nhv/2)
        )
        mo = maxord(a1)
        gamma_c = TI.cutord(sqrt(factor_out(H[1].v[1], 1)), mo)
        C = nv > 2 ? StaticArrays.sacollect(SMatrix{2,2,eltype(H[1].v)}, factor_out(H[1].v[i], j) for i in 3:4 for j in 3:4) : 0
        Ct = nv > 2 ? SA[C[2,2] -C[1,2]; -C[2,1] C[1,1]] : 0
        Vi = gamma_c*I + vcat(hcat(zero(C), -C), hcat(Ct, zero(Ct)))
        N = TI.cutord.(Vi*StaticArrays.sacollect(SMatrix{nhv,nhv,eltype(a1.v)}, factor_out(a1.v[i], j) for i in 1:nhv for j in 1:nhv), mo)
        beta = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, TI.cutord(N[i,i]^2, mo) for i in 1:2:Int(nhv))
        alpha = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, TI.cutord(-N[i,i+1]*N[i,i], mo) for i in 1:2:Int(nhv))
        return (; beta=beta, alpha=alpha, gamma_c=gamma_c, C=C)
      end
    end
  end
end
