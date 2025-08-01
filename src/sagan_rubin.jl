function compute_sagan_rubin(a1::DAMap{V0}, ::Val{linear}=Val{false}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix

  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  nhv = isodd(length(V0)) ? length(V0)-1 : length(V0)
  if linear
    let a1_mat = jacobian(a1, VARS_CPARAM), a1i_mat = inv(a1_mat), a1_mat_hv = jacobian(a1, HVARS)
      H = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nv/2))
      gamma_c = sqrt(H[1][1,1])
      C = nv > 2 ? StaticArrays.sacollect(SMatrix{2,2,eltype(H[1])}, H[1][i,j] for i in 3:4 for j in 3:4) : 0
      Ct = nv > 2 ? SA[C[2,2] -C[1,2]; -C[2,1] C[1,1]] : 0
      eta_full = StaticArrays.sacollect(SVector{nhv,eltype(H[1])}, H[end][i,nv] for i in 1:nhv)
      eta = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(eta_full)}, eta_full[i] for i in 1:2:nhv) 
      etap = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(eta_full)}, eta_full[i] for i in 2:2:nhv)
      Vi = gamma_c*I + vcat(hcat(zero(C), -C), hcat(Ct, zero(Ct))) #transpose(SDiagonal(Ct,-C)) #[gamma_c*I -C; Ct gamma_c*I]
      N = Vi*a1_mat_hv
      beta = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, N[i,i]^2 for i in 1:2:Int(nv/2))
      alpha = StaticArrays.sacollect(SVector{Int(nhv/2),eltype(N)}, -N[i+1,i]*N[i,i] for i in 1:2:Int(nv/2))
      return (; beta=beta, alpha=alpha, eta=eta, etap=etap, gamma_c=gamma_c, C=C)
    end
  else
    let
      error("Not implemented yet")
    end
  end
end
