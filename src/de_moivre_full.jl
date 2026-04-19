function compute_de_moivre_full(a1::DAMap{V0}, ::Val{linear}=Val{true}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1)
  if linear
    let a1_mat = jacobian(a1, HVARS), a1i_mat = inv(a1_mat)
      B = StaticArrays.sacollect(SVector{div(nhv, 2),typeof(a1_mat)}, a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:div(nhv, 2))
      K = StaticArrays.sacollect(SVector{div(nhv, 2),typeof(a1_mat)}, -j_mat(a1)*B[i] for i in 1:div(nhv, 2))
      E = StaticArrays.sacollect(SVector{div(nhv, 2),typeof(a1_mat)}, -B[i]*j_mat(a1) for i in 1:div(nhv, 2))
      H = StaticArrays.sacollect(SVector{div(nhv, 2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:div(nhv, 2))
      return (; H=H, B=B, E=E, K=K)
    end
  else
    let
      tmp = zero(a1)
      a1i = inv(a1)
      B = StaticArrays.sacollect(SVector{div(nhv, 2)}, begin
        setray!(tmp.v, v_matrix=jp_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:div(nhv, 2)
      )
      setray!(tmp.v, v_matrix=-j_mat(a1))
      K = StaticArrays.sacollect(SVector{div(nhv, 2)}, begin
        tmp∘B[i]
      end for i in 1:div(nhv, 2)
      )
      E = StaticArrays.sacollect(SVector{div(nhv, 2)}, begin
        B[i]∘tmp
      end for i in 1:div(nhv, 2)
      )
      H = StaticArrays.sacollect(SVector{div(nhv, 2)}, begin
        setray!(tmp.v, v_matrix=ip_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:div(nhv, 2)
      )
      return (; H=H, B=B, E=E, K=K)
    end
  end
end

function compute_de_moivre_full(a1::DAMap, ::Val{linear}=Val{true}()) where {linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1) 
  if linear
    let a1_mat = jacobian(a1, VARS_CPARAM), a1i_mat = inv(a1_mat)
      B = [a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:div(nhv, 2)]
      K = [-j_mat(a1)*B[i] for i in 1:div(nhv, 2)]
      E = [-B[i]*j_mat(a1) for i in 1:div(nhv, 2)]
      H = [a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:div(nhv, 2)]
      return (; H=H, B=B, E=E, K=K)
    end
  else
    let
      tmp = zero(a1)
      a1i = inv(a1)
      B = [begin
        setray!(tmp.v, v_matrix=jp_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:div(nhv, 2)
      ]
      setray!(tmp.v, v_matrix=-j_mat(a1))
      K = [begin
        tmp∘B[i]
      end for i in 1:div(nhv, 2)
      ]
      E = [begin
        B[i]∘tmp
      end for i in 1:div(nhv, 2)
      ]
      H = [begin
        setray!(tmp.v, v_matrix=ip_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:div(nhv, 2)
      ]
      return (; H=H, B=B, E=E, K=K)
    end
  end
end



