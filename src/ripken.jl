function compute_ripken(a1::DAMap{V0}, ::Val{linear}=Val{true}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1)
  if linear
    let a1_mat = jacobian(a1, HVARS), a1i_mat = inv(a1_mat)
      E = StaticArrays.sacollect(SVector{div(nhv, 2),typeof(a1_mat)}, -a1_mat*jp_mat(a1, i)*a1i_mat*j_mat(a1) for i in 1:div(nhv, 2))
      return E
    end
  else
    let
      tmp1 = zero(a1)
      tmp2 = zero(a1)
      a1i = inv(a1)
      setray!(tmp2.v, v_matrix=-j_mat(a1)) # note negative sign is here
      E = StaticArrays.sacollect(SVector{div(nhv, 2)}, begin
        setray!(tmp1.v, v_matrix=jp_mat(a1, i))
        a1∘tmp1∘a1i∘tmp2
      end for i in 1:div(nhv, 2)
      )
      return E
    end
  end
end

function compute_ripken(a1::DAMap, ::Val{linear}=Val{true}()) where {linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1) 
  if linear
    let a1_mat = jacobian(a1, VARS_CPARAM), a1i_mat = inv(a1_mat)
      E = [-a1_mat*jp_mat(a1, i)*a1i_mat*j_mat(a1) for i in 1:div(nhv, 2)]
      return (; E=E)
    end
  else
    let
      tmp1 = zero(a1)
      tmp2 = zero(a1)
      a1i = inv(a1)
      setray!(tmp2.v, v_matrix=-j_mat(a1)) # note negative sign is here
      E = [begin
        setray!(tmp1.v, v_matrix=jp_mat(a1, i))
        a1∘tmp1∘a1i∘tmp2
      end for i in 1:div(nhv, 2)
      ]
      return (; E=E)
    end
  end
end

# Linear:
function ripken_to_A(E::SVector{nd,<:SMatrix{nhv}}; symplectic_tol=1e-8) where {nd,nhv}
  A_blocks = StaticArrays.sacollect(SVector{nd}, reconstruct_A_block(E[i], i; symplectic_tol=symplectic_tol) for i in 1:nd)
  A = hcat(A_blocks...)
  if norm(A*j_mat(Val{nhv}())*transpose(A)-j_mat(Val{nhv}())) > symplectic_tol
    error("Unable to reconstruct A from Ripken functions: non-symplectic reconstruction is not implemented yet")
  end
  return A
end

# Nonlinear:
function ripken_to_A(E::SVector{<:DAMap})
  error("Not implemented")
end

# This was done with assistance of Claude
# B are the two columns of A defined by the given E_i
function reconstruct_A_block(E_i::SMatrix{nhv}, block_idx::Int; symplectic_tol=1e-8) where {nhv}
  # "Clean" row has the zero of B_i; next row gives second column.
  r1 = 2 * block_idx - 1      # 1-indexed row with the zero
  r2 = 2 * block_idx

  # --- First column of B_i, read off row r1 of E_i ---
  beta = E_i[r1, r1]
  if beta <= 0
    error("E[$(block_idx)][$(r1),$(r1)] = beta_$(block_idx)$(block_dix) = $(beta) is not positive")
  end

  pivot = sqrt(beta)   # = B_i[r1, 1] > 0  (sign convention)

  # First column of B:
  B1 = StaticArrays.sacollect(SVector{nhv}, E_i[r1, j] / pivot for j in 1:nhv)

  # B[r1, 1] = pivot  (consistent, automatic)
  # B[r1, 2] = 0      (forced by the zero constraint)

  # --- Second column of B_i, read off row r2 of E_i ---
  # E_i[r2, r2] = B[r2,1]^2 + B[r2,2]^2
  c_sq = E_i[r2, r2] - B1[r2]^2
  c_sq > 0|| throw(ArgumentError("E[$(block_idx)] inconsistent: computed c² = $(c_sq) < 0"))
  c_mag = sqrt(max(c_sq, zero(c_sq)))

  # Pick the sign so that B' J B == J_2 (not -J_2).
  best_sign = +1
  best_B    = hcat(B1, zero(B1))
  best_err  = Inf
  for sc in (+1, -1)
    c = sc * c_mag
    B2 = StaticArrays.sacollect(SVector{nhv, eltype(B1)}, begin
      if j == r1
        zero(c)
      elseif j == r2
        c
      elseif c_mag > symplectic_tol
        (E_i[r2, j] - B1[r2, 1] * B1[j, 1]) / c
      end
    end for j in 1:nhv)
    Btrial = hcat(B1, B2)
    err = maximum(abs.(Btrial' * j_mat(Val{nhv}()) * Btrial - j_mat(Val{2}())))
    if err < best_err
      best_err  = err
      best_sign = sc
      best_B    = Btrial
    end
  end
  return best_B
end