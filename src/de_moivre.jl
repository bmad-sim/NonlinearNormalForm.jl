function compute_de_moivre(a1::DAMap{V0}, ::Val{linear}=Val{true}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1)
  if linear
    let a1_mat = jacobian(a1, HVARS), a1i_mat = inv(a1_mat)
      B = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2))
      K = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, -j_mat(a1)*B[i] for i in 1:Int(nhv/2))
      E = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, -B[i]*j_mat(a1) for i in 1:Int(nhv/2))
      H = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2))
      return (; H=H, B=B, E=E, K=K)
    end
  else
    let
      tmp = zero(a1)
      a1i = inv(a1)
      B = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
        setray!(tmp.v, v_matrix=jp_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nhv/2)
      )
      setray!(tmp.v, v_matrix=-j_mat(a1))
      K = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
        tmp∘B[i]
      end for i in 1:Int(nhv/2)
      )
      E = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
        B[i]∘tmp
      end for i in 1:Int(nhv/2)
      )
      H = StaticArrays.sacollect(SVector{Int(nhv/2)}, begin
        setray!(tmp.v, v_matrix=ip_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nhv/2)
      )
      return (; H=H, B=B, E=E, K=K)
    end
  end
end

function compute_de_moivre(a1::DAMap, ::Val{linear}=Val{true}()) where {linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nhv = nhvars(a1) 
  if linear
    let a1_mat = jacobian(a1, VARS_CPARAM), a1i_mat = inv(a1_mat)
      B = [a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2)]
      K = [-j_mat(a1)*B[i] for i in 1:Int(nhv/2)]
      E = [-B[i]*j_mat(a1) for i in 1:Int(nhv/2)]
      H = [a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2)]
      return (; H=H, B=B, E=E, K=K)
    end
  else
    let
      tmp = zero(a1)
      a1i = inv(a1)
      B = [begin
        setray!(tmp.v, v_matrix=jp_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nhv/2)
      ]
      setray!(tmp.v, v_matrix=-j_mat(a1))
      K = [begin
        tmp∘B[i]
      end for i in 1:Int(nhv/2)
      ]
      E = [begin
        B[i]∘tmp
      end for i in 1:Int(nhv/2)
      ]
      H = [begin
        setray!(tmp.v, v_matrix=ip_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nhv/2)
      ]
      return (; H=H, B=B, E=E, K=K)
    end
  end
end

# =================================================================================== #
# Construct special matrices for lattice functions
function ip_mat(m::DAMap{V0}, i) where {V0<:StaticArray}
  nhv = nhvars(m) #isodd(length(V0)) ? length(V0)+1 : length(V0)
  return StaticArrays.sacollect(SMatrix{nhv,nhv,real(eltype(V0))}, 
  begin
    if col == 2*i-1 && row == 2*i-1
      1
    elseif col == 2*i && row == 2*i
      1
    else
      0
    end
  end for col in 1:nhv for row in 1:nhv)
end

function ip_mat(m::DAMap, i)
  nhv = nhvars(m)
  ip = zeros(real(eltype(m.v0)), nhv, nhv)
  ip[2*i-1, 2*i-1] = 1
  ip[2*i, 2*i] = 1
  return ip
end

function jp_mat(m::DAMap{V0}, i) where {V0<:StaticArray}
  nhv = nhvars(m) 
  return StaticArrays.sacollect(SMatrix{nhv,nhv,real(eltype(V0))}, 
  begin
    if col == 2*i && row == 2*i-1
      1
    elseif col == 2*i-1 && row == 2*i
      -1
    else
      0
    end
  end for col in 1:nhv for row in 1:nhv)
end

function jp_mat(m::DAMap, i)
  nhv = nhvars(m)
  jp = zeros(real(eltype(m.v0)), nhv, nhv)
  jp[2*i-1, 2*i] = 1
  jp[2*i, 2*i-1] = -1
  return jp
end

function j_mat(m::DAMap{V0}) where {V0<:StaticArray}
  nhv = nhvars(m) 
  return StaticArrays.sacollect(SMatrix{nhv,nhv,real(eltype(V0))}, 
  begin
    if fld1(col,2) != fld1(row,2) 
      0
    else # then we are in the block
      if mod1(row,2) == 1 && mod1(col,2) == 2 # First row second col is 1
        1
      elseif mod1(row,2) == 2 && mod1(col,2) == 1
        -1
      else
        0
      end
    end
  end for col in 1:nhv for row in 1:nhv)
end

function j_mat(m::DAMap)
  nhv = nhvars(m)
  j = zeros(real(eltype(m.v0)), nhv, nhv)
  for i in 1:Int(nhv/2)
    j[2*i-1, 2*i] = 1
    j[2*i, 2*i-1] = -1
  end
  return j
end




