struct DeMoivre{S}
  H::S
  B::S
  E::S
  K::S
end

function compute_de_moivre(a1::DAMap{V0}, ::Val{linear}=Val{false}()) where {V0<:StaticArray,linear}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix
  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  if linear
    let a1_mat = jacobian(a1, VARS_CPARAM), a1i_mat = inv(a1_mat)
      B = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:Int(nv/2))
      K = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, -j_mat(a1)*B[i] for i in 1:Int(nv/2))
      E = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, -B[i]*j_mat(a1) for i in 1:Int(nv/2))
      H = StaticArrays.sacollect(SVector{Int(nv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nv/2))
      return DeMoivre(H, B, E, K)
    end
  else
    let
      tmp = zero(a1)
      a1i = inv(a1)
      B = StaticArrays.sacollect(SVector{Int(nv/2)}, begin
        setray!(tmp.v, v_matrix=jp_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nv/2)
      )
      setray!(tmp.v, v_matrix=-j_mat(a1))
      K = StaticArrays.sacollect(SVector{Int(nv/2)}, begin
        tmp∘B[i]
      end for i in 1:Int(nv/2)
      )
      E = StaticArrays.sacollect(SVector{Int(nv/2)}, begin
        B[i]∘tmp
      end for i in 1:Int(nv/2)
      )
      H = StaticArrays.sacollect(SVector{Int(nv/2)}, begin
        setray!(tmp.v, v_matrix=ip_mat(a1, i))
        a1∘tmp∘a1i
      end for i in 1:Int(nv/2)
      )
      return DeMoivre(H, B, E, K)
    end
  end
end

function compute_de_moivre(a1::DAMap, ::Val{linear}=Val{false}()) where {linear}
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
      return DeMoivre(H, B, E, K)
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
      return DeMoivre(H, B, E, K)
    end
  end
end

# =================================================================================== #
# Construct special matrices for lattice functions
function ip_mat(m::DAMap{V0}, i) where {V0<:StaticArray}
  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  return StaticArrays.sacollect(SMatrix{nv,nv,real(eltype(V0))}, 
  begin
    if col == 2*i-1 && row == 2*i-1
      1
    elseif col == 2*i && row == 2*i
      1
    else
      0
    end
  end for col in 1:nv for row in 1:nv)
end

function ip_mat(m::DAMap, i)
  nv = nvars(m)
  if isodd(nv) 
    nv += 1
  end
  ip = zeros(real(eltype(m.v0)), nhv, nhv)
  ip[2*i-1, 2*i-1] = 1
  ip[2*i, 2*i] = 1
  return ip
end

function jp_mat(m::DAMap{V0}, i) where {V0<:StaticArray}
  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  coast = isodd(length(V0))
  coast_plane = coast && i == Int(nv/2)
  return StaticArrays.sacollect(SMatrix{nv,nv,real(eltype(V0))}, 
  begin
    if col == 2*i && row == 2*i-1
      1
    elseif col == 2*i-1 && row == 2*i #&& !coast_plane
      -1
    else
      0
    end
  end for col in 1:nv for row in 1:nv)
end

function jp_mat(m::DAMap, i)
  nv = nvars(m)
  if isodd(nv) 
    nv += 1
  end
  jp = zeros(real(eltype(m.v0)), nv, nv)
  jp[2*i-1, 2*i] = 1
  jp[2*i, 2*i-1] = -1
  return jp
end

function j_mat(m::DAMap{V0}) where {V0<:StaticArray}
  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  return StaticArrays.sacollect(SMatrix{nv,nv,real(eltype(V0))}, 
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
  end for col in 1:nv for row in 1:nv)
end

function j_mat(m::DAMap)
  nv = nvars(m)
  if isodd(nv) 
    nv += 1
  end
  j = zeros(real(eltype(m.v0)), nv, nv)
  for i in 1:Int(nv/2)
    j[2*i-1, 2*i] = 1
    j[2*i, 2*i-1] = -1
  end
  return j
end




