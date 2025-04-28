struct LatticeFunctions{S}
  H::S
  B::S
  E::S
  K::S
end


function compute_lattice_functions(a1::DAMap{V0}) where {V0<:StaticArray}
  # jp_mat[i] in FPP is J matrix restricted to i-th plane
  # ip_mat[i] in FPP is identity matrix restricted to i-th plane
  # jt_mat in FPP is symplectic s matrix

  let a1_mat = parametric_jacobian(a1), a1i_mat = iseven(length(V0)) ? inv(a1_mat) : parametric_jacobian(inv(a1))
    nhv = nhvars(a1) 
    B = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, a1_mat*jp_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2))
    K = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, -j_mat(a1)*B[i] for i in 1:Int(nhv/2))
    E = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, -B[i]*j_mat(a1) for i in 1:Int(nhv/2))
    H = StaticArrays.sacollect(SVector{Int(nhv/2),typeof(a1_mat)}, a1_mat*ip_mat(a1, i)*a1i_mat for i in 1:Int(nhv/2))
    return LatticeFunctions(H,B,E,K)
  end
  
#=
  if V0 <: StaticArray

  else

  end


  
  # Computing the de Moivre Lattice Functions
  # equivalent to Ripken ones
  # Coefficients of vector field matrices  B**2=-H
  #=
     do i = 1,3 ! jp_mat is symplectic J matrix (I call S matrix)
      f%B(i,:,:)=matmul(matmul(m%mat,jp_mat(i,:,:)),mi%mat)
     enddo

!!!! Coefficients of invariants
      do i = 1,3
      f%K(i,:,:)= -matmul(jt_mat,f%B(i,:,:))
     enddo
!!!! Coefficient of moments as functions of emittances
     do i = 1,3
      f%E(i,:,:)= -matmul(f%B(i,:,:),jt_mat)
     enddo
!!!!  Dispersive quantities containing zeta and eta for example
!!!!  as well as \gamma C of Teng-Edwards
     do i = 1,3
      f%H(i,:,:)=matmul(matmul(m%mat,ip_mat(i,:,:)),mi%mat)
     enddo
     =#


=#
end

function parametric_jacobian(a1::DAMap{V0,V}) where {V0<:StaticArray,V}
  if nparams(a1) == 0
    return jacobian(a1, HVARS) # scalar part
  else
    nhv = nhvars(a1) #isodd(length(V0)) ? length(V0)+1 : length(V0)
    return StaticArrays.sacollect(SMatrix{nhv,nhv,eltype(V)}, 
    begin
      t = zero(a1.v[row])
      TI.deriv!(t, a1.v[row], col)  
      t
    end for col in 1:nhv for row in 1:nhv)
  end
end

function parametric_jacobian(a1::DAMap) 
  if nparams(a1) == 0
    return jacobian(a1, HVARS) # scalar part
  else
    nhv = nhvars(a1)
    M = zeros(eltype(a1.v), nhv, nhv)
    for col in 1:nhv
      for row in 1:nhv
        t = zero(a1.v[row])
        TI.deriv!(t, a1.v[row], col)  
        M[row, col] = t
      end
    end
    return M
  end
end

# =================================================================================== #
# Construct special matrices for lattice functions
function ip_mat(m::DAMap{V0}, i) where {V0<:StaticArray}
  nhv = nhvars(m)
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
  if isodd(nhv)
    nhv += 1
  end
  j = zeros(real(eltype(m.v0)), nhv, nhv)
  for i in 1:Int(nhv/2)
    j[2*i-1, 2*i] = 1
    j[2*i, 2*i-1] = -1
  end
  return j
end




