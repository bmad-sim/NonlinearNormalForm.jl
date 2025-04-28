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

  a1_mat = parametric_jacobian(a1)

  #a1i_mat::(typeof(a1_mat)) = zero(a1_mat)
  #if iseven(length(V0))
  #  a1i_mat = inv(a1_mat)
  #else
    a1i_mat = parametric_jacobian(inv(a1))
  #end
  nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
  B = StaticArrays.sacollect(SVector{Int(nv/2)}, 
  begin
    a1_mat*jp_mat(a1, i)*a1i_mat
  end for i in 1:Int(nv/2))
    return B
  
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
    return jacobian(a1, VARS) # scalar part
  else
    nv = isodd(length(V0)) ? length(V0)+1 : length(V0)
    return StaticArrays.sacollect(SMatrix{nv,nv,eltype(V)}, 
    begin
      t = zero(a1.v[row])
      TI.deriv!(t, a1.v[row], col)  
      t
    end for col in 1:nv for row in 1:nv)
  end
end

function parametric_jacobian(a1::DAMap) 
  if nparams(a1) == 0
    return jacobian(a1, VARS) # scalar part
  else
    nv = nvars(a1)
    if isodd(nv)
      nv += 1
    end
    M = zeros(eltype(a1.v), nv, nv)
    for col in 1:nv
      for row in 1:nv
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
  ip = zeros(real(eltype(m.v0)), nv, nv)
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
    elseif col == 2*i-1 && row == 2*i && !coast_plane
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
  if !iscoasting(m) || i != Int(nv/2)
    jp[2*i, 2*i-1] = -1
  end
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




