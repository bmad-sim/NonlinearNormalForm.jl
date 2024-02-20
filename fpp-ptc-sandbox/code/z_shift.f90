program example
use c_tpsa
implicit none 
integer no,nv,k
type(c_taylor) f,g
 

no=6; nv= 4;    ! no: the order of the polynomial    nv: the number of variables   
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init


k=1
f=(2.d0.cmono.'023')+(1.d0.cmono.'011')+(4.d0.cmono.'004')+3.d0
g=f<=k     ! Shifts exponents downwards by k
  !  Creates (2.d0 x_1^2 x_2^3) +  (1.d0 + x_1 x_2) + (4.d0 x_2^4)  + 3.d0  <------------ comment added
call print(f,6)
call print(g,6)


k=2
f=(2.d0.cmono.'0032')+(1.d0.cmono.'0011')+(4.d0.cmono.'0040')+3.d0
g=f<=k     ! Shifts exponents downwards by k
 !  Creates (2.d0 x_1^3 x_2^2) +  (1.d0 + x_1 x_2) + (4.d0 x_1^4) +3.d0   <------------ comment added
call print(f,6)
call print(g,6)

!  This is not permissible
! FPP will return the exception 
! "trouble in dashift"
! because the fist monomial
! (2.d0.cmono.'0112') cannot be shifted by 2 notches

!k=2
!f=(2.d0.cmono.'0112')+(1.d0.mono.'0011')+(4.d0.cmono.'0040')+3.d0
!g=f<=k     




call kill(f,g)      ! must be destroyed
end program example