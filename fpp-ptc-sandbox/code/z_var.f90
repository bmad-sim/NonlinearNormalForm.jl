! These routines are obsolete and should be replaced by the operator  .mono.
! They originally overloaded Berz's DAVAR and should be avoided

program example
  use c_tpsa
  implicit none 
  integer no,nv,i,ii
  type(c_taylor) f
  complex(dp) r0, r(2)   ! or real(dp)
  
  no=4; nv= 2;    ! no: the order of the polynomial    nv: the number of variables   
  call c_init(no,nv)  ! initializes taylor series without maps
  
  call alloc(f)      ! must be constructed after init
  
  ! creates f=3.d0 + (one.cmono.i)
  i=1
  f=3.d0+dz_c(1)    ! Creates 3.d0 + x_1
  
  call print(f,6)
  
  
  ! creates f=r(1) + (r(2).cmono.2)
  r(1)=5.d0 ; 
  r(2)=6.d0
  f=r(1)+r(2)*dz_c(2)    ! Creates 5.d0 + (6.d0 x_2)
  
  call print(f,6)
  
  
  call kill(f)      ! must be destroyed
  end program example