program example
  use c_tpsa
  implicit none 
  integer no,nv,i
  type(c_taylor) f
  integer, allocatable :: j(:) ,k(:)
  
  no=4; nv= 2;    ! no: the order of the polynomial    nv: the number of variables   
  call c_init(no,nv)  ! initializes taylor series without maps
  
  call alloc(f)      ! must be constructed after init
  
  allocate(j(nv),k(nv)) 
  
  j=0
  j(1)=1;j(2)=2;
  
  k=0
  k(1)=3;k(2)=1;
  
  f=(2.d0.cmono.j) +  (4.d0.cmono.k )  !  Creates 2.d0 x_1 x_ 2 ^2 + 4.d0 x_1^3 x_ 2 
  
  call print(f,6)
  
  f=(2.d0.cmono.'12') + (3.d0.cmono.'11')    !  Creates 2.d0 x_1 x_ 2 ^2+3.d0x_1x_2
  
  call print(f,6)
  
  i=2
  f= (2.d0.cmono.i ) +  (3.d0.cmono.(i -1))    ! Creates 2.d0 x_2 +3.d0 x_1   <--------  changed
  
  call print(f,6)
  
  
  
  deallocate(j,k)
  call kill(f)      ! must be destroyed
  end program example