program example
  use c_tpsa
  implicit none 
  integer no,nv,i
  type(c_taylor) f,g
  
  
  no=4; nv= 2;  ! no: the order of the polynomial    nv: the number of variables
  call c_init(no,nv)  ! initializes taylor series without maps
  
  call alloc(f,g)      ! must be constructed after init
  
  i=3
  f=(2.d0.cmono.'11') + (3.d0.cmono.'21') + (4.d0.cmono.'04') + 5.d0  
    !  Creates 2.d0 x_1 x_ 2 + 3.d0x_1^2x_2 + 4.d0x_2^4 + 5.d0
  g=f.cut.i        !  Creates 2.d0 x_1 x_2 + 5.d0
  call print(f,6)
  call print(g,6)
  
  
  
  
  
  i=1
  f=(2.d0.cmono.'11') + (3.d0.cmono.'21') + (4.d0.cmono.'04') + 5.d0  
    !  Creates 2.d0 x_1 x_ 2 + 3.d0x_1^2x_2 + 4.d0x_2^4 + 5.d0
  g=f.cut.i        !  Creates 5.d0
  call print(f,6)
  call print(g,6)
  
  
  
  
  call kill(f,g)      ! must be destroyed
  end program example