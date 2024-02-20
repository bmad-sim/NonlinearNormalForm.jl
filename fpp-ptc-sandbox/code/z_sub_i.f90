program example
use c_tpsa
implicit none 
integer no,nv,i
type(c_taylor) f,g


no=4; nv= 2;  ! no: the order of the polynomial    nv: the number of variables
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init

i=4
f=(2.d0.cmono.'12') + (3.d0.cmono.'31') + (4.d0.cmono.'04') + 5.d0    
!  Creates 2.d0 x_1 x_ 2 ^2 + 3.d0 x_1^3x_2 + 4.d0 x_2^4 + 5.d0
g=f.sub.i        !  Creates 3.d0 x_1^3x_2 + 4.d0 x_2^4

call print(f,6)
call print(g,6)


call kill(f,g)      ! must be destroyed
end program example