program example
use c_tpsa
implicit none 
integer no,nv
type(c_taylor) f,g


no=4; nv=3;    ! no: the order of the polynomial    nv: the number of variables
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init


f=(2.d0.cmono.'12') + (3.d0.cmono.'11') + 4.d0  !  Creates 2.d0 x_1 x_ 2 ^2+3.d0x_1x_2 +4.d0
g=f.k.1       !     Creates 2.d0 x_2^2 + 3.d0 x_2

call print(f,6)
call print(g,6)



f=(2.d0.cmono.'121') + (3.d0.cmono.'103') + 4.d0  !  Creates 2.d0 x_1 x_ 2^2 x_3 +3.d0 x_1x_3^3 +4.d0
g=f.k.2       !     Creates 2.d0x_1 x_2  x_3 

call print(f,6)
call print(g,6)



call kill(f,g)      ! must be destroyed
end program example