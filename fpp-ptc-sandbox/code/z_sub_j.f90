program example
use c_tpsa
implicit none 
integer no,nv
real(dp) r
type(c_taylor) f,g

integer, allocatable :: j(:),k(:) 

no=4; nv= 2;    ! no: the order of the polynomial    nv: the number of variables   
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init

allocate(j(nv))
allocate(k(nv))
  
j=0
j(1)=1;j(2)=2;

k=0
k(1)=3;k(2)=1;


f=(2.d0.cmono.j)  + 5.d0   !  Creates 2.d0 x_1 x_ 2 ^2 + 5.d0
r=f.sub.j     ! Peeks coefficient   2.d0

call print( f,6)
write(6,*) r


f=(2.d0.cmono.j)  +   (4.d0.cmono.k) + 5.d0   !  Creates 2.d0 x_1 x_ 2 ^2  + 4.d0 x_1 ^3 x_ 2  + 5.d0
r=f.sub.k   ! Peeks coefficient   4.d0  

call print( f,6)
write(6,*) r



f=(2.d0.cmono.'12' ) + (4.d0.cmono.'13' ) + ( 3.d0.cmono.'01')  + 5.d0
    !  Creates 2.d0 x_1 x_ 2 ^2 + 4.d0 x_1 x_ 2 ^3 + 3.d0 x_ 2 + 5.d0
r=f.sub.'13'        ! Peeks coefficient   4.d0

call print( f,6)
write(6,*) r


deallocate(j)
call kill(f,g)      ! must be destroyed
end program example