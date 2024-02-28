program example
use c_tpsa
implicit none 
integer no,nv,n,i
type(c_taylor) f,g
integer, allocatable :: j(:),k(:),jj(:)

no=6; nv= 4;    ! no: the order of the polynomial    nv: the number of variables   
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init


n=2
allocate(j(nv),k(n),jj(nv)) 

j=0
j(1)=1;j(2)=2;j(3)=2;j(4)=1;

jj=0
jj(1)=1;jj(2)=1;jj(3)=0;jj(4)=3;


k=0
do i=1, n
     k(i)=j(i)
end do

f=(2.d0.cmono.j )+ (3.d0.cmono.jj ) + 5.d0  
     !  Creates 2.d0 x_1 x_ 2 ^2 x_3^2 x_ 4 + 3.d0 x_1 x_ 2  x_ 4^3 + 5.d0
g=f.par.k           !  Creates 2.d0 x_3^2 x_4 
call print(f,6)
call print(g,6)

deallocate(j,k,jj)



f=(2.d0.cmono.'1111')  + (4.d0.cmono.'211') +  (1.d0.cmono.'1112' )  + 5.d0
     ! Creates 2.d0 x_1 x_ 2 x_3 x_ 4 + 4.d0 x_1^2 x_ 2 x_3 + x_1 x_ 2 x_3 x_ 4^2+ 5.d0  
g=f.par.'11'      !  Creates 2.d0 x_3 x_4 + x_3 x_4^2

call print(f,6)
call print(g,6)



f=(2.d0.cmono.'2011')  +  (3.d0.cmono.'2012' ) + (6.d0.cmono.'001') + 5.d0
     !  Creates 2.d0 x_1^2  x_3 x_ 4 + 3.d0 x_1^2  x_3 x_ 4^2 + 6.d0 x_3 + 5.d0
g=f.par.'201'      !  Creates 2.d0 x_4 + 3.d0 x_4^2

call print(f,6)
call print(g,6)



call kill(f,g)      ! must be destroyed
end program example