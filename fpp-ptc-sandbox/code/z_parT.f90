program example
use c_tpsa
implicit none 
integer no,nv,n,i,k
type(c_taylor) f,g
type(sub_taylor) inf

no=6; nv= 4;    ! no: the order of the polynomial    nv: the number of variables   
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,g)      ! must be constructed after init


inf%j=0
inf%min=2
inf%max=3
inf%j(2)=1
inf%j(3)=1

f=(2.d0.cmono.'2101')  +  (3.d0.cmono.'2112' ) + (6.d0.cmono.'211') + 5.d0+(6.d0.cmono.'111') 
     !  Creates 2.d0 x_1^2  x_2 x_ 4 + 3.d0 x_1^2  x_2 x_3 x_ 4^2 + 6.d0 x_1^2 x_2 x_3 + 5.d0 + 6.d0 x_1 x_2 x_3 
g=f.parT.inf      !  Extracts  3.d0 x_1^2  x_2 x_3 x_ 4^2 + 6.d0 x_1^2 x_2 x_3 + 6.d0 x_1 x_2 x_3 
call print(f,6)
call print(g,6)


call kill(f,g)      ! must be destroyed
end program example