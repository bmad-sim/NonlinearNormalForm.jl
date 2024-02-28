program example
use c_TPSA
implicit none 
integer no,nv
type(c_taylor) f,df,idf,t


no=4; nv= 2;        ! no: the order of the polynomial    nv: the number of variables
call c_init(no,nv)  ! initializes taylor series without maps

call alloc(f,df,idf,t)      ! must be constructed after init

f=(2.d0.cmono.'12') + (3.d0.cmono.'11') + 4.d0  !  Creates 2.d0 x_1 x_ 2 ^2+3.d0 x_1x_2 +4.d0
df=f.d.2       !      df/dx_2    Creates 4.d0 x_1x_2 + 3.d0 x_1
idf=df.i.2       !      int df/dx_2    recreates 2.d0 x_1 x_ 2 ^2+3.d0 x_1x_2 WITHOUT 4.d0

call print(f,6)
call print(df,6)
call print(idf,6)

call kill(f,df,idf,t)      ! must be destroyed
end program example