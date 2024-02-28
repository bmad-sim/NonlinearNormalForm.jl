program example
use c_tpsa
implicit none 
integer no,nd,i,nv
type(c_taylor) f
type(c_damap) m,nl
type(c_vector_field) vf
type(c_universal_taylor) uf
integer, allocatable :: j(:) ,k(:)
real(dp) r1,r2
type(c_normal_form) normal
n_cai=-i_
no=5; nd= 1;    ! no: the order of the polynomial    nv: the number of variables   
use_quaternion=.true.
call c_init(no,nd,0)  ! initializes taylor series with maps
call alloc(f)      
call alloc(vf)      
call alloc(m,nl)
call alloc(normal)

r1=1.0_dp
r2=0.5d0

call c_TAYLOR_ran(f,r1,r2)
f=f-i_*aimag(f)  ! remove the imaginary part
f=f-(f.cut.2)   ! remove linear part
vf=getvectorfield(f)  ! producing a symplectic map
do i=1,3
 call c_TAYLOR_ran(vf%q%x(i),r1,r2)
 vf%q%x(i)=vf%q%x(i)-i_*aimag(vf%q%x(i))
enddo


m=exp(vf)

 call print(m)
 
write(6,*) " random map of order ",no
call print(m)
 
nl=m.sub.2
write(6,*) " extracting order 2 "
call print(nl)

nl=m.cut.3 
write(6,*) " cutting order 3 and above "
call print(nl)

nl=m.cut.(-3) 
write(6,*) " cutting order 3  for orbital and 2 for quaternion "
call print(nl)


call kill(f)      
call kill(vf)      
call kill(m,nl)
call kill(normal)
end program example