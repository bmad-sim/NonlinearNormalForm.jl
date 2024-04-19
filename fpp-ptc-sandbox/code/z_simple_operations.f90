program example
use c_tpsa
implicit none 
integer no,nd,i,nv
type(c_taylor) f
type(c_damap) m,mt,A,mtr,Identity
type(c_vector_field) vf,vfr,vft,vtemp
type(c_universal_taylor) uf
integer, allocatable :: j(:) ,k(:)
real(dp) r1,r2,sca
complex(dp) scac
n_cai=-i_
no=2; nd= 1;    ! no: the order of the polynomial    nv: the number of variables   
use_quaternion=.true.
call c_init(no,nd,0)  ! initializes taylor series with maps

call alloc(f)      
call alloc(vf)      
call alloc(vfr)      
call alloc(vft)
call alloc(vtemp)
      
call alloc(m,mt,A,mtr,Identity)
 
r1=1.0d0
r2=0.d0
sca=0.1d0

call c_TAYLOR_ran(f,r1,r2)
f=f-i_*aimag(f)  ! remove the imaginary part
f=f-(f.cut.2)   ! remove linear part
f=f*sca   !  to make map near the identity for log
vf=getvectorfield(f)  ! producing a symplectic map
do i=1,3
 call c_TAYLOR_ran(vf%q%x(i),r1,r2) ! putting spin
 vf%q%x(i)=vf%q%x(i)-i_*aimag(vf%q%x(i))
 vf%q%x(i)=vf%q%x(i) *sca !  to make map near the identity for log
enddo

 call c_TAYLOR_ran(f,r1,r2)
f=f-i_*aimag(f)  ! remove the imaginary part
f=f-(f.cut.2)   ! remove linear part
f=f*sca   !  to make map near the identity for log
vfr=getvectorfield(f)  ! producing a symplectic map
do i=1,3
 call c_TAYLOR_ran(vfr%q%x(i),r1,r2) ! putting spin
 vfr%q%x(i)=vfr%q%x(i)-i_*aimag(vfr%q%x(i))
 vfr%q%x(i)=vfr%q%x(i) *sca !  to make map near the identity for log
enddo


 
Identity=1

m=exp(vf,Identity)    !  m = exp(vf.grad)Identity    ! 

write(*,*) "m before mul:"
call print(m)

m = 2.0d0*m
write(*,*) "m after mul:"
call print(m)
stop

A=exp(vfr)   !  A = exp(vfr.grad)Identity   Identity is optional.


write(6,*) "Checking addition and subtraction on maps"
mtr=m+a

mtr=mtr-a

call print(m%v(1))
call print(mtr%v(1))
call print(m%q%x(1))
call print(mtr%q%x(1))

write(6,*) "Checking unary '-' " 
mtr=m+a
a=(-1.d0)*a
mtr=mtr+a

call print(m%v(1))
call print(mtr%v(1))
call print(m%q%x(1))
call print(mtr%q%x(1))

!scac=-1.d0
!mtr=m+a
!a=scac*a
!mtr=mtr+a

write(6,*) "Checking addition and subtraction on vector fields"

vfr=-vf


call print(vf%v(1))
call print(vfr%v(1))
call print(vf%q%x(1))
call print(vfr%q%x(1))

vfr=(-1.d0)*vf


call print(vf%v(1))
call print(vfr%v(1))
call print(vf%q%x(1))
call print(vfr%q%x(1))


scac=-i_
vfr=scac*vf


call print(vf%v(1))
call print(vfr%v(1))
call print(vf%q%x(1))
call print(vfr%q%x(1))

call kill(f)      
call kill(vf)      
call kill(vfr)      
call kill(vft)      
call kill(m,mt,A,mtr,Identity)

end program example