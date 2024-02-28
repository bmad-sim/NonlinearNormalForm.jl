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
A=exp(vfr)   !  A = exp(vfr.grad)Identity   Identity is optional.

mt=A**(-1)*m*A   ! similarity transformation

vft=exp_ad(vfr,vf) !   vft = exp(<vfr,>)vf

mtr=exp(vft)       ! A**(-1)*m*A = exp(vft.grad)I
 
mtr=mt*mtr**(-1)

call print(mtr)


!   transforming the vector field with A

vft=A*vf !
 
mtr=exp(vft)       ! A**(-1)*m*A = exp(vft.grad)I
 
mtr=mt*mtr**(-1)

call print(mtr)

!Lie bracket used for example in exp_ad
!vft=exp_ad(vfr,vf) !   vft = exp(<vfr,>)vf

write(6,*) " Cheap program reproducing exp_ad using .lb. "
vft=vf
vtemp=vf
r1=1
do i=1,100
r1=i 
vtemp=((1.d0/r1)*(vfr.lb.vtemp))
 vft=vtemp + vft
enddo
 
mtr=exp(vft)       ! A**(-1)*m*A = exp(vft.grad)I
 
mtr=mt*mtr**(-1)

call print(mtr)

call kill(f)      
call kill(vf)      
call kill(vfr)      
call kill(vft)      
call kill(m,mt,A,mtr,Identity)

end program example