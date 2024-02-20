program example
use c_tpsa
implicit none 
integer no,nd,i,nv
type(c_taylor) f
type(c_damap) m,m2
type(c_vector_field) vf,vfr
type(c_universal_taylor) uf
integer, allocatable :: j(:) ,k(:)
real(dp) r1,r2,sca
type(c_normal_form) normal
n_cai=-i_
no=2; nd= 1;    ! no: the order of the polynomial    nv: the number of variables   
use_quaternion=.true.
call c_init(no,nd,0)  ! initializes taylor series with maps
call alloc(f)      
call alloc(vf)      
call alloc(vfr)      
call alloc(m,m2)
call alloc(normal)

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

 
 

m=exp(vf)

write(6,*) " random vector field of order ",no
call print(vf)
 
vfr=vf.sub.2
write(6,*) " extracting order 2 "
call print(vfr)

vfr=vf.cut.3 
write(6,*) " cutting order 3 and above "
call print(vfr)

vfr=vf.cut.(-3) 
write(6,*) " cutting order 3  for orbital and 2 for quaternion "
call print(vfr)

write(6,*) " computing logarithm and comparing "
!  F=ln(M,H,epso,n,tpsa)  then   M=exp(F)I
! H is a guess for F is known
! epso = is a small number (optional) for convergence
! n = maximal number of iteration (optional, defaulted to 1000)
!  tpsa is true by default, otherwise it is a DA calculation 

vfr=ln(m)

vfr=vf-vfr
 
call print(vfr)

!!!! example of exponential which amounts to squaring a map
write(6,*) " squaring a map with Lie exponent and comparing "

m2=exp(vf,m)

m2=m2-(m**2)

call print(m2)

call kill(f)      
call kill(vf)      
call kill(m,m2)
call kill(normal)
end program example