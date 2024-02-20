program example
use c_tpsa
implicit none 
integer no,np,i
type(c_damap) map,id
type(c_taylor) h
type(c_vector_field) hf
type(c_taylor) k,m,x,p,pb
real(dp) time

no=4; np= 0;    ! no: the order of the polynomial    nv: the number of variables

! Here we use the Poisson bracket on a "real" map


CALL c_INIT_all(NO1=no,nd_or_nv=1,NP1=np,NDPT1 =0)  

call alloc(k,m,x,p,pb)      ! must be constructed after init
call alloc(h)      ! must be constructed after init
call alloc(map,id)      ! must be constructed after init
call alloc(hf)      ! must be constructed after init

x=1.d0.mono.1
p=1.d0.mono.2

time=0.01d0 ! seconds
k=2.d0
m=0.01d0

! Hamiltonian Lie Operator for a spring 

h= -time*(p**2/2/m + k*x**2/2)
hf=getvectorfield(h)
call print(hf,6)

!  The complex FPP has no Poisson Bracket operators
id=1

map=texp(hf,id)

call print(id,6)
call print(map,6)

pb=(id%v(1)).pb.(id%v(2))

write(6,*) " [x_0,p_0] = 1"
call print(pb,6)
 
! Preservation of the Poisson Bracket after a time of  0.01 
write(6,*) " [x_t,p_t] = 1"

pb=(map%v(1)).pb.(map%v(2))

call print(pb,6)
 

call kill(h)      ! must be constructed after init
call kill(map,id)      ! must be constructed after init
call kill(k,m,x,p,pb)      ! must be constructed after init

end program example