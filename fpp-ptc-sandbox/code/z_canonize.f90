program example
use c_tpsa
use gauss_dis
implicit none 
integer no,nd,i,nv,k,j(lnv),n
type(c_taylor) f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin
type(c_damap) m,m1,m2,m3,as,a0,a1,a2,a_cs
type(c_vector_field) vf1,vf2,vf3  !,h_rot,K_r,K_r_real
real(dp) r1,r2,sca,decrement
complex(dp) c1
type(c_normal_form) normal
type(c_ray) ray,ray2
type(c_lattice_function)   Lat_function
logical :: damp = .false.,fastcanonise=.true.
type(c_taylor), allocatable :: phase(:)
type(c_damap), allocatable :: mt(:)
real(dp)  ph(3),spintune(2),dampdec(3)
 

n_cai=-i_
no=3; nd= 3;    ! no: the order of the polynomial    nv: the number of variables   
c_lda_used=15000
use_quaternion=.true.
call c_init(no,nd,0)  ! initializes taylor series with maps
allocate(phase(nd))
 
call alloc(f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin)      ! must be constructed after init
call alloc(vf1,vf2,vf3)      
call alloc(m,m1,m2,m3,as,a0,a1,a2,a_cs)      ! 
call alloc(phase)
call alloc(nu_spin)

n=20
allocate(mt(n))
do i=1,n
 call alloc(mt(i))
enddo
 !default_fractional_tune_positive=.false.
  negative_synchrotron_tune=.false.

damp=.true.
fastcanonise=.true.
if(.not.fastcanonise) damp=.false.
do_damping=damp

 decrement=.001d0
 call alloc(normal)

r1=1.0_dp
r2=0.5d0

!
!write(6,'(A,/)') "creating a fifth order m with 2 parameters in 1-d-f"

do k=1,n
vf1=0

call c_TAYLOR_ran(f,r1,r2)
f=f-i_*aimag(f)  ! remove the imaginary part
f=f-(f.cut.2)
do i=1,nd
sca=0
do while(sca==0.d0) 
  call GRNF(sca,1.d0)
enddo
!sca=sign(sca,1.d0)+sca  
sca=abs(sca)*1000
 f=f+sca*(dz_c(2*i-1)**2+dz_c(2*i)**2)
enddo
 

f=-f*0.002d0
 
vf1=getvectorfield(f)  !  
if(damp) then
  do i=1,nd
     call GRNF(sca,decrement)
     sca=1.d0-sca
     vf1%v(2*i)=sca*vf1%v(2*i)
  enddo
endif
! putting spin
do i=1,3
 call c_TAYLOR_ran(vf1%q%x(i),r1,r2) ! putting spin
 vf1%q%x(i)=vf1%q%x(i)-i_*aimag(vf1%q%x(i))
 vf1%q%x(i)=5.d0*vf1%q%x(i) !   
enddo
 
mt(k)=exp(vf1)

enddo

m=1
do i=1,n
 m=mt(i)*m
enddo
 
call c_normal(m,normal,canonise=.false.,dospin=.true.,phase=phase,nu_spin=nu_spin)
do_damping=.true.
!call c_fast_canonise(normal%atot, a_cs)
!call print(a_cs)
!stop

call clean(phase,phase,prec=1.d-5)
call clean(nu_spin,nu_spin,prec=1.d-5)

write(6,'(A,3(1x,g23.16),/)') "The tune    is ",normal%tune(1:nd)
write(6,'(A,3(1x,g23.16),/)') "The damping is ",normal%damping(1:nd)
write(6,'(A,1x,g23.16,/)') "The spin tune  is ",normal%spin_tune 
write(6,'(a)') "Orbital tunes "
call print(phase)
write(6,'(a)') "Spin tune "
call print(nu_spin)

 

if(fastcanonise) then
 call c_fast_canonise(normal%atot,a_cs,dospin=.true.)
else
 call c_full_canonise(normal%atot,a_cs,as,a0,a1,a2)   
endif

do i=1,nd
phase(i)=0.d0
enddo
nu_spin=0.d0
ph=0
spintune=0
dampdec=0
 
! Tracking the canonised a_cs
!do k=1,n
a_cs = m*a_cs
!a_cs=mt(k)*a_cs
 if(fastcanonise) then
    call c_fast_canonise(a_cs,a_cs,phase=ph,damping=dampdec,spin_tune=spintune,dospin=.true.)  
  else
    call c_full_canonise(a_cs,a_cs,as,a0,a1,a2,phase=phase,nu_spin=nu_spin)
 endif
!enddo

call clean(phase,phase,prec=1.d-10)
call clean(nu_spin,nu_spin,prec=1.d-10)
write(6,'(a)') "Total Orbital Tunes (phase advances) "
if(fastcanonise) then
 write(6,*) ph(1:nd)
else
 call print(phase)
endif
write(6,'(a)') "Total Spin Tune (phase advance) "
if(fastcanonise) then
 write(6,*) spintune
else
 call print(nu_spin)
endif

if(fastcanonise) then
 write(6,'(a)') "Total damping "
 write(6,*)  dampdec(1:nd)
endif

write(6,'(A,/)') "  Fractional Results "
write(6,'(A,3(1x,g23.16),/)') "The tune    is ",normal%tune(1:nd)
write(6,'(A,3(1x,g23.16),/)') "The damping is ",normal%damping(1:nd)
write(6,'(A,2(1x,g23.16),/)') "The spin tune and spin angle are ",normal%spin_tune, &
normal%quaternion_angle/pi 

end program example