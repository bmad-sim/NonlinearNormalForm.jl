program normal_simple
use pointer_lattice
implicit none
type(c_damap) T,N
type(c_taylor) r2,cs_invariant
type(c_normal_form) Normal_form
real(dp) mat(2,2),alpha,beta,gamma
integer index_beta(2),index_gamma(2),index_alpha(2),mf
type(internal_state) state
mf = 6
! This program initialize FPP via PTC
! which is more in line with a BMAD usage
! but is obviously not necessary.
! one degree of freedom in PTC
! defined via an internal state
!call kanalnummer(mf," Results.tex")
! Set the phase space
! to 2-dim or 1-degree-of-freedom
state=only_2d0
! Initializes FPP to order 2
call init(state,2,0)
!call ptc_ini_no_append ! initializes PTC
call alloc(T); call alloc(N);
call alloc(r2);call alloc(cs_invariant);
call alloc(Normal_form)
mat(1,1)=1; mat(1,2)=1;mat(2,1)=-0.4_dp; mat(2,2)=0.6_dp;
! Matrix mat(2,2) put into a c_damap to create T
T=mat
call print(T)
call c_normal_new(T,Normal_form)
N=Normal_form%Atot**(-1)*T*Normal_form%Atot
! Creating x^2+p^2 which is the invariant of a rotation
r2=(1.0_dp.cmono.1)**2+(1.0_dp.cmono.2)**2
! Creating the invariant in the lab frame
cs_invariant=r2*Normal_form%Atot**(-1)
write(mf,*) " The map T"
call print(T,mf);write(mf,*);
write(mf,*) " The map A"
call print(Normal_form%Atot,mf);write(mf,*);
write(mf,*) " The map N : a rotation"
call print(N,mf);write(mf,*);
write(mf,*) " The invariant of the normal form ";write(mf,*);
call print(r2,mf);write(mf,*);
write(mf,*) " The Courant-Snyder Invariant";write(mf,*);
call print(cs_invariant,mf);write(mf,*);
!!! computing 1+alpha**2 and beta*gamma
index_gamma(1)=2; index_gamma(2)=0;
index_alpha(1)=1; index_alpha(2)=1;
index_beta(1)=0; index_beta(2)=2;
gamma=cs_invariant.sub.index_gamma
alpha=(cs_invariant.sub.index_alpha)/2
beta=cs_invariant.sub.index_beta
write(mf,'(a)',advance="no") "gamma, alpha, beta"
write(mf,"(3(1x,g13.6,1x))") gamma, alpha, beta
write(mf,'(a)',advance="no") "1.0_dp+alpha**2"
write(mf,*) 1.0_dp+alpha**2
write(mf,'(a)',advance="no") " beta*gamma"
write(mf,*) beta*gamma
write(mf,'(a)',advance="no") " The tune "
write(mf,*) Normal_form%tune(1)
! N is phasors basis (diagonal)
N=Ci_phasor()*N*C_phasor()
! Removes small numbers
call clean(N,N,prec=1.d-10)
write(mf,*) " The map N in phasors : diagonal";write(mf,*);
call print(N,mf);write(mf,*);
close(mf)
end program normal_simple