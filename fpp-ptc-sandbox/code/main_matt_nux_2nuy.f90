program Resonance
use pointer_lattice 
use gauss_dis
implicit none

! ptc variables 
 type(layout), pointer :: ptc_lattice
 type(fibre), pointer ::ptc_fibre
type(internal_state), target :: state
type(c_damap):: T,id,A,Nc,rot,u,uc,fr1,fr2,N,u_0,u_1,U_2,rot_c_inverse
type(c_normal_form) :: c_n
type(probe_8) :: ray_8,xs
type(probe)  :: ray_closed,ray,xs0,xst 
real(dp) closed_orbit(6) !,phi0,phase(3) !,om,g,ang(3),Jy
real(dp) jx(2),discriminant,k1,k0
real(dp) x1(6),x2(6),x1k(6),x2k(6),sigm(6),cut,par(2)
complex(dp) v 
integer i,k,n_arg,mf,mfoo,nturns,nmax0 ,j,mfr,no 
integer, allocatable :: jex(:)
character(100) lat_file
type(c_universal_taylor)   H_res 
type(c_ray) f1,f2
type(c_vector_field) fq,f_comoving,f_barber
type(c_universal_taylor), target ::   a_tot_inverse(6),a_tot(6)
complex(dp) ,allocatable ::f3_1(:),f3_1t(:),f3_3(:)
type(c_ray) fix1,fix2,fix1k,fix2k
real(dp), target ,allocatable ::  vr(:,:),vo(:,:),vi1(:),vi2(:) 
real(dp), target  :: mr
integer nlmin,nlmax,jl,je(6)
real(dp) s,circ,mu_tot,mr2 ,p_res,del, mu(3)
type(c_taylor) hh,phat(3),delorb
logical :: skew=.false.,mis=.false. ,removeit,recomputing_h_r
!complex(dp) vv(3,3),Ar
!real(dp) phir,r,alphar,delr,mphi1(2),mphi2(2)

 preffile= "C:\document\etienne_programs_new\programs_for_learning\park\pref3nux_als.txt"
 !c_lda_used =1500
 
 
call kanalnummer(mf,"C:\document\etienne_programs_new\programs_for_learning\suntao\results_rad.txt")

call ptc_ini_no_append ! initializes PTC

  
 call append_empty_layout(m_u)
 ptc_lattice=>m_u%start
 call build_lattice(ptc_lattice,mis,exact=.false.,thin=.true.,onecell=.true.) 
 
 state=only_4d0 

my_estate=>state

 

!!! next, use PTC normal form to compute intrinsic resonance strength
ptc_lattice=>m_u%start
ptc_fibre=>ptc_lattice%start
 
!call read_ptc_command('C:\document\etienne_programs_new\programs_for_learning\resonance orbital\fittune_resonance_nux_2nuy.txt')

no=3
write(6,*) "order >= 3 (try 7th order) "
read(5,*) no
if(no<3) no=3
!! Counting numbers of multipoles 

ptc_fibre=>ptc_lattice%start
 
use_quaternion=.true.
 n_cai=-i_
  
  
my_estate=>state

 ptc_fibre=>ptc_lattice%start
 
 

call init(state, no, 0)

call alloc(T,id,Nc,rot,a,u_0,u_1,U_2,N,rot_c_inverse)
call alloc(xs)
call alloc(ray_8)
call alloc(c_n)
call alloc(rot)

call alloc(hh,delorb)
call alloc(phat)
call alloc(fq,f_comoving,f_barber)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
closed_orbit=0

! calculation of closed orbit
call FIND_ORBIT_x(closed_orbit,state, 1.0e-7_dp, fibre1=ptc_fibre)

write(6,*) "closed orbit "
 write(6,format4)closed_orbit(1:6)
 
ray_closed=closed_orbit 

 

id=1
ray_8=id+ray_closed
 call propagate(ray_8,+state,fibre1=ptc_fibre)
T=ray_8

call kanalnummer(mfr,"map_res.txt")
 call print(t,mfr)
close(mfr)

 


call c_normal(T,c_n,phase=phat)


 
 
 
  c_n%positive=.false.
  c_n%nres=(c_%no+1)/3
  write(6,*) " # of resonance terms ",c_n%nres
!!!! Keeping a family in 
  do i=1, (c_%no+1)/3
    c_n%m(1,i)=i ! 1/3 order 
    c_n%m(2,i)=-i*2 ! 1/3 order 
  enddo

allocate(vr(c_%nd,c_%nd),vo(c_%nd,c_%nd),vi1(c_%nd),vi2(c_%nd))
vr=0
vo=0

vr(1,1:c_%nd) = c_n%m(1:c_%nd,1)
 write(6,*) vr(1,1:c_%nd)


call gramschmidt(vr)
!! vo is J_1 and J_2 in terms of Jr and Ja
vo=transpose(vr)
 


mr2=0
 do i=1,c_%nd
mr2=mr2+c_n%m(i,1)**2
enddo
 !write(6,*) " mr2 ",mr2
mr=sqrt(mr2)
 
  vr=mr*vr
 vo=(1.d0/mr)*vo
 write(6,*) "************** Vectors m and a **********************"
 write(6,*) vr(1,1:2)
 write(6,*) vr(2,1:2)
write(6,*) "******************************************************"

 
do i=1,c_%nd
mu(i)=phat(i)
enddo
 
call c_normal(T,c_n,canonize=.true.)

write(6,*) "************** tunes **********************"
write(6,*) c_n%tune(1:3)
write(6,*) c_n%spin_tune
write(6,*) "********************************************"

 !!!!!!!! for windows graphical interface : ignore!!!!!!!
  call alloc(a_tot_inverse,1,c_%NV,c_%nd2)
  call alloc(a_tot,1,c_%NV,c_%nd2)

id=1
do i=1,c_%nd2
id%v(i)=id%v(i)+closed_orbit(i)
enddo

a=id.o.c_n%atot

do i=1,c_%nd2
 a_tot(i)=a%v(i) 
enddo

a=c_n%atot**(-1)


id=1
do i=1,c_%nd2
id%v(i)=id%v(i)-closed_orbit(i)
enddo
a=a.o.id

do i=1,c_%nd2
 a_tot_inverse(i)= a%v(i)
enddo

my_euni_1=>a_tot_inverse
my_euni_2=>a_tot
my_evr=>vr
my_evo=>vo
my_evi1=>vi1
my_evi2=>vi2
 my_emr=>mr
  my_estate=>state
 
!!!!!!!!! end of windows graphical interface objects  !!!!!!!!!!!!!

!!!! Put map into a single resonance using A and later in phasors basis

N=c_n%atot**(-1)*T*c_n%atot
N=Ci_phasor()*N*c_phasor()

!!!! Computation of the co-moving rotations

   
 del= (vr(1,1)*mu(1)+vr(1,2)*mu(2))
 p_res=nint(del)
 del=p_res*twopi/mr**2


   do i=1,c_%nd
    f_comoving%v(2*i-1)=-i_*del* vr(1,i)  *dz_c(2*i-1)  
    f_comoving%v(2*i)=  i_*del*  vr(1,i)*dz_c(2*i)  
  enddo

 del= (vr(2,1)*mu(1)+vr(2,2)*mu(2))
  del=del*twopi/mr**2

   do i=1,c_%nd
    f_comoving%v(2*i-1)=-i_*del* vr(2,i)  *dz_c(2*i-1) +f_comoving%v(2*i-1)
    f_comoving%v(2*i)=  i_*del*  vr(2,i)*dz_c(2*i) + f_comoving%v(2*i)
  enddo

 

 rot_c_inverse= exp(-f_comoving)
 Nc=exp(-f_comoving,N) 
 fq=ln(Nc)

  
!!!!!!!!!!!!!  up to this point  the theory was like spin!!!!!!!!!!!    
 

  
Id=n*rot_c_inverse
U_1=rot_c_inverse*n
call c_full_norm_damap(U_1,k0)

u_2=U_1-id
call c_full_norm_damap(U_2,k1)
 write(6,*) "Should be zero 1",k1,k1/k0

 !!!!!!!!!!!!!!!  Hamiltonian vector field fq is computed as :H_r:
 fq=(-1.d0/twopi)*fq
 
  call d_field_for_demin(fq,H_res)
fq=(-twopi)*fq


call clean(H_res,H_res,prec=1.e-6_dp,relative=.true.)
call print(H_res)
 
 write(6,*) " End of Vector field "
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!    There is not clear way to go further as in spin except for one-dimensional
!!!!!!!!!!!!!    resonances: 3nux, 4 nux, etc....
!!!!!!!!!!!!!   It is the spectator aspect of spin that allows the easy computation of the ISF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
write(6,*) "recomputing_h_r "
read(5,*) recomputing_h_r

!!!!!!!!!!!!!!!   doing a style Barber calculation !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!  Removing linear terms in J_x : why not ? !!!!!!!!!!!!!!!!!!!!!!! 
 
 
 do i=1,c_%nd-1 

!  hh is the x-tune of H_r
 hh= (fq%v(2*i).k.(2*i)).sub.'0'
 
     f_barber%v(2*i-1)=-hh*dz_c(2*i-1)  
     f_barber%v(2*i)=hh*dz_c(2*i)  
    delorb=vr(1,i)*aimag(hh) 

  enddo
 
 i=c_%nd

 delorb= -delorb/vr(1,i)
 
    f_barber%v(2*i-1)=-i_*delorb  *dz_c(2*i-1)  
    f_barber%v(2*i)=  i_*delorb *dz_c(2*i)   

! f_barber is exactly on the resonance and thus commutes
! f_barber could habe beeb chosen using a contruction using vector "a"
! Vector(s) "a" allow for the constructions of nonlinear co-moving maps if one so desires
 
 if(recomputing_h_r) then

 f_comoving=f_comoving+f_barber


rot_c_inverse= exp(-f_comoving)

 
Nc= N*rot_c_inverse



fq=ln(Nc)
Id=n*rot_c_inverse
U_1=rot_c_inverse*n
call c_full_norm_damap(U_1,k0)

u_2=U_1-id
call c_full_norm_damap(U_2,k1)
 write(6,*) "Should be zero ",k1,k1/k0 
 
 else

 fq=fq-f_barber

endif

 
 !!!!!!!!!!!!!!!  Hamiltonian H_r is computed : linear J_x dependence
 
 call d_field_for_demin(fq,H_res)
 

call clean(H_res,H_res,prec=1.e-6_dp,relative=.true.)
call print(H_res)

 
close(mfoo)

!!!!!  This ray creates islands if plotted in the appropriate basis
call kanalnummer(mfoo,"f_as.txt")
    write(mfoo,*) 1,4,C_%NO
   write(mfoo,format6) -0.46057E-03,  0.14002E-05,  0.6461E-04 ,-0.6759E-04
close(mfoo)
 

 

1000 call ptc_end(graphics_maybe=1,flat_file=.false.)


 contains


subroutine  build_lattice(ALS,mis,error,exact,sl,thin,onecell)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: ALS
real(dp),optional :: error(6)
logical, optional :: exact,sl,thin,onecell
real(dp) :: alpha,lbend, cut, ksd, ksf,sig(6)
type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10 
type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,CAVM
type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1 
type(fibre) L27h,sfb
type(layout) :: sfline,sdline,sup1,supb,sup1s
logical  :: mis,thi=.false.,oneperiod
!-----------------------------------
if(present(thin)) thi=thin

call make_states(.true.)
exact_model = .false.;oneperiod = .false.
if(present(exact)) exact_model=exact
if(present(onecell)) oneperiod=onecell
call update_states
madlength = .false.

!old_integrator=-100
call set_mad(energy = 1.5d0, method = 2, step = 1)

madkind2 = matrix_kick_matrix
!madkind2 = drift_kick_drift


L1  = drift("L1 ",  2.832695d0);L2  = drift("L2 ",  0.45698d0);
L3  = drift("L3 ",  0.08902d0);L4  = drift("L4 ",  0.2155d0);
L5  = drift("L5 ",  0.219d0);L6  = drift("L6 ",  0.107078d0);
L7  = drift("L7 ",  0.105716d0);L8  = drift("L8 ",  0.135904d0);
L9  = drift("L9 ",  0.2156993d0);L10 = drift("L10",  0.089084d0);
L11= drift("L11",  0.235416d0);L12= drift("L12",  0.1245d0);
L13= drift("L13",  0.511844d0);L14= drift("L14",  0.1788541d0);
L15= drift("L15",  0.1788483d0);L16= drift("L16",  0.511849d0);
L17= drift("L17",  0.1245d0);L18= drift("L18",  0.235405d0);
L19= drift("L19",  0.089095d0);L20= drift("L20",  0.2157007d0);
L21= drift("L21",  0.177716d0);L22= drift("L22",  0.170981d0);
L23= drift("L23",  0.218997d0);L24 = drift ("L24",  0.215503d0);
L25 = drift ("L25",  0.0890187d0);L26 = drift ("L26",  0.45698d0);
L27 = drift ("L27",  2.832696d0);L27a  = drift (" L27a",  0.8596d0);
L27b  = drift (" L27b",  0.1524d0);L27c  = drift (" L27c",  0.04445d0);
L27d  = drift (" L27d",  1.776246d0);ds  = drift (" DS  ", 0.1015d0);
L27h = drift ("L27",  2.832696d0-0.1d0);
  
QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.13514880435077d0)
QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.46323712526101d0); 
QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  

!!! 1/2 mad-x value
ksf=-41.3355516397069748d0/20;
ksd=56.2564709584745489d0/20;

sf=sextupole ("sf",2.d0*0.1015d0, K2= ksf);
sd= sextupole("sd", 2.d0*0.1015d0, K2= ksd);
!if(skew) call add(sf,-3,1,300.d0)
 VC5=marker("vc5");
ALPHA=0.17453292519943295769236907684886d0;
 
LBEND=0.86621d0;
 
BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
 
CAVM=MARK("CAVM");
CAV=RFCAVITY("CAV",L=0.1000d0,VOLT=-10000.0d0,REV_FREQ=500.0d6)
 sfb=sextupole ("sfb",0.d0, K2= ksf*3.1d0);
call add(sfb,4,0,-25000.d0)
!if(thi) then
 sf=sextupole ("sf",0.d0, K2= ksf*0.203d0);
 sd= sextupole("sd", 0.d0, K2= ksd*0.203d0);
  sfline=(ds+sf+ds);
  sdline=(ds+sd+ds);
!else
! sfline=1*sf;
! sdline=1*sd;
!endif

SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27h+cav;


SUP1s=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27h+sfb+cav;
!           L26+VC5+L27+cavm;

SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND1+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cav;

if(oneperiod) then
 ALS = 3*sup1+sfb;  !11*sup1+supb;
else
 ALS = 11*sup1+supb;
endif
if(present(sl)) then
L1  = drift("L1 ",  2.832695d0);
 if( sl ) then
  Qf1 = QUADRUPOLE(" QF1 ",L=0.d0, K1= 0.01d0 ); L1  = drift("L1 ",L=0.1d0);
  ALS=L1+QF1;
 endif 
endif

ALS = .ring.ALS

call survey(ALS)


if(mis) then
 sig=1.d-5; cut=4.d0; 
 if(present(error)) sig=error
 call MESS_UP_ALIGNMENT(ALS,SIG,cut);
endif
end subroutine build_lattice


end program  



