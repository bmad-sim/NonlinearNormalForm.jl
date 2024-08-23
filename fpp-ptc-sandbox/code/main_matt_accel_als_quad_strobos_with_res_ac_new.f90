program spin_phase_advance_isf
use madx_ptc_module
use pointer_lattice
implicit none
! lesson 1
type(probe), target :: xs0,xs1,xst
type(probe_8) xs
type(layout), pointer :: als
integer mf,mf1,i,k,pos,no,kp,nturn,mfa,ks,kn,j
type(fibre),pointer:: p , fib
type(integration_node),pointer:: tt
type(internal_state),target :: state
real(dp)  prec,cut, closed(6),x(6), ph(4),stu(2),damp(3),thi,nu_mod,circ,epsr,f3r,f3i,emi,xturn,del,p_res
real(dp) spin1,spin2,mu(4),crossing_position,tunespin,pag0
complex(dp) epsb
logical first
logical :: mis=.false.,thin=.false.,recomputing_h_r
type(c_damap) id,Q,U_0,U_1,U_2,N,Nc,one_turn_map,rot_c_inverse,n_isf
type(c_taylor) phase(4),nu_spin,delspin,delorb
type(c_normal_form) c_n
type(c_spinor) isf
type(c_vector_field) fq,f_comoving,n_isvf
type(c_quaternion) qf
type(quaternion) q0,qr
TYPE(c_spinor) N_spin
TYPE(work) werk
real(dp) nus_in,nus_out,dnudn,dnus,pol_rat,dnuda,a_in,a_out,dnda,del_a,nus_shift,k0,k1
character(255) filename
real(dp), target ,allocatable ::  vr(:,:),vo(:,:),vecr(:) 
real(dp), target  :: mr  
     type(spinor) n_stroboscopic 
type (c_linear_map)  q_rot, q_as,q_cs,q_map,q_orb,q_u,q_u_inv
 
 
type(c_lattice_function) lf
 integer jj(lnv)
first=.true.;lmax = 10.d0;use_info = .true.;prec=1.d-7;thin=.false.

!ALWAYS_EXACT_PATCHING=.false.
use_quaternion=.true.
! check_excessive_cutting=.false.
call ptc_ini_no_append
 call append_empty_layout(m_u)

 als=>m_u%start
  n_cai=-i_

! call append_empty_layout(m_u)
!ALS=>m_u%start

call build_lattice_als(ALS,exact=.false.,thin=.true.,onecell=.FALSE.) 

 ! call read_lattice_append(M_U,'C:\document\etienne_programs_new\programs_for_learning\fu\ptc.lattice')
  
als=>m_u%start
call MAKE_node_LAYOUT(als)
call survey(als)
 
 

 call get_length(als,circ)
write(6,*) "circumference ",circ
! call read_ptc_command77("C:\document\etienne_programs_new\programs_for_learning\fu\AC_modulation_2.txt")

!!!!Fitting the tune to nu_x=0.334 !!!! 
 call kanalnummer(mf,file="AC_modulation.txt")
 write(mf,*)"select layout "
 write(mf,*)"          1 "
 write(mf,*)" MODULATE "
 write(mf,*)" QF2SPIN   1" 
 write(mf,*)"1.d0 0 0    "   !DC_ac,A_ac,theta_ac
 write(mf,*)"1.d0   2    "   ! D_ac,n_ac  
 write(mf,*)"1 0 0.0005d0     "    ! n d_bn(n) d_an(n)  0.0001d0 
 write(mf,*)"0  0 0 "
 write(mf,*)" return "
 close(mf)

 call read_ptc_command77("AC_modulation.txt")
 
  
 
state=only_4d0+spin0+modulation0

werk=als%start

p=>als%start  

write(6,*) p%ag/werk%gamma0I 

no=5
pag0=p%ag
pag0=pag0*1.0005d0
 p=>als%start

do i=1,als%n
  p%ag=pag0
p=>p%next
enddo
call get_length(als,circ)
write(6,*) "circumference ",circ 

 

!endif
 

 
 
 nu_mod=0.404074285438127d0
!!!! set a modulation clock !!!!!!
xs0%ac(1)%om=twopi*nu_mod/circ ! (B1) ! differs from the first edition   
xs0%ac(1)%x=0.d0 ;       ! (B2) ! differs from the first edition   
write(6,*) " Modulation tune in radians =",circ*xs0%ac(1)%om
xs0%use_q=.true.

xs1=xs0
xs1%ac(1)%x=0
xs1%ac(1)%x(1)=1.d0
xs1%ac(1)%om=twopi*nu_mod/circ ! (B1) ! differs from the first edition   
xs1%nac=1
write(6,*) xs1%ac(1)%x
my_eprobe=>xs1

write(6,*) xs0%nac
xs0%nac=1

 
closed=0
xs0=closed

 
p=>als%start    !%next%next 

werk=p
 closed=0
call FIND_ORBIT_x(closed,state, 1.0e-7_dp, fibre1=p)

write(6,*) check_stable
write(6,*) closed

 
!stop
!goto 1000

write(6,format6) closed
call init(state,no,0)
  
call alloc(id,Q,U_0,U_1,U_2,N,Nc,one_turn_map,rot_c_inverse,n_isf)
call alloc(c_n)
call alloc(xs)
call alloc(isf)
call alloc(phase)
call alloc(nu_spin)
call alloc(fq,f_comoving,n_isvf)
call alloc(qf)
call alloc(phase)
call alloc(nu_spin,delspin,delorb) 
call alloc(N_spin) 


xs0=closed
id=1
xs=xs0+id
 
 
call propagate(xs,state,fibre1=p)
one_turn_map=xs

!!  
c_n%M=0
c_n%NRES=1
c_n%M(3,1)=1
c_n%ms(1)=-1
c_n%positive=.false.
 
!!!!!!!!!!! Gram-Schmidt computation of perpendicular planes !!!!!!!!!!!!!

 allocate(vr(c_%nd+1,c_%nd+1),vo(c_%nd,c_%nd+1),vecr(c_%nd+1) )
vr=0
vo=0

vr(1,1:c_%nd) = c_n%m(1:c_%nd,1)
vr(1,c_%nd+1) = c_n%ms(1)
 
mr=(dot_product(vr(1,1:c_%nd+1),vr(1,1:c_%nd+1)) )**(0.5_dp)
write(6,*) mr
call gramschmidt(vr)

!! vo is J_1 and J_2 in terms of Jr and Ja

 vo=transpose(vr)

!!!!  The vectors are not normalized following the convention of section 5.4.3 of blue book
!!!!  or 5.3 of yellow book

 vr=mr*vr
 vo=vo/mr
pause 123
write(6,format4) vr(1,1:c_%nd+1) 
write(6,format4) vr(2,1:c_%nd+1) 
write(6,format4) vr(3,1:c_%nd+1) 
write(6,format4) vr(4,1:c_%nd+1) 
write(6,*)(dot_product(vr(1,1:c_%nd+1),vr(1,1:c_%nd+1)) )
write(6,*)(dot_product(vr(2,1:c_%nd+1),vr(2,1:c_%nd+1)) )
write(6,*)(dot_product(vr(3,1:c_%nd+1),vr(3,1:c_%nd+1)) )
write(6,*)(dot_product(vr(4,1:c_%nd+1),vr(4,1:c_%nd+1)) )
pause 678
 
!call c_normal(id,c_n, phase=phase, nu_spin=nu_spin,dospin=state%spin)
call c_normal(one_turn_map,c_n, phase=phase, nu_spin=nu_spin,dospin=state%spin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,c_%nd
 mu(i)=phase(i)  
enddo
mu(4)=  nu_spin  
write(6,*) " nus including spin "  
write(6,*) mu
 
tunespin= A_ELECTRON/werk%gamma0I-nint(A_ELECTRON/werk%gamma0I)
 
 
write(6,*) " x,y and spin tunes "
write(6,format4) c_n%tune(1:3), c_n%spin_tune
write(6,format4) mu
write(6,*) " term multipliying the quaternion  e_2 in units of pi "
write(6,*)  c_n%quaternion_angle/pi

 
!!!! co-moving  
 N=c_n%atot**(-1)*one_turn_map*c_n%atot
N=Ci_phasor()*N*c_phasor()
!!!! Computation of the inverse co-moving rotations

!!!!   Theory identical to transverse except for the obvious
!!!!   quaternion operator   :-J_spin: = -1/2 e_2 
!!!! notice that this is similar to J_x= (x^2+p^2)/2
!!!!  so e_2 is truly like the Courant-Snyder invariant
 

f_comoving=0

do k=1,c_%nd+1
del=0
 do ks=1,c_%nd+1
  del= vr(k,ks)*mu(ks)+del
 enddo
if(k==1) then
 p_res=nint(del)
 del=p_res*twopi/mr**2
else
 del=del*twopi/mr**2
endif

   do i=1,c_%nd
    f_comoving%v(2*i-1)=-i_*del* vr(k,i)  *dz_c(2*i-1)  + f_comoving%v(2*i-1)
    f_comoving%v(2*i)=  i_*del*  vr(k,i)*dz_c(2*i)  +   f_comoving%v(2*i)
  enddo
  f_comoving%q%x(2)=-vr(k,c_%nd+1)*del/2 + f_comoving%q%x(2)

enddo


write(6,*) vr(1,1)*mu(1)+vr(1,2)*mu(2)+vr(1,3)*mu(3)+c_n%ms(1)*mu(4)


 rot_c_inverse= exp(-f_comoving)

Nc= N*rot_c_inverse

 
  
fq=ln(Nc)

call clean(fq,fq,prec=1.d-4,relative=.true.)
call clean(fq,fq,prec=1.d-7 )

!call print(fq)


delorb=0.0_dp
do i=1,c_%nd
  delorb=vr(1,i)*aimag(fq%v(2*i).k.(2*i))+delorb
enddo 

!!! uncomment below and the ISF is slightly wrong as well as the prediction
!!! for the resonance crossing
!!!  since the Barber H_r is only correct in leading order

! delorb = delorb.sub.'0'     


!!!!  This is equation 6.54 of the blue book
recomputing_h_r=.true.
if(recomputing_h_r) then
  do i=1,c_%nd
    f_comoving%v(2*i-1)=fq%v(2*i-1) +f_comoving%v(2*i-1)
    f_comoving%v(2*i)= fq%v(2*i)  + f_comoving%v(2*i)
  enddo
 
 delspin= -delorb/c_n%ms(1)
 f_comoving%q%x(2)=-delspin/2 + f_comoving%q%x(2)
 
 
rot_c_inverse= exp(-f_comoving)

Nc= N*rot_c_inverse

 
fq=ln(Nc)
Id=n*rot_c_inverse
U_1=rot_c_inverse*n
call c_full_norm_damap(U_1,k0)

u_2=U_1-id
call c_full_norm_damap(U_2,k1)
 write(6,*) "Should be zero ",k1,k1/k0 
call clean(fq,fq,prec=1.d-10,relative=.true.)
 pause 7654

write(6,*) " Barber H_r : THE gadget "
call c_q0_to_qr(fq%q,qf)
 else
write(6,*) " Barber H_r : THE gadget "
 delspin= -delorb/c_n%ms(1)
call c_q0_to_qr(fq%q,qf)
 qf%x(2)= delspin/2 + qf%x(2)

endif

   
 

  

  call clean(qf,qf,prec=1.d-4 ,relative=.true.)
  call clean(qf,qf,prec=1.d-10  )
!call print(qf)


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! initial conditions in the laboratory frame (around closed orbit) !!!!!!!!!!!!!!!

x=0

x(c_%nd2-1)=1  ! clock arms set
x(c_%nd2)=0




do i=1,c_%nd2
u_1%v(i)=x(i)
enddo

 !!!! put Barber H_r in laboratory basis
!!!!! and evaluatec
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

U_0=c_n%atot*c_phasor()
U_2=U_0**(-1).o.U_1
 
do i=0,3
qf%x(i)=qf%x(i).o.u_2
enddo
 nus_shift= qf%x(2)/pi
 
epsb=qf%x(1)*2


epsr=abs(epsb)

epsr=epsr**2
write(6,*) "epsilon^2 ", epsr
write(6,*) "Tune Shift ", nus_shift
 
call c_qr_to_q0(qf,qf)

!!! Evaluating the ISF in the normalized frame
!!! via Barber's formulation (Exact result)

  damp(1)= qf%x(1)  
  damp(2)= qf%x(2)
  damp(3)= qf%x(3)   
 write(6,*) damp

 damp=damp/sqrt(damp(1)**2+damp(2)**2+damp(3)**2)
 
 if(.false.) then

n_isf=1
do i=1,3
 n_isf%q%x(i)=damp(i)
enddo

n_isf=c_n%atot*n_isf*c_n%atot**(-1)

!!!!  Moving the ISF into the laboratory variables
!!!!  as computed by stroboscopic average 


 n_isf%q=n_isf%q.o.u_1
do i=1,3
 damp(i)=n_isf%q%x(i)
enddo
else
n_isvf=0
do i=1,3
 n_isvf%q%x(i)=damp(i)
enddo

!n_isvf=c_n%atot*n_isvf*c_n%atot**(-1)
n_isvf= c_n%atot**(-1)*n_isvf

!!!!  Moving the ISF into the laboratory variables
!!!!  as computed by stroboscopic average 


 n_isvf%q=n_isvf%q.o.u_1
do i=1,3
 damp(i)=n_isvf%q%x(i)
enddo
endif
 damp=damp/sqrt(damp(1)**2+damp(2)**2+damp(3)**2)
  
if(damp(2)<0) damp=-damp
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use_quaternion=.true.
 write(6,format6) x
  


 
write(6,*) " pol_rat "
read(5,*) pol_rat
dnus=.02d0


nus_in= (A_ELECTRON/werk%gamma0I)- dnus/2.d0 
nus_out= (A_ELECTRON/werk%gamma0I) +dnus/2.d0
write(6,*) "nus_in,nus_out"

write(6,*) nus_in,nus_out
pause 666
a_in= nus_in*werk%gamma0I
a_out= nus_out*werk%gamma0I
write(6,*) a_in,a_out
write(6,*) nus_in,nus_out

 dnuda=(nus_out-nus_in)/(a_out-a_in)
 dnudn= -epsr/log((pol_rat+1.d0)/2.d0)/4.d0
write(6,*) dnuda,dnudn

 
dnda=dnuda/dnudn
xturn=dnda*(a_out-a_in)
cut=dnda*(a_out-a_in)


del_a=1.d0/dnda 
!nus_shift=nus_shift/del_a 
nturn=xturn
write(6,*)"# of turns ", xturn,nturn


x=0
xs0=closed+x 
xs0%ac(1)%om=twopi*nu_mod/circ ! (B1) ! differs from the first edition   
xs0%ac(1)%x=0.d0 ;       ! (B2) ! differs from the first edition   
xs0%ac(1)%x(1)=1.d0 ;       ! (B2) ! differs from the first edition 
xs0%ac(1)%x(2)=0.d0 ;       ! (B2) ! differs from the first edition 
xs0%q=1.d0
xs0%nac=1
 


xs=xs0

fib=>als%start
 
 do j=1,als%n
  fib%ag=fib%ag-(a_out-a_in)/2.d0
 fib=>fib%next
 enddo

call kanalnummer(mf,"plot.dat")


write(6,*) " Tune shift ", nus_shift 
crossing_position=nus_shift+tunespin
write(6,*) " position of crossing ", crossing_position

 first=.true.

p=>als%start
do i=1,nturn
fib=>als%start
 
 do j=1,als%n
  fib%ag=fib%ag+del_a
 fib=>fib%next
 enddo



call propagate(xs0,state,fibre1=p)
if(.not.check_stable) stop 777
if(mod(i,nturn/100)==0) then
!write(6,*) i, float(i)/nturn 

!write(6,format4) xs0%q%x(0:3)
  q0=2
qr=xs0%q
 q0=qr*q0*qr**(-1)
if(fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I)> crossing_position.and.first) then
first=.false.
Write(6,*) " Resonance crossing predicted around here "
endif
 write(mf,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I), q0%x(2)
 write(6,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I), q0%x(2)

 
endif
enddo
 
my_estate=>state
close(mf)

 p=>als%start

do i=1,als%n
  p%ag=pag0
p=>p%next
enddo

 
xs0=closed+x 
xs0%ac(1)%om=twopi*nu_mod/circ ! (B1) ! differs from the first edition   
xs0%ac(1)%x=0.d0 ;       ! (B2) ! differs from the first edition   
xs0%ac(1)%x(1)=1.d0 ;       ! (B2) ! differs from the first edition 
xs0%ac(1)%x(2)=0.d0 ;       ! (B2) ! differs from the first edition 
xs0%q=1.d0
xs0%nac=1
 


call stroboscopic_average(als,xs0,xst,1,state,20000,5000,n_stroboscopic,6)
if(n_stroboscopic%x(2)<0) n_stroboscopic%x=-n_stroboscopic%x

write(6,*) " numerical stroboscopic average"

write(6,*) n_stroboscopic%x

write(6,*) "We leave a resonance in the map computation "

write(6,*) damp(1:3)
call kill(c_n)
call alloc(c_n)
!!!!!!!!!!!!!!!!!!!!!!!!! 
call c_normal(one_turn_map,c_n, phase=phase, nu_spin=nu_spin,dospin=state%spin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! No resonances left in the map so the normalized ISF is e_2.
!!! two ways to get the isf 
n_isf=1
n_isf%q=2
n_isvf=0
n_isvf%q=2

n_isf=c_n%atot*n_isf*c_n%atot**(-1)
n_isvf=c_n%atot**(-1)*n_isvf
  
 n_isf%q=n_isvf%q.o.u_1
do i=1,3
 damp(i)=n_isf%q%x(i)
enddo
 damp=damp/sqrt(damp(1)**2+damp(2)**2+damp(3)**2)
 
if(damp(2)<0) damp=-damp

write(6,*) " leave no resonance resonance  "
write(6,*) damp
 
my_estate=>state


1000  call ptc_end(graphics_maybe=1,flat_file=.false.)

contains


subroutine  build_lattice_als(ALS,error,exact,sl,thin,onecell)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: ALS
real(dp),optional :: error(6)
logical, optional :: exact,sl,thin,onecell
real(dp) :: alpha,lbend, cut, ksd, ksf 
type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10 ,QF2spin
type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,CAVM
type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1 
type(layout) :: sfline,sdline,sup1,supb
logical(lp) :: thi=.false.,oneperiod
!-----------------------------------
if(present(thin)) thi=thin

call make_states(.true.)
exact_model = .false.;oneperiod = .false.
if(present(exact)) exact_model=exact
if(present(onecell)) oneperiod=onecell
call update_states
madlength = .false.

old_integrator=-10
call set_mad(energy = 1.5d0, method = 2, step = 1)

madkind2 = matrix_kick_matrix


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

QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+6.447435260914397D-03)
QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.3368D0-2.593018157427161D-02); 
QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  
QF2spin= QUADRUPOLE("QF2SPIN",0.344D0, K1= 2.2474D0)
!call add(QF2spin,4,0,1000000.d0)
!!! 1/2 mad-x value
ksf=-41.3355516397069748d0;
ksd=56.2564709584745489d0;

sf=sextupole ("sf",2.d0*0.1015d0, K2= ksf);
sd= sextupole("sd", 2.d0*0.1015d0, K2= ksd);


 VC5=marker("vc5");
ALPHA=0.17453292519943295769236907684886d0;
 
LBEND=0.86621d0;
 
BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
 
CAVM=MARK("CAVM");
CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1.0d0,REV_FREQ=500.0d6)

if(thi) then
 sf=sextupole ("sf",0.d0, K2= ksf*0.203d0);
 sd= sextupole("sd", 0.d0, K2= ksd*0.203d0);
!call add(sf,4,0,1000.d0)
!call add(sd,4,0,1000.d0)
  sfline=(ds+sf+ds);
  sdline=(ds+sd+ds);
else
 sfline=1*sf;
 sdline=1*sd;
endif

SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cavm;

SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND1+L21+&
           L22+QD2+L23+L24+QF2spin+L25+ &
           L26+VC5+L27+cav;

if(oneperiod) then
 ALS = sup1;  !11*sup1+supb;
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

 
end subroutine build_lattice_als



end program spin_phase_advance_isf




