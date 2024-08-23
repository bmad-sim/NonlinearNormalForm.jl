program spin_resonance_isf
use madx_ptc_module
use pointer_lattice
implicit none
! lesson 1
type(probe) xs0,xs1,xst
type(probe_8) xs
type(layout), pointer :: ring
integer mf,mf1,i,k,pos,no,kp,nturn,mfa,ks,kn,j
type(fibre),pointer:: p , fib
type(integration_node),pointer:: tt
type(internal_state),target :: state
real(dp)  prec,cut, closed(6),x(6), ph(4),stu(2),damp(3),thi,nu_mod,circ,epsr,f3r,f3i,emi,xturn
real(dp) mu(3),del,tunespin,crossing_position,k1,k0
complex(dp) epsb
logical first,recomputing_h_r
logical :: mis=.false.,thin=.false.
type(c_damap) id,one_turn_map,U_0,U_1,U_2,fp,N,Nc,n_isf,rot_c_inverse
type(c_taylor) phase(4),nu_spin,delspin,delorb
type(c_normal_form) c_n
type(c_spinor) isf
type(c_vector_field) fq,f_comoving,n_isvf
type(c_quaternion) qf
type(quaternion) q0,qr
TYPE(c_spinor) N_spin
TYPE(work) werk
real(dp) nus_in,nus_out,dnudn,dnus,pol_rat,dnuda,a_in,a_out,dnda,del_a,nus_shift,p_res
character(255) filename
real(dp), target ,allocatable ::  vr(:,:),vo(:,:) 
real(dp), target  :: mr 
    type(spinor) n_stroboscopic 
 integer jj(lnv)
 
first=.true.;lmax = 10.d0;use_info = .true.;prec=1.d-7;thin=.false.

 
use_quaternion=.true.
 
call ptc_ini_no_append
 call append_empty_layout(m_u)
 ring=>m_u%start
  n_cai=-i_
call build_fu(ring)

no=5

ring=>m_u%start
call MAKE_node_LAYOUT(ring)
call survey(ring)
 
!goto 1001
 
state=only_4d0+spin0 

werk=ring%start


call get_length(ring,circ)
write(6,*) "circumference ",circ 

p=>ring%start  
 

 p=>ring%start
tunespin= 0.159169892732279d0 
do i=1,ring%n
 p%ag=(2.d0+ tunespin)*werk%gamma0I
if(p%mag%name(1:2)=="B2") then
prec=i*0.01d0
call add(p,4,1,-4000.d0)
call add(p,2,1,prec)
endif
p=>p%next
enddo


closed=0
xs0=closed

 
p=>ring%start   

werk=p
 closed=0
call FIND_ORBIT_x(closed,state, 1.0e-7_dp, fibre1=p)

write(6,*) check_stable
write(6,format6) closed
 

write(6,format6) closed
call init(state,no,0)
 
 
call alloc(id,one_turn_map,U_0,U_1,U_2,fp,N,Nc,n_isf,rot_c_inverse)
call alloc(c_n)
call alloc(xs)
call alloc(isf)
call alloc(phase)
call alloc(fq,f_comoving,n_isvf)
call alloc(nu_spin,delspin,delorb)
call alloc(qf)
call alloc(N_spin) 
call alloc(isf) 

xs0=closed
id=1
xs=xs0+id
 
call propagate(xs,state,fibre1=p)
one_turn_map=xs

!!!!!!!!!!!Resonance left in the normalized  map !!!!!!!!!!!!!
 kp=-1
c_n%M=0
c_n%NRES=1
c_n%M(2,1)=1
c_n%ms(1)=kp
c_n%positive=.false.

!!!!!!!!!!! Gram-Schmidt computation of perpendicular planes !!!!!!!!!!!!!

 allocate(vr(c_%nd+1,c_%nd+1),vo(c_%nd,c_%nd+1)  )
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

!!!!!!!!!!!! Normalization  !!!!!!!!!!!!! 
call c_normal(one_turn_map,c_n, phase=phase, nu_spin=nu_spin,dospin=state%spin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,c_%nd
 mu(i)=phase(i)  
enddo
mu(3)=  nu_spin  
write(6,*) " nus including spin "  
write(6,*) mu
 
write(6,*) " x,y and spin tunes "
write(6,*) c_n%tune(1:2), c_n%spin_tune
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



 rot_c_inverse= exp(-f_comoving)

Nc= N*rot_c_inverse

fq=ln(Nc)

!!!!!!!!!!!!!  up to this point  the theory was like orbital!!!!!!!!!!!    
 

!!!!!!!!!!!!!   fq is H_r !!!!!!!!!!!    
!!!!!!!!!!!!!   However it has some rotations in the orbital !!!!!!!!!!!    
!!!!!!!!!!!!!   It is possible to move this into the spin !!!!!!!!!!!    
!!!!!!!!!!!!!   via another co-moving rotation  !!!!!!!!!!!    
!!!!!!!!!!!!!   via another co-moving rotation  !!!!!!!!!!!    
!!!!!!!!!!!!!   This is equation 6.54 explained in a sloppy way in the blue book !!!!!!!!!!
 
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
 

write(6,*) " Barber H_r : THE gadget "
call c_q0_to_qr(fq%q,qf)
 else
write(6,*) " Barber H_r : THE gadget "
 delspin= -delorb/c_n%ms(1)
call c_q0_to_qr(fq%q,qf)
 qf%x(2)= delspin/2 + qf%x(2)

endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! initial conditions in the laboratory frame (around closed orbit) !!!!!!!!!!!!!!!

x=0
x(3)=9.3676d-005
x(4)=-1.9044d-005
  
 x=x*50  
!x=x*10  

do i=1,4
u_1%v(i)=x(i)
enddo

 !!!! put Barber H_r in laboratory basis
!!!!! and evaluate
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
 
 
 

!!!!  Moving the ISF into the laboratory variables
!!!!  as computed by stroboscopic average 


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

 n_isf%q=n_isf%q.o.u_1
do i=1,3
 damp(i)=n_isf%q%x(i)
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

!!!! fake acceleration computed to achieve desired polarization
!!!  via Froissart-Stora

nus_in= (A_ELECTRON/werk%gamma0I)- dnus/2.d0 
nus_out= (A_ELECTRON/werk%gamma0I) +dnus/2.d0

a_in= nus_in*werk%gamma0I 
a_out= nus_out*werk%gamma0I 

 dnuda=(nus_out-nus_in)/(a_out-a_in)
 dnudn= -epsr/log((pol_rat+1.d0)/2.d0)/4.d0
 
dnda=dnuda/dnudn
xturn=dnda*(a_out-a_in)
cut=dnda*(a_out-a_in)


del_a=1.d0/dnda 
nturn=xturn

Write(6,*) " Number of turns during fake acceleration "
write(6,*)  nturn

xs0=closed+x 
xs0%q=1.d0
 
xs=xs0
 
fib=>ring%start
 
 do j=1,ring%n
  fib%ag=fib%ag-(a_out-a_in)/2.d0
 fib=>fib%next
 enddo

call kanalnummer(mf,"plot.dat")


write(6,*) " Tune shift ", nus_shift 
crossing_position=nus_shift+tunespin
write(6,*) " position of crossing ", crossing_position
first=.true.

p=>ring%start
do i=1,nturn
fib=>ring%start
 
 do j=1,ring%n
  fib%ag=fib%ag+del_a
 fib=>fib%next
 enddo

 

call propagate(xs0,state,fibre1=p)
if(.not.check_stable) stop 777
if(mod(i,nturn/100)==0) then
 
 
 if(use_quaternion) then
 q0=2
qr=xs0%q
 q0=qr*q0*qr**(-1)
 
if(fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I)> crossing_position.and.first) then
first=.false.
Write(6,*) " Resonance crossing predicted around here "
endif
! write(mf,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I, q0%x(2)
 write(mf,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I), q0%x(2)
 write(6,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I-int(fib%ag/werk%gamma0I), q0%x(2)

 
 endif
endif
enddo
 
my_estate=>state
close(mf)

1221 continue 
 p=>ring%start

do i=1,ring%n

  p%ag=(2.d0+ tunespin)*werk%gamma0I
p=>p%next
enddo

 
 
xs0=closed+x 
xs0%q=1.d0
call stroboscopic_average(ring,xs0,xst,1,state,300000,50000,n_stroboscopic,6)
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

n_isf=1
n_isf%q=2

n_isf=c_n%atot*n_isf*c_n%atot**(-1)

 n_isf%q=n_isf%q.o.u_1
do i=1,3
 damp(i)=n_isf%q%x(i)
enddo
 damp=damp/sqrt(damp(1)**2+damp(2)**2+damp(3)**2)
 
if(damp(2)<0) damp=-damp

write(6,*) " leave no resonance resonance  "
write(6,*) damp
 
my_estate=>state
close(mf)

 
 stop 333
  
2222 continue
1001  call ptc_end(graphics_maybe=1,flat_file=.false.)


contains

subroutine  build_fu(ring)
use madx_ptc_module
use pointer_lattice
implicit none
type(layout), target :: ring
real(dp) :: alpha,lbend, cut, ksd, ksf,sig(6)
type(fibre) L27h
type(fibre)  QZ1,b2,beg,ende

call make_states(.true.)
 
 exact_model=.false.
 
call update_states
madlength = .false.

!old_integrator=-100
call set_mad(p0c = 0.94999986256844970d0, method = 2, step = 1)

madkind2 = matrix_kick_matrix


ksd=0.d0
 QZ1=QUADRUPOLE("QZ1", 1.d-002, K1= -125.66366121076938d0)
Lbend=1.d0
alpha=pi/2
B2=sbend("B2", LBEND, ANGLE=ALPHA)


 

 
 beg=marker("start");
 ende=marker("end");

 
 

 
 ring =beg+4*(QZ1 + B2 +QZ1) +ende;
 

ring = .ring.ring

call survey(ring)
 
end subroutine build_fu

end program spin_resonance_isf




