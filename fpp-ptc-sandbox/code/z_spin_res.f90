program spin_res
  use madx_ptc_module
use pointer_lattice
implicit none
! lesson 1
type(probe) xs0,xs1,xst
type(probe_8) xs
type(layout), pointer :: als
integer mf,mf1,i,k,pos,no,kp,nturn,mfa,ks,kn,j
type(fibre),pointer:: p , fib
type(integration_node),pointer:: tt
type(internal_state),target :: state
real(dp)  prec,cut, closed(6),x(6), ph(4),stu(2),damp(3),thi,nu_mod,circ,epsr,f3r,f3i,emi,xturn
real(dp) spin1,spin2,periodicity
complex(dp) epsb
logical first,breaksym
logical :: mis=.false.,thin=.false.
type(c_damap) id,Q,U_0,U_1,U_2,fp,N,Nc
type(c_taylor) phase(4),nu_spin
type(c_normal_form) c_n
type(c_spinor) isf
type(c_vector_field) fq
type(c_quaternion) qf
type(quaternion) q0,qr
TYPE(c_spinor) N_spin
TYPE(work) werk
real(dp) nus_in,nus_out,dnudn,dnus,pol_rat,dnuda,a_in,a_out,dnda,del_a,nus_shift
character(255) filename
 
 
 integer jj(lnv)
first=.true.;lmax = 10.d0;use_info = .true.;prec=1.d-7;thin=.false.

!ALWAYS_EXACT_PATCHING=.false.
use_quaternion=.true.
! check_excessive_cutting=.false.
call ptc_ini_no_append
 call append_empty_layout(m_u)
 als=>m_u%start
  n_cai=-i_
call build_lu(als)


breaksym=.true.
 



!  call read_lattice_append(M_U,'C:\document\etienne_programs_new\programs_for_learning\fu\Fu_flat_new.txt')
periodicity=1.d0

als=>m_u%start
call MAKE_node_LAYOUT(als)
call survey(als)
 
!goto 1001
 
state=only_4d0+spin0 !+modulation0

werk=als%start


call get_length(als,circ)
write(6,*) "circumference ",circ 

p=>als%start  

if(breaksym) then
call move_to(als,p,"QF2SPIN")
call add(p,2,0,0.01d0)
call add(p,2,0,40.0d0)
!call add(p,2,0,-20.d0)
!call add(p,3,0,-100.d0)
endif
 
write(6,*) p%ag/werk%gamma0I/periodicity  
  
 p=>als%start

do i=1,als%n
 

 p%ag=periodicity*(2.17557d0)*werk%gamma0I
 if(breaksym) p%ag=periodicity*(2.d0+ 0.139870246560962d0)*werk%gamma0I
 

 

p=>p%next

enddo


closed=0
xs0=closed

 
p=>als%start    !%next%next 

werk=p
 closed=0
call FIND_ORBIT_x(closed,state, 1.0e-7_dp, fibre1=p)

write(6,*) check_stable
write(6,format6) closed
 
 
!stop
!goto 1000

write(6,format6) closed
call init(state,3,0)
 
 
call alloc(id,Q,U_0,U_1,U_2,fp,N,Nc)
call alloc(c_n)
call alloc(xs)
call alloc(isf)
call alloc(phase)
call alloc(nu_spin)
call alloc(fq)
call alloc(qf)
call alloc(N_spin) 


xs0=closed
id=1
xs=xs0+id
! call print(xs%ac(1),6 ) ! (F2) ! differs from the first edition   
 
call propagate(xs,state,fibre1=p)
id=xs

 kp=-1
c_n%M=0
c_n%NRES=1
c_n%M(2,1)=1
c_n%ms(1)=kp
c_n%positive=.false.
 

call c_normal(id,c_n, phase=phase, nu_spin=nu_spin,dospin=state%spin)
write(*,*) "yo:L"
N=c_n%atot**(-1)*id*c_n%atot
call print(N)
stop

write(6,*) c_n%tune(1:2)
write(6,*) c_n%spin_tune,c_n%quaternion_angle/pi
 !call print(phase(2))
! call print(nu_spin)
 
!call print(c_n%h_l)
 
!call kanalnummer(i,"spinmap.txt")
!call print(id,i)
!close(i)

!stop
 
N=c_n%atot**(-1)*id*c_n%atot
call clean(N,N,prec=1.d-7)

U_0=N
U_0%q=1.d0
  fq%q%x(2)=-c_n%tune(2)*pi 
 
 U_0=u_0*exp(fq)
 
Nc=Ci_phasor()*N*U_0**(-1)*c_phasor()
fq=ln(Nc)

call c_q0_to_qr(fq%q,qf)
 
call clean(qf,qf,prec=1.d-7)
 call print(qf)
 

  emi=1.d-5
 u_2=c_n%atot*c_phasor()
 
u_1=0
u_1%v(3)=emi
u_1%v(4)=emi

u_2=u_2.o.u_1

x=0
do i=1,4
x(i)=u_2%v(i)
enddo

 
do i=0,3
qf%x(i)=qf%x(i).o.u_1
enddo
nus_shift=2*qf%x(2)
write(6,*) "nus_shift ", nus_shift


epsb=qf%x(1)*2

epsr=abs(epsb)


epsr=epsr**2
write(6,*) "square ", epsr
!call print(qf)
!STOP
pause 888

use_quaternion=.true.
!use_quaternion=.false.

 
write(6,*) " pol_rat "
read(5,*) pol_rat
dnus=.02d0


nus_in= (A_ELECTRON/werk%gamma0I)/periodicity- dnus/2.d0 
nus_out= (A_ELECTRON/werk%gamma0I)/periodicity +dnus/2.d0

a_in= nus_in*werk%gamma0I *periodicity
a_out= nus_out*werk%gamma0I*periodicity 
write(6,*) a_in,a_out
write(6,*) nus_in,nus_out

 dnuda=(nus_out-nus_in)/(a_out-a_in)
 dnudn= -epsr/log((pol_rat+1.d0)/2.d0)/4.d0
write(6,*) dnuda,dnudn
 
dnda=dnuda/dnudn
xturn=dnda*(a_out-a_in)
cut=dnda*(a_out-a_in)


del_a=1.d0/dnda 
nus_shift=nus_shift/del_a 
nturn=xturn
write(6,*) xturn,nturn

pause 777
write(6,format6) x


xs0=closed+x 
xs0%q=1.d0
 
xs=xs0

fib=>als%start
 
 do j=1,als%n
  fib%ag=fib%ag-(a_out-a_in)/2.d0
 fib=>fib%next
 enddo

call kanalnummer(mf,"plot.dat")

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
 
 
 if(use_quaternion) then
 q0=2
qr=xs0%q
 q0=qr*q0*qr**(-1)
 
 write(mf,format6)float(i), float(i)/nturn ,fib%ag/werk%gamma0I/periodicity, q0%x(2)
 write(6,format2)float(i)/nturn,  q0%x(2)

 
 endif
endif
enddo
 
my_estate=>state
close(mf)

 
 
  
1001  call ptc_end(graphics_maybe=1,flat_file=.false.)


contains

subroutine  build_lu(ALS)
use madx_ptc_module
use pointer_lattice
implicit none


type(layout), target :: ALS
 
real(dp) :: alpha,lbend, cut, ksd, ksf,sig(6)

type(fibre) L27h
type(fibre) QF2SPIN,QZ1,b2,beg,ende
 
 

call make_states(.true.)
 
 exact_model=.false.
 
call update_states
madlength = .false.

!old_integrator=-100
call set_mad(p0c = 0.94999986256844970d0, method = 2, step = 1)

madkind2 = matrix_kick_matrix
!madkind2 = drift_kick_drift


ksd=0.d0
QF2SPIN=QUADRUPOLE("QF2SPIN", 1.d-002, K1= ksd)
QZ1=QUADRUPOLE("QZ1", 1.d-002, K1= -125.66366121076938d0)
Lbend=1.d0
alpha=pi/2
B2=sbend("B2", LBEND, ANGLE=ALPHA)


 

 
 beg=marker("start");
 ende=marker("end");

 
 

 
 ALS =beg+ QF2SPIN+4*(QZ1 + B2 +QZ1) +ende;
 

ALS = .ring.ALS

call survey(ALS)


! sig=1.d-5; cut=4.d0; 
! call MESS_UP_ALIGNMENT(ALS,SIG,cut);

end subroutine build_lu
end program spin_res
   
    