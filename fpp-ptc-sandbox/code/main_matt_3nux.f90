program Resonance
use pointer_lattice 
use gauss_dis
implicit none

! bmad variables 
!type (lat_struct) lat0

! ptc variables 
type (fibre), pointer :: ptc_fibre ,ptc_oct
type(layout), pointer :: ptc_lattice
type(internal_state), target :: state
type(c_damap):: T,id,A,N_c,rot,u,uc,fr1,fr2,n,rot_c_inverse,U_1,U_2
type(c_normal_form) :: c_n
type(probe_8) :: ray_8,xs
type(probe)  :: ray_closed,ray,xs0,xst 
real(dp) closed_orbit(6),r,phi0,phase(3),om,g,ang(3),Jy,k1,k0
real(dp) jx(2),discriminant,cnu3,mu(3),del
real(dp) x1(6),x2(6),x1k(6),x2k(6),sigm(6),cut 
complex(dp) gamma 
integer i,k,n_arg,mf,mfoo,nturns,nmax0 ,j,mfr,no,p_res
character(100) lat_file
type(c_universal_taylor)   H_res,F_res,H_bengtsson
type(c_ray) f1,f2
complex(dp) dx,c3nu,alpha1,alpha2,alpha3,alpha_xx_half,alphaoct,f21,ld
complex(dp) f40_2,f40,f40_1,g3h,g3,f31_1,f31_2,f40_oct,f31_oct,f31,sylee,bengtsson
logical use_vector_field 
type(c_universal_taylor), target ::   a_tot_inverse(6),a_tot(6)
integer jinv(2),npos 
real(dp),allocatable ::  bet(:),bet2(:),pha(:),spos(:)
complex(dp) ,allocatable ::f3_1(:),f3_1t(:),f3_3(:)
type(c_ray) fix1,fix2,fix1k,fix2k
type(c_vector_field)ft ,ft3,ft4,fq,f_comoving
 
integer  jl,je(6)
real(dp) s,circ,mu_tot,mr2 
type(c_taylor) phat(3),hh
logical :: skew=.false.,mis=.false.
real(dp), target ,allocatable ::  vr(:,:),vo(:,:),vi1(:),vi2(:) 
real(dp), target  :: mr
 
 preffile= "C:\document\etienne_programs_new\programs_for_learning\park\pref3nux_als.txt"
 
!c_lda_used =1500
 
call kanalnummer(mf,"C:\document\etienne_programs_new\programs_for_learning\suntao\results_rad.txt")

call ptc_ini_no_append ! initializes PTC
 call append_empty_layout(m_u)
 ptc_lattice=>m_u%start
 call build_lattice(ptc_lattice,mis,exact=.false.,thin=.true.,onecell=.true.) 
 
write(6,*) " 1-d-f or 2-d-f "
read(5,*) i
if(i==1) then
 state=only_2d0 
else
 state=only_4d0 
endif
my_estate=>state

 

!!! next, use PTC normal form to compute intrinsic resonance strength
ptc_lattice=>m_u%start
ptc_fibre=>ptc_lattice%start

no=3
write(6,*) "order >= 3 (try 7th order) "
read(5,*) no
if(no<3) no=3
!! Counting numbers of multipoles 

ptc_fibre=>ptc_lattice%start
npos=0
 
 
use_quaternion=.true.
 n_cai=-i_
ptc_oct=>ptc_fibre
 
  
my_estate=>state

 ptc_fibre=>ptc_lattice%start
 
do i=1, ptc_lattice%n 

if(ptc_fibre%mag%p%nmul>=3) then

 if(ptc_fibre%mag%bn(3)>0) then
   call add(ptc_fibre,4,1,0.0_dp)
  else
   call add(ptc_fibre,4,1,0.0_dp)
 endif
endif
 ptc_fibre=>ptc_fibre%next
enddo
 
call init(state, no, 0)

call alloc(T,id,N_c,rot,a,rot_c_inverse,n,U_1,U_2)
call alloc(xs)
call alloc(ray_8)
call alloc(c_n)
call alloc(rot)
call alloc(ft,ft3,ft4,fq,f_comoving)
call alloc(phat)
call alloc(hh)
 

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

 
  c_n%positive=.false.
  c_n%nres=(c_%no+1)/3
  write(6,*) " # of resonance terms ",c_n%nres
  do i=1, (c_%no+1)/3
   c_n%m(1,i)=3*i ! 1/3 order 
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
 write(6,*) " mr2 ",mr2
mr=sqrt(mr2)
 
  vr=mr*vr
 vo=(1.d0/mr)*vo
 write(6,*) "*************************************"
 if(c_%nd==1) then 
   write(6,*) vr(1,1)
   write(6,*) vo(1,1)
 else
 write(6,*) vr(1,1:2)
 write(6,*) vr(2,1:2)
 write(6,*) vo(1,1:2)
 write(6,*) vo(2,1:2)
endif

 

 
call c_normal(T,c_n,canonize=.true.,phase=phat)

write(6,*) c_n%tune(1:2)
 !!!!!!!!!!!!!!!
 

do i=1,c_%nd
mu(i)=phat(i)
enddo

 !!!! for windows graphical interface !!!!!!!


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

  
!!!!!!!!!!!!!!!!!!!!!!

 
write(6,*) " mr2  ",mr2
 

!!!! Computation of the inverse co-moving rotations

del=0
do i=1,c_%nd
del=del+vr(1,i)*mu(i) 
enddo

 p_res=nint(del)
 
 del=p_res*twopi/mr**2

  do i=1,c_%nd
  f_comoving%v(2*i-1)=-i_*del* vr(1,i)  *dz_c(2*i-1) 
  f_comoving%v(2*i)= i_*del*  vr(1,i)*dz_c(2*i)   
enddo

if(c_%nd>=2) then

 del= (vr(2,1)*mu(1)+vr(2,2)*mu(2))
  del=del*twopi/mr**2

   do i=1,c_%nd
    f_comoving%v(2*i-1)=-i_*del* vr(2,i)  *dz_c(2*i-1) + f_comoving%v(2*i-1)
    f_comoving%v(2*i)=  i_*del*  vr(2,i)*dz_c(2*i)     + f_comoving%v(2*i)
  enddo

endif
 N=c_n%atot**(-1)*T*c_n%atot
n=ci_phasor()*n*c_phasor()
 
rot_c_inverse=exp(-f_comoving)

 N_c=n*rot_c_inverse
  

  
Id=n*rot_c_inverse
U_1=rot_c_inverse*n
call c_full_norm_damap(U_1,k0)

u_2=U_1-id
call c_full_norm_damap(U_2,k1)
 write(6,*) "Should be zero ",k1,k1/k0

 
 fq=c_logf_spin(N_c)

   call clean(fq,fq,prec=1.e-7_dp,relative=.true.)

call print(fq)
 

 
 write(6,*) " End of Vector field "
 
 call d_field_for_demin(fq,H_res)
 
  call clean(H_res,H_res,prec=1.e-5_dp,relative=.true.)

 write(6,*); write(6,*) " -2*pi*H_r "; write(6,*)

call print(H_res)
 
 write(mf,*); write(mf,*) " -2*pi*H_r "; write(mf,*)

 call print(H_res,mf)


 !!!!!!!!!!!!!!   Computation of fixed points: analytical calculation to order octupoles
je=0
je(1)=3

gamma= h_res.sub.je(1:c_%nd2)
 
je=0;je(1)=2;je(2)=2;
alpha_xx_half= h_res.sub.je(1:c_%nd2)
je=0;je(1)=1;je(2)=1;
del=h_res.sub.je(1:c_%nd2)
 !  Using H_res 

  phi0=atan2(aimag(gamma),real(gamma) )
  g=sqrt(aimag(gamma)**2+real(gamma)**2)
  write(6,*) g*exp(i_*phi0)

discriminant=1.d0-8.d0*alpha_xx_half*del/9.d0/g**2
 
jx(1)=-(3.d0*g/4.d0/alpha_xx_half)*(1.d0-sqrt(discriminant))
jx(2)=(3.d0*g/4.d0/alpha_xx_half)*(1.d0+sqrt(discriminant))


 write(6,*) " alpha_xx_half",alpha_xx_half
 write(6,*) " del",del
 write(6,*) " g",g
 write(6,*) " discriminant",discriminant
write(6,*) jx(1),jx(2)
 

 id=0
 fix1=0
 fix2=0
 fix1%x(1)=jx(1)*exp(-i_*(phi0)/3.d0)
 fix1%x(2)=jx(1)*exp(i_*(phi0)/3.d0)
 fix2%x(1)=jx(2)*exp(-i_*(pi+phi0)/3.d0)
 fix2%x(2)=jx(2)*exp(i_*(pi+phi0)/3.d0)


phi0=atan2(aimag(gamma),real(gamma))
id=0

f1=0
f1= fix1 
f2=0
f2= fix2

 !!!!!!!!!!!!!!   Computation of fixed points to all orders using Newtow search !!!!!!!!!!!!!!

write(6,"(8(1x,g15.8,1x))") f1%x(1:4)  
write(6,"(8(1x,g15.8,1x))") f2%x(1:4) 
 call Newton_search(fq,f1) 
call Newton_search(fq,f2) 
write(6,"(8(1x,g15.8,1x))") f1%x(1:4) 
write(6,"(8(1x,g15.8,1x))") f2%x(1:4) 
 

 u%n=c_%nv
uc%n=c_%nv
fr1%n=c_%nv
fr2%n=c_%nv
call alloc(u,uc,fr1,fr2)
 

u=1

u=fq
 
uc=1
do i=1,c_%nv
uc%v(i)=f1%x(i)+uc%v(i)
enddo

 

fr1=u.o.uc
 
!do i=3,c_%nd2
do i=3,c_%nv
fr1%v(i)=1.e0_dp.cmono.i
enddo

 !call print(fr1)
 

 
fr1=fr1.oo.(-1)

uc=1
do i=1,2
uc%v(i)=0
enddo

fr1=fr1.o.uc
 
do i=1,c_%nv
fr1%v(i)=f1%x(i)+fr1%v(i)
enddo

 
u=1

u=fq


uc=1
do i=1,c_%nv
uc%v(i)=f2%x(i)+uc%v(i)
enddo

 

fr2=u.o.uc
 

!do i=3,c_%nd2
do i=3,c_%nv
fr2%v(i)=1.e0_dp.cmono.i
enddo

  

 
fr2=fr2.oo.(-1)

uc=1
do i=1,2
uc%v(i)=0
enddo

fr2=fr2.o.uc

!do i=1,c_%nd2
do i=1,c_%nv
fr2%v(i)=f2%x(i)+fr2%v(i)
enddo

 
call kill(u,uc)
call alloc(u,uc)
 
jy=0
if(c_%nd2==4) then
 write(6,*) "sqrt(J_y) (<=0 nothing done)  -> 1.5e-4 or smaller (resonances above) "
 read(5,*)  jy
endif
 

fix1=f1
fix2=f2

fix1=c_phasor().o.fix1
fix2=c_phasor().o.fix2


fix1k=(c_n%atot.sub.1).o.fix1
fix2k=(c_n%atot.sub.1).o.fix2


fix1=c_n%atot.o.fix1
fix2=c_n%atot.o.fix2



x1=0
x2=0
x1k=0
x2k=0

x1=fix1%x(1:6) + closed_orbit
x2=fix2%x(1:6) + closed_orbit
x1k=fix1k%x(1:6) + closed_orbit
x2k=fix2k%x(1:6) + closed_orbit

   write(mf,*) " Fixed points of resonance "
   write(mf,format4) x1(1:4)
   write(mf,format4) x2(1:4)
   write(mf,*) "  "
!!!!!!!!!!!!!!!!!!!!  Fully analytical calculation in 1-d-f !!!!!!!!!!!!!!!!!!!!!!!!!
write(6,*) "Fully analytical calculation  0 or 1"
read(5,*) i
if(i==0) goto 998
!!!!   Fully analytical calculation  !!!

ptc_fibre=>ptc_lattice%start
npos=0
do i=1,ptc_lattice%n
if(ptc_fibre%mag%p%nmul>=3) then
 npos=npos+1
endif
ptc_fibre=>ptc_fibre%next
enddo
 
allocate(pha(npos), bet(npos), bet2(npos),spos(npos),f3_1(npos),f3_1t(npos),f3_3(npos))
 
!!! Courant-Snyder phase advance and beta function calculation
ptc_fibre=>ptc_lattice%start
ray_closed=closed_orbit
 call c_fast_canonise(c_n%atot,uc)
ray_8=ray_closed+uc
ptc_fibre=>ptc_lattice%start
 phase=0
 npos=0
 spos=0
s=0
call GET_LENGTH(ptc_lattice,circ)
write(6,*) " circ =", circ
 
del=-(h_res.sub.[1,1,0,0])
om=c_n%tune(1)*twopi
 

do i=1,ptc_lattice%n
 call propagate(ray_8,state,fibre1=ptc_fibre,fibre2=ptc_fibre%next)

s=s+ptc_fibre%mag%p%ld
u=ray_8
xs0=ray_8

if(ptc_fibre%mag%p%nmul>=3) then
 npos=1+npos
 pha(npos)=phase(1)*twopi
 bet(npos)=((uc%v(1).sub.'10')**2+(uc%v(1).sub.'01')**2)

ld=1
if(ptc_fibre%mag%l/=0.d0) then
 ld=ptc_fibre%mag%l
endif
 bet2(npos)=(bet(npos))**2d0*ptc_fibre%mag%bn(4) *ld
 

 bet(npos)=(bet(npos))**1.5d0*ptc_fibre%mag%bn(3) *ld
 f3_1(npos)=-bet(npos)*exp(-i_*pha(npos))/sqrt(8.d0) 
 f3_1t(npos)=-0.5d0*f3_1(npos)*(1.d0+exp(-I_*om))/(1.d0-exp(-I_*om))
 f3_3(npos)=-bet(npos)*exp(-i_*3*pha(npos))/sqrt(8.d0)/3.d0
 spos(npos)=s
endif

 call c_fast_canonise(u,uc,phase=phase)
ray_8=xs0+uc  


 ptc_fibre=>ptc_fibre%next
enddo
 
 
mu_tot=phase(1)*twopi
mu_tot=c_n%tune(1)*twopi

 
f31_2=0
 do i=1,npos
do  j=i+1,npos
 f31_2=6.d0*i_*(CONJG(f3_1(i))*f3_3(j) - CONJG(f3_1(j))*f3_3(i) )/2.d0 +f31_2
enddo
enddo

f31_2=f31_2/(exp(-I_*2.d0*om)-1.d0)

f40_2=0
 do i=1,npos
do j=i+1,npos
f40_2=3.d0*i_*(f3_1(i)*f3_3(j)-f3_1(j)*f3_3(i))/2.d0 +f40_2
enddo
enddo
f40_2=f40_2/(exp(-I_*4.d0*om)-1.d0)

 
f21=0
do i=1,npos
f21=f21-f3_1(i)
enddo
g3=f21
f21=-(1.d0/(exp(-I_*om)-1.d0)) *f21

g3h=0
do i=1,npos
g3h=g3h+f3_3(i)
enddo

!T,id,A,N_c,rot,u,uc
write(6,*) g3
write(6,*) g3h

A=c_n%atot.cut.2

ID=A**(-1)*c_n%atot
call c_factor_map(c_n%atot,U,ft,dir=-1)
call clean(ft,ft,prec=1.d-7)
 
call clean(id,id,prec=1.d-7)
 

ft3=c_logf_spin(id)
ft3=ft3.cut.3
call clean(ft3,ft3,prec=1.d-7)
hh=getpb(ft3)
hh=hh*c_phasor()
call clean(hh,hh,prec=1.d-7)
 
ID=exp(-ft3,ID)
ft4=c_logf_spin(id)
ft4=ft4.cut.4
hh=getpb(ft4)
hh=hh*c_phasor()
 
  
 A=A*exp(ft3)
 id=A**(-1)*T*A
call c_factor_map(id,id,ft,dir=1)
ft=ft.cut.4
hh=getpb(ft)
hh=hh*c_phasor()
call clean(hh,hh,prec=1.d-7)

!call print(hh)


f40_1=0
 do i=1,npos
do j=1,npos
 f40_1=3.d0*i_*(f3_1t(i)*f3_3(j)) +f40_1
enddo
enddo
 
f40_1=f40_1/(exp(-I_*4.d0*om)-1.d0)

 





f31_1=0
 do i=1,npos
do  j=1,npos
 f31_1=6.d0*i_*(CONJG(f3_1t(i))*f3_3(j)    ) +f31_1   !+ CONJG(f3_1t(j)) * f3_3(i)
enddo
enddo
f31_1=f31_1/(exp(-I_*2.d0*om)-1.d0)


 
 




 
 
id=c_n%atot
id=id.cut.2
id=id**(-1)
id=id*c_n%atot
id=ci_phasor()*id*c_phasor()
ft=c_logf_spin(id)
 call d_field_for_demin(ft,F_res)

write(mf,*) " This is F_res"
call print(F_res,mf)

write(mf,*) " This is -2pi H_res"
call print(H_res,mf)


 


c3nu=0
do i=1,npos
c3nu=c3nu+bet(i)*exp(-I_*3*pha(i))
enddo
c3nu=(i_*del/(exp(-I_*3*del)-1.d0 ))/(2.d0*sqrt(2.d0))*c3nu

 
alpha1=0
alphaoct=0
f40_oct=0.d0
f31_oct=0.d0
do i=1,npos
alphaoct= -bet2(i)*3.d0/8.d0  +alphaoct
f40_oct=f40_oct-exp(-I_*4.d0*pha(i))*bet2(i)/16.d0
f31_oct=f31_oct-exp(-I_*2.d0*pha(i))*bet2(i)/4.d0
do j=1,npos
alpha1=bet(i)*bet(j)*cos(3.d0*(pha(j)-pha(i))) +alpha1
enddo
enddo
alpha1=alpha1*(1.d0/16.d0)*(1.d0/tan(3*del/2.d0)-3.d0*del/(1.d0-cos(3*del))) 
f40_oct=f40_oct/(exp(-I_*4.d0*om)-1.d0)
f31_oct=f31_oct/(exp(-I_*2.d0*om)-1.d0)

 f40=f40_1+f40_2+f40_oct
f31=f31_2+f31_1+f31_oct
alpha2=0
do i=1,npos
do j=1,npos
alpha2=bet(i)*bet(j)*cos((pha(j)-pha(i))) +alpha2
enddo
enddo
alpha2=alpha2*(3.d0/tan(om/2.d0)/16.d0)

alpha3=0
do i=1,npos
do j=i+1,npos
alpha3=bet(i)*bet(j)*(3*sin(pha(j)-pha(i)) +sin(3*(pha(j)-pha(i)))     ) +alpha3
enddo
enddo
alpha3= (18d0*4.d0/24.d0**2)*alpha3


close(mf)

mf=6

write(mf,* ) " FPP/Analytical  results result A_s ";write(6,*)

write(mf,*) " f21 ",f21,abs(f21)
write(mf,*) " ptc ",(F_res.sub.[2,1,0,0]),abs(F_res.sub.[2,1,0,0])
write(mf,*) " f40 ",f40,abs(f40)
write(mf,*) " ptc ", (F_res.sub.[4,0,0,0]),abs(F_res.sub.[4,0,0,0])
write(mf,*) " f31 ",f31,abs(f31)
write(mf,*) " ptc ",(F_res.sub.[3,1,0,0]),abs(F_res.sub.[3,1,0,0])

write(mf,* ) " FPP/Analytical  results result H_r ";write(6,*)
write(mf,*) " Analytical results result for Gamma "
write(mf,*) c3nu,abs(c3nu)
write(mf,*) " FPP Gamma"
write(mf,*) (h_res.sub.[3,0,0,0]),abs((h_res.sub.[3,0,0,0]))

 
write(mf,*) " Analytical  results for H_res -2pi*alpha_xx/2"
alpha1= (real(alpha1)+real(alpha2)+real(alpha3)+real(alphaoct))
write(mf,*) "Tune shift term = ", real(alpha1) 
 
write(mf,*) " FPP result"
write(mf,*) "Tune shift term = ", real((H_res.sub.[2,2,0,0]))








998  continue
goto 500

!!! Stuff for Forest's Windows Graphical interface
500 continue
call kanalnummer(mfoo,"f_ask.txt")
    write(mfoo,*) 2,c_%nd2,C_%NO
   write(mfoo,format6) x1k
 
    write(mfoo,format6) x2k
close(mfoo)
call kanalnummer(mfoo,"f_as.txt")
    write(mfoo,*) 2,c_%nd2,C_%NO
   write(mfoo,format6) x1
 
    write(mfoo,format6) x2
close(mfoo)

call kanalnummer(mfoo,"f1.txt")
    write(mfoo,*) 1,c_%nd2,C_%NO
   write(mfoo,format6) x1

close(mfoo)

call kanalnummer(mfoo,"f2.txt")
    write(mfoo,*) 1,c_%nd2,C_%NO
 
    write(mfoo,format6) x2
close(mfoo)

 
  call alloc(a_tot_inverse,1,c_%NV,c_%nd2)
  call alloc(a_tot,1,c_%NV,c_%nd2)

id=c_n%atot**(-1)
do i=1,c_%nd2
 a_tot(i)=c_n%atot%v(i)
 a_tot_inverse(i)=id%v(i)
enddo
 


999 my_euni_1=>a_tot_inverse
 my_euni_2=>a_tot
 my_estate=>state

 


1001 call kill(T,id,A,N_c,rot,u,uc)
call kill(xs)
call kill(ray_8)
call kill(c_n)
 

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
ksd=6.447435260914397D-03
!if(.not.skew)ksd=6.447435260914397D-02-2.5d-2
!if(skew) 
ksd=6.447435260914397D-02+3.1d-2
QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.06810987731199d0)
QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.51445064774657d0); 
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


subroutine Newton_search(h_co,f)
implicit none
type(c_vector_field)  f_rot,h_co
type(c_ray) f,f1h
real(dp) epsn0,epsn,epsnb
integer nmax0,ncut,i,k
logical potential_exit
type(c_damap) id,rot
complex(dp) dx

 ncut=c_%no+1 
 epsn0=abs(f%x(1))+abs(f%x(2))+abs(f%x(3))+abs(f%x(4))
 epsn0=epsn0*1.e-6_dp

call alloc(f_rot)
call alloc(id,rot)
 
nmax0=100

  f1h=f     !  (6a) 
  potential_exit=.false.
  epsnb=1.e38_dp
  do k=1,nmax0

     id=1   ! (6b)
     do i=1,c_%nd2
        id%v(i)=id%v(i)+f1h%x(i)   ! (6c)
     enddo
 
!  evaluating the vector field around the fixed point    
     do i=1,c_%nd2
        f_rot%v(i)=h_co%v(i).o.id   ! (6d)
     enddo
 
! putting in a c_damap because FPP can only inverse maps
   id=1
    do i=1,c_%nd2
      id%v(i)=f_rot%v(i)   
     enddo
 
 
!    F(x)=0, the x is the fixed point
!  .oo. does a TPSA inversion of the map, 
!  ie, taking into account the constant part
! This inversion is exact if and only 
! if the calculation is done to infinite order
! Therefore this is part of a Newton search

!rot=1
!do i=3,c_%nd2
! rot%v(i)=0.0_dp
!enddo
! removing all "vertical dependence "  (6f)
!id=id.o.rot   
! adding identity in the y-py planes to allow inversion
!do i=3,c_%nd2
! id%v(i)=1.0_dp.cmono.i ! (6g)  
!enddo
 

!doing a linear Newton search if ncut = 2
     id=id.cut.ncut   !  (6h) 
 do i=1,c_%nd2
   rot%v(i)=-(id%v(i).sub.'0')
   id%v(i)=id%v(i)+rot%v(i)
 enddo

     id=id.oo.(-1)    !  (6i)    Inversion
     id=id.o.rot
 

!   The fixed point is updated by taking the constant part of the map id.
     epsn=0
     do i=1, 2  !c_%nd2
     dx=(id%v(i).sub.'0')  !
    epsn=abs(dx)+epsn
 
        f1h%x(i)=f1h%x(i)+dx   ! (6j) updating the fixed point
     enddo
     f1h%x(3:4)=f%x(3:4)
 !if(.false.) then
       if(potential_exit) then  
        if(abs(dx)==0.0_dp.or.epsn>=epsnb) exit
       endif
      if(epsn<epsn0) then
       potential_exit=.true.
     endif
!endif
         epsnb=epsn
 
!     write(6,format8) (f1h%x(1:4))  ! check convergence

  enddo

if(k>nmax0) then
 write(6,*) " Search was not succesful in", nmax0 , " steps"
else 
 write(6,*) " Search was succesful in", k, " steps"
endif


f=f1h
call kill(f_rot)
call kill(id,rot)
 

end subroutine 

 


end program
