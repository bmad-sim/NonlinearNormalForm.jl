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
  type(c_damap):: T,id,A,N_c,rot,u,uc
  type(c_normal_form) :: normal_form
  type(probe_8) :: ray_8,xs
  type(probe)  :: ray_closed,ray,xs0,xst 
  real(dp) closed_orbit(6),r,phi0,phase(3),om,g,ang(3)
  real(dp) jx(2),discriminant
  real(dp) x1(6),x2(6),x1k(6),x2k(6),sigm(6),cut
  complex(dp) gamma 
  integer i,k,n_arg,mf,mfoo,nturns,nmax0 ,j,mfr,no,p_res
  character(100) lat_file
  type(c_universal_taylor)   H_res,F_res,H_bengtsson
  type(c_ray) f1,f2
  complex(dp) dx,del,c3nu,alpha1,alpha2,alpha3,alpha_xx_half,alphaoct,f21,ld,mu(3)
  complex(dp) f40_2,f40,f40_1,g3h,g3,f31_1,f31_2,f40_oct,f31_oct,f31,sylee,bengtsson
  logical use_vector_field,res
  type(c_universal_taylor), target ::   a_tot_inverse(6),a_tot(6)
  integer jinv(2),npos,mis
  real(dp),allocatable ::  bet(:),bet2(:),pha(:),spos(:)
  complex(dp) ,allocatable ::f3_1(:),f3_1t(:),f3_3(:)
  type(c_ray) fix1,fix2,fix1k,fix2k
  type(c_vector_field)ft ,ft3,ft4
  complex(dp) ,allocatable :: A_l(:)
  integer nlmin,nlmax,jl,je(6)
  real(dp) s,circ,mu_tot, alpha0,mr2 
  type(c_taylor) hh,phat(3)
  !!!  input parameters
  preffile= "C:\document\etienne_programs_new\programs_for_learning\park\pref3nux.txt"
  
  
  nlmin=-100
  nlmax=100
  allocate(A_l(nlmin:nlmax))
  
  A_l=0
   
  call kanalnummer(mf,"C:\document\etienne_programs_new\programs_for_learning\suntao\results_rad.txt")
  
  call ptc_ini_no_append ! initializes PTC
  
  !call read_lattice_append(M_U,'C:\document\etienne_programs_new\programs_for_learning\park\park_cell_01.txt')
  !call read_ptc_command('C:\document\etienne_programs_new\programs_for_learning\park\fittune_resonance_3nu.txt')
  call read_lattice_append(M_U,"C:\document\my_tex_papers\fpp_handbook\talk\cesr_3nux_ptc.lat")
  !call read_lattice_append(M_U,"C:\document\my_tex_papers\fpp_handbook\talk\cesr_3nux_ptc_exactmodeloff.lat")
  
  
  !state=radiation0 !+time0
  !state=time0
  state=only_4d0 
  my_estate=>state
    
  ! goto 1000
  
  !!! next, use PTC normal form to compute intrinsic resonance strength
  ptc_lattice=>m_u%start
  ptc_fibre=>ptc_lattice%start
  
  !write(6,*) "Misalignment "
   
  !call read_ptc_command('C:\document\etienne_programs_new\programs_for_learning\resonance orbital\fittune_resonance_3nu.txt')
  
  
  
  no=3
  write(6,*) "order >= 3 "
  read(5,*) no
  if(no<3) no=3
  !! Counting numbers of multipoles 
  
  ptc_fibre=>ptc_lattice%start
  npos=0
  
  
  
  !call in_bmad_units
  use_quaternion=.true.
   n_cai=-i_
  ptc_oct=>ptc_fibre
   
    
  my_estate=>state
  
   ptc_fibre=>ptc_lattice%start
   
   
  call init(state, no, 0)
  
  call alloc(T,id,A,N_c,rot,u,uc)
  call alloc(xs)
  call alloc(ray_8)
  call alloc(normal_form)
  call alloc(rot)
  call alloc(ft,ft3,ft4)
  call alloc(hh)
  call alloc(phat)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  closed_orbit=0
  
  ! calculation of closed orbit
  call FIND_ORBIT_x(closed_orbit,state, 1.0e-7_dp, fibre1=ptc_fibre)
  
  write(6,*) "closed orbit "
   write(6,format4)closed_orbit(1:6)
  
  ray_closed=closed_orbit 
  
   
  
  id=1
  ray_8=id+ray_closed
   call propagate(ray_8,state,fibre1=ptc_fibre)
  T=ray_8
  
  call kanalnummer(mfr,"map_res.txt")
   call print(t,mfr)
  close(mfr)
  
  call c_normal(T,normal_form,phase=phat)
  alpha0=phat(1).sub.'11'
  
  !call clean(phat,phat,prec=1.d-8)
  !call print(phat(1))
  
  
   
   
  res=.true.
   
    normal_form%positive=.false.
  if(res) normal_form%nres=(c_%no+1)/3
    write(6,*) " # of resonance terms ",normal_form%nres
    do i=1, (c_%no+1)/3
    if(res) normal_form%m(1,i)=3*i ! 1/3 order 
    enddo
   
  mr2=0
  if(res) then
   do i=1,c_%nd
  mr2=mr2+normal_form%m(i,1)**2
  enddo
   
  s=0 
  do i=1,c_%nd
  mu(i)=phat(i)
  s=s+normal_form%m(i,1)*mu(i)
  mu(i)=mu(i)*twopi
  enddo
  p_res=nint(s)
   endif
  
   
  call c_normal(T,normal_form,canonize=.true.)
  
  write(6,*) normal_form%tune(1:3)
  write(6,*) normal_form%spin_tune
   !!!!!!!!!!!!!!!
    call alloc(a_tot_inverse,1,c_%NV,c_%nd2)
    call alloc(a_tot,1,c_%NV,c_%nd2)
  
  id=1
  do i=1,c_%nd2
  id%v(i)=id%v(i)+closed_orbit(i)
  enddo
  
  a=id.o.normal_form%atot
  
  do i=1,c_%nd2
   a_tot(i)=a%v(i) 
  enddo
  
  a=normal_form%atot**(-1)
  
  
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
  
   
  if(.not.res)  goto 1000
   
  !!!!!!!!!!!!!!!!!!!!!!
   call clean(normal_form%h_nl,normal_form%h_nl,prec=1.e-10_dp)
  
  use_vector_field=.false.
  if(use_vector_field) then
   normal_form%h_l%v(1)=(-i_*twopi/3.0_dp)*dz_c(1)+normal_form%h_l%v(1)
   normal_form%h_l%v(2)=( i_*twopi/3.0_dp)*dz_c(2)+normal_form%h_l%v(2)
   normal_form%h_l%v(3)=0.0_dp
   normal_form%h_l%v(4)=0.0_dp
   N_c=exp(normal_form%h_l)
   N_c=exp(normal_form%h_nl,N_c)
  else
  write(6,*) " mr2  ",mr2
   
  
  del=0
  do i=1,c_%nd
  del=del+normal_form%m(i,1)*mu(i)/mr2
  enddo
  !!! removeal of the perpendicular terms
    do i=1,c_%nd
    ft%v(2*i-1)=-i_* del*normal_form%m(i,1)  *dz_c(2*i-1) 
    ft%v(2*i)= i_* del*normal_form%m(i,1)  *dz_c(2*i) 
    ft%v(2*i-1)=ft%v(2*i-1)- normal_form%h_l%v(2*i-1)
    ft%v(2*i)= ft%v(2*i)-normal_form%h_l%v(2*i)
  enddo
   N_c=normal_form%atot**(-1)*T*normal_form%atot
  n_c=ci_phasor()*n_c*c_phasor()
  call clean(n_c,n_c,prec=1.d-10)
  call print(n_c)
   N_c=exp(ft,N_c) 
  call clean(n_c,n_c,prec=1.d-10)
  pause 77
   call print(n_c)
  write(6,*) p_res
  pause 888
  !  removal of co-moving
    do i=1,c_%nd
    ft3%v(2*i-1)=i_* p_res*normal_form%m(i,1)*twopi/mr2  *dz_c(2*i-1) 
    ft3%v(2*i)= -i_* p_res*normal_form%m(i,1)*twopi/mr2   *dz_c(2*i) 
  enddo
   N_c=exp(ft3,N_c)
  call clean(n_c,n_c,prec=1.d-10)
  
   call print(n_c)
   
  endif
   
    
   
   normal_form%h=c_logf_spin(N_c)
    call clean(normal_form%h,normal_form%h,prec=1.e-7_dp)
  write(6,*) " Vector field "
  write(6,*) normal_form%tune(1)+1.d0/3.d0
  write(6,*) normal_form%damping(1) 
  normal_form%h_l=(-1.d0/twopi)*normal_form%h
  call print(normal_form%h%v(1))
  call print(normal_form%h%v(2))
  call print(normal_form%h_l%v(1))
  call print(normal_form%h_l%v(2))
   write(6,*) " End of Vector field "
   
   call d_field_for_demin(normal_form%h,H_res)
   
    call clean(H_res,H_res,prec=1.e-5_dp)
  
  call print(H_res)
   
   write(mf,*); write(mf,*) " -2*pi*H_r "; write(mf,*)
   call print(H_res,mf)
   
  je=0
  je(1)=3
  !gamma= c_get_coeff(h_res,[3,0,0,0])  
  gamma= h_res.sub.je(1:c_%nd2)
  !alpha_xx_half= c_get_coeff(h_res,[2,2,0,0])
  !del=c_get_coeff(h_res,[1,1,0,0])
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
   call print(H_res,6)
  
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
  
  
  call clean(normal_form%h,normal_form%h,prec=1.d-4)
  call print(normal_form%h)
  
  
  
  write(6,*) " mess-up f2 ?"
  !read(5,*) i
  i=0
  if(i==1) then
  write(6,"(8(1x,g15.8,1x))") f2%x(1:4)  
  
   read(5,*) x1(3:4)
   fix2%x(3)=x1(3)
   fix2%x(4)=x1(4)
  x1=0
  endif
  write(6,"(8(1x,g15.8,1x))") f1%x(1:4)  
  write(6,"(8(1x,g15.8,1x))") f2%x(1:4) 
  !goto 1000
  !stop
  call Newton_search(normal_form%h,f1) 
  write(6,"(8(1x,g15.8,1x))") f1%x(1:4) 
  
  call Newton_search(normal_form%h,f2) 
  write(6,"(8(1x,g15.8,1x))") f2%x(1:4) 
  
  fix1=f1
  fix2=f2
  
  fix1=c_phasor().o.fix1
  fix2=c_phasor().o.fix2
  
  
  fix1k=(normal_form%atot.sub.1).o.fix1
  fix2k=(normal_form%atot.sub.1).o.fix2
  
  
  fix1=normal_form%atot.o.fix1
  fix2=normal_form%atot.o.fix2
  
  
  
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
  
  write(6,*) "Fully analytical calculation"
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
   call c_fast_canonise(normal_form%atot,uc)
  ray_8=ray_closed+uc
  ptc_fibre=>ptc_lattice%start
   phase=0
   npos=0
   spos=0
  s=0
  call GET_LENGTH(ptc_lattice,circ)
  write(6,*) " circ =", circ
   
   
  
    
  del=-(h_res.sub.[1,1,0,0])
  om=normal_form%tune(1)*twopi
   
  
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
  
  ld=-1.d0/3.d0+(-3.264239087234111E-003)
  write(6,*) phase(1)
  phase(1)=nint(phase(1))+ld
  write(6,*) phase(1)
  write(6,*) om/twopi,-17.d0+phase(1)
  del=0.2050981907203069E-01
  !write(6,*) del
  pause 12
  !stop
  mu_tot=phase(1)*twopi
  mu_tot=normal_form%tune(1)*twopi
    
  
   
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
  
  A=normal_form%atot.cut.2
  
  ID=A**(-1)*normal_form%atot
  call c_factor_map(normal_form%atot,U,ft,dir=-1)
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
  
  
   
   
  
  
  
  
   
   
  id=normal_form%atot
  id=id.cut.2
  id=id**(-1)
  id=id*normal_form%atot
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
   do jl=nlmin,nlmax
   A_l(jl)=A_l(jl) - bet(i)*exp(-I_*3*(pha(i)-spos(i)*mu_tot/Circ) ) &
  *exp(-I_*jl*spos(i)*2.d0*pi/Circ)/3.d0/sqrt(8.d0)/Circ
   enddo
  enddo
  c3nu=(i_*del/(exp(-I_*3*del)-1.d0 ))/(2.d0*sqrt(2.d0))*c3nu
  
  
   
   
   
  
  
  
  
  
  !stop
  
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
  
  
  
  if(.true.) then
   a=normal_form%atot.cut.2
    a=a**(-1)*T*a
   a=ci_phasor()*a*c_phasor()
   id=a.cut.2
  
  a=a*id**(-1)
   
   ft=c_logf_spin(a)
   call clean(ft,ft,prec=1.d-6)
  write(mf,'(/,(a))') " L exp(:hh:)  (Bengtsson factorization)   "
  write(mf,'(/,(a))') " vector field   "
  
   call print(ft,mf)
   hh=cgetpb(ft)
   call d_field_for_demin(ft,H_bengtsson)
   !call print(hh,mf)
  write(mf,'(/,(a))') " Poisson Bracket version   "
   call print(H_res,mf)
  bengtsson=H_bengtsson.sub.[2,2,0,0]
  Write(mf,'(/,(a))') " Bengtsson H_2200 "
  Write(mf,*) bengtsson
  endif
  
  
  
  
  
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
  jl=43
  dx=(circ/twopi)*(exp(i_*(jl*twopi-3.d0*mu_tot))-1.d0)/i_/(jl-3.d0*mu_tot/twopi)
  sylee=dx*a_l(jl)
  write(mf,*) " SY Lee Gamma l=nint(3*total_tune) "
  write(mf,"(i4,3(1x,g23.16,1x))") jl,sylee, abs(sylee)
   
  write(mf,*) " Analytical  H_res -2pi*alpha_xx/2"
  alpha1= (real(alpha1)+real(alpha2)+real(alpha3)+real(alphaoct))
  write(mf,*) "Tune shift term = ", real(alpha1) 
   
  write(mf,*) " FPP result"
  write(mf,*) "Tune shift term = ", (H_res.sub.[2,2,0,0])
  write(mf,*) "Bengtsson h2200 = ", bengtsson 
  
  
  
  
  
  write(mf,*)" Fourier terms in 's'  ";write(mf,*)
  
  
  c3nu=0
  do jl=nlmin,nlmax
  
  dx=(circ/twopi)*(exp(i_*(jl*twopi-3.d0*mu_tot))-1.d0)/i_/(jl-3.d0*mu_tot/twopi)
  write(mf,"(i4,3(1x,g23.16,1x))") jl,dx*a_l(jl), abs(dx*a_l(jl))
  c3nu=c3nu+dx*a_l(jl)
  enddo
  
  write(mf,*) c3nu,abs(c3nu)
  
  close(mf)
  
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
  
  id=normal_form%atot**(-1)
  do i=1,c_%nd2
   a_tot(i)=normal_form%atot%v(i)
   a_tot_inverse(i)=id%v(i)
  enddo
   
  
  
  999 my_euni_1=>a_tot_inverse
   my_euni_2=>a_tot
   my_estate=>state
  
   
  
  
  1001 call kill(T,id,A,N_c,rot,u,uc)
  call kill(xs)
  call kill(ray_8)
  call kill(normal_form)
   
  
  1000 call ptc_end(graphics_maybe=1,flat_file=.false.)
  
  
   contains
  
  
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
   
  !call print(id)
  !pause 1231
   
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
  