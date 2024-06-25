program example
  use c_tpsa
  use gauss_dis
  implicit none 
  integer no,nd,i,nv,k,j(lnv),n,np,j1,ndpt,i1,i2,is1(3),is2(3),is3(3),is2n(3),is3n(3)  
  type(c_taylor) f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin
  type(c_damap), target :: m,m1,rot,as,a0,a1,a2,a_cs,R_TE,CS_TE
  type(c_damap), pointer :: o
  type(c_vector_field) vf1,vf2,vf3  !,h_rot,K_r,K_r_real
  real(dp) r1,r2,sca,decrement(6),radfluc,radx,emi(3),eigen_emi(3) !,emi_chao(3)
  complex(dp) c1
  type(c_normal_form) normal
  type(c_vector_field) vf
  type(c_ray) ray,ray2
  type(c_lattice_function)   Lat_function
  logical :: damp = .false.,fastcanonize=.true.,COSLIKE,t_e
  type(c_taylor), allocatable :: phase(:)
  type(c_damap), allocatable :: mt(:)
  real(dp)  ph(3),spintune(2),dampdec(3)
  type(c_lattice_function) cl
   logical putspin
  
  radfluc=0.0001d0
  remove_tune_shift=.false.
  putspin=.true.
  n_cai=-i_
  write(6,*) " If deltap/p0 is a canonical variable  enter 3"
  write(6,*) " else enter 2"
  !read(5,*) nd
  nd=3  ! coasting
  !if(nd/=2.and.nd/=3) stop 44
  if(nd==1) ndpt = 0
  if(nd==2) ndpt = 0
  if(nd==3) ndpt = 6  ! BMAD choice
  !if(nd==3) ndpt = 0  ! BMAD choice
  no=3;     ! no: the order of the polynomial    nv: the number of variables   
  np=0
  c_lda_used=1500
   call gaussian_seed(2463)
  use_quaternion=.true.
  call c_init(no,nd,np1=np,ndpt1=ndpt)  ! initializes taylor series with maps
  allocate(phase(nd))
   
  call alloc(f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin)      ! must be constructed after init
  call alloc(vf1,vf2,vf3)      
  call alloc(m,m1,as,a0,a1,a2,a_cs,rot,R_TE,CS_TE)      ! 
  call alloc(phase)
  call alloc(nu_spin)
  call alloc(vf)
  
  n=1
  allocate(mt(n))
  do i=1,n
   call alloc(mt(i))
  enddo
   !default_fractional_tune_positive=.false.
    negative_synchrotron_tune=.false.
  
   
   call alloc(normal)
  
  r1=1.0_dp
  r2=0.5d0
  
  
  do k=1,n
  vf1=0
  
  call c_TAYLOR_ran(f,r1,r2)
  f=f-i_*aimag(f)  ! remove the imaginary part
  f=f-(f.cut.2)
  f=-f*0.002d0
  do i=1,nd
  sca=.1342d0
  !do while(sca==0.d0) 
  !  call GRNF(sca,5.d0)
  !enddo
  !sca=sign(sca,1.d0)+sca  
  sca=abs(sca)/i
   f=f+sca*(dz_c(2*i-1)**2+dz_c(2*i)**2)
  enddo
   
  ! Removing time terms for coasting
  if(c_%ndpt==6) then
   m1=1
   m1%v(2*nd-1)=0.0_dp
   f=f*m1
  endif
  if(c_%ndpt==5) then
   m1=1
   m1%v(2*nd)=0.0_dp
   f=f*m1
  endif
  
   
  vf1=getvectorfield(f)  !  
   
  ! putting spin
  if(putspin) then
  do i=1,3 
   call c_TAYLOR_ran(vf1%q%x(i),r1,r2) ! putting spin
   vf1%q%x(i)=vf1%q%x(i)-i_*aimag(vf1%q%x(i))
   vf1%q%x(i)=0.3d0*vf1%q%x(i)*m1 !   
  
  enddo
   endif
  
  mt(k)=exp(vf1)
  
  goto 999
  do i1=1,c_%nd2
  do i2=i1,c_%nd2
   call GRNF(radx,4.d0)
   radx=radx*radfluc
  if(i1==i2) then
   mt(k)%e_ij(i1,i2)=abs(radx)
  else
   mt(k)%e_ij(i1,i2)=abs(radx)/10
   mt(k)%e_ij(i2,i1)=abs(radx)/10
  
  endif
  enddo
  enddo
  999 continue
  
  
  enddo
  
  m=1
  do i=1,n
   m=mt(i)*m
  enddo
  
  
  
  
  decrement(1)=0.99999993d0
  decrement(2)=0.9999997d0
  decrement(3)=0.9999996d0
  decrement(4)=0.99999995d0
  decrement(5)=0.99999993d0
  decrement(6)=0.99999999d0
  
  decrement(1)=0.89d0
  decrement(2)=0.87d0
  decrement(3)=0.86d0
  decrement(4)=0.895d0
  decrement(5)=0.83d0
  decrement(6)=0.89d0
  
  decrement(1)=0.9993d0
  decrement(2)=0.997d0
  decrement(3)=0.996d0
  decrement(4)=0.9995d0
  decrement(5)=0.9993d0
  decrement(6)=0.9999d0
  
  !remove_tune_shift=.true.
  do i=1,c_%nd2
  m%v(i)=m%v(i) !*decrement(i)
  enddo
  call print(m)
  stop
  
  lielib_print(4)=0
  call c_normal(m,normal,dospin=putspin,phase=phase,nu_spin=nu_spin)
  
  !call c_fast_canonise(normal%atot,a1,dospin=putspin)
  call c_fast_canonise_clean_julia(normal%atot,a1,dospin=putspin)
  
  call print(a1)
  stop
  
  normal%atot=a1
   write(6,*) normal%damping(1:3)
  m=normal%atot**(-1)*m*normal%atot
    m=ci_phasor()*m*c_phasor()
  
  call clean(m,m,prec=1.d-10)
  m=m.cut.2
  call print(m)
  stop
  
  m1= mt(1)*m*mt(1)**(-1)
  !m1=m1.cut.2
  
  normal%atot= mt(1)*normal%atot
  a2=normal%atot
  !a2=a2.cut.(2)
  a_cs=a2.cut.(2)
  call print(a_cs)
  
  call c_fast_canonise(a2,a0,dospin=putspin)
  
  call print(a0)
  
  m1=a2**(-1)*m*a2
  m1=a0**(-1)*m*a0
  m1=ci_phasor()*m1*c_phasor()
  call clean(m1,m1,prec=1.d-10)
  m1=m1.cut.2
  call print(m1)
  
  stop
  !a0=a0.cut.(-2)
  !normal%atot=normal%atot.cut.(-2)
   
  call c_fast_canonise_clean_julia(normal%atot,a1,dospin=putspin)
  
  do i=1,c_%nd2
   a0%v(i)= a0%v(i)+i_*a1%v(i)
  enddo
  do i=0,3
   a0%q%x(i)=  a0%q%x(i)+i_*a1%q%x(i)
  enddo
  call print(a0)
  do i=1,c_%nd2
   a1%v(i)= a1%v(i)+i_*a2%v(i)
  enddo
  do i=0,3
   a1%q%x(i)=  a1%q%x(i)+i_*a2%q%x(i)
  enddo
  !call print(a1)
  
  stop
  !call print(phase(1))
  
  !call clean(phase,phase,prec=1.d-5)
  !call clean(nu_spin,nu_spin,prec=1.d-5)
  
  write(6,'(A,3(1x,g23.16),/)') "The tune    is ",normal%tune(1:nd)
  write(6,'(A,3(1x,g23.16),/)') "The damping is ",normal%damping(1:nd)
  write(6,'(A,1x,g23.16,/)') "The spin tune  is ",normal%spin_tune 
   
   
    m=ci_phasor()*m*c_phasor()
  
   
   call clean(normal%h_l,normal%h_l,prec=1.d-13)
  
   !call print(normal%h_l)
   
   do i1=1,c_%nd2,-1
   do i2=i1,c_%nd2,-1
    write(6,*) normal%s_ij0(i1,i2)
   enddo
  enddo 
  
    call compute_lattice_functions(normal%atot,cl)
  
  call eigen_emit(normal%s_ij0,nd,eigen_emi)
  
  call emit(normal%s_ij0,nd,emi)
  write(6,*) " Chao emittances  "
  write(6,*)normal%emittance 
  write(6,*) "Radar eigen emittances  "
  write(6,*)eigen_emi(1:nd)
  write(6,*) "Radar eigen emittances analytic "
  write(6,*)emi(1:nd)
    write(6,*) "<x^2> correct and via emittances"
   write(6,*) normal%s_ij0(1,1)  
   write(6,*) cl%e(1,1,1)*normal%emittance(1)+cl%e(2,1,1)*normal%emittance(2)+cl%e(3,1,1)*normal%emittance(3)
   
   call sort(normal%emittance,nd,is1)
  
  
  call sort(eigen_emi,nd,is2,is2n)
   
  call sort(emi,nd,is3,is3n)
   
   write(6,*) " Chao emittances  "
  write(6,*)normal%emittance 
  write(6,*) is1
  write(6,*) "Radar eigen emittances  "
  write(6,*)eigen_emi(1:nd)
  write(6,*) is2
  
  write(6,*) "Radar eigen emittances analytic "
  write(6,*)emi(1:nd)
  write(6,*) is3
   
  write(6,*) "emittances from stochastic normal form"
  write(6,*) normal%emittance(1:nd)
  write(6,*) "Radar eigen emittances  "
  write(6,*)eigen_emi(1:nd)
  write(6,*) "Radar eigen emittances analytic "
  write(6,*)emi(1:nd)
  
   
   
   write(6,*) "<x^2> correct and via emittances"
   write(6,*) normal%s_ij0(1,1)  
   write(6,*) cl%e(is1(1),1,1)*normal%emittance(1)+cl%e(is1(2),1,1)*normal%emittance(2)+cl%e(is1(3),1,1)*normal%emittance(3)
   write(6,*) cl%e(is1(1),1,1)*eigen_emi(1)+cl%e(is1(2),1,1)*eigen_emi(2)+cl%e(is1(3),1,1)*eigen_emi(3)
   write(6,*) cl%e(is1(1),1,1)*emi(1)+cl%e(is1(2),1,1)*emi(2)+cl%e(is1(3),1,1)*emi(3)
  
  write(6,*) is1
  write(6,*) is2
  write(6,*) is3
  write(6,*) is3n
  
  contains
  
  subroutine sort(emi,nd,is,isi)
  implicit none 
  real(dp) emi(3), emic(3)
  integer nd,i, k,is(3),isc(3)
  integer,optional:: isi(3)
  
  is=[1,2,3]
  isc=[1,2,3]
  
   emic=emi
  do k=1,nd**2
  do i=1,nd-1
  if(emic(i)>emi(i+1)) then
  isc(i)=is(i+1)
  isc(i+1)=is(i)
  emic(i)=emi(i+1)
  emic(i+1)=emi(i)
  else 
  emic(i)=emi(i)
  emic(i+1)=emi(i+1)
  endif
  emi=emic
  is=isc
  enddo
  enddo
  
  if(present(isi)) then
  do i=1,3
  write(6,*) is(i),i
  isi(is(i))=i
  enddo
  endif
  
  end subroutine sort
  
  subroutine emit(sig,nd,emi)
  implicit none 
  complex(dp) sig(6,6)
  real(dp) emi(3),i2,i4,i6,sig0(6,6),sig2(6,6),dis,xj(6,6),p,q
  integer nd,i,j,k
  
  if(nd==1) then
   emi(1)=sqrt(abs(sig(1,1)*sig(2,2)-sig(1,2)**2))
   return
  endif
  
  if(nd==2) then
  do i=1,6
  do j=1,6
  sig0(i,j)=sig(i,j)
  enddo
  enddo
  
  xj=0
  
  do i=1,nd
  xj(2*i-1,2*i)=1
  xj(2*i,2*i-1)=-1
  enddo
   
  sig0=matmul(sig0,xj)
  
  sig0=matmul(sig0,sig0)
  i2=sig0(1,1)+sig0(2,2)+sig0(3,3)+sig0(4,4)
  
  sig0=matmul(sig0,sig0)
  i4=sig0(1,1)+sig0(2,2)+sig0(3,3)+sig0(4,4)
  i6=0
  !write(6,*) " i2,i4,i6 analytical "
  !write(6,*) i2,i4,i6
  
  dis=sqrt(abs(i4- i2**2/4))
   
  emi(1)=sqrt(abs((i2/2+dis)))/sqrt(2.d0)
  emi(2)=sqrt(abs((i2/2-dis)))/sqrt(2.d0)
   
  endif
  
  if(nd==3) then
  do i=1,6
  do j=1,6
  sig0(i,j)=sig(i,j)
  enddo
  enddo
  
  xj=0
  
  do i=1,nd
  xj(2*i-1,2*i)=1
  xj(2*i,2*i-1)=-1
  enddo
   
  sig0=matmul(sig0,xj)
  
  sig2=matmul(sig0,sig0)
  sig0=sig2
  i2=sig0(1,1)+sig0(2,2)+sig0(3,3)+sig0(4,4)+sig0(5,5)+sig0(6,6)
  
  sig0=matmul(sig0,sig0)
  i4=sig0(1,1)+sig0(2,2)+sig0(3,3)+sig0(4,4)+sig0(5,5)+sig0(6,6)
   
  sig0=matmul(sig0,sig2)
  i6=sig0(1,1)+sig0(2,2)+sig0(3,3)+sig0(4,4)+sig0(5,5)+sig0(6,6)
   
  
  !write(6,*) " i2,i4,i6 analytical "
  !write(6,*) i2,i4,i6
   
   
  p=-i4/4+i2**2/24
   
  q= -i2**3/108 +i2*i4/12 -i6/6;
   
  dis=(3.0_dp*q/2.0_dp/p)
   
  dis=acos(dis*sqrt(3.0_dp/-p))
   
  do k=1,3
  
   emi(k)=2.0_dp*sqrt(-p/3)*cos((dis-k*twopi)/3.0_dp)+I2/6
   emi(k)=sqrt(abs(emi(k)))
  enddo
  
   
  endif
  end subroutine emit
  
  subroutine eigen_emit(sig,nd,emi)
  implicit none 
  complex(dp) sig(:,:)
  real(dp) emi(3),sig0(ndim2t,ndim2t),xj(ndim2t,ndim2t),i2,i4,i6
  integer nd,i,j
  real(dp) reval(ndim2t),imval(ndim2t),vr(ndim2t,ndim2t),vi(ndim2t,ndim2t),vrt(ndim2t,ndim2t),vit(ndim2t,ndim2t)
  complex(dp) emic(6)
  sig0=0
  sig0(1:6,1:6)=sig(1:6,1:6)
  
  
  xj=0
  
  do i=1,nd
  xj(2*i-1,2*i)=1
  xj(2*i,2*i-1)=-1
  enddo
  sig0=matmul(sig0,xj)
   
  
  
    call c_eig6(sig0,reval,imval,vr,vi)
    
  
   
  do i=1,c_%nd2
  emic(i)=reval(i)+i_*imval(i)
  enddo
  do i=1,c_%nd
  emi(i)=abs(imval(2*i))
  enddo
  i2=0
  i4=0
  i6=0
  do i=1,c_%nd2
  i2=i2+emic(i)**2
  i4=i4+emic(i)**4
  i6=i6+emic(i)**6
  enddo
   
  !write(6,*) " i2,i4,i6 numerical "
  
  ! write(6,*) i2,i4,i6
  !pause 667
  !write(6,*)"sum ",2*(imval(2)**2+imval(4)**2)
  
  end subroutine eigen_emit
  
  subroutine invert_with_log(M,MI)
  implicit none
  type(c_damap) m,mi,m1,nl
  type(c_vector_field) f
  call alloc(m1,nl)
  
  call alloc(f)
  
  m1=m.cut.2
  m1=m1**(-1)
  
  nl=m1*m
  
  f=ln(nl)
  f=-f
  
  nl=exp(f)
  
  mi=nl*m1
  
  call kill(f)
  
  call kill(m1,nl)
  
  end subroutine invert_with_log
  
  subroutine invert_with_1minusx(M,MI,minimum_step)
  implicit none
  type(c_damap) m,mi,m1,nl,dl,id
  integer i,nt
  logical minimum_step
  call alloc(m1,nl,dl,id)
  
  nt=c_%no
  if(minimum_step)  nt=(c_%no +mod(c_%no,2))/2
  
  write(6,*) nt,c_%no
  pause 100
  m1=m.cut.2
  m1=m1**(-1)
  
  nl=m1*m
  dl=nl
   call print(nl)
  
  pause 1
  
  mi=1
  id=1
  do i=1,nt
  dl=nl-id
  dl=id-dl
  call clean(nl,nl,prec=1.d-10)
   call print(nl%v(1))
   call print(nl%q%x(1))
  write(6,*) i,nt
  pause 2
  
  nl=dl*nl
  
  mi=dl*mi
  
  enddo
  
  mi=mi*m1
  call kill(m1,nl,dl,id)
  
  end subroutine invert_with_1minusx
  
  
  
  subroutine c_normal1(xyso3,n,dospin,no_used,rot,phase,nu_spin,canonize)
  !#general:  normal
  !# This routine normalises the map xy
  !# xy = n%a_t**(-1)*r*n%a_t 
  !# The linear part of r is described in Chap.4 for the orbital part
  !# and in Chap.6 for the spin. The nonlinear parts are in Chap.5 and 6.
  !# Dospin must be set to .true. if spin is to be normalised.
  !# Resonances can be left in the map. Their number is in n%nres.
  !# They are nres resonances The kth resonance is n%m(i,k).Q_i+n%ms(k)=integer
  !# canonize=.true. Then it is put into courant-snyder form or anti- courant-snyder form
  !# depending on the logical  courant_snyder_teng_edwards=true or false. (See blue or yellow book)
  !#  The map in phasors is exp(n%H_l.grad) exp(n%H_nl.grad)
  !# if fully normalized into a rotation then the map is exp(n%h.grad)
  
      implicit none
      type(c_damap) , intent(inout) :: xyso3
      type(c_damap) m1,ri,nonl,a1,a2,mt,AS,xy,Nf,N_cut_2,N_nl
      type(c_normal_form), intent(inout) ::  n
      type(c_damap), optional :: rot
      type(c_taylor), optional :: phase(:),nu_spin
      type(taylor) c1,s1
      integer,optional :: no_used
      integer i,j,k,l,kr,not,ncoast
      integer, allocatable :: je(:)
      logical(lp) removeit,rad_in
      complex(dp) v,lam,egspin(3)
      complex(dp), allocatable :: eg(:)
      real(dp) norm,alpha,prec !,cx,sx
      logical(lp), optional :: dospin,canonize
      logical dospinr,change
      type(c_spinor) n0,nr
      type(c_quaternion) qnr
      integer mker, mkers,mdiss,mdis,ndptbmad 
         real(dp), allocatable :: da(:)
       integer nd2t,ndc2t,nd2,ndptb,nv,no
  
       nd2t=c_%nd2t
       ndc2t=c_%ndc2t
       nd2=c_%nd2
       ndptb=c_%ndptb
       nv=c_%nv
       no=c_%no
  
      call kill(n%ker)
      call kill(n%g)
      call alloc(n%ker)
      call alloc(n%g)
      n%g%dir=-1
      n%ker%dir=1
  
      if(lielib_print(13)/=0) then
       call kanalnummer(mker,"kernel.txt")
       call kanalnummer(mdis,"distortion.txt")
       call kanalnummer(mkers,"kernel_spin.txt")
       call kanalnummer(mdiss,"distortion_spin.txt")
      endif
  
      dospinr=.false.
      if(present(dospin)) then
       dospinr=dospin
       else
        if(force_spin_input_normal) then
          write(6,*) " your default forces you to include dospin in the input of c_normal"
          stop
        endif
      endif
  
  
  if(bmad_automatic) then
    if(nd2t+ndc2t/=6) then
     write(6,*) " nd2t , ndc2t ",nd2t,ndc2t
     write(6,*) " not BMAD on entrance, suspicious"
    endif
    ndptbmad=0
      alpha=abs(xyso3%v(6).sub.'000001')
      norm=full_abs(xyso3%v(6))
      alpha=abs(alpha-1.0_dp)+abs(norm-1.0_dp)
      if(alpha<1.d-12) then 
        ndptbmad=6
       call in_bmad_units
      endif
       alpha=abs(xyso3%v(5).sub.'000010')
      norm=full_abs(xyso3%v(5))
      alpha=abs(alpha-1.0_dp)+abs(norm-1.0_dp)
      if(alpha<1.d-12) then
         ndptbmad=5
       call in_ptc_units
      endif
    call c_bmad_reinit(ndptbmad)
   
  
    if(use_quaternion) then
      call c_full_norm_quaternion(xyso3%q,k,norm)
      if(k==-1) dospinr=.true.
    else
     call c_full_norm_spin(xyso3%s,k,norm)
      if(k==-1) dospinr=.true.
    endif
  endif
  
  if(spin_automatic) then
    dospinr=.false.
    if(use_quaternion) then
      call c_full_norm_quaternion(xyso3%q,k,norm)
      if(k==-1) dospinr=.true.
    else
     call c_full_norm_spin(xyso3%s,k,norm)
      if(k==-1) dospinr=.true.
    endif
  write(6,*)"dospin ", dospinr
  endif
  
  inside_normal=.true.
  !call c_count_da(i_alloc)
  !write(6,*)" entering c_normal ", i_alloc
      change=.false.
      not=no
      if(present(no_used)) then
        not=no_used  ! sometimes only linear stuff is needed
      else
         if(complex_extra_order==1.and.special_extra_order_1) not=not-1
      endif
  
      call alloc(xy);
      xy=xyso3
      if(use_quaternion_in_so3.and.(.not.use_quaternion.and.dospinr)) then
       call makequaternion(xy)
       use_quaternion=.true.
       change=.true.
      endif
      call alloc(m1);call alloc(nonl);call alloc(a1);call alloc(a2);call alloc(ri);
   
      allocate(je(nv))    
      allocate(eg(xyso3%n))
      allocate(da(c_%nd))
       da=0.0_dp
      
      m1=xy
  
      ! Brings the map to the parameter dependent fixed point
      ! including the coasting beam gymnastic: time-energy is canonical
      ! but energy is constant. (Momentum compaction, phase slip etc.. falls from there)
   ! etienne
   
    if(c_skip_gofix) then
    a1=1
  else
      call  c_gofix(m1,a1) 
    
  endif 
   
   
    m1=c_simil(a1,m1,-1)
   
   
  
      ! Does the the diagonalisation into a rotation
      call c_linear_a(m1,a2)  
   
   
    
  
      !!! Now the linear map is normalised
      m1=c_simil(a2,m1,-1)
    
   
      !!! We go into the phasors' basis
      ri=from_phasor(-1)
   
      m1=c_simil(ri,m1,1)
   
   
  !stop 999
      ri=(m1.sub.-1)**(-1) 
  
      ri%s=1  ! make spin identity
      ri%q=1.0_dp  ! make spin identity
   
   
  
      !!! The tunes are stored for the nonlinear normal form recursive algorithm
      do k=1,xy%n
        if(coast(k)) then
          eg(k)=1
         else
          je=0
          je(k)=1
          eg(k)=ri%v(k).sub.je
           if(mod(k,2)==0) then
            da(k/2)=log(abs(eg(k)))
           endif
         endif  
      enddo
   
      n%ker=0  ! In case reusing normal form
  
      do i=2,not
        if(lielib_print(13)/=0) then
          write(mdis,*) " **************************************** " 
          write(mdis,*) "Order ",i
          write(mker,*) " **************************************** " 
          write(mker,*) "Order ",i
        endif
  
        nonl=(m1*ri)
        nonl= exp_inv(n%ker,nonl)
        nonl=nonl.sub.i
  
   
  
        do k=1,xy%n
          if(lielib_print(13)/=0) then
            write(mdis,*) " **************************************** " 
            write(mdis,*) "field component ",k
            write(mker,*) " **************************************** " 
            write(mker,*) "field component ",k
          endif
  
          n%g%f(i)%v(k)=0.0_dp
          n%ker%f(i)%v(k)=0.0_dp
  
  
          j=1
  
          do while(.true.) 
  
             call  c_cycle(nonl%v(k),j,v ,je); if(j==0) exit;
             call check_kernel(k,xy%n,je,removeit)
  
             if(n%nres>0.and.removeit) then 
               do kr=1,n%nres
                 if(n%ms(kr)/=0) cycle  ! a spin resonance
                 call check_resonance(k,xy%n,je,kr,n%m,removeit)
                 if(.not.removeit) then
                   exit
                 endif
               enddo
             endif
  
            if(removeit) then
  
              lam=1.0_dp
              je(k)=je(k)-1
              do l=1,xy%n 
                if(coast(l)) cycle 
                lam=lam*eg(l)**je(l)
              enddo
  
              if(lielib_print(13)/=0) then
                   write(mdis,*) k
                   write(mdis,'(6(1x,i4))') je(1:nd2)
                   write(mdis,*) v
                   write(mdis,*) abs(v/(1-lam))
              endif
  
              je(k)=je(k)+1
  
              n%g%f(i)%v(k)=n%g%f(i)%v(k)+(v.cmono.je)/(1.0_dp-lam)
  
            else ! Put in the kernel
  
              if(lielib_print(13)/=0) then
                 je(k)=je(k)-1
                 write(mker,*) k
                 write(mker,'(6(1x,i4))') je(1:nd2)
                 write(mker,*) v
                 write(mker,*) abs(v/(1-lam))
                 je(k)=je(k)+1
              endif
                 n%ker%f(i)%v(k)=n%ker%f(i)%v(k)+(v.cmono.je)
              endif
  
          enddo  ! over monomial
        enddo  ! over vector index
  
        m1=c_simil(n%g%f(i),m1,-1)
  !call c_full_norm_vector_field(n%g%f(i),norm)
  !write(6,*) " old ",i,norm
      enddo
   
    !     if(dospinr)then
          do i=1,size(n%g%f)
           n%g%f(i)%q=0.0_dp   ! makes identity 2024.1.2
          enddo
          do i=1,size(n%ker%f)
           n%ker%f(i)%q=0.0_dp
          enddo
    !     endif
   
  
      n%a_t=a1*a2*from_phasor()*texp(n%g)*from_phasor(-1)
   !
      n%a1=a1
      n%a2=a2
  
   
  !!!!!   here we put the normalised linear part into the factored vector field
  !!!!!   not necessary but useful
      do k=1,xy%n
         if(.not.coast(k)) then    
          je=0
          je(k)=1     
          n%ker%f(1)%v(k)=n%ker%f(1)%v(k)-(log(eg(k)).cmono.je)
  
          if(mod(k,2)==1) then
              n%tune((k+1)/2)=aimag(log(eg(k)))/twopi
              n%damping((k+1)/2)=real(log(eg(k)))
              if(n%tune((k+1)/2)<0.and.n%positive) n%tune((k+1)/2)=n%tune((k+1)/2)+1.0_dp
              if(n%tune((k+1)/2)<-0.5_dp.and.(.not.n%positive)) n%tune((k+1)/2)=n%tune((k+1)/2)+1.0_dp
              if(n%tune((k+1)/2)> 0.5_dp.and.(.not.n%positive)) n%tune((k+1)/2)=n%tune((k+1)/2)-1.0_dp
  
          endif
         endif 
        enddo
  
          if(c_skip_gofix) then
           do k=1,xy%n
                    if(mod(k,2)==1) then
                       if(n%tune((k+1)/2)>0.50_dp) n%tune((k+1)/2)=n%tune((k+1)/2)-1.0_dp
                      endif
           enddo
          endif
          if(nd2t==6) then
             if(n%tune(3)>0.50_dp.and.negative_synchrotron_tune) n%tune(3)=n%tune(3)-1.0_dp
          endif 
  
         if(ndpt/=0) then
          je=0
          je(ndpt)=1
          lam=(ri%v(ndptb).sub.je) 
          n%ker%f(1)%v(ndptb)=n%ker%f(1)%v(ndptb)-(lam.cmono.je)
              if(mod(ndpt,2)==0) then
                n%tune(ndpt/2)=-lam
              else
                n%tune(ndptb/2)=-lam
              endif
         endif
  
  
      if(dospinr) then
  
  if(use_quaternion)then
        call c_full_norm_quaternion(m1%q,k,norm) 
  else
        call c_full_norm_spin(m1%s,k,norm)   
  endif
        if(k>=0) then
          dospinr=.false.
           if(use_quaternion)  then
             write(6,*) " no quaternion spin in map: dospin command ignored "
           else
              write(6,*) " no spin matrix in map: dospin command ignored "
          endif
       endif
      endif
  
  
      if(dospinr) then
        call alloc(n0) 
        call alloc(nr)
        call alloc(mt) 
        call alloc(AS) 
        call alloc(qnr)
        n%AS=1
   
  
  if(use_quaternion)then
  
        call c_normal_spin_linear_quaternion(m1,m1,n%AS,alpha)
  
        n%quaternion_angle=alpha/2.0_dp
        ri=1 ; ri%q=m1%q.sub.0 ; ! exp(theta_0 L_y)   (2)
  !      sx=sqrt(ri%q%x(1)**2+ri%q%x(2)**2+ri%q%x(3)**2)
  !      cx=ri%q%x(0)
  !write(6,*) alpha
  !      alpha=-(-2*atan2(sx,cx))
  !write(6,*) alpha
  !pause 723
        egspin(3)=cos(alpha)-i_*sin(alpha)
        egspin(2)=1.0_dp
        egspin(1)=cos(alpha)+i_*sin(alpha) 
  else
         call c_normal_spin_linear(m1,m1,n%AS,n0)  ! (1)
         ri=1 ; ri%s=m1%s.sub.0 ; ! exp(theta_0 L_y)   (2)
        egspin(3)=ri%s%s(1,1)-i_*ri%s%s(1,3)
        egspin(2)=1.0_dp
        egspin(1)=ri%s%s(1,1)+i_*ri%s%s(1,3)
  endif
   
   
  
  
  
   
        if(lielib_print(13)/=0) then
          write(mdiss,*) " eg(1:4),spin_def_tune" ,spin_def_tune
          write(mdiss,*)eg(1)
          write(mdiss,*)eg(2)
          write(mdiss,*)eg(3)
          write(mdiss,*)eg(4)
          write(mdiss,*) " egspin(1:3)" 
          write(mdiss,*)egspin(1)
          write(mdiss,*)egspin(2)
          write(mdiss,*)egspin(3)
        endif
        if(lielib_print(13)/=0) then
          write(mkers,*) " eg(1:4),spin_def_tune" ,spin_def_tune
          write(mkers,*)eg(1)
          write(mkers,*)eg(2)
          write(mkers,*)eg(3)
          write(mkers,*)eg(4)
          write(mkers,*) " egspin(1:3)" 
          write(mkers,*)egspin(1)
          write(mkers,*)egspin(2)
          write(mkers,*)egspin(3)
        endif
        
        !!! tune is taken from egspin(1) or egspin(3)   spin_def_tune= +/- 1
         n%spin_tune=aimag(log(egspin(2+spin_def_tune))/twopi)  
   
    
   
        ! because  exp(a L_y) x = x- a z + O(a**2)
         ri=ri**(-1) ! exp(-alpha_0 L_y)   (3)
  
  
  if(use_quaternion)then
         nonl=m1.sub.1 ; nonl%q=1.0_dp ;nonl=nonl**(-1)  ! R_0^-1      (4)  
  else
       nonl=m1.sub.1 ; nonl%s=1 ;nonl=nonl**(-1)  ! R_0^-1      (4)  
  endif
  !       nonl=m1.sub.1 ; nonl%s=1 ;nonl=nonl**(-1)  ! R_0^-1      (4)          
          
  
         do i=1,no    !+2
            if(lielib_print(13)/=0) then
              write(mdiss,*) " **************************************** " 
              write(mdiss,*) "Order ",i
              write(mkers,*) " **************************************** " 
              write(mkers,*) "Order ",i
            endif
    
   
           
            mt=m1*ri !  S*exp(-theta_0 L_y)    (5)
   
   
  if(use_quaternion)then
         n0=mt%q
  else
        call c_find_om_da(mt%s,n0)   ! (4)  
  endif
  
             call c_n0_to_nr(n0,n0)   ! n0 = > eigen-operator of spin   (7)
   
   
   
            n0=n0*nonl               !  no * R^-1      (8)
   
  
            nr=0
            
            do k=1,3
              if(lielib_print(13)/=0) then 
                write(mdiss,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " 
                write(mdiss,*) "Spin component ",k
                write(mdiss,*) " "
                write(mkers,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " 
                write(mkers,*) "Spin component ",k
                write(mkers,*) " "
              endif
              
              j=1
              do while(.true.) 
                call  c_cycle(n0%v(k),j,v ,je); if(j==0) exit;
                call check_kernel_spin1(k,xy%n,je,da,removeit)
                if(n%nres>0.and.removeit) then 
                  do kr=1,n%nres
                    call check_resonance_spin(k,xy%n,je,kr,n%ms,n%m,removeit)
                    if(.not.removeit) then
                      exit
                    endif
                  enddo
                endif
                  
                if(removeit) then
                  lam=egspin(k) 
                  do l=1,xy%n 
                    if(coast(l)) cycle 
                    lam=lam*eg(l)**je(l)
                  enddo
                 
                  if(lielib_print(13)/=0) then 
                    !do kr=1,nd2
    !je(kr)=-(-1)**kr*je(kr)
                    !enddo
                    write(mdiss,'(6(1x,i4))') je(1:nd2)
                    write(mdiss,*)lam
                    write(mdiss,*) v
                    write(mdiss,*) abs(v/(1-lam))
                    !do kr=1,nd2
    !je(kr)=-(-1)**kr*je(kr)
                    !enddo
                  endif
                nr%v(k)=nr%v(k) +(v.cmono.je)/(1.0_dp-lam)   ! (9)
              else
                if(lielib_print(13)/=0) then 
                  do kr=1,nd2
                    je(kr)=-(-1)**kr*je(kr)
                  enddo      
                  write(mkers,'(6(1x,i4))') je(1:nd2)
                  write(mkers,*) v
                  do kr=1,nd2
                    je(kr)=-(-1)**kr*je(kr)
                  enddo
                endif
              endif
            enddo ! cycle
          enddo ! k
   
  
          call c_nr_to_n0(nr,nr)  !   (10)
   
  
   
   
  
  if(use_quaternion)then
  qnr=nr
   AS=1 ; 
  AS%q=exp(qnr)
  else
        AS=1 ; AS%s=exp(nr)*AS%s 
  endif
  
  
   
  
          n%AS=n%AS*AS             ! (12)
   
  
   
          m1=c_simil(AS,m1,-1) 
    
  
         enddo
  
        n%AS=from_phasor()*n%AS*from_phasor(-1)
   
        n%AS=n%A_t*n%AS*n%a_t**(-1)
   
        call kill(AS) 
        call kill(mt) 
        call kill(n0) 
        call kill(nr) 
        call kill(qnr) 
      endif
        
      call c_check_rad(m1%e_ij,rad_in)
      if(rad_in) call c_normal_radiation(m1,n)
  
       ri=from_phasor()
      n%n=c_simil(ri,m1,1)
      n%Atot=n%as*n%a_t
  
    if(present(canonize)) then
     if(canonize) call c_full_canonise(n%atot,n%atot)
    endif
  
  
  
      if(present(rot)) then
        rot=n%Atot**(-1)*xy*n%Atot
      endif
      
      if(present(nu_spin)) nu_spin=0.0_dp
       if(present(phase))   then
       do i=1,size(phase)
           phase(i)=0.0_dp
        enddo
  endif
      if((present(phase).or.present(nu_spin)).and.old_phase_calculation) then
  
   
          
        if(present(rot)) then
          m1=rot
        else
          m1=n%Atot**(-1)*xy*n%Atot
        endif
  
  
  
  
            qphase=.false.
        call c_full_canonise(m1,a1,phase=phase,nu_spin=nu_spin)
        ! if(dospinr.and.present(nu_spin)) then
        !  if(real(nu_spin.sub.'0')<0) nu_spin=-nu_spin   ! 2018.11.01  to match phase advance
        ! endif
            qphase=qphasedef
      endif
  
  ! error because m1 was reutilized in present(rot)
  !   call c_check_rad(m1%e_ij,rad_in)
   !!   if(rad_in) call c_normal_radiation(m1,n)
   !   call c_check_rad(m1%e_ij,rad_in)
   !!   if(rad_in) call c_normal_radiation(m1,n)
  
      call kill(m1);call kill(nonl);call kill(a1);call kill(a2);call kill(ri);
  
        deallocate(eg)
        deallocate(da)
        deallocate(je)
  
      if(lielib_print(13)/=0) then
       close(mker)
       close(mdis)
       close(mdiss)
       close(mkers)
      endif
  
      if(change) then
       call makeso3(n%a1)
       call makeso3(n%a2)
       call makeso3(n%a_t)
       call makeso3(n%atot)
       call makeso3(n%as)
       call makeso3(n%n)
       n%a1%q=0.0_dp
       n%a2%q=0.0_dp
       n%a_t%q=0.0_dp
       n%atot%q=0.0_dp
       n%as%q=0.0_dp
       n%n%q=0.0_dp
       use_quaternion=.false.
      endif
     call kill(xy)
  !call c_count_da(i_alloc)
  !write(6,*) " exiting c_normal ",i_alloc
  inside_normal=.false.
  !!!! finding a Lie exponent
  
  !  if(use_quaternion.and.rf==0)  then
     call alloc(Nf,N_cut_2,N_nl )
     Nf=n%atot**(-1)*xyso3*n%atot
     N_cut_2=Nf.cut.(-2)
  
  ! creating the linear vector field in phasors variables
  ncoast=0
  if(c_%ndpt/=0) ncoast=1
  !!!! create a full vector field for N_cut_2
  do i=1,(c_%nd-ncoast)
   n%H_l%v(2*i-1)=-(i_*twopi*n%tune(i))*dz_c(2*i-1)-n%damping(i)*dz_c(2*i-1)
   n%H_l%v(2*i)=(i_*twopi*n%tune(i))*dz_c(2*i)-n%damping(i)*dz_c(2*i)
  enddo
  if(ncoast/=0) then
   n%H_l%v(c_%ndptb)=n%tune(c_%nd)*dz_c(c_%ndpt)
  endif
  if(dospinr) then
   n%H_l%q=0.0_dp
   n%H_l%q%x(2)=n%quaternion_angle 
  endif 
  !!!! Going into variables moving on a circle 
  n%H_l=ci_phasor()*n%H_l
  
  !!! Reverse-Dragt-Finn order for Lie map
  N_nl=N_cut_2**(-1)*nf
   
  n%H_nl=ln(N_nl)
  
  n%H_l=c_phasor()*n%H_l
  n%H_nl=c_phasor()*n%H_nl
  
   n%H=n%H_l+n%H_nl
     call kill(Nf,N_cut_2,N_nl )
  !  endif
  
  if(.not.old_phase_calculation) then
  if(present(phase)) then
  ! ncoast=0
   !if(c_%ndpt/=0) ncoast=1
   do i=1,c_%nd  !-ncoast
    phase(i)=aimag(n%h%v(2*i).k.(2*i))/twopi
   enddo
  
   if(c_%ndptb/=0) then
     phase(c_%nd)=n%h%v(c_%ndptb)   ! overwrites if necessary
   endif
  endif
  
  if(present(nu_spin))  nu_spin=-spin_def_tune*n%h%q%x(2)/pi
   
  endif
  
  
   end subroutine c_normal1 !_with_quaternion
  
  
  subroutine c_fast_canonise_clean_julia(u,u_c,phase,damping,q_cs,q_as,q_orb,q_rot,spin_tune ,dospin)
  implicit none
  type(c_damap), intent(inout) ::  u,u_c
  real(dp), optional, intent(inout) :: phase(:),damping(:)
  real(dp), optional, intent(inout) :: spin_tune(2)
  type(c_linear_map), optional :: q_cs,q_as,q_rot,q_orb   ! q_c is properly factorised
  real(dp) b(6,6),b0(6,6),ri(6,6),id(6,6),st(6,6),ang,damp(3),t,cphi,sphi,s(6,6),aq,daq
  type(c_linear_map) q ,qr,qc,qrot
  complex(dp) cri(6,6)
  integer i
  logical dos,rota
  logical, optional :: dospin
  type(c_damap) uct
  integer ndt,ndptb,ndpt,ndct
  
  ! nd = number dimensions
  ! ndt = number harmonic oscillators = 4 if coastings
  
  ndt=c_%nd2t/2
  ndptb=c_%ndptb
  ndpt=c_%ndpt 
  ndct=c_%ndc2t/2  
  write(6,*) "ndt,ndptb,ndpt,ndct"
  write(6,*) ndt,ndptb,ndpt,ndct
  pause 
  !   ndc2t=2*ndct  ! 2 if coasting, otherwise 0
  
  call alloc(uct)
   
  dos=.false.
  if(present(dospin)) then
    dos=dospin
  else
        if(force_spin_input_normal) then
          write(6,*) " your default forces you to include dospin in the input of c_fast_canonise"
          stop
        endif
  endif
  s=0
  b0=0
  id=0
  do i=1,nd
  b0(2*i-1,2*i-1)=1 ! I
  b0(2*i,2*i)=1     ! 
  s(2*i-1,2*i)=1 
  s(2*i,2*i-1)=-1 
  enddo
  if(present(q_rot) ) then 
  qrot=0  ! actually makes identity
  endif
  b=0
  id=b0
  ri=0
  b=u
   
  if(ndpt/=0)  call extract_a0_mat(b,b0)
  
  
   
     damp=1
   
   
  
        do i=1,ndt 
       ! if((i<=ndt)) then  !.or.(i>nd-rf)) then
          t=sqrt(b(2*i-1,2*i-1)**2+b(2*i-1,2*i)**2)
          cphi=b(2*i-1,2*i-1)/t
          sphi=b(2*i-1,2*i)/t
              if(sphi*b(2*i-1,2*i)+cphi*b(2*i-1,2*i-1)< 0) then    !2023.12.10    bug in China?
               cphi=-cphi
               sphi=-sphi
              endif
          
   
         ri(2*i-1,2*i-1)=cphi 
         ri(2*i,2*i)=cphi 
         ri(2*i-1,2*i)=-sphi  
         ri(2*i,2*i-1)=sphi  
     ! endif
   
   if(present(phase)) then
       ang=-atan2(sphi,cphi)
    phase(i)=phase(i)-ang/twopi
   endif
  
  
        enddo
  
  
  
  if(ndpt/=0) then
          ri(5,5)=1
          ri(6,6)=1
          ri(ndptb,ndpt)=- b(ndptb,ndpt)
          if(mod(ndpt,2)==0) then
           i=ndpt/2
          else
           i=ndptb/2
          endif
           if(present(phase))       phase(i)=phase(i)+b(ndptb,ndpt)
           if(present(q_rot) ) then 
            qrot%mat(ndptb,ndpt)=-ri(ndptb,ndpt)
            qrot%mat(5,5)=ri(5,5)
            qrot%mat(6,6)=ri(6,6)
           endif
  endif
  
  
         st=matmul(b,ri)
         st=matmul(b0,st)
  if(do_damping) then
   
  
  call canonize_damping(st,id,damp)
   
         st=matmul(st,id)
   
  if(present(damping)) then
   do i=1,ndt
    damping(i)=damping(i)+log(damp(i))
   enddo
  endif
  
  endif
  
  if(present(q_rot) ) then 
  if(ndpt/=0) then  !2023.12.10
   do i=1,2
    qrot%mat(2*i-1,2*i-1)=ri(2*i-1,2*i-1)/damp(i) 
    qrot%mat(2*i,2*i)=ri(2*i,2*i)/damp(i) 
    qrot%mat(2*i-1,2*i)=-ri(2*i-1,2*i)/damp(i) 
    qrot%mat(2*i,2*i-1)=-ri(2*i,2*i-1)/damp(i) 
   enddo
   else
   do i=1,3
    qrot%mat(2*i-1,2*i-1)=ri(2*i-1,2*i-1)/damp(i) 
    qrot%mat(2*i,2*i)=ri(2*i,2*i)/damp(i) 
    qrot%mat(2*i-1,2*i)=-ri(2*i-1,2*i)/damp(i) 
    qrot%mat(2*i,2*i-1)=-ri(2*i,2*i-1)/damp(i) 
   enddo
   endif
  endif
  
      uct=st
  if(use_quaternion.and.dos) then
  q=1
   q=u%q
   
   
   qr=1
   qr=0.0_dp
  
   
   aq=-atan2(real(q%q(2,0)),real(q%q(0,0)))
   
   qr%q(0,0)= cos(aq)
   qr%q(2,0)= sin(aq)
  
   daq=0
   if(ndpt/=0) then  
  
    daq=(q%q(0,ndpt)*qr%q(2,0)+q%q(2,ndpt)*qr%q(0,0))/(q%q(2,0)*qr%q(2,0)-q%q(0,0)*qr%q(0,0))
    qr%q(0,ndpt)=  -daq*qr%q(2,0)
    qr%q(2,ndpt)= daq*qr%q(0,0)
   endif
   qc=q*qr
   
   if(present(spin_tune).and.dos) then
    spin_tune(1)=spin_tune(1)+aq/pi   ! changed 2018.11.01
   endif
  cri=ri
  qc=qc*cri
   uct%q=qc 
  
  
  
  endif
   
   if(present(spin_tune).and.dos) then
    spin_tune(2)=spin_tune(2)+daq/pi   ! changed 2018.11.01
   endif
  
  
  
  
  
  u_c=uct
   call kill(uct)
  end subroutine c_fast_canonise_clean_julia
  
  end program example