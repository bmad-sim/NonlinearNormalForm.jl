program example
  use c_tpsa
  use gauss_dis
  implicit none 
  integer no,nd,i,nv,k,j(lnv),n,np,j1,ndpt,i1,i2
  type(c_taylor) f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin
  type(c_damap), target :: m,m1,rot,as,a0,a1,a2,a_cs,R_TE,CS_TE
  type(c_damap), pointer :: o
  type(c_vector_field) vf1,vf2,vf3  !,h_rot,K_r,K_r_real
  real(dp) r1,r2,sca,decrement(6),radfluc,radx
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
  putspin=.false.
  n_cai=-i_
  write(6,*) " If deltap/p0 is a canonical variable  enter 3"
  write(6,*) " else enter 2"
  !read(5,*) nd
  nd=3  ! coasting
  !if(nd/=2.and.nd/=3) stop 44
  if(nd==1) ndpt = 0
  if(nd==2) ndpt = 0
  if(nd==3) ndpt = 0  !6  ! BMAD choice
  no=3;     ! no: the order of the polynomial    nv: the number of variables   
  np=2
  c_lda_used=1500
   call gaussian_seed(121)
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
  do i=1,nd
  sca=0
  do while(sca==0.d0) 
    call GRNF(sca,1.d0)
  enddo
  !sca=sign(sca,1.d0)+sca  
  sca=abs(sca)*1000
   f=f+sca*(dz_c(2*i-1)**2+dz_c(2*i)**2)
  enddo
   
  ! Removing time terms for coasting
  if(c_%ndpt/=0) then
   m1=1
   m1%v(2*nd-1)=0.0_dp
   f=f*m1
  endif
  
  f=-f*0.002d0
   
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
  
  do i1=1,c_%nd2
  do i2=i1,c_%nd2
   call GRNF(radx,4.d0)
   radx=radx*radfluc
   mt(k)%e_ij(i1,i2)=abs(radx)
   mt(k)%e_ij(i2,i1)=abs(radx)
  enddo
  enddo
  
  
  enddo
  
  m=1
  do i=1,n
   m=mt(i)*m
  enddo
  
  decrement(1)=0.99999d0
  decrement(2)=0.99997d0
  decrement(3)=0.99996d0
  decrement(4)=0.999995d0
  decrement(5)=0.99993d0
  decrement(6)=0.99999d0
  
  decrement(1)=0.9d0
  decrement(2)=0.7d0
  decrement(3)=0.6d0
  decrement(4)=0.95d0
  decrement(5)=0.3d0
  decrement(6)=0.9d0
  
   
  do i=1,c_%nd2
  m%v(i)=m%v(i)*decrement(i)
  enddo
 m1=3.d0 + m

   do i=1,6
    write(6,format6) real(m1%e_ij(i,1:6) )
  enddo
  stop

! do i=1,1000
!  m1  = m1*m1
! enddo
 !call print(m1)
   
!do i=1,6
!  write(6,format6) real(m1%e_ij(i,1:6) )
!enddo
!stop

  
  call c_normal(m,normal,dospin=putspin,phase=phase,nu_spin=nu_spin)
  !m1 =ci_phasor()*normal%atot**(-1)*m*normal%atot*c_phasor()
     
  !m = ci_phasor()*normal%atot**(-1)*m*normal%atot*c_phasor()
  !do i=1,6
  !  write(6,format6) real(%e_ij(i,1:6))
  !enddo
  
  !write(*,*) "\n"
  !do i=1,6
  !  write(6,format6) imag(m1%e_ij(i,1:6))
 ! enddo
 ! stop
   


  do i=1,6
    write(6,format6) real(normal%s_ij0(i,1:6) )
  enddo
  stop

  !write(*,*) "emit_1 = ", normal%emittance(1)
  !write(*,*) "emit_2 = ", normal%emittance(2)
  !write(*,*) "emit_3 = ", normal%emittance(3)
  !stop
  !call print(phase(1))
  
  !call clean(phase,phase,prec=1.d-5)
  !call clean(nu_spin,nu_spin,prec=1.d-5)
  
  write(6,'(A,3(1x,g23.16),/)') "The tune    is ",normal%tune(1:nd)
  write(6,'(A,3(1x,g23.16),/)') "The damping is ",normal%damping(1:nd)
  write(6,'(A,1x,g23.16,/)') "The spin tune  is ",normal%spin_tune 

    !m=ci_phasor()*m*c_phasor()
  
   
   call clean(m,m,prec=1.d-13)

   call print(m)
   stop
  
   call print(normal%h_l)
   pause 123
   call clean(normal%h_nl, normal%h_nl, prec=1.d-13)
   call print(normal%h_nl)
   pause 124
   
   do i1=1,c_%nd2
   do i2=i1,c_%nd2
    write(6,*) normal%s_ij0(i1,i2)
   enddo
  enddo 
  
    call compute_lattice_functions(normal%atot,cl)
  
  write(6,*) cl%e(1,1,1)*normal%emittance(1)+cl%e(2,1,1)*normal%emittance(2)+cl%e(3,1,1)*normal%emittance(3)
  contains
  
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

  m1=m.cut.2
  m1=m1**(-1)
  
  nl=m1*m
  dl=nl
   call print(nl)
  
  
  
  mi=1
  id=1
  do i=1,nt
  dl=nl-id
  dl=id-dl
  call clean(nl,nl,prec=1.d-10)
   call print(nl%v(1))
   call print(nl%q%x(1))
  write(6,*) i,nt
  
  
  nl=dl*nl
  
  mi=dl*mi
  
  enddo
  
  mi=mi*m1
  call kill(m1,nl,dl,id)
  
  end subroutine invert_with_1minusx
  
  
  end program example