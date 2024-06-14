program example
  use c_tpsa
  use gauss_dis
  implicit none 
  integer no,nd,i,nv,k,j(lnv),n,np,j1,ndpt
  type(c_taylor) f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin
  type(c_damap), target :: m,m1,rot,as,a0,a1,a2,a_cs,R_TE,CS_TE
  type(c_damap), pointer :: o
  type(c_vector_field) vf1,vf2,vf3  !,h_rot,K_r,K_r_real
  real(dp) r1,r2,sca,decrement
  complex(dp) c1
  type(c_normal_form) normal
  type(c_ray) ray,ray2
  type(c_lattice_function)   Lat_function
  logical :: damp = .false.,fastcanonize=.true.,COSLIKE,t_e
  type(c_taylor), allocatable :: phase(:)
  type(c_damap), allocatable :: mt(:)
  real(dp)  ph(3),spintune(2),dampdec(3)
   
  
  
  n_cai=-i_
  write(6,*) " If deltap/p0 is a canonical variable  enter 3"
  write(6,*) " else enter 2"
  !read(5,*) nd
  nd=1
  !if(nd/=2.and.nd/=3) stop 44
  if(nd==1) ndpt = 0
  if(nd==2) ndpt = 0
  if(nd==3) ndpt = 6  ! BMAD choice
  no=3;     ! no: the order of the polynomial    nv: the number of variables   
  np=2
  c_lda_used=15000
  use_quaternion=.true.
  call c_init(no,nd,np1=np,ndpt1=ndpt)  ! initializes taylor series with maps
  allocate(phase(nd))
   
  call alloc(f,F_FLOQUET,F_FLOQUET_cs,courant_snyder,nu_spin)      ! must be constructed after init
  call alloc(vf1,vf2,vf3)      
  call alloc(m,m1,as,a0,a1,a2,a_cs,rot,R_TE,CS_TE)      ! 
  call alloc(phase)
  call alloc(nu_spin)
  
  n=20
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
  do i=1,3,-1
   call c_TAYLOR_ran(vf1%q%x(i),r1,r2) ! putting spin
   vf1%q%x(i)=vf1%q%x(i)-i_*aimag(vf1%q%x(i))
   vf1%q%x(i)=0.3d0*vf1%q%x(i)*m1 !   
  
  enddo
   
  mt(k)=exp(vf1)
  
  enddo
  
  m=1
  do i=1,n
   m=mt(i)*m
  enddo
   
   
  

  call c_normal_new_no_fac(m,normal,dospin=.false.,phase=phase,nu_spin=nu_spin,doberz=.true.)
  !write(*,*) "atot =========================="
  !call print(normal%atot)
  !stop
  !call print(phase(1))
  
  call clean(phase,phase,prec=1.d-5)
  call clean(nu_spin,nu_spin,prec=1.d-5)
  
  write(6,'(A,3(1x,g23.16),/)') "The tune    is ",normal%tune(1:nd)
  write(6,'(A,3(1x,g23.16),/)') "The damping is ",normal%damping(1:nd)
  write(6,'(A,1x,g23.16,/)') "The spin tune  is ",normal%spin_tune 
   
  
  m=normal%atot**(-1)
    
  call print(m)
  stop
  
  m=ci_phasor()*m*c_phasor()
  
  
   call clean(m,m,prec=1.d-13)

   
   
  
  end program example