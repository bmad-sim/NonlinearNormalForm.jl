program example
  use S_extend_poly
  use gauss_dis
  implicit none 
  integer no,nd,i,nv,k,j(lnv),n,np
  type(c_taylor)  F_FLOQUET,F_FLOQUET_cs,cs,nu_spin,nu_spin1
  type(c_damap) m,m1,m2,a_ori,as,a0,a1,a2,a_cs,id
  type(c_vector_field) vf1,vf2,vf3  !,h_rot,K_r,K_r_real
  real(dp) r1,r2,sca,decrement,closed_orbit(6),normb
  complex(dp) c1,mat(6,6)
  type(c_normal_form) normal
  type(c_ray) ray,ray2
  type(c_lattice_function)   Lat_function
  logical :: damp = .false.,fastcanonize=.true.
  type(c_taylor), allocatable :: phase(:),phase1(:)
  type(c_damap), allocatable :: mt(:)
  real(dp)  ph(3),spintune(2),dampdec(3)
  type(c_taylor) beta_ori,beta_lin,beta_wrong
  type(taylor) f
  type(real_8) z(6),z0(6),kick(50),k1,k2l
  type(c_damap) z0j,z1j
  real(dp) :: lrad = 0 
  type(probe_8) xs
  type(probe) xs0
  n_cai=-i_
  np=2
  no=2; nd= 3;    ! no: the order of the polynomial    nv: the number of variables   
  c_lda_used=15000
  use_quaternion=.true.
  call c_init_all(no,nd,np,ndpt1=0)  ! initializes taylor series with maps
  allocate(phase(nd),phase1(nd))
   
  call alloc(F_FLOQUET,F_FLOQUET_cs,cs,nu_spin,nu_spin1)      ! must be constructed after init
  call alloc(vf1,vf2,vf3)      
  call alloc(m,m1,m2,a_ori,as,a0,a1,a2,a_cs,id)      ! 
  call alloc(phase)
  call alloc(phase1)
  call alloc(nu_spin,beta_ori,beta_lin,beta_wrong)
  
  call alloc(z)
  call alloc(z0)
  call alloc(kick)
  call alloc(k1,k2l)
  call alloc(f)
  call alloc(normal)
  z0j%n=6
  z1j%n=6
  call alloc(z0j,z1j)
  call alloc(xs)
  closed_orbit= 0.001d0
   
  xs0=closed_orbit
  id=1

  ! +closed_orbit

  ! id is a DAMAP!!!!!   
  ! here the DAMap is "inserted" into the probe
  xs=id+xs0

  ! in this step, xs0 is added both to x0 of probe and scalar part of TPSA in x
    
 ! call print(x)
 ! stop

  ! ins


  
  k1=0.36d0
  k2l= 1.2d0
      ! xs0 probe gets promoted to probe_8 here
  !call print(xs)
  write(*,*) 'BEFORE: x0 elements:'
   do i = 1, size(xs%x0)
      call print(xs%x0(i))
   end do
   write(*,*) 'BEFORE: x elements:'
   do i = 1, size(xs%x)
      call print(xs%x(i))
   end do

  call track_ring(xs,xs,k1,k2l,kick)
  write(*,*) 'AFTER: x0 elements:'
  do i = 1, size(xs%x0)
     call print(xs%x0(i))
  end do
  write(*,*) 'AFTER: x elements:'
  do i = 1, size(xs%x)
     call print(xs%x(i))
  end do
  stop
  m = xs

  call print(m)
  !call print(xs)
!  call print(xs)

  !m = xs


  
   
  kick(1)=0.0d0! morph(1.d0.mono.7)  !+0.1d0
  kick(2)=0.0d0!morph(1.d0.mono.8)  !+.23d0
  ! call track_drift(xs, xs) ! (xs, xs, k1, kick(1))
  !m=xs
  !call print(m )
  !stop

  !call print(xs)
  print *, k1
  stop

  m1 = m**(-1)
  m2 = m.oo.(-1)
  print *, m%x0
  call print(m)
  stop
  call print(m1)
  call print(m2)
 !  call checksymp(m,normb)
  
 ! write(6,*) normb
    
   
  stop
  
   
  
   
  !call c_normal(m,normal,dospin=.true.,phase=phase,nu_spin=nu_spin)
  
  
   
   
   
  call c_normal(m,normal,phase=phase)  !,dospin=.false.,phase=phase)
  
  call c_full_factorise(normal%atot,as,a0,a1,a2,dir=-1) 
  
  call print(a0)
  
  call clean(phase,phase,prec=1.d-10)
  call print(phase)
  
  
  contains 
  
  
  subroutine remove_nonlinear(m)
  implicit none
  type(c_damap) m
  type(c_taylor) t
  integer i,k,n,it,i1
  integer,allocatable :: j(:)
  complex(dp) value
  
  
  call alloc(t)
  allocate(J(c_%nv))
  do k=1,m%n
     call c_taylor_cycle(m%v(k),size=n)
   t=0.0_dp
      do i=1,n
         call c_taylor_cycle(m%v(k),ii=i,value=value,j=j)
  it=0
       do i1=1,c_%nd2
          it=it+j(i1)
       enddo
  
       if (it==1.or.it==0) then
        t=t+(value.cmono.j)
       endif
   !      x=real(value)
    !     c_real=c_real+(x.cmono.j)
     enddo
     m%v(k)=t
  enddo
  
  deallocate(J)
  call kill(t)
  
  end subroutine
  
  subroutine remove_parameters(m)
  implicit none
  type(c_damap) m
  type(c_taylor) t
  integer i,k,n,it,i1
  integer,allocatable :: j(:)
  complex(dp) value
   
  
  call alloc(t)
  allocate(J(c_%nv))
  do k=1,m%n
     call c_taylor_cycle(m%v(k),size=n)
   t=0.0_dp
      do i=1,n
         call c_taylor_cycle(m%v(k),ii=i,value=value,j=j)
  it=0
       do i1=c_%nd2+1,c_%nv
          it=it+j(i1)
       enddo
  
       if (it<=0) then
        t=t+(value.cmono.j)
       endif
   !      x=real(value)
    !     c_real=c_real+(x.cmono.j)
     enddo
     m%v(k)=t
  enddo
  
  deallocate(J)
  call kill(t)
  
  end subroutine
  
  subroutine track_qf(x0,x, k1, hkick)
  implicit none 
  type(probe_8) x0,x,z
  type(real_8)  k1,hkick,h,L 
  real(dp) lbend
    call alloc(h,L)
    call alloc(z)
   
   lbend=0.1d0
  
  
    L = 0.5d0/(1.d0+x0%x(6))
    h=-L*((x0%x(2)**2+k1*x0%x(1)**2)+(x0%x(4)**2-k1*x0%x(3)**2))/(1.d0+x0%x(6))/2.d0
    z%x(1) =  cos(sqrt(k1)*L)*x0%x(1)           + 1.d0 /sqrt(k1)*sin(sqrt(k1)*L)*x0%x(2)
    z%x(2) =  -sqrt(k1)*sin(sqrt(k1)*L)*x0%x(1) + cos(sqrt(k1)*L)*x0%x(2) + hkick 
    z%x(3) =  cosh(sqrt(k1)*L)*x0%x(3)          + 1.d0  /sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(4)
    z%x(4) =  sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(3) + cosh(sqrt(k1)*L)*x0%x(4)
    z%x(5)=x0%x(5)+h  
    z%x(6)=x0%x(6)
    z%x(2) = z%x(2) +lbend*z%x(6)
    z%x(5)=  z%x(5)-lbend*z%x(1)
  
    x=z
    call kill(h,L)
    call kill(z)
  end subroutine
  
  subroutine track_qd(x0,x, k1, vkick)
  implicit none 
  type(probe_8) x0,x,z
  type(real_8)  k1,vkick,h,L
  
    call alloc(h,L)
    call alloc(z)
  
   
    L = 0.5d0/(1.d0+x0%x(6))
    h=-L*((x0%x(2)**2-k1*x0%x(1)**2)+(x0%x(4)**2+k1*x0%x(3)**2))/(1.d0+x0%x(6))/2.d0
    z%x(1) =  cosh(sqrt(k1)*L)*x0%x(1)  + 1.d0 /sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(2)
    z%x(2) =  sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(1) + cosh(sqrt(k1)*L)*x0%x(2) 
    z%x(3) =  cos(sqrt(k1)*L)*x0%x(3)+ 1.d0  /sqrt(k1)*sin(sqrt(k1)*L)*x0%x(4)
    z%x(4) =  -sqrt(k1)*sin(sqrt(k1)*L)*x0%x(3) + cos(sqrt(k1)*L)*x0%x(4) + vkick
    z%x(5)=x0%x(5)+h
    z%x(6)=x0%x(6)
  
    x=z
    call kill(h,L)
    call kill(z)
   
  end subroutine
  
  subroutine track_drift(x0,x)
  implicit none 
  type(probe_8) x0,x 
  real(dp) L
    
    L = 0.75d0
    x%x(1) =  x0%x(1)+x0%x(2)*L/(1.d0+x0%x(6))
    x%x(3) =  x0%x(3)+x0%x(4)*L/(1.d0+x0%x(6))
    x%x(2) =  x0%x(2) 
    x%x(4) =  x0%x(4) 
    x%x(5) =  x0%x(5) -L*((x0%x(2)**2)+(x0%x(4)**2))/(1.d0+x0%x(6))**2/2.d0
    x%x(6) =  x0%x(6) 
   
  end subroutine
  
  subroutine track_sextupole(x0,x, k2l)
  implicit none 
  type(real_8)  k2l
  type(probe_8) x0,x 
  
    x%x(2) =  x0%x(2) -k2l *(x0%x(1)**2 - x0%x(3)**2)
    x%x(4) =  x0%x(4) +k2l*2.d0*x0%x(1)*x0%x(3) 
    x%x(1) =  x0%x(1) 
    x%x(3) =  x0%x(3) 
    x%x(5) =  x0%x(5) 
    x%x(6) =  x0%x(6) 
   
  end subroutine
  
  
  subroutine track_fodo(z0,z, k1, k2l, kick)
  implicit none
  type(real_8) k2l, k1,kick,vkick
  integer i
  type(probe_8) z,z0 
  
  call alloc(vkick)
  vkick=0.d0
    call track_qf(z0,z, k1, kick)
     call track_sextupole(z,z0, k2l)
     call track_drift(z0,z)
     call track_qd(z,z0, k1, vkick)
     call track_sextupole(z0,z, -k2l)
     call track_drift(z,z0)
        z=z0
  call kill(vkick)
  end subroutine
  
  subroutine track_ring(z0,z, k1 , k2l, kick )
  implicit none
  type(real_8)  k2l, k1,kick(:)
  type(probe_8) z,z0
  integer i,j
  
    do i=1,2
   
       call track_fodo(z0,z, k1, k2l, kick(i))
  
    enddo
    if(c_%ndpt==0) call track_cav(z)
  
  end subroutine
  
  subroutine track_cav(z0)
  implicit none
  type(probe_8) z0
    z0%x(6) =  z0%x(6) + 0.0001d0*z0%x(5)
  
   end subroutine
  
  end program example