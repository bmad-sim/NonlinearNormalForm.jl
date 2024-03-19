program example
   use S_extend_poly
   use gauss_dis
   implicit none 
   integer no,nd,i,nv,k,j(lnv),n,np
   type(c_taylor)  F_FLOQUET,F_FLOQUET_cs,cs,nu_spin,nu_spin1
   type(c_damap) m,m1,m2,a_ori,as,a0,a1,a2,a_cs,id
   type(c_vector_field) vf1,vf2,vfr,h_rot1,h_rot2,h_rot0
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
   type(real_8) z(6),z0(6),kick(50),vkick(50),k1,k2l
   type(c_damap) z0j,z1j
   real(dp) :: lrad = 0 
   type(probe_8) xs
   type(probe) xs0
   real(dp), parameter:: a = 0.00115965218128d0
   real(dp), parameter:: gamma_0 = 40.5d0/a
   logical :: quaternion_on=.true.
   logical :: radiation_on = .false.
   character(255) filename
   n_cai=-i_
   np=2
   no=2; nd= 3;    ! no: the order of the polynomial    nv: the number of variables   
   c_lda_used=15000
   use_quaternion=.true.
    
   call c_init_all(no,nd,np,ndpt1=0)  ! initializes taylor series with maps
   allocate(phase(nd),phase1(nd))
    
   call alloc(F_FLOQUET,F_FLOQUET_cs,cs,nu_spin,nu_spin1)      ! must be constructed after init
   call alloc(vf1,vf2,vfr,h_rot1,h_rot2,h_rot0)      
   call alloc(m,m1,m2,a_ori,as,a0,a1,a2,a_cs,id)      ! 
   call alloc(phase)
   call alloc(phase1)
   call alloc(nu_spin,beta_ori,beta_lin,beta_wrong)
   
   call alloc(z)
   call alloc(z0)
   call alloc(kick)
   call alloc(vkick)
   call alloc(k1,k2l)
   call alloc(f)
   call alloc(normal)
   z0j%n=6
   z1j%n=6
   call alloc(z0j,z1j)
   call alloc(xs)
   
   
   k1=0.36d0
   k2l= 1.2d0
   !vkick(1)=morph(1.d0.mono.7)!+1.d-4
   !vkick(2)=morph(1.d0.mono.8)  !+.23d0
   c_%ndpt=0
    
   closed_orbit=0
   !call find_fix_point(closed_orbit, k1 , k2l, kick,vkick )
    
   xs0=closed_orbit
   id=1
   xs=id+xs0

   !call print(xs0)
    
    !.m = xs
    !call print(m)
   
   call track_ring(xs,xs,k1,k2l,kick,vkick)
      m=xs
      call print(m)
      stop


   filename="/Users/mgs255/.julia/dev/NonlinearNormalForm/src/mymap.txt"
   !filename="C:\document\my_tex_papers\julia\one_turn_map_no_quaternion.txt"
   
   call read_julia_map(filename,id)
    
   
   id=id.o.(m.oo.(-1))
   call print(id)
   stop
   !
   !id=id*(m**(-1))

    !call print(xs)
   
 
    call c_gofix(m,m1)
    call print(m1)
    stop


    !m%x0(1:6)=closed_orbit
   !call c_normal(m,normal,dospin=.true.,phase=phase,nu_spin=nu_spin)
   
    !call c_full_canonise(normal%atot,a0) 
   
    a1=(a0.cut.(-2))
    
    a2=a1**(-1)*a0
   
   call clean(a1,a1,prec=1.d-12)
   call clean(a2,a2,prec=1.d-12)
   call clean(normal%h,normal%h,prec=1.d-6)
   !vf2=ln(a2)
    
    
   !h_rot1=ci_phasor()*normal%h
   !h_rot1=a1*h_rot1
   !h_rot1=exp(vf2,h_rot1)
   !call print(h_rot1%v(1))
   
   !h_rot2=ci_phasor()*normal%h
   !h_rot2=a0*h_rot2
   !call print(h_rot2%v(1))
   
   !h_rot0=ci_phasor()*normal%h
   !h_rot0=a0*h_rot0
   !call print(h_rot0%v(1))
    
   !filename="C:\document\my_tex_papers\julia\one_turn_map.txt"
   !filename="C:\document\my_tex_papers\julia\one_turn_map_no_quaternion.txt"
   
   !call read_julia_map(filename,id)
    
   
   !id=id.o.(m.oo.(-1))
   !
   !id=id*(m**(-1))
   !call clean(id,id,prec=1.d-7)
   
   
   !do i=1,6
   !call print(id%v(i))
   !pause 1232
   !enddo
   
   !do i=0,3
   !!call print(id%q%x(i))
   !pause 1233
   !enddo
   
   !call clean(id,id,prec=1.d-10)
   !call print(id)
   !stop
     
   contains 
   
   
   
   
   subroutine track_qf(x0,x, k1, hkick)
   implicit none 
   type(probe_8) x0,x,z
   type(real_8)  k1,hkick,h,L 
   type(real_8) gx,sf,cf,chi,psi,zeta,alf,bet
   type(real_8)  kx,w_x,w_y,sx,cx,sy,cy,sig,xi   
   type(quaternion_8) q,Qkick
   real(dp) lbend,lrad(3)
     call alloc(h,L)
     call alloc(z)
     call alloc(gx,sf,cf,chi,psi,zeta,alf,bet)
     call alloc( kx,w_x,w_y,sx,cx,sy,cy,sig,xi)
     call alloc(q)
     call alloc(Qkick)
     
    
    lbend=0.1d0
   
   
     L = 0.5d0/(1.d0+x0%x(6))
     h=-L*((x0%x(2)**2+k1*x0%x(1)**2)+(x0%x(4)**2-k1*x0%x(3)**2))/(1.d0+x0%x(6))/2.d0
     z%x(1) =  cos(sqrt(k1)*L)*x0%x(1)           + 1.d0 /sqrt(k1)*sin(sqrt(k1)*L)*x0%x(2)
     z%x(2) =  -sqrt(k1)*sin(sqrt(k1)*L)*x0%x(1) + cos(sqrt(k1)*L)*x0%x(2) + hkick 
     z%x(3) =  cosh(sqrt(k1)*L)*x0%x(3)+1/sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(4)
     z%x(4) =  sqrt(k1)*sinh(sqrt(k1)*L)*x0%x(3) + cosh(sqrt(k1)*L)*x0%x(4)
     z%x(5)=x0%x(5)+h  
     z%x(6)=x0%x(6)
     z%x(2) = z%x(2) +lbend*z%x(6)
     z%x(5)=  z%x(5)-lbend*z%x(1)
      
   
    if(quaternion_on) then
     gx =   lbend/L
     sf =  sin(a*lbend*gamma_0/2.d0)
     cf =   cos(a*lbend*gamma_0/2.d0)
     chi = 1d0+a*gamma_0
     psi = gamma_0**2-1d0
     zeta = gamma_0-1d0
     alf =  2*(a**2*gamma_0**2*gx**2+k1)
     bet =  a*gx*k1*(gamma_0*chi-zeta)
     kx =   k1+gx**2
     w_x =  sqrt(kx)
     w_y =  sqrt(k1)
     sx =   sin(L*w_x)
     cx =   cos(L*w_x)
     sy =   sinh(L*w_y)
     cy =   cosh(L*w_y)
     sig =  w_y*(k1 + a*k1*gamma_0 + a**2*gx**2*zeta*gamma_0)
     xi =   w_y*(k1*chi + a**2*gx**2*zeta*gamma_0)
   
     q%x(0) =  cf - x0%x(1)*kx*chi/(2*w_x)*sx*sf - x0%x(2)*kx*chi/(2*w_x**2)*(1-cx)*sf + x0%x(6)*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*sf
     q%x(1) =  x0%x(3)*-1/alf*(bet*(1+cy)*sf + sig*sy*cf) + x0%x(4)*1/(w_y*alf)*(xi*(1-cy)*cf-bet*sy*sf)
     q%x(2) =  -sf - x0%x(1)*kx*chi/(2*w_x)*sx*cf - x0%x(2)*kx*chi/(2*w_x**2)*(1-cx)*cf + x0%x(6)*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*cf
     q%x(3) =  x0%x(3)/alf*(bet*(1-cy)*cf + sig*sy*sf) + x0%x(4)/(w_y*alf)*(xi*(1+cy)*sf-bet*sy*cf)
   
   gx=sqrt(q%x(0)**2+q%x(1)**2+q%x(2)**2+q%x(3)**2)
   do i=0,3
    q%x(i)=q%x(i)/gx
   enddo
   !call print(q)
   !pause 1
   !call print(x0%q)
   !pause 2
   z%q=q*x0%q
   !call print(z%q)
   !pause 3
   !cos(hkick*(1+a*gamma_0)/2),0, sin(hkick*(1+a*gamma_0)/2), 0]
   !(Q, p.Q, Q)
   Qkick=0.d0
   Qkick%x(0)=cos(hkick*(1+a*gamma_0)/2)
   Qkick%x(2)=sin(hkick*(1+a*gamma_0)/2)
   z%q=Qkick*z%q
   !call print(Qkick)
   !pause 4
   !call print(z%q)
   !pause 5
   endif
   if (radiation_on) then
     do i=1,3 
      z%x(2*i-1)=exp(lrad(i)*(1.d0+z%x(2*i-1)**2))* z%x(2*i-1)
   !   z%x(2*i-1)=exp(lrad(i))* z%x(2*i-1)
     enddo
    endif
    
   
     x=z
    
     call kill(h,L)
     call kill(z)
     call kill(gx,sf,cf,chi,psi,zeta,alf,bet)
     call kill( kx,w_x,w_y,sx,cx,sy,cy,sig,xi)
     call kill(q)
     call kill(Qkick)
   
   end subroutine
   
   subroutine track_qd(x0,x, k1, vkick)
   implicit none 
   type(probe_8) x0,x,z
   type(real_8)    k1,vkick,h,L
   type(real_8)  gx,sf,cf,chi,psi,zeta,alf,bet
   type(real_8) kx,w_x,w_y,sx,cx,sy,cy,sig,xi
   type(quaternion_8)  q,Qkick
     call alloc(gx,sf,cf,chi,psi,zeta,alf,bet)
     call alloc( kx,w_x,w_y,sx,cx,sy,cy,sig,xi)
     call alloc(q)
     call alloc(Qkick)
   
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
    
    if(quaternion_on) then
     gx =   0.d0
     sf = 0.d0
     cf =  1.d0
     chi = 1d0+a*gamma_0
     psi = gamma_0**2-1d0
     zeta = gamma_0-1d0
     alf =  2*k1
     bet =  0
     kx =   k1 
     w_x =  sqrt(kx)
     w_y =  w_x
     sx =   sinh(L*w_x)
     cx =   cosh(L*w_x)
     sy =   sin(L*w_y)
     cy =   cos(L*w_y)
     sig = w_y*k1*chi
     xi =  sig
   
       q%x(0) =  cf - x0%x(1)*kx*chi/(2*w_x)*sx*sf + x0%x(2)*kx*chi/(2*w_x**2)*(1-cx)*sf + x0%x(6)*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*sf
       q%x(1) =  x0%x(3)*-1/alf*(bet*(1+cy)*sf - sig*sy*cf) + x0%x(4)*1/(w_y*alf)*(xi*(1-cy)*cf-bet*sy*sf)
       q%x(2) =  -sf - x0%x(1)*kx*chi/(2*w_x)*sx*cf + x0%x(2)*kx*chi/(2*w_x**2)*(1-cx)*cf + x0%x(6)*gx/2*(chi*sx/w_x - a*L*psi/gamma_0)*cf
       q%x(3) =  x0%x(3)/alf*(bet*(1-cy)*cf - sig*sy*sf) + x0%x(4)/(w_y*alf)*(xi*(1+cy)*sf-bet*sy*cf)
    
   gx=sqrt(q%x(0)**2+q%x(1)**2+q%x(2)**2+q%x(3)**2)
   do i=0,3
    q%x(i)=q%x(i)/gx
   enddo
   
   z%q=q*x0%q
   !cos(hkick*(1+a*gamma_0)/2),0, sin(hkick*(1+a*gamma_0)/2), 0]
   !(Q, p.Q, Q)
   Qkick=0.d0
   Qkick%x(0)=cos(vkick*(1.d0+a*gamma_0)/2)
   Qkick%x(1)=sin(vkick*(1.d0+a*gamma_0)/2)
   z%q=Qkick*z%q
   endif
   
     x=z
    
     call kill(h,L)
     call kill(z)
      call kill(gx,sf,cf,chi,psi,zeta,alf,bet)
     call kill( kx,w_x,w_y,sx,cx,sy,cy,sig,xi)
     call kill(q)
     call kill(Qkick)
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
      x%q=x0%q
   
   end subroutine
   
   subroutine track_sextupole(x0,x, k2l)
   implicit none 
   type(real_8)  k2l,gx
   type(probe_8) x0,x 
   type(quaternion_8)  q
   call alloc(q)
   call alloc(gx)
   
     x%x(2) =  x0%x(2) -k2l *(x0%x(1)**2 - x0%x(3)**2)
     x%x(4) =  x0%x(4) +k2l*2.d0*x0%x(1)*x0%x(3) 
     x%x(1) =  x0%x(1) 
     x%x(3) =  x0%x(3) 
     x%x(5) =  x0%x(5) 
     x%x(6) =  x0%x(6) 
         x%q=x0%q
   
   if(quaternion_on) then
    q%x(0)=1.d0
    q%x(1)=-k2l*2.d0*x0%x(1)*x0%x(3)*(1.d0+a*gamma_0)/2.d0
    q%x(2)=-k2l*(x0%x(1)**2-x0%x(3)**2)*(1+a*gamma_0)/2
    q%x(3)=0
   
   gx=sqrt(q%x(0)**2+q%x(1)**2+q%x(2)**2+q%x(3)**2)
   do i=0,3
    q%x(i)=q%x(i)/gx
   enddo
    
   x%q=q*x0%q
    
   endif
   call kill(q)
   call kill(gx)
   
   end subroutine
   
   
   subroutine track_fodo(z0,z, k1, k2l, kick, vkick)
   implicit none
   type(real_8) k2l, k1,kick,vkick
   integer i
   type(probe_8) z,z0 
   
   
      call track_qf(z0,z, k1, kick)
      call track_sextupole(z,z0, k2l)
      call track_drift(z0,z)
      call track_qd(z,z0, k1, vkick)
      call track_sextupole(z0,z, -k2l)
      call track_drift(z,z0)
         z=z0
    
   end subroutine
   
   subroutine track_ring(z0,z, k1 , k2l, kick,vkick )
   implicit none
   type(real_8)  k2l, k1,kick(:),vkick(:)
   type(probe_8) z,z0
   integer i,j
   
     do i=1,50
    
        call track_fodo(z0,z, k1, k2l, kick(i), vkick(i))
        z0=z
     enddo
     if(c_%ndpt==0) call track_cav(z)
   
   end subroutine
   
   subroutine track_cav(z0)
   implicit none
   type(probe_8) z0
     z0%x(6) =  z0%x(6) + 0.0001d0*z0%x(5)
   
    end subroutine
   
   
   subroutine find_fix_point(fix, k1 , k2l, kick,vkick )
   implicit none
   type(real_8)  k2l, k1,kick(:),vkick(:)
   real(dp) fix(6)
   type(probe_8)  z,z0
   type(probe)   p0
   type(c_damap) id,m
   type(c_ray) ray
   integer k
   call alloc(z0)
   call alloc(z)
   call alloc(id,m)
   
   do k=1,10
   123 p0=fix
   id=1
   z0=p0+id
   call track_ring(z0,z, k1 , k2l, kick,vkick )
   m=z
   ray=0
   if(c_%ndpt/=0) then
    do i=1,4
     ray%x(i)=fix(i)-(m%v(i).sub.'0')
    enddo
   else
    do i=1,6
     ray%x(i)=fix(i)-(m%v(i).sub.'0')
   enddo
   endif
    
   
    
   m=id-m
   
   
    
   if(c_%ndpt/=0) then
   m%v(5)=1.d0.cmono.5
   m%v(6)=1.d0.cmono.6
   endif
   id=m**(-1)
   
   
   ray=id.o.ray
   
    
   if(c_%ndpt/=0) then
    do i=1,4
     fix(i)=fix(i)-ray%x(i)
    enddo
   else
    do i=1,6
     fix(i)=fix(i)-ray%x(i)
    enddo
   endif
   
   write(6,format6) fix(1:6)
   enddo
   !write(6,*)  " more "
   !read(5,*) i
   !if(i==1) goto 123
    
   call kill(id,m)
   call kill(z0)
   call kill(z)
   
   end subroutine
   
   
   subroutine read_julia_map(filename,id)
   implicit none
   integer i,mf,k 
   character(*) filename
   character(255) line,a5
   real(dp) x(6),cc
   integer, allocatable :: J(:)
   type(c_damap) id
   
   id=0
   
   allocate(J(1:c_%nv))
   
   call  kanalnummer(mf,filename)
    read(mf,'(a255)') line
    read(mf,'(a255)') line
    
   
   do i=1,6
    read(mf,'(a255)') line
    write(6,'(a255)') line
    read(line(6:40),*)  x(i)
    write(6,*)  x(i)
   enddo
   id%x0(1:6)=x(1:6)
   
    read(mf,'(a255)') line
    read(mf,'(a255)') line
    read(mf,'(a255)') line
    read(mf,'(a255)') line
   k=1
    do i=1,10000
    read(mf,'(a255)',end=1001) line
    write(6,'(a255)') line
   
   if(index(line,"---")/=0)  then
    k = k+1
   if(index(line,"endendend")/=0) then
    goto 1001 
   endif
   elseif(line(1:10)=='          ') then
    read(mf,'(a255)') line
   if(index(line,"Quaternion")/=0) then
    goto 1000 
   endif
    write(6,'(a255)') line
     read(mf,'(a255)') line
    write(6,'(a255)') line
     read(mf,'(a255)') line
    write(6,'(a255)') line
     read(mf,'(a255)') line
    write(6,'(a255)') line
   else
     read(line(7:34),*) cc
     read(line(41:100),*) j(1:6)
     read(line(69:100),*) j(7:c_%nv)
      id%v(k) = id%v(k) + (cc.cmono.j)
    
     write(6,*) cc
     write(6,'(10(1x,i4))') k,j
    
   
   endif
   
    
   enddo
   
   
   
   1000 continue
    read(mf,'(a255)') line
    write(6,'(a255)') line
    read(mf,'(a255)') line
    write(6,'(a255)') line
   
   k=0
    do i=1,10000
    read(mf,'(a255)',end=1001) line
    write(6,'(a255)') line
   
    if(index(line,"endendend")/=0) then
     goto 1001 
    endif
   
   if(index(line,"---")/=0)  then
    k = k+1
    if(k==4) goto 1001
   else
     read(line(7:34),*) cc
     read(line(41:100),*) j(1:6)
     read(line(69:100),*) j(7:c_%nv)
      id%q%x(k) = id%q%x(k) + (cc.cmono.j)
    
     write(6,*) cc
     write(6,'(10(1x,i4))') k,j
   
   
   endif
   
    
   enddo
   
   1001 deallocate(j)
   
   close(mf)
   end subroutine
   
   
   end program example