program matt
  use pointer_lattice
  use gauss_dis
  
  implicit none
  type(layout), pointer:: ring
  real(dp) prec,x_closed,x_closed1 
  real(dp) x0,x1,x2,x_closed_tpsa
  type(c_ray) f_ray
  type(c_damap)  one_turn_map_AP, identity_AP,two_turn_map_AP,lagrange_map_ap,go_to_orbit
  type(c_damap) map1,map2,map12, f_map,c_w0_inv,c_w0
  type(c_taylor) x
  integer no,nv,np
  integer i,j,k  !,nd1,ic(4)
  real(dp)  m(2,2)
  prec=1.d-10 ! for printing
  
  write(6,*) " give order no>=1 "
  read(5,*) no

  nV=1 ;  ! nv=7+12=19
  call c_init_all(no,nv)  
  map1%n=nV   ! size of the c_damap must be explicitely stated
  map2%n=nV   ! size of the c_damap must be explicitely stated
  map12%n=nV   ! size of the c_damap must be explicitely stated
  one_turn_map_AP%n=nV
  go_to_orbit%n=nV
  call alloc(x); call alloc(map1,map2,map12,go_to_orbit,one_turn_map_AP);
  
  x1= 0.015d0
  map1%x0(1)=x1
  x=map1%x0(1)+(1.d0.cmono.1)
  map1%v(1)=sin(x/2.d0)+0.3d0*sin(x)**2+0.05d0
  write(*,*) "Map 1:"
  call print(map1)

  x2=.02d0
  map2%x0(1)=x2
  x=map2%x0(1)+(1.d0.cmono.1)
  map2%v(1)=sin(0.3d0*x)+0.2d0*sin(x)**2+0.03d0
  write(*,*) "Map 2:"
  call print(map2)

  !!! multiplying the TPSA maps
  map12=map2.o.map1
  !!map12=map1*map2
  write(*,*) "Map 2 .o. Map 1:"
  call print(map12)
  
end program matt