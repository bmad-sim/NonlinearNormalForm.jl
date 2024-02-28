program tpsa_1d_map
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
call c_init_all(no,nv)   ! FPP command to initialize the TPSA
map1%n=nV   ! size of the c_damap must be explicitely stated
map2%n=nV   ! size of the c_damap must be explicitely stated
map12%n=nV   ! size of the c_damap must be explicitely stated
one_turn_map_AP%n=nV
go_to_orbit%n=nV
call alloc(x); call alloc(map1,map2,map12,go_to_orbit,one_turn_map_AP);
f_map%n=nv
c_w0_inv%n=nv
c_w0%n=nv
call alloc(f_map,c_w0_inv,c_w0)
 !!!! DA

x_closed=.05469119581164052d0
write(6,*) "Imitating PTC: tracking through  "

x=x_closed+(1.d0.cmono.1)
x=sin(x/2.d0)+0.3d0*sin(x)**2+0.05d0     
x_closed1=x	! recording orbit for future use
x=sin(0.3d0*x)+0.2d0*sin(x)**2+0.03d0

call print(x)

write(6,*) "Creating two DA maps : map1 and map2 "

! computing map1
x=x_closed+(1.d0.cmono.1)
x=sin(x/2.d0)+0.3d0*sin(x)**2+0.05d0  
map1%v(1)=x
! computing map2
x=x_closed1+(1.d0.cmono.1)
x=sin(0.3d0*x)+0.2d0*sin(x)**2+0.03d0
map2%v(1)=x
one_turn_map_AP = map2*map1
Write(6,*) " This is map1 "
call print(map1)
Write(6,*) " This is map2 "
call print(map2)
write(6,*) "one_turn_map_AP = map2*map1 "
call print(one_turn_map_AP)
!!! Saving DA coefficients to imitate a TRANSPORT-like code
m(1,1)=map1%v(1).sub.'1'
m(1,2)=map1%v(1).sub.'2'
m(2,1)=map2%v(1).sub.'1'
m(2,2)=map2%v(1).sub.'2'
write(6,*) "Coefficients using Taylor Map Multiplication "
write(6,*) m(1,1)*m(2,1), m(2,1)*m(1,2)+m(2,2)*m(1,1)**2

!  Creating TPSA maps around the design orbit 

x1= 0.015d0
map1%x0(1)=x1
x=map1%x0(1)+(1.d0.cmono.1)
map1%v(1)=sin(x/2.d0)+0.3d0*sin(x)**2+0.05d0

x2=.02d0
map2%x0(1)=x2
x=map2%x0(1)+(1.d0.cmono.1)
map2%v(1)=sin(0.3d0*x)+0.2d0*sin(x)**2+0.03d0

!!! multiplying the TPSA maps
 map12=map2.o.map1
  
 
x=map12%v(1)+i_*one_turn_map_AP%v(1)
 write(6,*)  "    The map map12 for one turn"
 call print(map12)
 write(6,*)  "    Comparing the TPSA one-turn map with the DA one-turn map"
 write(6,*)  "    TPSA coefficients      DA coefficient ",k
 call print(x)
  

write(6,*)
write(6,*)
write(6,*)
! TPSA maps re-expressed around the closed orbit !

write(6,*) " closed orbit computation "
write(6,*) map12%x0(1)
c_w0%x0=map12%x0
c_w0%v(1)=map12%v(1)-(1.d0.cmono.1)-map12%x0(1)
c_w0_inv=c_w0.oo.(-1)

f_map%x0=0
f_map%v(1)=0.d0
f_map=c_w0_inv.o.f_map
f_ray=0
f_ray=c_w0_inv.o.f_ray
call print(c_w0_inv)
call print(f_map)
write(6,*) "exact ", x_closed
x_closed_tpsa=f_map%v(1)
write(6,*) "TPSA  ", x_closed_tpsa
write(6,*) "f_ray ", f_ray%x(1)
 

!x_closed_tpsa=x_closed
go_to_orbit%v(1)= (1.d0.cmono.1)+x_closed_tpsa
go_to_orbit%x0=x_closed_tpsa 

map12=(go_to_orbit.o.map12).o.(go_to_orbit.oo.(-1))

call print(map12)

 

x=map12%v(1)+i_*one_turn_map_AP%v(1)
 write(6,*)  "    Maps around the TPSA closed orbit "
 write(6,*)  "    Comparing TPSA one-turn map around the TPSA fixed point "
write(6,*)  "    with the DA one-turn map"
 write(6,*)  "    TPSA coefficients      DA coefficient "

 call print(x)



end program tpsa_1d_map

