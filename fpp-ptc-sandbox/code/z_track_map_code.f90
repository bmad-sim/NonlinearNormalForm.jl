program tpsa_da_map
  !use madx_ptc_module
  use pointer_lattice
  use gauss_dis
  
  implicit none
  type(layout), pointer:: ring
  real(dp) prec  
  real(dp) k1 
  type(c_ray) f_ray
  type(c_normal_form) normal_form
   integer  no,np,neq,NV
  integer i,j,k, case_section
  type(probe_8) xs  
  !!!!!!!!!!!!!!!!!!!!!
  type(fibre), pointer :: p 
  type(probe) ray_TC 
  type(probe_8) ray_8_TC
  type(c_damap)  one_turn_map_AP, identity_AP,two_turn_map_AP,lagrange_map_ap
  real(dp)  closed_orbit(6),goal(9),entrance_orbit(6)
  type(internal_state),target :: bmad_state
  type(c_taylor) ct,phase(3),eq(9)  
  type(c_universal_taylor) tpos(9)  
  
  c_verbose=.false.
  prec=1.d-10 ! for printing
  use_info=.true.
  use_quaternion=.true.
   
  call ptc_ini_no_append
  
   
   call append_empty_layout(m_u)
   
  ring=>m_u%start
  
   call build_lattice_als(ring,mis=.true.,exact=.true., thin=.false. ,onecell=.true.)
  
  
  call make_node_layout(ring)
  call in_bmad_units()
  
  bmad_state=default+time0
  write(6,*) " Select case 1,2,3 or 4"
  write(6,*) " 1 : comparing DA maps vs tracking code on closed orbit "
  write(6,*) " 2 : one-turn x(1) as a function of 6 parameters using PTC "
  write(6,*) " 3 : fitting x(1:4) using 6 parameters and storing in universal taylors "
  write(6,*) " 4 : fitting x(1:4), 2 tunes, time slip, using 12 parameters "
  read(5,*) case_section
  
  
  
  select case(case_section)
  case(1)
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
  call make_it_knob(p%magp%bn(1),1)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QF2")   !  Locating QF2 in TC
  call make_it_knob(p%magp%an(1),2)  !  PP : making a Polymorph into a knob
  
  p=>ring%start
  closed_orbit=0
  call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
  write(6,*) " closed orbit "
  write(6,"(6(1x,g11.4))") closed_orbit
   
  no=1 ; np=2 ;  ! nv=6+2=8
  call init(bmad_state,no,np)   ! Berz's TPSA package is initialized
  
  call alloc(one_turn_map_AP, identity_AP,two_turn_map_AP )
  call alloc(ray_8_TC)
  
  p=>ring%start  !  Moving the tracking point to the start of the ring
  
  ray_TC=closed_orbit   !  For TC
  
  identity_AP=1  !  For AP
  
  ray_8_TC=ray_TC+identity_AP   ! Connect AP with TC
  write(6,*) " The initial x "
  call print(ray_8_TC%x(1))
  call propagate(ray_8_TC,+bmad_state,fibre1=p) ! TC of PTC proper
  
  one_turn_map_AP=ray_8_TC ! TC into AP : makes a map out of ray_8_TC
  
  write(6,*) "x_final for one turn "
  call print(one_turn_map_AP%v(1))
  
  two_turn_map_AP=one_turn_map_AP*one_turn_map_AP
  write(6,*) "x_final for two turns : squaring the map "
  call print(two_turn_map_AP%v(1))
  
  !  Tracking two turns using TC
  ray_8_TC=ray_TC+identity_AP   ! Connect AP with TC
  call propagate(ray_8_TC,+bmad_state,fibre1=p) ! TC of PTC proper
  call propagate(ray_8_TC,+bmad_state,fibre1=p) ! TC of PTC proper
  two_turn_map_AP=ray_8_TC ! TC into AP : makes a map out of ray_8_TC
  write(6,*) "x_final for two turns : tracking two turns "
  call print(two_turn_map_AP%v(1))
  call KILL_PARA(ring)
  
  case(2)
  
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
  call make_it_knob(p%magp%bn(1),1)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QF2")   ! Locating QF2 in TC
  call make_it_knob(p%magp%an(1),2)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD1")     !  Locating QD1 in TC
  call make_it_knob(p%magp%bn(1),3)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD2")     !  Locating QD2 in TC
  call make_it_knob(p%magp%an(1),4)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND")     !  Locating BEND in TC
  call make_it_knob(p%magp%bn(1),5)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND1")    !  Locating BEND1 in TC
  call make_it_knob(p%magp%an(1),6)  !  PP : making a Polymorph into a knob
  
   
   
  bmad_state=bmad_state+nocavity0 ! turns cavities into drifts
  
  p=>ring%start
  closed_orbit=0
  call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
  write(6,*) " closed orbit "
  write(6,"(6(1x,g11.4))") closed_orbit
   
  np=6     ! number of free parameters to solve these equations
  no=1; nV=NP ;  ! nv=6 
  call c_init_all(no,nv)   ! FPP command to initialize the TPSA
    
  call alloc(ct);call alloc(ray_8_TC);
   
  ray_TC=entrance_orbit   !  For TC
  ray_8_TC=ray_TC  ! Connect AP with TC: no identity added
  write(6,*) "The variable x = z_1^0 "
  call print(ray_8_TC%x(1))
  p=>ring%start
  call propagate(ray_8_TC,+bmad_state,fibre1=p) !  TC of PTC proper
  call print(ray_8_TC%x(1))
  
  
   call KILL_PARA(ring)
  
  case(3)
  
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
  call make_it_knob(p%magp%bn(1),1)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QF2")   ! Locating QF2 in TC
  call make_it_knob(p%magp%an(1),2)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD1")     !  Locating QD1 in TC
  call make_it_knob(p%magp%bn(1),3)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD2")     !  Locating QD2 in TC
  call make_it_knob(p%magp%an(1),4)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND")     !  Locating BEND in TC
  call make_it_knob(p%magp%bn(1),5)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND1")    !  Locating BEND1 in TC
  call make_it_knob(p%magp%an(1),6)  !  PP : making a Polymorph into a knob
  
   
   
  bmad_state=bmad_state+nocavity0 ! turns cavities into drifts
  
  p=>ring%start
  closed_orbit=0
  call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
  write(6,*) " closed orbit "
  write(6,"(6(1x,g11.4))") closed_orbit
   
   goal=0.d0
  entrance_orbit=0
  goal(1:6)=entrance_orbit
  
  do k=1,4    ! performing 4 iterations
  
     !!!!!! Part 1 : calling TC !!!!!!
     neq=4    ! number of equations to solve
     np=6     ! number of free parameters to solve these equations
     no=1; nV=NP ;  ! nv=6 
     call c_init_all(no,nv)   ! FPP command to initialize the TPSA
      
     call alloc(ct);call alloc(ray_8_TC);
     
     ray_TC=entrance_orbit   !  For TC
     ray_8_TC=ray_TC  ! Connect AP with TC: no identity added
  
     p=>ring%start
     call propagate(ray_8_TC,+bmad_state,fibre1=p) !  TC of PTC proper
     
     !!!!!! Part 2 : saving the relevant Taylor series  !!!!!!
  
     call ALLOC(tpos,1,NV,0) ! c_universal_taylor
     do i=1,4
      ct=ray_8_TC%x(i)%t   ! saving final transverse positions and momenta
      tpos(i)=ct
     enddo
     call kill(ct);call kill(ray_8_TC);
  
  !!!!!!!! Part 3 : creating a map of some sort and using it !!!!!!
  
  
     nV=NP+neq ;  ! nv=6+4=10
     call c_init_all(no,nv)   ! FPP command to initialize the TPSA
     call alloc(ct);call alloc(eq);
  
     lagrange_map_ap%n=nv
     call alloc(lagrange_map_ap)
     
     do i=1,4
      eq(i)=tpos(i)   ! putting the c_universal_taylor into c_taylor
     enddo
     
     !!!  Construction of Lagrange function
     do i=1,np
      lagrange_map_ap%v(i)=1.d0.cmono.i 
      do j=np+1,np+neq
  
       ct=((eq(j-np).d.i)*(1.d0.cmono.j))
       lagrange_map_ap%v(i)=lagrange_map_ap%v(i)+ct
      enddo         
     enddo 
  
     do i=np+1,np+neq
      lagrange_map_ap%v(i)=eq(i-np)-goal(i-np)       
     enddo 
     
     lagrange_map_ap=lagrange_map_ap.oo.(-1)
  
     f_ray=0
     f_ray=lagrange_map_ap.o.f_ray
     
     p=>ring%start
     call move_to(ring,p,"QF1")   !  Locating QF1 in TC
      k1=f_ray%x(1)   
      call add(p,1,1,k1)
     call move_to(ring,p,"QF2") 
      k1=f_ray%x(2)   
      call add(p,-1,1,k1)
     call move_to(ring,p,"QD1")   !  Locating QF1 in TC
      k1=f_ray%x(3)   
      call add(p,1,1,k1)
     call move_to(ring,p,"QD2")   !  Locating QF2 in TC
      k1=f_ray%x(4)   
      call add(p,-1,1,k1)
     call move_to(ring,p,"BEND")   !  Locating QF2 in TC
      k1=f_ray%x(5)  
      call add(p,1,1,k1)
     call move_to(ring,p,"BEND1")   !  Locating QF2 in TC
      k1=f_ray%x(6)   
      call add(p,-1,1,k1)
  
     p=>ring%start
  
     ray_TC=entrance_orbit   !  For TC
      
     call propagate(ray_TC,bmad_state,fibre1=p) !  TC of PTC proper
  
     write(6,*) " exit orbit " 
     write(6,"(4(1x,g11.4))") ray_TC%x(1:4)
     call kill(lagrange_map_ap);call kill(ct);call kill(eq);call kill(tpos);
  
  enddo
  
   call KILL_PARA(ring)
  case(4)
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
  call make_it_knob(p%magp%bn(2),1)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QF2")   ! Locating QF2 in TC
  call make_it_knob(p%magp%an(2),2)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD1")     !  Locating QD1 in TC
  call make_it_knob(p%magp%bn(2),3)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD2")     !  Locating QD2 in TC
  call make_it_knob(p%magp%an(2),4)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND")     !  Locating BEND in TC
  call make_it_knob(p%magp%bn(2),5)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND1")    !  Locating BEND1 in TC
  call make_it_knob(p%magp%an(2),6)  !  PP : making a Polymorph into a knob
  
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
  call make_it_knob(p%magp%bn(1),7)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QF2")   ! Locating QF2 in TC
  call make_it_knob(p%magp%an(1),8)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD1")     !  Locating QD1 in TC
  call make_it_knob(p%magp%bn(1),9)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"QD2")     !  Locating QD2 in TC
  call make_it_knob(p%magp%an(1),10)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND")     !  Locating BEND in TC
  call make_it_knob(p%magp%bn(1),11)  !  PP : making a Polymorph into a knob
  call move_to(ring,p,"BEND1")    !  Locating BEND1 in TC
  call make_it_knob(p%magp%an(1),12)  !  PP : making a Polymorph into a knob
  
  
  write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
  
  bmad_state=bmad_state+nocavity0
  
  p=>ring%start
  closed_orbit=0
  call find_orbit_x(closed_orbit,bmad_state, 1.0e-7_dp, fibre1=p)
  write(6,*) " closed orbit "
  write(6,"(6(1x,g11.4))") closed_orbit
   
  goal=0.d0
  entrance_orbit=closed_orbit
  goal(1)=0.28d0
  goal(2)=0.79d0
  goal(3)=(-1)**ndpt_bmad*2.5d-2
  goal(4:7)=entrance_orbit(1:4)
  
  do k=1,8
  
  !!!!!! Part 1 : calling TC !!!!!!
  neq=7     ! number of equations to solve
  np=12     ! number of free parameters to solve these equations
  no=2 ;    ! nv=6+12=18
  call init(bmad_state,no,np)   ! Berz's TPSA package is initialized
  nV=c_%nv   !     
  call alloc(ct)  
  call alloc(identity_AP,one_turn_map_AP);call alloc(normal_form);
  call alloc(ray_8_TC);call alloc(phase);call alloc(eq)
  p=>ring%start
     
  ray_TC=entrance_orbit   !  For TC
  identity_AP=1
  ray_8_TC=ray_TC + identity_AP ! Connect AP with TC: identity added
   
  p=>ring%start
  call propagate(ray_8_TC,+bmad_state,fibre1=p) !  TC of PTC proper
   
  one_turn_map_AP=ray_8_TC
  
  ! compute 2 tunes and the phase slip
  call c_normal(one_turn_map_AP,normal_form,phase=phase) 
  ray_TC=ray_8_TC
  
  write(6,*) 
  write(6,"(a27,3(1x,g17.10,1x))") "2 tunes and the phase slip",normal_form%tune(1:3)
  write(6,"(a13,4(1x,g17.10,1x))") "Closed orbit",ray_TC%x(1:4)
   identity_AP=0
  
  phase(3)=phase(3).d.(5+ndpt_bmad)    ! compute phase slip dbeta*t/ddelta 
   
  phase(1)=(phase(1).o.identity_AP)<=6  ! removing delta dependence
  phase(2)=(phase(2).o.identity_AP)<=6  ! removing delta dependence
  phase(3)=(phase(3).o.identity_AP)<=6  ! removing delta dependence
   
  
  eq(1)=phase(1)-goal(1)  ! tune x
  eq(2)=phase(2)-goal(2)  ! tune y
  eq(3)=phase(3)-goal(3)  ! phase slip  
   
  do i=1,4
  ct=ray_8_TC%x(i)%t
  eq(3+i)=((ct.o.identity_AP)<=6)-goal(3+i)
  enddo
  
  !!!!!! Part 2 : saving the relevant Taylor series  !!!!!!
     call ALLOC(tpos,1,NV,0)
  do i=1,7
   tpos(i)=eq(i)
  enddo
  
  
  call kill(identity_AP,one_turn_map_AP);call kill(normal_form);
  call kill(ray_8_TC);call kill(phase);call kill(eq);call kill(ct)  ;
  !!!!!! Part 3 : creating a map of some sort and using it !!!!!!
  
  neq=7    ! number of equations to solve
  np=12     ! number of free parameters to solve these equations
  no=no-1; nV=NP+neq ;  ! nv=7+12=19
  call c_init_all(no,nv)   ! FPP command to initialize the TPSA
  lagrange_map_ap%n=nV   ! size of the c_damap must be explicitely stated
  call alloc(ct);call alloc(eq);call alloc(lagrange_map_ap);
  do i=1,7
   eq(i)=tpos(i)
  enddo
  !!!  Construction of Lagrange function
  do i=1,np
   lagrange_map_ap%v(i)=1.d0.cmono.i 
   do j=np+1,np+neq 
    ct=(eq(j-np).d.i)*(1.d0.cmono.j)
     lagrange_map_ap%v(i)=lagrange_map_ap%v(i)+ct
   enddo         
  enddo 
   
  do i=np+1,np+neq 
   lagrange_map_ap%v(i)=eq(i-np) 
  enddo 
  
  lagrange_map_ap=lagrange_map_ap.oo.(-1)
  
  f_ray=0
  f_ray=lagrange_map_ap.o.f_ray
  
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
   k1=f_ray%x(1)   
   call add(p,2,1,k1)
  call move_to(ring,p,"QF2") 
   k1=f_ray%x(2)   
   call add(p,-2,1,k1)
  call move_to(ring,p,"QD1")   !  Locating QF1 in TC
   k1=f_ray%x(3)   
   call add(p,2,1,k1)
  call move_to(ring,p,"QD2")   !  Locating QF2 in TC
   k1=f_ray%x(4)   
   call add(p,-2,1,k1)
  call move_to(ring,p,"BEND")   !  Locating QF2 in TC
   k1=f_ray%x(5)   
   call add(p,2,1,k1)
  call move_to(ring,p,"BEND1")   !  Locating QF2 in TC
   k1=f_ray%x(6)   
   call add(p,-2,1,k1)
  
  p=>ring%start
  call move_to(ring,p,"QF1")   !  Locating QF1 in TC
   k1=f_ray%x(7)   
   call add(p,1,1,k1)
  call move_to(ring,p,"QF2")   ! Locating QF2 in TC
   k1=f_ray%x(8)   
   call add(p,-1,1,k1)
  call move_to(ring,p,"QD1")     !  Locating QD1 in TC
   k1=f_ray%x(9)   
   call add(p,1,1,k1)
  call move_to(ring,p,"QD2")     !  Locating QD2 in TC
   k1=f_ray%x(10)   
   call add(p,-1,1,k1)
  call move_to(ring,p,"BEND")     !  Locating BEND in TC
   k1=f_ray%x(11)   
   call add(p,1,1,k1)
  call move_to(ring,p,"BEND1")    !  Locating BEND1 in TC
   k1=f_ray%x(12)   
   call add(p,-1,1,k1)
  
  call kill(ct) ;call kill(eq);call kill(lagrange_map_ap);
  enddo
   
  case default
   stop 666
  end select
  
  call ptc_end(graphics_maybe=0,flat_file=.false.)
  contains
  
  
  
  subroutine  build_lattice_als(ALS,mis,error,exact,sl,thin,onecell)
  use madx_ptc_module
  use pointer_lattice
  implicit none
  
  type(layout), target :: ALS
  real(dp),optional :: error(6)
  logical, optional :: exact,sl,thin,onecell
  real(dp) :: alpha,lbend, cut, ksd, ksf 
  type(fibre)  L1,L2,L3,L4,L5,L6,L7,L8,L9,L10 
  type(fibre)  L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,CAVM
  type(fibre)  L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
  type(fibre)  QF1,QF2,QD1,QD2,QFA1,QFA2,sf,sd,cav,bend,vc5,bend1 
  type(layout) :: sfline,sdline,sup1,supb
  logical(lp) :: mis,thi=.false.,oneperiod
  real(dp) sig(6)
  !-----------------------------------
  if(present(thin)) thi=thin
  
  
  call make_states(.true.)
  exact_model = .false.;oneperiod = .false.
  if(present(exact)) exact_model=exact
  if(present(onecell)) oneperiod=onecell
  call update_states
  madlength = .false.
  
  
  call set_mad(energy = 1.5d0, method = 2, step = 8)
  
  !madkind2 = matrix_kick_matrix
  madkind2 = drift_kick_drift
  
  
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
  
  QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+6.447435260914397D-03)
  QF2 = QUADRUPOLE(" QF2 ",0.344D0, K1= 2.2474D0)
  QD1 = QUADRUPOLE(" QD1 ",0.187D0, K1= -2.3368D0-2.593018157427161D-02); 
  QD2 = QUADRUPOLE(" QD2 ",0.187D0, K1= -2.3368D0);  
  QFA1= QUADRUPOLE(" QFA1",0.448D0, K1= 2.8856D0);  
  QFA2= QUADRUPOLE(" QFA2",0.448D0, K1= 2.8856D0);  
  
  !!! 1/2 mad-x value
  ksf=-41.3355516397069748d0;
  ksd=56.2564709584745489d0;
  
  sf=sextupole ("sf",2.d0*0.1015d0, K2= ksf);
  sd= sextupole("sd", 2.d0*0.1015d0, K2= ksd);
  
   VC5=marker("vc5");
  ALPHA=0.17453292519943295769236907684886d0;
   
  LBEND=0.86621d0;
   
  BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
  BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
   
  CAVM=MARK("CAVM");
  N_CAV4_F=3
  CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1000.0d0,REV_FREQ=500.0d6)
  !CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1.0d0,REV_FREQ=500.0d6)
  write(6,*) cav%mag%c4%nf
  write(6,*) cav%mag%c4%f
  cav%mag%c4%f(1)=1
  cav%magp%c4%f(1)=cav%mag%c4%f(1)
  !cav%mag%c4%f(3)=-0.3d0*cav%mag%c4%f(1)
  !cav%magp%c4%f(3)=cav%mag%c4%f(3)
  
  if(thi) then
   sf=sextupole ("sf",0.d0, K2= ksf*0.203d0);
   sd= sextupole("sd", 0.d0, K2= ksd*0.203d0);
    sfline=(ds+sf+ds);
    sdline=(ds+sd+ds);
  else
   sfline=1*sf;
   sdline=1*sd;
  endif
  
  SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
             L11+QFA1+L12+sdline+L13+ &
             L14+BEND+L15+L16+sdline+L17+ &
             QFA2+L18+L19+sfline+L20+BEND+L21+&
             L22+QD2+L23+L24+QF2+L25+ &
             L26+VC5+L27+cavm;
  
  SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
             L11+QFA1+L12+sdline+L13+ &
             L14+BEND+L15+L16+sdline+L17+ &
             QFA2+L18+L19+sfline+L20+BEND1+L21+&
             L22+QD2+L23+L24+QF2+L25+ &
             L26+VC5+L27+cav;
  
  if(oneperiod) then
   ALS = supb;  !11*sup1+supb;
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
  
  call gaussian_seed(1000)
  if(mis) then
   sig=1.d-5; cut=4.d0; 
   if(present(error)) sig=error
   call MESS_UP_ALIGNMENT(ALS,SIG,cut);
  endif
  end subroutine build_lattice_als
  
  
  
  end program tpsa_da_map
  
  