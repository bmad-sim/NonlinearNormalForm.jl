program Guignard_normal_form_average_x
   use madx_ptc_module
   use pointer_lattice
   use c_TPSA
   implicit none
   
   
   type(layout), pointer:: ALS
   real(dp) prec
   real(dp) circ,dtheta,s,total_tunes(4),ts,xx,ph(3),damp(3),misa(6)
   type(internal_state),target :: state 
   logical(lp) :: mis=.false. 
   type(c_taylor) fonction,fonction_FLOQUET,phase(3)
   type(c_damap)  one_turn_map,id_s,U_c,U,A ,a_cs ,N_c 
   type(c_normal_form) normal_form,normal_form_no_fac
   type(c_vector_field) H_theta,f_non,H0_hat ,f_com
   integer :: pos =1
   integer i,map_order,mf,mf1,n_s,n_mode,km ,n
   type(probe) ray_closed
   type(probe_8) ray,ray_cs_twiss
   type(fibre), pointer :: p
   type(integration_node), pointer :: t
   character*120 :: dsc
   logical :: doit=.true. ,skew=.false.
   type(c_vector_field_fourier) H,F,K,f1
   real(dp) DX_AVERAGE_DCS,betax_1,theta
   complex(dp) coe
   !!!!!!!!!!!!!!!!!!!!!
   integer,target :: mybeg, M_res(-1:ndim)
   real(dp) , target :: mondelta,closed_orbit(6)  
   integer selcas
   
   use_quaternion=.true.
   c_verbose=.false.
   prec=1.d-6 ! for printing 
   use_info = .true.
   longprint=.false.
   c_lda_used=c_lda_used*30
   n_cai=-i_
   mis=.false.
   !state=nocavity0 +time0+spin0
   state=nocavity0 +time0+spin0
   !state=default0 +spin0  +radiation0
   state=only_4d0  +spin0
   skew=.false.
   !remove_tune_shift=.false.
   radfac=100  
   call in_bmad_units
   
   call ptc_ini_no_append
   call append_empty_layout(m_u)
   ALS=>m_u%start
   
   call build_lattice(ALS,mis,exact=.false.,thin=.false.,onecell=.false.) 
   
   
   
    !call read_ptc_command77("fit_tune.txt")  
   p=>als%start
   misa=0
   call move_to(als,p,"BEND")
   if(state%nocavity) then
   misa(4)=0.03d0
   misa(6)=0.03d0
   else
   misa(4)=0.003d0
   misa(6)=0.003d0
   endif
   misa=0.1d0*misa
   call MISALIGN_FIBRE(p,misa)
   
   courant_snyder_teng_edwards=.true.
   time_lie_choice=.true.
   
   p=>als%start
   
    call MAKE_NODE_LAYOUT(als)
    
   my_ering=>als
   my_estate=>state
   my_start=>mybeg
   mybeg=1
   mondelta=0
   closed_orbit=0
   MY_DELTA=>mondelta
   my_fix=>closed_orbit
   
   !call phase_advance_n(6)
   !write(6,format3) als%end%tm%lf%k(1,2,2)
   write(6,format3) als%end%t2%lf%phase
   write(6,format3) als%end%t2%lf%damping
   write(6,format3) als%end%t2%lf%spin
   write(6,format1) als%end%t2%lf%phase(1)*5
   
   !goto 1001
   call kanalnummer(mf,"C:\document\my_tex_papers\fpp_handbook\julia\lat.txt")
   
    
   
    
     map_order=3   
     skew=.true. 
    
   
   call init_all(state,map_order,0)
   
   call alloc(one_turn_map, id_s,U_c,A,U,a_cs,N_c ) 
   call alloc(normal_form); call alloc(ray); call alloc(normal_form_no_fac)
   call alloc(f_non);  call alloc(f_com)
   call alloc(H_theta); call alloc(H0_hat)
   call alloc(fonction,fonction_FLOQUET)
   call alloc(H);call alloc(F);call alloc(K);
   call alloc(phase);call alloc(ray_cs_twiss)
    
   DX_AVERAGE_DCS=0.d0
   closed_orbit=0.d0;                                            
   call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=1)   
    
   ray_closed=closed_orbit     
   id_s=1;   
   ray=id_s+ray_closed;          
                             
   call propagate(als,RAY,state,fibre1=1)   ! (4)
   
   one_turn_map=ray                         ! (5)      
    
    
   one_turn_map%x0=0
    
   call kanalnummer(i,"spinmap.txt")
   call print(one_turn_map,i)
   close(I)
   !stop
   call c_normal(one_turn_map,normal_form,dospin=state%spin,phase=phase)  ! (6)
   
   id_s=one_turn_map
   
   
   one_turn_map=ci_phasor()*normal_form%atot**(-1)*id_s*normal_form%atot*c_phasor()
   call clean(one_turn_map, one_turn_map, prec=1.e-7_dp)
   one_turn_map%e_ij = 0
   !call c_q0_to_qr(one_turn_map%q, one_turn_map%q) ! only necessary when leaving spin resonance in so we make i, k parts 0
   call print(one_turn_map)
write(*,*) "hi"
   stop
   
   
   id_s=ci_phasor()*one_turn_map*c_phasor()
   
   one_turn_map=id_s
   
   one_turn_map=(one_turn_map.cut.(-2))**(-1)*one_turn_map
   call  print(one_turn_map)
   
   pause 444
   H_theta=ln(one_turn_map)
   
   call clean(H_theta,H_theta,prec=1.d-9)
    call print(H_theta)
   goto 1001
   stop
   
   A=normal_form%atot 
    
   !id=1
   !id%
   !stop
   
    
   call  c_fast_canonise(normal_form%A_t ,U_c)       ! (9a)
   !call  c_full_canonise(normal_form%A_t ,U_c)       ! (9a)
    
    
         
   dtheta=twopi/n_s          ! (11)
                    ! 
    
    
   ray=U_c+ray_closed;        ! (12)
   
   p=>als%start               ! (13a)
   t=>p%t1                    ! (13b)
   
   do i=1,als%t%n
   
    if(mod(i,n_s/10)==0) then
      write(6,*) als%t%n-i, " steps remaining "
    endif
    
    call propagate(als,ray,state,node1=i,node2=i+1)   ! (14a)
    
     ray_closed=ray ! Saving orbit  ! (15a)
    
   U=ray !copying the tracking result {M}_{i i+1}o{A}_{i}} into a map U (15)
   
   call  c_fast_canonise(U ,U_c)       ! (9a)
   
   
   
    ray=ray_closed+U_c                     ! (21)
   
    t=>t%next
    p=>t%parent_fibre
   enddo
   
    close(mf)
   
   stop 1111
   
   1001 call ptc_end(graphics_maybe=1,flat_file=.false.)
   
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
   type(fibre) L27h
   type(layout) :: sfline,sdline,sup1,supb
   logical(lp) :: mis,thi=.false.,oneperiod
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
   if(.not.skew)ksd=6.447435260914397D-02-2.5d-2
   !if(skew) ksd=6.447435260914397D-02-2.1d-2
   QF1 = QUADRUPOLE(" QF1 ",0.344D0, K1= 2.2474D0+ksd)
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
   if(skew) call add(sf,-3,1,300.d0)
    VC5=marker("vc5");
   ALPHA=0.17453292519943295769236907684886d0;
    
   LBEND=0.86621d0;
    
   BEND = RBEND("BEND", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
   BEND1 = RBEND("BEND1", LBEND, ANGLE=ALPHA).q.(-0.778741d0)
    
   CAVM=MARK("CAVM");
   !CAV=RFCAVITY("CAV",L=0.1000d0,VOLT=-100.0d0,REV_FREQ=500.0d6)
   CAV=RFCAVITY("CAV",L=0.0000d0,VOLT=-1000.0d0,REV_FREQ=500.0d6)
   
   if(thi) then
    sf=sextupole ("sf",0.d0, K2= ksf*0.203d0);
    sd= sextupole("sd", 0.d0, K2= ksd*0.203d0);
     sfline=(ds+sf+ds);
     sdline=(ds+sd+ds);
   else
    sfline=1*sf;
    sdline=1*sd;
   endif
   
   if(oneperiod) then
   SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
              L11+QFA1+L12+sdline+L13+ &
              L14+BEND+L15+L16+sdline+L17+ &
              QFA2+L18+L19+sfline+L20+BEND+L21+&
              L22+QD2+L23+L24+QF2+L25+ &
              L26+VC5+L27+cav;
   else
   SUP1=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
              L11+QFA1+L12+sdline+L13+ &
              L14+BEND+L15+L16+sdline+L17+ &
              QFA2+L18+L19+sfline+L20+BEND+L21+&
              L22+QD2+L23+L24+QF2+L25+ &
              L26+VC5+L27+cavm;
   endif
   !           L26+VC5+L27+cavm;
   
   SUPb=L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
              L11+QFA1+L12+sdline+L13+ &
              L14+BEND+L15+L16+sdline+L17+ &
              QFA2+L18+L19+sfline+L20+BEND1+L21+&
              L22+QD2+L23+L24+QF2+L25+ &
              L26+VC5+L27+cav;
   
   if(oneperiod) then
    ALS = sup1;  !11*sup1+supb;
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
   
   subroutine phase_advance_n(mf)
   implicit none
   type(fibre),pointer:: f
   type(probe) xs0,xs1
   type(probe_8) xs
   type(c_damap) id
   type(c_normal_form) n
   type(integration_node), pointer :: t
   real(dp) phase(3), spin_tune(2),damping(3)
   integer i,mff
   integer,optional :: mf
   type(c_taylor) nu_spin
   use_quaternion=.true.
   mff=0
   if(present(mf)) mff=mf
   
   f=>my_ering%start
   do i=1,my_start-1
   f=>f%next
   enddo
    
   write(6,*) f%mag%name
   my_fix=0
   
   my_fix(ndpt_bmad+5)=my_delta
   
   call find_orbit_x(my_fix,my_estate,1.e-8_dp,fibre1=f) 
   
   if(.not.check_stable) then
     write(6,*) "Could not find closed orbit "
     write(6,*) " No calculation done "
    return
   endif
   
   call init(my_estate,1,0)
   
   call alloc(id)
   call alloc(xs)
   call alloc(n)
   call alloc(nu_spin)
   xs0=my_fix
   id=1
   xs=xs0+id
    call propagate(xs,my_estate,fibre1=f)
   id=xs
   
   call c_normal(id,n,dospin=my_estate%spin,nu_spin=nu_spin)
   
   if(mff/=0) then
    write(mff,*) " Linear A from c_normal "
    call print(n%atot,mff)
   endif
   
   
   
   
   call c_fast_canonise(n%atot,n%atot, dospin=my_estate%spin)
   
   if(mff/=0) then
    write(mff,*) " Linear A canonised "
    call print(n%atot,mff)
    write(mff,*) " end of Info from phase_advance "
   endif
   
   !if(mff/=0) then
    write(mff,*) " Info from map :tunes, damping, spin,quaternion_angle/pi"
    write(mff,*) n%tune(1:c_%nd)
    write(mff,*) n%damping(1:c_%nd)
    write(mff,*) n%spin_tune,n%quaternion_angle/pi
    write(mff,*) " end of Info from map "
   !endif
   phase=0
   spin_tune=0
   damping=0
   t=>f%t1
   
   xs=xs0+n%atot
   
   do i=1,my_ering%t%n
   
    
    call propagate(xs,my_estate,node1=t,node2=t%next)
   xs0=xs
   n%atot=xs
   
    f%tm%lf%symplectic=.not.my_estate%radiation
    
     call c_fast_canonise(n%atot,n%atot,phase=phase,damping=damping,spin_tune=spin_tune,dospin=my_estate%spin)
    
   xs=xs0+n%atot
   t=>t%next
   
     call compute_lattice_functions(n%atot,t%lf)
    t%lf%phase=phase
    t%lf%damping=damping
    t%lf%spin=spin_tune
    t%lf%fix=xs0%x
    
   if(i==my_ering%t%n) then
    t%lf%phase=0
    t%lf%damping=0
    t%lf%spin=0
   endif
   
   
   enddo
   
   write(6,*)  "   "
   write(6,*)  " Phase advance and fractional"
   write(6,format3) phase(1:c_%nd)
   write(6,format3) n%tune(1:c_%nd)
    
   write(6,*)  " damping advance "
   write(6,format3) damping
   write(6,*)  " spin advance and chromaticity "
   write(6,format2) spin_tune
   write(6,format1) n%spin_tune
    write(6,*)  " Closed orbit before and after "
    write(6,format6)  my_fix
    write(6,format6)  xs0%x
   write(6,*)  "   "
   
   write(6,*) "spin from normal form"
   call print(nu_spin)
   
   call kill(id)
   call kill(xs)
   call kill(n)
   call kill(nu_spin)
   
   end subroutine phase_advance_n
   
   
   
   end program Guignard_normal_form_average_x
   
    