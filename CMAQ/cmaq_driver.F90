subroutine cmaq_driver( )

  use cgrid_conv,      only: conv_cgrid, rev_cgrid
  use sedv_defn,       only: aero_sedi
  use mem_ijtabs,      only: istp, mrl_endl
  use module_cu_g3_aq, only: grell_aq_driver

  use cgrid_defn
  use cgrid_spcs
  use aero_data

#ifdef CMAQ_TIMING
  use misc_coms,  only: io6
#endif

  implicit none

#ifdef CMAQ_TIMING
  real :: t1,w1,w2,t2
#endif

  ! Photolysis rates

  if (istp == 1) call phot( )

  if (mrl_endl(istp) > 0) then

  ! Cloud/aqueous processes

#ifdef CMAQ_TIMING
     call walltime(w1)
     call cpu_time(t1)
#endif

     call grell_aq_driver( )

     call rescld         ( )

#ifdef CMAQ_TIMING
     call walltime(w2)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aqch [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species to densities

     call rev_cgrid( )

  ! Aerosol sedimentation

#ifdef CMAQ_TIMING
     call walltime(w1)
     call cpu_time(t1)
#endif

     call aero_sedi( )

#ifdef CMAQ_TIMING
     call walltime(w2)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Sedi [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Gas chemistry

     call chem( )

#ifdef CMAQ_TIMING
     call walltime(w1)
     call cpu_time(t1)
     write(io6,'(a,2f13.3)') '++++++++++ Chem [cpu,wall(sec)]: ',t1-t2,w1-w2
#endif

  ! Aerosol microphysics

     call aero( )

#ifdef CMAQ_TIMING
     call walltime(w2)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aero [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species back to concentrations

     call conv_cgrid ( )
  endif

end subroutine cmaq_driver
