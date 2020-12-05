subroutine cmaq_driver( )

  use cgrid_conv, only: conv_cgrid, rev_cgrid
  use sedv_defn,  only: aero_sedi
  use mem_ijtabs, only: istp, mrl_endl

#ifdef CMAQ_TIMING
  use misc_coms,  only: io6
#endif

  implicit none

  integer :: mrl

#ifdef CMAQ_TIMING
  real, external :: walltime
  real :: wtime_start,t1,w1,w2,t2
#endif

  mrl = mrl_endl(istp)

  ! Photolysis rates

  if (istp == 1) call phot( )

  if (mrl > 0) then

  ! Cloud/aqueous processes

#ifdef CMAQ_TIMING
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

     call rescld( mrl )

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aqch [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species to densities

     call rev_cgrid( mrl )

  ! Aerosol sedimentation

#ifdef CMAQ_TIMING
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

     call aero_sedi( mrl )

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Sedi [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Gas chemistry

#ifdef CMAQ_TIMING
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

     call chem( mrl )

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Chem [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Aerosol microphysics

#ifdef CMAQ_TIMING
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

     call aero( mrl )

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aero [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species back to concentrations

     call conv_cgrid ( mrl )
  endif

end subroutine cmaq_driver
