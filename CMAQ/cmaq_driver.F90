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

  integer :: mrl

#ifdef CMAQ_TIMING
  real, external :: walltime
  real :: wtime_start,t1,w1,w2,t2
#endif

  mrl = mrl_endl(istp)

  ! Photolysis rates

  if (istp == 1) call phot( )

  if (mrl > 0) then

!  write(*,'(A,10g13.5)') " ", cgrid(6,535,91), cgrid(6,535,n_gc_spc+60), cgrid(6,535,n_gc_spc+59), cgrid(6,535,n_gc_spc+62)

  ! Cloud/aqueous processes

#ifdef CMAQ_TIMING
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

!     write(*,'(A,10g13.5)') "a", cgrid(6,535,91), cgrid(6,535,n_gc_spc+60), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+62), cgrid(6,535,num_str:num_end)

     call grell_aq_driver( mrl )

!     write(*,'(A,10g13.5)') "b", cgrid(6,535,91), cgrid(6,535,n_gc_spc+60), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

     call rescld         ( mrl )

!     write(*,'(A,10g13.5)') "c", cgrid(6,535,91), cgrid(6,535,n_gc_spc+60), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aqch [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species to densities

     call rev_cgrid( mrl )

!     write(*,'(A,10g13.5)') "d", cgrid(6,535,91), cgrid(6,535,n_gc_spc+48), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

  ! Aerosol sedimentation

#ifdef CMAQ_TIMING
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)
     call cpu_time(t1)
#endif

     call aero_sedi( mrl )

!     write(*,'(A,10g13.5)') "e", cgrid(6,535,91), cgrid(6,535,n_gc_spc+48), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

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

!     write(*,'(A,10g13.5)') "f", cgrid(6,535,91), cgrid(6,535,n_gc_spc+48), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

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

!     write(*,'(A,10g13.5)') "g", cgrid(6,535,91), cgrid(6,535,n_gc_spc+48), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

#ifdef CMAQ_TIMING
     w2 = walltime(wtime_start)
     call cpu_time(t2)
     write(io6,'(a,2f13.3)') '++++++++++ Aero [cpu,wall(sec)]: ',t2-t1,w2-w1
#endif

  ! Convert aerosol species back to concentrations

     call conv_cgrid ( mrl )
  endif

!  write(*,'(A,10g13.5)') "h", cgrid(6,535,91), cgrid(6,535,n_gc_spc+48), cgrid(6,535,n_gc_spc+47), cgrid(6,535,n_gc_spc+50), cgrid(6,535,num_str:num_end)

end subroutine cmaq_driver
