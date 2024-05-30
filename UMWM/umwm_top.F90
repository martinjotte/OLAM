module umwm_top

  implicit none

contains

!===============================================================================

subroutine umwm_initialize (dtolam)

  use umwm_init,   only: nmlassign, init
  use umwm_stokes, only: stokes_drift
  use umwm_module, only: dt_olam, alloc_umwm, filltab_umwm

  real, intent(in) :: dtolam

  dt_olam = dtolam

  print '(a)', ' University of Miami Wave Model version 2.0.0'

  call nmlassign()          ! fill the namelist values by in-code assignment
  call alloc_umwm()         ! allocate arrays
  call filltab_umwm()       ! add UMWM array(s) to history I/O table
  call init()               ! initialize model variables
  call stokes_drift('init') ! initialize stokes drift arrays

end subroutine umwm_initialize

!===============================================================================

subroutine umwm_plot_prep()

  use umwm_physics,          only: umwm_diag
  use umwm_oforcing,         only: oforcing
  use umwm_stress,           only: stress
  use umwm_source_functions, only: sin_d12, sds_d12, snl_d12, s_ice

  ! Call a group of subroutines right before plotting fields at the start of
  ! INITIAL or after any HISTORY read to make available various quantities for
  ! plotting.  (The first 3 routines are also called on each timestep.)

  call oforcing()    ! update force fields from OLAM
  call sin_d12()     ! compute source input term ssin (uses ustar)
  call stress('atm') ! compute wind stress, drag coefficient, and ustar (uses ssin)
  call umwm_diag()   ! diagnose several wave properties

end subroutine umwm_plot_prep
 
!===============================================================================

subroutine umwm_step()

  use umwm_module
  use umwm_physics,          only: source, umwm_diag
  use umwm_advection,        only: propagation, refraction
  use umwm_oforcing,         only: oforcing
  use umwm_stokes,           only: stokes_drift
  use umwm_stress,           only: stress
  use umwm_source_functions, only: sin_d12, sds_d12, snl_d12, s_ice
  use mem_sea,               only: msea, omsea
  use mem_sfcg,              only: sfcg
  use misc_coms,             only: time_istp8
  use consts_coms,           only: r8, piu180

  real, external :: walltime
  real :: wtime_start,t1,w1,w2,t2
  integer :: i,o,p

  ! Carry out one time step of UMWM

     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call oforcing()    ! update force fields from OLAM

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime1 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call sin_d12()     ! compute source input term Sin

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime2 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call sds_d12()     ! compute source dissipation term Sds

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime3 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call snl_d12()     ! compute non-linear source term Snl

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime4 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call s_ice()       ! compute sea ice attenuation term Sice

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime5 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call source()      ! integrate source functions

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime6 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call propagation() ! compute advection and integrate

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime7 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  call refraction()  ! compute refraction and integrate

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime8 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  !$omp parallel do
  do i = 2,msea
     evs(:,:,i) = evsf(:,:,i) ! Copy updated wave variance sprctrum to e array
  enddo
  !$omp end parallel do

  call stress('atm') ! compute wind stress and drag coefficient

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime9 ',w2-w1
     wtime_start = walltime(0.)
     w1 = walltime(wtime_start)

  if (stokes) call stokes_drift()

     w2 = walltime(wtime_start)
     write(6,'(a,f13.3)') 'toptime10 ',w2-w1

  ! call stress('ocn') ! compute stress into ocean top and bottom (not needed for OLAM coupling)

  i = 91457 - omsea

  print*, ' '
  do p = 1,pm
     write(6,'(a,i2,37f7.1)') 'evs ',p,(evs(o,p,i),o=1,om)
  enddo
  print*, ' '
  do p = 1,pm
     write(6,'(a,i2,37f7.1)') 'evsk ',p,(evs(o,p,i)*kdk(o,i)*1.e3,o=1,om)
  enddo
  print*, ' '

end subroutine umwm_step

end module umwm_top
