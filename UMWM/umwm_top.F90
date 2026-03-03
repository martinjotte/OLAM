module umwm_top

contains

!===============================================================================

subroutine umwm_initialize (dtolam)

  use umwm_init,   only: nmlassign, init
  use umwm_stokes, only: stokes_init
  use umwm_module, only: dt_olam, alloc_umwm, filltab_umwm, stokes
  use misc_coms,   only: io6

  implicit none

  real, intent(in) :: dtolam

  dt_olam = dtolam

  write(io6,*)
  write(io6,*) 'University of Miami Wave Model version 2.0.0'
  write(io6,*)

  ! fill the namelist values by in-code assignment
  call nmlassign()

  ! allocate arrays
  call alloc_umwm()

  ! add UMWM array(s) to history I/O table
  call filltab_umwm()

  ! initialize model variables
  call init()

  ! initialize stokes drift arrays
  if (stokes) call stokes_init()

end subroutine umwm_initialize

!===============================================================================

subroutine umwm_hist_prep()

  use umwm_physics,          only: umwm_diag
  use umwm_oforcing,         only: oforcing
  use mem_sea,               only: msea, omsea
  use mem_sfcg,              only: itab_wsfc
  use mem_para,              only: myrank
  use umwm_stokes,           only: stokes_drift
  use umwm_module,           only: stokes, umwm

  implicit none

  integer :: i, iwsfc

  ! For restart or PLOTONLY runs, recompute some wave diagnostic fields that
  ! aren't saved in the HISTORY files.

  !$omp parallel do private(iwsfc)
  do i = 2, msea
     iwsfc = i + omsea

     ! skip cells that are not primary on this rank
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     call oforcing(i, iwsfc)

     ! skip cells where UMWM was not active
     if (.not. umwm%iactive(i)) cycle

     ! stokes drift u,v components at sea surface
     if (stokes) call stokes_drift(i)

     ! diagnose several wave properties
     call umwm_diag(i)

  enddo
  !$omp end parallel do

end subroutine umwm_hist_prep

!===============================================================================

subroutine umwm_step()

  use umwm_module !,           only: evs, evsf, stokes, om, pm, dta, dt_olam, fice_uth, umwm, swh, ssin, sdv
  use umwm_physics,          only: umwm_diag
  use umwm_advection,        only: propagation, refraction
  use umwm_oforcing,         only: oforcing
  use umwm_stokes,           only: stokes_drift
  use umwm_stress,           only: stress_atm
  use umwm_source_functions, only: sin_d12, sds_d12, snl_d12, s_ice, source
  use mem_sea,               only: sea, msea, omsea
  use mem_sfcg,              only: sfcg, itab_wsfc
  use misc_coms,             only: time_istp8, io6, iparallel
  use consts_coms,           only: r8, piu180
  use mem_para,              only: myrank
  use olam_mpi_sfc,          only: mpi_send_wsfc, mpi_recv_wsfc

  implicit none

  real, external :: walltime
  real           :: wtime_start, w1, w2
  integer        :: i, o, p, iwsfc
  real           :: dummy(pm,om), sds(pm,om), snl(pm,om), sdt(om), sice(om)

  ! Set timestep to value from OLAM

  dta = dt_olam

  ! if first time step, set time step to zero and integrate only the diagnostic part

  if (time_istp8 < 1.e-3_r8) dta = 0.

  ! Carry out one time step of UMWM

  wtime_start = walltime(0.)
  w1 = walltime(wtime_start)

  !$omp parallel do private(iwsfc,dummy,sds,snl,sdt,sice)
  do i = 2, msea
     iwsfc = i + omsea

     ! skip cells that are not primary on this rank
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! update atmospheric wind forcing from OLAM
     call oforcing(i,iwsfc)

     ! skip cells with sufficient seaice
     if (sea%seaicec(i) > fice_uth .or. sfcg%glatw(iwsfc) > 83.) cycle

     ! compute wind input source term Ssin
     call sin_d12(i,iwsfc)

     ! compute source dissipation term Sds
     call sds_d12(i,dummy,sds)

     ! compute non-linear source term Snl
     call snl_d12(i,iwsfc,snl,sdt,sds)

     ! compute sea ice attenuation term Sice
     call s_ice(i,sice)

     ! integrate source functions forward in time
     call source(i,dummy,sds,snl,sdt,sice)

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_wsfc(set="umwm", umwmvar=evs)
     call mpi_recv_wsfc(set="umwm", umwmvar=evs)
  endif

  !$omp parallel
  !$omp do private(iwsfc)
  do i = 2, msea
     iwsfc = i + omsea

     ! skip cells that are not primary on this rank
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! skip cells with sufficient seaice
     if (sea%seaicec(i) > fice_uth .or. sfcg%glatw(iwsfc) > 83.) cycle

     ! compute advection and integrate
     call propagation(i, iwsfc)

     ! compute refraction and integrate
     call refraction(i, iwsfc)

  enddo
  !$omp end do

  !$omp do private(iwsfc,o,p)
  do i = 2, msea
     iwsfc = i + omsea

     ! skip cells that are not primary on this rank
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     if (sea%seaicec(i) > fice_uth .or. sfcg%glatw(iwsfc) > 83.) then

        ! special for seaice
        umwm%iactive(i) = .false.
        evs     (:,:,i) = 0.

     else

        ! Copy updated wave variance spectrum to e array
        do o = 1, om
           do p = 1, pm
              evs(p,o,i) = evsf(p,o,i)
              if (evs(p,o,i) < 1.e-30) evs(p,o,i) = 0.
           enddo
        enddo

        umwm%iactive(i) = .true.

     endif

     ! stokes drift u,v components at sea surface
     if (stokes) call stokes_drift(i)

     ! compute wind stress and drag coefficient
     call stress_atm(i,iwsfc)

     ! diagnose other quantities (wave height, wave direction)
     call umwm_diag(i)

  enddo
  !$omp end do nowait
  !$omp end parallel

  w2 = walltime(wtime_start)
  write(io6,'(a,f0.3,a)') '     UMWM walltime = ', w2-w1, ' sec'

!!  i = 13710
!!
!!! write(*,*) 'OC: ', oc(i)
!!  write(*,*) 'Wave height, Wind speed and dir', swh(i), umwm%wspd(i), umwm%wdir(i)
!!
!!  write(*,'(2x,37g10.2)') (o,o=10,om)
!!
!!  do p = 1,pm
!!     write(*,'(a,i2,37g10.2)') 'evs ',p,(evs(p,o,i),o=10,om)
!!  enddo
!!
!!  write(*,*)
!!
!!  do p = 1,pm
!!     write(*,'(a,i2,37g10.2)') 'ssin',p,(ssin(p,o,i),o=10,om)
!!  enddo

end subroutine umwm_step

end module umwm_top
