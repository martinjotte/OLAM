module obs_nudge_mod

  use consts_coms, only: r8
  implicit none

  real(r8)             :: past_mass
  real(r8)             :: future_mass
  real(r8)             :: current_mass
  logical, allocatable :: skipw(:)
  integer              :: maxmrl = -1

  private :: r8

contains

!===============================================================================

subroutine nudge_prep_obs( iaction )

  use mem_nudge,   only: rho_obsp, theta_obsp, rrw_obsp, uzonal_obsp, umerid_obsp, &
                         rho_obsf, theta_obsf, rrw_obsf, uzonal_obsf, umerid_obsf, &
                         nudflag, nudnxp
  use consts_coms, only: r8
  use isan_coms,   only: o_rho, o_theta, o_rrw, o_uzonal, o_umerid
  use mem_basic,   only: rho
  use check_nan,   only: mass_sum_from_rho
  use misc_coms,   only: io6, iparallel
  use mem_grid,    only: mwa
  use mem_ijtabs,  only: jtab_w, itab_w, jtw_prog
  use oname_coms,  only: nl

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iaction
  integer             :: j, iw, ier

  ! Make sure we are doing obs/gridpoint nudging

  if (nudflag < 1) return
  if (nudflag > 0 .and. nudnxp > 0) return

  ! On first call, save the max MRL and store points that we will skip,
  ! If we are not skipping any points, compute total (dry) mass in atmosphere.

  if (iaction == 0 .and. nl%nud_preserve_total_mass) then
     maxmrl = maxval( itab_w(2:mwa)%mrlw )

#ifdef OLAM_MPI
     if (iparallel == 1) then
        call MPI_Allreduce( MPI_IN_PLACE, maxmrl, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier )
     endif
#endif

     if (maxmrl > nl%max_nud_mrl) then
        allocate( skipw(mwa) )

        !$omp parallel do private(iw)
        do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
           skipw(iw) = (itab_w(iw)%mrlw > nl%max_nud_mrl)
        enddo
        !$omp end parallel do

     else

        call mass_sum_from_rho( current_mass, rho, allnodes=.true. )

     endif

  endif

  ! Swap future data time into past data time if necessary

  if (iaction == 1) then
     call move_alloc(    rho_obsf,    rho_obsp )
     call move_alloc(  theta_obsf,  theta_obsp )
     call move_alloc(    rrw_obsf,    rrw_obsp )
     call move_alloc( uzonal_obsf, uzonal_obsp )
     call move_alloc( umerid_obsf, umerid_obsp )

     past_mass = future_mass
  endif

  ! Swap gridded data fields into future nudging arrays

  call move_alloc( o_rho   ,    rho_obsf )
  call move_alloc( o_theta ,  theta_obsf )
  call move_alloc( o_rrw   ,    rrw_obsf )
  call move_alloc( o_uzonal, uzonal_obsf )
  call move_alloc( o_umerid, umerid_obsf )

  call mass_sum_from_rho( future_mass, rho_obsf, allnodes=.true., mask=skipw )

end subroutine nudge_prep_obs

!===============================================================================

subroutine obs_nudge()

  use mem_nudge,   only: tnudcent,    rhot_nud,    &
                         rho_obsp,    rho_obsf,    &
                         theta_obsp,  theta_obsf,  &
                         rrw_obsp,    rrw_obsf,    &
                         uzonal_obsp, uzonal_obsf, &
                         umerid_obsp, umerid_obsf, nudflag, nudnxp

  use mem_basic,   only: rho, theta, rr_v, ue, ve, vxe, vye, vze
  use mem_grid,    only: mza, lpw, vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
  use misc_coms,   only: s1900_sim, io6
  use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog
  use mem_tend,    only: thilt, rr_wt, vmxet, vmyet, vmzet
  use isan_coms,   only: ifgfile, s1900_fg
  use oname_coms,  only: nl
  use consts_coms, only: r8

  use check_nan

  implicit none

! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data

  integer  :: j, iw, k
  real     :: umzonalt, ummeridt
  real     :: tp, tf, tnudi, tnudr, rho4
  real(r8) :: tpr, tfr, rho_obs
  real     :: theta_obs, rrv_obs, uzonal_obs, umerid_obs

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK
!
! integer, save :: icall = 0
! real, save, allocatable :: wtnud(:)
!
! if (icall /= 1) then
!    icall = 1
!
!    allocate (wtnud(mwa))
!
! ! Horizontal loop over T points
!
!    do iw = 2, mwa
!
! ! Default: Uniform nudging weight = 1
!
!       wtnud(iw) = 1.
!
! ! Sample code for modifying nudging weight
!
! ! Transform current IW point to polar stereographic coordinates using specified
! ! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
!
!       call e_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)
!
!       dist = sqrt(xw ** 2 + yw ** 2)
!
!       if (dist > 5000.e3) then
!          wtnud(iw) = 1.
!       elseif (dist < 3600.e3) then
!          wtnud(iw) = 0.
!       else
!          wtnud(iw) = ((dist - 3600.e3) / 1400.e3) ** 2
!       endif
!
!    enddo
!
! endif
!----------------------------------------------------------------------

  if (nudflag < 1) return
  if (nudflag > 0 .and. nudnxp > 0) return

  ! Time interpolation coefficients

  tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
            / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

  tp = 1. - tf

  tnudi = 1.0 / tnudcent

  if (nl%nud_preserve_total_mass) then

     if (maxmrl > nl%max_nud_mrl) then
        call mass_sum_from_rho( current_mass, rho, allnodes=.true., mask=skipw )
     endif

     tpr = (1._r8 - real(tf,r8)) * current_mass /   past_mass
     tfr =          real(tf,r8)  * current_mass / future_mass

  else

     tpr = (1._r8 - real(tf,r8))
     tfr =          real(tf,r8)

  endif

  ! Horizontal loop over W columns
  !----------------------------------------------------------------------
  !$omp parallel do private( iw,k,rho4,tnudr,umzonalt,ummeridt, &
  !$omp                      rho_obs, theta_obs, rrv_obs, uzonal_obs, umerid_obs )
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
  !---------------------------------------------------------------------

     ! Skip obs nudging if mrl of this column is greater than max nudging mrl

     if (itab_w(iw)%mrlw > nl%max_nud_mrl) cycle

     ! Inverse of nudging time scale
     ! If using spatial nudging mask defined above in this subroutine
     ! (will also need to add tnudi to list of OpenMP private vars):

     ! tnudi = wtnud(iw) / tnudcent

     ! Compute nudging tendencies, interpolating observational fields in time

     do k = lpw(iw), mza

        rho_obs    = tpr *    rho_obsp(k,iw) + tfr *    rho_obsf(k,iw)
        theta_obs  = tp  *  theta_obsp(k,iw) + tf  *  theta_obsf(k,iw)
        rrv_obs    = tp  *    rrw_obsp(k,iw) + tf  *    rrw_obsf(k,iw)
        uzonal_obs = tp  * uzonal_obsp(k,iw) + tf  * uzonal_obsf(k,iw)
        umerid_obs = tp  * umerid_obsp(k,iw) + tf  * umerid_obsf(k,iw)

        rho4  = rho(k,iw)
        tnudr = rho4 * tnudi

        ! Density nudging

        rhot_nud(k,iw) = tnudi * (rho_obs - rho(k,iw))

        ! Heat (rho Theta) nudging

        thilt(k,iw) = thilt(k,iw) &
                    + tnudi * (rho_obs * theta_obs - rho4 * theta(k,iw))

        ! Total water mixing ratio nudging

        rr_wt(k,iw) = rr_wt(k,iw) + tnudr * (rrv_obs - rr_v(k,iw))

        ! Momentum nudging (including density change)

        umzonalt = tnudr * (uzonal_obs - ue(k,iw))

        ummeridt = tnudr * (umerid_obs - ve(k,iw))

        vmxet(k,iw) = vmxet(k,iw) + vxn_ew(iw) * umzonalt + vxn_ns(iw) * ummeridt &
                    + vxe(k,iw) * rhot_nud(k,iw)

        vmyet(k,iw) = vmyet(k,iw) + vyn_ew(iw) * umzonalt + vyn_ns(iw) * ummeridt &
                    + vye(k,iw) * rhot_nud(k,iw)

        vmzet(k,iw) = vmzet(k,iw) + vzn_ns(iw) * ummeridt &
                    + vze(k,iw) * rhot_nud(k,iw)
     enddo

  enddo
  !$omp end parallel do

end subroutine obs_nudge

end module obs_nudge_mod
