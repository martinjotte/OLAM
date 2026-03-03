module obs_nudge_mod

contains

!===============================================================================

subroutine nudge_prep_obs( iaction )

  use mem_nudge,   only: rho_obsp, theta_obsp, rrv_obsp, uzonal_obsp, umerid_obsp, &
                         rho_obsf, theta_obsf, rrv_obsf, uzonal_obsf, umerid_obsf, &
                         nudflag, nudnxp
  use consts_coms, only: r8
  use isan_coms,   only: o_rho, o_theta, o_rrw, o_uzonal, o_umerid
  use mem_basic,   only: rho
  use check_nan,   only: mass_sum_from_rho
  use misc_coms,   only: io6
  use mem_grid,    only: mwa, mza, lpw
  use mem_ijtabs,  only: jtab_w, itab_w, jtw_prog
  use oname_coms,  only: nl

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iaction
  integer             :: j, iw, k, ier
  real                :: fact
  real(r8)            :: future_mass, current_mass

  ! Make sure we are doing obs/gridpoint nudging

  if (nudflag < 1) return
  if (nudflag > 0 .and. nudnxp > 0) return

  ! Swap future data time into past data time if necessary

  if (iaction == 1) then
     call move_alloc(    rho_obsf,    rho_obsp )
     call move_alloc(  theta_obsf,  theta_obsp )
     call move_alloc(    rrv_obsf,    rrv_obsp )
     call move_alloc( uzonal_obsf, uzonal_obsp )
     call move_alloc( umerid_obsf, umerid_obsp )
  endif

  ! Swap gridded data fields into future nudging arrays

  call move_alloc( o_rho   ,    rho_obsf )
  call move_alloc( o_theta ,  theta_obsf )
  call move_alloc( o_rrw   ,    rrv_obsf )
  call move_alloc( o_uzonal, uzonal_obsf )
  call move_alloc( o_umerid, umerid_obsf )

  ! Scale rho_obsf so that the analysis global mass sum equals the current
  ! simulation's atmospheric mass sum

  call mass_sum_from_rho( current_mass, rho,      allnodes=.true. )
  call mass_sum_from_rho(  future_mass, rho_obsf, allnodes=.true. )

  fact = real( current_mass / future_mass )

  !$omp parallel do private(iw,k)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw), mza
        rho_obsf(k,iw) = rho_obsf(k,iw) * fact
     enddo
  enddo
  !$omp end parallel do

end subroutine nudge_prep_obs

!===============================================================================

subroutine obs_nudge()

  use mem_nudge,   only: tnudcent,    rhot_nud,    &
                         rho_obsp,    rho_obsf,    &
                         theta_obsp,  theta_obsf,  &
                         rrv_obsp,    rrv_obsf,    &
                         uzonal_obsp, uzonal_obsf, &
                         umerid_obsp, umerid_obsf, &
                         tf, tp, nudflag, nudnxp, wmanud, vmanud, vovera
  use mem_basic,   only: rho, theta, thil, rr_v, ue, ve, vxe, vye, vze
  use mem_grid,    only: mza, lpv, lpw, arw, arw0, arv, dzm, &
                         vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns, volti
  use misc_coms,   only: s1900_sim, io6, iparallel
  use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog, itab_v, jtab_v, jtv_prog
  use mem_tend,    only: thilt, rr_wt, vmxet, vmyet, vmzet
  use isan_coms,   only: ifgfile, s1900_fg
  use oname_coms,  only: nl
  use consts_coms, only: r8
  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

  implicit none

! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data

  integer  :: j, iw, k, iv, iw1, iw2
  real     :: umzonalt, ummeridt
  real     :: tnudi, tnudr, dxi
  real(r8) :: rho_obs
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
! ! pole point location (pole point lat/lon = 4th & 5th arguments of ec_ps)
!
!       call ec_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)
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

  ! Horizontal loop over W columns
  !----------------------------------------------------------------------
  !$omp parallel do private( iw,k,tnudr,dxi,umzonalt,ummeridt, &
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

        rho_obs    = tp *    rho_obsp(k,iw) + tf *    rho_obsf(k,iw)
        theta_obs  = tp *  theta_obsp(k,iw) + tf *  theta_obsf(k,iw)
        rrv_obs    = tp *    rrv_obsp(k,iw) + tf *    rrv_obsf(k,iw)
        uzonal_obs = tp * uzonal_obsp(k,iw) + tf * uzonal_obsf(k,iw)
        umerid_obs = tp * umerid_obsp(k,iw) + tf * umerid_obsf(k,iw)

        tnudr = real(rho(k,iw)) * tnudi

        ! Density nudging

        rhot_nud(k,iw) = tnudi * (rho_obs - rho(k,iw))

        ! Theta nudging

        thilt(k,iw) = thilt(k,iw) + tnudr * (theta_obs - theta(k,iw))

        ! Total water mixing ratio nudging

        rr_wt(k,iw) = rr_wt(k,iw) + tnudr * (rrv_obs - rr_v(k,iw))

        ! Momentum nudgingq

        umzonalt = tnudr * (uzonal_obs - ue(k,iw))
        ummeridt = tnudr * (umerid_obs - ve(k,iw))

        vmxet(k,iw) = vmxet(k,iw) + vxn_ew(iw) * umzonalt + vxn_ns(iw) * ummeridt
        vmyet(k,iw) = vmyet(k,iw) + vyn_ew(iw) * umzonalt + vyn_ns(iw) * ummeridt
        vmzet(k,iw) = vmzet(k,iw) +                         vzn_ns(iw) * ummeridt

     enddo

     if (.not. nl%nud_preserve_total_mass) then

        ! Include effect of density change on momentum and heat
        do k = lpw(iw), mza
           vmxet(k,iw) = vmxet(k,iw) + vxe (k,iw) * rhot_nud(k,iw)
           vmyet(k,iw) = vmyet(k,iw) + vye (k,iw) * rhot_nud(k,iw)
           vmzet(k,iw) = vmzet(k,iw) + vze (k,iw) * rhot_nud(k,iw)
           thilt(k,iw) = thilt(k,iw) + thil(k,iw) * rhot_nud(k,iw)
        enddo

     else

        rhot_nud(lpw(iw)-1,iw) = 0._r8
        do k = lpw(iw), mza
           rhot_nud(k,iw) = rhot_nud(k,iw) * vovera(k,iw)
        enddo

        dxi = 1. / sqrt(arw0(iw))
        do k = lpw(iw), mza-1
           wmanud(k,iw) = real(rhot_nud(k+1,iw) - rhot_nud(k,iw)) * dzm(k) * dxi * arw(k,iw)
        enddo

     endif

  enddo
  !$omp end parallel do

  if (nl%nud_preserve_total_mass) then

     if (iparallel == 1) then
        call mpi_send_w(dvara1=rhot_nud)
        call mpi_recv_w(dvara1=rhot_nud)
     endif

     !$omp parallel do private(iv,iw1,iw2,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if ( itab_w(iw1)%mrlw > nl%max_nud_mrl .and. &
             itab_w(iw2)%mrlw > nl%max_nud_mrl ) cycle

        do k = lpv(iv), mza
           vmanud(k,iv) = real(rhot_nud(k,iw2) - rhot_nud(k,iw1)) * arv(k,iv)
        enddo

     enddo
     !$omp end parallel do

     if (iparallel == 1) then
        call mpi_send_v(rvara1=vmanud)
        call mpi_recv_v(rvara1=vmanud)
     endif

  endif

end subroutine obs_nudge

end module obs_nudge_mod
