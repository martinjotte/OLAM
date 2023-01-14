subroutine nudge_prep_obs(iaction, o_rho, o_theta, o_rrw, o_uzonal, o_umerid)

use mem_nudge,   only: rho_obsp, theta_obsp, rrw_obsp, &
                       rho_obsf, theta_obsf, rrw_obsf, &
                       uzonal_obsp, umerid_obsp,       &
                       uzonal_obsf, umerid_obsf
use mem_grid,    only: mza, mwa, lpw
use mem_ijtabs,  only: jtab_w, jtw_init
use consts_coms, only: r8

implicit none

integer, intent(in) :: iaction
real,    intent(in) :: o_rho   (mza,mwa)
real,    intent(in) :: o_theta (mza,mwa)
real,    intent(in) :: o_rrw   (mza,mwa)
real,    intent(in) :: o_uzonal(mza,mwa)
real,    intent(in) :: o_umerid(mza,mwa)
integer             :: j, iw, k

! Begin OpenMP parallel block
!$omp parallel

! Swap future data time into past data time if necessary.

if (iaction == 1) then

!----------------------------------------------------------------------
   !$omp do private(iw,k)
   do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------
      do k = lpw(iw), mza
            rho_obsp(k,iw) =    rho_obsf(k,iw)
          theta_obsp(k,iw) =  theta_obsf(k,iw)
            rrw_obsp(k,iw) =    rrw_obsf(k,iw)
         uzonal_obsp(k,iw) = uzonal_obsf(k,iw)
         umerid_obsp(k,iw) = umerid_obsf(k,iw)

      enddo
   enddo
   !$omp end do

endif

! If doing point-by-point (non-spectral) nudging, fill future observational
! arrays with simple memory copy

! Horizontal loop over W columns

!----------------------------------------------------------------------
!$omp do private(iw,k)
do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

   do k = lpw(iw), mza
         rho_obsf(k,iw) = o_rho   (k,iw)
       theta_obsf(k,iw) = o_theta (k,iw)
         rrw_obsf(k,iw) = o_rrw   (k,iw)
      uzonal_obsf(k,iw) = o_uzonal(k,iw)
      umerid_obsf(k,iw) = o_umerid(k,iw)
   enddo

enddo
!$omp end do

! End OpenMP parallel block
!$omp end parallel

end subroutine nudge_prep_obs

!===============================================================================

subroutine obs_nudge()

use mem_nudge, only:   tnudcent, rhot_nud,                 &
                        rho_obs, rho_obsp,    rho_obsf,    &
                      theta_obs, theta_obsp,  theta_obsf,  &
                        rrw_obs, rrw_obsp,    rrw_obsf,    &
                     uzonal_obs, uzonal_obsp, uzonal_obsf, &
                     umerid_obs, umerid_obsp, umerid_obsf

use mem_basic,   only: rho, theta, rr_w, ue, ve, vxe, vye, vze
use mem_grid,    only: mza, lpw, vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
use misc_coms,   only: s1900_sim
use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog
use mem_tend,    only: thilt, rr_wt, vmxet, vmyet, vmzet
use isan_coms,   only: ifgfile, s1900_fg
use oname_coms,  only: nl

implicit none

! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data

integer :: j, iw, k
real    :: umzonalt, ummeridt
real    :: tp, tf, tnudi, tnudr, rho4

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK
!
!integer, save :: icall = 0
!real, save, allocatable :: wtnud(:)
!
!if (icall /= 1) then
!   icall = 1
!
!   allocate (wtnud(mwa))
!
!! Horizontal loop over T points
!
!   do iw = 2, mwa
!
!! Default: Uniform nudging weight = 1
!
!      wtnud(iw) = 1.
!
!! Sample code for modifying nudging weight
!
!! Transform current IW point to polar stereographic coordinates using specified
!! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
!
!      call e_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)
!
!      dist = sqrt(xw ** 2 + yw ** 2)
!
!      if (dist > 5000.e3) then
!         wtnud(iw) = 1.
!      elseif (dist < 3600.e3) then
!         wtnud(iw) = 0.
!      else
!         wtnud(iw) = ((dist - 3600.e3) / 1400.e3) ** 2
!      endif
!
!   enddo
!
!endif
!----------------------------------------------------------------------

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf

! Horizontal loop over W columns
!----------------------------------------------------------------------
!$omp parallel do private( iw,k,tnudi,rho4,tnudr,umzonalt,ummeridt )
do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
!---------------------------------------------------------------------

! Skip obs nudging if mrl of this column is greater than max nudging mrl

   if (itab_w(iw)%mrlw > nl%max_nud_mrl) cycle

! Interpolate observational fields in time

   do k = lpw(iw), mza
         rho_obs(k,iw) = tp *    rho_obsp(k,iw) + tf *    rho_obsf(k,iw)
       theta_obs(k,iw) = tp *  theta_obsp(k,iw) + tf *  theta_obsf(k,iw)
         rrw_obs(k,iw) = tp *    rrw_obsp(k,iw) + tf *    rrw_obsf(k,iw)
      uzonal_obs(k,iw) = tp * uzonal_obsp(k,iw) + tf * uzonal_obsf(k,iw)
      umerid_obs(k,iw) = tp * umerid_obsp(k,iw) + tf * umerid_obsf(k,iw)
   enddo

! Inverse of nudging time scale
! Use spatial nudging mask defined above in this subroutine

!  tnudi = wtnud(iw) / tnudcent
   tnudi = 1.0 / tnudcent

! Compute nudging tendencies, interpolating observational fields in time

   do k = lpw(iw), mza

      rho4  = rho(k,iw)
      tnudr = rho4 * tnudi

      ! Density nudging

      rhot_nud(k,iw) = tnudi * (rho_obs(k,iw) - rho4)

      ! Heat (rho Theta) nudging

      thilt(k,iw) = thilt(k,iw) &
                  + tnudi * (rho_obs(k,iw) * theta_obs(k,iw) - rho4 * theta(k,iw))

      ! Total water mixing ratio nudging

      rr_wt(k,iw) = tnudr * (rrw_obs(k,iw) - rr_w(k,iw))

      ! Momentum nudging (including density change)

      umzonalt = tnudr * (uzonal_obs(k,iw) - ue(k,iw))

      ummeridt = tnudr * (umerid_obs(k,iw) - ve(k,iw))

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
