!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================

subroutine nudge_prep_obs(iaction, o_rho, o_theta, o_rrw, o_uzonal, o_umerid)

use mem_nudge,   only: rho_obsp, theta_obsp, rrw_obsp, &
                       rho_obsf, theta_obsf, rrw_obsf, &
                       uzonal_obsp, umerid_obsp,       &
                       uzonal_obsf, umerid_obsf

use mem_grid,    only: mza, mwa, lpw
use mem_ijtabs,  only: jtab_w, jtw_init
use consts_coms, only: r8

implicit none

integer,  intent(in) :: iaction
real(r8), intent(in) :: o_rho   (mza,mwa)
real,     intent(in) :: o_theta (mza,mwa)
real,     intent(in) :: o_rrw   (mza,mwa)
real,     intent(in) :: o_uzonal(mza,mwa)
real,     intent(in) :: o_umerid(mza,mwa)
integer              :: j,iw,k

! Begin OpenMP parallel block
!$omp parallel

! Swap future data time into past data time if necessary.

if (iaction == 1) then

!----------------------------------------------------------------------
   !$omp do private(iw,k)
   do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
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
do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
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

subroutine obs_nudge(mrl)

use mem_nudge, only:   tnudcent, rhot_nud,                 &
                        rho_obs, rho_obsp,    rho_obsf,    &
                      theta_obs, theta_obsp,  theta_obsf,  &
                        rrw_obs, rrw_obsp,    rrw_obsf,    &
                     uzonal_obs, uzonal_obsp, uzonal_obsf, &
                     umerid_obs, umerid_obsp, umerid_obsf

use mem_basic,   only: rho, theta, rr_w, vxe, vye, vze
use mem_grid,    only: mza, lpw, vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
use misc_coms,   only: s1900_sim, dtlm
use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog
use mem_tend,    only: thilt, rr_wt, vmxet, vmyet, vmzet
use isan_coms,   only: ifgfile, s1900_fg
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w
use oname_coms,  only: nl

implicit none

! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data

integer, intent(in) :: mrl

integer :: j, iw, k
real    :: umzonalt, ummeridt
real    :: uzonal, umerid
real    :: tp, tf, tnudi, tnudr, rho4, dti, rrw_nudget

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

! Check whether it is time to nudge

if (mrl < 1) return

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf

! Horizontal loop over W columns
!----------------------------------------------------------------------
!$omp parallel do private( iw,k,tnudi,dti,rho4,tnudr,uzonal,umerid,rrw_nudget, &
!$omp                      umzonalt,ummeridt )
do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
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

   dti = 1.0 / real(dtlm(itab_w(iw)%mrlw))

   do k = lpw(iw), mza

      rho4  = rho(k,iw)
      tnudr = rho4 * tnudi

      ! Reconstruct UZONAL and UMERID from VXE, VYE, VZE

      uzonal = vxe(k,iw) * vxn_ew(iw) + vye(k,iw) * vyn_ew(iw)

      umerid = vxe(k,iw) * vxn_ns(iw) + vye(k,iw) * vyn_ns(iw) &
             + vze(k,iw) * vzn_ns(iw)

      ! Density nudging

      rhot_nud(k,iw) = tnudi * (rho_obs(k,iw) - rho4)

      ! Heat (rho Theta) nudging

      thilt(k,iw) = thilt(k,iw) &
                  + tnudi * (rho_obs(k,iw) *  theta_obs(k,iw) - rho4 * theta(k,iw))

      ! Total water mixing ratio nudging

      rrw_nudget = tnudr * (rrw_obs(k,iw) - rr_w(k,iw))

      rrw_nudget = max(rrw_nudget, -.95 * max(rr_w(k,iw) * rho4 * dti + rr_wt(k,iw), 0.))

      rr_wt(k,iw) = rr_wt(k,iw) + rrw_nudget

      ! Momentum nudging (including density change)

      umzonalt = tnudr * (uzonal_obs(k,iw) - uzonal)

      ummeridt = tnudr * (umerid_obs(k,iw) - umerid)

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
