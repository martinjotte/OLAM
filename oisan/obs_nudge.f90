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

subroutine nudge_prep_obs(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid)

use mem_nudge,   only: rho_obsp, theta_obsp, shw_obsp, &
                       rho_obsf, theta_obsf, shw_obsf, &
                       uzonal_obsp, umerid_obsp,       &
                       uzonal_obsf, umerid_obsf

use mem_grid,    only: mza, mwa, lpw
use mem_ijtabs,  only: jtab_w, jtw_init
use consts_coms, only: r8

implicit none

integer,  intent(in) :: iaction
real(r8), intent(in) :: o_rho   (mza,mwa)
real,     intent(in) :: o_theta (mza,mwa)
real,     intent(in) :: o_shv   (mza,mwa)
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
            shw_obsp(k,iw) =    shw_obsf(k,iw)
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
         shw_obsf(k,iw) = o_shv   (k,iw)
      uzonal_obsf(k,iw) = o_uzonal(k,iw)
      umerid_obsf(k,iw) = o_umerid(k,iw)
   enddo

enddo
!$omp end do

! End OpenMP parallel block
!$omp end parallel

end subroutine nudge_prep_obs

!===========================================================================

subroutine obs_nudge(rhot)

use mem_nudge, only:   tnudcent,                                       &
                        rho_obs, rho_obsp,    rho_obsf,    &
                      theta_obs, theta_obsp,  theta_obsf,  &
                        shw_obs, shw_obsp,    shw_obsf,    &
                     uzonal_obs, uzonal_obsp, uzonal_obsf, &
                     umerid_obs, umerid_obsp, umerid_obsf

use mem_basic,   only: rho, theta, sh_w, vxe, vye, vze
use mem_grid,    only: mza, mwa, lpw, xew, yew, zew
use misc_coms,   only: s1900_sim
use mem_ijtabs,  only: istp, jtab_w, mrl_begl, jtw_prog
use consts_coms, only: eradi
use mem_tend,    only: thilt, sh_wt, vmxet, vmyet, vmzet
use isan_coms,   only: ifgfile, s1900_fg
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w

implicit none

! Nudge selected model fields (rho, thil, sh_w, vmc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: j, iw, k, mrl

real :: umzonalt(mza)
real :: ummeridt(mza)
real :: uzonal  (mza)
real :: umerid  (mza)

real :: tp,tf,tnudi,tnudirho
real :: raxis,raxisi,uvtr

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

mrl = mrl_begl(istp)
if (mrl < 1) return

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf

! Horizontal loop over W columns

!----------------------------------------------------------------------
!$omp parallel do private(iw,raxis,raxisi,k,uzonal,umerid,tnudi,tnudirho,
!$omp                     umzonalt,ummeridt,uvtr)
do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!---------------------------------------------------------------------

! Reconstruct UZONAL(k) and UMERID(k) from VXE, VYE, VZE

   raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
   raxisi = 1.0 / max(raxis, 1.e-12)

   if (raxis > 1.e3) then

      do k = lpw(iw), mza

         uzonal(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
         umerid(k) = vze(k,iw) * raxis * eradi &
            - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi

      enddo

   else
      uzonal(:) = 0.
      umerid(:) = 0.
   endif

! Interpolate observational fields in time

   do k = lpw(iw), mza
         rho_obs(k,iw) = tp *    rho_obsp(k,iw) + tf *    rho_obsf(k,iw)
       theta_obs(k,iw) = tp *  theta_obsp(k,iw) + tf *  theta_obsf(k,iw)
         shw_obs(k,iw) = tp *    shw_obsp(k,iw) + tf *    shw_obsf(k,iw)
      uzonal_obs(k,iw) = tp * uzonal_obsp(k,iw) + tf * uzonal_obsf(k,iw)
      umerid_obs(k,iw) = tp * umerid_obsp(k,iw) + tf * umerid_obsf(k,iw)
   enddo

! Inverse of nudging time scale
! Use spatial nudging mask defined above in this subroutine

!  tnudi = wtnud(iw) / tnudcent
   tnudi = 1.0 / tnudcent

! Compute nudging tendencies, interpolating observational fields in time

   do k = lpw(iw), mza

      tnudirho = tnudi * rho(k,iw)

          rhot(k,iw) = rhot(k,iw) &
                     + tnudi * (rho_obs(k,iw) - rho(k,iw))

         thilt(k,iw) = thilt(k,iw) &
                     + tnudirho * (theta_obs(k,iw) - theta(k,iw))
                 
         sh_wt(k,iw) = sh_wt(k,iw) &
                     + tnudirho * (shw_obs(k,iw) - sh_w(k,iw))

      umzonalt(k)    = tnudirho * (uzonal_obs(k,iw) - uzonal(k))

      ummeridt(k)    = tnudirho * (umerid_obs(k,iw) - umerid(k))

   enddo

   if (raxis > 1.e3) then
      do k = lpw(iw), mza
         uvtr = -ummeridt(k) * zew(iw) * eradi
         vmxet(k,iw) = vmxet(k,iw) + (-umzonalt(k) * yew(iw) + uvtr * xew(iw)) * raxisi
         vmyet(k,iw) = vmyet(k,iw) + ( umzonalt(k) * xew(iw) + uvtr * yew(iw)) * raxisi
         vmzet(k,iw) = vmzet(k,iw) +   ummeridt(k) * raxis * eradi 
      enddo
   endif

enddo
!$omp end parallel do

end subroutine obs_nudge
