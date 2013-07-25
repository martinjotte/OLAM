!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

subroutine nudge_prep_obs(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid)

use mem_nudge,   only: nudflag, nudnxp,                &
                       rho_obsp, theta_obsp, shw_obsp, &
                       rho_obsf, theta_obsf, shw_obsf, &
                       uzonal_obsp, umerid_obsp,       &
                       uzonal_obsf, umerid_obsf

use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use mem_ijtabs, only: jtab_w, jtw_init

implicit none

integer,      intent(in) :: iaction
real(kind=8), intent(in) :: o_rho   (mza,mwa)
real,         intent(in) :: o_theta (mza,mwa)
real,         intent(in) :: o_shv   (mza,mwa)
real,         intent(in) :: o_uzonal(mza,mwa)
real,         intent(in) :: o_umerid(mza,mwa)
integer                  :: j,iw,k

! Begin OpenMP parallel block
!$omp parallel

! Swap future data time into past data time if necessary.

if (iaction == 1) then

!----------------------------------------------------------------------
   !$omp do private(iw,k)
   do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------
      do k = 2, mza-1
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

call psub()
!----------------------------------------------------------------------
!$omp do private(iw,k)
do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------
   call qsub('W',iw)

   do k = 2, mza-1
         rho_obsf(k,iw) = o_rho   (k,iw)
       theta_obsf(k,iw) = o_theta (k,iw)
         shw_obsf(k,iw) = o_shv   (k,iw)
      uzonal_obsf(k,iw) = o_uzonal(k,iw)
      umerid_obsf(k,iw) = o_umerid(k,iw)
   enddo

enddo
!$omp end do
call rsub('Wbb',7)

! End OpenMP parallel block
!$omp end parallel

end subroutine nudge_prep_obs

!===========================================================================

subroutine obs_nudge(rhot)

use mem_nudge, only: nudnxp, tnudcent, mwnud,                          &
                        rho_sim,    rho_obs, rho_obsp,    rho_obsf,    &
                      theta_sim,  theta_obs, theta_obsp,  theta_obsf,  &
                        shw_sim,    shw_obs, shw_obsp,    shw_obsf,    &
                     uzonal_sim, uzonal_obs, uzonal_obsp, uzonal_obsf, &
                     umerid_sim, umerid_obs, umerid_obsp, umerid_obsf

use mem_basic,   only: uc, vc, rho, theta, sh_w, vxe, vye, vze
use mem_grid,    only: mza, mwa, lpu, lpv, lpw, &
                       xeu, yeu, zeu, xev, yev, zev, xew, yew, zew, &
                       unx, uny, unz, vnx, vny, vnz, volt, glatw, glonw
use misc_coms,   only: io6, time8, s1900_sim, meshtype, iparallel
use mem_ijtabs,  only: istp, jtab_u, jtab_v, jtab_w, itab_u, itab_v, itab_w, &
                       mrl_begl, jtu_prog, jtv_prog, jtw_prog
use consts_coms, only: eradi
use mem_tend,    only: umt, vmt, thilt, sh_wt
use isan_coms,   only: ifgfile, s1900_fg
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w

implicit none

! Nudge selected model fields (rho, thil, sh_w, umc, vmc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: j, iw, k, mrl, iu, iv, iw1, iw2

real :: umzonalt(mza,mwa)
real :: ummeridt(mza,mwa)
real :: uzonal  (mza)
real :: umerid  (mza)

real :: tp,tf,tnudi,tnudirho
real :: raxis,raxisi
real :: umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK

integer, save :: icall = 0
real, save, allocatable :: wtnud(:)
real :: xw,yw,dist

if (icall /= 1) then
   icall = 1
   
   allocate (wtnud(mwa))
   
! Horizontal loop over T points

   do iw = 2, mwa

! Default: Uniform nudging weight = 1

      wtnud(iw) = 1.

! Sample code for modifying nudging weight
      
! Transform current IW point to polar stereographic coordinates using specified
! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
   
!      call e_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)

!      dist = sqrt(xw ** 2 + yw ** 2)

!      if (dist > 5000.e3) then
!         wtnud(iw) = 1.
!      elseif (dist < 3600.e3) then
!         wtnud(iw) = 0.
!      else
!         wtnud(iw) = ((dist - 3600.e3) / 1400.e3) ** 2
!      endif

   enddo

endif
!----------------------------------------------------------------------

! Check whether it is time to nudge

mrl = mrl_begl(istp)
if (mrl < 1) return

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,raxis,raxisi,k,uzonal,umerid,tnudi,tnudirho)
do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!---------------------------------------------------------------------
call qsub('W',iw)

! Reconstruct UZONAL(k) and UMERID(k) from VXE, VYE, VZE

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = lpw(iw), mza-1

         uzonal(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
         umerid(k) = vze(k,iw) * raxis * eradi &
            - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi

      enddo

   else
      uzonal(:) = 0.
      umerid(:) = 0.
   endif

! Interpolate observational fields in time

   do k = lpw(iw), mza-1
         rho_obs(k,iw) = tp *    rho_obsp(k,iw) + tf *    rho_obsf(k,iw)
       theta_obs(k,iw) = tp *  theta_obsp(k,iw) + tf *  theta_obsf(k,iw)
         shw_obs(k,iw) = tp *    shw_obsp(k,iw) + tf *    shw_obsf(k,iw)
      uzonal_obs(k,iw) = tp * uzonal_obsp(k,iw) + tf * uzonal_obsf(k,iw)
      umerid_obs(k,iw) = tp * umerid_obsp(k,iw) + tf * umerid_obsf(k,iw)
   enddo

! Inverse of nudging time scale
! Use spatial nudging mask defined above in this subroutine

   tnudi = wtnud(iw) / tnudcent

! Compute nudging tendencies, interpolating observational fields in time

   do k = lpw(iw), mza-1

      tnudirho = tnudi * rho(k,iw)

          rhot(k,iw) = rhot(k,iw) &
                     + tnudi * (rho_obs(k,iw) - rho(k,iw))

         thilt(k,iw) = thilt(k,iw) &
                     + tnudirho * (theta_obs(k,iw) - theta(k,iw))
                 
         sh_wt(k,iw) = sh_wt(k,iw) &
                     + tnudirho * (shw_obs(k,iw) - sh_w(k,iw))

      umzonalt(k,iw) = tnudirho * (uzonal_obs(k,iw) - uzonal(k))

      ummeridt(k,iw) = tnudirho * (umerid_obs(k,iw) - umerid(k))
   enddo

enddo
!$omp end parallel do
call rsub('Wb',23)

if (iparallel == 1) then
   call mpi_send_w('V', vxe=umzonalt, vye=ummeridt, domrl=mrl)
   call mpi_recv_w('V', vxe=umzonalt, vye=ummeridt, domrl=mrl)
endif

if (meshtype == 1) then

! UMT

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iu,iw1,iw2,k, &
   !$omp                     raxis,raxisi,umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt)
   do j = 1,jtab_u(jtu_prog)%jend(mrl); iu = jtab_u(jtu_prog)%iu(j)
      iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then

         raxisi = 1. / raxis

! Average momentum tendencies (times volume) to U point and rotate at U point

         do k = lpu(iu), mza-1

            umgt = .5 * (umzonalt(k,iw1) * volt(k,iw1) &
                       + umzonalt(k,iw2) * volt(k,iw2))

            vmgt = .5 * (ummeridt(k,iw1) * volt(k,iw1) &
                       + ummeridt(k,iw2) * volt(k,iw2))

            uvmgrt = -vmgt * zeu(iu) * eradi  ! radially outward from axis

            uvmgxt = (-umgt * yeu(iu) + uvmgrt * xeu(iu)) * raxisi
            uvmgyt = ( umgt * xeu(iu) + uvmgrt * yeu(iu)) * raxisi
            uvmgzt =   vmgt * raxis * eradi 

            umt(k,iu) = umt(k,iu) &
                      + uvmgxt * unx(iu) + uvmgyt * uny(iu) + uvmgzt * unz(iu)

         enddo

      endif

   enddo
   !$omp end parallel do
   call rsub('U',13)

else

! VMT

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iv,iw1,iw2,k, &
   !$omp                     raxis,raxisi,umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt)
   do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
      iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   call qsub('V',iv)

      raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
      
         raxisi = 1. / raxis

! Average momentum tendencies to V point and rotate at V point

         do k = lpv(iv), mza-1
            umgt = .5 * (umzonalt(k,iw1) + umzonalt(k,iw2))

            vmgt = .5 * (ummeridt(k,iw1) + ummeridt(k,iw2))

            uvmgrt = -vmgt * zev(iv) * eradi  ! radially outward from axis

            uvmgxt = (-umgt * yev(iv) + uvmgrt * xev(iv)) * raxisi 
            uvmgyt = ( umgt * xev(iv) + uvmgrt * yev(iv)) * raxisi 
            uvmgzt =   vmgt * raxis * eradi 

            vmt(k,iv) = vmt(k,iv) &
                      + uvmgxt * vnx(iv) + uvmgyt * vny(iv) + uvmgzt * vnz(iv)

         enddo

      endif

   enddo
   !$omp end parallel do
   call rsub('V',13)

endif

end subroutine obs_nudge
