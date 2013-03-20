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
subroutine nudge_prep(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid)

use mem_nudge,   only: nudflag, nudnxp, mwnud, &
                       rho_obsp, theta_obsp, shw_obsp, &
                       rho_obsf, theta_obsf, shw_obsf, &
                       uzonal_obsp, umerid_obsp, &
                       uzonal_obsf, umerid_obsf

use mem_grid,    only: mza, mva, mwa, lpw, vnx, vny, vnz, xeu, yeu, zeu, &
                       xev, yev, zev, xew, yew, zew, volt
use misc_coms,   only: io6, iparallel, runtype, meshtype
use mem_ijtabs,  only: jtab_w, itab_w, jtw_init

implicit none

integer, intent(in) :: iaction

real(kind=8), intent(in) :: o_rho   (mza,mwa)
real,         intent(in) :: o_theta (mza,mwa)
real,         intent(in) :: o_shv   (mza,mwa)
real,         intent(in) :: o_uzonal(mza,mwa)
real,         intent(in) :: o_umerid(mza,mwa)

integer :: j,iw,k,ka,iu,iv,iw1,iw2,iup,ivp,iwnud,iwnud1
integer :: npoly,kb,jv
real    :: volwnudi

! Automatic arrays

real :: volwnud(mza,mwnud)

! Swap future data time into past data time if necessary.

if (iaction == 1) then

   do iwnud = 2,mwnud
      do k = 2,mza-1
            rho_obsp(k,iwnud) =    rho_obsf(k,iwnud)
          theta_obsp(k,iwnud) =  theta_obsf(k,iwnud)
            shw_obsp(k,iwnud) =    shw_obsf(k,iwnud)
         uzonal_obsp(k,iwnud) = uzonal_obsf(k,iwnud)
         umerid_obsp(k,iwnud) = umerid_obsf(k,iwnud)
      enddo
   enddo

endif

! If doing spectral nudging, zero out nudging polygon arrays and volume counter
! prior to summing

if (nudnxp > 0) then

   do iwnud = 2,mwnud
      do k = 2,mza-1
            rho_obsf(k,iwnud) = 0.
          theta_obsf(k,iwnud) = 0.
            shw_obsf(k,iwnud) = 0.
         uzonal_obsf(k,iwnud) = 0.
         umerid_obsf(k,iwnud) = 0.

             volwnud(k,iwnud) = 0.
      enddo
   enddo

endif

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------
call qsub('W',iw)

! If doing point-by-point (non-spectral) nudging, fill future observational
! arrays with simple memory copy

   if (nudnxp == 0) then

      do k = 2,mza-1
            rho_obsf(k,iw) = o_rho   (k,iw)
          theta_obsf(k,iw) = o_theta (k,iw)
            shw_obsf(k,iw) = o_shv   (k,iw)
         uzonal_obsf(k,iw) = o_uzonal(k,iw)
         umerid_obsf(k,iw) = o_umerid(k,iw)
      enddo

! If doing spectral nudging, sum data to nudging polygon arrays
! (Will need more work to do this in parallel!)

   else

      do k = kb,mza-1
             volwnud(k,iwnud1) =     volwnud(k,iwnud1) + volt(k,iw) 

            rho_obsf(k,iwnud1) =    rho_obsf(k,iwnud1) + o_rho   (k,iw) * volt(k,iw)
          theta_obsf(k,iwnud1) =  theta_obsf(k,iwnud1) + o_theta (k,iw) * volt(k,iw)
            shw_obsf(k,iwnud1) =    shw_obsf(k,iwnud1) + o_shv   (k,iw) * volt(k,iw)
         uzonal_obsf(k,iwnud1) = uzonal_obsf(k,iwnud1) + o_uzonal(k,iw) * volt(k,iw)
         umerid_obsf(k,iwnud1) = umerid_obsf(k,iwnud1) + o_umerid(k,iw) * volt(k,iw)
      enddo

   endif

enddo
call rsub('Wbb',7)

! If doing spectral nudging, normalize nudging point sums to get average values

if (nudnxp > 0) then

! Horizontal loop over nudging polygons

   do iwnud = 2,mwnud
      do k = 2,mza-1

                     volwnudi = 1. / max(1.,volwnud(k,iwnud))

            rho_obsf(k,iwnud) =    rho_obsf(k,iwnud) * volwnudi
          theta_obsf(k,iwnud) =  theta_obsf(k,iwnud) * volwnudi
            shw_obsf(k,iwnud) =    shw_obsf(k,iwnud) * volwnudi
         uzonal_obsf(k,iwnud) = uzonal_obsf(k,iwnud) * volwnudi
         umerid_obsf(k,iwnud) = umerid_obsf(k,iwnud) * volwnudi

      enddo
   enddo

endif

return
end subroutine nudge_prep

!===========================================================================

subroutine obs_nudge(rhot)

use mem_nudge, only: nudnxp, tnudcent, mwnud,                          &
                        rho_sim,    rho_obs,    rho_obsp,    rho_obsf, &
                      theta_sim,  theta_obs,  theta_obsp,  theta_obsf, &
                        shw_sim,    shw_obs,    shw_obsp,    shw_obsf, &
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

!$ use omp_lib

implicit none

! Nudge selected model fields (rho, thil, sh_w, umc, vmc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: iwnud,k,j,jv,iv,iw,iwnud1,iwnud2,iwnud3,iw1,iw2,iu,mrl,npoly,kb

real :: volwnud (mza,mwnud)
real :: umzonalt (mza,mwa)
real :: ummeridt (mza,mwa)
real :: uzonal(mza),umerid(mza)

real :: volwnudi,tp,tf,tnudi,tnudirho
real :: raxis,raxisi
real :: umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt
real :: fnud1,fnud2,fnud3

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK

integer, save :: icall = 0
real, save, allocatable :: wtnud(:)
real :: xw,yw,dist

if (icall /= 1) then
   icall = 1
   
   allocate (wtnud(mwa))
   
! Horizontal loop over T points

   do iw = 2,mwa

! Transform current IW point to polar stereographic coordinates using specified
! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
   
      call e_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)

      dist = sqrt(xw ** 2 + yw ** 2)
      
      if (dist > 5000.e3) then
         wtnud(iw) = 1.
      elseif (dist < 3600.e3) then
         wtnud(iw) = 0.
      else
         wtnud(iw) = ((dist - 3600.e3) / 1400.e3) ** 2
      endif

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

! If doing spectral nudging, zero out nudging polygon arrays and volume counter
! prior to summing

if (nudnxp > 0) then

   !$omp parallel do private (k)
   do iwnud = 2,mwnud
      do k = 2,mza-1
            rho_sim(k,iwnud) = 0.
          theta_sim(k,iwnud) = 0.
            shw_sim(k,iwnud) = 0.
         uzonal_sim(k,iwnud) = 0.
         umerid_sim(k,iwnud) = 0.
      
            volwnud(k,iwnud) = 0.
      enddo
   enddo
   !$omp end parallel do

endif

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------
call qsub('W',iw)

! Reconstruct UZONAL(k) and UMERID(k) from VXE, VYE, VZE

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb,mza-1
         uzonal(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
         umerid(k) = vze(k,iw) * raxis * eradi &
            - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi
      enddo

   else
      uzonal(:) = 0.
      umerid(:) = 0.
   endif

! If doing point-by-point (non-spectral) nudging, fill model value
! arrays with simple memory copy

   if (nudnxp == 0) then

      do k = 2,mza-1
            rho_sim(k,iw) = rho  (k,iw)
          theta_sim(k,iw) = theta(k,iw)
            shw_sim(k,iw) = sh_w (k,iw)
         uzonal_sim(k,iw) = uzonal(k)
         umerid_sim(k,iw) = umerid(k)
      enddo

! If doing spectral nudging, sum model fields and volume to nudging polygon arrays

   else

      do k = kb,mza-1
            volwnud(k,iwnud1) =    volwnud(k,iwnud1) + volt(k,iw) 

            rho_sim(k,iwnud1) =    rho_sim(k,iwnud1) + rho  (k,iw) * volt(k,iw)
          theta_sim(k,iwnud1) =  theta_sim(k,iwnud1) + theta(k,iw) * volt(k,iw)
            shw_sim(k,iwnud1) =    shw_sim(k,iwnud1) + sh_w (k,iw) * volt(k,iw)
         uzonal_sim(k,iwnud1) = uzonal_sim(k,iwnud1) + uzonal(k)   * volt(k,iw)
         umerid_sim(k,iwnud1) = umerid_sim(k,iwnud1) + umerid(k)   * volt(k,iw)
      enddo

   endif

enddo
call rsub('Wa',23)

! Horizontal loop over nudging polygons

!$omp parallel do private(k,volwnudi)
do iwnud = 2,mwnud

! If doing spectral nudging, normalize nudging point sums to get average values

   if (nudnxp > 0) then

! Vertical loop over nudging polygons

      do k = 2,mza-1

! Inverse volume of nudging cell

         volwnudi = 1. / max(1.,volwnud(k,iwnud))

! Normalize sum to obtain average model value for each polygon nudge point

            rho_sim(k,iwnud) =    rho_sim(k,iwnud) * volwnudi
          theta_sim(k,iwnud) =  theta_sim(k,iwnud) * volwnudi
            shw_sim(k,iwnud) =    shw_sim(k,iwnud) * volwnudi
         uzonal_sim(k,iwnud) = uzonal_sim(k,iwnud) * volwnudi
         umerid_sim(k,iwnud) = umerid_sim(k,iwnud) * volwnudi

      enddo
   endif

! Vertical loop over nudging polygons

      do k = 2,mza-1

! Interpolate observational fields in time

         rho_obs(k,iwnud) = tp *    rho_obsp(k,iwnud) + tf *    rho_obsf(k,iwnud)
       theta_obs(k,iwnud) = tp *  theta_obsp(k,iwnud) + tf *  theta_obsf(k,iwnud)
         shw_obs(k,iwnud) = tp *    shw_obsp(k,iwnud) + tf *    shw_obsf(k,iwnud)
      uzonal_obs(k,iwnud) = tp * uzonal_obsp(k,iwnud) + tf * uzonal_obsf(k,iwnud)
      umerid_obs(k,iwnud) = tp * umerid_obsp(k,iwnud) + tf * umerid_obsf(k,iwnud)

   enddo
enddo
!$omp end parallel do

! If doing point-by-point (non-spectral) nudging, compute tendencies using
! point-by-point information

if (nudnxp == 0) then

! Loop over all W columns

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iw,iwnud1,iwnud2,iwnud3,k,tnudi,fnud1,fnud2,fnud3)
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!---------------------------------------------------------------------
   call qsub('W',iw)

! Inverse of nudging time scale

! Default case - no spatial nudging mask

!      tnudi = 1. / tnudcent

! Use spatial nudging mask defined above in this subroutine

      tnudi = wtnud(iw) / tnudcent

      do k = lpw(iw),mza-1

         tnudirho = tnudi * rho(k,iw)

             rhot(k,iw) = rhot(k,iw) &
                        + tnudi * (rho_obs(k,iw) - rho_sim(k,iw))

            thilt(k,iw) = thilt(k,iw) &
                        + tnudirho * (theta_obs(k,iw) - theta_sim(k,iw))
                 
            sh_wt(k,iw) = sh_wt(k,iw) &
                        + tnudirho * (shw_obs(k,iw) - shw_sim(k,iw))

         umzonalt(k,iw) = tnudirho * (uzonal_obs(k,iw) - uzonal_sim(k,iw))

         ummeridt(k,iw) = tnudirho * (umerid_obs(k,iw) - umerid_sim(k,iw))

      enddo

   enddo
   !$omp end parallel do
   call rsub('Wb',23)

else

! Loop over all W columns, find 3 neighboring polygon points for each,
! and interpolate (obs - model) differences at each polygon point to the W point

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iw,iwnud1,iwnud2,iwnud3,k,tnudi,fnud1,fnud2,fnud3)
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
      iwnud1 = itab_w(iw)%iwnud(1);  fnud1 = itab_w(iw)%fnud(1)
      iwnud2 = itab_w(iw)%iwnud(2);  fnud2 = itab_w(iw)%fnud(2)
      iwnud3 = itab_w(iw)%iwnud(3);  fnud3 = itab_w(iw)%fnud(3)
!---------------------------------------------------------------------
   call qsub('W',iw)

! Inverse of nudging time scale

! Default case - no spatial nudging mask

!      tnudi = 1. / tnudcent

! Use spatial nudging mask defined above in this subroutine

      tnudi = wtnud(iw) / tnudcent

      do k = lpw(iw),mza-1

         tnudirho = tnudi * rho(k,iw)

          rhot(k,iw) = rhot(k,iw) + tnudi * ( &
                       fnud1 * (rho_obs(k,iwnud1) - rho_sim(k,iwnud1)) &
                     + fnud2 * (rho_obs(k,iwnud2) - rho_sim(k,iwnud2)) &
                     + fnud3 * (rho_obs(k,iwnud3) - rho_sim(k,iwnud3)) )

         thilt(k,iw) = thilt(k,iw) + tnudirho * ( &
                       fnud1 * (theta_obs(k,iwnud1) - theta_sim(k,iwnud1)) &
                     + fnud2 * (theta_obs(k,iwnud2) - theta_sim(k,iwnud2)) &
                     + fnud3 * (theta_obs(k,iwnud3) - theta_sim(k,iwnud3)) )
                 
         sh_wt(k,iw) = sh_wt(k,iw) + tnudirho * ( &
                       fnud1 * (shw_obs(k,iwnud1) - shw_sim(k,iwnud1)) &
                     + fnud2 * (shw_obs(k,iwnud2) - shw_sim(k,iwnud2)) &
                     + fnud3 * (shw_obs(k,iwnud3) - shw_sim(k,iwnud3)) )

      umzonalt(k,iw) = tnudirho * ( &
                       fnud1 * (uzonal_obs(k,iwnud1) - uzonal_sim(k,iwnud1)) &
                     + fnud2 * (uzonal_obs(k,iwnud2) - uzonal_sim(k,iwnud2)) &
                     + fnud3 * (uzonal_obs(k,iwnud3) - uzonal_sim(k,iwnud3)) )

      ummeridt(k,iw) = tnudirho * ( &
                       fnud1 * (umerid_obs(k,iwnud1) - umerid_sim(k,iwnud1)) &
                     + fnud2 * (umerid_obs(k,iwnud2) - umerid_sim(k,iwnud2)) &
                     + fnud3 * (umerid_obs(k,iwnud3) - umerid_sim(k,iwnud3)) )

      enddo

   enddo
   !$omp end parallel do
   call rsub('Wb',23)

endif

if (iparallel == 1) then
   call mpi_send_w('V', vxe=umzonalt, vye=ummeridt)
   call mpi_recv_w('V', vxe=umzonalt, vye=ummeridt)
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

         do k = lpu(iu),mza-1

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

         do k = lpv(iv),mza-1
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

return
end subroutine obs_nudge
