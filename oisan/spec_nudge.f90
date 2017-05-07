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
subroutine nudge_prep_spec(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid)

use mem_nudge,   only: mwnud, &
                       rho_obsp, theta_obsp, shw_obsp, &
                       rho_obsf, theta_obsf, shw_obsf, &
                       uzonal_obsp, umerid_obsp, &
                       uzonal_obsf, umerid_obsf

use mem_grid,    only: mza, mwa, volt
use misc_coms,   only: iparallel
use mem_ijtabs,  only: jtab_w, itab_w, jtw_init
use olam_mpi_atm,only: mpi_send_wnud, mpi_recv_wnud
use consts_coms, only: r8

implicit none

integer, intent(in) :: iaction

real(r8), intent(in) :: o_rho   (mza,mwa)
real,     intent(in) :: o_theta (mza,mwa)
real,     intent(in) :: o_shv   (mza,mwa)
real,     intent(in) :: o_uzonal(mza,mwa)
real,     intent(in) :: o_umerid(mza,mwa)

integer :: j,iw,k,iwnud,iwnud1
real    :: volwnudi

! Automatic arrays

real :: volwnud(mza,mwnud)

! Swap future data time into past data time if necessary.

if (iaction == 1) then

   do iwnud = 2, mwnud
      do k = 2, mza
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

do iwnud = 2, mwnud
   do k = 2, mza
         rho_obsf(k,iwnud) = 0.
       theta_obsf(k,iwnud) = 0.
         shw_obsf(k,iwnud) = 0.
      uzonal_obsf(k,iwnud) = 0.
      umerid_obsf(k,iwnud) = 0.

          volwnud(k,iwnud) = 0.
   enddo
enddo

! If doing spectral nudging, sum data to nudging polygon arrays

!----------------------------------------------------------------------
do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------

   do k = 2, mza
          volwnud(k,iwnud1) =     volwnud(k,iwnud1) + volt(k,iw) 
         rho_obsf(k,iwnud1) =    rho_obsf(k,iwnud1) + o_rho   (k,iw) * volt(k,iw)
       theta_obsf(k,iwnud1) =  theta_obsf(k,iwnud1) + o_theta (k,iw) * volt(k,iw)
         shw_obsf(k,iwnud1) =    shw_obsf(k,iwnud1) + o_shv   (k,iw) * volt(k,iw)
      uzonal_obsf(k,iwnud1) = uzonal_obsf(k,iwnud1) + o_uzonal(k,iw) * volt(k,iw)
      umerid_obsf(k,iwnud1) = umerid_obsf(k,iwnud1) + o_umerid(k,iw) * volt(k,iw)
   enddo

enddo

! MPI SEND/RECV of nudging arrays

if (iparallel == 1) call mpi_send_wnud(rvara1=volwnud,    rvara2=rho_obsf, &
                                       rvara3=theta_obsf, rvara4=shw_obsf, &
                                       rvara5=uzonal_obsf,rvara6=umerid_obsf)

if (iparallel == 1) call mpi_recv_wnud(rvara1=volwnud,    rvara2=rho_obsf, &
                                       rvara3=theta_obsf, rvara4=shw_obsf, &
                                       rvara5=uzonal_obsf,rvara6=umerid_obsf)

! If doing spectral nudging, normalize nudging point sums to get average values

! Horizontal loop over nudging polygons

!$omp parallel do private(k,volwnudi)
do iwnud = 2, mwnud
   do k = 2, mza

                  volwnudi = 1. / max(1.,volwnud(k,iwnud))

         rho_obsf(k,iwnud) =    rho_obsf(k,iwnud) * volwnudi
       theta_obsf(k,iwnud) =  theta_obsf(k,iwnud) * volwnudi
         shw_obsf(k,iwnud) =    shw_obsf(k,iwnud) * volwnudi
      uzonal_obsf(k,iwnud) = uzonal_obsf(k,iwnud) * volwnudi
      umerid_obsf(k,iwnud) = umerid_obsf(k,iwnud) * volwnudi

   enddo
enddo
!$omp end parallel do

end subroutine nudge_prep_spec

!===========================================================================

subroutine spec_nudge(rhot)

use mem_nudge, only:   tnudcent,      mwnud,                           &
                        rho_sim,    rho_obs,    rho_obsp,    rho_obsf, &
                      theta_sim,  theta_obs,  theta_obsp,  theta_obsf, &
                        shw_sim,    shw_obs,    shw_obsp,    shw_obsf, &
                     uzonal_sim, uzonal_obs, uzonal_obsp, uzonal_obsf, &
                     umerid_sim, umerid_obs, umerid_obsp, umerid_obsf

use mem_basic,   only: rho, theta, sh_w, vxe, vye, vze
use mem_grid,    only: mza, mwa, lpw, xew, yew, zew, volt
use misc_coms,   only: s1900_sim, iparallel, dtlm
use mem_ijtabs,  only: istp, jtab_w, itab_w, mrl_begl, jtv_prog, jtw_prog
use consts_coms, only: eradi
use mem_tend,    only: thilt, sh_wt, vmxet, vmyet, vmzet
use isan_coms,   only: ifgfile, s1900_fg
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_wnud, mpi_recv_wnud

implicit none

! Nudge selected model fields (rho, thil, sh_w, vmc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: iwnud,k,j,iw,iwnud1,iwnud2,iwnud3,mrl,kb

real :: volwnud (mza,mwnud)
real :: umzonalt(mza)
real :: ummeridt(mza)
real :: uzonal  (mza)
real :: umerid  (mza)

real :: volwnudi,tp,tf,tnudi,tnudirho
real :: raxis,raxisi,uvtr
real :: fnud1,fnud2,fnud3

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK
!
!integer, save :: icall = 0
!real, save, allocatable :: wtnud(:)
!real :: xw,yw,dist
!
!if (icall /= 1) then
!   icall = 1
!   
!   allocate (wtnud(mwa))
!   
!! Horizontal loop over T points
!
!   do iw = 2,mwa
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

! If doing spectral nudging, zero out nudging polygon arrays and volume counter
! prior to summing

!$omp parallel do private (k)
do iwnud = 2,mwnud
   do k = 2,mza
         rho_sim(k,iwnud) = 0.
       theta_sim(k,iwnud) = 0.
         shw_sim(k,iwnud) = 0.
      uzonal_sim(k,iwnud) = 0.
      umerid_sim(k,iwnud) = 0.
      
         volwnud(k,iwnud) = 0.
   enddo
enddo
!$omp end parallel do

! Horizontal loop over W columns

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------

! Reconstruct UZONAL(k) and UMERID(k) from VXE, VYE, VZE

   kb = lpw(iw)
   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb, mza

         uzonal(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
         umerid(k) = vze(k,iw) * raxis * eradi &
            - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi

      enddo

   else
      uzonal(:) = 0.
      umerid(:) = 0.
   endif

! If doing spectral nudging, sum model fields and volume to nudging polygon arrays
   do k = kb, mza
         volwnud(k,iwnud1) =    volwnud(k,iwnud1) + volt(k,iw) 

         rho_sim(k,iwnud1) =    rho_sim(k,iwnud1) + rho  (k,iw) * volt(k,iw)
       theta_sim(k,iwnud1) =  theta_sim(k,iwnud1) + theta(k,iw) * volt(k,iw)
         shw_sim(k,iwnud1) =    shw_sim(k,iwnud1) + sh_w (k,iw) * volt(k,iw)
      uzonal_sim(k,iwnud1) = uzonal_sim(k,iwnud1) + uzonal(k)   * volt(k,iw)
      umerid_sim(k,iwnud1) = umerid_sim(k,iwnud1) + umerid(k)   * volt(k,iw)
   enddo

enddo

! MPI SEND/RECV of nudging arrays

if (iparallel == 1) call mpi_send_wnud(rvara1=volwnud,   rvara2=rho_sim, &
                                       rvara3=theta_sim, rvara4=shw_sim, &
                                       rvara5=uzonal_sim,rvara6=umerid_sim)

if (iparallel == 1) call mpi_recv_wnud(rvara1=volwnud,   rvara2=rho_sim, &
                                       rvara3=theta_sim, rvara4=shw_sim, &
                                       rvara5=uzonal_sim,rvara6=umerid_sim)

! Horizontal loop over nudging polygons

!$omp parallel do private(k,volwnudi)
do iwnud = 2,mwnud

! If doing spectral nudging, normalize nudging point sums to get average values

! Vertical loop over nudging polygons

   do k = 2,mza

! Inverse volume of nudging cell

                 volwnudi = 1. / max(1.,volwnud(k,iwnud))

         rho_sim(k,iwnud) =    rho_sim(k,iwnud) * volwnudi
       theta_sim(k,iwnud) =  theta_sim(k,iwnud) * volwnudi
         shw_sim(k,iwnud) =    shw_sim(k,iwnud) * volwnudi
      uzonal_sim(k,iwnud) = uzonal_sim(k,iwnud) * volwnudi
      umerid_sim(k,iwnud) = umerid_sim(k,iwnud) * volwnudi

   enddo

! Vertical loop over nudging polygons

   do k = 2, mza

! Interpolate observational fields in time

         rho_obs(k,iwnud) = tp *    rho_obsp(k,iwnud) + tf *    rho_obsf(k,iwnud)
       theta_obs(k,iwnud) = tp *  theta_obsp(k,iwnud) + tf *  theta_obsf(k,iwnud)
         shw_obs(k,iwnud) = tp *    shw_obsp(k,iwnud) + tf *    shw_obsf(k,iwnud)
      uzonal_obs(k,iwnud) = tp * uzonal_obsp(k,iwnud) + tf * uzonal_obsf(k,iwnud)
      umerid_obs(k,iwnud) = tp * umerid_obsp(k,iwnud) + tf * umerid_obsf(k,iwnud)
   enddo
      
enddo
!$omp end parallel do

! Loop over all W columns, find 3 neighboring polygon points for each,
! and interpolate (obs - model) differences at each polygon point to the W point

!----------------------------------------------------------------------
!$omp parallel do private(iw,iwnud1,iwnud2,iwnud3,k,tnudi,tnudirho, &
!$omp                     fnud1,fnud2,fnud3,umzonalt,ummeridt,uvtr, &
!$omp                     raxis,raxisi)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1);  fnud1 = itab_w(iw)%fnud(1)
   iwnud2 = itab_w(iw)%iwnud(2);  fnud2 = itab_w(iw)%fnud(2)
   iwnud3 = itab_w(iw)%iwnud(3);  fnud3 = itab_w(iw)%fnud(3)
!---------------------------------------------------------------------

! Inverse of nudging time scale
! Use spatial nudging mask defined above in this subroutine

!  tnudi = wtnud(iw) / tnudcent
   tnudi = 1.0 / tnudcent

   do k = lpw(iw), mza

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

      umzonalt(k)    = tnudirho * ( &
                       fnud1 * (uzonal_obs(k,iwnud1) - uzonal_sim(k,iwnud1)) &
                     + fnud2 * (uzonal_obs(k,iwnud2) - uzonal_sim(k,iwnud2)) &
                     + fnud3 * (uzonal_obs(k,iwnud3) - uzonal_sim(k,iwnud3)) )

      ummeridt(k)    = tnudirho * ( &
                       fnud1 * (umerid_obs(k,iwnud1) - umerid_sim(k,iwnud1)) &
                     + fnud2 * (umerid_obs(k,iwnud2) - umerid_sim(k,iwnud2)) &
                     + fnud3 * (umerid_obs(k,iwnud3) - umerid_sim(k,iwnud3)) )

      ! With nudging, scalar mass is no longer conserved. Make sure the long timestep
      ! tendency does not drive humidity negative when nudging is included:

      sh_wt(k,iw) = max( sh_wt(k,iw), &
           -0.999 * max(sh_w(k,iw), 0.0) * real( rho(k,iw) / dtlm(itab_w(iw)%mrlw) ) )

   enddo

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = lpw(iw), mza
         uvtr = -ummeridt(k) * zew(iw) * eradi
         vmxet(k,iw) = vmxet(k,iw) + (-umzonalt(k) * yew(iw) + uvtr * xew(iw)) * raxisi
         vmyet(k,iw) = vmyet(k,iw) + ( umzonalt(k) * xew(iw) + uvtr * yew(iw)) * raxisi
         vmzet(k,iw) = vmzet(k,iw) +   ummeridt(k) * raxis * eradi 
      enddo
   endif

enddo
!$omp end parallel do

end subroutine spec_nudge
