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
subroutine fldsisan(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid, o_uvc)

use mem_nudge,   only: nudflag, nudnxp, mwnud, &
                       rho_obsp, theta_obsp, shw_obsp, uzonal_obsp, &
                       umerid_obsp, rho_obsf, theta_obsf, shw_obsf, &
                       uzonal_obsf, umerid_obsf
use mem_basic,   only: umc, ump, uc, vmc, vmp, vc, thil, sh_w, sh_v, wmc, wc, &
                       theta, rho, press
use mem_grid,    only: mza, mua, mva, mwa, lcu, lcv, lpw, zm, zt, &
                       unx, uny, unz, vnx, vny, vnz, xeu, yeu, zeu, &
                       xev, yev, zev, xew, yew, zew, aru, arv, volt
use misc_coms,   only: io6, deltax, iparallel, runtype, &
                       meshtype
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_u, jtab_v, jtab_w, itab_u, itab_v, itab_w
use consts_coms, only: pc1, rdry, rvap, cpocv, erad, eradi
use massflux,    only: diagnose_uc, diagnose_vc

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                        mpi_send_u, mpi_recv_u, &
                        mpi_send_v, mpi_recv_v 

implicit none

integer, intent(in)  :: iaction

real(kind=8), intent(in) :: o_rho   (mza,mwa)
real,         intent(in) :: o_theta (mza,mwa)
real,         intent(in) :: o_shv   (mza,mwa)
real,         intent(in) :: o_uzonal(mza,mwa)
real,         intent(in) :: o_umerid(mza,mwa)
real,         intent(in) :: o_uvc   (mza,mva)

integer :: j,iw,k,ka,iu,iv,iw1,iw2,iup,ivp,iwnud,iwnud1
integer :: npoly,kb,jv

real :: rcloud,temp,rvls,uu,vv,angle
real :: alph_p,farv2,raxis,raxisi
real :: volwnudi

! Automatic arrays

real :: volwnud(mza,mwnud)
real :: ouzonal(mza),oumerid(mza)
real :: vxe(mza),vye(mza),vze(mza)

! init_hurricane is reset in subroutine hurricane, if calls are uncommented below

integer :: init_hurricane = 0

! Check whether initializing model

if (iaction == 0 .and. runtype == 'INITIAL') then

! If initializing the model, fill the main model arrays
! and initialize related arrays

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
!---------------------------------------------------------------------
   call qsub('W',iw)
      do k = 1,mza

         rho  (k,iw) = o_rho(k,iw)
         theta(k,iw) = o_theta(k,iw)
         sh_v (k,iw) = o_shv(k,iw)

         thil(k,iw) = theta(k,iw)
         sh_w(k,iw) = sh_v(k,iw)   ! no condensate considered here
         if (level > 1) then
            sh_c(k,iw) = 0.        ! no condensate considered here
         endif

! alph_p is like alpha_press except that the (theta/thil) factor is excluded
! here (it's equal to 1)

         alph_p = pc1 * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)) &
                ** cpocv

         press(k,iw) = alph_p * (rho(k,iw) * thil(k,iw)) ** cpocv

         wc(k,iw) = 0.
         wmc(k,iw) = 0.

      enddo
   enddo
   call rsub('Wb',7)

!-------------------------------------------------------------------------------
!  call hurricane_init('T',init_hurricane) ! returns init_hurricane = 0, 1, or 2
!-------------------------------------------------------------------------------

   if (iparallel == 1) then
      call mpi_send_w('I')  ! Send W group
      call mpi_recv_w('I')  ! Recv W group
   endif

   if (meshtype == 1) then

! If triangular mesh, initialize UMC, UC

      call psub()
!----------------------------------------------------------------------
      do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)
         iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
      call qsub('U',iu)

         do k = 1,mza
            uc(k,iu) = o_uvc(k,iu)
            umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
         enddo

      enddo
   
      call rsub('Ub',7)

! Horizontal loop over U points

      call psub()
!----------------------------------------------------------------------
      do iu = 2,mua
!----------------------------------------------------------------------
      call qsub('U',iu) !QQQQQ

! For below-ground points, set UC to LCU value.

         ka = lcu(iu)
         uc(1:ka-1,iu) = uc(ka,iu)

      enddo
      call rsub('Uc',0)

! MPI parallel send/recv of U group

      if (iparallel == 1) then
         call mpi_send_u('I')
         call mpi_recv_u('I')
      endif

! Set UMP to UMC

      ump(:,:) = umc(:,:)

! Diagnose VC
   
      call diagnose_vc()

!-------------------------------------------------------------------------------
!  call hurricane_init('U',init_hurricane) ! returns init_hurricane = 0, 1, or 2

! If initial wind has been altered by hurricane initialization, repeat
! parallel send/receive of wind field and repeat VC diagnosis

      if (init_hurricane > 0) then

         if (iparallel == 1) then
            call mpi_send_u('I')
            call mpi_recv_u('I')
         endif

         ump(:,:) = umc(:,:)

         call diagnose_vc()

      endif
!-------------------------------------------------------------------------------

   else

! If using hexagonal mesh, initialize VMC, VC

      call psub()
!----------------------------------------------------------------------
      do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)
         iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
      call qsub('V',iv)

         do k = 1,mza
            vc(k,iv) = o_uvc(k,iv)
            vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
         enddo

      enddo
   
      call rsub('Vb',7)

! Horizontal loop over V points

      call psub()
!----------------------------------------------------------------------
      do iv = 2,mva
!----------------------------------------------------------------------
      call qsub('V',iv)

! For below-ground points, set VC to LCV value.

         ka = lcv(iv)
         vc(1:ka-1,iv) = vc(ka,iv)

      enddo
      call rsub('Vc',0)
   
! MPI parallel send/recv of V group

      if (iparallel == 1) then
         call mpi_send_v('I')
         call mpi_recv_v('I')
      endif

! Set VMP to VMC

      vmp(:,:) = vmc(:,:)

! Diagnose UC
   
      call diagnose_uc()

!-------------------------------------------------------------------------------
!  call hurricane_init('V',init_hurricane) ! returns init_hurricane = 0, 1, or 2

      if (init_hurricane > 0) then

! If initial wind has been altered by hurricane initialization, repeat
! parallel send/receive of wind field and UC diagnosis

         if (iparallel == 1) then
            call mpi_send_v('I')
            call mpi_recv_v('I')
         endif

         vmp(:,:) = vmc(:,:)

         call diagnose_uc()

      endif
!-------------------------------------------------------------------------------

   endif

! Print out initial state column from column 2

   iw = 2

write(io6,*)' '
write(io6,*)'========================================================================'
write(io6,*)'                    OLAM INITIAL STATE COLUMN (vari)'
write(io6,*)'========================================================================'
write(io6,*)'   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*)'========================================================================'
write(io6,*)' '

   do k = mza-1,2,-1
      write(io6, '(f10.2,1x,9(''-------''))') zm(k)
      write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)') &
          k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
   enddo

   write(io6, '(f10.2,1x,9(''-------''))') zm(1)
   write(io6,*)' '

endif

! Return if nudging flag is not activated

if (nudflag < 1 .or. nudnxp < 1) return

! If we got here, nudging will be done in this model run, so fill nudging arrays.

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

! Zero out nudging polygon arrays and volume counter prior to summing

do iwnud = 2,mwnud
   do k = 2,mza-1
          volwnud(k,iwnud) = 0.

         rho_obsf(k,iwnud) = 0.
       theta_obsf(k,iwnud) = 0.
         shw_obsf(k,iwnud) = 0.
      uzonal_obsf(k,iwnud) = 0.
      umerid_obsf(k,iwnud) = 0.
   enddo
enddo

! Sum data to nudging polygon arrays

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------
call qsub('W',iw)

! Reconstruct OUZONAL(k) and OUMERID(k) from O_UVC instead of using 
! O_UZONAL(k,iw) and O_UMERID(k,iw) that were interpolated directly from obs

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)
   
   vxe(:) = 0.
   vye(:) = 0.
   vze(:) = 0.

   do jv = 1,npoly

      if (meshtype == 1) then

         iv = itab_w(iw)%iu(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + itab_w(iw)%vxu(jv) * o_uvc(k,iv)
            vye(k) = vye(k) + itab_w(iw)%vyu(jv) * o_uvc(k,iv)
            vze(k) = vze(k) + itab_w(iw)%vzu(jv) * o_uvc(k,iv)
         enddo

      else

         iv = itab_w(iw)%iv(jv)
         farv2 = 2. * itab_w(iw)%farv(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + farv2 * o_uvc(k,iv) * vnx(iv)
            vye(k) = vye(k) + farv2 * o_uvc(k,iv) * vny(iv)
            vze(k) = vze(k) + farv2 * o_uvc(k,iv) * vnz(iv)
         enddo

      endif

   enddo

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb,mza-1
         ouzonal(k) = (vye(k) * xew(iw) - vxe(k) * yew(iw)) * raxisi
         oumerid(k) = vze(k) * raxis * eradi &
            - (vxe(k) * xew(iw) + vye(k) * yew(iw)) * zew(iw) * raxisi * eradi
      enddo

   else
      ouzonal(:) = 0.
      oumerid(:) = 0.
   endif

! With observed values now interpolated to model T point, add these values to
! sum at local nudging point

   do k = kb,mza-1
          volwnud(k,iwnud1) =     volwnud(k,iwnud1) + volt(k,iw) 

         rho_obsf(k,iwnud1) =    rho_obsf(k,iwnud1) + o_rho  (k,iw) * volt(k,iw)
       theta_obsf(k,iwnud1) =  theta_obsf(k,iwnud1) + o_theta(k,iw) * volt(k,iw)
         shw_obsf(k,iwnud1) =    shw_obsf(k,iwnud1) + o_shv  (k,iw) * volt(k,iw)
      uzonal_obsf(k,iwnud1) = uzonal_obsf(k,iwnud1) + ouzonal(k)    * volt(k,iw)
      umerid_obsf(k,iwnud1) = umerid_obsf(k,iwnud1) + oumerid(k)    * volt(k,iw)
   enddo

enddo
call rsub('Wbb',7)

! Normalize nudging point sums to get average values

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

return
end subroutine fldsisan

!===========================================================================

subroutine obs_nudge(rhot)

use mem_nudge,   only: tnudcent, mwnud,                                  &
                       rho_sim, rho_obs, rho_obsp, rho_obsf,             &
                       theta_sim, theta_obs, theta_obsp, theta_obsf,     &
                       uzonal_sim, uzonal_obs, uzonal_obsp, uzonal_obsf, &
                       umerid_sim, umerid_obs, umerid_obsp, umerid_obsf, &
                       shw_sim, shw_obs, shw_obsp, shw_obsf
use mem_basic,   only: uc, vc, rho, theta, sh_w
use mem_grid,    only: mza, mwa, lpu, lpv, lpw, &
                       xeu, yeu, zeu, xev, yev, zev, xew, yew, zew, &
                       unx, uny, unz, vnx, vny, vnz, volt, glatw, glonw
use misc_coms,   only: io6, time8, s1900_sim, meshtype
use mem_ijtabs,  only: istp, jtab_u, jtab_v, jtab_w, itab_u, itab_v, itab_w, &
                       mrl_begl
use consts_coms, only: erad, eradi
use mem_tend,    only: umt, vmt, thilt, sh_wt
use isan_coms,   only: ifgfile, s1900_fg

!$ use omp_lib

implicit none

! Nudge selected model fields (rho, thil, sh_w, umc, vmc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: iwnud,k,j,jv,iv,iw,iwnud1,iwnud2,iwnud3,iw1,iw2,iu,mrl,npoly,kb

real :: volwnud (mza,mwnud)
real :: umzonalt (mza,mwa)
real :: ummeridt (mza,mwa)
real :: vxe(mza),vye(mza),vze(mza)
real :: uzonal(mza),umerid(mza)

real :: volwnudi,tp,tf,tnudi
real :: raxis,raxisi
real :: umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt,farv2
real :: fnud1,fnud2,fnud3

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK

!integer, save :: icall = 0
!real, save, allocatable :: wtnud(:)
!real :: xw,yw,dist

!if (icall /= 1) then
!   icall = 1
   
!   allocate (wtnud(mwa))
   
! Horizontal loop over T points

!   do iw = 2,mwa

! Transform current IW point to polar stereographic coordinates using specified
! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
   
!      call e_ps(xew(iw),yew(iw),zew(iw),36.,-120.,xw,yw)

!      dist = sqrt(xw ** 2 + yw ** 2)
      
!      if (dist > 4000.e3) then
!         wtnud(iw) = 1.
!      elseif (dist < 3600.e3) then
!         wtnud(iw) = 0.
!      else
!         wtnud(iw) = (dist - 3600.e3) / 400.e3
!      endif

!   enddo

!endif
!----------------------------------------------------------------------

! Check whether it is time to nudge

mrl = mrl_begl(istp)
if (mrl < 1) return

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf

! Zero out polygon arrays for model fields and volume counter prior to summing

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

! Sum model values to polygon arrays

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(23)%jend(mrl); iw = jtab_w(23)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1)
!---------------------------------------------------------------------
call qsub('W',iw)

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

   vxe(:) = 0.
   vye(:) = 0.
   vze(:) = 0.

   do jv = 1,npoly

      if (meshtype == 1) then
      
         iv = itab_w(iw)%iu(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + itab_w(iw)%vxu(jv) * uc(k,iv)
            vye(k) = vye(k) + itab_w(iw)%vyu(jv) * uc(k,iv)
            vze(k) = vze(k) + itab_w(iw)%vzu(jv) * uc(k,iv)
         enddo

      else

         iv = itab_w(iw)%iv(jv)
         farv2 = 2. * itab_w(iw)%farv(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + farv2 * vc(k,iv) * vnx(iv)
            vye(k) = vye(k) + farv2 * vc(k,iv) * vny(iv)
            vze(k) = vze(k) + farv2 * vc(k,iv) * vnz(iv)
         enddo

      endif

   enddo

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb,mza-1
         uzonal(k) = (vye(k) * xew(iw) - vxe(k) * yew(iw)) * raxisi
         umerid(k) = vze(k) * raxis * eradi &
            - (vxe(k) * xew(iw) + vye(k) * yew(iw)) * zew(iw) * raxisi * eradi
      enddo

   else
      uzonal(:) = 0.
      umerid(:) = 0.
   endif

! Sum model fields and volume to polygon arrays

   do k = kb,mza-1
         volwnud(k,iwnud1) =    volwnud(k,iwnud1) + volt(k,iw) 

         rho_sim(k,iwnud1) =    rho_sim(k,iwnud1) + rho  (k,iw) * volt(k,iw)
       theta_sim(k,iwnud1) =  theta_sim(k,iwnud1) + theta(k,iw) * volt(k,iw)
         shw_sim(k,iwnud1) =    shw_sim(k,iwnud1) + sh_w (k,iw) * volt(k,iw)
      uzonal_sim(k,iwnud1) = uzonal_sim(k,iwnud1) + uzonal(k)   * volt(k,iw)
      umerid_sim(k,iwnud1) = umerid_sim(k,iwnud1) + umerid(k)   * volt(k,iw)
   enddo

enddo
call rsub('Wa',23)

! Horizontal loop over nudging polygons

!$omp parallel do private(k,volwnudi)
do iwnud = 2,mwnud

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

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,iwnud1,iwnud2,iwnud3,k,tnudi,fnud1,fnud2,fnud3)
do j = 1,jtab_w(23)%jend(mrl); iw = jtab_w(23)%iw(j)
   iwnud1 = itab_w(iw)%iwnud(1);  fnud1 = itab_w(iw)%fnud(1)
   iwnud2 = itab_w(iw)%iwnud(2);  fnud2 = itab_w(iw)%fnud(2)
   iwnud3 = itab_w(iw)%iwnud(3);  fnud3 = itab_w(iw)%fnud(3)
!---------------------------------------------------------------------
call qsub('W',iw)

! Inverse of nudging time scale

! Default case - no spatial nudging mask

   tnudi = 1. / tnudcent

! Use spatial nudging mask defined above in this subroutine

!   tnudi = wtnud(iw) / tnudcent

   do k = lpw(iw),mza-1

       rhot(k,iw) = rhot(k,iw) + tnudi * ( &
                    fnud1 * (rho_obs(k,iwnud1) - rho_sim(k,iwnud1)) &
                  + fnud2 * (rho_obs(k,iwnud2) - rho_sim(k,iwnud2)) &
                  + fnud3 * (rho_obs(k,iwnud3) - rho_sim(k,iwnud3)) )

      thilt(k,iw) = thilt(k,iw) + tnudi * rho(k,iw) * ( &
                    fnud1 * (theta_obs(k,iwnud1) - theta_sim(k,iwnud1)) &
                  + fnud2 * (theta_obs(k,iwnud2) - theta_sim(k,iwnud2)) &
                  + fnud3 * (theta_obs(k,iwnud3) - theta_sim(k,iwnud3)) )
                 
      sh_wt(k,iw) = sh_wt(k,iw) + tnudi * rho(k,iw) * ( &
                    fnud1 * (shw_obs(k,iwnud1) - shw_sim(k,iwnud1)) &
                  + fnud2 * (shw_obs(k,iwnud2) - shw_sim(k,iwnud2)) &
                  + fnud3 * (shw_obs(k,iwnud3) - shw_sim(k,iwnud3)) )

   umzonalt(k,iw) = tnudi * rho(k,iw) * ( &
                    fnud1 * (uzonal_obs(k,iwnud1) - uzonal_sim(k,iwnud1)) &
                  + fnud2 * (uzonal_obs(k,iwnud2) - uzonal_sim(k,iwnud2)) &
                  + fnud3 * (uzonal_obs(k,iwnud3) - uzonal_sim(k,iwnud3)) )

   ummeridt(k,iw) = tnudi * rho(k,iw) * ( &
                    fnud1 * (umerid_obs(k,iwnud1) - umerid_sim(k,iwnud1)) &
                  + fnud2 * (umerid_obs(k,iwnud2) - umerid_sim(k,iwnud2)) &
                  + fnud3 * (umerid_obs(k,iwnud3) - umerid_sim(k,iwnud3)) )

   enddo

enddo
!$omp end parallel do
call rsub('Wb',23)

if (meshtype == 1) then

! UMT

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iu,iw1,iw2,k, &
   !$omp                     raxis,raxisi,umgt,vmgt,uvmgrt,uvmgxt,uvmgyt,uvmgzt)
   do j = 1,jtab_u(13)%jend(mrl); iu = jtab_u(13)%iu(j)
      iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      if (iw1 < 2) iw1 = iw2
      if (iw2 < 2) iw2 = iw1

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
   do j = 1,jtab_v(13)%jend(mrl); iv = jtab_v(13)%iv(j)
      iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   call qsub('V',iv)

      if (iw1 < 2) iw1 = iw2
      if (iw2 < 2) iw2 = iw1

      raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
      
         raxisi = 1. / raxis

! Average momentum tendencies to V point and rotate at V point

         do k = lpv(iv),mza-1
            umgt = .5 * (umzonalt(k,iw1) + umzonalt(k,iw2))

            vmgt = .5 * (ummeridt(k,iw1) + ummeridt(k,iw2))

            uvmgrt = -vmgt * zev(iv) / erad  ! radially outward from axis

            uvmgxt = (-umgt * yev(iv) + uvmgrt * xev(iv)) * raxisi 
            uvmgyt = ( umgt * xev(iv) + uvmgrt * yev(iv)) * raxisi 
            uvmgzt =   vmgt * raxis / erad 

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
