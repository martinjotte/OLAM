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
subroutine cuparm_driver(rhot)

use mem_grid,         only: mwa, mza, lpw, arw0, lpv
use module_cu_g3,     only: grell_driver
use module_cu_gf,     only: gf_driver
use module_cu_kfeta,  only: cuparm_kfeta, kf_lutab
use module_cu_tiedtke,only: cuparm_tiedtke
use module_cu_emanuel,only: cuparm_emanuel
use misc_coms,        only: io6, time_istp8, time_istp8p, nqparm, confrq, &
                            dtlong, initial, itime1, iparallel
use mem_ijtabs,       only: itab_w, jtab_w, mrl_begl, istp, mrls, jtw_prog, jtw_wadj
use mem_cuparm,       only: thsrc, rtsrc, aconpr, conprr, vxsrc, vysrc, vzsrc
use mem_tend,         only: thilt, sh_wt, vmxet, vmyet, vmzet
use mem_basic,        only: rho
use consts_coms,      only: r8
use olam_mpi_atm,     only: mpi_send_w, mpi_recv_w

implicit none

real, intent(inout) :: rhot(mza,mwa)
real(r8), parameter :: cptime = 0.0
integer             :: j, iw, k, mrl, mrlw
integer, save       :: init_kf = 0

real    :: dthmax, ftcon, dtlong4, confrq4, confrq4i
integer :: iwqmax, kqmax, km, km1

integer, save :: ncall = 0
integer :: npoly, jwn, iwn, iv, nblocked

real :: dx, ratio, total
real, save, allocatable :: tkeep(:), tsend(:)

real :: thsrc_distrib(mza,mwa)

real, parameter :: ths_min = 0.0
!real, parameter :: ths_min = 25.0 / 86400.

logical :: uvmix

!--------------------------------------------
! THSRC(k,iw) will contain the heating from parameterized convection in the IW
! vertical column.  We choose to apply only a fraction, tkeep, to the IW column
! itself and to divide the remainder equally among the nearest neighbors of IW,
! with each getting a fraction, tsend, of the total thsrc.  We first assign the
! ratio (tsend / tkeep) as a function of horizontal grid spacing dx.

if (ncall /= 1) then
   ncall = 1

   allocate (tkeep(mwa), tsend(mwa))

   do j = 1,jtab_w(jtw_wadj)%jend(1); iw = jtab_w(jtw_wadj)%iw(j) ! jend(1) for mrl = 1
      npoly = itab_w(iw)%npoly

      dx = sqrt(arw0(iw))

      if (dx < 100.e3) then
         ratio = 0.3 * dx / 100.e3  ! Linear increase from 0 to 100 km
      elseif (dx < 200.e3) then
         ratio = 0.3 + 0.2 * (dx - 100.e3) / 100.e3 ! linear increase from 100 to 200 km
      elseif (dx < 400.e3) then
         ratio = 0.5 - 0.5 * (dx - 200.e3) / 200.e3 ! linear decrease from 200 to 400 km
      else
         ratio = 0.  ! zero above 400 km
      endif

      total = 1. + ratio * real(npoly)

      tkeep(iw) = 1. / total
      tsend(iw) = ratio * tkeep(iw)

   enddo

endif
!--------------------------------------------

! For KF_eta parameterization, initialize scheme if needed

if ( any(nqparm(1:mrls) == 3) ) then
   if (init_kf == 0) then
      init_kf = 1
      call kf_lutab()
   endif
endif

uvmix = (allocated(vxsrc) .and. allocated(vysrc) .and. allocated(vzsrc))

! If model run has been initialized from observational data, avoid cumulus
! parameterization during an initial period to allow gravity waves to settle.

if ((initial == 2) .and. (time_istp8 < cptime)) return

! Check whether it is time to update cumulus parameterization tendencies

if ((istp == 1) .and. (mod(time_istp8p, confrq) < dtlong)) then

! Print message that cumulus parameterization is being computed

   write(io6, '(A,f0.2,A)') &
        ' Convective tendencies updated at ', time_istp8/3600., &
        ' hrs into simulation'

   dtlong4  = real(dtlong)
   confrq4  = real(confrq)
   confrq4i = 1. / confrq4

   thsrc_distrib(:,:) = 0.0
   thsrc        (:,:) = 0.0
   rtsrc        (:,:) = 0.0
   conprr         (:) = 0.0

   if (uvmix) then
      vxsrc(:,:) = 0.0
      vysrc(:,:) = 0.0
      vzsrc(:,:) = 0.0
   endif

! Loop over all IW grid cells where cumulus parameterization may be done

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(iw,mrlw,k,km,km1) 
   do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j) ! jend(1) for mrl = 1
!----------------------------------------------------------------------
   call qsub('W',iw)

! MRL for current IW column

      mrlw = itab_w(iw)%mrlw

! Select cumulus convection scheme based on MRL of current IW column
      
      if (nqparm(mrlw) == 1) then
   
! Tiedtke convective parameterization

         km = mza - lpw(iw) ! km  = # of T levels in cuparm_tiedtke
         km1 = km + 1       ! km1 = # of W levels in cuparm_tiedtke

         call cuparm_tiedtke(iw,km,km1,dtlong4,confrq4,confrq4i)

      elseif (nqparm(mrlw) == 2) then
   
! Grell deep and shallow convection

         call grell_driver( iw, dtlong4 )

      elseif (nqparm(mrlw) == 3) then
   
! Kain-Fritsch deep convection

         call cuparm_kfeta(iw,dtlong4)

      elseif (nqparm(mrlw) == 4) then
   
! Emanuel convective parameterization

         call cuparm_emanuel(iw,dtlong4)

      elseif (nqparm(mrlw) == 5) then
   
! Grell-Freitas deep and shallow convection

         call gf_driver( iw, dtlong4 )

      endif

! Distribute any deep convective heating among local and neighboring cells

      if (nqparm(mrlw) /= 0 .and. conprr(iw) > 1.e-16) then
         do k = lpw(iw), mza-1

            ! heating above ths_min will be a candidate to be distributed
            thsrc_distrib(k,iw) = max( thsrc(k,iw) - ths_min, 0.0 ) 

            ! the rest always stays in the local cell
            thsrc(k,iw) = thsrc(k,iw) - thsrc_distrib(k,iw)
         enddo
      endif

   enddo
   !$omp end parallel do
   call rsub('Wa',15)

! Perform MPI send of thsrc_distrib

   if (iparallel == 1) call mpi_send_w('V', vxe=thsrc_distrib, domrl=1)

! Print maximum heating rate in current parallel sub-domain

   dthmax = 0.
   iwqmax = 0
   kqmax  = 0

   ! this will need to be changed to work with OpenMP
   do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
      do k = lpw(iw), mza-1
         ftcon = thsrc(k,iw) + thsrc_distrib(k,iw)
         if (ftcon > dthmax) then
            dthmax = ftcon
           iwqmax = iw
            kqmax = k
         endif
      enddo
   enddo

   write(io6, '(A,I0,A,I0,A,F0.3,A,F0.4)') " MAX CONVECTIVE HEATING RATE AT IW=",  &
        iwqmax, " K=", kqmax, " IS ", dthmax*86400., " K/DAY"

! Perform MPI receive of thsrc_distrib

   if (iparallel == 1) call mpi_recv_w('V', vxe=thsrc_distrib, domrl=1)

! Distribute convective heating among neighboring cells

   !$omp parallel do private(iw,npoly,k,nblocked,jwn,iwn,iv) 
   do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

      npoly = itab_w(iw)%npoly

      do k = lpw(iw), mza-1
         nblocked = count( lpw( itab_w(iw)%iw(1:npoly)  ) > k )
         thsrc(k,iw) = thsrc(k,iw) + (tkeep(iw)+nblocked*tsend(iw)) * thsrc_distrib(k,iw)
      enddo

      do jwn = 1, npoly
         iwn = itab_w(iw)%iw(jwn)
         iv  = itab_w(iw)%iv(jwn)
         do k = lpv(iv), mza-1
            thsrc(k,iw) = thsrc(k,iw) + tsend(iwn) * thsrc_distrib(k,iwn)
         enddo
      enddo

   enddo
   !$omp end parallel do

endif

! Add current value of convective tendencies to thilt and sh_wt arrays
! every long timestep (whether cumulus parameterization is updated 
! this timestep or not)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   aconpr(iw) = aconpr(iw) + conprr(iw) * dtlong

   if ( any(nqparm(1:mrls) > 0) ) then

      do k = lpw(iw), mza-1
         thilt(k,iw) = thilt(k,iw) + thsrc(k,iw) * rho(k,iw)

         sh_wt(k,iw) = sh_wt(k,iw) + rtsrc(k,iw) * rho(k,iw)

         rhot (k,iw) = rhot (k,iw) + rtsrc(k,iw) * rho(k,iw)
      enddo

      if (uvmix) then
         do k = lpw(iw), mza-1
            vmxet(k,iw) = vmxet(k,iw) + vxsrc(k,iw) * rho(k,iw)

            vmyet(k,iw) = vmyet(k,iw) + vysrc(k,iw) * rho(k,iw)

            vmzet(k,iw) = vmzet(k,iw) + vzsrc(k,iw) * rho(k,iw)
         enddo
      endif

   endif

enddo
!$omp end parallel do
endif
call rsub('Wb',15)

return
end subroutine cuparm_driver
