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
subroutine sea_init_atm()

use mem_sea,     only: sea, itabg_ws, itab_ws
use sea_coms,    only: mws, iupdsst, s1900_sst, isstfile, nsstfiles, dt_sea,  &
                       iupdseaice, s1900_seaice, iseaicefile, nseaicefiles
                      
use mem_sflux,   only: mseaflux, seaflux, jseaflux
use mem_basic,   only: press, rho, theta, sh_v
use misc_coms,   only: io6, time8, s1900_sim, iparallel, runtype
use mem_ijtabs,  only: itabg_w
use mem_para,    only: myrank
                       
use consts_coms, only: p00i, rocp

implicit none

integer :: isf
integer :: iw
integer :: kw
integer :: iws
integer :: mrl
integer :: j

real :: airtemp
real :: timefac_sst
real :: timefac_seaice
real :: arf_sea

real, external :: rhovsl

! This subroutine fills the primary LEAF3 arrays which depend on current
! atmospheric conditions.

! Time interpolation factor for updating SST and SEA ICE

timefac_sst = 0.
timefac_seaice = 0.

if (iupdsst == 1 .and. nsstfiles > 1) then
   timefac_sst = (s1900_sim             - s1900_sst(isstfile))  &
               / (s1900_sst(isstfile+1) - s1900_sst(isstfile))
endif

if (iupdseaice == 1 .and. nseaicefiles > 1) then
   timefac_seaice = (s1900_sim                   - s1900_seaice(iseaicefile)) &
                  / (s1900_seaice(iseaicefile+1) - s1900_seaice(iseaicefile))
endif

if (runtype /= "INITIAL") return

! Loop over all SEAFLUX cells that are EVALUATED on this rank, and transfer 
! atmospheric properties to each, with weighting according to arf_sea

do j = 1,jseaflux(1)%jend(1)
   isf = jseaflux(1)%iseaflux(j)

   iw = seaflux(isf)%iw  ! global index

! If run is parallel, convert iw to local domain

   if (iparallel == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif

   kw      = seaflux(isf)%kw
   arf_sea = seaflux(isf)%arf_sfc

   airtemp = theta(kw,iw) * (p00i * press(kw,iw)) ** rocp

   seaflux(isf)%rhos    = arf_sea * rho(kw,iw)
   seaflux(isf)%airtemp = arf_sea * airtemp
   seaflux(isf)%airshv  = arf_sea * sh_v(kw,iw)   
enddo   

! Do parallel send/recv of ATM properties in seaflux cells

if (iparallel == 1) then
   mrl = 1

   call mpi_send_wsf('A',mrl)
   call mpi_recv_wsf('A',mrl)
endif

! Loop over all SEA cells and zero properties before summation

do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

   if (iparallel == 1 .and. itab_ws(iws)%irank /= myrank) cycle

   sea%rhos(iws) = 0.
   sea%can_temp(iws) = 0.
   sea%can_shv(iws) = 0.

enddo

! Loop over all SEAFLUX cells that are APPLIED on this rank, and sum 
! atmospheric properties to corresponding SEA cell

do j = 1,jseaflux(2)%jend(1)
   isf = jseaflux(2)%iseaflux(j)

   iws = seaflux(isf)%iwls ! global index

! If run is parallel, convert iws to local domain

   if (iparallel == 1) then
      iws = itabg_ws(iws)%iws_myrank
   endif

   sea%rhos(iws)     = sea%rhos(iws)     + seaflux(isf)%rhos
   sea%can_temp(iws) = sea%can_temp(iws) + seaflux(isf)%airtemp
   sea%can_shv(iws)  = sea%can_shv(iws)  + seaflux(isf)%airshv
enddo   

! Loop over all SEA cells to initialize remaining fields

do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

   if (iparallel == 1 .and. itab_ws(iws)%irank /= myrank) cycle

   sea%seatc(iws) = sea%seatp(iws)  &
                  + (sea%seatf(iws) - sea%seatp(iws)) * timefac_sst

   sea%seaicec(iws) = sea%seaicep(iws)  &
                    + (sea%seaicef(iws) - sea%seaicep(iws)) * timefac_seaice

   if (nint(sea%seaicec(iws)) == 0) then
      sea%rough(iws) = .001
   else
      sea%rough(iws) = .0005
   endif

   sea%surface_ssh(iws) = rhovsl(sea%seatc(iws)-273.15) / sea%rhos(iws)
   sea%can_depth(iws) = 20. * max(1.,.025 * dt_sea)

enddo

! Do parallel send/recv of SEA fields

if (iparallel == 1) then
   call mpi_send_ws('A')
   call mpi_recv_ws('A')
endif

return
end subroutine sea_init_atm
