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
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles, nzi
  use mem_sflux,   only: mseaflux, seaflux, jseaflux
  use mem_basic,   only: rho, tair, sh_v
  use misc_coms,   only: io6, time8, s1900_sim, iparallel, isubdomain, runtype
  use mem_ijtabs,  only: itabg_w
  use mem_para,    only: myrank
  use consts_coms, only: t00

  implicit none

  integer :: isf
  integer :: iw
  integer :: kw
  integer :: iws
  integer :: j

  real :: timefac_sst
  real :: timefac_seaice
  real :: arf_sea
  real :: dum1, dum2

  real, external :: rhovsl, rhovsil

! This subroutine fills the primary LEAF3 arrays which depend on current
! atmospheric conditions.

! Time interpolation factor for updating SST and SEA ICE

  timefac_sst    = 0.0
  timefac_seaice = 0.0

  if (iupdsst == 1 .and. nsstfiles > 1) then
     timefac_sst = (s1900_sim           - s1900_sst(isstfile-1))  &
                 / (s1900_sst(isstfile) - s1900_sst(isstfile-1))
  endif

  if (iupdseaice == 1 .and. nseaicefiles > 1) then
     timefac_seaice = (s1900_sim                 - s1900_seaice(iseaicefile-1)) &
                    / (s1900_seaice(iseaicefile) - s1900_seaice(iseaicefile-1))
  endif

  if (runtype /= "INITIAL") return

! Loop over all SEAFLUX cells that are EVALUATED on this rank, and transfer 
! atmospheric properties to each, with weighting according to arf_sea

  !$omp parallel do private (isf,iw,kw,arf_sea)
  do j = 1,jseaflux(1)%jend(1)
     isf = jseaflux(1)%iseaflux(j)

     iw = seaflux(isf)%iw  ! global index

     ! If run is parallel, convert iw to local domain

     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     kw      = seaflux(isf)%kw
     arf_sea = seaflux(isf)%arf_sfc

     seaflux(isf)%rhos    = arf_sea * rho(kw,iw)
     seaflux(isf)%airtemp = arf_sea * tair(kw,iw)
     seaflux(isf)%airshv  = arf_sea * sh_v(kw,iw)   
  enddo
  !$omp end parallel do

! Do parallel send of ATM properties in seaflux cells

  if (iparallel == 1) then
     call mpi_send_wsf('A', 1)
  endif

! Loop over all SEA cells and zero properties before summation

  !$omp parallel do
  do iws = 2, mws

     ! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK
     if (isubdomain == 1 .and. itab_ws(iws)%irank /= myrank) cycle

     sea%rhos(iws)     = 0.0
     sea%can_temp(iws) = 0.0
     sea%can_shv(iws)  = 0.0
  enddo
  !$omp end parallel do

! Do parallel recv of ATM properties in seaflux cells

  if (iparallel == 1) then
     call mpi_recv_wsf('A',1)
  endif

! Loop over all SEAFLUX cells that are APPLIED on this rank, and sum 
! atmospheric properties to corresponding SEA cell

  !$omp parallel do private (isf,iws)
  do j = 1,jseaflux(2)%jend(1)
     isf = jseaflux(2)%iseaflux(j)

     iws = seaflux(isf)%iwls ! global index

     ! If run is parallel, convert iws to local domain

     if (isubdomain == 1) then
        iws = itabg_ws(iws)%iws_myrank
     endif

     sea%rhos(iws)     = sea%rhos(iws)     + seaflux(isf)%rhos
     sea%can_temp(iws) = sea%can_temp(iws) + seaflux(isf)%airtemp
     sea%can_shv(iws)  = sea%can_shv(iws)  + seaflux(isf)%airshv
  enddo
  !$omp end parallel do

! Loop over all SEA cells to initialize remaining fields

  !$omp parallel do private (dum1,dum2)
  do iws = 2, mws

     ! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

     if (isubdomain == 1 .and. itab_ws(iws)%irank /= myrank) cycle

     sea%seatc(iws) = sea%seatp(iws)  &
                    + (sea%seatf(iws) - sea%seatp(iws)) * timefac_sst

     sea%seaicec(iws) = sea%seaicep(iws)  &
                      + (sea%seaicef(iws) - sea%seaicep(iws)) * timefac_seaice

     sea%can_depth(iws) = 20. * max(1.,.025 * dt_sea)

     sea%nlev_seaice(iws) = 0

     ! Seawater quantities

     sea%sea_rough  (iws) = .001
     sea%seacan_temp(iws) = sea%can_temp(iws)
     sea%seacan_shv (iws) = sea%can_shv(iws)
     sea%sea_sfc_ssh(iws) = rhovsl(sea%seatc(iws)-t00) / sea%rhos(iws)

     ! Seaice quantities
   
     call prep_seaice(sea%seatc              (iws), &
                      sea%seaicec            (iws), &
                      sea%seacan_temp        (iws), &
                      sea%icecan_temp        (iws), &
                      sea%seaice_energy(1:nzi,iws), &
                      sea%seaice_tempk (1:nzi,iws), &
                      sea%nlev_seaice        (iws), &
                      sea%ice_albedo         (iws), &
                      sea%ice_rlongup        (iws), &
                      dum1                        , &
                      dum2                        , &
                      0.0                         , &
                      0.0                         , &
                      sea%ice_rough          (iws), &
                      sea%seacan_shv         (iws), &
                      sea%icecan_shv         (iws), &
                      sea%sea_ustar          (iws), &
                      sea%ice_ustar          (iws), &
                      sea%sea_ggaer          (iws), &
                      sea%ice_ggaer          (iws), &
                      sea%ice_sxfer_t        (iws), &
                      sea%ice_sxfer_r        (iws)  )

     if (sea%nlev_seaice(iws) > 0) then

        sea%ice_sfc_ssh(iws) = rhovsil(sea%seaice_tempk(sea%nlev_seaice(iws),iws)) &
                             / sea%rhos(iws)

        sea%rough      (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_rough  (iws) + &
                                      sea%seaicec(iws)  * sea%ice_rough  (iws)

        sea%can_temp   (iws) = (1.0 - sea%seaicec(iws)) * sea%seacan_temp(iws) + &
                                      sea%seaicec(iws)  * sea%icecan_temp(iws)

        sea%can_shv    (iws) = (1.0 - sea%seaicec(iws)) * sea%seacan_shv (iws) + &
                                      sea%seaicec(iws)  * sea%icecan_shv (iws)

        sea%surface_ssh(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfc_ssh(iws) + &
                                      sea%seaicec(iws)  * sea%ice_sfc_ssh(iws)
 
     else

        sea%rough      (iws) = sea%sea_rough  (iws)
        sea%can_temp   (iws) = sea%seacan_temp(iws)
        sea%can_shv    (iws) = sea%seacan_shv (iws)
        sea%surface_ssh(iws) = sea%sea_sfc_ssh(iws)

     endif

  enddo
  !$omp end parallel do

! Do parallel send/recv of SEA fields

  if (iparallel == 1) then
     call mpi_send_ws('A')
     call mpi_recv_ws('A')
  endif

end subroutine sea_init_atm
