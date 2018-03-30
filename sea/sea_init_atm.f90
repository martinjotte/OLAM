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
subroutine sea_init_atm()

  use mem_sea,     only: sea, itab_ws
  use sea_coms,    only: mws, iupdsst, s1900_sst, isstfile, nsstfiles, dt_sea,  &
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles, nzi
  use mem_basic,   only: rho, press, vxe, vye, vze, tair, sh_v, theta
  use misc_coms,   only: io6, time8, s1900_sim, iparallel, isubdomain, runtype
  use mem_ijtabs,  only: itabg_w
  use consts_coms, only: t00, p00i, grav, rocp
  use mem_grid,    only: dzt_bot
  use therm_lib,   only: rhovsl, rhovsil

  implicit none

  integer :: iw
  integer :: kw
  integer :: iws

  real :: timefac_sst
  real :: timefac_seaice
  real :: dum1, dum2
  real :: prss

  ! Initialize sea quantities that do not depend on atmospheric conditions

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

  !$omp parallel do
  do iws = 2,mws

  ! Initialize sea temperature, sea ice, and canopy depth

     sea%seatc(iws) = sea%seatp(iws)  &
                    + (sea%seatf(iws) - sea%seatp(iws)) * timefac_sst

     sea%seaicec(iws) = sea%seaicep(iws)  &
                      + (sea%seaicef(iws) - sea%seaicep(iws)) * timefac_seaice

     sea%can_depth(iws) = 20. * max(1.,.025 * dt_sea)

  enddo  ! iws
  !$omp end parallel do

  ! End of initialization that does not depend on atmospheric conditions

  if (runtype /= "INITIAL") return

  ! Initialize sea quantities that depend on atmospheric conditions

  !$omp parallel do private (iw,kw,prss,dum1,dum2)
  do iws = 2,mws
     
     iw = itab_ws(iws)%iw  ! global index
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     kw = itab_ws(iws)%kw

! Transfer atmospheric properties to sea cells

     prss             = press(kw,iw) + dzt_bot(kw) * rho(kw,iw) * grav
!    sea%cantemp(iws) = tair(kw,iw)
     sea%cantemp(iws) = theta(kw,iw) * (prss * p00i)**rocp
     sea%canshv (iws) = sh_v(kw,iw)
     sea%ustar  (iws) = 0.1
     sea%ggaer  (iws) = 0.0
     sea%wthv   (iws) = 0.0

     sea%nlev_seaice(iws) = 0

     ! Seawater quantities

     sea%sea_rough  (iws) = .001
     sea%sea_cantemp(iws) = sea%cantemp(iws)
     sea%sea_canshv (iws) = sea%canshv(iws)
     sea%sea_sfc_ssh(iws) = rhovsl(sea%seatc(iws)-t00) / rho(kw,iw)
     sea%sea_ustar  (iws) = 0.1
     sea%sea_ggaer  (iws) = 0.0
     sea%sea_wthv   (iws) = 0.0

     ! Seaice quantities
   
     call prep_seaice(sea%seatc              (iws), &
                      sea%seaicec            (iws), &
                      sea%sea_cantemp        (iws), &
                      sea%ice_cantemp        (iws), &
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
                      sea%sea_canshv         (iws), &
                      sea%ice_canshv         (iws), &
                      sea%sea_ustar          (iws), &
                      sea%ice_ustar          (iws), &
                      sea%sea_ggaer          (iws), &
                      sea%ice_ggaer          (iws), &
                      sea%sea_wthv           (iws), &
                      sea%ice_wthv           (iws), &
                      sea%ice_sxfer_t        (iws), &
                      sea%ice_sxfer_r        (iws)  )     

     if (sea%nlev_seaice(iws) > 0) then

        sea%ice_sfc_ssh(iws) = rhovsil(sea%seaice_tempk(sea%nlev_seaice(iws),iws)) &
                             / rho(kw,iw)

        sea%rough      (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_rough  (iws) + &
                                      sea%seaicec(iws)  * sea%ice_rough  (iws)

        sea%cantemp    (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_cantemp(iws) + &
                                      sea%seaicec(iws)  * sea%ice_cantemp(iws)

        sea%canshv     (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_canshv (iws) + &
                                      sea%seaicec(iws)  * sea%ice_canshv (iws)

        sea%surface_ssh(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfc_ssh(iws) + &
                                      sea%seaicec(iws)  * sea%ice_sfc_ssh(iws)
 
     else

        sea%rough      (iws) = sea%sea_rough  (iws)
        sea%cantemp    (iws) = sea%sea_cantemp(iws)
        sea%canshv     (iws) = sea%sea_canshv (iws)
        sea%surface_ssh(iws) = sea%sea_sfc_ssh(iws)

     endif

  enddo
  !$omp end parallel do

end subroutine sea_init_atm
