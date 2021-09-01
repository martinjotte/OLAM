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

  use mem_sea,     only: sea, msea, omsea
  use sea_coms,    only: iupdsst, s1900_sst, isstfile, nsstfiles, dt_sea,  &
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles, nzi
  use mem_basic,   only: rho, press, rr_v, theta
  use mem_micro,   only: rr_c
  use misc_coms,   only: s1900_sim, isubdomain, runtype
  use mem_ijtabs,  only: itabg_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use consts_coms, only: t00, p00i, rocp
  use therm_lib,   only: rhovsl, rhovsil
  use mem_para,    only: myrank

  implicit none

  integer :: iw
  integer :: kw
  integer :: isea, iwsfc, j

  real :: timefac_sst
  real :: timefac_seaice
  real :: dum1, dum2
  real :: vels, psfc

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

  !$omp parallel do private (iwsfc)
  do isea = 2,msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Initialize sea depth, sea temperature, sea ice, and canopy depth

     sea%wdepth(isea) = sfcg%topw(iwsfc) - sfcg%bathym(iwsfc)

     sea%seatc(isea) = sea%seatp(isea)  &
                     + (sea%seatf(isea) - sea%seatp(isea)) * timefac_sst

     sea%seaicec(isea) = sea%seaicep(isea)  &
                       + (sea%seaicef(isea) - sea%seaicep(isea)) * timefac_seaice

     sfcg%can_depth(iwsfc) = 20. * max(1.,.025 * dt_sea)
  enddo
  !$omp end parallel do

  ! End of initialization that does not depend on atmospheric conditions

  if (runtype /= "INITIAL") return

  !$omp parallel do private (iwsfc,dum1,dum2)
  do isea = 2,msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply initial atmospheric properties to sea "canopy"

     sfcg%cantemp(iwsfc) = sfcg%airtheta(iwsfc) * (sfcg%prss(iwsfc) * p00i)**rocp
     sfcg%canrrv (iwsfc) = sfcg%airrrv(iwsfc)
     sfcg%ustar  (iwsfc) = 0.1
     sfcg%ggaer  (iwsfc) = 0.0
     sfcg%wthv   (iwsfc) = 0.0

     sea%nlev_seaice(isea) = 0

     ! Seawater quantities

     sea%sea_rough   (isea) = .001
     sea%sea_cantemp (isea) = sfcg%cantemp(iwsfc)
     sea%sea_canrrv  (isea) = sfcg%canrrv(iwsfc)
     sea%sea_sfc_srrv(isea) = rhovsl(sea%seatc(isea)-t00) / sfcg%rhos(iwsfc)
     sea%sea_ustar   (isea) = 0.1
     sea%sea_ggaer   (isea) = 0.0
     sea%sea_wthv    (isea) = 0.0

     ! Seaice quantities

     call prep_seaice(sea%seatc              (isea), &
                      sea%seaicec            (isea), &
                      sea%sea_cantemp        (isea), &
                      sea%ice_cantemp        (isea), &
                      sea%seaice_energy(1:nzi,isea), &
                      sea%seaice_tempk (1:nzi,isea), &
                      sea%nlev_seaice        (isea), &
                      sea%ice_albedo         (isea), &
                      sea%ice_rlongup        (isea), &
                      dum1                         , &
                      dum2                         , &
                      0.0                          , &
                      0.0                          , &
                      sea%ice_rough          (isea), &
                      sea%sea_canrrv         (isea), &
                      sea%ice_canrrv         (isea), &
                      sea%sea_ustar          (isea), &
                      sea%ice_ustar          (isea), &
                      sea%sea_ggaer          (isea), &
                      sea%ice_ggaer          (isea), &
                      sea%sea_wthv           (isea), &
                      sea%ice_wthv           (isea), &
                      sea%ice_sxfer_t        (isea), &
                      sea%ice_sxfer_r        (isea)  )

     if (sea%nlev_seaice(isea) > 0) then

        sea%ice_sfc_srrv(isea) = rhovsil(sea%seaice_tempk(sea%nlev_seaice(isea),isea)) &
                               / sfcg%rhos(iwsfc)

        sfcg%rough    (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_rough  (isea) + &
                                       sea%seaicec(isea)  * sea%ice_rough  (isea)

        sfcg%cantemp  (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_cantemp(isea) + &
                                       sea%seaicec(isea)  * sea%ice_cantemp(isea)

        sfcg%canrrv   (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_canrrv (isea) + &
                                       sea%seaicec(isea)  * sea%ice_canrrv (isea)

        sea%surface_srrv(isea) = (1.0 - sea%seaicec(isea)) * sea%sea_sfc_srrv(isea) + &
                                        sea%seaicec(isea)  * sea%ice_sfc_srrv(isea)
 
     else

        sfcg%rough      (iwsfc) = sea%sea_rough   (isea)
        sfcg%cantemp    (iwsfc) = sea%sea_cantemp (isea)
        sfcg%canrrv     (iwsfc) = sea%sea_canrrv  (isea)
        sea%surface_srrv(isea)  = sea%sea_sfc_srrv(isea)

     endif

  enddo
  !$omp end parallel do

end subroutine sea_init_atm
