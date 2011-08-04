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
subroutine seacells()

use sea_coms,  only: mws, iupdsst, s1900_sst, isstfile, nsstfiles,  &
                     iupdseaice, s1900_seaice, iseaicefile, nseaicefiles

use mem_sea,   only: sea, itab_ws
use misc_coms, only: io6, time8, s1900_sim, iparallel
use mem_para,  only: myrank

!$ use omp_lib

implicit none

! Local variables

integer :: iws      ! sea cell loop counter

real :: timefac_sst   ! fraction of elapsed time from past to future SST obs
real :: timefac_seaice   ! fraction of elapsed time from past to future SEA ICE obs

! Time interpolation factors for updating SST and SEA ICE

timefac_sst = 0.
timefac_seaice = 0.

if (iupdsst == 1 .and. nsstfiles > 1) then
   timefac_sst = (s1900_sim           - s1900_sst(isstfile-1))  &
               / (s1900_sst(isstfile) - s1900_sst(isstfile-1))
endif

if (iupdseaice == 1 .and. nseaicefiles > 1) then
   timefac_seaice = (s1900_sim                 - s1900_seaice(iseaicefile-1)) &
                  / (s1900_seaice(iseaicefile) - s1900_seaice(iseaicefile-1))
endif

! Loop over ALL SEA CELLS

!$omp parallel do
do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

   if (iparallel == 1 .and. itab_ws(iws)%irank /= myrank) cycle

! Update SEA fields

   call seacell(iws                 ,  &
                timefac_sst         ,  &
                timefac_seaice      ,  &
                sea%rhos       (iws),  &
                sea%ustar      (iws),  &
                sea%sxfer_t    (iws),  &
                sea%sxfer_r    (iws),  &
                sea%can_depth  (iws),  &
                sea%seatp      (iws),  &
                sea%seatf      (iws),  &
                sea%seatc      (iws),  &
                sea%seaicep    (iws),  &
                sea%seaicef    (iws),  &
                sea%seaicec    (iws),  &
                sea%can_temp   (iws),  &
                sea%can_shv    (iws),  &
                sea%surface_ssh(iws),  &
                sea%rough      (iws)   )

! Zero out SEA%SXFER_T(IWS) and SEA%SXFER_R(IWS) now that they have 
! been applied to the canopy (but save previous values for plotting)

   sea%sxfer_tsav(iws) = sea%sxfer_t(iws)
   sea%sxfer_rsav(iws) = sea%sxfer_r(iws)

   sea%sxfer_t(iws) = 0.
   sea%sxfer_r(iws) = 0.
                   
enddo
!$omp end parallel do

return
end subroutine seacells

!===============================================================================

subroutine seacell(iws,timefac_sst,timefac_seaice,rhos,ustar,sxfer_t,sxfer_r, &
                   can_depth,seatp,seatf,seatc,seaicep,seaicef,seaicec,       &
                   can_temp,can_shv,surface_ssh,rough                         )

use sea_coms,    only: dt_sea
use consts_coms, only: cp, grav
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iws     ! current sea cell index

real, intent(in) :: timefac_sst  ! frac of elapsed time from past to future SST obs
real, intent(in) :: timefac_seaice  ! frac of elapsed time from past to future  &
                                    ! SEA ICE obs
real, intent(in) :: rhos         ! air density [kg/m^3]
real, intent(in) :: ustar        ! friction velocity [m/s]
real, intent(in) :: sxfer_t      ! can_air to atm heat xfer this step [kg_air K/m^2]
real, intent(in) :: sxfer_r      ! can_air to atm vapor xfer this step (kg_vap/m^2]
real, intent(in) :: can_depth    ! "canopy" depth for heat and vap capacity [m]
real, intent(in) :: seatp        ! past sea temp (obs time) [K]
real, intent(in) :: seatf        ! future sea temp (obs time) [K]
real, intent(in) :: seaicep         ! past sea ice cover (obs time) [0-1]
real, intent(in) :: seaicef         ! future sea ice cover (obs time) [0-1]

real, intent(out  ) :: seatc     ! current sea temp [K]
real, intent(out  ) :: seaicec      ! current sea ice cover [0-1]
real, intent(inout) :: can_temp  ! "canopy" air temp [K]
real, intent(inout) :: can_shv   ! "canopy" air vapor spec hum [kg_vap/kg_air]

real, intent(out) :: surface_ssh ! sea surface sat spec hum [kg_vap/kg_air] 
real, intent(out) :: rough       ! sea cell roughess height [m] 

! Local parameter

real, parameter :: z0fac_water = .016 / grav  ! factor for sea roughness height

! Local variables

real :: rdi     ! stomatal conductance [m/s]
real :: hxfergc ! heat xfer from soil/sfcwater to can_air this step [J/m^2]
real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]
real :: wxfergc ! vapor xfer from soil/sfcwater to can_air this step [kg_vap/m^2]

real, external :: rhovsil  ! function to compute saturation vapor density

! Update SEATC and sea ice

seatc   = seatp   + timefac_sst * (seatf - seatp)
seaicec = seaicep + timefac_seaice * (seaicef - seaicep)

! Evaluate surface saturation specific humidity

surface_ssh = rhovsil(seatc-273.15) / rhos

! Update temperature and vapor specific humidity of "canopy" from
! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity derived from Garratt (1992).

rdi = .2 * ustar
hxfergc = dt_sea * cp * rhos * rdi * (seatc - can_temp)
wxfergc = dt_sea *      rhos * rdi * (surface_ssh - can_shv)

hxferca = cp * sxfer_t  ! sxfer_t and sxfer_r already incorporate dt_sea

can_temp = can_temp + (hxfergc - hxferca) / (can_depth * rhos * cp)
can_shv  = can_shv  + (wxfergc - sxfer_r) / (can_depth * rhos)             

! Evaluate sea cell roughness height

if (nint(seaicec) == 0)then
!   rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)
   rough = 10. * exp(-10. / ustar ** .333333) ! Davis et al. (2008)
   rough = max(.125e-6,min(2.85e-3,rough))    ! Davis et al. (2008)
else
   rough = 5.0e-4 ! [m]; taken from Los Alamos sea ice model.
endif

return
end subroutine seacell
