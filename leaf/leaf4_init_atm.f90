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

subroutine leaf4_init_atm()

use mem_leaf,    only: land, itabg_wl, itab_wl

use leaf_coms,   only: mwl, nzg, nzs, &
                       veg_ht, soilcp, slmstr, slmsts, slcpd, soil_rough, &
                       iupdndvi, s1900_ndvi, indvifile, nndvifiles, &
                       dt_leaf, isoilstateinit, iwatertabflg, watertab_db, &
                       water_frac_ph0, water_def_ph0, slpott, slpots, slbs, slzt
                      
use mem_sflux,    only: landflux, jlandflux
use mem_basic,    only: rho, tair, sh_v
use misc_coms,    only: io6, time8, s1900_sim, iparallel, runtype, initial
use mem_ijtabs,   only: itabg_w
use consts_coms,  only: cliq, cice, alli, cliq1000, cice1000, alli1000
use mem_para,     only: myrank
use leaf4_canopy, only: vegndvi
use leaf_db,      only: leaf_database_read

implicit none

integer :: k
integer :: ntext
integer :: iw
integer :: kw
integer :: ilf
integer :: iwl
integer :: leaf_class
integer :: nlsw1 ! maximum of (1,nlev_sfcwater)
integer :: mrl
integer :: j

real :: airtemp
real :: timefac_ndvi
real :: arf_land
real :: wq, wq_added
real :: headp_phi  ! Inverse of derivative of hydraulic pressure head
                   ! wrt soil water
real :: psi

! automatic arrays

real :: soil_tempc(nzg,mwl)  ! initial soil temperature (C)
real :: fracliq(nzg,mwl)     ! initial soil liquid fraction (0-1)
real :: wtd(mwl)             ! watertable depth from database

real, external :: rhovsl

! This subroutine fills the primary LEAF4 arrays that depend on current
! atmospheric conditions.

! Time interpolation factor for updating NDVI

timefac_ndvi = 0.

if (iupdndvi == 1 .and. nndvifiles > 1) then
   timefac_ndvi = (s1900_sim               - s1900_ndvi(indvifile))  &
                / (s1900_ndvi(indvifile+1) - s1900_ndvi(indvifile))
endif

if (runtype /= "INITIAL") return

! If iwatertabflg = 1, read water table depth from database; put into wtd array

   if (iwatertabflg == 1) then

      call leaf_database_read(mwl, &
           land%glatw,             &
           land%glonw,             &
           watertab_db,            &
           watertab_db,            &
           'wtd',                  &
           datp=wtd                )

   endif

! Loop over all LANDFLUX cells that are EVALUATED on this rank, and transfer
! atmospheric properties to each, with weighting according to arf_land

do j = 1,jlandflux(1)%jend(1)
   ilf = jlandflux(1)%ilandflux(j)

   iw = landflux(ilf)%iw  ! global index

! If run is parallel, convert iw to local domain

   if (iparallel == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif

   kw       = landflux(ilf)%kw
   arf_land = landflux(ilf)%arf_sfc

   landflux(ilf)%rhos    = arf_land * rho (kw,iw)
   landflux(ilf)%airtemp = arf_land * tair(kw,iw)
   landflux(ilf)%airshv  = arf_land * sh_v(kw,iw)
enddo

! Do parallel send/recv of ATM properties in landflux cells

if (iparallel == 1) then
   mrl = 1

   call mpi_send_wlf('A',mrl)
   call mpi_recv_wlf('A',mrl)
endif

! Loop over all LAND cells and zero properties before summation

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   land%rhos(iwl) = 0.
   land%can_temp(iwl) = 0.
   land%can_shv(iwl) = 0.

enddo

! Loop over all LANDFLUX cells that are APPLIED on this rank, and sum
! atmospheric properties to corresponding LAND cell.

do j = 1,jlandflux(2)%jend(1)
   ilf = jlandflux(2)%ilandflux(j)

   iwl = landflux(ilf)%iwls ! global index

! If run is parallel, convert iwl to local domain.

   if (iparallel == 1) then
      iwl = itabg_wl(iwl)%iwl_myrank
   endif

   land%rhos    (iwl) = land%rhos    (iwl) + landflux(ilf)%rhos
   land%can_temp(iwl) = land%can_temp(iwl) + landflux(ilf)%airtemp
   land%can_shv (iwl) = land%can_shv (iwl) + landflux(ilf)%airshv

enddo

! Loop over all LAND cells and perform default initialization of remaining fields

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

! Set vegetation parameters

   leaf_class = land%leaf_class(iwl)

   land%rough      (iwl) = soil_rough
   land%veg_rough  (iwl) = .13 * veg_ht(leaf_class)
   land%veg_height (iwl) = veg_ht(leaf_class)   
   land%stom_resist(iwl) = 1.e6
   land%veg_temp   (iwl) = land%can_temp(iwl)
   land%veg_water  (iwl) = 0.

! For now, choose heat/vapor capacities for stability based on timestep   

   land%can_depth(iwl) = 20. * max(1.,.025 * dt_leaf)
   land%hcapveg  (iwl) = 3.e4 * max(1.,.025 * dt_leaf)

! Initialize vegetation TAI, LAI, fractional area, albedo, and roughness

   call vegndvi(iwl,                    &
                leaf_class           ,  &
                timefac_ndvi         ,  &
                land%veg_height  (iwl), &
                land%veg_ndvip   (iwl), &
                land%veg_ndvif   (iwl), &
                land%veg_ndvic   (iwl), &
                land%veg_tai     (iwl), &
                land%veg_lai     (iwl), &
                land%veg_fracarea(iwl), &
                land%veg_albedo  (iwl), &
                land%veg_rough   (iwl)  )

! Initialize hydraulic head at bottom of soil grid...

   if (leaf_class == 17 .or. leaf_class == 20) then

! For bog, marsh, wetland areas, use head0 = +10 cm; this takes precedence
! over using watertable database
 
      land%head0(iwl) = 0.1

   elseif (iwatertabflg == 1 .and. wtd(iwl) > -1.e-6) then

! If area is not bog, marsh, or wetland, and if iwatertabflg = 1, use
! head0 = -watertable depth from database, unless it is "missing" (negative)

      land%head0(iwl) = -wtd(iwl)

   elseif (leaf_class == 3) then

! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is desert

      land%head0(iwl) = -100.0  ! Desert

   elseif (leaf_class == 10) then

! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is semi-desert

      land%head0(iwl) = -30.0  ! Semi-desert

   else

! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is none of the above

      land%head0(iwl) = -10.0

   endif

! Default initialization of sfcwater_mass, soil_tempc, and soil_water

   land%sfcwater_mass  (1:nzs,iwl) = 0.
   land%sfcwater_energy(1:nzs,iwl) = 0.
   land%sfcwater_depth (1:nzs,iwl) = 0.

   soil_tempc(1:nzg,iwl) = land%can_temp(iwl) - 273.15
   fracliq(1:nzg,iwl) = 1.0

! Loop over soil layers

   do k = 1,nzg
      ntext = land%ntext_soil(k,iwl)

! If initial head exceeds slpott, estimate soil_water based on total
! head (molecular plus hydraulic)       

      if (land%head0(iwl) - slzt(k) > slpott(ntext)) then
      
         headp_phi = (water_def_ph0 * slmsts(ntext)) / (10. - slpott(ntext))

! First estimate using slpott in place of psi

         land%soil_water(k,iwl) = water_frac_ph0 * slmsts(ntext) &
            + (land%head0(iwl) - slzt(k) - slpott(ntext)) * headp_phi 

! Evaluate molecular water potential from soil water estimate

         psi = slpots(ntext) &
             * (slmsts(ntext) / land%soil_water(k,iwl)) ** slbs(ntext)

! Correction using psi

         land%soil_water(k,iwl) = water_frac_ph0 * slmsts(ntext) &
            + (land%head0(iwl) - slzt(k) - psi) * headp_phi
            
      else

! If initial head does not exceed slpott, estimate soil_water based on
! molecular head alone       

         land%soil_water(k,iwl) = max(soilcp(ntext),slmsts(ntext) &
            * (slpots(ntext) / (land%head0(iwl) - slzt(k))) ** (1./slbs(ntext)))

      endif

!! Old_method: Default soil moisture initialization from namelist slmstr(k)
!!     land%soil_water(k,iwl) = max(soilcp(ntext),slmstr(k) * slmsts(ntext))
      
   enddo

enddo

! Overwrite the default soil initialization with observed data if specified

if (isoilstateinit == 1) then

   ! read sfcwater mass, soil_tempc, and soil_water from the 2 X 2.5 degree
   ! NCEP/NCAR reanalysis in netcdf format, obtained from
   ! ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/

   call read_soil_moist_temp(soil_tempc)

elseif (isoilstateinit == 2 .and. initial == 2) then

   ! read sfcwater mass, soil_tempc, and soil_water saved in the initial
   ! degribbed analysis files

   call read_soil_analysis(soil_tempc)

endif

!--------------------------------------------------------------------------------
! ADD A METHOD HERE TO INITIALIZE FRACTION OF SOIL WATER THAT IS LIQUID (FRACLIQ)
! BASED ON MODEL OBSERVED AND/OR MODEL CLIMATOLOGY
!--------------------------------------------------------------------------------

! Loop over all LAND cells

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

! Leaf classes 17 and 20 represent persistent wetlands (bogs, marshes, fens,
! swamps).  Initialize these areas with saturated soil, and with 0.1 m of standing
! surface water (sfcwater) added to whatever is already present (e.g., from obs).

   if (land%leaf_class(iwl) == 17 .or. land%leaf_class(iwl) == 20) then

! Vertical loop over soil layers: saturate soil layers

      do k = 1,nzg
         ntext = land%ntext_soil(k,iwl)
         land%soil_water(k,iwl) = slmsts(ntext)
      enddo

! Since sfcwater_energy has units of J/kg, first convert to J/m^2 before adding
! wetland sfcwater.
      
      wq = land%sfcwater_mass(1,iwl) * land%sfcwater_energy(1,iwl)

! Add wetland sfcwater mass and depth

      land%sfcwater_mass(1,iwl) = land%sfcwater_mass(1,iwl)  &
                                + 100.  ! 100 kg/m^2 equivalent to 0.1 m

      land%sfcwater_depth(1,iwl) = land%sfcwater_depth(1,iwl)  &
                                 + .1   ! 0.1 m added depth

! Add wetland sfcwater energy, which is assumed to have energy of liquid
! water at canopy air temperature, which could be below freezing.

      wq_added = 100.  &  ! 100 kg/m^2 added mass
               * ((land%can_temp(iwl) - 273.15) * cliq + alli) ! J/kg of added
                                                               ! liquid water
! Diagnose new sfcwater energy

      land%sfcwater_energy(1,iwl) = (wq + wq_added) / land%sfcwater_mass(1,iwl)
      
! Set head0 to maintain 10 cm of surface water

      land%head0(iwl) = .10

   endif
   
! If leaf_class of this IWL land cell is ice cap or glacier, assume that 
! all 'soil' water is in ice phase and that soil_tempc is at or below 0.
! (Soil will be replaced by firn model in the future.)

   if (land%leaf_class(iwl) == 2) then

      soil_tempc(1:nzg,iwl) = min(0.,soil_tempc(1:nzg,iwl))
      fracliq(1:nzg,iwl) = 0.

   endif

! Initialize soil energy [J/m^3] from given soil textural class, temperature, 
! total water content, and liquid fraction (as opposed to ice fraction) of the 
! soil water that is present.  

   do k = 1,nzg
      ntext = land%ntext_soil(k,iwl)

      if (soil_tempc(k,iwl) > 0.) then

         land%soil_energy(k,iwl)                                      &
             = soil_tempc(k,iwl) * slcpd(ntext)                       &
             + soil_tempc(k,iwl) * land%soil_water(k,iwl) * cliq1000  &
             + fracliq(k,iwl)    * land%soil_water(k,iwl) * alli1000
             
      else
      
         land%soil_energy(k,iwl)                                       &
            = soil_tempc(k,iwl) * slcpd(ntext)                         &
            + soil_tempc(k,iwl) * land%soil_water(k,iwl) * cice1000    &
            + fracliq(k,iwl)    * land%soil_water(k,iwl) * alli1000
             
      endif
   enddo

! Determine active number of surface water levels

   do k = 1,nzs
      if (land%sfcwater_mass(k,iwl) > 1.e-3) then
         land%nlev_sfcwater(iwl) = k
      endif
   enddo

! Initialize ground (soil) and surface vapor specific humidity

   nlsw1 = max(1,land%nlev_sfcwater(iwl))
   
   call grndvap(iwl,                             &
                land%nlev_sfcwater        (iwl), &
                land%ntext_soil       (nzg,iwl), &
                land%soil_water       (nzg,iwl), &
                land%soil_energy      (nzg,iwl), &
                land%sfcwater_energy(nlsw1,iwl), &
                land%rhos                 (iwl), &
                land%can_shv              (iwl), &
                land%surface_ssh          (iwl), &
                land%ground_shv           (iwl)  )

enddo

! Do parallel send/recv of LEAF fields

if (iparallel == 1) then
  call mpi_send_wl('A')
  call mpi_recv_wl('A')
endif

return
end subroutine leaf4_init_atm

!=====================================================================

subroutine read_soil_moist_temp(soil_tempc)

  use misc_coms, only: io6, iyear1, imonth1, idate1, itime1
  use leaf_coms, only: nzg, mwl, slz, soilcp, slmsts, slcpd, soilstate_db
  use mem_leaf,  only: land, itab_wl
  use consts_coms, only: pio180, piu180, erad, cliq1000, alli1000, cice,   &
       cice1000

  implicit none

  real, intent(inout) :: soil_tempc(nzg,mwl)

  integer :: ilatd
  integer :: ilond
  integer :: i
  real, dimension(144*73) :: soilt1 ! soil temperature, 0-10 cm
  real, dimension(144*73) :: soilt2 ! soil temperature, 10-200cm
  real, dimension(144*73) :: soilw1 ! soil water, 0-10 cm
  real, dimension(144*73) :: soilw2 ! soil water, 10-200cm
  real, dimension(144*73) :: snow ! snow depth [kg/m2]
  integer, dimension(144*73) :: miss_flag
  real :: latd
  real :: lond
  real :: best_dist
  integer :: ilatdd
  integer :: ilondd
  real :: latdd
  real :: londd
  integer :: ii
  real :: dist
  integer :: iwl
  real :: glat
  real :: glon
  real :: snowdens
  integer :: k
  integer :: ntext

  ! Acquire soil temperature, moisture, and snow depth from reanalysis.
  ! Input data is on a 2.5 x 2.5 degree grid.  Latitudes start at 90N, 
  ! longitudes at 1.25 E.

  call cdc_stw(iyear1, imonth1, idate1, itime1,   &
       trim(soilstate_db)//char(0), soilt1, soilt2, soilw1, soilw2, snow)

  ! Determine which grid cells have missing data.  

  do ilatd = 1,73
     do ilond = 1,144
        i = 144 * (ilatd-1) + ilond

        ! Make sure there are values for soil temp, moist.  There may
        ! be nans in NCEP.

        if (soilt1(i) == soilt1(i) .and. soilt2(i) == soilt2(i) .and.  &
            soilw1(i) == soilw1(i) .and. soilw2(i) == soilw2(i)) then

           if (soilt1(i) > 0.0 .and. soilt2(i) > 0.0 .and.  &
               soilw1(i) > 0.0 .and. soilw2(i) > 0.0) then
              miss_flag(i) = 0
           else
              miss_flag(i) = 1
           endif
        else
           miss_flag(i) = 1
        endif

     enddo
  enddo

  ! Now fill missing data with nearest neighbor

  do ilatd = 1,73
     latd = pio180 * (92.5 - ilatd *2.5)
     do ilond = 1, 144
        lond = pio180 * (2.5 * ilond - 1.25)
        i = 144 * (ilatd-1) + ilond 
        if(miss_flag(i) == 1)then
           ! Find nearest neighbor with viable data
           best_dist = 1.0e30
           do ilatdd = 1, 73
              latdd = pio180 * (92.5 - ilatdd *2.5)
              do ilondd = 1, 144
                 londd = pio180 * (2.5 * ilondd - 1.25)
                 ii = 144 * (ilatdd-1) + ilondd 
                 ! Make sure we are comparing to a viable datum
                 if(miss_flag(ii) == 0)then
                    dist = sqrt((cos(latdd)*cos(londd)-cos(latd)*  &
                         cos(lond))**2+(cos(latdd)*sin(londd)-cos(latd)* &
                         sin(lond))**2+(sin(latd)-sin(latdd))**2)
                    if(dist < best_dist)then
                       best_dist = dist
                       soilt1(i) = soilt1(ii)
                       soilt2(i) = soilt2(ii)
                       soilw1(i) = soilw1(ii)
                       soilw2(i) = soilw2(ii)
                    endif
                 endif
              enddo
           enddo
        endif
     enddo
  enddo
  
  ! Now, fill the land% arrays.
  do ilatd = 1,73
     do ilond = 1,144
        i = 144 * (ilatd-1) + ilond 

        ! Loop over land points
        do iwl = 2,mwl
           
           ! Land point lat, lon

           glat = land%glatw(iwl)
           glon = land%glonw(iwl)
           
           if (glon < 0.0) glon = glon + 360.0
           
           ! Find reanalysis point corresponding to this land point

           if (glat >= 0.0) then
              ilatdd = nint((90.0 - glat)/2.5) + 1
           else
              ilatdd = 73 - nint((90.0 - abs(glat))/2.5)
           endif
           ilondd = int(glon/2.5) + 1
           
           ! If there we are at the right point, fill the array

           if (ilatdd == ilatd .and. ilondd == ilond) then
              
              if (snow(i) > 1.0e-3 .and. snow(i) == snow(i)) then
                 land%sfcwater_mass(1,iwl) = snow(i)
                 land%sfcwater_energy(1,iwl) = min(0.,   &
                      (land%can_temp(iwl) - 273.15) * cice) 
                 ! snow density calculation comes from CLM3.0 documentation 
                 ! which is based on Anderson 1975 NWS Technical Doc # 19 
                 snowdens = 50.0
                 if (land%can_temp(iwl) > 258.15) snowdens =   &
                      50.0 + 1.5 * (land%can_temp(iwl) - 258.15)**1.5
                 land%sfcwater_depth(1,iwl) = land%sfcwater_mass(1,iwl) / snowdens
              endif

              ! Vertical loop over soil layers

              do k = 1,nzg

                 ntext = land%ntext_soil(k,iwl)

                 ! Determine if this layer corresponds to the first 
                 ! or second 'data' layer.

                 if (abs(slz(k)) < 0.1) then
                    soil_tempc(k,iwl) = soilt1(i) - 273.15
                    land%soil_water(k,iwl) = max(soilcp(ntext),   &
                         soilw1(i) * slmsts(ntext))
                 else
                    soil_tempc(k,iwl) = soilt2(i) - 273.15
                    land%soil_water(k,iwl) = max(soilcp(ntext),   &
                         soilw2(i) * slmsts(ntext))
                 endif

              enddo
              
           endif
           
        enddo

     enddo
  enddo

  return
end subroutine read_soil_moist_temp
