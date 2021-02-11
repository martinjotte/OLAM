subroutine read_soil_moist_temp(soil_tempc)

  use misc_coms, only: io6, iyear1, imonth1, idate1, itime1, isubdomain
  use leaf_coms, only: soilstate_db
  use mem_land,  only: land, mland, omland, nzg, slz, wresid_vg, wsat_vg
  use mem_sfcg,   only: sfcg, itab_wsfc
  use consts_coms, only: pio180, piu180, erad, cliq1000, alli1000, cice,   &
       cice1000
  use mem_para,   only: myrank

  implicit none

  real, intent(inout) :: soil_tempc(nzg,mland)

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
  integer :: iland, iwsfc
  real :: glat
  real :: glon
  real :: snowdens
  integer :: k

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
        do iland = 2,mland
           iwsfc = iland + omland
           
           ! Skip this cell if running in parallel and cell rank is not MYRANK
           if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

           ! Land point lat, lon

           glat = sfcg%glatw(iwsfc)
           glon = sfcg%glonw(iwsfc)
           
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
                 land%sfcwater_mass(1,iland) = snow(i)
                 land%sfcwater_energy(1,iland) = min(0.,   &
                      (sfcg%cantemp(iwsfc) - 273.15) * cice) 
                 ! snow density calculation comes from CLM3.0 documentation 
                 ! which is based on Anderson 1975 NWS Technical Doc # 19 
                 if (sfcg%cantemp(iwsfc) > 258.15) then
                    snowdens = 50.0 + 1.5 * (sfcg%cantemp(iwsfc) - 258.15)**1.5
                 else
                    snowdens = 50.0
                 endif
                 land%sfcwater_depth(1,iland) = land%sfcwater_mass(1,iland) / snowdens
              endif

              ! Vertical loop over soil layers

              do k = 1,nzg

                 ! Determine if this layer corresponds to the first 
                 ! or second 'data' layer.

                 if (abs(slz(k)) < 0.1) then
                    soil_tempc(k,iland) = soilt1(i) - 273.15
                    land%soil_water(k,iland) = soilw1(i)
                 else
                    soil_tempc(k,iland) = soilt2(i) - 273.15
                    land%soil_water(k,iland) = soilw2(i)
                 endif

                 land%soil_water(k,iland) = max(land%wresid_vg(k,iland)), &
                 land%soil_water(k,iland) * land%wsat_vg(k,iland)

              enddo
              
           endif
           
        enddo

     enddo
  enddo

end subroutine read_soil_moist_temp
