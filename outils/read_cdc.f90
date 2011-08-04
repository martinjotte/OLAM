subroutine cdc_stw(iyear1, imonth1, idate1, itime1,   &
       soilstate_db, soilt1, soilt2, soilw1, soilw2, snow)

implicit none

integer, intent(in) :: iyear1,imonth1,idate1,itime1

character(*), intent(in) :: soilstate_db

real, intent(out) :: soilt1(144,73)
real, intent(out) :: soilt2(144,73)
real, intent(out) :: soilw1(144,73)
real, intent(out) :: soilw2(144,73)
real, intent(out) :: snow  (144,73)

soilt1(:,:) = 273.15
soilt2(:,:) = 273.15

soilw1(:,:) = .5
soilw2(:,:) = .5

snow(:,:) = 0. 

return
end subroutine cdc_stw
