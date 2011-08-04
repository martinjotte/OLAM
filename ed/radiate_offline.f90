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
subroutine radiate_offline()
  
  use mem_leaf, only: land, first_site
  use ed_structure_defs
  use misc_coms, only: current_time, radfrq, dtlong
  use canopy_radiation_coms, only: visible_fraction_dir, visible_fraction_dif

  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: cp
  real :: total_beam_radiation
  real :: total_diffuse_radiation

  ! Check whether it is time to update radiative fluxes and heating rates
  
  if (mod(current_time%time + .001d0,dble(radfrq)) < dtlong) then
     
     ! Print message that radiative transfer is being computed
     
!     write(*,'(a,i)')'Radiation tendencies updated at UTC time =',  &
!          nint(current_time%time)
     
     ! Compute solar zenith angle [land%cosz]
     
     call zen_offline()

     ! Set the pointer for the first ED site, and then loop over sites
     cs => first_site
     do while(associated(cs))

        ! Compute the visible fraction of diffuse and beam radiation
        ! needed by the radiative transfer routine.
        total_beam_radiation = cs%metinput%nir_beam + cs%metinput%par_beam
        total_diffuse_radiation = cs%metinput%nir_diffuse +   &
             cs%metinput%par_diffuse
        if(total_beam_radiation > 0.0)then
           visible_fraction_dir = cs%metinput%par_beam / total_beam_radiation
        else
           visible_fraction_dir = 0.5
        endif
        if(total_diffuse_radiation > 0.0)then
           visible_fraction_dif = cs%metinput%par_diffuse /   &
                total_diffuse_radiation 
        else
           visible_fraction_dif = 0.5
        endif

        ! Define rshort, rshort_diffuse
        land%rshort(cs%iland) = total_beam_radiation + total_diffuse_radiation
        land%rshort_diffuse(cs%iland) = total_diffuse_radiation

        ! Loop over subgrid-scale patches.
        cp => cs%oldest_patch
        do while(associated(cp))

           ! Get unnormalized radiative transfer information.
           call sfcrad_ed(land%cosz(cs%iland), cp, cp%cohort_count)

           ! Normalize the absorbed radiations.
           call scale_ed_radiation(cp)
           cp => cp%younger
        enddo

        ! Update all albedos and rlongup.
        call ed2land_radiation(cs)

        cs => cs%next_site
     enddo

  endif

  return
end subroutine radiate_offline


!******************************************************************************

subroutine zen_offline()

  use misc_coms, only: current_time
  use consts_coms, only: pio180
  use mem_leaf, only: land, first_site
  use ed_structure_defs

  implicit none
  
  integer :: jday
  integer, external :: julday
  real :: declin
  real :: sdec
  real :: cdec
  real :: d0
  real :: d02
  real :: solfac
  real :: dayhr
  type(site), pointer :: cs
  real :: radlat
  real :: cslcsd
  real :: snlsnd
  real :: dayhrr
  real :: hrangl

  jday  = julday(current_time%month, current_time%date, current_time%year)

  ! sdec - sine of declination, cdec - cosine of declination
  
  declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pio180
  sdec = sin(declin)
  cdec = cos(declin)

  ! Find the factor, solfac, to multiply the solar constant to correct
  ! for Earth's varying distance to the sun.
  
  d0 = 6.2831853 * float(jday-1) / 365.
  d02 = d0 * 2.
  solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
       + 0.000719 * cos(d02) + 0.000077 * sin(d02)

  ! Find the hour angle, then get cosine of zenith angle.
  
  dayhr = current_time%time / 3600.

  cs => first_site
  do while(associated(cs))
     radlat = cs%lat * pio180
     if (radlat == declin) radlat = radlat + 1.e-5
     cslcsd = cos(radlat) * cdec
     snlsnd = sin(radlat) * sdec
     dayhrr = mod(dayhr+cs%lon/15.+24.,24.)
     hrangl = 15. * (dayhrr - 12.) * pio180
     land%cosz(cs%iland) = snlsnd + cslcsd * cos(hrangl)
     cs => cs%next_site
  enddo

  return
end subroutine zen_offline

