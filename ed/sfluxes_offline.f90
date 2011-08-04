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
subroutine sfluxes_offline()

  use ed_structure_defs
  use mem_leaf, only: land, first_site
  use consts_coms, only: p00i, rocp, cp, r8
  use misc_coms, only: dtlm

  implicit none

  type(site),  pointer :: cs
  type(patch), pointer :: ed_patch
  real :: vkmsfc
  real :: sfluxt
  real :: sfluxr
  real :: ustar0
  real :: thetaatm
  real :: thetacan
  real :: hgtmin
  real(r8) :: rhos

  cs => first_site
  do while(associated(cs))

     land%ustar(cs%iland) = 0.0
     rhos = land%rhos(cs%iland)

     ed_patch => cs%oldest_patch
     do while(associated(ed_patch))
        hgtmin = max(ed_patch%rough + 1.0e-3,  &
             max(50.0, cs%metinput%geoht))
!             max(ed_patch%veg_height, cs%metinput%geoht))
        call stars(hgtmin,         &
             ed_patch%rough  ,     &
             land%vels(cs%iland),  &
             rhos,                 &
             cs%metinput%atm_tmp,  &
             cs%metinput%atm_shv,  &
             ed_patch%can_temp,    &
             ed_patch%can_shv ,    &
             vkmsfc,               &
             sfluxt,               &
             sfluxr,               &
             ustar0                )
        
        ed_patch%sxfer_t = dtlm(1) * sfluxt
        ed_patch%sxfer_r = dtlm(1) * sfluxr
        ed_patch%ustar = ustar0

        land%sxfer_t(cs%iland) = land%sxfer_t(cs%iland)  &
                               + dtlm(1) * sfluxt * ed_patch%area
        land%sxfer_r(cs%iland) = land%sxfer_r(cs%iland)  &
                               + dtlm(1) * sfluxr * ed_patch%area

        land%ustar(cs%iland) = land%ustar(cs%iland) + ustar0 * ed_patch%area

        ed_patch => ed_patch%younger
     enddo

     cs => cs%next_site
  enddo

  return
end subroutine sfluxes_offline
