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
Module mem_rayf

! Memory for Rayleigh friction layer

real :: rayf_zmin
real :: rayf_distim
real :: rayf_expon
real :: rayfw_zmin
real :: rayfw_distim
real :: rayfw_expon  

real, allocatable :: rayf_cof(:)
real, allocatable :: rayf_cofw(:)

Contains

   subroutine rayf_init(mza,zm,zt)

! Initialize Rayleigh friction vertical profile coefficients

   implicit none
   
   integer, intent(in) :: mza
   real, intent(in) :: zm(mza),zt(mza)

   integer :: k
   real :: distimi,distimwi
   
   allocate(rayf_cof(mza),rayf_cofw(mza))

! RAYF coefficient for THIL and UMC  

   rayf_cof(1:mza) = 0.

   if (rayf_distim > 1.e-6) then
      distimi = 1. / rayf_distim
      do k = 2,mza-1
         if (zt(k) > rayf_zmin) then
            rayf_cof(k) = distimi   &
               * ((zt(k) - rayf_zmin) / (zm(mza-1) - rayf_zmin)) ** rayf_expon
         endif
      enddo
   endif

! RAYF coefficient for WMC

   rayf_cofw(1:mza) = 0.

   if (rayfw_distim > 1.e-6) then
      distimwi = 1. / rayfw_distim
      do k = 2,mza-2
         if (zm(k) > rayfw_zmin) then
            rayf_cofw(k) = distimwi   &
               * ((zm(k) - rayfw_zmin) / (zm(mza-1) - rayfw_zmin)) ** rayfw_expon
         endif
      enddo
   endif

  return
  end subroutine rayf_init

end module mem_rayf
