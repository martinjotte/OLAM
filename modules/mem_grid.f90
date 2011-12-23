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
Module mem_grid

   integer :: &  ! N values for full-domain reference on any process

      nza           &  ! Vertical number of all points
     ,nsw_max       &  ! Max # vert atm levels in IW column with sfc flux
     ,nma,nua,nva,nwa  ! Horiz number of all M,U,V,W points in full domain

   integer :: &

      mza           &  ! Vertical number of all points
     ,mma,mua,mva,mwa  ! Horiz number of all M,U,V,W points in sub-process

   integer, allocatable, dimension(:) ::  &

      lpm,lpu,lpv,lpw   &  ! Lowest prognosed M,U,V,W/T
     ,lcu,lcv           &  ! Lowest nonzero control volume for U,V
     ,lsw                  ! number of W/T levels in contact with surface

   real, allocatable, dimension(:) ::  &

      zm    ,zt        &  ! Z coordinate at M,T point
     ,dzm   ,dzt       &  ! Delta Z at M,T point
     ,dzim  ,dzit      &  ! Delta Z inverse (1/dz) at M,T point
     
     ,zfacm ,zfact     &  ! expansion factor of delta_x with height at M,T point
     ,zfacim,zfacit    &  ! inverse of zfacm, zfact

     ,xem ,yem, zem    &  ! XE,YE,ZE coordinates of M point
     ,xeu ,yeu, zeu    &  ! XE,YE,ZE coordinates of U point
     ,xev ,yev, zev    &  ! XE,YE,ZE coordinates of V point
     ,xew ,yew, zew    &  ! XE,YE,ZE coordinates of W point

     ,unx, uny, unz    &  ! U face normal unit vector components
     ,vnx, vny, vnz    &  ! V face normal unit vector components
     ,wnx, wny, wnz    &  ! W face normal unit vector components

     ,dnu              &  ! dxy across U face
     ,dniu             &  ! (1/dxy) across U face for turb U XY flx

     ,dnv              &  ! dxy across V face for PGF
     ,dniv             &  ! (1/dxy) across V face for turb V XY flx

     ,arm0             &  ! Area of IM triangle at earth surface
     ,arw0             &  ! Area of IW polygon at earth surface

     ,glatm ,glonm     &  ! Latitude/longitude at M point
     ,glatu ,glonu     &  ! Latitude/longitude at U point
     ,glatv ,glonv     &  ! Latitude/longitude at V point
     ,glatw ,glonw     &  ! Latitude/longitude at W point

     ,topm, topw          ! Topography height at M,W point

   real, allocatable, dimension(:,:) ::  &

      aru, arv, arw        &  ! Aperture area of U,V,W face
     ,volui, volvi ,volwi     ! (1/Volume) of U,V,W cell

   real(kind=8), allocatable, dimension(:,:) ::  &

      volt                 &  ! Volume of T cell
     ,volti                   ! (1/Volume) of T cell

   integer :: impent(12) ! Scratch array for storing 12 pentagonal IM indices

   integer, parameter :: nrows = 5
   integer :: mrows

Contains

!===============================================================================

   subroutine alloc_gridz()

   implicit none

   allocate (zm    (mza));  zm    (:) = 0.
   allocate (zt    (mza));  zt    (:) = 0.
   allocate (dzm   (mza));  dzm   (:) = 0.
   allocate (dzt   (mza));  dzt   (:) = 0.
   allocate (dzim  (mza));  dzim  (:) = 0.
   allocate (dzit  (mza));  dzit  (:) = 0.
   allocate (zfacm (mza));  zfacm (:) = 0.
   allocate (zfacim(mza));  zfacim(:) = 0.
   allocate (zfact (mza));  zfact (:) = 0.
   allocate (zfacit(mza));  zfacit(:) = 0.

   return
   end subroutine alloc_gridz
   
!===============================================================================

   subroutine alloc_xyzem(lma)

   implicit none
   
   integer, intent(in) :: lma

   allocate (xem(lma));  xem(1:lma) = 0.
   allocate (yem(lma));  yem(1:lma) = 0.
   allocate (zem(lma));  zem(1:lma) = 0.

   return
   end subroutine alloc_xyzem
   
!===============================================================================

   subroutine alloc_xyzew(lwa)

   implicit none

   integer, intent(in) :: lwa

   allocate (xew(lwa));  xew(1:lwa) = 0.
   allocate (yew(lwa));  yew(1:lwa) = 0.
   allocate (zew(lwa));  zew(1:lwa) = 0.

   return
   end subroutine alloc_xyzew
   
!===============================================================================

   subroutine alloc_grid1(meshtype, lma, lua, lva, lwa)

   use misc_coms, only: io6

   implicit none
   
   integer, intent(in) :: meshtype, lma, lua, lva, lwa
   
! Allocate and initialize arrays (xem, yem, zem are already allocated)

   allocate (lsw(lwa));  lsw(1:lwa) = 0

   if (meshtype == 1) then

      allocate (xeu(lua));  xeu(1:lua) = 0.
      allocate (yeu(lua));  yeu(1:lua) = 0.
      allocate (zeu(lua));  zeu(1:lua) = 0.

      allocate (glatu(lua));  glatu(1:lua) = 0.
      allocate (glonu(lua));  glonu(1:lua) = 0.

   elseif (meshtype == 2) then
   
      allocate (xev(lva));  xev(1:lva) = 0.
      allocate (yev(lva));  yev(1:lva) = 0.
      allocate (zev(lva));  zev(1:lva) = 0.

      allocate (glatv(lva));  glatv(1:lva) = 0.
      allocate (glonv(lva));  glonv(1:lva) = 0.
      
   endif

   allocate (unx(lva));  unx(1:lva) = 0.
   allocate (uny(lva));  uny(1:lva) = 0.
   allocate (unz(lva));  unz(1:lva) = 0.

   allocate (vnx(lva));  vnx(1:lva) = 0.
   allocate (vny(lva));  vny(1:lva) = 0.
   allocate (vnz(lva));  vnz(1:lva) = 0.

   allocate (wnx(lwa));  wnx(1:lwa) = 0.
   allocate (wny(lwa));  wny(1:lwa) = 0.
   allocate (wnz(lwa));  wnz(1:lwa) = 0.

   allocate (dnu  (lva));  dnu (1:lva) = 0.
   allocate (dniu (lva));  dniu(1:lva) = 0.

   allocate (dnv  (lva));  dnv (1:lva) = 0.
   allocate (dniv (lva));  dniv(1:lva) = 0.

   allocate  (arw0(lwa));   arw0(1:lwa) = 0.
   allocate  (topw(lwa));   topw(1:lwa) = 0.
   allocate (glatw(lwa));  glatw(1:lwa) = 0.
   allocate (glonw(lwa));  glonw(1:lwa) = 0.

   allocate  (arm0(lma));   arm0(1:lma) = 0.
   allocate  (topm(lma));   topm(1:lma) = 0.
   allocate (glatm(lma));  glatm(1:lma) = 0.
   allocate (glonm(lma));  glonm(1:lma) = 0.

   write(io6,*) 'finishing alloc_grid1'
            
   return
  
   end subroutine alloc_grid1

!===============================================================================

   subroutine alloc_grid2(meshtype, lma, lua, lva, lwa)

   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: meshtype, lma, lua, lva, lwa
   
! Allocate  and initialize arrays

   write(io6,*) 'alloc_grid2 ',lma,lua,lva,lwa

   if (meshtype == 1) then

      allocate (lpu(lua));  lpu(1:lua) = 0
      allocate (lcu(lua));  lcu(1:lua) = 0

      allocate (aru  (mza,lua));  aru  (1:mza,1:lua) = 0.
      allocate (volui(mza,lua));  volui(1:mza,1:lua) = 0.

   elseif (meshtype == 2) then

      allocate (lpv(lva));  lpv(1:lva) = 0
      allocate (lcv(lva));  lcv(1:lva) = 0

      allocate (aru  (mza,lva));  aru  (1:mza,1:lva) = 0.
      allocate (arv  (mza,lva));  arv  (1:mza,1:lva) = 0.
      allocate (volvi(mza,lva));  volvi(1:mza,1:lva) = 0.

   endif

   allocate (lpm(lma));  lpm(1:lma) = 0  ! In vtables
   allocate (lpw(lwa));  lpw(1:lwa) = 0  ! In vtables

   allocate (arw  (mza,lwa));  arw  (1:mza,1:lwa) = 0.

   allocate (volt (mza,lwa));  volt (1:mza,1:lwa) = 0.
   allocate (volti(mza,lwa));  volti(1:mza,1:lwa) = 0.
   allocate (volwi(mza,lwa));  volwi(1:mza,1:lwa) = 0.
            
   write(io6,*) 'finishing alloc_grid2'
            
   return
  
   end subroutine alloc_grid2

End Module mem_grid

