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
Module mem_grid

   use consts_coms, only: r8
   implicit none

   private :: r8

   integer :: &  ! N values for full-domain reference on any process

      nza           &  ! Vertical number of all points
     ,nsw_max       &  ! Max # vert atm levels in IW column with sfc flux
     ,nma,nua,nva,nwa  ! Horiz number of all M,U,V,W points in full domain

   integer :: &

      mza           &  ! Vertical number of all points
     ,mma,mua,mva,mwa  ! Horiz number of all M,U,V,W points in sub-process

   integer, allocatable, dimension(:) ::  &

      lpm,lpv,lpw   &  ! Lowest prognosed M,V,W/T
     ,lsw              ! number of W/T levels in contact with surface

   real, allocatable, dimension(:) ::  &

      zm    ,zt        &  ! Z coordinate at M,T point
     ,dzm   ,dzt       &  ! Delta Z at M,T point
     ,dzim  ,dzit      &  ! Delta Z inverse (1/dz) at M,T point
     
     ,zfacm ,zfact     &  ! expansion factor of delta_x with height at M,T point
     ,zfacim,zfacit    &  ! inverse of zfacm, zfact

     ,xem ,yem, zem    &  ! XE,YE,ZE coordinates of M point
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
     ,glatv ,glonv     &  ! Latitude/longitude at V point
     ,glatw ,glonw     &  ! Latitude/longitude at W point

     ,topm, topw          ! Topography height at M,W point

   real, allocatable, dimension(:,:) ::  &

      arv, arw            ! Aperture area of V,W face

   real(r8), allocatable, dimension(:,:) ::  &

      volt, volti         ! Volume of T cell; (1/Volume) of T cell

   integer :: impent(12)  ! Scratch array for storing 12 pentagonal IM indices

   integer, parameter :: nrows = 5
   integer :: mrows

! "Derived" variables computed at beginning of integration rather than
! at the MAKEGRID stage

   real, allocatable, dimension(:) ::  &

        wnxo2, wnyo2, wnzo2,  & ! W-face unit normals divided by 2
        vnxo2, vnyo2, vnzo2,  & ! V-face unit normals divided by 2
        
        dzt_top,              & ! distance between ZM(k) and ZT(k)
        dzt_bot,              & ! distance between ZT(k) and ZM(k-1)

        zwgt_top, zwgt_bot      ! weights for interpolating T levels to W

   real, allocatable, dimension(:,:) ::  &

        gxps_coef, gyps_coef    ! combined weights for grad_t3d

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

   subroutine alloc_grid1(lma, lva, lwa)

   use misc_coms, only: io6

   implicit none
   
   integer, intent(in) :: lma, lva, lwa
   
! Allocate and initialize arrays (xem, yem, zem are already allocated)

   allocate (lsw(lwa));  lsw(1:lwa) = 0

   allocate (xev(lva));  xev(1:lva) = 0.
   allocate (yev(lva));  yev(1:lva) = 0.
   allocate (zev(lva));  zev(1:lva) = 0.

   allocate (glatv(lva));  glatv(1:lva) = 0.
   allocate (glonv(lva));  glonv(1:lva) = 0.
      
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

   subroutine alloc_grid2(lma, lva, lwa)

   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: lma, lva, lwa
   
! Allocate  and initialize arrays

   write(io6,*) 'alloc_grid2 ',lma,lva,lwa

   allocate (lpv(lva));  lpv(1:lva) = 0

   allocate (arv  (mza,lva));  arv  (1:mza,1:lva) = 0.

   allocate (lpm(lma));  lpm(1:lma) = 0  ! In vtables
   allocate (lpw(lwa));  lpw(1:lwa) = 0  ! In vtables

   allocate (arw  (mza,lwa));  arw  (1:mza,1:lwa) = 0.

   allocate (volt (mza,lwa));  volt (1:mza,1:lwa) = 0.
   allocate (volti(mza,lwa));  volti(1:mza,1:lwa) = 0.
            
   write(io6,*) 'finishing alloc_grid2'
            
   return
  
   end subroutine alloc_grid2

!===============================================================================

   subroutine alloc_grid_other()
     
     use consts_coms, only: r8
     use mem_ijtabs,  only: itab_w
     implicit none
     
     integer :: iw, iv, k, n1, n2

     ! This routine allocates and defines grid arrays that were not computed
     ! during the MAKEGRID stage
     
     allocate(dzt_top(mza))
     allocate(dzt_bot(mza))

     dzt_top(1) = zm(1) - zt(1)
     dzt_bot(1) = dzt(1) - dzt_top(1)

     ! Loop over T levels
     do k = 2, mza
        dzt_top(k) = zm(k) - zt(k)
        dzt_bot(k) = zt(k) - zm(k-1)
     enddo

     allocate(zwgt_top(mza))
     allocate(zwgt_bot(mza))

     ! Loop over W levels
     do k = 1, mza-1
        zwgt_top(k) = dzt_top(k)   * dzim(k)
        zwgt_bot(k) = dzt_bot(k+1) * dzim(k)
     enddo
     
     zwgt_top(mza) = zwgt_top(mza-1)
     zwgt_bot(mza) = zwgt_bot(mza-1)

     allocate(wnxo2(mwa))
     allocate(wnyo2(mwa))
     allocate(wnzo2(mwa))
     
     wnxo2(1) = 0.0
     wnyo2(1) = 0.0
     wnzo2(1) = 0.0

     !$omp parallel do
     do iw = 2, mwa
        wnxo2(iw) = wnx(iw) * 0.5
        wnyo2(iw) = wny(iw) * 0.5
        wnzo2(iw) = wnz(iw) * 0.5
     enddo
     !$omp end parallel do

     allocate(vnxo2(mva))
     allocate(vnyo2(mva))
     allocate(vnzo2(mva))

     vnxo2(1) = 0.0
     vnyo2(1) = 0.0
     vnzo2(1) = 0.0

     !$omp parallel do
     do iv = 2, mva
        vnxo2(iv) = vnx(iv) * 0.5
        vnyo2(iv) = vny(iv) * 0.5
        vnzo2(iv) = vnz(iv) * 0.5
     enddo
     !$omp end parallel do

     allocate(gxps_coef(mwa,7))
     allocate(gyps_coef(mwa,7))

     !$omp parallel do
     do iw = 2, mwa
        do n1 = 1, itab_w(iw)%npoly

           if (n1 == 1) then
              n2 = itab_w(iw)%npoly
           else
              n2 = n1 - 1
           endif

           gxps_coef(iw,n1) = itab_w(iw)%gxps1(n1) + itab_w(iw)%gxps2(n2)
           gyps_coef(iw,n1) = itab_w(iw)%gyps1(n1) + itab_w(iw)%gyps2(n2)

        enddo
     enddo

   end subroutine alloc_grid_other

End Module mem_grid
