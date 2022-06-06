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
     ,nve2_max      &  ! Max # underground v[xyz]e2 levels in IW column
     ,nma,nua,nva,nwa  ! Horiz number of all M,U,V,W points in full domain

   integer :: &

      mza           &  ! Vertical number of all points
     ,mma,mua,mva,mwa  ! Horiz number of all M,U,V,W points in sub-process

   integer, allocatable, dimension(:) ::  &

      lpm,lpv,lpw   &  ! Lowest prognosed M,V,W/T
     ,lsw           &  ! number of W/T levels in contact with surface
     ,lve2             ! number of v[xyz]e2 levels in contact with surface

   real, allocatable, dimension(:) ::  &

      zm    ,zt        &  ! Z coordinate at M,T point
     ,dzm   ,dzt       &  ! Delta Z at M,T point
     ,dzim  ,dzit      &  ! Delta Z inverse (1/dz) at M,T point

     ,zfacm ,zfact     &  ! expansion factor of delta_x with height at M,T point
     ,zfacim,zfacit    &  ! inverse of zfacm, zfact
     ,zfacm2,zfacim2   &  ! expansion factor of arw with height, and its inverse

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

   real, allocatable :: gravm(:) ! gravity at M levels
   real, allocatable :: gravt(:) ! gravity at T levels

   real, allocatable, dimension(:,:) ::  &

      arv, arw            ! Aperture area of V,W face

   real(r8), allocatable, dimension(:,:) ::  &

      volt                ! Volume of T cell

   real, allocatable, dimension(:,:) :: &
        volti, &          ! 1 / (Volume of T cell)
        volwi, &          ! 1 / (Volumes of adjacent T cells)
        volvi             ! 1 / (Volumes of adjacent T cells)

   integer :: impent(12)  ! Scratch array for storing 12 pentagonal IM indices

   integer, parameter :: nrows = 5
   integer :: mrows

! "Derived" variables computed at beginning of integration rather than
! at the MAKEGRID stage

   real, allocatable, dimension(:) ::  &

        wnxo2, wnyo2, wnzo2,  & ! W-face unit normals divided by 2
        vnxo2, vnyo2, vnzo2,  & ! V-face unit normals divided by 2
        unxo2, unyo2, unzo2,  & ! U-face unit normals divided by 2

        dzt_top,              & ! distance between ZM(k) and ZT(k)
        dzt_bot,              & ! distance between ZT(k) and ZM(k-1)

        dzit_top,             & ! distance between ZM(k) and ZT(k)
        dzit_bot,             & ! distance between ZT(k) and ZM(k-1)

        zwgt_top, zwgt_bot,   & ! weights for interpolating T levels to W
        dzto2,    dzto4,      & ! dzt(k)    / 2, dzt(k)    / 4
        dztsqo2,  dztsqo4,    & ! dzt(k)**2 / 2, dzt(k)**2 / 4
        dztsqo6,  dztsqo12,   & ! dzt(k)**2 / 6, dzt(k)**2 / 12
        dzimsq,               & ! dzim(k)**2

        voa0,                 & ! ratio of cell volume to bottom area w/o terrain

        gdz_belo, gdz_abov,   & ! weights for hydrostatic integration

        gdzim,                & ! gravm / dzm

        arw0i,                & ! 1 / arw0
        dnivo2                  ! 1/(2dxy) across V face

   ! double precision weights for interpolating T levels to W

   real(r8), allocatable :: zwgt_top8(:), zwgt_bot8(:)
   real(r8), allocatable :: gdz_belo8(:), gdz_abov8(:)
   real(r8), allocatable :: gdz_wgtm8(:), gdz_wgtp8(:)
   real,     allocatable :: gdz_wgtm (:), gdz_wgtp (:)

   real, allocatable, dimension(:,:) ::  &

        gxps_coef, gyps_coef    ! combined weights for grad_t2d

   real, allocatable, dimension(:) ::  &

        vxn_ew, vyn_ew, vzn_ew, & ! unit normals of zonal (east-west) direction in earth cartesian coordinates
        vxn_ns, vyn_ns, vzn_ns, & ! unit normals of meridional (north-south) direction in earth cartesian coordinates
        vcn_ew, vcn_ns            ! components of zonal and merdional vectors in the direction of VC

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
    allocate (zfacm2 (mza)); zfacm2 (:) = 0.
    allocate (zfacim2(mza)); zfacim2(:) = 0.

    allocate (gravm (mza));  gravm (:) = 0.
    allocate (gravt (mza));  gravt (:) = 0.

  end subroutine alloc_gridz

!===============================================================================

  subroutine alloc_xyzem(lma)

    implicit none

    integer, intent(in) :: lma

    allocate (xem(lma));  xem(1:lma) = 0.
    allocate (yem(lma));  yem(1:lma) = 0.
    allocate (zem(lma));  zem(1:lma) = 0.

  end subroutine alloc_xyzem

!===============================================================================

  subroutine alloc_xyzew(lwa)

    implicit none

    integer, intent(in) :: lwa

    allocate (xew(lwa));  xew(1:lwa) = 0.
    allocate (yew(lwa));  yew(1:lwa) = 0.
    allocate (zew(lwa));  zew(1:lwa) = 0.

  end subroutine alloc_xyzew

!===============================================================================

  subroutine alloc_grid1(lma, lva, lwa)

    use misc_coms, only: io6
    implicit none

    integer, intent(in) :: lma, lva, lwa

! Allocate and initialize arrays (xem, yem, zem are already allocated)

    allocate (lsw (lwa));  lsw (1:lwa) = 0
    allocate (lve2(lwa));  lve2(1:lwa) = 0

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

    allocate (vxn_ew(lwa)) ; vxn_ew = 0.
    allocate (vyn_ew(lwa)) ; vyn_ew = 0.
    allocate (vzn_ew(lwa)) ; vzn_ew = 0.

    allocate (vxn_ns(lwa)) ; vxn_ns = 0.
    allocate (vyn_ns(lwa)) ; vyn_ns = 0.
    allocate (vzn_ns(lwa)) ; vzn_ns = 0.

    allocate (vcn_ew(lva)) ; vcn_ew = 0.
    allocate (vcn_ns(lva)) ; vcn_ns = 0.

    write(io6,*) 'finishing alloc_grid1'

  end subroutine alloc_grid1

!===============================================================================

   subroutine alloc_grid2(lma, lva, lwa)

     use consts_coms, only: r8
     use misc_coms,   only: io6

     implicit none

     integer, intent(in) :: lma, lva, lwa

! Allocate  and initialize arrays

     write(io6,*) 'alloc_grid2 ',lma,lva,lwa

     allocate (lpv(lva)); lpv(1:lva) = 0
     allocate (lpm(lma)); lpm(1:lma) = 0
     allocate (lpw(lwa)); lpw(1:lwa) = 0

     allocate (arv  (mza,lva));  arv  (1:mza,1:lva) = 0.
     allocate (arw  (mza,lwa));  arw  (1:mza,1:lwa) = 0.
     allocate (volt (mza,lwa));  volt (1:mza,1:lwa) = 0._r8

     write(io6,*) 'finishing alloc_grid2'

   end subroutine alloc_grid2

!===============================================================================

   subroutine alloc_gridz_other()

     ! This routine allocates and defines vertical grid arrays that were
     ! not computed during the MAKEGRID stage

     implicit none

     integer :: k

     ! Allocate and define 1D variables defined at T levels

     allocate(dzt_top(mza))
     allocate(dzt_bot(mza))

     dzt_top(1) = zm(1) - zt(1)
     dzt_bot(1) = dzt(1) - dzt_top(1)

     ! Loop over T levels
     do k = 2, mza
        dzt_top(k) = zm(k) - zt(k)
        dzt_bot(k) = zt(k) - zm(k-1)
     enddo

     allocate(dzit_top(mza))
     allocate(dzit_bot(mza))

     do k = 1, mza
        dzit_top(k) = 1.0 / dzt_top(k)
        dzit_bot(k) = 1.0 / dzt_bot(k)
     enddo

     ! Allocate and define 1D variables defined at W levels

     allocate(gdzim(mza))
     gdzim = gravm * dzim

     ! weights for averaging T variables to W levels
     allocate(zwgt_top(mza))
     allocate(zwgt_bot(mza))

     allocate(zwgt_top8(mza))
     allocate(zwgt_bot8(mza))

     ! weights for hydrostatic integration at W level

     allocate(gdz_belo(mza))
     allocate(gdz_abov(mza))

     allocate(gdz_belo8(mza))
     allocate(gdz_abov8(mza))

     allocate(gdz_wgtm8(mza))
     allocate(gdz_wgtp8(mza))

     allocate(gdz_wgtm(mza))
     allocate(gdz_wgtp(mza))

     ! Loop over W levels
     do k = 1, mza-1
        zwgt_top8(k) = dzt_bot(k+1) * dzim(k)
        zwgt_bot8(k) = dzt_top(k)   * dzim(k)

        gdz_abov8(k) = dzt_bot(k+1) * gravm(k)
        gdz_belo8(k) = dzt_top(k)   * gravm(k)
     enddo

     zwgt_top8(mza) = zwgt_top8(mza-1)
     zwgt_bot8(mza) = zwgt_bot8(mza-1)

     gdz_abov8(mza) = gdz_abov8(mza-1)
     gdz_belo8(mza) = gdz_belo8(mza-1)

     do k = 1, mza
        gdz_wgtp8(k) = zwgt_top8(k) * gravm(k)
        gdz_wgtm8(k) = zwgt_bot8(k) * gravm(k)

        gdz_wgtp(k) = real(gdz_wgtp8(k))
        gdz_wgtm(k) = real(gdz_wgtm8(k))

        gdz_abov(k) = real(gdz_abov8(k))
        gdz_belo(k) = real(gdz_belo8(k))

        zwgt_top(k) = real(zwgt_top8(k))
        zwgt_bot(k) = real(zwgt_bot8(k))
     enddo

     allocate(dzto2   (mza))
     allocate(dzto4   (mza))
     allocate(dztsqo2 (mza))
     allocate(dztsqo4 (mza))
     allocate(dztsqo6 (mza))
     allocate(dztsqo12(mza))
     allocate(dzimsq  (mza))
     allocate(voa0    (mza))

     do k = 1, mza
        dzto2   (k) = dzt    (k) * 0.50
        dzto4   (k) = dzt    (k) * 0.25
        dztsqo2 (k) = dzto2  (k) * dzt(k)
        dztsqo4 (k) = dztsqo2(k) * 0.5
        dztsqo6 (k) = dztsqo2(k) / 3.0
        dztsqo12(k) = dztsqo6(k) * 0.5
        dzimsq  (k) = dzim   (k) * dzim(k)
     enddo

     do k = 2, mza
        voa0(k) = dzt(k) * zfact(k)**2 * zfacim2(k-1)
     enddo
     voa0(1) = voa0(2)

   end subroutine alloc_gridz_other

!===============================================================================

   subroutine alloc_grid_other()

     ! This routine allocates and defines grid arrays that were not computed
     ! during the MAKEGRID stage

     use consts_coms, only: r8, eradi
     use mem_ijtabs,  only: itab_w, itab_v, jtab_v, jtv_wadj
     use misc_coms,   only: mdomain

     implicit none

     integer :: iw, j, iv, iw1, iw2, k, n1, n2, kbot
     real    :: raxis, vxn_ewv, vyn_ewv, vxn_nsv, vyn_nsv, vzn_nsv

     ! Allocate and define variables defined at V faces

     allocate(vnxo2     (mva))
     allocate(vnyo2     (mva))
     allocate(vnzo2     (mva))
     allocate(unxo2     (mva))
     allocate(unyo2     (mva))
     allocate(unzo2     (mva))
     allocate(dnivo2    (mva))
     allocate(volvi (mza,mva))

     vnxo2  (1) = 0.0
     vnyo2  (1) = 0.0
     vnzo2  (1) = 0.0
     unxo2  (1) = 0.0
     unyo2  (1) = 0.0
     unzo2  (1) = 0.0
     dnivo2 (1) = 0.0
     volvi(:,1) = 0.0

     allocate(wnxo2    (mwa))
     allocate(wnyo2    (mwa))
     allocate(wnzo2    (mwa))
     allocate(arw0i    (mwa))
     allocate(volti(mza,mwa))
     allocate(volwi(mza,mwa))
     allocate(gxps_coef(7,mwa))
     allocate(gyps_coef(7,mwa))

     wnxo2(1) = 0.0
     wnyo2(1) = 0.0
     wnzo2(1) = 0.0

     volti(:,1) = 0.0
     volwi(:,1) = 0.0

     !$omp parallel
     !$omp do
     do iv = 2, mva
        vnxo2 (iv) = vnx (iv) * 0.5
        vnyo2 (iv) = vny (iv) * 0.5
        vnzo2 (iv) = vnz (iv) * 0.5
        unxo2 (iv) = unx (iv) * 0.5
        unyo2 (iv) = uny (iv) * 0.5
        unzo2 (iv) = unz (iv) * 0.5
        dnivo2(iv) = dniv(iv) * 0.5
     enddo
     !omp end do

     !$omp do private(iv,iw1,iw2,k,vxn_ewv,vyn_ewv,vxn_nsv,vyn_nsv,vzn_nsv,kbot)
     do j = 1,jtab_v(jtv_wadj)%jend(1); iv = jtab_v(jtv_wadj)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

!       kbot = max(lpw(iw1),lpw(iw2))
        kbot = min(lpw(iw1),lpw(iw2))

!       do k = lpv(iv), mza
        do k = kbot, mza
           volvi(k,iv) = real( 1.0_r8 / (volt(k,iw1) + volt(k,iw2)) )
        enddo
        volvi(1:kbot-1,iv) = 0.0

        if (mdomain > 1) then

           vcn_ew(iv) = vnx(iv)
           vcn_ns(iv) = vny(iv)

        else

           raxis = sqrt(xev(iv)**2 + yev(iv)**2)
           if (raxis > 1.e3) then

              vxn_ewv = -yev(iv) / raxis
              vyn_ewv =  xev(iv) / raxis

              vxn_nsv = -xev(iv) * zev(iv) * eradi / raxis
              vyn_nsv = -yev(iv) * zev(iv) * eradi / raxis
              vzn_nsv =  raxis * eradi

              vcn_ew(iv) = vxn_ewv * vnx(iv) + vyn_ewv * vny(iv)
              vcn_ns(iv) = vxn_nsv * vnx(iv) + vyn_nsv * vny(iv) + vzn_nsv * vnz(iv)

           else

              vcn_ew(iv) = 0.
              vcn_ns(iv) = 0.

           endif
        endif

     enddo
     !$omp end do nowait

     !$omp do private(n1,n2,raxis)
     do iw = 2, mwa

        wnxo2(iw) = wnx(iw) * 0.5
        wnyo2(iw) = wny(iw) * 0.5
        wnzo2(iw) = wnz(iw) * 0.5

        volti(1:lpw(iw)-1,iw) = 0.0
        volti(lpw(iw):mza,iw) = real( 1.0_r8 / volt(lpw(iw):mza,iw) )

        volwi(1:lpw(iw)-1,iw) = 0.0
        volwi(mza,iw)         = 0.0

        volwi(lpw(iw):mza-1,iw) = real( 1.0_r8 / ( volt(lpw(iw)  :mza-1,iw) &
                                                 + volt(lpw(iw)+1:mza  ,iw) ) )
        arw0i(iw) = 1.0 / arw0(iw)

        do n1 = 1, itab_w(iw)%npoly
           if (n1 == 1) then
              n2 = itab_w(iw)%npoly
           else
              n2 = n1 - 1
           endif

           gxps_coef(n1,iw) = itab_w(iw)%gxps1(n1) + itab_w(iw)%gxps2(n2)
           gyps_coef(n2,iw) = itab_w(iw)%gyps1(n1) + itab_w(iw)%gyps2(n2)
        enddo

        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

        if (mdomain < 2 .and. raxis > 1.e3) then

           vxn_ew(iw) = -yew(iw) / raxis
           vyn_ew(iw) =  xew(iw) / raxis
           vzn_ew(iw) =  0.0

           vxn_ns(iw) = -xew(iw) * zew(iw) * eradi / raxis
           vyn_ns(iw) = -yew(iw) * zew(iw) * eradi / raxis
           vzn_ns(iw) =  raxis * eradi

        else
           vxn_ew(iw) = 1.0
           vyn_ew(iw) = 0.0
           vzn_ew(iw) = 0.0

           vxn_ns(iw) = 0.0
           vyn_ns(iw) = 1.0
           vzn_ns(iw) = 0.0
        endif

     enddo
     !$omp end do
     !$omp end parallel

   end subroutine alloc_grid_other

End Module mem_grid
