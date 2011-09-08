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

subroutine makesfc()

! This subroutine generates LAND and SEA files for OLAM runs. 

  use mem_grid,    only: nza, nma, nua, nwa, xem, yem, zem, xew, yew, zew,  &
                         arw0, glatw, glonw, arm0, glatm, glonm, zm, topm

  use misc_coms,   only: io6, runtype, itopoflg, rinit

  use sea_coms,    only: nms, nus, nws, seafile

  use leaf_coms,   only: nml, nul, nwl, nzg, nzs, nslcon, nvgcon,  &
                         isfcl, ivegflg, landusefile,  &
                         isoilflg, soil_database, veg_database

  use mem_ijtabs,  only: itab_md, itab_ud, itab_wd

  use consts_coms, only: erad, piu180

  use leaf_db,     only: leaf_database_read

  use max_dims,    only: maxnlspoly

  use mem_mksfc,   only: itab_uls_vars, itab_wls_vars

  use mem_leaf,    only: itab_ul, itab_wl, land, alloc_land_grid

  use mem_sea,     only: itab_us, itab_ws, sea, alloc_sea_grid

  implicit none

  integer :: k
  integer :: ndims, idims(1)
  integer :: im, iu, iw, iw1, iw2
  integer :: ius, iul, npoly
  integer :: imin, imax, kmin, kmax, nku
  integer :: im1ls, im2ls, iw1ls, iw2ls
  integer :: nmls, nuls, nwls
  integer :: iuls, iwls, imls, jmls
  integer :: iws, ims
  integer :: iwl, iml
  integer :: idatq, datsoil
  integer :: nzg0
  integer :: nws0
  integer :: nwl0
  integer :: nmls0
  integer :: ivegflg0, isoilflg0
  integer :: im1,im2,im3
  integer :: ju, nkw, iu1,iu2,iu3
  integer :: iuls_lo, iuls_hi, iu1ls, iu2ls, iu3ls
  integer :: iu1ls_mlo,iu1ls_mhi,iu2ls_mlo,iu2ls_mhi,iu3ls_mlo,iu3ls_mhi
  integer :: kku, kku1, kku2, kku3, kkw, kkw1, kkw2
  integer :: im1k, im2k, im3k

  real :: expansion, raxis
  real :: s1,s2,s3,hper
  real :: fracdist, fracz, ef
  real :: baryx, baryy, baryx1, baryx2, baryx3, baryy1, baryy2, baryy3
  real :: area, area1, area2, area3
  real :: wnx1, wny1, wnz1
  real :: wnx2, wny2, wnz2
  real :: wnx3, wny3, wnz3
  real :: xp(maxnlspoly), yp(maxnlspoly)
  real :: xemsfc(maxnlspoly), yemsfc(maxnlspoly), zemsfc(maxnlspoly)

! automatic arrays

  integer :: k1topm(nma),k2topm(nma)

  integer :: kmin_u(nua)
  integer :: kmax_u(nua)
  integer :: kadd_u(nua)

  integer :: kmin_w(nwa)
  integer :: kmax_w(nwa)
  integer :: kadd_w(nwa)

  integer :: iuls0_u(nua)
  integer :: imls0_u(nua)

  integer :: iuls0_w(nwa)
  integer :: iwls0_w(nwa)

  type(itab_uls_vars), allocatable :: mksfc_itab_uls(:)
  type(itab_wls_vars), allocatable :: mksfc_itab_wls(:)

! LOCAL ARRAYS

  integer, allocatable :: mpt_sea   (:) ! flag/counter of M pt for sea cell
  integer, allocatable :: mpt_land  (:) ! flag/counter of M pt for land cell
  integer, allocatable :: upt_sea   (:) ! flag/counter of U pt for sea cell
  integer, allocatable :: upt_land  (:) ! flag/counter of U pt for land cell
  integer, allocatable :: wpt_sea   (:) ! Counter of W pt for sea cell
  integer, allocatable :: wpt_land  (:) ! Counter of W pt for land cell
  integer, allocatable :: idatp     (:) ! integer array for storing database data

  integer, allocatable :: leaf_class(:) ! leaf (vegetation) class
  real,    allocatable :: areals    (:) ! land/sea cell area [m^2]
  real,    allocatable :: xewls     (:) ! earth x coord of land/sea grid W pts [m]
  real,    allocatable :: yewls     (:) ! earth y coord of land/sea grid W pts [m]
  real,    allocatable :: zewls     (:) ! earth z coord of land/sea grid W pts [m]
  real,    allocatable :: glatwls   (:) ! latitude of sfc cell W points
  real,    allocatable :: glonwls   (:) ! longitude of sfc cell W points
  real,    allocatable :: wnxls     (:) ! norm unit vec x comp of land/sea cell [m] 
  real,    allocatable :: wnyls     (:) ! norm unit vec y comp of land/sea cell [m] 
  real,    allocatable :: wnzls     (:) ! norm unit vec z comp of land/sea cell [m] 

  real,    allocatable :: xemls     (:) ! earth x coord of land/sea grid M pts [m]
  real,    allocatable :: yemls     (:) ! earth y coord of land/sea grid M pts [m]
  real,    allocatable :: zemls     (:) ! earth z coord of land/sea grid M pts [m]
  real,    allocatable :: zmls      (:) ! topo height of land/sea grid M pts [m]
  real,    allocatable :: glatmls   (:) ! latitude of sfc cell W points
  real,    allocatable :: glonmls   (:) ! longitude of sfc cell W points

! Subroutine makesfc is called from subroutine gridinit after OLAM has 
! initialized its Delaunay grid but before initializing its Voronoi grid.

! Allocate and fill XEW, YEW, ZEW, GLATW, GLONW, ARM0, GLATM, GLONM, TOPM 
! for Delaunay grid since they are not yet initialized 
! (and will not otherwise be used on the Delaunay grid)

! (ASSUME THAT MDOMAIN = 0.)

  allocate (xew  (nwa)) ; xew   = rinit
  allocate (yew  (nwa)) ; yew   = rinit
  allocate (zew  (nwa)) ; zew   = rinit
  allocate (glatw(nwa)) ; glatw = rinit
  allocate (glonw(nwa)) ; glonw = rinit
  allocate (arm0 (nma)) ; arm0  = rinit
  allocate (glatm(nma)) ; glatm = rinit
  allocate (glonm(nma)) ; glonm = rinit
  allocate (topm (nma)) ; topm  = rinit

! Initialize automatic arrays

  k1topm(:) = 0
  k2topm(:) = 0

  kmin_u(:) = 0
  kmax_u(:) = 0
  kadd_u(:) = 0

  kmin_w(:) = 0
  kmax_w(:) = 0
  kadd_w(:) = 0

  iuls0_u(:) = 0
  imls0_u(:) = 0

  iuls0_w(:) = 0
  iwls0_w(:) = 0

! Loop over all Delaunay W points

  do iw = 2,nwa

     im1 = itab_wd(iw)%im(1)
     im2 = itab_wd(iw)%im(2)
     im3 = itab_wd(iw)%im(3)

! OPTION 1: barycentric coordinates for IW point

     xew(iw) = (xem(im1) + xem(im2) + xem(im3)) / 3.
     yew(iw) = (yem(im1) + yem(im2) + yem(im3)) / 3.
     zew(iw) = (zem(im1) + zem(im2) + zem(im3)) / 3.

! Push W point coordinates out to earth radius

     expansion = erad / sqrt(xew(iw) ** 2  &
          + yew(iw) ** 2  &
          + zew(iw) ** 2  )

     xew(iw) = xew(iw) * expansion
     yew(iw) = yew(iw) * expansion
     zew(iw) = zew(iw) * expansion

     raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

     glatw(iw) = atan2(zew(iw),raxis)   * piu180
     glonw(iw) = atan2(yew(iw),xew(iw)) * piu180

  enddo

! Loop over all Delaunay M points

  do im = 2,nma

     raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
     glatm(im) = atan2(zem(im),raxis)   * piu180
     glonm(im) = atan2(yem(im),xem(im)) * piu180

! Set area of dual cell surrounding M to zero prior to summation over triangles

     arm0(im) = 0.

! Loop over all U neighbors of M

     do ju = 1,itab_md(im)%npoly

        iu = itab_md(im)%iu(ju)

        iw1 = itab_ud(iu)%iw(1)
        iw2 = itab_ud(iu)%iw(2)

! Contribution to dual-cell area around M point of triangle M-IW1-IW2

        s1 = sqrt( (xew(iw1) - xew(iw2))**2  &
             + (yew(iw1) - yew(iw2))**2  &
             + (zew(iw1) - zew(iw2))**2  )

        s2 = sqrt( (xem(im) - xew(iw2))**2  &
             + (yem(im) - yew(iw2))**2  &
             + (zem(im) - zew(iw2))**2  )

        s3 = sqrt( (xew(iw1) - xem(im))**2  &
             + (yew(iw1) - yem(im))**2  &
             + (zew(iw1) - zem(im))**2  )

        hper = .5 * (s1 + s2 + s3)  ! half perimeter of triangle

        arm0(im) = arm0(im)  &
             + sqrt(hper * (hper - s1) * (hper - s2) * (hper - s3))

     enddo

  enddo

! Read topography to triangle M points

  if (itopoflg == 1) then
! Read TOPM from dataset
     call topo_database_read(nma,xem,yem,zem,topm,arm0,glatm,glonm)
  else
! Initialize TOPM from default value (may customize here)
     call topo_init(nma,topm,glatm,glonm,xem,yem,zem)
  endif

! TOPM not allowed to be lower than zm(1)

  topm(2:nma) = max(topm(2:nma),zm(1))

! Loop over all Delaunay M points

  do im = 2,nma

! Loop over KM levels

     do k = 1,nza-2

! Check whether topm(im) is at or above current KM level

        if (zm(k) <= topm(im) .and. zm(k+1) >= topm(im)) then

           k1topm(im) = k
           k2topm(im) = k

! Check topm(im) within height range of KT level

           fracz = (topm(im) - zm(k)) / (zm(k+1) - zm(k))

! If topm(im) is in lowest or highest 1% of KT level, move it to the limit.      

           if (fracz < .01) then
              topm(im) = zm(k)
              k1topm(im) = k - 1
           elseif (fracz > .99) then
              topm(im) = zm(k+1)
              k2topm(im) = k + 1
           endif

           exit

        endif

     enddo

  enddo

! Initialize land/sea surface cell counters; use existing M points only

  nmls = nma
  nuls = 1
  nwls = 1

! Loop over all Delaunay U points

  do iu = 2,nua

     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

! Compute and store number of M and U points to be added along U edge

     kmin_u(iu) = min(k2topm(im1),k2topm(im2))
     kmax_u(iu) = max(k1topm(im1),k1topm(im2))

     if (kmax_u(iu) < kmin_u(iu)) kmax_u(iu) = kmin_u(iu)

     kadd_u(iu) = kmax_u(iu) - kmin_u(iu)

! Index of lowest IULS point and NULS increment for this U edge

     iuls0_u(iu) = nuls + 1
     nuls = nuls + kadd_u(iu) + 1

     if (kadd_u(iu) > 0) then

! Index of lowest IMLS point and NMLS increment for this U edge

        imls0_u(iu) = nmls + 1
        nmls = nmls + kadd_u(iu)

     endif

  enddo

! Loop over all Delaunay W points

  do iw = 2,nwa

     im1 = itab_wd(iw)%im(1)
     im2 = itab_wd(iw)%im(2)
     im3 = itab_wd(iw)%im(3)

! Compute and store number of U and W points to be added inside W triangle

     kmin_w(iw) = min(k2topm(im1),k2topm(im2),k2topm(im3))
     kmax_w(iw) = max(k1topm(im1),k1topm(im2),k1topm(im3))

     if (kmax_w(iw) < kmin_w(iw)) kmax_w(iw) = kmin_w(iw)

     kadd_w(iw) = kmax_w(iw) - kmin_w(iw)

! Index of lowest IWLS point and NWLS increment for this W triangle

     iwls0_w(iw) = nwls + 1
     nwls = nwls + kadd_w(iw) + 1

     if (kadd_w(iw) > 0) then

! Index of lowest IULS point and NULS increment for this W triangle

        iuls0_w(iw) = nuls + 1
        nuls = nuls + kadd_w(iw)

     endif

  enddo

! Allocate mksfc arrays for combined land/sea grid

  allocate (mksfc_itab_uls(nuls))
  allocate (mksfc_itab_wls(nwls))

! Allocate local arrays for combined land/sea grid

  allocate (mpt_sea   (nmls)) ; mpt_sea    = 0
  allocate (mpt_land  (nmls)) ; mpt_land   = 0
  allocate (upt_sea   (nuls)) ; upt_sea    = 0
  allocate (upt_land  (nuls)) ; upt_land   = 0
  allocate (wpt_sea   (nwls)) ; wpt_sea    = 0
  allocate (wpt_land  (nwls)) ; wpt_land   = 0
  allocate (idatp     (nwls)) ; idatp      = 0

  allocate (leaf_class(nwls)) ; leaf_class = 0
  allocate (areals    (nwls)) ; areals     = rinit
  allocate (xewls     (nwls)) ; xewls      = rinit
  allocate (yewls     (nwls)) ; yewls      = rinit
  allocate (zewls     (nwls)) ; zewls      = rinit
  allocate (glatwls   (nwls)) ; glatwls    = rinit
  allocate (glonwls   (nwls)) ; glonwls    = rinit
  allocate (wnxls     (nwls)) ; wnxls      = rinit
  allocate (wnyls     (nwls)) ; wnyls      = rinit
  allocate (wnzls     (nwls)) ; wnzls      = rinit

  allocate (xemls     (nmls)) ; xemls      = rinit
  allocate (yemls     (nmls)) ; yemls      = rinit
  allocate (zemls     (nmls)) ; zemls      = rinit
  allocate (zmls      (nmls)) ; zmls       = rinit
  allocate (glatmls   (nmls)) ; glatmls    = rinit
  allocate (glonmls   (nmls)) ; glonmls    = rinit

! Fill land/sea M point values that correspond with existing ATM M points

  xemls(1:nma) = xem(1:nma)
  yemls(1:nma) = yem(1:nma)
  zemls(1:nma) = zem(1:nma)
  zmls (1:nma) = topm(1:nma)

  do imls = 1, nma
     raxis = sqrt( xemls(imls)**2 + yemls(imls)**2 )
     glatmls(imls) = atan2( zemls(imls), raxis      ) * piu180
     glonmls(imls) = atan2( yemls(imls), xemls(imls)) * piu180
  enddo

! Loop over all Delaunay U points

  do iu = 2,nua

     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     if (topm(im1) <= topm(im2)) then
        imin = im1
        imax = im2
        iw1 = itab_ud(iu)%iw(1)
        iw2 = itab_ud(iu)%iw(2)
     else
        imin = im2
        imax = im1
        iw1 = itab_ud(iu)%iw(2)
        iw2 = itab_ud(iu)%iw(1)
     endif

! Loop over active K levels for this IU (those that intersect topography)

     do k = kmin_u(iu),kmax_u(iu)

        kku = k - kmin_u(iu)

        kkw1 = k - kmin_w(iw1)
        kkw2 = k - kmin_w(iw2)

        iuls = iuls0_u(iu) + kku

        iw1ls = iwls0_w(iw1) + kkw1
        iw2ls = iwls0_w(iw2) + kkw2

! default M neighbors of IULS

        mksfc_itab_uls(iuls)%im(1) = imin
        mksfc_itab_uls(iuls)%im(2) = imax

! corrections if more than 1 level

        if (k > kmin_u(iu)) then
           imls = imls0_u(iu) + kku - 1
           mksfc_itab_uls(iuls)%im(1) = imls
        endif

        if (k < kmax_u(iu)) then
           imls = imls0_u(iu) + kku
           mksfc_itab_uls(iuls)%im(2) = imls
        endif

        mksfc_itab_uls(iuls)%iw(1) = iw1ls
        mksfc_itab_uls(iuls)%iw(2) = iw2ls

        if (k > kmin_u(iu)) then
           imls = imls0_u(iu) + kku - 1

           fracdist = (zm(k) - topm(imin)) / (topm(imax) - topm(imin))

           xemls(imls) = xem(imin) + (xem(imax) - xem(imin)) * fracdist
           yemls(imls) = yem(imin) + (yem(imax) - yem(imin)) * fracdist
           zemls(imls) = zem(imin) + (zem(imax) - zem(imin)) * fracdist

           zmls(imls) = zm(k)

! Expand MLS point coordinates out to earth radius

           expansion = erad / sqrt(xemls(imls) ** 2  &
                + yemls(imls) ** 2  &
                + zemls(imls) ** 2  )

           xemls(imls) = xemls(imls) * expansion
           yemls(imls) = yemls(imls) * expansion
           zemls(imls) = zemls(imls) * expansion

           raxis = sqrt( xemls(imls)**2 + yemls(imls)**2 )
           glatmls(imls) = atan2( zemls(imls), raxis      ) * piu180
           glonmls(imls) = atan2( yemls(imls), xemls(imls)) * piu180

        endif

     enddo

  enddo

! Loop over all W points

  do iw = 2,nwa

     im1 = itab_wd(iw)%im(1)
     im2 = itab_wd(iw)%im(2)
     im3 = itab_wd(iw)%im(3)

     if     (topm(im1) <= topm(im2) .and. topm(im1) <= topm(im3)) then

        iu1 = itab_wd(iw)%iu(1)
        iu2 = itab_wd(iw)%iu(2)
        iu3 = itab_wd(iw)%iu(3)

     elseif (topm(im2) <= topm(im1) .and. topm(im2) <= topm(im3)) then

        im1 = itab_wd(iw)%im(2)
        im2 = itab_wd(iw)%im(3)
        im3 = itab_wd(iw)%im(1)

        iu1 = itab_wd(iw)%iu(2)
        iu2 = itab_wd(iw)%iu(3)
        iu3 = itab_wd(iw)%iu(1)

     else

        im1 = itab_wd(iw)%im(3)
        im2 = itab_wd(iw)%im(1)
        im3 = itab_wd(iw)%im(2)

        iu1 = itab_wd(iw)%iu(3)
        iu2 = itab_wd(iw)%iu(1)
        iu3 = itab_wd(iw)%iu(2)

     endif

! Loop over active K levels for this IW (those that intersect topography)

     do k = kmin_w(iw),kmax_w(iw)

! Get kk index for W, U1, U2, U3

        kkw = k - kmin_w(iw)

        kku1 = k - kmin_u(iu1)
        kku2 = k - kmin_u(iu2)
        kku3 = k - kmin_u(iu3)

        iwls = iwls0_w(iw) + kkw

! Initialize horizontal indices to zero

        iuls_lo = 0
        iuls_hi = 0

        iu1ls = 0
        iu2ls = 0
        iu3ls = 0

        iu1ls_mlo = 0
        iu1ls_mhi = 0

        iu2ls_mlo = 0
        iu2ls_mhi = 0

        iu3ls_mlo = 0
        iu3ls_mhi = 0

        im1k = 0
        im2k = 0
        im3k = 0

! Set horizontal indices to nonzero values if relevant

        if (k == k1topm(im1) .or. k == k2topm(im1)) im1k = im1
        if (k == k1topm(im2) .or. k == k2topm(im2)) im2k = im2
        if (k == k1topm(im3) .or. k == k2topm(im3)) im3k = im3

        if (kkw >= 1)         iuls_lo = iuls0_w(iw) + kkw - 1
        if (kkw < kadd_w(iw)) iuls_hi = iuls0_w(iw) + kkw

        if (kku1 >= 0 .and. kku1 <= kadd_u(iu1)) iu1ls = iuls0_u(iu1) + kku1
        if (kku2 >= 0 .and. kku2 <= kadd_u(iu2)) iu2ls = iuls0_u(iu2) + kku2
        if (kku3 >= 0 .and. kku3 <= kadd_u(iu3)) iu3ls = iuls0_u(iu3) + kku3

        if (iu1ls == 0 .and. im2k > 0 .and. im3k > 0) iu1ls = iuls0_u(iu1)
        if (iu2ls == 0 .and. im3k > 0 .and. im1k > 0) iu2ls = iuls0_u(iu2)
        if (iu3ls == 0 .and. im1k > 0 .and. im2k > 0) iu3ls = iuls0_u(iu3)

! Get M endpoints for IULS1 (not used if kadd_u(iu1) = 0)

        if (kku1 >= 0 .and. kku1 <= kadd_u(iu1)) then

           iu1ls_mlo = mksfc_itab_uls(iu1ls)%im(1)
           iu1ls_mhi = mksfc_itab_uls(iu1ls)%im(2)

        endif

! Get M endpoints for IULS2 (not used if kadd_u(iu2) = 0)

        if (kku2 >= 0 .and. kku2 <= kadd_u(iu2)) then

           iu2ls_mlo = mksfc_itab_uls(iu2ls)%im(1)
           iu2ls_mhi = mksfc_itab_uls(iu2ls)%im(2)

        endif

! Get M endpoints for IULS3 (not used if kadd_u(iu3) = 0)

        if (kku3 >= 0 .and. kku3 <= kadd_u(iu3)) then

           iu3ls_mlo = mksfc_itab_uls(iu3ls)%im(1)
           iu3ls_mhi = mksfc_itab_uls(iu3ls)%im(2)

        endif

        if (im1k > 0 .and. im2k > 0 .and. im3k > 0) then  ! 4

           mksfc_itab_wls(iwls)%im(1) = im1k
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = im2k
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = im3k
           mksfc_itab_wls(iwls)%iu(3) = iu2ls

           mksfc_itab_wls(iwls)%npoly = 3

        elseif (im1k > 0 .and. im2k > 0 .and. iu1ls > 0 .and. iu2ls > 0) then  ! 8

           mksfc_itab_wls(iwls)%im(1) = im1k
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = im2k
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iuls_hi
           mksfc_itab_wls(iwls)%im(4) = iu2ls_mhi
           mksfc_itab_wls(iwls)%iu(4) = iu2ls

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu2ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu1ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        elseif (im2k > 0 .and. im3k > 0 .and. iu2ls > 0 .and. iu3ls > 0) then  ! 7

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = im2k
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = im3k
           mksfc_itab_wls(iwls)%iu(3) = iu2ls
           mksfc_itab_wls(iwls)%im(4) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(4) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        elseif (im3k > 0 .and. im1k > 0 .and. iu3ls > 0 .and. iu1ls > 0) then  ! 9

           mksfc_itab_wls(iwls)%im(1) = im1k
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = iu3ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu1ls
           mksfc_itab_wls(iwls)%im(4) = im3k
           mksfc_itab_wls(iwls)%iu(4) = iu2ls

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu1ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu3ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        elseif (im1k > 0 .and. iu1ls == 0) then  ! 1

           mksfc_itab_wls(iwls)%im(1) = im1k
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = iu3ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu2ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu2ls

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu2ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu3ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 3

        elseif (im3k > 0 .and. iu3ls == 0) then  ! 2

           mksfc_itab_wls(iwls)%im(1) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iuls_lo
           mksfc_itab_wls(iwls)%im(2) = iu1ls_mlo
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = im3k
           mksfc_itab_wls(iwls)%iu(3) = iu2ls

!------------------------------------------------------
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 3

        elseif (im2k > 0 .and. iu2ls == 0) then  ! 3

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = im2k
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mlo
           mksfc_itab_wls(iwls)%iu(3) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 3

        elseif (im3k > 0 .and. iu1ls > 0 .and. iu2ls > 0 .and. iu3ls > 0) then  ! 5

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = iu3ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu1ls
           mksfc_itab_wls(iwls)%im(4) = im3k
           mksfc_itab_wls(iwls)%iu(4) = iu2ls
           mksfc_itab_wls(iwls)%im(5) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(5) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu1ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu3ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 5

        elseif (im2k > 0 .and. iu1ls > 0 .and. iu2ls > 0 .and. iu3ls > 0) then  ! 6

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = im2k
           mksfc_itab_wls(iwls)%iu(2) = iu1ls
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iuls_hi
           mksfc_itab_wls(iwls)%im(4) = iu2ls_mhi
           mksfc_itab_wls(iwls)%iu(4) = iu2ls
           mksfc_itab_wls(iwls)%im(5) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(5) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu2ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu1ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 5

        elseif (iu2ls > 0 .and. iu3ls > 0) then  ! 10

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = iu3ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu2ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu2ls
           mksfc_itab_wls(iwls)%im(4) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(4) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu2ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu3ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        elseif (iu1ls > 0 .and. iu2ls > 0) then  ! 11

           mksfc_itab_wls(iwls)%im(1) = iu1ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu1ls
           mksfc_itab_wls(iwls)%im(2) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu2ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu2ls
           mksfc_itab_wls(iwls)%im(4) = iu2ls_mlo
           mksfc_itab_wls(iwls)%iu(4) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu2ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu1ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        elseif (iu3ls > 0 .and. iu1ls > 0) then  ! 12

           mksfc_itab_wls(iwls)%im(1) = iu3ls_mlo
           mksfc_itab_wls(iwls)%iu(1) = iu3ls
           mksfc_itab_wls(iwls)%im(2) = iu3ls_mhi
           mksfc_itab_wls(iwls)%iu(2) = iuls_hi
           mksfc_itab_wls(iwls)%im(3) = iu1ls_mhi
           mksfc_itab_wls(iwls)%iu(3) = iu1ls
           mksfc_itab_wls(iwls)%im(4) = iu1ls_mlo
           mksfc_itab_wls(iwls)%iu(4) = iuls_lo

!------------------------------------------------------
           mksfc_itab_uls(iuls_hi)%im(1) = iu1ls_mhi
           mksfc_itab_uls(iuls_hi)%im(2) = iu3ls_mhi

           mksfc_itab_uls(iuls_hi)%iw(1) = iwls
           mksfc_itab_uls(iuls_lo)%iw(2) = iwls
!------------------------------------------------------

           mksfc_itab_wls(iwls)%npoly = 4

        endif

! Loop over land/sea polygon vertices/edges

        npoly = mksfc_itab_wls(iwls)%npoly

        do jmls = 1,npoly

           imls = mksfc_itab_wls(iwls)%im(jmls)
           iuls = mksfc_itab_wls(iwls)%iu(jmls)

! Evaluate x,y coordinates of current M point on polar stereographic plane
! tangent at IW

           call e_ps(xemls(imls), yemls(imls), zemls(imls),  &
                     glatw(iw), glonw(iw), xp(jmls), yp(jmls))

! Expand earth relative surface coordinates out to terrain height
! for later computing unit normals for sloping surface cells

           ef = (erad + zmls(imls)) / erad
           xemsfc(jmls) = ef * xemls(imls)
           yemsfc(jmls) = ef * yemls(imls)
           zemsfc(jmls) = ef * zemls(imls)

        enddo

! Find barycenter and unit normal of current polygon

        if (npoly == 3) then

           baryx = (xp(1) + xp(2) + xp(3)) / 3.      
           baryy = (yp(1) + yp(2) + yp(3)) / 3.      

           area = abs(.5 * (xp(1) * (yp(2) - yp(3))  &
                + xp(2) * (yp(3) - yp(1))  &
                + xp(3) * (yp(1) - yp(2))))

           call unit_normal ( xemsfc(1),   yemsfc(1),   zemsfc(1),  &
                              xemsfc(2),   yemsfc(2),   zemsfc(2),  &
                              xemsfc(3),   yemsfc(3),   zemsfc(3),  &
                              wnxls(iwls), wnyls(iwls), wnzls(iwls) )
           
        elseif (npoly == 4) then

! 2 sub-triangles

           baryx1 = (xp(1) + xp(2) + xp(3)) / 3.      
           baryy1 = (yp(1) + yp(2) + yp(3)) / 3.      

           baryx2 = (xp(1) + xp(3) + xp(4)) / 3.
           baryy2 = (yp(1) + yp(3) + yp(4)) / 3.

           area1 = abs(.5 * (xp(1) * (yp(2) - yp(3))  &
                + xp(2) * (yp(3) - yp(1))  &
                + xp(3) * (yp(1) - yp(2))))

           area2 = abs(.5 * (xp(1) * (yp(3) - yp(4))  &
                + xp(3) * (yp(4) - yp(1))  &
                + xp(4) * (yp(1) - yp(3))))

           area = area1 + area2

           baryx = (baryx1 * area1 + baryx2 * area2) / area                 
           baryy = (baryy1 * area1 + baryy2 * area2) / area         

           call unit_normal ( xemsfc(1), yemsfc(1), zemsfc(1), &
                              xemsfc(2), yemsfc(2), zemsfc(2), &
                              xemsfc(3), yemsfc(3), zemsfc(3), &
                              wnx1,      wny1,      wnz1       )

           call unit_normal ( xemsfc(1), yemsfc(1), zemsfc(1), &
                              xemsfc(3), yemsfc(3), zemsfc(3), &
                              xemsfc(4), yemsfc(4), zemsfc(4), &
                              wnx2,      wny2,      wnz2       )

           wnxls(iwls) = (wnx1 * area1 + wnx2 * area2) / area
           wnyls(iwls) = (wny1 * area1 + wny2 * area2) / area
           wnzls(iwls) = (wnz1 * area1 + wnz2 * area2) / area

        elseif (npoly == 5) then

! 3 sub-triangles

           baryx1 = (xp(1) + xp(2) + xp(3)) / 3.      
           baryy1 = (yp(1) + yp(2) + yp(3)) / 3.      

           baryx2 = (xp(1) + xp(3) + xp(4)) / 3.
           baryy2 = (yp(1) + yp(3) + yp(4)) / 3.

           baryx3 = (xp(1) + xp(4) + xp(5)) / 3.
           baryy3 = (yp(1) + yp(4) + yp(5)) / 3.

           area1 = abs(.5 * (xp(1) * (yp(2) - yp(3))  &
                + xp(2) * (yp(3) - yp(1))  &
                + xp(3) * (yp(1) - yp(2))))

           area2 = abs(.5 * (xp(1) * (yp(3) - yp(4))  &
                + xp(3) * (yp(4) - yp(1))  &
                + xp(4) * (yp(1) - yp(3))))

           area3 = abs(.5 * (xp(1) * (yp(4) - yp(5))  &
                + xp(4) * (yp(5) - yp(1))  &
                + xp(5) * (yp(1) - yp(4))))

           area = area1 + area2 + area3

           baryx = (baryx1 * area1 + baryx2 * area2 + baryx3 * area3) / area      
           baryy = (baryy1 * area1 + baryy2 * area2 + baryy3 * area3) / area                 

           call unit_normal ( xemsfc(1), yemsfc(1), zemsfc(1), &
                              xemsfc(2), yemsfc(2), zemsfc(2), &
                              xemsfc(3), yemsfc(3), zemsfc(3), &
                              wnx1,      wny1,      wnz1       )

           call unit_normal ( xemsfc(1), yemsfc(1), zemsfc(1), &
                              xemsfc(3), yemsfc(3), zemsfc(3), &
                              xemsfc(4), yemsfc(4), zemsfc(4), &
                              wnx2,      wny2,      wnz2       )

           call unit_normal ( xemsfc(1), yemsfc(1), zemsfc(1), &
                              xemsfc(4), yemsfc(4), zemsfc(4), &
                              xemsfc(5), yemsfc(5), zemsfc(5), &
                              wnx3,      wny3,      wnz3       )
           
           wnxls(iwls) = (wnx1 * area1 + wnx2 * area2 + wnx3 * area3) / area
           wnyls(iwls) = (wny1 * area1 + wny2 * area2 + wny3 * area3) / area
           wnzls(iwls) = (wnz1 * area1 + wnz2 * area2 + wnz3 * area3) / area

        endif

! Store surface cell area

        areals(iwls) = area

! Transform barycenter from ps to earth coordinates      

        call ps_e ( xewls(iwls), yewls(iwls), zewls(iwls), &
                    glatw(iw),   glonw(iw),                &
                    baryx,       baryy                     )

! Expand WLS point coordinates out to earth radius

        expansion = erad / &
             sqrt( xewls(iwls)**2 + yewls(iwls)**2 + zewls(iwls)**2  )

        xewls(iwls) = xewls(iwls) * expansion
        yewls(iwls) = yewls(iwls) * expansion
        zewls(iwls) = zewls(iwls) * expansion

        raxis = sqrt( xewls(iwls)**2 + yewls(iwls)**2 )
        glatwls(iwls) = atan2( zewls(iwls), raxis      ) * piu180
        glonwls(iwls) = atan2( yewls(iwls), xewls(iwls)) * piu180

     enddo ! k

  enddo ! iw

! Deallocate atmosphere grid arrays defined on Delaunay grid

  deallocate (xew, yew, zew, glatw, glonw)
  deallocate (arm0, glatm, glonm)

! Initialize mksfc land/sea flags

  do imls = 1,nmls
     mpt_sea (imls) = 0
     mpt_land(imls) = 0
  enddo

  do iuls = 2,nuls
     upt_sea (iuls) = 0
     upt_land(iuls) = 0
  enddo

  do iwls = 2,nwls
     wpt_sea (iwls) = 0
     wpt_land(iwls) = 0
  enddo

  if (ivegflg == 2) then

! If ivegflg == 2, fill sea/land cell areas and IW values, plus
! leaf_class for land cells, from default value defined in OLAMIN.
! User customization can be done here.

     leaf_class(2:nwls) = nvgcon

  else

! If ivegflg == 1, fill sea/land cells with leaf_class from leaf database

     call leaf_database_read(nwls,             &
          glatwls(:), &
          glonwls(:), &
          veg_database,     &
          veg_database,     &
          'leaf_class',     &
          idatp=idatp       )

     do iwls = 2,nwls
        call datp_datq( idatp(iwls), idatq)
        leaf_class(iwls) = idatq
     enddo

  endif

! Initialize global sea/land cell counters to 1

  nws = 1
  nwl = 1

! Loop over all land/sea cells

  do iwls = 2,nwls

! If area of land/sea cell is below threshold, skip cell

     if (areals(iwls) < 1.e2) cycle

! Count up individual land & sea cells and assign indices to them.
! Flag M and U points for association with land and/or sea cells.

     if (leaf_class(iwls) <= 1) then
        nws = nws + 1

        wpt_sea(iwls) = nws

        do jmls = 1,mksfc_itab_wls(iwls)%npoly
           imls = mksfc_itab_wls(iwls)%im(jmls)
           iuls = mksfc_itab_wls(iwls)%iu(jmls)

           mpt_sea(imls) = 1
           upt_sea(iuls) = 1
        enddo

     else  
        nwl = nwl + 1

        wpt_land(iwls) = nwl

        do jmls = 1,mksfc_itab_wls(iwls)%npoly
           imls = mksfc_itab_wls(iwls)%im(jmls)
           iuls = mksfc_itab_wls(iwls)%iu(jmls)

           mpt_land(imls) = 1
           upt_land(iuls) = 1
        enddo
     endif


  enddo

! Count up land & sea M and U pts and assign indices to them.

  nms = 1
  nml = 1

  do imls = 2,nmls

     if (mpt_sea(imls) == 1) then
        nms = nms + 1
        mpt_sea(imls) = nms
     endif

     if (mpt_land(imls) == 1) then
        nml = nml + 1
        mpt_land(imls) = nml
     endif
  enddo

  nus = 1
  nul = 1

  do iuls = 2,nuls

     if (upt_sea (iuls) == 1) then
        nus = nus + 1
        upt_sea (iuls) = nus
     endif

     if (upt_land(iuls) == 1) then
        nul = nul + 1
        upt_land(iuls) = nul
     endif
  enddo

! Now that final nws, nwl, nus, nul, nms, and nml values are determined, 
! allocate separate sea and land cell arrays.  Copy data from temporary 
! arrays to these.

  call alloc_sea_grid (nms, nus, nws)
  call alloc_land_grid(nml, nul, nwl, nzg)

! Copy earth coordinate values from combined land + sea arrays to individual
! land & sea arrays.

  do imls = 2,nmls

     if (mpt_sea(imls) > 1) then
        ims = mpt_sea(imls)

        sea%xem  (ims) = xemls  (imls)
        sea%yem  (ims) = yemls  (imls)
        sea%zem  (ims) = zemls  (imls)
        sea%zm   (ims) = zmls   (imls)
        sea%glatm(ims) = glatmls(imls)
        sea%glonm(ims) = glonmls(imls)

     endif

     if (mpt_land(imls) > 1) then
        iml = mpt_land(imls)

        land%xem  (iml) = xemls(imls)
        land%yem  (iml) = yemls(imls)
        land%zem  (iml) = zemls(imls)
        land%zm   (iml) = zmls (imls)
        land%glatm(iml) = glatmls(imls)
        land%glonm(iml) = glonmls(imls)

     endif

  enddo

! Copy data from combined land + sea arrays to individual land & sea arrays.

  do iwls = 2,nwls

     if (leaf_class(iwls) <= 1) then

        iws = wpt_sea(iwls)

! If sea cell has been skipped due to small size, iws = 0; skip cell

        if (iws < 2) cycle
        
        sea%leaf_class(iws) = leaf_class(iwls)
        sea%area      (iws) = areals    (iwls)
        sea%xew       (iws) = xewls     (iwls)
        sea%yew       (iws) = yewls     (iwls)
        sea%zew       (iws) = zewls     (iwls)
        sea%glatw     (iws) = glatwls   (iwls)
        sea%glonw     (iws) = glonwls   (iwls)
        
        itab_ws(iws)%npoly = mksfc_itab_wls(iwls)%npoly

        do jmls = 1, itab_ws(iws)%npoly
           imls = mksfc_itab_wls(iwls)%im(jmls)
           iuls = mksfc_itab_wls(iwls)%iu(jmls)

           itab_ws(iws)%im(jmls) = mpt_sea(imls)
           itab_ws(iws)%iu(jmls) = upt_sea(iuls)
        enddo

     else

        iwl = wpt_land(iwls)

! If land cell has been skipped due to small size, iwl = 0; skip cell

        if (iwl < 2) cycle

        land%leaf_class(iwl) = leaf_class(iwls)
        land%area      (iwl) = areals    (iwls)
        land%xew       (iwl) = xewls     (iwls)
        land%yew       (iwl) = yewls     (iwls)
        land%zew       (iwl) = zewls     (iwls)
        land%glatw     (iwl) = glatwls   (iwls)
        land%glonw     (iwl) = glonwls   (iwls)
        land%wnx       (iwl) = wnxls     (iwls)
        land%wny       (iwl) = wnyls     (iwls)
        land%wnz       (iwl) = wnzls     (iwls)

        itab_wl(iwl)%npoly = mksfc_itab_wls(iwls)%npoly

        do jmls = 1, itab_wl(iwl)%npoly
           imls = mksfc_itab_wls(iwls)%im(jmls)
           iuls = mksfc_itab_wls(iwls)%iu(jmls)

           itab_wl(iwl)%im(jmls) = mpt_land(imls)
           itab_wl(iwl)%iu(jmls) = upt_land(iuls)
        enddo

     endif

  enddo

! Copy U pt neighbor data from combined land + sea arrays to individual 
! land & sea arrays.

  do iuls = 2,nuls

     im1ls = mksfc_itab_uls(iuls)%im(1)
     im2ls = mksfc_itab_uls(iuls)%im(2)
     iw1ls = mksfc_itab_uls(iuls)%iw(1)
     iw2ls = mksfc_itab_uls(iuls)%iw(2)

     ius = upt_sea(iuls)
     iul = upt_land(iuls)

     if (ius > 1) then
        itab_us(ius)%im(1) = mpt_sea(im1ls)
        itab_us(ius)%im(2) = mpt_sea(im2ls)
        itab_us(ius)%iw(1) = max(1,wpt_sea(iw1ls))
        itab_us(ius)%iw(2) = max(1,wpt_sea(iw2ls))
     endif

     if (iul > 1) then
        itab_ul(iul)%im(1) = mpt_land(im1ls)
        itab_ul(iul)%im(2) = mpt_land(im2ls)
        itab_ul(iul)%iw(1) = max(1,wpt_land(iw1ls))
        itab_ul(iul)%iw(2) = max(1,wpt_land(iw2ls))
     endif

  enddo

! Initialize soil textural class

  if (isoilflg == 2) then

! If soilflg == 2, fill land cells with default horizontally homogeneous
! soil textural class value defined in OLAMIN.
! User customization can be done here.

     land%ntext_soil(1:nzg,2:nwl) = nslcon

  else

! If soilflg == 2, read soil textural class from database

     call leaf_database_read(nwl, &
          land%glatw,             &
          land%glonw,             &
          soil_database,          &
          soil_database,          &
          'soil_text',            &
          idatp=idatp             )

! Loop over all land cells (already defined and filled with leaf_class)

     do iwl = 2,nwl
        call datp_datsoil(idatp(iwl), datsoil)

! For now, assign single-level FAO textural class to all soil layers.

        do k = 1, nzg
           land%ntext_soil(k,iwl) = datsoil
        enddo

     enddo

  endif

! Write land file (contains LEAF CLASS, SOIL TEXTURAL CLASS, and grid info)   

  write(io6,'(/,a)') 'calling landfile_write'

  call landfile_write()

! Write sea file (contains LEAF CLASS and grid info)   

  write(io6,'(/,a)') 'calling seafile_write'

  call seafile_write()

  write(io6,'(/,a)') 'called seafile_write'

! Deallocate mksfc arrays

  deallocate (mksfc_itab_uls)
  deallocate (mksfc_itab_wls)

  deallocate (mpt_sea, mpt_land)
  deallocate (upt_sea, upt_land)
  deallocate (wpt_sea, wpt_land)
  deallocate (idatp)

  deallocate (leaf_class, areals)
  deallocate (xewls, yewls, zewls)
  deallocate (glatwls, glonwls)
  deallocate (wnxls, wnyls, wnzls)

  deallocate (xemls, yemls, zemls, zmls, glatmls, glonmls)

  return
end subroutine makesfc

!==========================================================================

subroutine datp_datq(idatp,idatq)

! This subroutine maps the input idatp classes to a smaller set idatq
! which represents the full set of LEAF-2 or LEAF-3 classes for which 
! LSP values are
! defined.

implicit none

integer, intent(in)  :: idatp
integer, intent(out) :: idatq

integer :: catb(0:100)

!  Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes
!  (see leaf3_document).

!-------------------------------------------!
data catb/ 0,                             & !
          19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  0
          10, 2,17, 1, 0,12,13,14,18, 4,  & ! 10
           4, 4,14,14, 6, 6, 4, 7, 7,15,  & ! 20
          15, 6, 7, 7,15,16,16,16,16, 8,  & ! 30
           8, 8,18,17,17,12,12, 7,10, 3,  & ! 40
          10,10,11,14,18,18,18,18,13, 6,  & ! 50
           5, 4,11,12, 0, 0, 0, 0, 3, 2,  & ! 60
           3,20, 0,17,17,17, 4,14, 7, 3,  & ! 70
           3, 3, 3, 3, 3, 3, 8,12, 7, 6,  & ! 80
          18,15,15,15, 4, 5, 0, 0, 0, 0   / ! 90  ! 97 & 98 not used
!-------------------------------------------!     ! 99 is Goode Homolosine 
                                                  !    empty space
!          1  2  3  4  5  6  7  8  9 10           ! 100 is missing data
                                                  ! Map all of these to ocean (idatq=0)
idatq = catb(idatp)

return
end subroutine datp_datq

!==========================================================================

subroutine datp_datsoil(idatp,idatsoil)

! This subroutine maps the input idatp soil classes to a smaller set idatsoil
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

implicit none

integer, intent(in) :: idatp
integer, intent(out) :: idatsoil

integer :: catb(0:133)

! (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our 
! own missing value designated 133, to the USDA soil textural classes.  FAO 
! classes [0] (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not 
! used in original FAO dataset), [132] (water), and [133] (our missing value) 
! are all mapped to a default class of sandy clay loam in case they happen to
! correspond to a land surface area in the landuse dataset that RAMS uses to
! define land area.  We wrote missing value class 133 to the RAMS FAO files
! whenever a negative value (which is out of range of defined values) was
! found in the original FAO dataset, which occurred in about 2.6% of the
! pixels.  For the remaining FAO classes, a cross reference table to Zobler 
! soil texture classes that was provided, plus our own cross referencing table
! from Zobler to USDA classes listed below, provides the mapping from FAO to 
! USDA.  In this mapping, we use only organic USDA classes and omit nonorganic
! classes since most soils do contain organic matter, and organic content 
! information is not provided in the Zobler classes.

!  Zobler Class              USDA Class

!  1  Coarse                 2  Loamy sand
!  2  Medium                 4  Silt loam
!  3  Fine                   8  Clay loam
!  4  Coarse-medium          3  Sandy loam
!  5  Coarse-fine            6  Sandy clay loam
!  6  Medium-fine            7  Silty clay loam
!  7  Coarse-medium-fine     6  Sandy clay loam
!  8  Organic matter         5  Loam

!                            1  Sand (not used)
!                            9  Sandy clay (not used)
!                           10  Silty clay (not used)
!                           11  Clay (not used)
!                           12  Peat (not used)

!-------------------------------------------!
data catb/ 6,                             & !
           6, 4, 4, 7, 7, 8, 6, 4, 4, 4,  & !   0
           7, 4, 4, 4, 8, 4, 8, 4, 4, 8,  & !  10
           4, 2, 4, 4, 4, 4, 6, 8, 8, 8,  & !  20
           4, 8, 8, 2, 6, 4, 7, 4, 4, 3,  & !  30
           4, 6, 7, 4, 4, 4, 4, 4, 4, 4,  & !  40
           4, 4, 4, 4, 4, 4, 2, 4, 4, 2,  & !  50
           4, 3, 4, 2, 7, 6, 4, 4, 6, 8,  & !  60
           8, 7, 2, 5, 4, 5, 6, 6, 4, 2,  & !  70
           2, 2, 4, 6, 2, 2, 2, 2, 2, 4,  & !  80
           2, 2, 2, 4, 2, 4, 3, 6, 2, 7,  & !  90
           4, 4, 4, 8, 8, 8, 3, 7, 4, 4,  & ! 100
           4, 3, 6, 4, 2, 4, 4, 4, 2, 2,  & ! 110
           2, 4, 6, 4, 4, 7, 7, 6, 3, 2,  & ! 120
           2, 6, 6                        / ! 130
!-------------------------------------------!
!          1  2  3  4  5  6  7  8  9 10

idatsoil = catb(idatp)

return
end subroutine datp_datsoil

!==========================================================================

subroutine landfile_write()

use max_dims,   only: maxnlspoly, pathlen
use leaf_coms,  only: nzg, nml, nul, nwl, landusefile, slz
use mem_leaf,   only: land, itab_ul, itab_wl
use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
use misc_coms,  only: io6, ngrids, nxp

implicit none

integer :: iclobber1 = 1
integer :: ndims, idims(2)
integer :: iwl
integer :: iscr(maxnlspoly,nwl)

character(pathlen) :: flnm

! Write land file header information

flnm=trim(landusefile)//'-S'//'.h5'

call shdf5_open(flnm,'W',iclobber1)

ndims = 1
idims(1) = 1

call shdf5_orec(ndims, idims, 'nml'   , ivars=nml)
call shdf5_orec(ndims, idims, 'nul'   , ivars=nul)
call shdf5_orec(ndims, idims, 'nwl'   , ivars=nwl)
call shdf5_orec(ndims, idims, 'nzg'   , ivars=nzg)
call shdf5_orec(ndims, idims, 'ngrids', ivars=ngrids)
call shdf5_orec(ndims, idims, 'nxp'   , ivars=nxp)

! Write arrays to land file

ndims = 1
idims(1) = nzg

call shdf5_orec(ndims, idims, 'slz', rvara=slz)

idims(1) = nwl

call shdf5_orec(ndims, idims, 'leaf_class', ivara=land%leaf_class)
call shdf5_orec(ndims, idims, 'land_area' , rvara=land%area)
call shdf5_orec(ndims, idims, 'xewl'      , rvara=land%xew)
call shdf5_orec(ndims, idims, 'yewl'      , rvara=land%yew)
call shdf5_orec(ndims, idims, 'zewl'      , rvara=land%zew)
call shdf5_orec(ndims, idims, 'glatwl'    , rvara=land%glatw)
call shdf5_orec(ndims, idims, 'glonwl'    , rvara=land%glonw)
call shdf5_orec(ndims, idims, 'wnxl'      , rvara=land%wnx)
call shdf5_orec(ndims, idims, 'wnyl'      , rvara=land%wny)
call shdf5_orec(ndims, idims, 'wnzl'      , rvara=land%wnz)
call shdf5_orec(ndims, idims, 'nlpoly'    , ivara=itab_wl(:)%npoly)

idims(1) = nul

call shdf5_orec(ndims, idims, 'itab_ul%im1', ivara=itab_ul(:)%im(1))
call shdf5_orec(ndims, idims, 'itab_ul%im2', ivara=itab_ul(:)%im(2))
call shdf5_orec(ndims, idims, 'itab_ul%iw1', ivara=itab_ul(:)%iw(1))
call shdf5_orec(ndims, idims, 'itab_ul%iw2', ivara=itab_ul(:)%iw(2))

idims(1) = nml

call shdf5_orec(ndims, idims, 'xeml'  , rvara=land%xem)
call shdf5_orec(ndims, idims, 'yeml'  , rvara=land%yem)
call shdf5_orec(ndims, idims, 'zeml'  , rvara=land%zem)
call shdf5_orec(ndims, idims, 'zml'   , rvara=land%zm)
call shdf5_orec(ndims, idims, 'glatml', rvara=land%glatm)
call shdf5_orec(ndims, idims, 'glonml', rvara=land%glonm)

ndims = 2
idims(1) = nzg
idims(2) = nwl

call shdf5_orec(ndims, idims, 'ntext_soil', ivara=land%ntext_soil)

ndims = 2
idims(1) = maxnlspoly
idims(2) = nwl

! Copy itab_wl%im to scratch array for output

do iwl = 1,nwl
   iscr(1:maxnlspoly,iwl) = itab_wl(iwl)%im(1:maxnlspoly)
enddo

call shdf5_orec(ndims,idims,'itab_wl%im',ivara=iscr)

! Copy itab_wl%iu to scratch array for output

do iwl = 1,nwl
   iscr(1:maxnlspoly,iwl) = itab_wl(iwl)%iu(1:maxnlspoly)
enddo

call shdf5_orec(ndims,idims,'itab_wl%iu',ivara=iscr)

call shdf5_close()

return
end subroutine landfile_write

!==========================================================================

subroutine seafile_write()

use max_dims,   only: maxnlspoly, pathlen
use sea_coms,   only: nms, nus, nws, seafile
use mem_sea,    only: sea, itab_us, itab_ws
use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
use misc_coms,  only: io6, ngrids, nxp

implicit none

integer :: iclobber1 = 1
integer :: ndims, idims(2)
integer :: iws
integer :: iscr(maxnlspoly,nws)

character(pathlen) :: flnm

! Write sea file header information

flnm=trim(seafile)//'-S'//'.h5'

call shdf5_open(flnm,'W',iclobber1)

ndims = 1
idims(1) = 1

call shdf5_orec(ndims, idims, 'nms'   , ivars=nms)
call shdf5_orec(ndims, idims, 'nus'   , ivars=nus)
call shdf5_orec(ndims, idims, 'nws'   , ivars=nws)
call shdf5_orec(ndims, idims, 'ngrids', ivars=ngrids)
call shdf5_orec(ndims, idims, 'nxp'   , ivars=nxp)

! Write arrays to sea file

ndims = 1
idims(1) = nws

call shdf5_orec(ndims, idims, 'leaf_class', ivara=sea%leaf_class)
call shdf5_orec(ndims, idims, 'sea_area'  , rvara=sea%area)
call shdf5_orec(ndims, idims, 'xews'      , rvara=sea%xew)
call shdf5_orec(ndims, idims, 'yews'      , rvara=sea%yew)
call shdf5_orec(ndims, idims, 'zews'      , rvara=sea%zew)
call shdf5_orec(ndims, idims, 'glatws'    , rvara=sea%glatw)
call shdf5_orec(ndims, idims, 'glonws'    , rvara=sea%glonw)
call shdf5_orec(ndims, idims, 'nspoly'    , ivara=itab_ws(:)%npoly)

idims(1) = nus

call shdf5_orec(ndims, idims, 'itab_us%im1', ivara=itab_us(:)%im(1))
call shdf5_orec(ndims, idims, 'itab_us%im2', ivara=itab_us(:)%im(2))
call shdf5_orec(ndims, idims, 'itab_us%iw1', ivara=itab_us(:)%iw(1))
call shdf5_orec(ndims, idims, 'itab_us%iw2', ivara=itab_us(:)%iw(2))

idims(1) = nms

call shdf5_orec(ndims, idims, 'xems'  , rvara=sea%xem)
call shdf5_orec(ndims, idims, 'yems'  , rvara=sea%yem)
call shdf5_orec(ndims, idims, 'zems'  , rvara=sea%zem)
call shdf5_orec(ndims, idims, 'zms'   , rvara=sea%zm)
call shdf5_orec(ndims, idims, 'glatms', rvara=sea%glatm)
call shdf5_orec(ndims, idims, 'glonms', rvara=sea%glonm)

ndims = 2
idims(1) = maxnlspoly
idims(2) = nws

! Copy itab_ws%im to temporary array for output

do iws = 1,nws
   iscr(1:maxnlspoly,iws) = itab_ws(iws)%im(1:maxnlspoly)
enddo

call shdf5_orec(ndims, idims, 'itab_ws%im', ivara=iscr)

! Copy itab_ws%iu to scratch array for output

do iws = 1,nws
   iscr(1:maxnlspoly,iws) = itab_ws(iws)%iu(1:maxnlspoly)
enddo

call shdf5_orec(ndims,idims,'itab_ws%iu',ivara=iscr)

call shdf5_close()

return
end subroutine seafile_write

!==========================================================================

subroutine landfile_read()

use max_dims,   only: maxnlspoly, pathlen
use leaf_coms,  only: nzg, nml, nul, nwl, mml, mul, mwl, landusefile, slz
use mem_leaf,   only: land, itab_ul, itab_wl, alloc_land_grid
use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
use misc_coms,  only: io6

implicit none

integer              :: ndims, idims(2)
character(pathlen)   :: flnm
logical              :: there
integer              :: iwl
integer, allocatable :: iscr(:,:)

!-------------------------------------------------------------------------------
! STEP 1: Open land file and read 4 array dimensions
!-------------------------------------------------------------------------------

flnm = trim(landusefile)//'-S'//'.h5'

write(io6,*) 'Checking leaf file ',trim(flnm)

inquire(file=flnm, exist=there)

if (.not. there) then
   write(io6,*) 'Land file was not found - stopping run'
   stop 'stop: no land file'
endif

call shdf5_open(flnm,'R')

ndims = 1
idims(1) = 1

call shdf5_irec(ndims, idims, 'nml', ivars=nml)
call shdf5_irec(ndims, idims, 'nul', ivars=nul)
call shdf5_irec(ndims, idims, 'nwl', ivars=nwl)
call shdf5_irec(ndims, idims, 'nzg', ivars=nzg)

write(io6, '(/,a)')   '==============================================='
write(io6, '(a)')     'Reading from land file:'
write(io6, '(a,4i8)') '  nml, nul, nwl, nzg = ', nml, nul, nwl, nzg
write(io6, '(a,/)')   '==============================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate land grid arrays
!-------------------------------------------------------------------------------

call alloc_land_grid(nml,nul,nwl,nzg)
allocate (iscr(maxnlspoly,nwl))

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from land file
!-------------------------------------------------------------------------------

ndims = 1
idims(1) = nzg

call shdf5_irec(ndims, idims, 'slz', rvara=slz)

! Write arrays to land file

idims(1) = nwl

call shdf5_irec(ndims, idims, 'leaf_class', ivara=land%leaf_class)
call shdf5_irec(ndims, idims, 'land_area' , rvara=land%area)
call shdf5_irec(ndims, idims, 'xewl'      , rvara=land%xew)
call shdf5_irec(ndims, idims, 'yewl'      , rvara=land%yew)
call shdf5_irec(ndims, idims, 'zewl'      , rvara=land%zew)
call shdf5_irec(ndims, idims, 'glatwl'    , rvara=land%glatw)
call shdf5_irec(ndims, idims, 'glonwl'    , rvara=land%glonw)
call shdf5_irec(ndims, idims, 'wnxl'      , rvara=land%wnx)
call shdf5_irec(ndims, idims, 'wnyl'      , rvara=land%wny)
call shdf5_irec(ndims, idims, 'wnzl'      , rvara=land%wnz)
call shdf5_irec(ndims ,idims, 'nlpoly'    , ivara=itab_wl(:)%npoly)

idims(1) = nul

call shdf5_irec(ndims, idims, 'itab_ul%im1', ivara=itab_ul(:)%im(1))
call shdf5_irec(ndims, idims, 'itab_ul%im2', ivara=itab_ul(:)%im(2))
call shdf5_irec(ndims, idims, 'itab_ul%iw1', ivara=itab_ul(:)%iw(1))
call shdf5_irec(ndims, idims, 'itab_ul%iw2', ivara=itab_ul(:)%iw(2))

idims(1) = nml

call shdf5_irec(ndims, idims, 'xeml',   rvara=land%xem)
call shdf5_irec(ndims, idims, 'yeml',   rvara=land%yem)
call shdf5_irec(ndims, idims, 'zeml',   rvara=land%zem)
call shdf5_irec(ndims, idims, 'zml' ,   rvara=land%zm)
call shdf5_irec(ndims, idims, 'glatml', rvara=land%glatm)
call shdf5_irec(ndims, idims, 'glonml', rvara=land%glonm)

ndims = 2
idims(1) = nzg
idims(2) = nwl

call shdf5_irec(ndims, idims, 'ntext_soil', ivara=land%ntext_soil)

ndims = 2
idims(1) = maxnlspoly
idims(2) = nwl

call shdf5_irec(ndims,idims,'itab_wl%im',ivara=iscr)

! Copy input scratch array to itab_wl%im

do iwl = 1,nwl
   itab_wl(iwl)%im(1:maxnlspoly) = iscr(1:maxnlspoly,iwl)
enddo

call shdf5_irec(ndims,idims,'itab_wl%iu',ivara=iscr)

! Copy input scratch array to itab_wl%iu

do iwl = 1,nwl
   itab_wl(iwl)%iu(1:maxnlspoly) = iscr(1:maxnlspoly,iwl)
enddo

call shdf5_close()
deallocate(iscr)

mml = nml
mul = nul
mwl = nwl

return
end subroutine landfile_read

!==========================================================================

subroutine seafile_read()

use max_dims,   only: maxnlspoly, pathlen
use sea_coms,   only: nms, nus, nws, mms, mus, mws, seafile
use mem_sea,    only: sea, itab_us, itab_ws, alloc_sea_grid
use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
use misc_coms,  only: io6

implicit none

integer              :: ndims, idims(2)
character(pathlen)   :: flnm
logical              :: there
integer              :: iws
integer, allocatable :: iscr(:,:)

!-------------------------------------------------------------------------------
! STEP 1: Open SEAFILE and read 3 array dimensions
!-------------------------------------------------------------------------------

flnm = trim(seafile)//'-S'//'.h5'

write(io6,*) 'Checking sea file ', trim(flnm)

inquire(file=flnm, exist=there)

if (.not. there) then
   write(io6,*) 'SEA file was not found - stopping run'
   stop 'stop: no sea file'
endif

call shdf5_open(flnm,'R')

ndims = 1
idims(1) = 1

call shdf5_irec(ndims, idims, 'nms', ivars=nms)
call shdf5_irec(ndims, idims, 'nus', ivars=nus)
call shdf5_irec(ndims, idims, 'nws', ivars=nws)

write(io6, '(/,a)')   '====================================='
write(io6, '(a)')     'Reading from sea file:'
write(io6, '(a,3i8)') '  nms, nus, nws = ', nms, nus, nws
write(io6, '(a,/)')   '====================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate SEA grid arrays
!-------------------------------------------------------------------------------

call alloc_sea_grid(nms,nus,nws)
allocate (iscr(maxnlspoly,nws))

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from sea file
!-------------------------------------------------------------------------------

ndims = 1
idims(1) = nws

call shdf5_irec(ndims, idims, 'leaf_class', ivara=sea%leaf_class)
call shdf5_irec(ndims, idims, 'sea_area'  , rvara=sea%area)
call shdf5_irec(ndims, idims, 'xews'      , rvara=sea%xew)
call shdf5_irec(ndims, idims, 'yews'      , rvara=sea%yew)
call shdf5_irec(ndims, idims, 'zews'      , rvara=sea%zew)
call shdf5_irec(ndims, idims, 'glatws'    , rvara=sea%glatw)
call shdf5_irec(ndims, idims, 'glonws'    , rvara=sea%glonw)
call shdf5_irec(ndims, idims, 'nspoly'    , ivara=itab_ws(:)%npoly)

idims(1) = nus

call shdf5_irec(ndims, idims, 'itab_us%im1', ivara=itab_us(:)%im(1))
call shdf5_irec(ndims, idims, 'itab_us%im2', ivara=itab_us(:)%im(2))
call shdf5_irec(ndims, idims, 'itab_us%iw1', ivara=itab_us(:)%iw(1))
call shdf5_irec(ndims, idims, 'itab_us%iw2', ivara=itab_us(:)%iw(2))

idims(1) = nms

call shdf5_irec(ndims, idims, 'xems'  , rvara=sea%xem)
call shdf5_irec(ndims, idims, 'yems'  , rvara=sea%yem)
call shdf5_irec(ndims, idims, 'zems'  , rvara=sea%zem)
call shdf5_irec(ndims, idims, 'zms'   , rvara=sea%zm)
call shdf5_irec(ndims, idims, 'glatms', rvara=sea%glatm)
call shdf5_irec(ndims, idims, 'glonms', rvara=sea%glonm)

ndims = 2
idims(1) = maxnlspoly
idims(2) = nws

iscr(:,:) = 0

call shdf5_irec(ndims, idims, 'itab_ws%im', ivara=iscr)

! Copy temporary scratch array to itab_ws%im

do iws = 1,nws
   itab_ws(iws)%im(1:maxnlspoly) = iscr(1:maxnlspoly,iws)
enddo

call shdf5_irec(ndims, idims, 'itab_ws%iu', ivara=iscr)

! Copy temporary scratch array to itab_ws%iu

do iws = 1,nws
   itab_ws(iws)%iu(1:maxnlspoly) = iscr(1:maxnlspoly,iws)
enddo

call shdf5_close()
deallocate(iscr)

mms = nms
mus = nus
mws = nws

return
end subroutine seafile_read
