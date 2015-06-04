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
subroutine makesfc2()

! This subroutine generates LAND and SEA files for OLAM runs. 

  use mem_grid,    only: nza, nma, nwa, zm, &
                         xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                         arw0, glatw, glonw, glatm, glonm, topm, topw

  use mem_ijtabs,  only: itab_w, jtab_w, jtw_grid, jtw_lbcp

  use mem_grid,    only: xem, yem, zem, xev, xew, yew, zew, dnu, dnv, arw0

  use misc_coms,   only: io6, itopoflg, rinit, topo_database

  use sea_coms,    only: nws, seafile, iseagrid

  use leaf_coms,   only: nwl, nzg, nzs, nslcon, nvgcon, &
                         isfcl, ivegflg, landusefile, isoilflg, &
                         soil_database, veg_database

  use consts_coms, only: erad, piu180

  use land_db,     only: land_database_read

  use max_dims,    only: maxnlspoly

  use mem_mksfc,   only: itab_wls_vars, landsea_grid_vars

  use mem_leaf,    only: itab_wl, land, alloc_land_grid

  use mem_sea,     only: itab_ws, sea, alloc_sea_grid

  implicit none

  integer :: k, im, iw, im1, im2, iq, j, jm, jv, jm1, jm2, iv, ipat
  integer :: ndims, idims(1)
  integer :: npoly
  integer :: nwls, nwls_old, nwls_est, nwls_inc, ifsize
  integer :: iwls, imls, jmls
  integer :: iws, ims
  integer :: iwl, iml
  integer :: leafclass, datsoil

  real :: expansion, raxis
  real :: fracdist, fracz, ef, efw, efm1, efm2
  real :: pnx, pny, pnz, dvm1, dvm2, topv
  real :: wnx1, wny1, wnz1
  real :: wnx2, wny2, wnz2
  real :: wnx3, wny3, wnz3
  real :: xp(maxnlspoly), yp(maxnlspoly)
  real :: area, xc, yc

  integer :: nqa, jp1, jp2
  real :: topwmin, topwmax

  real :: z1,z2,z3,x1,x2,x3,y1,y2,y3
  real :: xm1, ym1, xm2, ym2, xem1, yem1, zem1
  real :: topc, xec, yec, zec
  real :: xmpat(5,nza), ympat(5,nza)
  integer :: kwpat(nza),lpoly(nza),npat

  real, allocatable :: qlat(:),qlon(:),topq(:)
  real, allocatable :: xeq(:),yeq(:),zeq(:)

  type(itab_wls_vars), allocatable :: itab_wls(:), itab_wls_temp(:)
  type(landsea_grid_vars), allocatable :: landsea_grid(:), landsea_grid_temp(:)

  integer :: ii,iii

! Estimate required array size and increment, and allocate arrays

  nwls_est = 50 * nwa
  nwls_inc = 10 * nwa

! Allocate mksfc arrays for combined land/sea grid

  allocate (itab_wls  (nwls_est))
  allocate (landsea_grid(nwls_est))

! Initialize topography at M and W points

  nqa = nma + nwa - 1

  allocate (qlat(nqa), qlon(nqa), topq(nqa))
  allocate ( xeq(nqa),  yeq(nqa),  zeq(nqa))

  qlat(1:nma) = glatm(1:nma)
  qlon(1:nma) = glonm(1:nma)

  xeq(1:nma) = xem(1:nma)
  yeq(1:nma) = yem(1:nma)
  zeq(1:nma) = zem(1:nma)

  qlat(nma+1:nqa) = glatw(2:nma)
  qlon(nma+1:nqa) = glonw(2:nma)

  xeq(nma+1:nqa) = xew(2:nma)
  yeq(nma+1:nqa) = yew(2:nma)
  zeq(nma+1:nqa) = zew(2:nma)

! Interpolate TOPO from database or initialize from customizable topo_init

  if (itopoflg == 1) then
     call land_database_read(nqa, qlat, qlon, &
        topo_database, topo_database, 'topo', datq=topq)
  else
     call topo_init(nqa,topq,qlat,qlon,xeq,yeq,zeq)
  endif

! Prevent TOPQ lower than lowest model level zm(1)

  topq(2:nqa) = max(topq(2:nqa),zm(1))

! Loop over all TOPQ points

  !$omp parallel do private(k,fracz)
  do iq = 2,nqa

! Loop over KM levels

     do k = 1,nza-1

! Check whether topq(iq) is at or above current KM level

        if (zm(k) <= topq(iq) .and. zm(k+1) >= topq(iq)) then

! Check topq(iq) within height range of KT level

           fracz = (topq(iq) - zm(k)) / (zm(k+1) - zm(k))

! If topq(iq) is in lowest or highest 2% of KT level, move it to the limit.      

           if (fracz < .02) then
              topq(iq) = zm(k)
           elseif (fracz > .98) then
              topq(iq) = zm(k+1)
           endif

           exit

        endif

     enddo

  enddo
  !$omp end parallel do

! Fill individual topm and topw values from topq array

  topm(2:nma) = topq(2:nma)
  topw(2:nwa) = topq(nma+1:nqa)

  deallocate (qlat, qlon, topq)
  deallocate (xeq, yeq, zeq)

! Initialize land/sea surface index

  nwls = 1

! Loop over Hex W points

  do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
     npoly = itab_w(iw)%npoly

! Constrain TOPW to lie within the max/min range of adjacent TOPM values

     topwmin = 1.e9
     topwmax = -1.e9

     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        if (topwmin > topm(im)) topwmin = topm(im)
        if (topwmax < topm(im)) topwmax = topm(im)
     enddo

     topw(iw) = max(topwmin,min(topwmax,topw(iw)))

     nwls_old = nwls

! Loop over all polygon edges

     do jv = 1,npoly
        jm1 = jv - 1
        if (jm1 == 0) jm1 = npoly
        jm2 = jv

        iv  = itab_w(iw)%iv(jv)
        im1 = itab_w(iw)%im(jm1)
        im2 = itab_w(iw)%im(jm2)

! Each triangular sector of a hexagonal (or pentagonal or heptagonal) grid
! cell, formed by the IW point and two adjacent IM points, has topographic
! height specified on each of its vertices, so the triangle is treated as a
! planar topographic surface.  Find lines on this surface where it is
! intersected by each zm level of the atmospheric grid, and construct polygonal
! areas bounded by the lines of intersection and triangular sector edges (using
! the algorithm in cont3f).  Each of these polygons defines a separate cell of
! the land/sea grid, and its area is later summed (in subroutine ctrlvols) with
! other cells to compute arw(k,iw) and volt(k,iw) for the atmosphere cell.

! Expand earth relative surface coordinates out to terrain height
! for computing unit normals for sloping surface cells

        efw  = (erad + topw(iw )) / erad
        efm1 = (erad + topm(im1)) / erad
        efm2 = (erad + topm(im2)) / erad

        call unit_normal ( xew(iw )*efw , yew(iw )*efw , zew(iw )*efw , &
                           xem(im1)*efm1, yem(im1)*efm1, zem(im1)*efm1, &
                           xem(im2)*efm2, yem(im2)*efm2, zem(im2)*efm2, &
                           pnx, pny, pnz )

! For cont3f algorithm, use local coordinate system with IW point at origin
! and IV edge being a line of constant positive y.

        dvm1 = sqrt((xev(iv) - xem(im1))**2 &
             +      (yev(iv) - yem(im1))**2 &
             +      (zev(iv) - zem(im1))**2)

        dvm2 = sqrt((xev(iv) - xem(im2))**2 &
             +      (yev(iv) - yem(im2))**2 &
             +      (zev(iv) - zem(im2))**2)

        topv = (dvm1 * topm(im2) + dvm2 * topm(im1)) / dnu(iv)

        x1  = 0.
        x2  = dvm1
        x3  = -dvm2
        y1  = 0.
        y2  = 0.5 * dnv(iv)
        y3  = 0.5 * dnv(iv)
        z1  = topw(iw)
        z2  = topm(im1)
        z3  = topm(im2)

        if     (z1 >= z2 .and. z2 >= z3) then
           call cont3sfc(z1,z2,z3,x1,x2,x3,y1,y2,y3,xmpat,ympat,kwpat,lpoly,npat,iw)
        elseif (z1 >= z3 .and. z3 >= z2) then
           call cont3sfc(z1,z3,z2,x1,x3,x2,y1,y3,y2,xmpat,ympat,kwpat,lpoly,npat,iw)
           call reverse_polygon(xmpat,ympat,lpoly,npat)
        elseif (z2 >= z1 .and. z1 >= z3) then
           call cont3sfc(z2,z1,z3,x2,x1,x3,y2,y1,y3,xmpat,ympat,kwpat,lpoly,npat,iw)
           call reverse_polygon(xmpat,ympat,lpoly,npat)
        elseif (z2 >= z3 .and. z3 >= z1) then
           call cont3sfc(z2,z3,z1,x2,x3,x1,y2,y3,y1,xmpat,ympat,kwpat,lpoly,npat,iw)
        elseif (z3 >= z1 .and. z1 >= z2) then
           call cont3sfc(z3,z1,z2,x3,x1,x2,y3,y1,y2,xmpat,ympat,kwpat,lpoly,npat,iw)
        elseif (z3 >= z2 .and. z2 >= z1) then
           call cont3sfc(z3,z2,z1,x3,x2,x1,y3,y2,y1,xmpat,ympat,kwpat,lpoly,npat,iw)
           call reverse_polygon(xmpat,ympat,lpoly,npat)
        endif

! Catalogue new patches that have been found in this sector, 
! and find their areas

! Loop over new patches that have been found in this sector

        do ipat = 1,npat

! Increment global land/sea surface index

           nwls = nwls + 1

! Allocate more space if necessary

           ifsize = size(itab_wls)
           if (nwls > ifsize) then
              allocate (itab_wls_temp(ifsize+nwls_inc))
              itab_wls_temp(1:ifsize) = itab_wls(1:ifsize)
              call move_alloc(itab_wls_temp, itab_wls)

              allocate (landsea_grid_temp(ifsize+nwls_inc))
              landsea_grid_temp(1:ifsize) = landsea_grid(1:ifsize)
              call move_alloc(landsea_grid_temp, landsea_grid)
           endif

! Compute polygon area and centroid for current ipat

           area = 0.
           xc = 0.
           yc = 0.

           do jp1 = 1,lpoly(ipat)
              jp2 = jp1 + 1
              if (jp2 > lpoly(ipat)) jp2 = 1

              xm1 = xmpat(jp1,ipat)
              ym1 = ympat(jp1,ipat)

              xm2 = xmpat(jp2,ipat)
              ym2 = ympat(jp2,ipat)

              area = area + 0.5 * (xm1 * ym2 - xm2 * ym1)

              xc = xc + (xm1 + xm2) * (xm1 * ym2 - xm2 * ym1)
              yc = yc + (ym1 + ym2) * (xm1 * ym2 - xm2 * ym1)

! Transform M points to earth coordinates

              xem1 = xew(iw) + xm1 * (xem(im1) - xem(im2)) / dnu(iv)  &
                             + ym1 * (xev(iv ) - xew(iw )) / (0.5 * dnv(iv))
              yem1 = yew(iw) + xm1 * (yem(im1) - yem(im2)) / dnu(iv)  &
                             + ym1 * (yev(iv ) - yew(iw )) / (0.5 * dnv(iv))
              zem1 = zew(iw) + xm1 * (zem(im1) - zem(im2)) / dnu(iv)  &
                             + ym1 * (zev(iv ) - zew(iw )) / (0.5 * dnv(iv))

              expansion = erad / sqrt(xem1**2 + yem1**2 + zem1**2)

              itab_wls(nwls)%xem(jp1) = xem1 * expansion
              itab_wls(nwls)%yem(jp1) = yem1 * expansion
              itab_wls(nwls)%zem(jp1) = zem1 * expansion

           enddo

           xc = xc / (6. * area)
           yc = yc / (6. * area)

! Interpolate earth coordinates and topo height to centroid

           xec = xew(iw) + xc * (xem(im1) - xem(im2)) / dnu(iv) &
                         + yc * (xev(iv)  - xew(iw )) / (0.5 * dnv(iv))
           yec = yew(iw) + xc * (yem(im1) - yem(im2)) / dnu(iv) &
                         + yc * (yev(iv)  - yew(iw )) / (0.5 * dnv(iv))
           zec = zew(iw) + xc * (zem(im1) - zem(im2)) / dnu(iv) &
                         + yc * (zev(iv)  - zew(iw )) / (0.5 * dnv(iv))

           topc = topw(iw) + xc * (topm(im1) - topm(im2)) / dnu(iv) &
                           + yc * (topv - topw(iw )) / (0.5 * dnv(iv))

           expansion = erad / sqrt(xec**2 + yec**2 + zec**2)

           xec = xec * expansion
           yec = yec * expansion
           zec = zec * expansion

! Compute latitude and longitude of centroid; catalogue sfc cell variables

           raxis = sqrt(xec**2 + yec**2)

           landsea_grid(nwls)%area = area

           landsea_grid(nwls)%glatw = atan2(zec,raxis) * piu180
           landsea_grid(nwls)%glonw = atan2(yec,xec) * piu180

           landsea_grid(nwls)%xew = xec
           landsea_grid(nwls)%yew = yec
           landsea_grid(nwls)%zew = zec

           landsea_grid(nwls)%topw = topc

           landsea_grid(nwls)%wnx = pnx
           landsea_grid(nwls)%wny = pny
           landsea_grid(nwls)%wnz = pnz

           itab_wls(nwls)%iw    = iw
           itab_wls(nwls)%kw    = kwpat(ipat)
           itab_wls(nwls)%npoly = lpoly(ipat)

           itab_wls(nwls)%arf_iw = area / arw0(iw)

        enddo  ! ipat

     enddo   ! jv

  enddo   ! iw

! Initialize leaf_class


  if (ivegflg == 2) then

! If ivegflg == 2, fill sea/land cell areas and IW values, plus
! leaf_class for land cells, from default value defined in OLAMIN.
! User customization can be done here.

     landsea_grid(2:nwls)%leaf_class = nvgcon

  else

! If ivegflg == 1, fill sea/land cells with leaf_class from leaf database

     call land_database_read(nwls, &
          landsea_grid(:)%glatw, &
          landsea_grid(:)%glonw, &
          veg_database, &
          veg_database, &
          'leaf_class', &
          idatq=landsea_grid(:)%idatq)




     do iwls = 2,nwls
        call oge_leafclass(landsea_grid(iwls)%idatq, leafclass)
        landsea_grid(iwls)%leaf_class = leafclass

     enddo

  endif

! Initialize global sea/land cell counters to 1

  nws = 1
  nwl = 1

! Loop over all land/sea cells

  do iwls = 2,nwls

! If area of land/sea cell is below threshold, skip cell

     if (landsea_grid(iwls)%area < 1.e2) cycle

! Count up individual land & sea cells and assign indices to them.
! Flag M and U points for association with land and/or sea cells.

     if (landsea_grid(iwls)%leaf_class <= 1) then
        nws = nws + 1
        landsea_grid(iwls)%wpt_sea = nws
     else  
        nwl = nwl + 1
        landsea_grid(iwls)%wpt_land = nwl
     endif

  enddo

! Now that final nws and nwl values are determined, 
! allocate separate sea and land cell arrays.  Copy data from temporary 
! arrays to these.

  call alloc_sea_grid (nws)
  call alloc_land_grid(nwl, nzg)

! Copy earth coordinate values from combined land + sea arrays to individual
! land & sea arrays.

! Copy data from combined land + sea arrays to individual land & sea arrays.

  do iwls = 2,nwls

     if (landsea_grid(iwls)%leaf_class <= 1) then

        iws = landsea_grid(iwls)%wpt_sea

! If sea cell has been skipped due to small size, iws = 0; skip cell

        if (iws < 2) cycle
        
        sea%area      (iws) = landsea_grid(iwls)%area
        sea%glatw     (iws) = landsea_grid(iwls)%glatw
        sea%glonw     (iws) = landsea_grid(iwls)%glonw
        sea%xew       (iws) = landsea_grid(iwls)%xew
        sea%yew       (iws) = landsea_grid(iwls)%yew
        sea%zew       (iws) = landsea_grid(iwls)%zew
        sea%topw      (iws) = landsea_grid(iwls)%topw
        sea%leaf_class(iws) = landsea_grid(iwls)%leaf_class
        
        itab_ws(iws)%iw     = itab_wls(iwls)%iw
        itab_ws(iws)%kw     = itab_wls(iwls)%kw
        itab_ws(iws)%npoly  = itab_wls(iwls)%npoly
        itab_ws(iws)%arf_iw = itab_wls(iwls)%arf_iw
        itab_ws(iws)%arf_kw = itab_wls(iwls)%arf_kw

        itab_ws(iws)%xem(1:maxnlspoly) = itab_wls(iwls)%xem(1:maxnlspoly)
        itab_ws(iws)%yem(1:maxnlspoly) = itab_wls(iwls)%yem(1:maxnlspoly)
        itab_ws(iws)%zem(1:maxnlspoly) = itab_wls(iwls)%zem(1:maxnlspoly)

     else

        iwl = landsea_grid(iwls)%wpt_land

! If land cell has been skipped due to small size, iwl = 0; skip cell

        if (iwl < 2) cycle

        land%area      (iwl) = landsea_grid(iwls)%area
        land%glatw     (iwl) = landsea_grid(iwls)%glatw
        land%glonw     (iwl) = landsea_grid(iwls)%glonw
        land%xew       (iwl) = landsea_grid(iwls)%xew
        land%yew       (iwl) = landsea_grid(iwls)%yew
        land%zew       (iwl) = landsea_grid(iwls)%zew
        land%topw      (iwl) = landsea_grid(iwls)%topw
        land%wnx       (iwl) = landsea_grid(iwls)%wnx
        land%wny       (iwl) = landsea_grid(iwls)%wny
        land%wnz       (iwl) = landsea_grid(iwls)%wnz
        land%leaf_class(iwl) = landsea_grid(iwls)%leaf_class

        itab_wl(iwl)%iw     = itab_wls(iwls)%iw
        itab_wl(iwl)%kw     = itab_wls(iwls)%kw
        itab_wl(iwl)%npoly  = itab_wls(iwls)%npoly
        itab_wl(iwl)%arf_iw = itab_wls(iwls)%arf_iw
        itab_wl(iwl)%arf_kw = itab_wls(iwls)%arf_kw

        itab_wl(iwl)%xem(1:maxnlspoly) = itab_wls(iwls)%xem(1:maxnlspoly)
        itab_wl(iwl)%yem(1:maxnlspoly) = itab_wls(iwls)%yem(1:maxnlspoly)
        itab_wl(iwl)%zem(1:maxnlspoly) = itab_wls(iwls)%zem(1:maxnlspoly)

     endif

  enddo

! Initialize soil textural class

  if (isoilflg == 2) then

! If soilflg == 2, fill land cells with default horizontally homogeneous
! soil textural class value defined in OLAMIN.
! User customization can be done here.

     land%ntext_soil(1:nzg,1)     = 0
     land%ntext_soil(1:nzg,2:nwl) = nslcon

  else

! If soilflg == 1, read soil textural class from database

     call land_database_read(nwl,  &
          land%glatw(:),           &
          land%glonw(:),           &
          soil_database,           &
          soil_database,           &
          'soil_text',             &
          idatq=landsea_grid(:)%idatq)

! Loop over all land cells (already defined and filled with leaf_class)

     do iwl = 2,nwl
        call fao_usda(landsea_grid(iwl)%idatq, datsoil)

! For now, assign single-level FAO textural class to all soil layers.

        do k = 1, nzg
           land%ntext_soil(k,iwl) = datsoil
        enddo

     enddo

  endif

! If iseagrid = 1, combine sea grid cells that share a
! sea-only atm grid column if all do so at the same kw level.

  if (iseagrid == 1) call combine_sea_cells()

! Deallocate mksfc arrays

  deallocate (itab_wls, landsea_grid)

end subroutine makesfc2

!===============================================================================

subroutine topo_init(nqa,topq,glatq,glonq,xeq,yeq,zeq)

use misc_coms,   only: io6, deltax
use consts_coms, only: pi1, pio180

implicit none

integer, intent(in) :: nqa
real, intent(in) :: glatq(nqa),glonq(nqa),xeq(nqa),yeq(nqa),zeq(nqa)

real, intent(out) :: topq(nqa)

integer :: iq

real :: hfwid
real :: hgt
real :: hfwid2

real :: r, r0

! Fill the TOPQ array with a default value of 0 or modify it as desired.
! If itopoflg is set to 1, these values will be overridden in the call to
! topo_database, which inputs a standard OLAM topography dataset.

topq(:) = 0.

hfwid = 10000.

r0 = pi1 / 9.

! dudhia expts
! hfwid = 5. * deltax

! hgt = 405.
! hgt = 1012.
! end dudhia expts

hfwid2 = hfwid**2

do iq = 2,nqa
   topq(iq) = 0.

!   topq(iq) = 200. * mod(iq,4)

! SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   topq(iq) = max(0.,hgt * hfwid2 / (hfwid2 + xeq(iq)**2) - 1.)
!   write(io6,*) 'topq ',iq,xeq(iq),topq(iq)
! TOPQ = 0 AT LARGE DISTANCE FROM HILL
! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SPECIAL WM5 EXPT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find lat/lon of current M point

!   r = sqrt((glonq(iq) * pio180 + 0.5 * pi1) ** 2 &
!          + (glatq(iq) * pio180 - pi1 / 6.) ** 2)

!   topq(iq) = max(0., 2000. * (1. - r / r0))

!   print*, 'topoinit ',iq,r,r0,topq(iq)

! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine topo_init

!===============================================================================

subroutine cont3sfc(z1,z2,z3,x1,x2,x3,y1,y2,y3,xmpat,ympat,kwpat,lpoly,npat,iw)

use mem_grid, only: nza, zm

implicit none

integer, intent(in) :: iw
real, intent(in) :: z1,z2,z3,x1,x2,x3,y1,y2,y3
real, intent(inout) :: xmpat(5,nza), ympat(5,nza)
integer, intent(out) :: kwpat(nza),lpoly(nza),npat

integer :: km, iflag
real :: xo1, xo2, yo1, yo2 ,contlev    ! "old" values of xmpat1,xmpat2,ympat1,ympat2

! Determine lowest zm level that exceeds z3, the lowest triangle vertex

km = 1
do while (z3 >= zm(km) .and. km < nza-1)
   km = km + 1
enddo
contlev = zm(km)
npat = 1

! If all 3 points are in the same contour interval

if (contlev >= z1 .or. contlev <= z3) then

   xmpat(1,npat) = x1
   ympat(1,npat) = y1
   xmpat(2,npat) = x2
   ympat(2,npat) = y2
   xmpat(3,npat) = x3
   ympat(3,npat) = y3

   kwpat(npat) = km

   lpoly(npat) = 3

   return
endif

! Draw all contour lines of value less than z1

iflag = 0

do while (contlev < z1)

   if (iflag > 0) then
      xo1 = xmpat(1,npat-1)
      yo1 = ympat(1,npat-1)
      xo2 = xmpat(2,npat-1)
      yo2 = ympat(2,npat-1)
   endif

   if ( abs(z1-z3) > 1.e-25 ) then
      xmpat(1,npat) = x3 + (x1-x3) * (contlev-z3) / (z1-z3)
      ympat(1,npat) = y3 + (y1-y3) * (contlev-z3) / (z1-z3)
   else
      xmpat(1,npat) = x3
      ympat(1,npat) = y3
   endif

   if (contlev < z2) then

      if ( abs(z2-z3) > 1.e-25 ) then
         xmpat(2,npat) = x3 + (x2-x3) * (contlev-z3) / (z2-z3)
         ympat(2,npat) = y3 + (y2-y3) * (contlev-z3) / (z2-z3)
      else
         xmpat(2,npat) = x3
         ympat(2,npat) = y3
      endif

      if (iflag == 0) then  ! lowest contour interval: z3 is a node
         xmpat(3,npat) = x3   
         ympat(3,npat) = y3

         kwpat(npat) = km
         lpoly(npat) = 3
      else
         xmpat(3,npat) = xo2
         ympat(3,npat) = yo2
         xmpat(4,npat) = xo1
         ympat(4,npat) = yo1

         kwpat(npat) = km
         lpoly(npat) = 4
      endif

      if (km == nza-1) then  ! highest model level: z1 and z2 are nodes
                           ! (should not actually happen)
         xmpat(3,npat) = x2
         ympat(3,npat) = y2
         xmpat(4,npat) = x1         
         ympat(4,npat) = y1

         kwpat(npat) = km
         lpoly(npat) = 4

         return
      elseif (zm(km+1) >= z1) then  ! highest contour interval:
                                    ! z1 and z2 are nodes
         xmpat(1,npat+1) = xmpat(2,npat)
         ympat(1,npat+1) = ympat(2,npat)
         xmpat(2,npat+1) = xmpat(1,npat)
         ympat(2,npat+1) = ympat(1,npat)

         xmpat(3,npat+1) = x1
         ympat(3,npat+1) = y1
         xmpat(4,npat+1) = x2         
         ympat(4,npat+1) = y2

         kwpat(npat+1) = km + 1
         lpoly(npat+1) = 4

         npat = npat + 1

         return
      endif
      
      iflag = 1
         
   else   
      
      if ( abs(z1-z2) > 1.e-25 ) then
         xmpat(2,npat) = x2 + (x1-x2) * (contlev-z2) / (z1-z2)
         ympat(2,npat) = y2 + (y1-y2) * (contlev-z2) / (z1-z2)
      else
         xmpat(2,npat) = x2
         ympat(2,npat) = y2
      endif
      
      if (iflag == 0) then  ! lowest contour color: z2 and z3 are nodes
         xmpat(3,npat) = x2
         ympat(3,npat) = y2
         xmpat(4,npat) = x3        
         ympat(4,npat) = y3         

         kwpat(npat) = km
         lpoly(npat) = 4
      elseif (iflag == 1) then ! switching from contlev < z2 to contlev > z2
         xmpat(3,npat) = x2
         ympat(3,npat) = y2
         xmpat(4,npat) = xo2
         ympat(4,npat) = yo2
         xmpat(5,npat) = xo1         
         ympat(5,npat) = yo1

         kwpat(npat) = km
         lpoly(npat) = 5
      else                     ! keeping to contlev > z2
         xmpat(3,npat) = xo2
         ympat(3,npat) = yo2
         xmpat(4,npat) = xo1
         ympat(4,npat) = yo1

         kwpat(npat) = km
         lpoly(npat) = 4
      endif
      
      if (km == nza-1) then  ! highest model level: z1 is a node 
         xmpat(3,npat) = x1      ! (should not actually happen)
         ympat(3,npat) = y1

         kwpat(npat) = km
         lpoly(npat) = 3

         return
      elseif (zm(km+1) >= z1) then  ! highest contour interval: z1 is a node
         xmpat(1,npat+1) = xmpat(2,npat)
         ympat(1,npat+1) = ympat(2,npat)

         xmpat(2,npat+1) = xmpat(1,npat)
         ympat(2,npat+1) = ympat(1,npat)

         xmpat(3,npat+1) = x1
         ympat(3,npat+1) = y1

         kwpat(npat+1) = km + 1
         lpoly(npat+1) = 3

         npat = npat + 1

         return
      endif

      iflag = 2
         
   endif

   npat = npat + 1
   km = km + 1
   contlev = zm(km)

enddo

end subroutine cont3sfc

!==========================================================================

subroutine reverse_polygon(xmpat,ympat,lpoly,npat)

use mem_grid, only: nza

implicit none

integer, intent(in) :: npat, lpoly(nza)

real, intent(inout) :: xmpat(5,nza), ympat(5,nza)

integer :: ipat, j, jc
real :: store

do ipat = 1,npat
   do j = 1,lpoly(ipat) / 2
      jc = lpoly(ipat) + 1 - j

      store = xmpat(j,ipat)
      xmpat(j,ipat) = xmpat(jc,ipat)
      xmpat(jc,ipat) = store

      store = ympat(j,ipat)
      ympat(j,ipat) = ympat(jc,ipat)
      ympat(jc,ipat) = store
   enddo
enddo

end subroutine reverse_polygon

!==========================================================================

subroutine oge_leafclass(ioge,leafclass)

! This subroutine maps the input ioge classes to a smaller set leafclass
! which represents the full set of LEAF-2 or LEAF-3 classes for which 
! LSP values are defined.

  implicit none

  integer, intent(in)  :: ioge
  integer, intent(out) :: leafclass

  integer :: catb(0:101)

! Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes
! (see leaf3_document).

!--------------------------------------------!
  data catb/ 0,                            & !
            19, 8, 4, 5, 6, 7, 9, 3,11,16, & !  0
            10, 2,17, 1, 0,12,13,14,18, 4, & ! 10
             4, 4,14,14, 6, 6, 4, 7, 7,15, & ! 20
            15, 6, 7, 7,15,16,16,16,16, 8, & ! 30
             8, 8,18,17,17,12,12, 7,10, 3, & ! 40
            10,10,11,14,18,18,18,18,13, 6, & ! 50
             5, 4,11,12, 0, 0, 0, 0, 3, 2, & ! 60
             3,20, 0,17,17,17, 4,14, 7, 3, & ! 70
             3, 3, 3, 3, 3, 3, 8,12, 7, 6, & ! 80
            18,15,15,15, 4, 5, 0, 0, 0, 0, & ! 90  ! 97 & 98 not used
            21                             / !
!--------------------------------------------!     ! 99 is Goode Homolosine 
                                                   !    empty space
!            1  2  3  4  5  6  7  8  9 10          ! 100 is missing data
                                                   ! Map all of these to ocean
                                                   ! (leafclass=0)
  leafclass = catb(ioge)

end subroutine oge_leafclass

!==========================================================================

subroutine fao_usda(ifao,iusda)

! This subroutine maps the input ifao soil classes to a smaller set iusda
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

  implicit none

  integer, intent(in) :: ifao
  integer, intent(out) :: iusda

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

!---------------------------------------------!
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
!---------------------------------------------!
!            1  2  3  4  5  6  7  8  9 10

  iusda = catb(ifao)

end subroutine fao_usda

!============================================================================

subroutine combine_sea_cells()

  use mem_sea,    only: sea, sea_vars, itab_ws
  use mem_ijtabs, only: itab_w, mrls
  use leaf_coms,  only: nwl
  use sea_coms,   only: nws
  use misc_coms,  only: io6
  use mem_mksfc,  only: itab_wls_vars
  use mem_leaf,   only: itab_wl
  use mem_sea,    only: itab_ws, sea
  use mem_grid,   only: nwa, xem, yem, zem, xew, yew, zew, &
                        glatw, glonw, arw0, topw

  implicit none

  type(itab_wls_vars), allocatable :: ltab_ws(:)

  type(sea_vars) :: sea_t

  integer :: jws, iws, iw, kw, jpoly, im, iwl

  integer, allocatable :: iwflag(:)

! Loop over all land cells, and flag the atmosphere grid IW columns that
! are attached to them

  allocate (iwflag(nwa))

  iwflag(:) = 0

  do iwl = 2,nwl
     iw = itab_wl(iwl)%iw

     if (iwflag(iw) == 0) iwflag(iw) = -1
  enddo

! Loop over all sea cells, and for those attached to an atmosphere grid
! column that has no land cell attached, compare sea cell kw levels
! within the same atmosphere grid column.  If any differ, flag atmosphere
! grid column with -1 value.

  do iws = 2,nws
     iw = itab_ws(iws)%iw
     kw = itab_ws(iws)%kw

     if (iwflag(iw) == 0) then
        iwflag(iw) = kw
     elseif (iwflag(iw) > 0 .and. iwflag(iw) /= kw) then
        iwflag(iw) = -1
     endif
  enddo

! Loop over all sea cells and count the number (jws) to be retained

  jws = 1

  do iws = 2,nws
     iw = itab_ws(iws)%iw

     if (iwflag(iw) == -1) then

! This sea cell shares an atmospheric grid column with one or more land
! cells and/or with sea cells having a different kw value, so increment jws
! to retain this sea cell as is

        jws = jws + 1

     elseif (iwflag(iw) /= -2) then

! This sea cell shares an atmospheric grid column only with other sea
! cells, and all have the same kw value, so increment jws only on first
! encounter with atmospheric cell, and reflag atmosphere cell (with -2 value)
! as having been visited

        jws = jws + 1
        iwflag(iw) = -2

     endif

  enddo

! Allocate temporary sea arrays to new size for permanent arrays

  allocate (ltab_ws(jws))

  allocate (sea_t%area      (jws))
  allocate (sea_t%glatw     (jws))
  allocate (sea_t%glonw     (jws))
  allocate (sea_t%xew       (jws))
  allocate (sea_t%yew       (jws))
  allocate (sea_t%zew       (jws))
  allocate (sea_t%topw      (jws))
  allocate (sea_t%leaf_class(jws))

! Loop again over all sea cells, using same logic for incrementing jws,
! but this time fill temporary sea array data values

  jws = 1

  do iws = 2,nws
     iw = itab_ws(iws)%iw

     if (iwflag(iw) == -1) then

! This sea cell shares an atmospheric grid column with one or more land
! cells and/or with sea cells having a different kw value, so increment jws
! to retain this sea cell as is

        jws = jws + 1

! Copy sea values to jws index in temporary arrays

        ltab_ws(jws)%iw     = itab_ws(iws)%iw
        ltab_ws(jws)%kw     = itab_ws(iws)%kw
        ltab_ws(jws)%npoly  = itab_ws(iws)%npoly
        ltab_ws(jws)%arf_iw = itab_ws(iws)%arf_iw
        ltab_ws(jws)%arf_kw = itab_ws(iws)%arf_kw

        do jpoly = 1,itab_ws(iws)%npoly
           ltab_ws(jws)%xem(jpoly) = itab_ws(iws)%xem(jpoly)
           ltab_ws(jws)%yem(jpoly) = itab_ws(iws)%yem(jpoly)
           ltab_ws(jws)%zem(jpoly) = itab_ws(iws)%zem(jpoly)
        enddo

! Copy sea values to jws index in temporary arrays

        sea_t%area (jws) = sea%area(iws)
        sea_t%glatw(jws) = sea%glatw(iws)
        sea_t%glonw(jws) = sea%glonw(iws)
        sea_t%xew  (jws) = sea%xew(iws)
        sea_t%yew  (jws) = sea%yew(iws)
        sea_t%zew  (jws) = sea%zew(iws)
        sea_t%topw (jws) = sea%topw(iws)

        sea_t%leaf_class(jws) = sea%leaf_class(iws)

     elseif (iwflag(iw) /= -3) then

! This sea cell shares an atmospheric grid column only with other sea
! cells, and all have the same kw value, so increment jws only on first
! encounter with atmospheric cell, and reflag atmosphere cell (with -3 value
! this time because -2 was used previously) as having been visited

        jws = jws + 1
        iwflag(iw) = -3

! Copy sea values to jws index in temporary arrays

        ltab_ws(jws)%iw     = itab_ws(iws)%iw
        ltab_ws(jws)%kw     = itab_ws(iws)%kw
        ltab_ws(jws)%npoly  = itab_w(iw)%npoly
        ltab_ws(jws)%arf_iw = 1.0
        ltab_ws(jws)%arf_kw = 1.0

        do jpoly = 1,itab_w(iw)%npoly
           im = itab_w(iw)%im(jpoly)

           ltab_ws(jws)%xem(jpoly) = xem(im)
           ltab_ws(jws)%yem(jpoly) = yem(im)
           ltab_ws(jws)%zem(jpoly) = zem(im)
        enddo

! Copy sea values to jws index in temporary arrays

        sea_t%area (jws) = arw0(iw)
        sea_t%glatw(jws) = glatw(iw)
        sea_t%glonw(jws) = glonw(iw)
        sea_t%xew  (jws) = xew(iw)
        sea_t%yew  (jws) = yew(iw)
        sea_t%zew  (jws) = zew(iw)
        sea_t%topw (jws) = topw(iw)

        sea_t%leaf_class(jws) = sea%leaf_class(iws)

     endif

  enddo

  nws = jws

! Move temporary arrays to original sea arrays

  call move_alloc(ltab_ws, itab_ws)

  call move_alloc(sea_t%area,       sea%area)
  call move_alloc(sea_t%glatw,      sea%glatw)
  call move_alloc(sea_t%glonw,      sea%glonw)
  call move_alloc(sea_t%xew,        sea%xew)
  call move_alloc(sea_t%yew,        sea%yew)
  call move_alloc(sea_t%zew,        sea%zew)
  call move_alloc(sea_t%topw,       sea%topw)
  call move_alloc(sea_t%leaf_class, sea%leaf_class)

end subroutine combine_sea_cells

!==========================================================================

subroutine landfile_write()

  use max_dims,   only: maxnlspoly, pathlen
  use leaf_coms,  only: nzg, nwl, landusefile, slz
  use mem_leaf,   only: land, itab_wl
  use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer :: iclobber1 = 1
  integer :: ndims, idims(2)
  integer :: iwl
  real    :: rscr(maxnlspoly,nwl)

  character(pathlen) :: flnm

! Open landfile

  flnm = trim(landusefile)//'.h5'

  call shdf5_open(flnm,'W',iclobber1)

! Write land grid quantities to landfile

  ndims = 1
  idims(1) = 1

  call shdf5_orec(ndims, idims, 'nwl'   , ivars=nwl)
  call shdf5_orec(ndims, idims, 'nzg'   , ivars=nzg)

  ndims = 1
  idims(1) = nzg

  call shdf5_orec(ndims, idims, 'slz', rvara=slz)

  idims(1) = nwl

  call shdf5_orec(ndims, idims, 'iw'        , ivara=itab_wl(:)%iw)
  call shdf5_orec(ndims, idims, 'kw'        , ivara=itab_wl(:)%kw)
  call shdf5_orec(ndims, idims, 'nlpoly'    , ivara=itab_wl(:)%npoly)
  call shdf5_orec(ndims, idims, 'arf_iw'    , rvara=itab_wl(:)%arf_iw)
  call shdf5_orec(ndims, idims, 'arf_kw'    , rvara=itab_wl(:)%arf_kw)

  call shdf5_orec(ndims, idims, 'land_area' , rvara=land%area)
  call shdf5_orec(ndims, idims, 'glatwl'    , rvara=land%glatw)
  call shdf5_orec(ndims, idims, 'glonwl'    , rvara=land%glonw)
  call shdf5_orec(ndims, idims, 'xewl'      , rvara=land%xew)
  call shdf5_orec(ndims, idims, 'yewl'      , rvara=land%yew)
  call shdf5_orec(ndims, idims, 'zewl'      , rvara=land%zew)
  call shdf5_orec(ndims, idims, 'topwl'     , rvara=land%topw)
  call shdf5_orec(ndims, idims, 'wnxl'      , rvara=land%wnx)
  call shdf5_orec(ndims, idims, 'wnyl'      , rvara=land%wny)
  call shdf5_orec(ndims, idims, 'wnzl'      , rvara=land%wnz)
  call shdf5_orec(ndims, idims, 'leaf_class', ivara=land%leaf_class)

  ndims = 2
  idims(1) = nzg
  idims(2) = nwl

  call shdf5_orec(ndims, idims, 'ntext_soil', ivara=land%ntext_soil)

  ndims = 2
  idims(1) = maxnlspoly
  idims(2) = nwl

! Copy itab_wl arrays to scratch array rscr for output

  do iwl = 1,nwl
     rscr(1:maxnlspoly,iwl) = itab_wl(iwl)%xem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims,idims,'itab_wl%xem',rvara=rscr)

  do iwl = 1,nwl
     rscr(1:maxnlspoly,iwl) = itab_wl(iwl)%yem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims,idims,'itab_wl%yem',rvara=rscr)

  do iwl = 1,nwl
     rscr(1:maxnlspoly,iwl) = itab_wl(iwl)%zem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims,idims,'itab_wl%zem',rvara=rscr)

  call shdf5_close()

end subroutine landfile_write

!==========================================================================

subroutine seafile_write()

  use max_dims,   only: maxnlspoly, pathlen
  use sea_coms,   only: nws, seafile
  use mem_sea,    only: sea, itab_ws
  use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer :: iclobber1 = 1
  integer :: ndims, idims(2)
  integer :: iws
  real    :: rscr(maxnlspoly,nws)

  character(pathlen) :: flnm

! Open seafile

  flnm=trim(seafile)//'.h5'

  call shdf5_open(flnm,'W',iclobber1)

! Write sea grid quantities to seafile

  ndims = 1
  idims(1) = 1

  call shdf5_orec(ndims, idims, 'nws'   , ivars=nws)

  ndims = 1
  idims(1) = nws

  call shdf5_orec(ndims, idims, 'iw'        , ivara=itab_ws(:)%iw)
  call shdf5_orec(ndims, idims, 'kw'        , ivara=itab_ws(:)%kw)
  call shdf5_orec(ndims, idims, 'nspoly'    , ivara=itab_ws(:)%npoly)
  call shdf5_orec(ndims, idims, 'arf_iw'    , rvara=itab_ws(:)%arf_iw)
  call shdf5_orec(ndims, idims, 'arf_kw'    , rvara=itab_ws(:)%arf_kw)

  call shdf5_orec(ndims, idims, 'sea_area'  , rvara=sea%area)
  call shdf5_orec(ndims, idims, 'glatws'    , rvara=sea%glatw)
  call shdf5_orec(ndims, idims, 'glonws'    , rvara=sea%glonw)
  call shdf5_orec(ndims, idims, 'xews'      , rvara=sea%xew)
  call shdf5_orec(ndims, idims, 'yews'      , rvara=sea%yew)
  call shdf5_orec(ndims, idims, 'zews'      , rvara=sea%zew)
  call shdf5_orec(ndims, idims, 'topws'     , rvara=sea%topw)
  call shdf5_orec(ndims, idims, 'leaf_class', ivara=sea%leaf_class)

  ndims = 2
  idims(1) = maxnlspoly
  idims(2) = nws

! Copy itab_ws arrays to temporary array rscr for output

  do iws = 1,nws
     rscr(1:maxnlspoly,iws) = itab_ws(iws)%xem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims, idims, 'itab_ws%xem', rvara=rscr)

  do iws = 1,nws
     rscr(1:maxnlspoly,iws) = itab_ws(iws)%yem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims, idims, 'itab_ws%yem', rvara=rscr)

  do iws = 1,nws
     rscr(1:maxnlspoly,iws) = itab_ws(iws)%zem(1:maxnlspoly)
  enddo

  call shdf5_orec(ndims, idims, 'itab_ws%zem', rvara=rscr)

  call shdf5_close()

end subroutine seafile_write

!==========================================================================

subroutine landfile_read()

  use max_dims,   only: maxnlspoly, pathlen
  use leaf_coms,  only: nzg, nwl, mwl, landusefile, slz
  use mem_leaf,   only: land, itab_wl, alloc_land_grid
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: iwl
  real, allocatable  :: rscr(:,:)

!-------------------------------------------------------------------------------
! STEP 1: Open landfile and read 2 array dimensions
!-------------------------------------------------------------------------------

  flnm = trim(landusefile)//'.h5'

  write(io6,*) 'Checking landfile ',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'Landfile was not found - stopping run'
     stop 'stop: no landfile'
  endif

  call shdf5_open(flnm,'R')

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'nwl', ivars=nwl)
  call shdf5_irec(ndims, idims, 'nzg', ivars=nzg)

  write(io6, '(/,a)')   '==============================================='
  write(io6, '(a)')     'Reading from landfile:'
  write(io6, '(a,4i8)') '  nwl, nzg = ', nwl, nzg
  write(io6, '(a,/)')   '==============================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate land grid arrays
!-------------------------------------------------------------------------------

  call alloc_land_grid(nwl,nzg)
  allocate (rscr(maxnlspoly,nwl))

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from landfile
!-------------------------------------------------------------------------------

  ndims = 1
  idims(1) = nzg

  call shdf5_irec(ndims, idims, 'slz', rvara=slz)

! Read arrays from landfile

  idims(1) = nwl

  call shdf5_irec(ndims, idims, 'iw'        , ivara=itab_wl(:)%iw)
  call shdf5_irec(ndims, idims, 'kw'        , ivara=itab_wl(:)%kw)
  call shdf5_irec(ndims, idims, 'nlpoly'    , ivara=itab_wl(:)%npoly)
  call shdf5_irec(ndims, idims, 'arf_iw'    , rvara=itab_wl(:)%arf_iw)
  call shdf5_irec(ndims, idims, 'arf_kw'    , rvara=itab_wl(:)%arf_kw)

  call shdf5_irec(ndims, idims, 'land_area' , rvara=land%area)
  call shdf5_irec(ndims, idims, 'glatwl'    , rvara=land%glatw)
  call shdf5_irec(ndims, idims, 'glonwl'    , rvara=land%glonw)
  call shdf5_irec(ndims, idims, 'xewl'      , rvara=land%xew)
  call shdf5_irec(ndims, idims, 'yewl'      , rvara=land%yew)
  call shdf5_irec(ndims, idims, 'zewl'      , rvara=land%zew)
  call shdf5_irec(ndims, idims, 'topwl'     , rvara=land%topw)
  call shdf5_irec(ndims, idims, 'wnxl'      , rvara=land%wnx)
  call shdf5_irec(ndims, idims, 'wnyl'      , rvara=land%wny)
  call shdf5_irec(ndims, idims, 'wnzl'      , rvara=land%wnz)
  call shdf5_irec(ndims, idims, 'leaf_class', ivara=land%leaf_class)

  ndims = 2
  idims(1) = nzg
  idims(2) = nwl

  call shdf5_irec(ndims, idims, 'ntext_soil', ivara=land%ntext_soil)

  ndims = 2
  idims(1) = maxnlspoly
  idims(2) = nwl

  call shdf5_irec(ndims,idims,'itab_wl%xem',rvara=rscr)

  do iwl = 1,nwl
     itab_wl(iwl)%xem(1:maxnlspoly) = rscr(1:maxnlspoly,iwl)
  enddo

  call shdf5_irec(ndims,idims,'itab_wl%yem',rvara=rscr)

  do iwl = 1,nwl
     itab_wl(iwl)%yem(1:maxnlspoly) = rscr(1:maxnlspoly,iwl)
  enddo

  call shdf5_irec(ndims,idims,'itab_wl%zem',rvara=rscr)

  do iwl = 1,nwl
     itab_wl(iwl)%zem(1:maxnlspoly) = rscr(1:maxnlspoly,iwl)
  enddo

  call shdf5_close()
  deallocate(rscr)

  mwl = nwl

end subroutine landfile_read

!==========================================================================

subroutine seafile_read()

  use max_dims,   only: maxnlspoly, pathlen
  use sea_coms,   only: nws, mws, seafile, iseagrid
  use mem_sea,    only: sea, itab_ws, alloc_sea_grid
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: iws
  real, allocatable  :: rscr(:,:)

!-------------------------------------------------------------------------------
! STEP 1: Open SEAFILE and read 1 array dimension
!-------------------------------------------------------------------------------

  flnm = trim(seafile)//'.h5'

  write(io6,*) 'Checking seafile ', trim(flnm)

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'SEAFILE was not found - stopping run'
     stop 'stop: no seafile'
  endif

  call shdf5_open(flnm,'R')

  ndims = 1
  idims(1) = 1

  call shdf5_irec(ndims, idims, 'nws', ivars=nws)

  write(io6, '(/,a)')   '====================================='
  write(io6, '(a)')     'Reading from seafile:'
  write(io6, '(a,3i8)') '  nws = ', nws
  write(io6, '(a,/)')   '====================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate SEA grid arrays
!-------------------------------------------------------------------------------

  call alloc_sea_grid(nws)

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from seafile
!-------------------------------------------------------------------------------

  ndims = 1
  idims(1) = nws

  call shdf5_irec(ndims, idims, 'iw'        , ivara=itab_ws(:)%iw)
  call shdf5_irec(ndims, idims, 'kw'        , ivara=itab_ws(:)%kw)
  call shdf5_irec(ndims, idims, 'nspoly'    , ivara=itab_ws(:)%npoly)
  call shdf5_irec(ndims, idims, 'arf_iw'    , rvara=itab_ws(:)%arf_iw)
  call shdf5_irec(ndims, idims, 'arf_kw'    , rvara=itab_ws(:)%arf_kw)

  call shdf5_irec(ndims, idims, 'sea_area'  , rvara=sea%area)
  call shdf5_irec(ndims, idims, 'glatws'    , rvara=sea%glatw)
  call shdf5_irec(ndims, idims, 'glonws'    , rvara=sea%glonw)
  call shdf5_irec(ndims, idims, 'xews'      , rvara=sea%xew)
  call shdf5_irec(ndims, idims, 'yews'      , rvara=sea%yew)
  call shdf5_irec(ndims, idims, 'zews'      , rvara=sea%zew)
  call shdf5_irec(ndims, idims, 'topws'     , rvara=sea%topw)
  call shdf5_irec(ndims, idims, 'leaf_class', ivara=sea%leaf_class)

  allocate (rscr(maxnlspoly,nws))

  ndims = 2
  idims(1) = maxnlspoly
  idims(2) = nws

  rscr(:,:) = 0

  call shdf5_irec(ndims, idims, 'itab_ws%xem', rvara=rscr)

  do iws = 1,nws
     itab_ws(iws)%xem(1:maxnlspoly) = rscr(1:maxnlspoly,iws)
  enddo

  call shdf5_irec(ndims, idims, 'itab_ws%yem', rvara=rscr)

  do iws = 1,nws
     itab_ws(iws)%yem(1:maxnlspoly) = rscr(1:maxnlspoly,iws)
  enddo

  call shdf5_irec(ndims, idims, 'itab_ws%zem', rvara=rscr)

  do iws = 1,nws
     itab_ws(iws)%zem(1:maxnlspoly) = rscr(1:maxnlspoly,iws)
  enddo

  deallocate(rscr)

  call shdf5_close()

  mws = nws

end subroutine seafile_read

!==========================================================================

subroutine landfile_read_oldgrid()

  use max_dims,   only: maxnlspoly, pathlen
  use leaf_coms,  only: landusefile
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6
  use mem_addgrid,only: nwl_og, nzg_og, xewl_og, yewl_og, zewl_og, ntext_soil_og

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: iwl_og

  flnm = trim(landusefile)//'-OG'//'.h5'

  write(io6,*) 'Checking OLD landfile ',trim(flnm)

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'OLD Landfile was not found - stopping run'
     stop 'stop: no old landfile'
  endif

  call shdf5_open(flnm,'R')

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'nwl', ivars=nwl_og)
  call shdf5_irec(ndims, idims, 'nzg', ivars=nzg_og)

  allocate (xewl_og(nwl_og))
  allocate (yewl_og(nwl_og))
  allocate (zewl_og(nwl_og))

  allocate (ntext_soil_og(nzg_og,nwl_og))

  ndims = 1
  idims(1) = nwl_og

  call shdf5_irec(ndims, idims, 'xewl', rvara=xewl_og)
  call shdf5_irec(ndims, idims, 'yewl', rvara=yewl_og)
  call shdf5_irec(ndims, idims, 'zewl', rvara=zewl_og)

  ndims = 2
  idims(1) = nzg_og
  idims(2) = nwl_og

  call shdf5_irec(ndims, idims, 'ntext_soil', ivara=ntext_soil_og)

  call shdf5_close()

end subroutine landfile_read_oldgrid

!==========================================================================

subroutine seafile_read_oldgrid()

  use max_dims,   only: maxnlspoly, pathlen
  use sea_coms,   only: seafile
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6
  use mem_addgrid,only: nws_og, xews_og, yews_og, zews_og

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: iws

  flnm = trim(seafile)//'-OG'//'.h5'

  write(io6,*) 'Checking OLD seafile ', trim(flnm)

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'OLD SEAFILE was not found - stopping run'
     stop 'stop: no old seafile'
  endif

  call shdf5_open(flnm,'R')

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'nws', ivars=nws_og)

  allocate (xews_og(nws_og))
  allocate (yews_og(nws_og))
  allocate (zews_og(nws_og))

  idims(1) = nws_og

  call shdf5_irec(ndims, idims, 'xews', rvara=xews_og)
  call shdf5_irec(ndims, idims, 'yews', rvara=yews_og)
  call shdf5_irec(ndims, idims, 'zews', rvara=zews_og)

  call shdf5_close()

end subroutine seafile_read_oldgrid
