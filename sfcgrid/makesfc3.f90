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
subroutine makesfc3()

! This subroutine generates the SURFACE GRID FILE for OLAM runs. 

  use mem_grid,    only: nza, nma, nwa, zm, &
                         xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                         arw0, glatw, glonw, glatm, glonm, topm, topw, dnu, dnv

  use mem_ijtabs,  only: itab_w, jtab_w, jtw_grid, jtw_lbcp

  use misc_coms,   only: io6, itopoflg, rinit, topo_database, mdomain

  use leaf_coms,   only: dt_leaf, nzs, nvgcon, isfcl, ivegflg, isoilflg, &
                         soil_database, veg_database

  use consts_coms, only: erad, eradi, piu180, r8

  use land_db,     only: land_database_read

  use soilgrids_db, only: soilgrids_read

  use max_dims,    only: maxnlspoly

  use mem_sfcg,    only: nsfcgrids, sfcgfile, nmsfc, nvsfc, nwsfc, mmsfc, mvsfc, mwsfc, &
                         itab_msfc, itab_vsfc, itab_wsfc, sfcg, &
                         itab_msfc_vars, itab_vsfc_vars, itab_wsfc_vars, resize_sfcgrid

  use mem_land,    only: nland, onland, itab_land, land, nzg, landgrid_dztop, &
                         landgrid_depth, slz, dslz, dslzo2, dslzi, slzt, &
                         alloc_landcol

  use mem_lake,    only: nlake, onlake, itab_lake, lake 
  use mem_sea,     only: nsea,  onsea,  itab_sea,  sea
  use oname_coms,  only: nl

  implicit none

  integer :: k, im, iw, j
  integer :: nqa, iqa, kw_sea, kw_land
  integer :: npoly, iland, inew
  integer :: nmsfc_inc, nvsfc_inc, nwsfc_inc
  integer :: imsfc, ivsfc, iwsfc
  integer :: leafclass
  integer :: nasfc, sumleafclassm

  integer, parameter :: nsgdata = 1 ! number of SoilGrids datasets to read

  integer, allocatable :: miw(:), miv(:,:), mim(:,:,:)

  real, allocatable :: qlat(:),qlon(:),topq(:)
  real, allocatable :: xeq(:),yeq(:),zeq(:)
  integer, allocatable :: ioge(:), leaf_class(:)

  type(itab_msfc_vars), allocatable :: itab_msfc_temp(:)
  type(itab_vsfc_vars), allocatable :: itab_vsfc_temp(:)
  type(itab_wsfc_vars), allocatable :: itab_wsfc_temp(:)

  real(r8), allocatable :: tot_area(:)
  real(r8) :: area_tot

  integer, allocatable :: leaf_classm(:), iwnew(:)

  integer :: jj,iwn
  real :: topwmin, topwmax, fracz, topr

  integer :: iter
  real :: dz, srati

  integer :: nzga, kk
  real :: thick

  ! At this point in the process of surface grid initialization, the surface
  ! grid has been copied in Delaunay form (triangle mesh) from the atmospheric
  ! grid (by a call to subroutine copy_sfc_tri_grid).  This copy may, but need
  ! not, include all local mesh refinements applied to the atmosphere grid.

  ! If any additional local mesh refinements are to be made on the surface grid
  ! independent of those on the atmospheric grid, perform them next.

  if (nsfcgrids > 0) call spawn_nest_sfc()

  ! Generate Voronoi form of surface grid and compute its geometric properties.

  call voronoi_sfc()
  call grid_geometry_hex_sfc()

  ! Topography height and surface type (ocean, lake, land) will be filled not
  ! only for Voronoi (hexagon) cell centers (W points), but also for the
  ! vertices (M points).  In preparation, allocate and define arrays of spatial
  ! coordinates that are composites of both sets of points (M + W). 

  nqa = nmsfc + nwsfc - 1

  allocate (qlat(nqa), qlon(nqa), topq(nqa))
  allocate ( xeq(nqa),  yeq(nqa),  zeq(nqa))
  allocate (ioge(nqa), leaf_class(nqa))
  allocate (leaf_classm(nmsfc))

  qlat(1:nmsfc) = sfcg%glatm(1:nmsfc)
  qlon(1:nmsfc) = sfcg%glonm(1:nmsfc)
  xeq (1:nmsfc) = sfcg%xem  (1:nmsfc)
  yeq (1:nmsfc) = sfcg%yem  (1:nmsfc)
  zeq (1:nmsfc) = sfcg%zem  (1:nmsfc)

  qlat(nmsfc+1:nqa) = sfcg%glatw(2:nwsfc)
  qlon(nmsfc+1:nqa) = sfcg%glonw(2:nwsfc)
  xeq (nmsfc+1:nqa) = sfcg%xew  (2:nwsfc)
  yeq (nmsfc+1:nqa) = sfcg%yew  (2:nwsfc)
  zeq (nmsfc+1:nqa) = sfcg%zew  (2:nwsfc)

  ! Interpolate or assign topography height to (M + W) points of surface grid

  if (itopoflg == 1) then  ! from database
     call land_database_read(nqa, qlat, qlon, &
         topo_database, topo_database, 'topo', datq=topq)
  else
     call topo_init(nqa,topq,qlat,qlon,xeq,yeq,zeq)
  endif

  ! Read or assign surface type to (M + W) points of surface grid

  if (ivegflg == 1) then  ! from database
     call land_database_read(nqa, qlat, qlon, &
          veg_database, veg_database, 'leaf_class', idatq=ioge)

     do iqa = 2,nqa
        call oge_leafclass(ioge(iqa), leaf_class(iqa))
     enddo
  else  ! from user-defined values in namelist file
     leaf_class(2:nqa) = nvgcon
  endif

  ! Loop over all (M + W) points of surface grid

  do iqa = 2,nqa

     ! Prevent TOPQ lower than lowest model level zm(1)

     topq(iqa) = max(topq(iqa),zm(1))

     ! If point is ocean but has nonzero topography height, print message
     ! and reset topography to zero.

     if (leaf_class(iqa) == 0 .and. topq(iqa) /= 0.) then
        print*, 'topq height of ocean point is nonzero: resetting to zero '
        topq(iqa) = 0.
     endif

  enddo

  ! Copy topography height and surface type from temporary composite arrays to
  ! individual W and M arrays of surface grid, and deallocate composite arrays.

  sfcg%topm(2:nmsfc) = topq(2:nmsfc)
  sfcg%topw(2:nwsfc) = topq(nmsfc+1:nqa)

  sfcg%ioge(2:nwsfc) = ioge(nmsfc+1:nqa)

       leaf_classm(2:nmsfc) = leaf_class(2:nmsfc)
  sfcg%leaf_class (2:nwsfc) = leaf_class(nmsfc+1:nqa)

  deallocate (qlat, qlon, topq)
  deallocate (xeq, yeq, zeq)
  deallocate (ioge, leaf_class)

  ! Allocate arrays to store list of M points during subdivision of surface
  ! grid hex cells

  allocate (miw(nwsfc))      ; miw(:)     = 0
  allocate (miv(nvsfc,nza))  ; miv(:,:)   = 0
  allocate (mim(nwsfc,7,nza)); mim(:,:,:) = 0

  ! Set amount of extra array space in case cell subdivision is done

  nmsfc_inc = 10 * nmsfc
  nvsfc_inc = 10 * nvsfc
  nwsfc_inc = 10 * nwsfc

  ! Initialize grid size counters at current values

  mmsfc = nmsfc
  mvsfc = nvsfc
  mwsfc = nwsfc

  ! Vertical index for sea (ocean) cells

  kw_sea = 2
  do while(zm(kw_sea) < 1.) ! 1-meter threshold
     kw_sea = kw_sea + 1
  enddo

  ! Initialize ATM grid topography values to large negative ("missing") value

  topw(:) = -2.e3
  topm(:) = -2.e3

  ! Loop over surface grid W cells

  do iwsfc = 2,nwsfc

     npoly = itab_wsfc(iwsfc)%npoly

     if (itab_wsfc(iwsfc)%ivoronoi >= 2) then  ! (iwatm = 1)

        ! This surface cell is flagged to remain a Voronoi cell (i.e., to not
        ! be subdivided).  Determine overlay with atmospheric grid.

        call sfc_atm_hex_overlay(iwsfc)

        if (sfcg%leaf_class(iwsfc) == 0) then

           itab_wsfc(iwsfc)%kwatm(1: ) = kw_sea

           sfcg%wnx(iwsfc) = sfcg%xew(iwsfc) * eradi
           sfcg%wny(iwsfc) = sfcg%yew(iwsfc) * eradi
           sfcg%wnz(iwsfc) = sfcg%zew(iwsfc) * eradi

        else

           sfcg%wnx(iwsfc) = sfcg%xew(iwsfc) * eradi
           sfcg%wny(iwsfc) = sfcg%yew(iwsfc) * eradi
           sfcg%wnz(iwsfc) = sfcg%zew(iwsfc) * eradi

           kw_land = 2
           do while(zm(kw_land) < 1. + sfcg%topw(iwsfc)) ! 1-meter threshold
              kw_land = kw_land + 1
           enddo
           itab_wsfc(iwsfc)%kwatm = kw_land

        endif

     else

        ! This cell is not flagged to remain a Voronoi cell, which means that
        ! in its present Voronoi form, it should exactly match an atmospheric
        ! grid column in both index value and location (until subdivision, if
        ! any, is done in subroutine sfc_grid_topocut).  Check here for exact
        ! coincidence.

        if (abs(sfcg%xew(iwsfc) - xew(iwsfc)) > 1. .or. &
            abs(sfcg%yew(iwsfc) - yew(iwsfc)) > 1. .or. &
            abs(sfcg%zew(iwsfc) - zew(iwsfc)) > 1.) then

           print*, 'Sfc grid non-Voronoi cell does not match ATM cell ',iwsfc
           print*, 'sfc2 ',itab_wsfc(iwsfc)%npoly,sfcg%glatw(iwsfc),sfcg%glonw(iwsfc)
           stop 'stop 101 '
        endif

        sumleafclassm = 0
        do j = 1,npoly
           imsfc = itab_wsfc(iwsfc)%imn(j)
           iwn = itab_wsfc(iwsfc)%iwn(j)

           if (abs(sfcg%xem(imsfc) - xem(imsfc)) > 2. .or. &
               abs(sfcg%yem(imsfc) - yem(imsfc)) > 2. .or. &
               abs(sfcg%zem(imsfc) - zem(imsfc)) > 2.) then

              print*, 'Sfc grid non-Voronoi cell does not match ATM cell ',iwsfc,j
              print*, 'sfc3 ',npoly,sfcg%glatw(iwsfc),sfcg%glonw(iwsfc)
              print*, 'sfc4 ',abs(sfcg%xem(imsfc) - xem(imsfc)), &
                              abs(sfcg%yem(imsfc) - yem(imsfc)), &
                              abs(sfcg%zem(imsfc) - zem(imsfc))

              do jj = 1,npoly
                 iwn = itab_wsfc(iwsfc)%iwn(jj)
                 print*, 'sfc5 ',itab_wsfc(iwn)%npoly
              enddo

              stop 'stop 102 '
           endif

           sumleafclassm = sumleafclassm + leaf_classm(imsfc)
        enddo

        ! Perform topography adjustment for the W and M points of this surface
        ! cell prior to topocutting

        topwmin = 1.e9
        topwmax = -1.e9

        do j = 0,npoly
           if (j == 0) then
              topr = sfcg%topw(iwsfc)
           else
              imsfc = itab_wsfc(iwsfc)%imn(j)
              topr = sfcg%topm(imsfc)
           endif

           ! Loop over KM levels

           do k = 1,nza-1

              ! Check whether topr is at or above current KM level

              if (zm(k) <= topr .and. zm(k+1) >= topr) then

                 ! Check topr within height range of KT level

                 fracz = (topr - zm(k)) / (zm(k+1) - zm(k))

                 ! If topr is in lowest or highest 2% of KT level, move it to the limit.      

                 if (fracz < .02) then
                    topr = zm(k)
                 elseif (fracz > .98) then
                    topr = zm(k+1)
                 endif

                 exit

              endif

           enddo

           ! Copy adjusted topography back to SFC grid W or M point, and copy
           ! M point topography to coincident ATM grid M point

           if (j == 0) then
              sfcg%topw(iwsfc) = topr
           else
              sfcg%topm(imsfc) = topr
              topm(imsfc) = sfcg%topm(imsfc)

              if (topwmin > sfcg%topm(imsfc)) topwmin = sfcg%topm(imsfc)
              if (topwmax < sfcg%topm(imsfc)) topwmax = sfcg%topm(imsfc)
           endif

        enddo

        ! Constrain SFC grid TOPW by max/min range of adjacent TOPM values, and
        ! copy to coincident ATM grid W point

        sfcg%topw(iwsfc) = max(topwmin,min(topwmax,sfcg%topw(iwsfc)))
        topw(iwsfc) = sfcg%topw(iwsfc)

        ! Check surface type of this surface grid cell

        if (sfcg%leaf_class(iwsfc) == 0 .and. sumleafclassm == 0) then

           ! Surface grid W cell and all surrounding M points are ocean.  Cell
           ! will not be subdivided.  Set ivoronoi flag = 1 to indicate that
           ! no topocutting is done for this cell.

           itab_wsfc(iwsfc)%kwatm(1) = kw_sea
           itab_wsfc(iwsfc)%arc(1) = sfcg%area(iwsfc)
           itab_wsfc(iwsfc)%ivoronoi = 1

        else 

           ! Surface grid W cell is not flagged to remain a Voronoi cell and
           ! it is not all ocean, so cell will be subdivided, based in part on
           ! how the topographic surface intersects the horizontal layers of
           ! the atmospheric grid.  The subdivision of the surface grid follows
           ! the lines of this intersection such that each surface grid cell
           ! interfaces with only a single atmospheric cell (single by virtue
           ! of both its horizontal and vertical index).  In preparation for
           ! the new surface grid cells generated by this process, allocate
           ! more space in the surface grid arrays.

           if (mmsfc > size(sfcg%xem) - 100  .or. &
               mvsfc > size(sfcg%xev) - 100  .or. &
               mwsfc > size(sfcg%xew) - 100) then

              allocate (itab_msfc_temp(mmsfc + nmsfc_inc))
              itab_msfc_temp(1:mmsfc) = itab_msfc(1:mmsfc)
              call move_alloc(itab_msfc_temp, itab_msfc)

              allocate (itab_vsfc_temp(mvsfc + nvsfc_inc))
              itab_vsfc_temp(1:mvsfc) = itab_vsfc(1:mvsfc)
              call move_alloc(itab_vsfc_temp, itab_vsfc)

              allocate (itab_wsfc_temp(mwsfc + nwsfc_inc))
              itab_wsfc_temp(1:mwsfc) = itab_wsfc(1:mwsfc)
              call move_alloc(itab_wsfc_temp, itab_wsfc)

              call resize_sfcgrid(mmsfc + nmsfc_inc, mvsfc + nvsfc_inc, mwsfc + nwsfc_inc)
           endif

           ! Proceed with subdivision by atmospheric grid cut cells.

           call sfc_grid_topocut(iwsfc,nwsfc,nvsfc,mwsfc,mmsfc,miw,miv,mim)

        endif

     endif

  enddo

  ! For any surface grid cells that were added or altered in sfc_grid_topocut,
  ! fill oge and leaf_class values

  if (mwsfc > nwsfc) then

     if (ivegflg == 2) then

        ! If ivegflg == 2, fill sea/land cell areas and IW values, plus
        ! leaf_class for land cells, from default value defined in OLAMIN.
        ! User customization can be done here.

        do iwsfc = 2,mwsfc
           if (itab_wsfc(iwsfc)%ivoronoi == 0) then
              sfcg%leaf_class(iwsfc) = nvgcon
           endif
        enddo

     else

        allocate (qlat(mwsfc), qlon(mwsfc))
        allocate (ioge(mwsfc), leaf_class(mwsfc))
        nqa = 1

        do iwsfc = 2,mwsfc
           if (itab_wsfc(iwsfc)%ivoronoi == 0) then
              nqa = nqa + 1
              qlat(nqa) = sfcg%glatw(iwsfc)
              qlon(nqa) = sfcg%glonw(iwsfc)
           endif
        enddo

        ! If ivegflg == 1, fill sea/land cells with leaf_class from leaf database

        call land_database_read(nqa, &
             qlat(1:nqa), &
             qlon(1:nqa), &
             veg_database, &
             veg_database, &
             'leaf_class', &
             idatq=ioge(1:nqa))

        iqa = 1
        do iwsfc = 2,mwsfc
           if (itab_wsfc(iwsfc)%ivoronoi == 0) then
              iqa = iqa + 1
              call oge_leafclass(ioge(iqa), leafclass)
              sfcg%ioge(iwsfc) = ioge(iqa)
              sfcg%leaf_class(iwsfc) = leafclass
           endif
        enddo

        deallocate (qlat, qlon, ioge, leaf_class)
     endif

  endif

  ! call plot_fields(0)  ! Only for LEAF_CLASS

  ! Check all ATM grid topography values to make sure they have been set

  do im = 2,nma
     if (topm(im) < -1.e3) then
        print*, 'topm not set ',im,topm(im),glatm(im),glonm(im)
        stop 'stop topm '
     endif
  enddo

  do iw = 2,nwa
     if (topw(iw) < -1.e3) then
        print*, 'topw not set ',iw,topw(iw),glatw(iw),glonw(iw)
        stop 'stop topw '
     endif
  enddo

  ! Now, all surface grid cells have been generated.  Resize surface arrays in
  ! order to eliminate unused space.  At the same time, reorder surface grid W
  ! cell indices so that land cells are first, deep lake cells are second, and
  ! sea cells are third.

  nmsfc = mmsfc
  nvsfc = mvsfc
  nwsfc = mwsfc

  nasfc = 1
  nland = 1
  nlake = 1
  nsea  = 1

  allocate (iwnew(nwsfc))
  iwnew(1) = 1

  ! Loop over all surface cells and search for those that are land or shallow lake cells

  do iwsfc = 2,nwsfc
     if (sfcg%leaf_class(iwsfc) >= 2) then
        nasfc = nasfc + 1
        nland = nland + 1
        iwnew(iwsfc) = nasfc
     endif
  enddo

  ! Loop over all surface cells and search for those that are deep lake cells

  do iwsfc = 2,nwsfc
     if (sfcg%leaf_class(iwsfc) == 1) then
        nasfc = nasfc + 1
        nlake = nlake + 1
        iwnew(iwsfc) = nasfc
     endif
  enddo

  ! Loop over all surface cells and search for those that are sea cells

  do iwsfc = 2,nwsfc
     if (sfcg%leaf_class(iwsfc) == 0) then
        nasfc = nasfc + 1
        nsea = nsea + 1
        iwnew(iwsfc) = nasfc
     endif
  enddo

  onland = 0
  onlake = nland - 1
  onsea  = nland + nlake - 2

  print*, 'nmsfc ', nmsfc
  print*, 'nvsfc ', nvsfc
  print*, 'nwsfc ', nwsfc

  print*, 'nasfc ', nasfc
  print*, 'nland ', nland
  print*, 'nlake ', nlake
  print*, 'nsea  ', nsea

  print*, 'onland ', onland
  print*, 'onlake ', onlake
  print*, 'onsea  ', onsea

  if (nasfc /= nwsfc) stop 'stop: nasfc /= nwsfc'

  ! Loop through SFC grid M points in old order (same as new order)

  allocate (itab_msfc_temp(nmsfc))
  do imsfc = 2,nmsfc
     itab_msfc_temp(imsfc) = itab_msfc(imsfc)
     itab_msfc_temp(imsfc)%iwn(1) = iwnew(itab_msfc(imsfc)%iwn(1))
     itab_msfc_temp(imsfc)%iwn(2) = iwnew(itab_msfc(imsfc)%iwn(2))
     itab_msfc_temp(imsfc)%iwn(3) = iwnew(itab_msfc(imsfc)%iwn(3))
  enddo
  call move_alloc(itab_msfc_temp, itab_msfc)

  ! Loop through SFC grid V points in old order (same as new order)

  allocate (itab_vsfc_temp(nvsfc))
  do ivsfc = 2,nvsfc
     itab_vsfc_temp(ivsfc) = itab_vsfc(ivsfc)
     itab_vsfc_temp(ivsfc)%iwn(1) = iwnew(itab_vsfc(ivsfc)%iwn(1))
     itab_vsfc_temp(ivsfc)%iwn(2) = iwnew(itab_vsfc(ivsfc)%iwn(2))
  enddo
  call move_alloc(itab_vsfc_temp, itab_vsfc)

  ! Loop through SFC grid W points in old order

  allocate (itab_wsfc_temp(nwsfc))
  do iwsfc = 2,nwsfc
     inew = iwnew(iwsfc)

     itab_wsfc_temp(inew) = itab_wsfc(iwsfc)

     ! Convert IWN neighbor indices of inew W point if inew cell is Voronoi

     if (itab_wsfc(iwsfc)%ivoronoi == 3) then
        do j = 1,itab_wsfc(iwsfc)%npoly
           itab_wsfc_temp(inew)%iwn(j) = iwnew(itab_wsfc(iwsfc)%iwn(j))
        enddo
     endif

  enddo
  call move_alloc(itab_wsfc_temp, itab_wsfc)

  call resize_sfcgrid(nmsfc, nvsfc, nwsfc, iwnew=iwnew)

  deallocate(iwnew)

  allocate(tot_area(nwa))
  tot_area(:) = 0.0_r8

  ! Loop over all surface cells

  do iwsfc = 2,nwsfc

     ! Loop over ATM coupling areas of current surface cell

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw = itab_wsfc(iwsfc)%iwatm(j)
        tot_area(iw) = tot_area(iw) + itab_wsfc(iwsfc)%arc(j)
     enddo
  enddo

  ! Scale surface cell and ATM coupling areas slightly so that coupling areas sum
  ! almost exactly to arw0 within an ATM column
  
  do iwsfc = 2,nwsfc
     area_tot = 0.

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw  = itab_wsfc(iwsfc)%iwatm(j)

        itab_wsfc(iwsfc)%arc(j) = itab_wsfc(iwsfc)%arc(j) * arw0(iw) / tot_area(iw)

        area_tot = area_tot + itab_wsfc(iwsfc)%arc(j)
     enddo

     sfcg%area(iwsfc) = area_tot
  enddo

  deallocate(tot_area)

  ! Initialize soil grid depth and vertical spacing arrays

  ! Soil vertical grid spacing arrays

  call alloc_landcol()

  ! Iteration to find land grid vertical stretch rate

  write(io6,'(/,a,/)') 'Iterations for land grid vertical stretch ratio'

  ! If this is long-term groundwater spin-up simulation with long time steps,
  ! soil layers must (generally) be thicker than 0.5 m for stability.  At the
  ! same time, it is desirable for the thicker/deeper soil layers to be identical
  ! with those of the subsequent simulation that will be initialized with the
  ! spun-up groundwater, and for the thinner soil layers in the subsequent
  ! simulation to be groupable together into sets that correspond to individual
  ! layers in the present spin-up simulation.  To achieve this, set nzga 
  ! (an alternate value of nzg) to be the nzg value in the subsequent 
  ! simulation so that the soil depths of that simulation can be replicated
  ! here, and then group the thinner layers together here to be consistent
  ! with nzg (which is less than nzga) that is used for this spin-up simulation.

  ! In order to determine the correct value of nzg to specify in OLAMIN for the
  ! spin-up simulation, first look at the soil layer depths in the subsequent
  ! simulation and determine how the thin layers would be grouped into each
  ! layer of the spin-up simulation.

  if (nl%igw_spinup == 1) then
     nzga = 25
  else
     nzga = nzg
  endif

  srati = 0.5
  do iter = 1,20
     srati = 1. / (landgrid_depth * (1. - srati) / (srati * landgrid_dztop) + 1.)**(1./real(nzga))
     print*, 'iter, stretch ratio ',iter,1./srati
  enddo

  ! Compute soil grid levels

  dz = landgrid_dztop
  slz(nzg+1) = 0.
  thick = 0.

  k = nzg + 1

  do kk = nzga,1,-1

     thick = thick + dz

     if (nl%igw_spinup /= 1 .or. thick > 0.5) then
        k = k - 1 
        slz(k) = slz(k+1) - thick

        dslz  (k) = slz(k+1) - slz(k)
        dslzo2(k) = .5 * dslz(k)
        dslzi (k) = 1. / dslz(k)
        slzt  (k) = .5 * (slz(k) + slz(k+1))

        thick = 0.
     endif

     dz = dz / srati
  enddo

  call landgrid_print()

  ! Initialize soil static properties

  allocate (land%usdatext            (nland)) ; land%usdatext         = 0
  allocate (land%z_bedrock           (nland)) ; land%z_bedrock        = 0.
  allocate (land%gpp                 (nland)) ; land%gpp              = 0.
  allocate (land%glhymps_ksat        (nland)) ; land%glhymps_ksat     = 0.
  allocate (land%glhymps_ksat_pfr    (nland)) ; land%glhymps_ksat_pfr = 0.
  allocate (land%glhymps_poros       (nland)) ; land%glhymps_poros    = 0.
  allocate (land%sand            (nzg,nland)) ; land%sand             = 0.
  allocate (land%clay            (nzg,nland)) ; land%clay             = 0.
  allocate (land%silt            (nzg,nland)) ; land%silt             = 0.
  allocate (land%organ           (nzg,nland)) ; land%organ            = 0.
  allocate (land%bulkdens_drysoil(nzg,nland)) ; land%bulkdens_drysoil = 0.
  allocate (land%pH_soil         (nzg,nland)) ; land%pH_soil          = 0.
  allocate (land%cec_soil        (nzg,nland)) ; land%cec_soil         = 0.

  ! Always read FAO soil data and fill single-level usdatext array.  This may be needed
  ! permanently to fill in holes in the SoilGrids maps.

  allocate (qlat(nqa), qlon(nqa), topq(nqa))

  call land_database_read(nland, &
       sfcg%glatw(1:nland),      &
       sfcg%glonw(1:nland),      &
       soil_database,            &
       soil_database,            &
       'soil_text',              &
       idatq=sfcg%ioge(1:nland))

  ! Loop over all land cells (already defined and filled with leaf_class)

  do iland = 2,nland

     ! Convert from FAO soil type to USDA soil textural class

     call fao_usda(sfcg%ioge(iland), land%usdatext(iland))

     ! Customization of soil composition can be done inside subroutine
     ! usda_composition by modifying usdatext(iland) and/or how it is used

     call usda_composition(iland, land%usdatext(iland))

  enddo

  if (isoilflg == 1) then

     ! Read SoilGrids soil composition from SoilGrids data files

     call soilgrids_read()

     ! Read GLHYMPS soil permeability and porosity database files

     print*, 'calling glhymps_read '

     call glhymps_read()

  endif

print*, 'finished soilgrids_read '

! REMAINING TASKS:
! (1) Method for computing unit normal vectors for sfc hex cells over land
! (2) Check to make sure all "flux cell" area is accounted for


! RELATED TASKS:
! (1) 3D soil water
! (2) MPI parallelism for soil
! (3) Deep lake model
! (4) Shallow lake model
! (5) 

end subroutine makesfc3

!===============================================================================

subroutine topo_init(nqa,topq,glatq,glonq,xeq,yeq,zeq)

  ! Subroutine topo_init fills the TOPQ array with topography height.  TOPQ and
  ! its corresponding horizontal coordinates (GLATQ, GLONQ) or (XEQ, YEQ, ZEQ)
  ! are composite arrays that include both M and W horizontal stagger points
  ! in the OLAM hexagonal grid.

  use misc_coms,   only: io6, deltax
  use consts_coms, only: pi1, pio180, grav

  use oname_coms, only: nl

  use dcmip_initial_conditions_test_1_2_3, only: &
     test1_advection_orography, &
     test2_steady_state_mountain, &
     test2_schaer_mountain

  implicit none

  integer, intent(in) :: nqa
  real, intent(in) :: glatq(nqa),glonq(nqa),xeq(nqa),yeq(nqa),zeq(nqa)

  real, intent(out) :: topq(nqa)

  integer :: iq

  real :: hfwid
  real :: hgt
  real :: hfwid2
  real :: hoffset

  real :: r, r0

  !-------------------------------------------------------------------
  ! Variables for NCAR DCMIP 2012 TEST CASES

  real(8) :: zm0, rhom0, u0, v0, wm0
  real(8) :: lon,lat,p,t,phis,ps,q,q1,q2,q3,q4
  real(8) :: hyam, hybm, gc
  real(8) :: time0   = 0.0d0
  integer :: zcoords = 1.0d0
  integer :: cfv = 0
  integer :: shear = 0
  logical :: hybrid_eta = .false.
  !-------------------------------------------------------------------

  ! By default, the TOPQ array is filled with a value of 0.

  topq(:) = 0.

  ! The remainder of this subroutine is for re-defining topography, and a few
  ! commonly used options are given below.  However, if itopoflg is set to 1 in
  ! OLAMIN, topography values set in this subroutine will be overridden in
  ! subroutine topo_database, which interpolates topography to the OLAM grid
  ! from a standard topography dataset.

  ! Horizontal loop over all land topography points

  do iq = 2,nqa

     !=============================================================================
     ! WITCH OF AGNESI MOUNTAIN
     !=============================================================================

     ! Half-width and height commonly used for "10-meter-mountain" test

     ! hfwid = 10000.
     ! hgt = 10.0
     ! hoffset = 0.

     ! Half-width and height used for Dudhia simulations

     ! hfwid = 5. * deltax ! Value for Dudhia simulations
     ! hgt = 405.          ! Value for Dudhia simulations
     ! hgt = 1012.         ! Value for Dudhia simulations
     ! hoffset = 1.

     ! hfwid2 = hfwid**2

     ! The following form (with a nonzero hoffset value subtracted) gives zero
     ! topography at large distance from mountain top.

     ! topq(iq) = max(0.,hgt * hfwid2 / (hfwid2 + xeq(iq)**2) - hoffset) 

     ! write(io6,*) 'topq ',iq,xeq(iq),topq(iq)

     !=============================================================================
     ! SHALLOW WATER TEST CASE 5
     !=============================================================================

     if (nl%test_case == 5) then
        r0 = pi1 / 9.

        r = sqrt((glonq(iq) * pio180 + 0.5 * pi1) ** 2 &
          + (glatq(iq) * pio180 - pi1 / 6.) ** 2)

        topq(iq) = max(0., 2000. * (1. - r / r0))

        print*, 'topoinit ',iq,r,r0,topq(iq)
     endif

     !=============================================================================
     ! NCAR DCMIP 2012 TEST CASES
     !=============================================================================

     if (nl%test_case ==  13 .or. &
         nl%test_case == 200 .or. &
         nl%test_case == 201 .or. &
         nl%test_case ==  21 .or. &
         nl%test_case ==  22) then

        lon = pio180 * glonq(iq)
        lat = pio180 * glatq(iq)

        if (glonq(iq) < 0.) lon = pio180 * (glonq(iq) + 360.)

        p = 1.0d0
        zm0 = 0.0d0

        if (nl%test_case == 13) then

           !===================================================================
           ! DCMIP-2012 TEST CASE 13 - HORIZONTAL ADVECTION OF THIN
           ! CLOUD-LIKE TRACERS IN THE PRESENCE OF OROGRAPHY
           !===================================================================

           call test1_advection_orography(lon,lat,p,zm0,zcoords, &
              cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wm0, &
              t,phis,ps,rhom0,q,q1,q2,q3,q4)

           topq(iq) = phis / grav

        elseif (nl%test_case == 200 .or. &
                nl%test_case == 201) then

           !==================================================================
           ! DCMIP-2012 TEST CASES 200, 201 - STEADY STATE ATMOSPHERE AT REST
           ! IN THE PRESENCE OF OROGRAPHY
           !==================================================================

           call test2_steady_state_mountain(lon,lat,p,zm0,zcoords, &
              hybrid_eta,hyam,hybm,u0,v0,wm0, &
              t,phis,ps,rhom0,q)

           topq(iq) = phis / grav

        elseif (nl%test_case == 21 .or. &
                nl%test_case == 22) then

           !==================================================================
           ! DCMIP-2012 Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves
           ! over a Schaer-type Mountain
           !==================================================================

           call test2_schaer_mountain(lon,lat,p,zm0,zcoords, &
              hybrid_eta,hyam,hybm,shear,u0,v0,wm0, &
              t,phis,ps,rhom0,q)

           topq(iq) = phis / grav

        endif

        print*, 'topodefn ',iq,topq(iq)

     endif

  enddo

end subroutine topo_init

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

  !   1 Urban
  !   2 Low Sparse Grassland
  !   3 Coniferous Forest
  !   4 Deciduous Conifer Forest
  !   5 Deciduous Broadleaf Forest
  !   6 Evergreen Broadleaf Forests
  !   7 Tall Grasses and Shrubs
  !   8 Bare Desert
  !   9 Upland Tundra
  !  10 Irrigated Grassland
  !  11 Semi Desert
  !  12 Glacier Ice
  !  13 Wooded Wet Swamp
  !  14 Inland Water
  !  15 Sea Water
  !  16 Shrub Evergreen
  !  17 Shrub Deciduous
  !  18 Mixed Forest and Field
  !  19 Evergreen Forest and Fields
  !  20 Cool Rain Forest
  !  21 Conifer Boreal Forest
  !  22 Cool Conifer Forest
  !  23 Cool Mixed Forest
  !  24 Mixed Forest
  !  25 Cool Broadleaf Forest
  !  26 Deciduous Broadleaf Forest
  !  27 Conifer Forest
  !  28 Montane Tropical Forests
  !  29 Seasonal Tropical Forest
  !  30 Cool Crops and Towns
  !  31 Crops and Town
  !  32 Dry Tropical Woods
  !  33 Tropical Rainforest
  !  34 Tropical Degraded Forest
  !  35 Corn and Beans Cropland
  !  36 Rice Paddy and Field
  !  37 Hot Irrigated Cropland
  !  38 Cool Irrigated Cropland
  !  39 Cold Irrigated Cropland
  !  40 Cool Grasses and Shrubs
  !  41 Hot and Mild Grasses and Shrubs
  !  42 Cold Grassland
  !  43 Savanna (Woods)
  !  44 Mire, Bog, Fen
  !  45 Marsh Wetland
  !  46 Mediterranean Scrub
  !  47 Dry Woody Scrub
  !  48 Dry Evergreen Woods
  !  49 Volcanic Rock
  !  50 Sand Desert
  !  51 Semi Desert Shrubs
  !  52 Semi Desert Sage
  !  53 Barren Tundra
  !  54 Cool Southern Hemisphere Mixed Forests
  !  55 Cool Fields and Woods
  !  56 Forest and Field
  !  57 Cool Forest and Field
  !  58 Fields and Woody Savanna
  !  59 Succulent and Thorn Scrub
  !  60 Small Leaf Mixed Woods
  !  61 Deciduous and Mixed Boreal Forest
  !  62 Narrow Conifers
  !  63 Wooded Tundra
  !  64 Heath Scrub
  !  65 Coastal Wetland, NW
  !  66 Coastal Wetland, NE
  !  67 Coastal Wetland, SE
  !  68 Coastal Wetland, SW
  !  69 Polar and Alpine Desert
  !  70 Glacier Rock
  !  71 Salt Playas
  !  72 Mangrove
  !  73 Water and Island Fringe
  !  74 Land, Water, and Shore (see Note 1)
  !  75 Land and Water, Rivers (see Note 1)
  !  76 Crop and Water Mixtures
  !  77 Southern Hemisphere Conifers
  !  78 Southern Hemisphere Mixed Forest
  !  79 Wet Sclerophylic Forest
  !  80 Coastline Fringe
  !  81 Beaches and Dunes
  !  82 Sparse Dunes and Ridges
  !  83 Bare Coastal Dunes
  !  84 Residual Dunes and Beaches
  !  85 Compound Coastlines
  !  86 Rocky Cliffs and Slopes
  !  87 Sandy Grassland and Shrubs
  !  88 Bamboo
  !  89 Moist Eucalyptus
  !  90 Rain Green Tropical Forest
  !  91 Woody Savanna
  !  92 Broadleaf Crops
  !  93 Grass Crops
  !  94 Crops, Grass, Shrubs
  !  95 Evergreen Tree Crop
  !  96 Deciduous Tree Crop
  !  99 Interrupted Areas (Goodes Homolosine Projection)
  ! 100 Missing Data

  !------------------------------------------!
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
            21                             / !     ! 99 is Goode Homolosine
  !------------------------------------------!     !    empty space
  !          1  2  3  4  5  6  7  8  9 10          ! 100 is missing data
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

  iusda = catb(ifao)

end subroutine fao_usda

!==========================================================================

subroutine usda_composition(iland,usdatext)

  use mem_land, only: land, nzg, slzt

  implicit none

  integer, intent(in) :: iland
  integer, intent(in) :: usdatext

  integer :: k

  ! This subroutine fills soil composition parameters based on USDA soil textural class

  ! (This option is being replaced by utilization of SoilGrids soil composition data,
  ! but will be retained for the interim in order to allow intercomparison and testing.
  ! It may prove necessary to retain this subroutine permanently as a means to fill
  ! in missing data in the SoilGrids composition datasets.)

  ! van Genuchten parameters are not applied here; they are included in the following
  ! parameter statement only to keep it identical with subroutine soil_ptf.

  real, parameter :: soilparms4(10,12) = reshape( (/ &
  !------------------------------------------------------------------------------------------------------
  ! sand   clay   silt  organ  wresid_vg wsat_vg  ksat_vg  alpha_vg  en_vg  slwilt      USDA SOIL CLASS
  !                                                 (m/s)    (1/m)                       # AND NAME
  !------------------------------------------------------------------------------------------------------
    .92,   .03,   .05,   .00,   .045,    .43,    .825e-4,    14.5,   2.68,  .070,   & !  1 sand
    .82,   .06,   .12,   .01,   .057,    .41,    .405e-4,    12.4,   2.28,  .075,   & !  2 loamy sand
    .58,   .10,   .32,   .02,   .065,    .41,    .123e-4,     7.5,   1.89,  .114,   & !  3 sandy loam
    .17,   .13,   .70,   .03,   .067,    .45,    .125e-5,     2.0,   1.41,  .179,   & !  4 silt loam
    .43,   .18,   .39,   .05,   .078,    .43,    .289e-5,     3.6,   1.56,  .155,   & !  5 loam
    .58,   .27,   .15,   .04,   .100,    .39,    .364e-5,     5.9,   1.48,  .175,   & !  6 sandy clay loam
    .10,   .34,   .56,   .06,   .089,    .43,    .194e-6,     1.0,   1.23,  .218,   & !  7 silty clay loam
    .32,   .34,   .34,   .07,   .095,    .41,    .722e-6,     1.9,   1.31,  .250,   & !  8 clay loam
    .52,   .42,   .06,   .08,   .100,    .38,    .333e-6,     2.7,   1.23,  .219,   & !  9 sandy clay
    .06,   .47,   .47,   .09,   .070,    .36,    .556e-7,      .5,   1.09,  .283,   & ! 10 silty clay
    .22,   .58,   .20,   .10,   .068,    .38,    .556e-6,      .8,   1.09,  .286,   & ! 11 clay 
    .10,   .06,   .84,   .05,   .034,    .46,    .694e-6,     1.6,   1.37,  .200/), & ! 12 silt
     (/10,12/) )

  land%z_bedrock       (iland) = -100.
  land%gpp             (iland) = 0.
  land%glhymps_ksat    (iland) = 0.
  land%glhymps_ksat_pfr(iland) = 0.
  land%glhymps_poros   (iland) = 0.

  ! For this option, assign single-level FAO textural class to all soil layers.

  do k = 1,nzg

     if (slzt(k) > land%z_bedrock(iland)) then
        land%sand            (k,iland) = soilparms4(1,usdatext)
        land%clay            (k,iland) = soilparms4(2,usdatext)
        land%silt            (k,iland) = soilparms4(3,usdatext)
        land%organ           (k,iland) = soilparms4(4,usdatext)
        land%bulkdens_drysoil(k,iland) = 2700. * (1. - soilparms4(6,usdatext))
     else
        land%sand            (k,iland) = -0.001 ! Negative sand indicates bedrock
        land%clay            (k,iland) = 0.
        land%silt            (k,iland) = 0.
        land%organ           (k,iland) = 0.
        land%bulkdens_drysoil(k,iland) = 2700. ! [kg/m^3]
     endif

     land%pH_soil (k,iland) = 7.
     land%cec_soil(k,iland) = 0.

  enddo

end subroutine usda_composition

!==========================================================================

subroutine sfcgfile_write()

  ! Write sfcg quantities to sfcgfile

  use max_dims,   only: maxnlspoly, pathlen
  use mem_sfcg,   only: nmsfc, nvsfc, nwsfc, itab_msfc, itab_vsfc, itab_wsfc, sfcg, sfcgfile
  use mem_land,   only: nland, onland, land, nzg, &
                        slz, dslz, dslzo2, dslzi, slzt
  use mem_lake,   only: nlake, onlake
  use mem_sea,    only: nsea, onsea
  use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer :: iclobber1 = 1
  integer :: ndims, idims(2)
  integer :: im, iv, iw

  real,    allocatable :: rscr1(:)
  integer, allocatable :: iscr1(:)
  real,    allocatable :: rscr2(:,:)
  integer, allocatable :: iscr2(:,:)

  integer, pointer :: i1pointer(:)

  character(pathlen) :: flnm

  ! Open sfcgfile

  flnm = trim(sfcgfile)//'.h5'

  call shdf5_open(flnm,'W',iclobber1)

  ! Write the grid dimensions

  ndims = 1
  idims(1) = 1

  call shdf5_orec(ndims, idims, 'nmsfc'  , ivars=nmsfc)
  call shdf5_orec(ndims, idims, 'nvsfc'  , ivars=nvsfc)
  call shdf5_orec(ndims, idims, 'nwsfc'  , ivars=nwsfc)
  call shdf5_orec(ndims, idims, 'nland'  , ivars=nland)
  call shdf5_orec(ndims, idims, 'nlake'  , ivars=nlake)
  call shdf5_orec(ndims, idims, 'nsea'   , ivars=nsea)
  call shdf5_orec(ndims, idims, 'onland' , ivars=onland)
  call shdf5_orec(ndims, idims, 'onlake' , ivars=onlake)
  call shdf5_orec(ndims, idims, 'onsea'  , ivars=onsea)
  call shdf5_orec(ndims, idims, 'nzg'    , ivars=nzg)

  ! Write ITAB_M ARRAYS

  ndims    = 2
  idims(1) = 3
  idims(2) = nmsfc

  allocate (iscr2(3,nmsfc))
  do im = 1,nmsfc
     iscr2(1:3,im) = itab_msfc(im)%ivn(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_msfc%ivn',ivar2=iscr2)

  do im = 1,nmsfc
     iscr2(1:3,im) = itab_msfc(im)%iwn(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_msfc%iwn',ivar2=iscr2)
  deallocate(iscr2)

  ! Write ITAB_V ARRAYS

  ndims    = 2
  idims(1) = 2
  idims(2) = nvsfc

  allocate (iscr2(2,nvsfc))
  do iv = 1,nvsfc
     iscr2(1:2,iv) = itab_vsfc(iv)%imn(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%imn',ivar2=iscr2)

  do iv = 1,nvsfc
     iscr2(1:2,iv) = itab_vsfc(iv)%iwn(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%iwn',ivar2=iscr2)
  deallocate(iscr2)

  ! Write ITAB_WSFC SCALARS

  ndims    = 1
  idims(1) = nwsfc

  allocate (iscr1(nwsfc))
  do iw = 1,nwsfc
     iscr1(iw) = itab_wsfc(iw)%ivoronoi
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%ivoronoi' ,ivar1=iscr1)

  do iw = 1,nwsfc
     iscr1(iw) = itab_wsfc(iw)%npoly
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%npoly'    ,ivar1=iscr1)

  do iw = 1,nwsfc
     iscr1(iw) = itab_wsfc(iw)%nwatm
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%nwatm'    ,ivar1=iscr1)
  deallocate(iscr1)

  ! Write ITAB_WSFC ARRAYS

  ndims = 2
  idims(1) = 7
  idims(2) = nwsfc

  allocate (iscr2(7,nwsfc))
  do iw = 1,nwsfc
     iscr2(1:7,iw) = itab_wsfc(iw)%imn(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%imn',ivar2=iscr2)

  do iw = 1,nwsfc
     iscr2(1:7,iw) = itab_wsfc(iw)%ivn(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%ivn',ivar2=iscr2)

  do iw = 1,nwsfc
     iscr2(1:7,iw) = itab_wsfc(iw)%iwn(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%iwn',ivar2=iscr2)
  deallocate(iscr2)

  idims(1) = 8

  allocate (iscr2(8,nwsfc))
  do iw = 1,nwsfc
     iscr2(1:8,iw) = itab_wsfc(iw)%iwatm(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%iwatm',ivar2=iscr2)

  do iw = 1,nwsfc
     iscr2(1:8,iw) = itab_wsfc(iw)%kwatm(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%kwatm',ivar2=iscr2)
  deallocate(iscr2)

  allocate (rscr2(8,nwsfc))
  do iw = 1,nwsfc
     rscr2(1:8,iw) = itab_wsfc(iw)%arc(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%arc',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:8,iw) = itab_wsfc(iw)%arcoarsfc(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%arcoarsfc',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:8,iw) = itab_wsfc(iw)%arcoariw(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%arcoariw',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:8,iw) = itab_wsfc(iw)%arcoarkw(1:8)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%arcoarkw',rvar2=rscr2)
  deallocate(rscr2)

  ! Write grid structure variables

  ndims = 1
  idims(1) = nzg+1

  call shdf5_orec(ndims, idims, 'slz', rvar1=slz)

  idims(1) = nzg

  call shdf5_orec(ndims, idims, 'dslz'   , rvar1=dslz)
  call shdf5_orec(ndims, idims, 'dslzo2' , rvar1=dslzo2)
  call shdf5_orec(ndims, idims, 'dslzi'  , rvar1=dslzi)
  call shdf5_orec(ndims, idims, 'slzt'   , rvar1=slzt)

  ndims    = 1
  idims(1) = nmsfc

  call shdf5_orec(ndims, idims, 'xem'       , rvar1=sfcg%xem)
  call shdf5_orec(ndims, idims, 'yem'       , rvar1=sfcg%yem)
  call shdf5_orec(ndims, idims, 'zem'       , rvar1=sfcg%zem)
  call shdf5_orec(ndims, idims, 'glatm'     , rvar1=sfcg%glatm)
  call shdf5_orec(ndims, idims, 'glonm'     , rvar1=sfcg%glonm)
  call shdf5_orec(ndims, idims, 'topm'      , rvar1=sfcg%topm)

  idims(1) = nvsfc

  call shdf5_orec(ndims, idims, 'xev'       , rvar1=sfcg%xev)
  call shdf5_orec(ndims, idims, 'yev'       , rvar1=sfcg%yev)
  call shdf5_orec(ndims, idims, 'zev'       , rvar1=sfcg%zev)
  call shdf5_orec(ndims, idims, 'dnu'       , rvar1=sfcg%dnu)
  call shdf5_orec(ndims, idims, 'dniu'      , rvar1=sfcg%dniu)
  call shdf5_orec(ndims, idims, 'dnv'       , rvar1=sfcg%dnv)
  call shdf5_orec(ndims, idims, 'dniv'      , rvar1=sfcg%dniv)

  idims(1) = nwsfc

  call shdf5_orec(ndims, idims, 'area'      , rvar1=sfcg%area)
  call shdf5_orec(ndims, idims, 'xew'       , rvar1=sfcg%xew)
  call shdf5_orec(ndims, idims, 'yew'       , rvar1=sfcg%yew)
  call shdf5_orec(ndims, idims, 'zew'       , rvar1=sfcg%zew)
  call shdf5_orec(ndims, idims, 'glatw'     , rvar1=sfcg%glatw)
  call shdf5_orec(ndims, idims, 'glonw'     , rvar1=sfcg%glonw)
  call shdf5_orec(ndims, idims, 'topw'      , rvar1=sfcg%topw)
  call shdf5_orec(ndims, idims, 'wnx'       , rvar1=sfcg%wnx)
  call shdf5_orec(ndims, idims, 'wny'       , rvar1=sfcg%wny)
  call shdf5_orec(ndims, idims, 'wnz'       , rvar1=sfcg%wnz)
  call shdf5_orec(ndims, idims, 'dzt_bot'   , rvar1=sfcg%dzt_bot)

  ! Write permanent grid cell data

  call shdf5_orec(ndims, idims, 'leaf_class', ivar1=sfcg%leaf_class)
  call shdf5_orec(ndims, idims, 'oge'       , ivar1=sfcg%ioge)

  idims(1) = nland

  call shdf5_orec(ndims, idims, 'usdatext'        , ivar1=land%usdatext)
  call shdf5_orec(ndims, idims, 'z_bedrock'       , rvar1=land%z_bedrock)
  call shdf5_orec(ndims, idims, 'gpp'             , rvar1=land%gpp)
  call shdf5_orec(ndims, idims, 'glhymps_ksat'    , rvar1=land%glhymps_ksat)
  call shdf5_orec(ndims, idims, 'glhymps_ksat_pfr', rvar1=land%glhymps_ksat_pfr)
  call shdf5_orec(ndims, idims, 'glhymps_poros'   , rvar1=land%glhymps_poros)

  ndims    = 2
  idims(1) = nzg
  idims(2) = nland

  call shdf5_orec(ndims, idims, 'sand'            , rvar2=land%sand)
  call shdf5_orec(ndims, idims, 'clay'            , rvar2=land%clay)
  call shdf5_orec(ndims, idims, 'silt'            , rvar2=land%silt)
  call shdf5_orec(ndims, idims, 'organ'           , rvar2=land%organ)
  call shdf5_orec(ndims, idims, 'bulkdens_drysoil', rvar2=land%bulkdens_drysoil)
  call shdf5_orec(ndims, idims, 'pH_soil'         , rvar2=land%pH_soil)
  call shdf5_orec(ndims, idims, 'cec_soil'        , rvar2=land%cec_soil)

  call shdf5_close()

end subroutine sfcgfile_write

!==========================================================================

subroutine sfcgfile_read_pd()

  use max_dims,   only: pathlen
  use mem_sfcg,   only: nmsfc, mmsfc, nvsfc, mvsfc, nwsfc, mwsfc, &
                        itab_msfc_pd, itab_vsfc_pd, itab_wsfc_pd, sfcgfile
  use mem_land,   only: nland, onland, nzg
  use mem_lake,   only: nlake, onlake
  use mem_sea,    only: nsea, onsea
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  integer            :: imsfc, ivsfc, iwsfc
  character(pathlen) :: flnm
  logical            :: there
  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  flnm = trim(sfcgfile)//'.h5'

  write(io6,*) 'Checking sfcgfile ',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'sfcgfile was not found - stopping run'
     stop 'stop: no sfcgfile'
  endif

  call shdf5_open(flnm,'R')

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'nmsfc' , ivars=nmsfc)
  call shdf5_irec(ndims, idims, 'nvsfc' , ivars=nvsfc)
  call shdf5_irec(ndims, idims, 'nwsfc' , ivars=nwsfc)
  call shdf5_irec(ndims, idims, 'nland' , ivars=nland)
  call shdf5_irec(ndims, idims, 'nlake' , ivars=nlake)
  call shdf5_irec(ndims, idims, 'nsea'  , ivars=nsea)
  call shdf5_irec(ndims, idims, 'onland', ivars=onland)
  call shdf5_irec(ndims, idims, 'onlake', ivars=onlake)
  call shdf5_irec(ndims, idims, 'onsea' , ivars=onsea)
  call shdf5_irec(ndims, idims, 'nzg'   , ivars=nzg)

  write(io6, '(/,a)')    '==============================================='
  write(io6, '(a)')      'Reading from sfcgfile:'
  write(io6, '(2(a,i0))')'  nwsfc = ', nwsfc, ', nzg = ', nzg
  write(io6, '(a,/)')    '==============================================='

  ndims = 2
  idims(1) = 8
  idims(2) = nwsfc

  allocate (itab_msfc_pd(nmsfc))
  allocate (itab_vsfc_pd(nvsfc))
  allocate (itab_wsfc_pd(nwsfc))

  ! Read ITAB_MSFC ARRAYS

  ndims    = 2
  idims(1) = 3
  idims(2) = nmsfc

  allocate (iscr2(3,nmsfc))
  call shdf5_irec(ndims,idims,'itab_msfc%ivn',ivar2=iscr2)
  do imsfc = 1,nmsfc
     itab_msfc_pd(imsfc)%ivn(1:3) = iscr2(1:3,imsfc)
  enddo

  call shdf5_irec(ndims,idims,'itab_msfc%iwn',ivar2=iscr2)
  do imsfc = 1,nmsfc
     itab_msfc_pd(imsfc)%iwn(1:3) = iscr2(1:3,imsfc)
  enddo
  deallocate(iscr2)

  ! Read ITAB_VSFC ARRAYS

  ndims    = 2
  idims(1) = 2
  idims(2) = nvsfc

  allocate (iscr2(2,nvsfc))
  call shdf5_irec(ndims,idims,'itab_vsfc%imn',ivar2=iscr2)
  do ivsfc = 1,nvsfc
     itab_vsfc_pd(ivsfc)%imn(1:2) = iscr2(1:2,ivsfc)
  enddo

  call shdf5_irec(ndims,idims,'itab_vsfc%iwn',ivar2=iscr2)
  do ivsfc = 1,nvsfc
     itab_vsfc_pd(ivsfc)%iwn(1:2) = iscr2(1:2,ivsfc)
  enddo
  deallocate(iscr2)

  ! Read ITAB_WSFC SCALARS

  ndims    = 1
  idims(1) = nwsfc

  allocate (iscr1(nwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%ivoronoi'  ,ivar1=iscr1)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%ivoronoi = iscr1(iwsfc)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%nwatm'     ,ivar1=iscr1)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%nwatm = iscr1(iwsfc)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%npoly'     ,ivar1=iscr1)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%npoly = iscr1(iwsfc)
  enddo

  call shdf5_irec(ndims,idims,'leaf_class'          ,ivar1=iscr1)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%leaf_class = iscr1(iwsfc)
  enddo
  deallocate(iscr1)

  ! Read ITAB_WSFC ARRAYS

  ndims = 2
  idims(1) = 7
  idims(2) = nwsfc

  allocate (iscr2(7,nwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%imn',ivar2=iscr2)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%imn(1:7) = iscr2(1:7,iwsfc)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%ivn',ivar2=iscr2)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%ivn(1:7) = iscr2(1:7,iwsfc) 
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%iwn',ivar2=iscr2)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%iwn(1:7) = iscr2(1:7,iwsfc)
  enddo
  deallocate(iscr2)

  idims(1) = 8

  allocate (iscr2(8,nwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%iwatm',ivar2=iscr2)
  do iwsfc = 1,nwsfc
     itab_wsfc_pd(iwsfc)%iwatm(1:8) = iscr2(1:8,iwsfc)
  enddo

  call shdf5_close()

  mmsfc = nmsfc
  mvsfc = nvsfc
  mwsfc = nwsfc

end subroutine sfcgfile_read_pd

!==========================================================================

subroutine sfcgfile_read()

  use max_dims,   only: maxnlspoly, pathlen
  use mem_sfcg,   only: sfcg, itab_msfc, itab_vsfc, itab_wsfc, &
                        sfcgfile, mmsfc, mvsfc, mwsfc
  use mem_land,   only: land, itab_land, mland, nzg, slz, dslz, dslzo2, &
                        dslzi, slzt, alloc_landcol
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: im, iv, iw

  ! Scratch arrays for copying input

  real,    allocatable :: rscr2(:,:)
  integer, allocatable :: iscr2(:,:)

  ! Pointers to the global index of the local point

  integer :: lgmsfc(mmsfc)
  integer :: lgvsfc(mvsfc)
  integer :: lgwsfc(mwsfc)
  integer :: lgland(mland)

  lgmsfc = itab_msfc(1:mmsfc)%imglobe
  lgvsfc = itab_vsfc(1:mvsfc)%ivglobe
  lgwsfc = itab_wsfc(1:mwsfc)%iwglobe
  lgland = itab_land(1:mland)%iwglobe

  ! Check if sfcgfile exists

  flnm = trim(sfcgfile)//'.h5'

  write(io6,*) 'Checking sfcgfile ',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'sfcgfile was not found - stopping run'
     stop 'stop: no sfcgfile'
  endif

  call shdf5_open(flnm,'R')

 ! Read ITAB_WSFC ARRAYS

  ndims = 2
  idims(1) = 8
  idims(2) = mwsfc

  allocate (iscr2(8,mwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%kwatm',ivar2=iscr2, points=lgwsfc)

  do iw = 1,mwsfc
     itab_wsfc(iw)%kwatm(1:8) = iscr2(1:8,iw)
  enddo

  deallocate(iscr2)

  allocate (rscr2(8,mwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%arc',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arc(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoarsfc',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoarsfc(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoariw',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoariw(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoarkw',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoarkw(1:8) = rscr2(1:8,iw)
  enddo
  deallocate(rscr2)

  ! Read grid structure variables

  call alloc_landcol()

  ndims = 1
  idims(1) = nzg+1

  call shdf5_irec(ndims, idims, 'slz'    , rvar1=slz)

  idims(1) = nzg

  call shdf5_irec(ndims, idims, 'dslz'   , rvar1=dslz)
  call shdf5_irec(ndims, idims, 'dslzo2' , rvar1=dslzo2)
  call shdf5_irec(ndims, idims, 'dslzi'  , rvar1=dslzi)
  call shdf5_irec(ndims, idims, 'slzt'   , rvar1=slzt)

  ndims    = 1
  idims(1) = mmsfc

  call shdf5_irec(ndims, idims, 'xem'       , rvar1=sfcg%xem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'yem'       , rvar1=sfcg%yem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'zem'       , rvar1=sfcg%zem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'glatm'     , rvar1=sfcg%glatm, points=lgmsfc)
  call shdf5_irec(ndims, idims, 'glonm'     , rvar1=sfcg%glonm, points=lgmsfc)
  call shdf5_irec(ndims, idims, 'topm'      , rvar1=sfcg%topm,  points=lgmsfc)

  idims(1) = mvsfc

  call shdf5_irec(ndims, idims, 'xev'       , rvar1=sfcg%xev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'yev'       , rvar1=sfcg%yev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'zev'       , rvar1=sfcg%zev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dnu'       , rvar1=sfcg%dnu,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dniu'      , rvar1=sfcg%dniu, points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dnv'       , rvar1=sfcg%dnv,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dniv'      , rvar1=sfcg%dniv, points=lgvsfc)

  idims(1) = mwsfc

  call shdf5_irec(ndims, idims, 'area'      , rvar1=sfcg%area,    points=lgwsfc)
  call shdf5_irec(ndims, idims, 'xew'       , rvar1=sfcg%xew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'yew'       , rvar1=sfcg%yew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'zew'       , rvar1=sfcg%zew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'glatw'     , rvar1=sfcg%glatw,   points=lgwsfc)
  call shdf5_irec(ndims, idims, 'glonw'     , rvar1=sfcg%glonw,   points=lgwsfc)
  call shdf5_irec(ndims, idims, 'topw'      , rvar1=sfcg%topw,    points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wnx'       , rvar1=sfcg%wnx,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wny'       , rvar1=sfcg%wny,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wnz'       , rvar1=sfcg%wnz,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'dzt_bot'   , rvar1=sfcg%dzt_bot, points=lgwsfc)

  ! Read permanent grid cell data

  call shdf5_irec(ndims, idims, 'leaf_class', ivar1=sfcg%leaf_class, points=lgwsfc)
  call shdf5_irec(ndims, idims, 'oge'       , ivar1=sfcg%ioge,       points=lgwsfc)

  allocate (land%usdatext            (mland)) ; land%usdatext         = 0
  allocate (land%z_bedrock           (mland)) ; land%z_bedrock        = 0.
  allocate (land%gpp                 (mland)) ; land%gpp              = 0.
  allocate (land%glhymps_ksat        (mland)) ; land%glhymps_ksat     = 0.
  allocate (land%glhymps_ksat_pfr    (mland)) ; land%glhymps_ksat_pfr = 0.
  allocate (land%glhymps_poros       (mland)) ; land%glhymps_poros    = 0.
  allocate (land%sand            (nzg,mland)) ; land%sand             = 0.
  allocate (land%clay            (nzg,mland)) ; land%clay             = 0.
  allocate (land%silt            (nzg,mland)) ; land%silt             = 0.
  allocate (land%organ           (nzg,mland)) ; land%organ            = 0.
  allocate (land%bulkdens_drysoil(nzg,mland)) ; land%bulkdens_drysoil = 0.
  allocate (land%pH_soil         (nzg,mland)) ; land%pH_soil          = 0.
  allocate (land%cec_soil        (nzg,mland)) ; land%cec_soil         = 0.

  idims(1) = mland

  call shdf5_irec(ndims, idims, 'usdatext'        , ivar1=land%usdatext,        points=lgland)
  call shdf5_irec(ndims, idims, 'z_bedrock'       , rvar1=land%z_bedrock,       points=lgland)
  call shdf5_irec(ndims, idims, 'gpp'             , rvar1=land%gpp,             points=lgland)
  call shdf5_irec(ndims, idims, 'glhymps_ksat'    , rvar1=land%glhymps_ksat,    points=lgland)
  call shdf5_irec(ndims, idims, 'glhymps_ksat_pfr', rvar1=land%glhymps_ksat_pfr,points=lgland)
  call shdf5_irec(ndims, idims, 'glhymps_poros'   , rvar1=land%glhymps_poros,   points=lgland)

  ndims    = 2
  idims(1) = nzg
  idims(2) = mland

  call shdf5_irec(ndims, idims, 'sand'            , rvar2=land%sand,             points=lgland)
  call shdf5_irec(ndims, idims, 'clay'            , rvar2=land%clay,             points=lgland)
  call shdf5_irec(ndims, idims, 'silt'            , rvar2=land%silt,             points=lgland)
  call shdf5_irec(ndims, idims, 'organ'           , rvar2=land%organ,            points=lgland)
  call shdf5_irec(ndims, idims, 'bulkdens_drysoil', rvar2=land%bulkdens_drysoil, points=lgland)
  call shdf5_irec(ndims, idims, 'pH_soil'         , rvar2=land%pH_soil,          points=lgland)
  call shdf5_irec(ndims, idims, 'cec_soil'        , rvar2=land%cec_soil,         points=lgland)

  call shdf5_close()

end subroutine sfcgfile_read

!===============================================================================

subroutine landgrid_print()

  use mem_land,  only: nland, onland, itab_land, land, nzg, landgrid_dztop, &
                       landgrid_depth, slz, dslz, dslzo2, dslzi, slzt, &
                       alloc_landcol
  use misc_coms, only: io6

  implicit none

  integer :: k

  ! Print land grid vertical geometry 

  write(io6,'(/,a)') '=========================================='
  write(io6,'(a)'  ) '        LAND GRID VERTICAL GEOMETRY'
  write(io6,'(a)'  ) '=========================================='
  write(io6,'(a)'  ) '   k     slz(m)   k    dslz(m)   slzt(m)  '
  write(io6,'(a,/)') '=========================================='

  do k = nzg+1,1,-1
     if (k == nzg+1 .or. k == 1) then
        write(io6,11) k, slz(k)
     else
        write(io6,12) k, slz(k)
     endif
     if (k > 1) write(io6,13) k-1, dslz(k-1), slzt(k-1)
  enddo

  write(io6,'(//)')

  11 format (i4,f10.3,1x,3('======'))
  12 format (i4,f10.3,1x,3('------'))
  13 format (15x,i4,3f10.3)

end subroutine landgrid_print

