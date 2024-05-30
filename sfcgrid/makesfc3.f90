subroutine makesfc3()

! This subroutine generates the SURFACE GRID FILE for OLAM runs.

  use mem_grid,    only: zm, nwa, nva, nma, xew, yew, zew

  use misc_coms,   only: io6, itopoflg, topo_database, bathym_database, &
                         ibathflg, runtype, mdomain, topodb_cutoff

  use leaf_coms,   only: nvgcon, ivegflg, isoilflg, soil_database, veg_database

  use consts_coms, only: eradi, piu180, r8

  use land_db,     only: land_database_read

  use soilgrids_db,only: soilgrids_read

  use max_dims,    only: maxnlspoly

  use mem_ijtabs,  only: itab_m, itab_v, itab_w

  use mem_sfcg,    only: nsfcgrids, nmsfc, nvsfc, nwsfc, sfcg, &
                         itab_msfc, itab_vsfc, itab_wsfc, &
                         itab_msfc_vars, itab_vsfc_vars, itab_wsfc_vars, &
                         nswmzons, nswmzonll, swmzonrad, swmzonlat, swmzonlon, &
                         alloc_sfcgrid1

  use mem_land,    only: nland, onland, land, nzg, &
                         slz, dslz, dslzo2, dslzi, slzt, &
                         landgrid_dztop, landgrid_depth, alloc_landcol

  use mem_lake,    only: nlake, onlake

  use mem_sea,     only: nsea, onsea, sea, &
                         npomzons, npomzonll, pomzonrad, pomzonlat, pomzonlon

  use pom2k1d,     only: nzpom, pom, yy, alloc_pomgrid, pom_levels

  use oname_coms,  only: nl

  use mem_sfcnud,  only: nzg_nl, nzg_sp, kspm

  use sea_swm,     only: depthmax_swe

  implicit none

  integer :: k, iw, j, iv, im, np, kpom
  integer :: iland, inew, iswmzon, minside, ipomzon, isea
  integer :: imsfc, ivsfc, iwsfc
  integer :: nasfc

  integer, parameter :: nsgdata = 1 ! number of SoilGrids datasets to read

  real,    allocatable :: rscr(:)
  integer, allocatable :: iscr(:)

  type(itab_wsfc_vars), allocatable :: itab_wsfc_temp(:)

  real, allocatable :: slz_temp   (:)
  real, allocatable :: dslz_temp  (:)
  real, allocatable :: dslzo2_temp(:)
  real, allocatable :: dslzi_temp (:)
  real, allocatable :: slzt_temp  (:)

  integer, allocatable :: iwnew (:)

  integer :: iter
  real :: dz, srati, area_cutoff

  integer :: kk, koff
  real :: thick

  ! At this point in the process of surface grid initialization, the surface
  ! grid has been copied in Delaunay form (triangle mesh) from the atmospheric
  ! grid (by a call to subroutine copy_sfc_tri_grid).  This copy may, but need
  ! not, include all local mesh refinements applied to the atmosphere grid.

  ! If any additional local mesh refinements are to be made on the surface grid
  ! independent of those on the atmospheric grid, perform them next.

  if (nsfcgrids > 0) call spawn_nest(.false.)

  if (runtype == 'MAKEGRID_PLOT') return

  ! Generate Voronoi form of surface grid and compute its geometric properties.

  if (mdomain /= 4) then

     call voronoi_sfc()

  else

     ! Special for cartesian hex grid

     nmsfc = nma
     nvsfc = nva
     nwsfc = nwa

     call alloc_sfcgrid1(nmsfc, nvsfc, nwsfc)

     sfcg%xew = xew
     sfcg%yew = yew
     sfcg%zew = zew

     do iv = 2,nvsfc
        itab_vsfc(iv)%imn(1:2)  = itab_v(iv)%im(1:2)
        itab_vsfc(iv)%iwn(1:2)  = itab_v(iv)%iw(1:2)
     enddo

     do im = 2,nmsfc
        itab_msfc(im)%ivn(1:3) = itab_m(im)%iv(1:3)
        itab_msfc(im)%iwn(1:3) = itab_m(im)%iw(1:3)
     enddo

     do iw = 2,nwsfc
        np = itab_w(iw)%npoly

        itab_wsfc(iw)%npoly      = np
        itab_wsfc(iw)%imn (1:np) = itab_w(iw)%im  (1:np)
        itab_wsfc(iw)%ivn (1:np) = itab_w(iw)%iv  (1:np)
        itab_wsfc(iw)%iwn (1:np) = itab_w(iw)%iw  (1:np)
        itab_wsfc(iw)%dirv(1:np) = itab_w(iw)%dirv(1:np)
     enddo

  endif

  call grid_geometry_hex_sfc()

  ! Interpolate or assign topography height to W points of surface grid

  if (itopoflg == 1) then  ! from database

     call land_database_read(nwsfc, sfcg%glatw, sfcg%glonw, &
          topo_database(1), topo_database(1), 'topo', datq=sfcg%topw)

     area_cutoff = topodb_cutoff ** 2

     if (any(sfcg%area(2:) < area_cutoff)) then
        call land_database_read(nwsfc, sfcg%glatw, sfcg%glonw, &
             topo_database(2), topo_database(2), 'topo2', &
             datq=sfcg%topw, area=sfcg%area)
     endif

  else

     call topo_init(nwsfc, sfcg%topw, sfcg%glatw, sfcg%glonw, &
                    sfcg%xew, sfcg%yew,sfcg%zew)

  endif

  ! Read or assign surface type to W points of surface grid

  if (ivegflg == 1) then  ! from database

     call land_database_read(nwsfc, sfcg%glatw, sfcg%glonw, &
          veg_database, veg_database, 'leaf_class', idatq=sfcg%ioge)

     do iwsfc = 2, nwsfc
        call oge_leafclass(sfcg%ioge(iwsfc), sfcg%leaf_class(iwsfc))
     enddo

  else  ! from user-defined values in namelist file

     sfcg%leaf_class(2:nwsfc) = nvgcon

  endif

  ! Loop over all W points of surface grid

  do iwsfc = 2, nwsfc

     ! Prevent TOPQ lower than lowest model level zm(1)

     sfcg%topw(iwsfc) = max(sfcg%topw(iwsfc),zm(1))

     ! If point is ocean but has nonzero topography height,
     ! reset topography to zero.

     if (sfcg%leaf_class(iwsfc) == 0) sfcg%topw(iwsfc) = 0.

  enddo

  ! Now that land use has been assigned, Reorder surface grid W cell indices so
  ! that land cells are first, deep lake cells are second, and sea cells are third.

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

  do imsfc = 2, nmsfc
     itab_msfc(imsfc)%iwn(1) = iwnew( itab_msfc(imsfc)%iwn(1) )
     itab_msfc(imsfc)%iwn(2) = iwnew( itab_msfc(imsfc)%iwn(2) )
     itab_msfc(imsfc)%iwn(3) = iwnew( itab_msfc(imsfc)%iwn(3) )
  enddo

  ! Loop through SFC grid V points in old order (same as new order)

  do ivsfc = 2,nvsfc
     itab_vsfc(ivsfc)%iwn(1) = iwnew( itab_vsfc(ivsfc)%iwn(1) )
     itab_vsfc(ivsfc)%iwn(2) = iwnew( itab_vsfc(ivsfc)%iwn(2) )
  enddo

  ! Store previous W itabs in a temporary array

  call move_alloc(itab_wsfc, itab_wsfc_temp)

  allocate( itab_wsfc(nwsfc) )

  ! Loop through SFC grid W points in old order

  do iwsfc = 2, nwsfc
     inew = iwnew(iwsfc)

     itab_wsfc(inew) = itab_wsfc_temp(iwsfc)

     ! Convert IWN neighbor indices of inew W points
     do j = 1, itab_wsfc_temp(iwsfc)%npoly
        itab_wsfc(inew)%iwn(j) = iwnew( itab_wsfc_temp(iwsfc)%iwn(j) )
     enddo
  enddo

  deallocate(itab_wsfc_temp)

  ! Reorder the sfcg W arrays that have been previously computed

  allocate(rscr(nwsfc))

  rscr = sfcg%area
  sfcg%area(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%glatw
  sfcg%glatw(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%glonw
  sfcg%glonw(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%xew
  sfcg%xew(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%yew
  sfcg%yew(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%zew
  sfcg%zew(iwnew(1:nwsfc)) = rscr

  rscr = sfcg%topw
  sfcg%topw(iwnew(1:nwsfc)) = rscr

  deallocate(rscr)
  allocate(iscr(nwsfc))

  iscr = sfcg%leaf_class
  sfcg%leaf_class(iwnew(1:nwsfc)) = iscr

  iscr = sfcg%ioge
  sfcg%ioge(iwnew(1:nwsfc)) = iscr

  deallocate(iscr)
  deallocate(iwnew)

  do iwsfc = 2, nwsfc
     sfcg%wnx(iwsfc) = sfcg%xew(iwsfc) * eradi
     sfcg%wny(iwsfc) = sfcg%yew(iwsfc) * eradi
     sfcg%wnz(iwsfc) = sfcg%zew(iwsfc) * eradi
  enddo

  ! Allocate and initialize soil grid depth and vertical spacing arrays

  call alloc_landcol()

  allocate (kspm(nzg))

  ! Iteration to find land grid vertical stretch rate

  write(io6,'(/,a,/)') 'Iterations for land grid vertical stretch ratio'

  srati = 0.5
  do iter = 1,20
     srati = 1. / (landgrid_depth * (1. - srati) / (srati * landgrid_dztop) + 1.)**(1./real(nzg))
     print*, 'iter, stretch ratio ',iter,1./srati
  enddo

  ! If this is a long-term groundwater spin-up simulation with long time steps,
  ! soil layers must (generally) be thicker than 0.5 m for stability.  At the
  ! same time, it is desirable for the thicker/deeper soil layers to be identical
  ! with those of the subsequent simulation that will be initialized with the
  ! spun-up groundwater, and for the thinner soil layers in the subsequent
  ! simulation to be groupable together into sets that correspond to individual
  ! layers in the present spin-up simulation.  To achieve this, the following
  ! code computes soil layer thicknesses and depths the same as in the subsequent
  ! simulation, but then combines the thinner layers into groups that are each
  ! at least 0.5 m thick.  THEREFORE, THE USER MUST LEAVE NZG IN OLAMIN SET TO
  ! THE SAME VALUE THAT WILL BE USED IN THE SUBSEQUENT SIMULATION.  The code
  ! sets nzg to a new reduced value that reflects the grouping of thinner
  ! layers, but saves the original OLAMIN value in nzg_nl.  It also makes a
  ! mapping table that relates the vertical indexes of the soil layers between
  ! the grouped and ungrouped sets, and writes that table, along with nzg and
  ! nzg_nl to the output file that contains the spun-up groundwater
  ! and temperature values.  In the subsequent simulation that is initialized
  ! from the spun-up soil, the table is used to map the spun-up soil layers to
  ! the thinner and more numerous layers used in the subsequent simulation.

  ! Save nzg value read from OLAMIN in nzg_nl and, provisionally, in nzg_sp

  nzg_nl = nzg
  nzg_sp = nzg

  ! Compute soil grid levels

  dz = landgrid_dztop
  slz(nzg+1) = 0.
  thick = 0.

  k = nzg        ! k counts over grouped layers (if igw_spinup = 1)

  do kk = nzg,1,-1   ! kk counts over original ungrouped layers

     ! Map kk soil layers into k soil layers

     kspm(kk) = k

     ! print*, 'kspm ',k,kk,kspm(kk)

     thick = thick + dz

     if (nl%igw_spinup /= 1 .or. thick > 0.5) then
        slz(k) = slz(k+1) - thick

        dslz  (k) = slz(k+1) - slz(k)
        dslzo2(k) = .5 * dslz(k)
        dslzi (k) = 1. / dslz(k)
        slzt  (k) = .5 * (slz(k) + slz(k+1))

        thick = 0.
        k = k - 1
     endif

     dz = dz / srati

  enddo

  ! If this is a long-term groundwater spin-up simulation and soil layers have
  ! been grouped into thicker layers, the k index for the layer thicknesses and
  ! depth assigned above did not terminate at 1.  Therefore, shift the indexes
  ! here, and reduce nzg accordingly.

  koff = k
  if (koff > 0) then

     ! Computed reduced nzg value for grouped soil layers

     nzg_sp = nzg - koff
     nzg = nzg_sp

     kspm(1:nzg_nl) = kspm(1:nzg_nl) - koff

     print*, 'final nzg, nzg_sp, nzg_nl ',nzg, nzg_sp, nzg_nl
     print*, ' '
     print*, 'final kspm(1:nzg_nl) ',kspm(1:nzg_nl)

     ! Move soil vertical grid spacing arrays to temporary space

     call move_alloc (slz   , slz_temp)
     call move_alloc (dslz  , dslz_temp)
     call move_alloc (dslzo2, dslzo2_temp)
     call move_alloc (dslzi , dslzi_temp)
     call move_alloc (slzt  , slzt_temp)

     ! Reallocate soil vertical grid spacing arrays at reduced nzg size

     call alloc_landcol()

     ! Fill soil vertical grid spacing arrays with index-shifted values

     slz(nzg+1) = 0.

     do k = 1,nzg
           slz(k) =    slz_temp(k + koff)
          dslz(k) =   dslz_temp(k + koff)
        dslzo2(k) = dslzo2_temp(k + koff)
         dslzi(k) =  dslzi_temp(k + koff)
          slzt(k) =   slzt_temp(k + koff)
     enddo

     deallocate (slz_temp, dslz_temp, dslzo2_temp, dslzi_temp, slzt_temp)

  endif

  call landgrid_print()

  ! Fill bathymetry data

  if (ibathflg == 1) then  ! from database

     call land_database_read( nlake+nsea-1, sfcg%glatw(nland:), sfcg%glonw(nland:), &
          bathym_database, bathym_database, 'etopo1', datq=sfcg%bathym(nland:) )

     ! For land cells, set bathym equal to topw
     sfcg%bathym(2:nland) = sfcg%topw(2:nland)

     ! Prevent bathym from exceeding (topw - 1.0) for sea and lake cells
     sfcg%bathym(nland+1:) = min(sfcg%bathym(nland+1:), sfcg%topw(nland+1:) - 1.0)

  else

     sfcg%bathym(2       :nland  ) = sfcg%topw(2       :nland  )
     sfcg%bathym(2+onlake:1+onsea) = sfcg%topw(2+onlake:1+onsea) -  10.0
     sfcg%bathym(2+onsea :nwsfc  ) = sfcg%topw(2+onsea :nwsfc  ) - 100.0

  endif

  ! Set logical flag for IWSFC cells that use Shallow Water Model (SWM).
  ! Require that bathym depth be no greater than depthmax_swe

  do iswmzon = 1,nswmzons
     do iwsfc = 2,nwsfc
        call ngr_area(iswmzon,minside,sfcg%xew(iwsfc),sfcg%yew(iwsfc),sfcg%zew(iwsfc), &
                      nswmzonll, swmzonrad, swmzonlat, swmzonlon)

        if (minside == 1 .and. sfcg%bathym(iwsfc) > depthmax_swe) sfcg%swm_active(iwsfc) = .true.
     enddo
  enddo

  ! Set logical flag for SEA cells that use POM1D.
  ! Require that water be deeper than 30 m.

  allocate(sea%pom_active(nsea)); sea%pom_active = .false.

  do ipomzon = 1,npomzons
     do isea = 2,nsea
        iwsfc = isea + onsea

        call ngr_area(ipomzon,minside,sfcg%xew(iwsfc),sfcg%yew(iwsfc), &
                      sfcg%zew(iwsfc), npomzonll, pomzonrad, pomzonlat, pomzonlon)

        if ( minside == 1 .and. (.not. sfcg%swm_active(iwsfc)) .and. &
             sfcg%bathym(iwsfc) < -30. ) sea%pom_active(isea) = .true.

     enddo
  enddo

  ! Allocate POM1D grid arrays and define vertical levels

  call alloc_pomgrid(nsea)

  call pom_levels()

  ! Fill POM1D horizontally-dependent number of layers.

  do isea = 2,nsea
     iwsfc = isea + onsea

     if (sea%pom_active(isea)) then
        do kpom = 1,nzpom
           if (sfcg%bathym(iwsfc) < yy(kpom)) then
              pom%kba(isea) = kpom
           else
              exit
           endif
        enddo
     endif
  enddo

  ! Initialize soil static properties

  allocate (land%usdatext            (nland)) ; land%usdatext         = 0
  allocate (land%z_bedrock           (nland)) ; land%z_bedrock        = 0.
  allocate (land%gpp                 (nland)) ; land%gpp              = 0.
  allocate (land%glhymps_ksat        (nland)) ; land%glhymps_ksat     = 0.
! allocate (land%glhymps_ksat_pfr    (nland)) ; land%glhymps_ksat_pfr = 0.
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

! allocate (qlat(nqa), qlon(nqa), topq(nqa))

  if (isoilflg == 3) then

     do iland = 2, nland
        land%usdatext(iland) = nl%isoiltext

        ! Customization of soil composition can be done inside subroutine
        ! usda_composition by modifying usdatext(iland) and/or how it is used
        call usda_composition(iland, land%usdatext(iland))
     enddo

  else

     allocate(iscr(nland))

     call land_database_read(nland, &
          sfcg%glatw(1:nland),      &
          sfcg%glonw(1:nland),      &
          soil_database,            &
          soil_database,            &
          'soil_text',              &
          idatq=iscr)

     ! Loop over all land cells (already defined and filled with leaf_class)

     do iland = 2,nland

        ! Convert from FAO soil type to USDA soil textural class

        call fao_usda(iscr(iland), land%usdatext(iland))

     enddo

     do iland = 2,nland

        ! Customization of soil composition can be done inside subroutine
        ! usda_composition by modifying usdatext(iland) and/or how it is used

        call usda_composition(iland, land%usdatext(iland))

     enddo

     deallocate(iscr)

  endif

  if (isoilflg == 1) then

     ! Read SoilGrids soil composition from SoilGrids data files

     call soilgrids_read()

     ! Read GLHYMPS soil permeability and porosity database files

     print*, 'calling glhymps_read '

     call glhymps_read()

  endif

  print*, 'finished soilgrids_read '

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

  integer, intent( in) :: nqa
  real,    intent( in) :: glatq(nqa),glonq(nqa),xeq(nqa),yeq(nqa),zeq(nqa)
  real,    intent(out) :: topq(nqa)

  integer :: iq

! WITCH OF AGNESI MOUNTAIN
! real :: hfwid
! real :: hgt
! real :: hfwid2
! real :: hoffset

  real :: r, r0

  !-------------------------------------------------------------------
  ! Variables for NCAR DCMIP 2012 TEST CASES

  real(8) :: zm0, rhom0, u0, v0, wm0
  real(8) :: lon,lat,p,t,phis,ps,q,q1,q2,q3,q4
  real(8) :: hyam, hybm, gc
  integer :: zcoords = 1
  integer :: cfv = 0
  integer :: shear = 0
  logical :: hybrid_eta = .false.
  !-------------------------------------------------------------------

  ! By default, the TOPQ array is filled with a value of 0.

  topq(:) = 0.

  ! The remainder of this subroutine is defining topography in the special case
  ! where itopflg = 2.  A few occasionally used options are given below.

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

        ! print*, 'topoinit ',iq,r,r0,topq(iq)
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

        ! print*, 'topodefn ',iq,topq(iq)

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

  use mem_land,   only: land, nzg, slzt
  use oname_coms, only: nl

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

  land%z_bedrock       (iland) = nl%zbedrock
  land%gpp             (iland) = 0.
  land%glhymps_ksat    (iland) = nl%gnd_ksat
! land%glhymps_ksat_pfr(iland) = 0.
  land%glhymps_poros   (iland) = nl%gnd_poros

  ! For this option, assign single-level FAO textural class to all soil layers.

  do k = 1, nzg

     land%sand            (k,iland) = soilparms4(1,usdatext)
     land%clay            (k,iland) = soilparms4(2,usdatext)
     land%silt            (k,iland) = soilparms4(3,usdatext)
     land%organ           (k,iland) = soilparms4(4,usdatext)
     land%bulkdens_drysoil(k,iland) = 2700. * (1. - soilparms4(6,usdatext))
     land%pH_soil         (k,iland) = 7.
     land%cec_soil        (k,iland) = 30.

  enddo

end subroutine usda_composition

!==========================================================================

subroutine sfcgfile_write()

  ! Write sfcg quantities to sfcgfile

  use max_dims,   only: maxnlspoly, pathlen
  use mem_sfcg,   only: nmsfc, nvsfc, nwsfc, itab_msfc, itab_vsfc, itab_wsfc, &
                        sfcg, sfcgfile
  use mem_land,   only: nland, onland, land, nzg, &
                        slz, dslz, dslzo2, dslzi, slzt
  use mem_lake,   only: nlake, onlake
  use mem_sea,    only: nsea, onsea, sea
  use pom2k1d,    only: nzpom, pom, y, yy, dy, dyy
  use mem_sfcnud, only: nzg_nl, nzg_sp, kspm
  use oname_coms, only: nl

  use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close

  implicit none

  integer :: iclobber1 = 1
  integer :: ndims, idims(2)
  integer :: im, iv, iw

  integer, allocatable :: iscr1(:)
  real,    allocatable :: rscr2(:,:)
  integer, allocatable :: iscr2(:,:)

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
  call shdf5_orec(ndims, idims, 'nzpom'  , ivars=nzpom)

  ! Write nzg_nl, nzg_sp, and kspm to the sfcgrid file only if this MAKEGRID
  ! run is for a groundwater spin-up simulation.

  if (nl%igw_spinup == 1) then
     call shdf5_orec(ndims, idims, 'nzg_nl' , ivars=nzg_nl)
     call shdf5_orec(ndims, idims, 'nzg_sp' , ivars=nzg_sp)

     idims(1) = nzg_nl

     call shdf5_orec(ndims, idims, 'kspm'   , ivar1=kspm)
  endif

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

  do im = 1,nmsfc
     iscr2(1:3,im) = itab_msfc(im)%imn(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_msfc%imn',ivar2=iscr2)
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

  allocate (rscr2(2,nvsfc))
  do iv = 1,nvsfc
     rscr2(1:2,iv) = itab_vsfc(iv)%cosv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%cosv',rvar2=rscr2)

  do iv = 1,nvsfc
     rscr2(1:2,iv) = itab_vsfc(iv)%sinv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%sinv',rvar2=rscr2)

  do iv = 1,nvsfc
     rscr2(1:2,iv) = itab_vsfc(iv)%dxps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%dxps',rvar2=rscr2)

  do iv = 1,nvsfc
     rscr2(1:2,iv) = itab_vsfc(iv)%dyps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_vsfc%dyps',rvar2=rscr2)
  deallocate(rscr2)

  ! Write ITAB_WSFC SCALARS

  ndims    = 1
  idims(1) = nwsfc

  allocate (iscr1(nwsfc))
  do iw = 1,nwsfc
     iscr1(iw) = itab_wsfc(iw)%npoly
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%npoly'    ,ivar1=iscr1)

!TR  do iw = 1,nwsfc
!TR     iscr1(iw) = itab_wsfc(iw)%nwatm
!TR  enddo
!TR  call shdf5_orec(ndims,idims,'itab_wsfc%nwatm'    ,ivar1=iscr1)
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

  allocate (rscr2(7,nwsfc))
  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%dirv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%dirv',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%farm(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%farm',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%farv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%farv',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%gxps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%gxps1',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%gyps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%gyps1',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%gxps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%gxps2',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%gyps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%gyps2',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%ecvec_vx(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%ecvec_vx',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%ecvec_vy(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%ecvec_vy',rvar2=rscr2)

  do iw = 1,nwsfc
     rscr2(1:7,iw) = itab_wsfc(iw)%ecvec_vz(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_wsfc%ecvec_vz',rvar2=rscr2)
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

  idims(1) = nzpom

  call shdf5_orec(ndims, idims, 'y'   , rvar1=y)
  call shdf5_orec(ndims, idims, 'yy'  , rvar1=yy)
  call shdf5_orec(ndims, idims, 'dy'  , rvar1=dy)
  call shdf5_orec(ndims, idims, 'dyy' , rvar1=dyy)

  ndims    = 1
  idims(1) = nmsfc

  call shdf5_orec(ndims, idims, 'xem'       , rvar1=sfcg%xem)
  call shdf5_orec(ndims, idims, 'yem'       , rvar1=sfcg%yem)
  call shdf5_orec(ndims, idims, 'zem'       , rvar1=sfcg%zem)
  call shdf5_orec(ndims, idims, 'glatm'     , rvar1=sfcg%glatm)
  call shdf5_orec(ndims, idims, 'glonm'     , rvar1=sfcg%glonm)
  call shdf5_orec(ndims, idims, 'arm0'      , rvar1=sfcg%arm0)
! call shdf5_orec(ndims, idims, 'topm'      , rvar1=sfcg%topm)

  idims(1) = nvsfc

  call shdf5_orec(ndims, idims, 'xev'       , rvar1=sfcg%xev)
  call shdf5_orec(ndims, idims, 'yev'       , rvar1=sfcg%yev)
  call shdf5_orec(ndims, idims, 'zev'       , rvar1=sfcg%zev)
  call shdf5_orec(ndims, idims, 'dnu'       , rvar1=sfcg%dnu)
  call shdf5_orec(ndims, idims, 'dniu'      , rvar1=sfcg%dniu)
  call shdf5_orec(ndims, idims, 'dnv'       , rvar1=sfcg%dnv)
  call shdf5_orec(ndims, idims, 'dniv'      , rvar1=sfcg%dniv)
  call shdf5_orec(ndims, idims, 'unx'       , rvar1=sfcg%unx)
  call shdf5_orec(ndims, idims, 'uny'       , rvar1=sfcg%uny)
  call shdf5_orec(ndims, idims, 'unz'       , rvar1=sfcg%unz)
  call shdf5_orec(ndims, idims, 'vnx'       , rvar1=sfcg%vnx)
  call shdf5_orec(ndims, idims, 'vny'       , rvar1=sfcg%vny)
  call shdf5_orec(ndims, idims, 'vnz'       , rvar1=sfcg%vnz)

  idims(1) = nwsfc

  call shdf5_orec(ndims, idims, 'area'      , rvar1=sfcg%area)
  call shdf5_orec(ndims, idims, 'xew'       , rvar1=sfcg%xew)
  call shdf5_orec(ndims, idims, 'yew'       , rvar1=sfcg%yew)
  call shdf5_orec(ndims, idims, 'zew'       , rvar1=sfcg%zew)
  call shdf5_orec(ndims, idims, 'glatw'     , rvar1=sfcg%glatw)
  call shdf5_orec(ndims, idims, 'glonw'     , rvar1=sfcg%glonw)
  call shdf5_orec(ndims, idims, 'topw'      , rvar1=sfcg%topw)
  call shdf5_orec(ndims, idims, 'bathym'    , rvar1=sfcg%bathym)
  call shdf5_orec(ndims, idims, 'wnx'       , rvar1=sfcg%wnx)
  call shdf5_orec(ndims, idims, 'wny'       , rvar1=sfcg%wny)
  call shdf5_orec(ndims, idims, 'wnz'       , rvar1=sfcg%wnz)

  ! Write permanent grid cell data

  call shdf5_orec(ndims, idims, 'leaf_class', ivar1=sfcg%leaf_class)
  call shdf5_orec(ndims, idims, 'oge'       , ivar1=sfcg%ioge)
  call shdf5_orec(ndims, idims, 'swm_active', lvar1=sfcg%swm_active)

  idims(1) = nland

  call shdf5_orec(ndims, idims, 'usdatext'        , ivar1=land%usdatext)
  call shdf5_orec(ndims, idims, 'z_bedrock'       , rvar1=land%z_bedrock)
  call shdf5_orec(ndims, idims, 'gpp'             , rvar1=land%gpp)
  call shdf5_orec(ndims, idims, 'glhymps_ksat'    , rvar1=land%glhymps_ksat)
! call shdf5_orec(ndims, idims, 'glhymps_ksat_pfr', rvar1=land%glhymps_ksat_pfr)
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

  ndims    = 1
  idims(1) = nsea

  call shdf5_orec(ndims, idims, 'pom_active'      , lvar1=sea%pom_active)
  call shdf5_orec(ndims, idims, 'pom%kba'         , ivar1=pom%kba)

  call shdf5_close()

end subroutine sfcgfile_write

!==========================================================================

subroutine sfcgfile_read_pd()

  use max_dims,   only: pathlen
  use mem_sfcg,   only: nmsfc, mmsfc, nvsfc, mvsfc, nwsfc, mwsfc, &
                        itab_msfc_pd, itab_vsfc_pd, itab_wsfc_pd, &
                        itabg_msfc, itabg_vsfc, itabg_wsfc, sfcgfile
  use mem_land,   only: nland, mland, onland, omland, nzg
  use mem_lake,   only: nlake, mlake, onlake, omlake
  use mem_sea,    only: nsea, msea, onsea, omsea
  use pom2k1d,    only: nzpom

!UNN  use mem_sfcnud, only: nzg_nl, nzg_sp, kspm

  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  integer            :: imsfc, ivsfc, iwsfc
  integer            :: nzg0, nzpom0, ierr
  character(pathlen) :: flnm
  logical            :: there
  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  ierr = 0

  flnm = trim(sfcgfile)//'.h5'

  write(io6,*) 'sfcgfile_read_pd - checking sfcgfile (1) ',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'sfcgfile was not found - stopping run'
     stop 'stop: no sfcgfile'
  endif

  call shdf5_open(flnm,'R')

  ndims    = 1
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
  call shdf5_irec(ndims, idims, 'nzg'   , ivars=nzg0)
  call shdf5_irec(ndims, idims, 'nzpom' , ivars=nzpom0)

  if (nzg0   /= nzg  ) ierr = 1
  if (nzpom0 /= nzpom) ierr = 1

  if (ierr == 1) then

     write(io6,*) 'SFCGFILE mismatch with OLAMIN namelist: Stopping model run'
     write(io6,*) 'Values: gridfile, namelist'
     write(io6,*) '-----------------------------------------------'
     write(io6,*)              'nzg:      ', nzg0  , nzg
     write(io6,*)              'nzpom:    ', nzpom0, nzpom
     write(io6,*) '-----------------------------------------------'

     stop 'stop - surface gridfile mismatch'

  endif

  ! Copy grid dimensions

  mwsfc  = nwsfc
  mmsfc  = nmsfc
  mvsfc  = nvsfc
  mland  = nland
  mlake  = nlake
  msea   = nsea
  omland = onland
  omlake = onlake
  omsea  = onsea

  ! In this subroutine, do not read nzg_nl, nzg_sp, or kspm from the sfcgfile.
  ! When needed for groundwater initialization, they are read by another
  ! subroutine from a sfcgfile that was written by a separate simulation.

  write(io6, '(/,a)')    '==============================================='
  write(io6, '(a)')      'Reading from sfcgfile:'
  write(io6, '(2(a,i0))')'  nwsfc = ', nwsfc, ', nzg = ', nzg
  write(io6, '(a,/)')    '==============================================='

  ! Allocate permanent itabg data structures

  allocate (itab_msfc_pd(nmsfc))
  allocate (itab_vsfc_pd(nvsfc))
  allocate (itab_wsfc_pd(nwsfc))

  ! Allocate permanent itabg data structures

  allocate(itabg_msfc(nmsfc))
  allocate(itabg_vsfc(nvsfc))
  allocate(itabg_wsfc(nwsfc))

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

  call shdf5_irec(ndims,idims,'itab_msfc%imn',ivar2=iscr2)
  do imsfc = 1,nmsfc
     itab_msfc_pd(imsfc)%imn(1:3) = iscr2(1:3,imsfc)
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
!TR  call shdf5_irec(ndims,idims,'itab_wsfc%nwatm'     ,ivar1=iscr1)
!TR  do iwsfc = 1,nwsfc
!TR     itab_wsfc_pd(iwsfc)%nwatm = iscr1(iwsfc)
!TR  enddo

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

!TR  idims(1) = 8

!TR  allocate (iscr2(8,nwsfc))
!TR  call shdf5_irec(ndims,idims,'itab_wsfc%iwatm',ivar2=iscr2)
!TR  do iwsfc = 1,nwsfc
!TR     itab_wsfc_pd(iwsfc)%iwatm(1:8) = iscr2(1:8,iwsfc)
!TR  enddo

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
  use mem_sea,    only: sea, itab_sea, msea
  use pom2k1d,    only: nzpom, pom, y, yy, dy, dyy
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  character(pathlen) :: flnm
  logical            :: there
  integer            :: iw, iv
  character(2)       :: type

  ! Scratch arrays for copying input

  real, allocatable :: rscr2(:,:)

  ! Pointers to the global index of the local point

  integer :: lgmsfc(mmsfc)
  integer :: lgvsfc(mvsfc)
  integer :: lgwsfc(mwsfc)
  integer :: lgland(mland)
  integer :: lgsea (msea)

  lgmsfc = itab_msfc(1:mmsfc)%imglobe
  lgvsfc = itab_vsfc(1:mvsfc)%ivglobe
  lgwsfc = itab_wsfc(1:mwsfc)%iwglobe
  lgland = itab_land(1:mland)%iwglobe
  lgsea  = itab_sea (1:msea )%iwglobe

  ! Check if sfcgfile exists

  flnm = trim(sfcgfile)//'.h5'

  write(io6,*) 'sfcgfile_read - checking sfcgfile (2)',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'sfcgfile was not found - stopping run'
     stop 'stop: no sfcgfile'
  endif

  call shdf5_open(flnm,'R')

  ! Read ITAB_VSFC ARRAYS

  ndims    = 2
  idims(1) = 2
  idims(2) = mvsfc
  type     = 'CV'

  allocate (rscr2(2,mvsfc))
  call shdf5_irec(ndims,idims,'itab_vsfc%cosv',rvar2=rscr2,  points=lgvsfc)
  do iv = 1,mvsfc
     itab_vsfc(iv)%cosv(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_vsfc%sinv',rvar2=rscr2,  points=lgvsfc)
  do iv = 1,mvsfc
     itab_vsfc(iv)%sinv(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_vsfc%dxps',rvar2=rscr2,  points=lgvsfc)
  do iv = 1,mvsfc
     itab_vsfc(iv)%dxps(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_vsfc%dyps',rvar2=rscr2,  points=lgvsfc)
  do iv = 1,mvsfc
     itab_vsfc(iv)%dyps(1:2) = rscr2(1:2,iv)
  enddo
  deallocate(rscr2)

 ! Read ITAB_WSFC ARRAYS

  ndims    = 2
  idims(1) = 7
  idims(2) = mwsfc
  type     = 'CW'

  allocate (rscr2(7,mwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%dirv',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%dirv(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%farm',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%farm(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%farv',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%farv(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%gxps1',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%gxps1(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%gyps1',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%gyps1(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%gxps2',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%gxps2(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%gyps2',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%gyps2(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%ecvec_vx',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%ecvec_vx(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%ecvec_vy',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%ecvec_vy(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%ecvec_vz',rvar2=rscr2, points=lgwsfc)
  do iw = 1,mwsfc
     itab_wsfc(iw)%ecvec_vz(1:7) = rscr2(1:7,iw)
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

  idims(1) = nzpom

  call shdf5_irec(ndims, idims, 'y'   , rvar1=y)
  call shdf5_irec(ndims, idims, 'yy'  , rvar1=yy)
  call shdf5_irec(ndims, idims, 'dy'  , rvar1=dy)
  call shdf5_irec(ndims, idims, 'dyy' , rvar1=dyy)

  ndims    = 1
  idims(1) = mmsfc
  type     = 'CM'

  call shdf5_irec(ndims, idims, 'xem'       , rvar1=sfcg%xem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'yem'       , rvar1=sfcg%yem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'zem'       , rvar1=sfcg%zem,   points=lgmsfc)
  call shdf5_irec(ndims, idims, 'glatm'     , rvar1=sfcg%glatm, points=lgmsfc)
  call shdf5_irec(ndims, idims, 'glonm'     , rvar1=sfcg%glonm, points=lgmsfc)
  call shdf5_irec(ndims, idims, 'arm0'      , rvar1=sfcg%arm0,  points=lgmsfc)
! call shdf5_irec(ndims, idims, 'topm'      , rvar1=sfcg%topm,  points=lgmsfc)

  idims(1) = mvsfc
  type     = 'CV'

  call shdf5_irec(ndims, idims, 'xev'       , rvar1=sfcg%xev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'yev'       , rvar1=sfcg%yev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'zev'       , rvar1=sfcg%zev,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dnu'       , rvar1=sfcg%dnu,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dniu'      , rvar1=sfcg%dniu, points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dnv'       , rvar1=sfcg%dnv,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'dniv'      , rvar1=sfcg%dniv, points=lgvsfc)
  call shdf5_irec(ndims, idims, 'unx'       , rvar1=sfcg%unx,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'uny'       , rvar1=sfcg%uny,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'unz'       , rvar1=sfcg%unz,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'vnx'       , rvar1=sfcg%vnx,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'vny'       , rvar1=sfcg%vny,  points=lgvsfc)
  call shdf5_irec(ndims, idims, 'vnz'       , rvar1=sfcg%vnz,  points=lgvsfc)

  idims(1) = mwsfc
  type     = 'CW'

  call shdf5_irec(ndims, idims, 'area'      , rvar1=sfcg%area,    points=lgwsfc)
  call shdf5_irec(ndims, idims, 'xew'       , rvar1=sfcg%xew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'yew'       , rvar1=sfcg%yew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'zew'       , rvar1=sfcg%zew,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'glatw'     , rvar1=sfcg%glatw,   points=lgwsfc)
  call shdf5_irec(ndims, idims, 'glonw'     , rvar1=sfcg%glonw,   points=lgwsfc)
  call shdf5_irec(ndims, idims, 'topw'      , rvar1=sfcg%topw,    points=lgwsfc)
  call shdf5_irec(ndims, idims, 'bathym'    , rvar1=sfcg%bathym,  points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wnx'       , rvar1=sfcg%wnx,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wny'       , rvar1=sfcg%wny,     points=lgwsfc)
  call shdf5_irec(ndims, idims, 'wnz'       , rvar1=sfcg%wnz,     points=lgwsfc)

  ! Read permanent grid cell data

! call shdf5_irec(ndims, idims, 'leaf_class', ivar1=sfcg%leaf_class, points=lgwsfc)
  call shdf5_irec(ndims, idims, 'oge'       , ivar1=sfcg%ioge,       points=lgwsfc)
  call shdf5_irec(ndims, idims, 'swm_active', lvar1=sfcg%swm_active, points=lgwsfc)

  allocate (land%usdatext            (mland)) ; land%usdatext         = 0
  allocate (land%z_bedrock           (mland)) ; land%z_bedrock        = 0.
  allocate (land%gpp                 (mland)) ; land%gpp              = 0.
  allocate (land%glhymps_ksat        (mland)) ; land%glhymps_ksat     = 0.
! allocate (land%glhymps_ksat_pfr    (mland)) ; land%glhymps_ksat_pfr = 0.
  allocate (land%glhymps_poros       (mland)) ; land%glhymps_poros    = 0.
  allocate (land%sand            (nzg,mland)) ; land%sand             = 0.
  allocate (land%clay            (nzg,mland)) ; land%clay             = 0.
  allocate (land%silt            (nzg,mland)) ; land%silt             = 0.
  allocate (land%organ           (nzg,mland)) ; land%organ            = 0.
  allocate (land%bulkdens_drysoil(nzg,mland)) ; land%bulkdens_drysoil = 0.
  allocate (land%pH_soil         (nzg,mland)) ; land%pH_soil          = 0.
  allocate (land%cec_soil        (nzg,mland)) ; land%cec_soil         = 0.

  idims(1) = mland
  type     = 'LW'

  call shdf5_irec(ndims, idims, 'usdatext'        , ivar1=land%usdatext,        points=lgland)
  call shdf5_irec(ndims, idims, 'z_bedrock'       , rvar1=land%z_bedrock,       points=lgland)
  call shdf5_irec(ndims, idims, 'gpp'             , rvar1=land%gpp,             points=lgland)
  call shdf5_irec(ndims, idims, 'glhymps_ksat'    , rvar1=land%glhymps_ksat,    points=lgland)
! call shdf5_irec(ndims, idims, 'glhymps_ksat_pfr', rvar1=land%glhymps_ksat_pfr,points=lgland)
  call shdf5_irec(ndims, idims, 'glhymps_poros'   , rvar1=land%glhymps_poros,   points=lgland)

  ndims    = 2
  idims(1) = nzg
  idims(2) = mland
  type     = 'LW'

  call shdf5_irec(ndims, idims, 'sand'            , rvar2=land%sand,             points=lgland)
  call shdf5_irec(ndims, idims, 'clay'            , rvar2=land%clay,             points=lgland)
  call shdf5_irec(ndims, idims, 'silt'            , rvar2=land%silt,             points=lgland)
  call shdf5_irec(ndims, idims, 'organ'           , rvar2=land%organ,            points=lgland)
  call shdf5_irec(ndims, idims, 'bulkdens_drysoil', rvar2=land%bulkdens_drysoil, points=lgland)
  call shdf5_irec(ndims, idims, 'pH_soil'         , rvar2=land%pH_soil,          points=lgland)
  call shdf5_irec(ndims, idims, 'cec_soil'        , rvar2=land%cec_soil,         points=lgland)

  ndims    = 1
  idims(1) = msea
  type     = 'SW'

  call shdf5_irec(ndims, idims, 'pom_active', lvar1=sea%pom_active, points=lgsea)
  call shdf5_irec(ndims, idims, 'pom%kba'   , ivar1=pom%kba,        points=lgsea)

  call shdf5_close()

end subroutine sfcgfile_read

!==========================================================================

subroutine sfcgfile_read_makeregrid()

  use max_dims,   only: pathlen
  use mem_sfcg,   only: nwsfc, mwsfc, nmsfc, mmsfc, itab_wsfc, sfcg, sfcgfile
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
  use misc_coms,  only: io6

  implicit none

  integer            :: ndims, idims(2)
  integer            :: iwsfc
  character(pathlen) :: flnm
  logical            :: there
  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  flnm = trim(sfcgfile)//'.h5'

  write(io6,*) 'sfcgfile_read_makeregrid - checking sfcgfile (1) ',flnm

  inquire(file=flnm, exist=there)

  if (.not. there) then
     write(io6,*) 'sfcgfile was not found - stopping run'
     stop 'stop: no sfcgfile'
  endif

  call shdf5_open(flnm,'R')

  ndims    = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'nwsfc' , ivars=nwsfc)
  call shdf5_irec(ndims, idims, 'nmsfc' , ivars=nmsfc)

  ! Copy grid dimensions

  mwsfc  = nwsfc
  mmsfc  = nmsfc

  write(io6, '(/,a)')    '==============================================='
  write(io6, '(a)')      'Reading from sfcgfile:'
  write(io6, '(2(a,i0))')'  nwsfc = ', nwsfc
  write(io6, '(a,/)')    '==============================================='

  ! Allocate itabw data structure

  allocate (itab_wsfc(nwsfc))

  ! Read ITAB_WSFC SCALARS

  ndims    = 1
  idims(1) = nwsfc

  allocate (iscr1(nwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%npoly'     ,ivar1=iscr1)
  do iwsfc = 1,nwsfc
     itab_wsfc(iwsfc)%npoly = iscr1(iwsfc)
  enddo
  deallocate(iscr1)

  ! Read ITAB_WSFC ARRAYS

  ndims = 2
  idims(1) = 7
  idims(2) = nwsfc

  allocate (iscr2(7,nwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%imn',ivar2=iscr2)
  do iwsfc = 1,nwsfc
     itab_wsfc(iwsfc)%imn(1:7) = iscr2(1:7,iwsfc)
  enddo
  deallocate(iscr2)

  ! Allocate sfcg arrays

  allocate (sfcg%xem(nmsfc))
  allocate (sfcg%yem(nmsfc))
  allocate (sfcg%zem(nmsfc))

  allocate (sfcg%xew  (nwsfc))
  allocate (sfcg%yew  (nwsfc))
  allocate (sfcg%zew  (nwsfc))
  allocate (sfcg%glatw(nwsfc))
  allocate (sfcg%glonw(nwsfc))
  allocate (sfcg%topw (nwsfc))
  allocate (sfcg%area (nwsfc))

  allocate (sfcg%leaf_class (nwsfc))

  ndims    = 1
  idims(1) = nmsfc

  call shdf5_irec(ndims, idims, 'xem'        , rvar1=sfcg%xem)
  call shdf5_irec(ndims, idims, 'yem'        , rvar1=sfcg%yem)
  call shdf5_irec(ndims, idims, 'zem'        , rvar1=sfcg%zem)

  idims(1) = nwsfc

  call shdf5_irec(ndims, idims, 'area'       , rvar1=sfcg%area)
  call shdf5_irec(ndims, idims, 'xew'        , rvar1=sfcg%xew)
  call shdf5_irec(ndims, idims, 'yew'        , rvar1=sfcg%yew)
  call shdf5_irec(ndims, idims, 'zew'        , rvar1=sfcg%zew)
  call shdf5_irec(ndims, idims, 'glatw'      , rvar1=sfcg%glatw)
  call shdf5_irec(ndims, idims, 'glonw'      , rvar1=sfcg%glonw)
  call shdf5_irec(ndims, idims, 'topw'       , rvar1=sfcg%topw)

  call shdf5_irec(ndims, idims, 'leaf_class' , ivar1=sfcg%leaf_class)

  call shdf5_close()

end subroutine sfcgfile_read_makeregrid

!===============================================================================

subroutine landgrid_print()

  use mem_land,  only: nzg, slz, dslz, slzt
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
