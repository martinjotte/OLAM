subroutine leaf4_init_atm()

  use leaf_coms,    only: nzs, veg_ht, soil_rough, snow_rough, &
                          iupdndvi, s1900_ndvi, indvifile, nndvifiles, &
                          dt_leaf, isoilstateinit, iwatertabflg, watertab_db, &
                          wcap_min
  use mem_land,     only: mland, land, omland, nzg, slzt
  use leaf4_surface,only: grndvap
  use misc_coms,    only: s1900_sim, iparallel, runtype, initial
  use mem_sfcg,     only: itab_wsfc, sfcg
  use consts_coms,  only: cliq, cice, alli, cliq1000, cice1000, alli1000, &
                          grav, p00i, rocp, t00
  use mem_para,     only: myrank
  use leaf4_canopy, only: vegndvi
  use land_db,      only: land_database_read
  use leaf4_soil,   only: soil_pot2wat, soil_wat2pot
  use oname_coms,   only: nl

  implicit none

  integer :: k
  integer :: iland, iwsfc
  integer :: leaf_class
  integer :: nlsw1 ! maximum of (1,nlev_sfcwater)

  real :: timefac_ndvi
  real :: wq, wq_added, snowdens
  real :: psi, psi_slope

  real :: soil_tempc(nzg,mland)  ! initial soil temperature (C)
  real :: fracliq(nzg,mland)     ! initial soil liquid fraction (0-1)
  real :: wtd(mland)             ! watertable depth from database

  real, external :: rhovsl

  ! Leaf quantities that get initialized at the start of any model run

  timefac_ndvi = 0.

  if (iupdndvi == 1 .and. nndvifiles > 1) then
     timefac_ndvi = (s1900_sim               - s1900_ndvi(indvifile))  &
                  / (s1900_ndvi(indvifile+1) - s1900_ndvi(indvifile))
  endif

  ! Allocate and fill temporary lat/lon arrays for land points only

  ! Subgrid orographic roughness

  land%slope_fact(1) = 1.0

  if (nl%iorogslopeflg == 1 .and. len_trim(nl%orog_slope_db) > 0) then

     ! assume omland = 0
     call land_database_read(mland, &
          sfcg%glatw,               &
          sfcg%glonw,               &
          nl%orog_slope_db,         &
          nl%orog_slope_db,         &
          'orog',                   &
          datq=wtd                  )

     !$omp parallel do
     do iland = 2, mland
        land%slope_fact(iland) = 1.0 + 2.0 * wtd(iland)
     enddo
     !$omp end parallel do

  else

     !$omp parallel do
     do iland = 2, mland
        land%slope_fact(iland) = 1.0
     enddo
     !$omp end parallel do

  endif

  ! If iwatertabflg = 1, read water table depth from database; put into wtd array

  if (iwatertabflg == 1) then

     ! assume omland = 0
     call land_database_read(mland, &
          sfcg%glatw,               &
          sfcg%glonw,               &
          watertab_db,              &
          watertab_db,              &
          'wtd',                    &
          datq=wtd                  )

  endif

  !$omp parallel do private (iwsfc,leaf_class)
  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Set vegetation parameters

     leaf_class = sfcg%leaf_class(iwsfc)

     land%veg_rough  (iland) = .13 * veg_ht(leaf_class)
     sfcg%rough      (iwsfc) = max(soil_rough, land%veg_rough(iland))
     land%veg_height (iland) = veg_ht(leaf_class)
     land%stom_resist(iland) = 1.e6

     ! For now, choose heat/vapor capacities for stability based on timestep

     sfcg%can_depth(iwsfc) = 20. * max(1.,.030 * dt_leaf)
     land%hcapveg  (iland) = 3.e4 * max(1.,.025 * dt_leaf)

     ! Initialize vegetation TAI, LAI, fractional area, albedo, and roughness

     call vegndvi(iland, iwsfc, timefac_ndvi, &
                  sfcg%leaf_class  (iwsfc), &
                  land%veg_height  (iland), &
                  land%veg_rough   (iland), &
                  land%veg_ndvip   (iland), &
                  land%veg_ndvif   (iland), &
                  land%veg_ndvic   (iland), &
                  land%veg_tai     (iland), &
                  land%veg_lai     (iland), &
                  land%veg_fracarea(iland), &
                  land%veg_albedo  (iland)  )

     ! Initialize hydraulic head at bottom of soil grid...

     if (leaf_class == 17 .or. leaf_class == 20) then

        ! For bog, marsh, wetland areas, use head0 = +10 cm; this takes precedence
        ! over using watertable database

        land%head0(iland) = 0.1

     elseif (iwatertabflg == 1 .and. wtd(iland) > -1.e-6) then

        ! If area is not bog, marsh, or wetland, and if iwatertabflg = 1, use
        ! head0 = -watertable depth from database, unless it is "missing" (negative)

        land%head0(iland) = -wtd(iland)

     elseif (leaf_class == 3) then

        ! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is desert

        land%head0(iland) = -100.0  ! Desert

     elseif (leaf_class == 10) then

        ! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is semi-desert

        land%head0(iland) = -30.0  ! Semi-desert

     else

        ! If iwatertabflg /= 1 (or wtd is "missing") and landuse class is none of the above

        land%head0(iland) = -10.0

     endif

  enddo  ! iland

  ! End of leaf quantities that are initialized on any run

  if (runtype /= "INITIAL") return

  ! Leaf quantities that are initialized only on 'INITIAL' run

  !$omp parallel do private (iwsfc, k, psi)
  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply initial atmospheric properties to land canopy

     sfcg%cantemp  (iwsfc) = sfcg%airtheta(iwsfc) * (sfcg%prss(iwsfc) * p00i)**rocp
     sfcg%canrrv   (iwsfc) = sfcg%airrrv(iwsfc)
     land%veg_temp (iland) = sfcg%cantemp(iwsfc)
     land%veg_water(iland) = 0.
     sfcg%ustar    (iwsfc) = 0.1
     sfcg%sfluxt   (iwsfc) = 0.
     sfcg%sfluxr   (iwsfc) = 0.
     sfcg%wthv     (iwsfc) = 0.

     land%veg_energy(iland) = (land%veg_temp(iland) - 273.15) * land%hcapveg(iland)

     ! Default initialization of sfcwater_mass, soil_tempc, and soil_water

     land%sfcwater_mass  (1:nzs,iland) = 0.
     land%sfcwater_energy(1:nzs,iland) = 0.
     land%sfcwater_depth (1:nzs,iland) = 0.

     soil_tempc(1:nzg,iland) = sfcg%cantemp(iwsfc) - 273.15
     fracliq(1:nzg,iland) = 1.0

     ! Loop over soil layers

     do k = 1,nzg

        ! Diagnose pressure head and soil water from given head0 and slzt

        psi = land%head0(iland) - slzt(k)

        call soil_pot2wat(psi, land%wresid_vg(k,iland), land%wsat_vg(k,iland), &
                          land%alpha_vg(k,iland), land%en_vg(k,iland), &
                          land%soil_water(k,iland))

     enddo

  enddo  ! iland
  !$omp end parallel do

  ! Overwrite the default soil initialization with observed data if specified

  ! read sfcwater mass, soil_tempc, and soil_water from the 2 X 2.5 degree
  ! NCEP/NCAR reanalysis in netcdf format, obtained from
  ! ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/
  ! call read_soil_moist_temp(soil_tempc)

  if (isoilstateinit > 0 .and. initial == 2) then

     ! read sfcwater mass, soil_tempc, and soil_water saved in the initial
     ! degribbed analysis files

     call read_soil_analysis(soil_tempc)

  endif

  ! Loop over all LAND cells

  !$omp parallel do private(iwsfc, k, wq, wq_added, snowdens, nlsw1, psi)
  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     !--------------------------------------------------------------------------------
     ! ADD A METHOD HERE TO INITIALIZE FRACTION OF SOIL WATER THAT IS LIQUID (FRACLIQ)
     ! BASED ON MODEL OBSERVED AND/OR MODEL CLIMATOLOGY
     !--------------------------------------------------------------------------------

     ! If leaf_class of this iland land cell is ice cap or glacier, assume that
     ! all 'soil' water is in ice phase and that soil_tempc is at or below 0.
     ! (Soil will be replaced by firn model in the future.)

     if (sfcg%leaf_class(iwsfc) == 2) then

        soil_tempc(1:nzg,iland) = min(0.,soil_tempc(1:nzg,iland))
        fracliq(1:nzg,iland) = 0.

     endif

     ! Initialize soil energy [J/m^3] from given soil textural class, temperature,
     ! total water content, and liquid fraction (as opposed to ice fraction) of the
     ! soil water that is present.

     do k = 1,nzg

        if (soil_tempc(k,iland) > 0.) then

           land%soil_energy(k,iland)                                        &
               = soil_tempc(k,iland) * land%specifheat_drysoil(k,iland)     &
               + soil_tempc(k,iland) * land%soil_water(k,iland) * cliq1000  &
               + fracliq(k,iland)    * land%soil_water(k,iland) * alli1000

        else

           land%soil_energy(k,iland)                                       &
              = soil_tempc(k,iland) * land%specifheat_drysoil(k,iland)     &
              + soil_tempc(k,iland) * land%soil_water(k,iland) * cice1000  &
              + fracliq(k,iland)    * land%soil_water(k,iland) * alli1000

        endif
     enddo

     ! Leaf classes 17 and 20 represent persistent wetlands (bogs, marshes,
     ! fens, swamps).  Initialize these areas with 0.1 m of standing surface
     ! water (sfcwater) added to whatever is already present (e.g., from obs).

     if (sfcg%leaf_class(iwsfc) == 17 .or. sfcg%leaf_class(iwsfc) == 20) then

        ! Since sfcwater_energy has units of J/kg, first convert to J/m^2 before adding
        ! wetland sfcwater.

        wq = land%sfcwater_mass(1,iland) * land%sfcwater_energy(1,iland)

        ! Add wetland sfcwater mass and depth

        land%sfcwater_mass(1,iland) = land%sfcwater_mass(1,iland)  &
                                    + 100.  ! 100 kg/m^2 equivalent to 0.1 m depth
        land%sfcwater_depth(1,iland) = land%sfcwater_depth(1,iland)  &
                                     + .1   ! 0.1 m added depth

        ! Add wetland sfcwater energy, which is assumed to have energy of liquid
        ! water at canopy air temperature, which could be below freezing.

        wq_added = 100.  &  ! 100 kg/m^2 added mass
                 * ((sfcg%cantemp(iwsfc) - 273.15) * cliq + alli) ! J/kg of added
                                                                ! liquid water
        ! Diagnose new sfcwater energy

        land%sfcwater_energy(1,iland) = (wq + wq_added) / land%sfcwater_mass(1,iland)

     endif

     ! Leaf classes 2 represents persistent ice/snow (icecap, glacier).
     ! Initialize these areas with at least 100 kg/m^2 of frozen surface
     ! water (sfcwater).

     if (sfcg%leaf_class(iwsfc) == 2) then

        ! Add wetland sfcwater mass
        land%sfcwater_mass(1,iland) = min(land%sfcwater_mass(1,iland), 100.)

        ! Snow is frozen at canopy temperature
        land%sfcwater_energy(1,iland) = min(0., (sfcg%cantemp(iwsfc) - t00) * cice)

        ! Snow density calculation comes from CLM3.0 documentation,
        ! which is based on Anderson 1975 NWS Technical Doc # 19
        if (sfcg%cantemp(iwsfc) > 258.15) then
           snowdens = 50.0 + 1.5 * (sfcg%cantemp(iwsfc) - 258.15)**1.5
        else
           snowdens = 50.0
        endif
        land%sfcwater_depth(1,iland) = land%sfcwater_mass(1,iland) / snowdens

     endif

     ! Determine active number of surface water levels

     land%nlev_sfcwater(iland) = 0

     do k = 1,nzs
        if (land%sfcwater_mass(k,iland) < wcap_min) then
           land%sfcwater_mass  (k:nzs,iland) = 0.
           land%sfcwater_energy(k:nzs,iland) = 0.
           land%sfcwater_depth (k:nzs,iland) = 0.
           exit
        else
           land%nlev_sfcwater(iland) = k
        endif
     enddo

     ! Initialize ground (soil) and surface vapor specific humidity

     nlsw1 = max(1,land%nlev_sfcwater(iland))

     call soil_wat2pot(nzg, iland, land%soil_water(nzg,iland), land%wresid_vg(nzg,iland), &
                       land%wsat_vg(nzg,iland), land%alpha_vg(nzg,iland), &
                       land%en_vg(nzg,iland), psi, psi_slope)

     call grndvap(iland,                              &
                  sfcg%rhos                  (iwsfc), &
                  sfcg%canrrv                (iwsfc), &
                  land%nlev_sfcwater         (iland), &
                  land%surface_srrv          (iland), &
                  land%ground_rrv            (iland), &
                  land%sfcwater_energy (nlsw1,iland), &
                  land%soil_water        (nzg,iland), &
                  land%soil_energy       (nzg,iland), &
                  psi                               , &
                  land%specifheat_drysoil(nzg,iland)  )

     ! Initialize snowfac

     land%snowfac(iland) = 0.
     do k = 1,land%nlev_sfcwater(iland)
        land%snowfac(iland) = land%snowfac(iland) + land%sfcwater_depth(k,iland)
     enddo
     land%snowfac(iland) = land%snowfac(iland) / max(.001,land%veg_height(iland))
     if (land%snowfac(iland) > 0.9) land%snowfac(iland) = 1.0

     sfcg%rough(iwsfc) = max(snow_rough, &
        max(soil_rough, land%veg_rough(iland)) * (1. - land%snowfac(iland)))

     ! Diagnose head1 based on sfcwater mass

     sfcg%head1(iwsfc) = .001 * sum(land%sfcwater_mass(1:land%nlev_sfcwater(iland),iland))

  enddo

end subroutine leaf4_init_atm
