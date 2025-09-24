subroutine landcells(iland, timefac_ndvi, head_slope, soil_watfrac)

  use leaf_coms,     only: wcap_min, soil_rough, snow_rough, &
                           thermcond_dry_organic, thermcond_sat_organic, &
                           thermcond_bedrock, thermcond_liq, thermcond_ice, &
                           thermcond_firn
  use mem_land,      only: land, omland, nzg, nzs, dslz
  use misc_coms,     only: iparallel, time_istp8
  use consts_coms,   only: grav, p00i, rocp, eps_virt, cpi, r8
  use mem_sfcg,      only: itab_wsfc, sfcg
  use leaf4_canopy,  only: canopy, vegndvi, fast_canopy
  use leaf4_surface, only: skncomp_diagnose, sfcwater, remove_runoff
  use leaf4_soil,    only: soil, soil_wat2khyd
  use leaf4_plot,    only: leaf_plot
  use therm_lib,     only: qwtk, qtk
  use mem_para,      only: myrank
  use oname_coms,    only: nl
  use mem_sfcnud,    only: sfcwat_nud, sfctemp_nud, fracliq_nud

  implicit none

  integer, intent(in)    :: iland
  real,    intent(in)    :: timefac_ndvi
  real,    intent(in)    :: head_slope  (nzg+1)
  real,    intent(inout) :: soil_watfrac(nzg+1)

  ! Local arrays

  real :: sfcwater_tempk  (nzs)   ! surface water temperature [K]
  real :: sfcwater_fracliq(nzs)   ! fraction of sfc water in liquid phase
  real :: sfcwater_epm2   (nzs)   ! sfcwater energy per m^2 [J/m^2]

  real :: soil_tempk      (nzg)   ! soil temperature [K]
  real :: soil_fracliq    (nzg)   ! fraction of soil moisture in liquid phase
  real :: thermcond_soil  (nzg)   ! soil thermal conductivity [W/(K m)]
  real :: dheight         (nzg)   ! change in water height (of a T cell) from lateral water fluxes [m]
  real :: energyin        (nzg)   ! change in energy (of T cell) from lateral water fluxes [J/m^2]

  ! Local variables

  integer :: iwsfc

  real :: canexneri, cantheta, canthetav
  real :: airthetav, ufree

  integer :: k      ! vertical index over soil layers
  integer :: ktrans ! vertical index of soil layer supplying transpiration

  real :: transp  ! transpiration xfer this LEAF timestep [kg/m^2]
  real :: wfree1  ! sfcwater bottom layer free water mass [kg/m^2]
  real :: qwfree1 ! sfcwater bottom layer free water energy [J/m^2]
  real :: dwfree1 ! sfcwater bottom layer free water depth [m]

  real :: thermcond_dry_mineral
  real :: thermcond_sat_mineral

  real :: mineral
  real :: bulkdens_mineral
  real :: thermcond_drysoil
  real :: thermcond_sat_solids
  real :: thermcond_sat_soil
  real :: kersten_liq, kersten_ice
  real :: kersten

  integer, parameter :: iland_print = 0

  real, parameter :: onethird = 1./3.

  iwsfc = iland - omland

  ! Skip this cell if running in parallel and cell rank is not MYRANK

  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  ! Initialize local sfcwater arrays to zero

  sfcwater_tempk  (:) = 0.
  sfcwater_fracliq(:) = 0.
  sfcwater_epm2   (:) = 0.

  ! Diagnose current sfcwater properties and skinlayer composition based on those
  ! properties.  Make adjustments to sfcwater and top of soil as appropriate.

  call skncomp_diagnose(iland, iwsfc,                     &
                        land%skncomp           (  iland), &
                        land%sfcwater_mass     (:,iland), &
                        land%sfcwater_energy   (:,iland), &
                        land%sfcwater_depth    (:,iland), &
                        land%soil_water        (:,iland), &
                        land%soil_energy       (:,iland), &
                        land%specifheat_drysoil(:,iland)  )

!!  if ( land%sfcwater_mass(1,iland) > wcap_min ) &
!!       call qtk(land%sfcwater_energy(1,iland), sfcwater_tempk(1), sfcwater_fracliq(1))
!!
!!  if ( land%skncomp(iland) == 2 ) &
!!       call qtk(land%sfcwater_energy(2,iland), sfcwater_tempk(2), sfcwater_fracliq(2))
!!
!!  if (itab_wsfc(iwsfc)%iwglobe == 5326) then
!!     write(*,*) sfcg%glatw(iwsfc), sfcg%glonw(iwsfc)
!!     write(*,*) land%skncomp           (  iland)
!!     write(*,*) sfcwater_tempk(:), sfcwater_fracliq(:)
!!     write(*,*) land%sfcwater_mass     (:,iland)
!!  endif

  ! Diagnose soil temperature and liquid fraction

  do k = 1, nzg
     call qwtk(land%soil_energy(k,iland),land%soil_water(k,iland)*1.e3, &
               land%specifheat_drysoil(k,iland),soil_tempk(k),soil_fracliq(k))
!!     if (itab_wsfc(iwsfc)%iwglobe == 5326) write(*,*) k, soil_tempk(k),soil_fracliq(k), dslz(k)
  enddo

  ! Diagnose thermal conductivity

  if (sfcg%leaf_class(iwsfc) == 2) then

     ! Case for firn/glacier locations
     do k = 1, nzg
        thermcond_soil(k) = thermcond_firn
     enddo

  else

     ! Case for bedrock
     do k = 1, land%k_bedrock(iland)
        thermcond_soil(k) = thermcond_bedrock  ! Assuming that porosity is accounted for.
     enddo

     ! Case for soil
     do k = land%k_bedrock(iland), nzg

        ! Bulk density of mineral part of soil

        ! (Note: bulkdens_mineral uses wsat_vg which in turn uses bulkdens_drysoil, which
        ! includes both mineral and organic components.  Is there a more direct way
        ! to get bulkdens_mineral from bulkdens_drysoil?  One idea is:

        ! bulkdens_mineral = land%bulkdens_drysoil(k,iland) * (1. - land%organ(k,iland))

        ! This formula would presumably be more accurate for peat soil, for example,
        ! but would it ultimately give better values for thermcond_dry_mineral and
        ! thermcond_drysoil?

        bulkdens_mineral = 2700. * (1. - land%wsat_vg(k,iland))

        ! Thermal conductivity for dry soil material

        mineral = 1. - land%organ(k,iland)

        thermcond_dry_mineral = (0.135 * bulkdens_mineral + 64.7) &
                              / (2700. - 0.947 * bulkdens_mineral)

        thermcond_drysoil = mineral             * thermcond_dry_mineral &
                          + land%organ(k,iland) * thermcond_dry_organic

        ! Thermal conductivity for saturated soil material

!       thermcond_sat_mineral = (8.80 * land%sand(k,iland) + 2.92 * land%clay(k,iland)) &
!                             / (land%sand(k,iland) + land%clay(k,iland))
        thermcond_sat_mineral = 8.80 * land%sand(k,iland) + 2.92 * (1. - land%sand(k,iland))

        thermcond_sat_solids = mineral             * thermcond_sat_mineral &
                             + land%organ(k,iland) * thermcond_sat_organic

        ! For soil case, use porosity and fractions of sand, clay, and organic matter

        ! Compute thermcond_sat_soil as geometric mean of solids, liquid, and ice contributions

        thermcond_sat_soil = thermcond_sat_solids**(1.-land%wsat_vg(k,iland)) &
                           * thermcond_liq**(land%wsat_vg(k,iland)*soil_fracliq(k)) &
                           * thermcond_ice**(land%wsat_vg(k,iland)*(1.-soil_fracliq(k)))

        kersten_liq = 1. + log10( max(.1, soil_watfrac(k)) )
        kersten_ice = soil_watfrac(k)
        kersten = kersten_liq * soil_fracliq(k) + kersten_ice * (1. - soil_fracliq(k))

        thermcond_soil(k) = kersten  * thermcond_sat_soil &
                    + (1. - kersten) * thermcond_drysoil
     enddo

  endif

  if (nl%igw_spinup == 1) then

     ! With igw_spinup = 1, call subroutine fast_canopy to nudge surface water
     ! mass and energy

     ktrans = nzg
     transp = 0.

     call fast_canopy(iland, iwsfc,                     &
                      sfcwat_nud             (iwsfc),   &
                      sfctemp_nud            (iwsfc),   &
                      fracliq_nud            (iwsfc),   &
                      land%sfcwater_mass     (1,iland), &
                      land%sfcwater_energy   (1,iland), &
                      land%sfcwater_depth    (1,iland), &
                      sfcwater_tempk         (1),       &
                      sfcwater_fracliq       (1),       &
                      sfcwater_epm2          (1),       &
                      land%soil_water        (:,iland), &
                      land%soil_energy       (:,iland), &
                      land%specifheat_drysoil(:,iland), &
                      soil_tempk             (:),       &
                      soil_fracliq           (:)        )

  else

     ! With igw_spinup /= 1, call subroutines vegndvi and canopy for standard
     ! atmosphere-surface interaction

     ! Update vegetation TAI, LAI, fractional area, albedo, and roughness

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

     ! Compute roughness length based on vegetation and snow.

     sfcg%rough(iwsfc) = max(max(soil_rough, land%veg_rough(iland)) &
                       * (1. - land%snowfac(iland)), snow_rough)

     ! Evaluate surface layer exchange coefficients vkmsfc and vkhsfc

     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))
     canexneri = 1. / sfcg%canexner(iwsfc)

     cantheta  = sfcg%cantemp(iwsfc) * canexneri
     canthetav = cantheta * (1.0 + eps_virt * sfcg%canrrv(iwsfc))

     ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sfcg%wthv(iwsfc),0.0) / airthetav) ** onethird

     if (time_istp8 < 1.e-3_r8) ufree = max(ufree, 2.0)

     call stars(sfcg%dzt_bot (iwsfc), &
                sfcg%rough   (iwsfc), &
                sfcg%vels    (iwsfc), &
                sfcg%rhos    (iwsfc), &
                ufree               , &
                airthetav           , &
                canthetav           , &
                sfcg%vkmsfc  (iwsfc), &
                sfcg%vkhsfc  (iwsfc), &
                sfcg%ustar   (iwsfc), &
                sfcg%ggaer   (iwsfc)  )

     if (nl%iorogslopeflg > 1) then
        sfcg%vkmsfc(iwsfc) = sfcg%vkmsfc(iwsfc) * land%slope_fact(iland)
     endif

     ! Evaluate turbulent exchanges of heat and moisture between vegetation and canopy air
     ! and also between soil or snow surface and canopy air.  Evaluate transfer of
     ! precipitation moisture and heat to vegetation and shed from vegetation to surface.
     ! Update vegetation and canopy air temperatures resulting from these
     ! plus radiative fluxes.

     call canopy(iland, iwsfc, ktrans, transp,     &
                 sfcg%leaf_class        (  iwsfc), &
                 sfcg%can_depth         (  iwsfc), &
                 sfcg%rhos              (  iwsfc), &
                 sfcg%vels              (  iwsfc), &
                 sfcg%ustar             (  iwsfc), &
                 sfcg%vkhsfc            (  iwsfc), &
                 sfcg%sfluxt            (  iwsfc), &
                 sfcg%sfluxr            (  iwsfc), &
                 sfcg%pcpg              (  iwsfc), &
                 sfcg%qpcpg             (  iwsfc), &
                 sfcg%dpcpg             (  iwsfc), &
                 sfcg%rshort            (  iwsfc), &
                 sfcg%cantemp           (  iwsfc), &
                 sfcg%canrrv            (  iwsfc), &
                 sfcg%glatw             (  iwsfc), &
                 sfcg%glonw             (  iwsfc), &
                 sfcg%airtheta          (  iwsfc), &
                 sfcg%airrrv            (  iwsfc), &
                 sfcg%canexner          (  iwsfc), &
                 land%snowfac           (  iland), &
                 land%vf                (  iland), &
                 land%stom_resist       (  iland), &
                 land%veg_height        (  iland), &
                 land%veg_rough         (  iland), &
                 land%veg_tai           (  iland), &
                 land%veg_lai           (  iland), &
                 land%hcapveg           (  iland), &
                 land%rshort_s          (  iland), &
                 land%rshort_v          (  iland), &
                 land%rlong_s           (  iland), &
                 land%rlong_v           (  iland), &
                 land%veg_water         (  iland), &
                 land%veg_energy        (  iland), &
                 land%veg_temp          (  iland), &
                 land%skncomp           (  iland), &
                 land%sfcwater_mass     (:,iland), &
                 land%sfcwater_energy   (:,iland), &
                 land%sfcwater_depth    (:,iland), &
                 sfcwater_tempk         (:),       &
                 sfcwater_fracliq       (:),       &
                 sfcwater_epm2          (:),       &
                 soil_watfrac           (:),       &
                 land%soil_water        (:,iland), &
                 land%soil_energy       (:,iland), &
                 land%wsat_vg           (:,iland), &
                 land%wresid_vg         (:,iland), &
                 land%soilfldcap        (  iland), &
                 land%ksat_vg           (:,iland), &
                 land%specifheat_drysoil(:,iland), &
                 land%head              (:,iland), &
                 head_slope             (:      ), &
                 soil_tempk             (:      ), &
                 soil_fracliq           (:      )  )

  endif

! Original calculation of wthv for sfluxt units [kg_dry K m^-2 s^-1]

!  sfcg%wthv(iwsfc) = ( sfcg%sfluxt(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
!       + sfcg%sfluxr(iwsfc) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

  ! New calculation of wthv for sfluxt units [W m^-2]

  sfcg%wthv(iwsfc) = ( sfcg%sfluxt(iwsfc) * cpi * canexneri * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
                   + sfcg%sfluxr(iwsfc) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

  ! CALL SFCWATER:
  !  1. Compute soil and sfcwater heat conductivities
  !  2. Compute surface heat flux for top layer of soil or sfcwater
  !  3. Evaluate internal and bottom fluxes for sfcwater
  !  4. Update sfcwater layer energies due to heat flux and solar radiation
  !  5. Evaluate melting and percolation of liquid through sfcwater layers

  call sfcwater(iland, iwsfc, wfree1, qwfree1, dwfree1, &
                sfcg%head1             (  iwsfc), &
                land%skncomp           (  iland), &
                land%sfcwater_mass     (:,iland), &
                land%sfcwater_energy   (:,iland), &
                land%sfcwater_depth    (:,iland), &
                sfcwater_tempk         (:      ), &
                sfcwater_fracliq       (:      ), &
                sfcwater_epm2          (:      ), &
                land%soil_water        (:,iland), &
                land%soil_energy       (:,iland), &
                land%specifheat_drysoil(:,iland), &
                soil_tempk             (:      ), &
                thermcond_soil         (:      )  )

  call soil(iland, iwsfc, ktrans, transp, wfree1, qwfree1, &
            sfcg%leaf_class     (  iwsfc), &
            sfcg%glatw          (  iwsfc), &
            sfcg%glonw          (  iwsfc), &
            sfcg%head1          (  iwsfc), &
            land%head0          (  iland), &
            land%gpp            (  iland), &
            land%skncomp        (  iland), &
            land%sfcwater_mass  (1,iland), &
            land%sfcwater_energy(1,iland), &
            land%sfcwater_depth (1,iland), &
            sfcwater_epm2       (1      ), &
            land%soil_water     (:,iland), &
            land%soil_energy    (:,iland), &
            land%wresid_vg      (:,iland), &
            land%wsat_vg        (:,iland), &
            land%ksat_vg        (:,iland), &
            land%lambda_vg      (:,iland), &
            land%en_vg          (:,iland), &
            soil_watfrac        (:      ), &
            head_slope          (:      ), &
            land%head           (:,iland), &
            soil_tempk          (:      ), &
            soil_fracliq        (:      ), &
            thermcond_soil      (:      ), &
            dheight             (:      ), & ! included here only to pass to leaf_plot
            energyin            (:      )  )

  !-----------------------------------------------------------------------------
  ! TEMPORARY UNTIL FULL LEAF-HYDRO MODEL IS COMPLETED WITH STREAM/RIVER RUNOFF:
  ! Simple representation of runoff

  if (land%sfcwater_mass(1,iland) >= wcap_min) then

     call remove_runoff(iland, iwsfc,                  &
                        sfcg%leaf_class       (iwsfc), &
!needed?                land%skncomp          (iland), &
                        land%sfcwater_mass  (1,iland), &
                        land%sfcwater_energy(1,iland), &
                        land%sfcwater_depth (1,iland), &
                        sfcg%runoff           (iwsfc)  )

  endif

  !-----------------------------------------------------------------------------

  ! Call land patch plot routine for selected iland values.

  if (iland == iland_print) call leaf_plot(               &
        iland,                                            &
        linit            = 1,                             &
        lframe           = 1,                             &
        ktrans           = ktrans,                        &
        sfcwater_tempk   = sfcwater_tempk      (:      ), &
        sfcwater_fracliq = sfcwater_fracliq    (:      ), &
        soil_tempk       = soil_tempk          (:      ), &
        soil_fracliq     = soil_fracliq        (:      ), &
        leaf_class       = sfcg%leaf_class     (  iwsfc), &
        can_depth        = sfcg%can_depth      (  iwsfc), &
        rhos             = sfcg%rhos           (  iwsfc), &
        vels             = sfcg%vels           (  iwsfc), &
        prss             = sfcg%prss           (  iwsfc), &
        pcpg             = sfcg%pcpg           (  iwsfc), &
        qpcpg            = sfcg%qpcpg          (  iwsfc), &
        dpcpg            = sfcg%dpcpg          (  iwsfc), &
        ustar            = sfcg%ustar          (  iwsfc), &
        cantemp          = sfcg%cantemp        (  iwsfc), &
        canrrv           = sfcg%canrrv         (  iwsfc), &
        head1            = sfcg%head1          (  iwsfc), &
        head0            = land%head0          (  iland), &
        rshort           = sfcg%rshort         (  iwsfc), &
        rshort_s         = land%rshort_s       (  iland), &
        rshort_v         = land%rshort_v       (  iland), &
        rlong_s          = land%rlong_s        (  iland), &
        rlong_v          = land%rlong_v        (  iland), &
        veg_height       = land%veg_height     (  iland), &
        veg_rough        = land%veg_rough      (  iland), &
        veg_tai          = land%veg_tai        (  iland), &
        veg_lai          = land%veg_lai        (  iland), &
        hcapveg          = land%hcapveg        (  iland), &
        snowfac          = land%snowfac        (  iland), &
        vf               = land%vf             (  iland), &
        veg_water        = land%veg_water      (  iland), &
        veg_temp         = land%veg_temp       (  iland), &
        transp           = transp                       , &
        stom_resist      = land%stom_resist    (  iland), &
        veg_fracarea     = land%veg_fracarea   (  iland), &
        skncomp          = land%skncomp        (  iland), &
        sfcwater_mass    = land%sfcwater_mass  (:,iland), &
        sfcwater_energy  = land%sfcwater_energy(:,iland), &
        sfcwater_depth   = land%sfcwater_depth (:,iland), &
        soil_water       = land%soil_water     (:,iland), &
        soil_energy      = land%soil_energy    (:,iland), &
        head             = land%head           (:,iland)  )

end subroutine landcells
