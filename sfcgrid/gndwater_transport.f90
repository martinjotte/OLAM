subroutine gndwater_transport(head_slope, soil_watfrac)

  use mem_sfcg,      only: mvsfc
  use mem_land,      only: land, mland, nzg
  use mem_lake,      only: lake, mlake
  use misc_coms,     only: iparallel
  use leaf4_canopy,  only: canopy, vegndvi, fast_canopy
  use olam_mpi_sfc,  only: mpi_send_vsfc, mpi_recv_vsfc

  implicit none

  real, intent(in)    :: head_slope  (nzg,mland)
  real, intent(inout) :: soil_watfrac(nzg,mland)

  ! Local variables

  integer :: iland, ilake
  integer :: ivsfc

  real, allocatable :: watflux   (:,:) ! [m/s]
  real, allocatable :: energyflux(:,:)

  allocate(watflux   (nzg,mvsfc))
  allocate(energyflux(nzg,mvsfc))

  !$omp parallel do
  do ivsfc = 2,mvsfc

     watflux   (:,ivsfc) = 0.
     energyflux(:,ivsfc) = 0.

     ! Evaluate water and energy fluxes across V face, with flux direction
     ! fdefined by flux = (head(v%iw1) - head(v%iw2)) / resist

     call comp_horiz_grndwatfluxes( ivsfc, soil_watfrac, &
                                    watflux(:,ivsfc), energyflux(:,ivsfc) )
  enddo
  !$omp end parallel do

  ! MPI SEND/RECV watflux, energyflux

  if (iparallel == 1) then
     call mpi_send_vsfc(watflux=watflux,energyflux=energyflux)
     call mpi_recv_vsfc(watflux=watflux,energyflux=energyflux)
  endif

  !-----------------------------------------------------------------------
  ! Loop over ALL LAKE CELLS
  !-----------------------------------------------------------------------

  !$omp parallel
  !$omp do
  do ilake = 2, mlake

     call apply_lake_fluxes( ilake, watflux, energyflux, &
                             lake%depth(ilake), lake%lake_energy(ilake) )

  enddo
  !$omp end do nowait

  !-----------------------------------------------------------------------
  ! Loop over ALL LAND CELLS
  !-----------------------------------------------------------------------

  !$omp do
  do iland = 2, mland

     call apply_land_fluxes( iland, watflux, energyflux,                          &
                             head_slope      (:,iland), land%soil_water(:,iland), &
                             land%soil_energy(:,iland), land%head      (:,iland), &
                             soil_watfrac    (:,iland), land%wsat_vg   (:,iland), &
                             land%wresid_vg  (:,iland)                            )

  enddo
  !$omp end do nowait
  !$omp end parallel

end subroutine gndwater_transport

!=========================================================================

subroutine comp_horiz_grndwatfluxes( ivsfc, soil_watfrac, watflux, energyflux )

  use mem_land,    only: land, omland, mland, nzg, slzt, dslz
  use mem_lake,    only: lake, omlake, timescale_lake
  use mem_sea,     only: sea, omsea
  use leaf_coms,   only: z_root
  use mem_sfcg,    only: itab_vsfc, sfcg
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use consts_coms, only: cliq1000, alli1000
  use sea_swm,     only: depthmin_flux
  use therm_lib,   only: qwtk, qtk
  use leaf4_soil,  only: soil_wat2khyd

  implicit none

  integer, intent(in)    :: ivsfc
  real,    intent(inout) :: watflux     (nzg)
  real,    intent(inout) :: energyflux  (nzg)
  real,    intent(in)    :: soil_watfrac(nzg,mland)

  ! Local arrays

  real :: soil_tempk       (nzg) ! soil temperature [K]
  real :: soil_fracliq     (nzg) ! fraction of soil moisture in liquid phase
  real :: hydresist1       (nzg)
  real :: hydresist2       (nzg)

  ! Local variables

  integer :: iw1, iw2, iland1, iland2, ilake1, ilake2, isea1, isea2
  integer :: k

  real :: dnvo2, soil_watfracv
  real :: khyd1, khyd2
  real :: laketemp, fracliq

  iw1 = itab_vsfc(ivsfc)%iwn(1)
  iw2 = itab_vsfc(ivsfc)%iwn(2)

  ! Skip this vsfc face if neighbor WSFC cell is absent
  if (iw1 < 2 .or. iw2 < 2) return

  ! Skip this vsfc face if running in parallel and cell rank is not MYRANK
  if (iparallel == 1 .and. itab_vsfc(ivsfc)%irank /= myrank) return

  ! Compute half distance along topographic slope between adjacent sfcg w points
  dnvo2 = 0.5 * sqrt(sfcg%dnv(ivsfc)**2 + (sfcg%topw(iw1) - sfcg%topw(iw2))**2)

  if (sfcg%leaf_class(iw1) >= 2 .and. sfcg%leaf_class(iw2) >= 2) then

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !LAND-LAND
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     ! Fluxes at current IVSFC location are between two land cells

     iland1 = iw1 - omland
     iland2 = iw2 - omland

     ! Loop vertically over soil levels
     do k = 1,nzg

        ! Compute soil water fraction at ivsfc and hydraulic resistance on
        ! iw1 and iw2 sides of ivsfc face

        soil_watfracv = .5 * (soil_watfrac(k,iland1) + soil_watfrac(k,iland2))

        call soil_wat2khyd(soil_watfracv, land%ksat_vg(k,iland1), land%lambda_vg(k,iland1), &
                           land%en_vg(k,iland1), khyd1, land%gpp(iland1), slzt(k),          &
                           z_root(sfcg%leaf_class(iw1)), land%head(k,iland1))

        call soil_wat2khyd(soil_watfracv, land%ksat_vg(k,iland2), land%lambda_vg(k,iland2), &
                           land%en_vg(k,iland2), khyd2, land%gpp(iland2), slzt(k),          &
                           z_root(sfcg%leaf_class(iw2)), land%head(k,iland2))

        ! Compute hydraulic resistances between iw1 and iw2 cells

        hydresist1(k) = dnvo2 / khyd1
        hydresist2(k) = dnvo2 / khyd2

        ! Compute horizontal mass and energy fluxes across this cell face,
        ! taking into account the topographic gradient since it is not
        ! included in head(:,:).  Modulate fluxes if water is partially or
        ! completely frozen.

        watflux(k) = (land%head(k,iland1) + sfcg%topw(iw1) - land%head(k,iland2) - sfcg%topw(iw2)) &
                   / (hydresist1(k) + hydresist2(k)) ! [m/s]

        if (watflux(k) > 0.) then

           call qwtk(land%soil_energy(k,iland1),land%soil_water(k,iland1)*1.e3, &
                     land%specifheat_drysoil(k,iland1),soil_tempk(k),soil_fracliq(k))

           watflux(k) = watflux(k) * soil_fracliq(k)

           energyflux(k) = watflux(k) &
                         * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

        else

           call qwtk(land%soil_energy(k,iland2),land%soil_water(k,iland2)*1.e3, &
                     land%specifheat_drysoil(k,iland2),soil_tempk(k),soil_fracliq(k))

           watflux(k) = watflux(k) * soil_fracliq(k)

           energyflux(k) = watflux(k) &
                         * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]
        endif

     enddo

  elseif (sfcg%leaf_class(iw1) >= 2) then

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !LAND-SEA or LAND-LAKE
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     ! Fluxes at current IVSFC location are between a land cell (iw1) and
     ! a water cell (iw2) - either sea or lake.  Hydraulic resistance is
     ! from land side only with assumed saturation at interface.  Head is
     ! the same at all levels in water cells (relative to surface elevation).

     iland1 = iw1 - omland

     ! Loop over vertical levels

     do k = 1,nzg

        ! Compute soil water fraction at ivsfc and hydraulic resistance on
        ! iw1 side of ivsfc face

        soil_watfracv = .5 * (soil_watfrac(k,iland1) + 1.0)

        ! Compute hydraulic resistances between iw1 and iw2 cells

        call soil_wat2khyd(soil_watfracv, land%ksat_vg(k,iland1), land%lambda_vg(k,iland1), &
                           land%en_vg(k,iland1), khyd1, land%gpp(iland1), slzt(k),          &
                           z_root(sfcg%leaf_class(iw1)), land%head(k,iland1))

        ! Compute hydraulic resistances between iw1 and iw2 cells

        hydresist1(k) = dnvo2 / khyd1

        ! Compute horizontal mass and energy fluxes across this cell face,
        ! taking into account the topographic gradient since it is not
        ! included in head(:,:).  Modulate fluxes if water is partially
        ! or completely frozen.

        watflux(k) = (land%head(k,iland1) + sfcg%topw(iw1) - sfcg%head1(iw2) - sfcg%topw(iw2)) &
                   / hydresist1(k) ! [m/s]

        if (watflux(k) > 0.) then

           call qwtk(land%soil_energy(k,iland1),land%soil_water(k,iland1)*1.e3, &
                     land%specifheat_drysoil(k,iland1),soil_tempk(k),soil_fracliq(k))

           watflux(k) = watflux(k) * soil_fracliq(k)

           energyflux(k) = watflux(k) &
                         * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

        else

           if (sfcg%leaf_class(iw2) == 0) then
              isea2 = iw2 - omsea
              energyflux(k) = watflux(k) &
                            * (cliq1000 * (sea%seatc(isea2) - 273.15) + alli1000) ! [J/(m^2 s)]
           else

              ilake2 = iw2 - omlake

              ! Zero flux out of lake if depth is below minimum
              if (lake%depth(ilake2) < depthmin_flux) then
                 watflux(k) = 0.
              endif

              call qtk(lake%lake_energy(ilake2), laketemp, fracliq)

              energyflux(k) = watflux(k) &
                            * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
           endif

        endif

     enddo ! k

  elseif (sfcg%leaf_class(iw2) >= 2) then

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !SEA-LAND or LAKE-LAND
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     ! Fluxes at current IVSFC location are between a land cell (iw2) and
     ! a water cell (iw1) - either sea or lake.  Hydraulic resistance is
     ! from land side only with assumed saturation at interface.  Head is
     ! the same at all levels in water cells (relative to surface elevation).

     iland2 = iw2 - omland

     ! Loop over vertical levels

     do k = 1,nzg

        ! Compute soil water fraction at ivsfc and hydraulic resistance on
        ! iw2 side of ivsfc face

        soil_watfracv = .5 * (soil_watfrac(k,iland2) + 1.0)

        call soil_wat2khyd(soil_watfracv, land%ksat_vg(k,iland2), land%lambda_vg(k,iland2), &
                           land%en_vg(k,iland2), khyd2, land%gpp(iland2), slzt(k),          &
                           z_root(sfcg%leaf_class(iw2)), land%head(k,iland2))

        ! Compute hydraulic resistances between iw1 and iw2 cells

        hydresist2(k) = dnvo2 / khyd2

        ! Compute horizontal mass and energy fluxes across this cell face,
        ! taking into account the topographic gradient since it is not
        ! included in head(:,:).  Modulate fluxes if water is partially
        ! or completely frozen.

        watflux(k) = (sfcg%head1(iw1) + sfcg%topw(iw1) - land%head(k,iland2) - sfcg%topw(iw2)) &
                   / hydresist2(k) ! [m/s]

        if (watflux(k) > 0.) then

           if (sfcg%leaf_class(iw1) == 0) then
              isea1 = iw1 - omsea
              energyflux(k) = watflux(k) &
                            * (cliq1000 * (sea%seatc(isea1) - 273.15) + alli1000) ! [J/(m^2 s)]
           else
              ilake1 = iw1 - omlake

              ! Zero flux out of lake if depth is below minimum
              if (lake%depth(ilake1) < depthmin_flux) then
                 watflux(k) = 0.
              endif

              call qtk(lake%lake_energy(ilake1), laketemp, fracliq)

              energyflux(k) = watflux(k) &
                            * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
           endif

        else

           call qwtk(land%soil_energy(k,iland2),land%soil_water(k,iland2)*1.e3, &
                     land%specifheat_drysoil(k,iland2),soil_tempk(k),soil_fracliq(k))

           watflux(k) = watflux(k) * soil_fracliq(k)

           energyflux(k) = watflux(k) &
                         * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

        endif

     enddo ! k

  elseif (sfcg%leaf_class(iw1) == 1 .and. sfcg%leaf_class(iw2) == 1) then

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !LAKE-LAKE
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     ! Fluxes at current IVSFC location are between two lake cells.  If they
     ! have different topographic heights, assume that they are unconnected
     ! and skip this IVN point.

     if (abs(sfcg%topw(iw1) - sfcg%topw(iw2)) > 0.1) return

     ilake1 = iw1 - omlake
     ilake2 = iw2 - omlake

     ! Given their (nearly) identical topographic heights, the two lake
     ! cells are assumed to be interconnected.  Compute a water flux at
     ! ivsfc that is proportional to the difference in head1 values between
     ! the lake cells so that water levels tend toward equilibration.
     ! ALTHOUGH MASS CONSERVING, THIS IS STRICTLY A RELAXATION TERM, NOT A
     ! COMPUTATION OF DYNAMICS BASED ON FORCE, DRAG, AND/OR MOMENTUM.

     watflux(nzg) = (sfcg%head1(iw1) + sfcg%topw(iw1) - sfcg%head1(iw2) - sfcg%topw(iw2)) &
                  * min(sfcg%area(iw1),sfcg%area(iw2)) &
                  / (sfcg%dnu(ivsfc) * dslz(nzg) * timescale_lake)  ! [m/s]

     if (watflux(nzg) > 0.) then

        ! Zero flux out of lake cell if depth is below minimum
        if (lake%depth(ilake1) < depthmin_flux) then
           watflux(nzg) = 0.
        endif

        call qtk(lake%lake_energy(ilake1), laketemp, fracliq)

        energyflux(nzg) = watflux(nzg) &
                        * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
     else

        ! Zero flux out of lake cell if depth is below minimum
        if (lake%depth(ilake2) < depthmin_flux) then
           watflux(nzg) = 0.
        endif

        call qtk(lake%lake_energy(ilake2), laketemp, fracliq)

        energyflux(nzg) = watflux(nzg) &
                        * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]

     endif

  endif ! sfcg%leaf_class(iw1 & iw2)

end subroutine comp_horiz_grndwatfluxes

!=============================================================================

subroutine apply_lake_fluxes( ilake, watflux,    energyflux, &
                                     lake_depth, lake_energy )

  use mem_land,  only: nzg, dslz
  use mem_lake,  only: omlake
  use mem_sfcg,  only: mvsfc, itab_wsfc, itab_vsfc, sfcg
  use mem_para,  only: myrank
  use misc_coms, only: iparallel
  use leaf_coms, only: wcap_min, dt_leaf

  implicit none

  integer, intent(in)    :: ilake
  real,    intent(in)    :: watflux   (nzg,mvsfc)
  real,    intent(in)    :: energyflux(nzg,mvsfc)

  real,    intent(inout) :: lake_depth
  real,    intent(inout) :: lake_energy

  integer :: iwsfc, j, ivn, iwn, k
  real    :: dheight
  real    :: energyin
  real    :: energy_per_m2
  real    :: dtoa
  real    :: dirv
  real    :: fconv
  real    :: factor

  iwsfc = ilake + omlake

  ! Skip this cell if running in parallel and cell rank is not MYRANK
  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  dheight  = 0.
  energyin = 0.

  dtoa = dt_leaf / sfcg%area(iwsfc) ! [s/m^s]

  ! Loop over all lateral faces of this ilake cell

  do j = 1,itab_wsfc(iwsfc)%npoly
     ivn = itab_wsfc(iwsfc)%ivn(j)
     iwn = itab_wsfc(iwsfc)%iwn(j)

     if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
        dirv = 1.0
     else
        dirv = -1.0
     endif

     fconv = dirv * dtoa * sfcg%dnu(ivn) ! [s/m]

     if (sfcg%leaf_class(iwn) >= 2) then

        ! Neighbor cell is land; loop over its vertical levels and sum
        ! mass and energy fluxes.

        do k = 1,nzg
           factor   = fconv * dslz(k)                       ! [s]
           dheight  = dheight  + factor * watflux   (k,ivn) ! [m]
           energyin = energyin + factor * energyflux(k,ivn) ! [J/m^2]
        enddo

     elseif (sfcg%leaf_class(iwn) == 1) then

        ! Neighbor cell is lake; only top-level flux is nonzero.

        factor   = fconv * dslz(k)                         ! [s]
        dheight  = dheight  + factor * watflux   (nzg,ivn) ! [m]
        energyin = energyin + factor * energyflux(nzg,ivn) ! [J/m^2]

     endif

  enddo ! j

  ! Apply height and energy changes to cell

  energy_per_m2 = lake_energy * lake_depth * 1000.

  lake_depth  = lake_depth + dheight

  lake_energy = (energy_per_m2 + energyin) &
                          / (max(wcap_min, lake_depth) * 1000.) ! water density = 1000 kg/m^3

end subroutine apply_lake_fluxes

!=============================================================================

subroutine apply_land_fluxes( iland, watflux,      energyflux,  head_slope, &
                                     soil_water,   soil_energy, head,       &
                                     soil_watfrac, wsat_vg,     wresid_vg   )

  use mem_land,  only: nzg, omland
  use mem_sfcg,  only: mvsfc, itab_wsfc, itab_vsfc, sfcg
  use mem_para,  only: myrank
  use misc_coms, only: iparallel
  use leaf_coms, only: dt_leaf

  implicit none

  integer, intent(in)    :: iland
  real,    intent(in)    :: watflux     (nzg,mvsfc)
  real,    intent(in)    :: energyflux  (nzg,mvsfc)
  real,    intent(in)    :: head_slope  (nzg)
  real,    intent(in)    :: wsat_vg     (nzg)
  real,    intent(in)    :: wresid_vg   (nzg)

  real,    intent(inout) :: soil_water  (nzg)
  real,    intent(inout) :: soil_energy (nzg)
  real,    intent(inout) :: head        (nzg)
  real,    intent(inout) :: soil_watfrac(nzg)

  real    :: dheight (nzg)
  real    :: energyin(nzg)

  real    :: dtoa, dirv, factor
  integer :: iwsfc, j, ivn, k

  iwsfc = iland + omland

  ! Skip this cell if running in parallel and cell rank is not MYRANK
  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  ! Apply horizontal groundwater mass and energy fluxes

  dheight (:) = 0.
  energyin(:) = 0.

  dtoa = dt_leaf / sfcg%area(iwsfc) ! [s/m^s]

  ! Loop over all lateral faces of this iland cell

  do j = 1,itab_wsfc(iwsfc)%npoly
     ivn = itab_wsfc(iwsfc)%ivn(j)

     if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
        dirv = 1.0
     else
        dirv = -1.0
     endif

     factor = dirv * dtoa * sfcg%dnu(ivn) ! [s/m]

     ! Loop over vertical levels; sum mass and energy fluxes over j faces

     do k = 1, nzg
        dheight (k) = dheight (k) + factor * watflux   (k,ivn) ! [Vol/Vol]
        energyin(k) = energyin(k) + factor * energyflux(k,ivn) ! [J/m^3]
     enddo

  enddo

  ! Apply mass and energy changes to cell, and adjust hydraulic head due to
  ! water mass change according to constant slope fit for current timestep

  do k = 1, nzg
     factor =  wsat_vg(k) - wresid_vg(k)

     soil_water  (k) = soil_water  (k) + dheight (k)
     soil_energy (k) = soil_energy (k) + energyin(k)
     head        (k) = head        (k) + dheight (k) * head_slope(k)
     soil_watfrac(k) = soil_watfrac(k) + dheight (k) / factor

     soil_watfrac(k) = max(0., min(1., soil_watfrac(k)))
  enddo

end subroutine apply_land_fluxes

!=============================================================================
