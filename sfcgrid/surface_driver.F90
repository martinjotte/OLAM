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

subroutine surface_driver()

  use leaf_coms,     only: nzs, iupdndvi, s1900_ndvi, dt_leaf, &
                           indvifile, nndvifiles, soil_rough, snow_rough, z_root
  use mem_ijtabs,    only: itabg_w, itab_w
  use mem_sfcg,      only: itab_wsfc, itab_vsfc, sfcg, mvsfc, mwsfc, &
                           jtab_wsfc_swm
  use mem_land,      only: land, mland, omland, nzg, slz, slzt, dslz
  use mem_lake,      only: lake, mlake, omlake, timescale_lake
  use mem_sea,       only: sea, msea, omsea
  use misc_coms,     only: io6, s1900_sim, isubdomain, time8p, iparallel
  use mem_basic,     only: rho, press, theta, tair, vxe, vye, vze, rr_v
  use mem_micro,     only: rr_c
  use consts_coms,   only: grav, cliq1000, alli1000, erad, eradi
  use leaf4_canopy,  only: canopy, vegndvi, fast_canopy
  use leaf4_surface, only: sfcwater, sfcwater_adjust, remove_runoff, grndvap
  use leaf4_soil,    only: soil, soil_wat2khyd, head_column
  use leaf4_plot,    only: leaf_plot
  use leaf_coms,     only: wcap_min
  use sea_swm,       only: swm_grad2d, swm_hflux, swm_progw, swm_progv, &
                           swm_diagvel, depthmin_flux
  use therm_lib,     only: qwtk, qtk
  use mem_para,      only: myrank
  use olam_mpi_sfc,  only: mpi_send_wsfc, mpi_recv_wsfc, mpi_send_vsfc, mpi_recv_vsfc
  use oname_coms,    only: nl
  use mem_sfcnud,    only: sfcwat_nud, sfctemp_nud, fracliq_nud

  implicit none

  ! Local arrays

  real :: sfcwater_tempk   (nzs) ! surface water temperature [K]
  real :: sfcwater_fracliq (nzs) ! fraction of sfc water in liquid phase
  real :: soil_tempk       (nzg) ! soil temperature [K]
  real :: soil_fracliq     (nzg) ! fraction of soil moisture in liquid phase
  real :: energy_per_m2    (nzs) ! sfcwater or lake energy [J/m^2]
  real :: hydresist1       (nzg)
  real :: hydresist2       (nzg)
  real :: dheight          (nzg) ! change in water height (of a T cell)
                                 ! from lateral water fluxes [m]
  real :: energyin         (nzg) ! change in energy (of T cell) from lateral water fluxes [J/m^2]

  real :: head_slope  (nzg,mland)
  real :: soil_watfrac(nzg,mland)
  real :: watflux     (nzg,mvsfc) ! [m/s]
  real :: energyflux  (nzg,mvsfc)

  ! Local variables

  integer :: iland, ilandn, j, iwsfc, ilake, isea
  integer :: ivn, iwn, ivsfc
  integer :: iw1, iw2, iland1, iland2, ilake1, ilake2, isea1, isea2

  integer :: k     ! vertical index over soil layers
  integer :: nlsw1 ! maximum of (1,land%nlev_sfcwater(iland))
  integer :: icomb ! implicit heat balance flag [0=no, 1=yes]

  integer :: ktrans ! vertical index of soil layer supplying transpiration

  real :: transp  ! transpiration xfer this LEAF timestep [kg/m^2]
  real :: wfree1  ! sfcwater bottom layer free water mass [kg/m^2]
  real :: qwfree1 ! sfcwater bottom layer free water energy [J/m^2]
  real :: dwfree1 ! sfcwater bottom layer free water depth [m]

  real :: timefac_ndvi  ! fraction of elapsed time from past to future NDVI obs
  real :: dnvo2, dirv, soil_watfracv

  real :: khyd1, khyd2
  real :: factor

  real, parameter :: thermcond_dry_organic = 0.05  ! [W/(m K)] 
  real, parameter :: thermcond_sat_organic = 0.25  ! [W/(m K)] 
  real, parameter :: thermcond_bedrock     = 3.0   ! [W/(m K)] 
  real, parameter :: thermcond_liq         = 0.6   ! [W/(m K)] 
  real, parameter :: thermcond_ice         = 2.29  ! [W/(m K)] 

  real :: thermcond_dry_mineral
  real :: thermcond_sat_mineral

  real :: mineral
  real :: bulkdens_mineral
  real :: thermcond_drysoil
  real :: thermcond_sat_solids
  real :: thermcond_sat_soil
  real :: kersten_liq, kersten_ice
  real :: kersten
  real :: thermcond_soil(nzg)

  integer, parameter :: iland_print = 0

  real, external :: rhovsil ! function to compute sat spec hum over liquid or ice

  real :: laketemp, fracliq

  watflux   (:,:) = 0.
  energyflux(:,:) = 0.

  ! Loop over all land cells (in this subdomain)

  !$omp parallel do private(iwsfc)
  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Compute total head (relative to local topographic datum) based on soil
     ! moisture, pressure head, and depth in soil

     call head_column(nzg, iland,                &
                      slzt            (:),       &
                      land%soil_water (:,iland), &
                      land%wresid_vg  (:,iland), &
                      land%wsat_vg    (:,iland), &
                      land%alpha_vg   (:,iland), &
                      land%en_vg      (:,iland), &
                      land%head       (:,iland), &
                      head_slope      (:,iland), &
                      soil_watfrac    (:,iland)  )

  enddo ! iland
  !$omp end parallel do

  ! Evaluate horizontal gradients of temp, VXE, VYE, and VZE for ocean cells
  ! that are updated with Shallow Water Model (SWM).

  if (nl%igw_spinup /= 1) call swm_grad2d()

  if (iparallel == 1) then

     ! MPI SEND/RECV land%head, soil_watfrac, sea%gxps_vxe, sea%gyps_vxe,
     !               sea%gxps_vye, sea%gyps_vye, sea%gxps_vze, sea%gyps_vze

     call mpi_send_wsfc(set='head_swm_grad',soil_watfrac=soil_watfrac)
     call mpi_recv_wsfc(set='head_swm_grad',soil_watfrac=soil_watfrac)

  endif

  ! Loop over all SFC grid V faces in this subdomain, except where neighbor
  ! WSFC cell is absent.  Evaluate water and energy fluxes across V face, with
  ! flux direction defined by flux = (head(v%iw1) - head(v%iw2)) / resist

  !$omp parallel private(soil_tempk, soil_fracliq, hydresist1, hydresist2)
  !$omp do private(iw1, iw2, iland1, iland2, k, ilake1, ilake2, isea1, isea2, &
  !$omp            dnvo2, soil_watfracv, khyd1, khyd2)
  do ivsfc = 2,mvsfc
     iw1 = itab_vsfc(ivsfc)%iwn(1)
     iw2 = itab_vsfc(ivsfc)%iwn(2)

     if (iw1 < 2 .or. iw2 < 2) cycle

     ! Skip this vsfc face if running in parallel and cell rank is not MYRANK

     if (isubdomain == 1 .and. itab_vsfc(ivsfc)%irank /= myrank) cycle

     ! Compute half distance along topographic slope between adjacent sfcg w points

     dnvo2 = 0.5 * sqrt(sfcg%dnv(ivsfc)**2 + (sfcg%topw(iw1) - sfcg%topw(iw2))**2)

     if (sfcg%leaf_class(iw1) >= 2 .and. sfcg%leaf_class(iw2) >= 2) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!LAND-LAND
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

           watflux(k,ivsfc) = (land%head(k,iland1) + sfcg%topw(iw1) - land%head(k,iland2) - sfcg%topw(iw2)) &
                            / (hydresist1(k) + hydresist2(k)) ! [m/s]

           if (watflux(k,ivsfc) > 0.) then

              call qwtk(land%soil_energy(k,iland1),land%soil_water(k,iland1)*1.e3, &
                        land%specifheat_drysoil(k,iland1),soil_tempk(k),soil_fracliq(k))

              watflux(k,ivsfc) = watflux(k,ivsfc) * soil_fracliq(k)

              energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                  * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

           else

              call qwtk(land%soil_energy(k,iland2),land%soil_water(k,iland2)*1.e3, &
                        land%specifheat_drysoil(k,iland2),soil_tempk(k),soil_fracliq(k))

              watflux(k,ivsfc) = watflux(k,ivsfc) * soil_fracliq(k)

              energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                  * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]
           endif

        enddo

     elseif (sfcg%leaf_class(iw1) >= 2) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!LAND-SEA or LAND-LAKE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

           watflux(k,ivsfc) = (land%head(k,iland1) + sfcg%topw(iw1) - sfcg%head1(iw2) - sfcg%topw(iw2)) &
                            / hydresist1(k) ! [m/s]

           if (watflux(k,ivsfc) > 0.) then

              call qwtk(land%soil_energy(k,iland1),land%soil_water(k,iland1)*1.e3, &
                        land%specifheat_drysoil(k,iland1),soil_tempk(k),soil_fracliq(k))

              watflux(k,ivsfc) = watflux(k,ivsfc) * soil_fracliq(k)

              energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                  * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

           else

              if (sfcg%leaf_class(iw2) == 0) then
                 isea2 = iw2 - omsea
                 energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                     * (cliq1000 * (sea%seatc(isea2) - 273.15) + alli1000) ! [J/(m^2 s)]
              else

                 ilake2 = iw2 - omlake

                 ! Zero flux out of lake if depth is below minimum
                 if (lake%depth(ilake2) < depthmin_flux) then
                    watflux(k,ivsfc) = 0.
                 endif

                 call qtk(lake%lake_energy(ilake2), laketemp, fracliq)

                 energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                     * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
              endif

           endif

        enddo ! k

     elseif (sfcg%leaf_class(iw2) >= 2) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!SEA-LAND or LAKE-LAND
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

           watflux(k,ivsfc) = (sfcg%head1(iw1) + sfcg%topw(iw1) - land%head(k,iland2) - sfcg%topw(iw2)) &
                            / hydresist2(k) ! [m/s]

           if (watflux(k,ivsfc) > 0.) then

              if (sfcg%leaf_class(iw1) == 0) then
                 isea1 = iw1 - omsea
                 energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                     * (cliq1000 * (sea%seatc(isea1) - 273.15) + alli1000) ! [J/(m^2 s)]
              else
                 ilake1 = iw1 - omlake

                 ! Zero flux out of lake if depth is below minimum
                 if (lake%depth(ilake1) < depthmin_flux) then
                    watflux(k,ivsfc) = 0.
                 endif

                 call qtk(lake%lake_energy(ilake1), laketemp, fracliq)

                 energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                     * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
              endif

           else

              call qwtk(land%soil_energy(k,iland2),land%soil_water(k,iland2)*1.e3, &
                        land%specifheat_drysoil(k,iland2),soil_tempk(k),soil_fracliq(k))

              watflux(k,ivsfc) = watflux(k,ivsfc) * soil_fracliq(k)

              energyflux(k,ivsfc) = watflux(k,ivsfc) &
                                  * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000) ! [J/(m^2 s)]

           endif

        enddo ! k

     elseif (sfcg%leaf_class(iw1) == 1 .and. sfcg%leaf_class(iw2) == 1) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!LAKE-LAKE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! Fluxes at current IVSFC location are between two lake cells.  If they
        ! have different topographic heights, assume that they are unconnected
        ! and skip this IVN point.

        if (abs(sfcg%topw(iw1) - sfcg%topw(iw2)) > 0.1) cycle

        ilake1 = iw1 - omlake
        ilake2 = iw2 - omlake

        ! Given their (nearly) identical topographic heights, the two lake
        ! cells are assumed to be interconnected.  Compute a water flux at
        ! ivsfc that is proportional to the difference in head1 values between
        ! the lake cells so that water levels tend toward equilibration.
        ! ALTHOUGH MASS CONSERVING, THIS IS STRICTLY A RELAXATION TERM, NOT A
        ! COMPUTATION OF DYNAMICS BASED ON FORCE, DRAG, AND/OR MOMENTUM.

        watflux(nzg,ivsfc) = (sfcg%head1(iw1) + sfcg%topw(iw1) - sfcg%head1(iw2) - sfcg%topw(iw2)) &
                           * min(sfcg%area(iw1),sfcg%area(iw2)) &
                           / (sfcg%dnu(ivsfc) * dslz(nzg) * timescale_lake)  ! [m/s]

        if (watflux(nzg,ivsfc) > 0.) then

           ! Zero flux out of lake cell if depth is below minimum
           if (lake%depth(ilake1) < depthmin_flux) then
              watflux(nzg,ivsfc) = 0.
           endif

           call qtk(lake%lake_energy(ilake1), laketemp, fracliq)

           energyflux(nzg,ivsfc) = watflux(nzg,ivsfc) &
                                 * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]
        else

           ! Zero flux out of lake cell if depth is below minimum
           if (lake%depth(ilake2) < depthmin_flux) then
              watflux(nzg,ivsfc) = 0.
           endif

           call qtk(lake%lake_energy(ilake2), laketemp, fracliq)

           energyflux(nzg,ivsfc) = watflux(nzg,ivsfc) &
                                 * (cliq1000 * (laketemp - 273.15) + alli1000) ! [J/(m^2 s)]

        endif

     endif ! sfcg%leaf_class(iw1 & iw2)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!LAND-LAND, LAND-LAKE, LAND-SEA, LAKE-LAND, SEA-LAND, SEA-SEA surface fluxes in SWM-ACTIVE regions
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     if (sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) then ! Includes sea, lake, land cells
        if (nl%igw_spinup /= 1) call swm_hflux(ivsfc, iw1, iw2)
     endif

  enddo ! ivsfc
  !$omp end do
  !$omp end parallel

  if (iparallel == 1) then

     ! MPI SEND/RECV watflux, energyflux, sfcg%hflux_wat, sfcg%hflux_enr,
     !               sfcg%hflux_vxe, sfcg%hflux_vye, sfcg%hflux_vze,
     !               sfcg%vmp, sfcg%vmc, sfcg%vc

     call mpi_send_vsfc(watflux=watflux,energyflux=energyflux)
     call mpi_recv_vsfc(watflux=watflux,energyflux=energyflux)

  endif

  ! Time interpolation factor for updating NDVI

  timefac_ndvi = 0.

  if (iupdndvi == 1 .and. nndvifiles > 1) then
     timefac_ndvi = (s1900_sim             - s1900_ndvi(indvifile-1))  &
                  / (s1900_ndvi(indvifile) - s1900_ndvi(indvifile-1))
  endif

  !----------------------------------------------------------------------------------
  ! Loop over SEA CELLS that use SWM; prognose WDEPTH, HEAD1, & W-cell horiz momentum
  !----------------------------------------------------------------------------------

  if (nl%igw_spinup /= 1) then

     !$omp parallel do private(iwsfc) 
     do j = 1,jtab_wsfc_swm%jend; iwsfc = jtab_wsfc_swm%iwsfc(j)
        call swm_progw(iwsfc)
     enddo
     !$omp end parallel do 

     ! MPI SEND/RECV of quantities needed for prog_v

     if (iparallel == 1) then
        call mpi_send_wsfc(set='swm_progw')
        call mpi_recv_wsfc(set='swm_progw')
     endif

     ! UPDATE SFCG%VMC AND SFCG%VC for sea cells that use SWM

     call swm_progv()

     ! MPI send/recv of SFCG%VMC, SFCG%VC

     if (iparallel == 1) then
        call mpi_send_vsfc()
        call mpi_recv_vsfc()
     endif

     ! Compute earth cartesian velocities (VXE, VYE, VZE) for sea cells that use SWM

     if (nl%igw_spinup /= 1) call swm_diagvel()

     ! Parallel send-recieve of SWM Earth Cartesian velocities

     if (iparallel == 1) then
        call mpi_send_wsfc(set='swm_diagvel')
        call mpi_recv_wsfc(set='swm_diagvel')
     endif

  endif

  ! Update surface temperature, ice, and canopy variables in all sea cells.
  ! Update POM1D column variables in sea cells that are pom_active and not swm_active.

  call seacells()

  ! NOTE: If the call to plot_pom gets uncommented, subroutine plot_pom needs
  ! to be customized so that the 4 IWSFC/ISEA points that it plots are actually
  ! in the pom_active region of the simulation that you are running.

  ! print*, 'calling plot_pom'

  ! call plot_pom()

  ! print*, 'returned from plot_pom'

  !-----------------------------------------------------------------------
  ! Loop over ALL LAKE CELLS
  !-----------------------------------------------------------------------

  !$omp parallel private(dheight, energyin)
  !$omp do private(iwsfc, j, ivn, iwn, dirv, k, factor)
  do ilake = 2,mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     dheight (:) = 0.
     energyin(:) = 0.

     ! Loop over all lateral faces of this ilake cell

     do j = 1,itab_wsfc(iwsfc)%npoly
        ivn = itab_wsfc(iwsfc)%ivn(j)
        iwn = itab_wsfc(iwsfc)%iwn(j)

        if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
           dirv = 1.0
        else
           dirv = -1.0
        endif

        if (sfcg%leaf_class(iwn) >= 2) then

           ! Neighbor cell is land; loop over its vertical levels and sum
           ! mass and energy fluxes.  Accumulate sum only in nzg level.

           do k = 1,nzg
              factor = dirv * dt_leaf * sfcg%dnu(ivn) * dslz(k) / sfcg%area(iwsfc) ! [s] to here

              dheight (nzg) = dheight (nzg) + factor * watflux   (k,ivn) ! [m]
              energyin(nzg) = energyin(nzg) + factor * energyflux(k,ivn) ! [J/m^2]

           enddo

           if (nl%igw_spinup /= 1) then

              ! Add mass and energy contributions from soil sfcwater in neighbor cell

              factor = dirv * dt_leaf / sfcg%area(iwsfc) ! [s/m^2]

              dheight (nzg) = dheight (nzg) + factor * sfcg%hflux_wat(ivn) ! [m] or [m^3/m^2]
              energyin(nzg) = energyin(nzg) + factor * sfcg%hflux_enr(ivn) ! [J/m^2]

           endif

        elseif (sfcg%leaf_class(iwn) == 1) then

           ! Neighbor cell is lake; only top-level flux is nonzero.

           factor = dirv * dt_leaf * sfcg%dnu(ivn) * dslz(nzg) / sfcg%area(iwsfc) ! [s/m]

           dheight (nzg) = dheight (nzg) + factor * watflux   (nzg,ivn) ! [m]
           energyin(nzg) = energyin(nzg) + factor * energyflux(nzg,ivn) ! [J/m^2]

        endif

     enddo ! j

     ! Add height and energy changes to cell

     energy_per_m2(1) = lake%lake_energy(ilake) * lake%depth(ilake) * 1000.

     lake%depth(ilake) = lake%depth(ilake) + dheight(nzg)

     lake%lake_energy(ilake) = (energy_per_m2(1) + energyin(nzg)) &
                             / (max(wcap_min,lake%depth(ilake)) * 1000.) ! water density = 1000 kg/m^3

     sfcg%head1(iwsfc) = lake%depth(ilake) + sfcg%bathym(iwsfc) - sfcg%topw(iwsfc)
  enddo
  !$omp end do
  !$omp end parallel

  call lakecells()

  !-----------------------------------------------------------------------
  ! Loop over ALL LAND CELLS
  !-----------------------------------------------------------------------

  !$omp parallel private(dheight, energyin, soil_tempk, soil_fracliq, &
  !$omp                  thermcond_soil, sfcwater_tempk, sfcwater_fracliq, &
  !$omp                  energy_per_m2)

  !$omp do private(iwsfc, j, ivn, dirv, k, factor, icomb, nlsw1, &
  !$omp            bulkdens_mineral, mineral, thermcond_dry_mineral, &
  !$omp            thermcond_drysoil, thermcond_sat_mineral, &
  !$omp            thermcond_sat_solids, thermcond_sat_soil, &
  !$omp            kersten_liq, kersten_ice, kersten, &
  !$omp            ktrans, transp, wfree1, qwfree1, dwfree1)
  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply horizontal groundwater mass and energy fluxes

     dheight (:) = 0.
     energyin(:) = 0.

     energy_per_m2(1) = land%sfcwater_energy(1,iland) * land%sfcwater_mass(1,iland)

     ! Loop over all lateral faces of this iland cell

     do j = 1,itab_wsfc(iwsfc)%npoly
        ivn = itab_wsfc(iwsfc)%ivn(j)

        if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
           dirv = 1.0
        else
           dirv = -1.0
        endif

        ! Loop over vertical levels; sum mass and energy fluxes over j faces

        do k = 1,nzg
           factor = dirv * dt_leaf * sfcg%dnu(ivn) * dslz(k) / sfcg%area(iwsfc) ! [s]

           dheight (k) = dheight (k) + factor * watflux   (k,ivn) ! [m]
           energyin(k) = energyin(k) + factor * energyflux(k,ivn) ! [J/m^2]
        enddo ! k

        if (nl%igw_spinup /= 1) then

           ! Add mass and energy contributions from overland flow from neighbor cell

           factor = dirv * dt_leaf / sfcg%area(iwsfc) ! [s/m^2]

           land%sfcwater_depth(1,iland)  = land%sfcwater_depth(1,iland)  &
                                         + factor * sfcg%hflux_wat(ivn) ! [m] or [m^3/m^2]

           land%sfcwater_mass(1,iland)   = land%sfcwater_mass(1,iland)   &
                                         + factor * sfcg%hflux_wat(ivn) * 1000. ! [kg/m^2]

           energy_per_m2(1)              = energy_per_m2(1) & 
                                         + factor * sfcg%hflux_enr(ivn) ! [J/m^2]

        endif

        land%sfcwater_energy(1,iland) = energy_per_m2(1) &
                                      / max(wcap_min,land%sfcwater_mass(1,iland))

     enddo ! j, ivn

     ! Reset relationship between sfcwater presence and nlev_sfcwater setting

     if (land%sfcwater_mass(1,iland) > wcap_min .and. land%nlev_sfcwater(iland) == 0) then
        land%nlev_sfcwater(iland) = 1
        land%sfcwater_energy(1,iland) = energy_per_m2(1) / land%sfcwater_mass(1,iland)
     elseif (land%sfcwater_mass(1,iland) < wcap_min .and. land%nlev_sfcwater(iland) > 0) then
        land%nlev_sfcwater(iland) = 0
        land%sfcwater_mass(1,iland) = 0.
        land%sfcwater_energy(1,iland) = 0.
     endif

     ! Add mass and energy changes to cell, and adjust hydraulic head due to
     ! water mass change according to constant slope fit for current timestep

     do k = 1,nzg
        land%soil_water (k,iland) = land%soil_water (k,iland) + dheight(k)
        land%soil_energy(k,iland) = land%soil_energy(k,iland) + energyin(k)

        land%head(k,iland) = land%head(k,iland) + head_slope(k,iland) * dheight(k)
     enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! START OF OLD LEAF4_LANDCELL SECTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     ! Initialize variables for this land cell

     icomb = 1  ! Default initialization
     nlsw1 = max(land%nlev_sfcwater(iland),1)

     ! Diagnose soil temperature, liquid fraction, and thermal conductivity

     do k = 1,nzg
        call qwtk(land%soil_energy(k,iland),land%soil_water(k,iland)*1.e3, &
           land%specifheat_drysoil(k,iland),soil_tempk(k),soil_fracliq(k))

        ! Check whether bedrock or soil (negative ph_soil has been set as a flag for bedrock)

        if (land%ph_soil(k,iland) < 0.) then

           thermcond_soil(k) = thermcond_bedrock  ! Assuming that porosity is accounted for

        else

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

           thermcond_sat_mineral = (8.80 * land%sand(k,iland) + 2.92 * land%clay(k,iland)) &
                                 / (land%sand(k,iland) + land%clay(k,iland))

           thermcond_sat_solids = mineral             * thermcond_sat_mineral &
                                + land%organ(k,iland) * thermcond_sat_organic 
           !--------------------------------------------------------------------------------

           ! For soil case, use porosity and fractions of sand, clay, and organic matter

           ! Compute thermcond_sat_soil as geometric mean of solids, liquid, and ice contributions

           thermcond_sat_soil = thermcond_sat_solids**(1.-land%wsat_vg(k,iland)) &
                              * thermcond_liq**(land%wsat_vg(k,iland)*soil_fracliq(k)) &
                              * thermcond_ice**(land%wsat_vg(k,iland)*(1.-soil_fracliq(k)))

           kersten_liq = max(0.,log10(soil_watfrac(k,iland)) + 1.)
           kersten_ice = soil_watfrac(k,iland)
           kersten = kersten_liq * soil_fracliq(k) + kersten_ice * (1. - soil_fracliq(k))

           thermcond_soil(k) = kersten * thermcond_sat_soil &
                       + (1. - kersten) * thermcond_drysoil

        endif

     enddo

     ! Loop over all possible sfcwater (a.k.a. snowcover) layers

     do k = 1,nzs

        ! Check (temporary) for consistency between nlev_sfcwater and 
        ! sfcwater_mass.  If inconsistency found, print message and stop.

        if (land%sfcwater_mass(k,iland) > wcap_min .and. k > land%nlev_sfcwater(iland)) then
           write(io6,*) 'sfcwater mass is present in layer above nlev_sfcwater'
           write(io6,*) 'iland,k,glatw,glonw = ',iland,k, sfcg%glatw(iwsfc), sfcg%glonw(iwsfc)
           write(io6,*) 'nlev_sfcwater,sfcwater_mass(k) = ', &
                        land%nlev_sfcwater(iland),land%sfcwater_mass(k,iland)
           stop 'stop1 leaf4'
        endif

        if (land%sfcwater_mass(k,iland) < wcap_min .and. k <= land%nlev_sfcwater(iland)) then
           write(io6,*) 'sfcwater mass is absent in layer at or below nlev_sfcwater'
           write(io6,*) 'iland,k,glatw,glonw =',iland,k, sfcg%glatw(iwsfc), sfcg%glonw(iwsfc)
           write(io6,*) 'nlev_sfcwater,sfcwater_mass(k) = ', &
                        land%nlev_sfcwater(iland),land%sfcwater_mass(k,iland)
           stop 'stop2 leaf4'
        endif

        if (k <= land%nlev_sfcwater(iland)) then

           ! Diagnose sfcwater temperature, liquid fraction, and energy per m^2

           call qtk(land%sfcwater_energy(k,iland),sfcwater_tempk(k),sfcwater_fracliq(k))

           energy_per_m2(k) = land%sfcwater_energy(k,iland) * land%sfcwater_mass(k,iland)

        else

           ! Assign default values above active sfcwater layers

           land%sfcwater_mass  (k,iland) = 0.
           land%sfcwater_energy(k,iland) = 0.
           land%sfcwater_depth (k,iland) = 0.

           sfcwater_tempk  (k) = 273.15
           sfcwater_fracliq(k) = 0.
           energy_per_m2   (k) = 0.

        endif
     enddo

     if (nl%igw_spinup == 1) then

        ! With igw_spinup = 1, call subroutine fast_canopy to nudge surface water
        ! mass and energy

        ktrans = nzg
        transp = 0.

        call fast_canopy(iland, iwsfc, nlsw1, icomb,           &
                         sfcwat_nud                   (iwsfc), &
                         sfctemp_nud                  (iwsfc), &
                         fracliq_nud                  (iwsfc), &
                         land%sfcwater_mass     (nlsw1,iland), &
                         land%sfcwater_energy   (nlsw1,iland), &
                         land%sfcwater_depth    (nlsw1,iland), &
                         energy_per_m2          (nlsw1),       &
                         sfcwater_tempk         (nlsw1),       &
                         sfcwater_fracliq       (nlsw1),       &
                         land%soil_water        (1:nzg,iland), &
                         land%soil_energy       (1:nzg,iland), &
                         land%specifheat_drysoil(1:nzg,iland), &
                         soil_tempk             (1:nzg),       &
                         soil_fracliq           (1:nzg)        )

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

        ! Evaluate turbulent exchanges of heat and moisture between vegetation and canopy air
        ! and also between soil or snow surface and canopy air.  Evaluate transfer of
        ! precipitation moisture and heat to vegetation and shed from vegetation to surface.
        ! Update vegetation and canopy air temperatures resulting from these 
        ! plus radiative fluxes.

        call canopy(iland, iwsfc, nlsw1, icomb, ktrans, transp, &
                    sfcg%leaf_class           (iwsfc), &
                    sfcg%can_depth            (iwsfc), &
                    sfcg%rhos                 (iwsfc), &
                    sfcg%vels                 (iwsfc), &
                    sfcg%ustar                (iwsfc), &
                    sfcg%sxfer_t              (iwsfc), &
                    sfcg%sxfer_r              (iwsfc), &
                    sfcg%pcpg                 (iwsfc), &
                    sfcg%qpcpg                (iwsfc), &
                    sfcg%dpcpg                (iwsfc), &
                    sfcg%rshort               (iwsfc), &
                    sfcg%cantemp              (iwsfc), &
                    sfcg%canrrv               (iwsfc), &
                    sfcg%glatw                (iwsfc), &
                    sfcg%glonw                (iwsfc), & 
                    land%snowfac              (iland), &
                    land%vf                   (iland), &
                    land%stom_resist          (iland), &
                    land%veg_height           (iland), &
                    land%veg_rough            (iland), &
                    land%veg_tai              (iland), &
                    land%veg_lai              (iland), &
                    land%hcapveg              (iland), &
                    land%rshort_v             (iland), &
                    land%rshort_g             (iland), &
                    land%rlong_v              (iland), &
                    land%rlong_s              (iland), &
                    land%rlong_g              (iland), &
                    land%veg_water            (iland), &
                    land%veg_energy           (iland), &
                    land%veg_temp             (iland), &   
                    land%sfcwater_mass  (nlsw1,iland), &
                    land%sfcwater_energy(nlsw1,iland), &
                    land%sfcwater_depth (nlsw1,iland), &
                    energy_per_m2       (nlsw1),       &
                    sfcwater_tempk      (nlsw1),       &
                    sfcwater_fracliq    (nlsw1),       &
                    land%soil_water        (1:nzg,iland), &
                    land%soil_energy       (1:nzg,iland), &
                    land%wsat_vg           (1:nzg,iland), &
                    land%ksat_vg           (1:nzg,iland), &
                    land%specifheat_drysoil(1:nzg,iland), &
                    land%head              (1:nzg,iland), &
                    head_slope             (1:nzg,iland), &
                    soil_tempk             (1:nzg),       &
                    soil_fracliq           (1:nzg)        )

     endif

     ! CALL SFCWATER:
     !  1. Compute soil and sfcwater heat conductivities
     !  2. Compute surface heat flux for top layer of soil or sfcwater
     !  3. Evaluate internal and bottom fluxes for sfcwater
     !  4. Update sfcwater layer energies due to heat flux and solar radiation
     !  5. Evaluate melting and percolation of liquid through sfcwater layers

     call sfcwater(iland, iwsfc, icomb, wfree1, qwfree1, dwfree1, &
                   sfcg%head1                   (iwsfc), &
                   land%nlev_sfcwater           (iland), &
                   land%sfcwater_mass     (1:nzs,iland), &
                   land%sfcwater_energy   (1:nzs,iland), &
                   land%sfcwater_depth    (1:nzs,iland), &
                   sfcwater_tempk         (1:nzs),       &
                   sfcwater_fracliq       (1:nzs),       &
                   energy_per_m2          (1:nzs),       &
                   land%soil_water        (1:nzg,iland), &
                   land%soil_energy       (1:nzg,iland), &
                   land%specifheat_drysoil(1:nzg,iland), &
                   soil_watfrac           (1:nzg,iland), &
                   soil_tempk             (1:nzg),       &
                   soil_fracliq           (1:nzg),       &
                   thermcond_soil         (1:nzg)        )

     call soil(iland, iwsfc, ktrans, transp, wfree1, qwfree1, &
               sfcg%leaf_class          (iwsfc), &
               sfcg%glatw               (iwsfc), &
               sfcg%glonw               (iwsfc), &
               sfcg%head1               (iwsfc), &
               land%head0               (iland), &
               land%nlev_sfcwater       (iland), &
               land%gpp                 (iland), &
               land%sfcwater_mass     (1,iland), &
               land%sfcwater_energy   (1,iland), &
               land%sfcwater_depth    (1,iland), &
               energy_per_m2          (1),       &
               land%soil_water    (1:nzg,iland), &
               land%soil_energy   (1:nzg,iland), &
               land%wresid_vg     (1:nzg,iland), &
               land%wsat_vg       (1:nzg,iland), &
               land%ksat_vg       (1:nzg,iland), &
               land%lambda_vg     (1:nzg,iland), &
               land%en_vg         (1:nzg,iland), &
               land%alpha_vg      (1:nzg,iland), &
               soil_watfrac       (1:nzg,iland), &
               head_slope         (1:nzg,iland), &
               land%head          (1:nzg,iland), &
               soil_tempk         (1:nzg),       &
               soil_fracliq       (1:nzg),       &
               dheight            (1:nzg),       & ! included here only to pass to leaf_plot
               energyin           (1:nzg),       & ! included here only to pass to leaf_plot
               thermcond_soil     (1:nzg)        )

     ! Inventory sfcwater layer(s) and adjust thicknesses if more than one

     if (nzs > 1) then
        call sfcwater_adjust(iland, iwsfc,                      &
                             sfcg%glatw                (iwsfc), &
                             sfcg%glonw                (iwsfc), &
                             land%nlev_sfcwater        (iland), &
                             land%sfcwater_mass  (1:nzs,iland), &
                             land%sfcwater_energy(1:nzs,iland), &
                             land%sfcwater_depth (1:nzs,iland), &
                             energy_per_m2       (1:nzs)        )
     endif

     ! Compute surface and ground vap mxrat for next timestep; put into surface_srrv 
     ! and ground_rrv.

     nlsw1 = max(land%nlev_sfcwater(iland),1)

     call grndvap(iland,                              &
                  sfcg%rhos                  (iwsfc), &
                  sfcg%canrrv                (iwsfc), &
                  land%nlev_sfcwater         (iland), &
                  land%surface_srrv          (iland), &
                  land%ground_rrv            (iland), &
                  land%sfcwater_energy (nlsw1,iland), &
                  land%soil_water        (nzg,iland), &
                  land%soil_energy       (nzg,iland), &
                  land%head              (nzg,iland), &
                  land%specifheat_drysoil(nzg,iland)  )

     !-----------------------------------------------------------------------------
     ! TEMPORARY UNTIL FULL LEAF-HYDRO MODEL IS COMPLETED WITH STREAM/RIVER RUNOFF:
     ! Simple representation of runoff

     if (land%nlev_sfcwater(iland) >= 1) then

        call remove_runoff(iland, iwsfc,                  &
                           sfcg%leaf_class       (iwsfc), &
                           land%sfcwater_mass  (1,iland), &
                           land%sfcwater_energy(1,iland), &
                           land%sfcwater_depth (1,iland), &
                           sfcg%runoff           (iwsfc)  )

     endif

     !-----------------------------------------------------------------------------

     ! Call land patch plot routine for selected iland values.

     if (iland == iland_print) call leaf_plot(                &
        iland,                                                &
        land%nlev_sfcwater(iland),                            &
        linit            = 1,                                 &
        lframe           = 1,                                 &
        ktrans           = ktrans,                            &
        sfcwater_tempk   = sfcwater_tempk      (1:nzs),       &
        sfcwater_fracliq = sfcwater_fracliq    (1:nzs),       &
        soil_tempk       = soil_tempk          (1:nzg),       &
        soil_fracliq     = soil_fracliq        (1:nzg),       &
        leaf_class       = sfcg%leaf_class           (iwsfc), &
        can_depth        = sfcg%can_depth            (iwsfc), &
        rhos             = sfcg%rhos                 (iwsfc), &
        vels             = sfcg%vels                 (iwsfc), &
        prss             = sfcg%prss                 (iwsfc), &
        pcpg             = sfcg%pcpg                 (iwsfc), &
        qpcpg            = sfcg%qpcpg                (iwsfc), &
        dpcpg            = sfcg%dpcpg                (iwsfc), &
        sxfer_t          = sfcg%sxfer_t              (iwsfc), &
        sxfer_r          = sfcg%sxfer_r              (iwsfc), &
        ustar            = sfcg%ustar                (iwsfc), &
        cantemp          = sfcg%cantemp              (iwsfc), &
        canrrv           = sfcg%canrrv               (iwsfc), &
        rshort           = sfcg%rshort               (iwsfc), &
        head1            = sfcg%head1                (iwsfc), &
        head0            = land%head0                (iland), &
        rshort_v         = land%rshort_v             (iland), &
        rshort_g         = land%rshort_g             (iland), &
        rlong_v          = land%rlong_v              (iland), &
        rlong_s          = land%rlong_s              (iland), &
        rlong_g          = land%rlong_g              (iland), &
        veg_height       = land%veg_height           (iland), &
        veg_rough        = land%veg_rough            (iland), &
        veg_tai          = land%veg_tai              (iland), &
        veg_lai          = land%veg_lai              (iland), &
        hcapveg          = land%hcapveg              (iland), &
        snowfac          = land%snowfac              (iland), &
        vf               = land%vf                   (iland), &
        surface_srrv     = land%surface_srrv         (iland), &
        ground_rrv       = land%ground_rrv           (iland), &
        veg_water        = land%veg_water            (iland), &
        veg_temp         = land%veg_temp             (iland), &
        transp           = transp                           , &
        stom_resist      = land%stom_resist          (iland), &
        veg_fracarea     = land%veg_fracarea         (iland), &
        sfcwater_mass    = land%sfcwater_mass  (1:nzs,iland), &
        sfcwater_energy  = land%sfcwater_energy(1:nzs,iland), &
        sfcwater_depth   = land%sfcwater_depth (1:nzs,iland), &
        soil_water       = land%soil_water     (1:nzg,iland), &
        soil_energy      = land%soil_energy    (1:nzg,iland), &
        head             = land%head           (1:nzg,iland)  )

  enddo  ! iland
  !$omp end do
  !$omp end parallel

  ! MPI SEND/RECV of quantities updated in seacells, lakecells, land cells

  if (iparallel == 1) then
     call mpi_send_wsfc(set='sfc_driv_end')
     call mpi_recv_wsfc(set='sfc_driv_end')
  endif

end subroutine surface_driver

!=========================================================================

subroutine sfcg_avgatm()

  use mem_basic,    only: rho, press, theta, tair, rr_v, vxe, vye, vze
  use mem_micro,    only: rr_c
  use mem_co2,      only: co2flag, rr_co2
  use misc_coms,    only: isubdomain, iparallel
  use consts_coms,  only: grav
  use mem_ijtabs,   only: itab_w
  use mem_grid,     only: gdz_abov8
  use mem_para,     only: myrank
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc
  use mem_sfcg,     only: mwsfc, sfcg, itab_wsfc
  use mem_sea,      only: sea, omsea

  implicit none

  integer :: iwsfc, j, iw, kw, isea
  real :: vels, psfc

  ! Average atmospheric quantities to each SFC grid cell location

  ! Loop over all SFC grid cells

  !$omp parallel do private(j,iw,kw,vels,psfc)
  do iwsfc = 2, mwsfc

     ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     sfcg%vels    (iwsfc) = 0.
     sfcg%prss    (iwsfc) = 0.
     sfcg%rhos    (iwsfc) = 0.
     sfcg%airtemp (iwsfc) = 0.
     sfcg%airtheta(iwsfc) = 0.
     sfcg%airrrv  (iwsfc) = 0.

     if (co2flag /= 0) then
        sfcg%airco2(iwsfc) = 0.
     endif

     ! Loop over all ATM grid cells that couple to this SFC grid cell

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
        kw = itab_wsfc(iwsfc)%kwatm(j)

        vels = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
        psfc = press(kw,iw) + gdz_abov8(kw-1) * rho(kw,iw)  ! hydrostatic eqn.

        sfcg%vels    (iwsfc) = sfcg%vels    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * vels
        sfcg%prss    (iwsfc) = sfcg%prss    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * psfc
        sfcg%rhos    (iwsfc) = sfcg%rhos    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rho  (kw,iw)
        sfcg%airtemp (iwsfc) = sfcg%airtemp (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * tair (kw,iw)
        sfcg%airtheta(iwsfc) = sfcg%airtheta(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * theta(kw,iw)

        if (allocated(rr_c)) then
           sfcg%airrrv(iwsfc) = sfcg%airrrv(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * (rr_v(kw,iw) + rr_c(kw,iw))
        else
           sfcg%airrrv(iwsfc) = sfcg%airrrv(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rr_v(kw,iw)
        endif

        if (co2flag /= 0) then
           sfcg%airco2(iwsfc) = sfcg%airco2(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rr_co2(kw,iw)
        endif

        ! Surface stress components for sea cells

        if (sfcg%leaf_class(iwsfc) == 0) then
           isea = iwsfc - omsea

           if (j == 1) then
              sea%windxe(isea) = itab_wsfc(iwsfc)%arcoarsfc(j) * vxe(kw,iw)
              sea%windye(isea) = itab_wsfc(iwsfc)%arcoarsfc(j) * vye(kw,iw)
              sea%windze(isea) = itab_wsfc(iwsfc)%arcoarsfc(j) * vze(kw,iw)
           else
              sea%windxe(isea) = sea%windxe(isea) + itab_wsfc(iwsfc)%arcoarsfc(j) * vxe(kw,iw)
              sea%windye(isea) = sea%windye(isea) + itab_wsfc(iwsfc)%arcoarsfc(j) * vye(kw,iw)
              sea%windze(isea) = sea%windze(isea) + itab_wsfc(iwsfc)%arcoarsfc(j) * vze(kw,iw)
           endif
        endif

     enddo

  enddo
  !$omp end parallel do

  ! MPI send/recv sfcg%vels, sfcg%prss, sfcg%rhos, sfcg%airtemp,
  ! sfcg%airtheta, sfcg%airrrv, sea%windxe, sea%windye, sea%windze 

  if (iparallel == 1) then
     call mpi_send_wsfc(set='avgatm')
     call mpi_recv_wsfc(set='avgatm')
  endif

end subroutine sfcg_avgatm

