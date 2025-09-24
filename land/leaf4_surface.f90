Module leaf4_surface

Contains

  subroutine skncomp_diagnose(iland, iwsfc, skncomp, sfcwater_mass, sfcwater_energy, &
                              sfcwater_depth, soil_water, soil_energy, specifheat_drysoil)

  use mem_land,   only: nzg, nzs
  use leaf_coms,  only: wcap_min
  use therm_lib,  only: qtk, qwtk
  use oname_coms, only: nl

  implicit none

  integer, intent(in   ) :: iland   ! land cell horizontal index []
  integer, intent(in   ) :: iwsfc   ! sfcgrid cell horizontal index []
  integer, intent(inout) :: skncomp ! skinlayer composition type (1, 2, or 3)

  real, intent(inout) :: sfcwater_mass   (nzs) ! sfcwater mass [kg/m^2]
  real, intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
  real, intent(inout) :: sfcwater_depth  (nzs) ! surface water energy [J/kg]

  real, intent(inout) :: soil_water        (nzg) ! soil water [water_vol/total_vol]
  real, intent(inout) :: soil_energy       (nzg) ! soil internal energy [J/m^3]
  real, intent(in   ) :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]

  ! Local variables

  real :: sfcwater_tempk     ! sfcwater(1) temperature [K]
  real :: sfcwater_fracliq   ! fraction of sfcwater(1) in liquid phase
  real :: sfcwater_epm2(nzs) ! sfcwater energy per m^2 [J/m^2]

  real :: wxfer  ! water mass transfer between sfcwater(1) and sfcwater(2) [kg/m^2]
  real :: qwxfer ! water energy transfer between sfcwater(1) and sfcwater(2) [J/m^2]
  real :: dwxfer ! water thickness transfer between sfcwater(1) and sfcwater(2) [m]

  ! This subroutine defines the composition of the land surface skin layer based
  ! on the amount of sfcwater present, if any, and the liquid vs ice content of
  ! sfcwater.  The skinlayer is the uppermost portion of sfcwater and/or soil,
  ! and is usually assigned a thickness small enough that its temperature relaxes
  ! in less than 10 minutes in response to changes in energy fluxes.

  ! The following skncomp values define the 4 options for skinlayer composition:

  ! skncomp = 0:  USED IF AND ONLY IF NL%IGW_SPINUP = 1, WHICH IS FOR LONG-TERM
  !               GROUNDWATER SPIN-UP SIMULATIONS.  SFCWATER(2) AND SOIL(NZG+1) ARE
  !               ALWAYS INACTIVE, AND SFCWATER(1) AND SOIL(NZG) ARE ALWAYS
  !               THERMALLY COMBINED.

  ! skncomp = 1:  Sfcwater is either absent or is present only in small quantity.
  !               The skinlayer is comprised of the sfcwater and the top centimeter
  !               (or so) of the soil.  Both are thermally combined so that they
  !               share the same temperature and fracliq values.  The sfcwater
  !               component has vertical index (1), and the soil component has
  !               vertical index (nzg+1).  The (nzg+1) soil layer is segmented
  !               from the (nzg) layer, and the effective thickness of the (nzg)
  !               soil layer is correspondingly reduced from its nominal value
  !               so that total mass and energy are conserved.

  ! skncomp = 2:  Enough sfcwater is present to fully comprise the skinlayer, but
  !               not enough for sfcwater to be divided into two layers.  The
  !               skinlayer is sfcwater(k=1).  Soil layer nzg+1 is not active as
  !               an independent layer, but instead is reincorporated within
  !               the nzg layer, which assumes its nominal thickness.

  ! skncomp = 3:  Enough sfcwater is present for division into two layers, and the
  !               skinlayer is the upper layer, sfcwater(k=2).  Soil layer nzg+1 is
  !               not active as an independent layer, but instead is reincorporated
  !               within the nzg layer, which assumes its nominal thickness.

  ! If current sfcwater conditions warrant a change from the previous skncomp value,
  ! reapportion physical components of the skinlayer as appropriate and update skncomp.

  !-----------------------------------------------------------------------------------

  ! If this run is performing a long-term groundwater spin-up (igw_spinup = 1),
  ! zero out values for inactive sfcwater(2) layer, set inactive soil(nzg+1)
  ! layer water and energy to (nzg) values, set skncomp = 0, and return.

  if (nl%igw_spinup == 1) then
     sfcwater_mass  (2) = 0.
     sfcwater_energy(2) = 0.
     sfcwater_depth (2) = 0.

     skncomp = 0

     return
  endif

! If sfcwater_mass is less than threshold, set it to zero and return.

  if (sum(sfcwater_mass(:)) < wcap_min) then

     sfcwater_mass  (:) = 0.
     sfcwater_energy(:) = 0.
     sfcwater_depth (:) = 0.

     skncomp = 1

     return
  endif

  ! Diagnose sfcwater_tempk and sfcwater_fracliq for sfcwater layer (1)

  call qtk(sfcwater_energy(1),sfcwater_tempk,sfcwater_fracliq)

!!  if (sum(sfcwater_mass(:)) < 5.0) then
!!
!!     ! Under this low sfcwater_mass condition, the skinlayer must include the
!!     ! nzg+1 soil layer, so the nzg+1 layer must be segmented from the nzg
!!     ! layer.  If skncomp was /= 1 upon entering this subroutine, perform the
!!     ! separation here.  Skncomp will be set to 1 below, which signals that
!!     ! the effectve thickness of the nzg layer is reduced from its nominal
!!     ! value, which conserves total mass and energy.
!!
!!     ! If skncomp was previously = 3, combine both layers of sfcwater into
!!     ! a single layer.  Otherwise, just diagnose sfcwater_epm2 for layer (1).
!!
!!     if (skncomp == 3) then
!!        sfcwater_epm2  (1) = sum(sfcwater_mass (:) * sfcwater_energy(:))
!!        sfcwater_mass  (1) = sum(sfcwater_mass (:))
!!        sfcwater_depth (1) = sum(sfcwater_depth(:))
!!        sfcwater_energy(1) = sfcwater_epm2(1) / max(wcap_min,sfcwater_mass(1))
!!     else
!!        sfcwater_epm2(1) = sfcwater_mass(1) * sfcwater_energy(1)
!!     endif
!!
!!     ! Set sfcwater mass, energy, and depth in layer (2) to zero.
!!
!!     sfcwater_mass  (2) = 0.
!!     sfcwater_energy(2) = 0.
!!     sfcwater_depth (2) = 0.
!!
!!     ! Thermally combine sfcwater(1) and soil(nzg)
!!
!!     call sfcwater_soil_comb(iland, iwsfc,                         &
!!                             soil_water(nzg),soil_energy(nzg), &
!!                             specifheat_drysoil(nzg),              &
!!                             sfcwater_mass(1),sfcwater_epm2(1),    &
!!                             sfcwater_tempk,sfcwater_fracliq       )
!!
!!     skncomp = 1

  if (sum(sfcwater_mass(:)) < 12.0 .or. sfcwater_fracliq > 0.3 .or. &
       (sfcwater_fracliq  > 0.2 .and. skncomp /= 2))         then
!            (sfcwater_fracliq  > 0.2 .and. skncomp /= 3))         then

     ! (Use of skncomp in the above elseif statement provides hysteresis that
     ! prevents oscillating transition between skncomp values 1 and 2.)

     ! If skncomp was previously = 2, combine both layers of sfcwater into
     ! a single layer.

!!     if (iwsfc == 6238) then
!!        write(*,*) sfcwater_mass(:)
!!        write(*,*) sfcwater_energy(:)
!!     endif

     if (skncomp == 2) then
        sfcwater_epm2  (1) = sum(sfcwater_mass (:) * sfcwater_energy(:))
        sfcwater_mass  (1) = sum(sfcwater_mass (:))
        sfcwater_depth (1) = sum(sfcwater_depth(:))
        sfcwater_energy(1) = sfcwater_epm2(1) / max(wcap_min,sfcwater_mass(1))
     else
        sfcwater_epm2  (1) = sfcwater_mass(1) * sfcwater_energy(1)
     endif

     ! Set sfcwater mass, energy, and depth in layer (2) to zero.

     sfcwater_mass  (2) = 0.
     sfcwater_energy(2) = 0.
     sfcwater_depth (2) = 0.

     skncomp = 1

  else

     ! Transfer mass and energy between sfcwater layers (1) and (2) when necessary in
     ! order to maintain mass in layer (2) within a prescribed range.  If skncomp was
     ! previously < 2, then this operation newly separates layer (2) from layer (1).

     if (sfcwater_mass(2) < 4.5) then

        wxfer  = 5. - sfcwater_mass(2)
        qwxfer = wxfer * sfcwater_energy(1)
        dwxfer = wxfer * sfcwater_depth (1) / sfcwater_mass(1)

        sfcwater_epm2  (2) = sfcwater_mass  (2) * sfcwater_energy(2) + qwxfer
        sfcwater_mass  (2) = sfcwater_mass  (2) + wxfer
        sfcwater_depth (2) = sfcwater_depth (2) + dwxfer
        sfcwater_energy(2) = sfcwater_epm2  (2) / sfcwater_mass(2)

        sfcwater_epm2  (1) = sfcwater_mass  (1) * sfcwater_energy(1) - qwxfer
        sfcwater_mass  (1) = sfcwater_mass  (1) - wxfer
        sfcwater_depth (1) = sfcwater_depth (1) - dwxfer
        sfcwater_energy(1) = sfcwater_epm2  (1) / sfcwater_mass(1)

     elseif (sfcwater_mass(2) > 5.5) then

        wxfer  = 5. - sfcwater_mass(2)
        qwxfer = wxfer * sfcwater_energy(2)
        dwxfer = wxfer * sfcwater_depth (2) / sfcwater_mass(2)

        sfcwater_epm2  (2) = sfcwater_mass  (2) * sfcwater_energy(2) + qwxfer
        sfcwater_mass  (2) = sfcwater_mass  (2) + wxfer
        sfcwater_depth (2) = sfcwater_depth (2) + dwxfer
        sfcwater_energy(2) = sfcwater_epm2  (2) / sfcwater_mass(2)

        sfcwater_epm2  (1) = sfcwater_mass  (1) * sfcwater_energy(1) - qwxfer
        sfcwater_mass  (1) = sfcwater_mass  (1) - wxfer
        sfcwater_depth (1) = sfcwater_depth (1) - dwxfer
        sfcwater_energy(1) = sfcwater_epm2  (1) / sfcwater_mass(1)

     else

        sfcwater_epm2  (1) = sfcwater_mass  (1) * sfcwater_energy(1)

     endif

     skncomp = 2

  endif

  ! Thermally combine sfcwater(1) and soil(nzg)

  call sfcwater_soil_comb(iland, iwsfc,                       &
                          soil_water(nzg), soil_energy(nzg),  &
                          specifheat_drysoil(nzg),            &
                          sfcwater_mass(1), sfcwater_epm2(1), &
                          sfcwater_tempk, sfcwater_fracliq    )

  sfcwater_energy(1) = sfcwater_epm2(1) / max(wcap_min,sfcwater_mass(1))

  end subroutine skncomp_diagnose

  !=============================================================================

  subroutine sfcwater(iland, iwsfc, wfree1, qwfree1, dwfree1, head1,           &
                      skncomp, sfcwater_mass, sfcwater_energy, sfcwater_depth, &
                      sfcwater_tempk, sfcwater_fracliq, sfcwater_epm2,         &
                      soil_water, soil_energy, specifheat_drysoil,             &
                      soil_tempk, thermcond_soil )

  use leaf_coms,   only: dt_leaf, wcap_min
  use mem_land,    only: nzg, nzs, dslzi, dslzo2
  use consts_coms, only: alvi, cice, cliq, alli
  use leaf4_plot,  only: leaf_plot
  use therm_lib,   only: qwtk

  implicit none

  integer, intent(in)    :: iland         ! current land cell index
  integer, intent(in)    :: iwsfc         ! current sfcg cell index
  real,    intent(out)   :: wfree1        ! free liquid in lowest sfcwater layer [kg/m^2]
  real,    intent(out)   :: qwfree1       ! energy carried by wfree1 [J/m^2]
  real,    intent(out)   :: dwfree1       ! depth carried by wfree1 [m]
  real,    intent(inout) :: head1         ! Top boundary hydraulic head for soil model [m]
  integer, intent(in)    :: skncomp       ! surface skinlayer composition type []

  real,    intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
  real,    intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
  real,    intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
  real,    intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
  real,    intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
  real,    intent(inout) :: sfcwater_epm2   (nzs) ! sfcwater energy [J/m^2]

  real,    intent(inout) :: soil_water        (nzg) ! soil water [water_vol/total_vol]
  real,    intent(inout) :: soil_energy       (nzg) ! soil internal energy [J/m^3]
  real,    intent(in)    :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]
  real,    intent(in)    :: soil_tempk        (nzg) ! soil temperature [K]
  real,    intent(in)    :: thermcond_soil    (nzg) ! soil thermal conductivity [W/(K m)]

  ! Local variables

  real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]

  real :: soil_rfactor ! soil thermal resistance [K m^2/W]
  real :: hxferss   ! sensible heat transfer from sfcwater(1) to sfcwater(2) [J/m^2]
  real :: hxfergs   ! sensible heat from soil(nzg or nzg+1) to sfcwater(1) [J/m^2]
  real :: snden     ! sfcwater density [kg/m^3]
  real :: wfree     ! free liquid in sfcwater layer that can percolate out [kg/m^2]
  real :: qwfree    ! energy carried by wfree [J/m^2]
  real :: dwfree    ! depth carried by wfree [m]
  real :: fracstep  ! ratio of leaf timestep to snow density exponential decay time
  real :: tempk     ! Kelvin temperature of snowcover, limited to max of 273.15 [K]

  ! Local parameters

  integer, parameter :: iland_print = 0

  ! Initialize free water values to zero

  wfree1   = 0.
  qwfree1  = 0.
  dwfree1  = 0.
  head1    = 0.

  ! Take inventory of sfcwater_mass now that exchange with canopy is
  ! complete.  If sfcwater_mass is less than threshold, set it to zero and return.

  if (sum(sfcwater_mass(:)) < wcap_min) then

     sfcwater_mass  (:) = 0.
     sfcwater_energy(:) = 0.
     sfcwater_depth (:) = 0.
     sfcwater_epm2  (:) = 0.

     sfcwater_tempk  (:) = soil_tempk(nzg)
     sfcwater_fracliq(:) = 0.

     return
  endif

  ! If we got here, some sfcwater is present.  If implicit heat balance between
  ! sfcwater and soil(nzg) was NOT done in subroutine canopy, compute heat
  ! transfer between sfcwater(1) and soil(nzg).

  if (skncomp > 1) then

     ! Compute snow heat resistance times HALF layer depth (sfcwater_rfactor).
     ! Sfcwater_tempk should be correctly balanced value at this point, so
     ! sfcwater_rfactor should have correct value.  Formula applies to snow,
     ! so limit temperature to no greater than 273.15.

     call qwtk(sfcwater_epm2(1),sfcwater_mass(1),10.,sfcwater_tempk(1),sfcwater_fracliq(1))

     snden = sfcwater_mass(1) / sfcwater_depth(1)
     tempk = min(273.15,sfcwater_tempk(1))

     sfcwater_rfactor(1) = .5 * sfcwater_depth(1) &
        / (1.093e-3 * exp(.028 * tempk) *         &
        (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

     ! If second sfcwater layer is active, compute rfactor for it also, evaluate
     ! sensible heat transfer between 2 sfcwater layers, and apply the heat
     ! transfer to both layers

!    if (skncomp == 3) then

        call qwtk(sfcwater_epm2(2),sfcwater_mass(2),10.,sfcwater_tempk(2),sfcwater_fracliq(2))

        snden = sfcwater_mass(2) / sfcwater_depth(2)
        tempk = min(273.15,sfcwater_tempk(2))

        sfcwater_rfactor(2) = .5 * sfcwater_depth(2) &
           / (1.093e-3 * exp(.028 * tempk) *         &
             (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

        hxferss = dt_leaf * (sfcwater_tempk(1) - sfcwater_tempk(2)) &
                / (sfcwater_rfactor(1) + sfcwater_rfactor(2))

        sfcwater_epm2(1) = sfcwater_epm2(1) - hxferss
        sfcwater_epm2(2) = sfcwater_epm2(2) + hxferss

!!        if (iwsfc == 2024) then
!!           write(*,*) "sfcwater3: ", sfcwater_tempk(2), sfcwater_tempk(1), hxferss
!!        endif

!    endif

!!     ! Compute heat resistance of top HALF of top soil layer (soil_rfactor).
!!
!!     soil_rfactor = dslzo2(nzg) / thermcond_soil(nzg)
!!
!!     ! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
!!     ! layer, and apply to both.
!!
!!     hxfergs = dt_leaf * (soil_tempk(nzg) - sfcwater_tempk(1))   &
!!             / (soil_rfactor + sfcwater_rfactor(1))
!!
!!     if (iwsfc == 2024) then
!!        write(*,*) "sfcwater: ", sfcwater_tempk(1), soil_tempk(nzg), hxfergs
!!     endif
!!
!!     sfcwater_epm2(1) = sfcwater_epm2(1) + hxfergs
!!
!!     soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)
!!
!!  endif

  ! Now that new mass and energy contributions have been applied, process
  ! existing sfcwater layers beginning with top layer

!!  if (skncomp == 2) then

     ! Rediagnose sfcwater temperature and liquid fraction.  Use qwtk instead of qtk because
     ! sfcwater_epm2 is used instead of sfcwater_energy.  ("dryhcap" = 10 is very small value)

     call qwtk(sfcwater_epm2(2), sfcwater_mass(2), 10., sfcwater_tempk(2), sfcwater_fracliq(2))

     ! Diagnose sfcwater density.  Make sure sfcwater_depth is not too near zero.

     sfcwater_depth(2) = max(sfcwater_depth(2), .001 * sfcwater_mass(2))
     snden = sfcwater_mass(2) / sfcwater_depth(2)

     ! Assume that as snow ages on ground, its density difference with a limiting
     ! maximum value (currently set to 400 kg/m^3) decays exponentially (with a
     ! decay time currently set to about 3 weeks).  If sfcwater density is less
     ! than this limiting value, apply the density increase for this timestep.

     ! This formulation and decay constants are very crude approximations to a few
     ! widely variable observations of snowpack density and are intended only as a
     ! rough representation of the tendency for snowcover to compress with time.
     ! A better formulation that accounts for several environmental factors would
     ! be desirable here.

     if (snden < 400.) then
        fracstep = .5e-6 * dt_leaf  ! .5e-6 is inverse decay time scale
        snden = snden * (1. - fracstep) + 400. * fracstep
        sfcwater_depth(2) = sfcwater_mass(2) / snden
     endif

     ! If liquid exists in current sfcwater layer, any low-density ice structure
     ! tends to collapse.  Increase density accordingly using simple linear relation.

     if (snden < 1.e3 * sfcwater_fracliq(2)) then
        snden = 1.e3 * sfcwater_fracliq(2)
        sfcwater_depth(2) = sfcwater_mass(2) / snden
     endif

     ! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

     if (sfcwater_fracliq(2) > .10) then

        wfree = sfcwater_mass(2) * (sfcwater_fracliq(2) - .10) / .90

        ! Evaluate energy and depth contained in wfree (which is in liquid phase)

        qwfree = wfree * (cliq * (sfcwater_tempk(2) - 273.15) + alli)
        dwfree = wfree * .001

        ! Drain mass, energy, and depth from this layer and add to layer below.
        ! Check if essentially all of sfcwater_mass(k) will drain from layer

        if (wfree > .999 * sfcwater_mass(2)) then

           ! All sfcwater_mass(2) drains from layer.  Set layer quantities to zero to
           ! avoid truncation error.

           sfcwater_mass  (2) = 0.
           sfcwater_epm2  (2) = 0.
           sfcwater_depth (2) = 0.

        else

           ! Not all sfcwater_mass(2) drains from layer.  Drain mass, energy, and depth
           ! of free water out of current layer

           sfcwater_mass  (2) = sfcwater_mass (2) - wfree
           sfcwater_epm2  (2) = sfcwater_epm2 (2) - qwfree
           sfcwater_depth (2) = sfcwater_depth(2) - dwfree

        endif

        sfcwater_mass (1) = sfcwater_mass (1) + wfree
        sfcwater_epm2 (1) = sfcwater_epm2 (1) + qwfree
        sfcwater_depth(1) = sfcwater_depth(1) + dwfree

     endif

     ! Now that all contributions to sfcwater_mass(2) and sfcwater_epm2(2) are
     ! complete, rediagnose sfcwater_energy(2)

     sfcwater_energy(2) = sfcwater_epm2(2) / max(wcap_min,sfcwater_mass(2))

  endif

  ! Now, do sfcwater layer (1) as done for (2) above; comments omitted this time

!  if (skncomp > 1) then
!     call qwtk(sfcwater_epm2(1), sfcwater_mass(1), 10., sfcwater_tempk(1), sfcwater_fracliq(1))
!  else
     call sfcwater_soil_comb(iland, iwsfc, soil_water(nzg), soil_energy(nzg), &
                             specifheat_drysoil(nzg), sfcwater_mass(1), sfcwater_epm2(1), &
                             sfcwater_tempk(1), sfcwater_fracliq(1))
!  endif

  sfcwater_energy(1) = sfcwater_epm2(1) / max(wcap_min,sfcwater_mass(1))

  sfcwater_depth(1) = max(sfcwater_depth(1), .001 * sfcwater_mass(1))
  snden = sfcwater_mass(1) / sfcwater_depth(1)

  if (snden < 400.) then
     fracstep = .5e-6 * dt_leaf  ! .5e-6 is inverse decay time scale
     snden = snden * (1. - fracstep) + 400. * fracstep
     sfcwater_depth(1) = sfcwater_mass(1) / snden
  endif

  if (snden < 1.e3 * sfcwater_fracliq(1)) then
     snden = 1.e3 * sfcwater_fracliq(1)
     sfcwater_depth(1) = sfcwater_mass(1) / snden
  endif

  if (sfcwater_fracliq(1) > .10) then
     wfree = sfcwater_mass(1) * (sfcwater_fracliq(1) - .10) / .90

     qwfree = wfree * (cliq * (sfcwater_tempk(1) - 273.15) + alli)
     dwfree = wfree * .001

     ! Copy any free water in sfcwater layer (1) to output arrays

     wfree1  = wfree
     qwfree1 = qwfree
     dwfree1 = dwfree
  endif

  ! Sfcwater mass, depth, epm2, and energy in layer (1) are not updated here, but
  ! instead in subroutine soil after it is determined how much free water (wfree1)
  ! can actually drain from sfcwater(1) and enter the top soil layer.

  ! 11 Feb 2015: Define head1 from total sfcwater_mass (but leaf4_soil will
  ! still allow no more than wfree1 to actually enter soil on current timestep

  ! July 2017: With LEAF as a 3D groundwater model, limit on head1 has been removed.

  head1 = .001 * sum(sfcwater_mass(1:nzs))

  if (iland == iland_print) then

     call leaf_plot(iland,                             &
                    linit           = 1,               &
                    lframe          = 1,               &
                    sfcwater_mass   = sfcwater_mass,   &
                    sfcwater_energy = sfcwater_energy, &
                    sfcwater_epm2   = sfcwater_epm2,   &
                    sfcwater_depth  = sfcwater_depth   )

  endif

  end subroutine sfcwater

!===============================================================================

  subroutine sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy,       &
                                specifheat_drysoil, sfcwater_mass, sfcwater_epm2, &
                                sfcwater_tempk, sfcwater_fracliq                  )

  use mem_land,    only: nzg, dslz, dslzi
  use consts_coms, only: cice, cliq, alli
  use therm_lib,   only: qwtk

  implicit none

  integer, intent(in) :: iland              ! current land cell index
  integer, intent(in) :: iwsfc              ! current sfcg cell index

  real, intent(in)    :: soil_water         ! soil water [water_vol/total_vol]
  real, intent(inout) :: soil_energy        ! soil internal energy [J/m^3]
  real, intent(in)    :: specifheat_drysoil ! specific heat of dry soil [J/(m^3 K)]
  real, intent(in)    :: sfcwater_mass      ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_epm2      ! sfcwater energy [J/m^2]
  real, intent(out)   :: sfcwater_tempk     ! surface water temp [K]
  real, intent(out)   :: sfcwater_fracliq   ! fraction of sfc water in liq phase

  ! Local variables

  real :: w_comb       ! (sfcwater + soil) water mass [kg/m^2]
  real :: qw_comb      ! (sfcwater + soil) energy [J/m^2]
  real :: tempk_comb   ! (sfcwater + soil) Kelvin temp [K]
  real :: fracliq_comb ! (sfcwater + soil) frac of water in liq phase
  real :: hcapsoil     ! soil heat capacity [J/(m^2 K)]
  real :: flmin        ! lower bound on sfcwater_fracliq in balance with soil
  real :: flmax        ! upper bound on sfcwater_fracliq in balance with soil

  ! Combined sfcwater and soil water mass per square meter

  w_comb = sfcwater_mass + soil_water * 1.e3 * dslz(nzg)

  ! Combined sfcwater and soil energy per square meter

  qw_comb = sfcwater_epm2 + soil_energy * dslz(nzg)

  ! Soil heat capacity per square meter

  hcapsoil = specifheat_drysoil * dslz(nzg)

  ! Diagnose equilibrium temperature and fractional liquid/ice water phases

  call qwtk(qw_comb,w_comb,hcapsoil,tempk_comb,fracliq_comb)

!!  if (iwsfc == 6238) then
!!     write(*,*) tempk_comb, fracliq_comb, sfcwater_epm2, soil_energy * dslz(nzg)
!!  endif

  ! Diagnose new energy value for sfcwater based on qw_comb value.

  if (qw_comb < 0.) then

     ! Case of equilibrium temperature below 0 deg C.  Sfcwater_fracliq = 0.

     sfcwater_fracliq = 0.
     sfcwater_tempk   = tempk_comb
     sfcwater_epm2    = sfcwater_mass * cice * (tempk_comb - 273.15)

  elseif (qw_comb > w_comb * alli) then

     ! Case of equilibrium temperature above 0 deg C.  Sfcwater_fracliq = 1.

     sfcwater_fracliq = 1.
     sfcwater_tempk   = tempk_comb
     sfcwater_epm2    = sfcwater_mass * (cliq * (tempk_comb - 273.15) + alli)

  else

     ! Equilibrium temperature is 0 deg C.  If sfcwater_mass is near zero,
     ! assume identical values for sfcwater_fracliq and soil_fracliq.

     if (sfcwater_mass < 1.e-6) then

        sfcwater_fracliq = fracliq_comb
        sfcwater_tempk   = 273.15
        sfcwater_epm2    = sfcwater_mass * sfcwater_fracliq * alli

     else

        ! If sfcwater_mass is larger, determine separate values for
        ! sfcwater_fracliq and soil_fracliq using constraint that the sum
        ! of (mass * fracliq) over both components is (w_comb * fracliq_comb).

        ! Lower bound on sfcwater_fracliq: case with soil_water all liquid

        flmin = (fracliq_comb * w_comb - soil_water * 1.e3 * dslz(nzg)) &
              / sfcwater_mass

        ! Upper bound on sfcwater_fracliq: case with soil_water all ice

        flmax = fracliq_comb * w_comb / sfcwater_mass

        ! New sfcwater_fracliq value becomes closest value within bounds to old value

        sfcwater_fracliq = max(0.,flmin,min(1.0,flmax,sfcwater_fracliq))
        sfcwater_tempk   = 273.15
        sfcwater_epm2    = sfcwater_mass * sfcwater_fracliq * alli

     endif

  endif

  ! New energy value for soil is combined energy minus new sfcwater energy

  soil_energy = (qw_comb - sfcwater_epm2) * dslzi(nzg)

  end subroutine sfcwater_soil_comb

!===============================================================================

  subroutine remove_runoff(iland, iwsfc, leaf_class, sfcwater_mass, sfcwater_energy, &
                           sfcwater_depth, runoff)

  use leaf_coms,   only: dt_leaf
  use consts_coms, only: alli, cliq
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in) :: iland
  integer, intent(in) :: iwsfc
  integer, intent(in) :: leaf_class

  real, intent(inout) :: sfcwater_mass
  real, intent(inout) :: sfcwater_energy
  real, intent(inout) :: sfcwater_depth
  real, intent(inout) :: runoff

  ! Local variables

  real :: sfcwater_tempk
  real :: sfcwater_fracliq
  real :: sfcwater_epm2
  real :: sfcwater_mass_thresh
  real :: wfree
  real :: qrunoff, drunoff

  real, parameter :: runoff_time = 86400.  ! time scale for surface-water runoff [s]

  ! Set a threshold value of sfcwater_mass for runoff to occur that is based on
  ! leaf_class

  if (leaf_class == 17 .or. leaf_class == 20) then
     sfcwater_mass_thresh = 1000.  ! 1000 kg/m^2 equivalent to 1.0 m depth
  else
     sfcwater_mass_thresh = 1.     ! 1 kg/m^2 equivalent to 1.0 mm depth
  endif

  if (sfcwater_mass <= sfcwater_mass_thresh) then
     runoff = 0.0
     return
  endif

  call qtk(sfcwater_energy, sfcwater_tempk, sfcwater_fracliq)

  ! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

  if (sfcwater_fracliq > .10) then
     wfree = sfcwater_mass * (sfcwater_fracliq - .10) / .90

     ! Fraction of wfree removed as runoff this time step
     ! (In future development, runoff_time should be function of topography)

     runoff = wfree * dt_leaf / runoff_time

     ! Evaluate energy and depth contained in runoff (which is in liquid phase)

     qrunoff = runoff * (cliq * (sfcwater_tempk - 273.15) +  alli)
     drunoff = runoff * 0.001

     ! Subtract runoff from sfcwater mass, energy, and depth

     sfcwater_epm2 = sfcwater_mass * sfcwater_energy
     sfcwater_mass = sfcwater_mass - runoff
     sfcwater_energy = (sfcwater_epm2 - qrunoff) / sfcwater_mass
     sfcwater_depth = sfcwater_depth - drunoff
  endif

  end subroutine remove_runoff

!===============================================================================

  subroutine grndvap(iland, rhos, canshv, surface_ssh, ground_shv, skncomp, &
                     sfcwater_mass, sfcwater_energy, soil_water, soil_energy, &
                     head, specifheat_drysoil, wresid, sfldcap)

    use consts_coms, only: grav, rvap
    use therm_lib,   only: qtk, rhovsil, qwtk
    use leaf_coms,   only: wcap_min

    implicit none

    integer, intent(in)  :: iland              ! current land cell number
    real,    intent(in)  :: rhos               ! air density [kg/m^3]
    real,    intent(in)  :: canshv             ! canopy vapor spec hum [kg_vap/kg_air]
    real,    intent(out) :: surface_ssh        ! surface (saturation) spec hum [kg_vap/kg_air]
    real,    intent(out) :: ground_shv         ! ground equilibrium spec hum [kg_vap/kg_air]
    integer, intent(in)  :: skncomp            ! surface skinlayer composition type (1, 2, or 3)
    real,    intent(in)  :: sfcwater_mass(:)   ! [kg/m^2]
    real,    intent(in)  :: sfcwater_energy(:) ! [J/kg]
    real,    intent(in)  :: soil_water         ! soil water content [vol_water/vol_tot]
    real,    intent(in)  :: soil_energy        ! [J/m^3]
    real,    intent(in)  :: head               ! hydraulic head of top soil layer [m]
    real,    intent(in)  :: specifheat_drysoil
    real,    intent(in)  :: wresid             ! residual soil water content [vol_water/vol_tot]
    real,    intent(in)  :: sfldcap            ! soil water content at field capacity [vol_water/vol_tot]

    ! Local variables

    real :: tempk     ! surface water temp [K]
    real :: fracliq   ! fraction of surface water in liquid phase
    real :: sfc_rhovs ! ground sfc saturation vapor density [kg_vap/m^3]
    real :: alpha     ! "alpha" term in Lee and Pielke (1992)
    real :: beta      ! "beta" term in Lee and Pielke (1992)

    if (sfcwater_mass(1) < wcap_min) then

       ! Without sfcwater, sfc_rhovs is the saturation vapor density at the
       ! temperature of the soil surface, and is used for computing dew/frost
       ! formation only.  skncomp = 1 in this case.

       call qwtk(soil_energy,soil_water*1.e3,specifheat_drysoil,tempk,fracliq)
       call grndvap_ab(iland,tempk,soil_water,wresid,sfldcap,head,alpha,beta)

       sfc_rhovs = rhovsil(tempk-273.15)
       surface_ssh = sfc_rhovs / rhos
       ground_shv  = surface_ssh * alpha * beta + (1. - beta) * canshv

    else

       ! With surface water (or snowcover) present, sfc_rhovs is the saturation
       ! vapor density of the water/snow surface and is used for computing
       ! dew/frost formation and sfcwater evaporation.

       if (skncomp < 3) then
          call qtk(sfcwater_energy(1),tempk,fracliq)
       else
          call qtk(sfcwater_energy(2),tempk,fracliq)
       endif

       sfc_rhovs   = rhovsil(tempk-273.15)
       surface_ssh = sfc_rhovs / rhos
       ground_shv  = surface_ssh

    endif

  end subroutine grndvap

!===============================================================================

  subroutine grndvap_ab(iland,tempk,soil_water,wresid,sfldcap,head,alpha,beta)

    use consts_coms, only: grav, rvap, pi1

    implicit none

    integer, intent(in) :: iland       ! current land cell number

    real, intent(in)  :: tempk       ! soil temperature [K]
    real, intent(in)  :: soil_water  ! soil water content [vol_water/vol_tot]
    real, intent(in)  :: wresid      ! residual soil water content [vol_water/vol_tot]
    real, intent(in)  :: sfldcap     ! soil water content at field capacity [vol_water/vol_tot]
    real, intent(in)  :: head        ! hydraulic head [m]
    real, intent(out) :: alpha       ! "alpha" term in Lee and Pielke (1992)
    real, intent(out) :: beta        ! "beta" term in Lee and Pielke (1992)

    ! Local parameter

    real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

    if (head >= 0.) then
       alpha = 0.999
    else
       alpha = min( 0.999, exp( max(-80., gorvap * head / tempk) ) )
    endif

    if (soil_water >= sfldcap) then
       beta = 1.0
    elseif (soil_water <= wresid) then
       beta = 0.0
    else
       ! Original Lee & Pielke:
       ! beta = .25 * (1. - cos(soil_water / sfldcap(nts) * pi1) )**2
       !
       ! Noilhan & Planton:
       ! beta = .50 * (1. - cos(soil_water / sfldcap(nts) * pi1) )
       !
       ! My modification:
       beta = .50 * (1. - cos( (soil_water-wresid) / (sfldcap-wresid) * pi1 ) )
    endif

  end subroutine grndvap_ab

!===============================================================================

  subroutine sfcrad_land(iland, iwsfc, leaf_class,                                &
                         rshort, rlong, slong, vlong, gnd_albedo, gnd_emiss,      &
                         vf, veg_albedo, rshort_s, rlong_s, rshort_v, rlong_v     )

  use leaf_coms, only: emisv

  implicit none

  integer, intent(in)  :: iland         ! horizontal index of current land cell
  integer, intent(in)  :: iwsfc         ! horizontal index of current sfcgrid cell
  integer, intent(in)  :: leaf_class    ! leaf class
  real,    intent(in)  :: rshort        ! downward surface incident s/w rad flux [W/m^2]
  real,    intent(in)  :: rlong         ! downward surface incident l/w rad flux [W/m^2]
  real,    intent(in)  :: slong         ! gnd outgoing l/w rad flux [W/m^2]
  real,    intent(in)  :: vlong         ! veg outgoing l/w rad flux [W/m^2]
  real,    intent(in)  :: gnd_albedo    ! ground albedo
  real,    intent(in)  :: gnd_emiss     ! ground emiss
  real,    intent(in)  :: vf            ! fractional coverage of non-buried part of veg
  real,    intent(in)  :: veg_albedo    ! veg albedo
  real,    intent(out) :: rshort_s      ! s/w net rad flux to surface [W/m^2]
  real,    intent(out) :: rshort_v      ! s/w net rad flux to veg [W/m^2]
  real,    intent(out) :: rlong_s       ! l/w net rad flux to surface [W/m^2]
  real,    intent(out) :: rlong_v       ! l/w net rad flux to veg [W/m^2]

  ! Local variables

  real :: rlonga_v           ! longwave radiative flux from atm  to veg [W/m^2]
  real :: rlonga_s           ! longwave radiative flux from atm  to sfc [W/m^2]
  real :: rlongv_a           ! longwave radiative flux from veg  to atm [W/m^2]
  real :: rlongv_s           ! longwave radiative flux from veg  to sfc [W/m^2]
  real :: rlongs_a           ! longwave radiative flux from snow to atm [W/m^2]
  real :: rlongs_v           ! longwave radiative flux from snow to veg [W/m^2]
  real :: vfc                ! 1 - vf []
  real :: veg_emiss          ! vegetation emissivity

  ! This routine is called twice (for each land cell) by the radiation
  ! parameterization driver, performing exactly the same operations on both calls.
  ! All that is used from the first call are net surface albedo and upward
  ! longwave radiative flux for each land cell.  The radiation parameterization
  ! carries out atmospheric radiative transfer computations following this first
  ! call using the albedo and upward longwave flux.  The second call to this
  ! subroutine is made after the atmospheric radiative fluxes are computed.
  ! This call provides net radiative fluxes to vegetation, snowcover, and soil in
  ! each land cell, plus functions of snowcover, all of which are used in leaf4.

  ! Change for LEAF4: rshort_s always gets all shortwave that is absorbed by the
  ! surface, whether it is sfcwater, soil, or a combination. Subroutine leaf4_canopy
  ! further sorts out where rshort_s should go.

  vfc      = 1. - vf
  rshort_s = rshort * vfc * (1. - gnd_albedo)
  rshort_v = rshort * vf  * (1. - veg_albedo + vfc * gnd_albedo)
 !rshort_a = rshort * albedo_beam

  ! Longwave radiation calculations

  veg_emiss = emisv(leaf_class)

  rlonga_v = rlong * vf  * (veg_emiss + vfc * (1. - gnd_emiss))
  rlonga_s = rlong * vfc * gnd_emiss
  rlongv_s = vlong * vf  * gnd_emiss
  rlongv_a = vlong * vf  * (2. - gnd_emiss - vf + gnd_emiss * vf)
  rlongs_v = slong * vf  * veg_emiss
  rlongs_a = slong * vfc

  rlong_s = rlonga_s - rlongs_a + rlongv_s - rlongs_v
  rlong_v = rlonga_v - rlongv_a + rlongs_v - rlongv_s

  end subroutine sfcrad_land

!===============================================================================

  subroutine sfcrad_prep(iland, iwsfc, leaf_class, wnxl, wnyl, wnzl,              &
                         skncomp, sfcwater_mass, sfcwater_energy, sfcwater_depth, &
                         veg_tempk, veg_fracarea, veg_height, veg_albedo,         &
                         soil_energy, soil_water, wresid_vg, wsat_vg,             &
                         specifheat_drysoil, sand, snowfac, vf, cosz, rlongup,    &
                         rlong_albedo, albedo_beam, gnd_albedo, gnd_emiss,        &
                         slong, vlong                                             )

  use leaf_coms,   only: emisv, emisg, emisw, wcap_min
  use mem_land,    only: nzs
  use consts_coms, only: stefan
  use mem_radiate, only: sunx, suny, sunz
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in)    :: iland         ! horizontal index of current land cell
  integer, intent(in)    :: iwsfc         ! horizontal index of current sfcgrid cell
  integer, intent(in)    :: leaf_class    ! leaf class
  real,    intent(in)    :: wnxl          ! land cell norm unit vec x comp [m]
  real,    intent(in)    :: wnyl          ! land cell norm unit vec y comp [m]
  real,    intent(in)    :: wnzl          ! land cell norm unit vec z comp [m]
  integer, intent(in)    :: skncomp       ! skinlayer composition type (1, 2, or 3)
  real,    intent(in)    :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
  real,    intent(in)    :: sfcwater_energy(nzs) ! surface water internal energy [J/kg]
  real,    intent(in)    :: sfcwater_depth (nzs) ! surface water depth [m]
  real,    intent(in)    :: veg_tempk     ! veg temp [K]
  real,    intent(in)    :: veg_fracarea  ! veg fractional area coverage
  real,    intent(in)    :: veg_height    ! veg height [m]
  real,    intent(in)    :: veg_albedo    ! veg albedo
  real,    intent(inout) :: soil_energy   ! soil internal energy [J/m^3]
  real,    intent(in)    :: soil_water    ! soil water content [vol_water/vol_tot]
  real,    intent(in)    :: wresid_vg
  real,    intent(in)    :: wsat_vg
  real,    intent(in)    :: specifheat_drysoil
  real,    intent(in)    :: sand

  real,    intent(out)   :: snowfac       ! fractional veg burial by snowcover
  real,    intent(out)   :: vf            ! fractional coverage of non-buried part of veg
  real,    intent(out)   :: cosz          ! solar zenith angle for land cell
  real,    intent(out)   :: rlongup       ! upward sfc l/w rad flux [W/m^2]
  real,    intent(out)   :: rlong_albedo  ! albedo for upward sfc l/w
  real,    intent(out)   :: albedo_beam   ! land cell albedo (beam)
  real,    intent(out)   :: gnd_albedo    ! ground albedo
  real,    intent(out)   :: gnd_emiss     ! ground emiss
  real,    intent(out)   :: slong         ! surface l/w rad emission [W/m^2]
  real,    intent(out)   :: vlong         ! veg l/w rad emission [W/m^2]

  ! Local variables

  real :: sfcwater_epm       ! sfcwater energy per m^2 [J/m^2]
  real :: skn_tempk          ! surface skinlayer temp [K]
  real :: skn_fracliq        ! fraction of surface skinlayer water in liquid phase []
  real :: rlongv_a           ! longwave radiative flux from veg  to atm [W/m^2]
  real :: rlongs_a           ! longwave radiative flux from snow to atm [W/m^2]
  real :: soil_watfrac_ul    ! soil water fraction (unlimited) []
  real :: soil_watfrac       ! soil water fraction (limited) []
  real :: sfcwater_fracarea  ! fractional coverage of surface water []
  real :: vfc                ! 1 - vf []
  real :: fc50               ! minimum of soil water fraction and 50% []
  real :: alg                ! soil surface albedo []
  real :: alw                ! surface water albedo [] (liquid and frozen)
  real :: veg_emiss          ! veg emissivity []

  ! This routine is called by the radiation scheme to compute the net surface
  ! albedo and upward longwave radiative flux for each land cell.

  ! Compute solar incidence angle for land cells (accounts for local topo slope)

  cosz = wnxl * sunx + wnyl * suny + wnzl * sunz

  ! Evaluate surface skin temperature and albedo depending on skncomp value

  if (skncomp == 1) then
     sfcwater_epm = sfcwater_mass(1) * sfcwater_energy(1)
     call qtk(sfcwater_energy(1),skn_tempk,skn_fracliq)
     call sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy, &
                             specifheat_drysoil, sfcwater_mass(1), sfcwater_epm, &
                             skn_tempk, skn_fracliq)
  elseif (skncomp == 2) then
     call qtk(sfcwater_energy(1),skn_tempk,skn_fracliq)
  else
     call qtk(sfcwater_energy(2),skn_tempk,skn_fracliq)
  endif

  ! Bare ground (soil surface) albedo

  soil_watfrac_ul = (soil_water - wresid_vg) / (wsat_vg - wresid_vg)
  soil_watfrac = min(1.0,max(0.001,soil_watfrac_ul))
  fc50 = min(.50, soil_watfrac)

  if (leaf_class == 2) then
     alg = .80  ! Firn/glacier albedo
  elseif (sand > 0.7) then
     alg = .41 - .34 * fc50  ! sandy soil albedo
  else
     alg = .31 - .34 * fc50  ! non-sandy soil albedo
  endif

  if (sfcwater_mass(1) < wcap_min) then  ! Case with no surface water

     snowfac    = 0.
     gnd_emiss  = emisg
     gnd_albedo = alg

  else                                   ! Case with surface water

     ! Surface water albedo (frozen and/or liquid)

     if (leaf_class == 2) then
        alw = .80                    ! Firn/glacier albedo
     else
        alw = .5 - .36 * skn_fracliq ! alw = .5 for all ice, .14 for all liquid
     endif

     sfcwater_fracarea = sqrt(min(1.,0.5 * sfcwater_mass(1)))
     snowfac           = sum(sfcwater_depth(:)) / max(.001,veg_height)
     if (snowfac > .9) snowfac = 1.

     gnd_emiss  = sfcwater_fracarea * emisw + (1. - sfcwater_fracarea) * emisg
     gnd_albedo = sfcwater_fracarea * alw   + (1. - sfcwater_fracarea) * alg

  endif

  ! Total shortwave albedo

  vf  = veg_fracarea * (1. - snowfac)
  vfc = 1. - vf

  albedo_beam = vf * veg_albedo + vfc * vfc * gnd_albedo

  ! Longwave albedo (1. - emissivity)

  veg_emiss    = emisv(leaf_class)
  rlong_albedo = (vf * (1. - veg_emiss) + vfc * vfc * (1. - gnd_emiss))

  ! Outgoing longwave

  vlong = veg_emiss * stefan * veg_tempk**4
  slong = gnd_emiss * stefan * skn_tempk**4

  rlongv_a = vlong * vf * (2. - gnd_emiss - vf + gnd_emiss * vf)
  rlongs_a = slong * vfc

  rlongup  = rlongv_a + rlongs_a

  end subroutine sfcrad_prep

!===============================================================================

  subroutine sfcrad_rlongup(iland, iwsfc, leaf_class, skncomp, sfcwater_mass, &
                            sfcwater_energy, sfcwater_depth, veg_tempk,       &
                            soil_energy, soil_water, specifheat_drysoil,      &
                            gnd_emiss, vf, rlongup, slong, vlong)

  use leaf_coms,   only: emisv
  use mem_land,    only: nzs
  use consts_coms, only: stefan
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in)    :: iland         ! horizontal index of current land cell
  integer, intent(in)    :: iwsfc         ! horizontal index of current sfcgrid cell
  integer, intent(in)    :: leaf_class    ! leaf class
  integer, intent(in)    :: skncomp       ! skinlayer composition type (1, 2, or 3)
  real,    intent(in)    :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
  real,    intent(in)    :: sfcwater_energy(nzs) ! surface water internal energy [J/kg]
  real,    intent(in)    :: sfcwater_depth (nzs) ! surface water depth [m]
  real,    intent(in)    :: veg_tempk     ! veg temp [K]
  real,    intent(inout) :: soil_energy   ! soil internal energy [J/m^3]
  real,    intent(in)    :: soil_water    ! soil water content [vol_water/vol_tot]
  real,    intent(in)    :: specifheat_drysoil
  real,    intent(in)    :: gnd_emiss     ! ground emiss
  real,    intent(in)    :: vf            ! fractional coverage of non-buried part of veg
  real,    intent(out)   :: rlongup       ! upward sfc l/w rad flux [W/m^2]
  real,    intent(out)   :: slong         ! surface l/w rad emission [W/m^2]
  real,    intent(out)   :: vlong         ! veg l/w rad emission [W/m^2]

  ! Local variables

  real :: sfcwater_epm       ! sfcwater energy per m^2 [J/m^2]
  real :: skn_tempk          ! surface skinlayer temp [K]
  real :: skn_fracliq        ! fraction of surface skinlayer water in liquid phase []
  real :: rlongv_a           ! longwave radiative flux from veg  to atm [W/m^2]
  real :: rlongs_a           ! longwave radiative flux from snow to atm [W/m^2]
  real :: veg_emiss          ! veg emissivity []

  ! This routine is called by the radiation scheme to compute the net surface
  ! upward longwave radiative flux and skin temperature for each land cell.
  ! The surface emissivity is held constant between full radiation updates.

  ! Evaluate surface skin temperature and albedo depending on skncomp value

  if (skncomp == 1) then
     sfcwater_epm = sfcwater_mass(1) * sfcwater_energy(1)
     call qtk(sfcwater_energy(1),skn_tempk,skn_fracliq)
     call sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy, &
                             specifheat_drysoil, sfcwater_mass(1), sfcwater_epm, &
                             skn_tempk, skn_fracliq)
  elseif (skncomp == 2) then
     call qtk(sfcwater_energy(1),skn_tempk,skn_fracliq)
  else
     call qtk(sfcwater_energy(2),skn_tempk,skn_fracliq)
  endif

  ! Evaluate surface outgoing radiation

  veg_emiss = emisv(leaf_class)

  vlong = veg_emiss * stefan * veg_tempk**4
  slong = gnd_emiss * stefan * skn_tempk**4

  rlongv_a = vlong * vf * (2. - gnd_emiss - vf + gnd_emiss * vf)
  rlongs_a = slong * (1. - vf)

  rlongup  = rlongv_a + rlongs_a

  end subroutine sfcrad_rlongup


End Module leaf4_surface
