module sea_swm

  implicit none

  real, parameter :: depthmin_flux = 0.1   ! Minimum depth [m] for sfcwater flow
                                           ! across V edge of sfcgrid

  real, parameter :: depthmin_swe  = 1.0   ! Minimum depth [m] for application
                                           ! full SW equations

  real, parameter :: depthmax_swe  = 100.0

  real, parameter :: mc_land = 0.10       ! Manning coefficient (relatively rough surface)
  real, parameter :: mc_sea  = 0.03       ! Manning coefficient (relatively smooth surface)
  real, parameter :: mc_sea2 = mc_sea**2
  real, parameter :: minusonethird = -1./3.

  integer :: niter_swm
  real    :: dt_swm

  ! Definitions related to height/depth

  ! sfcg%topo(iw)     Topographic elevation [m] of land or water surface
  !                   relative to mean sea level (MSL).
  !
  ! sfcg%bathym(iw)   Topographic elevation [m] of land surface or bottom
  !                   surface of a lake, ocean, or permanent large river
  !                   relative to mean sea level.  sfcg%bathym is equal to
  !                   sfcg%topo for land areas, and is less than sfcg%topo
  !                   for permanent water bodies.
  !
  ! sea%swmdepth(iw)         Depth [m] of ocean water in SWM.
  ! land%sfcwater_depth(iw)  Depth [m] of surface water over land.
  ! lake%_depth(iw)          Depth [m] of lake water
  !
  ! sfcg%head1(iw)    Elevation of water surface [m] relative to local
  !                   topography sfcg%topm(iw).  In all locations,
  !                   sfcg%head1(iw) + sfcg%topo1(iw) = water_depth(iw) + sfcg%bathym(iw).
  !                   In land locations only, it is also true that
  !                   sfcg%head1(iw) = land%sfcwater_depth(iw).

contains

!===============================================================================

subroutine swm_driver()

  use misc_coms,    only: iparallel
  use mem_sfcg,     only: jtab_vsfc_swm, sfcg, jtab_wsfc_swm, itab_wsfc, itab_vsfc, &
                          nswmzons, mvsfc
  use leaf_coms,    only: dt_leaf, wcap_min
  use sea_coms,     only: dt_sea
  use mem_land,     only: land, mland, omland
  use mem_lake,     only: lake, mlake, omlake
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc, mpi_send_vsfc, mpi_recv_vsfc
  use mem_para,     only: myrank
  use oname_coms,   only: nl
  use umwm_module,  only: umwmflg
  use umwm_top,     only: umwm_step

  implicit none

  integer :: iter_swm, ivsfc, iw1, iw2, j, iwsfc, ilake, iland, ivn, iwn, isea
  real    :: factor, dheight, energyin, lake_epm2, dirv

  real, allocatable :: vmts(:)

  ! Could put call to umwm_step inside of iter_swm loop if necessary to use smaller timestep in umwm

  if (umwmflg == 1) call umwm_step()

  if (nswmzons < 1) return

  allocate(vmts(mvsfc)) ; vmts = 0.0

  call vort_damp_swm(vmts)

  call divh_damp_swm(vmts, dt_swm) ! This call is outside iter_swm loop: Should it use dt_sea?

  ! START MAIN SWM LOOP

  do iter_swm = 1,niter_swm

     ! Evaluate horizontal gradients of temp, VXE, VYE, and VZE for ocean cells
     ! that are updated with Shallow Water Model (SWM).

     if (nl%igw_spinup /= 1) call swm_grad2d()

     if (iparallel == 1) then

        ! MPI SEND/RECV sea%gxps_vxe, sea%gyps_vxe, sea%gxps_vye,
        !               sea%gyps_vye, sea%gxps_vze, sea%gyps_vze

        call mpi_send_wsfc(set='swm_grad')
        call mpi_recv_wsfc(set='swm_grad')

     endif

     ! Evaluate horizontal fluxes of temp, VXE, VYE, and VZE for all cells
     ! that are updated with Shallow Water Model (SWM).
     !$omp parallel do private(ivsfc, iw1, iw2)
     do j = 1,jtab_vsfc_swm%jend
        ivsfc = jtab_vsfc_swm%ivsfc(j)

        ! Skip this vsfc face if running in parallel and cell rank is not MYRANK

        if (itab_vsfc(ivsfc)%irank /= myrank) cycle

        iw1 = itab_vsfc(ivsfc)%iwn(1)
        iw2 = itab_vsfc(ivsfc)%iwn(2)

        ! LAND-LAND, LAND-LAKE, LAND-SEA, LAKE-LAND, SEA-LAND, SEA-SEA horiz SWM fluxes

        call swm_hflux(ivsfc, iw1, iw2)

     enddo ! ivsfc
     !$omp end parallel do

     if (iparallel == 1) then

        ! MPI SEND/RECV sfcg%hflux_wat, sfcg%hflux_enr, sfcg%hflux_vxe,
        !               sfcg%hflux_vye, sfcg%hflux_vze, sfcg%vmp, sfcg%vmc, sfcg%vc

        call mpi_send_vsfc(set='swm_hflux')
        call mpi_recv_vsfc(set='swm_hflux')

     endif

     ! Loop over SEA CELLS that use SWM; prognose WDEPTH, HEAD1, & W-cell horiz momentum

     !$omp parallel do private(iwsfc)
     do j = 1,jtab_wsfc_swm%jend; iwsfc = jtab_wsfc_swm%iwsfc(j)
        if (itab_wsfc(iwsfc)%irank /= myrank) cycle

        if (sfcg%swm_active(iwsfc)) then
           call swm_progw(iwsfc)      ! Active SWM cells
        else
           call swm_progw_lbc(iwsfc)  ! Boundary SWM cells
        endif
     enddo
     !$omp end parallel do

     ! MPI SEND/RECV of quantities needed for prog_v

     if (iparallel == 1) then
        call mpi_send_wsfc(set='swm_progw')
        call mpi_recv_wsfc(set='swm_progw')
     endif

     ! UPDATE SFCG%VMC AND SFCG%VC for sea cells that use SWM

     call swm_progv(vmts)

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

     ! Change in lake height from SWM fluxes

     !$omp parallel
     !$omp do private(iwsfc, dheight, energyin, j, ivn, iwn, dirv, factor, &
     !$omp            lake_epm2)
     do ilake = 2,mlake
        iwsfc = ilake + omlake

        ! Skip this cell if running in parallel and cell rank is not MYRANK

        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle
        if (.not. sfcg%swm_active(iwsfc)) cycle

        dheight  = 0.
        energyin = 0.

        ! Loop over all lateral faces of this ilake cell

        do j = 1,itab_wsfc(iwsfc)%npoly
           ivn = itab_wsfc(iwsfc)%ivn(j)
           iwn = itab_wsfc(iwsfc)%iwn(j)

           if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
              dirv = 1.0
           else
              dirv = -1.0
           endif

!           if (sfcg%leaf_class(iwn) >= 2) then

              ! Neighbor cell is land

              ! Add mass and energy contributions from soil sfcwater in neighbor cell

              factor = dirv * dt_swm / sfcg%area(iwsfc) ! [s/m^2]

              dheight  = dheight  + factor * sfcg%hflux_wat(ivn) ! [m] or [m^3/m^2]
              energyin = energyin + factor * sfcg%hflux_enr(ivn) ! [J/m^2]

!           endif

        enddo ! j

        ! Add height and energy changes to cell

        lake_epm2 = lake%lake_energy(ilake) * lake%depth(ilake) * 1000.

        lake%depth(ilake) = lake%depth(ilake) + dheight

        lake%lake_energy(ilake) = (lake_epm2 + energyin) &
                                / (max(wcap_min,lake%depth(ilake)) * 1000.) ! water density = 1000 kg/m^3

        sfcg%head1(iwsfc) = lake%depth(ilake) + sfcg%bathym(iwsfc) - sfcg%topw(iwsfc)
     enddo
     !$omp end do

     ! Change in LAND sfcwater height from SWM fluxes

     !$omp do private(iwsfc,lake_epm2,j,ivn,dirv,factor) schedule(guided)
     do iland = 2,mland
        iwsfc = iland + omland

        ! Skip if not swm active
        if (.not. sfcg%swm_active(iwsfc)) cycle

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        ! Apply horizontal groundwater mass and energy fluxes

        lake_epm2 = land%sfcwater_energy(1,iland) * land%sfcwater_mass(1,iland)

        ! Loop over all lateral faces of this iland cell

        do j = 1,itab_wsfc(iwsfc)%npoly
           ivn = itab_wsfc(iwsfc)%ivn(j)

           if (iwsfc == itab_vsfc(ivn)%iwn(2)) then
              dirv = 1.0
           else
              dirv = -1.0
           endif

           if (nl%igw_spinup /= 1) then

              ! Add mass and energy contributions from overland flow from neighbor cell

              factor = dirv * dt_swm / sfcg%area(iwsfc) ! [s/m^2]

              land%sfcwater_depth(1,iland) = land%sfcwater_depth(1,iland)  &
                                           + factor * sfcg%hflux_wat(ivn) ! [m] or [m^3/m^2]

              land%sfcwater_mass(1,iland)  = land%sfcwater_mass(1,iland)   &
                                           + factor * sfcg%hflux_wat(ivn) * 1000. ! [kg/m^2]

              lake_epm2                    = lake_epm2 &
                                           + factor * sfcg%hflux_enr(ivn) ! [J/m^2]

           endif

        enddo ! j, ivn

        land%sfcwater_energy(1,iland) = lake_epm2 &
                                      / max(wcap_min,land%sfcwater_mass(1,iland))

        if (land%sfcwater_mass  (1,iland) < wcap_min) then
            land%sfcwater_mass  (:,iland) = 0.
            land%sfcwater_energy(:,iland) = 0.
            land%sfcwater_depth (:,iland) = 0.
        endif

     enddo
     !$omp end do
     !$omp end parallel

  enddo ! iter_swm

end subroutine swm_driver

!===============================================================================

subroutine swm_grad2d()

  use mem_sfcg, only: itab_wsfc, jtab_wsfc_swm, sfcg
  use mem_sea,  only: sea, omsea
  use mem_para, only: myrank

  implicit none

  integer :: j, iw, isea, npoly, jv, iv, iwn, isean, jm1, jm2

  real :: gxps_coef, gyps_coef
  real :: dvxe, dvye, dvze

  !$omp parallel do private(iw,isea,npoly,jv,iv,iwn,isean,jm1,jm2, &
  !$omp                     gxps_coef,gyps_coef,dvxe,dvye,dvze)
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     isea = iw - omsea

     if (itab_wsfc(iw)%irank /= myrank) cycle

     sea%gxps_vxe(isea) = 0.
     sea%gxps_vye(isea) = 0.
     sea%gxps_vze(isea) = 0.

     sea%gyps_vxe(isea) = 0.
     sea%gyps_vye(isea) = 0.
     sea%gyps_vze(isea) = 0.

     if (sea%swmdepth(isea) < depthmin_swe) cycle

     npoly = itab_wsfc(iw)%npoly

     ! Loop over V/W neighbors of this W cell

     do jv = 1, npoly
        iv  = itab_wsfc(iw)%ivn(jv)
        iwn = itab_wsfc(iw)%iwn(jv)

        jm2 = jv
        if (jm2 == 1) then
           jm1 = npoly
        else
           jm1 = jm2 - 1
        endif

        gxps_coef = itab_wsfc(iw)%gxps1(jm2) + itab_wsfc(iw)%gxps2(jm1)
        gyps_coef = itab_wsfc(iw)%gyps1(jm2) + itab_wsfc(iw)%gyps2(jm1)

        dvxe = 0.  ! ASSUME ZERO GRADIENT IF IWN NOT FILLED SWM CELL
        dvye = 0.
        dvze = 0.

        if ( sfcg%leaf_class(iwn) == 0 .and. &
             sfcg%swm_active(iwn) ) then

           isean = iwn - omsea

           dvxe = sea%vxe(isean) - sea%vxe(isea)
           dvye = sea%vye(isean) - sea%vye(isea)
           dvze = sea%vze(isean) - sea%vze(isea)
        endif

        sea%gxps_vxe(isea) = sea%gxps_vxe(isea) + gxps_coef * dvxe
        sea%gxps_vye(isea) = sea%gxps_vye(isea) + gxps_coef * dvye
        sea%gxps_vze(isea) = sea%gxps_vze(isea) + gxps_coef * dvze

        sea%gyps_vxe(isea) = sea%gyps_vxe(isea) + gyps_coef * dvxe
        sea%gyps_vye(isea) = sea%gyps_vye(isea) + gyps_coef * dvye
        sea%gyps_vze(isea) = sea%gyps_vze(isea) + gyps_coef * dvze

     enddo

  enddo
  !$omp end parallel do

end subroutine swm_grad2d

!=========================================================================

subroutine swm_hflux(iv, iw1, iw2)

  use mem_sfcg,    only: sfcg, itab_vsfc
  use mem_sea,     only: sea, omsea
  use mem_lake,    only: lake, omlake
  use mem_land,    only: land, omland
  use consts_coms, only: cliq1000, alli1000, grav
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in) :: iv, iw1, iw2

  integer :: isea1, isea2, iland1, iland2, ilake1, ilake2

  real :: cosv1, sinv1, cosv2, sinv2
  real :: dxps_v, dyps_v

  real :: ufacev, vmcfman, vmcfpgf, vmcf, vcf, vdepth, slope, mc1, mc2
  real :: wdepth1, wdepth2, tempk1, tempk2, fracliq1, fracliq2

  sfcg%hflux_wat(iv) = 0.
  sfcg%hflux_enr(iv) = 0.
  sfcg%hflux_vxe(iv) = 0.
  sfcg%hflux_vye(iv) = 0.
  sfcg%hflux_vze(iv) = 0.

  ! No ocean-to-lake or lake-to-lake fluxes (lake-to-lake done elsewhere),
  ! and no fluxes to/from non-swm_active cells unless they are sea cells

!  if ((sfcg%leaf_class(iw1) <  2  .and. sfcg%leaf_class(iw2) == 1) .or. &
!      (sfcg%leaf_class(iw1) == 1  .and. sfcg%leaf_class(iw2) <  2) .or. &
!      (.not. sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) >  0) .or. &
!      (.not. sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) >  0)) then

  if ( ( sfcg%leaf_class(iw1) == 1 .and. sfcg%leaf_class(iw2) == 1 &
              .and. (abs(sfcg%topw(iw2) - sfcg%topw(iw1)) < .02) ) .or. &
      (.not. sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) >  0) .or. &
      (.not. sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) >  0)) then

     sfcg%vmp(iv) = 0.
     sfcg%vmc(iv) = 0.
     sfcg%vc (iv) = 0.

     return
  endif

  if (sfcg%leaf_class(iw1) == 0) then
     isea1 = iw1 - omsea
     wdepth1 = sea%swmdepth(isea1)
     mc1 = mc_sea
     tempk1 = sea%seatc(isea1)
     fracliq1 = 1.0
  elseif (sfcg%leaf_class(iw1) == 1) then
     ilake1 = iw1 - omlake
!    wdepth1 = lake%depth(ilake1)
     wdepth1 = lake%depth(ilake1) - sfcg%topw(iw1) + sfcg%bathym(iw1)
     mc1 = mc_sea
     call qtk(lake%lake_energy(ilake1), tempk1, fracliq1)
  else
     iland1 = iw1 - omland
     wdepth1 = land%sfcwater_depth(1,iland1)
     mc1 = mc_land
     call qtk(land%sfcwater_energy(1,iland1), tempk1, fracliq1)
  endif

  if (sfcg%leaf_class(iw2) == 0) then
     isea2 = iw2 - omsea
     wdepth2 = sea%swmdepth(isea2)
     mc2 = mc_sea
     tempk2 = sea%seatc(isea2)
     fracliq2 = 1.0
  elseif (sfcg%leaf_class(iw2) == 1) then
     ilake2 = iw2 - omlake
!    wdepth2 = lake%depth(ilake2)
     wdepth2 = lake%depth(ilake2) - sfcg%topw(iw2) + sfcg%bathym(iw2)
     mc2 = mc_sea
     call qtk(lake%lake_energy(ilake2), tempk2, fracliq2)
  else
     iland2 = iw2 - omland
     wdepth2 = land%sfcwater_depth(1,iland2)
     mc2 = mc_land
     call qtk(land%sfcwater_energy(1,iland2), tempk2, fracliq2)
  endif

  vdepth = 0.5 * (wdepth1 + wdepth2)
  slope = sfcg%dniv(iv) * (sfcg%head1(iw1) + sfcg%topw(iw1) &
                         - sfcg%head1(iw2) - sfcg%topw(iw2))

  if (slope >= 0. .and. wdepth1 <= depthmin_flux .or. &
      slope <= 0. .and. wdepth2 <= depthmin_flux) then

     ! If donor cell wdepth < depthmin_flux, prevent flux across IV

     sfcg%vmp(iv) = 0.
     sfcg%vmc(iv) = 0.
     sfcg%vc (iv) = 0.

     return

  elseif (wdepth1 < depthmin_swe   .or. &
          wdepth2 < depthmin_swe   .or. &
          sfcg%leaf_class(iw1) > 1 .or. &
          sfcg%leaf_class(iw2) > 1) then

     ! Donor cell has sufficient water depth for outflow.  If at least one cell
     ! does not have sufficient water for full SWE model or is land, maximum flux
     ! across IV is based on Manning's equation applied to steady state flow over
     ! a sloping rough surface.  However, also take into account the inertia of
     ! the flux and the maximum amount that it can be changed from its previous
     ! value by the PGF.

     vmcfman = 2. / (mc1 + mc2) * vdepth**(5./3.) * sign(sqrt(abs(slope)), slope)

     vmcfpgf = sfcg%vmc(iv) + dt_swm * vdepth * grav * slope

     if (slope > 0.) then
        vmcf = max(0.,min(vmcfman, vmcfpgf))
     else
        vmcf = min(0.,max(vmcfman, vmcfpgf))
     endif

     sfcg%vmp(iv) = vmcf
     sfcg%vmc(iv) = vmcf
     sfcg%vc (iv) = vmcf / max(depthmin_flux,vdepth)

  else

     ! Both adjacent cells are sea cells with depth > depthmin_swe.

!!     if (sfcg%swm_active(iw1) .and. &
!!         sfcg%swm_active(iw2)) then

        ! If both cells are swm_active, use full SWE model, which
        ! extrapolates VM to time T + 1/2.

        vmcf         = 1.5 * sfcg%vmc(iv) - 0.5 * sfcg%vmp(iv)
        sfcg%vmp(iv) = sfcg%vmc(iv)

!!     else
!!
!!        ! If one cell is not swm_active, prognose vmc solely from pgf (defined in progv)
!!
!!        vmcf = sfcg%vmc(iv) + dt_swm * slope * grav * vdepth
!!
!!        sfcg%vmc(iv) = vmcf
!!        sfcg%vmp(iv) = vmcf
!!        sfcg%vc(iv)  = vmcf / vdepth
!!
!!     endif

  endif

  vcf = vmcf / max(depthmin_swe,vdepth)
  sfcg%hflux_wat(iv) = vmcf * sfcg%dnu(iv)

  ! Compute X, Y displacement components for V face relative to T point

  if (vcf > 0.0) then
     sfcg%hflux_enr(iv) = sfcg%hflux_wat(iv) &
                        * (cliq1000 * (tempk1 - 273.15) + alli1000) ! [J/s]

     if (sfcg%leaf_class(iw1) == 0 .and. sfcg%swm_active(iw1)) then

        if (sfcg%leaf_class(iw2) == 0) then
           ufacev = .5 * (sfcg%unx(iv) * (sea%vxe(isea1) + sea%vxe(isea2)) &
                        + sfcg%uny(iv) * (sea%vye(isea1) + sea%vye(isea2)) &
                        + sfcg%unz(iv) * (sea%vze(isea1) + sea%vze(isea2)))
        else
           ufacev = .5 * (sfcg%unx(iv) * sea%vxe(isea1) &
                        + sfcg%uny(iv) * sea%vye(isea1) &
                        + sfcg%unz(iv) * sea%vze(isea1))
        endif

        cosv1 = itab_vsfc(iv)%cosv(1)
        sinv1 = itab_vsfc(iv)%sinv(1)

        dxps_v = -0.5 * dt_swm * (vcf * cosv1 - ufacev * sinv1) + itab_vsfc(iv)%dxps(1)
        dyps_v = -0.5 * dt_swm * (vcf * sinv1 + ufacev * cosv1) + itab_vsfc(iv)%dyps(1)

        sfcg%hflux_vxe(iv) = sfcg%hflux_wat(iv) * (sea%vxe(isea1) &
                                   + dxps_v * sea%gxps_vxe(isea1) &
                                   + dyps_v * sea%gyps_vxe(isea1))

        sfcg%hflux_vye(iv) = sfcg%hflux_wat(iv) * (sea%vye(isea1) &
                                   + dxps_v * sea%gxps_vye(isea1) &
                                   + dyps_v * sea%gyps_vye(isea1))

        sfcg%hflux_vze(iv) = sfcg%hflux_wat(iv) * (sea%vze(isea1) &
                                   + dxps_v * sea%gxps_vze(isea1) &
                                   + dyps_v * sea%gyps_vze(isea1))

     endif

  else
     sfcg%hflux_enr(iv) = sfcg%hflux_wat(iv) &
                        * (cliq1000 * (tempk2 - 273.15) + alli1000) ! [J/s]

     if (sfcg%leaf_class(iw2) == 0 .and. sfcg%swm_active(iw2)) then

        if (sfcg%leaf_class(iw1) == 0) then
           ufacev = .5 * (sfcg%unx(iv) * (sea%vxe(isea1) + sea%vxe(isea2)) &
                        + sfcg%uny(iv) * (sea%vye(isea1) + sea%vye(isea2)) &
                        + sfcg%unz(iv) * (sea%vze(isea1) + sea%vze(isea2)))
        else
           ufacev = .5 * (sfcg%unx(iv) * sea%vxe(isea2) &
                        + sfcg%uny(iv) * sea%vye(isea2) &
                        + sfcg%unz(iv) * sea%vze(isea2))
        endif

        cosv2 = itab_vsfc(iv)%cosv(2)
        sinv2 = itab_vsfc(iv)%sinv(2)

        dxps_v = -0.5 * dt_swm * (vcf * cosv2 - ufacev * sinv2) + itab_vsfc(iv)%dxps(2)
        dyps_v = -0.5 * dt_swm * (vcf * sinv2 + ufacev * cosv2) + itab_vsfc(iv)%dyps(2)

        sfcg%hflux_vxe(iv) = sfcg%hflux_wat(iv) * (sea%vxe(isea2) &
                                   + dxps_v * sea%gxps_vxe(isea2) &
                                   + dyps_v * sea%gyps_vxe(isea2))

        sfcg%hflux_vye(iv) = sfcg%hflux_wat(iv) * (sea%vye(isea2) &
                                   + dxps_v * sea%gxps_vye(isea2) &
                                   + dyps_v * sea%gyps_vye(isea2))

        sfcg%hflux_vze(iv) = sfcg%hflux_wat(iv) * (sea%vze(isea2) &
                                   + dxps_v * sea%gxps_vze(isea2) &
                                   + dyps_v * sea%gyps_vze(isea2))

     endif

  endif

end subroutine swm_hflux

!=========================================================================

subroutine swm_progw(iw)

  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, omsea
  use consts_coms, only: omega2, grav
  use misc_coms,   only: time8

  implicit none

  integer, intent(in) :: iw

  integer :: npoly, jv, iv, isea

  real :: areai, dirv
  real :: wdepthi, tide_tend

  real :: wd_tend   ! Water depth tendency [m/s]
  real :: delex_wd  ! Change in water depth [m]

  real :: vmxe1, vmye1, vmze1
  real :: speed         ! Water flow speed [m/s]
  real :: cdtop, cdbot  ! Top and bottom drag coefficients [m/s]

  real, parameter :: fp = .80  ! head1 forward bias coefficient
  real, parameter :: rhow = 1000. ! Density of liquid water [kg/m^3]

  real, parameter :: tide_amp    = 0.0     ! [m] for testing only
  real, parameter :: tide_period = 14400.  ! [s] for testing only

  isea = iw - omsea

  areai = 1.0 / sfcg%area(iw)
  wd_tend = 0.

  ! Current T cell depth-weighted velocity

  vmxe1 = sea%swmdepth(isea) * sea%vxe(isea)
  vmye1 = sea%swmdepth(isea) * sea%vye(isea)
  vmze1 = sea%swmdepth(isea) * sea%vze(isea)

  ! Evaluate momentum tendencies from top and bottom drag forces and Coriolis force

  ! CDTOP is the drag coefficient between wind and water at the top water surface
  ! and is based on vkmsfc computed in subroutine stars for surface wind stress.

  ! CDBOT is derived from the Manning formula for water flow over a sloping surface,
  ! combined with an equation for the downslope component of the gravity force that
  ! such a flow would have in the steady state, and with the form in which CDBOT is
  ! used with units of [m/s] in the shallow water model.  Algebraic manipulation
  ! yields CDBOT in terms of gravity, the Manning coefficient, the flow depth, and
  ! the flow speed.  This eliminates the slope from the formula, yielding a form
  ! that is assumed to be applicable to the more general inertial flow in the SWM.

  speed = sqrt(sea%vxe(isea)**2 + sea%vye(isea)**2 + sea%vze(isea)**2)

  cdtop = sfcg%vkmsfc(iw) / (sfcg%dzt_bot(iw) * rhow)
  cdbot = grav * mc_sea2 * speed * ( max(depthmin_flux,sea%swmdepth(isea)) )**minusonethird

  sea%vmxet(isea) = cdtop * (sea%windxe(isea) - sea%vxe(isea)) - cdbot * sea%vxe(isea) + omega2 * vmye1
  sea%vmyet(isea) = cdtop * (sea%windye(isea) - sea%vye(isea)) - cdbot * sea%vye(isea) - omega2 * vmxe1
  sea%vmzet(isea) = cdtop * (sea%windze(isea) - sea%vze(isea)) - cdbot * sea%vze(isea)

  npoly = itab_wsfc(iw)%npoly
  do jv = 1,npoly
     iv = itab_wsfc(iw)%ivn(jv)
     dirv = itab_wsfc(iw)%dirv(jv)

     if (jv == 1) then
        sea%vmxet_area(isea) = sfcg%hflux_vxe(iv) * dirv
        sea%vmyet_area(isea) = sfcg%hflux_vye(iv) * dirv
        sea%vmzet_area(isea) = sfcg%hflux_vze(iv) * dirv
     else
        sea%vmxet_area(isea) = sea%vmxet_area(isea) + sfcg%hflux_vxe(iv) * dirv
        sea%vmyet_area(isea) = sea%vmyet_area(isea) + sfcg%hflux_vye(iv) * dirv
        sea%vmzet_area(isea) = sea%vmzet_area(isea) + sfcg%hflux_vze(iv) * dirv
     endif

     wd_tend = wd_tend + sfcg%hflux_wat(iv) * dirv * areai

  enddo

  ! Explicit tendency for height and velocity

  tide_tend = (tide_amp / tide_period) * cos(real(time8) / tide_period)

  delex_wd = dt_swm * (wd_tend + tide_tend)

  ! Update head1 at forward-biased intermediate time t + fp for horiz pgf

! sfcg%head1(iw) = max(0., sea%swmdepth(isea) + fp * delex_wd) + sfcg%bathym(iw)
  sfcg%head1(iw) = max(0., sea%swmdepth(isea) + fp * delex_wd) - min(depthmax_swe, -sfcg%bathym(iw))

  ! Update wdepth from (t) to (t+1)

  sea%swmdepth(isea) = max(0., sea%swmdepth(isea) + delex_wd)

  ! IN CASE this IW cell has any closed V faces, estimate velocity in T cells
  ! at (t+1) by prognostic method.  This will be used only in cut cells to
  ! project onto closed faces for re-diagnosis of vxe,vye,vze

  wdepthi = 1.0 / max(depthmin_flux,sea%swmdepth(isea))

  sea%vxe1(isea) = (vmxe1 + dt_swm * (sea%vmxet(isea) + sea%vmxet_area(isea) * areai)) * wdepthi
  sea%vye1(isea) = (vmye1 + dt_swm * (sea%vmyet(isea) + sea%vmyet_area(isea) * areai)) * wdepthi
  sea%vze1(isea) = (vmze1 + dt_swm * (sea%vmzet(isea) + sea%vmzet_area(isea) * areai)) * wdepthi

end subroutine swm_progw

!=========================================================================

subroutine swm_progw_lbc(iw)

  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, omsea
  use consts_coms, only: omega2, grav
  use misc_coms,   only: time8

  implicit none

  integer, intent(in) :: iw

  integer :: npoly, jv, ivn, iwn, isea

  real :: areai, dirv
  real :: wdepthi, tide_tend, depth0, dtfact

  real :: wd_tend   ! Water depth tendency [m/s]
  real :: delex_wd  ! Change in water depth [m]

  real :: vmxe1, vmye1, vmze1
  real :: speed         ! Water flow speed [m/s]
  real :: cdtop, cdbot  ! Top and bottom drag coefficients [m/s]

  real, parameter :: fp = .80  ! head1 forward bias coefficient
! real, parameter :: rhow = 1000. ! Density of liquid water [kg/m^3]

!  real, parameter :: tide_amp    = 0.0     ! [m] for testing only
!  real, parameter :: tide_period = 14400.  ! [s] for testing only

  real, parameter :: relax = 1000.0

  isea = iw - omsea

  areai = 1.0 / sfcg%area(iw)
  wd_tend = 0.

  ! Flux tendencies with fully prognosed neighbors

  npoly = itab_wsfc(iw)%npoly
  do jv = 1,npoly

     iwn = itab_wsfc(iw)%iwn(jv)
     if (sfcg%swm_active(iwn)) then

        ivn = itab_wsfc(iw)%ivn(jv)
        dirv = itab_wsfc(iw)%dirv(jv)
        wd_tend = wd_tend + sfcg%hflux_wat(ivn) * dirv

     endif
  enddo

  wd_tend = wd_tend * areai

  ! Relaxation term

  depth0 = min(depthmax_swe, sfcg%topw(iw) - sfcg%bathym(iw))
  dtfact = min(0.5, depth0 / relax)

  wd_tend = wd_tend - dtfact * (sea%swmdepth(isea) - depth0)

  ! Depth change this timestep

  delex_wd = dt_swm * wd_tend

  ! Update head1 at forward-biased intermediate time t + fp for horiz pgf

! sfcg%head1(iw) = max(0., sea%swmdepth(isea) + fp * delex_wd) + sfcg%bathym(iw))
  sfcg%head1(iw) = max(0., sea%swmdepth(isea) + fp * delex_wd) - min(depthmax_swe, -sfcg%bathym(iw))

  ! Update wdepth from (t) to (t+1)

  sea%swmdepth(isea) = max(0., sea%swmdepth(isea) + delex_wd)

end subroutine swm_progw_lbc

!============================================================================

subroutine swm_progv(vmts)

  use mem_sfcg,    only: sfcg, jtab_vsfc_swm, jtab_msfc_swm, jtab_wsfc_swm, &
                         itab_vsfc, itab_msfc, itab_wsfc, mmsfc, mvsfc
  use mem_sea,     only: sea, msea, omsea
  use misc_coms,   only: iparallel
  use consts_coms, only: grav
  use olam_mpi_sfc,only: mpi_send_wsfc, mpi_recv_wsfc, mpi_send_vsfc, mpi_recv_vsfc
  use mem_para,    only: myrank

  implicit none

  real, intent(in) :: vmts(mvsfc)

  integer          :: iv, iw, iw1, iw2, j, jv, im, isea, isea1, isea2, iwn, isean
  integer          :: im1, im2
  real             :: vdepth, slope, pgf

  !$omp parallel do private(iv,iw1,iw2,isea1,isea2,im1,im2,vdepth,pgf)
  do j = 1,jtab_vsfc_swm%jend
     iv = jtab_vsfc_swm%ivsfc(j)

     if (itab_vsfc(iv)%irank /= myrank) cycle

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     ! Skip this IV point unless both adjacent IW cells are sea cells

     if (sfcg%leaf_class(iw1) > 0 .or. sfcg%leaf_class(iw2) > 0) cycle

     isea1 = iw1 - omsea
     isea2 = iw2 - omsea

     if (sea%swmdepth(isea1) < depthmin_swe .or. &
         sea%swmdepth(isea2) < depthmin_swe) then

        sfcg%vmp(iv) = 0.
        sfcg%vmc(iv) = 0.
        sfcg%vc(iv)  = 0.

        cycle
     endif

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     vdepth = 0.5 * (sea%swmdepth(isea1) + sea%swmdepth(isea2))

     ! The following "pressure gradient force" term has a form convenient for the
     ! shallow water equations with units of [m^2/s^2], which is acceleration
     ! times water depth.  When multiplied by the timestep [s], the result has units
     ! of [m^2/s] representing velocity times water depth, which is the same as vmc.
     ! vc has units of [m/s].

     slope = sfcg%dniv(iv) * (sfcg%head1(iw1) + sfcg%topw(iw1) &
                            - sfcg%head1(iw2) - sfcg%topw(iw2))

     pgf = vdepth * grav * slope

     ! Update VMC

     if (.not. sfcg%swm_active(iw2)) then

        sfcg%vmc(iv) = sfcg%vmc(iv) + dt_swm * (vmts(iv) + pgf &

             + sfcg%vnx(iv) * sea%vmxet(isea1) &
             + sfcg%vny(iv) * sea%vmyet(isea1) &
             + sfcg%vnz(iv) * sea%vmzet(isea1) &

             + ( sfcg%vnx(iv) * sea%vmxet_area(isea1) &
               + sfcg%vny(iv) * sea%vmyet_area(isea1) &
               + sfcg%vnz(iv) * sea%vmzet_area(isea1) &
               ) / sfcg%area(iw1) )

     elseif (.not. sfcg%swm_active(iw1)) then

        sfcg%vmc(iv) = sfcg%vmc(iv) + dt_swm * (vmts(iv) + pgf &

             + sfcg%vnx(iv) * sea%vmxet(isea2) &
             + sfcg%vny(iv) * sea%vmyet(isea2) &
             + sfcg%vnz(iv) * sea%vmzet(isea2) &

             + ( sfcg%vnx(iv) * sea%vmxet_area(isea2) &
               + sfcg%vny(iv) * sea%vmyet_area(isea2) &
               + sfcg%vnz(iv) * sea%vmzet_area(isea2) &
               ) / sfcg%area(iw2) )

     else

        sfcg%vmc(iv) = sfcg%vmc(iv) + dt_swm * (vmts(iv) + pgf &

             + ( sfcg%vnx(iv) * (sea%vmxet(isea1) + sea%vmxet(isea2)) &
             + sfcg%vny(iv) * (sea%vmyet(isea1) + sea%vmyet(isea2)) &
             + sfcg%vnz(iv) * (sea%vmzet(isea1) + sea%vmzet(isea2)) ) * 0.5 &

             + ( sfcg%vnx(iv) * (sea%vmxet_area(isea1) + sea%vmxet_area(isea2)) &
               + sfcg%vny(iv) * (sea%vmyet_area(isea1) + sea%vmyet_area(isea2)) &
               + sfcg%vnz(iv) * (sea%vmzet_area(isea1) + sea%vmzet_area(isea2)) &
               ) / (sfcg%area(iw1) + sfcg%area(iw2)) )

     endif

     sfcg%vc(iv) = sfcg%vmc(iv) / max(depthmin_flux,vdepth)

  enddo

end subroutine swm_progv

!===============================================================================

subroutine swm_diagvel()

  use mem_sfcg,    only: sfcg, jtab_wsfc_swm, itab_wsfc, nswmzons
  use mem_sea,     only: sea, omsea
  use mem_para,    only: myrank

  implicit none

  integer :: j, iw, isea, npoly, jv, iv, iwn, isean
  real    :: vcwall

  if (nswmzons < 1) return

  !$omp parallel do private(iw, isea, npoly, jv, iv, iwn, isean, vcwall)
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     isea = iw - omsea

     if (itab_wsfc(iw)%irank /= myrank) cycle

     sea%vxe(isea) = 0.
     sea%vye(isea) = 0.
     sea%vze(isea) = 0.

     if (sea%swmdepth(isea) < depthmin_swe) cycle

     if (.not. sfcg%swm_active(iw)) cycle

     ! Loop over V neighbors of this W cell

     npoly = itab_wsfc(iw)%npoly

     do jv = 1, npoly
        iv    = itab_wsfc(iw)%ivn(jv)
        iwn   = itab_wsfc(iw)%iwn(jv)

        if (sfcg%leaf_class(iwn) == 0 .and. sfcg%swm_active(iwn)) then
           isean = iwn - omsea

           if (sea%swmdepth(isean) >= depthmin_swe) then

              sea%vxe(isea) = sea%vxe(isea) + itab_wsfc(iw)%ecvec_vx(jv) * sfcg%vc(iv)
              sea%vye(isea) = sea%vye(isea) + itab_wsfc(iw)%ecvec_vy(jv) * sfcg%vc(iv)
              sea%vze(isea) = sea%vze(isea) + itab_wsfc(iw)%ecvec_vz(jv) * sfcg%vc(iv)

           else

              vcwall = sfcg%vnx(iv) * sea%vxe1(isea) &
                     + sfcg%vny(iv) * sea%vye1(isea) &
                     + sfcg%vnz(iv) * sea%vze1(isea)

              sea%vxe(isea) = sea%vxe(isea) + itab_wsfc(iw)%ecvec_vx(jv) * vcwall
              sea%vye(isea) = sea%vye(isea) + itab_wsfc(iw)%ecvec_vy(jv) * vcwall
              sea%vze(isea) = sea%vze(isea) + itab_wsfc(iw)%ecvec_vz(jv) * vcwall

           endif

        else

           vcwall = sfcg%vnx(iv) * sea%vxe1(isea) &
                  + sfcg%vny(iv) * sea%vye1(isea) &
                  + sfcg%vnz(iv) * sea%vze1(isea)

           sea%vxe(isea) = sea%vxe(isea) + itab_wsfc(iw)%ecvec_vx(jv) * vcwall
           sea%vye(isea) = sea%vye(isea) + itab_wsfc(iw)%ecvec_vy(jv) * vcwall
           sea%vze(isea) = sea%vze(isea) + itab_wsfc(iw)%ecvec_vz(jv) * vcwall

        endif

     enddo

  enddo
  !$omp end parallel do

end subroutine swm_diagvel

!===============================================================================

subroutine swm_init()

  use mem_sfcg, only: sfcg, jtab_wsfc_swm, itab_wsfc, nswmzons
  use mem_sea,  only: sea, omsea

  implicit none

  integer :: j, iwsfc, isea, npoly, jv, iv
  real    :: tide

  if (nswmzons < 1) return

  tide = 0.

  !$omp parallel do private(iwsfc,isea,npoly,jv,iv)
  do j = 1,jtab_wsfc_swm%jend
     iwsfc = jtab_wsfc_swm%iwsfc(j)
     isea = iwsfc - omsea

     sea%vxe1(isea) = 0.
     sea%vye1(isea) = 0.
     sea%vze1(isea) = 0.

!    sea%swmdepth(isea) = max(depthmin_flux, tide - sfcg%bathym(iwsfc))
     sea%swmdepth(isea) = tide + min(depthmax_swe, -sfcg%bathym(iwsfc))

     sfcg%head1(iwsfc) = sea%swmdepth(isea) - min(depthmax_swe, -sfcg%bathym(iwsfc))

     if (sea%swmdepth(isea) < depthmin_swe) cycle

     ! Loop over adjacent V faces, regardless of whether prognosed or underground,
     ! to diagnose initial vxe1,vye2,vze1 from initial vc.

     npoly = itab_wsfc(iwsfc)%npoly
     do jv = 1, npoly
        iv = itab_wsfc(iwsfc)%ivn(jv)

        sea%vxe1(isea) = sea%vxe1(isea) + itab_wsfc(iwsfc)%ecvec_vx(jv) * sfcg%vc(iv)
        sea%vye1(isea) = sea%vye1(isea) + itab_wsfc(iwsfc)%ecvec_vy(jv) * sfcg%vc(iv)
        sea%vze1(isea) = sea%vze1(isea) + itab_wsfc(iwsfc)%ecvec_vz(jv) * sfcg%vc(iv)
     enddo

  enddo
  !$omp end parallel do

end subroutine swm_init

!===============================================================================

subroutine vort_damp_swm(vmts)

  use olam_mpi_sfc, only: mpi_send_msfc, mpi_recv_msfc
  use mem_sfcg,     only: sfcg, jtab_vsfc_swm, jtab_msfc_swm, &
                          itab_vsfc, itab_msfc, itab_wsfc, mmsfc, mvsfc
  use mem_sea,      only: sea, omsea
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
  use oname_coms,   only: nl

  implicit none

  real, intent(inout) :: vmts(mvsfc)

  integer                 :: iv, iw1, iw2, k, j, jv, im, kb
  integer                 :: iv1, iv2, iv3
  integer                 :: im1, im2, im3, isea1, isea2
  real                    :: del_vort(mmsfc)
  real                    :: arm0i, cv, vbar, dnui, vdepth
  real, allocatable, save :: cm(:)
  logical,           save :: firstime = .true.

  real, parameter         :: onethird = 1./3.
  real, parameter         :: twothird = 2./3.

  if (nl%akmin_vort < 1.e-7) return

  if (firstime) then
     firstime = .false.

     allocate(cm(mmsfc))

     cv = nl%akmin_vort  * 0.075

     !$omp parallel do private(im,iv1,iv2,iv3)
     do j = 1,jtab_msfc_swm%jend; im = jtab_msfc_swm%imsfc(j)

        if (itab_msfc(im)%irank /= myrank) cycle

        iv1  = itab_msfc(im)%ivn(1)
        iv2  = itab_msfc(im)%ivn(2)
        iv3  = itab_msfc(im)%ivn(3)

        cm(im) = cv * sfcg%arm0(im)**twothird * sfcg%dnu(iv1) * sfcg%dnu(iv2) * sfcg%dnu(iv3) / &
                 ( sfcg%dnu(iv1) *sfcg%dnu(iv2) *sfcg%dnv(iv3) &
                 + sfcg%dnu(iv1) *sfcg%dnv(iv2) *sfcg%dnu(iv3) &
                 + sfcg%dnv(iv1) *sfcg%dnu(iv2) *sfcg%dnu(iv3) )

     enddo
     !$omp end parallel do

  endif  ! firstime

  sfcg%vort(:) = 0.

  !$omp parallel do private(im,jv,iv)
  do j = 1,jtab_msfc_swm%jend; im = jtab_msfc_swm%imsfc(j)

     if (itab_msfc(im)%irank /= myrank) cycle

!    if ( any( sea%swmdepth( itab_msfc(im)%iwn(1:3) ) < depthmin_swe ) ) cycle

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1, 3
        iv = itab_msfc(im)%ivn(jv)

        if (itab_vsfc(iv)%imn(2) == im) then

           sfcg%vort(im) = sfcg%vort(im) + sfcg%vc(iv) * sfcg%dnv(iv)

        else

           sfcg%vort(im) = sfcg%vort(im) - sfcg%vc(iv) * sfcg%dnv(iv)

        endif
     enddo

     ! Convert circulation to relative vertical vorticity at M
     ! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
     ! divide their quotient by zfact)

     sfcg%vort(im) = sfcg%vort(im) / sfcg%arm0(im)

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_msfc(vort=sfcg%vort)
     call mpi_recv_msfc(vort=sfcg%vort)
  endif

  del_vort(:) = 0.0

  !$omp parallel do private(im,im1,im2,im3,vbar)
  do j = 1,jtab_msfc_swm%jend; im = jtab_msfc_swm%imsfc(j)

     if (itab_msfc(im)%irank /= myrank) cycle

!    if ( any( sea%swmdepth( itab_msfc(im)%iwn(1:3) ) < depthmin_swe ) ) cycle

     im1  = itab_msfc(im)%imn(1)
     im2  = itab_msfc(im)%imn(2)
     im3  = itab_msfc(im)%imn(3)

     vbar = onethird * (sfcg%vort(im1) + sfcg%vort(im2) + sfcg%vort(im3))
     del_vort(im) = cm(im) * (vbar - sfcg%vort(im))

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_msfc(vort=del_vort)
     call mpi_recv_msfc(vort=del_vort)
  endif

  !$omp parallel do private(iv,iw1,iw2,isea1,isea2,im1,im2)
  do j = 1,jtab_vsfc_swm%jend; iv = jtab_vsfc_swm%ivsfc(j)

     if (itab_vsfc(iv)%irank /= myrank) cycle

     iw1 = itab_vsfc(iv)%iwn(1); iw2 = itab_vsfc(iv)%iwn(2)

    ! Skip this IV point unless both adjacent IW cells are sea cells and swm_active

     if (sfcg%leaf_class(iw1) > 0 .or. sfcg%leaf_class(iw2) > 0) cycle
     if (.not. (sfcg%swm_active(iw1) .and. sfcg%swm_active(iw2))) cycle

     isea1 = iw1 - omsea
     isea2 = iw2 - omsea

     if ( sea%swmdepth(isea1) < depthmin_swe .or. &
          sea%swmdepth(isea2) < depthmin_swe ) cycle

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     vdepth = 0.5 * (sea%swmdepth(isea1) + sea%swmdepth(isea2))

     ! Horizontal filter for vertical vorticity

     vmts(iv) = vmts(iv) + vdepth * (del_vort(im2) - del_vort(im1)) / sfcg%dnu(iv)

  enddo
  !$omp end parallel do

end subroutine vort_damp_swm

!===============================================================================

subroutine divh_damp_swm(vmts, dtsm)

  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc
  use mem_sfcg,     only: sfcg, jtab_vsfc_swm, jtab_wsfc_swm, &
                          itab_vsfc, itab_wsfc, mwsfc, mvsfc
  use mem_sea,      only: omsea, sea
  use misc_coms,    only: iparallel
  use oname_coms,   only: nl
  use mem_para,     only: myrank

  implicit none

  real,    intent(inout)  :: vmts(mvsfc)
  real,    intent(in)     :: dtsm

  integer                 :: iv, iw, iw1, iw2, k, j, jv, ivn, iwn
  integer                 :: isea, isean, isea1, isea2
  real                    :: fact1, fact2
  real                    :: div2d(mwsfc)
  real                    :: del2d(mwsfc)

  real, allocatable, save :: fdiv(:), lplfct(:,:)
  logical,           save :: firstime = .true.

  if (nl%divh_damp_fact < 1.e-7) return

  if (firstime) then
     firstime = .false.

     allocate(fdiv(mvsfc)) ; fdiv = 0.0

     fact2 = nl%divh_damp_fact**2 / real(dtsm)

     do j = 1,jtab_vsfc_swm%jend; iv = jtab_vsfc_swm%jend
        if (itab_vsfc(iv)%irank /= myrank) cycle
        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)
        fdiv(iv) = -fact2 * min(sfcg%area(iw1),sfcg%area(iw2))**2 / sfcg%dnv(iv)
     enddo

     allocate(lplfct(7,mwsfc)) ; lplfct(:,:) = 0.0

     do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
        if (itab_wsfc(iw)%irank /= myrank) cycle
        do jv = 1, itab_wsfc(iw)%npoly
           ivn = itab_wsfc(iw)%ivn(jv)
           lplfct(jv,iw) = sfcg%dnu(ivn) / (sfcg%dnv(ivn) * sfcg%area(iw))
        enddo
     enddo

  endif  ! firstime

  div2d = 0.0

  !$omp parallel do private(iw,isea,jv,iv)
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     if (itab_wsfc(iw)%irank /= myrank) cycle
     isea = iw - omsea

     do jv = 1, itab_wsfc(iw)%npoly
        iv = itab_wsfc(iw)%ivn(jv)
        div2d(iw) = div2d(iw) - itab_wsfc(iw)%dirv(jv) * sfcg%dnu(iv) * sfcg%vmc(iv)
     enddo

     div2d(iw) = div2d(iw) / (sfcg%area(iw) * max(depthmin_flux,sea%swmdepth(isea)))

  enddo
  !$omp end parallel do

  ! MPI send/recv of div2d
  if (iparallel == 1) then
     call mpi_send_wsfc(set='swm_div2d_ex', div2d_ex=div2d)
     call mpi_recv_wsfc(set='swm_div2d_ex', div2d_ex=div2d)
  endif

  del2d = 0.0

  !$omp parallel do private(iw,jv,ivn,iwn,isea,isean)
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     if (itab_wsfc(iw)%irank /= myrank) cycle
     isea = iw - omsea
     if (sea%swmdepth(isea) < depthmin_swe) cycle

     do jv = 1, itab_wsfc(iw)%npoly
        ivn = itab_wsfc(iw)%ivn(jv)
        iwn = itab_wsfc(iw)%iwn(jv)
        isean = iwn - omsea
        if (sfcg%swm_active(iwn) .and. (sea%swmdepth(isean) > depthmin_swe)) then
           del2d(iw) = del2d(iw) + (div2d(iwn) - div2d(iw)) * lplfct(jv,iw)
        endif
     enddo

  enddo
  !$end parallel do

  ! MPI send/recv of del2d
  if (iparallel == 1) then
     call mpi_send_wsfc(set='swm_div2d_ex', div2d_ex=del2d)
     call mpi_recv_wsfc(set='swm_div2d_ex', div2d_ex=del2d)
  endif

  !$omp parallel do private(iv,iw1,iw2,isea1,isea2)
  do j = 1,jtab_vsfc_swm%jend; iv = jtab_vsfc_swm%ivsfc(j)

     if (itab_vsfc(iv)%irank /= myrank) cycle

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     ! Skip this IV point unless both adjacent IW cells are sea cells and swm_active

     if (sfcg%leaf_class(iw1) > 0 .or. sfcg%leaf_class(iw2) > 0) cycle
     if (.not. (sfcg%swm_active(iw1) .and. sfcg%swm_active(iw2))) cycle

     isea1 = iw1 - omsea
     isea2 = iw2 - omsea

     if ( sea%swmdepth(isea1) < depthmin_swe .or. &
          sea%swmdepth(isea2) < depthmin_swe ) cycle

     ! Damp (remove) laplacian of horizontal divergence

     vmts(iv) = vmts(iv) + fdiv(iv) * (del2d(iw2) - del2d(iw1))

  enddo
  !$omp end parallel do

end subroutine divh_damp_swm

end module sea_swm
