subroutine surface_driver()

  use leaf_coms,     only: iupdndvi, s1900_ndvi, indvifile, nndvifiles
  use mem_ijtabs,    only: itab_w, jtab_w, jtw_prog
  use mem_basic,     only: theta, rr_v
  use mem_sfcg,      only: itab_wsfc, sfcg
  use mem_land,      only: land, mland, omland, nzg, slzt
  use mem_lake,      only: mlake
  use mem_sea,       only: msea
  use sea_coms,      only: iupdsst, s1900_sst, isstfile, nsstfiles, &
                           iupdseaice, s1900_seaice, iseaicefile, nseaicefiles
  use misc_coms,     only: s1900_sim, iparallel
  use mem_grid,      only: lpw, volti
  use consts_coms,   only: grav, eps_virt, vonk, cp
  use leaf4_soil,    only: head_column
  use sea_swm,       only: swm_driver
  use mem_turb,      only: akm_sfc, vkm_sfc, ustar, sfluxt, sfluxr, &
                           wstar, wtv0, pblh, moli, ustar_k, wtv0_k
  use mem_tend,      only: thilt, rr_wt
  use mem_para,      only: myrank
  use olam_mpi_sfc,  only: mpi_send_wsfc, mpi_recv_wsfc
  use oname_coms,    only: nl
  use landcells_mod, only: landcells

  implicit none

  ! Local arrays

  real :: head_slope  (nzg,mland)
  real :: soil_watfrac(nzg,mland)

  ! Local variables

  integer :: iwsfc, iland, ilake, isea
  integer :: j, iw, ks, kw, ka
  integer :: jsfc, jasfc

  real :: timefac_ndvi   ! fraction of elapsed time from past to future NDVI obs
  real :: timefac_sst    ! fraction of elapsed time from past to future SST obs
  real :: timefac_seaice ! fraction of elapsed time from past to future SEAICE obs

  real, parameter :: onethird = 1./3.

!!-----------------------------------------------------------------------------
!!  ADVANCE SWM BY DTLONG (POSSIBLY USING SUB-CYCLING)
!!-----------------------------------------------------------------------------

  if (nl%igw_spinup /= 1) call swm_driver()

!!-----------------------------------------------------------------------------
!!  LOOP OVER ALL LAND CELLS (IN THIS SUBDOMAIN) TO COMPUTE HEAD, HEAD_SLOPE,
!!  AND SOIL_WATFRAC
!!-----------------------------------------------------------------------------

  !$omp parallel do private(iwsfc)
  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Compute total head (relative to local topographic datum) based on soil
     ! moisture, pressure head, and depth in soil

     call head_column(nzg, iland,                &
                      slzt            (:),       &
                      land%soil_water (:,iland), &
                      land%wfrac_low  (:,iland), &
                      land%wresid_vg  (:,iland), &
                      land%wsat_vg    (:,iland), &
                      land%alpha_vg   (:,iland), &
                      land%en_vg      (:,iland), &
                      land%head       (:,iland), &
                      head_slope      (:,iland), &
                      soil_watfrac    (:,iland)  )

  enddo
  !$omp end parallel do

!!-----------------------------------------------------------------------------
!!  MPI SEND/RECV OF LAND%HEAD AND SOIL_WATFRAC FOR A PARALLEL RUN WITH
!!  HORIZONTAL GROUND WATER TRANSPORT
!!-----------------------------------------------------------------------------

  if (iparallel == 1 .and. nl%ihoriz_gndwater_transport == 1) then
     call mpi_send_wsfc(set='head',soil_watfrac=soil_watfrac)
     call mpi_recv_wsfc(set='head',soil_watfrac=soil_watfrac)
  endif

!!-----------------------------------------------------------------------------
!!  CALL DRIVER ROUTINE TO PERFORM HORIZONTAL WATER TRANSPORT
!!-----------------------------------------------------------------------------

  if (nl%ihoriz_gndwater_transport == 1) then
     call gndwater_transport(head_slope, soil_watfrac)
  endif

!!-----------------------------------------------------------------------------
!! Time interpolation factor for updating NDVI, SST, and SEA ICE
!!-----------------------------------------------------------------------------

  timefac_ndvi   = 0.0
  timefac_sst    = 0.0
  timefac_seaice = 0.0

  if (iupdndvi == 1 .and. nndvifiles > 1) then
     timefac_ndvi = (s1900_sim             - s1900_ndvi(indvifile-1))  &
                  / (s1900_ndvi(indvifile) - s1900_ndvi(indvifile-1))
  endif

  if (iupdsst == 1 .and. nsstfiles > 1) then
     timefac_sst = (s1900_sim           - s1900_sst(isstfile-1))  &
                 / (s1900_sst(isstfile) - s1900_sst(isstfile-1))
  endif

  if (iupdseaice == 1 .and. nseaicefiles > 1) then
     timefac_seaice = (s1900_sim                 - s1900_seaice(iseaicefile-1)) &
                    / (s1900_seaice(iseaicefile) - s1900_seaice(iseaicefile-1))
  endif

!!-----------------------------------------------------------------------------
!! Update surface temperature, ice, and canopy variables in all sea cells.
!! Update POM1D column variables in sea cells that are pom_active.
!!-----------------------------------------------------------------------------

  !$omp parallel
  !$omp do
  do isea = 2, msea

     call seacells(isea, timefac_sst,  timefac_seaice)

  enddo
  !$omp end do nowait

!!-----------------------------------------------------------------------------
!! Loop over ALL LAKE CELLS
!!-----------------------------------------------------------------------------

  !$omp do
  do ilake = 2, mlake

     call lakecells(ilake)

  enddo
  !$omp end do nowait

!!-----------------------------------------------------------------------------
!! Loop over ALL LAND CELLS
!!-----------------------------------------------------------------------------

  !$omp do
  do iland = 2, mland

     call landcells(iland, timefac_ndvi, head_slope(:,iland), soil_watfrac(:,iland))

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! NOTE: If the call to plot_pom gets uncommented, subroutine plot_pom needs
  ! to be customized so that the 4 IWSFC/ISEA points that it plots are actually
  ! in the pom_active region of the simulation that you are running.

  ! print*, 'calling plot_pom'
  ! call plot_pom()
  ! print*, 'returned from plot_pom'

!!-----------------------------------------------------------------------------
!! MPI SEND/RECV of quantities updated in seacells, lakecells, land cells
!!-----------------------------------------------------------------------------

  if (iparallel == 1) then
     call mpi_send_wsfc(set='sfc_driv_end')
     call mpi_recv_wsfc(set='sfc_driv_end')
  endif

  ! Reset to zero atmospheric surface exchange variables prior to summing over
  ! the contributing surface cells for each atmospheric column.

  ustar   = 0.
  wtv0    = 0.
  sfluxt  = 0.
  sfluxr  = 0.
  vkm_sfc = 0.
  akm_sfc = 0.
  ustar_k = 0.
  wtv0_k  = 0.

  ! Loop over ATM grid columns that are primary in this subdomain

  !$omp parallel do private(iw, jsfc, iwsfc, jasfc, kw, ka, ks)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ka = lpw(iw)

     ! Loop over all SFC grid cells that couple to this ATM grid column

     do jsfc = 1,itab_w(iw)%jsfc2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)
        ks = kw - ka + 1

        ! Values averaged over all surface overlaps of an atm column
        ustar  (iw) = ustar  (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%ustar (iwsfc)
        wtv0   (iw) = wtv0   (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%wthv  (iwsfc)
        sfluxt (iw) = sfluxt (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%sfluxt(iwsfc)
        sfluxr (iw) = sfluxr (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%sfluxr(iwsfc)
        vkm_sfc(iw) = vkm_sfc(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%vkmsfc(iwsfc)

        ! Values averaged by vertical level of an atm column
        ustar_k(ks,iw) = ustar_k(ks,iw) + itab_wsfc(iwsfc)%arcoarkw(jasfc) * sfcg%ustar(iwsfc)
        wtv0_k (ks,iw) = wtv0_k (ks,iw) + itab_wsfc(iwsfc)%arcoarkw(jasfc) * sfcg%wthv (iwsfc)

        ! Values summed to tendency arrays
        akm_sfc(ks,iw) = akm_sfc(ks,iw) + itab_wsfc(iwsfc)%arc(jasfc) * sfcg%vkmsfc(iwsfc)

        ! New calculation of thilt for sfluxt units [W m^-2]
        thilt  (kw,iw) = thilt  (kw,iw) + itab_wsfc(iwsfc)%arc(jasfc) * sfcg%sfluxt(iwsfc) &
                                          / (cp * sfcg%canexner(iwsfc)) * volti(kw,iw)
        rr_wt  (kw,iw) = rr_wt  (kw,iw) + itab_wsfc(iwsfc)%arc(jasfc) * sfcg%sfluxr(iwsfc) * volti(kw,iw)
!       rr_co2t(kw,iw) = rr_co2t(kw,iw) + itab_wsfc(iwsfc)%arc(jasfc) * sfcg%sfluxc(iwsfc) * volti(kw,iw)

     enddo

     ! Compute some derived surface quantities

     moli(iw) = - grav * vonk * wtv0(iw) /  &
          ( ustar(iw)**3 * theta(ka,iw) * (1.0 + eps_virt * rr_v(ka,iw)) )

     if (wtv0(iw) > 0.0) then
        wstar(iw) = (grav * pblh(iw) * wtv0(iw) / theta(ka,iw)) ** onethird
     else
        wstar(iw) = 0.0
     endif

  enddo
  !$omp end parallel do

  ! No MPI send/recv communication of ATM turbulent fluxes are required for
  ! averaging to SFC grid cells because the fluxes are directly computed on
  ! SFC grid cells.

end subroutine surface_driver

!=============================================================================

subroutine sfcg_avgatm()

  use mem_basic,    only: rho, press, theta, tair, rr_v, vxe, vye, vze
  use consts_coms,  only: rocp, p00i
  use mem_micro,    only: rr_c
  use mem_co2,      only: co2flag, rr_co2
  use misc_coms,    only: iparallel
  use mem_grid,     only: gdz_abov8
  use mem_para,     only: myrank
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc
  use mem_sfcg,     only: mwsfc, sfcg, itab_wsfc
  use mem_sea,      only: sea, omsea
  use mem_turb,     only: pblh

  implicit none

  integer :: iwsfc, j, iw, kw, isea
  real :: vels, psfc

  ! Average atmospheric quantities to each SFC grid cell location

  ! Loop over all SFC grid cells

  !$omp parallel do private(j,iw,kw,vels,psfc,isea)
  do iwsfc = 2, mwsfc

     ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     sfcg%vels    (iwsfc) = 0.
     sfcg%prss    (iwsfc) = 0.
     sfcg%rhos    (iwsfc) = 0.
     sfcg%airtemp (iwsfc) = 0.
     sfcg%airtheta(iwsfc) = 0.
     sfcg%airrrv  (iwsfc) = 0.
     sfcg%pblh    (iwsfc) = 0.

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
        sfcg%pblh    (iwsfc) = sfcg%pblh    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * pblh    (iw)

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
     sfcg%canexner(iwsfc)  = (sfcg%prss(iwsfc) * p00i) ** rocp

  enddo
  !$omp end parallel do

  ! MPI send/recv sfcg%vels, sfcg%prss, sfcg%rhos, sfcg%airtemp,
  ! sfcg%airtheta, sfcg%airrrv, sea%windxe, sea%windye, sea%windze

  if (iparallel == 1) then
     call mpi_send_wsfc(set='avgatm')
     call mpi_recv_wsfc(set='avgatm')
  endif

end subroutine sfcg_avgatm
