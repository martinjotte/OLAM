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

module sea_swm

  implicit none

  real, parameter :: depthmin_flux = 0.1  ! Minimum depth [m] for sfcwater flow
                                          ! across V edge of sfcgrid

  real, parameter :: depthmin_swe  = 1.0  ! Minimum depth [m] for application
                                          ! full SW equations

  real, parameter :: mc_land = 0.10       ! Manning coefficient (relatively rough surface)
  real, parameter :: mc_sea  = 0.03       ! Manning coefficient (relatively smooth surface)
  real, parameter :: mc_sea2 = mc_sea**2
  real, parameter :: minusonethird = -1./3.

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
  ! sea%wdepth(iw)           Depth [m] of ocean water in SWM.
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

subroutine swm_grad2d()

  use mem_sfcg, only: mwsfc, itab_wsfc, jtab_wsfc_swm, sfcg
  use mem_sea,  only: sea, msea, omsea

  implicit none

  integer :: j, iw, isea, npoly, jv, iv, iwn, isean, jm1, jm2

  real :: gxps_coef, gyps_coef
  real :: dvxe, dvye, dvze

  !$omp parallel do private(iw,isea,npoly,jv,iv,iwn,isean,jm1,jm2, &
  !$omp                     gxps_coef,gyps_coef,dvxe,dvye,dvze)
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     isea = iw - omsea

     sea%gxps_vxe(isea) = 0.
     sea%gxps_vye(isea) = 0.
     sea%gxps_vze(isea) = 0.

     sea%gyps_vxe(isea) = 0.
     sea%gyps_vye(isea) = 0.
     sea%gyps_vze(isea) = 0.

     if (sea%wdepth(isea) < depthmin_swe) cycle

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

        if (sfcg%leaf_class(iwn) == 0 .and. &
            sfcg%swm_active(iwn)) then

           isean = iwn - omsea

           if (sea%wdepth(isean) > depthmin_swe) then
              dvxe = sea%vxe(isean) - sea%vxe(isea)
              dvye = sea%vye(isean) - sea%vye(isea)
              dvze = sea%vze(isean) - sea%vze(isea)
           endif
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
  use misc_coms,   only: dtlm, dtsm, mdomain
  use consts_coms, only: eradi, cliq1000, alli1000, grav
  use oname_coms,  only: nl
  use leaf_coms,   only: dt_leaf
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in) :: iv, iw1, iw2

  integer :: isea1, isea2, iland1, iland2, ilake1, ilake2

  real :: cosv1, sinv1, cosv2, sinv2
  real :: dxps_v, dyps_v

  real :: ufacev, vmcf, vcf, vdepth, slope, mc1, mc2
  real :: wdepth1, wdepth2, tempk1, tempk2, fracliq1, fracliq2

  sfcg%hflux_wat(iv) = 0.
  sfcg%hflux_enr(iv) = 0.
  sfcg%hflux_vxe(iv) = 0.
  sfcg%hflux_vye(iv) = 0.
  sfcg%hflux_vze(iv) = 0.

  ! No ocean-to-lake or lake-to-lake fluxes (lake-to-lake done elsewhere),
  ! and no fluxes to/from non-swm_active cells unless they are sea cells

  if ((sfcg%leaf_class(iw1) <  2  .and. sfcg%leaf_class(iw2) == 1) .or. &
      (sfcg%leaf_class(iw1) == 1  .and. sfcg%leaf_class(iw2) <  2) .or. &
      (.not. sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) >  0) .or. &
      (.not. sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) >  0)) then

     sfcg%vmp(iv) = 0.
     sfcg%vmc(iv) = 0.
     sfcg%vc (iv) = 0.

     return
  endif

  if (sfcg%leaf_class(iw1) == 0) then
     isea1 = iw1 - omsea
     wdepth1 = sea%wdepth(isea1)
     mc1 = mc_sea
     tempk1 = sea%seatc(isea1)
     fracliq1 = 1.0
  elseif (sfcg%leaf_class(iw1) == 1) then
     ilake1 = iw1 - omlake
     wdepth1 = lake%depth(ilake1)
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
     wdepth2 = sea%wdepth(isea2)
     mc2 = mc_sea
     tempk2 = sea%seatc(isea2)
     fracliq2 = 1.0
  elseif (sfcg%leaf_class(iw2) == 1) then
     ilake2 = iw2 - omlake
     wdepth2 = lake%depth(ilake2)
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
     ! does not have sufficient water for full SWE model or is land, set flux
     ! across IV based on Manning's equation applied to steady state flow over
     ! a sloping rough surface.

     vmcf = 2. / (mc1 + mc2) * vdepth**(5./3.) * sign(sqrt(abs(slope)), slope)
     sfcg%vmp(iv) = vmcf
     sfcg%vmc(iv) = vmcf
     sfcg%vc (iv) = vmcf / max(depthmin_flux,vdepth)

  else

     ! Both adjacent cells are sea cells with depth > depthmin_swe.

     if (sfcg%swm_active(iw1) .and. &
         sfcg%swm_active(iw2)) then

        ! If both cells are swm_active, use full SWE model, which
        ! extrapolates VM to time T + 1/2.

        vmcf         = 1.5 * sfcg%vmc(iv) - 0.5 * sfcg%vmp(iv)
        sfcg%vmp(iv) = sfcg%vmc(iv)

     else

        ! If one cell is not swm_active, prognose vmc solely from pgf (defined in progv)

        vmcf = sfcg%vmc(iv) + dt_leaf * slope * grav * vdepth

        sfcg%vmc(iv) = vmcf
        sfcg%vmp(iv) = vmcf
        sfcg%vc(iv)  = vmcf / vdepth

     endif

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

        dxps_v = -0.5 * dt_leaf * (vcf * cosv1 - ufacev * sinv1) + itab_vsfc(iv)%dxps(1)
        dyps_v = -0.5 * dt_leaf * (vcf * sinv1 + ufacev * cosv1) + itab_vsfc(iv)%dyps(1)

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

        dxps_v = -0.5 * dt_leaf * (vcf * cosv2 - ufacev * sinv2) + itab_vsfc(iv)%dxps(2)
        dyps_v = -0.5 * dt_leaf * (vcf * sinv2 + ufacev * cosv2) + itab_vsfc(iv)%dyps(2)

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

  use mem_sfcg,    only: itab_wsfc, sfcg, itab_vsfc
  use mem_sea,     only: sea, msea, omsea
  use consts_coms, only: omega2, grav
  use misc_coms,   only: time8
  use leaf_coms,   only: dt_leaf

  use mem_ijtabs,  only: itab_w

  implicit none

  integer, intent(in) :: iw

  integer :: npoly, jv, iv, isea

  real :: areai, dirv, dirv_areai
  real :: vmt1, fact, wdepthi, tide_tend

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

  vmxe1 = sea%wdepth(isea) * sea%vxe(isea)
  vmye1 = sea%wdepth(isea) * sea%vye(isea)
  vmze1 = sea%wdepth(isea) * sea%vze(isea)

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
  cdbot = grav * mc_sea2 * speed * ( max(depthmin_flux,sea%wdepth(isea)) )**minusonethird

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

  delex_wd = dt_leaf * (wd_tend + tide_tend)

  ! Update wdepth from (t) to (t+1)

  sea%wdepth(isea) = max(0.,sea%wdepth(isea) + delex_wd)

  ! Update head1 at forward-biased intermediate time t + fp for horiz pgf

  sfcg%head1(iw) = sea%wdepth(isea) + fp * delex_wd + sfcg%bathym(iw)

  ! IN CASE this IW cell has any closed V faces, estimate velocity in T cells
  ! at (t+1) by prognostic method.  This will be used only in cut cells to
  ! project onto closed faces for re-diagnosis of vxe,vye,vze

  wdepthi = 1.0 / max(depthmin_flux,sea%wdepth(isea))

  sea%vxe1(isea) = (vmxe1 + dt_leaf * (sea%vmxet(isea) + sea%vmxet_area(isea) * areai)) * wdepthi
  sea%vye1(isea) = (vmye1 + dt_leaf * (sea%vmyet(isea) + sea%vmyet_area(isea) * areai)) * wdepthi
  sea%vze1(isea) = (vmze1 + dt_leaf * (sea%vmzet(isea) + sea%vmzet_area(isea) * areai)) * wdepthi

end subroutine swm_progw

!============================================================================

subroutine swm_progv()

  use mem_sfcg,    only: sfcg, jtab_vsfc_swm, jtab_msfc_swm, jtab_wsfc_swm, &
                         itab_vsfc, itab_msfc, itab_wsfc, mmsfc, mvsfc, mwsfc
  use mem_sea,     only: sea, msea, omsea
  use misc_coms,   only: iparallel
  use consts_coms, only: grav
  use leaf_coms,   only: dt_leaf
  use olam_mpi_sfc,only: mpi_send_wsfc, mpi_recv_wsfc, mpi_send_vsfc, mpi_recv_vsfc

  implicit none

  integer                 :: iv, iw, iw1, iw2, j, jv, im, isea, isea1, isea2, iwn, isean
  integer                 :: im1, im2, im3, im4, im5, im6
  real                    :: vortp_ex(mmsfc)
  real                    :: vc_ex   (mvsfc)
  real                    :: div2d_ex(msea)
  logical,           save :: firstime = .true.
  real, allocatable, save :: c1(:), c2(:)
  real,              save :: dti_leaf
  real                    :: cd, dnusum, f1, f2
  real                    :: vdepth, pgf, vmt_vortdamp, vmt_divdamp

  real, parameter         :: onethird = 1./3.
  real, parameter         :: vort_damp_swm = 0.01
  real, parameter         :: akmin_div2d_swm = 0.01

  ! Initialize c1 & c2 arrays on first call

  if (firstime) then
     firstime = .false.

     dti_leaf = 1.0 / dt_leaf

     allocate(c1(mvsfc))
     allocate(c2(mvsfc))

     cd = akmin_div2d_swm * 0.075
     
     !$omp parallel do private(iv,iw1,iw2)
     do j = 1,jtab_vsfc_swm%jend; iv = jtab_vsfc_swm%ivsfc(j)

        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)

        c1(iv) = cd * (sfcg%area(iw1)**onethird)**2 * sfcg%dniv(iv)
        c2(iv) = cd * (sfcg%area(iw2)**onethird)**2 * sfcg%dniv(iv)
     enddo
     !$omp end parallel do

  endif  ! firstime

  ! Compute excess V velocity and set vorticity = 0 at V endpoints

  vc_ex(:) = 0.  ! Could move this to progw for more selective zeroing

  !$omp parallel do private(iv, iw1, iw2, isea1, isea2)
  do j = 1,jtab_vsfc_swm%jend
     iv = jtab_vsfc_swm%ivsfc(j)

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     isea1 = iw1 - omsea
     isea2 = iw2 - omsea

     vortp_ex(im1) = 0.
     vortp_ex(im2) = 0.

     if (sea%wdepth(isea1) > depthmin_swe .and. &
         sea%wdepth(isea2) > depthmin_swe) then

        vc_ex(iv) = sfcg%vc(iv) &
                  - 0.5 * (sfcg%vnx(iv) * (sea%vxe(isea1) + sea%vxe(isea2)) &
                         + sfcg%vny(iv) * (sea%vye(isea1) + sea%vye(isea2)) &
                         + sfcg%vnz(iv) * (sea%vze(isea1) + sea%vze(isea2)))

     endif
  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_vsfc(vc_ex=vc_ex)
     call mpi_recv_vsfc(vc_ex=vc_ex)
  endif

  ! Compute vorticity for vorticity damping at all swm-active msfc points in
  ! this subdomain (vortp_ex does not require mpi send/recv)
  !

  !$omp parallel
  !$omp do private(im,jv,iv,iw,isea)
  do j = 1,jtab_msfc_swm%jend; im = jtab_msfc_swm%imsfc(j)

     ! Loop over V neighbors to evaluate circulation around M (at time T)

     do jv = 1, 3
        iv = itab_msfc(im)%ivn(jv)
        iw = itab_msfc(im)%iwn(jv)

        isea = iw - omsea

        if (sea%wdepth(isea) < depthmin_swe) then
           vortp_ex(im) = 0.
           exit
        endif

        if (itab_vsfc(iv)%imn(2) == im) then
           vortp_ex(im) = vortp_ex(im) + sfcg%dnv(iv) * vc_ex(iv)
        else
           vortp_ex(im) = vortp_ex(im) - sfcg%dnv(iv) * vc_ex(iv)
        endif
     enddo

  enddo
  !$omp end do

  ! Compute divergence for divergence damping

  !$omp do private(iw,isea,jv,iv,iwn,isean)
  do j = 1, jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     isea = iw - omsea

     div2d_ex(isea) = 0.0

     if (sea%wdepth(isea) < depthmin_swe) cycle

     do jv = 1, itab_wsfc(iw)%npoly
        iv    = itab_wsfc(iw)%ivn(jv)
        iwn   = itab_wsfc(iw)%iwn(jv)

        if (sfcg%leaf_class(iwn) == 0 .and. &
            sfcg%swm_active(iwn)) then

           isean = iwn - omsea

           if (sea%wdepth(isean) >= depthmin_swe) then
              div2d_ex(isea) = div2d_ex(isea) - itab_wsfc(iw)%dirv(jv) &
                             * vc_ex(iv) * sfcg%dnu(iv) / sfcg%area(iw)
           endif

        endif

     enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! MPI send/recv of div2d

  if (iparallel == 1) then
     call mpi_send_wsfc(set='swm_div2d_ex', div2d_ex=div2d_ex)
     call mpi_recv_wsfc(set='swm_div2d_ex', div2d_ex=div2d_ex)
  endif

  !$omp parallel do private(iv,iw1,iw2,isea1,isea2,im1,im2,vdepth, &
  !$omp                     vmt_vortdamp,vmt_divdamp,pgf)
  do j = 1,jtab_vsfc_swm%jend
     iv = jtab_vsfc_swm%ivsfc(j)

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     isea1 = iw1 - omsea
     isea2 = iw2 - omsea

     if (sea%wdepth(isea1) < depthmin_swe .or. &
         sea%wdepth(isea2) < depthmin_swe) then

        sfcg%vmp(iv) = 0.
        sfcg%vmc(iv) = 0.
        sfcg%vc(iv)  = 0.

        cycle
     endif

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     vdepth = 0.5 * (sea%wdepth(isea1) + sea%wdepth(isea2))

     vmt_vortdamp = vdepth * vort_damp_swm * sfcg%dniv(iv) * dti_leaf &
                  * (vortp_ex(im1) - vortp_ex(im2))

     vmt_divdamp = sea%wdepth(isea2) * c2(iv) * div2d_ex(isea2) &
                 - sea%wdepth(isea1) * c1(iv) * div2d_ex(isea1)

     ! The following "pressure gradient force" term has a form convenient for the
     ! shallow water equations with units of [m^2/s^2], which is acceleration
     ! times water depth.  When multiplied by the timestep [s], the result has units
     ! of [m^2/s] representing velocity times water depth, which is the same as vmc.
     ! vc has units of [m/s].

     pgf = vdepth * grav * sfcg%dniv(iv) &
         * (sfcg%head1(iw1) + sfcg%topw(iw1) - sfcg%head1(iw2) - sfcg%topw(iw2))

        ! Update VMC

     sfcg%vmc(iv) = sfcg%vmc(iv) + dt_leaf * (vmt_vortdamp + vmt_divdamp + pgf &

              + (sfcg%vnx(iv) * (sea%vmxet(isea1) + sea%vmxet(isea2))  &
              +  sfcg%vny(iv) * (sea%vmyet(isea1) + sea%vmyet(isea2))  &
              +  sfcg%vnz(iv) * (sea%vmzet(isea1) + sea%vmzet(isea2))) * 0.5 &

              + (sfcg%vnx(iv) * (sea%vmxet_area(isea1) + sea%vmxet_area(isea2))  &
              +  sfcg%vny(iv) * (sea%vmyet_area(isea1) + sea%vmyet_area(isea2))  &
              +  sfcg%vnz(iv) * (sea%vmzet_area(isea1) + sea%vmzet_area(isea2))) &
              / (sfcg%area(iw1) + sfcg%area(iw2)))

     sfcg%vc(iv) = sfcg%vmc(iv) / max(depthmin_flux,vdepth)

  enddo

end subroutine swm_progv

!===============================================================================

subroutine swm_diagvel()

  use mem_sfcg,    only: sfcg, jtab_wsfc_swm, itab_vsfc, itab_wsfc
  use mem_sea,     only: sea, omsea
  use misc_coms,   only: iparallel

  implicit none

  integer :: j, iw, isea, npoly, jv, iv, iwn, isean
  real :: vcwall

  !$omp parallel do private(iw, isea, npoly, jv, iv, iwn, isean, vcwall) 
  do j = 1,jtab_wsfc_swm%jend; iw = jtab_wsfc_swm%iwsfc(j)
     isea = iw - omsea

     sea%vxe(isea) = 0.
     sea%vye(isea) = 0.
     sea%vze(isea) = 0.

     if (sea%wdepth(isea) < depthmin_swe) cycle

     ! Loop over V neighbors of this W cell

     npoly = itab_wsfc(iw)%npoly

     do jv = 1, npoly
        iv    = itab_wsfc(iw)%ivn(jv)
        iwn   = itab_wsfc(iw)%iwn(jv)

        if (sfcg%leaf_class(iwn) == 0 .and. sfcg%swm_active(iwn)) then
           isean = iwn - omsea

           if (sea%wdepth(isean) >= depthmin_swe) then

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

  use mem_sfcg, only: sfcg, jtab_wsfc_swm, itab_vsfc, itab_wsfc
  use mem_sea,  only: sea, omsea

  implicit none

  integer :: j, iwsfc, isea, npoly, jv, iv
  real :: tide

  tide = 0.

  !$omp parallel do private(iwsfc,isea,npoly,jv,iv) 
  do j = 1,jtab_wsfc_swm%jend
     iwsfc = jtab_wsfc_swm%iwsfc(j)
     isea = iwsfc - omsea

     sea%vxe1(isea) = 0.
     sea%vye1(isea) = 0.
     sea%vze1(isea) = 0.

     sea%wdepth(isea) = max(depthmin_flux, tide - sfcg%bathym(iwsfc))

     sfcg%head1(iwsfc) = sea%wdepth(isea) + sfcg%bathym(iwsfc)

     if (sea%wdepth(isea) < depthmin_swe) cycle

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

end module sea_swm
