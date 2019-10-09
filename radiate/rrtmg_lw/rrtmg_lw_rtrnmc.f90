!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:25 $
!
      module rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

      use parkind,  only: im => kind_im, rb => kind_rb, cldmin
      use parrrtm,  only: mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, delwaveg, ngb, ngs, ngc, ngs1
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl

      implicit none

      real(kind=rb), parameter :: wtdiff = 0.5_rb
      real(kind=rb), parameter :: rec_6  = -1._rb / 6._rb
      real(kind=rb), parameter :: fluxfacw = fluxfac * wtdiff
      real(kind=rb), parameter :: delwavef(ngptlw) = delwaveg * fluxfacw

      private
      public :: rtrnmc_noclr

      contains

!-----------------------------------------------------------------------------

      subroutine rtrnmc_noclr( &
                        nlayers, nsfc, semiss, &
                        cldfr, taucmc, planklayb, dplankdnb, dplankupb, &
                        planksfcb, dplnksfcb, frac_sfck, &
                        pwvcm, fracs, taut, &
                        totuflux, totdflux, totuflux_sfc, totdflux_sfc )

!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that
!  use of the emissivity angle for the flux integration can cause errors of
!  1 to 4 W/m2 within cloudy layers.
!  Clouds are treated with the McICA stochastic approach and maximum-random
!  cloud overlap.
!  This subroutine also provides the optional capability to calculate
!  the derivative of upward flux respect to surface temperature using
!  the pre-tabulated derivative of the Planck function with respect to
!  temperature integrated over each spectral band.
!***************************************************************************

! ------- Declarations -------

      implicit none

! ----- Input -----
      integer(kind=im), intent(in) :: nlayers         ! total number of layers

      integer(kind=im), intent(in) :: nsfc            ! total number of surface layers

! Atmosphere
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=rb), intent(in) :: semiss(nbndlw,nsfc)! lw surface emissivity

      real(kind=rb), intent(in) :: planklayb(nbndlw,nlayers)
      real(kind=rb), intent(in) :: dplankdnb(nbndlw,nlayers)
      real(kind=rb), intent(in) :: dplankupb(nbndlw,nlayers)
      real(kind=rb), intent(in) :: planksfcb(nbndlw,nsfc)
      real(kind=rb), intent(in) :: dplnksfcb(nbndlw,nsfc)

      real(kind=rb), intent(in) :: fracs(ngptlw,nlayers)

      real(kind=rb), intent(inout) :: taut(ngptlw,nlayers) ! gaseous + aerosol optical depths

      real(kind=rb), intent(in) :: frac_sfck(nsfc)

! Clouds
      real(kind=rb), intent(in) :: cldfr        (nlayers) ! layer cloud fraction [mcica]
      real(kind=rb), intent(in) :: taucmc(ngptlw,nlayers) ! layer cloud optical depth [mcica]

! ----- Output -----
      real(kind=rb), intent(out) :: totuflux(nlayers)  ! upward longwave flux (w/m2)
      real(kind=rb), intent(out) :: totdflux(nlayers)  ! downward longwave flux (w/m2)
      real(kind=rb), intent(out) :: totuflux_sfc(nsfc) ! upward longwave flux (w/m2)
      real(kind=rb), intent(out) :: totdflux_sfc(nsfc) ! downward longwave flux (w/m2)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=rb) :: urad(ngptlw,nlayers)
      real(kind=rb) :: drad(ngptlw,nlayers)

      real(kind=rb) :: drad_sfc(ngptlw,1:nsfc)
      real(kind=rb) :: urad_sfc(ngptlw,1:nsfc)

      real(kind=rb) :: reflect(nbndlw,nsfc)
      real(kind=rb) :: transtot(ngptlw,nlayers)
      real(kind=rb) :: secdif, secdiff(ngptlw)  ! secant of diffusivity angle
      real(kind=rb) :: tfactot, atot, odtot, tf0, rad0

      real(kind=rb) :: bbutot    (ngptlw,nlayers)
      real(kind=rb) :: bbdtot    (ngptlw,nlayers)
      real(kind=rb) :: bbdtot_sfc(ngptlw,nsfc   )

      integer(kind=im) :: ig, iband, lev     ! loop indices

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdclr                       ! gas-only planck function for downward rt
!    bbuclr                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance
!    radclru                      ! spectrally summed clear sky upward radiance
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance
!    radclrd                      ! spectrally summed clear sky downward radiance
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer
!    d_radlu_dt                   ! spectrally summed upward radiance
!    d_radclru_dt                 ! spectrally summed clear sky upward radiance
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)
!    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature
!    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature

! This secant and weight corresponds to the standard diffusivity
! angle.  This initial value is redefined below for some bands.

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.

      real(kind=rb), parameter :: a0(nbndlw) =         &
           (/-1.660_rb, -1.550_rb, -1.58_rb, -1.66_rb, &
             -1.540_rb, -1.454_rb, -1.89_rb, -1.33_rb, &
             -1.668_rb, -1.660_rb, -1.66_rb, -1.66_rb, &
             -1.660_rb, -1.660_rb, -1.66_rb, -1.66_rb /)

      real(kind=rb), parameter :: a1(nbndlw) =         &
           (/ 0.000_rb,  0.250_rb,  0.22_rb,  0.00_rb, &
              0.130_rb,  0.446_rb, -0.10_rb,  0.40_rb, &
             -0.006_rb,  0.000_rb,  0.00_rb,  0.00_rb, &
              0.000_rb,  0.000_rb,  0.00_rb,  0.00_rb /)

      real(kind=rb), parameter :: a2(nbndlw) =          &
           (/ 0.000_rb, -12.00_rb, -11.7_rb,  0.000_rb, &
             -0.720_rb, -0.243_rb,  0.19_rb, -0.062_rb, &
              0.414_rb,  0.000_rb,  0.00_rb,  0.000_rb, &
              0.000_rb,  0.000_rb,  0.00_rb,  0.000_rb /)

      do iband = 1, nbndlw
         if (iband == 1 .or. iband == 4 .or. iband > 9) then
            secdif = -1.666
         else
            secdif = a0(iband) - a1(iband)*exp(a2(iband)*pwvcm)
            if (secdif .lt. -1.80_rb) secdif = -1.80_rb
            if (secdif .gt. -1.50_rb) secdif = -1.50_rb
         endif
         do ig = ngs1(iband), ngs(iband)
            secdiff(ig) = secdif
         enddo
      enddo

      do lev = 1, nsfc
         do iband = 1,nbndlw
            reflect(iband,lev) = 1.0 - semiss(iband,lev)
         enddo
      enddo

      do lev = 1, nlayers

         ! Compute optical depths/transmissivities/planck functions

         if (cldfr(lev) > cldmin) then
            do ig = 1, ngptlw
               taut(ig,lev) = taut(ig,lev) + taucmc(ig,lev)
            enddo
         endif

         do ig = 1, ngptlw
            iband = ngb(ig)

            odtot            = secdiff(ig) * taut(ig,lev)
            transtot(ig,lev) = exp(odtot)
            atot             = fracs(ig,lev) - fracs(ig,lev) * transtot(ig,lev)

            if (odtot > -0.1) then
               tfactot = rec_6 * odtot * atot
            else
               tfactot = 2._rb * (fracs(ig,lev) + atot / odtot) - atot
            endif

            tf0            = atot * planklayb(iband,lev)
            bbdtot(ig,lev) = tf0 + tfactot * dplankdnb(iband,lev)
            bbutot(ig,lev) = tf0 + tfactot * dplankupb(iband,lev)

            ! special for shaved cells:
            if (lev > 1 .and. lev <= nsfc) then
               bbdtot_sfc(ig,lev) = tf0 + tfactot * dplnksfcb(iband,lev)
            endif

         enddo

      enddo

      bbdtot_sfc(1:ngptlw,1) = bbdtot(1:ngptlw,1)

! Downward radiative transfer loop

      drad(1:ngptlw,nlayers  ) = 0.0
      drad(1:ngptlw,nlayers-1) = bbdtot(1:ngptlw,nlayers)

      do lev = nlayers-1, 2, -1
         !dir$ ivdep
         do ig = 1, ngptlw
            drad(ig,lev-1) = drad(ig,lev) * transtot(ig,lev) + bbdtot(ig,lev)
         enddo
      enddo

      ! special at surface and with shaved cells
      do lev = nsfc, 1, -1
         do ig = 1, ngptlw
            drad_sfc(ig,lev) = drad(ig,lev) * transtot(ig,lev) + bbdtot_sfc(ig,lev)
         enddo
      enddo

!  Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when
!  they are defined in subroutine setcoef.

      do ig = 1, ngptlw
         iband = ngb(ig)
         rad0           = fracs(ig,1) * planksfcb(iband,1)
         urad_sfc(ig,1) = rad0 + reflect(iband,1) * drad_sfc(ig,1)
         urad    (ig,1) = urad_sfc(ig,1) * transtot(ig,1) + bbutot(ig,1)
      enddo

      ! special for shaved cells:
      do lev = 2, nsfc

         !dir$ ivdep
         do ig = 1, ngptlw
            iband = ngb(ig)

            rad0             = fracs(ig,lev) * planksfcb(iband,lev)

            urad_sfc(ig,lev) = rad0 + reflect(iband,lev) * drad_sfc(ig,lev)

            urad    (ig,lev) = ( urad_sfc(ig,lev  ) *     frac_sfck(lev)    &
                               + urad    (ig,lev-1) * (1.-frac_sfck(lev)) ) &
                              * transtot(ig,lev) + bbutot(ig,lev)
         enddo
      enddo

! Upward radiative transfer loop

      do lev = nsfc+1, nlayers
         !dir$ ivdep
         do ig = 1, ngptlw
            urad(ig,lev) = urad(ig,lev-1) * transtot(ig,lev) + bbutot(ig,lev)
         enddo
      enddo

! Calculate total upward and downward fluxes

      do lev = 1, nlayers
         totuflux(lev) = sum( urad(:,lev) * delwavef(:) )
         totdflux(lev) = sum( drad(:,lev) * delwavef(:) )
      enddo

      do lev = 1, nsfc
         totuflux_sfc(lev) = sum( urad_sfc(:,lev) * delwavef(:) )
         totdflux_sfc(lev) = sum( drad_sfc(:,lev) * delwavef(:) )
      enddo

      end subroutine rtrnmc_noclr

      end module rrtmg_lw_rtrnmc
