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

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngb, ngs, ngc
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrrtc, hnamrtc

      implicit none

      contains

!-----------------------------------------------------------------------------

      subroutine rtrnmc(nlayers, nsfc, pz, semiss, ncbands, &
                        cldfmc, taucmc, planklay, planklev, &
                        planksfc, frac_sfck, &
                        pwvcm, fracs, taut, &
                        totuflux, totdflux, &
                        totuclfl, totdclfl, &
                        totuflux_sfc, totdflux_sfc, &
                        totuclfl_sfc, totdclfl_sfc )

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
      real(kind=rb), intent(in) :: pz(0:nlayers+1)    ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=rb), intent(in) :: semiss(nsfc,nbndlw)! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: planklay(nlayers+1,nbndlw)
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=rb), intent(in) :: planklev(0:nlayers+1,nbndlw)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(in) :: planksfc(nsfc,nbndlw)   ! 
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: fracs(nlayers+1,ngptlw)
                                                      !    Dimensions: (nlayers,ngptw)
      real(kind=rb), intent(in) :: taut(nlayers+1,ngptlw)! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=rb), intent(in) :: frac_sfck(nsfc)

! Clouds
      integer(kind=im), intent(in) :: ncbands         ! number of cloud spectral bands
      real(kind=rb), intent(in) :: cldfmc(ngptlw,nlayers+1) ! layer cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: taucmc(ngptlw,nlayers+1)        ! layer cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)

! ----- Output -----
      real(kind=rb), intent(out) :: totuflux(nlayers+1)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux(nlayers+1)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: fnet(0:nlayers+1) ! net longwave flux (w/m2)
!                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: htr(0:nlayers+1)  ! longwave heating rate (k/day)
!                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl(nlayers+1)! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl(nlayers+1)! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: fnetc(0:nlayers+1)! clear sky net longwave flux (w/m2)
!                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: htrc(0:nlayers+1) ! clear sky longwave heating rate (k/day)
!                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuflux_sfc(nsfc)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux_sfc(nsfc)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl_sfc(nsfc)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl_sfc(nsfc)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
! ----- Local -----
! Declarations for radiative transfer
!     real(kind=rb) :: abscld(nlayers,ngptlw)
!     real(kind=rb) :: atot(nlayers)
!     real(kind=rb) :: atrans(nlayers)
      real(kind=rb) :: bbuclr(nlayers)
      real(kind=rb) :: bbutot(nlayers)
      real(kind=rb) :: clrurad(0:nlayers)
      real(kind=rb) :: clrdrad(0:nlayers)
!     real(kind=rb) :: efclfrac(nlayers,ngptlw)
      real(kind=rb) :: uflux(0:nlayers)
      real(kind=rb) :: dflux(0:nlayers)
      real(kind=rb) :: urad(0:nlayers)
      real(kind=rb) :: drad(0:nlayers)
      real(kind=rb) :: uclfl(0:nlayers)
      real(kind=rb) :: dclfl(0:nlayers)
!     real(kind=rb) :: odcld(nlayers,ngptlw)

      real(kind=rb) :: drad_sfc(0:nsfc)
      real(kind=rb) :: urad_sfc(0:nsfc)

      real(kind=rb) :: clrdrad_sfc(0:nsfc)
      real(kind=rb) :: clrurad_sfc(0:nsfc)

      real(kind=rb) :: transclr(nlayers)
      real(kind=rb) :: transtot(nlayers)
 
      real(kind=rb) :: secdiff(nbndlw)                 ! secant of diffusivity angle
      real(kind=rb) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=rb) :: odepth, odtot, odepth_rec, odtot_rec !, gassrc
      real(kind=rb) :: tblind, tfactot, tfacgas, transc, tausfac, aclr, atot
      real(kind=rb) :: radlu, radclru, dplanksfc
      real(kind=rb) :: reflect(nsfc,nbndlw)
      real(kind=rb) :: plankbnd(nsfc,nbndlw)

      real(kind=rb) :: rad0(0:nsfc)
      real(kind=rb) :: bbdclr(nlayers), bbdtot(nlayers)
      real(kind=rb) :: bbdclr_sfc(nsfc), bbdtot_sfc(nsfc)

      integer(kind=im) :: iband, lev, l      ! loop indices
      integer(kind=im) :: igc                          ! g-point interval counter
      integer(kind=im) :: iclddn                       ! flag for cloud in down path
      integer(kind=im) :: ittot, itgas, itr            ! lookup table indices

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
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

      real(kind=rb), parameter :: wtdiff = 0.5_rb
      real(kind=rb), parameter :: rec_6  = 0.166667_rb
      real(kind=rb), parameter :: fluxfacw = fluxfac * wtdiff

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.

      real(kind=rb), parameter :: a0(nbndlw) =         &
           (/ 1.660_rb,  1.550_rb,  1.58_rb,  1.66_rb, &
              1.540_rb,  1.454_rb,  1.89_rb,  1.33_rb, &
              1.668_rb,  1.660_rb,  1.66_rb,  1.66_rb, &
              1.660_rb,  1.660_rb,  1.66_rb,  1.66_rb /)

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

!     hvrrtc = '$Revision: 1.7 $'

      do iband = 1,nbndlw
         if (iband.eq.1 .or. iband.eq.4 .or. iband.ge.10) then
            secdiff(iband) = 1.66_rb
         else
            secdiff(iband) = a0(iband) + a1(iband)*exp(a2(iband)*pwvcm)
            if (secdiff(iband) .gt. 1.80_rb) secdiff(iband) = 1.80_rb
            if (secdiff(iband) .lt. 1.50_rb) secdiff(iband) = 1.50_rb
         endif
      enddo

      do lev = 1, nlayers+1
         totuflux(lev) = 0.0_rb
         totdflux(lev) = 0.0_rb
         totuclfl(lev) = 0.0_rb
         totdclfl(lev) = 0.0_rb
      enddo

      do lev = 1, nsfc
         totuflux_sfc(lev) = 0.0_rb
         totdflux_sfc(lev) = 0.0_rb
         totuclfl_sfc(lev) = 0.0_rb
         totdclfl_sfc(lev) = 0.0_rb
      enddo

      drad   (nlayers) = 0.0_rb
      clrdrad(nlayers) = 0.0_rb

      do iband = 1,nbndlw
         do lev = 1, nsfc
            reflect (lev,iband)  = 1.0 - semiss(lev,iband)
            plankbnd(lev,iband) = semiss(lev,iband) * planksfc(lev,iband)
         enddo
      enddo


! Loop over g-channels
      do igc = 1, ngptlw
         iband = ngb(igc)

! Compute optical depths/transmissivities/planck functions

         do lev = 1, nlayers
            plfrac = fracs(lev,igc)
            blay = planklay(lev,iband)
            dplankup = planklev(lev,iband) - blay
            dplankdn = planklev(lev-1,iband) - blay

            !  Clear layer properties

            odepth = max(secdiff(iband) * taut(lev,igc), 0.0_rb)

            ! If not using lookup tables:
            !
            !transclr(lev) = max( exp(-odepth), 1.e-20 )
            !aclr          = 1._rb - transclr(lev)
            !if (odepth < 0.3) then
            !   tausfac  = rec_6*odepth
            !else
            !   tausfac  = 1._rb - 2._rb * ( 1._rb / odepth - transclr(lev) / aclr )
            !endif
               
            if (odepth .le. 0.06_rb) then
               transclr(lev) = 1.0 - odepth + 0.5_rb*odepth*odepth
               tausfac       = rec_6*odepth
            else
               tblind        = odepth/(bpade+odepth)
               itr           = tblint*tblind+0.5_rb
               transclr(lev) = exp_tbl(itr)
               tausfac       = tfn_tbl(itr)
            endif
            
            aclr        = 1.0 - transclr(lev)
            bbdclr(lev) = plfrac * (blay + tausfac * dplankdn) * aclr
            bbuclr(lev) = plfrac * (blay + tausfac * dplankup) * aclr

            ! special for shaved cells:
            if (lev <= nsfc) then
               dplanksfc       = planksfc(lev,iband) - blay
               bbdclr_sfc(lev) = plfrac * (blay + tausfac * dplanksfc) * aclr
            endif

            !  Cloudy layer properties

            if (cldfmc(igc,lev) > 0.99_rb) then

               odtot = odepth + secdiff(iband) * taucmc(igc,lev)

               ! If not using lookup tables:
               !
               !transtot(lev) = max( exp(-odtot), 1.e-20 )
               !atot          = 1._rb - transtot(lev)
               !if (odepth < 0.3) then
               !   tfactot  = rec_6*odtot
               !else
               !   tfactot  = 1._rb - 2._rb * ( 1._rb / odtot - transtot(lev) / atot )
               !endif

               if (odtot .lt. 0.06_rb) then
                  transtot(lev) = 1.0 - odtot + 0.5_rb*odtot*odtot
                  tfactot       = rec_6*odtot
               else
                  tblind = odtot/(bpade+odtot)
                  ittot = tblint*tblind + 0.5_rb
                  transtot(lev) = exp_tbl(ittot)
                  tfactot = tfn_tbl(ittot)
               endif

               atot        = 1.0 - transtot(lev)
               bbdtot(lev) = plfrac * (blay + tfactot * dplankdn) * atot
               bbutot(lev) = plfrac * (blay + tfactot * dplankup) * atot

               ! special for shaved cells:
               if (lev <= nsfc) then
                  bbdtot_sfc(lev) = plfrac * (blay + tfactot * dplanksfc) * atot
               endif

             else

                transtot(lev) = transclr(lev)
                bbdtot  (lev) = bbdclr  (lev)
                bbutot  (lev) = bbuclr  (lev)

                ! special for shaved cells:
                if (lev <= nsfc) then
                   bbdtot_sfc(lev) = bbdclr_sfc(lev)
                endif

             endif

         enddo

! Radiative transfer starts here.

         iclddn = 0

! Downward radiative transfer loop

         do lev = nlayers, 1, -1

            if (cldfmc(igc,lev) > 0.99_rb) iclddn = 1

! Total Sky
            drad(lev-1) = drad(lev) * transtot(lev) + bbdtot(lev)

!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.

            if (iclddn == 1) then
               clrdrad(lev-1) = clrdrad(lev) * transclr(lev) + bbdclr(lev)
            else
               clrdrad(lev-1) = drad(lev-1)
            endif

         enddo

         ! special for shaved cells:
         do lev = nsfc, 1, -1
            drad_sfc(lev-1) = drad(lev) * transtot(lev) + bbdtot_sfc(lev)

            if (iclddn == 1) then
               clrdrad_sfc(lev-1) = clrdrad(lev) * transclr(lev) + bbdclr_sfc(lev)
            else
               clrdrad_sfc(lev-1) = drad_sfc(lev-1)
            endif
         enddo

!  Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
!  they are defined in subroutine setcoef. 

         rad0(0) = fracs(1,igc) * plankbnd(1,iband)

!  Add in specular reflection of surface downward radiance.

         urad_sfc(0) = rad0(0) + reflect(1,iband) * drad_sfc(0)
         urad    (0) = urad_sfc(0)

         ! special for shaved cells
         if (nsfc > 1) then
            do lev = 1, nsfc-1

               rad0(lev) = fracs(lev+1,igc) * plankbnd(lev+1,iband)

               urad_sfc(lev) = rad0(lev) + reflect(lev+1,iband) * drad_sfc(lev)

               urad(lev) = ( urad_sfc(lev-1) *     frac_sfck(lev)  &
                           + urad    (lev-1) * (1.-frac_sfck(lev)) ) * transtot(lev) + bbutot(lev)
            enddo
            
            urad(nsfc) = ( urad_sfc(nsfc-1) *     frac_sfck(nsfc)  &
                         + urad    (nsfc-1) * (1.-frac_sfck(nsfc)) ) * transtot(nsfc) + bbutot(nsfc)
         else
            urad(nsfc) = urad(nsfc-1) * transtot(nsfc) + bbutot(nsfc)
         endif

! Upward radiative transfer loop, total sky

         do lev = nsfc+1, nlayers
            urad(lev) = urad(lev-1) * transtot(lev) + bbutot(lev)
         enddo

! Now do clear sky

         if (iclddn == 1) then

            clrurad_sfc(0) = rad0(0) + reflect(1,iband) * clrdrad_sfc(0)
            clrurad    (0) = clrurad_sfc(0)

            ! special for shaved cells
            if (nsfc > 1) then
               do lev = 1, nsfc-1

                  clrurad_sfc(lev) = rad0(lev) + reflect(lev+1,iband) * clrdrad_sfc(lev)

                  clrurad(lev) = ( clrurad_sfc(lev-1) *     frac_sfck(lev)  &
                                 + clrurad    (lev-1) * (1.-frac_sfck(lev)) ) * transclr(lev) + bbuclr(lev)
               enddo

               clrurad(nsfc) = ( clrurad_sfc(nsfc-1) *     frac_sfck(nsfc)  &
                               + clrurad    (nsfc-1) * (1.-frac_sfck(nsfc)) ) * transclr(nsfc) + bbuclr(nsfc)
            else
               clrurad(nsfc) = clrurad(nsfc-1) * transclr(nsfc) + bbuclr(nsfc)
            endif

            ! Upward radiative transfer loop, clear sky
            do lev = nsfc+1, nlayers
               clrurad(lev) = clrurad(lev-1) * transclr(lev) + bbuclr(lev)
            enddo

         else

            clrurad_sfc(0:nsfc-1)  = urad_sfc(0:nsfc-1)
            clrurad    (0:nlayers) = urad    (0:nlayers)

         endif

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.

         do lev = 1, nlayers+1
            totuflux(lev) = totuflux(lev) +    urad(lev-1) * delwave(iband)
            totdflux(lev) = totdflux(lev) +    drad(lev-1) * delwave(iband)
            totuclfl(lev) = totuclfl(lev) + clrurad(lev-1) * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + clrdrad(lev-1) * delwave(iband)
         enddo

         do lev = 1, nsfc
            totuflux_sfc(lev) = totuflux_sfc(lev) +    urad_sfc(lev-1) * delwave(iband)
            totdflux_sfc(lev) = totdflux_sfc(lev) +    drad_sfc(lev-1) * delwave(iband)
            totuclfl_sfc(lev) = totuclfl_sfc(lev) + clrurad_sfc(lev-1) * delwave(iband)
            totdclfl_sfc(lev) = totdclfl_sfc(lev) + clrdrad_sfc(lev-1) * delwave(iband)
         enddo

! End spectral g-point loop
      enddo

! Calculate fluxes at model levels

      do lev = 1, nlayers+1
         totuflux(lev) = totuflux(lev) * fluxfacw
         totdflux(lev) = totdflux(lev) * fluxfacw
         totuclfl(lev) = totuclfl(lev) * fluxfacw
         totdclfl(lev) = totdclfl(lev) * fluxfacw
      enddo

      do lev = 1, nsfc
         totuflux_sfc(lev) = totuflux_sfc(lev) * fluxfacw
         totdflux_sfc(lev) = totdflux_sfc(lev) * fluxfacw

         totuclfl_sfc(lev) = totuclfl_sfc(lev) * fluxfacw
         totdclfl_sfc(lev) = totdclfl_sfc(lev) * fluxfacw
      enddo
      
! Calculate heating rates at model layers

!!      do l = 0, nlayers-1
!!         lev = l + 1
!!         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
!!         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
!!      enddo

! Set heating rate to zero in top layer

!!      htr (nlayers) = 0.0_rb
!!      htrc(nlayers) = 0.0_rb

      end subroutine rtrnmc

!-----------------------------------------------------------------------------

      subroutine rtrnmc_noclr( &
                        nlayers, nsfc, pz, semiss, ncbands, &
                        cldfmc, taucmc, planklay, planklev, &
                        planksfc, frac_sfck, &
                        pwvcm, fracs, taut, &
                        totuflux, totdflux, &
                        totuclfl, totdclfl, &
                        totuflux_sfc, totdflux_sfc, &
                        totuclfl_sfc, totdclfl_sfc )

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
      real(kind=rb), intent(in) :: pz(0:nlayers+1)    ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=rb), intent(in) :: semiss(nsfc,nbndlw)! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: planklay(nlayers+1,nbndlw)
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=rb), intent(in) :: planklev(0:nlayers+1,nbndlw)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(in) :: planksfc(nsfc,nbndlw)
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: fracs(nlayers+1,ngptlw)
                                                      !    Dimensions: (nlayers,ngptw)
      real(kind=rb), intent(in) :: taut(nlayers+1,ngptlw)! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=rb), intent(in) :: frac_sfck(nsfc)

! Clouds
      integer(kind=im), intent(in) :: ncbands         ! number of cloud spectral bands
      real(kind=rb), intent(in) :: cldfmc(ngptlw,nlayers+1) ! layer cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: taucmc(ngptlw,nlayers+1)        ! layer cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)

! ----- Output -----
      real(kind=rb), intent(out) :: totuflux(nlayers+1)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux(nlayers+1)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: fnet(0:nlayers+1) ! net longwave flux (w/m2)
!                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: htr(0:nlayers+1)  ! longwave heating rate (k/day)
!                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl(nlayers+1)! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl(nlayers+1)! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: fnetc(0:nlayers+1)! clear sky net longwave flux (w/m2)
!                                                      !    Dimensions: (0:nlayers)
!      real(kind=rb), intent(out) :: htrc(0:nlayers+1) ! clear sky longwave heating rate (k/day)
!                                                      !    Dimensions: (0:nlayers)

      real(kind=rb), intent(out) :: totuflux_sfc(nsfc)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux_sfc(nsfc)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl_sfc(nsfc)! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl_sfc(nsfc)! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
! ----- Local -----
! Declarations for radiative transfer
!     real(kind=rb) :: abscld(nlayers,ngptlw)
!     real(kind=rb) :: atot(nlayers)
!     real(kind=rb) :: atrans(nlayers)
      real(kind=rb) :: bbuclr(nlayers)
      real(kind=rb) :: bbutot(nlayers)
      real(kind=rb) :: clrurad(0:nlayers)
      real(kind=rb) :: clrdrad(0:nlayers)
!     real(kind=rb) :: efclfrac(nlayers,ngptlw)
      real(kind=rb) :: uflux(0:nlayers)
      real(kind=rb) :: dflux(0:nlayers)
      real(kind=rb) :: urad(0:nlayers)
      real(kind=rb) :: drad(0:nlayers)
      real(kind=rb) :: uclfl(0:nlayers)
      real(kind=rb) :: dclfl(0:nlayers)
!     real(kind=rb) :: odcld(nlayers,ngptlw)

      real(kind=rb) :: drad_sfc(0:nsfc)
      real(kind=rb) :: urad_sfc(0:nsfc)

      real(kind=rb) :: transclr(nlayers)
      real(kind=rb) :: transtot(nlayers)
 
      real(kind=rb) :: secdiff(nbndlw)                 ! secant of diffusivity angle
      real(kind=rb) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=rb) :: odepth, odtot, odepth_rec, odtot_rec !, gassrc
      real(kind=rb) :: tblind, tfactot, tfacgas, transc, tausfac, aclr, atot
      real(kind=rb) :: rad0, radlu, radclru, dplanksfc
      real(kind=rb) :: reflect (nsfc,nbndlw)
      real(kind=rb) :: plankbnd(nsfc,nbndlw)

      real(kind=rb) :: bbdtot(nlayers)
      real(kind=rb) :: bbdtot_sfc(nsfc)

      integer(kind=im) :: iband, lev, l      ! loop indices
      integer(kind=im) :: igc                          ! g-point interval counter
      integer(kind=im) :: iclddn                       ! flag for cloud in down path
      integer(kind=im) :: ittot, itgas, itr            ! lookup table indices

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
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

      real(kind=rb), parameter :: wtdiff = 0.5_rb
      real(kind=rb), parameter :: rec_6  = 0.166667_rb
      real(kind=rb), parameter :: fluxfacw = fluxfac * wtdiff

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.

      real(kind=rb), parameter :: a0(nbndlw) =         &
           (/ 1.660_rb,  1.550_rb,  1.58_rb,  1.66_rb, &
              1.540_rb,  1.454_rb,  1.89_rb,  1.33_rb, &
              1.668_rb,  1.660_rb,  1.66_rb,  1.66_rb, &
              1.660_rb,  1.660_rb,  1.66_rb,  1.66_rb /)

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

!     hvrrtc = '$Revision: 1.7 $'

      do iband = 1,nbndlw
         if (iband.eq.1 .or. iband.eq.4 .or. iband.ge.10) then
            secdiff(iband) = 1.66_rb
         else
            secdiff(iband) = a0(iband) + a1(iband)*exp(a2(iband)*pwvcm)
            if (secdiff(iband) .gt. 1.80_rb) secdiff(iband) = 1.80_rb
            if (secdiff(iband) .lt. 1.50_rb) secdiff(iband) = 1.50_rb
         endif
      enddo

      do lev = 1, nlayers+1
         totuflux(lev) = 0.0_rb
         totdflux(lev) = 0.0_rb
         totuclfl(lev) = 0.0_rb
         totdclfl(lev) = 0.0_rb
      enddo

      do lev = 1, nsfc
         totuflux_sfc(lev) = 0.0_rb
         totdflux_sfc(lev) = 0.0_rb
         totuclfl_sfc(lev) = 0.0_rb
         totdclfl_sfc(lev) = 0.0_rb
      enddo

      drad   (nlayers) = 0.0_rb
      clrdrad(nlayers) = 0.0_rb

      do iband = 1,nbndlw
         do lev = 1, nsfc
            reflect (lev,iband)  = 1.0 - semiss(lev,iband)
            plankbnd(lev,iband) = semiss(lev,iband) * planksfc(lev,iband)
         enddo
      enddo

! Loop over g-channels
      do igc = 1, ngptlw
         iband = ngb(igc)

! Compute optical depths/transmissivities/planck functions

         do lev = 1, nlayers
            plfrac = fracs(lev,igc)
            blay = planklay(lev,iband)
            dplankup = planklev(lev,iband) - blay
            dplankdn = planklev(lev-1,iband) - blay

            !  Total sky properties

            odtot = max( secdiff(iband) * (taut(lev,igc) + taucmc(igc,lev)), 0.0_rb)

            ! If not using lookup tables:
            !
            !transtot(lev) = max( exp(-odtot), 1.e-20 )
            !atot          = 1._rb - transtot(lev)
            !if (odepth < 0.3) then
            !   tfactot  = rec_6*odtot
            !else
            !   tfactot  = 1._rb - 2._rb * ( 1._rb / odtot - transtot(lev) / atot )
            !endif

            if (odtot .lt. 0.06_rb) then
               transtot(lev) = 1.0 - odtot + 0.5_rb*odtot*odtot
               tfactot       = rec_6*odtot
            else
               tblind = odtot/(bpade+odtot)
               ittot = tblint*tblind + 0.5_rb
               transtot(lev) = exp_tbl(ittot)
               tfactot = tfn_tbl(ittot)
            endif

            atot        = 1.0 - transtot(lev)
            bbdtot(lev) = plfrac * (blay + tfactot * dplankdn) * atot
            bbutot(lev) = plfrac * (blay + tfactot * dplankup) * atot

            ! special for shaved cells:
            if (lev <= nsfc) then
               dplanksfc       = planksfc(lev,iband) - blay
               bbdtot_sfc(lev) = plfrac * (blay + tfactot * dplanksfc) * atot
            endif
         enddo

! Radiative transfer starts here.

! Downward radiative transfer loop, total sky only.

         do lev = nlayers, 1, -1
            drad(lev-1) = drad(lev) * transtot(lev) + bbdtot(lev)
         enddo

         ! special for shaved cells:
         do lev = nsfc, 1, -1
            drad_sfc(lev-1) = drad(lev) * transtot(lev) + bbdtot_sfc(lev)
         enddo

!  Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
!  they are defined in subroutine setcoef. 

         rad0 = fracs(1,igc) * plankbnd(1,iband)

!  Add in specular reflection of surface downward radiance.

         urad_sfc(0) = rad0 + reflect(1,iband) * drad_sfc(0)
         urad    (0) = urad_sfc(0)

         ! special for shaved cells
         if (nsfc > 1) then
            do lev = 1, nsfc-1

               rad0 = fracs(lev+1,igc) * plankbnd(lev+1,iband)

               urad_sfc(lev) = rad0 + reflect(lev+1,iband) * drad_sfc(lev)

               urad(lev) = ( urad_sfc(lev-1) *     frac_sfck(lev)  &
                           + urad    (lev-1) * (1.-frac_sfck(lev)) ) * transtot(lev) + bbutot(lev)
            enddo
            
            urad(nsfc) = ( urad_sfc(nsfc-1) *     frac_sfck(nsfc)  &
                         + urad    (nsfc-1) * (1.-frac_sfck(nsfc)) ) * transtot(nsfc) + bbutot(nsfc)
         else
            urad(nsfc) = urad(nsfc-1) * transtot(nsfc) + bbutot(nsfc)
         endif

! Upward radiative transfer loop, total sky

         do lev = nsfc+1, nlayers
            urad(lev) = urad(lev-1) * transtot(lev) + bbutot(lev)
         enddo

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.

         do lev = 1, nlayers+1
            totuflux(lev) = totuflux(lev) + urad(lev-1) * delwave(iband)
            totdflux(lev) = totdflux(lev) + drad(lev-1) * delwave(iband)
         enddo

         do lev = 1, nsfc
            totuflux_sfc(lev) = totuflux_sfc(lev) + urad_sfc(lev-1) * delwave(iband)
            totdflux_sfc(lev) = totdflux_sfc(lev) + drad_sfc(lev-1) * delwave(iband)
         enddo

! End spectral g-point loop
      enddo

! Calculate fluxes at model levels

      do lev = 1, nlayers+1
         totuflux(lev) = totuflux(lev) * fluxfacw
         totdflux(lev) = totdflux(lev) * fluxfacw
!        fnet    (lev-1) = totuflux(lev) - totdflux(lev)
!        fnetc   (lev-1) = 0.0
      enddo

!!      if (nsfc > 2) then
!!
!!      write(io6,*) nsfc
!!      write(io6,'(20g13.5)') totuflux_sfc(1:nsfc) * fluxfacw
!!      write(io6,'(20g13.5)') totdflux_sfc(1:nsfc) * fluxfacw
!!      write(io6,*)
!!      write(io6,'(20g13.5)') totuflux(0:4)
!!      write(io6,'(20g13.5)') totdflux(0:4)
!!      write(io6,*)
!!      write(io6,'(20g13.5)') planklay(1:nsfc,1)
!!      write(io6,'(20g13.5)') planklev(0:nsfc-1,1)
!!      write(io6,'(20g13.5)') planksfc(1:nsfc,1)
!!      write(io6,'(20g13.5)') plankbnd(1:nsfc,1)
!!      write(io6,*)
!!      write(io6,'(20g13.5)') planklay(1:nsfc,16)
!!      write(io6,'(20g13.5)') planklev(0:nsfc-1,16)
!!      write(io6,'(20g13.5)') planksfc(1:nsfc,16)
!!      write(io6,'(20g13.5)') plankbnd(1:nsfc,16)
!!      write(io6,*)
!!      write(io6,'(20g13.5)') reflect(1:nsfc,1)
!!      write(io6,'(20g13.5)') semiss(1:nsfc,1)
!!      write(io6,*)
!!      write(io6,'(20g13.5)') reflect(1:nsfc,16)
!!      write(io6,'(20g13.5)') semiss(1:nsfc,16)
!!
!!      write(io6,*)
!!      write(io6,'(20g13.5)') frac_sfck(1:nsfc)
!!      write(io6,*)
!!      write(io6,*)
!!      
!!
!!      endif
!!!     stop
!!
      do lev = 1, nsfc
         totuflux_sfc(lev) = totuflux_sfc(lev) * fluxfacw
         totdflux_sfc(lev) = totdflux_sfc(lev) * fluxfacw
      enddo

! Calculate heating rates at model layers

!!      do l = 0, nlayers-1
!!         lev = l + 1
!!         htr (l) = heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
!!         htrc(l) = 0.0
!!      enddo

! Set heating rate to zero in top layer

!!      htr (nlayers) = 0.0_rb
!!      htrc(nlayers) = 0.0_rb

      end subroutine rtrnmc_noclr


      end module rrtmg_lw_rtrnmc
