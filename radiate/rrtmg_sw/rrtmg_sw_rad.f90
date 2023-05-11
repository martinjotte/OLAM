!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $

       module rrtmg_sw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------

      use parkind, only: im => kind_im, rb => kind_rb

! public interfaces/functions/subroutines

      private
      public :: rrtmg_sw, earth_sun

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_sw(                                  &
             nlay    ,iaer    ,nsfc    ,frac_sfck,          &
             pavel   ,tavel   ,coldry  ,                    &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr ,    &
             asdir   ,asdif   ,aldir   ,aldif   ,           &
             coszen  ,adjes   ,dyofyr  ,scon    ,           &
             inflag  ,iceflag ,liqflag ,cldfrac ,           &
             taucmc  ,ssacmc  ,asmcmc  ,                    &
             ciwpmc  ,clwpmc  ,reicmc  ,relqmc ,            &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   ,           &
             swuflx  ,swdflx  ,                             &
             swuflx_sfc, swdflx_sfc, swdflx_dir_sfc,        &
             zbbfu_sfc, zbbfd_sfc, zbbfddir_sfc             )

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization
!  area, since this has to be called only once.
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa.
!    c) calls CLDPRMC_SW to set cloud optical depth for McICA based
!       on input cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for
!       the radiative transfer algorithm
!    e) calls SPCVMC to call the two-stream model that in turn
!       calls TAUMOL to calculate gaseous optical depths for each
!       of the 16 spectral bands and to perform the radiative transfer
!       using McICA, the Monte-Carlo Independent Column Approximation,
!       to represent sub-grid scale cloud variability
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al.,
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky or partial cloud conditions.
!
! This call to RRTMG_SW must be preceeded by a call to the module
!     mcica_subcol_gen_sw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngptsw) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprmc.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);
!       cloud optical properties are calculated by cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs.
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag.
!     Dec 2008: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parkind,          only: im => kind_im, rb => kind_rb
      use parrrsw,          only: nbndsw, ngptsw, naerec, rrsw_scon
      use rrsw_aer,         only: rsrtaua, rsrpiza, rsrasya
      use rrsw_con,         only: zepzen
      use rrtmg_sw_taumol,  only: taumol_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvmc,  only: spcvmc_sw_noclr

! ------- Declarations

      implicit none

! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol)

      integer(kind=im), intent(in) :: nlay            ! Number of model layers

      integer(kind=im), intent(in) :: iaer            ! Aerosol option flag
                                                      !    0: No aerosol
                                                      !    6: ECMWF method
                                                      !    10:Input aerosol optical
                                                      !       properties
      integer(kind=im), intent(in) :: nsfc

      real(kind=rb), intent(in) :: frac_sfck(nsfc)

      real(kind=rb), intent(in) :: pavel(nlay)          ! layer pressures (mb)

      real(kind=rb), intent(in) :: tavel(nlay)          ! layer temperatures (K)

      real(kind=rb), intent(in) :: coldry(nlay)       ! dry air column amount

      real(kind=rb), intent(in) :: h2ovmr(nlay)        ! H2O volume mixing ratio

      real(kind=rb), intent(in) :: o3vmr(nlay)         ! O3 volume mixing ratio

      real(kind=rb), intent(in) :: co2vmr(nlay)        ! CO2 volume mixing ratio

      real(kind=rb), intent(in) :: ch4vmr(nlay)        ! Methane volume mixing ratio

      real(kind=rb), intent(in) :: o2vmr(nlay)         ! Oxygen volume mixing ratio

      real(kind=rb), intent(in) :: asdir(nsfc)        ! UV/vis surface albedo direct rad

      real(kind=rb), intent(in) :: aldir(nsfc)        ! Near-IR surface albedo direct rad

      real(kind=rb), intent(in) :: asdif(nsfc)        ! UV/vis surface albedo: diffuse rad

      real(kind=rb), intent(in) :: aldif(nsfc)        ! Near-IR surface albedo: diffuse rad

      integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                      !  distance if adjflx not provided)
      real(kind=rb), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance

      real(kind=rb), intent(in) :: coszen             ! Cosine of solar zenith angle

      real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)

      real(kind=rb), intent(inout) :: tauaer(nbndsw,nlay)  ! Aerosol optical depth (iaer=10 only)
                                                           ! (non-delta scaled)
      real(kind=rb), intent(inout) :: ssaaer(nbndsw,nlay)  ! Aerosol single scattering albedo (iaer=10 only)
                                                           ! (non-delta scaled)
      real(kind=rb), intent(inout) :: asmaer(nbndsw,nlay)  ! Aerosol asymmetry parameter (iaer=10 only)
                                                           ! (non-delta scaled)
      real(kind=rb), intent(in) :: ecaer(naerec,nlay)      ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                           ! (non-delta scaled)
! ----- Output -----

      real(kind=rb), intent(out) :: swuflx(nlay)    ! Total sky shortwave upward flux (W/m2)

      real(kind=rb), intent(out) :: swdflx(nlay)    ! Total sky shortwave downward flux (W/m2)

      real(kind=rb), intent(out) :: swuflx_sfc(nsfc)  ! Total sky shortwave upward flux (W/m2)

      real(kind=rb), intent(out) :: swdflx_sfc(nsfc)  ! Total sky shortwave downward flux (W/m2)

      real(kind=rb), intent(out) :: swdflx_dir_sfc(nsfc)  ! Total sky shortwave downward flux (W/m2)

      real(kind=rb), intent(out) :: zbbfu_sfc(nsfc,nbndsw)    ! Upward shortwave flux by band (w/m2)
      real(kind=rb), intent(out) :: zbbfd_sfc(nsfc,nbndsw)    ! Downward shortwave flux by band (w/m2)
      real(kind=rb), intent(out) :: zbbfddir_sfc(nsfc,nbndsw) ! Downward direct shortwave flux by band (w/m2)

! ----- Local -----

      integer(kind=im) :: i                   ! layer loop index                       ! jk
      integer(kind=im) :: ib                  ! band loop index                        ! jsw
      integer(kind=im) :: ia                  ! indices

! Atmosphere
      real(kind=rb) :: cossza                 ! Cosine of solar zenith angle
      real(kind=rb) :: adjflx
      real(kind=rb) :: adjflux                ! adjustment for current Earth/Sun distance
      real(kind=rb) :: solvar                 ! solar constant scaling factor from rrtmg_sw
                                              !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=rb) :: albdir(nsfc,nbndsw)    ! surface albedo, direct          ! zalbp
      real(kind=rb) :: albdif(nsfc,nbndsw)    ! surface albedo, diffuse         ! zalbd

! Atmosphere - setcoef
      integer(kind=im) :: laytrop             ! tropopause layer index
      integer(kind=im) :: jp(nlay)          !
      integer(kind=im) :: jt(nlay)          !
      integer(kind=im) :: jt1(nlay)         !

      real(kind=rb) :: colh2o(nlay)         ! column amount (h2o)
      real(kind=rb) :: colco2(nlay)         ! column amount (co2)
      real(kind=rb) :: colo3(nlay)          ! column amount (o3)
      real(kind=rb) :: colch4(nlay)         ! column amount (ch4)
      real(kind=rb) :: colo2(nlay)          ! column amount (o2)
      real(kind=rb) :: colmol(nlay)         ! column amount

      integer(kind=im) :: indself(nlay)
      integer(kind=im) :: indfor(nlay)

      real(kind=rb) :: selffac(nlay)
      real(kind=rb) :: selffrac(nlay)
      real(kind=rb) :: forfac(nlay)
      real(kind=rb) :: forfrac(nlay)

      real(kind=rb) :: fac00(nlay), fac01(nlay), &
                       fac10(nlay), fac11(nlay)

! Atmosphere/clouds - cldprop
      integer(kind=im), intent(in) :: inflag              ! flag for cloud property method
      integer(kind=im), intent(in) :: iceflag             ! flag for ice cloud properties
      integer(kind=im), intent(in) :: liqflag             ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb), intent(in) :: ciwpmc(ngptsw,nlay)    ! in-cloud ice water path [mcica]
      real(kind=rb), intent(in) :: clwpmc(ngptsw,nlay)    ! in-cloud liquid water path [mcica]
      real(kind=rb), intent(in) :: relqmc(nlay)           ! liquid particle effective radius (microns)
      real(kind=rb), intent(in) :: reicmc(nlay)           ! ice particle effective size (microns)

! Atmosphere/clouds/aerosol - spcvrt,spcvmc

      real(kind=rb), intent(in)    :: cldfrac(nlay)         ! layer cloud fraction
      real(kind=rb), intent(inout) :: taucmc(ngptsw,nlay)   ! cloud optical depth [mcica]
      real(kind=rb), intent(inout) :: asmcmc(ngptsw,nlay)   ! cloud asymmetry parameter [mcica]
      real(kind=rb), intent(inout) :: ssacmc(ngptsw,nlay)   ! cloud single scattering albedo [mcica]

! Arrays from rrtmg_sw_taumol routines

      real(kind=rb) :: ztaug(ngptsw,nlay)
      real(kind=rb) :: zsflxzen(ngptsw)

      real(kind=rb) :: tau, ssa, asm

      if (dyofyr == 0) then
         adjflx = adjes
      else
         adjflx = earth_sun(dyofyr)
      endif

      solvar = scon / rrsw_scon
      adjflux = adjflx * solvar

!  For cloudy atmosphere, use cldprmc to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprmc.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprmc.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprmc.

      if (inflag /= 0) then
         call cldprmc_sw(nlay, inflag, iceflag, liqflag, cldfrac, &
                         ciwpmc, clwpmc, reicmc, relqmc, &
                         taucmc, ssacmc, asmcmc)
      endif

! Calculate coefficients for the temperature and pressure dependence of the
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.


      call setcoef_sw(nlay, pavel, tavel, coldry, &
                      h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, &
                      laytrop, jp, jt, jt1, &
                      colch4, colco2, colh2o, colmol, &
                      colo2, colo3, fac00, fac01, fac10, fac11, &
                      selffac, selffrac, indself, forfac, forfrac, indfor)

! Cosine of the solar zenith angle
! Prevent using value of zero; ideally, SW model is not called from host model when sun
! is below horizon

      cossza = max(coszen, zepzen)

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
      do i = 1, nsfc
         do ib=1,9
            albdir(i,ib) = aldir(i)
            albdif(i,ib) = aldif(i)
         enddo
         albdir(i,nbndsw) = aldir(i)
         albdif(i,nbndsw) = aldif(i)
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(i,ib) = asdir(i)
            albdif(i,ib) = asdif(i)
         enddo
      enddo

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer),
! or set manually here for each aerosol and layer.
      if (iaer == 6) then

         tauaer = 0._rb
         asmaer = 0._rb
         ssaaer = 0._rb

         do i = 1, nlay
            do ib = 1, nbndsw
               do ia = 1, naerec
                  tau = rsrtaua(ia,ib) * ecaer(ia,i)
                  ssa = rsrpiza(ia,ib) * tau
                  asm = rsrasya(ia,ib) * ssa

                  tauaer(ib,i) = tauaer(ib,i) + tau
                  ssaaer(ib,i) = ssaaer(ib,i) + ssa
                  asmaer(ib,i) = asmaer(ib,i) + asm
               enddo
            enddo
         enddo

      endif

! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      call taumol_sw(nlay, &
                     colh2o, colco2, colch4, colo2, colo3, &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     zsflxzen, ztaug)

! Call the 2-stream radiation transfer model

      call spcvmc_sw_noclr &
                 (nlay, iaer, &
                  nsfc, frac_sfck, colmol, &
                  albdif, albdir, &
                  ztaug, zsflxzen, &
                  taucmc, asmcmc, ssacmc, cldfrac, &
                  tauaer, asmaer, ssaaer, cossza, adjflux, &
                  swuflx, swdflx, &
                  zbbfd_sfc, zbbfu_sfc, zbbfddir_sfc)

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

      do i = 1, nsfc
         swuflx_sfc    (i) = sum(zbbfu_sfc   (i,1:nbndsw))
         swdflx_sfc    (i) = sum(zbbfd_sfc   (i,1:nbndsw))
         swdflx_dir_sfc(i) = sum(zbbfddir_sfc(i,1:nbndsw))
      enddo

      end subroutine rrtmg_sw


!*************************************************************************
      real(kind=rb) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

      use rrsw_con, only: pi

      implicit none

      integer(kind=im), intent(in) :: idn

      real(kind=rb) :: gamma
      integer       :: idnn

      idnn = max(1,min(366,idn))-1

      gamma = 2._rb*pi*real(idnn,rb)/365._rb

      ! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_rb + .034221_rb * cos(gamma) + .001289_rb * sin(gamma) + &
                   .000719_rb * cos(2._rb*gamma) + .000077_rb * sin(2._rb*gamma)

      end function earth_sun


      end module rrtmg_sw_rad
