!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rad.f90,v $
!     author:    $Author: miacono $
!     revision:  $Revision: 1.12 $
!     created:   $Date: 2011/04/08 20:25:01 $

       module rrtmg_lw_rad

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
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            *
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  miacono@aer.com                            *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! public interfaces/functions/subroutines

      private
      public :: rrtmg_lw

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_lw &
            (nlay    ,nsfc    ,iaer    ,                            &
             pavel   ,tavel   ,tbound  ,frac_sfck, coldry,          &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr   , &
             covmr   ,cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,semiss  , &
             pwvcm   ,inflag  ,iceflag ,liqflag ,cldfr   ,          &
             taucmc  ,ciwpmc  ,clwpmc  ,reicmc  ,relqmc  ,taua    , &
             uflx    ,dflx    ,uflx_sfc                             )

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once.
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa.
!    b) calls CLDPRMC to set cloud optical depth for McICA based
!       on input cloud properties
!    c) calls SETCOEF to calculate various quantities needed for
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo
!       Independent Column Approximation, to represent sub-grid scale
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al.,
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!     Two random number generators are available for use when imca = 1.
!     This is chosen by setting flag irnd on input to mcica_subcol_gen_lw.
!     1) KISSVEC (irnd = 0)
!     2) Mersenne-Twister (irnd = 1)
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprmc.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);
!       cloud optical properties are calculated by cldprmc or cldprmc based
!       on input settings of iceflglw and liqflglw.  Ice particle size provided
!       must be appropriately defined for the ice parameterization selected.
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input
!     flag iaer; see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available.
!
! The optional calculation of the change in upward flux as a function of surface
! temperature is available (controlled by input flag idrv).  This can be utilized
! to approximate adjustments to the upward flux profile caused only by a change in
! surface temperature between full radiation calls.  This feature uses the pre-
! calculated derivative of the Planck function with respect to surface temperature.
!
!    1) Normal forward calculation for the input profile (idrv=0)
!    2) Normal forward calculation with optional calculation of the change
!       in upward flux as a function of surface temperature for clear sky
!       and total sky flux.  Flux partial derivatives are provided in arrays
!       duflx_dt and duflxc_dt for total and clear sky.  (idrv=1)
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced
! set of g-points for application to GCMs.
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.
!-- Added capability to calculate derivative of upward flux wrt surface temperature.
!     Nov 2009: M. J. Iacono, E. J. Mlawer, AER, Inc.

! --------- Modules ----------

      use parkind,          only: im => kind_im, rb => kind_rb
      use rrtmg_lw_cldprmc, only: cldprmc
      use rrtmg_lw_rtrnmc,  only: rtrnmc_noclr
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol,  only: taumol
      use parrrtm,          only: nbndlw, ngptlw, maxxsec
      use rrlw_wvn,         only: ngb

      implicit none

! ------- Declarations -------

! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol)

      integer(kind=im), intent(in) :: nlay            ! Number of model layers
      integer(kind=im), intent(in) :: nsfc            ! Number of layers that intersect surface
      integer(kind=im), intent(in) :: iaer            ! aerosol option flag

      real(kind=rb), intent(in) :: pavel(nlay)        ! layer pressures (mb)
      real(kind=rb), intent(in) :: tavel(nlay)        ! layer temperatures (K)
      real(kind=rb), intent(in) :: tbound(nsfc)       ! surface temperature (K)
      real(kind=rb), intent(in) :: coldry(nlay)       ! dry air column density (mol/cm2)
      real(kind=rb), intent(in) :: frac_sfck(nsfc)
      real(kind=rb), intent(in) :: h2ovmr(nlay)        ! H2O volume mixing ratio
      real(kind=rb), intent(in) :: o3vmr(nlay)         ! O3 volume mixing ratio
      real(kind=rb), intent(in) :: co2vmr(nlay)        ! CO2 volume mixing ratio
      real(kind=rb), intent(in) :: ch4vmr(nlay)        ! Methane volume mixing ratio
      real(kind=rb), intent(in) :: covmr(nlay)         ! CO volume mixing ratio
      real(kind=rb), intent(in) :: n2ovmr(nlay)        ! Nitrous oxide volume mixing ratio
      real(kind=rb), intent(in) :: o2vmr(nlay)         ! Oxygen volume mixing ratio
      real(kind=rb), intent(in) :: cfc11vmr(nlay)      ! CFC11 volume mixing ratio
      real(kind=rb), intent(in) :: cfc12vmr(nlay)      ! CFC12 volume mixing ratio
      real(kind=rb), intent(in) :: cfc22vmr(nlay)      ! CFC22 volume mixing ratio
      real(kind=rb), intent(in) :: ccl4vmr(nlay)       ! CCL4 volume mixing ratio
      real(kind=rb), intent(in) :: semiss(nbndlw,nsfc) ! Surface emissivity
      real(kind=rb), intent(in) :: cldfr(nlay)         ! cloud fraction [mcica]
      real(kind=rb), intent(in) :: pwvcm               ! precipitable water vapor (cm)

      integer(kind=im), intent(in)  :: inflag          ! flag for cloud property method
      integer(kind=im), intent(in)  :: iceflag         ! flag for ice cloud properties
      integer(kind=im), intent(in)  :: liqflag         ! flag for liquid cloud properties

      real(kind=rb), intent(in)   :: ciwpmc(ngptlw,nlay)  ! in-cloud ice water path [mcica]
      real(kind=rb), intent(in)   :: clwpmc(ngptlw,nlay)  ! in-cloud liquid water path [mcica]
      real(kind=rb), intent(in)   :: relqmc(nlay)         ! liquid particle effective radius (microns)
      real(kind=rb), intent(in)   :: reicmc(nlay)         ! ice particle effective size (microns)
      real(kind=rb), intent(inout):: taucmc(ngptlw,nlay)  ! in-cloud optical depth [mcica]
      real(kind=rb), intent(in)   :: taua(nbndlw,nlay)    ! aerosol optical depth

! ----- Output -----

      real(kind=rb), intent(out) :: uflx(nlay+1)      ! Total sky longwave upward flux (W/m2)
      real(kind=rb), intent(out) :: dflx(nlay+1)      ! Total sky longwave downward flux (W/m2)
!     real(kind=rb), intent(out) :: uflxc(nlay)     ! Clear sky longwave upward flux (W/m2)
!     real(kind=rb), intent(out) :: dflxc(nlay)     ! Clear sky longwave downward flux (W/m2)
      real(kind=rb), intent(out) :: uflx_sfc(nsfc)  ! Total sky longwave upward flux (W/m2)
!     real(kind=rb), intent(out) :: dflx_sfc(nsfc)  ! Total sky longwave downward flux (W/m2)
!     real(kind=rb), intent(out) :: uflxc_sfc(nsfc) ! Clear sky longwave upward flux (W/m2)
!     real(kind=rb), intent(out) :: dflxc_sfc(nsfc) ! Clear sky longwave downward flux (W/m2)

! ----- Local -----

! Control
      integer(kind=im) :: k                  ! layer loop index
      integer(kind=im) :: ig                 ! g-point loop index

! Atmosphere
      real(kind=rb) :: wx(nlay,maxxsec)      ! cross-section amounts (mol/cm-2)
      real(kind=rb) :: fracs(ngptlw,nlay)    !
      real(kind=rb) :: taug (ngptlw,nlay)    ! gaseous optical depths

! Atmosphere - setcoef
      integer(kind=im) :: laytrop            ! tropopause layer index
      integer(kind=im) :: jp(nlay)           ! lookup table index
      integer(kind=im) :: jt(nlay)           ! lookup table inde
      integer(kind=im) :: jt1(nlay)          ! lookup table ind

      real(kind=rb) :: planklay(nbndlw,nlay) !
      real(kind=rb) :: plankbnd(nbndlw,nsfc) !

      real(kind=rb) :: colh2o(nlay)          ! column amount (h2o)
      real(kind=rb) :: colco2(nlay)          ! column amount (co2)
      real(kind=rb) :: colo3 (nlay)          ! column amount (o3)
      real(kind=rb) :: coln2o(nlay)          ! column amount (n2o)
      real(kind=rb) :: colco (nlay)          ! column amount (co)
      real(kind=rb) :: colch4(nlay)          ! column amount (ch4)
      real(kind=rb) :: colo2 (nlay)          ! column amount (o2)
      real(kind=rb) :: colbrd(nlay)          ! column amount (broadening gases)

      integer(kind=im) :: indself(nlay)
      integer(kind=im) :: indfor(nlay)
      integer(kind=im) :: indminor(nlay)

      real(kind=rb) :: selffac(nlay)
      real(kind=rb) :: selffrac(nlay)
      real(kind=rb) :: forfac(nlay)
      real(kind=rb) :: forfrac(nlay)
      real(kind=rb) :: minorfrac(nlay)
      real(kind=rb) :: scaleminor(nlay)
      real(kind=rb) :: scaleminorn2(nlay)

      real(kind=rb) :: fac00(nlay), fac01(nlay), &
                       fac10(nlay), fac11(nlay)

      real(kind=rb) :: rat_h2oco2(nlay), rat_h2oco2_1(nlay), &
                       rat_h2oo3 (nlay), rat_h2oo3_1 (nlay), &
                       rat_h2on2o(nlay), rat_h2on2o_1(nlay), &
                       rat_h2och4(nlay), rat_h2och4_1(nlay), &
                       rat_n2oco2(nlay), rat_n2oco2_1(nlay), &
                       rat_o3co2 (nlay), rat_o3co2_1 (nlay)

!  For cloudy atmosphere, use cldprmc to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprmc.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprmc.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprmc.

      if (inflag == 2) then

         call cldprmc(nlay, inflag, iceflag, liqflag, cldfr, &
                      ciwpmc, clwpmc, reicmc, relqmc, taucmc)
      endif

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres.

      call setcoef(nlay, nsfc, pavel, tavel, tbound, semiss, &
                   coldry, wx, h2ovmr, o3vmr, co2vmr, ch4vmr, covmr,&
                   n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
                   laytrop, jp, jt, jt1, planklay, plankbnd, &
                   colh2o, colco2, colo3, coln2o, &
                   colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, &
                   rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                   rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                   rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                   selffac, selffrac, indself, forfac, forfrac, indfor, &
                   minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for
!  each longwave spectral band.

      call taumol(nlay, pavel, wx, coldry, &
                  laytrop, jp, jt, jt1, &
                  colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                  colbrd, fac00, fac01, fac10, fac11, &
                  rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                  rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                  rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                  selffac, selffrac, indself, forfac, forfrac, indfor, &
                  minorfrac, scaleminor, scaleminorn2, indminor, &
                  fracs, taug, n2ovmr, co2vmr)

! Combine gaseous and aerosol optical depths, if aerosol active

      if (iaer .eq. 10) then
         do k = 1, nlay
            do ig = 1, ngptlw
               taug(ig,k) = taug(ig,k) + taua(ngb(ig),k)
            enddo
         enddo
      endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.

      call rtrnmc_noclr(nlay, nsfc, semiss, &
                        cldfr, taucmc, planklay, &
                        plankbnd, frac_sfck, &
                        pwvcm, fracs, taug, &
                        uflx, dflx, uflx_sfc)

      end subroutine rrtmg_lw

      end module rrtmg_lw_rad
