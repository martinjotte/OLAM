!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $

      module rrtmg_sw_spcvmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband, jpb1, jpb2
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspc, hnamspc, hvrrft, hnamrft
      use rrsw_wvn, only : ngc, ngs
      use rrtmg_sw_taumol, only: taumol_sw

      implicit none

      contains

! ---------------------------------------------------------------------------
      subroutine spcvmc_sw &
            (nlayers, istart, iend, icpr, idelm, iout, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pcldfmc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
             ptaua, pasya, pomga, prmu0, coldry, wkl, adjflux, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             zbbfd, zbbfu, zbbcd, zbbcu, &
             zbbfddir, zbbcddir)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker and McICA, the Monte-Carlo
!          Independent Column Approximation, for the representation of 
!          sub-grid cloud variability (i.e. cloud overlap).
!
! Interface:  *spcvmc_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Modified to use McICA: MJIacono, AER, Nov 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers
      integer(kind=im), intent(in) :: istart
      integer(kind=im), intent(in) :: iend
      integer(kind=im), intent(in) :: icpr
      integer(kind=im), intent(in) :: idelm   ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]
      integer(kind=im), intent(in) :: iout
      integer(kind=im), intent(in) :: laytrop
      integer(kind=im), intent(in) :: layswtch
      integer(kind=im), intent(in) :: laylow

      integer(kind=im), intent(in) :: indfor(:)
                                                               !   Dimensions: (nlayers)
      integer(kind=im), intent(in) :: indself(:)
                                                               !   Dimensions: (nlayers)
      integer(kind=im), intent(in) :: jp(:)
                                                               !   Dimensions: (nlayers)
      integer(kind=im), intent(in) :: jt(:)
                                                               !   Dimensions: (nlayers)
      integer(kind=im), intent(in) :: jt1(:)
                                                               !   Dimensions: (nlayers)

      real(kind=rb), intent(in) :: pavel(:)                    ! layer pressure (hPa, mb) 
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: tavel(:)                    ! layer temperature (K)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: pz(0:)                      ! level (interface) pressure (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: tz(0:)                      ! level temperatures (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: tbound                      ! surface temperature (K)
      real(kind=rb), intent(in) :: wkl(:,:)                    ! molecular amounts (mol/cm2) 
                                                               !   Dimensions: (mxmol,nlayers)
      real(kind=rb), intent(in) :: coldry(:)                   ! dry air column density (mol/cm2)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: colmol(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: adjflux(:)                  ! Earth/Sun distance adjustment
                                                               !   Dimensions: (jpband)

      real(kind=rb), intent(in) :: palbd(:)                    ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      real(kind=rb), intent(in) :: palbp(:)                    ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)
      real(kind=rb), intent(in) :: prmu0                       ! cosine of solar zenith angle
      real(kind=rb), intent(in) :: pcldfmc(:,:)                ! cloud fraction [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real(kind=rb), intent(in) :: ptaucmc(:,:)                ! cloud optical depth [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real(kind=rb), intent(in) :: pasycmc(:,:)                ! cloud asymmetry parameter [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real(kind=rb), intent(in) :: pomgcmc(:,:)                ! cloud single scattering albedo [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real(kind=rb), intent(in) :: ptaormc(:,:)                ! cloud optical depth, non-delta scaled [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real(kind=rb), intent(in) :: ptaua(:,:)                  ! aerosol optical depth
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=rb), intent(in) :: pasya(:,:)                  ! aerosol asymmetry parameter
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=rb), intent(in) :: pomga(:,:)                  ! aerosol single scattering albedo
                                                               !   Dimensions: (nlayers,nbndsw)

      real(kind=rb), intent(in) :: colh2o(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: colco2(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: colch4(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: co2mult(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: colo3(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: colo2(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: coln2o(:)
                                                               !   Dimensions: (nlayers)

      real(kind=rb), intent(in) :: forfac(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: forfrac(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: selffac(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: selffrac(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: fac00(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: fac01(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: fac10(:)
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: fac11(:)
                                                               !   Dimensions: (nlayers)

! ------- Output -------
                                                               !   All Dimensions: (nlayers+1)
!!      real(kind=rb), intent(out) :: pbbcd(:)
!!      real(kind=rb), intent(out) :: pbbcu(:)
!!      real(kind=rb), intent(out) :: pbbfd(:)
!!      real(kind=rb), intent(out) :: pbbfu(:)
!!      real(kind=rb), intent(out) :: pbbfddir(:)
!!      real(kind=rb), intent(out) :: pbbcddir(:)
!!
!!      real(kind=rb), intent(out) :: puvcd(:)
!!      real(kind=rb), intent(out) :: puvfd(:)
!!      real(kind=rb), intent(out) :: puvcddir(:)
!!      real(kind=rb), intent(out) :: puvfddir(:)
!!
!!      real(kind=rb), intent(out) :: pnicd(:)
!!      real(kind=rb), intent(out) :: pnifd(:)
!!      real(kind=rb), intent(out) :: pnicddir(:)
!!      real(kind=rb), intent(out) :: pnifddir(:)

! Output - inactive                                            !   All Dimensions: (nlayers+1)
!      real(kind=rb), intent(out) :: puvcu(:)
!      real(kind=rb), intent(out) :: puvfu(:)
!      real(kind=rb), intent(out) :: pnicu(:)
!      real(kind=rb), intent(out) :: pnifu(:)
!      real(kind=rb), intent(out) :: pvscd(:)
!      real(kind=rb), intent(out) :: pvscu(:)
!      real(kind=rb), intent(out) :: pvsfd(:)
!      real(kind=rb), intent(out) :: pvsfu(:)

! ------- Local -------

      logical :: lrtchkclr(nlayers),lrtchkcld(nlayers)

      integer(kind=im)  :: klev
      integer(kind=im) :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer(kind=im) :: iw, jb, jg, jl, jk
!      integer(kind=im), parameter :: nuv = ?? 
!      integer(kind=im), parameter :: nvs = ?? 
      integer(kind=im) :: itind

      real(kind=rb) :: tblind, ze1
      real(kind=rb) :: zclear, zcloud
      real(kind=rb) :: zdbt(nlayers+1), zdbt_nodel(nlayers+1)
      real(kind=rb) :: zgc(nlayers), zgcc(nlayers), zgco(nlayers)
      real(kind=rb) :: zomc(nlayers), zomcc(nlayers), zomco(nlayers)
      real(kind=rb) :: zrdnd(nlayers+1), zrdndc(nlayers+1)
      real(kind=rb) :: zref(nlayers+1), zrefc(nlayers+1), zrefo(nlayers+1)
      real(kind=rb) :: zrefd(nlayers+1), zrefdc(nlayers+1), zrefdo(nlayers+1)
      real(kind=rb) :: zrup(nlayers+1), zrupd(nlayers+1)
      real(kind=rb) :: zrupc(nlayers+1), zrupdc(nlayers+1)
      real(kind=rb) :: zs1(nlayers+1)
      real(kind=rb) :: ztauc(nlayers), ztauo(nlayers)
      real(kind=rb) :: ztdn(nlayers+1), ztdnd(nlayers+1), ztdbt(nlayers+1)
      real(kind=rb) :: ztoc(nlayers), ztor(nlayers)
      real(kind=rb) :: ztra(nlayers+1), ztrac(nlayers+1), ztrao(nlayers+1)
      real(kind=rb) :: ztrad(nlayers+1), ztradc(nlayers+1), ztrado(nlayers+1)
      real(kind=rb) :: zdbtc(nlayers+1), ztdbtc(nlayers+1)
      real(kind=rb) :: zincflx(ngptsw), zdbtc_nodel(nlayers+1) 
      real(kind=rb) :: ztdbt_nodel(nlayers+1), ztdbtc_nodel(nlayers+1)

      real(kind=rb) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=rb) :: zwf, tauorig, ap

      real(kind=rb) :: ptXpo(nlayers+1), ptXpoXpa(nlayers+1)
!     real(kind=rb) :: zincflux                                   ! inactive

      real(kind=rb)            :: prmu0m1
      real(kind=rb), parameter :: repclc = 1.e-12_rb

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=rb) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=rb) :: zsflxzen(16)
      real(kind=rb) :: ztaug(nlayers,ngptsw), ztaur(nlayers,ngptsw)
      real(kind=rb) :: zsflxzen(ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=rb) :: zcd(nlayers+1,ngptsw), zcu(nlayers+1,ngptsw)
      real(kind=rb) :: zfd(nlayers+1,ngptsw), zfu(nlayers+1,ngptsw)

! Inactive arrays
     real(kind=rb), intent(out) :: zbbcd   (:,:), zbbcu   (:,:)
     real(kind=rb), intent(out) :: zbbfd   (:,:), zbbfu   (:,:)
     real(kind=rb), intent(out) :: zbbfddir(:,:), zbbcddir(:,:)

! ------------------------------------------------------------------

! Initializations

!     ib1 = istart
!     ib2 = iend
      klev = nlayers
      iw = 0
!      repclc = 1.e-12_rb
!      zincflux = 0.0_rb

!!      do jk=1,klev+1
!!         pbbcd(jk)=0._rb
!!         pbbcu(jk)=0._rb
!!         pbbfd(jk)=0._rb
!!         pbbfu(jk)=0._rb
!!         pbbcddir(jk)=0._rb
!!         pbbfddir(jk)=0._rb
!!         puvcd(jk)=0._rb
!!         puvfd(jk)=0._rb
!!         puvcddir(jk)=0._rb
!!         puvfddir(jk)=0._rb
!!         pnicd(jk)=0._rb
!!         pnifd(jk)=0._rb
!!         pnicddir(jk)=0._rb
!!         pnifddir(jk)=0._rb
!!      enddo

      do ibm = 1, nbndsw
         do jk = 1, klev + 1
            zbbcd   (jk,ibm) =0._rb
            zbbcu   (jk,ibm) =0._rb
            zbbfd   (jk,ibm) =0._rb
            zbbfu   (jk,ibm) =0._rb
            zbbcddir(jk,ibm) =0._rb
            zbbfddir(jk,ibm) =0._rb
         enddo
      enddo
      
      prmu0m1 = 1.0_rb / prmu0

      lrtchkclr(1:klev)=.true.

! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      call taumol_sw(klev, &
                     colh2o, colco2, colch4, colo2, colo3, colmol, &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     zsflxzen, ztaug, ztaur)

! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      do jb = jpb1, jpb2
         ibm = jb-15
!        igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
!        if (iout.gt.0.and.ibm.ge.2) iw = ngs(ibm-1)

!        do jk=1,klev+1
!           zbbcd(jk)=0.0_rb
!           zbbcu(jk)=0.0_rb
!           zbbfd(jk)=0.0_rb
!           zbbfu(jk)=0.0_rb
!        enddo

         ap = adjflux(jb) * prmu0

         do jk=1,klev
            ikl=klev+1-jk
            ptXpo   (ikl) = ptaua(ikl,ibm) * pomga(ikl,ibm)
            ptXpoXpa(ikl) = pasya(ikl,ibm) * ptXpo(ikl)
         enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1, ngc(ibm)
            iw = iw+1

! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
            zincflx(iw) = ap * zsflxzen(iw)
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ztdbtc(1)=1.0_rb
            ztdbtc_nodel(1)=1.0_rb
!   Surface values
            zdbtc(klev+1) =0.0_rb
            ztrac(klev+1) =0.0_rb
            ztradc(klev+1)=0.0_rb
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)
           
! Cloudy-sky    
!   Surface values
            ztrao(klev+1) =0.0_rb
            ztrado(klev+1)=0.0_rb
            zrefo(klev+1) =palbp(ibm)
            zrefdo(klev+1)=palbd(ibm)
           
! Total sky    
!   TOA direct beam    
            ztdbt(1)=1.0_rb
            ztdbt_nodel(1)=1.0_rb
!   Surface values
            zdbt(klev+1) =0.0_rb
            ztra(klev+1) =0.0_rb
            ztrad(klev+1)=0.0_rb
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)
    
! Top of layer loop
            do jk=1,klev

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here

               ikl=klev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
!              lrtchkclr(jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
!              lrtchkcld(jk)=.false.
               lrtchkcld(jk)=(pcldfmc(ikl,iw) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_rb
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaucmc(ikl,iw)
!               zomco(jk) = ptaucmc(ikl,iw) * pomgcmc(ikl,iw) + ztaur(ikl,iw)
!               zgco (jk) = (ptaucmc(ikl,iw) * pomgcmc(ikl,iw) * pasycmc(ikl,iw) + &
!                           ztaur(ikl,iw) * 0.0001_rb) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm)
               zomcc(jk) = ztaur(ikl,iw) + ptXpo(ikl)
               zgcc(jk)  = ptXpoXpa(ikl) / zomcc(jk)
               zomcc(jk) = zomcc(jk) / ztauc(jk)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for unscaled direct beam calculation
!!               if (idelm .eq. 0) then
!!!     
!!                  zclear = 1.0_rb - pcldfmc(ikl,iw)
!!                  zcloud = pcldfmc(ikl,iw)
!!
!!! Clear
!!!                   zdbtmc = exp(-ztauc(jk) / prmu0)
!!
!!! Use exponential lookup table for transmittance, or expansion of exponential for low tau
!!                  ze1 = ztauc(jk) / prmu0
!!                  if (ze1 .le. od_lo) then
!!                     zdbtmc = 1._rb - ze1 + 0.5_rb * ze1 * ze1
!!                  else 
!!                     tblind = ze1 / (bpade + ze1)
!!                     itind = tblint * tblind + 0.5_rb
!!                     zdbtmc = exp_tbl(itind)
!!                  endif
!!
!!                  zdbtc_nodel(jk) = zdbtmc
!!                  ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)
!!
!!! Clear + Cloud
!!                  tauorig = ztauc(jk) + ptaormc(ikl,iw)
!!!                   zdbtmo = exp(-tauorig / prmu0)
!!
!!! Use exponential lookup table for transmittance, or expansion of exponential for low tau
!!                  ze1 = tauorig / prmu0
!!                  if (ze1 .le. od_lo) then
!!                     zdbtmo = 1._rb - ze1 + 0.5_rb * ze1 * ze1
!!                  else
!!                     tblind = ze1 / (bpade + ze1)
!!                     itind = tblint * tblind + 0.5_rb
!!                     zdbtmo = exp_tbl(itind)
!!                  endif
!!
!!                  zdbt_nodel(jk) = zclear*zdbtmc + zcloud*zdbtmo
!!                  ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)
!!
!!               endif
!   /\/\/\ Above code only needed for unscaled direct beam calculation


! Delta scaling - clear   
               zf = zgcc(jk) * zgcc(jk)
               zwf = zomcc(jk) * zf
               ztauc(jk) = (1.0_rb - zwf) * ztauc(jk)
               zomcc(jk) = (zomcc(jk) - zwf) / (1.0_rb - zwf)
               zgcc (jk) = (zgcc(jk) - zf) / (1.0_rb - zf)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
!               if (icpr .ge. 1) then
               if (lrtchkcld(jk)) then
                  ztauo(jk) = ztauc(jk) + ptaucmc(ikl,iw)
                  zomco(jk) = ztauc(jk) * zomcc(jk) + ptaucmc(ikl,iw) * pomgcmc(ikl,iw) 
                  zgco (jk) = (ptaucmc(ikl,iw) * pomgcmc(ikl,iw) * pasycmc(ikl,iw) + &
                              ztauc(jk) * zomcc(jk) * zgcc(jk)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)
               endif

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
!               elseif (icpr .eq. 0) then
!                  ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm) + ptaucmc(ikl,iw)
!                  zomco(jk) = ptaua(ikl,ibm) * pomga(ikl,ibm) + ptaucmc(ikl,iw) * pomgcmc(ikl,iw) + &
!                              ztaur(ikl,iw) * 1.0_rb
!                  zgco (jk) = (ptaucmc(ikl,iw) * pomgcmc(ikl,iw) * pasycmc(ikl,iw) + &
!                              ptaua(ikl,ibm)*pomga(ikl,ibm)*pasya(ikl,ibm)) / zomco(jk)
!                  zomco(jk) = zomco(jk) / ztauo(jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
!                  zf = zgco(jk) * zgco(jk)
!                  zwf = zomco(jk) * zf
!                  ztauo(jk) = (1._rb - zwf) * ztauo(jk)
!                  zomco(jk) = (zomco(jk) - zwf) / (1.0_rb - zwf)
!                  zgco (jk) = (zgco(jk) - zf) / (1.0_rb - zf)
!               endif 

! End of layer loop
            enddo    

! Clear sky reflectivities
            call reftra_sw (klev, &
                            lrtchkclr, zgcc, prmu0, prmu0m1, ztauc, zomcc, &
                            zrefc, zrefdc, ztrac, ztradc)

! Total sky reflectivities      
            call reftra_sw (klev, &
                            lrtchkcld, zgco, prmu0, prmu0m1, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado)

            do jk=1,klev

               if (lrtchkcld(jk)) then
! Combine clear and cloudy contributions for total sky
                  ikl = klev+1-jk 
                  zclear = 1.0_rb - pcldfmc(ikl,iw)
                  zcloud = pcldfmc(ikl,iw)

                  zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
                  zrefd(jk)= zclear*zrefdc(jk) + zcloud*zrefdo(jk)
                  ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
                  ztrad(jk)= zclear*ztradc(jk) + zcloud*ztrado(jk)
               else
                  zref(jk) = zrefc(jk)
                  zrefd(jk)= zrefdc(jk)
                  ztra(jk) = ztrac(jk)
                  ztrad(jk)= ztradc(jk)
               endif

! Direct beam transmittance        

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) * prmu0m1
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._rb - ze1 + 0.5_rb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_rb
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (lrtchkcld(jk)) then

                  ze1 = ztauo(jk) * prmu0m1
                  if (ze1 .le. od_lo) then
                     zdbtmo = 1._rb - ze1 + 0.5_rb * ze1 * ze1
                  else
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_rb
                     zdbtmo = exp_tbl(itind)
                  endif

                  zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
                  ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)

               else
                  
                  zdbt (jk)   = zdbtmc
                  ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)

               endif
        
            enddo
                 
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw(klev, iw, &
                           zrefc, zrefdc, ztrac, ztradc, &
                           zdbtc, zrdndc, zrupc, zrupdc, ztdbtc, &
                           zcd, zcu)
      
! Vertical quadrature for cloudy fluxes

            if (any(lrtchkcld(1:klev))) then

               call vrtqdr_sw(klev, iw, &
                              zref, zrefd, ztra, ztrad, &
                              zdbt, zrdnd, zrup, zrupd, ztdbt, &
                              zfd, zfu)
            else
               
               zfd(1:klev+1,iw) = zcd(1:klev+1,iw)
               zfu(1:klev+1,iw) = zcu(1:klev+1,iw)

            endif

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
               zbbfu   (ikl,ibm) = zbbfu   (ikl,ibm) + zincflx(iw)*zfu   (jk,iw)  
               zbbfd   (ikl,ibm) = zbbfd   (ikl,ibm) + zincflx(iw)*zfd   (jk,iw)
               zbbcu   (ikl,ibm) = zbbcu   (ikl,ibm) + zincflx(iw)*zcu   (jk,iw)
               zbbcd   (ikl,ibm) = zbbcd   (ikl,ibm) + zincflx(iw)*zcd   (jk,iw)
               zbbfddir(ikl,ibm) = zbbfddir(ikl,ibm) + zincflx(iw)*ztdbt (jk)
               zbbcddir(ikl,ibm) = zbbcddir(ikl,ibm) + zincflx(iw)*ztdbtc(jk)

! Accumulate spectral fluxes over whole spectrum  
!!               pbbfu(ikl) = pbbfu(ikl) + zincflx(iw)*zfu(jk,iw)
!!               pbbfd(ikl) = pbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
!!               pbbcu(ikl) = pbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
!!               pbbcd(ikl) = pbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
!!               if (idelm .eq. 0) then 
!!                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!!                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
!!               elseif (idelm .eq. 1) then
!!                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt(jk)
!!                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc(jk)
!!               endif

! Accumulate direct fluxes for UV/visible bands
!!               if (ibm >= 10 .and. ibm <= 13) then
!!                  puvcd(ikl) = puvcd(ikl) + zincflx(iw)*zcd(jk,iw)
!!                  puvfd(ikl) = puvfd(ikl) + zincflx(iw)*zfd(jk,iw)
!!                  if (idelm .eq. 0) then 
!!                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!!                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
!!                  elseif (idelm .eq. 1) then
!!                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt(jk)
!!                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc(jk)
!!                  endif
! Accumulate direct fluxes for near-IR bands
!!               else if (ibm == 14 .or. ibm <= 9) then  
!!                  pnicd(ikl) = pnicd(ikl) + zincflx(iw)*zcd(jk,iw)
!!                  pnifd(ikl) = pnifd(ikl) + zincflx(iw)*zfd(jk,iw)
!!                  if (idelm .eq. 0) then 
!!                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!!                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
!!                  elseif (idelm .eq. 1) then
!!                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt(jk)
!!                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc(jk)
!!                  endif
!!               endif

            enddo

! End loop on jg, g-point interval
         enddo             

! End loop on jb, spectral band
      enddo                    

    contains

      subroutine reftra_sw(nlayers, lrtchk, pgg, prmuz, prmuzm1, ptau, pw, &
                           pref, prefd, ptra, ptrad)

! --------------------------------------------------------------------
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                         ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (nlayers)

      real(kind=rb), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: ptau(:)                     ! optical depth
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: pw(:)                       ! single scattering albedo 
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: prmuz                       ! cosine of solar zenith angle
      real(kind=rb), intent(in) :: prmuzm1                     ! cosine of solar zenith angle

! ------- Output -------

      real(kind=rb), intent(inout) :: pref(:)                  ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prefd(:)                 ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptra(:)                  ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptrad(:)                 ! diffuse beam transmissivity
                                                               !   Dimensions: (nlayers+1)

! ------- Local -------

      integer(kind=im) :: jk, jl, kmodts
      integer(kind=im) :: itind

      real(kind=rb) :: tblind
      real(kind=rb) :: za, za1, za2
      real(kind=rb) :: zbeta, zdend, zdenr, zdent
      real(kind=rb) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=rb) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=rb) :: zr1, zr2, zr3, zr4, zr5
      real(kind=rb) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=rb) :: zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=rb) :: zw, zwo

      real(kind=rb), parameter :: eps = 1.e-08_rb
      real(kind=rb), parameter :: zsr3 = sqrt(3._rb)
      real(kind=rb), parameter :: zwcrit = 0.9999995_rb
!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision: 11661 $'

!     zsr3=sqrt(3._rb)
!     zwcrit=0.9999995_rb
!     kmodts=2

      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._rb
            ptra(jk) =1._rb
            prefd(jk)=0._rb
            ptrad(jk)=1._rb
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)  

! General two-stream expressions

            zg3= 3._rb * zg
!            if (kmodts == 1) then
!               zgamma1= (7._rb - zw * (4._rb + zg3)) * 0.25_rb
!               zgamma2=-(1._rb - zw * (4._rb - zg3)) * 0.25_rb
!               zgamma3= (2._rb - zg3 * prmuz ) * 0.25_rb
!            else if (kmodts == 2) then  
               zgamma1= (8._rb - zw * (5._rb + zg3)) * 0.25_rb
               zgamma2=  3._rb *(zw * (1._rb - zg )) * 0.25_rb
               zgamma3= (2._rb - zg3 * prmuz ) * 0.25_rb
!            else if (kmodts == 3) then  
!               zgamma1= zsr3 * (2._rb - zw * (1._rb + zg)) * 0.5_rb
!               zgamma2= zsr3 * zw * (1._rb - zg ) * 0.5_rb
!               zgamma3= (1._rb - zsr3 * zg * prmuz ) * 0.5_rb
!            end if
            zgamma4= 1._rb - zgamma3

! Recompute original s.s.a. to test for conservative solution

!           zwo= zw / (1._rb - (1._rb - zw) * (zg / (1._rb - zg))**2)
            zwo= zw * (1._rb - zg)**2 / ( 1._rb - 2._rb*zg + zw*zg**2 )

            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 * prmuzm1 , 500._rb)
!               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  ze2 = 1._rb - ze1 + 0.5_rb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_rb
                  ze2 = exp_tbl(itind)
               endif
!

               pref(jk) = (zgt - za1 * (1._rb - ze2)) / (1._rb + zgt)
               ptra(jk) = 1._rb - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._rb + zgt)
               ptrad(jk) = 1._rb - prefd(jk)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0_rb) then 
                  pref(jk) = 0.0_rb
                  ptra(jk) = 1.0_rb
                  prefd(jk) = 0.0_rb
                  ptrad(jk) = 1.0_rb
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._rb + zrp
               zrm1 = 1._rb - zrp
               zrk2 = 2._rb * zrk
               zrpp = 1._rb - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5

! mji - reformulated code to avoid potential floating point exceptions
!               zbeta = - zr5 / zr4
               zbeta = (zgamma1 - zrk) / zrkg
!!
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._rb)
               ze2 = min ( zto1 * prmuzm1 , 500._rb)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
!              zep1 = exp( ze1 )
!              zem1 = 1._rb / zep1
!              zep2 = exp( ze2 )
!              zem2 = 1._rb / zep2
!
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  zem1 = 1._rb - ze1 + 0.5_rb * ze1 * ze1
                  zep1 = 1._rb / zem1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_rb
                  zem1 = exp_tbl(itind)
                  zep1 = 1._rb / zem1
               endif

               if (ze2 .le. od_lo) then 
                  zem2 = 1._rb - ze2 + 0.5_rb * ze2 * ze2
                  zep2 = 1._rb / zem2
               else
                  tblind = ze2 / (bpade + ze2)
                  itind = tblint * tblind + 0.5_rb
                  zem2 = exp_tbl(itind)
                  zep2 = 1._rb / zem2
               endif

! collimated beam

! mji - reformulated code to avoid potential floating point exceptions
!               zdenr = zr4*zep1 + zr5*zem1
!               pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
!               zdent = zt4*zep1 + zt5*zem1
!               ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent

               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(jk) = eps
                  ptra(jk) = zem2
               else 
                  pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif
!!

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._rb / ( (1._rb - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._rb - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

         endif         

      enddo    

      end subroutine reftra_sw


      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer(kind=im), intent (in) :: klev                   ! number of model layers
      integer(kind=im), intent (in) :: kw                     ! g-point index

      real(kind=rb), intent(in) :: pref(:)                    ! direct beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: prefd(:)                   ! diffuse beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptra(:)                    ! direct beam transmissivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptrad(:)                   ! diffuse beam transmissivity
                                                              !   Dimensions: (nlayers+1)

      real(kind=rb), intent(in) :: pdbt(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptdbt(:)
                                                              !   Dimensions: (nlayers+1)

      real(kind=rb), intent(inout) :: prdnd(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prup(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prupd(:)
                                                              !   Dimensions: (nlayers+1)

! Output
      real(kind=rb), intent(out) :: pfd(:,:)                  ! downwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real(kind=rb), intent(out) :: pfu(:,:)                  ! upwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle

! Local

      integer(kind=im) :: ikp, ikx, jk

      real(kind=rb) :: zreflect
      real(kind=rb) :: ztdn(klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      zreflect = 1._rb / (1._rb - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._rb / (1._rb -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._rb
      prdnd(1) = 0._rb
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._rb / (1._rb - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._rb / (1._rb - prdnd(jk) * prupd(jk))
         pfu(jk,kw) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(jk,kw) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

      end subroutine vrtqdr_sw

      end subroutine spcvmc_sw

      end module rrtmg_sw_spcvmc

