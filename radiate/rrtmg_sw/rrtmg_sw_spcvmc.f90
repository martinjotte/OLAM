!     Path:      $Source$
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

      private
      public :: spcvmc_sw_noclr

      contains

! ---------------------------------------------------------------------------
      subroutine spcvmc_sw_noclr &
            (nlayers, iaer, &
             nsfc, frac_sfck, colmol, &
             palbd, palbp, &
             ztaug, zsflxzen, &
             ztauc, zasyc, zomgc, cldf, &
             ptaua, pasya, pomga, &
             prmu0, adjflux, &
             swuflx, swdflx, &
             zbbfd_sfc, zbbfu_sfc, zbbfddir_sfc)
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

      use parkind,  only: im => kind_im, rb => kind_rb, cldmin
      use parrrsw,  only: nbndsw, ngptsw, mxmol, jpband, jpb1, jpb2
      use rrsw_wvn, only: ngc, ngs, ngs1, raylt

      implicit none

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers
      integer(kind=im), intent(in) :: iaer
      integer(kind=im), intent(in) :: nsfc

      real(kind=rb), intent(in) :: frac_sfck(nsfc)
      real(kind=rb), intent(in) :: palbd(nsfc,nbndsw)      ! surface albedo (diffuse)
      real(kind=rb), intent(in) :: palbp(nsfc,nbndsw)      ! surface albedo (direct)

      real(kind=rb), intent(in) :: ztauc(ngptsw,nlayers)   ! cloud optical depth [mcica]
      real(kind=rb), intent(in) :: zasyc(ngptsw,nlayers)   ! cloud asymmetry parameter [mcica]
      real(kind=rb), intent(in) :: zomgc(ngptsw,nlayers)   ! cloud single scattering albedo [mcica]

      real(kind=rb), intent(in) :: cldf(nlayers)           ! cloud fraction
      real(kind=rb), intent(in) :: colmol(nlayers)

      real(kind=rb), intent(in) :: ptaua(nbndsw,nlayers)   ! aerosol optical depth
      real(kind=rb), intent(in) :: pasya(nbndsw,nlayers)   ! aerosol asymmetry parameter
      real(kind=rb), intent(in) :: pomga(nbndsw,nlayers)   ! aerosol single scattering albedo

      real(kind=rb), intent(in) :: prmu0                   ! cosine of solar zenith angle
      real(kind=rb), intent(in) :: adjflux                 ! Earth/Sun distance adjustment

! Arrays from rrtmg_sw_taumol routines

      real(kind=rb), intent(in) :: ztaug(ngptsw,nlayers)
      real(kind=rb), intent(in) :: zsflxzen(ngptsw)

! ------- Output -------

      real(kind=rb), intent(out) :: swuflx(nlayers)
      real(kind=rb), intent(out) :: swdflx(nlayers)

      real(kind=rb), intent(out) :: zbbfd_sfc   (nsfc,nbndsw), zbbfu_sfc(nsfc,nbndsw)
      real(kind=rb), intent(out) :: zbbfddir_sfc(nsfc,nbndsw)

! ------- Local -------

      integer(kind=im) :: ibm, ikl
      integer(kind=im) :: ig, jb, jk, js

      real(kind=rb) :: zref (ngptsw,nlayers+1)  ! direct albedo
      real(kind=rb) :: zrefd(ngptsw,nlayers+1)  ! diffuse albedo
      real(kind=rb) :: ztra (ngptsw,nlayers+1)  ! direct transmittance
      real(kind=rb) :: ztrad(ngptsw,nlayers+1)  ! diffuse transmittance
      real(kind=rb) :: zdbt (ngptsw,nlayers+1)  ! direct beam transmittance
      real(kind=rb) :: ztdbt(ngptsw,nlayers+1)  ! total direct beam transmittance

      real(kind=rb) :: ztauo(ngptsw)
      real(kind=rb) :: zomco(ngptsw)
      real(kind=rb) :: zggco(ngptsw)

      real(kind=rb) :: zincflx(ngptsw)

      real(kind=rb) :: zf, zwf, ap
      real(kind=rb) :: prmu0m1, zom, ztaur

      real(kind=rb) :: ztaua(ngptsw,nlayers)   ! aerosol optical depth
      real(kind=rb) :: zasya(ngptsw,nlayers)   ! aerosol asymmetry parameter
      real(kind=rb) :: zomga(ngptsw,nlayers)   ! aerosol single scattering albedo

      real(kind=rb), parameter :: onem = 1.0 - 3.0 * epsilon(1.0)

      ! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=rb) :: zfd(ngptsw,nlayers+1)
      real(kind=rb) :: zfu(ngptsw,nlayers+1)

      real(kind=rb) :: zfd_sfc(ngptsw,nlayers-nsfc+2:nlayers+1)
      real(kind=rb) :: zfu_sfc(ngptsw,nlayers-nsfc+2:nlayers+1)

      real(kind=rb) :: zrup_sfc (ngptsw,nlayers-nsfc+2:nlayers+1)
      real(kind=rb) :: zrupd_sfc(ngptsw,nlayers-nsfc+2:nlayers+1)

      prmu0m1 = 1.0_rb / prmu0

      ap = adjflux * prmu0

      do ig = 1, ngptsw
         ztdbt(ig,1)=1.0_rb
         zincflx(ig) = ap * zsflxzen(ig)
      enddo

      do jb = jpb1, jpb2
         ibm = jb-15
         do jk = nlayers-nsfc+2, nlayers+1
            js = nlayers+2-jk
            do ig = ngs1(ibm), ngs(ibm)
               zrup_sfc (ig,jk) = palbp(js,ibm)
               zrupd_sfc(ig,jk) = palbd(js,ibm)
            enddo
         enddo

         if (iaer > 0) then
            do jk=1,nlayers
               do ig = ngs1(ibm), ngs(ibm)
                  ztaua(ig,jk) = ptaua(ibm,jk)
                  zasya(ig,jk) = pasya(ibm,jk)
                  zomga(ig,jk) = pomga(ibm,jk)
               enddo
            enddo
         endif

      enddo

      ! Note: two-stream calculations proceed from top to bottom;
      ! RRTMG_SW quantities are given bottom to top and are reversed here

      do jk=1,nlayers
         ikl=nlayers+1-jk

         ! Total sky optical parameters

         if (iaer <= 0) then

            if (cldf(ikl) <= cldmin) then

               do ig = 1, ngptsw
                  ztaur = colmol(ikl) * raylt(ig)

                  ztauo(ig) = ztaur + ztaug(ig,ikl)
                  zomco(ig) = min(ztaur / ztauo(ig), onem)
                  zggco(ig) = 0._rb
               enddo

            else

               do ig = 1, ngptsw
                  ztaur = colmol(ikl) * raylt(ig)

                  ztauo(ig) = ztaur + ztaug(ig,ikl) + ztauc(ig,ikl)
                  zom       = ztaur + zomgc(ig,ikl)
                  zggco(ig) = zasyc(ig,ikl) / zom
                  zomco(ig) = min(zom / ztauo(ig), onem)
               enddo
            endif

         else

            if (cldf(ikl) <= cldmin) then
               do ig = 1, ngptsw
                  ztaur = colmol(ikl) * raylt(ig)

                  ztauo(ig) = ztaur + ztaug(ig,ikl) + ztaua(ig,ikl)
                  zom       = ztaur + zomga(ig,ikl)
                  zggco(ig) = zasya(ig,ikl) / zom
                  zomco(ig) = min(zom / ztauo(ig), onem)
               enddo
            else
               do ig = 1, ngptsw
                  ztaur = colmol(ikl) * raylt(ig)

                  ztauo(ig) =  ztaur + ztaug(ig,ikl) + ztaua(ig,ikl) + ztauc(ig,ikl)
                  zom       =  ztaur + zomga(ig,ikl) + zomgc(ig,ikl)
                  zggco(ig) = (zasya(ig,ikl) + zasyc(ig,ikl)) / zom
                  zomco(ig) = min(zom / ztauo(ig), onem)
               enddo
            endif

         endif

         ! Delta scaling - clouds and aerosols

         if (iaer > 0 .or. cldf(ikl) > cldmin) then
            do ig = 1, ngptsw
               zf        = zggco(ig) * zggco(ig)
               zwf       = zomco(ig) * zf
               ztauo(ig) = (1._rb - zwf) * ztauo(ig)
               zomco(ig) = (zomco(ig) - zwf) / (1.0_rb - zwf)
               zggco(ig) = (zggco(ig) - zf ) / (1.0_rb - zf)
            enddo
         endif

         ! Reflectivities

         call reftra_sw (jk, nlayers, zggco, prmu0, prmu0m1, ztauo, &
                         zomco, zref, zrefd, ztra, ztrad, zdbt)

         ! End of layer loop
      enddo

      ! Direct beam transmittance

      do jk = 2, nlayers+1
         !dir$ ivdep
         do ig = 1, ngptsw
            ztdbt(ig,jk) = zdbt(ig,jk-1) * ztdbt(ig,jk-1)
         enddo
      enddo

      ! Vertical quadrature for total-sky fluxes

      call vrtqdr_sw(nlayers, nsfc, zrup_sfc, zrupd_sfc, frac_sfck, &
                     zref, zrefd, ztra, ztrad, zdbt, ztdbt, &
                     zfd, zfu, zfd_sfc, zfu_sfc, zincflx)

      ! Upwelling and downwelling fluxes at levels
      ! Two-stream calculations go from top to bottom;
      ! layer indexing is reversed to go bottom to top for output arrays

      do ikl = 1, nlayers
         jk = nlayers+1-ikl
         swuflx(ikl) = sum( zfu(1:ngptsw,jk) )
         swdflx(ikl) = sum( zfd(1:ngptsw,jk) )
      enddo

      ! Special at surface and with shaved cells

      do ikl = 1, nsfc
         jk = nlayers+2-ikl

         do ig = 1, ngptsw
            ztdbt(ig,jk) = ztdbt(ig,jk) * zincflx(ig)
         enddo

         do ibm = 1, nbndsw
            zbbfu_sfc   (ikl,ibm) = sum( zfu_sfc( ngs1(ibm):ngs(ibm), jk ) )
            zbbfd_sfc   (ikl,ibm) = sum( zfd_sfc( ngs1(ibm):ngs(ibm), jk ) )
            zbbfddir_sfc(ikl,ibm) = sum( ztdbt  ( ngs1(ibm):ngs(ibm), jk ) )
         enddo
      enddo

    end subroutine spcvmc_sw_noclr


! --------------------------------------------------------------------
    subroutine reftra_sw(jk, nlayers, pgg, prmuz, prmuzm1, ptau, pw, &
                         pref, prefd, ptra, ptrad, pdbt)
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

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ngptsw

      implicit none

! ------- Input -------

      integer(kind=im), intent(in) :: jk
      integer(kind=im), intent(in) :: nlayers

      real(kind=rb), intent(in) :: pgg(ngptsw)                ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: ptau(ngptsw)               ! optical depth
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: pw(ngptsw)                 ! single scattering albedo
                                                               !   Dimensions: (nlayers)
      real(kind=rb), intent(in) :: prmuz                       ! cosine of solar zenith angle
      real(kind=rb), intent(in) :: prmuzm1                     ! 1 / cosine of solar zenith angle

! ------- Output -------

      real(kind=rb), intent(inout) :: pref(ngptsw,nlayers+1)          ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prefd(ngptsw,nlayers+1)           ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptra(ngptsw,nlayers+1)          ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: ptrad(ngptsw,nlayers+1)         ! diffuse beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: pdbt(ngptsw,nlayers+1)          ! total transmissivity
                                                               !   Dimensions: (nlayers+1)
! ------- Local -------

      real(kind=rb) :: za1, za2
      real(kind=rb) :: zdend, zdenr, denom
      real(kind=rb) :: zem1, zem2, zemm, zemm1, zemp1
      real(kind=rb) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4
      real(kind=rb) :: zr1, zr2, zr3
      real(kind=rb) :: zrk, zrk2, zrm1, zrp, zrp1
      real(kind=rb) :: zt1, zt2, zt3, zto1, zrpm1
      real(kind=rb) :: zw

      integer(kind=im) :: ig

      real(kind=rb), parameter :: zsr3 = 0.5_rb * sqrt(3._rb)
      real(kind=rb), parameter :: one  = nearest(1._rb, -1._rb)

      do ig = 1, ngptsw

         zto1 = -ptau(ig)
         zw   =  pw  (ig)
         zg   =  pgg (ig)

! General two-stream expressions

         ! kmodts == 1
!        zg3 = 3._rb * zg
!        zgamma1 = (7._rb - zw * (4._rb + zg3)) * 0.25_rb
!        zgamma2 =-(1._rb - zw * (4._rb - zg3)) * 0.25_rb
!        zgamma3 = (2._rb - zg3 * prmuz ) * 0.25_rb

         ! kmodts == 2
         zg3 = 0.75_rb * zg
         zgamma1 = 2._rb - zw * (1.25_rb + zg3)
         zgamma2 = zw * (0.75_rb - zg3)
         zgamma3 = 0.50_rb - zg3 * prmuz

         ! kmodts == 3
!        zgamma1 = zsr3 * (2._rb - zw * (1._rb + zg))
!        zgamma2 = zsr3 * zw * (1._rb - zg )
!        zgamma3 = (0.5_rb - zsr3 * zg * prmuz )

         zgamma4 = 1._rb - zgamma3

! Non-conservative scattering

         za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
         za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
         zrk = sqrt ( zgamma1**2 - zgamma2**2 )
!        zrp = zrk * prmuz

!        Avoid zrp too close to 1
         zrpm1 = zrk * prmuz - 1._rb
         zrp   = 1._rb + sign( max(2.5e-7, abs(zrpm1)), zrpm1 )

         zrp1 = 1._rb + zrp
         zrm1 = 1._rb - zrp

         zrk2 = 2._rb * zrk

         zr1  = zrm1 * (za2 + zrk * zgamma3)
         zr2  = zrp1 * (za2 - zrk * zgamma3)
         zr3  = zrk2 * (zgamma3 - za2 * prmuz )

         zt1  = zrp1 * (za1 + zrk * zgamma4)
         zt2  = zrm1 * (za1 - zrk * zgamma4)
         zt3  = zrk2 * (zgamma4 + za1 * prmuz )

         zem1 = exp( zto1 * zrk)
         zem2 = exp( zto1 * prmuzm1 )

         zemm  = zem1 * zem1
         zemm1 = 1._rb - zemm
         zemp1 = 1._rb + zemm

         pdbt(ig,jk) = zem2

! collimated beam

         zdenr = zw / (zrp1 * zrm1 * (zrk * zemp1 + zgamma1 * zemm1))

         pref(ig,jk) = (zr1 - zr2*zemm - zr3*zem2*zem1) * zdenr
         pref(ig,jk) = max(0., min(1., pref(ig,jk)))

         ptra(ig,jk) = zem2 - (zem2*zt1 - zem2*zt2*zemm - zt3*zem1) * zdenr
         ptra(ig,jk) = max(0., min(1.-pref(ig,jk), ptra(ig,jk)))

! diffuse beam

         zdend = 1.0 / (zrk2 + zemm1 * (zgamma1 - zrk))

         prefd(ig,jk) = zgamma2 * zemm1 * zdend
         ptrad(ig,jk) = zrk2 * zem1 * zdend

      enddo

      end subroutine reftra_sw

! --------------------------------------------------------------------------

      subroutine vrtqdr_sw(nlay, nsfc, prup_sfc, prupd_sfc, frac_sfck, &
                           pref, prefd, ptra, ptrad, pdbt, ptdbt, &
                           pfd, pfu, pfd_sfc, pfu_sfc, pincflx)

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

        use parkind, only: im => kind_im, rb => kind_rb
        use parrrsw, only: ngptsw

        implicit none

! ------- Declarations -------

! Input

      integer(kind=im), intent (in) :: nlay                ! number of model layers

      integer(kind=im), intent (in) :: nsfc

      real(kind=rb), intent(in) :: prup_sfc (ngptsw,nlay-nsfc+2:nlay+1)

      real(kind=rb), intent(in) :: prupd_sfc(ngptsw,nlay-nsfc+2:nlay+1)

      real(kind=rb), intent(in) :: frac_sfck(nsfc)

      real(kind=rb), intent(in) :: pref(ngptsw,nlay+1)     ! direct beam reflectivity

      real(kind=rb), intent(in) :: prefd(ngptsw,nlay+1)    ! diffuse beam reflectivity

      real(kind=rb), intent(in) :: ptra(ngptsw,nlay+1)     ! direct beam transmissivity

      real(kind=rb), intent(in) :: ptrad(ngptsw,nlay+1)    ! diffuse beam transmissivity

      real(kind=rb), intent(in) :: pdbt(ngptsw,nlay+1)     ! layer mean direct transmittance

      real(kind=rb), intent(in) :: ptdbt(ngptsw,nlay+1)    ! total direct transmittance

      real(kind=rb), intent(in) :: pincflx(ngptsw)

! Output
      real(kind=rb), intent(out) :: pfd(ngptsw,nlay+1)     ! downwelling flux (W/m2) unadjusted for
                                                           ! earth/sun distance or zenith angle
      real(kind=rb), intent(out) :: pfu(ngptsw,nlay+1)     ! upwelling flux (W/m2) unadjusted for
                                                           ! earth/sun distance or zenith angle
      real(kind=rb), intent(out) :: pfd_sfc(ngptsw,nlay-nsfc+2:nlay+1)

      real(kind=rb), intent(out) :: pfu_sfc(ngptsw,nlay-nsfc+2:nlay+1)

! Local
      real(kind=rb) :: prdn (ngptsw,nlay+1)
      real(kind=rb) :: prdnd(ngptsw,nlay+1)

      real(kind=rb) :: prup (ngptsw,nlay+1)
      real(kind=rb) :: prupd(ngptsw,nlay+1)

      integer(kind=im) :: ikp, jk, js, ig

      real(kind=rb) :: zreflect, prup_tot(ngptsw), prupd_tot(ngptsw)

!-----------------------------------------------------------------------------

! Link lowest layers with surface

      do ig = 1, ngptsw
         prup (ig,nlay+1) = prup_sfc (ig,nlay+1)
         prupd(ig,nlay+1) = prupd_sfc(ig,nlay+1)

         prup_tot (ig) = prup_sfc (ig,nlay+1)
         prupd_tot(ig) = prupd_sfc(ig,nlay+1)
      enddo

      do js = 2, nsfc+1
         jk = nlay-js+2
         ikp = jk + 1

         do ig = 1, ngptsw
            zreflect = ptrad(ig,jk) / (1._rb - prupd_tot(ig) * prefd(ig,jk))

            prup(ig,jk) = pref(ig,jk) + &
                        ((ptra(ig,jk) - pdbt(ig,jk)) * prupd_tot(ig) + &
                          pdbt(ig,jk) * prup_tot(ig)) * zreflect

            prupd(ig,jk) = prefd(ig,jk) + ptrad(ig,jk) * prupd_tot(ig) * zreflect
         enddo

         if (js <= nsfc) then
            do ig = 1, ngptsw
               prup_tot (ig) = prup_sfc (ig,jk) * frac_sfck(js) + prup (ig,jk) * (1.0 - frac_sfck(js))
               prupd_tot(ig) = prupd_sfc(ig,jk) * frac_sfck(js) + prupd(ig,jk) * (1.0 - frac_sfck(js))
            enddo
         endif
      enddo

! Pass from bottom to top

      do jk = nlay-nsfc, 1, -1
         ikp = jk+1
         do ig = 1, ngptsw
            zreflect = ptrad(ig,jk) / (1._rb - prupd(ig,ikp) * prefd(ig,jk))

            prup(ig,jk) = pref(ig,jk) + &
                        ((ptra(ig,jk) - pdbt(ig,jk)) * prupd(ig,ikp) + &
                          pdbt(ig,jk) * prup(ig,ikp)) * zreflect

            prupd(ig,jk) = prefd(ig,jk) + ptrad(ig,jk) * prupd(ig,ikp) * zreflect
         enddo
      enddo

! Upper boundary conditions

      do ig = 1, ngptsw
         prdn(ig,1) = 1.0_rb
         prdn(ig,2) = ptra(ig,1)

         prdnd(ig,1) = 0.0_rb
         prdnd(ig,2) = prefd(ig,1)
      enddo

! Pass from top to bottom

      do jk = 2,nlay
         ikp = jk+1
         !dir$ ivdep
         do ig = 1, ngptsw
            zreflect = ptrad(ig,jk) / (1._rb - prefd(ig,jk) * prdnd(ig,jk))

            prdn(ig,ikp) = ptdbt(ig,jk) * ptra(ig,jk) + &
                          (ptdbt(ig,jk) * pref(ig,jk) * prdnd(ig,jk) + &
                           prdn(ig,jk) - ptdbt(ig,jk)) * zreflect

            prdnd(ig,ikp) = prefd(ig,jk) + ptrad(ig,jk) * prdnd(ig,jk) * zreflect
         enddo
      enddo

! Up and down-welling fluxes at levels

      do jk = 1,nlay+1
         !dir$ ivdep
         do ig = 1, ngptsw
            zreflect = pincflx(ig) / (1._rb - prdnd(ig,jk) * prupd(ig,jk))

            pfu(ig,jk) = (ptdbt(ig,jk) * prup(ig,jk) + &
                         (prdn(ig,jk) - ptdbt(ig,jk)) * prupd(ig,jk)) * zreflect

            pfd(ig,jk) = ptdbt(ig,jk) * pincflx(ig) + (prdn(ig,jk) - ptdbt(ig,jk) + &
                         ptdbt(ig,jk) * prup(ig,jk) * prdnd(ig,jk)) * zreflect
         enddo
      enddo

      if (nsfc > 1) then
         do jk = nlay-nsfc+2, nlay
            do ig = 1, ngptsw
               zreflect = pincflx(ig) / (1._rb - prdnd(ig,jk) * prupd_sfc(ig,jk))
               pfu_sfc(ig,jk) = (ptdbt(ig,jk) * prup_sfc(ig,jk) + &
                                (prdn(ig,jk) - ptdbt(ig,jk)) * prupd_sfc(ig,jk)) * zreflect
               pfd_sfc(ig,jk) = ptdbt(ig,jk) * pincflx(ig) + (prdn(ig,jk) - ptdbt(ig,jk)+ &
                                ptdbt(ig,jk) * prup_sfc(ig,jk) * prdnd(ig,jk)) * zreflect
            enddo
         enddo
      endif

      do ig = 1, ngptsw
         pfu_sfc(ig,nlay+1) = pfu(ig,nlay+1)
         pfd_sfc(ig,nlay+1) = pfd(ig,nlay+1)
      enddo

      end subroutine vrtqdr_sw


      end module rrtmg_sw_spcvmc
