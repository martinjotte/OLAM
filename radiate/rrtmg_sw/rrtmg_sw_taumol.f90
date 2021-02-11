!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $

      module rrtmg_sw_taumol

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
      public :: taumol_sw

      contains

!----------------------------------------------------------------------------
      subroutine taumol_sw(nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           selffac, selffrac, indself, forfac, forfrac, indfor,&
                           sfluxzen, taug)
!----------------------------------------------------------------------------

! ******************************************************************************
! *                                                                            *
! *                 Optical depths developed for the                           *
! *                                                                            *
! *               RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
! *                                                                            *
! *                                                                            *
! *           ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
! *                       131 HARTWELL AVENUE                                  *
! *                       LEXINGTON, MA 02421                                  *
! *                                                                            *
! *                                                                            *
! *                          ELI J. MLAWER                                     *
! *                        JENNIFER DELAMERE                                   *
! *                        STEVEN J. TAUBMAN                                   *
! *                        SHEPARD A. CLOUGH                                   *
! *                                                                            *
! *                                                                            *
! *                                                                            *
! *                                                                            *
! *                      email:  mlawer@aer.com                                *
! *                      email:  jdelamer@aer.com                              *
! *                                                                            *
! *       The authors wish to acknowledge the contributions of the             *
! *       following people:  Patrick D. Brown, Michael J. Iacono,              *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
! *                                                                            *
! ******************************************************************************
! *    TAUMOL                                                                  *
! *                                                                            *
! *    This file contains the subroutines TAUGBn (where n goes from            *
! *    1 to 28).  TAUGBn calculates the optical depths and Planck fractions    *
! *    per g-value and layer for band n.                                       *
! *                                                                            *
! * Output:  optical depths (unitless)                                         *
! *          fractions needed to compute Planck functions at every layer       *
! *              and g-value                                                   *
! *                                                                            *
! *    COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
! *    COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
! *                                                                            *
! * Input                                                                      *
! *                                                                            *
! *    PARAMETER (MG=16, MXLAY=203, NBANDS=14)                                 *
! *                                                                            *
! *    COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
! *    COMMON /PRECISE/  ONEMINUS                                              *
! *    COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
! *   &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
! *    COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
! *   &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
! *   &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
! *   &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
! *    COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
! *   &                  FAC10(MXLAY),FAC11(MXLAY)                             *
! *    COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
! *    COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
! *                                                                            *
! *    Description:                                                            *
! *    NG(IBAND) - number of g-values in band IBAND                            *
! *    NSPA(IBAND) - for the lower atmosphere, the number of reference         *
! *                  atmospheres that are stored for band IBAND per            *
! *                  pressure level and temperature.  Each of these            *
! *                  atmospheres has different relative amounts of the         *
! *                  key species for the band (i.e. different binary           *
! *                  species parameters).                                      *
! *    NSPB(IBAND) - same for upper atmosphere                                 *
! *    ONEMINUS - since problems are caused in some cases by interpolation     *
! *               parameters equal to or greater than 1, for these cases       *
! *               these parameters are set to this value, slightly < 1.        *
! *    PAVEL - layer pressures (mb)                                            *
! *    TAVEL - layer temperatures (degrees K)                                  *
! *    PZ - level pressures (mb)                                               *
! *    TZ - level temperatures (degrees K)                                     *
! *    LAYTROP - layer at which switch is made from one combination of         *
! *              key species to another                                        *
! *    COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
! *              vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
! *              respectively (molecules/cm**2)                                *
! *    CO2MULT - for bands in which carbon dioxide is implemented as a         *
! *              trace species, this is the factor used to multiply the        *
! *              band's average CO2 absorption coefficient to get the added    *
! *              contribution to the optical depth relative to 355 ppm.        *
! *    FACij(LAY) - for layer LAY, these are factors that are needed to        *
! *                 compute the interpolation factors that multiply the        *
! *                 appropriate reference k-values.  A value of 0 (1) for      *
! *                 i,j indicates that the corresponding factor multiplies     *
! *                 reference k-value for the lower (higher) of the two        *
! *                 appropriate temperatures, and altitudes, respectively.     *
! *    JP - the index of the lower (in altitude) of the two appropriate        *
! *         reference pressure levels needed for interpolation                 *
! *    JT, JT1 - the indices of the lower of the two appropriate reference     *
! *              temperatures needed for interpolation (for pressure           *
! *              levels JP and JP+1, respectively)                             *
! *    SELFFAC - scale factor needed to water vapor self-continuum, equals     *
! *              (water vapor density)/(atmospheric density at 296K and        *
! *              1013 mb)                                                      *
! *    SELFFRAC - factor needed for temperature interpolation of reference     *
! *               water vapor self-continuum data                              *
! *    INDSELF - index of the lower of the two appropriate reference           *
! *              temperatures needed for the self-continuum interpolation      *
! *                                                                            *
! * Data input                                                                 *
! *    COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
! *       (note:  n is the band number)                                        *
! *                                                                            *
! *    Description:                                                            *
! *    KA - k-values for low reference atmospheres (no water vapor             *
! *         self-continuum) (units: cm**2/molecule)                            *
! *    KB - k-values for high reference atmospheres (all sources)              *
! *         (units: cm**2/molecule)                                            *
! *    SELFREF - k-values for water vapor self-continuum for reference         *
! *              atmospheres (used below LAYTROP)                              *
! *              (units: cm**2/molecule)                                       *
! *                                                                            *
! *    DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
! *    EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
! *                                                                            *
! *****************************************************************************
!
! Modifications
!
! Revised: Adapted to F90 coding, J.-J.Morcrette, ECMWF, Feb 2003
! Revised: Modified for g-point reduction, MJIacono, AER, Dec 2003
! Revised: Reformatted for consistency with rrtmg_lw, MJIacono, AER, Jul 2006
!
! ------- Declarations -------

      use parkind,  only: im => kind_im, rb => kind_rb
      use parrrsw,  only: ngptsw
      use rrsw_con, only: oneminus
      use rrsw_wvn, only: nspa, nspb

      implicit none

! ----- Input -----

      integer(kind=im), intent(in) :: nlayers            ! total number of layers

      integer(kind=im), intent(in) :: laytrop            ! tropopause layer index

      integer(kind=im), intent(in) :: jp(nlayers)

      integer(kind=im), intent(in) :: jt(nlayers)

      integer(kind=im), intent(in) :: jt1(nlayers)

      real(kind=rb), intent(in) :: colh2o(nlayers)       ! column amount (h2o)

      real(kind=rb), intent(in) :: colco2(nlayers)       ! column amount (co2)

      real(kind=rb), intent(in) :: colo3(nlayers)        ! column amount (o3)

      real(kind=rb), intent(in) :: colch4(nlayers)       ! column amount (ch4)

      real(kind=rb), intent(in) :: colo2(nlayers)        ! column amount (o2)

      integer(kind=im), intent(in) :: indself(nlayers)

      integer(kind=im), intent(in) :: indfor(nlayers)

      real(kind=rb), intent(in) :: selffac(nlayers)

      real(kind=rb), intent(in) :: selffrac(nlayers)

      real(kind=rb), intent(in) :: forfac(nlayers)

      real(kind=rb), intent(in) :: forfrac(nlayers)

      real(kind=rb), intent(in) :: &
                       fac00(nlayers), fac01(nlayers), &
                       fac10(nlayers), fac11(nlayers)

! ----- Output -----

      real(kind=rb), intent(out) :: sfluxzen(ngptsw)     ! solar source function

      real(kind=rb), intent(out) :: taug(ngptsw,nlayers) ! gaseous optical depth

!     hvrtau = '$Revision: 11661 $'

      real(kind=rb) :: selffac0(nlayers), selffac1(nlayers)
      real(kind=rb) :: forfac0 (nlayers), forfac1 (nlayers)
      real(kind=rb) :: fac00h(nlayers)
      real(kind=rb) :: fac10h(nlayers)
      real(kind=rb) :: fac01h(nlayers)
      real(kind=rb) :: fac11h(nlayers)

      integer(kind=im) :: lay

! pre-compute some quantities

      do lay = 1, nlayers
         forfac0(lay) = colh2o(lay) * forfac(lay) * (1.0 - forfrac(lay))
         forfac1(lay) = colh2o(lay) * forfac(lay) * forfrac(lay)

         fac00h(lay) = fac00(lay) * colh2o(lay)
         fac10h(lay) = fac10(lay) * colh2o(lay)
         fac01h(lay) = fac01(lay) * colh2o(lay)
         fac11h(lay) = fac11(lay) * colh2o(lay)
      enddo

      do lay = 1, laytrop
         selffac0(lay) = colh2o(lay) * selffac(lay) * (1.0 - selffrac(lay))
         selffac1(lay) = colh2o(lay) * selffac(lay) * selffrac(lay)
      enddo

! Calculate gaseous optical depth and planck fractions for each spectral band.

      call taumol16
      call taumol17
      call taumol18
      call taumol19
      call taumol20
      call taumol21
      call taumol22
      call taumol23
      call taumol24
      call taumol25
      call taumol26
      call taumol27
      call taumol28
      call taumol29

!-------------
      contains
!-------------

!----------------------------------------------------------------------------
      subroutine taumol16
!----------------------------------------------------------------------------
!
!     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng16
      use rrsw_kg16, only: absa, absb, forref, selfref, sfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real(kind=rb), parameter :: strrat1 = 252.131_rb

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         speccomb = colh2o(lay) + strrat1*colch4(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(16) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(16) + js

         !dir$ vector always
         do ig = 1, ng16
            taug(ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(16) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(16) + 1

         !dir$ vector always
         do ig = 1, ng16
            taug(ig,lay) = colch4(lay) * &
                (fac00(lay) * absb(ig,ind0  ) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1  ) + &
                 fac11(lay) * absb(ig,ind1+1))
         enddo
      enddo

      do ig = 1, ng16
         sfluxzen(ig) = sfluxref(ig)
      enddo

      end subroutine taumol16

!----------------------------------------------------------------------------
      subroutine taumol17
!----------------------------------------------------------------------------
!
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng17, ngs16
      use rrsw_kg17, only: absa, absb, forref, selfref, sfluxref, &
                           dsfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 0.364641_rb
      integer(kind=im), parameter :: layreffr = 30

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         speccomb = colh2o(lay) + strrat*colco2(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(17) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(17) + js

         !dir$ vector always
         do ig = 1, ng17
            taug(ngs16+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo
      enddo

      laysolfr = nlayers

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         if (jp(lay-1) .lt. layreffr .and. jp(lay) .ge. layreffr) &
            laysolfr = lay
         speccomb = colh2o(lay) + strrat*colco2(lay)
         specmult = 4._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(17) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(17) + js

         !dir$ vector always
         do ig = 1, ng17
            taug(ngs16+ig,lay) = speccomb * &
                (fac000 * absb(ig,ind0) + &
                 fac100 * absb(ig,ind0+1) + &
                 fac010 * absb(ig,ind0+5) + &
                 fac110 * absb(ig,ind0+6) + &
                 fac001 * absb(ig,ind1) + &
                 fac101 * absb(ig,ind1+1) + &
                 fac011 * absb(ig,ind1+5) + &
                 fac111 * absb(ig,ind1+6)) + &

                 forfac0(lay) * forref(ig,3) + &
                 forfac1(lay) * forref(ig,4)
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng17
               sfluxzen(ngs16+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

      end subroutine taumol17

!----------------------------------------------------------------------------
      subroutine taumol18
!----------------------------------------------------------------------------
!
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng18, ngs17
      use rrsw_kg18, only: absa, absb, forref, selfref, sfluxref, &
                           rayl, dsfluxref

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 38.9589_rb
      integer(kind=im), parameter :: layreffr = 6

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

      laysolfr = laytrop

! Lower atmosphere loop
      do lay = 1, laytrop
         if (jp(lay) .lt. layreffr .and. jp(lay+1) .ge. layreffr) &
            laysolfr = min(lay+1,laytrop)
         speccomb = colh2o(lay) + strrat*colch4(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(18) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(18) + js

         !dir$ vector always
         do ig = 1, ng18
            taug(ngs17+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng18
               sfluxzen(ngs17+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(18) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(18) + 1

         !dir$ vector always
         do ig = 1, ng18
            taug(ngs17+ig,lay) = colch4(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))
         enddo
       enddo

       end subroutine taumol18

!----------------------------------------------------------------------------
      subroutine taumol19
!----------------------------------------------------------------------------
!
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng19, ngs18
      use rrsw_kg19, only: absa, absb, forref, selfref, sfluxref, &
                           dsfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 5.49281_rb
      integer(kind=im), parameter :: layreffr = 3

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

      laysolfr = laytrop

! Lower atmosphere loop
      do lay = 1, laytrop
         if (jp(lay) .lt. layreffr .and. jp(lay+1) .ge. layreffr) &
            laysolfr = min(lay+1,laytrop)
         speccomb = colh2o(lay) + strrat*colco2(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(19) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(19) + js

         !dir$ vector always
         do ig = 1 , ng19
            taug(ngs18+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng19
               sfluxzen(ngs18+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(19) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(19) + 1

         !dir$ vector always
         do ig = 1 , ng19
            taug(ngs18+ig,lay) = colco2(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))
         enddo
      enddo

      end subroutine taumol19

!----------------------------------------------------------------------------
      subroutine taumol20
!----------------------------------------------------------------------------
!
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng20, ngs19
      use rrsw_kg20, only: absa, absb, forref, selfref, sfluxref, &
                           absch4, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult, specparm

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(20) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(20) + 1

         !dir$ vector always
         do ig = 1, ng20
            taug(ngs19+ig,lay) = &
                 fac00h(lay) * absa(ig,ind0) + &
                 fac10h(lay) * absa(ig,ind0+1) + &
                 fac01h(lay) * absa(ig,ind1) + &
                 fac11h(lay) * absa(ig,ind1+1) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1) + &

                 colch4(lay) * absch4(ig)
         enddo
      enddo

      do ig = 1, ng20
         sfluxzen(ngs19+ig) = sfluxref(ig)
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(20) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(20) + 1

         !dir$ vector always
         do ig = 1, ng20
            taug(ngs19+ig,lay) = &
                 fac00h(lay) * absb(ig,ind0) + &
                 fac10h(lay) * absb(ig,ind0+1) + &
                 fac01h(lay) * absb(ig,ind1) + &
                 fac11h(lay) * absb(ig,ind1+1) + &

                 forfac0(lay) * forref(ig,3) + &
                 forfac1(lay) * forref(ig,4) + &

                 colch4(lay) * absch4(ig)
         enddo
      enddo

      end subroutine taumol20

!----------------------------------------------------------------------------
      subroutine taumol21
!----------------------------------------------------------------------------
!
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng21, ngs20
      use rrsw_kg21, only: absa, absb, forref, selfref, sfluxref, &
                           rayl, dsfluxref

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 0.0045321_rb
      integer(kind=im), parameter :: layreffr = 8

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

      laysolfr = laytrop

! Lower atmosphere loop
      do lay = 1, laytrop
         if (jp(lay) .lt. layreffr .and. jp(lay+1) .ge. layreffr) &
            laysolfr = min(lay+1,laytrop)
         speccomb = colh2o(lay) + strrat*colco2(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb , oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(21) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(21) + js

         !dir$ vector always
         do ig = 1, ng21
            taug(ngs20+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng21
               sfluxzen(ngs20+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         speccomb = colh2o(lay) + strrat*colco2(lay)
         specmult = 4._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(21) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(21) + js

         !dir$ vector always
         do ig = 1, ng21
            taug(ngs20+ig,lay) = speccomb * &
                (fac000 * absb(ig,ind0) + &
                 fac100 * absb(ig,ind0+1) + &
                 fac010 * absb(ig,ind0+5) + &
                 fac110 * absb(ig,ind0+6) + &
                 fac001 * absb(ig,ind1) + &
                 fac101 * absb(ig,ind1+1) + &
                 fac011 * absb(ig,ind1+5) + &
                 fac111 * absb(ig,ind1+6)) + &

                 forfac0(lay) * forref(ig,3) + &
                 forfac1(lay) * forref(ig,4)
         enddo
      enddo

      end subroutine taumol21

!----------------------------------------------------------------------------
      subroutine taumol22
!----------------------------------------------------------------------------
!
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng22, ngs21
      use rrsw_kg22, only: absa, absb, forref, selfref, sfluxref, &
                           rayl, dsfluxref

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real(   kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult, o2cont

! The following factor is the ratio of total O2 band intensity (lines
! and Mate continuum) to O2 band intensity (line only).  It is needed
! to adjust the optical depths since the k's include only lines.

      real   (kind=rb), parameter :: o2adj = 1.6_rb
      real   (kind=rb), parameter :: strrat = 0.022708_rb
      integer(kind=im), parameter :: layreffr = 2

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

      laysolfr = laytrop

! Lower atmosphere loop
      do lay = 1, laytrop
         if (jp(lay) .lt. layreffr .and. jp(lay+1) .ge. layreffr) &
            laysolfr = min(lay+1,laytrop)
         o2cont = 4.35e-4_rb*colo2(lay)/(350.0_rb*2.0_rb)
         speccomb = colh2o(lay) + o2adj*strrat*colo2(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(22) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(22) + js

         !dir$ vector always
         do ig = 1, ng22
            taug(ngs21+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1) + &

                 o2cont
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng22
               sfluxzen(ngs21+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         o2cont = 4.35e-4_rb*colo2(lay)/(350.0_rb*2.0_rb)
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(22) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(22) + 1

         !dir$ vector always
         do ig = 1, ng22
            taug(ngs21+ig,lay) = colo2(lay) * o2adj * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) + &
                 o2cont
         enddo
      enddo

      end subroutine taumol22

!----------------------------------------------------------------------------
      subroutine taumol23
!----------------------------------------------------------------------------
!
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng23, ngs22
      use rrsw_kg23, only: absa, forref, selfref, sfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult, specparm

! Average Giver et al. correction factor for this band.

      real(kind=rb), parameter :: givfac = 1.029_rb

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(23) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(23) + 1

         !dir$ vector always
         do ig = 1, ng23
            taug(ngs22+ig,lay) = givfac * &
                (fac00h(lay) * absa(ig,ind0) + &
                 fac10h(lay) * absa(ig,ind0+1) + &
                 fac01h(lay) * absa(ig,ind1) + &
                 fac11h(lay) * absa(ig,ind1+1)) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo
      enddo

      do ig = 1, ng23
         sfluxzen(ngs22+ig) = sfluxref(ig)
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         do ig = 1, ng23
            taug(ngs22+ig,lay) = 0._rb
         enddo
      enddo

      end subroutine taumol23

!----------------------------------------------------------------------------
      subroutine taumol24
!----------------------------------------------------------------------------
!
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng24, ngs23
      use rrsw_kg24, only: absa, absb, forref, selfref, sfluxref, &
                           dsfluxref, abso3a, abso3b, rayla, raylb, drayla

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 0.124692_rb
      integer(kind=im), parameter :: layreffr = 1

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

      laysolfr = laytrop

! Lower atmosphere loop
      do lay = 1, laytrop
         if (jp(lay) .lt. layreffr .and. jp(lay+1) .ge. layreffr) &
            laysolfr = min(lay+1,laytrop)
         speccomb = colh2o(lay) + strrat*colo2(lay)
         specmult = 8._rb * min(colh2o(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(24) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(24) + js

         !dir$ vector always
         do ig = 1, ng24
            taug(ngs23+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10)) + &

                 colo3(lay) * abso3a(ig) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1)
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng24
               sfluxzen(ngs23+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(24) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(24) + 1

         !dir$ vector always
         do ig = 1, ng24
            taug(ngs23+ig,lay) = colo2(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) + &
                 colo3(lay) * abso3b(ig)
         enddo
      enddo

      end subroutine taumol24

!----------------------------------------------------------------------------
      subroutine taumol25
!----------------------------------------------------------------------------
!
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng25, ngs24
      use rrsw_kg25, only: absa, sfluxref, abso3a, abso3b, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult, specparm

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(25) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(25) + 1

         !dir$ vector always
         do ig = 1, ng25
            taug(ngs24+ig,lay) = &
                 fac00h(lay) * absa(ig,ind0) + &
                 fac10h(lay) * absa(ig,ind0+1) + &
                 fac01h(lay) * absa(ig,ind1) + &
                 fac11h(lay) * absa(ig,ind1+1) + &
                 colo3(lay) * abso3a(ig)
         enddo
      enddo

      do ig = 1, ng25
         sfluxzen(ngs24+ig) = sfluxref(ig)
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         do ig = 1, ng25
            taug(ngs24+ig,lay) = colo3(lay) * abso3b(ig)
         enddo
      enddo

      end subroutine taumol25

!----------------------------------------------------------------------------
      subroutine taumol26
!----------------------------------------------------------------------------
!
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng26, ngs25
      use rrsw_kg26, only: sfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, lay

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower/Upper atmosphere loop

      do lay = 1, nlayers
         do ig = 1, ng26
            taug(ngs25+ig,lay) = 0._rb
         enddo
      enddo

      do ig = 1, ng26
         sfluxzen(ngs25+ig) = sfluxref(ig)
      enddo

      end subroutine taumol26

!----------------------------------------------------------------------------
      subroutine taumol27
!----------------------------------------------------------------------------
!
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng27, ngs26
      use rrsw_kg27, only: absa, absb, sfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, lay

      real(kind=rb), parameter :: scalekur = 50.15_rb/48.37_rb

! Kurucz solar source function
! The values in sfluxref were obtained using the "low resolution"
! version of the Kurucz solar source function.  For unknown reasons,
! the total irradiance in this band differs from the corresponding
! total in the "high-resolution" version of the Kurucz function.
! Therefore, these values are scaled below by the factor SCALEKUR.

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(27) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(27) + 1

         !dir$ vector always
         do ig = 1, ng27
            taug(ngs26+ig,lay) = colo3(lay) * &
                (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(27) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(27) + 1

         !dir$ vector always
         do ig = 1, ng27
            taug(ngs26+ig,lay) = colo3(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))
         enddo
      enddo

      do ig = 1, ng27
         sfluxzen(ngs26+ig) = scalekur * sfluxref(ig)
      enddo

      end subroutine taumol27

!----------------------------------------------------------------------------
      subroutine taumol28
!----------------------------------------------------------------------------
!
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng28, ngs27
      use rrsw_kg28, only: absa, absb, sfluxref, dsfluxref, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay, laysolfr, isp
      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult

      real   (kind=rb), parameter :: strrat = 6.67029e-07_rb
      integer(kind=im), parameter :: layreffr = 58

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         speccomb = colo3(lay) + strrat*colo2(lay)
         specmult = 8._rb * min(colo3(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(28) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(28) + js

         !dir$ vector always
         do ig = 1, ng28
            taug(ngs27+ig,lay) = speccomb * &
                (fac000 * absa(ig,ind0) + &
                 fac100 * absa(ig,ind0+1) + &
                 fac010 * absa(ig,ind0+9) + &
                 fac110 * absa(ig,ind0+10) + &
                 fac001 * absa(ig,ind1) + &
                 fac101 * absa(ig,ind1+1) + &
                 fac011 * absa(ig,ind1+9) + &
                 fac111 * absa(ig,ind1+10))
         enddo
      enddo

      laysolfr = nlayers

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         if (jp(lay-1) .lt. layreffr .and. jp(lay) .ge. layreffr) &
            laysolfr = lay
         speccomb = colo3(lay) + strrat*colo2(lay)
         specmult = 4._rb * min(colo3(lay)/speccomb, oneminus)
         isp = int(specmult)
         js = 1 + isp
         fs = specmult - isp
         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs) * fac01(lay)
         fac011 = (1._rb - fs) * fac11(lay)
         fac101 = fs * fac01(lay)
         fac111 = fs * fac11(lay)
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(28) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(28) + js

         !dir$ vector always
         do ig = 1, ng28
            taug(ngs27+ig,lay) = speccomb * &
                (fac000 * absb(ig,ind0) + &
                 fac100 * absb(ig,ind0+1) + &
                 fac010 * absb(ig,ind0+5) + &
                 fac110 * absb(ig,ind0+6) + &
                 fac001 * absb(ig,ind1) + &
                 fac101 * absb(ig,ind1+1) + &
                 fac011 * absb(ig,ind1+5) + &
                 fac111 * absb(ig,ind1+6))
         enddo

         if (lay .eq. laysolfr) then
            !dir$ vector always
            do ig = 1, ng28
               sfluxzen(ngs27+ig) = sfluxref(ig,js) + fs * dsfluxref(ig,js)
            enddo
         endif
      enddo

      end subroutine taumol28

!----------------------------------------------------------------------------
      subroutine taumol29
!----------------------------------------------------------------------------
!
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)
!
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrsw,   only: ng29, ngs28
      use rrsw_kg29, only: absa, absb, forref, selfref, &
                           sfluxref, absh2o, absco2, rayl

      implicit none

! ------- Declarations -------

      integer(kind=im) :: ig, ind0, ind1, js, lay

      real   (kind=rb) :: fac000, fac001, fac010, fac011, fac100, fac101, &
                          fac110, fac111, fs, speccomb, specmult, specparm

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(29) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(29) + 1

         !dir$ vector always
         do ig = 1, ng29
            taug(ngs28+ig,lay) = &
                 fac00h(lay) * absa(ig,ind0) + &
                 fac10h(lay) * absa(ig,ind0+1) + &
                 fac01h(lay) * absa(ig,ind1) + &
                 fac11h(lay) * absa(ig,ind1+1) + &

                 selffac0(lay) * selfref(ig,indself(lay)) + &
                 selffac1(lay) * selfref(ig,indself(lay)+1) + &

                 forfac0(lay) * forref(ig,indfor(lay)) + &
                 forfac1(lay) * forref(ig,indfor(lay)+1) + &

                  colco2(lay) * absco2(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(29) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(29) + 1

         !dir$ vector always
         do ig = 1, ng29
            taug(ngs28+ig,lay) = colco2(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + colh2o(lay) * absh2o(ig)
         enddo
      enddo

      do ig = 1, ng29
         sfluxzen(ngs28+ig) = sfluxref(ig)
      enddo

      end subroutine taumol29

      end subroutine taumol_sw

      end module rrtmg_sw_taumol
