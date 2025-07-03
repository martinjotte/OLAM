!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_taumol.f90,v $
!     author:    $Author: miacono $
!     revision:  $Revision: 1.8 $
!     created:   $Date: 2011/04/08 20:25:01 $
!
      module rrtmg_lw_taumol

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

      private
      public :: taumol

      contains

!----------------------------------------------------------------------------
      subroutine taumol(nlayers, pavel, wx, coldry, &
                        laytrop, jp, jt, jt1, &
                        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                        colbrd, fac00, fac01, fac10, fac11, &
                        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                        rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                        rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                        selffac, selffrac, indself, forfac, forfrac, indfor, &
                        minorfrac, scaleminor, scaleminorn2, indminor, &
                        fracs, taug, n2ovmr, co2vmr)
!----------------------------------------------------------------------------

! *******************************************************************************
! *                                                                             *
! *                  Optical depths developed for the                           *
! *                                                                             *
! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
! *                                                                             *
! *                                                                             *
! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
! *                        131 HARTWELL AVENUE                                  *
! *                        LEXINGTON, MA 02421                                  *
! *                                                                             *
! *                                                                             *
! *                           ELI J. MLAWER                                     *
! *                         JENNIFER DELAMERE                                   *
! *                         STEVEN J. TAUBMAN                                   *
! *                         SHEPARD A. CLOUGH                                   *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                       email:  mlawer@aer.com                                *
! *                       email:  jdelamer@aer.com                              *
! *                                                                             *
! *        The authors wish to acknowledge the contributions of the             *
! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *
! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
! *                                                                             *
! *******************************************************************************
! *                                                                             *
! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
! *                                                                             *
! *******************************************************************************
! *     TAUMOL                                                                  *
! *                                                                             *
! *     This file contains the subroutines TAUGBn (where n goes from            *
! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
! *     per g-value and layer for band n.                                       *
! *                                                                             *
! *  Output:  optical depths (unitless)                                         *
! *           fractions needed to compute Planck functions at every layer       *
! *               and g-value                                                   *
! *                                                                             *
! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
! *     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
! *                                                                             *
! *  Input                                                                      *
! *                                                                             *
! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
! *     COMMON /PRECISE/  ONEMINUS                                              *
! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
! *     COMMON /PROFDATA/ LAYTROP,                                              *
! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
! *    &                  COLO2(MXLAY)
! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
! *    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
! *                                                                             *
! *     Description:                                                            *
! *     NG(IBAND) - number of g-values in band IBAND                            *
! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
! *                   atmospheres that are stored for band IBAND per            *
! *                   pressure level and temperature.  Each of these            *
! *                   atmospheres has different relative amounts of the         *
! *                   key species for the band (i.e. different binary           *
! *                   species parameters).                                      *
! *     NSPB(IBAND) - same for upper atmosphere                                 *
! *     ONEMINUS - since problems are caused in some cases by interpolation     *
! *                parameters equal to or greater than 1, for these cases       *
! *                these parameters are set to this value, slightly < 1.        *
! *     PAVEL - layer pressures (mb)                                            *
! *     TAVEL - layer temperatures (degrees K)                                  *
! *     PZ - level pressures (mb)                                               *
! *     TZ - level temperatures (degrees K)                                     *
! *     LAYTROP - layer at which switch is made from one combination of         *
! *               key species to another                                        *
! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
! *               respectively (molecules/cm**2)                                *
! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
! *                  compute the interpolation factors that multiply the        *
! *                  appropriate reference k-values.  A value of 0 (1) for      *
! *                  i,j indicates that the corresponding factor multiplies     *
! *                  reference k-value for the lower (higher) of the two        *
! *                  appropriate temperatures, and altitudes, respectively.     *
! *     JP - the index of the lower (in altitude) of the two appropriate        *
! *          reference pressure levels needed for interpolation                 *
! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
! *               temperatures needed for interpolation (for pressure           *
! *               levels JP and JP+1, respectively)                             *
! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
! *               (water vapor density)/(atmospheric density at 296K and        *
! *               1013 mb)                                                      *
! *     SELFFRAC - factor needed for temperature interpolation of reference     *
! *                water vapor self-continuum data                              *
! *     INDSELF - index of the lower of the two appropriate reference           *
! *               temperatures needed for the self-continuum interpolation      *
! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
! *     FORFRAC - factor needed for temperature interpolation of reference      *
! *                water vapor foreign-continuum data                           *
! *     INDFOR  - index of the lower of the two appropriate reference           *
! *               temperatures needed for the foreign-continuum interpolation   *
! *                                                                             *
! *  Data input                                                                 *
! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
! *        (note:  n is the band number,'MGAS' is the species name of the minor *
! *         gas)                                                                *
! *                                                                             *
! *     Description:                                                            *
! *     KA - k-values for low reference atmospheres (key-species only)          *
! *          (units: cm**2/molecule)                                            *
! *     KB - k-values for high reference atmospheres (key-species only)         *
! *          (units: cm**2/molecule)                                            *
! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
! *          (units: cm**2/molecule)                                            *
! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
! *          (units: cm**2/molecule)                                            *
! *     SELFREF - k-values for water vapor self-continuum for reference         *
! *               atmospheres (used below LAYTROP)                              *
! *               (units: cm**2/molecule)                                       *
! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
! *               atmospheres (used below/above LAYTROP)                        *
! *               (units: cm**2/molecule)                                       *
! *                                                                             *
! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
! *                                                                             *
!*******************************************************************************

! ------- Declarations -------

      use parkind,  only: im => kind_im, rb => kind_rb
      use rrlw_ref, only: chi_mls4, chi_mls4i
      use parrrtm,  only: maxxsec, ngptlw
      use rrlw_con, only: oneminus
      use rrlw_wvn, only: nspa, nspb

      implicit none

! ----- Input -----

      integer(kind=im), intent(in) :: nlayers         ! total number of layers

      real(kind=rb), intent(in) :: pavel(nlayers)     ! layer pressures (mb)

      real(kind=rb), intent(in) :: wx(nlayers,maxxsec)! cross-section amounts (mol/cm2)

      real(kind=rb), intent(in) :: coldry(nlayers)    ! column amount (dry air)

      integer(kind=im), intent(in) :: laytrop         ! tropopause layer index

      integer(kind=im), intent(in) :: jp(nlayers)

      integer(kind=im), intent(in) :: jt(nlayers)

      integer(kind=im), intent(in) :: jt1(nlayers)

      real(kind=rb), intent(in) :: colh2o(nlayers)    ! column amount (h2o)

      real(kind=rb), intent(in) :: colco2(nlayers)    ! column amount (co2)

      real(kind=rb), intent(in) :: colo3(nlayers)     ! column amount (o3)

      real(kind=rb), intent(in) :: coln2o(nlayers)    ! column amount (n2o)

      real(kind=rb), intent(in) :: colco(nlayers)     ! column amount (co)

      real(kind=rb), intent(in) :: colch4(nlayers)    ! column amount (ch4)

      real(kind=rb), intent(in) :: colo2(nlayers)     ! column amount (o2)

      real(kind=rb), intent(in) :: colbrd(nlayers)    ! column amount (broadening gases)

      integer(kind=im), intent(in) :: indself(nlayers)

      integer(kind=im), intent(in) :: indfor(nlayers)

      real(kind=rb), intent(in) :: selffac(nlayers)

      real(kind=rb), intent(in) :: selffrac(nlayers)

      real(kind=rb), intent(in) :: forfac(nlayers)

      real(kind=rb), intent(in) :: forfrac(nlayers)

      integer(kind=im), intent(in) :: indminor(nlayers)

      real(kind=rb), intent(in) :: minorfrac(nlayers)

      real(kind=rb), intent(in) :: scaleminor(nlayers)

      real(kind=rb), intent(in) :: scaleminorn2(nlayers)

      real(kind=rb), intent(in) :: &
                       fac00(nlayers), fac01(nlayers), &
                       fac10(nlayers), fac11(nlayers)
      real(kind=rb), intent(in) :: &
                       rat_h2oco2(nlayers),rat_h2oco2_1(nlayers), &
                       rat_h2oo3(nlayers),rat_h2oo3_1(nlayers),   &
                       rat_h2on2o(nlayers),rat_h2on2o_1(nlayers), &
                       rat_h2och4(nlayers),rat_h2och4_1(nlayers), &
                       rat_n2oco2(nlayers),rat_n2oco2_1(nlayers), &
                       rat_o3co2(nlayers),rat_o3co2_1(nlayers)

      real(kind=rb), intent(in) :: n2ovmr(nlayers), co2vmr(nlayers)

! ----- Output -----

      real(kind=rb), intent(out) :: fracs(ngptlw,nlayers)        ! planck fractions
      real(kind=rb), intent(out) :: taug(ngptlw,nlayers)         ! gaseous optical depth

! ----- Local -----

      real(kind=rb) :: selffac0(nlayers), selffac1(nlayers)
      real(kind=rb) :: forfac0 (nlayers), forfac1 (nlayers)
      real(kind=rb) :: adjcoln2o(nlayers), adjfac, ratn2o

      integer(kind=im) :: lay

! pre-compute some quantities

      do lay = 1, laytrop
         selffac0(lay) = selffac(lay) * (1.0 - selffrac(lay))
         selffac1(lay) = selffac(lay) * selffrac(lay)
      enddo

      do lay = 1, nlayers
         forfac0(lay) = forfac(lay) * (1.0 - forfrac(lay))
         forfac1(lay) = forfac(lay) * forfrac(lay)

         !  In atmospheres where the amount of N2O is too great to be considered
         !  a minor species, adjust the column amount of N2O by an empirical
         !  factor to obtain the proper contribution.

         ratn2o = n2ovmr(lay) * chi_mls4i(jp(lay)+1)
         if (ratn2o .gt. 1.5_rb) then
            adjfac = 0.5_rb+(ratn2o-0.5_rb)**0.65_rb
            adjcoln2o(lay) = adjfac*chi_mls4(jp(lay)+1)*coldry(lay)*1.e-20_rb
         else
            adjcoln2o(lay) = coln2o(lay)
         endif
      enddo

! Calculate gaseous optical depth and planck fractions for each spectral band.

      call taugb1
      call taugb2
      call taugb3
      call taugb4
      call taugb5
      call taugb6
      call taugb7
      call taugb8
      call taugb9
      call taugb10
      call taugb11
      call taugb12
      call taugb13
      call taugb14
      call taugb15
      call taugb16

      contains

!----------------------------------------------------------------------------
      subroutine taugb1
!----------------------------------------------------------------------------

! ------- Modifications -------
!  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
!
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1:
!           10-250 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng1
      use rrlw_kg01, only : fracrefa, fracrefb, absa, absb, &
                            ka_mn2, kb_mn2, selfref, forref, &
                            dka_mn2, dkb_mn2

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: corradj, scalen2, tauself, taufor, taun2

! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

! Compute the optical depth by interpolating in ln(pressure) and
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(1) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(1) + 1

         if (pavel(lay) .lt. 250._rb) then
            corradj = 1._rb - 0.15_rb * (250._rb-pavel(lay)) / 154.4_rb
         else
            corradj =  1.
         endif

         scalen2 = colbrd(lay) * scaleminorn2(lay)

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng1
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            taun2 = scalen2 * (ka_mn2(ig,indminor(lay)) + minorfrac(lay) * dka_mn2(ig,indminor(lay)))
            taug(ig,lay) = corradj * (colh2o(lay) * &
                (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1)) &
                 + tauself + taufor + taun2)
             fracs(ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(1) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(1) + 1

         corradj =  1._rb - 0.15_rb * (pavel(lay) / 95.6_rb)

         scalen2 = colbrd(lay) * scaleminorn2(lay)

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng1
            taufor = forfac0(lay) * forref(ig,3) + &
                     forfac1(lay) * forref(ig,4)
            taun2 = scalen2 * (kb_mn2(ig,indminor(lay)) + minorfrac(lay) * dkb_mn2(ig,indminor(lay)))
            taug(ig,lay) = corradj * (colh2o(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + taufor + taun2)
            fracs(ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb1

!----------------------------------------------------------------------------
      subroutine taugb2
!----------------------------------------------------------------------------
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2:
!           250 - 500 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng2, ngs1
      use rrlw_kg02, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: corradj, tauself, taufor

! Compute the optical depth by interpolating in ln(pressure) and
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(2) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(2) + 1

         corradj = 1._rb - .05_rb * (pavel(lay) - 100._rb) / 900._rb

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng2
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            taug(ngs1+ig,lay) = corradj * (colh2o(lay) * &
                (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1)) &
                 + tauself + taufor)
            fracs(ngs1+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(2) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(2) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng2
            taufor = forfac0(lay) * forref(ig,3) + &
                     forfac1(lay) * forref(ig,4)
            taug(ngs1+ig,lay) = colh2o(lay) * &
                (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + taufor
            fracs(ngs1+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb2

!----------------------------------------------------------------------------
      subroutine taugb3
!----------------------------------------------------------------------------
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm,   only: ng3, ngs2
      use rrlw_ref,  only: chi_mls1_2
      use rrlw_kg03, only: fracrefa, fracrefb, absa, absb, &
                           ka_mn2o, kb_mn2o, selfref, forref, &
                           dfracrefa, dfracrefb, dka_mn2o, dkb_mn2o

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmn2o, jpl
      integer(kind=im) :: is, is1, imn2o, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
                       fmn2o, fmn2omf
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_planck_b
      real(kind=rb) :: refrat_m_a
      real(kind=rb) :: refrat_m_b

!  P = 212.725 mb
      refrat_planck_a = chi_mls1_2(9)

!  P = 95.58 mb
      refrat_planck_b = chi_mls1_2(13)

!  P = 706.270mb
      refrat_m_a = chi_mls1_2(3)

!  P = 95.58 mb
      refrat_m_b = chi_mls1_2(13)


! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the water vapor
! self-continuum and foreign continuum is interpolated (in temperature)
! separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mn2o = colh2o(lay) + refrat_m_a*colco2(lay)
         specparm_mn2o = min(colh2o(lay)/speccomb_mn2o, oneminus)
         specmult_mn2o = 8._rb*specparm_mn2o
         imn2o = int(specmult_mn2o)
         jmn2o = 1 + imn2o
         fmn2o = specmult_mn2o - imn2o
         fmn2omf = minorfrac(lay)*fmn2o

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(3) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(3) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1_rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng3
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            n2om1 = ka_mn2o(ig,jmn2o,indminor(lay)  ) + fmn2o * dka_mn2o(ig,jmn2o,indminor(lay)  )
            n2om2 = ka_mn2o(ig,jmn2o,indminor(lay)+1) + fmn2o * dka_mn2o(ig,jmn2o,indminor(lay)+1)
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) +  &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs2+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o(lay)*absn2o

            fracs(ngs2+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colh2o(lay) + rat_h2oco2(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 4._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 4._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs1) * fac01(lay)
         fac011 = (1._rb - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_mn2o = colh2o(lay) + refrat_m_b*colco2(lay)
         specparm_mn2o = min(colh2o(lay)/speccomb_mn2o, oneminus)
         specmult_mn2o = 4._rb*specparm_mn2o
         imn2o = int(specmult_mn2o)
         jmn2o = 1 + imn2o
         fmn2o = specmult_mn2o - imn2o
         fmn2omf = minorfrac(lay)*fmn2o

         speccomb_planck = colh2o(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 4._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(3) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(3) + js1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng3
            taufor = forfac0(lay) * forref(ig,3) + &
                     forfac1(lay) * forref(ig,4)

            n2om1 = kb_mn2o(ig,jmn2o,indminor(lay)  ) + fmn2o * dkb_mn2o(ig,jmn2o,indminor(lay))
            n2om2 = kb_mn2o(ig,jmn2o,indminor(lay)+1) + fmn2o * dkb_mn2o(ig,jmn2o,indminor(lay)+1)
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

            taug(ngs2+ig,lay) = speccomb * &
                (fac000 * absb(ig,ind0) + &
                fac100 * absb(ig,ind0+1) + &
                fac010 * absb(ig,ind0+5) + &
                fac110 * absb(ig,ind0+6)) &
                + speccomb1 * &
                (fac001 * absb(ig,ind1) +  &
                fac101 * absb(ig,ind1+1) + &
                fac011 * absb(ig,ind1+5) + &
                fac111 * absb(ig,ind1+6))  &
                + taufor &
                + adjcoln2o(lay)*absn2o

            fracs(ngs2+ig,lay) = fracrefb(ig,jpl) + fpl * dfracrefb(ig,jpl)
         enddo
      enddo

      end subroutine taugb3

!----------------------------------------------------------------------------
      subroutine taugb4
!----------------------------------------------------------------------------
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng4, ngs3
      use rrlw_ref, only : chi_mls1_2, chi_mls3_2
      use rrlw_kg04, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref, dfracrefa, dfracrefb

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jpl
      integer(kind=im) :: is, is1, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor
      real(kind=rb) :: tau_major, tau_major1

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_planck_b

! P =   142.5940 mb
      refrat_planck_a = chi_mls1_2(11)

! P = 95.58350 mb
      refrat_planck_b = chi_mls3_2(13)

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated (in temperature)
! separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(4) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(4) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng4
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) +  &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs3+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor

            fracs(ngs3+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colo3(lay) + rat_o3co2(lay)
         specparm = min(colo3(lay)/speccomb, oneminus)
         specmult = 4._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colo3(lay) + rat_o3co2_1(lay)
         specparm1 = min(colo3(lay)/speccomb1, oneminus)
         specmult1 = 4._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs1) * fac01(lay)
         fac011 = (1._rb - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = min(colo3(lay)/speccomb_planck, oneminus)
         specmult_planck = 4._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl= 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(4) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(4) + js1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng4
            taug(ngs3+ig,lay) =  speccomb * &
                (fac000 * absb(ig,ind0) + &
                fac100 * absb(ig,ind0+1) + &
                fac010 * absb(ig,ind0+5) + &
                fac110 * absb(ig,ind0+6)) &
                + speccomb1 * &
                (fac001 * absb(ig,ind1) +  &
                fac101 * absb(ig,ind1+1) + &
                fac011 * absb(ig,ind1+5) + &
                fac111 * absb(ig,ind1+6))

            fracs(ngs3+ig,lay) = fracrefb(ig,jpl) + fpl * dfracrefb(ig,jpl)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for co2.  Revised to apply weighting for g-point reduction in this band.

!!         taug(ngs3+8,lay)=taug(ngs3+8,lay)*0.92
!!         taug(ngs3+9,lay)=taug(ngs3+9,lay)*0.88
!!         taug(ngs3+10,lay)=taug(ngs3+10,lay)*1.07
!!         taug(ngs3+11,lay)=taug(ngs3+11,lay)*1.1
!!         taug(ngs3+12,lay)=taug(ngs3+12,lay)*0.99
!!         taug(ngs3+13,lay)=taug(ngs3+13,lay)*0.88
!!         taug(ngs3+14,lay)=taug(ngs3+14,lay)*0.943

      enddo

      end subroutine taugb4

!----------------------------------------------------------------------------
      subroutine taugb5
!----------------------------------------------------------------------------
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng5, ngs4
      use rrlw_ref, only : chi_mls1_2, chi_mls3_2
      use rrlw_kg05, only : fracrefa, fracrefb, absa, absb, &
                            ka_mo3, selfref, forref, ccl4, &
                            dfracrefa, dfracrefb

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmo3, jpl
      integer(kind=im) :: is, is1, imo3, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor, o3m1, o3m2, abso3
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_planck_b
      real(kind=rb) :: refrat_m_a

! P = 473.420 mb
      refrat_planck_a = chi_mls1_2(5)

! P = 0.2369 mb
      refrat_planck_b = chi_mls3_2(43)

! P = 317.3480
      refrat_m_a = chi_mls1_2(7)

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the
! water vapor self-continuum and foreign continuum is
! interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mo3 = colh2o(lay) + refrat_m_a*colco2(lay)
         specparm_mo3 = min(colh2o(lay)/speccomb_mo3, oneminus)
         specmult_mo3 = 8._rb*specparm_mo3
         imo3 = int(specmult_mo3)
         jmo3 = 1 + imo3
         fmo3 =specmult_mo3 - imo3

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl =specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(5) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(5) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng5
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            o3m1 = ka_mo3(ig,jmo3,indminor(lay)) + fmo3 * &
                 (ka_mo3(ig,jmo3+1,indminor(lay))-ka_mo3(ig,jmo3,indminor(lay)))
            o3m2 = ka_mo3(ig,jmo3,indminor(lay)+1) + fmo3 * &
                 (ka_mo3(ig,jmo3+1,indminor(lay)+1)-ka_mo3(ig,jmo3,indminor(lay)+1))
            abso3 = o3m1 + minorfrac(lay)*(o3m2-o3m1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs4+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + abso3*colo3(lay) &
                 + wx(lay,1) * ccl4(ig)

            fracs(ngs4+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colo3(lay) + rat_o3co2(lay)
         specparm = min(colo3(lay)/speccomb, oneminus)
         specmult = 4._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colo3(lay) + rat_o3co2_1(lay)
         specparm1 = min(colo3(lay)/speccomb1, oneminus)
         specmult1 = 4._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         fac000 = (1._rb - fs) * fac00(lay)
         fac010 = (1._rb - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._rb - fs1) * fac01(lay)
         fac011 = (1._rb - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = min(colo3(lay)/speccomb_planck, oneminus)
         specmult_planck = 4._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(5) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(5) + js1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng5
            taug(ngs4+ig,lay) = speccomb * &
                (fac000 * absb(ig,ind0) + &
                fac100 * absb(ig,ind0+1) + &
                fac010 * absb(ig,ind0+5) + &
                fac110 * absb(ig,ind0+6)) &
                + speccomb1 * &
                (fac001 * absb(ig,ind1) + &
                fac101 * absb(ig,ind1+1) + &
                fac011 * absb(ig,ind1+5) + &
                fac111 * absb(ig,ind1+6))  &
                + wx(lay,1) * ccl4(ig)

            fracs(ngs4+ig,lay) = fracrefb(ig,jpl) + fpl * dfracrefb(ig,jpl)
         enddo
      enddo

      end subroutine taugb5

!----------------------------------------------------------------------------
      subroutine taugb6
!----------------------------------------------------------------------------
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng6, ngs5
      use rrlw_ref, only : chi_mls2, chi_mls2i
      use rrlw_kg06, only : fracrefa, absa, ka_mco2, dka_mco2, &
                            selfref, forref, cfc11adj, cfc12

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: ratco2, adjfac, adjcolco2
      real(kind=rb) :: tauself, taufor, absco2

! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

! Compute the optical depth by interpolating in ln(pressure) and
! temperature. The water vapor self-continuum and foreign continuum
! is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

! In atmospheres where the amount of CO2 is too great to be considered
! a minor species, adjust the column amount of CO2 by an empirical factor
! to obtain the proper contribution.
         ratco2 = co2vmr(lay)*chi_mls2i(jp(lay)+1)
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 2.0_rb+(ratco2-2.0_rb)**0.77_rb
            adjcolco2 = adjfac*chi_mls2(jp(lay)+1)*coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(6) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(6) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng6
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            absco2 = ka_mco2(ig,indminor(lay)) + minorfrac(lay) * dka_mco2(ig,indminor(lay))
            taug(ngs5+ig,lay) = colh2o(lay) * &
                (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) +  &
                 fac11(lay) * absa(ig,ind1+1))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(lay,2) * cfc11adj(ig) &
                 + wx(lay,3) * cfc12(ig)
            fracs(ngs5+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
! Nothing important goes on above laytrop in this band.
      do lay = laytrop+1, nlayers

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng6
            taug(ngs5+ig,lay) = wx(lay,2) * cfc11adj(ig) &
                              + wx(lay,3) * cfc12(ig)
            fracs(ngs5+ig,lay) = fracrefa(ig)
         enddo
      enddo

      end subroutine taugb6

!----------------------------------------------------------------------------
      subroutine taugb7
!----------------------------------------------------------------------------
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng7, ngs6
      use rrlw_ref, only : chi_mls1_3, chi_mls2, chi_mls2i
      use rrlw_kg07, only : fracrefa, fracrefb, dfracrefa, absa, absb, forref, &
                            ka_mco2, kb_mco2, selfref, dka_mco2, dkb_mco2

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmco2, jpl
      integer(kind=im) :: is, is1, imco2, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor, co2m1, co2m2, absco2
      real(kind=rb) :: ratco2, adjfac, adjcolco2
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_m_a

! P = 706.2620 mb
      refrat_planck_a = chi_mls1_3(3)

! P = 706.2720 mb
      refrat_m_a = chi_mls1_3(3)

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oo3(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oo3_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mco2 = colh2o(lay) + refrat_m_a*colo3(lay)
         specparm_mco2 = min(colh2o(lay)/speccomb_mco2, oneminus)
         specmult_mco2 = 8._rb*specparm_mco2
         imco2 = int(specmult_mco2)
         jmco2 = 1 + imco2
         fmco2 = specmult_mco2 - imco2

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         ratco2 = co2vmr(lay)*chi_mls2i(jp(lay)+1)
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 3.0_rb+(ratco2-3.0_rb)**0.79_rb
            adjcolco2 = adjfac*chi_mls2(jp(lay)+1)*coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         speccomb_planck = colh2o(lay)+refrat_planck_a*colo3(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl= 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(7) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(7) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng7
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            co2m1 = ka_mco2(ig,jmco2,indminor(lay)  ) + fmco2 * dka_mco2(ig,jmco2,indminor(lay))
            co2m2 = ka_mco2(ig,jmco2,indminor(lay)+1) + fmco2 * dka_mco2(ig,jmco2,indminor(lay)+1)
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) +  &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs6+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2

            fracs(ngs6+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         ratco2 = co2vmr(lay)*chi_mls2i(jp(lay)+1)
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 2.0_rb+(ratco2-2.0_rb)**0.79_rb
            adjcolco2 = adjfac*chi_mls2(jp(lay)+1)*coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(7) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(7) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng7
            absco2 = kb_mco2(ig,indminor(lay)) + minorfrac(lay) * dkb_mco2(ig,indminor(lay))
            taug(ngs6+ig,lay) = colo3(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + adjcolco2 * absco2
            fracs(ngs6+ig,lay) = fracrefb(ig)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for o3.  Revised to apply weighting for g-point reduction in this band.

!!         taug(ngs6+6,lay)=taug(ngs6+6,lay)*0.92_rb
!!         taug(ngs6+7,lay)=taug(ngs6+7,lay)*0.88_rb
!!         taug(ngs6+8,lay)=taug(ngs6+8,lay)*1.07_rb
!!         taug(ngs6+9,lay)=taug(ngs6+9,lay)*1.1_rb
!!         taug(ngs6+10,lay)=taug(ngs6+10,lay)*0.99_rb
!!         taug(ngs6+11,lay)=taug(ngs6+11,lay)*0.855_rb

      enddo

      end subroutine taugb7

!----------------------------------------------------------------------------
      subroutine taugb8
!----------------------------------------------------------------------------
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng8, ngs7
      use rrlw_ref, only : chi_mls2, chi_mls2i
      use rrlw_kg08, only : fracrefa, fracrefb, absa, absb, &
                            ka_mco2, ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, &
                            dka_mco2, dka_mn2o, dka_mo3, dkb_mco2, dkb_mn2o, &
                            selfref, forref, cfc12, cfc22adj

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: tauself, taufor, absco2, abso3, absn2o
      real(kind=rb) :: ratco2, adjfac, adjcolco2

! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the water vapor
! self-continuum and foreign continuum is interpolated (in temperature)
! separately.

! Lower atmosphere loop
      do lay = 1, laytrop

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         ratco2 = co2vmr(lay)*chi_mls2i(jp(lay)+1)
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 2.0_rb+(ratco2-2.0_rb)**0.65_rb
            adjcolco2 = adjfac*chi_mls2(jp(lay)+1)*coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(8) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(8) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng8
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            absco2 = ka_mco2(ig,indminor(lay)) + minorfrac(lay) * dka_mco2(ig,indminor(lay))

            abso3 = ka_mo3(ig,indminor(lay)) + minorfrac(lay) * dka_mo3(ig,indminor(lay))

            absn2o = ka_mn2o(ig,indminor(lay)) + minorfrac(lay) * dka_mn2o(ig,indminor(lay))

            taug(ngs7+ig,lay) = colh2o(lay) * &
                 (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) +  &
                 fac11(lay) * absa(ig,ind1+1)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(lay) * abso3 &
                 + coln2o(lay) * absn2o &
                 + wx(lay,3) * cfc12(ig) &
                 + wx(lay,4) * cfc22adj(ig)

            fracs(ngs7+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         ratco2 = co2vmr(lay)*chi_mls2i(jp(lay)+1)
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 2.0_rb+(ratco2-2.0_rb)**0.65_rb
            adjcolco2 = adjfac*chi_mls2(jp(lay)+1) * coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(8) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(8) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng8
            absco2 = kb_mco2(ig,indminor(lay)) + minorfrac(lay) * dkb_mco2(ig,indminor(lay))
            absn2o = kb_mn2o(ig,indminor(lay)) + minorfrac(lay) * dkb_mn2o(ig,indminor(lay))
            taug(ngs7+ig,lay) = colo3(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + adjcolco2*absco2 &
                 + coln2o(lay)*absn2o &
                 + wx(lay,3) * cfc12(ig) &
                 + wx(lay,4) * cfc22adj(ig)
            fracs(ngs7+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb8

!----------------------------------------------------------------------------
      subroutine taugb9
!----------------------------------------------------------------------------
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng9, ngs8
      use rrlw_ref, only : chi_mls1_6, chi_mls4, chi_mls4i
      use rrlw_kg09, only : fracrefa, fracrefb, absa, absb, dfracrefa, &
                            ka_mn2o, kb_mn2o, selfref, forref, &
                            dka_mn2o, dkb_mn2o

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmn2o, jpl
      integer(kind=im) :: is, is1, imn2o, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_m_a

! P = 212 mb
      refrat_planck_a = chi_mls1_6(9)

! P = 706.272 mb
      refrat_m_a = chi_mls1_6(3)

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2och4(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2och4_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mn2o = colh2o(lay) + refrat_m_a*colch4(lay)
         specparm_mn2o = min(colh2o(lay)/speccomb_mn2o, oneminus)
         specmult_mn2o = 8._rb*specparm_mn2o
         imn2o = int(specmult_mn2o)
         jmn2o = 1 + imn2o
         fmn2o = specmult_mn2o - imn2o

         speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(9) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(9) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng9
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            n2om1 = ka_mn2o(ig,jmn2o,indminor(lay)  ) + fmn2o * dka_mn2o(ig,jmn2o,indminor(lay))
            n2om2 = ka_mn2o(ig,jmn2o,indminor(lay)+1) + fmn2o * dka_mn2o(ig,jmn2o,indminor(lay)+1)
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs8+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o(lay)*absn2o

            fracs(ngs8+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(9) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(9) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng9
            absn2o = kb_mn2o(ig,indminor(lay)) + minorfrac(lay) * dkb_mn2o(ig,indminor(lay))

            taug(ngs8+ig,lay) = colch4(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) +  &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + adjcoln2o(lay)*absn2o

            fracs(ngs8+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb9

!----------------------------------------------------------------------------
      subroutine taugb10
!----------------------------------------------------------------------------
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng10, ngs9
      use rrlw_kg10, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: tauself, taufor

! Compute the optical depth by interpolating in ln(pressure) and
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(10) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(10) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng10
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            taug(ngs9+ig,lay) = colh2o(lay) * &
                 (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1))  &
                 + tauself + taufor
            fracs(ngs9+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(10) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(10) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng10
            taufor = forfac0(lay) * forref(ig,3) + &
                     forfac1(lay) * forref(ig,4)
            taug(ngs9+ig,lay) = colh2o(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) +  &
                 fac11(lay) * absb(ig,ind1+1)) &
                 + taufor
            fracs(ngs9+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb10

!----------------------------------------------------------------------------
      subroutine taugb11
!----------------------------------------------------------------------------
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng11, ngs10
      use rrlw_kg11, only : fracrefa, fracrefb, absa, absb, &
                            ka_mo2, kb_mo2, selfref, forref, &
                            dka_mo2, dkb_mo2

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: scaleo2, tauself, taufor, tauo2

! Minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

! Compute the optical depth by interpolating in ln(pressure) and
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(11) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(11) + 1
         scaleo2 = colo2(lay)*scaleminor(lay)

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng11
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            tauo2 =  scaleo2 * (ka_mo2(ig,indminor(lay)) + minorfrac(lay) * dka_mo2(ig,indminor(lay)))
            taug(ngs10+ig,lay) = colh2o(lay) * &
                 (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1)) &
                 + tauself + taufor &
                 + tauo2
            fracs(ngs10+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(11) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(11) + 1
         scaleo2 = colo2(lay)*scaleminor(lay)

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng11
            taufor = forfac0(lay) * forref(ig,3) + &
                     forfac1(lay) * forref(ig,4)
            tauo2 =  scaleo2 * (kb_mo2(ig,indminor(lay)) + minorfrac(lay) * dkb_mo2(ig,indminor(lay)))
            taug(ngs10+ig,lay) = colh2o(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))  &
                 + taufor &
                 + tauo2
            fracs(ngs10+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb11

!----------------------------------------------------------------------------
      subroutine taugb12
!----------------------------------------------------------------------------
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng12, ngs11
      use rrlw_ref, only : chi_mls1_2
      use rrlw_kg12, only : fracrefa, absa, dfracrefa, selfref, forref

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jpl
      integer(kind=im) :: is, is1, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor
      real(kind=rb) :: tau_major, tau_major1

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

      real(kind=rb) :: refrat_planck_a

! P =   174.164 mb
      refrat_planck_a = chi_mls1_2(10)

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum adn foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl= 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(12) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(12) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng12
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs11+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor

            fracs(ngs11+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop

      do lay = laytrop+1, nlayers
         do ig = 1, ng12
            taug(ngs11+ig,lay) = 0.0_rb
            fracs(ngs11+ig,lay) = 0.0_rb
         enddo
      enddo

      end subroutine taugb12

!----------------------------------------------------------------------------
      subroutine taugb13
!----------------------------------------------------------------------------
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng13, ngs12
      use rrlw_ref, only : chi_mls1_4
      use rrlw_kg13, only : fracrefa, fracrefb, absa, dfracrefa, &
                            ka_mco2, ka_mco, kb_mo3, selfref, forref, &
                            dka_mco2, dka_mco, dkb_mo3

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmco2, jmco, jpl
      integer(kind=im) :: is, is1, imco2, imco, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=rb) :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor, co2m1, co2m2, absco2
      real(kind=rb) :: com1, com2, absco, abso3
      real(kind=rb) :: ratco2, adjfac, adjcolco2
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_m_a
      real(kind=rb) :: refrat_m_a3

! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls1_4(5)

! P = 1053. (Level 1)
      refrat_m_a = chi_mls1_4(1)

! P = 706. (Level 3)
      refrat_m_a3 = chi_mls1_4(3)

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2on2o(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2on2o_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mco2 = colh2o(lay) + refrat_m_a*coln2o(lay)
         specparm_mco2 = min(colh2o(lay)/speccomb_mco2, oneminus)
         specmult_mco2 = 8._rb*specparm_mco2
         imco2 = int(specmult_mco2)
         jmco2 = 1 + imco2
         fmco2 = specmult_mco2 - imco2

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         ratco2 = co2vmr(lay)/3.55e-4_rb
         if (ratco2 .gt. 3.0_rb) then
            adjfac = 2.0_rb+(ratco2-2.0_rb)**0.68_rb
            adjcolco2 = adjfac*3.55e-4*coldry(lay)*1.e-20_rb
         else
            adjcolco2 = colco2(lay)
         endif

         speccomb_mco = colh2o(lay) + refrat_m_a3*coln2o(lay)
         specparm_mco = min(colh2o(lay)/speccomb_mco, oneminus)
         specmult_mco = 8._rb*specparm_mco
         imco = int(specmult_mco)
         jmco = 1 + imco
         fmco = specmult_mco - imco

         speccomb_planck = colh2o(lay)+refrat_planck_a*coln2o(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(13) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(13) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng13
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            co2m1 = ka_mco2(ig,jmco2,indminor(lay)  ) + fmco2 * dka_mco2(ig,jmco2,indminor(lay))
            co2m2 = ka_mco2(ig,jmco2,indminor(lay)+1) + fmco2 * dka_mco2(ig,jmco2,indminor(lay)+1)
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)

            com1 = ka_mco(ig,jmco,indminor(lay)  ) + fmco * dka_mco(ig,jmco,indminor(lay))
            com2 = ka_mco(ig,jmco,indminor(lay)+1) + fmco * dka_mco(ig,jmco,indminor(lay)+1)
            absco = com1 + minorfrac(lay) * (com2 - com1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs12+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(lay)*absco

            fracs(ngs12+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng13
            abso3 = kb_mo3(ig,indminor(lay)) + minorfrac(lay) * dkb_mo3(ig,indminor(lay))
            taug(ngs12+ig,lay) = colo3(lay)*abso3
            fracs(ngs12+ig,lay) =  fracrefb(ig)
         enddo
      enddo

      end subroutine taugb13

!----------------------------------------------------------------------------
      subroutine taugb14
!----------------------------------------------------------------------------
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng14, ngs13
      use rrlw_kg14, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      real(kind=rb) :: tauself, taufor

! Compute the optical depth by interpolating in ln(pressure) and
! temperature.  Below laytrop, the water vapor self-continuum
! and foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(14) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(14) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng14
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)
            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)
            taug(ngs13+ig,lay) = colco2(lay) * &
                 (fac00(lay) * absa(ig,ind0) + &
                 fac10(lay) * absa(ig,ind0+1) + &
                 fac01(lay) * absa(ig,ind1) + &
                 fac11(lay) * absa(ig,ind1+1)) &
                 + tauself + taufor
            fracs(ngs13+ig,lay) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(14) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(14) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng14
            taug(ngs13+ig,lay) = colco2(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))
            fracs(ngs13+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb14

!----------------------------------------------------------------------------
      subroutine taugb15
!----------------------------------------------------------------------------
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng15, ngs14
      use rrlw_ref, only : chi_mls4_2
      use rrlw_kg15, only : fracrefa, absa, ka_mn2, selfref, forref, &
                            dfracrefa, dka_mn2

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jmn2, jpl
      integer(kind=im) :: is, is1, imn2, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: scalen2, tauself, taufor, n2m1, n2m2, taun2
      real(kind=rb) :: tau_major, tau_major1

! Minor gas mapping level :
!     Lower - Nitrogen Continuum, P = 1053., T = 294.

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

      real(kind=rb) :: refrat_planck_a
      real(kind=rb) :: refrat_m_a

! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls4_2(1)

! P = 1053.
      refrat_m_a = chi_mls4_2(1)

! Compute the optical depth by interpolating in ln(pressure),
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = coln2o(lay) + rat_n2oco2(lay)
         specparm = min(coln2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = coln2o(lay) + rat_n2oco2_1(lay)
         specparm1 = min(coln2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_mn2 = coln2o(lay) + refrat_m_a*colco2(lay)
         specparm_mn2 = min(coln2o(lay)/speccomb_mn2, oneminus)
         specmult_mn2 = 8._rb*specparm_mn2
         imn2 = int(specmult_mn2)
         jmn2 = 1 + imn2
         fmn2 = specmult_mn2 - imn2

         speccomb_planck = coln2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = min(coln2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl= 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(15) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(15) + js1

         scalen2 = colbrd(lay)*scaleminor(lay)

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng15
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            n2m1 = ka_mn2(ig,jmn2,indminor(lay)  ) + fmn2 * dka_mn2(ig,jmn2,indminor(lay))
            n2m2 = ka_mn2(ig,jmn2,indminor(lay)+1) + fmn2 * dka_mn2(ig,jmn2,indminor(lay)+1)
            taun2 = scalen2 * (n2m1 + minorfrac(lay) * (n2m2 - n2m1))

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs14+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + taun2

            fracs(ngs14+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop

      do lay = laytrop+1, nlayers
         do ig = 1, ng15
            taug (ngs14+ig,lay) = 0.0_rb
            fracs(ngs14+ig,lay) = 0.0_rb
         enddo
      enddo

      end subroutine taugb15

!----------------------------------------------------------------------------
      subroutine taugb16
!----------------------------------------------------------------------------
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!----------------------------------------------------------------------------

! ------- Modules -------

      use parrrtm, only : ng16, ngs15
      use rrlw_ref, only : chi_mls1_6
      use rrlw_kg16, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref, dfracrefa

! ------- Declarations -------

! Local
      integer(kind=im) :: lay, ind0, ind1, inds, indf, ig
      integer(kind=im) :: js, js1, jpl
      integer(kind=im) :: is, is1, ipl
      real(kind=rb) :: speccomb, specparm, specmult, fs
      real(kind=rb) :: speccomb1, specparm1, specmult1, fs1
      real(kind=rb) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=rb) :: p, p4, fk0, fk1, fk2
      real(kind=rb) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=rb) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=rb) :: tauself, taufor
      real(kind=rb) :: tau_major, tau_major1

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

      real(kind=rb) :: refrat_planck_a

! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls1_6(6)

! Compute the optical depth by interpolating in ln(pressure),
! temperature,and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated
! (in temperature) separately.

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2och4(lay)
         specparm = min(colh2o(lay)/speccomb, oneminus)
         specmult = 8._rb*specparm
         is = int(specmult)
         js = 1 + is
         fs = specmult - is

         speccomb1 = colh2o(lay) + rat_h2och4_1(lay)
         specparm1 = min(colh2o(lay)/speccomb1, oneminus)
         specmult1 = 8._rb*specparm1
         is1 = int(specmult1)
         js1 = 1 + is1
         fs1 = specmult1 - is1

         speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
         specparm_planck = min(colh2o(lay)/speccomb_planck, oneminus)
         specmult_planck = 8._rb*specparm_planck
         ipl = int(specmult_planck)
         jpl = 1 + ipl
         fpl = specmult_planck - ipl

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(16) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(16) + js1

         if (specparm .lt. 0.125_rb) then
            p = fs - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_rb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._rb - fs) * fac00(lay)
            fac010 = (1._rb - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_rb) then
            p = fs1 - 1._rb
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_rb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._rb - p - 2.0_rb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._rb - fs1) * fac01(lay)
            fac011 = (1._rb - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng16
            tauself = selffac0(lay) * selfref(ig,indself(lay)) + &
                      selffac1(lay) * selfref(ig,indself(lay)+1)

            taufor = forfac0(lay) * forref(ig,indfor(lay)) + &
                     forfac1(lay) * forref(ig,indfor(lay)+1)

            if (specparm .lt. 0.125_rb) then
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac200 * absa(ig,ind0+2) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10) + &
                    fac210 * absa(ig,ind0+11))
            else if (specparm .gt. 0.875_rb) then
               tau_major = speccomb * &
                    (fac200 * absa(ig,ind0-1) + &
                    fac100 * absa(ig,ind0) + &
                    fac000 * absa(ig,ind0+1) + &
                    fac210 * absa(ig,ind0+8) + &
                    fac110 * absa(ig,ind0+9) + &
                    fac010 * absa(ig,ind0+10))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ig,ind0) + &
                    fac100 * absa(ig,ind0+1) + &
                    fac010 * absa(ig,ind0+9) + &
                    fac110 * absa(ig,ind0+10))
            endif

            if (specparm1 .lt. 0.125_rb) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac201 * absa(ig,ind1+2) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10) + &
                    fac211 * absa(ig,ind1+11))
            else if (specparm1 .gt. 0.875_rb) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ig,ind1-1) + &
                    fac101 * absa(ig,ind1) + &
                    fac001 * absa(ig,ind1+1) + &
                    fac211 * absa(ig,ind1+8) + &
                    fac111 * absa(ig,ind1+9) + &
                    fac011 * absa(ig,ind1+10))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ig,ind1) + &
                    fac101 * absa(ig,ind1+1) + &
                    fac011 * absa(ig,ind1+9) + &
                    fac111 * absa(ig,ind1+10))
            endif

            taug(ngs15+ig,lay) = tau_major + tau_major1 &
                 + tauself + taufor

            fracs(ngs15+ig,lay) = fracrefa(ig,jpl) + fpl * dfracrefa(ig,jpl)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(16) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(16) + 1

         !dir$ vector always
         !dir$ ivdep
         do ig = 1, ng16
            taug(ngs15+ig,lay) = colch4(lay) * &
                 (fac00(lay) * absb(ig,ind0) + &
                 fac10(lay) * absb(ig,ind0+1) + &
                 fac01(lay) * absb(ig,ind1) + &
                 fac11(lay) * absb(ig,ind1+1))
            fracs(ngs15+ig,lay) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb16

      end subroutine taumol

      end module rrtmg_lw_taumol
