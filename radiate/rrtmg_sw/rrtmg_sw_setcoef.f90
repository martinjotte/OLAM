!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $

      module rrtmg_sw_setcoef

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
      public :: setcoef_sw

      contains

!----------------------------------------------------------------------------
      subroutine setcoef_sw(nlayers, pavel, tavel, coldry, &
                            h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, &
                            laytrop, jp, jt, jt1, &
                            colch4, colco2, colh2o, colmol, &
                            colo2, colo3, fac00, fac01, fac10, fac11, &
                            selffac, selffrac, indself, forfac, forfrac, indfor)
!----------------------------------------------------------------------------
!
! Purpose:  For a given atmosphere, calculate the indices and
! fractions related to the pressure and temperature interpolations.

! Modifications:
! Original: J. Delamere, AER, Inc. (version 2.5, 02/04/01)
! Revised: Rewritten and adapted to ECMWF F90, JJMorcrette 030224
! Revised: For uniform rrtmg formatting, MJIacono, Jul 2006

! ------ Declarations -------

      use parkind,  only: im => kind_im, rb => kind_rb
      use rrsw_ref, only: preflog, tref

      implicit none

! ----- Input -----
      integer(kind=im), intent(in) :: nlayers         ! total number of layers

      real(kind=rb), intent(in) :: pavel(nlayers)     ! layer pressures (mb)

      real(kind=rb), intent(in) :: tavel(nlayers)     ! layer temperatures (K)

      real(kind=rb), intent(in) :: coldry(nlayers)    ! dry air column density (mol/cm2)

      real(kind=rb), intent(in) :: h2ovmr(nlayers)    ! H2O volume mixing ratio

      real(kind=rb), intent(in) :: o3vmr(nlayers)     ! O3 volume mixing ratio

      real(kind=rb), intent(in) :: co2vmr(nlayers)    ! CO2 volume mixing ratio

      real(kind=rb), intent(in) :: ch4vmr(nlayers)    ! Methane volume mixing ratio

      real(kind=rb), intent(in) :: o2vmr(nlayers)     ! Oxygen volume mixing ratio

! ----- Output -----
      integer(kind=im), intent(out) :: laytrop        ! tropopause layer index

      integer(kind=im), intent(out) :: jp(nlayers)

      integer(kind=im), intent(out) :: jt(nlayers)

      integer(kind=im), intent(out) :: jt1(nlayers)

      real(kind=rb), intent(out) :: colh2o(nlayers)   ! column amount (h2o)

      real(kind=rb), intent(out) :: colco2(nlayers)   ! column amount (co2)

      real(kind=rb), intent(out) :: colo3(nlayers)    ! column amount (o3)

      real(kind=rb), intent(out) :: colch4(nlayers)   ! column amount (ch4)

      real(kind=rb), intent(out) :: colo2(nlayers)    ! column amount (o2)

      real(kind=rb), intent(out) :: colmol(nlayers)

      integer(kind=im), intent(out) :: indself(nlayers)

      integer(kind=im), intent(out) :: indfor(nlayers)

      real(kind=rb), intent(out) :: selffac(nlayers)

      real(kind=rb), intent(out) :: selffrac(nlayers)

      real(kind=rb), intent(out) :: forfac(nlayers)

      real(kind=rb), intent(out) :: forfrac(nlayers)

      real(kind=rb), intent(out) :: fac00(nlayers), fac01(nlayers), &
                                    fac10(nlayers), fac11(nlayers)

! ----- Local -----

      integer(kind=im) :: lay
      integer(kind=im) :: jp1

      real(kind=rb) :: plog
      real(kind=rb) :: fp
      real(kind=rb) :: tt
      real(kind=rb) :: ft
      real(kind=rb) :: ft1
      real(kind=rb) :: scalefac
      real(kind=rb) :: factor
      real(kind=rb) :: compfp
      real(kind=rb) :: cdf

      real(kind=rb), parameter :: stpfac = 296._rb/1013._rb

      do lay = 1, nlayers

! Find the two reference pressures on either side of the
! layer pressure.  Store them in JP and JP1.  Store in FP the
! fraction of the difference (in ln(pressure)) between these
! two values that the layer pressure lies.

         plog = log(pavel(lay))
         jp(lay) = min(58,max(1,int(36._rb - 5._rb*(plog+0.04_rb))))
         fp = 5._rb * (preflog(jp(lay)) - plog)
         jp1 = jp(lay) + 1

! Determine, for each reference pressure (JP and JP1), which
! reference temperature (these are different for each
! reference pressure) is nearest the layer temperature but does
! not exceed it.  Store these indices in JT and JT1, resp.
! Store in FT (resp. FT1) the fraction of the way between JT
! (JT1) and the next highest reference temperature that the
! layer temperature falls.

         tt      = (tavel(lay) - tref(jp(lay))) / 15._rb
         jt(lay) = min(4,max(1,int(3._rb + tt)))
         ft      = tt - real((jt(lay)-3),kind=rb)

         tt       = (tavel(lay) - tref(jp1)) / 15._rb
         jt1(lay) = min(4,max(1,int(3._rb + tt)))
         ft1      = tt - real((jt1(lay)-3),kind=rb)

         scalefac    = pavel(lay) * stpfac / tavel(lay)
         forfac(lay) = scalefac / (1._rb + h2ovmr(lay))

         cdf = 1.e-20_rb * coldry(lay)

         colh2o(lay) = cdf * h2ovmr(lay)
         colco2(lay) = cdf * co2vmr(lay)
         colo3 (lay) = cdf * o3vmr(lay)
         colch4(lay) = cdf * ch4vmr(lay)
         colo2 (lay) = cdf * o2vmr(lay)
         colmol(lay) = cdf + colh2o(lay)

! We have now isolated the layer ln pressure and temperature,
! between two reference pressures and two reference temperatures
! (for each reference pressure).  We multiply the pressure
! fraction FP with the appropriate temperature fractions to get
! the factors that will be needed for the interpolation that yields
! the optical depths (performed in routines TAUGBn for band n).

         compfp = 1._rb - fp
         fac10(lay) = compfp * ft
         fac00(lay) = compfp * (1._rb - ft)
         fac11(lay) = fp * ft1
         fac01(lay) = fp * (1._rb - ft1)

      enddo

!  Find tropopause level

      laytrop = 1
      do lay = 1, nlayers
         if (jp(lay) >= 13) exit
      enddo
      laytrop = min(lay-1, nlayers)

!  If the pressure is more than ~100mb, perform a different
!  set of species interpolations.

      do lay = 1, laytrop

! Set up factors needed to separately include the water vapor
! foreign-continuum in the calculation of absorption coefficient.

         factor = (332.0_rb-tavel(lay))/36.0_rb
         indfor(lay) = min(2, max(1, int(factor)))
         forfrac(lay) = factor - real(indfor(lay),kind=rb)

! Set up factors needed to separately include the water vapor
! self-continuum in the calculation of absorption coefficient.

         selffac(lay) = h2ovmr(lay) * forfac(lay)
         factor = (tavel(lay)-188.0_rb)/7.2_rb
         indself(lay) = min(9, max(1, int(factor)-7))
         selffrac(lay) = factor - real((indself(lay) + 7),kind=rb)

      enddo

 !  Above laytrop.

      do lay = laytrop + 1, nlayers

! Set up factors needed to separately include the water vapor
! foreign-continuum in the calculation of absorption coefficient.

         factor = (tavel(lay)-188.0_rb)/36.0_rb
         forfrac(lay) = factor - 1.0_rb

      enddo

      end subroutine setcoef_sw

      end module rrtmg_sw_setcoef
