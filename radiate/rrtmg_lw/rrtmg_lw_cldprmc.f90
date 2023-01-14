!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_cldprmc.f90,v $
!     author:    $Author: miacono $
!     revision:  $Revision: 1.9 $
!     created:   $Date: 2011/04/08 20:25:00 $
!
      module rrtmg_lw_cldprmc

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
      public :: cldprmc

      contains

! ------------------------------------------------------------------------------
      subroutine cldprmc(nlayers, inflag, iceflag, liqflag, cldf, &
                         ciwpmc, clwpmc, reicmc, relqmc, taucmc)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

      use parkind,  only: im => kind_im, rb => kind_rb, cldmin
      use parrrtm,  only: ngptlw, nbndlw
      use rrlw_cld, only: absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_wvn, only: ngb

      implicit none

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: inflag          ! see definitions
      integer(kind=im), intent(in) :: iceflag         ! see definitions
      integer(kind=im), intent(in) :: liqflag         ! see definitions

      real(kind=rb), intent(in) :: cldf(nlayers)      ! mean cloud fraction

      real(kind=rb), intent(in) :: ciwpmc(ngptlw,nlayers)   ! cloud ice water path [mcica]
      real(kind=rb), intent(in) :: clwpmc(ngptlw,nlayers)   ! cloud liquid water path [mcica]

      real(kind=rb), intent(in) :: relqmc(nlayers)    ! liquid particle effective radius (microns)
      real(kind=rb), intent(in) :: reicmc(nlayers)    ! ice particle effective radius (microns)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
! ------- Output -------

      real(kind=rb), intent(inout) :: taucmc(ngptlw,nlayers)     ! cloud optical depth [mcica]

! ------- Local -------

      integer(kind=im) :: lay                         ! Layer index
      integer(kind=im) :: ib                          ! spectral band index
      integer(kind=im) :: ig                          ! g-point interval index
      integer(kind=im) :: index

      real(kind=rb) :: abscoice(ngptlw)               ! ice absorption coefficients
      real(kind=rb) :: abscoliq(ngptlw)               ! liquid absorption coefficients
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)
      real(kind=rb) :: factor                         !
      real(kind=rb) :: fint                           !
      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical
!                  depth is computed as in CCM2.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CCM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent
!                     as possible:
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input
!                     and the optical depths due to water clouds are computed
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a
!                     range of effective radii by an averaging procedure
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption
!                     coefficients for the input effective radius.

      integer(kind=im), parameter :: icb(nbndlw) = (/ 1,2,3,3,3,4,4,4,5,5,5,5,5,5,5,5 /)

      if (inflag /= 2 ) return

      ! Main layer loop
      do lay = 1, nlayers

         if (cldf(lay) >= cldmin) then

            ! Calculation of absorption coefficients due to ice clouds.
            if (reicmc(lay) > 0.01) then

               if (iceflag .le. 0) then

                  radice = max(10.0_rb, reicmc(lay))
                  abscoice(1:ngptlw) = absice0(1) + absice0(2) / radice

               elseif (iceflag .eq. 1) then

                  radice = max(13.0_rb, min(130._rb, reicmc(lay)))

                  do ig = 1, ngptlw
                     ib     = icb(ngb(ig))
                     abscoice(ig) = absice1(1,ib) + absice1(2,ib) / radice
                  enddo

               elseif (iceflag .eq. 2) then

                  radice = max(5.0_rb, min(131._rb, reicmc(lay)))
                  factor = (radice - 2._rb)/3._rb
                  index  = min(int(factor), 42)
                  fint   = factor - real(index)

                  do ig = 1, ngptlw
                     ib     = ngb(ig)
                     abscoice(ig) = absice2(index,ib) + fint * &
                          (absice2(index+1,ib) - (absice2(index,ib)))
                  enddo

               elseif (iceflag .ge. 3) then

                  radice = max(5.0_rb, min(140._rb, reicmc(lay)))
                  factor = (radice - 2._rb)/3._rb
                  index  = min(int(factor), 45)
                  fint   = factor - real(index)

                  do ig = 1, ngptlw
                     ib     = ngb(ig)
                     abscoice(ig) = absice3(index,ib) + fint * &
                          (absice3(index+1,ib) - (absice3(index,ib)))
                  enddo

               endif

            else

               abscoice(1:ngptlw) = 0.0

            endif

            ! Calculation of absorption coefficients due to water clouds.

            if (relqmc(lay) > 0.01) then

               if (liqflag .le. 0) then

                  abscoliq(1:ngptlw) = absliq0

               elseif (liqflag .ge. 1) then

                  radliq = max(2.5_rb, min(60._rb, relqmc(lay))) - 1.5_rb
                  index  = max(1, min(57, int(radliq)))
                  fint = radliq - real(index)

                  do ig = 1, ngptlw
                     ib = ngb(ig)
                     abscoliq(ig) = absliq1(index,ib) + fint * &
                          (absliq1(index+1,ib) - (absliq1(index,ib)))
                  enddo

               endif

            else

               abscoliq(1:ngptlw) = 0.0

            endif

            do ig = 1, ngptlw

               taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) &
                              + clwpmc(ig,lay) * abscoliq(ig)
            enddo

         endif
      enddo

      end subroutine cldprmc

      end module rrtmg_lw_cldprmc
