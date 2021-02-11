!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $

      module rrtmg_sw_cldprmc

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
      public :: cldprmc_sw

      contains

! ----------------------------------------------------------------------------
      subroutine cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldf, &
                            ciwpmc, clwpmc, reicmc, relqmc, &
                            taucmc, ssacmc, asmcmc)
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer
! and g-point interval for use by the McICA method.
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=1,2,3 are available;
! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented.

      use parkind,  only : im => kind_im, rb => kind_rb, cldmin
      use parrrsw,  only : ngptsw, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2, ngbo

      implicit none

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers             ! total number of layers
      integer(kind=im), intent(in) :: inflag              ! see definitions
      integer(kind=im), intent(in) :: iceflag             ! see definitions
      integer(kind=im), intent(in) :: liqflag             ! see definitions

      real(kind=rb), intent(in) :: cldf(nlayers)          ! mean cloud fraction

      real(kind=rb), intent(in) :: ciwpmc(ngptsw,nlayers) ! cloud ice water path [mcica]

      real(kind=rb), intent(in) :: clwpmc(ngptsw,nlayers) ! cloud liquid water path [mcica]

      real(kind=rb), intent(in) :: relqmc(nlayers)        ! cloud liquid particle effective radius (microns)

      real(kind=rb), intent(in) :: reicmc(nlayers)        ! cloud ice particle effective radius (microns)

                                            ! specific definition of reicmc depends on setting of iceflag:
                                            ! iceflag = 0: (inactive)
                                            !
                                            ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                            !              r_ec range is limited to 13.0 to 130.0 microns
                                            ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                            !              r_k range is limited to 5.0 to 131.0 microns
                                            ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                            !              dge range is limited to 5.0 to 140.0 microns
                                            !              [dge = 1.0315 * r_ec]
! ------- Output -------

      real(kind=rb), intent(out) :: taucmc(ngptsw,nlayers) ! cloud optical depth (non-delta scaled)

      real(kind=rb), intent(out) :: ssacmc(ngptsw,nlayers) ! single scattering albedo (non-delta scaled)

      real(kind=rb), intent(out) :: asmcmc(ngptsw,nlayers) ! asymmetry parameter (non-delta scaled)

! ------- Local -------

      integer(kind=im) :: ib, lay, index, icx, ig

      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)
      real(kind=rb) :: factor
      real(kind=rb) :: fint
      real(kind=rb) :: scatice, tauice, scatliq, tauliq
      real(kind=rb) :: extcoice(ngptsw), gice(ngptsw), ssacoice(ngptsw)
      real(kind=rb) :: extcoliq(ngptsw), gliq(ngptsw), ssacoliq(ngptsw)

      if (inflag /= 2) return

      ! Main layer loop
      do lay = 1, nlayers

         if (cldf(lay) >= cldmin) then

            ! Calculation of absorption coefficients due to ice clouds.
            if (reicmc(lay) > 0.01) then

! (iceflag = 1):
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat unjustified for large ice particles
               if (iceflag .le. 1) then

                  radice = max(13._rb, min(130.0_rb, reicmc(lay)))

                  do ig = 1, ngptsw
                     ib = ngbo(ig)
                     if (wavenum2(ib) .gt. 1.43e04_rb) then
                        icx = 1
                     elseif (wavenum2(ib) .gt. 7.7e03_rb) then
                        icx = 2
                     elseif (wavenum2(ib) .gt. 5.3e03_rb) then
                        icx = 3
                     elseif (wavenum2(ib) .gt. 4.0e03_rb) then
                        icx = 4
                     elseif (wavenum2(ib) .ge. 2.5e03_rb) then
                        icx = 5
                     endif
                     extcoice(ig) = max(0., abari(icx) + bbari(icx)/radice)
                     ssacoice(ig) = max(0., min(1., 1._rb - cbari(icx) - dbari(icx) * radice))
                     gice    (ig) = max(0., min(1., ebari(icx) + fbari(icx) * radice))
                  enddo

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then

                  radice = max(5.0_rb, min(131._rb, reicmc(lay)))
                  factor = (radice - 2._rb)/3._rb
                  index = min(42,int(factor))
                  fint = factor - real(index,kind=rb)

                  do ig = 1, ngptsw
                     ib = ngbo(ig)
                     extcoice (ig) = max(0., extice2(index,ib) + fint * &
                                     (extice2(index+1,ib) -  extice2(index,ib)))
                     ssacoice (ig) = max(0., min(1., ssaice2(index,ib) + fint * &
                                     (ssaice2(index+1,ib) -  ssaice2(index,ib))))
                     gice     (ig) = max(0., min(1., asyice2(index,ib) + fint * &
                                     (asyice2(index+1,ib) -  asyice2(index,ib))))
                  enddo

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               else

                  radice = max(5.0_rb, min(140._rb, reicmc(lay)))
                  factor = (radice - 2._rb)/3._rb
                  index = min(45,int(factor))
                  fint = factor - real(index,kind=rb)

                  do ig = 1, ngptsw
                     ib = ngbo(ig)
                     extcoice(ig) = max(0., extice3(index,ib) + fint * &
                                    (extice3(index+1,ib) - extice3(index,ib)))
                     ssacoice(ig) = max(0., min(1., ssaice3(index,ib) + fint * &
                                    (ssaice3(index+1,ib) - ssaice3(index,ib))))
                     gice    (ig) = max(0., min(1., asyice3(index,ib) + fint * &
                                    (asyice3(index+1,ib) - asyice3(index,ib))))
                  enddo

               endif

            else

               do ig = 1, ngptsw
                  extcoice(ig) = 0.0
                  ssacoice(ig) = 0.0
                  gice    (ig) = 0.0
               enddo

            endif

            ! Calculation of absorption coefficients due to water clouds.
            if (reicmc(lay) > 0.01) then

               radliq = max(2.5_rb, min(60._rb, relqmc(lay))) - 1.5_rb
               index = max(1, min(57, int(radliq)))
               fint = radliq - real(index,kind=rb)

               do ig = 1, ngptsw
                  ib = ngbo(ig)
                  extcoliq(ig) = max(0., extliq1(index,ib) + fint * &
                                 (extliq1(index+1,ib) - extliq1(index,ib)))
                  ssacoliq(ig) = max(0.0, min(1.0, ssaliq1(index,ib) + fint * &
                                 (ssaliq1(index+1,ib) - ssaliq1(index,ib))))
                  gliq    (ig) = max(0., min(1.0, asyliq1(index,ib) + fint * &
                                 (asyliq1(index+1,ib) - asyliq1(index,ib))))
               enddo

            else

               do ig = 1, ngptsw
                  extcoliq(ig) = 0.0_rb
                  ssacoliq(ig) = 0.0_rb
                  gliq    (ig) = 0.0_rb
               enddo

            endif

            do ig = 1, ngptsw
               tauliq = clwpmc(ig,lay) * extcoliq(ig)
               tauice = ciwpmc(ig,lay) * extcoice(ig)

               scatliq = ssacoliq(ig) * tauliq
               scatice = ssacoice(ig) * tauice

               taucmc(ig,lay) = tauliq + tauice
               ssacmc(ig,lay) = scatliq + scatice
               asmcmc(ig,lay) = scatliq * gliq(ig) + scatice * gice(ig)
            enddo

         endif

         ! End layer loop
      enddo

      end subroutine cldprmc_sw

      end module rrtmg_sw_cldprmc
