      module rrlw_wvn

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : nbndlw, mg, ngptlw, maxinpx, maxxsec

      implicit none

      private :: im, rb, nbndlw, mg, ngptlw, maxinpx, maxxsec, i

!------------------------------------------------------------------
! rrtmg_lw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: For the lower atmosphere, the number of reference
!                    atmospheres that are stored for each spectral band
!                    per pressure level and temperature.  Each of these
!                    atmospheres has different relative amounts of the
!                    key species for the band (i.e. different binary
!                    species parameters).
! nspb   :  integer: Same as nspa for the upper atmosphere
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
! totplnk:  real   : Integrated Planck value for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if
!                    individual band output requested
!totplnkderiv: real: Integrated Planck function derivative with respect
!                    to temperature for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16deriv:real: Integrated Planck function derivative with respect
!                    to temperature for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if
!                    individual band output requested
!
! ngc    :  integer: The number of new g-intervals in each band
! ngs    :  integer: The cumulative sum of new g-intervals for each band
! ngm    :  integer: The index of each new g-interval relative to the
!                    original 16 g-intervals in each band
! ngn    :  integer: The number of original g-intervals that are
!                    combined to make each new g-intervals in each band
! ngb    :  integer: The band index for each new g-interval
! wt     :  real   : RRTM weights for the original 16 g-intervals
! rwgt   :  real   : Weights for combining original 16 g-intervals
!                    (256 total) into reduced set of g-intervals
!                    (140 total)
! nxmol  :  integer: Number of cross-section molecules
! ixindx :  integer: Flag for active cross-sections in calculation
!------------------------------------------------------------------

      integer(kind=im), parameter :: ng(nbndlw) = &
           (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)

      integer(kind=im), parameter :: nspa(nbndlw) = &
           (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)

      integer(kind=im), parameter :: nspb(nbndlw) = &
           (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

      real(kind=rb), parameter :: wavenum1(nbndlw) = &
           (/ 10._rb, 350._rb, 500._rb, 630._rb, 700._rb, 820._rb, &
             980._rb,1080._rb,1180._rb,1390._rb,1480._rb,1800._rb, &
             2080._rb,2250._rb,2380._rb,2600._rb/)

      real(kind=rb), parameter :: wavenum2(nbndlw) = &
           (/350._rb, 500._rb, 630._rb, 700._rb, 820._rb, 980._rb, &
            1080._rb,1180._rb,1390._rb,1480._rb,1800._rb,2080._rb, &
            2250._rb,2380._rb,2600._rb,3250._rb/)

      real(kind=rb), parameter :: delwave(nbndlw) = &
           (/ 340._rb, 150._rb, 130._rb,  70._rb, 120._rb, 160._rb, &
              100._rb, 100._rb, 210._rb,  90._rb, 320._rb, 280._rb, &
              170._rb, 130._rb, 220._rb, 650._rb/)

      integer :: i

      real(kind=rb), parameter :: delwaveg(ngptlw) = &
           (/ (340._rb, i=1,10), &
              (150._rb, i=1,12), &
              (130._rb, i=1,16), &
              ( 70._rb, i=1,14), &
              (120._rb, i=1,16), &
              (160._rb, i=1, 8), &
              (100._rb, i=1,12), &
              (100._rb, i=1, 8), &
              (210._rb, i=1,12), &
              ( 90._rb, i=1, 6), &
              (320._rb, i=1, 8), &
              (280._rb, i=1, 8), &
              (170._rb, i=1, 4), &
              (130._rb, i=1, 2), &
              (220._rb, i=1, 2), &
              (650._rb, i=1, 2) /)

      real(kind=rb), allocatable ::  totplnk(:,:)
      real(kind=rb), allocatable :: dtotplnk(:,:)

      real(kind=rb), allocatable ::  totplnkderiv(:,:)
      real(kind=rb), allocatable :: dtotplnkderiv(:,:)

      integer(kind=im), parameter :: ngc(nbndlw) = &
           (/10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)

      integer(kind=im), parameter :: nga(nbndlw) = &
           (/1,11,23,39,53,69,77,89,97,109,115,123,131,135,137,139/)

      integer(kind=im), parameter :: ngs(nbndlw) = &
           (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)

      integer(kind=im), parameter :: ngn(ngptlw) = &
           (/1,1,2,2,2,2,2,2,1,1, &                       ! band 1
             1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
             1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
             1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
             1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
             2,2,2,2,2,2,2,2, &                           ! band 6
             2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
             2,2,2,2,2,2,2,2, &                           ! band 8
             1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
             2,2,2,2,4,4, &                               ! band 10
             1,1,2,2,2,2,3,3, &                           ! band 11
             1,1,1,1,2,2,4,4, &                           ! band 12
             3,3,4,6, &                                   ! band 13
             8,8, &                                       ! band 14
             8,8, &                                       ! band 15
             4,12/)                                       ! band 16

      integer(kind=im), parameter :: ngb(ngptlw) = &
           (/1,1,1,1,1,1,1,1,1,1, &                       ! band 1
             2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
             3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
             4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
             5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
             6,6,6,6,6,6,6,6, &                           ! band 6
             7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
             8,8,8,8,8,8,8,8, &                           ! band 8
             9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
             10,10,10,10,10,10, &                         ! band 10
             11,11,11,11,11,11,11,11, &                   ! band 11
             12,12,12,12,12,12,12,12, &                   ! band 12
             13,13,13,13, &                               ! band 13
             14,14, &                                     ! band 14
             15,15, &                                     ! band 15
             16,16/)                                      ! band 16

      integer(kind=im), parameter :: ngm(nbndlw*mg) = &
           (/1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
             1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
             1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
             1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
             1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
             1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
             1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
             1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
             1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
             1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
             1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
             1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
             1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
             1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
             1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
             1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16

      real(kind=rb), parameter :: wt(mg) = &
           (/ 0.1527534276_rb, 0.1491729617_rb, 0.1420961469_rb, &
              0.1316886544_rb, 0.1181945205_rb, 0.1019300893_rb, &
              0.0832767040_rb, 0.0626720116_rb, 0.0424925000_rb, &
              0.0046269894_rb, 0.0038279891_rb, 0.0030260086_rb, &
              0.0022199750_rb, 0.0014140010_rb, 0.0005330000_rb, &
              0.0000750000_rb/)

      real(kind=rb) :: rwgt(nbndlw*mg)

      integer(kind=im), parameter :: nxmol = maxxsec
      integer(kind=im), parameter :: ixindx(maxxsec) = (/1, 2, 3, 4/)

      end module rrlw_wvn
