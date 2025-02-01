      module rrsw_wvn

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : nbndsw, mg, ngptsw, jpb1, jpb2

      implicit none

      private :: im, rb, nbndsw, mg, ngptsw, jpb1, jpb2

!------------------------------------------------------------------
! rrtmg_sw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer:
! nspb   :  integer:
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
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
!                    (224 total) into reduced set of g-intervals
!                    (112 total)
!------------------------------------------------------------------

      integer(kind=im), parameter :: ng(jpb1:jpb2) = &
           (/16,16,16,16,16,16,16,16,16,16,16,16,16,16/)

      integer(kind=im), parameter :: nspa(jpb1:jpb2) = &
           (/9,9,9,9,1,9,9,1,9,1,0,1,9,1/)

      integer(kind=im), parameter :: nspb(jpb1:jpb2) = &
           (/1,5,1,1,1,5,1,0,1,0,0,1,5,1/)

      real(kind=rb), parameter :: wavenum1(jpb1:jpb2) =                                     &
                   (/ 2600._rb, 3250._rb, 4000._rb, 4650._rb, 5150._rb, 6150._rb, 7700._rb, &
                      8050._rb,12850._rb,16000._rb,22650._rb,29000._rb,38000._rb,  820._rb /)

      real(kind=rb), parameter :: wavenum2(jpb1:jpb2) =                                     &
                   (/ 3250._rb, 4000._rb, 4650._rb, 5150._rb, 6150._rb, 7700._rb, 8050._rb, &
                     12850._rb,16000._rb,22650._rb,29000._rb,38000._rb,50000._rb, 2600._rb /)

      real(kind=rb), parameter :: delwave(jpb1:jpb2) =                                      &
                   (/  650._rb,  750._rb,  650._rb,  500._rb, 1000._rb, 1550._rb,  350._rb, &
                      4800._rb, 3150._rb, 6650._rb, 6350._rb, 9000._rb,12000._rb, 1780._rb /)

      integer(kind=im), parameter :: ngc(nbndsw) = &
           (/ 6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8, 6,12 /)

      integer(kind=im), parameter :: nga(nbndsw) = &
           (/ 1,7,19,27,35,45,55,57,67,75,81,87,95,101 /)

      integer(kind=im), parameter :: ngs(nbndsw) = &
           (/ 6,18,26,34,44,54,56,66,74,80,86,94,100,112 /)

      integer(kind=im), parameter :: ngn(ngptsw) = &
           (/ 2,2,2,2,4,4, &                               ! band 16
              1,1,1,1,1,2,1,2,1,2,1,2, &                   ! band 17
              1,1,1,1,2,2,4,4, &                           ! band 18
              1,1,1,1,2,2,4,4, &                           ! band 19
              1,1,1,1,1,1,1,1,2,6, &                       ! band 20
              1,1,1,1,1,1,1,1,2,6, &                       ! band 21
              8,8, &                                       ! band 22
              2,2,1,1,1,1,1,1,2,4, &                       ! band 23
              2,2,2,2,2,2,2,2, &                           ! band 24
              1,1,2,2,4,6, &                               ! band 25
              1,1,2,2,4,6, &                               ! band 26
              1,1,1,1,1,1,4,6, &                           ! band 27
              1,1,2,2,4,6, &                               ! band 28
              1,1,1,1,2,2,2,2,1,1,1,1 /)                   ! band 29

      integer(kind=im), parameter :: ngbo(ngptsw) = &
           (/ 16,16,16,16,16,16, &                         ! band 16
              17,17,17,17,17,17,17,17,17,17,17,17, &       ! band 17
              18,18,18,18,18,18,18,18, &                   ! band 18
              19,19,19,19,19,19,19,19, &                   ! band 19
              20,20,20,20,20,20,20,20,20,20, &             ! band 20
              21,21,21,21,21,21,21,21,21,21, &             ! band 21
              22,22, &                                     ! band 22
              23,23,23,23,23,23,23,23,23,23, &             ! band 23
              24,24,24,24,24,24,24,24, &                   ! band 24
              25,25,25,25,25,25, &                         ! band 25
              26,26,26,26,26,26, &                         ! band 26
              27,27,27,27,27,27,27,27, &                   ! band 27
              28,28,28,28,28,28, &                         ! band 28
              29,29,29,29,29,29,29,29,29,29,29,29 /)       ! band 29

      integer(kind=im), parameter :: ngb(ngptsw) = ngbo - 15

      integer(kind=im), parameter :: ngm(nbndsw*mg) = &
           (/ 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 16
              1,2,3,4,5,6,6,7,8,8,9,10,10,11,12,12, &      ! band 17
              1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 18
              1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 19
              1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 20
              1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 21
              1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 22
              1,1,2,2,3,4,5,6,7,8,9,9,10,10,10,10, &       ! band 23
              1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 24
              1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 25
              1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 26
              1,2,3,4,5,6,7,7,7,7,8,8,8,8,8,8, &           ! band 27
              1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 28
              1,2,3,4,5,5,6,6,7,7,8,8,9,10,11,12 /)        ! band 29

      real(kind=rb), parameter :: wt(mg) = &
           (/ 0.1527534276_rb, 0.1491729617_rb, 0.1420961469_rb, &
              0.1316886544_rb, 0.1181945205_rb, 0.1019300893_rb, &
              0.0832767040_rb, 0.0626720116_rb, 0.0424925000_rb, &
              0.0046269894_rb, 0.0038279891_rb, 0.0030260086_rb, &
              0.0022199750_rb, 0.0014140010_rb, 0.0005330000_rb, &
              0.0000750000_rb /)

      real(kind=rb) :: rwgt(nbndsw*mg)

      real(kind=rb) :: raylt(ngptsw)

      end module rrsw_wvn
