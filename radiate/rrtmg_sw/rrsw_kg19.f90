      module rrsw_kg19

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ng19

      implicit none
      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 19
! band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real
! kbo     : real
! selfrefo: real
! forrefo : real
!sfluxrefo: real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no19 = 16

      real(kind=rb) :: kao(9,5,13,no19)
      real(kind=rb) :: kbo(5,13:59,no19)
      real(kind=rb) :: selfrefo(10,no19), forrefo(3,no19)
      real(kind=rb) :: sfluxrefo(no19,9)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 19
! band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real
! kb      : real
! absa    : real
! absb    : real
! selfref : real
! forref  : real
! sfluxref: real
!-----------------------------------------------------------------

      real(kind=rb) :: absa(ng19,585)
      real(kind=rb) :: absb(ng19,235)
      real(kind=rb) :: selfref(ng19,10), forref(ng19,3)
      real(kind=rb) :: sfluxref(ng19,9), dsfluxref(ng19,9)

      end module rrsw_kg19

