      module rrsw_kg29

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ng29

      implicit none
      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 29
! band 29:  820-2600 cm-1 (low - h2o; high - co2)
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
! absh2oo : real
! absco2o : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no29 = 16

      real(kind=rb) :: kao(5,13,no29)
      real(kind=rb) :: kbo(5,13:59,no29)
      real(kind=rb) :: selfrefo(10,no29), forrefo(4,no29)
      real(kind=rb) :: sfluxrefo(no29)
      real(kind=rb) :: absh2oo(no29), absco2o(no29)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 29
! band 29:  820-2600 cm-1 (low - h2o; high - co2)
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
! selfref : real
! forref  : real
! sfluxref: real
! absh2o  : real
! absco2  : real
!-----------------------------------------------------------------

      real(kind=rb) :: absa(ng29,65)
      real(kind=rb) :: absb(ng29,235)

      real(kind=rb) :: selfref(ng29,10), forref(ng29,4)
      real(kind=rb) :: sfluxref(ng29)
      real(kind=rb) :: absh2o(ng29), absco2(ng29)

      end module rrsw_kg29

