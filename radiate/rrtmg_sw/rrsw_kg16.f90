      module rrsw_kg16

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ng16

      implicit none
      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      integer(kind=im), parameter :: no16 = 16

      real(kind=rb) :: kao(9,5,13,no16)
      real(kind=rb) :: kbo(5,13:59,no16)
      real(kind=rb) :: selfrefo(10,no16), forrefo(3,no16)
      real(kind=rb) :: sfluxrefo(no16)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      real(kind=rb) :: absa(ng16,585)
      real(kind=rb) :: absb(ng16,235)
      real(kind=rb) :: selfref(ng16,10), forref(ng16,3)
      real(kind=rb) :: sfluxref(ng16)

      end module rrsw_kg16

