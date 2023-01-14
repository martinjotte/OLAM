      module rrsw_kg20

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ng20

      implicit none
      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
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
! absch4o : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no20 = 16

      real(kind=rb) :: kao(5,13,no20)
      real(kind=rb) :: kbo(5,13:59,no20)
      real(kind=rb) :: selfrefo(10,no20), forrefo(4,no20)
      real(kind=rb) :: sfluxrefo(no20)
      real(kind=rb) :: absch4o(no20)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
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
! absch4  : real
!-----------------------------------------------------------------

      real(kind=rb) :: absa(ng20,65)
      real(kind=rb) :: absb(ng20,235)
      real(kind=rb) :: selfref(ng20,10), forref(ng20,4)
      real(kind=rb) :: sfluxref(ng20)
      real(kind=rb) :: absch4(ng20)

      end module rrsw_kg20

