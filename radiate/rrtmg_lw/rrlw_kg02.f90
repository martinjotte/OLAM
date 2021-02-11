      module rrlw_kg02

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real
!fracrefbo: real
! kao     : real
! kbo     : real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no2 = 16

      real(kind=rb) :: fracrefao(no2)   , fracrefbo(no2)
      real(kind=rb) :: kao(5,13,no2)
      real(kind=rb) :: kbo(5,13:59,no2)
      real(kind=rb) :: selfrefo(10,no2) , forrefo(4,no2)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real
!fracrefb : real
! ka      : real
! kb      : real
! absa    : real
! absb    : real
! selfref : real
! forref  : real
!
! refparam: real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng2 = 12

      real(kind=rb) :: fracrefa(ng2)  , fracrefb(ng2)
      real(kind=rb) :: absa(ng2,65)
      real(kind=rb) :: absb(ng2,235)
      real(kind=rb) :: selfref(ng2,10), forref(ng2,4)

      real(kind=rb) :: refparam(13)

      end module rrlw_kg02


