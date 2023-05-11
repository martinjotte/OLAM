      module rrlw_kg14

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
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

      integer(kind=im), parameter :: no14 = 16

      real(kind=rb) , dimension(no14) :: fracrefao
      real(kind=rb) , dimension(no14) :: fracrefbo

      real(kind=rb) :: kao(5,13,no14)
      real(kind=rb) :: kbo(5,13:59,no14)
      real(kind=rb) :: selfrefo(10,no14)
      real(kind=rb) :: forrefo(4,no14)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
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
! selfref : real
! forref  : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng14 = 2

      real(kind=rb) , dimension(ng14) :: fracrefa
      real(kind=rb) , dimension(ng14) :: fracrefb

      real(kind=rb) :: absa(ng14,65)
      real(kind=rb) :: absb(ng14,235)
      real(kind=rb) :: selfref(ng14,10)
      real(kind=rb) :: forref(ng14,4)

      end module rrlw_kg14
