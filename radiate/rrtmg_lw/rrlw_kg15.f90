      module rrlw_kg15

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real
! kao     : real
! kao_mn2 : real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no15 = 16

      real(kind=rb) :: fracrefao(no15,9)
      real(kind=rb) :: kao(9,5,13,no15)
      real(kind=rb) :: kao_mn2(9,19,no15)
      real(kind=rb) :: selfrefo(10,no15)
      real(kind=rb) :: forrefo(4,no15)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real
! ka      : real
! ka_mn2  : real
! selfref : real
! forref  : real
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng15 = 2

      real(kind=rb) :: fracrefa(ng15,9)
      real(kind=rb) :: dfracrefa(ng15,9)
      real(kind=rb) :: absa(ng15,585)
      real(kind=rb) :: ka_mn2(ng15,9,19)
      real(kind=rb) :: dka_mn2(ng15,9,19)
      real(kind=rb) :: selfref(ng15,10)
      real(kind=rb) :: forref(ng15,4)

      end module rrlw_kg15
