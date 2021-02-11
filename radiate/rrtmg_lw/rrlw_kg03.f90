      module rrlw_kg03

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
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
! kao_mn2o: real
! kbo_mn2o: real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no3 = 16

      real(kind=rb) :: fracrefao(no3,9) ,fracrefbo(no3,5)
      real(kind=rb) :: kao(9,5,13,no3)
      real(kind=rb) :: kbo(5,5,13:59,no3)
      real(kind=rb) :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
      real(kind=rb) :: selfrefo(10,no3)
      real(kind=rb) :: forrefo(4,no3)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
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
! ka_mn2o : real
! kb_mn2o : real
! selfref : real
! forref  : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng3 = 16

      real(kind=rb) :: fracrefa(ng3,9) ,fracrefb(ng3,5)
      real(kind=rb) :: dfracrefa(ng3,9),dfracrefb(ng3,5)
      real(kind=rb) :: absa(ng3,585)
      real(kind=rb) :: absb(ng3,1175)
      real(kind=rb) :: ka_mn2o(ng3,9,19), kb_mn2o(ng3,5,19)
      real(kind=rb) :: dka_mn2o(ng3,9,19), dkb_mn2o(ng3,5,19)
      real(kind=rb) :: selfref(ng3,10)
      real(kind=rb) :: forref(ng3,4)

      end module rrlw_kg03


