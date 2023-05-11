      module rrlw_kg07

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! kao_mco2: real
! kbo_mco2: real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no7 = 16

      real(kind=rb) , dimension(no7) :: fracrefbo
      real(kind=rb) :: fracrefao(no7,9)
      real(kind=rb) :: kao(9,5,13,no7)
      real(kind=rb) :: kbo(5,13:59,no7)
      real(kind=rb) :: kao_mco2(9,19,no7)
      real(kind=rb) :: kbo_mco2(19,no7)
      real(kind=rb) :: selfrefo(10,no7)
      real(kind=rb) :: forrefo(4,no7)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! ka_mco2 : real
! kb_mco2 : real
! selfref : real
! forref  : real
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng7 = 12

      real(kind=rb) , dimension(ng7) :: fracrefb
      real(kind=rb) :: fracrefa(ng7,9)
      real(kind=rb) :: dfracrefa(ng7,9)
      real(kind=rb) :: absa(ng7,585)
      real(kind=rb) :: absb(ng7,235)
      real(kind=rb) :: ka_mco2(ng7,9,19)
      real(kind=rb) :: kb_mco2(ng7,19)
      real(kind=rb) :: dka_mco2(ng7,9,19)
      real(kind=rb) :: dkb_mco2(ng7,19)
      real(kind=rb) :: selfref(ng7,10)
      real(kind=rb) :: forref(ng7,4)

      end module rrlw_kg07
