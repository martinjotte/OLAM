      module rrlw_kg08

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! kao_mn2o: real
! kbo_mn2o: real
! kao_mo3 : real
! selfrefo: real
! forrefo : real
! cfc12o  : real
!cfc22adjo: real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no8 = 16

      real(kind=rb) , dimension(no8) :: fracrefao
      real(kind=rb) , dimension(no8) :: fracrefbo
      real(kind=rb) , dimension(no8) :: cfc12o
      real(kind=rb) , dimension(no8) :: cfc22adjo

      real(kind=rb) :: kao(5,13,no8)
      real(kind=rb) :: kao_mco2(19,no8)
      real(kind=rb) :: kao_mn2o(19,no8)
      real(kind=rb) :: kao_mo3(19,no8)
      real(kind=rb) :: kbo(5,13:59,no8)
      real(kind=rb) :: kbo_mco2(19,no8)
      real(kind=rb) :: kbo_mn2o(19,no8)
      real(kind=rb) :: selfrefo(10,no8)
      real(kind=rb) :: forrefo(4,no8)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! ka_mn2o : real
! kb_mn2o : real
! ka_mo3  : real
! selfref : real
! forref  : real
! cfc12   : real
! cfc22adj: real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng8 = 8

      real(kind=rb) , dimension(ng8) :: fracrefa
      real(kind=rb) , dimension(ng8) :: fracrefb
      real(kind=rb) , dimension(ng8) :: cfc12
      real(kind=rb) , dimension(ng8) :: cfc22adj

      real(kind=rb) :: absa(ng8,65)
      real(kind=rb) :: absb(ng8,235)
      real(kind=rb) :: ka_mco2(ng8,19)
      real(kind=rb) :: ka_mn2o(ng8,19)
      real(kind=rb) :: ka_mo3(ng8,19)
      real(kind=rb) :: kb_mco2(ng8,19)
      real(kind=rb) :: kb_mn2o(ng8,19)
      real(kind=rb) :: dka_mco2(ng8,19)
      real(kind=rb) :: dka_mn2o(ng8,19)
      real(kind=rb) :: dka_mo3(ng8,19)
      real(kind=rb) :: dkb_mco2(ng8,19)
      real(kind=rb) :: dkb_mn2o(ng8,19)
      real(kind=rb) :: selfref(ng8,10)
      real(kind=rb) :: forref(ng8,4)

      end module rrlw_kg08

