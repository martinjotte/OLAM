      module rrlw_kg01

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
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
! kao_mn2 : real
! kbo_mn2 : real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no1 = 16

      real(kind=rb) :: fracrefao(no1)  , fracrefbo(no1)
      real(kind=rb) :: kao(5,13,no1)
      real(kind=rb) :: kbo(5,13:59,no1)
      real(kind=rb) :: kao_mn2(19,no1) , kbo_mn2(19,no1)
      real(kind=rb) :: selfrefo(10,no1), forrefo(4,no1)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
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
! ka_mn2  : real
! kb_mn2  : real
! selfref : real
! forref  : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng1 = 10

      real(kind=rb) :: fracrefa(ng1)  , fracrefb(ng1)
      real(kind=rb) :: absa(ng1,65)
      real(kind=rb) :: absb(ng1,235)
      real(kind=rb) :: ka_mn2(ng1,19) , kb_mn2(ng1,19)
      real(kind=rb) :: dka_mn2(ng1,19), dkb_mn2(ng1,19)
      real(kind=rb) :: selfref(ng1,10), forref(ng1,4)

      end module rrlw_kg01
