      module rrlw_kg04

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
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

      integer(kind=im), parameter :: no4 = 16

      real(kind=rb) :: fracrefao(no4,9)  ,fracrefbo(no4,5)
      real(kind=rb) :: kao(9,5,13,no4)
      real(kind=rb) :: kbo(5,5,13:59,no4)
      real(kind=rb) :: selfrefo(10,no4)  ,forrefo(4,no4)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! absa    : real
! absb    : real
!fracrefa : real
!fracrefb : real
! ka      : real
! kb      : real
! selfref : real
! forref  : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng4  = 14

      real(kind=rb) :: fracrefa(ng4,9) , fracrefb(ng4,5)
      real(kind=rb) :: dfracrefa(ng4,9), dfracrefb(ng4,5)
      real(kind=rb) :: absa(ng4,585)
      real(kind=rb) :: absb(ng4,1175)
      real(kind=rb) :: selfref(ng4,10)  ,forref(ng4,4)

      end module rrlw_kg04
