      module rrlw_kg16

      use parkind, only: im => kind_im, rb => kind_rb

      implicit none

      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
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
! kbo     : real
! selfrefo: real
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no16 = 16

      real(kind=rb) , dimension(no16) :: fracrefbo

      real(kind=rb) :: fracrefao(no16,9)
      real(kind=rb) :: kao(9,5,13,no16)
      real(kind=rb) :: kbo(5,13:59,no16)
      real(kind=rb) :: selfrefo(10,no16)
      real(kind=rb) :: forrefo(4,no16)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
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
! kb      : real
! selfref : real
! forref  : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng16 = 2

      real(kind=rb) , dimension(ng16) :: fracrefb

      real(kind=rb) :: fracrefa(ng16,9)
      real(kind=rb) :: dfracrefa(ng16,9)
      real(kind=rb) :: absa(ng16,585)
      real(kind=rb) :: absb(ng16,235)
      real(kind=rb) :: selfref(ng16,10)
      real(kind=rb) :: forref(ng16,4)

      end module rrlw_kg16

