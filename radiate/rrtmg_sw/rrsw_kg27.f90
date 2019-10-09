      module rrsw_kg27

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: ng27

      implicit none
      private :: im, rb

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
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
!sfluxrefo: real     
! raylo   : real     
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no27 = 16

      real(kind=rb) :: kao(5,13,no27)
      real(kind=rb) :: kbo(5,13:59,no27)
      real(kind=rb) :: sfluxrefo(no27)
      real(kind=rb) :: raylo(no27)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
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
! sfluxref: real     
! rayl    : real     
!-----------------------------------------------------------------

      real(kind=rb) :: absa(ng27,65)
      real(kind=rb) :: absb(ng27,235)
      real(kind=rb) :: sfluxref(ng27)
      real(kind=rb) :: rayl(ng27)

      end module rrsw_kg27

