      module rrsw_con

      use parkind, only : im => kind_im, rb => kind_rb

      implicit none

!------------------------------------------------------------------
! rrtmg_sw constants

! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! fluxfac:  real   : radiance to flux conversion factor
! heatfac:  real   : flux to heating rate conversion factor
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
! grav   :  real   : acceleration of gravity
! planck :  real   : planck constant
! boltz  :  real   : boltzmann constant
! clight :  real   : speed of light
! avogad :  real   : avogadro constant
! alosmt :  real   : loschmidt constant
! gascon :  real   : molar gas constant
! radcn1 :  real   : first radiation constant
! radcn2 :  real   : second radiation constant
! sbcnst :  real   : stefan-boltzmann constant
!  secdy :  real   : seconds per day
!------------------------------------------------------------------

!     real(kind=rb) :: fluxfac, heatfac

      real(kind=rb), parameter :: zepsec = 1.e-06_rb
      real(kind=rb), parameter :: zepzen = 1.e-10_rb
      real(kind=rb), parameter :: oneminus = 1.0_rb - zepsec
      real(kind=rb), parameter :: pi = 2._rb * asin(1._rb)

      real(kind=rb), parameter :: grav = 9.8066_rb           ! (m s-2)
      real(kind=rb), parameter :: planck = 6.62606876e-27_rb ! (ergs s; g cm2 s-1)
      real(kind=rb), parameter :: boltz = 1.3806503e-16_rb   ! (ergs K-1; g cm2 s-2 K-1)
      real(kind=rb), parameter :: clight = 2.99792458e+10_rb ! (cm s-1)
      real(kind=rb), parameter :: avogad = 6.02214199e+23_rb ! (mol-1)
      real(kind=rb), parameter :: alosmt = 2.6867775e+19_rb  ! (cm-3)
      real(kind=rb), parameter :: gascon = 8.31447200e+07_rb ! (ergs mol-1 K-1)
      real(kind=rb), parameter :: radcn1 = 1.191042772e-12_rb! (W cm2 sr-1)
      real(kind=rb), parameter :: radcn2 = 1.4387752_rb      ! (cm K)
      real(kind=rb), parameter :: sbcnst = 5.670400e-04_rb   ! (W cm-2 K-4)
      real(kind=rb), parameter :: secdy = 8.6400e4_rb        ! (s d-1)

      end module rrsw_con

