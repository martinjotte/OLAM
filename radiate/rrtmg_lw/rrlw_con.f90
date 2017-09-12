      module rrlw_con

      use parkind, only : rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw constants

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

      real(kind=rb), parameter :: grav = 9.8066_rb           ! Acceleration of gravity
                                                             ! (m s-2)
      real(kind=rb), parameter :: planck = 6.62606876e-27_rb ! Planck constant
                                                             ! (ergs s; g cm2 s-1)
      real(kind=rb), parameter :: boltz = 1.3806503e-16_rb   ! Boltzmann constant
                                                             ! (ergs K-1; g cm2 s-2 K-1)
      real(kind=rb), parameter :: clight = 2.99792458e+10_rb ! Speed of light in a vacuum  
                                                             ! (cm s-1)
      real(kind=rb), parameter :: avogad = 6.02214199e+23_rb ! Avogadro constant
                                                             ! (mol-1)
      real(kind=rb), parameter :: alosmt = 2.6867775e+19_rb  ! Loschmidt constant
                                                             ! (cm-3)
      real(kind=rb), parameter :: gascon = 8.31447200e+07_rb ! Molar gas constant
                                                             ! (ergs mol-1 K-1)
      real(kind=rb), parameter :: radcn1 = 1.191042722e-12_rb! First radiation constant
                                                             ! (W cm2 sr-1)
      real(kind=rb), parameter :: radcn2 = 1.4387752_rb      ! Second radiation constant
                                                             ! (cm K)
      real(kind=rb), parameter :: sbcnst = 5.670400e-04_rb   ! Stefan-Boltzmann constant
                                                             ! (W cm-2 K-4)
      real(kind=rb), parameter :: secdy = 8.6400e4_rb        ! Number of seconds per day
                                                             ! (s d-1)
      real(kind=rb), parameter :: oneminus = 1._rb - 1.e-6_rb

      real(kind=rb), parameter :: pi = 2._rb * asin(1._rb)

      real(kind=rb), parameter :: fluxfac = pi * 2.e4_rb     ! orig:   fluxfac = pi * 2.d4  

      real(kind=rb), parameter :: heatfac = grav * secdy / (1004._rb * 1.e2_rb)

      end module rrlw_con

