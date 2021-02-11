      module rrsw_aer

      use parkind, only: im => kind_im, rb => kind_rb
      use parrrsw, only: nbndsw, naerec

      implicit none

!------------------------------------------------------------------
! rrtmg_sw aerosol optical properties
!
!  Data derived from six ECMWF aerosol types and defined for
!  the rrtmg_sw spectral intervals
!
! Initial: J.-J. Morcrette, ECMWF, mar2003
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------
!
!-- The six ECMWF aerosol types are respectively:
!
!  1/ continental average                 2/ maritime
!  3/ desert                              4/ urban
!  5/ volcanic active                     6/ stratospheric background
!
! computed from Hess and Koepke (con, mar, des, urb)
!          from Bonnel et al.   (vol, str)
!
! rrtmg_sw 14 spectral intervals (microns):
!  3.846 -  3.077
!  3.077 -  2.500
!  2.500 -  2.150
!  2.150 -  1.942
!  1.942 -  1.626
!  1.626 -  1.299
!  1.299 -  1.242
!  1.242 -  0.7782
!  0.7782-  0.6250
!  0.6250-  0.4415
!  0.4415-  0.3448
!  0.3448-  0.2632
!  0.2632-  0.2000
! 12.195 -  3.846
!
!------------------------------------------------------------------
!
!  name     type     purpose
! -----   : ----   : ----------------------------------------------
! rsrtaua : real   : ratio of average optical thickness in 
!                    spectral band to that at 0.55 micron
! rsrpiza : real   : average single scattering albedo (unitless)
! rsrasya : real   : average asymmetry parameter (unitless)
!------------------------------------------------------------------

      real(kind=rb) :: rsrtaua(naerec,nbndsw)
      real(kind=rb) :: rsrpiza(naerec,nbndsw)
      real(kind=rb) :: rsrasya(naerec,nbndsw)

      contains

      subroutine swaerpr

      implicit none

      rsrtaua = reshape( shape = (/ naerec,nbndsw /), source = (/               &
        0.10849_rb, 0.66699_rb, 0.65255_rb, 0.11600_rb, 0.06529_rb, 0.04468_rb, &
        0.10849_rb, 0.66699_rb, 0.65255_rb, 0.11600_rb, 0.06529_rb, 0.04468_rb, &
        0.20543_rb, 0.84642_rb, 0.84958_rb, 0.21673_rb, 0.28270_rb, 0.10915_rb, &
        0.20543_rb, 0.84642_rb, 0.84958_rb, 0.21673_rb, 0.28270_rb, 0.10915_rb, &
        0.20543_rb, 0.84642_rb, 0.84958_rb, 0.21673_rb, 0.28270_rb, 0.10915_rb, &
        0.20543_rb, 0.84642_rb, 0.84958_rb, 0.21673_rb, 0.28270_rb, 0.10915_rb, &
        0.20543_rb, 0.84642_rb, 0.84958_rb, 0.21673_rb, 0.28270_rb, 0.10915_rb, &
        0.52838_rb, 0.93285_rb, 0.93449_rb, 0.53078_rb, 0.67148_rb, 0.46608_rb, &
        0.52838_rb, 0.93285_rb, 0.93449_rb, 0.53078_rb, 0.67148_rb, 0.46608_rb, &
        1.69446_rb, 1.11855_rb, 1.09212_rb, 1.72145_rb, 1.03858_rb, 1.12044_rb, &
        1.69446_rb, 1.11855_rb, 1.09212_rb, 1.72145_rb, 1.03858_rb, 1.12044_rb, &
        1.69446_rb, 1.11855_rb, 1.09212_rb, 1.72145_rb, 1.03858_rb, 1.12044_rb, &
        1.69446_rb, 1.11855_rb, 1.09212_rb, 1.72145_rb, 1.03858_rb, 1.12044_rb, &
        0.10849_rb, 0.66699_rb, 0.65255_rb, 0.11600_rb, 0.06529_rb, 0.04468_rb /))

      rsrpiza = reshape( shape = (/ naerec,nbndsw /), source = (/                     &
        .5230504_rb, .7868518_rb, .8531531_rb, .4048149_rb, .8748231_rb, .2355667_rb, &
        .5230504_rb, .7868518_rb, .8531531_rb, .4048149_rb, .8748231_rb, .2355667_rb, &
        .8287144_rb, .9949396_rb, .9279543_rb, .6765051_rb, .9467578_rb, .9955938_rb, &
        .8287144_rb, .9949396_rb, .9279543_rb, .6765051_rb, .9467578_rb, .9955938_rb, &
        .8287144_rb, .9949396_rb, .9279543_rb, .6765051_rb, .9467578_rb, .9955938_rb, &
        .8287144_rb, .9949396_rb, .9279543_rb, .6765051_rb, .9467578_rb, .9955938_rb, &
        .8287144_rb, .9949396_rb, .9279543_rb, .6765051_rb, .9467578_rb, .9955938_rb, &
        .8970131_rb, .9984940_rb, .9245594_rb, .7768385_rb, .9532763_rb, .9999999_rb, &
        .8970131_rb, .9984940_rb, .9245594_rb, .7768385_rb, .9532763_rb, .9999999_rb, &
        .9148907_rb, .9956173_rb, .7504584_rb, .8131335_rb, .9401905_rb, .9999999_rb, &
        .9148907_rb, .9956173_rb, .7504584_rb, .8131335_rb, .9401905_rb, .9999999_rb, &
        .9148907_rb, .9956173_rb, .7504584_rb, .8131335_rb, .9401905_rb, .9999999_rb, &
        .9148907_rb, .9956173_rb, .7504584_rb, .8131335_rb, .9401905_rb, .9999999_rb, &
        .5230504_rb, .7868518_rb, .8531531_rb, .4048149_rb, .8748231_rb, .2355667_rb /))

      rsrasya = reshape( shape = (/ naerec,nbndsw /), source = (/                     &
        0.700610_rb, 0.818871_rb, 0.702399_rb, 0.689886_rb, .4629866_rb, .1907639_rb, &
        0.700610_rb, 0.818871_rb, 0.702399_rb, 0.689886_rb, .4629866_rb, .1907639_rb, &
        0.636342_rb, 0.802467_rb, 0.691305_rb, 0.627497_rb, .6105750_rb, .4760794_rb, &
        0.636342_rb, 0.802467_rb, 0.691305_rb, 0.627497_rb, .6105750_rb, .4760794_rb, &
        0.636342_rb, 0.802467_rb, 0.691305_rb, 0.627497_rb, .6105750_rb, .4760794_rb, &
        0.636342_rb, 0.802467_rb, 0.691305_rb, 0.627497_rb, .6105750_rb, .4760794_rb, &
        0.636342_rb, 0.802467_rb, 0.691305_rb, 0.627497_rb, .6105750_rb, .4760794_rb, &
        0.668431_rb, 0.788530_rb, 0.698682_rb, 0.657422_rb, .6735182_rb, .6519706_rb, &
        0.668431_rb, 0.788530_rb, 0.698682_rb, 0.657422_rb, .6735182_rb, .6519706_rb, &
        0.729019_rb, 0.803129_rb, 0.784592_rb, 0.712208_rb, .7008249_rb, .7270548_rb, &
        0.729019_rb, 0.803129_rb, 0.784592_rb, 0.712208_rb, .7008249_rb, .7270548_rb, &
        0.729019_rb, 0.803129_rb, 0.784592_rb, 0.712208_rb, .7008249_rb, .7270548_rb, &
        0.729019_rb, 0.803129_rb, 0.784592_rb, 0.712208_rb, .7008249_rb, .7270548_rb, &
        0.700610_rb, 0.818871_rb, 0.702399_rb, 0.689886_rb, .4629866_rb, .1907639_rb /))

      end subroutine swaerpr

      end module rrsw_aer

