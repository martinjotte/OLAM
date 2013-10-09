subroutine rrtmg_raddriv(iw, ka, nrad, koff, rlong_previous)

  use mem_grid,    only: mza, zm, zt, glatw, glonw

  use mem_basic,   only: rho, press, theta, tair, sh_v

  use misc_coms,   only: io6, iswrtyp, ilwrtyp, time8

  use consts_coms, only: rocp, p00i, stefan, cp, eps_virt, eps_vapi

  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         albedt_beam, albedt_diffuse, rshort_diffuse, &
                         rlong_albedo

  use micro_coms,  only: ncat

  use rrtmg_sw_rad_nomcica, only: rrtmg_sw_nomcica
  use rrtmg_sw_rad,         only: rrtmg_sw

  implicit none

  integer, intent(in) :: iw
  integer, intent(in) :: ka
  integer, intent(in) :: nrad
  integer, intent(in) :: koff
  real,    intent(in) :: rlong_previous
  
  integer, parameter :: ncol   = 1
  integer, parameter :: nflglw = 0
  
  real :: plev(ncol, nrad+1)
  real :: tlev(ncol, nrad+1)

  real :: play    (ncol, nrad)
  real :: tlay    (ncol, nrad)
  real :: h2ovmr  (ncol, nrad)
  real :: co2vmr  (ncol, nrad)
  real :: o3vmr   (ncol, nrad)
  real :: o2vmr   (ncol, nrad)
  real :: ch4vmr  (ncol, nrad)
  real :: n2ovmr  (ncol, nrad)
  real :: cfc11vmr(ncol, nrad)
  real :: cfc12vmr(ncol, nrad)
  real :: cfc22vmr(ncol, nrad)
  real :: ccl4vmr (ncol, nrad)

  integer :: k, krad

! Set trace gas volume mixing ratios, 2005 values, IPCC (2007)

  real, parameter :: co2   = 379.e-6  ! carbon dioxide (379 ppmv)
  real, parameter :: ch4   = 1774.e-9 ! methane (1774 ppbv)
  real, parameter :: n2o   = 319.e-9  ! nitrous oxide (319 ppbv)
  real, parameter :: cfc11 = 0.251e-9 ! cfc-11 (251 ppt)
  real, parameter :: cfc12 = 0.538e-9 ! cfc-12 (538 ppt)
  real, parameter :: cfc22 = 0.169e-9 ! cfc-22 (169 ppt)
  real, parameter :: ccl4  = 0.093e-9 ! ccl4 (93 ppt)

! Set oxygen volume mixing ratio (for o2mmr=0.23143)

  real, parameter :: o2 = 0.209488

  do k = ka, mza-1
     krad = k - koff

   ! h2ovmr  (ncol,krad) = qv * eps_vapi
     h2ovmr  (ncol,krad) = sh_v(k,iw) * eps_vapi / (1.0 + sh_v(k,iw) * eps_virt)

  enddo

  do krad = 1, nrad
     co2vmr  (ncol,krad) = co2
     o2vmr   (ncol,krad) = o2
     ch4vmr  (ncol,krad) = ch4
     n2ovmr  (ncol,krad) = n2o
     cfc11vmr(ncol,krad) = cfc11
     cfc12vmr(ncol,krad) = cfc12
     cfc22vmr(ncol,krad) = cfc22
     ccl4vmr (ncol,krad) = ccl4
  enddo

end subroutine rrtmg_raddriv
