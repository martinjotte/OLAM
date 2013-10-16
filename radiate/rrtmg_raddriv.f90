!!module rrtmg_drivermod
!!  implicit none
!!
!!  integer, parameter :: ncol   = 1
!!  integer, parameter :: icld   = 2
!!  integer, parameter :: dyofyr = 0
!!  integer, parameter :: inflg  = 0
!!  integer, parameter :: iceflg = 0
!!  integer, parameter :: liqflg = 0
!!  
!!contains
!!
!!  subroutine alloc_rrtmg_columns
!!
!!    integer, intent(in) :: 
!!
!!

subroutine rrtmg_raddriv(iw, ka, nrad, koff)

  use mem_grid,    only: mza, zm, zt, glatw, glonw, dzt, dzim

  use mem_basic,   only: rho, press, theta, tair, sh_v

  use misc_coms,   only: io6, iswrtyp, ilwrtyp, time8

  use consts_coms, only: stefan, eps_virt, eps_vapi, grav, solar, cp

  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         albedt_beam, albedt_diffuse, rshort_diffuse, &
                         rlong_albedo, solfac

  use micro_coms,  only: ncat

  use parrrtm,              only: nbndlw
  use parrrsw,              only: nbndsw
  use rrtmg_sw_rad_nomcica, only: rrtmg_sw_nomcica
  use rrtmg_sw_rad,         only: rrtmg_sw
  use rrtmg_lw_rad_nomcica, only: rrtmg_lw_nomcica
  use rrtmg_lw_rad,         only: rrtmg_lw

  implicit none

  integer, intent(in) :: iw
  integer, intent(in) :: ka
  integer, intent(in) :: nrad
  integer, intent(in) :: koff
  
  integer, parameter :: ncol   = 1
  integer, parameter :: icld   = 2
! integer, parameter :: icld   = 0
  integer, parameter :: dyofyr = 0
  integer, parameter :: inflg  = 0
  integer, parameter :: iceflg = 0
  integer, parameter :: liqflg = 0

  integer, parameter :: idrv   = 0

  real :: coszen(ncol)
  real :: tsfc  (ncol)
  real :: asdir (ncol)
  real :: aldir (ncol)
  real :: asdif (ncol)
  real :: aldif (ncol)
  real :: emis  (ncol, nbndlw)

  real :: zsfc
  real :: emiss
  real :: p1, p2

  real :: plev(ncol, nrad+1)
  real :: tlev(ncol, nrad+1)

  real :: play    (ncol, nrad)
  real :: tlay    (ncol, nrad)
  real :: h2ovmr  (ncol, nrad)
  real :: cldfr   (ncol, nrad)
  real :: cicewp  (ncol, nrad)
  real :: cliqwp  (ncol, nrad)
  real :: reice   (ncol, nrad)
  real :: reliq   (ncol, nrad)

  real :: co2vmr  (ncol, nrad)
  real :: o3vmr   (ncol, nrad)
  real :: o2vmr   (ncol, nrad)
  real :: ch4vmr  (ncol, nrad)
  real :: n2ovmr  (ncol, nrad)
  real :: cfc11vmr(ncol, nrad)
  real :: cfc12vmr(ncol, nrad)
  real :: cfc22vmr(ncol, nrad)
  real :: ccl4vmr (ncol, nrad)

  real :: taucldl(nbndlw, ncol, nrad)
  real :: tauclds(nbndsw, ncol, nrad)
  real :: ssaclds(nbndsw, ncol, nrad)
  real :: asmclds(nbndsw, ncol, nrad)
  real :: fsfclds(nbndsw, ncol, nrad)

  real :: dl (nrad)
  real :: zml(nrad)
  real :: ztl(nrad)
  real :: exl(nrad)
  real :: dzl(nrad)

  real :: swuflx (ncol, nrad+1)
  real :: swdflx (ncol, nrad+1)
  real :: swhr   (ncol, nrad  )
  real :: swuflxc(ncol, nrad+1)
  real :: swdflxc(ncol, nrad+1)
  real :: swhrc  (ncol, nrad  )

  real :: lwuflx (ncol, nrad+1) 
  real :: lwdflx (ncol, nrad+1)
  real :: lwhr   (ncol, nrad  )
  real :: lwuflxc(ncol, nrad+1)
  real :: lwdflxc(ncol, nrad+1)
  real :: lwhrc  (ncol, nrad  )

  real :: duflx_dt (ncol, nrad+1)
  real :: duflxc_dt(ncol, nrad+1)

  real :: tauaerl(ncol, nrad, nbndlw)
  real :: tauaers(ncol, nrad, nbndsw)
  real :: ssaaers(ncol, nrad, nbndsw)
  real :: asmaers(ncol, nrad, nbndsw)
  real :: ecaer  (ncol, nrad, nbndsw)

  integer :: k, krad, icloud

! Set gas volume mixing ratios, 2005 values, IPCC (2007)

  real, parameter :: co2   = 379.e-6  ! carbon dioxide (379 ppmv)
  real, parameter :: ch4   = 1774.e-9 ! methane (1774 ppbv)
  real, parameter :: n2o   = 319.e-9  ! nitrous oxide (319 ppbv)
  real, parameter :: cfc11 = 0.251e-9 ! cfc-11 (251 ppt)
  real, parameter :: cfc12 = 0.538e-9 ! cfc-12 (538 ppt)
  real, parameter :: cfc22 = 0.169e-9 ! cfc-22 (169 ppt)
  real, parameter :: ccl4  = 0.093e-9 ! ccl4 (93 ppt)
  real, parameter :: o2    = 0.209488 ! oxygen, for o2mmr=0.23143

! Molecular weight of dry air / ozone

  real, parameter :: amdo3 = 0.603428

! Set some surface values needed by RRTMg

  emiss = 1.0 - rlong_albedo(iw)
  zsfc  = zm(ka-1)

  asdir(ncol) = albedt_beam(iw)
  aldir(ncol) = albedt_beam(iw)
  asdif(ncol) = albedt_diffuse(iw)
  aldif(ncol) = albedt_diffuse(iw)

  coszen(ncol) = cosz(iw)
  emis(ncol,:) = emiss
  tsfc(ncol)   = (rlongup(iw) / emiss / stefan) ** 0.25

! Copy column values from model to radiation memory space

  do k = ka, mza-1
     krad = k - koff

     play  (1,krad) = press(k,iw)
     tlay  (1,krad) = tair (k,iw)
     h2ovmr(1,krad) = sh_v (k,iw)
     dl      (krad) = rho  (k,iw)
     zml     (krad) = zm   (k)
     ztl     (krad) = zt   (k)
     exl     (krad) = theta(k,iw) / tair(k,iw)
     dzl     (krad) = dzt  (k)
  enddo

! Fill ozone column and any extra radiation layers at model top

  call rad_mclat(iw,nrad,koff,glatw(iw),dl,play,h2ovmr,tlay,o3vmr,zml,ztl,dzl)

! Gases not defined in OLAM

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

! Convert water vapor to molar mixing ratio, pressure to mb, and
! ozone to molar mixing ratio

  do krad = 1, nrad
     h2ovmr(ncol,krad) = h2ovmr(ncol,krad) * eps_vapi
     play  (ncol,krad) = play  (ncol,krad) * 0.01
     o3vmr (ncol,krad) = o3vmr (ncol,krad) * amdo3
  enddo

! surface pressure and temperature

  plev(ncol,1) = play(ncol,1) + (ztl(1) - zsfc) * dl(1) * grav * 0.01
  tlev(ncol,1) = tsfc(ncol)

! pressure and temperature at intermediate levels

  do krad = 2, nrad
     p1 = play(ncol,krad-1) + (ztl(krad-1)-zml(krad-1)) * dl(krad-1) * grav * 0.01
     p2 = play(ncol,krad)   + (ztl(krad)  -zml(krad-1)) * dl(krad)   * grav * 0.01
     plev(ncol,krad) = 0.5 * ( p1 + p2 )
     tlev(ncol,krad) = 0.5 * ( tlay(ncol,krad-1) + tlay(ncol,krad) )
  enddo
  
! pressure and temperature at top level

  plev(ncol,nrad+1) = play(ncol,nrad) + (ztl(nrad) - zml(nrad)) * dl(nrad) * grav * 0.01
  tlev(ncol,nrad+1) = tlay(ncol,nrad)

! NO CLOUDS FOR NOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  cldfr    (ncol,:) = 0.0
  taucldl(:,ncol,:) = 0.0
  tauclds(:,ncol,:) = 0.0
  ssaclds(:,ncol,:) = 0.0
  asmclds(:,ncol,:) = 0.0
  fsfclds(:,ncol,:) = asmclds(:,ncol,:) * asmclds(:,ncol,:)
  cicewp   (ncol,:) = 0.0
  cliqwp   (ncol,:) = 0.0
  reice    (ncol,:) = 0.0
  reliq    (ncol,:) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! No Aerosols for now

  if (iswrtyp == 2 .and. cosz(iw) > 0.03) then

     icloud = icld

!!     swuflx  = 0.0
!!     swdflx  = 0.0
!!     swhr    = 0.0
!!     swuflxc = 0.0
!!     swdflxc = 0.0
!!     swhrc   = 0.0


     call rrtmg_sw_nomcica( ncol   , nrad   , icloud ,                          &
                            play   , plev   , tlay   , tlev   , tsfc   ,        &
                            h2ovmr , o3vmr  , co2vmr , ch4vmr , n2ovmr , o2vmr, &
                            asdir  , asdif  , aldif  , aldif  ,                 &
                            coszen , solfac , dyofyr , solar  ,                 &   
                            inflg  , iceflg , liqflg , cldfr  ,                 &
                            tauclds, ssaclds, asmclds, fsfclds,                 &
                            cicewp , cliqwp , reice  , reliq  ,                 &
                            tauaers, ssaaers, asmaers, ecaer  ,                 &
                            swuflx , swdflx , swhr   , swuflxc, swdflxc, swhrc  )

     rshort        (iw) = swdflx(1,1)
!!   rshort_diffuse(iw) = flx_diff
     rshort_top    (iw) = swdflx(1,nrad)
     rshortup_top  (iw) = swuflx(1,nrad)
     albedt        (iw) = swuflx(1,1) / swdflx(1,1)

     do k = ka, mza-1
        krad = k - koff
        fthrd_sw(k,iw) = swhr(1,krad) * exl(krad) / 86400.0
     enddo

  endif

  if (ilwrtyp == 2) then

     icloud = 1

     call rrtmg_lw_nomcica( ncol    , nrad    , icloud  , idrv    ,                 &
                            play    , plev    , tlay    , tlev    , tsfc  ,         &
                            h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr, o2vmr,  &
                            cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis  ,         &
                            inflg   , iceflg  , liqflg  , cldfr   ,                 &
                            taucldl , cicewp  , cliqwp  , reice   , reliq ,         &
                            tauaerl , &
                            lwuflx  , lwdflx  , lwhr    , lwuflxc , lwdflxc, lwhrc, &
                            duflx_dt, duflxc_dt                                     )

     rlong(iw)       = lwdflx(1,1)
     rlongup(iw)     = lwuflx(1,1)
     rlongup_top(iw) = lwuflx(1,nrad)

     do k = ka, mza-1
        krad = k - koff
        fthrd_lw(k,iw) = lwhr(1,krad) * exl(krad) / 86400.0
     enddo

  endif

end subroutine rrtmg_raddriv
