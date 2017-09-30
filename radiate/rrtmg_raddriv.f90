subroutine rrtmg_raddriv(iw, ka, nrad, koff, nsfc, &
                         rlongup_ks, rlong_albedo_ks, albedt_ks, albedt_diffuse_ks )

  use mem_grid,    only: mza, zm, zt, glatw, dzt, dzit
  use mem_basic,   only: rho, press, theta, tair, sh_v
  use misc_coms,   only: iswrtyp, ilwrtyp, time8, dtlm, io6
  use consts_coms, only: stefan, eps_virt, eps_vapi, grav, solar, cp, pi1, &
                         t00, r8, cpi
  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         albedt_beam, albedt_diffuse, rshort_diffuse, &
                         rlong_albedo, solfac, cloud_frac, rshort_clr, &
                         rshortup_clr, rshort_top_clr, rshortup_top_clr, &
                         rlong_clr, rlongup_clr, rlongup_top_clr, &
                         par, par_diffuse, uva, uvb, uvc, pbl_cld_forc, &
                         ppfd, ppfd_diffuse, mcica_seed, &
                         rlong_ks, rshort_ks, rshort_diffuse_ks, &
                         ppfd_ks, ppfd_diffuse_ks
  use micro_coms,  only: ncat, rxmin, emb0, reffcof, pwmasi, dmncof, jhabtab, &
                         emb2, zfactor_ccn
  use mem_cuparm,  only: kcutop, kcubot, qwcon, conprr, iactcu
  use rrtmg_cloud, only: cloud_props
  use mem_turb,    only: frac_land, kpblh, frac_sfc, frac_sfck
  use mem_mclat,   only: rad_mclat
  use mem_ijtabs,  only: itab_w
  use mem_micro,   only: cldnum
  use therm_lib,   only: rhovsl

  use parrrtm,             only: nbndlw, ngptlw
  use parrrsw,             only: nbndsw, ngptsw
  use rrtmg_sw_rad,        only: rrtmg_sw
  use rrtmg_lw_rad,        only: rrtmg_lw
  use mcica_subcol_gen_sw, only: mcica_subcol_sw
  use mcica_subcol_gen_lw, only: mcica_subcol_lw
  use rrsw_wvn,            only: ngb_sw => ngb
  use rrlw_wvn,            only: ngb_lw => ngb

  implicit none

  integer, intent(in) :: iw
  integer, intent(in) :: ka
  integer, intent(in) :: nrad 
  integer, intent(in) :: koff
  integer, intent(in) :: nsfc
 
  real, intent(in) :: rlongup_ks(nsfc)
  real, intent(in) :: rlong_albedo_ks(nsfc)
  real, intent(in) :: albedt_ks(nsfc)
  real, intent(in) :: albedt_diffuse_ks(nsfc)

  integer, parameter :: ncol   = 1
  integer, parameter :: icld   = 2
  integer, parameter :: iaer   = 0
  integer, parameter :: dyofyr = 0
  integer, parameter :: inflg  = 0
  integer, parameter :: iceflg = 0
  integer, parameter :: liqflg = 0

  integer, parameter :: idrv   = 0

  real, parameter :: rliqland  = 14.0
  real, parameter :: rliqocean =  8.0

  integer :: jhcat(mza,ncat)  ! hydrom category table with ice habits

  real :: rhov (mza)       ! vapor density [kg_vap/m^3]
  real :: rx   (mza,ncat)  ! hydrom bulk spec dens [kg_hyd/kg_air]
  real :: cx   (mza,ncat)  ! hydrom bulk number [num_hyd/kg_air]
  real :: emb  (mza,ncat)  ! hydrom mean particle mass [kg/particle]

  real :: coszen(ncol)

  real :: par_k(nsfc)
  real :: par_diffuse_k(nsfc)
  real :: ppfd_k(nsfc)
  real :: ppfd_diffuse_k(nsfc)
  
  real :: uva_k(nsfc)
  real :: uvb_k(nsfc)
  real :: uvc_k(nsfc)

  real :: asdir (nsfc)
  real :: aldir (nsfc)
  real :: asdif (nsfc)
  real :: aldif (nsfc)

  real :: tsfc  (nsfc)
  real :: emis  (nsfc, nbndlw)

  real :: zsfc
  real :: emiss
  real :: p1, p2, tc, rh
  real :: fland

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

  real :: cldfmcl(ngptsw, ncol, nrad)
  real :: clwpmcl(ngptsw, ncol, nrad)
  real :: ciwpmcl(ngptsw, ncol, nrad)

  real :: relqmcl(ncol, nrad)
  real :: reicmcl(ncol, nrad)

  real :: taucmcl(ngptsw, ncol, nrad)
  real :: ssacmcl(ngptsw, ncol, nrad)
  real :: asmcmcl(ngptsw, ncol, nrad)
  real :: fsfcmcl(ngptsw, ncol, nrad)
  
  real :: taucmcl_lw(ngptlw, ncol, nrad)
  real :: cldfmcl_lw(ngptlw, ncol, nrad)
  real :: clwpmcl_lw(ngptlw, ncol, nrad)
  real :: ciwpmcl_lw(ngptlw, ncol, nrad)

  real :: dl (nrad)
  real :: zml(nrad)
  real :: ztl(nrad)
! real :: exl(nrad)
  real :: dzl(nrad)

  real :: swuflxt(nrad+1)
  real :: swdflxt(nrad+1)
  real :: swuflxc(nrad+1)
  real :: swdflxc(nrad+1)

  real :: swuflxt_sfc    (nsfc)
  real :: swdflxt_sfc    (nsfc)
  real :: swdflxt_dir_sfc(nsfc)

  real :: swuflxc_sfc    (nsfc)
  real :: swdflxc_sfc    (nsfc)
  real :: swdflxc_dir_sfc(nsfc)

  real :: swuflxt_band_sfc    (nsfc,nbndsw)
  real :: swdflxt_band_sfc    (nsfc,nbndsw)
  real :: swdflxt_band_dir_sfc(nsfc,nbndsw)

  real :: swuflxc_band_sfc    (nsfc,nbndsw)
  real :: swdflxc_band_sfc    (nsfc,nbndsw)
  real :: swdflxc_band_dir_sfc(nsfc,nbndsw)

  real :: lwuflxt(nrad+1) 
  real :: lwdflxt(nrad+1)
  real :: lwuflxc(nrad+1)
  real :: lwdflxc(nrad+1)

  real :: lwuflxt_sfc(nsfc)
  real :: lwdflxt_sfc(nsfc)
  real :: lwuflxc_sfc(nsfc)
  real :: lwdflxc_sfc(nsfc)

  real :: tauaerl(ncol, nrad, nbndlw)
  real :: tauaers(ncol, nrad, nbndsw)
  real :: ssaaers(ncol, nrad, nbndsw)
  real :: asmaers(ncol, nrad, nbndsw)
  real :: ecaer  (ncol, nrad, nbndsw)

  integer :: k, ks, krad, icloud, iaeros, ib, ig, krad1, krad2
  integer :: iplon, irng, permuteseed, ns, nt, iseed
  integer :: mc, mcat, ih, l, ntim, ngbmsw, ngbmlw

  real :: r_ef, dmean, watp, twc, prate, rshort_dir, rshortup

  real :: flux_net(mza), flux_net_bot
  real :: frac(mza)
  logical :: dosnow

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

! Array kradcat maps RAMS/OLAM microphysics hydrometeor categories to those
! represented in the cloud optics tables according to the following numbering:

!     Lookup table category             OLAM Microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   drizzle                     2.  rain
!  3:   rain                        3.  pristine ice columns
!  4:   hail                        4.  snow columns
!  5:   aggregates                  5.  aggregates
!  6:   hollow columns              6.  graupel
!  7:   solid columns               7.  hail
!  8:   hexagonal plates            8.  drizzle
!  9:   rosettes                    9.  pristine ice hexagonal plates
!                                  10.  pristine ice dendrites
!                                  11.  pristine ice needles
!                                  12.  pristine ice rosettes
!                                  13.  snow hexagonal plates
!                                  14.  snow dendrites
!                                  15.  snow needles
!                                  16.  snow rosettes

                                      ! 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
  integer, parameter :: kradcat(16) = (/1,3,6,6,5,4,4,2,8, 8, 7, 9, 8, 8, 7, 9/)

! Set some surface values needed by RRTMg

  do ks = 1, nsfc
     emiss      = 1.0 - rlong_albedo_ks(ks)
     emis(ks,:) = emiss
     tsfc(ks)   = (rlongup_ks(ks) / emiss / stefan) ** 0.25

     asdir(ks) = albedt_ks(ks)
     aldir(ks) = albedt_ks(ks)
     asdif(ks) = albedt_diffuse_ks(ks)
     aldif(ks) = albedt_diffuse_ks(ks)
  enddo

  zsfc  = zm(ka-1)


  coszen(ncol) = cosz(iw)

  if (allocated(frac_land)) then
     fland = frac_land(iw)
  else
     fland = 1.0
  endif

! Copy column values from model to radiation memory space

  do k = ka, mza
     krad = k - koff

     rhov(k) = max(0.,sh_v(k,iw)) * rho(k,iw)

     play  (1,krad) = press(k,iw)
     tlay  (1,krad) = tair (k,iw)
     dl      (krad) = rho  (k,iw)
!    exl     (krad) = theta(k,iw) / tair(k,iw)
     h2ovmr(1,krad) = rhov (k)
     zml     (krad) = zm   (k)
     ztl     (krad) = zt   (k)
     dzl     (krad) = dzt  (k)
  enddo

! Fill ozone column and any extra radiation layers at model top

  call rad_mclat(iw,nrad,koff,glatw(iw),dl,play,h2ovmr,tlay,o3vmr,zml,ztl,dzl)

! Gases not defined in OLAM

  do krad = 1, nrad
     co2vmr  (1,krad) = co2
     o2vmr   (1,krad) = o2
     ch4vmr  (1,krad) = ch4
     n2ovmr  (1,krad) = n2o
     cfc11vmr(1,krad) = cfc11
     cfc12vmr(1,krad) = cfc12
     cfc22vmr(1,krad) = cfc22
     ccl4vmr (1,krad) = ccl4
  enddo

! Convert water vapor from density to molar mixing ratio, pressure to mb,
! and ozone from density to molar mixing ratio

  do krad = 1, nrad
     h2ovmr(1,krad) = h2ovmr(1,krad) * eps_vapi / dl(krad)
     play  (1,krad) = play  (1,krad) * 0.01
     o3vmr (1,krad) = o3vmr (1,krad) * amdo3 / dl(krad)
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

! initialize cloud properties to 0

  cldfr    (ncol,:) = 0.0
  taucldl(:,ncol,:) = 0.0
  tauclds(:,ncol,:) = 0.0
  ssaclds(:,ncol,:) = 0.0
  asmclds(:,ncol,:) = 0.0
  fsfclds(:,ncol,:) = 0.0
  cicewp   (ncol,:) = 0.0  ! not used when we specify cloud optical properties (iceflg=0)
  cliqwp   (ncol,:) = 0.0  ! not used when we specify cloud optical properties (liqflg=0)
  reice    (ncol,:) = 0.0  ! not used when we specify cloud optical properties (iceflg=0)
  reliq    (ncol,:) = 0.0  ! not used when we specify cloud optical properties (liqflg=0)

! Fill arrays rx, cx, and emb with hydrometeor properties

  call cloudprep_rad(iw,ka,mcat,jhcat,rhov,rx,cx,emb)

! Compute fractional cloudiness for the resolved microphysics moisture fields. 
! The cloud fraction estimated here will only be applied later in this routine
! if there are any resolved hydrometeors or subgrid cumulus

  call get_cloud_frac(iw, ka, frac)

! Get optical properties of resolved clouds

  do mc = 1, mcat

     do k = ka, mza
        krad = k - koff

        if (rx(k,mc) >= rxmin(mc) .and. emb(k,mc) >= emb0(mc)) then

           ! Set lower bound on frac(k) because there is condensate
           frac(k)       = max(frac(k),0.1)
           cldfr(1,krad) = frac(k)

           ! lookup table category
           ih = jhcat(k,mc)
           l  = kradcat(ih)

           ! effective radius in microns
           r_ef = 1.e6 * reffcof(ih) * emb(k,mc) ** pwmasi(ih)

           ! ice or liquid water path in g/m^2
           watp = rx(k,mc) * dl(krad) * 1000. * dzt(k)
           watp = watp / frac(k)

           call lookup_rrtmg_cld_optics( l, r_ef, watp, krad )

        endif
     enddo
  enddo

  ! Include the optical properties of subgrid convective clouds

  if (iactcu(iw) == 1) then

     do k = kcubot(iw), kcutop(iw)
        krad = k - koff

        ! Set lower bound on frac(k) because there is condensate
        frac(k)       = max(frac(k), 0.1)
        cldfr(1,krad) = frac(k)

        ! water path in g/m^2
        watp = qwcon(k,iw) * dl(krad) * 1000. * dzt(k)
        watp = watp / max( frac(k), 0.2 )

        tc = tair(k,iw) - t00

        if (tc > -10.0) then

           ! Add convective cloud water to cloud drops if warmer then -10C
           ih = 1
           l  = kradcat(ih)

           ! effective radius in microns
           r_ef = 1.e6 * reffcof(ih) &
                * ( qwcon(k,iw) / (cldnum(iw) * zfactor_ccn(k)) ) ** pwmasi(ih)

        else

           ! Add convective cloud water to pristine ice. Diagnose habit
           ! from temperature and humidity
           rh = min( 1., rhov(k) / rhovsl(tc) )
           ns = max( 1, nint(100. * rh) )
           nt = max( 1, min(31,-nint(tc)) )
           ih = jhabtab(nt,ns,1)
           l  = kradcat(ih)

           ! Mean maximum dimension of ice crystals as a function of T
           ! (see Kristjansson et al., 2000, JGR)
           dmean = 1030.7 * exp(0.05522*(tair(k,iw)-279.5))

           ! Convert mean diameter to an effective radius using the 
           ! microphysics power laws
           r_ef = reffcof(ih) / dmncof(ih) * dmean

        endif

        call lookup_rrtmg_cld_optics( l, r_ef, watp, krad )

     enddo

     ! Also include the optical properties of convective rain/snow

     if (conprr(iw) > 1.e-10) then

        ! Estimate rain/snow from parameterization used in EPA's CMAQ model

        prate = conprr(iw) * 3600.0          ! precip rate: mm / hr
        twc   = 0.06 * prate**0.846          ! tot wat: g / m3

        ! if entire cloud is below freezing, map to snow
        dosnow = all(tair(kcubot(iw):kcutop(iw),iw) < 273.)

        do k = ka, kcutop(iw)
           krad = k - koff

           tc = tair(k,iw) - t00

           ! Set lower bound on frac(k) because there is condensate
           frac(k)       = max( frac(k), 0.1)
           cldfr(1,krad) = frac(k)

           watp  = twc * dzt(k)                 ! water/ice path: g / m^2
           watp  = watp / max( frac(k), 0.2 )   ! scale by cloud fraction

           if (dosnow) then

              if (tc > 0.0) then
                 ! rain
                 mc = 2
                 ih = 2
              else
                 ! snow
                 mc = 4
                 rh = min( 1., rhov(k) / rhovsl(tc) )
                 ns = max( 1, nint(100. * rh) )
                 nt = max( 1, min(31,-nint(tc)) )
                 ih = jhabtab(nt,ns,2)
              endif

           else

              if (tc > -10.) then
                 ! rain
                 mc = 2
                 ih = 2
              else
                 ! hail
                 mc = 7
                 ih = 7
              endif

           endif

           ! cloud optics category
           l  = kradcat(ih)

           ! effective radius of convctive rain/snow
           r_ef = 1.e6 * reffcof(ih) * emb2(mc) ** pwmasi(ih)

           call lookup_rrtmg_cld_optics( l, r_ef, watp, krad )

        enddo

     endif  ! convective rain/snow

  endif  ! convective clouds

  ! Save cloud fraction in 3D variable for output or plotting

  cloud_frac(1:ka-1,iw) = 0.0

  do k = ka, mza
     krad = k - koff
     if (cldfr(1,krad) > 0.99) cldfr(1,krad) = 1.0
     if (cldfr(1,krad) < 0.01) cldfr(1,krad) = 0.0
     cloud_frac(k,iw) = cldfr(1,krad)
  enddo

  ! Compute the shortwave fluxes and heating rates

  if (iswrtyp > 0 .and. cosz(iw) > 0.03) then

     ! Combine optical properties

     do k = ka, mza
        krad = k - koff
        do ib = 1, nbndsw
           if (tauclds(ib,1,krad) > 1.e-12 .and. ssaclds(ib,1,krad) > 1.e-12) then
              asmclds(ib,1,krad) = asmclds(ib,1,krad) / ssaclds(ib,1,krad)
              ssaclds(ib,1,krad) = ssaclds(ib,1,krad) / tauclds(ib,1,krad)
              fsfclds(ib,1,krad) = asmclds(ib,1,krad) * asmclds(ib,1,krad)
           endif
        enddo
     enddo

     icloud = icld
     iaeros = iaer
     irng = 0
     permuteseed = 1 + nint( mod(time8/dtlm(1), 1.e8_r8) )

     if ( any( cloud_frac(ka:mza,iw) > 0.01 .and. cloud_frac(ka:mza,iw) < 0.99 ) ) then
     
        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_sw(nrad, icloud, permuteseed, mcica_seed(1:4,iw), irng, &
                             cldfr,  cicewp,  cliqwp, reice, reliq, &
                             tauclds, ssaclds, asmclds, fsfclds, &
                             cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                             taucmcl, ssacmcl, asmcmcl, fsfcmcl, inflg)
     else

        ! No fractional cloudiness (cloud fraction either 0 or 1)

        ngbmsw = ngb_sw(1) - 1

        do k = 1, nrad

           if (cldfr(1,k) < 0.5) then

              relqmcl(1,k) = 0.0
              reicmcl(1,k) = 0.0
              
              do ig = 1, ngptsw
                 cldfmcl(ig,1,k) = 0.0
                 clwpmcl(ig,1,k) = 0.0
                 ciwpmcl(ig,1,k) = 0.0
                 taucmcl(ig,1,k) = 0.0
                 ssacmcl(ig,1,k) = 1.0
                 asmcmcl(ig,1,k) = 0.0
                 fsfcmcl(ig,1,k) = 0.0
              enddo

           else

              relqmcl(1,k) = reliq(1,k)
              reicmcl(1,k) = reice(1,k)
              
              do ig = 1, ngptsw
                 ib = ngb_sw(ig) - ngbmsw

                 cldfmcl(ig,1,k) = 1.0
                 clwpmcl(ig,1,k) = cliqwp    (1,k)
                 ciwpmcl(ig,1,k) = cicewp    (1,k)
                 taucmcl(ig,1,k) = tauclds(ib,1,k)
                 ssacmcl(ig,1,k) = ssaclds(ib,1,k)
                 asmcmcl(ig,1,k) = asmclds(ib,1,k)
                 fsfcmcl(ig,1,k) = fsfclds(ib,1,k)
              enddo
              
           endif

        enddo

     endif

     call rrtmg_sw(ncol    ,nrad    ,icloud  ,iaeros  ,nsfc   ,frac_sfck(1:nsfc,iw), &
                   play    ,plev    ,tlay    ,tlev    ,tsfc   ,         &
                   h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr   ,&
                   asdir   ,asdif   ,aldir   ,aldif   ,                 &
                   coszen  ,solfac  ,dyofyr  ,solar   ,                 &
                   inflg   ,iceflg  ,liqflg  ,cldfmcl ,                 &
                   taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl ,                 &
                   ciwpmcl ,clwpmcl ,reicmcl ,relqmcl ,                 &
                   tauaers ,ssaaers ,asmaers ,ecaer   ,                 &
                   swuflxt ,swdflxt ,swuflxc ,swdflxc ,                 &
                   swuflxt_sfc, swdflxt_sfc,  swdflxt_dir_sfc,          &
                   swuflxc_sfc, swdflxc_sfc,  swdflxc_dir_sfc,          &
                   swuflxt_band_sfc, swdflxt_band_sfc, swdflxt_band_dir_sfc, &
                   swuflxc_band_sfc, swdflxc_band_sfc, swdflxc_band_dir_sfc )

     do krad = 1, mza - koff + 1
        flux_net(krad) = swuflxt(krad) - swdflxt(krad)
     enddo

     do k = ka, ka + nsfc - 1
        krad = k - koff
        flux_net_bot = (swuflxt_sfc(krad) - swdflxt_sfc(krad)) *        frac_sfck(krad,iw) &
                     + (swuflxt    (krad) - swdflxt    (krad)) * (1.0 - frac_sfck(krad,iw))
        fthrd_sw(k,iw) = (flux_net_bot - flux_net(krad+1)) * dzit(k) * cpi
     enddo

     do k = ka + nsfc, mza
        krad = k - koff
        fthrd_sw(k,iw) = (flux_net(krad) - flux_net(krad+1)) * dzit(k) * cpi
     enddo

     rshort        (iw) = surface_avg( swdflxt_sfc )
     rshortup           = surface_avg( swuflxt_sfc )
     rshort_dir         = surface_avg( swdflxt_dir_sfc )
     rshort_diffuse(iw) = rshort(iw) - rshort_dir
     albedt        (iw) = rshortup / rshort(iw)
     rshort_top    (iw) = swdflxt(nrad+1)
     rshortup_top  (iw) = swuflxt(nrad+1)

     rshort_clr      (iw) = surface_avg( swdflxc_sfc )
     rshortup_clr    (iw) = surface_avg( swuflxc_sfc )
     rshort_top_clr  (iw) = swdflxc(nrad+1)
     rshortup_top_clr(iw) = swuflxc(nrad+1)

     do krad = 1, nsfc

        par_k(krad) = 0.5268*swdflxt_band_sfc(krad, 9) &
                    +        swdflxt_band_sfc(krad,10) &
                    + 0.4724*swdflxt_band_sfc(krad,11)

        par_diffuse_k(krad) = par_k(krad) - ( 0.5268*swdflxt_band_dir_sfc(krad, 9) &
                                            +        swdflxt_band_dir_sfc(krad,10) &
                                            + 0.4724*swdflxt_band_dir_sfc(krad,11) )

        ppfd_k(krad) = 2.9156 * swdflxt_band_sfc(krad, 9) &
                     + 4.4568 * swdflxt_band_sfc(krad,10) &
                     + 1.6614 * swdflxt_band_sfc(krad,11)

        ppfd_diffuse_k(krad) = ppfd_k(krad) - ( 2.9156 * swdflxt_band_dir_sfc(krad, 9) &
                                              + 4.4568 * swdflxt_band_dir_sfc(krad,10) &
                                              + 1.6614 * swdflxt_band_dir_sfc(krad,11) )

        uva_k(krad) = 0.5276*swdflxt_band_sfc(krad,11) + 0.3932*swdflxt_band_sfc(krad,12)

        uvb_k(krad) = 0.3708*swdflxt_band_sfc(krad,12)

        uvc_k(krad) = 0.2360*swdflxt_band_sfc(krad,12) + swdflxt_band_sfc(krad,13)

     enddo

     par         (iw) = surface_avg( par_k )
     par_diffuse (iw) = surface_avg( par_diffuse_k )
     ppfd        (iw) = surface_avg( ppfd_k )
     ppfd_diffuse(iw) = surface_avg( ppfd_diffuse_k )
     uva         (iw) = surface_avg( uva_k )
     uvb         (iw) = surface_avg( uvb_k )
     uvc         (iw) = surface_avg( uvc_k )

     do krad = 1, nsfc

        rshort_ks        (krad,iw) = swdflxt_sfc(krad)
        rshort_diffuse_ks(krad,iw) = swdflxt_sfc(krad) - swdflxt_dir_sfc(krad)

!       par_ks         (krad,iw) = par_k(krad)
!       par_diffuse_ks (krad,iw) = par_diffuse_k(krad)
        ppfd_ks        (krad,iw) = ppfd_k(krad)
        ppfd_diffuse_ks(krad,iw) = ppfd_diffuse_k(krad)
        
     enddo

     krad1 = max(kpblh(iw) - 1, ka)
     krad2 = min(kpblh(iw) + 1, mza)

     pbl_cld_forc(iw) = - sum( fthrd_sw(krad1:krad2,iw) * dzt(krad1:krad2) ) / rho(kpblh(iw),iw)

  endif

  ! Compute the longwave fluxes and heating rates

  if (ilwrtyp > 0) then

     icloud = icld
     iaeros = iaer
     irng = 0
     permuteseed = 150 + nint( mod(time8/dtlm(1), 1.e8_r8) )

     if ( any( cloud_frac(ka:mza,iw) > 0.01 .and. cloud_frac(ka:mza,iw) < 0.99 ) ) then

        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_lw(nrad   , icloud    , permuteseed, mcica_seed(1:4,iw), irng , &
                             cldfr  , cicewp    , cliqwp     , reice             , reliq, &
                             taucldl, cldfmcl_lw, ciwpmcl_lw , clwpmcl_lw        ,        &
                             reicmcl, relqmcl   , taucmcl_lw , inflg                      )
     else

        ! No fractional cloudiness (cloud fraction either 0 or 1)

        ngbmlw = ngb_lw(1) - 1

        do k = 1, nrad

           if (cldfr(1,k) < 0.5) then

              relqmcl(1,k) = 0.0
              reicmcl(1,k) = 0.0
              
              do ig = 1, ngptlw
                 cldfmcl_lw(ig,1,k) = 0.0
                 clwpmcl_lw(ig,1,k) = 0.0
                 ciwpmcl_lw(ig,1,k) = 0.0
                 taucmcl_lw(ig,1,k) = 0.0
              enddo

           else

              relqmcl(1,k) = reliq(1,k)
              reicmcl(1,k) = reice(1,k)
              
              do ig = 1, ngptlw
                 ib = ngb_lw(ig) - ngbmlw

                 cldfmcl_lw(ig,1,k) = 1.0
                 clwpmcl_lw(ig,1,k) = cliqwp    (1,k)
                 ciwpmcl_lw(ig,1,k) = cicewp    (1,k)
                 taucmcl_lw(ig,1,k) = taucldl(ib,1,k)
              enddo

           endif

        enddo

     endif

     call rrtmg_lw(ncol       ,nrad       ,icloud     ,nsfc       ,iaeros  ,iw     , &
                   play       ,plev       ,tlay       ,tlev       ,tsfc    ,frac_sfck(1:nsfc,iw),&
                   h2ovmr     ,o3vmr      ,co2vmr     ,ch4vmr     ,n2ovmr  ,o2vmr  , &
                   cfc11vmr   ,cfc12vmr   ,cfc22vmr   ,ccl4vmr    ,emis    ,         &
                   inflg      ,iceflg     ,liqflg     ,cldfmcl_lw ,                  &
                   taucmcl_lw ,ciwpmcl_lw ,clwpmcl_lw ,reicmcl    ,relqmcl ,tauaerl, &
                   lwuflxt    ,lwdflxt    ,lwuflxc    ,lwdflxc    ,                  &
                   lwuflxt_sfc,lwdflxt_sfc,lwuflxc_sfc,lwdflxc_sfc                   )

     do krad = 1, mza - koff + 1
        flux_net(krad) = lwuflxt(krad) - lwdflxt(krad)
     enddo

     do k = ka, ka + nsfc - 1
        krad = k - koff
        flux_net_bot = (lwuflxt_sfc(krad) - lwdflxt_sfc(krad)) *        frac_sfck(krad,iw) &
                     + (lwuflxt    (krad) - lwdflxt    (krad)) * (1.0 - frac_sfck(krad,iw))
        fthrd_lw(k,iw) = (flux_net_bot - flux_net(krad+1)) * dzit(k) * cpi
     enddo

     do k = ka + nsfc, mza
        krad = k - koff
        fthrd_lw(k,iw) = (flux_net(krad) - flux_net(krad+1)) * dzit(k) * cpi
     enddo

     rlong      (iw) = surface_avg( lwdflxt_sfc )
     rlongup    (iw) = surface_avg( lwuflxt_sfc )
     rlongup_top(iw) = lwuflxt(nrad+1)

     rlong_clr      (iw) = surface_avg( lwdflxc_sfc )
     rlongup_clr    (iw) = surface_avg( lwuflxt_sfc )
     rlongup_top_clr(iw) = lwuflxc(nrad+1)

     do krad = 1, nsfc
        rlong_ks(krad,iw) = lwdflxt_sfc(krad)
     enddo

     krad1 = max(kpblh(iw) - 1, ka)
     krad2 = min(kpblh(iw) + 1, mza)

     pbl_cld_forc(iw) = pbl_cld_forc(iw) &
                      - sum( fthrd_lw(krad1:krad2,iw) * dzt(krad1:krad2) ) / rho(kpblh(iw),iw)

  endif


contains


  real function surface_avg( field )

    implicit none
    real, intent(in) :: field(nsfc)

    if (nsfc == 1) then
       surface_avg = field(1)
    else
       surface_avg = dot_product( field(1:nsfc), frac_sfc(1:nsfc,iw) )
    endif

  end function surface_avg


  subroutine lookup_rrtmg_cld_optics( l, r_ef, watp, krad )

    implicit none

    integer, intent(in) :: l     ! lookup table category
    real,    intent(in) :: r_ef  ! effective radius (um)
    real,    intent(in) :: watp  ! ice or liquid water path (g/m^2)
    integer, intent(in) :: krad

    integer :: num, index
    real    :: rstart, rend, rscale, reff
    real    :: fint0, fint1
    real    :: tau, ssa, asm

    num    = cloud_props(l)%num
    rstart = cloud_props(l)%start
    rend   = cloud_props(l)%end

    reff   = max(rstart, min(rend, r_ef))

    rscale = (reff - rstart) / cloud_props(l)%delr
    index  = max(1, min(num-1, int(rscale) + 1))
    fint1  = rscale - real(index-1)
    fint0  = 1.0 - fint1

    if (iswrtyp > 0 .and. cosz(iw) >= 0.03) then

       do ib = 1, nbndsw
          tau = ( fint0 * cloud_props(l)%extsw(ib, index  ) &
                + fint1 * cloud_props(l)%extsw(ib, index+1) ) * watp

          ssa = fint0 * cloud_props(l)%ssasw(ib, index  ) &
              + fint1 * cloud_props(l)%ssasw(ib, index+1)

          asm = fint0 * cloud_props(l)%asysw(ib, index  ) &
              + fint1 * cloud_props(l)%asysw(ib, index+1)

          tauclds(ib,1,krad) = tauclds(ib,1,krad) + tau
          ssaclds(ib,1,krad) = ssaclds(ib,1,krad) + tau * ssa
          asmclds(ib,1,krad) = asmclds(ib,1,krad) + tau * ssa * asm

       enddo

    endif

    if (ilwrtyp > 0) then

       do ib = 1, nbndlw
          tau = ( fint0 * cloud_props(l)%abslw(ib, index  ) &
                + fint1 * cloud_props(l)%abslw(ib, index+1) ) * watp

          taucldl(ib,1,krad) = taucldl(ib,1,krad) + tau
       enddo

    endif

  end subroutine lookup_rrtmg_cld_optics


end subroutine rrtmg_raddriv
