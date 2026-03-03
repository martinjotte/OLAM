module rrtmg_driver

contains

!==============================================================================

subroutine rrtmg_raddriv(iw, ka, nrad, koff, nsfc, &
                         rlongup_ks, rlong_albedo_ks, albedt_ks, albedt_diffuse_ks )

  use mem_grid,    only: mza, zm, zt, glatw, dzt, dzit
  use mem_basic,   only: rho, press, tair, rr_v
  use misc_coms,   only: iswrtyp, ilwrtyp, i_o3, i_co, i_ch4
  use consts_coms, only: stefan, eps_virt, eps_vapi, grav, solar, cp, pi1, &
                         t00, r8, cpi
  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz_rad, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         rshort_diffuse, solfac, cloud_frac, &
                         par, par_diffuse, uva, uvb, uvc, pbl_cld_forc, &
                         ppfd, ppfd_diffuse, mcica_seed, cosz_min, &
                         rlong_ks, rshort_ks, rshort_diffuse_ks, &
                         ppfd_ks, ppfd_diffuse_ks, dlong, ulong, ulong_sfc, &
                         Dulong_DT, Dulong_sfc_DT, Tskin_rad
  use micro_coms,  only: ncat, reffcof, pwmasi
  use rrtmg_cloud, only: cloud_props
  use mem_turb,    only: frac_land, kpblh, frac_sfc, frac_sfck
  use mem_mclat,   only: rad_mclat
  use therm_lib,   only: rhovsl
  use mem_co2,     only: rr_co2, i_co2, co2_initppm, co2_sh2ppm
  use var_tables,  only: scalar_tab

  use parkind,             only: cldmin, cldmax
  use parrrtm,             only: nbndlw, ngptlw
  use parrrsw,             only: nbndsw, ngptsw, naerec
  use rrtmg_sw_rad,        only: rrtmg_sw
  use rrtmg_lw_rad,        only: rrtmg_lw
  use mcica_subcol_gen_sw, only: mcica_subcol_sw
  use mcica_subcol_gen_lw, only: mcica_subcol_lw
  use rrsw_wvn,            only: nga_sw => nga, ngs_sw => ngs
  use rrlw_wvn,            only: nga_lw => nga, ngs_lw => ngs

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

  integer, parameter :: icld   = 2
  integer, parameter :: iaer   = 0
  integer, parameter :: dyofyr = 0
  integer, parameter :: inflg  = 0
  integer, parameter :: iceflg = 0
  integer, parameter :: liqflg = 0

  integer, parameter :: idrv   = 0

  real, parameter :: rliqland  = 14.0
  real, parameter :: rliqocean =  8.0

  integer :: jhcat(mza,ncat+1)  ! hydrom category table with ice habits

  real :: rhov (mza)       ! vapor density [kg_vap/m^3]
  real :: rx   (mza,ncat+1)  ! hydrom bulk spec dens [kg_hyd/kg_air]
  real :: emb  (mza,ncat+1)  ! hydrom mean particle mass [kg/particle]

  integer :: ktop(ncat+1)

  real :: coszen

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
  real :: emis  (nbndlw, nsfc)
  real :: emiss (nsfc)
  real :: fland

  real :: coldry  (nrad)  ! molecules dry air / cm^2
  real :: play    (nrad)
  real :: tlay    (nrad)
  real :: cldfr   (nrad)
  real :: cicewp  (nrad)
  real :: cliqwp  (nrad)
  real :: reice   (nrad)
  real :: reliq   (nrad)

  ! volume mixing ratios (moles / moles of dry air)
  real :: h2ovmr  (nrad)
  real :: co2vmr  (nrad)
  real :: o3vmr   (nrad)
  real :: o2vmr   (nrad)
  real :: ch4vmr  (nrad)
  real :: n2ovmr  (nrad)
  real :: covmr   (nrad)
  real :: cfc11vmr(nrad)
  real :: cfc12vmr(nrad)
  real :: cfc22vmr(nrad)
  real :: ccl4vmr (nrad)

  real :: taucldl(nbndlw, nrad)

  real :: tauclds(nbndsw, nrad)
  real :: ssaclds(nbndsw, nrad)
  real :: asmclds(nbndsw, nrad)

  real :: clwpmcs(ngptsw, nrad)
  real :: ciwpmcs(ngptsw, nrad)

  real :: taucmcs(ngptsw, nrad)
  real :: ssacmcs(ngptsw, nrad)
  real :: asmcmcs(ngptsw, nrad)

  real :: taucmcl(ngptlw, nrad)
  real :: clwpmcl(ngptlw, nrad)
  real :: ciwpmcl(ngptlw, nrad)

  real :: dl_wet(nrad)
  real :: dl_dry(nrad)
  real :: dl_vap(nrad)
  real :: zml   (nrad)
  real :: ztl   (nrad)
  real :: dzl   (nrad)

  real :: swuflxt(nrad)
  real :: swdflxt(nrad)
! real :: swuflxc(nrad)
! real :: swdflxc(nrad)

  real :: swuflxt_sfc    (nsfc)
  real :: swdflxt_sfc    (nsfc)
  real :: swdflxt_dir_sfc(nsfc)

! real :: swuflxc_sfc    (nsfc)
! real :: swdflxc_sfc    (nsfc)
! real :: swdflxc_dir_sfc(nsfc)

  real :: swuflxt_band_sfc    (nsfc,nbndsw)
  real :: swdflxt_band_sfc    (nsfc,nbndsw)
  real :: swdflxt_band_dir_sfc(nsfc,nbndsw)

! real :: swuflxc_band_sfc    (nsfc,nbndsw)
! real :: swdflxc_band_sfc    (nsfc,nbndsw)
! real :: swdflxc_band_dir_sfc(nsfc,nbndsw)

  real :: lwuflxt(nrad+1)
  real :: lwdflxt(nrad+1)
! real :: lwuflxc(nrad)
! real :: lwdflxc(nrad)

  real :: lwuflxt_sfc(nsfc)
! real :: lwdflxt_sfc(nsfc)
! real :: lwuflxc_sfc(nsfc)
! real :: lwdflxc_sfc(nsfc)

  real :: Dlwuflxt_DT(nrad+1)
  real :: Dlwuflxt_sfc_DT(nsfc)

  real :: tauaerl(nbndlw, nrad)
  real :: tauaers(nbndsw, nrad)
  real :: ssaaers(nbndsw, nrad)
  real :: asmaers(nbndsw, nrad)
  real :: ecaer  (naerec, nrad)

  integer :: k, ks, krad, icloud, iaeros, ib, ig, krad1, krad2
  integer :: ns, mc, mcat, ih, l, isfc_dT
  logical :: do_pblforc

  real :: r_ef, watp, rshort_dir, rshortup, pwvcm

  real :: flux_net(mza), flux_net_bot

! Set gas volume mixing ratios, 2005 values, IPCC (2007)

  real, parameter :: co2   = 379.e-6  ! carbon dioxide (379 ppmv)
  real, parameter :: ch4   = 1774.e-9 ! methane (1774 ppbv)
  real, parameter :: n2o   = 319.e-9  ! nitrous oxide (319 ppbv)
  real, parameter :: cfc11 = 0.251e-9 ! cfc-11 (251 ppt)
  real, parameter :: cfc12 = 0.538e-9 ! cfc-12 (538 ppt)
  real, parameter :: cfc22 = 0.169e-9 ! cfc-22 (169 ppt)
  real, parameter :: ccl4  = 0.093e-9 ! ccl4 (93 ppt)
  real, parameter :: o2    = 0.209488 ! oxygen, for o2mmr=0.23143
  real, parameter :: co    = 50.e-9   ! carbon monoxide (50 ppbv)

  real, parameter :: amd = 28.9660    ! Effective molecular weight of dry air (g/mol)
  real, parameter :: amw = 18.0160    ! Molecular weight of water vapor (g/mol)
  real, parameter :: amdryo3 = 0.603428
  real, parameter :: avogad = 6.022142e23 ! Avogadro constant
  real, parameter :: dn2col = avogad / amd / 10.

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
     emiss    (ks)    = 1.0 - rlong_albedo_ks(ks)
     tsfc     (ks)    = (rlongup_ks(ks) / emiss(ks) / stefan) ** 0.25
     Tskin_rad(ks,iw) = tsfc(ks)

     asdir(ks) = albedt_ks(ks)
     aldir(ks) = albedt_ks(ks)
     asdif(ks) = albedt_diffuse_ks(ks)
     aldif(ks) = albedt_diffuse_ks(ks)
  enddo

! Set emissivity constant across all bands for now

  do ks = 1, nsfc
     do ns = 1, nbndlw
        emis(ns,ks) = emiss(ks)
     enddo
  enddo

  coszen = cosz_rad(iw)

  if (allocated(frac_land)) then
     fland = frac_land(iw)
  else
     fland = 1.0
  endif

  do_pblforc = .false.
  if (kpblh(iw) > ka) then
     do_pblforc = any( cloud_frac(kpblh(iw)-1:kpblh(iw)+1,iw) > 0.05 )
  endif

! Copy column values from model to radiation memory space

  do k = ka, mza
     krad = k - koff

     dl_dry(krad) = rho(k,iw)
     rhov  (k)    = max(0., rr_v(k,iw) * dl_dry(krad))
     dl_vap(krad) = rhov(k)
     dl_wet(krad) = dl_dry(krad) + dl_vap(krad)

     play  (krad) = press(k,iw)
     tlay  (krad) = tair (k,iw)
     zml   (krad) = zm   (k)
     ztl   (krad) = zt   (k)
     dzl   (krad) = dzt  (k)
  enddo

! Fill ozone column and any extra radiation layers at model top

  call rad_mclat(iw,nrad,koff,glatw(iw),dl_wet,play,dl_vap,tlay,o3vmr,zml,ztl,dzl)

  ! dry density above model top
  do krad = mza-koff+1, nrad
     dl_dry(krad) = dl_wet(krad) - dl_vap(krad)
 enddo

! Convert water vapor and ozone from density to (dry) molar mixing ratio
! and pressure to mb.

  pwvcm = 0.0

  do krad = 1, nrad
     pwvcm        = pwvcm + dl_vap(krad) * dzl(krad)
     coldry(krad) = dl_dry(krad) * dzl(krad) * dn2col
     h2ovmr(krad) = dl_vap(krad) * eps_vapi / dl_dry(krad)
     o3vmr (krad) = o3vmr (krad) * amdryo3  / dl_dry(krad)
     play  (krad) = play  (krad) * 0.01
  enddo

  pwvcm = pwvcm * 1.e-1

! Gases not defined in OLAM

  do krad = 1, nrad
     co2vmr  (krad) = co2_initppm * 1.e-6
     o2vmr   (krad) = o2
     ch4vmr  (krad) = ch4
     n2ovmr  (krad) = n2o
     cfc11vmr(krad) = cfc11
     cfc12vmr(krad) = cfc12
     cfc22vmr(krad) = cfc22
     ccl4vmr (krad) = ccl4
     covmr   (krad) = co
  enddo

! Use CO2 if prognosed. Units are assumed to be kg/kg

  if (i_co2 > 0) then
     do k = ka, mza
        krad = k - koff
        co2vmr(krad) = rr_co2(k,iw) * co2_sh2ppm * 1.e-6
     enddo
  endif

! Use CMAQ ozone if prognosed. Units are assumed to be PPMV

  if (i_o3 > 0) then
     do k = ka, mza
        krad = k - koff
        o3vmr(krad) = scalar_tab(i_o3)%var_p(k,iw) * 1.e-6
     enddo
  endif

! Use CMAQ carbon monoxide if prognosed. Units are assumed to be PPMV

  if (i_co > 0) then
     do k = ka, mza
        krad = k - koff
        covmr(krad) = scalar_tab(i_co)%var_p(k,iw) * 1.e-6
     enddo
  endif

! Use CMAQ methane if prognosed. Units are assumed to be PPMV

  if (i_ch4 > 0) then
     do k = ka, mza
        krad = k - koff
        ch4vmr(krad) = scalar_tab(i_ch4)%var_p(k,iw) * 1.e-6
     enddo
  endif

! initialize cloud properties to 0

  cldfr = 0.0

  if (inflg == 0) then
     taucldl = 0.0
     tauclds = 0.0
     ssaclds = 0.0
     asmclds = 0.0
  else
     cicewp  = 0.0
     cliqwp  = 0.0
     reice   = 0.0
     reliq   = 0.0
  endif

! Fill arrays rx, jhcat, and emb with hydrometeor properties

  call cloudprep_rad(iw,ka,mcat,jhcat,rhov,rx,emb,ktop)

! Get optical properties of resolved and subgrid clouds

  do mc = 1, mcat

     do k = ka, ktop(mc)
        krad = k - koff

        if (rx(k,mc) > 1.e-10) then

           ! Set lower bound on cldfrac because there is condensate
           cldfr(krad) = max(cloud_frac(k,iw), 0.1)

           ! lookup table category
           ih = jhcat(k,mc)
           l  = kradcat(ih)

           ! effective radius in microns
           r_ef = 1.e6 * reffcof(ih) * emb(k,mc) ** pwmasi(ih)

           ! ice or liquid water path in g/m^2
           watp = rx(k,mc) * dl_dry(krad) * 1000. * dzt(k)
           watp = watp / cldfr(krad)

           call lookup_rrtmg_cld_optics( l, r_ef, watp, krad )

        endif
     enddo
  enddo

  ! Compute the shortwave fluxes and heating rates

  if (iswrtyp > 0 .and. coszen > cosz_min) then

     icloud = icld
     iaeros = iaer

     if ( any( cldfr(1:mza-koff) >= cldmin .and. cldfr(1:mza-koff) <= cldmax ) ) then

        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_sw(nrad,    icloud,  mcica_seed(1:4,iw),        &
                             cldfr,   cicewp,  cliqwp,  &
                             tauclds, ssaclds, asmclds,          &
                             ciwpmcs, clwpmcs, &
                             taucmcs, ssacmcs, asmcmcs, inflg    )
     else

        ! No fractional cloudiness (cloud fraction either 0 or 1)

        if (inflg == 0) then

           do k = 1, mza-koff
              if (cldfr(k) > 0.5) then
                 do ib = 1, nbndsw
                    do ig = nga_sw(ib), ngs_sw(ib)
                       taucmcs(ig,k) = tauclds(ib,k)
                       ssacmcs(ig,k) = ssaclds(ib,k)
                       asmcmcs(ig,k) = asmclds(ib,k)
                    enddo
                 enddo
              endif
           enddo

        else

           do k = 1, mza-koff
              if (cldfr(k) > 0.5) then
                 do ig = 1, ngptsw
                    clwpmcs(ig,k) = cliqwp(k)
                    ciwpmcs(ig,k) = cicewp(k)
                 enddo
              endif
           enddo

        endif

     endif

     call rrtmg_sw(nrad    ,iaeros  ,nsfc    ,frac_sfck(1:nsfc,iw),     &
                   play    ,tlay    ,coldry  , &
                   h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr,           &
                   asdir   ,asdif   ,aldir   ,aldif   ,                 &
                   coszen  ,solfac  ,dyofyr  ,solar   ,                 &
                   inflg   ,iceflg  ,liqflg  ,cldfr   ,                 &
                   taucmcs ,ssacmcs ,asmcmcs ,                          &
                   ciwpmcs ,clwpmcs ,reice   ,reliq   ,                 &
                   tauaers ,ssaaers ,asmaers ,ecaer   ,                 &
                   swuflxt ,swdflxt ,                                   &
                   swuflxt_sfc, swdflxt_sfc,  swdflxt_dir_sfc,          &
                   swuflxt_band_sfc, swdflxt_band_sfc, swdflxt_band_dir_sfc)

     do krad = 1, mza - koff
        k = krad + koff
        flux_net(k) = swuflxt(krad) - swdflxt(krad)
     enddo
     flux_net(ka-1) = 0.0

     do k = ka, ka + nsfc - 1
        krad = k - koff
        flux_net_bot =        frac_sfck(krad,iw)  * (swuflxt_sfc(krad) - swdflxt_sfc(krad)) &
                     + (1.0 - frac_sfck(krad,iw)) * flux_net(k-1)
        fthrd_sw(k,iw) = (flux_net_bot - flux_net(k)) * dzit(k) * cpi
     enddo

     do k = ka + nsfc, mza
        krad = k - koff
        fthrd_sw(k,iw) = (flux_net(k-1) - flux_net(k)) * dzit(k) * cpi
     enddo

     rshort        (iw) = surface_avg( swdflxt_sfc )
     rshortup           = surface_avg( swuflxt_sfc )
     rshort_dir         = surface_avg( swdflxt_dir_sfc )
     rshort_diffuse(iw) = rshort(iw) - rshort_dir
     albedt        (iw) = rshortup / rshort(iw)
     rshort_top    (iw) = swdflxt(nrad)
     rshortup_top  (iw) = swuflxt(nrad)

!    if (nl%iclrsky /= 0) then
!       rshort_clr      (iw) = surface_avg( swdflxc_sfc )
!       rshortup_clr    (iw) = surface_avg( swuflxc_sfc )
!       rshort_top_clr  (iw) = swdflxc(nrad)
!       rshortup_top_clr(iw) = swuflxc(nrad)
!    endif

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

  endif

  ! Compute the longwave fluxes and heating rates

  if (ilwrtyp > 0) then

     icloud  = icld
     iaeros  = iaer

     if (fland > .05) then
        isfc_dt = 1
     else
        isfc_dt = 0
        dlwuflxt_dT    (:) = 0.
        dlwuflxt_sfc_dT(:) = 0.
     endif

     if ( any( cldfr(1:mza-koff) >= cldmin .and. cldfr(1:mza-koff) <= cldmax ) ) then

        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_lw(nrad   , icloud , mcica_seed(1:4,iw), cldfr, &
                             cicewp , cliqwp , taucldl           ,        &
                             ciwpmcl, clwpmcl, taucmcl           , inflg  )
     else

        ! No fractional cloudiness (cloud fraction either 0 or 1)

        if (inflg == 0) then

           do k = 1, nrad
              if (cldfr(k) > 0.5) then
                 do ib = 1, nbndlw
                    do ig = nga_lw(ib), ngs_lw(ib)
                       taucmcl(ig,k) = taucldl(ib,k)
                    enddo
                 enddo
              endif
           enddo

        else

           do k = 1, nrad
              if (cldfr(k) > 0.5) then
                 do ig = 1, ngptlw
                    clwpmcl(ig,k) = cliqwp(k)
                    ciwpmcl(ig,k) = cicewp(k)
                 enddo
              endif
           enddo

        endif

     endif

     call rrtmg_lw(nrad   , nsfc    , iaeros     , isfc_dT    ,                   &
                   play   , tlay    , tsfc       , frac_sfck(:,iw)     , coldry , &
                   h2ovmr , o3vmr   , co2vmr     , ch4vmr     , n2ovmr , o2vmr  , &
                   covmr  , cfc11vmr, cfc12vmr   , cfc22vmr   , ccl4vmr, emis   , &
                   pwvcm  , inflg   , iceflg     , liqflg     , cldfr  ,          &
                   taucmcl, ciwpmcl , clwpmcl    , reice      , reliq  , tauaerl, &
                   lwuflxt, lwdflxt , lwuflxt_sfc, Dlwuflxt_DT, Dlwuflxt_sfc_DT )

     do k = ka-1, mza
        krad = k - koff + 1

        dlong    (k,iw) = lwdflxt    (krad)
        ulong    (k,iw) = lwuflxt    (krad)
        Dulong_DT(k,iw) = Dlwuflxt_dT(krad)
     enddo

     do krad = 1, nsfc
         ulong_sfc   (krad,iw) =  lwuflxt_sfc   (krad)
        Dulong_sfc_DT(krad,iw) = Dlwuflxt_sfc_DT(krad)
     enddo

     if (nsfc == 1) then

        do k = ka-1, mza
           flux_net(k) = ulong(k,iw) - dlong(k,iw)
        enddo

        do k = ka, mza
           fthrd_lw(k,iw) = (flux_net(k-1) - flux_net(k)) * dzit(k) * cpi
        enddo

     else

        fthrd_lw(ka,iw) = ( ulong(ka-1,iw) - dlong(ka-1,iw) &
                          - ulong(ka  ,iw) + dlong(ka  ,iw) ) * dzit(ka) * cpi

        do ks = 2, nsfc
           k = ka + ks - 1
           fthrd_lw(k,iw) = ( ulong(k-1,iw)    * (1. - frac_sfck(ks,iw)) &
                            + ulong_sfc(ks,iw) *       frac_sfck(ks,iw)  &
                            - dlong(k-1,iw) - ulong(k,iw) + dlong(k,iw) ) * dzit(k) * cpi
        enddo

        do k = ka + nsfc - 1, mza
           flux_net(k) = ulong(k,iw) - dlong(k,iw)
        enddo

        do k = ka + nsfc, mza
           fthrd_lw(k,iw) = (flux_net(k-1) - flux_net(k)) * dzit(k) * cpi
        enddo

     endif

     rlong      (iw) = surface_avg( lwdflxt )
     rlongup    (iw) = surface_avg( lwuflxt_sfc )
     rlongup_top(iw) = lwuflxt(nrad)

!    if (nl%iclrsky /= 0) then
!       rlong_clr      (iw) = surface_avg( lwdflxc_sfc )
!       rlongup_clr    (iw) = surface_avg( lwuflxt_sfc )
!       rlongup_top_clr(iw) = lwuflxc(nrad)
!    endif

     do krad = 1, nsfc
        rlong_ks(krad,iw) = lwdflxt(krad)
     enddo

     if (do_pblforc) then
        krad1 = max(kpblh(iw) - 1, ka)
        krad2 = min(kpblh(iw) + 1, mza)

        do k = krad1, krad2
           if (cloud_frac(k,iw) > 0.05) then
              pbl_cld_forc(iw) = pbl_cld_forc(iw) - &
                   min(fthrd_lw(k,iw) + fthrd_sw(k,iw), 0.) * dzt(k) / rho(k,iw)
           endif
        enddo

        pbl_cld_forc(iw) = max(0., pbl_cld_forc(iw) - 1.e-5 * dzt(kpblh(iw)))
     endif

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

    if (iswrtyp > 0 .and. coszen > cosz_min) then

       do ib = 1, nbndsw
          tau = ( fint0 * cloud_props(l)%extsw(ib, index  ) &
                + fint1 * cloud_props(l)%extsw(ib, index+1) ) * watp

          ssa = fint0 * cloud_props(l)%ssasw(ib, index  ) &
              + fint1 * cloud_props(l)%ssasw(ib, index+1)

          asm = fint0 * cloud_props(l)%asysw(ib, index  ) &
              + fint1 * cloud_props(l)%asysw(ib, index+1)

          tauclds(ib,krad) = tauclds(ib,krad) + tau
          ssaclds(ib,krad) = ssaclds(ib,krad) + tau * ssa
          asmclds(ib,krad) = asmclds(ib,krad) + tau * ssa * asm

       enddo

    endif

    if (ilwrtyp > 0) then

       do ib = 1, nbndlw
          tau = ( fint0 * cloud_props(l)%abslw(ib, index  ) &
                + fint1 * cloud_props(l)%abslw(ib, index+1) ) * watp

          taucldl(ib,krad) = taucldl(ib,krad) + tau
       enddo

    endif

  end subroutine lookup_rrtmg_cld_optics


end subroutine rrtmg_raddriv

!==============================================================================

subroutine update_rrtmg_lw_intermediate( iw, ka, lsw )

  use mem_ijtabs,  only: itab_w
  use mem_grid,    only: mza, dzit
  use mem_radiate, only: ulong, dlong, ulong_sfc, Dulong_DT, Dulong_sfc_DT, &
                         Tskin_rad, fthrd_lw, rlong_ks
  use mem_sfcg,    only: sfcg, itab_wsfc
  use consts_coms, only: stefan, cpi
  use mem_turb,    only: frac_sfck, frac_sfc

  implicit none

  integer, intent(in) :: iw
  integer, intent(in) :: ka
  integer, intent(in) :: lsw

  integer :: jsfc, iwsfc, jasfc, k, kw, ks
  real    :: wtk, emiss
  real    :: rlongup_ks(lsw), rlong_albedo_ks(lsw)
  real    :: DTsfc(lsw), Dulong, Dulong_sfc, Dulong_toa, dfact
  real    :: ulong_cur_sfc(lsw)
  real    :: ulong_cur(mza)
  real    :: dlong_cur(mza)
  real    :: flux_net(mza)

  real, parameter :: gamma_down = 0.2

  ! Loop over SFC grid cells beneath this ATM grid column and sum
  ! SFC grid radiative properties to this ATM column

  if (itab_w(iw)%jsfc2 == 1) then

     iwsfc = itab_w(iw)%iwsfc(1)
     rlongup_ks     (1) = sfcg%rlongup     (iwsfc)
     rlong_albedo_ks(1) = sfcg%rlong_albedo(iwsfc)

  else

     rlongup_ks     (:) = 0.
     rlong_albedo_ks(:) = 0.

     do jsfc = 1, itab_w(iw)%jsfc2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)
        wtk   = itab_wsfc(iwsfc)%arcoarkw(jasfc)
        kw    = itab_wsfc(iwsfc)%kwatm(jasfc)
        ks    = kw - ka + 1

        rlongup_ks     (ks) = rlongup_ks     (ks) + wtk * sfcg%rlongup     (iwsfc)
        rlong_albedo_ks(ks) = rlong_albedo_ks(ks) + wtk * sfcg%rlong_albedo(iwsfc)
     enddo

  endif

  ! Change in upward longwave at surface overlaps
  do ks = 1, lsw
     emiss             = 1. - rlong_albedo_ks(ks)
     DTsfc        (ks) = (rlongup_ks(ks) / (emiss * stefan))**0.25 - Tskin_rad(ks,iw)
     ulong_cur_sfc(ks) = ulong_sfc(ks,iw) + DTsfc(ks) * Dulong_sfc_DT(ks,iw)
  enddo

  ! Change in upward longwave at lowest model level
  ulong_cur(ka)   = ulong(ka,iw) + DTsfc(1) * Dulong_DT(ka,iw)
  ulong_cur(ka-1) = ulong_cur_sfc(1)

  ! Special for shaved cells
  do ks = 2, lsw
     k = ka + ks - 1

     ! redefine DTsfc to be an average over cell overlaps at and below ks
     DTsfc(ks) = DTsfc(ks)   *       frac_sfck(ks,iw)  &
               + DTsfc(ks-1) * (1. - frac_sfck(ks,iw))

     ulong_cur(k) = ulong(k,iw) + DTsfc(ks) * Dulong_DT(k,iw)
  enddo

  ! Above any surface overlaps
  do k = ka + lsw , mza
     ulong_cur(k) = ulong(k,iw) + DTsfc(lsw) * Dulong_DT(k,iw)
  enddo

  ! Change in downward longwave is a factor (gamma_down, suggested  0.2) of
  ! the change in upward longwave at the surface, scaled to go to zero at TOA.
  ! (from Hogan and Bozzo, 2015, as used by ECMWF)

  if (gamma_down > 1.e-6 .and. abs(DTsfc(lsw)) > 1.e-6) then

     if (lsw == 1) then
        Dulong_sfc =  ulong_cur_sfc(1) - ulong_sfc(1,iw)
     else
        Dulong_sfc = sum( frac_sfc(1:lsw,iw) * (ulong_cur_sfc(1:lsw) - ulong_sfc(1:lsw,iw)) )
     endif

     Dulong_toa = ulong_cur(mza) - ulong(mza,iw)

     if ( abs(Dulong_sfc - Dulong_toa) < 1.e-6) then

        ! Avoid divide by zero
        do k = ka-1, mza
           dlong_cur(k) = dlong(k,iw) + gamma_down * (ulong_cur(k) - ulong(k,iw))
        enddo

     else

        ! Scale dlong to go to zero at TOA
        dfact = gamma_down * Dulong_sfc / (Dulong_sfc - Dulong_toa)

        do k = ka-1, mza
           Dulong       = ulong_cur(k) - ulong(k,iw)
           dlong_cur(k) = dlong(k,iw) + dfact * (Dulong - Dulong_toa)
        enddo

     endif

  else

     dlong_cur(ka-1:mza) = dlong(ka-1:mza,iw)

  endif

  ! Update longwave heating rate

  if (lsw == 1) then

     do k = ka-1, mza
        flux_net(k) = ulong_cur(k) - dlong_cur(k)
     enddo

     do k = ka, mza
        fthrd_lw(k,iw) = (flux_net(k-1) - flux_net(k)) * dzit(k) * cpi
     enddo

  else

     fthrd_lw(ka,iw) = ( ulong_cur(ka-1) - dlong_cur(ka-1) &
                       - ulong_cur(ka  ) + dlong_cur(ka  ) ) * dzit(ka) * cpi

     do ks = 2, lsw
        k = ka + ks - 1
        fthrd_lw(k,iw) = ( ulong_cur(k-1)    * (1. - frac_sfck(ks,iw)) &
                         + ulong_cur_sfc(ks) *       frac_sfck(ks,iw)  &
                         - dlong_cur(k-1) - ulong_cur(k) + dlong_cur(k) ) * dzit(k) * cpi
     enddo

     do k = ka + lsw - 1, mza
        flux_net(k) = ulong_cur(k) - dlong_cur(k)
     enddo

     do k = ka + lsw, mza
        fthrd_lw(k,iw) = (flux_net(k-1) - flux_net(k)) * dzit(k) * cpi
     enddo

  endif

  ! Update downward longwve for the surface computations

  do ks = 1, lsw
     k = ka + ks - 2
     rlong_ks(ks,iw) = dlong_cur(k)
  enddo

end subroutine update_rrtmg_lw_intermediate

!==============================================================================

end module rrtmg_driver
