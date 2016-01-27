subroutine rrtmg_raddriv(iw, ka, nrad, koff)

  use mem_grid,    only: mza, zm, zt, glatw, glonw, dzt, dzim
  use mem_basic,   only: rho, press, theta, tair, sh_v, sh_w, thil, wc
  use misc_coms,   only: io6, iswrtyp, ilwrtyp, time8, nqparm, dtlm, &
                         icfrac, cfracrh1, cfracrh2, cfraccup, do_chem
  use consts_coms, only: stefan, eps_virt, eps_vapi, grav, solar, cp, pi1, t00
  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         albedt_beam, albedt_diffuse, rshort_diffuse, &
                         rlong_albedo, solfac, cloud_frac, rshort_clr, &
                         rshortup_clr, rshort_top_clr, rshortup_top_clr, &
                         rlong_clr, rlongup_clr, rlongup_top_clr, &
                         par, par_diffuse, uva, uvb, uvc
  use micro_coms,  only: ncat, rxmin, emb0, reffcof, pwmasi, dmncof, jhabtab
  use mem_ijtabs,  only: itab_w
  use mem_cuparm,  only: kcutop, kcubot, cbmf, qwcon, conprr
  use oname_coms,  only: nl
  use rrtmg_cloud, only: cloud_props
  use mem_turb,    only: frac_land, pblh, kpblh, wtv0, hkm
  use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h
  use mem_para,    only: myrank
  use mem_ijtabs,  only: itab_w
  use mem_mclat,   only: rad_mclat
  use cgrid_defn,  only: acflux_dir_dn_tot, acflux_dif_dn_tot, acflux_dif_up_tot, &
                         acflux_dir_dn_clr, acflux_dif_dn_clr, acflux_dif_up_clr

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
  real :: tsfc  (ncol)
  real :: asdir (ncol)
  real :: aldir (ncol)
  real :: asdif (ncol)
  real :: aldif (ncol)
  real :: emis  (ncol, nbndlw)

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
  real :: exl(nrad)
  real :: dzl(nrad)

  real :: swuflx (ncol, nrad+1)
  real :: swdflx (ncol, nrad+1)
  real :: swhr   (ncol, nrad  )
  real :: swuflxc(ncol, nrad+1)
  real :: swdflxc(ncol, nrad+1)
  real :: swhrc  (ncol, nrad  )

  real :: swuflxt_band    (nrad+1,nbndsw)
  real :: swdflxt_band    (nrad+1,nbndsw)
  real :: swdflxt_band_dir(nrad+1,nbndsw)

  real :: swuflxc_band    (nrad+1,nbndsw)
  real :: swdflxc_band    (nrad+1,nbndsw)
  real :: swdflxc_band_dir(nrad+1,nbndsw)

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

  integer :: k, krad, icloud, iaeros, index, ib, ig, n
  integer :: iplon, irng, permuteseed, ns, nt
  integer :: mc, mcat, ih, l, num, ntim, ngbmsw, ngbmlw
  logical :: iconv

  real :: tau, ssa, asm
  real :: r_ef, dmean, watp, rstart, rend, rscale, fint0, fint1

  real :: frac(mza)

  real, external :: rhovsl, rhovsi

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

  integer, parameter :: kradcat(16) = (/1,8,6,6,5,4,4,2,8,8,7,9,8,8,7,9/)

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
     exl     (krad) = theta(k,iw) / tair(k,iw)
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
! if there are any resolved hydrometeors.

  call get_cloud_frac(iw, ka, frac)

! Get optical properties of resolved clouds

  do mc = 1, mcat

     do k = ka, mza
        krad = k - koff

        if (rx(k,mc) >= rxmin(mc) .and. emb(k,mc) >= emb0(mc)) then

! Set lower bound on frac(k) because there is condensate

           frac(k) = max(frac(k),0.1)

           cldfr(1,krad) = frac(k)

           ih = jhcat(k,mc)
           l  = kradcat(ih)

           num    = cloud_props(l)%num
           rstart = cloud_props(l)%start
           rend   = cloud_props(l)%end

           r_ef = 1.e6 * reffcof(ih) * emb(k,mc) ** pwmasi(ih)
           r_ef = max(rstart, min(rend, r_ef))

           ! ice or liquid water path in g/m^2
           watp = rx(k,mc) * dl(krad) * 1000. * dzt(k)
           watp = watp / frac(k)

           rscale = (r_ef - rstart) / cloud_props(l)%delr
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
        endif
     enddo
  enddo

  ! Now include the optical properties of subgrid convective clouds

  if (iconv) then

     do k = kcubot(iw), kcutop(iw)
        krad = k - koff

        if (frac(k) > 1.e-10) then
              
           cldfr(1,krad) = frac(k)

           ! water path in g/m^2
           watp = max(qwcon(k,iw),1.e-5) * rho(k,iw) * 1000. * dzt(k)
           tc   = tair(k,iw) - t00
              
           if (tc > -10.0) then

              ! Add convective cloud water to cloud drops if warmer then 10C
              l = kradcat(1)

              ! Hardwire droplet effective radius to 14 over land 
              ! and 8 over sea following CAM physics
              r_ef = rliqland + (rliqocean-rliqland) * fland

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

           num    = cloud_props(l)%num
           rstart = cloud_props(l)%start
           rend   = cloud_props(l)%end
   
           r_ef = max(rstart, min(rend, r_ef))

           rscale = (r_ef - rstart) / cloud_props(l)%delr
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

        endif
     enddo
  endif

  ! Save cloud fraction in 3D variable for output or plotting

  cloud_frac(1:ka-1,iw) = 0.0

  do k = ka, mza
     krad = k - koff
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
     iplon = 1
     irng = 0
     permuteseed = 1
     ntim = nint(time8/dtlm(1))

     if ( any( cloud_frac(ka:mza,iw) >= 0.001 .and. cloud_frac(ka:mza,iw) <= 0.999 ) ) then
     
        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_sw(iw, ntim, iplon, ncol, nrad, icloud, &
                             permuteseed, irng, play, &
                             cldfr, cicewp, cliqwp, reice, reliq,  &
                             tauclds, ssaclds, asmclds, fsfclds, &
                             cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                             taucmcl, ssacmcl, asmcmcl, fsfcmcl)
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

     call rrtmg_sw(ncol    ,nrad    ,icloud  ,iaeros  ,                 &
                   play    ,plev    ,tlay    ,tlev    ,tsfc   ,         &
                   h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr   ,&
                   asdir   ,asdif   ,aldir   ,aldif   ,                 &
                   coszen  ,solfac  ,dyofyr  ,solar   ,                 &
                   inflg   ,iceflg  ,liqflg  ,cldfmcl ,                 &
                   taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl ,                 &
                   ciwpmcl ,clwpmcl ,reicmcl ,relqmcl ,                 &
                   tauaers ,ssaaers ,asmaers ,ecaer   ,                 &
                   swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc  ,&
                   swuflxt_band     ,swdflxt_band     ,swuflxc_band    ,& 
                   swdflxc_band     ,swdflxt_band_dir ,swdflxc_band_dir)

     rshort        (iw) = swdflx(1,1)
     rshort_diffuse(iw) = swdflx(1,1) - sum(swdflxt_band_dir(1,1:nbndsw))
     rshort_top    (iw) = swdflx(1,nrad+1)
     rshortup_top  (iw) = swuflx(1,nrad+1)
     albedt        (iw) = swuflx(1,1) / swdflx(1,1)

     rshort_clr      (iw) = swdflxc(1,1)
     rshortup_clr    (iw) = swuflxc(1,1)
     rshort_top_clr  (iw) = swdflxc(1,nrad+1)
     rshortup_top_clr(iw) = swuflxc(1,nrad+1)

     par(iw) = 0.5268*swdflxt_band(1, 9) + swdflxt_band(1,10) &
             + 0.4724*swdflxt_band(1,11)

     par_diffuse(iw) = par(iw) - ( 0.5268*swdflxt_band_dir(1, 9) &
                                 +        swdflxt_band_dir(1,10) &
                                 + 0.4724*swdflxt_band_dir(1,11) )

     uva(iw) = 0.5276*swdflxt_band(1,11) + 0.3932*swdflxt_band(1,12)

     uvb(iw) = 0.3708*swdflxt_band(1,12)

     uvc(iw) = 0.2360*swdflxt_band(1,12) + swdflxt_band(1,13)

     do k = ka, mza
        krad = k - koff
        fthrd_sw(k,iw) = swhr(1,krad) * exl(krad) / 86400.0
     enddo

     if (do_chem == 1) then
        call rrtmg_to_cmaq()
     endif

  else

     if (do_chem == 1) then
        do n = 1, 7
           acflux_dir_dn_tot(:,:,iw) = 0.0
           acflux_dif_dn_tot(:,:,iw) = 0.0
           acflux_dif_up_tot(:,:,iw) = 0.0

           acflux_dir_dn_clr(:,:,iw) = 0.0
           acflux_dif_dn_clr(:,:,iw) = 0.0
           acflux_dif_up_clr(:,:,iw) = 0.0
        enddo
     endif

  endif

  ! Compute the longwave fluxes and heating rates

  if (ilwrtyp > 0) then

     icloud = icld
     iaeros = iaer
     iplon = 1
     irng = 0
     permuteseed = 150

     if ( any( cloud_frac(ka:mza,iw) >= 0.001 .and. cloud_frac(ka:mza,iw) <= 0.999 ) ) then
     
        ! Subgrid (fractional) cloudiness present

        call mcica_subcol_lw(iw     , ntim       , iplon     , ncol      , nrad  , &
                             icloud , permuteseed, irng      , play      ,         &
                             cldfr  , cicewp     , cliqwp    , reice     , reliq , &
                             taucldl, cldfmcl_lw , ciwpmcl_lw, clwpmcl_lw,         &
                             reicmcl, relqmcl    , taucmcl_lw                      )
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
                 clwpmcl_lw(ig,1,k) = cliqwp       (1,k)
                 ciwpmcl_lw(ig,1,k) = cicewp       (1,k)
                 taucmcl_lw(ig,1,k) = taucmcl_lw(ib,1,k)
              enddo

           endif

        enddo

     endif

     call rrtmg_lw(ncol       ,nrad       ,icloud     ,idrv       ,                  &
                   play       ,plev       ,tlay       ,tlev       ,tsfc    ,         & 
                   h2ovmr     ,o3vmr      ,co2vmr     ,ch4vmr     ,n2ovmr  ,o2vmr  , &
                   cfc11vmr   ,cfc12vmr   ,cfc22vmr   ,ccl4vmr    ,emis    ,         &
                   inflg      ,iceflg     ,liqflg     ,cldfmcl_lw ,                  &
                   taucmcl_lw ,ciwpmcl_lw ,clwpmcl_lw ,reicmcl    ,relqmcl ,tauaerl, &
                   lwuflx     ,lwdflx     ,lwhr       ,lwuflxc    ,lwdflxc ,lwhrc  , &
                   duflx_dt   ,duflxc_dt                                             )

     rlong      (iw) = lwdflx(1,1)
     rlongup    (iw) = lwuflx(1,1)
     rlongup_top(iw) = lwuflx(1,nrad+1)

     rlong_clr      (iw) = lwdflxc(1,1)
     rlongup_clr    (iw) = lwuflxc(1,1)
     rlongup_top_clr(iw) = lwuflxc(1,nrad+1)

     do k = ka, mza
        krad = k - koff
        fthrd_lw(k,iw) = lwhr(1,krad) * exl(krad) / 86400.0
     enddo

  endif


contains


  subroutine rrtmg_to_cmaq()

    implicit none

    real,    parameter :: mu1 = 2.0

    integer, parameter :: nb(6) = (/ 12, 12, 12, 12, 12, 11 /)
    real,    parameter :: cc(6) = (/ 0.080617, 0.111306, 0.064996, 0.107632, 0.388907, 0.663493 /)

    integer, parameter :: n7(4) = (/ 8, 9, 10, 11 /)
    real,    parameter :: c7(4) = (/ 0.231458, 1.000000, 1.000000, 0.336507 /)
    
    real               :: mu
    integer            :: k, n
    real               :: irr_dir_dn_tot(7)
    real               :: irr_dif_dn_tot(7)
    real               :: irr_dif_up_tot(7)
    real               :: irr_dir_dn_clr(7)
    real               :: irr_dif_dn_clr(7)
    real               :: irr_dif_up_clr(7)

    mu = 1.0 / cosz(iw)

    do k = 1, mza-koff

       do n = 1, 6
          irr_dir_dn_tot(n) = cc(n) * swdflxt_band_dir(k,nb(n))
          irr_dif_dn_tot(n) = cc(n) * (swdflxt_band(k,nb(n)) - swdflxt_band_dir(k,nb(n)))
          irr_dif_up_tot(n) = cc(n) * swuflxt_band(k,nb(n))

          irr_dir_dn_clr(n) = cc(n) * swdflxc_band_dir(k,nb(n))
          irr_dif_dn_clr(n) = cc(n) * (swdflxc_band(k,nb(n)) - swdflxc_band_dir(k,nb(n)))
          irr_dif_up_clr(n) = cc(n) * swuflxc_band(k,nb(n))
       enddo

       irr_dir_dn_tot(7) = sum( c7(:) * swdflxt_band_dir(k,n7(:)) )
       irr_dif_dn_tot(7) = sum( c7(:) * (swdflxt_band(k,n7(:)) - swdflxt_band_dir(k,n7(:))) )
       irr_dif_up_tot(7) = sum( c7(:) * swuflxt_band(k,n7(:)) )

       irr_dir_dn_clr(7) = sum( c7(:) * swdflxc_band_dir(k,n7(:)) )
       irr_dif_dn_clr(7) = sum( c7(:) * (swdflxc_band(k,n7(:)) - swdflxc_band_dir(k,n7(:))) )
       irr_dif_up_clr(7) = sum( c7(:) * swuflxc_band(k,n7(:)) )

       do n = 1, 7
          acflux_dir_dn_tot(k,n,iw) = irr_dir_dn_tot(n) * mu
          acflux_dif_dn_tot(k,n,iw) = irr_dif_dn_tot(n) * mu1
          acflux_dif_up_tot(k,n,iw) = irr_dif_up_tot(n) * mu1

          acflux_dir_dn_clr(k,n,iw) = irr_dir_dn_clr(n) * mu
          acflux_dif_dn_clr(k,n,iw) = irr_dif_dn_clr(n) * mu1
          acflux_dif_up_clr(k,n,iw) = irr_dif_up_clr(n) * mu1
       enddo

    enddo

  end subroutine rrtmg_to_cmaq


end subroutine rrtmg_raddriv
