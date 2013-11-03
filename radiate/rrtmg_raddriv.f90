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

  use consts_coms, only: stefan, eps_virt, eps_vapi, grav, solar, cp, pi1

  use mem_radiate, only: rshort, rlong, fthrd_lw, rlongup, cosz, albedt, &
                         rshort_top, rshortup_top, rlongup_top, fthrd_sw, &
                         albedt_beam, albedt_diffuse, rshort_diffuse, &
                         rlong_albedo, solfac

  use micro_coms,  only: ncat, rxmin, emb0, reffcof, pwmasi, dnfac

  use parrrtm,              only: nbndlw
  use parrrsw,              only: nbndsw, jpb1, jpb2
  use rrtmg_sw_rad_nomcica, only: rrtmg_sw_nomcica
  use rrtmg_sw_rad,         only: rrtmg_sw
  use rrtmg_lw_rad_nomcica, only: rrtmg_lw_nomcica
  use rrtmg_lw_rad,         only: rrtmg_lw
  use rrsw_cld,             only: extliq1, ssaliq1, asyliq1, extice2, ssaice2, asyice2
  use rrlw_cld,             only: absliq1, absice2

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
  
  integer :: jhcat(mza,ncat)  ! hydrom category table with ice habits

  real :: rhov (mza)       ! vapor density [kg_vap/m^3]
  real :: rx   (mza,ncat)  ! hydrom bulk spec dens [kg_hyd/kg_air]
  real :: cx   (mza,ncat)  ! hydrom bulk number [num_hyd/kg_air]
  real :: emb  (mza,ncat)  ! hydrom mean particle mass [kg/particle]

  integer :: mc, mcat, ih
  real    :: r_ef, watp


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

  integer :: k, krad, icloud, index, ib, jb
  real    :: lwc, tau, ssa, asm, fint, ext, abs

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

! Convert water vapor to molar mixing ratio, pressure to mb, and
! ozone to molar mixing ratio

  do krad = 1, nrad
     h2ovmr(1,krad) = h2ovmr(1,krad) * eps_vapi / dl(krad)
     play  (1,krad) = play  (1,krad) * 0.01
     o3vmr (1,krad) = o3vmr (1,krad) * amdo3
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
  fsfclds(:,ncol,:) = 0.0
  cicewp   (ncol,:) = 0.0
  cliqwp   (ncol,:) = 0.0
  reice    (ncol,:) = 0.0
  reliq    (ncol,:) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fill arrays rx, cx, and emb with hydrometeor properties

  call cloudprep_rad(iw,ka,mcat,jhcat,rhov,rx,cx,emb)

  do mc = 1, 8

     do k = ka, mza-1
        krad = k - koff

        if (rx(k,mc) >= rxmin(mc) .and. emb(k,mc) >= emb0(mc)) then

           ih = jhcat(k,mc)
           cldfr(1,krad) = 1.0
           r_ef = 1.e6 * reffcof(ih) * emb(k,mc) ** pwmasi(ih)
           watp = rx(k,mc) * dl(krad) * 1000. * dzt(k)

           if (any( ih == (/ 1, 2, 8 /) )) then

              ! Spherical water droplet model for solar

              if (iswrtyp == 2 .and. cosz(iw) >= 0.03) then
                 do ib = 1, nbndsw
                    jb = ib + jpb1 - 1

                    if (r_ef <= 2.5) then
                       ext = extliq1(1,jb)
                       ssa = ssaliq1(1,jb)
                       asm = asyliq1(1,jb)
                    elseif (r_ef >= 59.5) then
                       ext = extliq1(58,jb)
                       ssa = ssaliq1(58,jb)
                       asm = asyliq1(58,jb)
                    else
                       index = max(1, min( int(r_ef - 1.5), 57) )
                       fint  = r_ef - 1.5 - real(index)

                       ext = extliq1(index,jb) + &
                             fint * (extliq1(index+1,jb) - extliq1(index,jb))

                       ssa = ssaliq1(index,jb) + &
                             fint * (ssaliq1(index+1,jb) - ssaliq1(index,jb))

                       asm = asyliq1(index,jb) + &
                             fint * (asyliq1(index+1,jb) - asyliq1(index,jb))
                    endif

                    tau = watp * ext
                    tauclds(ib,1,krad) = tauclds(ib,1,krad) + tau
                    ssaclds(ib,1,krad) = ssaclds(ib,1,krad) + tau * ssa
                    asmclds(ib,1,krad) = asmclds(ib,1,krad) + tau * ssa * asm

                 enddo
              endif

              ! Spherical water droplet model for longwave

              if (ilwrtyp == 2) then
                 do ib = 1, nbndlw
                 
                    if (r_ef <= 2.5) then
                       abs = absliq1(1,ib)
                    elseif (r_ef >= 59.5) then
                       abs = absliq1(58,ib)
                    else
                       index = max(1, min( int(r_ef - 1.5), 57) )
                       fint  = r_ef - 1.5 - real(index)
                       abs = absliq1(index,ib) + &
                             fint * (absliq1(index+1,ib) -absliq1(index,ib))
                    endif

                    taucldl(ib,1,krad) = taucldl(ib,1,krad) + watp * abs

                 enddo
              endif

           else

              ! Spherical ice model for solar

              if (iswrtyp == 2 .and. cosz(iw) >= 0.03) then
                 do ib = 1, nbndsw
                    jb = ib + jpb1 - 1

                    if (r_ef <= 5.0) then
                       ext = extice2(1,jb)
                       ssa = ssaice2(1,jb)
                       asm = asyice2(1,jb)
                    elseif (r_ef >= 131) then
                       ext = extice2(43,jb)
                       ssa = ssaice2(43,jb)
                       asm = asyice2(43,jb)
                    else
                       index = max(1, min( int((r_ef - 2.0)/3.0), 42) )
                       fint  = (r_ef - 2.0)/3.0 - real(index)

                       ext = extice2(index,jb) + &
                             fint * (extice2(index+1,jb) - extice2(index,jb))

                       ssa = ssaice2(index,jb) + &
                             fint * (ssaice2(index+1,jb) - ssaice2(index,jb))

                       asm = asyice2(index,jb) + &
                             fint * (asyice2(index+1,jb) - asyice2(index,jb))
                    endif
                 
                    tau = watp * ext
                    tauclds(ib,1,krad) = tauclds(ib,1,krad) + tau
                    ssaclds(ib,1,krad) = ssaclds(ib,1,krad) + tau * ssa
                    asmclds(ib,1,krad) = asmclds(ib,1,krad) + tau * ssa * asm

                 enddo
              endif

              if (ilwrtyp == 2) then
                 do ib = 1, nbndlw

                    if (r_ef <= 5.0) then
                       abs = absice2(1,ib)
                    elseif (r_ef >= 131) then
                       abs = absice2(43,ib)
                    else
                       index = max(1, min( int((r_ef - 2.0)/3.0), 42) )
                       fint  = (r_ef - 2.0)/3.0 - real(index)
                       abs = absice2(index,ib) + &
                             fint * (absice2(index+1,ib) - absice2(index,ib))
                    endif

                    taucldl(ib,1,krad) = taucldl(ib,1,krad) + watp * abs

                 enddo
              endif

           endif

        endif
     enddo
  enddo

  if (iswrtyp == 2 .and. cosz(iw) > 0.03) then

     ! Combine optical properties

     do k = ka, mza-1
        krad = k - koff
        do ib = 1, nbndsw
           if (tauclds(ib,1,krad) > 1.e-15 .and. ssaclds(ib,1,krad) > 1.e-15) then
              asmclds(ib,1,krad) = asmclds(ib,1,krad) / ssaclds(ib,1,krad)
              ssaclds(ib,1,krad) = ssaclds(ib,1,krad) / tauclds(ib,1,krad)
              fsfclds(ib,1,krad) = asmclds(ib,1,krad) * asmclds(ib,1,krad)
           endif
        enddo
     enddo
          
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
