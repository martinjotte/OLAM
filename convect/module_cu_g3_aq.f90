MODULE module_cu_g3_aq

  private
  public :: grell_aq_driver

CONTAINS

  subroutine grell_aq_driver ( mrl )

    use mem_ijtabs, only: jtab_w, jtw_prog

    implicit none

    integer, intent(in) :: mrl
    integer             :: j, iw

    !$omp parallel do private(iw) schedule(guided)
    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       call grell_aqmix( iw )

    enddo
    !$omp end parallel do

  end subroutine grell_aq_driver




  subroutine grell_aqmix(iw)

    use mem_cuparm,  only: kcutop, kcubot, kudbot, cbmf, iactcu, &
                           cddf, kddtop, kddmax, kstabi, qwcon, &
                           cu_pwa, cu_pev
    use mem_basic,   only: rho, press, tair
    use mem_grid,    only: lpw, arw, arw0, arw0i, volt
    use aq_data,     only: convf, convb
    use cgrid_defn,  only: cgrid
    use cgrid_spcs,  only: nspcsd
    use mem_radiate, only: cosz
    use mem_co2,     only: rr_co2, i_co2, co2_sh2ppm, co2_ppm2sh, co2_initppm
    use misc_coms,   only: dtlong, nqparm
    use mem_ijtabs,  only: itab_w

    implicit none

    integer, intent(in) :: iw

    real, allocatable :: dens(:), area(:), voa(:), qrc(:), pw(:), to(:), po(:), pwd(:)
    real, allocatable :: tr(:,:), dtr(:,:), tr1(:,:)

    real              :: a0, xmb, edt, xmbi, dtl
    integer           :: ktf, ktop, kbcon, k22, kstab, jmin, kdet
    integer           :: ka, k, kc, ns, nt, itype
    logical           :: idd, dark

    ! exit if no convection
    if (iactcu(iw) < 1 .or. cbmf(iw) < 1.e-7) return

    ! exit if we are not using Grell convection
    if (nqparm( itab_w(iw)%mrlw ) /= 2) return

    ka = lpw(iw)

    ktop  = kcutop(iw) - ka + 1
    kbcon = kcubot(iw) - ka + 1
    k22   = kudbot(iw) - ka + 1
    kstab = kstabi(iw) - ka + 1
    xmb   = cbmf(iw)
    xmbi  = 1. / xmb
    dtl   = dtlong

    if (iactcu(iw) == 1 .and. cddf(iw) >= 1.e-7) then
       idd  = .true.
       edt  = cddf(iw) * xmbi
       jmin = kddtop(iw) - ka + 1
       kdet = kddmax(iw) - ka + 1
    else
       idd = .false.
       edt   = 0.0
       jmin  = 1
       kdet  = 1
    endif

    ktf = ktop + 2
    nt  = nspcsd + 1

    allocate(dens(ktf))
    allocate(area(ktf))
    allocate(voa (ktf))
    allocate(qrc (ktf))
    allocate(pw  (ktf))
    allocate(pwd (ktf))
    allocate(po  (ktf))
    allocate(to  (ktf))
    allocate(tr1 (ktf,nt))

    allocate( tr (nt,ktf))
    allocate(dtr (nt,ktf))

    do kc = 1, ktf
       k = kc + ka - 1
       dens(kc) = rho(k,iw)
       po  (kc) = press(k,iw)
       to  (kc) = tair(k,iw)
       area(kc) = arw(k-1,iw)
       voa (kc) = volt(k,iw) * arw0i(iw)
    enddo
    area(1) = area(2)

    qrc = 0.0
    pw  = 0.0
    pwd = 0.0

    if     (iactcu(iw) == 1 .or. iactcu(iw) == 2) then
       itype = 1
    elseif (iactcu(iw) == 3) then
       itype = 2
    else
       return
    endif

    do kc = kbcon, ktop
       k = kc + ka - 1
       qrc(kc+1) = qwcon(k,iw)
       if (itype == 1) pw(kc+1) = cu_pwa(k,iw) * xmbi
    enddo

    if (idd) then
       do kc = 1, jmin
          k = kc + ka - 1
          pwd(kc) = cu_pev(k,iw) * xmbi
       enddo
    endif

    ! load CMAQ scalars
    do ns = 1, nspcsd
       do kc = 1, ktf
          k = kc + ka - 1
          tr1(kc,ns) = max( cgrid( k, iw, ns ) * convf( ns ),  1.e-36 )
       enddo
    enddo

    ! load CO2
    if (i_co2 > 0) then
       ns = nspcsd + 1
       do kc = 1, ktf
          k = kc + ka - 1
          tr1(kc,ns) = rr_co2(k,iw) * co2_sh2ppm * 1.e-6  ! convert to mol/molV
       enddo
    else
       tr1(1:ktf,ns) = co2_initppm * 1.e-6  ! convert to mol/molV
    endif

    tr = transpose(tr1)

    a0   = arw0(iw)
    dark = (cosz(iw) < 1.e-3)

    call cup_tracer_mix_aq(tr,dtr,itype,idd,dens,po,to,voa,area,a0,xmb,k22,kbcon,ktop, &
         kstab,edt,jmin,kdet,qrc,pw,pwd,dark,ktf,nt,iw)

    tr1 = transpose(dtr)

    do ns = 1, nspcsd
       do kc = 1, ktop
          k = kc + ka - 1
          cgrid( k, iw, ns ) = max( cgrid( k, iw, ns ) + dtl * tr1( kc, ns ) * convb( ns ), 1.e-30 )
       enddo
    enddo

    if (i_co2 > 0) then
       ns = nspcsd + 1
       do kc = 1, ktf
          k = kc + ka - 1
          rr_co2(k,iw) = rr_co2(k,iw) + dtl * tr1(kc,ns) * co2_ppm2sh * 1.e6
       enddo
    endif

  end subroutine grell_aqmix


  subroutine cup_tracer_mix_aq(tr,dtr,itype,idd,dens,po,to,dz_cup,arw,arw0,xmb,k22,kbcon, &
             ktop,kstabi,edt,jmin,kdet,qrc,pw,pwd,dark,ktf,nt,iw)

    use const_data, only: inv_mwairkg
    use module_cu_g3

    implicit none

    real,    intent(in)  :: xmb, edt, arw0
    integer, intent(in)  :: k22, kbcon, ktop, kstabi, jmin, kdet, ktf
    integer, intent(in)  :: nt, itype, iw
    real,    intent(in)  :: dz_cup(ktf), arw(ktf), dens(ktf), po(ktf), to(ktf)
    real,    intent(in)  :: qrc(ktf), pw(ktf), pwd(ktf)
    real,    intent(in)  :: tr(nt,ktf)
    real,    intent(out) :: dtr(nt,ktf)
    logical, intent(in)  :: idd, dark

    integer              :: ierr, k, kbot, kfull
    real, dimension(ktf) :: cdd, cd, etd, en, zu, zd, pcpflx
    real                 :: tr_cup(nt,ktf), tc(nt,ktf), tcd(nt,ktf)
    real                 :: tr_pw(nt,ktf)
    real, dimension(ktf) :: frac_act  ! fraction of aerosol activated in clouds
    real, dimension(ktf) :: wtbar, wcbar, fracliq
    real                 :: entr_rate, detr_rate
    real                 :: hp_rain, arw0m, rad
    real                 :: dtr_scav(nt,ktf), remov(nt), tr_flx(nt)
    real, dimension(ktf) :: facu1, facu2, facd1, facd2

    ierr     = 0
    frac_act = 0.
    pcpflx   = 0.
    dtr_scav = 0.
    wtbar    = 0.
    wcbar    = 0.

    do k = 1, ktop
       fracliq(k+1) = max(0.0, min(1.0, (to(k) - 238.15) / 35.))
    enddo

    ! fracliq monotonically increasing (no refreezing of drops)
    do k = ktop, 2, -1
       fracliq(k) = max(fracliq(k+1), fracliq(k))
    enddo
    fracliq(1) = fracliq(2)
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
    if (itype == 1) then

       rad = sqrt(arw0)
       if (rad < radius_deep) then
          entr_rate = 0.24 / rad
          detr_rate = detr_deep + (entr_rate - entr_deep)
       else
          entr_rate  = entr_deep
          detr_rate  = detr_deep
       endif

    elseif (itype == 2) then

       entr_rate = entr_shal
       detr_rate = detr_shal

    else

       return

    endif
!
!--- Special for shaved cells
!
    arw0m = nearest(arw0,-1.)
    do k = 1, ktf-1
       if (arw(k) >= arw0m) exit
    enddo
    kfull = k
!
!--- Compute entrainment and detrainment
!
    call cup_up_en_cd(en,cd,facu1,facu2,dz_cup,arw,entr_rate,detr_rate,kfull,ierr,ktf)
!
!--- Increase detrainment in stable layers for deep/mid convection
!
    if (itype == 1 .and. kstabi <= ktop) then
       do k = kstabi, ktop
          cd   (k) = min(cd(k-1)/dz_cup(k-2)+.15*entr_rate,1.5*entr_rate)*dz_cup(k-1)
          facu2(k) = en(k)  / (1.0 + en(k) - 0.5*cd(k))
          facu1(k) = 1.0 - facu2(k)
       enddo
    endif
!
!--- Normalized updraft mass flux profile
!
    call cup_up_nms(zu,cd,en,ktop,ierr,k22,facu1,facu2,arw,arw0,ktf)
!
!--- Tracer concentration at cloud (mid-point) levels
!
    call cup_env_clev_tr(tr,tr_cup,ktf,nt)

    if (itype == 2) then
!
!--- Tracer concentration in updraft for shallow clouds. Includes aqueous chem
!    and Aitken scavenging but no rainout/washout
!
       call cup_up_tracer_shal(k22,dz_cup,cd,zu,facu1,facu2,tr_cup,tr,tc, &
                               fracliq,wtbar,wcbar,frac_act,kbcon,ktop,qrc, &
                               dens,po,to,ierr,dark,ktf,nt,iw)

    else
!
!--- Tracer concentration in updraft for shallow clouds. Includes aqueous chem,
!    Aitken scavenging, and rainout
!
       call cup_up_tracer(k22,dz_cup,cd,zu,facu1,facu2,tr_cup,tr,tc,fracliq, &
                          wtbar,wcbar,frac_act,hp_rain,kbcon,ktop,qrc,pw, &
                          tr_pw,dens,po,to,ierr,dark,ktf,nt,iw)
!
!--- Precipatation washout above level of downdraft
!
       kbot = 1
       if (idd) kbot = jmin+1

       do k = ktop, 1, -1
          pcpflx(k) = pcpflx(k+1) + pw(k+1) + pwd(k)
       enddo

       tr_flx = tr_pw(:,ktop)

       do k = ktop-1, kbot, -1
          call precip_scav( tr(:,k), dtr_scav(:,k), remov, to(k), dens(k), &
               po(k), dz_cup(k), frac_act(k+1), hp_rain, pcpflx(k+1), tr_flx, &
               fracliq(k+1), wcbar(k+1) )

          if (k >= kbcon) then
             tr_flx = tr_flx + remov + tr_pw(:,k)
          else
             tr_flx = tr_flx + remov
          endif
       enddo

    endif

    if (idd) then
!
!--- Downdraft mass flux and entrainment/detrainment rates
!
       call cup_dd_nms(zd,dz_cup,cdd,etd,facd1,facd2,entr_rate, &
                       jmin,ierr,arw,arw0,kdet,kfull,ktf)
!
!--- Scale downdrafts
!
       do k = 2, jmin
          zd(k) = zd(k) * edt
       enddo
!
!--- Tracer concentration in downdraft. Includes evaporation and washout effects
!
       call cup_dd_tracer(tcd,dz_cup,facd1,facd2,to,dens,po,cdd,entr_rate,jmin, &
                          kbcon,ierr,tr,tr_cup,tr_pw,zd,pwd,pcpflx,tr_flx, &
                          frac_act,fracliq,wcbar,hp_rain,nt,ktf)
!
!--- Compute tendencies with downdrafts
!
       call cup_dellas_tr(ierr,dz_cup,dtr,zu,tr,tr_cup,tc,dens, &
                           ktop,k22,cd,en,tcd,zd,cdd,etd,jmin,nt,ktf)

   else
!
!--- Compute tendencies without downdrafts
!
      call cup_dellas_tr_nodd(ierr,dz_cup,dtr,zu,tr,tr_cup,tc,dens, &
                              ktop,k22,cd,en,nt,ktf)
      tr_flx = 0.0

   endif

   ! Include washout above downdraft in tendencies

   if (itype == 1) then
      do k = kbot, ktop-1
         dtr(:,k) = dtr(:,k) + dtr_scav(:,k)
      enddo
   endif

   dtr = dtr * xmb
   if (idd) tr_flx = tr_flx * xmb

!      ipr = 14
!    ipr = 14
!     ipr = n_gc_spc + 51
!     ipr = n_gc_spc + 1
!    ipr = n_gc_spc + n_ae_spc + 1
!    ipr = 14
!   ipr = 171

!      ipr = n_gc_spc + n_ae_spc + 2
!    ipr = 9
!      ipr = 4

!      write(*,*) xmb, k22, ktop, jmin
!      write(*,*) idd, itype
!      do k = ktop+1, 1, -1
!         write(*,'(I4,10g13.6)') k, tr(ipr,k), tc(ipr,k), tcd(ipr,k), dtr(ipr,k), tr_pw(ipr,k), dtr_scav(ipr,k)*dens(k)*dz_cup(k) * inv_mwairkg, pwd(k)
!      enddo
!      write(*,*)
!      write(*,*) sum(dtr(ipr,1:ktop) * dens(1:ktop) * dz_cup(1:ktop) * inv_mwairkg), tr_flx(ipr) !&

!           xmb*sum(tr_pwd(ipr,:))+xmb*sum(tr_pw(ipr,:)) !-xmb*sum(tr_scav(ipr,:)*dens(1:ktop) * dz_cup(1:ktop))
!      stop

  end subroutine cup_tracer_mix_aq



   SUBROUTINE cup_up_tracer(k22,dz_cup,cd,zu,fact1,fact2,tr_cup,tr,tc,fracliq, &
                            wtbar,wcbar,frac_act,hp_rain,kbcon,ktop,qrc,pw, &
                            tr_pw,dens,po,to,ierr,dark,ktf,nt,iw)

     use consts_coms, only: rdry
     use const_data,  only: inv_mwairkg

     implicit none
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        k22,ktop,kbcon,ktf,nt,ierr,iw
     real,    dimension (nt,ktf)                                       &
        ,intent (in   )                   ::                           &
        tr,tr_cup
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        cd,fact1,fact2,qrc,pw,dens,po,to,dz_cup,zu,fracliq
     logical, intent(in) :: dark
!
! input and output
!
     real                                         &
        ,intent (inout)                   ::                           &
        hp_rain
     real,    dimension (nt,ktf)                                       &
        ,intent (out  )                   ::                           &
         tc, tr_pw
     real,    dimension (ktf)                                          &
        ,intent (inout)                   ::                           &
         frac_act,wtbar,wcbar
!
!  local variables in this routine
!
     integer                              ::                           &
        k,n

     real :: airm(ktf), ewash(ktf)
     real :: washr, dt, hplus, pdry, qtot, pwat
     real :: prate, removac, fact3

!--- moist static energy inside cloud
!
     tc      = 0.
     tr_pw   = 0.
     hp_rain = 0.

     if(ierr.ne.0)return

     do n = 1, nt
        tc(n,k22) = tr_cup(n,k22)
     enddo

     do k = kbcon+1, ktop+1

        ! mass flux kg_air/m^2/sec
        fact3 = zu(k) + 0.5 * zu(k-1) * cd(k)

        prate    = pw(k)            ! precipitation rate [kg_w/m^2/sec]
        pwat     = prate / fact3    ! rain water [kg_w/kg_air]
        qtot     = qrc(k) + pwat    ! total water mixing ratio [kg_w/kg_air]

        washr    = pwat / qtot      ! washout ratio (fraction of water removed)
        ewash(k) = 1.0 - washr      ! fraction of water that remains in updraft

        wtbar(k) = qtot * dens(k-1)      ! tot water content [kg/m^3]
        wcbar(k) = wtbar(k) * fracliq(k) ! liq water content [kg_w/m^3]

        airm (k) = fact3 * inv_mwairkg

        ! fraction of aerosol activated in mixed clouds, from
        ! Verheggen et al. 2007 JGR
        if (fracliq(k) < .999) then
           frac_act(k) = 0.051 + 0.949 * exp(6.*(fracliq(k) - 1.0))
        else
           frac_act(k) = 1.0
        endif

        ! limit aerosol activation for small water content?
        if (wtbar(k) < 1.e-5) then
           frac_act(k) = min(frac_act(k), sqrt(wtbar(k) / 1.e-5))
        endif

     enddo

     do k = k22+1, ktop+1

        do n = 1, nt
           tc(n,k) = tc(n,k-1) * fact1(k) + tr(n,k-1) * fact2(k)
        enddo

        if (k > kbcon .and. wtbar(k) > 1.e-9) then

           dt = dz_cup(k-1) / zu(k-1)

           ! IN-CLOUD AITKEN SCAVENGING;
           ! SCAVENGED MASS ADDED TO ACCUMULATION MODE, NO DEPOSITION

           call aitken_incloud_scav( tc(:,k), wtbar(k), dens(k-1), to(k-1), &
                                     po(k-1), dt )

           if (wcbar(k) > 1.e-6) then

              pdry = dens(k-1) * to(k-1) * rdry
              removac = 0.0

              ! AQUEOUS CHEMISTRY
              call aq_map( k, iw, wcbar(k), to(k-1), pdry, dt, tc(:,k), airm(k), &
                           frac_act(k), tr_pw(:,k-1), hplus, removac, ewash(k), dark )

              hp_rain = hp_rain + hplus * pw(k)

              ! IN-CLOUD RAINOUT OF GASES NOT INVOLVED IN AQUEOUS (Henry's Law)
              call gas_incloud_rainout_skipaq( ewash(k), wcbar(k), to(k-1), hplus, &
                                               airm(k), tc(:,k), tr_pw(:,k-1) )

           else

              ! Typical H+ value for clean rain
              hplus   = 2.e-6
              hp_rain = hp_rain + hplus * pw(k)

              ! IN-CLOUD RAINOUT OF GASES (Henry's Law)
              call gas_incloud_rainout( ewash(k), wcbar(k), to(k-1), hplus, &
                                        airm(k), tc(:,k), tr_pw(:,k-1) )

           endif

           ! IN-CLOUD RAINOUT / SNOWOUT OF AEROSOLS
           call aero_incloud_rainout( ewash(k), airm(k), frac_act(k), &
                                      tc(:,k), tr_pw(:,k-1) )
        endif

     enddo

     hp_rain = hp_rain / max(1.e-30, sum(pw(kbcon+1:ktop+1)))

   END SUBROUTINE cup_up_tracer



   SUBROUTINE cup_up_tracer_shal(k22,dz_cup,cd,zu,fact1,fact2,tr_cup,tr,tc,fracliq, &
              wtbar,wcbar,frac_act,kbcon,ktop,qrc,dens,po,to,ierr,dark,ktf,nt,iw)

     use consts_coms, only: rdry
     use const_data,  only: inv_mwairkg

     IMPLICIT NONE
!
!  on input
!
     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,nt,ierr,iw
     real,    dimension (nt,ktf)                                       &
        ,intent (in   )                   ::                           &
         tr,tr_cup
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
         cd,fact1,fact2,qrc,dens,po,to,dz_cup,zu,fracliq
     integer                                      &
        ,intent (in   )                   ::                           &
        k22,ktop,kbcon
!
! input and output
!
     real,    dimension (nt,ktf)                                       &
        ,intent (out  )                   ::                           &
         tc
     real,    dimension (ktf)                                          &
        ,intent (inout)                   ::                           &
         frac_act,wtbar,wcbar
     logical, intent(in) :: dark
!
!  local variables in this routine
!
     integer :: k,n
     real    :: remov(nt), airm(ktf)
     real    :: removac, fact3, hplus, pdry, dt

     tc = 0.
     if (ierr.ne.0) return

     do n = 1, nt
        tc(n,k22) = tr_cup(n,k22)
     enddo

     do k = kbcon+1, ktop+1

        ! mass flux kg_air/m^2/sec
        fact3    = zu(k) + 0.5 * zu(k-1) * cd(k)
        airm(k)  = fact3 * inv_mwairkg

        wtbar(k) = qrc(k) * dens(k-1)    ! tot water content [kg/m^3]
        wcbar(k) = wtbar(k) * fracliq(k) ! liq water content [kg_w/m^3]

        ! fraction of aerosol activated in mixed clouds, from
        ! Verheggen et al. 2007 JGR
        if (fracliq(k) < .999) then
           frac_act(k) = 0.051 + 0.949 * exp(6.*(fracliq(k) - 1.0))
        else
           frac_act(k) = 1.0
        endif

        ! limit aerosol activation for small water content?
        if (wtbar(k) < 1.e-5) then
           frac_act(k) = min(frac_act(k), sqrt(wtbar(k) / 1.e-5))
        endif

     enddo

     do k = k22+1, ktop+1

        do n = 1, nt
           tc(n,k) = tc(n,k-1) * fact1(k) + tr(n,k-1) * fact2(k)
        enddo

        if (k > kbcon .and. wtbar(k) > 1.e-9) then

           dt = dz_cup(k-1) / zu(k-1)

           ! IN-CLOUD AITKEN SCAVENGING;
           ! SCAVENGED MASS ADDED TO ACCUMULATION MODE, NO DEPOSITION

           if (wtbar(k) > 1.e-6) then
              call aitken_incloud_scav( tc(:,k), wtbar(k), dens(k-1), to(k-1), &
                                        po(k-1), dt )
           endif

           ! AQUEOUS CHEMISTRY

           if (wcbar(k) > 1.e-6) then
              pdry    = dens(k-1) * to(k-1) * rdry
              removac = 0.0
              call aq_map( k, iw, wcbar(k), to(k-1), pdry, dt, tc(:,k), airm(k), &
                           frac_act(k), remov, hplus, removac, 1.0, dark )
           endif
        endif

     enddo

   END SUBROUTINE cup_up_tracer_shal



   SUBROUTINE CUP_ENV_CLEV_TR(TR,TR_CUP,KTF,NT)
     implicit none

     integer, intent(in)  :: ktf, nt
     real,    intent(in)  :: tr    (nt,ktf)
     real,    intent(out) :: tr_cup(nt,ktf)
     integer              :: k, n

     do n = 1, nt
        tr_cup(n,1) = tr(n,1)
     enddo

     do k = 2, ktf
        do n = 1, nt
           tr_cup(n,k) = 0.5 * (tr(n,k-1) + tr(n,k))
        enddo
     enddo

   END SUBROUTINE CUP_ENV_CLEV_TR



   SUBROUTINE cup_dellas_tr(ierr,dz_cup,della,zu,tr,tr_cup,tc,dens,  &
                            ktop,k22,cd,en,tcd,zd,cdd,etd,jmin,ntr,ktf)

     implicit none

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,ntr
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ntr,ktf)                              &
        ,intent (inout  )                 ::                           &
        della
     real,    dimension (ktf)                              &
        ,intent (in  )                    ::                            &
        dz_cup,cd,zu,en,dens,zd,cdd,etd
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop,k22,jmin
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (ntr,ktf)                              &
        ,intent (in   )                   ::                           &
        tr,tr_cup,tc,tcd
!
!  local variables in this routine
!
     integer :: k, n
     real, dimension(ktf) :: rhodzm
     real :: subin, subdn, entra, detup, detdo

     della = 0.

     IF (ierr.ne.0) return

     do k = 1, ktop
        rhodzm(k) = 1.0 / (dens(k) * dz_cup(k))
    enddo

    do n = 1, ntr
       della(n,1) = zd(2) * (.5*(tcd(n,1)+tcd(n,2)) - tr_cup(n,2)) * rhodzm(1)
    enddo
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
     do k = 2, ktop

        subin = zu(k+1) - zd(k+1)
        subdn = zu(k)   - zd(k)
        entra = en(k+1) * zu(k) + etd(k) * zd(k+1)
        detup = cd(k+1) * zu(k)
        detdo = cdd(k)  * zd(k+1)

        do n = 1, ntr
           della(n,k) = ( detup * .5*(tc(n,k+1)+tc(n,k)) &
                        + detdo * .5*(tcd(n,k+1)+tcd(n,k)) &
                        - entra * tr(n,k) &
                        + subin * tr_cup(n,k+1) &
                        - subdn * tr_cup(n,k) &
                        ) * rhodzm(k)
        enddo
     enddo

     do n = 1, ntr

        ! Special at jmin
        della(n,jmin) = della(n,jmin) - zd(jmin) * tcd(n,jmin+1) * rhodzm(jmin)

        ! Special at k22
        della(n,k22-1) = della(n,k22-1) - zu(k22) * tc(n,k22) * rhodzm(k22-1)

     enddo

   END SUBROUTINE cup_dellas_tr


   SUBROUTINE cup_dellas_tr_nodd(ierr,dz_cup,della,zu,tr,tr_cup,tc,dens,  &
                                 ktop,k22,cd,en,ntr,ktf)

     implicit none

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,ntr
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ntr,ktf)                              &
        ,intent (inout  )                 ::                           &
        della
     real,    dimension (ktf)                              &
        ,intent (in  )                    ::                            &
        dz_cup,cd,zu,en,dens
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop,k22
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (ntr,ktf)                              &
        ,intent (in   )                   ::                           &
        tr,tr_cup,tc
!
!  local variables in this routine
!
     integer :: k, n
     real, dimension(ktf) :: rhodzm
     real :: subin, subdn, entra, detup

     della = 0.

     IF (ierr.ne.0) return

     do k = k22, ktop
        rhodzm(k) = 1.0 / (dens(k) * dz_cup(k))
     enddo
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
     do k = k22, ktop

        subin = zu(k+1)
        subdn = zu(k)
        entra = en(k+1) * zu(k)
        detup = cd(k+1) * zu(k)

        do n = 1, ntr
           della(n,k) = ( detup * .5*(tc(n,k+1)+tc(n,k)) &
                        - entra * tr(n,k) &
                        + subin * tr_cup(n,k+1) &
                        - subdn * tr_cup(n,k) &
                        ) * rhodzm(k)
        enddo
     enddo

   END SUBROUTINE cup_dellas_tr_nodd



   SUBROUTINE cup_dd_tracer( tcd,dz_cup,fact1,fact2,to,dens,po,cdd,entr,jmin, &
                             kbcon,ierr,tr,tr_cup,tr_pw,zd,pwd,pcpflx,tr_flx, &
                             f_act,f_liq,wcbar,hplus,ntr,ktf )

    use cgrid_spcs, only: ae_strt, ae_fini
    use aq_data,    only: lsrggas, ln2o5
    use const_data, only: mwairkg

    implicit none
!
!  on input
!

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,ntr
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (ntr,ktf)                              &
        ,intent (in   )                   ::                           &
         tr,tr_cup
     real,    dimension (ntr,ktf)                              &
        ,intent (inout)                   ::                           &
         tr_pw
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        dz_cup,fact1,fact2,cdd,zd,pwd,pcpflx,to,dens,po,f_act,f_liq,wcbar
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer                                      &
        ,intent (in   )                   ::                           &
        jmin,kbcon
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (ntr,ktf)                              &
        ,intent (out  )                   ::                           &
        tcd
     real,    dimension (ntr)                                  &
        ,intent (inout)                   ::                           &
        tr_flx
     real                                                      &
        ,intent (in   )                   ::                           &
        hplus
!
!  local variables in this routine
!
     integer                              ::                           &
        k,n

     real, dimension(jmin) :: fact3
     real, dimension(ntr) :: remov, dtr
     real :: pr_flx, pfact

     integer :: n_n2o5
     logical :: doentrain

     tcd = 0.0

     IF (ierr.ne.0) return

     n_n2o5 = lsrggas( ln2o5 )

     doentrain = (entr > 1.e-20)

     do k = 1, jmin
        fact3(k) = mwairkg / (zd(k) + 0.5 * zd(k+1) * cdd(k))
     enddo

     ! top boundary condition
     tcd(:,jmin+1) = tr_cup(:,jmin)

     do k = jmin, 1, -1

        ! Downdraft solution

        if ( (.not. doentrain) .or. (k == 1) .or. (k == jmin) ) then
           tcd(:,k) = tcd(:,k+1)
        else
           tcd(:,k) = tcd(:,k+1) * fact1(k) + tr(:,k) * fact2(k)
        endif

        ! Evaporation of rain releasing gas/aerosols

        if (pcpflx(k+1) > 1.e-15 .and. pwd(k) < -1.e-15) then
           pfact = min(pwd(k) / pcpflx(k+1), 0.0)

           if (pfact > -0.99) then

              do n = 1, ntr
                 dtr(n) = tr_flx(n) * pfact
                 if (n == n_n2o5 .or. (n >= ae_strt .and. n <= ae_fini)) then
                    dtr(n) = dtr(n) * 0.5
                 endif
              enddo

           else

              dtr = -tr_flx

           endif

           do n = 1, ntr
              tcd(n,k)  = tcd(n,k)  - dtr(n) * fact3(k)
              tr_flx(n) = tr_flx(n) + dtr(n)
           enddo

        endif

        ! Washout from precip falling from above

        pr_flx = pcpflx(k+1) + pwd(k)

        call precip_scav( tcd(:,k), dtr, remov, to(k), dens(k), po(k), &
                          dz_cup(k), f_act(k+1), hplus, pr_flx, tr_flx, &
                          f_liq(k+1), wcbar(k+1) )

        tcd(:,k) = tcd(:,k) - remov * fact3(k)

        if (k >= kbcon) then
           tr_flx = tr_flx + remov + tr_pw(:,k)
        else
           tr_flx = tr_flx + remov
        endif

     enddo

   END SUBROUTINE cup_dd_tracer


end MODULE module_cu_g3_aq
