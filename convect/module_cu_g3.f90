MODULE module_cu_g3

  use misc_coms,   only: io6
  use consts_coms, only: cp, alvl, rvap, grav, rocp

  implicit none

  real,    parameter :: xmbmax = 1.0
  real,    parameter :: xl     = alvl
  real,    parameter :: xlm    = 1.0 / alvl
  real,    parameter :: rv     = rvap
  real,    parameter :: g      = grav
  real,    parameter :: tcrit  = 258.0
  real,    parameter :: gamfac = xl * xl / (cp * rv)
  real,    parameter :: gocp   = g / cp
  real,    parameter :: one3   = 1.0 / 3.0
  real,    parameter :: cpm    = 1.0 / cp
  real,    parameter :: lambau = 2.0

  real,    parameter :: radius_deep = 10000.
  real,    parameter :: entr_deep = 0.24 / radius_deep
  real,    parameter :: detr_deep = 0.02 * entr_deep

  real,    parameter :: entr_uv = 0.24 / 150.
  real,    parameter :: detr_uv = entr_uv
  real,    parameter :: uvfac   = 0.5

  real,    parameter :: radius_shal =   200.
  real,    parameter :: entr_shal = 0.24 / radius_shal
  real,    parameter :: detr_shal = entr_shal

  private
  public :: entr_deep, detr_deep, radius_deep, entr_shal, detr_shal, &
            grell_driver, cup_up_nms, cup_up_en_cd, cup_dd_nms

CONTAINS

  SUBROUTINE grell_driver(iw, dtlong)

     use mem_turb,    only: kpblh, frac_sfc, wstar, ustar, pblh
     use mem_ijtabs,  only: itab_w
     use misc_coms,   only: nqparm
     use consts_coms, only: p00i, rocp, gravi
     use mem_grid,    only: mza, lpw, arw0, zm, zt, lsw, gdz_abov8, arw, &
                            volt, arw0i
     use mem_radiate, only: pbl_cld_forc
     use mem_basic,   only: theta, tair, press, rho, rr_v, ue, ve
     use mem_cuparm,  only: thsrc, rtsrc, conprr, kcutop, kcubot, cbmf, &
                            qwcon, iactcu, kudbot, cu_pwa, cu_pev, cddf, &
                            kddtop, kddmax, kddbot, kstabi, umsrc, vmsrc
     implicit none

     integer, intent(in) :: iw
     real,    intent(in) :: dtlong

     real    :: outqc  (mza) ! output deep conv cloud water tendency
     real    :: z1           ! surface elevation
     real    :: tn     (mza) ! input forced temperature
     real    :: qo     (mza) ! input forced mixing ratio
     real    :: zo     (mza) ! height
     real    :: po     (mza) ! input forced pressure
     real    :: dens   (mza) ! air density
     real    :: pre          ! output precipitation rate
     real    :: outt   (mza) ! output deep conv temp tendency
     real    :: outq   (mza) ! output deep conv water vapor tendency
     real    :: psur         ! input surface pressure
     real    :: outu   (mza) ! output deep conv u tendency
     real    :: outv   (mza) ! output deep conv v tendency
     real    :: us     (mza) ! input zonal wind
     real    :: vs     (mza) ! input meridional wind
!    real    :: tshall (mza)   ! input PBL forced temperature for shallow convection
!    real    :: qshall (mza)   ! input PBL forced water vapor for shallow convection
     integer :: kpbl           ! layer corresponding to PBL height
     real    :: outts  (mza)   ! shallow conv temp tendency
     real    :: outqs  (mza)   ! shallow conv water vapor tendency
     real    :: tscl_kf        ! shallow convection timescale
     integer :: k23            ! shallow convection updraft source level
     integer :: kbcon3         ! shallow convection LCL
     integer :: ktop3          ! shallow convection cloud top
     real    :: xmb3           ! shallow convection mass flux
     integer :: k22            ! deep convectin updraft source level
     integer :: kbcon          ! deep convection LCL
     integer :: ktop           ! deep convection cloud top
     integer :: kstab
     integer :: jmin           ! downdraft originating level
     integer :: kdet           ! downdraft detraining level
     real    :: xmb            ! deep convection mass flux
     real    :: cupclw (mza)   ! deep convection cloud water
     real    :: edt_out
     real    :: cupclws(mza)   ! shallow convection cloud water
     real    :: zu     (mza)   ! normalized updraft mass flux
     real    :: zd     (mza)   ! normalized downdraft mass flus
     real    :: area   (mza)
     real    :: voa    (mza)
     integer :: ishallow_g3    ! flag = 1 activates shallow convection
     integer :: ideep_g3       ! flag = 1 activates deep convection
     integer :: ktf            ! ending   k level
     integer :: ierr           ! deep convection activation flag
     integer :: ierrd          ! deep convection downdraft activation flag
     integer :: ierrs          ! shallow convection activation flag
     real    :: a0
     real    :: wm3

     real :: pcpflx(mza)
     real :: pwa   (mza)
     real :: pev   (mza)

     ! local variables
     integer :: k, ka, kc, ks
     real    :: exner(mza)
     real    :: wc3

     pwa = 0.
     pev = 0.
     ka  = lpw(iw)

     ! Go no higher than 50mb for convective calculations to prevent any
     ! problems when esat gets near ambient pressure
     do k = ka, mza-1
        if (press(k,iw) < 50.e2) exit
     enddo
     ktf = k - ka + 1    ! number of ATM levels to process

     ! Level of PBL height
     kpbl= min(kpblh(iw) - ka + 1, ktf-1)

     ! A. Betts for shallow convection: suggestion for the KF timescale < DELTAX / 25 m/s
     tscl_kf = sqrt(arw0(iw)) / 25.

     psur  = 0.01*(press(ka,iw) + gdz_abov8(ka-1) * rho(ka,iw))
     z1    = zm(ka-1)
     a0    = arw0(iw)

     wc3   = grav * pblh(iw) * pbl_cld_forc(iw) / tair(kpbl,iw)
     wm3   = (wstar(iw)**3 + 6. * ustar(iw)**3 + wc3) ** one3

     ! Loop over T levels
     do kc = 1, ktf
        k  = kc + ka - 1

        exner(k) = tair(k,iw) / theta(k,iw)
        dens(kc) = real(rho(k,iw))
        area(kc) = min(arw(k-1,iw), a0)
        voa (kc) = volt(k,iw) * arw0i(iw)

        ! Current temp, water vapor, and pressure
        tn(kc) = max( tair(k,iw), 200.0 )
        qo(kc) = max( rr_v(k,iw), 1.e-8 )
        qo(kc) = qo(kc) / (1.0 + qo(kc))
        po(kc) = 0.01*press(k,iw)
        zo(kc) = zt(k)

        ! U and V winds
        us(kc) = ue(k,iw)
        vs(kc) = ve(k,iw)

!       ! Shallow convection uses current T and Q plus PBL tendencies
!       tshall(kc) = tn(kc) + fthpbl(k,iw) * dtlong * exner(k)
!       qshall(kc) = qo(kc) + fqtpbl(k,iw) * dtlong
!       qshall(kc) = qshall(kc) / (1.0 + qshall(kc))
     enddo

     area(1) = area(2)

     if (nqparm( itab_w(iw)%mrlw ) == 5) then
        ideep_g3 = -1
     else
        ideep_g3 = 1
     endif

     ishallow_g3 = 1
     ktop  = 0
     ktop3 = 0

     call CUP_enss_3d(iw,OUTQC,Z1,ideep_g3,                            &
              TN,QO,PO,zo,PRE,OUTT,OUTQ,OUTU,OUTV,PSUR,US,VS,          &
              kpbl,outts,outqs,tscl_kf,                                &
              k23,kbcon3,ktop3,xmb3,wm3,                               &
              k22,xmb,jmin,kdet,zu,zd,                                 &
              kbcon,ktop,kstab,cupclw,pwa,pev,                         &
              edt_out,dens,voa,area,a0,                                &
              ishallow_g3,cupclws,ktf,ierr,ierrd,ierrs                 )

     if (ierr==0) then

        ! Deep convecton is active

        kcutop(iw) = ktop  + ka - 1
        kcubot(iw) = kbcon + ka - 1
        kudbot(iw) = k22   + ka - 1
        kstabi(iw) = kstab + ka - 1
        cbmf  (iw) = xmb

        ! Save type of convection in iactcu

        if (ierrd==0) then
           iactcu(iw) = 1
        else
           iactcu(iw) = 2
        endif

        if (ierrd==0) then
           kddtop(iw) = jmin + ka - 1
           kddmax(iw) = kdet + ka - 1
           kddbot(iw) = ka
           cddf  (iw) = edt_out * xmb
        endif

        do kc = 1, ktop
           k  = kc + ka - 1

           ! Total water tendency
           rtsrc(k,iw) = outq(kc) * dens(kc)

           ! Heat (temperature) tendency
           thsrc(k,iw) = outt(kc) * dens(kc)

           ! Cloud water
           qwcon(k,iw) = cupclw(kc)

           ! Momentum tendency
           umsrc(k,iw) = outu(kc) * dens(kc)
           vmsrc(k,iw) = outv(kc) * dens(kc)

           ! Precipitation formation rate (kg/m^2/s)
           if (allocated(cu_pwa)) then
              cu_pwa(k,iw) = pwa(kc)
           endif

           ! Precipitation evaporation rate (kg/m^2/s)
           if (allocated(cu_pev) .and. ierrd==0) then
              cu_pev(k,iw) = pev(kc)
           endif

        enddo

        ! Convective precipitation flux
        pcpflx( ktop+1 ) = 0.0
        do kc = ktop, 1, -1
           pcpflx(kc) = pcpflx(kc+1) + pwa(kc) + pev(kc)
        enddo

        ! Surface precipitation
        do ks = 1, min(lsw(iw), ktop)
           conprr(iw) = conprr(iw) + frac_sfc(ks,iw) * pcpflx(ks)
        enddo

     else if (ishallow_g3 == 1 .and. ierrs == 0) then

        ! Shallow convection is active

        kcutop(iw) = ktop3  + ka - 1
        kcubot(iw) = kbcon3 + ka - 1
        kudbot(iw) = k23    + ka - 1

        iactcu(iw) = 3
        cbmf  (iw) = xmb3

        do kc = k23-1, ktop3
           k  = kc + ka - 1
           rtsrc(k,iw) = outqs(kc) * dens(kc)
           thsrc(k,iw) = outts(kc) * dens(kc)
           qwcon(k,iw) = cupclws(kc)
        enddo

     endif

   END SUBROUTINE grell_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE CUP_enss_3d(iw,OUTQC,Z1,ideep,                           &
              TN,QO,P,z,PRE,OUTT,OUTQ,OUTU,OUTV,PSUR,US,VS,            &
              kpbl,outts,outqs,tscl_kf,                                &
              k23,kbcon3,ktop3,xmb3,wm3,                               &
              k22,xmb,jmin,kdet,zu,zd,                                 &
              kbcon,ktop,kstabi,cupclw,pwa_out,pev_out,                &
              edt_out,dens,voa,arw,arw0,                               &
              ishallow_g3,cupclws,ktf,ierr_out,ierrd_out,ierrs_out     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,iw
     integer, intent (in   )              ::                           &
        ishallow_g3, ideep
     real                                      &
        ,intent (out )                    ::                           &
               edt_out
     real                                      &
        ,intent (in   )                   ::                           &
               wm3
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (ktf)                              &
        ,intent (inout  )                   ::                         &
        OUTT,OUTQ,OUTQC,cupclw,outts,outqs,                &
        cupclws,zu,zd,pwa_out,pev_out,outu,outv
     real                                      &
        ,intent (out  )                   ::                           &
        pre,xmb3,xmb
     integer                                   &
        ,intent (out  )                   ::                           &
        k22,kbcon,ktop,kstabi,k23,kbcon3,ktop3,jmin,kdet,              &
        ierr_out,ierrd_out,ierrs_out
     integer                                   &
        ,intent (in   )                   ::                            &
        kpbl
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        P,z,US,VS,tn,dens,voa,arw
     real,    dimension (ktf)                              &
        ,intent (inout)                   ::                           &
         QO
     real                                         &
        ,intent (in   )                   ::                           &
        Z1,PSUR,arw0
     real                                                              &
        ,intent (in   )                   ::                           &
        tscl_kf
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! entr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! entr_rate = entrainment rate
     real,    dimension (ktf) ::                                       &
!        he3,hes3,qes3,                                                 &
!        qes3_cup,q3_cup,he3_cup,hes3_cup,t3_cup,                       &
        xhe3,xhes3,xqes3,xt3,xq3,                                      &
        xqes3_cup,xq3_cup,xhe3_cup,xhes3_cup,                          &
        xt3_cup,xhc3,                                                  &
        qc3,hc3,qrc3,zu3,cd3,en3,DELLAH3,DELLAQ3,DELLAT3

     real,    dimension (ktf) ::                                       &
        heo,heso,qeso,                                                 &
        xhe,xhes,xqes,xt,xq,                                           &

        z_cup,p_cup,                                                   &
        qeso_cup,qo_cup,heo_cup,heso_cup,dz,dz_cup,                    &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xt_cup,                       &
        qco,pwdo,pwo,hcdo,qcdo,hco,qrco,xhc,                           &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,en,cdd,etd,DELLAH,DELLAQ,DELLAT,DELLAQC,rhodz_i,            &
        facu1,facu2,facd1,facd2,cdu,enu,facu1u,facu2u,facd1u,facd2u,   &
        cddu,etdu,ue_cup,ve_cup,uc,vc,ucd,vcd,dellau,dellav

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)

     real ::                                                           &
       aa3_0,aa3,xaa3,                                                 &
       edto,AA1,XAA0,                                                  &
       PWAVO,PWEVO,cap_max,                                            &
       cap_max_increment,cap_max3
     integer ::                                                        &
       kzdown,ierr,ierr2,ierr3,KBMAX,ierr5,ierrd

     integer                              ::                           &
       K,kbot,kfull
     real                                 ::                           &
      day,entr_rate,rad,                                               &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,dh

     logical :: keep_going, pbl_linked
     real :: xkshal,arw0m,detr_rate

      day=86400.
      xmb=0.
      xmb3=0.
!
!--- maximum depth (mb) of capping inversion that convection is
!--- allowed to overcome (smaller cap_max means less convection)
!
      cap_max = 25.0
      cap_max_increment=25.
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
      rad = sqrt(arw0)
      if (rad < radius_deep) then
         entr_rate = 0.24 / rad
         detr_rate = detr_deep + (entr_rate - entr_deep)
      else
         entr_rate  = entr_deep
         detr_rate  = detr_deep
      endif
!
!--- initial detrainmentrates
!
      cupclw=0.
      cupclws=0.
      cdd=0.
      etd=0.
      hcdo=0.
      dellaqc=0.
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax=.9
      edtmin=.1
!
!--- minimum depth (m) precipitating clouds must have
!
!     depth_min=750.
!     depth_min=1000.
!     depth_min=1250.
      depth_min=1500.

      kbmax=3
      jmin=0
      kdet=0
      aa3_0=0.
      aa1=0.
      aa3=0.
      edto=0.
      IERR=0
      IERR2=0
      IERR3=0
      IERRD=0
!
!---  Flag to disable deep convection for this grid cell
!
      if (ideep <= 0) then
         ierr=1
      endif
!
!---  Special for shaved cells
!
      arw0m = nearest(arw0,-1.)
      do k = 1, ktf-1
         if (arw(k) >= arw0m) exit
      enddo
      kfull = k

      call cup_up_en_cd(en,cd,facu1,facu2,voa,arw,entr_rate,detr_rate,kfull,ierr,ktf)
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax = 3500.
!     zkbmax = 4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1250.
!
!--- calculate moist static energy, heights, qes
!
      call cup_env_heights(z,p,z1,psur,p_cup,z_cup,dz,dz_cup,ktf)

      do k=3,ktf-1
         if(z_cup(k).gt.zkbmax+z1)exit
      enddo
      kbmax=max(k,kpbl+1)

      call cup_env(z,qeso,heo,heso,tn,qo,p,ierr,ktf)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(tn,qeso,qo,heo,heso,qeso_cup,qo_cup, &
           heo_cup,heso_cup,tn_cup,ierr,1,ktf,ktf)

      call cup_env_clev_uv(us,vs,ue_cup,ve_cup,ierr,1,ktf,ktf)
!
!----  DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      call cup_maximi(heo_cup,4,kbmax,k22,ierr,ktf)
      if (k22 >= kbmax) ierr=2
!
!---- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(k22,kbcon,heo_cup,heso_cup,ierr, &
           kbmax,ktop,p_cup,cap_max,ktf)

      ! ktop is the tentative cloud top without entrainment.
      ! It will be recomputed later and may be slightly lower

      IF (ierr.eq.0) THEN
         IF ((z_cup(KTOP)-z_cup(KBCON)<depth_min) .or. (ktop - kbcon < 3)) ierr=6
      ENDIF
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,kbcon,ktop+1,kstabi,ierr,ktf)

      IF(ierr.eq.0)THEN
         if (kstabi <= ktop) then
            do k = kstabi, ktop
               cd   (k) = min(cd(k-1)/voa(k-2)+.15*entr_rate,1.5*entr_rate)*voa(k-1)
               facu2(k) = en(k)  / (1.0 + en(k) - 0.5*cd(k))
               facu1(k) = 1.0 - facu2(k)
            enddo
         ENDIF
      ENDIF
!
!--- normalized updraft mass flux profile
!
      call cup_up_nms(zu,cd,en,ktop,ierr,k22,facu1,facu2,arw,arw0,ktf)
!
!--- calculate incloud moist static energy and ktop
!
      call cup_up_he(k22,heo_cup,hco,ktop,ierr,heo,facu1,facu2,ktf)

      call cup_up_ktop(cd,en,hco,zu,kbcon,ktop,ierr,heso_cup,facu1,facu2,ktf)

      ! Deep/mid convection must have sufficient depth
      IF (ierr.eq.0) THEN
         IF (ktop-kbcon<3 .or. z_cup(KTOP)-z_cup(KBCON)<depth_min) ierr=6
      ENDIF
!
!--- increased entrainment/detrainment for momentum
!
      if (ierr.eq.0) then

         cdu(k22) = cd(k22)
         enu(k22) = en(k22)

         facu2u(k22) = facu2(k22)
         facu1u(k22) = facu1(k22)

         cdu(ktop+1) = cd(ktop+1)
         enu(ktop+1) = en(ktop+1)

         facu2u(ktop+1) = facu2(ktop+1)
         facu1u(ktop+1) = facu1(ktop+1)

         cdu(1:k22-1) = 0.
         enu(1:k22-1) = 0.

         do k = k22+1, ktop
            cdu   (k) = cd(k) + detr_uv * voa(k-1)
            enu   (k) = en(k) + entr_uv * voa(k-1)
            facu2u(k) = enu(k)  / (1.0 + enu(k) - 0.5*cdu(k))
            facu1u(k) = 1.0 - facu2u(k)
         enddo

      endif

      call cup_up_uv(k22,ue_cup,ve_cup,uc,vc,ktop,ierr,us,vs,facu1u,facu2u,ktf)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa1,dz,zu,hco,heso_cup,qeso_cup,tn_cup,kbcon,ktop,ierr,ktf)

      if(ierr.eq.0)then
         if (aa1 <= 15.0) ierr=17
      endif
!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture(ierr,voa,qco,qrco,pwo,pwavo,kbcon,ktop,cd, &
           hco,heso_cup,qo,qeso_cup,zu,facu1,facu2,k22,qo_cup,tn_cup,ktf)

      if (ierr.eq.0) then
         if (pwavo < 1.e-9) ierr = 18
      endif

      if (ierr.eq.0) then
         do k = kbcon, ktop
            cupclw(k) = qrco(k+1)
         enddo
      endif
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      if(ierr.eq.0)then
         zktop=(z_cup(ktop)-z1)*.6
         zktop=min(zktop,zcutdown)+z1
         do k=1,ktf-3
            if(z_cup(k).gt.zktop)exit
         enddo
         kzdown=k

         call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr,ktf)
      endif
!
!--- check whether downdraft would have negative buoyancy, if there where
!--- no entrainment/detrainment
!
      ierrd = ierr
      IF(ierrd.eq.0)THEN

         keep_going = .TRUE.
         do while ( keep_going )
            keep_going = .FALSE.

            ! find first negatively buoyant layer below jmin
            if (heso_cup(jmin).ge.heso_cup(jmin-1)) then
               do k=jmin-1,4,-1
                  if (heso_cup(k).lt.heso_cup(k-1)) exit
               enddo
               jmin = k
            endif

            if (jmin.lt.4) then
               ierrd = 9
               exit
            endif

            dh=0.
            do k=jmin-1,1,-1
               dh=dh+dz_cup(k)*(heso_cup(jmin)-heso_cup(k))
               if(dh.gt.0.)then
                  jmin=jmin-1
                  if (jmin .gt. 3) then
                     keep_going = .TRUE.
                     exit
                  else
                     ierrd = 9
                     exit
                  endif
               endif
            enddo
         enddo

      ENDIF
!
!--- level where detrainment for downdraft starts
!
      if (ierrd.eq.0) then
         if (z_cup(jmin) > z_detr+z1) then
            do k=jmin-1,2,-1
               if(z_cup(k) <= z_detr+z1) exit
            enddo
            kdet = k
         else
            kdet = jmin - 1
         endif
      endif
!
!--- normalized downdraft mass flux profile, also work on downdraft entrainment
!--- and detrainment in this routine
!
      call cup_dd_nms(zd,voa,cdd,etd,facd1,facd2,entr_rate, &
                      jmin,ierrd,arw,arw0,kdet,kfull,ktf)
!
!--- downdraft moist static energy
!
      call cup_dd_he(heso_cup,hcdo,facd1,facd2,entr_rate,jmin,ierrd,heo,ktf)
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture(zd,hcdo,heso_cup,qcdo,qeso_cup,tn_cup, &
           pwdo,qo_cup,cdd,facd1,facd2,entr_rate,jmin,ierrd,pwevo,qo,ktf)
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierrd,us,vs,z,ktop,kbcon,edto,p,pwavo, &
           pwevo,edtmax,edtmin,ktf)

      if (ierrd.eq.0) then
         edt_out = edto
         do k = 1, jmin
            zd  (k) = zd  (k) * edto
            pwdo(k) = pwdo(k) * edto
         enddo
         pwevo = pwevo * edto
      else
         zd   = 0.0
         edto = 0.0
         jmin = 0
         kdet = 0
      endif
!
!--- increased entrainment/detrainment for momentum
!
      if (ierrd.eq.0) then

         cddu(1) = cdd(1)
         etdu(1) = etd(1)

         facd2u(1) = facd2(1)
         facd1u(1) = facd1(1)

         cddu(jmin) = cdd(jmin)
         etdu(jmin) = etd(jmin)

         facd2u(jmin) = facd2(jmin)
         facd1u(jmin) = facd1(jmin)

         cddu(jmin+1:ktop) = 0.
         etdu(jmin+1:ktop) = 0.

         if (entr_uv > 1.e-20) then
            do k = 2, jmin-1
               cddu  (k) = cdd(k) + detr_uv * voa(k)
               etdu  (k) = etd(k) + entr_uv * voa(k)
               facd2u(k) = etdu(k)  / (1.0 + etdu(k) - 0.5*cddu(k))
               facd1u(k) = 1.0 - facd2u(k)
            enddo
         endif

         call cup_dd_uv(ue_cup,ve_cup,ucd,vcd,facd1u,facd2u,entr_uv, &
                        jmin,ierr,us,vs,ktf)
      endif
!
!--- change per unit mass that a model cloud would modify the environment
!
      if (ierrd == 0) then

         if (ierr.eq.0) then
            kbot = 1
            do k = kbot, ktop
               rhodz_i(k) = 1.0 / (dens(k) * voa(k))
            enddo
         endif

         call cup_dellas_3d(ierr,rhodz_i,hcdo,zd,cdd, &
              etd,heo,dellah,zu,en,                   &
              cd,hco,ktop,k22,jmin,heo_cup,ktf        )

         call cup_dellas_3d(ierr,rhodz_i,qcdo,zd,cdd, &
              etd,qo,dellaq,zu,en,                    &
              cd,qco,ktop,k22,jmin,qo_cup,ktf         )

         call cup_dellas_3d(ierr,rhodz_i,ucd,zd,cddu, &
              etdu,us,dellau,zu,enu,                  &
              cdu,uc,ktop,k22,jmin,ue_cup,ktf         )

         call cup_dellas_3d(ierr,rhodz_i,vcd,zd,cddu, &
              etdu,vs,dellav,zu,enu,                  &
              cdu,vc,ktop,k22,jmin,ve_cup,ktf         )

      else

         if (ierr.eq.0) then
            kbot = k22
            do k = kbot, ktop
               rhodz_i(k) = 1.0 / (dens(k) * voa(k))
            enddo
         endif

         call cup_dellas_nodd(ierr,rhodz_i,heo,dellah,zu, &
              cd,hco,ktop,k22,en,heo_cup,ktf)

         call cup_dellas_nodd(ierr,rhodz_i,qo,dellaq,zu, &
              cd,qco,ktop,k22,en,qo_cup,ktf)

         call cup_dellas_nodd(ierr,rhodz_i,us,dellau,zu, &
              cdu,uc,ktop,k22,enu,ue_cup,ktf)

         call cup_dellas_nodd(ierr,rhodz_i,vs,dellav,zu, &
              cdu,vc,ktop,k22,enu,ve_cup,ktf)

      endif
!
!-- take out cloud liquid water for detrainment
!
      dellaqc = 0.
      if(ierr.eq.0)then
         do k = kbcon, ktop
            dellaqc(k) = cd(k+1) * zu(k) * 0.5*(qrco(k)+qrco(k+1)) * rhodz_i(k)
         enddo
      endif
!
!--- using dellas, calculate changed environmental profiles
!
      dellat=0.
      if(ierr.eq.0)then

         do k = kbot, ktop
            DELLAT(K) = (DELLAH(K)-xl*DELLAQ(K)) * cpm
         enddo

         do k = 1, ktop
            XHE(K) = DELLAH(K) + HEO(K)
            XQ (K) = DELLAQ(K) + QO (K)
            XT (K) = DELLAT(K) + TN (K)
            IF(XQ(K).LT.1.E-08)XQ(K)=1.E-08
         enddo
      endif
!
!--- calculate moist static energy, qes
!
      call cup_env2(xqes,xhe,xhes,xt,xq,p,k22-1,ktop,ierr,ktf)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xt_cup,ierr,k22,ktop,ktf)
!
!--- moist static energy inside cloud
!
      call cup_up_he(k22,xhe_cup,xhc,ktop,ierr,xhe,facu1,facu2,ktf)
!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,dz,zu,xhc,xhes_cup,xqes_cup,xt_cup,kbcon,ktop,ierr,ktf)

      call cup_forcing_ens_3d(aa1,xaa0,ierr,xmb)
!
!--- FEEDBACK
!
      call cup_output_ens_3d(ierr,dellat,dellaq,        &
           dellaqc,outt,outq,outqc,pwavo,pwevo,         &
           pwa_out,pev_out,pre,xmb,ktop,pwo,pwdo,       &
           dellau,dellav,outu,outv,ktf)

      PRE=MAX(PRE,0.)

      ierr_out  = ierr
      ierrd_out = ierr + ierrd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    NEXT section for shallow convection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ierrs_out = 999
      xmb3      = 0.0
      ierr5     = 0

      outts = 0.
      outqs = 0.

      if (ishallow_g3 .eq. 0 .or. ierr .eq. 0) return

      cap_max3 = 25.

!     ! increase cap_max for unstable stratification
!     cap_max3 = cap_max3 + cap_max_increment * min(wstar,1.0)

      call cup_up_en_cd(en3,cd3,facu1,facu2,voa,arw,entr_shal,detr_shal,kfull,ierr5,ktf)

      CALL cup_MAXIMI(HEO_CUP,3,kbmax,K23,ierr5,ktf)

      ! Iterate upward from K23 to find if parcel reaches LCL
      call cup_kbcon(k23,kbcon3,heo_cup,heso_cup,ierr5, &
           kbmax,ktop3,p_cup,cap_max3,ktf)

      ! Does parcel originate from PBL?
      pbl_linked = .false.
      if (ierr5==0) then
         if (k23.le.kpbl .and. kpbl.gt.4) pbl_linked = .true.
      endif

      call cup_up_nms(zu3,cd3,en3,ktop3,ierr5,k23,facu1,facu2,arw,arw0,ktf)

      call cup_up_he(k23,heo_cup,hc3,ktop3,ierr5,heo,facu1,facu2,ktf)

      call cup_up_ktop(cd3,en3,hc3,zu3,kbcon3,ktop3,ierr5, &
           heso_cup,facu1,facu2,ktf)

      call cup_up_moisture_shall(ierr5,qc3,qrc3,kbcon3,ktop3,hc3, &
           heso_cup,tn_cup,qo,qo_cup,qeso_cup,facu1,facu2,k23,ktf)

      if (ierr5.eq.0) then
         do k=kbcon3,ktop3
            cupclws(k)=qrc3(k+1)
         enddo
      endif

      IF (.NOT. PBL_LINKED) THEN

         ! Calculate cloud work functions
         call cup_up_aa0(aa3,dz,zu3,hc3,heso_cup,qeso_cup,tn_cup,kbcon3,ktop3,ierr5,ktf)

         if(ierr5.eq.0)then
            if (aa3 <= 2.0) ierr5=17
         endif

      ENDIF
!
!---  now what is necessary for feedbacks
!
      if (ierr5.eq.0) then
         do k = k23, ktop3
            rhodz_i(k) = 1.0 / (dens(k) * voa(k))
         enddo
      endif

      call cup_dellas_nodd(ierr5,rhodz_i,heo,dellah3,zu3, &
           cd3,hc3,ktop3,k23,en3,heo_cup,ktf)

      call cup_dellas_nodd(ierr5,rhodz_i,qo,dellaq3,zu3, &
           cd3,qc3,ktop3,k23,en3,qo_cup,ktf)

      dellat3 = 0.
      if(ierr5.eq.0)then
         do k=k23,ktop3
            DELLAT3(K)=(DELLAH3(K)-xl*DELLAQ3(K)) * cpm
         enddo
      endif

      IF (.NOT. PBL_LINKED) THEN

         if (ierr5.eq.0) then
            do k=k23-1,ktop3
               XHE3(K) = DELLAH3(K) + HEO(K)
               XQ3 (K) = DELLAQ3(K) + QO (K)
               XT3 (K) = DELLAT3(K) + TN (K)
               XQ3 (K) = MAX(XQ3(K),1.E-08)
            enddo
         endif
!
!--- calculate moist static energy, heights, qes
!
         call cup_env2(xqes3,xhe3,xhes3,xt3,xq3,p,k23-1,ktop3,ierr5,ktf)
!
!--- environmental values on cloud levels
!
         call cup_env_clev(xt3,xqes3,xq3,xhe3,xhes3,xqes3_cup,xq3_cup, &
              xhe3_cup,xhes3_cup,xt3_cup,ierr5,k23,ktop3,ktf)
!
!--- moist static energy inside cloud
!
         call cup_up_he(k23,xhe3_cup,xhc3,ktop3,ierr5,xhe3,facu1,facu2,ktf)
!
!--- updated cloud work function
!
         call cup_up_aa0(xaa3,dz_cup,zu3,xhc3,xhes3_cup,xqes3_cup,xt3_cup,kbcon3,&
              ktop3,ierr5,ktf)
!
!--- shallow forcing
!
         if (ierr5.eq.0) then
            xkshal = max(aa3-xaa3, 1.e-2)
            xmb3   = max(aa3/(xkshal*tscl_KF), 0.)
         endif

      else  ! not PBL_LINKED

         if (ierr5.eq.0) then
            xmb3 = max(0.03 * wm3 * dens(kbcon3), 0.0)
         endif

      endif

      if(ierr5.eq.0) then
         if(xmb3.le.0.)ierr5=22
      endif

      ierrs_out = ierr5

      if(ierr5.ne.0)then

         k23=0
         kbcon3=0
         ktop3=0

      else

         xmb3 = min(xmb3,0.1)

!
! got the mass flux, sanity check, first for heating rates
!
!!         trash=0.
!!         do k=23,ktop3
!!            trash=max(trash,86400.*dellat3(k)*xmb3)
!!         enddo
!!         if(trash.gt.150.)xmb3=xmb3*150./trash
!
! sanity check on moisture tendencies: do not allow anything that may allow neg tendencies
!
!!         do k=23,ktop3
!!            trash=q(k)+dellaq3(k)*xmb3*dtime
!!            if(trash.lt.1.e-12)then
!!               ! max allowable tendency over tendency that would lead to too small mix ratios
!!               trash = (1.e-12-q(k)) &
!!                     / (dellaq3(k)*xmb3*dtime)
!!               trash=max(0.,trash)
!!               trash=min(1.,trash)
!!               xmb3=trash*xmb3
!!            endif
!!         enddo
!
! final tendencies
!
         do k=k23,ktop3
            outts(k)=dellat3(k)*xmb3
            outqs(k)=dellaq3(k)*xmb3
         enddo

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     done shallow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE CUP_enss_3d



   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,ktf                         )

     IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        us,vs,z,p
     real                                                      &
        ,intent (out  )                   ::                           &
        edt
     real                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev
     real                                                              &
        ,intent (in   )                   ::                           &
        edtmax,edtmin
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     integer :: kk
     real    :: pef,pefb,prezk,zkbc,vws,vshear,edtmn2
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
     edt=0.

     IF (ierr.eq.0) then

!--- calculate an average wind shear over the depth of the cloud
        vws = 0.
        do kk = kbcon, ktop
           vws = vws + &
              ( abs((us(kk+1)-us(kk))/(z(kk+1)-z(kk))) &
              + abs((vs(kk+1)-vs(kk))/(z(kk+1)-z(kk)))) * &
              (p(kk) - p(kk+1))
        enddo
        vshear = 1.e3 * vws / (p(kbcon) - p(ktop+1))

        pef=1.591-.639*VSHEAR+.0953*VSHEAR**2 &
            -.00496*VSHEAR**3
        if(pef.gt.1.)pef=1.
        if(pef.lt.0.)pef=0.
!
!--- cloud base precip efficiency
!
        zkbc=z(kbcon)*3.281e-3
        prezk=.02
        if(zkbc.gt.3.)then
           prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
                *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
        endif
        if(zkbc.gt.25.)then
           prezk=2.4
        endif
        pefb=1./(1.+prezk)
        if(pefb.gt.1.)pefb=1.
        if(pefb.lt.0.)pefb=0.

        EDT=1.-.5*(pefb+pef)
!--- edt here is 1-precipeff!

        edt=edt*pwav/max(-pwev,1.e-12)

        edtmn2 = min(edtmin, .9999 * pwav / max(-pwev,1.e-12))
        IF(EDT.GT.edtmax)EDT=edtmax
        IF(EDT.LT.edtmn2)EDT=edtmn2
     endif

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_he(hes_cup,hcd,facd1,facd2,entr,jmin,ierr,he,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! cdd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        he,hes_cup,facd1,facd2
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        hcd
!
!  local variables in this routine
!
     integer                              ::                           &
        k

     hcd = 0.
     IF (ierr.ne.0) return

     if (entr <= 1.e-20) then

        do k = 1, jmin+1
           hcd(k) = hes_cup(jmin)
        enddo

     else

        hcd(jmin+1) = hes_cup(jmin)
        hcd(jmin)   = hes_cup(jmin)

        do k = jmin-1, 2, -1
           hcd(k) = facd1(k) * hcd(k+1) + facd2(k) * he(k)
        enddo

        hcd(1) = hcd(2)

     endif

   END SUBROUTINE cup_dd_he


   SUBROUTINE cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup,t_cup, &
              pwd,q_cup,cdd,facd1,facd2,entr,jmin,ierr,pwev,q,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
  ! cdd= detrainment function
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! entr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate

     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        zd,hes_cup,hcd,qes_cup,q_cup,t_cup,cdd,q,facd1,facd2
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        qcd,pwd
     real                                      &
        ,intent (out  )                   ::                           &
        pwev
!
!  local variables in this routine
!
     integer                              ::                           &
        k
     real                                 ::                           &
        dh,dqeva,gamma
     logical                              ::                           &
        noentrain
     real, dimension (jmin)               ::                           &
        qrcd,facd3

     pwev=0.
     qcd =0.
     pwd =0.

     IF(ierr.ne.0)return

     k = jmin
     qcd (k)  = q_cup(k)
     qrcd(k)  = qes_cup(k)
     pwd (k)  = min(0., (qcd(k) - qrcd(k)) * zd(k))
     pwev     = pwd (k)
     qcd(k+1) = qcd (k)
     qcd(k)   = qrcd(k)

     noentrain = (entr < 1.e-20)

     do k = 1, jmin-1
        dh       = max(hcd(k) - hes_cup(k), 0.)
        gamma    = gamfac * qes_cup(k) / (t_cup(k) * t_cup(k))
        qrcd (k) = qes_cup(k) + dh * gamma * xlm / (1.+gamma)
        facd3(k) = zd(k) + 0.5 * zd(k+1) * cdd(k)
     enddo

     do k = jmin-1,1,-1
        if (noentrain .or. k == 1) then
           qcd(k) = qcd(k+1)
        else
           qcd(k) = qcd(k+1) * facd1(k) + q(k) * facd2(k)
        endif

        dqeva  = min(qcd(k)-qrcd(k), 0.0)
        qcd(k) = qcd(k) - dqeva

        pwd(k) = dqeva * facd3(k)
        pwev   = pwev + pwd(k)
     enddo

   END SUBROUTINE cup_dd_moisture



   SUBROUTINE cup_dd_uv(ue_cup,ve_cup,ucd,vcd,facd1,facd2,entr, &
                        jmin,ierr,us,vs,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! cdd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        us,vs,ue_cup,ve_cup,facd1,facd2
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        ucd,vcd
!
!  local variables in this routine
!

     integer                              ::                           &
        k

     ucd = 0.
     vcd = 0.
     IF (ierr.ne.0) return

     if (entr <= 1.e-20) then

        do k = 1, jmin+1
           ucd(k) = ue_cup(jmin)
           vcd(k) = ve_cup(jmin)
        enddo

     else

        ucd(jmin+1) = ue_cup(jmin)
        ucd(jmin)   = ue_cup(jmin)

        vcd(jmin+1) = ve_cup(jmin)
        vcd(jmin)   = ve_cup(jmin)

        do k = jmin-1, 2, -1
           ucd(k) = ucd(k+1) * facd1(k) + us(k) * facd2(k)
           vcd(k) = vcd(k+1) * facd1(k) + vs(k) * facd2(k)
        enddo

        ucd(1) = ucd(2)
        vcd(1) = vcd(2)

     endif

   END SUBROUTINE cup_dd_uv



   SUBROUTINE cup_dd_nms(zd,dz_cup,cdd,etd,facd1,facd2,entr, &
                         jmin,ierr,arw,arw0,kdet,kfull,ktf)

     IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  ! z_cup = height of cloud model level
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining

     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        dz_cup,arw
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,arw0
     integer                                                           &
        ,intent (in   )                   ::                           &
        jmin,kdet,kfull
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer                                                           &
        ,intent (inout)                   ::                           &
        ierr

   ! zd is the normalized downdraft mass flux
   ! cdd is the downdraft detrainment function

     real,    dimension (ktf)                                          &
        ,intent (out  )                   ::                           &
        zd
     real,    dimension (ktf)                                          &
        ,intent (inout)                   ::                           &
        cdd, etd, facd1, facd2
!
!  local variables in this routine
!
     integer                              ::                           &
        ki
     real                                 ::                           &
        fact(ktf), dzsum(ktf)

     zd  = 0.
     cdd = 0.
     etd = 0.

     facd1 = 1.0
     facd2 = 0.0

     IF(ierr.ne.0)return

     ! special for shaved cells:
     zd(jmin)=arw(jmin) / arw0

     if (entr > 1.e-20) then
        do ki = 2, jmin-1
           etd(ki) = entr * dz_cup(ki)
           cdd(ki) = entr * dz_cup(ki)
        enddo
     endif
!
!--- integrate downward, specify detrainment(cdd)!
!
     dzsum(1) = dz_cup(1)
     do ki = 2, jmin
        dzsum(ki) = dzsum(ki-1) + dz_cup(ki)
     enddo

     do ki = 2, min(jmin-1, max(kfull,kdet))
        cdd(ki) = cdd(ki) + (1.- dzsum(ki-1) / dzsum(ki))
     enddo

     do ki = 2, jmin-1
        fact(ki) = 1. + etd(ki) - cdd(ki)
     enddo

     do ki = jmin-1, 2, -1
        zd(ki) = zd(ki+1) * fact(ki)
     enddo

     ! special with shaved cells
!!     do ki = 2, kfull
!!        zd(ki) = zd(ki) * min(arw(ki),arw0) / arw0
!!     enddo
!!
!!     do ki = 2, jmin-1
!!        cdd(ki) = entr + (zd(ki+1) - zd(ki)) / (zd(ki+1) * dz_cup(ki))
!!     enddo

     ! special at bottom
     zd (1) = 0.0
     cdd(1) = 1.0

     if (entr > 1.e-20) then
        do ki = 2, jmin-1
           facd2(ki) = etd(ki) / (1.0 + etd(ki) - 0.5 * cdd(ki))
           facd1(ki) = 1.0 - facd2(ki)
        enddo
     endif

   END SUBROUTINE cup_dd_nms


   SUBROUTINE cup_dellas_3d(ierr,rhodz_i,hcd,zd,cdd, &
              etd,he,della,zu,en,                    &
              cd,hc,ktop,k22,jmin,he_cup,ktf         )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (inout  )                   ::                           &
        della
     real,    dimension (ktf)                              &
        ,intent (in  )                   ::                            &
        rhodz_i,hcd,zd,cdd,etd,he,hc,cd,zu,he_cup,en
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop,k22,jmin
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     integer :: k
     real :: entdo,subin,detdo,entup,                              &
             detup,subdn

     della = 0.

     IF (ierr.ne.0) return

     della(1) = zd(2) * (.5*(HCD(1)+HCD(2)) - he_cup(2)) * rhodz_i(1)
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
     DO k = 2, ktop
        subin = zu(k+1) - zd(k+1)
        subdn = zu(k)   - zd(k)

        entup = en(k+1) * zu(k)
        detup = cd(k+1) * zu(k)

        detdo = cdd(k) * zd(k+1)
        entdo = etd(k) * zd(k+1)

        della(k) = ( detup * .5*(HC(K+1)+HC(K)) &
                   + detdo * .5*(HCD(K+1)+HCD(K)) &
                   - (entup + entdo) * he(k) &
                   + subin * he_cup(k+1) &
                   - subdn * he_cup(k) &
                   ) * rhodz_i(k)
     ENDDO

     ! Special at jmin
     della(jmin) = della(jmin) - zd(jmin) * hcd(jmin+1) * rhodz_i(jmin)

     ! Special at k22
     della(k22-1) = della(k22-1) - zu(k22) * hc(k22) * rhodz_i(k22-1)

   END SUBROUTINE cup_dellas_3d


   SUBROUTINE cup_dellas_nodd(ierr,rhodz_i,he,della,zu, &
              cd,hc,ktop,k22,en,he_cup,ktf)

     IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (inout  )                   ::                           &
        della
     real,    dimension (ktf)                              &
        ,intent (in  )                   ::                            &
        rhodz_i,he,hc,cd,en,zu,he_cup
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop,k22
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     integer :: k
     real :: subin,entup,detup,subdn

     della = 0.
     IF (ierr.ne.0) return
!
! no downdrafts for shallow convection
!
     DO k = k22, ktop
        subin = zu(k+1)
        subdn = zu(k)

        entup = en(k+1) * zu(k)
        detup = cd(k+1) * zu(k)
!
!--- CHANGE DUE TO SUBSIDENCE AND ENTRAINMENT
!
        della(k) = ( detup * .5*(hc(k+1)+hc(k)) &
                   - entup * he(k) &
                   + subin * he_cup(k+1) &
                   - subdn * he_cup(k) &
                   ) * rhodz_i(k)
     ENDDO

   END SUBROUTINE cup_dellas_nodd


   SUBROUTINE cup_env_heights(z,p,z1,psur,p_cup,z_cup,dz,dz_cup,ktf)

     IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation

     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        z,p
     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        p_cup,z_cup,dz_cup,dz
     real                                      &
        ,intent (in   )                   ::                           &
        psur,z1
!
!  local variables in this routine
!
     integer                              ::                           &
       k

     do k=2,ktf
        z_cup(k)=.5*(z(k-1)+z(k))
        p_cup(k)=.5*(p(k-1)+p(k))
     enddo

     z_cup(1)=z1
     p_cup(1)=psur

     do k = 1, ktf-1
        dz_cup(k) = z_cup(k+1) - z_cup(k)
     enddo

     do k = 2, ktf
        dz(k) = z(k) - z(k-1)
     enddo
     dz(1) = 2. * (z(1) - z_cup(1))

   END SUBROUTINE cup_env_heights


   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,ierr,ktf)

     use consts_coms, only: t00, eps_vap
     use therm_lib,   only: eslf, esif

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,ierr
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  !
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
        p,t,z
     real,    dimension (ktf)                                          &
        ,intent (inout)                   ::                           &
        he,hes,q,qes
!
!  local variables in this routine
!
     integer                              ::                           &
       k
     real                                 ::                           &
       e, h1

     if(ierr.ne.0) return

     do k=1,ktf
        if ( t(k) > tcrit ) then
           e = eslf( t(k) - t00 )
        else
           e = esif( t(k) - t00 )
        endif

        qes(k) = max( eps_vap * e / (100.*p(k) - (1.0-eps_vap)*e), 1.e-8)
        q  (k) = min( q(k), qes(k) )
     enddo
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
     h1     = g*z(1) + cp*t(1) + xl*q(1)
     he (1) = 0.0
     hes(1) = xl * (qes(1) - q(1))

     do k = 2, ktf
        he (k) = g*z(k) + cp*t(k) + xl*q(k) - h1
        hes(k) = he(k) + xl*(qes(k) - q(k))
     enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env2(qes,he,hes,t,q,p,ks,ke,ierr,ktf)

     use consts_coms, only: t00, eps_vap
     use therm_lib,   only: eslf, esif

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        p,t,he
     real,    dimension (ktf)                              &
        ,intent (inout  )                 ::                           &
        hes,q,qes
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        ks,ke
!
!  local variables in this routine
!
     integer                              ::                           &
       k

     real :: e

      if(ierr.ne.0) return

      do k = ks, ke
        if ( t(k) > tcrit ) then
           e = eslf( t(k) - t00 )
        else
           e = esif( t(k) - t00 )
        endif

        qes(k) = max( eps_vap * e / (100.*p(k) - (1.0-eps_vap)*e), 1.e-8)
        q  (k) = min( q(k), qes(k) )

        hes(k) = he(k) + xl*(qes(k) - q(k))
     enddo

   END SUBROUTINE cup_env2


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,qes_cup,q_cup, &
              he_cup,hes_cup,t_cup,ierr,ks,ke,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  !
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,t
     real,    dimension (ktf)                                          &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,t_cup
     integer                                                           &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        ks,ke
!
!  local variables in this routine
!
     integer                              ::                           &
       k,kbot

     if(ierr.ne.0)return

     kbot = max(2,ks)
     do k = kbot, ke
        qes_cup(k) = 0.5 * (qes(k-1) + qes(k))
        q_cup  (k) = 0.5 * (q  (k-1) + q  (k))
        hes_cup(k) = 0.5 * (hes(k-1) + hes(k))
        he_cup (k) = 0.5 * (he (k-1) + he (k))
        t_cup  (k) = 0.5 * (t  (k-1) + t  (k))
     enddo

     if (ks == 1) then
        qes_cup(1) = qes(1)
        q_cup  (1) = q  (1)
        hes_cup(1) = hes(1)
        he_cup (1) = he (1)
        t_cup  (1) = t  (1)
     endif

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_env_clev_uv(us,vs,ue_cup,ve_cup,ierr,ks,ke,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  !
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
        us,vs
     real,    dimension (ktf)                                          &
        ,intent (out  )                   ::                           &
        ue_cup,ve_cup
     integer                                                           &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        ks,ke
!
!  local variables in this routine
!
     integer                              ::                           &
       k,kbot

     if(ierr.ne.0)return

     kbot = max(2,ks)
     do k = kbot, ke
        ue_cup(k) = 0.5 * (us(k-1) + us(k))
        ve_cup(k) = 0.5 * (vs(k-1) + vs(k))
     enddo

     if (ks == 1) then
        ue_cup(1) = us(1)
        ve_cup(1) = vs(1)
     endif

   END SUBROUTINE cup_env_clev_uv


   SUBROUTINE cup_forcing_ens_3d(aa1,xaa0,ierr,xmb)

   IMPLICIT NONE

  ! aa1     = cloud work function
  ! xaa0    = cloud work function with cloud effects

     real                                                              &
        ,intent (out  )                   ::                           &
        xmb
     real                                                              &
        ,intent (in   )                   ::                           &
        aa1
     real                                                              &
        ,intent (in   )                   ::                           &
        xaa0
     integer                                                           &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     real                                 ::                           &
       xk

     real, parameter :: tau_kf = 60.*35.

     IF (ierr.ne.0) return
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
     if (aa1 <= xaa0) then

        ierr = 13
        xmb = 0.0

     else

        xk  = max(aa1-xaa0, 0.1)
        xmb = AA1 / (tau_kf*xk)

     endif

   END SUBROUTINE cup_forcing_ens_3d


   SUBROUTINE cup_kbcon(k22,kbcon,he_cup,hes_cup,ierr, &
              kbmax,ktop,p_cup,cap_max,ktf)

   IMPLICIT NONE

   ! only local wrf dimensions are need as of now in this routine
     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real                                      &
        ,intent (inout)                   ::                           &
        cap_max
     integer                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr,ktop
     integer                                                           &
        ,intent (in   )                   ::                           &
        kbmax
!
!  local variables in this routine
!
     integer                              ::                           &
        k
     real                                 ::                           &
        hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
     kbcon=1
     IF (ierr.ne.0) return

     k22loop: do while (k22 < kbmax)
        hetest = he_cup(k22)

        do k = k22+1, kbmax+2
           if (HETEST.GT.HES_cup(k)) exit
        enddo

        if (k .gt. kbmax+2) exit k22loop

        if (p_cup(k22) - p_cup(k) <= cap_max) then
           kbcon = k
           exit k22loop
        endif

        k22 = k22 + 1
     enddo k22loop

     ! if we are here and kbcon was not set, then no cloud base was found
     if (kbcon == 1) ierr=3

     ! estimate convective cloud top without entrainment
     if (ierr == 0) then
        do k = kbcon+1, ktf-1
           if (hetest < hes_cup(k)) exit
        enddo
        ktop = k-1
     endif

     if (ktop == kbcon) ierr=55
     if (ktop == ktf-1) ierr=5

   END SUBROUTINE cup_kbcon




   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,ktf)

   IMPLICIT NONE

  ! only local wrf dimensions are need as of now in this routine
     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
         array
     integer                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer                                      &
        ,intent (out  )                   ::                           &
         maxx

     maxx=ks
     if(ierr.eq.0 .and. ke.gt.ks)then
        maxx = ks + maxloc( array(ks:ke), 1 ) - 1
     endif

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
         array
     integer                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer                                      &
        ,intent (out  )                   ::                           &
         kt

     KT=KS
     if(ierr.eq.0 .and. kend.gt.ks) then
        kt = ks + minloc( array(ks:kend), 1 ) - 1
     endif

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_output_ens_3d(ierr,dellat,dellaq,dellaqc,        &
              outtem,outq,outqc,pwav,pwev,pwa_out,pev_out,pre,xmb, &
              ktop,pw,pwd,dellau,dellav,outu,outv,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,pwa_out,pev_out,outu,outv
     real                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real                                      &
        ,intent (in  )                  ::                           &
        pwav,pwev
     real,    dimension (ktf)                     &
        ,intent (in   )                   ::                           &
        dellat,dellaqc,dellaq,pw,pwd,dellau,dellav
     integer                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     integer                              ::                           &
        k

     outtem=0.
     outq=0.
     outqc=0.
     pwa_out=0.
     pev_out=0.
     outu=0.
     outv=0.
     pre=0.

     if(ierr.ne.0)return

     !mjo limit mass flux
     xmb = min(xmb, xmbmax)

     if(xmb.le.0.)then
        ierr=13
        xmb=0.
        return
     endif
!
!-- now do feedback
!
     do k = 1, ktop
        OUTTEM (K) = xmb * dellat(k)
        OUTQ   (K) = xmb * dellaq(k)
        OUTQC  (K) = xmb * dellaqc(k)
        PWA_OUT(K) = xmb * pw (k)
        PEV_OUT(K) = xmb * pwd(k)
        outu   (k) = xmb * dellau(k) * uvfac
        outv   (k) = xmb * dellav(k) * uvfac
     enddo

     pre = xmb * (pwav + pwev)

   END SUBROUTINE cup_output_ens_3d


   SUBROUTINE cup_up_aa0(aa0,dz,zu,hc,hes_cup,qes_cup,t_cup,kbcon,ktop,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf

  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        dz,zu,hc,hes_cup,qes_cup,t_cup
     integer                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!
     integer                                      &
        ,intent (inout)                   ::                           &
        ierr
     real                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!
     integer                              ::                           &
        k
     real                                 ::                           &
        dby

     aa0 = 0.0

     IF (ierr.ne.0) return

     do k = kbcon+1, ktop
        dby = max( hc(k-1) - hes_cup(k-1), 0.0 )

        aa0 = aa0 + gocp * zu(k) * dz(k) * dby * t_cup(k) &
                  / (t_cup(k)*t_cup(k) + gamfac * qes_cup(k))
     enddo

   END SUBROUTINE cup_up_aa0


   SUBROUTINE cup_up_he(k22,he_cup,hc,ktop,ierr,he,facu1,facu2,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf

  ! hc = cloud moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! cd= detrainment function
  ! dz_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,facu1,facu2
     integer                                      &
        ,intent (in   )                   ::                           &
        k22,ktop
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        hc
!
!  local variables in this routine
!
     integer                              ::                           &
        k
!
!--- moist static energy inside cloud
!
     hc = 0.
     if (ierr.ne.0) return

     hc(k22) = he_cup(k22)

     do k = k22+1, ktop+1
        hc(k) = hc(k-1) * facu1(k) + he(k-1) * facu2(k)
     enddo

   END SUBROUTINE cup_up_he



   SUBROUTINE cup_up_ktop(cd,en,hc,zu,kbcon,ktop,ierr, &
              hes_cup,facu1,facu2,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf

  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! cd= detrainment function
  ! entr = entrainment rate
  !
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
        hes_cup
     integer                                                           &
        ,intent (in   )                   ::                           &
        kbcon
!
! input and output
!
     integer                                                           &
        ,intent (inout)                   ::                           &
        ierr,ktop

     real,    dimension (ktf)                                          &
        ,intent (inout)                   ::                           &
        hc,cd,zu,en,facu1,facu2
!
!  local variables in this routine
!
     integer                              ::                           &
        k

     if(ierr.ne.0)return
!
!--- Search if real cloud top is less than ktop
!
     do k = kbcon+1, ktop
        if (hc(k) < hes_cup(k)) exit
     enddo
     ktop = k - 1

     if (ktop == kbcon) ierr=55
!
!--- Redefine values at new ktop
!
     if (ierr.eq.0) then
        hc   (ktop+1) = hc(ktop)
        cd   (ktop+1) = 1.0
        zu   (ktop+1) = 0.0
        en   (ktop+1) = 0.0
        facu1(ktop+1) = 1.0
        facu2(ktop+1) = 0.0
     endif

   END SUBROUTINE cup_up_ktop



   SUBROUTINE cup_up_moisture(ierr,dz_cup,qc,qrc,pw,pwav,kbcon,ktop,cd, &
              hc,hes_cup,q,qes_cup,zu,facu1,facu2,k22,qe_cup,t_cup,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! entr_rate = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        q,zu,qe_cup,qes_cup,dz_cup,cd,facu1,facu2,hc,hes_cup,t_cup
  ! entr= entrainment rate
     integer                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!
   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw
     real                                      &
        ,intent (out  )                   ::                           &
        pwav
!
!  local variables in this routine
!
     integer                              ::                           &
        k
     real                                 ::                           &
        clw_all,pwat,gamma,dby
     real,    dimension (ktf)             ::                           &
        qsat,facu3,frain

     ! precipitation efficiency
     real, parameter :: c0 = 0.002

     pwav=0.
     pw = 0.
     qc = 0.
     qrc= 0.

     if(ierr.ne.0)return

     do k = kbcon+1, ktop+1
        DBY      = MAX(HC(K) - HES_CUP(K), 0.)
        GAMMA    = GAMFAC * QES_CUP(K) / (T_CUP(K)*T_CUP(K))

        ! Saturation in cloud, this is what is allowed to be in it
        QSAT (K) = QES_CUP(K) + DBY * GAMMA * XLM / (1. + GAMMA)

        ! Fraction of available water rained-out
        FRAIN(K) = 1.0 / (1.+C0*DZ_CUP(K-1)*0.5*(ZU(K)+ZU(K-1)))

        ! Mass flux to convert mixing ratio to updraft fluxes
        FACU3(K) = ZU(K) + 0.5 * ZU(K-1) * CD(K)
     enddo

     qc(k22) = qe_cup(k22)

     ! Plume equation without condesation below cloud base
     do k = k22+1, kbcon
        qc(k) = qc(k-1) * facu1(k) + q(k-1) * facu2(k)
     enddo

     ! Within cloud
     do k = kbcon+1, ktop+1

        ! Steady state plume equation, for what could be in cloud
        ! without condensation
        qc(k) = qc(k-1) * facu1(k) + q(k-1) * facu2(k)

        ! Liquid water content in cloud
        CLW_ALL = MAX( QC(K)-QSAT(K), 0.0 )

        ! Liquid water content in cloud after rainout
        QRC(K) = CLW_ALL * FRAIN(K)

        ! Condensation
        PWAT    = CLW_ALL - QRC(K)
        QC(K)   = QC(K) - PWAT
        PW(K-1) = PWAT * FACU3(K)

        ! Integrated normalized condensate
        PWAV = PWAV + PW(K-1)

     enddo

   END SUBROUTINE cup_up_moisture


   SUBROUTINE cup_up_moisture_shall(ierr,qc,qrc,kbcon,ktop,hc,  &
              hes_cup,t_cup,q,qe_cup,qes_cup,facu1,facu2,k22,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! gamma_cup = gamma on model cloud levels
  ! entr_rate = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        hc,hes_cup,t_cup,q,qe_cup,qes_cup,facu1,facu2
  ! entr= entrainment rate
     integer                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!
   ! ierr error value, maybe modified in this routine

     integer                                      &
        ,intent (inout)                   ::                           &
        ierr

   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout

     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        qc,qrc
!
!  local variables in this routine
!
     integer                              ::                           &
        k
     real                                 ::                           &
        dby, gamma, qsat

     qc = 0.
     qrc= 0.

     if(ierr.ne.0)return

     qc(k22) = qe_cup(k22)

     ! Steady state plume equation
     DO K = K22+1, KTOP+1
        QC(K) = QC(K-1) * FACU1(K) + Q(K-1) * FACU2(K)
     ENDDO

     ! Liquid water content in cloud
     DO K = KBCON+1, KTOP+1
        DBY    = MAX(HC(K) - HES_CUP(K), 0.)
        GAMMA  = GAMFAC * QES_CUP(K) / (T_CUP(K)*T_CUP(K))
        QSAT   = QES_CUP(K) + DBY * GAMMA * XLM / (1. + GAMMA)

        QRC(K) = MAX(QC(K) - QSAT, 0.0)
     ENDDO

   END SUBROUTINE cup_up_moisture_shall



   SUBROUTINE cup_up_uv(k22,ue_cup,ve_cup,uc,vc,ktop,ierr,us,vs,facu1,facu2,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf

  ! hc = cloud moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! cd= detrainment function
  ! dz_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        facu1,facu2
     real,    dimension (ktf)                              &
        ,intent (in   )                   ::                           &
        ue_cup,ve_cup,us,vs
     integer                                               &
        ,intent (in   )                   ::                           &
        k22,ktop
     integer                                                &
        ,intent (in   )                   ::                           &
        ierr
!
! input and output
!
     real,    dimension (ktf)                              &
        ,intent (out  )                   ::                           &
        uc,vc
!
!  local variables in this routine
!
     integer                              ::                           &
        k
!
!--- moist static energy inside cloud
!
     uc = 0.0
     vc = 0.0
     if (ierr.ne.0) return

     uc(k22) = ue_cup(k22)
     vc(k22) = ve_cup(k22)

     do k = k22+1, ktop+1
        uc(k) = uc(k-1) * facu1(k) + us(k-1) * facu2(k)
        vc(k) = vc(k-1) * facu1(k) + vs(k-1) * facu2(k)
     enddo

   END SUBROUTINE cup_up_uv


   SUBROUTINE cup_up_nms(zu,cd,en,ktop,ierr,k22,facu1,facu2,arw,arw0,ktf)

   IMPLICIT NONE
!
!  on input
!
   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,k22,ktop
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
         arw
     real                                                              &
        ,intent (in   )                   ::                           &
         arw0
!
! input and output
!
     integer                                                           &
        ,intent (inout)                   ::                           &
         ierr
     real,    dimension (ktf)                                          &
        ,intent (inout  )                 ::                           &
         zu,cd,en,facu1,facu2
!
!  local variables in this routine
!
     integer                              ::                           &
         k

     IF(ierr.ne.0)return
!
! do normalized mass budget
!
     en(1:k22) = 0.0
     cd(1:k22) = 0.0
     zu(1:k22-1) = 0.0

     ! special for shaved cells
     zu(k22) = arw(k22) / arw0

     do k = k22+1, ktop
        zu(k) = zu(k-1) * (1.0 + en(k) - cd(k))
     enddo

     en   (ktop+1) = 0.0
     cd   (ktop+1) = 1.0
     zu   (ktop+1) = 0.0
     facu1(ktop+1) = 1.0
     facu2(ktop+1) = 0.0

   END SUBROUTINE cup_up_nms



   SUBROUTINE cup_up_en_cd(en,cd,facu1,facu2,dz_cup,arw,entr_rate,detr_rate, &
                           kfull,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!
     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kfull,ierr
     real,    dimension (ktf)                                          &
        ,intent (in   )                   ::                           &
         arw,dz_cup
     real                                                              &
        ,intent (in   )                   ::                           &
         entr_rate,detr_rate
!
! input and output
!
     real,    dimension (ktf)                                          &
        ,intent (out  )                   ::                           &
         cd,en,facu1,facu2
!
!  local variables in this routine
!
     integer                              ::                           &
         k

     if (ierr.ne.0) return

     ! special for shaved cells
     do k = 2, kfull
        en(k) = max(arw(k) / arw(k-1) - 1.0, entr_rate*dz_cup(k-1))
     enddo

     do k = kfull+1, ktf
        en(k) = entr_rate*dz_cup(k-1)
     enddo

     do k = 2, ktf
        cd   (k) = detr_rate * dz_cup(k-1)
        facu2(k) = en(k) / (1.0 + en(k) - 0.5 * cd(k))
        facu1(k) = 1.0 - facu2(k)
     enddo

   END SUBROUTINE cup_up_en_cd

END MODULE module_cu_g3
