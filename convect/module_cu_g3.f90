MODULE module_cu_g3

  use mem_para,  only: myrank
  use misc_coms, only: io6
  implicit none

  private :: myrank, io6

  integer, parameter, private :: iens    = 1
  integer, parameter, private :: maxens  = 3
  integer, parameter, private :: maxens2 = 3
  integer, parameter, private :: maxens3 = 16
  integer, parameter, private :: ensdim  = maxens*maxens2*maxens3

  integer, parameter, private :: its = 1, ite = 1, itf = 1
  integer, parameter, private :: jts = 1, jte = 1, jtf = 1

  real,    parameter, private :: xmbmax = 1.0

CONTAINS

  SUBROUTINE grell_driver(iw, dtlong)

     use mem_turb,    only: kpblh, frac_land, fthpbl, fqtpbl
     use consts_coms, only: cp, alvl, grav, p00i, rocp, rvap, eradi, &
                            gravi, alvlocp, r8
     use mem_radiate, only: rshort, fthrd_sw, fthrd_lw
     use mem_grid,    only: mza, lpw, arw0, zm, zt, xew, yew, zew, &
                            dzt, arw, lpv, arv, volt, volti
     use mem_ijtabs,  only: itab_w
     use mem_basic,   only: wmc, vmc, theta, tair, press, rho, sh_v, &
                            vxe, vye, vze
     use mem_cuparm,  only: thsrc, rtsrc, conprr, kcutop, kcubot, cbmf, &
                            qwcon, iactcu, kddtop, cddf, rdsrc

     implicit none

     integer, intent(in)  :: iw
     real,    intent(in)  :: dtlong

     real    :: outqc  (1,mza) ! output deep conv cloud water tendency
     integer :: j              ! dummy horizontal index
     real    :: aaeq   (1)     ! turn off convection if negative
     real    :: t      (1,mza) ! input temperature
     real    :: q      (1,mza) ! input mixing ratio
     real    :: z1     (1)     ! surface elevation
     real    :: sub_mas(1,mza)
     real    :: tn     (1,mza) ! input forced temperature
     real    :: qo     (1,mza) ! input forced mixing ratio
     real    :: po     (1,mza) ! input forced pressure
     real    :: pre    (1)     ! output precipitation rate
     real    :: p      (1,mza) ! input pressure
     real    :: outt   (1,mza) ! output deep conv temp tendency
     real    :: outq   (1,mza) ! output deep conv water vapor tendency
     real    :: dtime          ! time step
     integer :: ktau           ! random # seed, not used
     real    :: tkmax  (1)     ! not used
     real    :: psur   (1)     ! input surface pressure
     real    :: us     (1,mza) ! input zonal wind
     real    :: vs     (1,mza) ! input meridional wind
     real    :: tcrit          ! temperature for cloud water phase
     real    :: tx     (1,mza,8) ! input temperature ensembles 
     real    :: qx     (1,mza,8) ! input water vapor ensembles 
     real    :: tshall (1,mza)   ! input PBL forced temperature for shallow convection
     real    :: qshall (1,mza)   ! input PBL forced water vapor for shallow convection
     integer :: kpbl   (1)       ! layer corresponding to PBL height
     real    :: dhdt   (1,mza)   ! moist static energy tendency
     real    :: outts  (1,mza)   ! shallow conv temp tendency
     real    :: outqs  (1,mza)   ! shallow conv water vapor tendency
     real    :: tscl_kf          ! shallow convection timescale
     integer :: k23    (1)       ! shallow convection updraft source level
     integer :: kbcon3 (1)       ! shallow convection LCL
     integer :: ktop3  (1)       ! shallow convection cloud top
     real    :: xmb3   (1)       ! shallow convection mass flux
     real    :: mconv  (1,    8) ! column moisture convergence
     real    :: omeg   (1,mza,8) ! vertical velocity in pressure coordinates
     real    :: qce    (1,mza)   ! Moisture tendency due to precipitation only
                                 ! (kg water removed from grid cell / m^2 / sec )
     real    :: APR_GR (1,1)     ! precip (mm/hr) for different closures/caps
     real    :: APR_W  (1,1)
     real    :: APR_MC (1,1)
     real    :: APR_ST (1,1)
     real    :: APR_AS (1,1)
     real    :: APR_CAPMA(1,1)
     real    :: APR_CAPME(1,1)
     real    :: APR_CAPMI(1,1)

     integer :: k22    (1)       ! deep convectin updraft source level
     integer :: kbcon  (1)       ! deep convection LCL
     integer :: ktop   (1)       ! deep convection cloud top
     integer :: jmin   (1)       ! downdraft originating level
     real    :: xmb    (1)       ! deep convection mass flux
     real    :: cupclw (1,mza)   ! deep convection cloud water
     real    :: xf_ens (1,1,ensdim) ! mass flux ensembles
     real    :: pr_ens (1,1,ensdim) ! precipitation ensembles
     real    :: xland  (1)       ! 1 = land, 2 = water
     real    :: gsw    (1)       ! surface incoming radiation
     real    :: edt_out(1,1)
     real    :: subt   (1,mza)   ! temp tendency due to subsidence
     real    :: subq   (1,mza)   ! water vapor tendency due to subsidence
     real    :: cupclws(1,mza)   ! shallow convection cloud water
     real    :: xl               ! latent heat of vaporization
     real    :: rv               ! gas constant for water vapor
     real    :: cpd              ! heat capacity of air
     real    :: g                ! gravitational acceleration
     integer :: ichoice          ! flag if only want one closure (usually set to zero)
     integer :: ipr,jpr          ! unused
     integer :: ens4             ! number of horizontal ensembles
     integer :: high_resolution  ! flag = 1 activates  model changes for high res runs
     integer :: ishallow_g3      ! flag = 1 activates shallow convection

     integer :: kts              ! starting k level
     integer :: ktf              ! ending   k level
     integer :: kte              ! max size of k dimensions

     ! local variables
     integer :: k, ka, kc, npoly, n, iwn, jv, iv
     real    :: raxis, raxisi, omega_ave, dirv, flx
     real    :: exner(mza)
     real    :: vflux(mza), vflux_the(mza), vflux_vap(mza)
     real    :: hflux, hflux_the, hflux_vap
     real    :: fqvadv, fthadv
     real(r8):: qsum, qav, tsum, tav
     real    :: dtemp, dsh_v, fens4m
     real    :: v, vp, vm

     ka    = lpw(iw)
     npoly = itab_w(iw)%npoly
     kpbl  = kpblh(iw) - ka + 1

     ktf  = mza + 1 - ka  ! number of ATM levels to process
     kts  = 1
     kte  = mza  ! array dimension
     ens4 = npoly + 1

     fens4m = 1.0 / real(ens4)

     ! Go no higher than 50mb for convective calculations to prevent any
     ! problems when esat gets near ambient pressure
     do k = ka, mza-1
        if (press(k,iw) < 50.e2) then
           ktf = k - ka + 1
           exit
        endif
     enddo

     ! A. Betts for shallow convection: suggestion for the KF timescale < DELTAX / 25 m/s
     tscl_kf = sqrt(arw0(iw)) / 25.

     dtime = dtlong
     gsw   = rshort(iw)
     psur  = 0.01 * (press(ka,iw) + (zt(ka)-zm(ka-1))*rho(ka,iw)*grav)
     z1    = zm(ka-1)
     
     ! one if by land, two if by sea
     if (allocated(frac_land)) then
        xland = 2.0 - frac_land(iw)
     else
        xland = 1.0
     endif

     raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
     raxisi = 1.0 / max(raxis, 1.e-12)

     ! Vertical advective theta and water vapor fluxes (W levels)

     do k = ka, mza-1
        v  = arw(k,iw) * wmc(k,iw)
        vm = max(v, 0.0)
        vp = min(v, 0.0)

        ! upwinded fluxes
        vflux    (k) = v
        vflux_vap(k) = vm * sh_v (k,iw) + vp * sh_v (k+1,iw)
        vflux_the(k) = vm * theta(k,iw) + vp * theta(k+1,iw)
     enddo

     vflux    (ka-1)  = 0.
     vflux    (mza) = 0.
     vflux_the(ka-1)  = 0.
     vflux_the(mza) = 0.
     vflux_vap(ka-1)  = 0.
     vflux_vap(mza) = 0.

     ! Loop over T levels

     !dir$ simd
     do kc = 1, ktf
        k  = kc + ka - 1

        exner(k) = tair(k,iw) / theta(k,iw)

        ! Current temp, water vapor, and pressure

        t(1,kc)  = max( tair(k,iw), 200.0 )
        q(1,kc)  = max( sh_v(k,iw), 1.e-8 )
        p(1,kc)  = 0.01 * press(k,iw)

        ! Horizontal advective mass and water vapor fluxes

        hflux     = 0.
        hflux_the = 0.
        hflux_vap = 0.
      
        !dir$ loop count max=7, min=5, avg=6
        do jv = 1, npoly
           iv  = itab_w(iw)%iv(jv)

           if (k >= lpv(iv)) then
              iwn  = itab_w(iw)%iw(jv)
              dirv = itab_w(iw)%dirv(jv)

              flx   = dirv * vmc(k,iv) * arv(k,iv)
              hflux = hflux + flx

              ! upwinded
              if (flx >= 0.0) then
                 hflux_vap = hflux_vap + flx * sh_v (k,iwn)
                 hflux_the = hflux_the + flx * theta(k,iwn)
              else
                 hflux_vap = hflux_vap + flx * sh_v (k,iw)
                 hflux_the = hflux_the + flx * theta(k,iw)
              endif
           endif
        enddo

        fthadv = ((vflux_the(k-1) - vflux_the(k) + hflux_the) &
               - (vflux(k-1) - vflux(k) + hflux) * theta(k,iw)) &
               / (volt(k,iw) * rho(k,iw))

        fqvadv = ((vflux_vap(k-1) - vflux_vap(k) + hflux_vap) &
               - (vflux(k-1) - vflux(k) + hflux) * sh_v(k,iw)) &
               / (volt(k,iw) * rho(k,iw))

        ! "forced" temp, water vapor, and pressure
        dtemp = (fthrd_lw(k,iw) + fthrd_sw(k,iw) + &
                 fthpbl(k,iw) + fthadv) * dtlong * exner(k)
        dsh_v = (fqtpbl(k,iw) + fqvadv) * dtlong

        tn(1,kc) = tair(k,iw) + dtemp
        qo(1,kc) = sh_v(k,iw) + dsh_v
        tn(1,kc) = max( tn(1,kc), 200.0 )
        qo(1,kc) = max( qo(1,kc), 1.e-8 )
        po(1,kc) = p(1,kc)

        ! U and V winds
        if (raxis > 1.e3) then
           us(1,kc) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi

           vs(1,kc) = vze(k,iw) * raxis * eradi  &
                - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * eradi * raxisi
        else
           us(1,kc) = vxe(k,iw)
           vs(1,kc) = vye(k,iw)
        endif

        ! Horizontal ensemble members (1st ensemble is current column)
        tx  (1,kc,1) = tn(1,kc)
        qx  (1,kc,1) = qo(1,kc)
        omeg(1,kc,1) = -grav * wmc(k,iw)

        ! Shallow convection uses current T and Q plus PBL tendencies
        tshall(1,kc) = t(1,kc) + fthpbl(k,iw) * dtlong * exner(k)
        qshall(1,kc) = q(1,kc) + fqtpbl(k,iw) * dtlong
        dhdt  (1,kc) = cp * exner(k) * fthpbl(k,iw) + alvl * fqtpbl(k,iw)
        
        ! Set area ensembles of omega, T, and Q
        !dir$ loop count max=7, min=5, avg=6
        do n = 1, npoly

           if (k >= lpw(itab_w(iw)%iw(n))) then
              iwn = itab_w(iw)%iw(n)
           else
              iwn = iw
           endif
        
           tx  (1,kc,n+1) = max( tair(k,iwn) + dtemp, 200.0 )
           qx  (1,kc,n+1) = max( sh_v(k,iwn) + dsh_v, 1.e-8 )
           omeg(1,kc,n+1) = -grav * wmc(k,iwn)
        enddo
     enddo

     ! Compute area ensemble moisture convergence based on average omega

     mconv(1,:) = 0.0

     do kc = 1, ktf
        k  = kc + ka - 1
        omega_ave = sum(omeg(1,kc,1:ens4)) * fens4m
        mconv(1,1) = mconv(1,1) + omega_ave * (sh_v(k+1,iw) - sh_v(k,iw)) * gravi

        do n = 1, npoly
           iwn = itab_w(iw)%iw(n)

           if (k >= lpw(iwn)) then
              mconv(1,n+1) = mconv(1,n+1) + omega_ave * (sh_v(k+1,iwn) - sh_v(k,iwn)) * gravi
           endif
        enddo
     enddo

     mconv(1,1:ens4) = max(mconv(1,1:ens4), 0.00)

     ishallow_g3 = 1
     tcrit = 258.0
     aaeq = 1.
     j = 1
     ichoice = 0
     edt_out = 0.0
     high_resolution = 0
!    if (arw0(iw) < 1.e8) high_resolution = 1

     APR_GR    = 0.0
     APR_W     = 0.0
     APR_MC    = 0.0
     APR_ST    = 0.0
     APR_AS    = 0.0
     APR_CAPMA = 0.0
     APR_CAPME = 0.0
     APR_CAPMI = 0.0

     ipr = iw
     jpr = 0

     xl  = alvl
     rv  = rvap
     cpd = cp
     g   = grav

     outt (1,:) = 0.0
     outts(1,:) = 0.0
     outq (1,:) = 0.0
     outqs(1,:) = 0.0
     outqc(1,:) = 0.0
     subt (1,:) = 0.0
     subq (1,:) = 0.0
     qce  (1,:) = 0.0
     pre  (1)   = 0.0

     k23 = 0
     kbcon3 = 0
     ktop3 = 0
     xmb3 = 0.0
     ktau = 0
     tkmax = 0.0
     k22 = 0
     kbcon = 0
     ktop = 0
     xmb = 0.0
     cupclw = 0.0
     cupclws = 0.0
     xf_ens = 0.0
     pr_ens = 0.0
     sub_mas = 0.0
     jmin = 0

     call  CUP_enss_3d(iw,OUTQC,J,AAEQ,T,Q,Z1,sub_mas,                 &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,tkmax,PSUR,US,VS,    &
              TCRIT,tx,qx,                                             &
              tshall,qshall,kpbl,dhdt,outts,outqs,tscl_kf,             &
              k23,kbcon3,ktop3,xmb3,                                   &
              mconv,omeg,k22,xmb,jmin,                                 &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                       &
              APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,cupclw,         &
              xf_ens,pr_ens,xland,gsw,edt_out,subt,subq,               &
              xl,rv,cpd,g,ichoice,ipr,jpr,ens4,high_resolution,        &
              qce,ishallow_g3,cupclws,ktf,kts,kte                      )

     if (pre(1) > 1.e-16) then

        ! Deep convecton is active

        kcutop(iw) = ktop (1) + ka - 1
        kddtop(iw) = jmin(1) + ka - 1
        kcubot(iw) = kbcon(1) + ka - 1
        iactcu(iw) = 1
        cbmf  (iw) = xmb(1)
        cddf  (iw) = edt_out(1,1) * xmb(1)
        conprr(iw) = pre(1)

        do kc = 1, ktf
           k  = kc + ka - 1

           ! Total water tendency
           rtsrc(k,iw) = (outq(1,kc) + subq(1,kc) + 0.*outqc(1,kc)) * rho(k,iw)

           ! Any cloud condensate is evaporated since we do not feed back
           ! to resolved microphysics
           thsrc(k,iw) = (outt(1,kc) + subt(1,kc) - alvlocp * outqc(1,kc)) * rho(k,iw)

           ! Cloud water
           qwcon(k,iw) = cupclw(1,kc)

           ! Density change (Water removed)
           rdsrc(k,iw) = -qce(1,kc) * arw0(iw) * volti(k,iw)
        enddo

     else if (ishallow_g3 == 1 .and. kbcon3(1) > 0 .and. ktop3(1) >= kbcon3(1)) then

        ! Shallow convection is active

        kcutop(iw) = ktop3 (1) + ka - 1
        kcubot(iw) = kbcon3(1) + ka - 1
        iactcu(iw) = 1
        cbmf  (iw) = xmb3(1)

        do kc = 1, ktf
           k  = kc + ka - 1
           rtsrc(k,iw) = outqs(1,kc) * rho(k,iw)
           thsrc(k,iw) = outts(1,kc) * rho(k,iw)
           qwcon(k,iw) = cupclws(1,kc)
        enddo

     endif

   END SUBROUTINE grell_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE CUP_enss_3d(iw,OUTQC,J,AAEQ,T,Q,Z1,sub_mas,              &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,tkmax,PSUR,US,VS,    &
              TCRIT,tx,qx,                                             &
              tshall,qshall,kpbl,dhdt,outts,outqs,tscl_kf,             &
              k23,kbcon3,ktop3,xmb3,                                   &
              mconv,                                                   &
              omeg,k22,xmb,jmin,                                       &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                       &
              APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,cupclw,         &
              xf_ens,pr_ens,xland,gsw,edt_out,subt,subq,               &
              xl,rv,cp,g,ichoice,ipr,jpr,ens4,high_resolution,         &
              qce,ishallow_g3,cupclws,ktf,kts,kte                      )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        iw,ktf,ktau,                                              &
        kts,kte,ipr,jpr,ens4,high_resolution
     integer, intent (in   )              ::                           &
        j,ishallow_g3,ichoice     
  !
  ! 
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        xf_ens,pr_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout )                  ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI,edt_out
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                         &
        DHDT,OUTT,OUTQ,OUTQC,subt,subq,sub_mas,cupclw,outts,outqs,cupclws,qce
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb3,xmb
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        k22,kbcon,ktop,k23,kbcon3,ktop3,jmin
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        T,PO,P,US,VS,tn,tshall
     real,    dimension (its:ite,kts:kte,1:ens4)                       &
        ,intent (inout   )                   ::                           &
        omeg,tx,qx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO,qshall
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        Z1,PSUR,AAEQ,tkmax,xland
     real, dimension (its:ite,1:ens4)                                         &
        ,intent (in   )                   ::                           &
        mconv

       
       real                                                            &
        ,intent (in   )                   ::                           &
        dtime,tcrit,xl,cp,rv,g,tscl_kf

!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens,subt_ens,subq_ens
!
!
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
  ! mentr_rate = entrainment rate
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
  ! xf_ens = mass flux ensembles
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! mentr_rate = entrainment rate
     real,    dimension (its:ite,kts:kte) ::                           &
        he3,hes3,qes3,z3,zdo3,zu3_0,hc3_0,dby3_0,                      &
        qes3_cup,q3_cup,he3_cup,hes3_cup,z3_cup,gamma3_cup,t3_cup,     &
        xhe3,xhes3,xqes3,xz3,xt3,xq3,                                  &
        xqes3_cup,xq3_cup,xhe3_cup,xhes3_cup,xz3_cup,xgamma3_cup,      &
        xt3_cup,                                                       &
        xdby3,xqc3,xhc3,xqrc3,xzu3,                                    &
        dby3,qc3,pw3,hc3,qrc3,zu3,cd3,DELLAH3,DELLAQ3,                 &
        dsubt3,dsubq3,DELLAT3,DELLAQC3

     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,                                                  &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,          &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,scr1,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       aa3_0,aa3,hkb3,qkb3,pwav3,bu3,xaa3,xhkb3,                       &
       hkb3_0,edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,aad,XHKB,edt3,                                             &
       XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                    &
       PWEVO,BU,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,cap_max3
     real,    dimension (its:ite,1:ens4) ::                            &
        axx
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,kstabi,kstabm,K22x,jmin3,kdet3,                     &
       KBCONx,ierr,ierr2,ierr3,KBMAX,ierr5,ierr5_0 

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K
     real                                 ::                           &
      day,dz,mbdt,mbdt_s,entr_rate,radius,mentr_rate,mentrd_rate,      &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,            &
      dh,cap_maxs,trash,entr_rate3,mentr_rate3

     integer :: jmini
     logical :: keep_going
     real xff_shal(9),blqe,xkshal

      day=86400.
      do i=its,itf
        xmb3(i)=0.
        closure_n(i)=16.
        xland1(i)=1.
        if(xland(i).gt.1.5)xland1(i)=0.
!       cap_max_increment(i)=50.
        cap_max_increment(i)=25.
      enddo
!
!--- specify entrainmentrate and detrainmentrate
!
      if(iens.le.4)then
      radius=14000.-float(iens)*2000.
      else
      radius=10000.
      endif
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
      entr_rate =.2/radius
      entr_rate3=.2/300.
!
!--- entrainment of mass
!
      mentrd_rate=0.
      mentr_rate=entr_rate
      mentr_rate3=entr_rate3
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        cupclw(i,k)=0.
        cupclws(i,k)=0.
        cd(i,k)=0.01*entr_rate
        cd3(i,k)=entr_rate3
        cdd(i,k)=0.
        zdo3(i,k)=0.
        hcdo(i,k)=0.
        qrcdo(i,k)=0.
        dellaqc(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax=1.
      edtmin=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=750.
!
!--- maximum depth (mb) of capping inversion that convection is
!--- allowed to overcome (smaller cap_max means less convection)
!
!     cap_maxs=125.
!     cap_maxs=75.
      cap_maxs=50.

      DO i=its,itf
        kbmax(i)=1
        jmin3(i)=0
        kdet3(i)=0
        aa0(i)=0.
        aa3_0(i)=0.
        aa1(i)=0.
        aa3(i)=0.
        aad(i)=0.
        edt(i)=0.
        edt3(i)=0.
        kstabm(i)=ktf-1
        IERR(i)=0
        IERR2(i)=0
        IERR3(i)=0
        IERR5(i)=0
        IERR5_0(i)=0
      enddo

      do i=its,itf
          cap_max(i)=cap_maxs
          cap_max3(i)=25.
          if ((gsw(i,j).lt.1.0.and.xland1(i)>0.5).or.high_resolution.eq.1) then
             cap_max(i)=max(cap_max(i)-25.,25.)
          endif
      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1250.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(iw,z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           ktf,kts,kte)

      call cup_env(iw,zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      do i=its,itf
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!     if(ierr(i).eq.0)then
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
!     endif
      enddo
!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22,ierr, &
           ktf,kts,kte)
       DO 36 i=its,itf
         IF(ierr(I).eq.0)THEN
         IF(K22(I).GE.KBMAX(i))ierr(i)=2
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
           ierr,kbmax,po_cup,cap_max, &
           ktf,kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,  &
           ktf,kts,kte)
      do i=its,itf
      IF(ierr(I).eq.0)THEN
        if(kstabm(i)-1.gt.kstabi(i))then
           do k=kstabi(i),kstabm(i)-1
             cd(i,k)=cd(i,k-1)+.15*entr_rate
             if(cd(i,k).gt.1.0*entr_rate)cd(i,k)=1.0*entr_rate
           enddo
        ENDIF
      ENDIF
      ENDDO
!
!--- calculate incloud moist static energy
!
      call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbcon,ierr,dby,he,hes_cup,'deep', &
           ktf,kts,kte)
      call cup_up_he(k22,hkbo,zo_cup,cd,mentr_rate,heo_cup,hco, &
           kbcon,ierr,dbyo,heo,heso_cup,'deep', &
           ktf,kts,kte)
!
!--- DETERMINE CLOUD TOP - KTOP
!
      call cup_ktop(1,dbyo,kbcon,ktop,ierr, &
           ktf,kts,kte)
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,kte
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr, &
           ktf,kts,kte)
      DO 100 i=its,ite
      IF(ierr(I).eq.0)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
      jmini = jmin(i)
      keep_going = .TRUE.
      do while ( keep_going )
        keep_going = .FALSE.
        if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
        if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
        ki = jmini
        hcdo(i,ki)=heso_cup(i,ki)
        DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
        dh=0.
        do k=ki-1,1,-1
          hcdo(i,k)=heso_cup(i,jmini)
          DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
          dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
          if(dh.gt.0.)then
            jmini=jmini-1
            if ( jmini .gt. 3 ) then
              keep_going = .TRUE.
            else
              ierr(i) = 9
              exit
            endif
          endif
        enddo
      enddo
      jmin(i) = jmini 
      if ( jmini .le. 3 ) then
        ierr(i)=4
      endif
      ENDIF
100   continue
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
      IF(ierr(I).eq.0)THEN
      IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
            ierr(i)=6
      endif
      endif
      enddo

!
!c--- normalized updraft mass flux profile
!
      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           ktf,kts,kte)
      call cup_up_nms(zuo,zo_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           ktf,kts,kte)
!
!c--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      call cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr, &
           0,kdet,z1,                 &
           ktf,kts,kte)
      call cup_dd_nms(zdo,zo_cup,cdd,mentrd_rate,jmin,ierr, &
           1,kdet,z1,                 &
           ktf,kts,kte)
!
!--- downdraft moist static energy
!
      call cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate, &
           jmin,ierr,he,dbyd,he_cup,  &
           ktf,kts,kte)
      call cup_dd_he(heso_cup,zdo,hcdo,zo_cup,cdd,mentrd_rate, &
           jmin,ierr,heo,dbydo,he_cup,&
           ktf,kts,kte)
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_3d(zd,hcd,hes_cup,qcd,qes_cup, &
           pwd,q_cup,z_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           pwev,bu,qrcd,q,he,t_cup,2,xl,high_resolution, &
           ktf,kts,kte)
      call cup_dd_moisture_3d(zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,cdd,mentrd_rate,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,xl,high_resolution, &
           ktf,kts,kte)
!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture('deep',ierr,z_cup,qc,qrc,pw,pwav, &
           kbcon,ktop,cd,dby,mentr_rate,clw_all,      &
           q,GAMMA_cup,zu,qes_cup,k22,q_cup,xl, &
           ktf,kts,kte)
      do k=kts,ktf
      do i=its,itf
         cupclw(i,k)=qrc(i,k)
      enddo
      enddo
      call cup_up_moisture('deep',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           kbcon,ktop,cd,dbyo,mentr_rate,clw_all, &
           qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,&
           ktf,kts,kte)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
            if (aa1(i) < 1.e-12) then
               ierr(i)=17
           endif
         endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    NEXT section for shallow convection
!
      if (ishallow_g3.eq.1 .and. ierr(1) .ne. 0) then

      call cup_env(iw,z3,qes3,he3,hes3,tshall,qshall,po,z1, &
           psur,ierr5,tcrit,0,xl,cp,   &
           ktf,kts,kte)
      call cup_env_clev(tshall,qes3,qshall,he3,hes3,z3,po,qes3_cup,q3_cup, &
           he3_cup,hes3_cup,z3_cup,po_cup,gamma3_cup,t3_cup,psur,  &
           ierr5,z1,xl,rv,cp,          &
           ktf,kts,kte)
      CALL cup_MAXIMI(HE3_CUP,1,kbmax,K23,ierr5, &
           ktf,kts,kte)
       DO i=its,itf
         if(kpbl(i).gt.5)cap_max3(i)=po_cup(i,kpbl(i))
         IF(ierr5(I).eq.0)THEN
         IF(K23(I).Gt.Kbmax(i))ierr5(i)=2
         if(kpbl(i).gt.5)k23(i)=kpbl(i)
         endif
         ierr5_0(i)=ierr5(i)
       ENDDO
      call cup_kbcon(cap_max_increment,5,k23,kbcon3,he3_cup,hes3_cup, &
           ierr5,kbmax,po_cup,cap_max3, &
!          ierr5,kpbl,po_cup,cap_max3, &
           ktf,kts,kte)
      call cup_up_he(k23,hkb3,z3_cup,cd3,mentr_rate3,he3_cup,hc3, &
           kbcon3,ierr5,dby3,he3,hes3_cup,'shallow', &
           ktf,kts,kte)
      call cup_up_he(k23,hkb3_0,z_cup,cd3,mentr_rate3,he_cup,hc3_0, &
           kbcon3,ierr5,dby3_0,he,hes_cup,'shallow', &
           ktf,kts,kte)
      call cup_ktop(1,dby3,kbcon3,ktop3,ierr5, &
           ktf,kts,kte)
      call cup_up_nms(zu3,z3_cup,mentr_rate3,cd3,kbcon3,ktop3,    &
           ierr5,k23, &
           ktf,kts,kte)
      call cup_up_nms(zu3_0,z_cup,mentr_rate3,cd3,kbcon3,ktop3,    &
           ierr5,k23, &
           ktf,kts,kte)
!
! first calculate aa3_0_cup
!
      call cup_up_aa0(aa3_0,z,zu3_0,dby3_0,GAMMA3_CUP,t_cup, &
           kbcon3,ktop3,ierr5,           &
           ktf,kts,kte)
!
!  now what is necessary for aa3 and feedbacks
!
      call cup_up_moisture('shallow',ierr5,z3_cup,qc3,qrc3,pw3,pwav3, &
           kbcon3,ktop3,cd3,dby3,mentr_rate3,clw_all, &
           qshall,GAMMA3_cup,zu3,qes3_cup,k23,q3_cup,xl,&
           ktf,kts,kte)
      do k=kts,ktf
      do i=its,itf
         cupclws(i,k)=qrc3(i,k)
      enddo
      enddo
      call cup_up_aa0(aa3,z3,zu3,dby3,GAMMA3_CUP,t3_cup, &
           kbcon3,ktop3,ierr5,           &
           ktf,kts,kte)
!     do i=its,itf
!        if(ierr5(i).eq.0)then
!          if(aa3(i).eq.0.)then
!              ierr5(i)=17
!          endif
!        endif
!     enddo
!     call cup_dellabot('shallow',iw,ipr,jpr,q3_cup,ierr5,z3_cup,po,qrcdo,edto, &
!          zdo,cdd,q3,dellaq3,dsubq,j,mentrd_rate,z3,g,&
!          ktf,kts,kte)

      call cup_dellas_3d(iw,ierr5,z3_cup,po_cup,hcdo,edt3,zdo3,cdd,    &
           he3,dellah3,dsubt3,j,mentrd_rate,zu3,g,                     &
           cd3,hc3,ktop3,k23,kbcon3,mentr_rate3,jmin,he3_cup,kdet, &
           k23,ipr,jpr,'shallow',0,                                 &
           ktf,kts,kte)

      call cup_dellas_3d(iw,ierr5,z3_cup,po_cup,qrcdo,edt3,zdo3,cdd, &
           qshall,dellaq3,dsubq3,j,mentrd_rate,zu3,g, &
           cd3,qc3,ktop3,k23,kbcon3,mentr_rate3,jmin,q3_cup,kdet, &
           k23,ipr,jpr,'shallow',0,               &
           ktf,kts,kte    )
              mbdt_s=1.e-1*mbdt_ens(1)
              do k=kts,ktf
              do i=its,itf
                 dellat3(i,k)=0.
                 if(ierr5(i).eq.0)then
                    trash=dsubt3(i,k)
                    XHE3(I,K)=(dsubt3(i,k)+DELLAH3(I,K))*MBDT_S+HE3(I,K)
                    XQ3(I,K)=(dsubq3(i,k)+DELLAQ3(I,K))*MBDT_S+QSHALL(I,K)
                    DELLAT3(I,K)=(1./cp)*(DELLAH3(I,K)-xl*DELLAQ3(I,K))
                    dSUBT3(I,K)=(1./cp)*(dsubt3(i,k)-xl*dsubq3(i,k))
                    XT3(I,K)= (DELLAT3(I,K)+dsubt3(i,k))*MBDT_S+TSHALL(I,K)
                    IF(XQ3(I,K).LT.1.E-08)XQ3(I,K)=1.E-08
                 ENDIF
              enddo
              enddo
      do i=its,itf
      if(ierr5(i).eq.0)then
      XHE3(I,ktf)=HE3(I,ktf)
      XQ3(I,ktf)=QSHALL(I,ktf)
      XT3(I,ktf)=TSHALL(I,ktf)
      IF(XQ3(I,ktf).LT.1.E-08)XQ3(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(iw,xz3,xqes3,xhe3,xhes3,xt3,xq3,po,z1, &
           psur,ierr5,tcrit,2,xl,cp,   &
           ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt3,xqes3,xq3,xhe3,xhes3,xz3,po,xqes3_cup,xq3_cup, &
           xhe3_cup,xhes3_cup,xz3_cup,po_cup,gamma3_cup,xt3_cup,psur,   &
           ierr5,z1,xl,rv,cp,          &
           ktf,kts,kte)
!     
!     
!**************************** static control
!     
!--- moist static energy inside cloud
!
      do i=its,itf
        if(ierr5(i).eq.0)then
          xhkb3(i)=xhe3(i,k23(i))
        endif
      enddo
      call cup_up_he(k23,xhkb3,xz3_cup,cd3,mentr_rate3,xhe3_cup,xhc3, &
           kbcon3,ierr5,xdby3,xhe3,xhes3_cup,'shallow', &
           ktf,kts,kte)
!          
!c--- normalized mass flux profile and CWF
!          
      call cup_up_nms(xzu3,xz3_cup,mentr_rate3,cd3,kbcon3,ktop3,ierr5,k23, &
           ktf,kts,kte)
      call cup_up_aa0(xaa3,xz3,xzu3,xdby3,GAMMA3_CUP,xt3_cup, &
           kbcon3,ktop3,ierr5,           &
           ktf,kts,kte)
!
! now for shallow forcing
!
       do i=its,itf
        xmb3(i)=0.
        xff_shal(1:9)=0.
        if(ierr5(i).eq.0)then
          xkshal=(xaa3(i)-aa3(i))/mbdt_s
          if(xkshal.ge.0.)xkshal=+1.e6
          if(xkshal.gt.-1.e-4 .and. xkshal.lt.0.)xkshal=-1.e-4
          xff_shal(1)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          xff_shal(2)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          xff_shal(3)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          if(aa3_0(i).le.0)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
          if(aa3(i)-aa3_0(i).le.0.)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
! boundary layer QE (from Saulo Freitas)
          blqe=0.
          trash=0.
          if(k23(i).lt.kpbl(i)+1)then
             do k=1,kbcon3(i)-1
                blqe=blqe+100.*dhdt(i,k)*(p_cup(i,k)-p_cup(i,k+1))/g
             enddo
             trash=max((hc3(i,kbcon3(i))-he_cup(i,kbcon3(i))),1.e1)
             xff_shal(7)=max(0.,blqe/trash)
             xff_shal(7)=min(0.1,xff_shal(7))
          else
             xff_shal(7)=0.
          endif
          if((xkshal.lt.-1.1e-04) .and.  &
             ((aa3(i)-aa3_0(i).gt.0.) .or. (xff_shal(7).gt.0)))then
          xff_shal(4)=max(0.,-aa3(i)/(xkshal*tscl_KF))
          xff_shal(4)=min(0.1,xff_shal(4))
          xff_shal(5)=xff_shal(4)
          xff_shal(6)=xff_shal(4)
          else
           xff_shal(4)=0.
           xff_shal(5)=0.
           xff_shal(6)=0.
          endif
!         write(0,888)'i0=',i,j,kpbl(i),blqe,xff_shal(7)
888       format(a3,3(1x,i3),2e12.4)
          xff_shal(8)= xff_shal(7)
          xff_shal(9)= xff_shal(7)
          do k=1,9
           xmb3(i)=xmb3(i)+xff_shal(k)
          enddo
          xmb3(i)=min(.1,xmb3(i)/9.)
          if(xmb3(i).eq.0.)ierr5(i)=22
          if(xmb3(i).lt.0.)then
             ierr5(i)=21
          endif
        endif
        if(ierr5(i).ne.0)then
           k23(i)=0
           kbcon3(i)=0
           ktop3(i)=0
           xmb3(i)=0
           do k=kts,ktf
              outts(i,k)=0.
              outqs(i,k)=0.
           enddo
        else if(ierr5(i).eq.0)then
!
! got the mass flux, sanity check, first for heating rates
!
          trash=0.
          do k=2,ktop3(i)
           trash=max(trash,86400.*(dsubt3(i,k)+dellat3(i,k))*xmb3(i))
          enddo
          if(trash.gt.150.)xmb3(i)=xmb3(i)*150./trash
!
! sanity check on moisture tendencies: do not allow anything that may allow neg tendencies
!
          do k=2,ktop3(i)
           trash=q(i,k)+(dsubq3(i,k)+dellaq3(i,k))*xmb3(i)*dtime
          if(trash.lt.1.e-12)then
! max allowable tendency over tendency that would lead to too small mix ratios
!
            trash=((1.e-12-q(i,k))/dtime)                   &
                  /((dsubq3(i,k)+dellaq3(i,k))*xmb3(i))
            trash=max(0.,trash)
            trash=min(1.,trash)
            xmb3(i)=trash*xmb3(i)
          endif
          enddo
! 
! final tendencies
!
          do k=2,ktop3(i)
           outts(i,k)=(dsubt3(i,k)+dellat3(i,k))*xmb3(i)
           outqs(i,k)=(dsubq3(i,k)+dellaq3(i,k))*xmb3(i)
          enddo
        endif
       enddo

! done shallow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call cup_axx(tcrit,kbmax,z1,p,psur,xl,rv,cp,tx,qx,axx,ierr,    &
           cap_max,cap_max_increment,entr_rate,mentr_rate,&
           j,ktf,kts,kte,ens4)

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwevo,edtmax,edtmin,edtc, &
           ktf,kts,kte)
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=edtc(i,iedt)
         edto(i)=edtc(i,iedt)
         edtx(i)=edtc(i,iedt)
         edt_out(i,j)=edtc(i,2)
         if(high_resolution.eq.1)then
            edt(i)=edtc(i,3)
            edto(i)=edtc(i,3)
            edtx(i)=edtc(i,3)
            edt_out(i,j)=edtc(i,3)
         endif
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
      do i=its,itf
        aad(i)=0.
      enddo
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      call cup_dellabot('deep',iw,ipr,jpr,heo_cup,ierr,zo_cup,po,hcdo,edto, &
           zdo,cdd,heo,dellah,dsubt,j,mentrd_rate,zo,g, &
           ktf,kts,kte)
      call cup_dellabot('deep',iw,ipr,jpr,qo_cup,ierr,zo_cup,po,qrcdo,edto, &
           zdo,cdd,qo,dellaq,dsubq,j,mentrd_rate,zo,g,&
           ktf,kts,kte)
!
!--- 2. everywhere else
!
      call cup_dellas_3d(iw,ierr,zo_cup,po_cup,hcdo,edto,zdo,cdd,    &
           heo,dellah,dsubt,j,mentrd_rate,zuo,g,                     &
           cd,hco,ktop,k22,kbcon,mentr_rate,jmin,heo_cup,kdet, &
           k22,ipr,jpr,'deep',high_resolution,                                 &
           ktf,kts,kte)
!
!-- take out cloud liquid water for detrainment
!
      do k=kts,ktf-1
      do i=its,itf
       scr1(i,k)=0.
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
         scr1(i,k)=qco(i,k)-qrco(i,k)
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k)-po_cup(i,k+1))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*cd(i,k)*dz*zuo(i,k) &
                        *.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k)-po_cup(i,k+1))
         endif
       endif
      enddo
      enddo
      call cup_dellas_3d(iw,ierr,zo_cup,po_cup,qrcdo,edto,zdo,cdd, &
           qo,dellaq,dsubq,j,mentrd_rate,zuo,g, &
           cd,qco,ktop,k22,kbcon,mentr_rate,jmin,qo_cup,kdet, &
           k22,ipr,jpr,'deep',high_resolution,               &
           ktf,kts,kte    )
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(2)
      do i=its,itf
      xaa0_ens(i,1)=0.
      xaa0_ens(i,2)=0.
      xaa0_ens(i,3)=0.
      enddo

      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            trash=dsubt(i,k)
            XHE(I,K)=(dsubt(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ(I,K)=(dsubq(i,k)+DELLAQ(I,K))*MBDT+QO(I,K)
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            dSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k))*MBDT+TN(I,K)
            IF(XQ(I,K).LT.1.E-08)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LT.1.E-08)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(iw,xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,2,xl,cp,   &
           ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
      do i=its,itf
        if(ierr(i).eq.0)then
          xhkb(i)=xhe(i,k22(i))
        endif
      enddo
      call cup_up_he(k22,xhkb,xz_cup,cd,mentr_rate,xhe_cup,xhc, &
           kbcon,ierr,xdby,xhe,xhes_cup,'deep', &
           ktf,kts,kte)
!
!c--- normalized mass flux profile
!
      call cup_up_nms(xzu,xz_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           ktf,kts,kte)
!
!--- moisture downdraft
!
      call cup_dd_nms(xzd,xz_cup,cdd,mentrd_rate,jmin,ierr, &
           1,kdet,z1,                 &
           ktf,kts,kte)
      call cup_dd_he(xhes_cup,xzd,xhcd,xz_cup,cdd,mentrd_rate, &
           jmin,ierr,xhe,dbyd,xhe_cup,&
           ktf,kts,kte)
      call cup_dd_moisture_3d(xzd,xhcd,xhes_cup,xqcd,xqes_cup, &
           xpwd,xq_cup,xz_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           xpwev,bu,xqrcd,xq,xhe,xt_cup,3,xl,high_resolution, &
           ktf,kts,kte)

!
!------- MOISTURE updraft
!
      call cup_up_moisture('deep',ierr,xz_cup,xqc,xqrc,xpw,xpwav, &
           kbcon,ktop,cd,xdby,mentr_rate,clw_all, &
           xq,GAMMA_cup,xzu,xqes_cup,k22,xq_cup,xl, &
           ktf,kts,kte)
!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
      do 200 nens=1,maxens
      do i=its,itf 
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
                                 +edto(i)*pwdo(i,k)             &
                                    +pwo(i,k) 
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
                                 +.5*edto(i)*pwdo(i,k)          &
                                 +  pwo(i,k)
                 else
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)+edto(i)*pwdo(i,k)
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,j,nall+7).lt.1.e-6)then
            ierr(i)=18
            do nens3=1,maxens3
               pr_ens(i,j,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,j,nall+nens3).lt.1.e-4)then
            pr_ens(i,j,nall+nens3)=0.
           endif
         enddo
         endif
      enddo
 200  continue
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero
!
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
           ktf,kts,kte)
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
      enddo
      call cup_kbcon(cap_max_increment,2,k22x,kbconx,heo_cup, &
           heso_cup,ierr2,kbmax,po_cup,cap_max, &
           ktf,kts,kte)
      call cup_kbcon(cap_max_increment,3,k22x,kbconx,heo_cup, &
           heso_cup,ierr3,kbmax,po_cup,cap_max, &
           ktf,kts,kte)
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,                 &
           iedt,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           ichoice,edt_out,     &
           high_resolution,ktf,kts,kte,ens4,ktau,ipr)
!
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           subt_ens(i,k,iedt)=dsubt(i,k)
           subq_ens(i,k,iedt)=dsubq(i,k)
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edt(i)*pwdo(i,k)
        else 
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!
      call cup_output_ens_3d(xf_ens,ierr,dellat_ens,dellaq_ens, &
           dellaqc_ens,subt_ens,subq_ens,subt,subq,outt,     &
           outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
           j,'deep',ierr2,ierr3,         &
           pr_ens,qce,                    &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                &
           APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,   &
           ktf,kts,kte)
      k=1
      do i=its,itf
          if(ierr(i).eq.0.and.ierr5(i).eq.0.and.kbcon(i).lt.ktop3(i)+1)then
!            write(0,*)'both ier and ier5=0 at i,j=',i,j,kbcon(i),ktop3(i)
             if(high_resolution.eq.1)then
                outts(i,kts:kte)=0.
                outqs(i,kts:kte)=0.
             endif
          elseif (ierr5(i).eq.0)then
!            write(0,*)'ier5=0 at i,j=',i,j,k23(i),ktop3(i)
          endif

           PRE(I)=MAX(PRE(I),0.)
!           if(i.eq.ipr.and.j.eq.jpr)then
!             write(0,*)'i,j,pre(i),aa0(i),aa1(i)'
!             write(0,*)i,j,pre(i),aa0(i)
!           endif
      enddo
!
!---------------------------done------------------------------
!
!      do i=its,itf
!        if(ierr(i).eq.0)then
!       if(i.eq.ipr.and.j.eq.jpr)then
!         write(0,*)'on output, pre =',pre(i),its,itf,kts,ktf
!         do k=kts,ktf
!           write(0,*)z(i,k),outt(i,k)*86400.,subt(i,k)*86400.
!         enddo
!         write(0,*)i,j,(axx(i,k),k=1,ens4)
!       endif
!       endif
!      enddo
!     print *,'ierr(i) = ',ierr(i),pre(i)

   END SUBROUTINE CUP_enss_3d


   SUBROUTINE cup_dd_aa0(edt,ierr,aa0,jmin,gamma_cup,t_cup, &
              hcd,hes_cup,z,zd,                             &
              ktf,kts,kte                    )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! aa0 cloud work function for downdraft
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! hcd = moist static energy in downdraft
  ! edt = epsilon
  ! zd normalized downdraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zd,gamma_cup,t_cup,hes_cup,hcd
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,kk
     real                                 ::                           &
        dz
!
       do i=its,itf
        aa0(i)=0.
       enddo
!
!??    DO k=kts,kte-1
       DO k=kts,ktf-1
       do i=its,itf
         IF(ierr(I).eq.0.and.k.lt.jmin(i))then
         KK=JMIN(I)-K
!
!--- ORIGINAL
!
         DZ=(Z(I,KK)-Z(I,KK+1))
         AA0(I)=AA0(I)+zd(i,kk)*EDT(I)*DZ*(9.81/(1004.*T_cup(I,KK))) &
            *((hcd(i,kk)-hes_cup(i,kk))/(1.+GAMMA_cup(i,kk)))
         endif
      enddo
      enddo

   END SUBROUTINE CUP_dd_aa0


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,edtc,               &
              ktf,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        us,vs,z,p
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev
     real                                                              &
        ,intent (in   )                   ::                           &
        edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc,vws
     real,    dimension (its:ite)         ::                           &
      vshear

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo

       do i=its,itf
       IF (ierr(i).ne.0) cycle
       vws = 0.
       do kk = kbcon(i), ktop(i)
             vws = vws + &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +    abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
       end do
       vshear(i) = 1.e3 * vws / (p(i,kbcon(i)) - p(i,ktop(i)+1))
       end do

       do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.1.)pef=1.
            if(pef.lt.0.)pef=0.
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
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
            EDT(I)=1.-.5*(pefb+pef)
!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr,       &
              jmin,ierr,he,dby,he_cup,                       &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cdd,zd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hcd,dby
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dz


      do k=kts+1,ktf
      do i=its,itf
      dby(i,k)=0.
      IF(ierr(I).eq.0)then
         hcd(i,k)=hes_cup(i,k)
      endif
      enddo
      enddo
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      hcd(i,k)=hes_cup(i,k)
      dby(i,k)=hcd(i,jmin(i))-hes_cup(i,k)
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         HCD(i,Ki)=(HCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*HE(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
         dby(i,ki)=HCD(i,Ki)-hes_cup(i,ki)
      enddo
!
      endif
!--- end loop over i
100    continue


   END SUBROUTINE cup_dd_he


   SUBROUTINE cup_dd_moisture_3d(zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,high_resolution,           &
              ktf,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte,high_resolution
  ! cdd= detrainment function 
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,cdd,gamma_cup,q,he 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      if(high_resolution.eq.1)qcd(i,k)=.5*(qes_cup(i,k)+q_cup(i,k))
      qrcd(i,k)=qes_cup(i,k)
      pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
      pwev(i)=pwev(i)+pwd(i,jmin(i))
      qcd(i,k)=qes_cup(i,k)
!
      DH=HCD(I,k)-HES_cup(I,K)
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*q(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
!
!--- to be negatively buoyant, hcd should be smaller than hes!
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)dqeva=0.
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
      enddo
!
!--- end loop over i
       if(pwev(I).ge.0.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
       if(BU(I).ge.0.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_3d


   SUBROUTINE cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr,        &
              itest,kdet,z1,                                 &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte
  ! z_cup = height of cloud model level
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! itest = flag to whether to calculate cdd
  
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        z1 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr 
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin,kdet
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
                                                                 ierr
   ! zd is the normalized downdraft mass flux
   ! cdd is the downdraft detrainmen function

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
                                                             zd
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
                                                             cdd
!
!  local variables in this routine
!

     integer                              ::                           &
                                                  i,k,ki
     real                                 ::                           &
                                            a,perc,dz

!
!--- perc is the percentage of mass left when hitting the ground
!
      perc=.03

      do k=kts,ktf
      do i=its,itf
         zd(i,k)=0.
         if(itest.eq.0)cdd(i,k)=0.
      enddo
      enddo
      a=1.-perc
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      zd(i,jmin(i))=1.
!
!--- integrate downward, specify detrainment(cdd)!
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         if(ki.le.kdet(i).and.itest.eq.0)then
           cdd(i,ki)=entr+(1.- (a*(z_cup(i,ki)-z1(i)) &
                     +perc*(z_cup(i,kdet(i))-z1(i)) ) &
                         /(a*(z_cup(i,ki+1)-z1(i)) &
                      +perc*(z_cup(i,kdet(i))-z1(i))))/dz
         endif
         zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki))*dz)
      enddo
!
      endif
!--- end loop over i
100    continue

   END SUBROUTINE cup_dd_nms


   SUBROUTINE cup_dellabot(name,iw,ipr,jpr,he_cup,ierr,z_cup,p_cup,  &
              hcd,edt,zd,cdd,he,della,subs,j,mentrd_rate,z,g,     &
              ktf,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        iw,ktf,kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr
      character (*), intent (in)        ::                           &
       name
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della,subs
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,z,he_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

      integer i
      real detdo,detdo1,detdo2,entdo,dp,dz,subin,                      &
      totmas
!
!
!     if(name.eq.'shallow')then
!        edt(:)=0.
!        cdd(:,:)=0.
!     endif
      do 100 i=its,itf
      della(i,1)=0.
      subs(i,1)=0.
      if(ierr(i).ne.0)go to 100
      dz=z_cup(i,2)-z_cup(i,1)
      DP=100.*(p_cup(i,1)-P_cup(i,2))
      detdo1=edt(i)*zd(i,2)*CDD(i,1)*DZ
      detdo2=edt(i)*zd(i,1)
      entdo=edt(i)*zd(i,2)*mentrd_rate*dz
      subin=-EDT(I)*zd(i,2)
      detdo=detdo1+detdo2-entdo+subin
      DELLA(I,1)=(detdo1*.5*(HCD(i,1)+HCD(i,2)) &
                 +detdo2*hcd(i,1) &
                 +subin*he_cup(i,2) &
                 -entdo*he(i,1))*g/dp
      SUBS(I,1)=0.
!      if(i.eq.ipr.and.j.eq.jpr)then
!       write(0,*)'db1',della(i,1),subs(i,1),subin,entdo
!       write(0,*)'db2',detdo1,detdo2,detdo1+detdo2-entdo+subin
!      endif
 100  CONTINUE

   END SUBROUTINE cup_dellabot


   SUBROUTINE cup_dellas_3d(iw,ierr,z_cup,p_cup,hcd,edt,zd,cdd,              &
              he,della,subs,j,mentrd_rate,zu,g,                             &
              cd,hc,ktop,k22,kbcon,mentr_rate,jmin,he_cup,kdet,kpbl,   &
              ipr,jpr,name,high_res,                                            &
              ktf,kts,kte                               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        iw,ktf,kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr,high_res
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della,subs
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,hc,cd,zu,he_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate,mentr_rate
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22,jmin,kdet,kpbl
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character (*), intent (in)        ::                           &
       name
!
!  local variables in this routine
!

      integer i,k,kstart
      real entdo,dp,dz,subin,detdo,entup,                            &
      detup,subdown,entdoj,entupk,detupk,totmas
!
      i=ipr
      kstart=kts+1
      if(name.eq.'shallow')kstart=kts
       DO K=kstart,ktf
       do i=its,itf
          della(i,k)=0.
          subs(i,k)=0.
       enddo
       enddo
!
! no downdrafts for shallow convection
!
       DO 100 k=kts+1,ktf-1
       DO 100 i=its,ite
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         if(k.lt.k22(i)-1.and.name.eq.'shallow')GO TO 100
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
         DZ=Z_cup(I,K+1)-Z_cup(I,K)
         detdo=edt(i)*CDD(i,K)*DZ*ZD(i,k+1)
         entdo=edt(i)*mentrd_rate*dz*zd(i,k+1)
!3d        subin=zu(i,k+1)-zd(i,k+1)*edt(i)
         subin=-zd(i,k+1)*edt(i)
         entup=0.
         detup=0.
         if(k.ge.kbcon(i).and.k.lt.ktop(i))then
            entup=mentr_rate*dz*zu(i,k)
            detup=CD(i,K+1)*DZ*ZU(i,k)
         endif
!3d        subdown=(zu(i,k)-zd(i,k)*edt(i))
         subdown=-zd(i,k)*edt(i)
         entdoj=0.
         entupk=0.
         detupk=0.
!
         if(k.eq.jmin(i))then
         entdoj=edt(i)*zd(i,k)
         endif

         if(k.eq.k22(i)-1)then
         entupk=zu(i,kpbl(i))
         subin=zu(i,k+1)-zd(i,k+1)*edt(i)
         if(high_res.eq.1)subin=-zd(i,k+1)*edt(i)
!        subin=-zd(i,k+1)*edt(i)
         endif

         if(k.gt.kdet(i))then
            detdo=0.
         endif

         if(k.eq.ktop(i)-0)then
         detupk=zu(i,ktop(i))
         subin=0.
!
! this subsidene for ktop now in subs term!
!        subdown=zu(i,k)
         subdown=0.
         endif
         if(k.lt.kbcon(i))then
            detup=0.
         endif
!C
!C--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
!C
         totmas=subin-subdown+detup-entup-entdo+ &
                 detdo-entupk-entdoj+detupk
!         if(j.eq.jpr.and.i.eq.ipr)print *,'k,totmas,sui,sud = ',k,
!     1   totmas,subin,subdown
!         if(j.eq.jpr.and.i.eq.ipr)print *,'updr stuff = ',detup,
!     1      entup,entupk,detupk
!         if(j.eq.jpr.and.i.eq.ipr)print *,'dddr stuff = ',entdo,
!     1      detdo,entdoj
         if(abs(totmas).gt.1.e-6)then
!           print *,'*********************',i,j,k,totmas,name
!           print *,kpbl(i),k22(i),kbcon(i),ktop(i)
!c          print *,'updr stuff = ',subin,
!c    1      subdown,detup,entup,entupk,detupk
!c          print *,'dddr stuff = ',entdo,
!c    1      detdo,entdoj
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
         endif
         dp=100.*(p_cup(i,k-1)-p_cup(i,k))
         della(i,k)=(detup*.5*(HC(i,K+1)+HC(i,K)) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entup*he(i,k) &
                    -entdo*he(i,k) &
                    +subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detupk*(hc(i,ktop(i))-he_cup(i,ktop(i)))    &
                    -entupk*he_cup(i,k22(i)) &
                    -entdoj*he_cup(i,jmin(i)) &
                     )*g/dp

           if(high_res.eq.1)then
! the first term includes entr and detr into/from updraft as well as (entup-detup)*he(i,k) from
!  neighbouring point, to make things mass consistent....
!            if(k.ge.k22(i))then
                della(i,k)=(                          &
                    detup*.5*(HC(i,K+1)+HC(i,K))-entup*he(i,k)+(entup-detup)*he(i,k) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entdo*he(i,k) &
                    +subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detupk*(hc(i,ktop(i))-he(i,ktop(i)))    &
                    -entdoj*he_cup(i,jmin(i)) &
                    -entupk*he_cup(i,k22(i))+entupk*he(i,k) &
                     )*g/dp
!             else if(k.eq.k22(i)-1)then
!                  della(i,k)=(-entupk*he_cup(i,k22(i))+entupk*he(i,k))*g/dp
           endif
!3d        subin=zu(i,k+1)-zd(i,k+1)*edt(i)
!
! updraft subsidence only
!
         if(k.ge.k22(i).and.k.lt.ktop(i))then
         subs(i,k)=(zu(i,k+1)*he_cup(i,k+1) &
                    -zu(i,k)*he_cup(i,k))*g/dp
!        else if(k.eq.ktop(i))then
!        subs(i,k)=-detupk*he_cup(i,k)*g/dp
         endif
!
! in igh res case, subsidence terms are for meighbouring points only. This has to be 
! done mass consistent with the della term
         if(high_res.eq.1)then
            if(k.ge.k22(i).and.k.lt.ktop(i))then
               subs(i,k)=(zu(i,k+1)*he_cup(i,k+1)-zu(i,k)*he_cup(i,k)-(entup-detup)*he(i,k))*g/dp
            else if(k.eq.ktop(i))then
               subs(i,k)=detupk*(he(i,ktop(i))-he_cup(i,ktop(i)))*g/dp
            else if(k.eq.k22(i)-1)then
               subs(i,k)=(entupk*he(i,k)-entupk*he_cup(i,k))*g/dp
         endif
         endif

 100  CONTINUE

   END SUBROUTINE cup_dellas_3d


   SUBROUTINE cup_env(iw,z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              ktf,kts,kte                     )

     use consts_coms, only: t00, eps_vap
     use therm_lib,   only: eslf, esif

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        iw,ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z,q
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
!     real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar

!     HT(1)=XL/CP
!     HT(2)=2.834E6/CP
!     BE(1)=.622*HT(1)/.286
!     AE(1)=BE(1)/273.+ALOG(610.71)
!     BE(2)=.622*HT(2)/.286
!     AE(2)=BE(2)/273.+ALOG(610.71)

      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
!        IPH=1
!        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
!       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
!       QES(I,K)=.622*E/(100.*P(I,K)-(1.-0.622)*E)

        if ( t(i,k) > tcrit ) then
           e = eslf( t(i,k) - t00 )
        else
           e = esif( t(i,k) - t00 )
        endif

        qes(i,k) = max( eps_vap*e / (100.*p(i,k) - (1.0-eps_vap)*e), 1.e-8)

!       IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
!       IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
        IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.ne.2)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1) = max(0.,Z1(I)) + ALOG(PSUR(I)/P(I,1)) * 287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K) = Z(I,K-1) + ALOG(P(I,K-1)/P(I,K)) * 287.*TVBAR/9.81
           endif
         enddo
         enddo
      else
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.eq.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              ktf,kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,iedt,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,      &
              icoic,edt_out,            &
              high_resolution,ktf,kts,kte,ens4,ktau,ipr               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte,ens4,high_resolution,ktau,ipr
     integer, intent (in   )              ::                           &
        j,iedt
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout   )                   ::                           &
        edt_out
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte,1:ens4)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,xland
     real,    dimension (its:ite,1:ens4)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character (*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3
     real                                 ::                           &
       fens4,a1,a_ave,a_max,a_min,xff0,xxx,xomg

     integer :: nall2,ixxx
     integer,  dimension (33) :: seed

     integer, parameter :: irandom = 0

       seed=0
       seed(2)=j
       seed(3)=ktau
       nens=0
       fens4=float(ens4)

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
!
!---
!
             if(name.eq.'deeps')then

                a_ave = sum(axx(1,1:ens4))/fens4
                a_ave = min(a_ave,aa1(i))
                a_ave = max(0.,a_ave)

                a_min = minval(axx(1,1:ens4))
                a_min = max(0.,a_min)

                a_max = maxval(axx(1,1:ens4))
                a_max = max(0.,a_max)

                do ne=1,16
                  xff_ens3(ne)=0.
                enddo

                xff0= (AA1(I)-AA0(I))/DTIME
                if(high_resolution.eq.1)xff0= (a_ave-AA0(I))/DTIME

                xff_ens3(1)=(AA1(I)-AA0(I))/dtime
                xff_ens3(2)=(AA1(I)-AA0(I))/dtime

                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(3)=(axx(i,ixxx)-AA0(I))/dtime
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(13)=(axx(i,ixxx)-AA0(I))/dtime
                else
                   xff_ens3(3) =(a_ave-AA0(I))/dtime
                   xff_ens3(13)=(a_min-AA0(I))/dtime
                   xff_ens3(16)=(a_max-AA0(I))/dtime
                endif

                if(high_resolution.eq.1)then
                   xff_ens3(1)=(a_ave-AA0(I))/dtime
                   xff_ens3(2)=(a_ave-AA0(I))/dtime
                   xff_ens3(3)=(a_ave-AA0(I))/dtime
                   xff_ens3(13)=(a_ave-AA0(I))/dtime
                   xff_ens3(16)=(a_ave-AA0(I))/dtime
                endif
!   
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(14)=0.
                do ne=1,ens4
                  xff_ens3(14)=xff_ens3(14)-omeg(i,k22(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(14).lt.0.)xff_ens3(14)=0.
                xff_ens3(5)=0.
                do ne=1,ens4
                  xff_ens3(5)=xff_ens3(5)-omeg(i,kbcon(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
!  
! minimum below kbcon
!
                if(high_resolution.eq.0)then
                   xff_ens3(4)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.lt.xff_ens3(4))xff_ens3(4)=xomg
                   enddo
                   enddo
                   if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
!
! max below kbcon
                   xff_ens3(6)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
                   enddo
                   if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(5)=min(xff_ens3(5),xff_ens3(14))
                   xff_ens3(4)=xff_ens3(5)
                   xff_ens3(6)=xff_ens3(5)
                endif
!
!--- more like Krishnamurti et al.; pick max, min,  and average values
!
                xff_ens3(7)=mconv(i,1)
                xff_ens3(8)=mconv(i,1)
                xff_ens3(9)=mconv(i,1)
                if(ens4.gt.1)then
                   do ne=2,ens4
                      if (mconv(i,ne).gt.xff_ens3(7))xff_ens3(7)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      if (mconv(i,ne).lt.xff_ens3(8))xff_ens3(8)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      xff_ens3(9)=xff_ens3(9)+mconv(i,ne)
                   enddo
                   xff_ens3(9)=xff_ens3(9)/fens4
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(7)=xff_ens3(9)
                   xff_ens3(8)=xff_ens3(9)
                   xff_ens3(15)=xff_ens3(9)
                endif
!
                if(high_resolution.eq.0)then
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(15)=mconv(i,ixxx)
                else
                   xff_ens3(15)=mconv(i,1)
                endif
                endif
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=AA0(I)/(60.*40.)
                xff_ens3(11)=AA1(I)/(60.*40.)
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AXX(I,ixxx)/(60.*40.)
                else
                   xff_ens3(12)=A_AVE/(60.*40.)
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(11)=A_AVE/(60.*40.)
                   xff_ens3(12)=A_AVE/(60.*40.)
                endif

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2)
                   if(xk(nens).le.0.0.and.xk(nens).gt.-1.e-2) &
                           xk(nens)=-1.e-2
                   if(xk(nens).gt.0.0.and.xk(nens).lt.1.e-2) &
                           xk(nens)=1.e-2
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
! over water, enforce small cap for some of the closures
!
                if(xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =0.
                      xff_ens3(2) =0.
                      xff_ens3(3) =0.
                      xff_ens3(10) =0.
                      xff_ens3(11) =0.
                      xff_ens3(12) =0.
                      xff_ens3(7) =0.
                      xff_ens3(8) =0.
                      xff_ens3(9) =0.
                      xff_ens3(13) =0.
                      xff_ens3(15) =0.
                 endif
                endif
!
! end water treatment
!

!
!--- special treatment for stability closures
!
                      if(xff0.gt.0. .and. xk(ne).lt.0.)then
                         if(xff_ens3(1).gt.0.)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))
                         if(xff_ens3(2).gt.0.)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne))
                         if(xff_ens3(3).gt.0.)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne))
                         if(xff_ens3(13).gt.0.)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                         if(xff_ens3(16).gt.0.)xf_ens(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne))
                      else
                         xf_ens(i,j,nall+1)=0.
                         xf_ens(i,j,nall+2)=0.
                         xf_ens(i,j,nall+3)=0.
                         xf_ens(i,j,nall+13)=0.
                         xf_ens(i,j,nall+16)=0.
                      endif
!
!--- following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4))
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5))
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6))
                         xf_ens(i,j,nall+14)=max(0.,xff_ens3(14))
                         a1=max(1.e-2,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+15))
                         xf_ens(i,j,nall+15)=max(0.,xff_ens3(15)/a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
                            xf_ens(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))
                            xf_ens(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))
                         else
                            xf_ens(i,j,nall+10)=0.
                            xf_ens(i,j,nall+11)=0.
                            xf_ens(i,j,nall+12)=0.
                         endif
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
                endif
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif

 350            continue
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
             enddo
          endif
 100   continue

!mjo    Make sure any individual closures are not going crazy
        xf_ens(:,:,:) = min(xf_ens(:,:,:), 2.*xmbmax)

   END SUBROUTINE cup_forcing_ens_3d


   SUBROUTINE cup_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              ierr,kbmax,p_cup,cap_max,                         &
              ktf,kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max,cap_inc
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        pbcdif,plus,hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)ierr(i)=3
!        if(iloop.lt.4)ierr(i)=997
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=HE_cup(I,K22(I))
      if(iloop.eq.5)then
       do k=1,k22(i)
         hetest=max(hetest,he_cup(i,k))
       enddo
      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))GO TO 31

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).eq.0)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=25.
      if(iloop.eq.5.and.cap_max(i).gt.25.)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ilo,dby,kbcon,ktop,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,kstop

       do i=its,itf
         maxx(i)=ks
         if(ierr(i).eq.0)then
            kstop = max( ks, ke(i) )
            maxx(i) = ks + maxloc( array(i,ks:kstop), 1 ) - 1
         endif
       enddo

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,kstop

       DO i=its,itf
         KT(I)=KS(I)
         if(ierr(i).eq.0)then
            kstop = max( ks(i), kend(i) )
            kt(i) = ks(i) + minloc( array(i,ks(i):kstop), 1 ) - 1
         endif
       enddo

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_output_ens_3d(xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,ierr2,ierr3,pr_ens,qce,             &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              ktf,kts,kte)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
     integer, intent (in   )              ::                           &
        j
  ! xf_ens = ensemble mass fluxes
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
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout)                   ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI 

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas,qce
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite,kts:kte,1:maxens2)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n
     real                                 ::                           &
        ddtes,dtt,dtq,dtqc,dtpw,tuning,clos_wei
     real                                 ::                           &
        dtts,dtqs
!    real,    dimension (its:ite)         ::                           &
!      xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm
     real, dimension (its:ite,5) :: xmb_w

!
      character (*), intent (in)        ::                           &
       name
!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

     tuning=0.
!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        subt(i,k)=0.
        subq(i,k)=0.
        sub_mas(i,k)=0.
        qce(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
!       xfac1(i)=0.
!       xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*maxens*maxens2*maxens3+1,iens*maxens*maxens2*maxens3
           if(pr_ens(i,j,n).le.0.)then
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!
       call massflx_stats(xf_ens,      &
            xmb_ave,xmb_std,xmb_cur,xmb_ske,j,ierr,1,    &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            ktf,kts,kte                   )
       xmb_w=0.
       call massflx_stats(pr_ens,  &
            pr_ave,pr_std,pr_cur,pr_ske,j,ierr,2,        &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            ktf,kts,kte                   )
!
!-- now do feedback
!
      ddtes=100.
      do i=its,itf
        if(ierr(i).eq.0)then
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
         xmb(i)=max(.1*xmb_ave(i),xmb_ave(i)-tuning*xmb_std(i))
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
           if (xland1(i).lt.0.5)xmb(i)=xmb(i)*clos_wei
           if(xmb(i).le.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
!mjo limit mass flux
           xmb(i) = min(xmb(i),xmbmax)
!          xfac1(i)=xmb(i)
!          xfac2(i)=xmb(i)

        endif
!       if(weight(1).lt.-100.)xfac1(i)=xmb_ave(i)
!       if(weight(1).lt.-100.)xfac2(i)=xmb_ave(i)
      ENDDO

      DO k=kts,ktf
      do i=its,itf
            dtt=0.
            dtts=0.
            dtq=0.
            dtqs=0.
            dtqc=0.
            dtpw=0.
        IF(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,maxens2
              dtt=dtt+dellat(i,k,n)
              dtts=dtts+subt_ens(i,k,n)
              dtq=dtq+dellaq(i,k,n)
              dtqs=dtqs+subq_ens(i,k,n)
              dtqc=dtqc+dellaqc(i,k,n)
              dtpw=dtpw+pw(i,k,n)
           enddo
           OUTTEM(I,K)=XMB(I)*dtt/float(maxens2)
           SUBT(I,K)=XMB(I)*dtts/float(maxens2)
           OUTQ(I,K)=XMB(I)*dtq/float(maxens2)
           SUBQ(I,K)=XMB(I)*dtqs/float(maxens2)
           OUTQC(I,K)=XMB(I)*dtqc/float(maxens2)
           qce(i,k)=XMB(I)*dtpw/float(maxens2)
           PRE(I)=PRE(I)+XMB(I)*dtpw/float(maxens2)
           sub_mas(i,k)=zu(i,k)*xmb(i)
        endif
      enddo
      enddo

!     do i=its,itf
!       if(ierr(i).eq.0)then
!       do k=(iens-1)*maxens*maxens2*maxens3+1,iens*maxens*maxens2*maxens3
!         xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
!       enddo
!       endif
!     ENDDO

   END SUBROUTINE cup_output_ens_3d


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0


   SUBROUTINE cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,     &
              kbcon,ierr,dby,he,hes_cup,name,                &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte
      character (*), intent (in)        ::                           &
       name
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hc,dby
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        hkb
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz
!
!--- moist static energy inside cloud
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
         hkb(i)=0.
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
          hkb(i)=he_cup(i,k22(i))
          if(name.eq.'shallow')then
             do k=1,k22(i)
               hkb(i)=max(hkb(i),he_cup(i,k))
             enddo
          endif
          do k=1,k22(i)
              hc(i,k)=he_cup(i,k)
          enddo
          do k=k22(i),kbcon(i)-1
              hc(i,k)=hkb(i)
          enddo
          k=kbcon(i)
          hc(i,k)=hkb(i)
          DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
        endif
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(k.gt.kbcon(i).and.ierr(i).eq.0)then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr* &
                DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
      enddo

      enddo

   END SUBROUTINE cup_up_he


   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,                  &
              q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,          &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        q,zu,gamma_cup,qe_cup,dby,qes_cup,z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        mentr_rate,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character (*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav
!
!  local variables in this routine
!

     integer                              ::                           &
        iall,i,k
     real                                 ::                           &
        dh,qrch,c0,dz
!
        iall=0
        c0=.002
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')c0=0.
        do i=its,itf
          pwav(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          qrc(i,k)=0.
        enddo
        enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=k22(i),kbcon(i)-1
        qc(i,k)=qe_cup(i,k22(i))
      enddo
      endif
      enddo

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Lt.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
        QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate* &
                DZ*Q(i,K-1))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------- LIQUID WATER CONTENT IN cloud after rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
        if(qrc(i,k).lt.0.)then
          qrc(i,k)=0.
        endif
!
!-------   3.Condensation
!
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
!
!----- set next level
!
         QC(I,K)=QRC(I,K)+qrch
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
 100     CONTINUE

   END SUBROUTINE cup_up_moisture


   SUBROUTINE cup_up_nms(zu,z_cup,entr,cd,kbcon,ktop,ierr,k22,  &
              ktf,kts,kte                        )

   IMPLICIT NONE

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kts,kte
  ! cd= detrainment function 
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
         entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
         ierr
   ! zu is the normalized mass flux

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
         zu
!
!  local variables in this routine
!

     integer                              ::                           &
         i,k
     real                                 ::                           &
         dz
!
!   initialize for this go around
!
       do k=kts,ktf
       do i=its,itf
         zu(i,k)=0.
       enddo
       enddo
!
! do normalized mass budget
!
       do i=its,itf
          IF(ierr(I).eq.0)then
             do k=k22(i),kbcon(i)
               zu(i,k)=1.
             enddo
             DO K=KBcon(i)+1,KTOP(i)
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k))*DZ)
             enddo
          endif
       enddo

   END SUBROUTINE cup_up_nms


   SUBROUTINE massflx_stats(xf_ens, &
              xt_ave,xt_std,xt_cur,xt_ske,j,ierr,itest,           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                  &
              APR_CAPMA,APR_CAPME,APR_CAPMI,                      &
              pr_gr,pr_w,pr_mc,pr_st,pr_as,                       &
              pr_capma,pr_capme,pr_capmi,                         &
              ktf,kts,kte)

   IMPLICIT NONE

   integer, intent (in   )              ::                                    &
                     j,itest
   INTEGER,      INTENT(IN   ) ::                                             &
                                  ktf,kts,kte


     real, dimension (its:ite)                                                &
         , intent(inout) ::                                                   &
           xt_ave,xt_cur,xt_std,xt_ske
     integer, dimension (its:ite), intent (in) ::                             &
           ierr
     real, dimension (its:ite,jts:jte,1:ensdim)                               &
         , intent(in   ) ::                                                   &
           xf_ens
     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                                 &
           APR_CAPMA,APR_CAPME,APR_CAPMI
     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           pr_gr,pr_w,pr_mc,pr_st,pr_as,                                      &
           pr_capma,pr_capme,pr_capmi

!
! local stuff
!
     real, dimension (its:ite , 1:maxens3 )       ::                          &
           x_ave,x_cur,x_std,x_ske
     real, dimension (its:ite , 1:maxens  )       ::                          &
           x_ave_cap

      integer :: i,k
      integer :: num,kk,num2,iedt
      real :: a3,a4

      num=ensdim/maxens3
      num2=ensdim/maxens
      if(itest.eq.1)then
      do i=its,ite
       pr_gr(i,j) =  0.
       pr_w(i,j) =  0.
       pr_mc(i,j) = 0.
       pr_st(i,j) = 0.
       pr_as(i,j) = 0.
       pr_capma(i,j) =  0.
       pr_capme(i,j) = 0.
       pr_capmi(i,j) = 0.
      enddo
      endif

      do k=1,maxens
      do i=its,ite
        x_ave_cap(i,k)=0.
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
      enddo
      enddo
      do i=its,ite
        xt_ave(i)=0.
        xt_std(i)=0.
        xt_ske(i)=0.
        xt_cur(i)=0.
      enddo
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)+xf_ens(i,j,maxens3*(kk-1)+k)
        endif
      enddo
      enddo
      enddo
      do iedt=1,maxens2
      do k=1,maxens
      do kk=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)                               &
            +xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
        endif
      enddo
      enddo
      enddo
      enddo
      do k=1,maxens
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)/float(num2)
        endif
      enddo
      enddo

      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)/float(num)
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)+x_ave(i,k)
        endif
      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)/float(maxens3)
        endif
      enddo
!
!--- now do std, skewness,curtosis
!
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!       print *,i,j,k,kk,x_std(i,k),xf_ens(i,j,maxens3*(kk-1)+k),x_ave(i,k)
        x_std(i,k)=x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
        x_ske(i,k)=x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
        x_cur(i,k)=x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
        endif
      enddo
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        xt_std(i)=xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
        xt_ske(i)=xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
        xt_cur(i)=xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
           x_std(i,k)=x_std(i,k)/float(num)
           a3=max(1.e-6,x_std(i,k))
           x_std(i,k)=sqrt(a3)
           a3=max(1.e-6,x_std(i,k)**3)
           a4=max(1.e-6,x_std(i,k)**4)
           x_ske(i,k)=x_ske(i,k)/float(num)/a3
           x_cur(i,k)=x_cur(i,k)/float(num)/a4
        endif
!       print*,'                               '
!       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
!       print*,'statistics for closure number ',k
!       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
!       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
!       print*,'                               '

      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0.and.xt_std(i).gt.0.)then
           xt_std(i)=xt_std(i)/float(maxens3)
           a3=max(1.e-6,xt_std(i))
           xt_std(i)=sqrt(a3)
           a3=max(1.e-6,xt_std(i)**3)
           a4=max(1.e-6,xt_std(i)**4)
           xt_ske(i)=xt_ske(i)/float(maxens3)/a3
           xt_cur(i)=xt_cur(i)/float(maxens3)/a4
!       print*,'                               '
!       print*,'Total ensemble independent statistics at i =',i
!       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
!       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
!       print*,'                               '
!
!  first go around: store massflx for different closures/caps
!
      if(itest.eq.1)then
       pr_gr(i,j) = .25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))
       pr_w(i,j) = .25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))
       pr_mc(i,j) = .25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))
       pr_st(i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
       pr_as(i,j) = x_ave(i,16)
       pr_capma(i,j) = x_ave_cap(i,1)
       pr_capme(i,j) = x_ave_cap(i,2)
       pr_capmi(i,j) = x_ave_cap(i,3)
!
!  second go around: store preciprates (mm/hour) for different closures/caps
!
        else if (itest.eq.2)then
       APR_GR(i,j)=.25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))*      &
                  3600.*pr_gr(i,j) +APR_GR(i,j)
       APR_W(i,j)=.25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))*       &
                  3600.*pr_w(i,j) +APR_W(i,j)
       APR_MC(i,j)=.25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))*      &
                  3600.*pr_mc(i,j) +APR_MC(i,j)
       APR_ST(i,j)=.333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))*   &
                  3600.*pr_st(i,j) +APR_ST(i,j)
       APR_AS(i,j)=x_ave(i,16)*                       &
                  3600.*pr_as(i,j) +APR_AS(i,j)
       APR_CAPMA(i,j) = x_ave_cap(i,1)*                          &
                  3600.*pr_capma(i,j) +APR_CAPMA(i,j)
       APR_CAPME(i,j) = x_ave_cap(i,2)*                          &
                  3600.*pr_capme(i,j) +APR_CAPME(i,j)
       APR_CAPMI(i,j) = x_ave_cap(i,3)*                          &
                  3600.*pr_capmi(i,j) +APR_CAPMI(i,j)
        endif
        endif
      enddo

   END SUBROUTINE massflx_stats

   SUBROUTINE cup_axx(tcrit,kbmax,z1,p,psur,xl,rv,cp,tx,qx,axx,ierr,    &
           cap_max,cap_max_increment,entr_rate,mentr_rate,&
           j,ktf,kts,kte,ens4)

     use consts_coms, only: t00, eps_vap
     use therm_lib,   only: eslf, esif

   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                                             &
                                  j,ktf,kts,kte,ens4
     real, dimension (its:ite,kts:kte,1:ens4)                                 &
         , intent(inout) ::                                                   &
           tx,qx
     real, dimension (its:ite,kts:kte)                                 &
         , intent(in) ::                                                   &
           p
     real, dimension (its:ite)                                 &
         , intent(in) ::                                                   &
           z1,psur,cap_max,cap_max_increment
     real, intent(in) ::                                                   &
           tcrit,xl,rv,cp,mentr_rate,entr_rate
     real, dimension (its:ite,1:ens4)                                 &
         , intent(out) ::                                                   &
           axx
     integer, dimension (its:ite), intent (in) ::                             &
           ierr,kbmax
     integer, dimension (its:ite) ::                             &
           ierrxx,k22xx,kbconxx,ktopxx,kstabm,kstabi
!     real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: e,tvbar
     integer n,i,k,iph
     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,                                                  &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        tn_cup,                                                        &
        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,cd

     real,    dimension (its:ite) ::                                   &
       AA0,HKB,QKB,          &
       PWAV,BU
      do n=1,ens4
      do i=its,ite
       axx(i,n)=0.
      enddo
      enddo
!     HT(1)=XL/CP
!     HT(2)=2.834E6/CP
!     BE(1)=.622*HT(1)/.286
!     AE(1)=BE(1)/273.+ALOG(610.71)
!     BE(2)=.622*HT(2)/.286
!     AE(2)=BE(2)/273.+ALOG(610.71)
!
!
     do 100 n=1,ens4

      do k=kts,ktf
      do i=its,itf
        cd(i,k)=0.1*entr_rate
      enddo
      enddo


      do i=its,itf
        ierrxx(i)=ierr(i)
        k22xx(i)=1
        kbconxx(i)=1
        ktopxx(i)=1
        kstabm(i)=ktf-1
      enddo
      DO k=kts,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then

!        IPH=1
!        IF(Tx(I,K,n).LE.TCRIT)IPH=2
!        E=EXP(AE(IPH)-BE(IPH)/TX(I,K,N))
!        QES(I,K)=.622*E/(100.*P(I,K)-(1.-0.622)*E)

        if ( tx(i,k,n) > tcrit ) then
           e = eslf( tx(i,k,n) - t00 )
        else
           e = esif( tx(i,k,n) - t00 )
        endif

        qes(i,k) = max( eps_vap*e / (100.*p(i,k) - (1.0-eps_vap)*e), 1.e-8)

!       IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(Qx(I,K,N).GT.QES(I,K))Qx(I,K,N)=QES(I,K)
        TV(I,K)=Tx(I,K,N)+.608*Qx(I,K,N)*Tx(I,K,N)
        endif
      enddo
      enddo
!
         do i=its,itf
           if(ierrxx(i).eq.0)then
             Z(I,KTS) = Z1(I) + ALOG(PSUR(I)/P(I,KTS)) * 287.*TV(I,KTS)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierrxx(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K) = Z(I,K-1) + ALOG(P(I,K-1)/P(I,K)) * 287.*TVBAR/9.81
           endif
         enddo
         enddo
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierrxx(i).eq.0)then
         HE(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*Qx(I,K,n)
         HES(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

! cup levels
!
      do k=kts+1,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(qx(i,k-1,n)+qx(i,k,n))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(tx(i,k-1,n)+tx(i,k,n))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=qx(i,1,n)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=tx(i,1,n)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HE_CUP,3,KBMAX,K22XX,ierrxx, &
           ktf,kts,kte)
       DO 36 i=its,itf
         IF(ierrxx(I).eq.0)THEN
         IF(K22xx(I).GE.KBMAX(i))ierrxx(i)=2
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(cap_max_increment,1,k22xx,kbconxx,he_cup,hes_cup, &
           ierrxx,kbmax,p_cup,cap_max, &
           ktf,kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEs_cup,Kbconxx,kstabm,kstabi,ierrxx,  &
           ktf,kts,kte)
      do i=its,itf
      IF(ierrxx(I).eq.0)THEN
        if(kstabm(i)-1.gt.kstabi(i))then
           do k=kstabi(i),kstabm(i)-1
             cd(i,k)=cd(i,k-1)+1.5*entr_rate
             if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
           enddo
        ENDIF
      ENDIF
      ENDDO
!
!--- calculate incloud moist static energy
!
      call cup_up_he(k22xx,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbconxx,ierrxx,dby,he,hes_cup,'deep', &
           ktf,kts,kte)

!--- DETERMINE CLOUD TOP - KTOP
!
      call cup_ktop(1,dby,kbconxx,ktopxx,ierrxx, &
           ktf,kts,kte)
!
!c--- normalized updraft mass flux profile
!
      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbconxx,ktopxx,ierrxx,k22xx, &
           ktf,kts,kte)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbconxx,ktopxx,ierrxx,           &
           ktf,kts,kte)
      do i=its,itf
       if(ierrxx(i).eq.0)axx(i,n)=aa0(i)
      enddo
100   continue
     END SUBROUTINE cup_axx

END MODULE module_cu_g3
