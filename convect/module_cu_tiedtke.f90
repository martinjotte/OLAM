!####################TIEDTKE SCHEME#########################
!  Taken from the IPRC iRAM - Yuqing Wang, University of Hawaii
!  Added by Chunxi Zhang and Yuqing Wang to WRF3.2, May, 2010
!
! Tiedtke (1989, MWR, 117, 1779-1800)
! Nordeng, T.E., (1995), CAPE closure and organized entrainment/detrainment
! Yuqing Wang et al. (2003,J. Climate, 16, 1721-1738) improved cloud top detrainment
!                    (2004, Mon. Wea. Rev., 132, 274-296), improvements for PBL clouds
!                    (2007,Mon. Wea. Rev., 135, 567-585), diurnal pcp cycle
!###########################################################
MODULE module_cu_tiedtke

  real, parameter, private :: t000 = 273.15
  real, parameter, private :: hgfr = 233.15
  real, parameter, private :: ALV = 2.5008E6
  real, parameter, private :: ALS = 2.8345E6
  real, parameter, private :: ALF = ALS-ALV
  real, parameter, private :: CPD = 1005.46
  real, parameter, private :: CPV = 1869.46
  real, parameter, private :: RCPD = 1.0/CPD
  real, parameter, private :: RHOH2O = 1.0E03
  real, parameter, private :: TMELT = 273.16
  real, parameter, private :: G = 9.806
  real, parameter, private :: ZRG = 1.0/G
  real, parameter, private :: RD = 287.05
  real, parameter, private :: RV = 461.51
  real, parameter, private :: C1ES = 610.78
  real, parameter, private :: C2ES = C1ES*RD/RV
  real, parameter, private :: C3LES = 17.269
  real, parameter, private :: C4LES = 35.86
  real, parameter, private :: C5LES = C3LES*(TMELT-C4LES)
  real, parameter, private :: C3IES = 21.875
  real, parameter, private :: C4IES = 7.66
  real, parameter, private :: C5IES = C3IES*(TMELT-C4IES)
  real, parameter, private :: VTMPC1 = RV/RD-1.0
  real, parameter, private :: VTMPC2 = CPV/CPD-1.0
  real, parameter, private :: CVDIFTS = 1.0
  real, parameter, private :: CEVAPCU1 = 1.93E-6*261.0*0.5/G
  real, parameter, private :: CEVAPCU2 = 1.E3/(38.3*0.293)

! Specify TUNABLE parameters for massflux scheme

  real, parameter, private :: ENTRPEN = 1.0E-4 ! avg entrain rate for penetrative conv
  real, parameter, private :: ENTRSCV = 3.0E-4 ! avg entrain rate for shallow conv
  real, parameter, private :: ENTRMID = 1.0E-4 ! avg entrain rate for midlev conv
  real, parameter, private :: ENTRDD  = 2.0E-4 ! avg entrain rate for downdrafts
  real, parameter, private :: CMFCTOP = 0.30   ! relative cloud massflux at level
                                               ! above nonbuoyancy level
  real, parameter, private :: CMFCMAX = 1.0   ! max massflux allowed for updrafts, etc
  real, parameter, private :: CMFCMIN = 1.E-10 ! min massflux value (for safety)
  real, parameter, private :: CMFDEPS = 0.30 ! fractional massflux for downdrafts at LFS
  real, parameter, private :: CPRCON = 1.1E-3/G ! coeffs for conversion from cld water
  real, parameter, private :: ZDNOPRC = 1.5E4 ! pressure depth below which no precip

  integer, parameter, private :: orgen = 1  ! Old organized entrain rate
! integer, parameter, private :: orgen = 2  ! New organized entrain rate

  integer, parameter, private :: nturben = 1 ! old deep turbulent entrain/detrain rate
! integer, parameter, private :: nturben = 2 ! New deep turbulent entrain/detrain rate

  integer, parameter, private :: cutrigger = 1 ! Old trigger function
! integer, parameter, private :: cutrigger = 2 ! New trigger function

  real, parameter, private :: RHC = 0.80, RHM = 1.0, ZBUO0 = 0.50
  real, parameter, private :: CRIRH = 0.70, fdbk = 0.0, ZTAU = 2200.0

  logical, parameter, private :: LMFPEN = .TRUE.
  logical, parameter, private :: LMFMID = .TRUE.
  logical, parameter, private :: LMFSCV = .TRUE.
  logical, parameter, private :: LMFDD  = .TRUE.
  logical, parameter, private :: LMFDUDV= .TRUE.

CONTAINS

   subroutine cuparm_tiedtke(iw,km,km1,dtlong4,confrq4,confrq4i)

   use mem_grid,    only: mza, lpv, lpw, zt, xew, yew, zew, arv, arw, arw0, &
                          volti, vxn_ns, vyn_ns, vzn_ns, vxn_ew, vyn_ew
   use mem_cuparm,  only: thsrc, rtsrc, conprr, vxsrc, vysrc, vzsrc, rdsrc, &
                          kcubot, kcutop, cbmf, qwcon, iactcu, cddf, kddtop
   use mem_basic,   only: theta, tair, press, rho, vxe, vye, vze, rr_v, &
                          vmc, wmc
   use mem_turb,    only: frac_land, sfluxt, sfluxr, fqtpbl
   use consts_coms, only: eradi, gravo2
   use mem_ijtabs,  only: itab_w
   use oname_coms,  only: nl

   implicit none

   integer, intent(in) :: iw, km, km1
   real,    intent(in) :: dtlong4, confrq4, confrq4i

   real :: u1  (km) ! zonal wind component
   real :: v1  (km) ! meridional wind component
   real :: t1  (km) ! temperature
   real :: q1  (km) ! water vapor mixing ratio [kg/kg]

   real :: ght (km) ! geopotential height [m] at OLAM T levels
   real :: omg (km) ! 'omega' vertical velocity at OLAM T levels
   real :: prst(km) ! pressure at OLAM T levels
   real :: qsat(km) ! saturated vapor pressure
   real :: sig1(km) ! sigma_p value

   real :: dTdt(km) ! Temperature tendency
   real :: dQdt(km) ! Water Vapor tendency
   real :: dCdt(km) ! Cloud water/ice tendency
   real :: dUdt(km) ! east-west wind tendency
   real :: dVdt(km) ! nort-south wind tendency

   real :: prsw(km1) ! pressure at OLAM W levels

   real :: ztu (km) ! cloud temperature
   real :: zqu (km) ! cloud specific humidity
   real :: zlu (km) ! cloud liquid water
   real :: zlde(km) ! cloud water detrained to environment
   real :: zmfu(km) ! upraft mass flux
   real :: zmfd(km) ! downdraft mass flux

   real :: zdmfup(km)
   real :: zdmfdp(km)
   real :: zdpmel(km)

   real :: prsfc    ! surface rainfall rate [kg/(m^2 s)]
   real :: pssfc    ! surface snowfall rate [kg/(m^2 s)]

   real :: paprc, paprsm, paprs, zrain, psrain
   real :: psevap, psheat, psdiss, psmelt

   integer :: lndj  ! 0 over land, 1 over water
   integer :: ktype ! Convective closure type; output from Tiedtke param

   real :: evap  ! vapor flux at surface 
   real :: hfx   ! sensible heat flux at surface [W/m^2]
   real :: rhosf ! air density at lpw level

   integer :: ictop, icbot ! cloud top and bottom
   logical :: iact         ! was convection active in this cell

   real :: dQdt_sav(km)
   real :: vflux(mza), vflux_vap(mza)

   integer :: ka, k, kt, jv, iv, iwn, npoly

   real :: hflux, hflux_vap, dirv, flx, fqvadv
   real :: gnpoly1

   ka = lpw(iw)
   npoly = itab_w(iw)%npoly
   gnpoly1 = gravo2 / real(npoly+1)

! Vertical advective mass and water vapor fluxes (W levels)

   do k = ka,mza-1
      vflux(k) = arw(k,iw) * wmc(k,iw)

      ! upwinded
      if (wmc(k,iw) >= 0.0) then
         vflux_vap(k) = vflux(k) * rr_v(k,iw)
      else
         vflux_vap(k) = vflux(k) * rr_v(k+1,iw)
      endif
      
      ! centered
      ! vflux_vap(k) = vflux(k) * 0.5 * (rr_v(k,iw) + rr_v(k+1,iw))
   enddo
   vflux(ka-1) = 0.
   vflux(mza) = 0.
   vflux_vap(ka-1) = 0.
   vflux_vap(mza) = 0.

! Loop over T levels

   do k = ka,mza
      kt = mza + 1 - k        

      ! Compute zonal and meridional wind components

      u1(kt) = vxe(k,iw) * vxn_ew(iw) + vye(k,iw) * vyn_ew(iw)
      v1(kt) = vxe(k,iw) * vxn_ns(iw) + vye(k,iw) * vyn_ns(iw) + vze(k,iw) * vzn_ns(iw)

      ! Horizontal advective mass and water vapor fluxes

      hflux = 0.
      hflux_vap = 0.
      
      do jv = 1, npoly
         iv   = itab_w(iw)%iv(jv)

         if (k >= lpv(iv)) then
            dirv = itab_w(iw)%dirv(jv)
            iwn  = itab_w(iw)%iw(jv)

            flx   = dirv * vmc(k,iv) * arv(k,iv)
            hflux = hflux + flx

            ! upwinded
            if (flx >= 0.0) then
               hflux_vap = hflux_vap + flx * rr_v(k,iwn)
            else
               hflux_vap = hflux_vap + flx * rr_v(k,iw)
            endif

            ! centered
            ! hflux_vap = hflux_vap + flx * 0.5 * (rr_v(k,iw) + rr_v(k,iwn))
         endif
      enddo

      t1  (kt) = tair(k,iw)
      ght (kt) = g * zt(k)

      prst(kt) = press(k,iw)

      if (k == mza) then
         prsw(kt) = 2. * press(k,iw) - press(k-1,iw)
      else
         prsw(kt) = 0.5 * (press(k,iw) + press(k+1,iw))
      endif

      sig1(kt) = (press(k,iw)  - press(mza,iw)) &
               / (press(ka,iw) - press(mza,iw))

      ! use Tiedtke scheme's saturation calcs for now

      qsat(kt) = min(0.5, tlucua(t1(kt)) / prst(kt))
      qsat(kt) = qsat(kt) / (1.0 - vtmpc1*qsat(kt))
      q1  (kt) = min(rr_v(k,iw), qsat(kt))

      ! Average vertical velocity from current and surrounding cells

      omg (kt) = -gnpoly1 * (wmc(k-1,iw) + wmc(k,iw) + sum(wmc(k-1:k,itab_w(iw)%iw(1:npoly))))

      ! Initialize tendencies to zero

      dTdt(kt) = 0.0
      dCdt(kt) = 0.0
      dUdt(kt) = 0.0
      dVdt(kt) = 0.0

      ! Tiedtke scheme requires large scale moisture tendency

      fqvadv = ((vflux_vap(k-1) - vflux_vap(k) + hflux_vap) &
               - (vflux(k-1) - vflux(k) + hflux) * rr_v(k,iw)) &
               * volti(k,iw) / real(rho(k,iw))

      dQdt    (kt) = fqtpbl(k,iw) + fqvadv
      dQdt_sav(kt) = dQdt(kt)

   enddo

   ! prsw(1) is press at zm(mza); prsw(km1) is press at zm(ka-1)

   prsw(1) = max(1.e-3,real(1.5 * press(mza,iw) - 0.5 * press(mza-1,iw)))
   prsw(km1) = 1.5 * press(ka,iw) - 0.5 * press(ka+1,iw)

   evap  = sfluxr(iw)
   hfx   = sfluxt(iw) * tair(ka,iw) / theta(ka,iw)
   rhosf = rho(ka,iw)

   if (allocated(frac_land)) then
      lndj = nint(1.0 - frac_land(iw))
   else
      lndj = 1
   endif

   iact   = .false.
   prsfc  = 0.0
   pssfc  = 0.0
   paprc  = 0.0
   paprsm = 0.0
   paprs  = 0.0
   zrain  = 0.0
   psrain = 0.0
   psevap = 0.0
   psheat = 0.0
   psdiss = 0.0
   psmelt = 0.0

   ! Tiedtke convective parameterization for one IW column

   CALL CUMASTR_NEW(                                   &
        km,       km1,      km-1,     t1,              &
        q1,       u1,       v1,       omg,     qsat,   &
        evap,     confrq4,  prst,     prsw,    ght,    &
        dTdt,     dQdt,     dUdt,     dVdt,    prsfc,  & 
        pssfc,    paprc,    paprsm,   paprs,   iact,   &
        ktype,    icbot,    ictop,    ztu,     zqu,    &
        zlu,      zlde,     zmfu,     zmfd,    zrain,  &
        psrain,   psevap,   psheat,   psdiss,  psmelt, &
        dCdt,     hfx,      rhosf,    sig1,    lndj,   &
        zdmfup,   zdmfdp,   zdpmel                     )

! Apply convective parameterization results to main model arrays

   if (iact .and. ictop < icbot .and. zmfu(icbot) > 1.1e-10) then

      ! convection extends 1 level above ictop
      ictop = ictop - 1

      kcutop(iw) = mza - ictop + 1
      kcubot(iw) = mza - icbot + 1
      iactcu(iw) = 1
      cbmf  (iw) = zmfu(icbot)

      ! precipitation rate
      conprr(iw) = max(prsfc+pssfc, 0.0)

      do k = kcutop(iw), kcubot(iw), -1
         kt = mza + 1 - k
         if (zmfd(kt) < -1.e-10) then
            kddtop(iw) = k
            cddf  (iw) = -zmfd(kt)
            exit
         endif
      enddo

      do k = ka, kcutop(iw)
         kt = mza + 1 - k

         ! subtract off the input humdity tendency
         dQdt(kt) = dQdt(kt) - dQdt_sav(kt)

         ! If convection created any ice or liquid, add it to the total
         ! water and evaporate it. Shouldn't be needed with fdbk=0 though.
         !if (dCdt(kt) > 1.e-18) then
         !   dQdt(kt) = dQdt(kt) + dCdt(kt)
         !   
         !   if (t1(kt) > tmelt) then
         !      dTdt(kt) = dTdt(kt) - alv * rcpd * dCdt(kt)
         !   else
         !      dTdt(kt) = dTdt(kt) - als * rcpd * dCdt(kt)
         !   endif
         !endif

         ! store convective heating and moisture rates and cloud water

         thsrc(k,iw) = dTdt(kt) * real(rho(k,iw))
         rtsrc(k,iw) = dQdt(kt) * real(rho(k,iw))

         qwcon(k,iw) = zlu(kt)

         ! density tendency (water removed per grid cell)
         rdsrc(k,iw) = -(zdmfup(kt) + zdmfdp(kt)) * arw0(iw) * volti(k,iw)

      enddo

      ! convective momentum transport

      if (nl%conv_uv_mix > 0) then
         do k = ka, kcutop(iw)
            kt = mza + 1 - k

            vxsrc(k,iw) = real(rho(k,iw)) * (dUdt(kt) * vxn_ew(iw) + dVdt(kt) * vxn_ns(iw))
            vysrc(k,iw) = real(rho(k,iw)) * (dUdt(kt) * vyn_ew(iw) + dVdt(kt) * vyn_ns(iw))
            vzsrc(k,iw) = real(rho(k,iw)) * (                        dVdt(kt) * vzn_ns(iw))

         enddo
      endif

   endif

   end subroutine cuparm_tiedtke


!------------This is the combined version for tiedtke---------------
!  In this module only the mass flux convection scheme of the ECMWF is included
!
!***********************************************************
      SUBROUTINE CUMASTR_NEW                             &
         (KLEV,     KLEVP1,   KLEVM1,   PTEN,            &
          PQEN,     PUEN,     PVEN,     PVERV,    PQSEN, &
          PQHFL,    ZTMST,    PAP,      PAPH,     PGEO,  &
          PTTE,     PQTE,     PVOM,     PVOL,     PRSFC, &
          PSSFC,    PAPRC,    PAPRSM,   PAPRS,    LDCUM, &
          KTYPE,    KCBOT,    KCTOP,    PTU,      PQU,   &
          PLU,      PLUDE,    PMFU,     PMFD,     PRAIN, &
          PSRAIN,   PSEVAP,   PSHEAT,   PSDISS,   PSMELT,& 
          PCTE,     PHHFL,    RHO,      sig1,     lndj,  &
          zdmfup,   zdmfdp,   zdpmel                     )
!
!***CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!***PURPOSE
!   -------
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!***INTERFACE.
!   ----------
!          *CUMASTR* IS CALLED FROM *MSSFLX*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
!        (USED IN SUBROUTINE CLOUD)
!***METHOD
!   ------
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!***EXTERNALS.
!   ----------
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!***SWITCHES.
!   --------
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!***
!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY
!                LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!***REFERENCE.
!   ----------
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!-----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KLEVM1
      REAL      ZTMST
      REAL      PSRAIN, PSEVAP, PSHEAT, PSDISS, PSMELT, ZCONS2
      INTEGER   JK,IKB
      REAL      ZQUMQE, ZDQMIN, ZMFMAX, ZALVDCP, ZQALV
      REAL      ZHSAT, ZGAM, ZZZ, ZHHAT, ZBI, ZRO, ZDZ, ZDHDZ, ZDEPTH
      REAL      ZFAC, ZRH, ZPBMPT, DEPT, ZHT, ZEPS
      INTEGER   ICUM, ITOPM2
      REAL     PTEN(KLEV),        PQEN(KLEV), &
              PUEN(KLEV),        PVEN(KLEV),  &
              PTTE(KLEV),        PQTE(KLEV),  &
              PVOM(KLEV),        PVOL(KLEV),  &
              PQSEN(KLEV),       PGEO(KLEV),  &
              PAP(KLEV),         PAPH(KLEVP1),& 
              PVERV(KLEV),       PQHFL,       &
              PHHFL,             RHO
      REAL     PTU(KLEV),        PQU(KLEV),   &
              PLU(KLEV),         PLUDE(KLEV), &
              PMFU(KLEV),        PMFD(KLEV),  &
              PAPRC,            PAPRS,      &
              PAPRSM,           PRAIN,      &
              PRSFC,            PSSFC
      REAL     ZTENH(KLEV),       ZQENH(KLEV),&
              ZGEOH(KLEV),       ZQSENH(KLEV),&
              ZTD(KLEV),         ZQD(KLEV),   &
              ZMFUS(KLEV),       ZMFDS(KLEV), &
              ZMFUQ(KLEV),       ZMFDQ(KLEV), &
              ZDMFUP(KLEV),      ZDMFDP(KLEV),& 
              ZMFUL(KLEV),       ZRFL,       &
              ZUU(KLEV),         ZVU(KLEV),   &
              ZUD(KLEV),         ZVD(KLEV)
      REAL     ZENTR,            ZHCBASE,   &
              ZMFUB,            ZMFUB1,     &
              ZDQPBL,           ZDQCV 
      REAL     ZSFL,             ZDPMEL(KLEV), &
              PCTE(KLEV),        ZCAPE,        &
              ZHEAT,            ZHHATT(KLEV),  &
              ZHMIN,            ZRELH
      REAL     sig1(KLEV)
      INTEGER  ILAB(KLEV),        IDTOP,   &
              ICTOP0,           ILWMIN    
      INTEGER  KCBOT,            KCTOP,   &
              KTYPE,            IHMIN,    &
              KTOP0,                  lndj
      LOGICAL  LDCUM
      LOGICAL  LODDRAF,          LLO1
      REAL     CRIRH1
!-------------------------------------------
!     1.    SPECIFY CONSTANTS AND PARAMETERS
!-------------------------------------------
  100 CONTINUE
      ZCONS2=1./(G*ZTMST)
!--------------------------------------------------------------
!*    2.    INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!--------------------------------------------------------------
  200 CONTINUE
      CALL CUINI &
         (KLEV,     KLEVP1,   KLEVM1,   PTEN,            &
          PQEN,     PQSEN,    PUEN,     PVEN,     PVERV, &
          PGEO,     PAPH,     ZGEOH,    ZTENH,    ZQENH,  &
          ZQSENH,   ILWMIN,   PTU,      PQU,      ZTD,   &
          ZQD,      ZUU,      ZVU,      ZUD,      ZVD,   &
          PMFU,     PMFD,     ZMFUS,    ZMFDS,    ZMFUQ, &
          ZMFDQ,    ZDMFUP,   ZDMFDP,   ZDPMEL,   PLU,  &
          PLUDE,    ILAB)
!----------------------------------
!*    3.0   CLOUD BASE CALCULATIONS
!----------------------------------
  300 CONTINUE
!*         (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!          -------------------------------------------
      CALL CUBASE &
         (KLEV,     KLEVP1,   KLEVM1,   ZTENH,           &
          ZQENH,    ZGEOH,    PAPH,     PTU,      PQU,   &
          PLU,      PUEN,     PVEN,     ZUU,      ZVU,   &
          LDCUM,    KCBOT,    ILAB)
!*          (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*              THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!               -----------------------------------------
       JK=1
       ZDQCV =PQTE(JK)*(PAPH(JK+1)-PAPH(JK))
       ZDQPBL=0.0
       IDTOP=0
       DO 320 JK=2,KLEV
       ZDQCV=ZDQCV+PQTE(JK)*(PAPH(JK+1)-PAPH(JK))
       IF(JK.GE.KCBOT) ZDQPBL=ZDQPBL+PQTE(JK)  &
                                    *(PAPH(JK+1)-PAPH(JK))
  320 CONTINUE

      if(cutrigger .eq. 1) then
         KTYPE=0
        IF(ZDQCV.GT.MAX(0.,1.1*PQHFL*G)) THEN
         KTYPE=1
        ELSE
         KTYPE=2
        ENDIF
      else if(cutrigger .eq. 2) then
         CALL CUTYPE  &
          (KLEV,     KLEVP1,   KLEVM1,                       &
          ZTENH,   ZQENH,       ZQSENH,    ZGEOH,     PAPH,  &
          RHO,     PHHFL,         PQHFL,    KTYPE,    lndj   )
      end if
!*         (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*             AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*             THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!              ------------------------------------------
!        if(ktype .ge. 1 ) then
!              write(6,*)"ktype=", KTYPE
!        end if

      IKB=KCBOT
      ZQUMQE=PQU(IKB)+PLU(IKB)-ZQENH(IKB)
      ZDQMIN=MAX(0.01*ZQENH(IKB),1.E-10)
      IF(ZDQPBL.GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM) THEN
         ZMFUB=ZDQPBL/(G*MAX(ZQUMQE,ZDQMIN))
      ELSE
         ZMFUB=0.01
         LDCUM=.FALSE.
      ENDIF
      ZMFMAX=(PAPH(IKB)-PAPH(IKB-1))*ZCONS2
      ZMFUB=MIN(ZMFUB,ZMFMAX)
!------------------------------------------------------
!*    4.0   DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!------------------------------------------------------
  400 CONTINUE
!*         (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*             CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*             FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
! -------------------------------------------------------------
      IKB=KCBOT
      ZHCBASE=CPD*PTU(IKB)+ZGEOH(IKB)+ALV*PQU(IKB)
      ICTOP0=KCBOT-1
      ZALVDCP=ALV/CPD
      ZQALV=1./ALV
      DO 420 JK=KLEVM1,3,-1
      ZHSAT=CPD*ZTENH(JK)+ZGEOH(JK)+ALV*ZQSENH(JK)
      ZGAM=C5LES*ZALVDCP*ZQSENH(JK)/  &
          ((1.-VTMPC1*ZQSENH(JK))*(ZTENH(JK)-C4LES)**2)
      ZZZ=CPD*ZTENH(JK)*0.608
      ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)* &
                 MAX(ZQSENH(JK)-ZQENH(JK),0.)
      ZHHATT(JK)=ZHHAT
      IF(JK.LT.ICTOP0.AND.ZHCBASE.GT.ZHHAT) ICTOP0=JK
  420 CONTINUE
      JK=KCBOT
      ZHSAT=CPD*ZTENH(JK)+ZGEOH(JK)+ALV*ZQSENH(JK)
      ZGAM=C5LES*ZALVDCP*ZQSENH(JK)/   &
          ((1.-VTMPC1*ZQSENH(JK))*(ZTENH(JK)-C4LES)**2)
      ZZZ=CPD*ZTENH(JK)*0.608
      ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)* &
                 MAX(ZQSENH(JK)-ZQENH(JK),0.)
      ZHHATT(JK)=ZHHAT
!
! Find lowest possible org. detrainment level
!
         ZHMIN = 0.
         IF( LDCUM.AND.KTYPE.EQ.1 ) THEN
            IHMIN = KCBOT
         ELSE
            IHMIN = -1
         END IF
!
      ZBI = 1./(25.*G)
      DO 450 JK = KLEV, 1, -1
      LLO1 = LDCUM.AND.KTYPE.EQ.1.AND.IHMIN.EQ.KCBOT
      IF (LLO1.AND.JK.LT.KCBOT.AND.JK.GE.ICTOP0) THEN
        IKB = KCBOT
        ZRO = RD*ZTENH(JK)/(G*PAPH(JK))
        ZDZ = (PAPH(JK)-PAPH(JK-1))*ZRO
        ZDHDZ=(CPD*(PTEN(JK-1)-PTEN(JK))+ALV*(PQEN(JK-1)-   &
          PQEN(JK))+(PGEO(JK-1)-PGEO(JK)))*G/(PGEO(JK-1)-PGEO(JK))
        ZDEPTH = ZGEOH(JK) - ZGEOH(IKB)
        ZFAC = SQRT(1.+ZDEPTH*ZBI)
        ZHMIN = ZHMIN + ZDHDZ*ZFAC*ZDZ
        ZRH = -ALV*(ZQSENH(JK)-ZQENH(JK))*ZFAC
        IF (ZHMIN.GT.ZRH) IHMIN = JK
      END IF
 450  CONTINUE 
      IF (LDCUM.AND.KTYPE.EQ.1) THEN
        IF (IHMIN.LT.ICTOP0) IHMIN = ICTOP0
      END IF
      IF(KTYPE.EQ.1) THEN
        ZENTR=ENTRPEN
      ELSE
        ZENTR=ENTRSCV
      ENDIF
      if(lndj.eq.1) ZENTR=ZENTR*1.05
!*         (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!----------------------------------------------------------
      CALL CUASC_NEW &
         (KLEV,     KLEVP1,   KLEVM1,   ZTENH,             &
          ZQENH,    PUEN,     PVEN,     PTEN,     PQEN,    &
          PQSEN,    PGEO,     ZGEOH,    PAP,      PAPH,    &
          PQTE,     PVERV,    ILWMIN,   LDCUM,    ZHCBASE, &
          KTYPE,    ILAB,     PTU,      PQU,      PLU,     &
          ZUU,      ZVU,      PMFU,     ZMFUB,    ZENTR,   &
          ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP,  &
          KCBOT,    KCTOP,    ICTOP0,   ICUM,     ZTMST,   &
          IHMIN,    ZHHATT,   ZQSENH)
      IF(ICUM.EQ.0) GO TO 1000
!*     (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!          CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!------------------------------------------------------------------
      ZPBMPT=PAPH(KCBOT)-PAPH(KCTOP)
      IF(LDCUM) ICTOP0=KCTOP
      IF(LDCUM.AND.KTYPE.EQ.1.AND.ZPBMPT.LT.ZDNOPRC) KTYPE=2
      IF(KTYPE.EQ.2) then
        ZENTR=ENTRSCV
        if(lndj.eq.1) ZENTR=ZENTR*1.05
      endif
      ZRFL=ZDMFUP(1)
      DO 490 JK=2,KLEV
          ZRFL=ZRFL+ZDMFUP(JK)
  490 CONTINUE
!-----------------------------------------
!*    5.0   CUMULUS DOWNDRAFT CALCULATIONS
!-----------------------------------------
  500 CONTINUE
      IF(LMFDD) THEN
!*      (A) DETERMINE LFS IN 'CUDLFS'
!--------------------------------------
         CALL CUDLFS &
         (KLEV,     KLEVP1,   ZTENH,    ZQENH,            &
          PUEN,     PVEN,     ZGEOH,    PAPH,     PTU,    &
          PQU,      ZUU,      ZVU,      LDCUM,    KCBOT,  &
          KCTOP,    ZMFUB,    ZRFL,     ZTD,      ZQD,    &
          ZUD,      ZVD,      PMFD,     ZMFDS,    ZMFDQ,  &
          ZDMFDP,   IDTOP,    LODDRAF)
!*     (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!------------------------------------------------------------
         CALL CUDDRAF &
         (KLEV,     KLEVP1,   ZTENH,    ZQENH,            &
          PUEN,     PVEN,     ZGEOH,    PAPH,     ZRFL,   &
          LODDRAF,  ZTD,      ZQD,      ZUD,      ZVD,    &
          PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP)
!*     (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!           DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!-----------------------------------------------------------
      END IF
!
!-- 5.1 Recalculate cloud base massflux from a cape closure
!       for deep convection (ktype=1) and by PBL equilibrium
!       taking downdrafts into account for shallow convection
!       (ktype=2)
!       implemented by Y. WANG based on ECHAM4 in Nov. 2001.
!
        ZHEAT=0.0
        ZCAPE=0.0
        ZRELH=0.0
        ZMFUB1=ZMFUB
!
      IF(LDCUM.AND.KTYPE.EQ.1) THEN
      do jk=KLEVM1,2,-1
      if(abs(paph(jk)*0.01 - 300) .lt. 50.) then
        KTOP0=MAX(jk,KCTOP)
        exit
      end if
      end do
!      KTOP0=MAX(12,KCTOP)
       DO JK=2,KLEV
       IF(JK.LE.KCBOT.AND.JK.GT.KCTOP) THEN
         ZRO=PAPH(JK)/(RD*ZTENH(JK))
         ZDZ=(PAPH(JK)-PAPH(JK-1))/(G*ZRO)
         ZHEAT=ZHEAT+((PTEN(JK-1)-PTEN(JK)   &
           +G*ZDZ/CPD)/ZTENH(JK)+0.608*(PQEN(JK-1)-  &
           PQEN(JK)))*(PMFU(JK)+PMFD(JK))*G/ZRO
         ZCAPE=ZCAPE+G*((PTU(JK)*(1.+.608*PQU(JK) &
           -PLU(JK)))/(ZTENH(JK)*(1.+.608*ZQENH(JK))) &
           -1.0)*ZDZ
       ENDIF
       IF(JK.LE.KCBOT.AND.JK.GT.KTOP0) THEN
         dept=(PAPH(JK)-PAPH(JK-1))/(PAPH(KCBOT)-  &
            PAPH(KTOP0))
         ZRELH=ZRELH+dept*PQEN(JK)/PQSEN(JK)
       ENDIF
       ENDDO
!
       
       if(cutrigger .eq. 1 ) then 
         IF(lndj.EQ.1) then
           CRIRH1=CRIRH*0.8
         ELSE
           CRIRH1=CRIRH
         ENDIF
       else
          CRIRH1=0.
       end if

       IF(ZRELH.GE.CRIRH1 .AND. ZCAPE .GT. 100.) THEN
         IKB=KCBOT
         ZHT=ZCAPE/(ZTAU*ZHEAT)
         ZMFUB1=MAX(ZMFUB*ZHT,0.01)
         ZMFMAX=(PAPH(IKB)-PAPH(IKB-1))*ZCONS2
         ZMFUB1=MIN(ZMFUB1,ZMFMAX)
       ELSE
         ZMFUB1=0.01
         ZMFUB=0.01
         LDCUM=.FALSE.
        ENDIF
       ENDIF
!
!*  5.2   RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!         DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!--------------------------------------------------------
        IF(KTYPE.NE.1) THEN
           IKB=KCBOT
           IF(PMFD(IKB).LT.0.0.AND.LODDRAF) THEN
              ZEPS=CMFDEPS
           ELSE
              ZEPS=0.
           ENDIF
           ZQUMQE=PQU(IKB)+PLU(IKB)-          &
                 ZEPS*ZQD(IKB)-(1.-ZEPS)*ZQENH(IKB)
           ZDQMIN=MAX(0.01*ZQENH(IKB),1.E-10)
           ZMFMAX=(PAPH(IKB)-PAPH(IKB-1))*ZCONS2
           IF(ZDQPBL.GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM &
             .AND.ZMFUB.LT.ZMFMAX) THEN
              ZMFUB1=ZDQPBL/(G*MAX(ZQUMQE,ZDQMIN))
           ELSE
              ZMFUB1=ZMFUB
           ENDIF
           LLO1=(KTYPE.EQ.2).AND.ABS(ZMFUB1  &
                -ZMFUB).LT.0.2*ZMFUB
           IF(.NOT.LLO1) ZMFUB1=ZMFUB
           ZMFUB1=MIN(ZMFUB1,ZMFMAX)
        END IF
        DO 530 JK=1,KLEV
        IF(LDCUM) THEN
           ZFAC=ZMFUB1/MAX(ZMFUB,1.E-10)
           PMFD(JK)=PMFD(JK)*ZFAC
           ZMFDS(JK)=ZMFDS(JK)*ZFAC
           ZMFDQ(JK)=ZMFDQ(JK)*ZFAC
           ZDMFDP(JK)=ZDMFDP(JK)*ZFAC
        ELSE
           PMFD(JK)=0.0
           ZMFDS(JK)=0.0
           ZMFDQ(JK)=0.0
           ZDMFDP(JK)=0.0
        ENDIF
  530   CONTINUE
           IF(LDCUM) THEN
              ZMFUB=ZMFUB1
           ELSE
              ZMFUB=0.0
           ENDIF
!
!---------------------------------------------------------------
!*    6.0      DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*             FOR PENETRATIVE CONVECTION (TYPE=1),
!*             FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*             AND FOR MID-LEVEL CONVECTION (TYPE=3).
!---------------------------------------------------------------
  600 CONTINUE
      CALL CUASC_NEW &
         (KLEV,     KLEVP1,   KLEVM1,   ZTENH,            &
          ZQENH,    PUEN,     PVEN,     PTEN,     PQEN,   &
          PQSEN,    PGEO,     ZGEOH,    PAP,      PAPH,   &
          PQTE,     PVERV,    ILWMIN,   LDCUM,    ZHCBASE,& 
          KTYPE,    ILAB,     PTU,      PQU,      PLU,    &
          ZUU,      ZVU,      PMFU,     ZMFUB,    ZENTR,  &
          ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP, &
          KCBOT,    KCTOP,    ICTOP0,   ICUM,     ZTMST,  &
          IHMIN,    ZHHATT,   ZQSENH)
!----------------------------------------------------------
!*    7.0      DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!----------------------------------------------------------
  700 CONTINUE
      CALL CUFLX &
         (KLEV,     KLEVP1,   PQEN,     PQSEN,            &
          ZTENH,    ZQENH,    PAPH,     ZGEOH,    KCBOT,  &
          KCTOP,    IDTOP,    KTYPE,    LODDRAF,  LDCUM,  &
          PMFU,     PMFD,     ZMFUS,    ZMFDS,    ZMFUQ,  &
          ZMFDQ,    ZMFUL,    PLUDE,    ZDMFUP,   ZDMFDP, &
          ZRFL,     PRAIN,    PTEN,     ZSFL,     ZDPMEL, &
          ITOPM2,   ZTMST,    sig1)
!----------------------------------------------------------------
!*    8.0      UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!----------------------------------------------------------------
  800 CONTINUE
      CALL CUDTDQ                                          &
         (KLEV,     KLEVP1,   ITOPM2,   PAPH,              &
          LDCUM,    PTEN,     PTTE,     PQTE,     ZMFUS,   &
          ZMFDS,    ZMFUQ,    ZMFDQ,    ZMFUL,    ZDMFUP,  &
          ZDMFDP,   ZTMST,    ZDPMEL,   PRAIN,    ZRFL,    &
          ZSFL,     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT,  &
          PRSFC,    PSSFC,    PAPRC,    PAPRSM,   PAPRS,   &
          PQEN,     PQSEN,    PLUDE,    PCTE)
!----------------------------------------------------------------
!*    9.0      UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!----------------------------------------------------------------
  900 CONTINUE
      IF(LMFDUDV) THEN
      CALL CUDUDV  &
         (KLEV,     KLEVP1,   ITOPM2,   KTYPE,             &
          KCBOT,    PAPH,     LDCUM,    PUEN,     PVEN,    &
          PVOM,     PVOL,     ZUU,      ZUD,      ZVU,     &
          ZVD,      PMFU,     PMFD,     PSDISS)
      END IF
 1000 CONTINUE
      RETURN
      END SUBROUTINE CUMASTR_NEW
!

!#############################################################
!
!             LEVEL 3 SUBROUTINEs
!
!#############################################################
!**********************************************
!       SUBROUTINE CUINI
!**********************************************
!
      SUBROUTINE CUINI                                    &
         (KLEV,     KLEVP1,   KLEVM1,   PTEN,             &
          PQEN,     PQSEN,    PUEN,     PVEN,     PVERV,  &
          PGEO,     PAPH,     PGEOH,    PTENH,    PQENH,  &
          PQSENH,   KLWMIN,   PTU,      PQU,      PTD,    &
          PQD,      PUU,      PVU,      PUD,      PVD,    &
          PMFU,     PMFD,     PMFUS,    PMFDS,    PMFUQ,  &
          PMFDQ,    PDMFUP,   PDMFDP,   PDPMEL,   PLU,    &
          PLUDE,    KLAB)
!      M.TIEDTKE         E.C.M.W.F.     12/89
!***PURPOSE
!   -------
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!***METHOD.
!  --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!***EXTERNALS
!   ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   klevm1
      INTEGER   JK,IK, ICALL
      REAL      ZDP, ZZS
      REAL     PTEN(KLEV),        PQEN(KLEV),    &
              PUEN(KLEV),        PVEN(KLEV),     &
              PQSEN(KLEV),       PVERV(KLEV),    &
              PGEO(KLEV),        PGEOH(KLEV),    &
              PAPH(KLEVP1),      PTENH(KLEV),    &
              PQENH(KLEV),       PQSENH(KLEV)
      REAL     PTU(KLEV),         PQU(KLEV),     &
              PTD(KLEV),         PQD(KLEV),      &
              PUU(KLEV),         PUD(KLEV),      &
              PVU(KLEV),         PVD(KLEV),      &
              PMFU(KLEV),        PMFD(KLEV),     &
              PMFUS(KLEV),       PMFDS(KLEV),    &
              PMFUQ(KLEV),       PMFDQ(KLEV),    &
              PDMFUP(KLEV),      PDMFDP(KLEV),   & 
              PLU(KLEV),         PLUDE(KLEV)
      REAL     ZWMAX,            ZPH,          &
              PDPMEL(KLEV)
      INTEGER  KLAB(KLEV),        KLWMIN
      LOGICAL  LOFLAG
!------------------------------------------------------------
!*    1.       SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*             ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*             FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
! -----------------------------------------------------------
  100 CONTINUE
      ZDP=0.5
      DO 130 JK=2,KLEV
      PGEOH(JK)=PGEO(JK)+(PGEO(JK-1)-PGEO(JK))*ZDP
      PTENH(JK)=(MAX(CPD*PTEN(JK-1)+PGEO(JK-1),   &
                  CPD*PTEN(JK)+PGEO(JK))-PGEOH(JK))*RCPD
      PQSENH(JK)=PQSEN(JK-1)
      ZPH=PAPH(JK)
      LOFLAG=.TRUE.
      IK=JK
      ICALL=0
      CALL CUADJTQ(KLEV,IK,ZPH,PTENH,PQSENH,LOFLAG,ICALL)
      PQENH(JK)=MIN(PQEN(JK-1),PQSEN(JK-1))    &
                 +(PQSENH(JK)-PQSEN(JK-1))
      PQENH(JK)=MAX(PQENH(JK),0.)
  130 CONTINUE
      PTENH(KLEV)=(CPD*PTEN(KLEV)+PGEO(KLEV)-   &
                     PGEOH(KLEV))*RCPD
      PQENH(KLEV)=PQEN(KLEV)
      PTENH(1)=PTEN(1)
      PQENH(1)=PQEN(1)
      PGEOH(1)=PGEO(1)
      KLWMIN=KLEV
      ZWMAX=0.
      DO 160 JK=KLEVM1,2,-1
      ZZS=MAX(CPD*PTENH(JK)+PGEOH(JK),   &
             CPD*PTENH(JK+1)+PGEOH(JK+1))
      PTENH(JK)=(ZZS-PGEOH(JK))*RCPD
  160 CONTINUE
      DO 190 JK=KLEV,3,-1
      IF(PVERV(JK).LT.ZWMAX) THEN
         ZWMAX=PVERV(JK)
         KLWMIN=JK
      END IF
  190 CONTINUE
!-----------------------------------------------------------
!*    2.0      INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!-----------------------------------------------------------
  200 CONTINUE
      DO 230 JK=1,KLEV
      IK=JK-1
      IF(JK.EQ.1) IK=1
      PTU(JK)=PTENH(JK)
      PTD(JK)=PTENH(JK)
      PQU(JK)=PQENH(JK)
      PQD(JK)=PQENH(JK)
      PLU(JK)=0.
      PUU(JK)=PUEN(IK)
      PUD(JK)=PUEN(IK)
      PVU(JK)=PVEN(IK)
      PVD(JK)=PVEN(IK)
      PMFU(JK)=0.
      PMFD(JK)=0.
      PMFUS(JK)=0.
      PMFDS(JK)=0.
      PMFUQ(JK)=0.
      PMFDQ(JK)=0.
      PDMFUP(JK)=0.
      PDMFDP(JK)=0.
      PDPMEL(JK)=0.
      PLUDE(JK)=0.
      KLAB(JK)=0
  230 CONTINUE
      RETURN
      END SUBROUTINE CUINI   

!**********************************************
!       SUBROUTINE CUBASE
!********************************************** 
      SUBROUTINE CUBASE &
         (KLEV,     KLEVP1,   KLEVM1,   PTENH,           &
          PQENH,    PGEOH,    PAPH,     PTU,      PQU,   &
          PLU,      PUEN,     PVEN,     PUU,      PVU,   &
          LDCUM,    KCBOT,    KLAB)
!      THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!      FOR CUMULUS PARAMETERIZATION
!      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF.  12/89
!***PURPOSE.
!   --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!***METHOD.
!  --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!***EXTERNALS
!   ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   klevm1
      INTEGER   JK,IS,IK,ICALL,IKB
      REAL      ZBUO,ZZ
      REAL     PTENH(KLEV),       PQENH(KLEV),  &
              PGEOH(KLEV),       PAPH(KLEVP1)
      REAL     PTU(KLEV),         PQU(KLEV),   &
              PLU(KLEV)
      REAL     PUEN(KLEV),        PVEN(KLEV),  &
              PUU(KLEV),         PVU(KLEV) 
      REAL     ZQOLD(KLEV),       ZPH
      INTEGER  KLAB(KLEV),        KCBOT
      LOGICAL  LDCUM,            LOFLAG
!***INPUT VARIABLES:
!       PTENH [ZTENH] - Environment Temperature on half levels. (CUINI)
!       PQENH [ZQENH] - Env. specific humidity on half levels. (CUINI)
!       PGEOH [ZGEOH] - Geopotential on half levels, (MSSFLX)
!       PAPH - Pressure of half levels. (MSSFLX)
!***VARIABLES MODIFIED BY CUBASE:
!       LDCUM - Logical denoting profiles. (CUBASE)
!       KTYPE - Convection type - 1: Penetrative  (CUMASTR)
!                                 2: Stratocumulus (CUMASTR)
!                                 3: Mid-level  (CUASC)
!       PTU - Cloud Temperature.
!       PQU - Cloud specific Humidity.
!       PLU - Cloud Liquid Water (Moisture condensed out)
!       KCBOT - Cloud Base Level. (CUBASE)
!       KLAB [ILAB] - Level Label - 1: Sub-cloud layer (CUBASE)
!------------------------------------------------
!     1.       INITIALIZE VALUES AT LIFTING LEVEL
!------------------------------------------------
  100 CONTINUE
        KLAB(KLEV)=1
        KCBOT=KLEVM1
        LDCUM=.FALSE.
        PUU(KLEV)=PUEN(KLEV)*(PAPH(KLEVP1)-PAPH(KLEV))
        PVU(KLEV)=PVEN(KLEV)*(PAPH(KLEVP1)-PAPH(KLEV))
!-------------------------------------------------------
!     2.0      DO ASCENT IN SUBCLOUD LAYER,
!              CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!              ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!              CHECK FOR BUOYANCY AND SET FLAGS
!-------------------------------------------------------
      DO 200 JK=1,KLEV
        ZQOLD(JK)=0.0
  200 CONTINUE
      DO 290 JK=KLEVM1,2,-1
        IS=0
          IF(KLAB(JK+1).EQ.1) THEN
             IS=IS+1
             LOFLAG=.TRUE.
          ELSE
             LOFLAG=.FALSE.
          ENDIF
          ZPH=PAPH(JK)
        IF(IS.EQ.0) GO TO 290
          IF(LOFLAG) THEN
             PQU(JK)=PQU(JK+1)
             PTU(JK)=(CPD*PTU(JK+1)+PGEOH(JK+1)  &
                       -PGEOH(JK))*RCPD
             ZBUO=PTU(JK)*(1.+VTMPC1*PQU(JK))-      &
                 PTENH(JK)*(1.+VTMPC1*PQENH(JK))+ZBUO0
             IF(ZBUO.GT.0.) KLAB(JK)=1
             ZQOLD(JK)=PQU(JK)
          END IF
        IK=JK
        ICALL=1
        CALL CUADJTQ(KLEV,IK,ZPH,PTU,PQU,LOFLAG,ICALL)
          IF(LOFLAG.AND.PQU(JK).NE.ZQOLD(JK)) THEN
             KLAB(JK)=2
             PLU(JK)=PLU(JK)+ZQOLD(JK)-PQU(JK)
             ZBUO=PTU(JK)*(1.+VTMPC1*PQU(JK))-      &
                 PTENH(JK)*(1.+VTMPC1*PQENH(JK))+ZBUO0
             IF(ZBUO.GT.0.) THEN
                KCBOT=JK
                LDCUM=.TRUE.
             END IF
          END IF
!             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
        IF(LMFDUDV) THEN
             IF(JK.GE.KCBOT) THEN
                PUU(KLEV)=PUU(KLEV)+           &
                          PUEN(JK)*(PAPH(JK+1)-PAPH(JK))
                PVU(KLEV)=PVU(KLEV)+           &
                          PVEN(JK)*(PAPH(JK+1)-PAPH(JK))
             END IF
        END IF
  290 CONTINUE
      IF(LMFDUDV) THEN
         IF(LDCUM) THEN
            IKB=KCBOT
            ZZ=1./(PAPH(KLEVP1)-PAPH(IKB))
            PUU(KLEV)=PUU(KLEV)*ZZ
            PVU(KLEV)=PVU(KLEV)*ZZ
         ELSE
            PUU(KLEV)=PUEN(KLEVM1)
            PVU(KLEV)=PVEN(KLEVM1)
         END IF
      END IF
      RETURN
      END SUBROUTINE CUBASE

!**********************************************
!       SUBROUTINE CUTYPE
!********************************************** 
      SUBROUTINE CUTYPE    &
        ( KLEV,     KLEVP1,   KLEVM1,                  &
          PTENH,   PQENH,     PQSENH,    PGEOH,   PAPH,&
          RHO,      HFX,         QFX,    KTYPE,   lndj   )
!      THIS ROUTINE CALCULATES CLOUD BASE and TOP
!      AND RETURN CLOUD TYPES
!      ZHANG & WANG      IPRC           12/2010
!***PURPOSE.
!   --------
!          TO PRODUCE CLOUD TYPE for CU-PARAMETERIZATIONS
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD TYPES AS FOLLOWS;
!                 KTYPE=1 FOR deep cumulus
!                 KTYPE=2 FOR shallow cumulus
!***METHOD.
!  --------
!          based on a simplified updraught equation
!            partial(Hup)/partial(z)=eta(H - Hup)
!            eta is the entrainment rate for test parcel
!            H stands for dry static energy or the total water specific humidity
!            references: Christian Jakob, 2003: A new subcloud model for mass-flux convection schemes
!                        influence on triggering, updraft properties, and model climate, Mon.Wea.Rev.
!                        131, 2765-2778
!            and
!                        IFS Documentation - Cy33r1 
!          
!***EXTERNALS
!   ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   klevm1
      INTEGER   JK,IS,IK,ICALL,IKB,LEVELS
      REAL     PTENH(KLEV),       PQENH(KLEV), &
                                       PQSENH(KLEV),&
               PGEOH(KLEV),       PAPH(KLEVP1)
      REAL     ZRELH
      REAL     QFX,RHO,HFX
      REAL     ZQOLD(KLEV),       ZPH
      INTEGER  KCTOP,KCBOT
      INTEGER  KTYPE,LCLFLAG
      LOGICAL  TOPFLAG,DEEPFLAG,MYFLAG

      REAL     part1, part2, root
      REAL     conw,deltT,deltQ
      REAL     eta,dz,coef
      REAL     dhen(KLEV), dh(KLEV),qh(KLEV)
      REAL      Tup(KLEV),Qup(KLEV),ql(KLEV)
      REAL       ww(KLEV),Kup(KLEV)
      REAL     Vtup(KLEV),Vten(KLEV),buoy(KLEV)

      INTEGER  lndj
      REAL     CRIRH1
!***INPUT VARIABLES:
!       PTENH [ZTENH] - Environment Temperature on half levels. (CUINI)
!       PQENH [ZQENH] - Env. specific humidity on half levels. (CUINI)
!       PGEOH [ZGEOH] - Geopotential on half levels, (MSSFLX)
!       PAPH - Pressure of half levels. (MSSFLX)
!       RHO  - Density of the lowest Model level
!       QFX  - net upward moisture flux at the surface (kg/m^2/s)
!       HFX  - net upward heat flux at the surface (W/m^2)
!***VARIABLES OUTPUT BY CUTYPE:
!       KTYPE - Convection type - 1: Penetrative  (CUMASTR)
!                                 2: Stratocumulus (CUMASTR)
!                                 3: Mid-level  (CUASC)
!--------------------------------------------------------------
        KCBOT=KLEVM1
        KCTOP=KLEVM1
        KTYPE=0
!-----------------------------------------------------------
! let's do test,and check the shallow convection first
! the first level is JK+1
! define deltaT and deltaQ
!-----------------------------------------------------------
      DO JK=1,KLEV
        ZQOLD(JK)=0.0
           ql(jk)=0.0  ! parcel liquid water
          Tup(jk)=0.0  ! parcel temperature
          Qup(jk)=0.0  ! parcel specific humidity
           dh(jk)=0.0  ! parcel dry static energy
           qh(jk)=0.0  ! parcel total water specific humidity
           ww(jk)=0.0  ! parcel vertical speed (m/s)
         dhen(jk)=0.0  ! environment dry static energy
          Kup(jk)=0.0  ! updraught kinetic energy for parcel
         Vtup(jk)=0.0  ! parcel virtual temperature considering water-loading
         Vten(jk)=0.0  ! environment virtual temperature
         buoy(jk)=0.0  ! parcel buoyancy
      END DO

         lclflag = 0  ! flag for the condensation level
         conw    = 0.0 ! convective-scale velocity,also used for the vertical speed at the first level
         myflag  = .true. ! just as input for cuadjqt subroutine
        topflag  = .false.! flag for whether the cloud top is found

! check the levels from lowest level to second top level
      do JK=KLEVM1,2,-1
          ZPH=PAPH(JK)

! define the variables at the first level      
      if(jk .eq. KLEVM1) then
        part1 = 1.5*0.4*pgeoh(jk+1)/(rho*ptenh(jk+1))
        part2 = hfx/cpd+0.61*ptenh(jk+1)*qfx
        root = 0.001-part1*part2
        if(root .gt. 0) then
          conw = 1.2*(root)**(1.0/3.0)
        else
          conw = -1.2*(-root)**(1.0/3.0)
        end if
        deltT = -1.5*hfx/(rho*cpd*conw)
        deltQ = -1.5*qfx/(rho*conw)

        Tup(jk+1) = ptenh(jk+1) + deltT
        Qup(jk+1) = pqenh(jk+1) + deltQ
         ql(jk+1) = 0.
         dh(jk+1) = pgeoh(jk+1) + Tup(jk+1)*cpd
         qh(jk+1) = pqenh(jk+1) + deltQ + ql(jk+1)
         ww(jk+1) = conw
      end if

! the next levels, we use the variables at the first level as initial values
      if(.not. topflag) then
        eta = 0.5*(0.55/(pgeoh(jk)*zrg)+1.0e-3)
        dz  = (pgeoh(jk)-pgeoh(jk+1))*zrg
        coef= eta*dz
        dhen(jk) = pgeoh(jk) + cpd*ptenh(jk)
        dh(jk) = (coef*dhen(jk) + dh(jk+1))/(1+coef)
        qh(jk) = (coef*pqenh(jk)+ qh(jk+1))/(1+coef)
        Tup(jk) = (dh(jk)-pgeoh(jk))*RCPD
        Qup(jk) = qh(jk) - ql(jk+1)
        zqold(jk) = Qup(jk)
      end if
! check if the parcel is saturated
      ik=jk
      icall=1
      call CUADJTQ(klev,ik,zph,Tup,Qup,myflag,icall)
        if( .not. topflag .and. zqold(jk) .ne. Qup(jk) ) then
          lclflag = lclflag + 1
          ql(jk) = ql(jk+1) + zqold(jk) - Qup(jk)
          dh(jk) = pgeoh(jk) + cpd*Tup(jk)
        end if

! compute the updraft speed
        if(.not. topflag)then
          Kup(jk+1) = 0.5*ww(jk+1)**2
          Vtup(jk) = Tup(jk)*(1.+VTMPC1*Qup(jk)-ql(jk))
          Vten(jk) = ptenh(jk)*(1.+VTMPC1*pqenh(jk))
          buoy(jk) = (Vtup(jk) - Vten(jk))/Vten(jk)*g
          Kup(jk)  = (Kup(jk+1) + 0.333*dz*buoy(jk))/ &
                        (1+2*2*eta*dz)
          if(Kup(jk) .gt. 0 ) then
             ww(jk) = sqrt(2*Kup(jk))
             if(lclflag .eq. 1 ) kcbot = jk
             if(jk .eq. 2) then
                kctop = jk
                topflag= .true.
             end if
          else
             ww(jk) = 0
             kctop = jk + 1
             topflag = .true.
          end if
         end if
      end do ! end all the levels

        if(paph(kcbot) - paph(kctop) .lt. ZDNOPRC .and. &
          paph(kcbot) - paph(kctop) .gt. 0 &
           .and. lclflag .gt. 0) then
           ktype = 2
         end if

!-----------------------------------------------------------
! Next, let's check the deep convection
! the first level is JK
! define deltaT and deltaQ
!----------------------------------------------------------
! we check the parcel starting level by level (from the second lowest level to the next 12th level,
! usually, the 12th level around 700 hPa for common eta levels)
      do levels=KLEVM1-1,KLEVM1-12,-1
      DO JK=1,KLEV
        ZQOLD(JK)=0.0
           ql(jk)=0.0  ! parcel liquid water
          Tup(jk)=0.0  ! parcel temperature
          Qup(jk)=0.0  ! parcel specific humidity
           dh(jk)=0.0  ! parcel dry static energy
           qh(jk)=0.0  ! parcel total water specific humidity
           ww(jk)=0.0  ! parcel vertical speed (m/s)
         dhen(jk)=0.0  ! environment dry static energy
          Kup(jk)=0.0  ! updraught kinetic energy for parcel
         Vtup(jk)=0.0  ! parcel virtual temperature considering water-loading
         Vten(jk)=0.0  ! environment virtual temperature
         buoy(jk)=0.0  ! parcel buoyancy
      END DO

         lclflag = 0  ! flag for the condensation level
         kctop = levels
         kcbot = levels
         myflag  = .true. ! just as input for cuadjqt subroutine
        topflag  = .false.! flag for whether the cloud top is found

! check the levels from lowest level to second top level
      do JK=levels,2,-1
          ZPH=PAPH(JK)

! define the variables at the first level      
      if(jk .eq. levels) then
        deltT = 0.2
        deltQ = 1.0e-4

        if(paph(KLEVM1-1)-paph(jk) .le. 6.e3) then
         ql(jk+1) = 0.
        Tup(jk+1) = 0.25*(ptenh(jk+1)+ptenh(jk)+ &
                             ptenh(jk-1)+ptenh(jk-2)) + &
                      deltT
        dh(jk+1) = 0.25*(pgeoh(jk+1)+pgeoh(jk)+ &
                            pgeoh(jk-1)+pgeoh(jk-2)) + &
                      Tup(jk+1)*cpd 
        qh(jk+1) = 0.25*(pqenh(jk+1)+pqenh(jk)+ &
                            pqenh(jk-1)+pqenh(jk-2))+ &
                      deltQ + ql(jk+1)
        Qup(jk+1) = qh(jk+1) - ql(jk+1)
        else
         ql(jk+1) = 0.
        Tup(jk+1) = ptenh(jk+1) + deltT
         dh(jk+1) = pgeoh(jk+1) + Tup(jk+1)*cpd
         qh(jk+1) = pqenh(jk+1) + deltQ
        Qup(jk+1) =    qh(jk+1) - ql(jk+1)
        end if
      ww(jk+1) = 1.0

      end if

! the next levels, we use the variables at the first level as initial values
      if(.not. topflag) then
        eta = 1.1e-4
        dz  = (pgeoh(jk)-pgeoh(jk+1))*zrg
        coef= eta*dz
        dhen(jk) = pgeoh(jk) + cpd*ptenh(jk)
        dh(jk) = (coef*dhen(jk) + dh(jk+1))/(1+coef)
        qh(jk) = (coef*pqenh(jk)+ qh(jk+1))/(1+coef)
        Tup(jk) = (dh(jk)-pgeoh(jk))*RCPD
        Qup(jk) = qh(jk) - ql(jk+1)
        zqold(jk) = Qup(jk)
      end if
! check if the parcel is saturated
      ik=jk
      icall=1
      call CUADJTQ(klev,ik,zph,Tup,Qup,myflag,icall)
        if( .not. topflag .and. zqold(jk) .ne. Qup(jk) ) then
          lclflag = lclflag + 1
          ql(jk) = ql(jk+1) + zqold(jk) - Qup(jk)
          dh(jk) = pgeoh(jk) + cpd*Tup(jk)
        end if

! compute the updraft speed
        if(.not. topflag)then
          Kup(jk+1) = 0.5*ww(jk+1)**2
          Vtup(jk) = Tup(jk)*(1.+VTMPC1*Qup(jk)-ql(jk))
          Vten(jk) = ptenh(jk)*(1.+VTMPC1*pqenh(jk))
          buoy(jk) = (Vtup(jk) - Vten(jk))/Vten(jk)*g
          Kup(jk)  = (Kup(jk+1) + 0.333*dz*buoy(jk))/ &
                        (1+2*2*eta*dz)
          if(Kup(jk) .gt. 0 ) then
             ww(jk) = sqrt(2*Kup(jk))
             if(lclflag .eq. 1 ) kcbot = jk
             if(jk .eq. 2) then
                kctop = jk
                topflag= .true.
             end if
          else
             ww(jk) = 0
             kctop = jk + 1
             topflag = .true.
          end if
         end if
      end do ! end all the levels

       if(paph(kcbot) - paph(kctop) .gt. ZDNOPRC .and. &
              lclflag .gt. 0 ) then
         ZRELH = 0.
         do jk=kcbot,kctop,-1
           ZRELH=ZRELH+ PQENH(JK)/PQSENH(JK)
         end do
         ZRELH = ZRELH/(kcbot-kctop+1)

         if(lndj .eq. 1) then
           CRIRH1 = CRIRH*0.8
         else
           CRIRH1 = CRIRH
         end if
         if(ZRELH .ge. CRIRH1) ktype  = 1
       end if

       end do ! end all cycles

      END SUBROUTINE CUTYPE

!
!**********************************************
!       SUBROUTINE CUASC_NEW
!********************************************** 
      SUBROUTINE CUASC_NEW &
         (KLEV,     KLEVP1,   KLEVM1,   PTENH,            &
          PQENH,    PUEN,     PVEN,     PTEN,     PQEN,   &
          PQSEN,    PGEO,     PGEOH,    PAP,      PAPH,   &
          PQTE,     PVERV,    KLWMIN,   LDCUM,    PHCBASE,& 
          KTYPE,    KLAB,     PTU,      PQU,      PLU,    &
          PUU,      PVU,      PMFU,     PMFUB,    PENTR,  &
          PMFUS,    PMFUQ,    PMFUL,    PLUDE,    PDMFUP, & 
          KCBOT,    KCTOP,    KCTOP0,   KCUM,     ZTMST,  &
          KHMIN,    PHHATT,   PQSENH)
!     THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!     FOR CUMULUS PARAMETERIZATION
!     M.TIEDTKE         E.C.M.W.F.     7/86 MODIF.  12/89
!     Y.WANG            IPRC           11/01 MODIF.
!***PURPOSE.
!   --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!***METHOD.
!  --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!***EXTERNALS
!   ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR_NEW*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!***REFERENCE
!   ---------
!          (TIEDTKE,1989)
!***INPUT VARIABLES:
!       PTENH [ZTENH] - Environ Temperature on half levels. (CUINI)
!       PQENH [ZQENH] - Env. specific humidity on half levels. (CUINI)
!       PUEN - Environment wind u-component. (MSSFLX)
!       PVEN - Environment wind v-component. (MSSFLX)
!       PTEN - Environment Temperature. (MSSFLX)
!       PQEN - Environment Specific Humidity. (MSSFLX)
!       PQSEN - Environment Saturation Specific Humidity. (MSSFLX)
!       PGEO - Geopotential. (MSSFLX)
!       PGEOH [ZGEOH] - Geopotential on half levels, (MSSFLX)
!       PAP - Pressure in Pa.  (MSSFLX)
!       PAPH - Pressure of half levels. (MSSFLX)
!       PQTE - Moisture convergence (Delta q/Delta t). (MSSFLX)
!       PVERV - Large Scale Vertical Velocity (Omega). (MSSFLX)
!       KLWMIN [ILWMIN] - Level of Minimum Omega. (CUINI)
!       KLAB [ILAB] - Level Label - 1: Sub-cloud layer.
!                                   2: Condensation Level (Cloud Base)
!       PMFUB [ZMFUB] - Updraft Mass Flux at Cloud Base. (CUMASTR)
!***VARIABLES MODIFIED BY CUASC:
!       LDCUM - Logical denoting profiles. (CUBASE)
!       KTYPE - Convection type - 1: Penetrative  (CUMASTR)
!                                 2: Stratocumulus (CUMASTR)
!                                 3: Mid-level  (CUASC)
!       PTU - Cloud Temperature.
!       PQU - Cloud specific Humidity.
!       PLU - Cloud Liquid Water (Moisture condensed out)
!       PUU [ZUU] - Cloud Momentum U-Component.
!       PVU [ZVU] - Cloud Momentum V-Component.
!       PMFU - Updraft Mass Flux.
!       PENTR [ZENTR] - Entrainment Rate. (CUMASTR ) (CUBASMC)
!       PMFUS [ZMFUS] - Updraft Flux of Dry Static Energy. (CUBASMC)
!       PMFUQ [ZMFUQ] - Updraft Flux of Specific Humidity.
!       PMFUL [ZMFUL] - Updraft Flux of Cloud Liquid Water.
!       PLUDE - Liquid Water Returned to Environment by Detrainment.
!       PDMFUP [ZMFUP] - FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS
!       KCBOT - Cloud Base Level. (CUBASE)
!       KCTOP -
!       KCTOP0 [ICTOP0] - Estimate of Cloud Top. (CUMASTR)
!       KCUM [ICUM] -
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   klevm1,kcum
      REAL      ZTMST,ZCONS2,ZDZ,ZDRODZ
      INTEGER   JK,IKB,IK,IS,IKT,ICALL
      REAL      ZMFMAX,ZFAC,ZMFTEST,ZDPRHO,ZMSE,ZNEVN,ZODMAX
      REAL      ZQEEN,ZSEEN,ZSCDE,ZGA,ZDT,ZSCOD
      REAL      ZQUDE,ZQCOD, ZMFUSK, ZMFUQK,ZMFULK
      REAL      ZBUO, ZPRCON, ZLNEW, ZZ, ZDMFEU, ZDMFDU
      REAL      ZBUOYZ,ZZDMF
      REAL     PTENH(KLEV),       PQENH(KLEV), &
              PUEN(KLEV),        PVEN(KLEV),   &
              PTEN(KLEV),        PQEN(KLEV),   &
              PGEO(KLEV),        PGEOH(KLEV),  &
              PAP(KLEV),         PAPH(KLEVP1), &
              PQSEN(KLEV),       PQTE(KLEV),   &
              PVERV(KLEV),       PQSENH(KLEV)  
      REAL     PTU(KLEV),         PQU(KLEV),   &
              PUU(KLEV),         PVU(KLEV),    &
              PMFU(KLEV),        ZPH,         &
              PMFUB,            PENTR,       &
              PMFUS(KLEV),       PMFUQ(KLEV),  &
              PLU(KLEV),         PLUDE(KLEV),  &
              PMFUL(KLEV),       PDMFUP(KLEV)
      REAL     ZDMFEN,           ZDMFDE,     &
              ZMFUU,            ZMFUV,       &
              ZPBASE,           ZQOLD,       &
              PHHATT(KLEV),      ZODETR(KLEV), &
              ZOENTR(KLEV),      ZBUOY
      REAL     PHCBASE
      INTEGER  KLWMIN,           KTYPE,      &
              KLAB(KLEV),        KCBOT,       &
              KCTOP,            KCTOP0,      &
              KHMIN
      LOGICAL LDCUM,            LOFLAG
      integer leveltop,levelbot
      real    tt,ttb
      real    zqsat, zqsatb
      real    fscale

!--------------------------------
!*    1.       SPECIFY PARAMETERS
!--------------------------------
  100 CONTINUE
      ZCONS2=1./(G*ZTMST)
!---------------------------------
!     2.        SET DEFAULT VALUES
!---------------------------------
  200 CONTINUE
        ZMFUU=0.
        ZMFUV=0.
        ZBUOY=0.
        IF(.NOT.LDCUM) KTYPE=0
      DO 230 JK=1,KLEV
          PLU(JK)=0.
          PMFU(JK)=0.
          PMFUS(JK)=0.
          PMFUQ(JK)=0.
          PMFUL(JK)=0.
          PLUDE(JK)=0.
          PDMFUP(JK)=0.
          ZOENTR(JK)=0.
          ZODETR(JK)=0.
          IF(.NOT.LDCUM.OR.KTYPE.EQ.3) KLAB(JK)=0
          IF(.NOT.LDCUM.AND.PAPH(JK).LT.4.E4) KCTOP0=JK
  230 CONTINUE
!------------------------------------------------
!     3.0      INITIALIZE VALUES AT LIFTING LEVEL
!------------------------------------------------
        KCTOP=KLEVM1
        IF(.NOT.LDCUM) THEN
           KCBOT=KLEVM1
           PMFUB=0.
           PQU(KLEV)=0.
        END IF
        PMFU(KLEV)=PMFUB
        PMFUS(KLEV)=PMFUB*(CPD*PTU(KLEV)+PGEOH(KLEV))
        PMFUQ(KLEV)=PMFUB*PQU(KLEV)
        IF(LMFDUDV) THEN
           ZMFUU=PMFUB*PUU(KLEV)
           ZMFUV=PMFUB*PVU(KLEV)
        END IF
!
!-- 3.1 Find organized entrainment at cloud base
!
      LDCUM=.FALSE.
      IF (KTYPE.EQ.1) THEN
       IKB = KCBOT
       if(orgen .eq. 1 ) then
! old scheme
       ZBUOY=G*((PTU(IKB)-PTENH(IKB))/PTENH(IKB)+ &
               0.608*(PQU(IKB)-PQENH(IKB)))
       IF (ZBUOY.GT.0.) THEN
        ZDZ = (PGEO(IKB-1)-PGEO(IKB))*ZRG
        ZDRODZ = -LOG(PTEN(IKB-1)/PTEN(IKB))/ZDZ -  &
                 G/(RD*PTENH(IKB))
        ZOENTR(IKB-1)=ZBUOY*0.5/(1.+ZBUOY*ZDZ) &
                +ZDRODZ
        ZOENTR(IKB-1) = MIN(ZOENTR(IKB-1),1.E-3)
        ZOENTR(IKB-1) = MAX(ZOENTR(IKB-1),0.)
       END IF
! New scheme
! Let's define the fscale
        else if(orgen .eq. 2 ) then
        tt = ptenh(ikb-1)
        zqsat = TLUCUA(tt)/paph(ikb-1)
        zqsat = zqsat/(1.-VTMPC1*zqsat)
        ttb = ptenh(ikb)
        zqsatb = TLUCUA(ttb)/paph(ikb)
        zqsatb = zqsatb/(1.-VTMPC1*zqsatb)
        fscale = (zqsat/zqsatb)**3
! end of defining the fscale
        zoentr(ikb-1) = 1.E-3*(1.3-PQEN(ikb-1)/PQSEN(ikb-1))*fscale
        zoentr(ikb-1) = MIN(zoentr(ikb-1),1.E-3)
        zoentr(ikb-1) = MAX(zoentr(ikb-1),0.)
       end if
      END IF
!
!-----------------------------------------------------------------
!     4.       DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!              BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!              BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!              THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!-----------------------------------------------------------------
  400 CONTINUE

! let's define the levels in which the middle level convection could be activated
      do jk=KLEVM1,2,-1
      if(abs(paph(jk)*0.01 - 250) .lt. 50.) then
        leveltop = jk
        exit
      end if
      end do
      leveltop = min(KLEV-15,leveltop)
      levelbot = KLEVM1 - 4
        
      DO 480 JK=KLEVM1,2,-1
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
! ---------------------------------------------------------------------
      IK=JK
      IF(LMFMID.AND.IK.LT.levelbot.AND.IK.GT.leveltop) THEN
      CALL CUBASMC  &
         (KLEV,     KLEVM1,   IK,      PTEN,            &
          PQEN,     PQSEN,    PUEN,     PVEN,    PVERV, &
          PGEO,     PGEOH,    LDCUM,    KTYPE,   KLAB,  &
          PMFU,     PMFUB,    PENTR,    KCBOT,   PTU,   &
          PQU,      PLU,      PUU,     PVU,      PMFUS, &
          PMFUQ,    PMFUL,    PDMFUP,  ZMFUU,    ZMFUV)
      ENDIF
      IS=0
        ZQOLD=0.0
        IS=IS+KLAB(JK+1)
        IF(KLAB(JK+1).EQ.0) KLAB(JK)=0
        LOFLAG=KLAB(JK+1).GT.0
        ZPH=PAPH(JK)
        IF(KTYPE.EQ.3.AND.JK.EQ.KCBOT) THEN
           ZMFMAX=(PAPH(JK)-PAPH(JK-1))*ZCONS2
           IF(PMFUB.GT.ZMFMAX) THEN
              ZFAC=ZMFMAX/PMFUB
              PMFU(JK+1)=PMFU(JK+1)*ZFAC
              PMFUS(JK+1)=PMFUS(JK+1)*ZFAC
              PMFUQ(JK+1)=PMFUQ(JK+1)*ZFAC
              ZMFUU=ZMFUU*ZFAC
              ZMFUV=ZMFUV*ZFAC
              PMFUB=ZMFMAX
           END IF
        END IF
      IF(IS.EQ.0) GO TO 480
!
!*     SPECIFY ENTRAINMENT RATES IN *CUENTR_NEW*
! -------------------------------------
      IK=JK
      CALL CUENTR_NEW &
         (KLEV,     KLEVP1,   IK,       PTENH,          &
          PAPH,     PAP,      PGEOH,    KLWMIN,   LDCUM,&
          KTYPE,    KCBOT,    KCTOP0,   ZPBASE,   PMFU, &
          PENTR,    ZDMFEN,   ZDMFDE,   ZODETR,   KHMIN)
!
!      DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
! -------------------------------------------------------
! Do adiabatic ascent for entraining/detraining plume
! the cloud ensemble entrains environmental values
! in turbulent detrainment cloud ensemble values are detrained
! in organized detrainment the dry static energy and
! moisture that are neutral compared to the
! environmental air are detrained
!
      IF(LOFLAG) THEN
        IF(JK.LT.KCBOT) THEN
         ZMFTEST=PMFU(JK+1)+ZDMFEN-ZDMFDE
         ZMFMAX=MIN(ZMFTEST,(PAPH(JK)-PAPH(JK-1))*ZCONS2)
         ZDMFEN=MAX(ZDMFEN-MAX(ZMFTEST-ZMFMAX,0.),0.)
        END IF
        ZDMFDE=MIN(ZDMFDE,0.75*PMFU(JK+1))
        PMFU(JK)=PMFU(JK+1)+ZDMFEN-ZDMFDE
        IF (JK.LT.kcbot) THEN
          zdprho = (pgeoh(jk)-pgeoh(jk+1))*zrg
          zoentr(jk) = zoentr(jk)*zdprho*pmfu(jk+1)
          zmftest = pmfu(jk) + zoentr(jk)-zodetr(jk)
          zmfmax = MIN(zmftest,(paph(jk)-paph(jk-1))*zcons2)
          zoentr(jk) = MAX(zoentr(jk)-MAX(zmftest-zmfmax,0.),0.)
        END IF
!
! limit organized detrainment to not allowing for too deep clouds
!
        IF (ktype.EQ.1.AND.jk.LT.kcbot.AND.jk.LE.khmin) THEN
          zmse = cpd*ptu(jk+1) + alv*pqu(jk+1) + pgeoh(jk+1)
          ikt = kctop0
          znevn=(pgeoh(ikt)-pgeoh(jk+1))*(zmse-phhatt(jk+1))*zrg
          IF (znevn.LE.0.) znevn = 1.
          zdprho = (pgeoh(jk)-pgeoh(jk+1))*zrg
          zodmax = ((phcbase-zmse)/znevn)*zdprho*pmfu(jk+1)
          zodmax = MAX(zodmax,0.)
          zodetr(jk) = MIN(zodetr(jk),zodmax)
        END IF
        zodetr(jk) = MIN(zodetr(jk),0.75*pmfu(jk))
        pmfu(jk) = pmfu(jk) + zoentr(jk) - zodetr(jk)
        ZQEEN=PQENH(JK+1)*ZDMFEN
        zqeen=zqeen + pqenh(jk+1)*zoentr(jk)
        ZSEEN=(CPD*PTENH(JK+1)+PGEOH(JK+1))*ZDMFEN
        zseen=zseen+(cpd*ptenh(jk+1)+pgeoh(jk+1))*  &
             zoentr(jk)
        ZSCDE=(CPD*PTU(JK+1)+PGEOH(JK+1))*ZDMFDE
! find moist static energy that give nonbuoyant air
        zga = alv*pqsenh(jk+1)/(rv*(ptenh(jk+1)**2))
        zdt = (plu(jk+1)-0.608*(pqsenh(jk+1)-pqenh(jk+1)))&
            /(1./ptenh(jk+1)+0.608*zga)
        zscod = cpd*ptenh(jk+1) + pgeoh(jk+1) + cpd*zdt
        zscde = zscde + zodetr(jk)*zscod
        zqude = pqu(jk+1)*zdmfde
        zqcod = pqsenh(jk+1) + zga*zdt
        zqude = zqude + zodetr(jk)*zqcod
        plude(jk) = plu(jk+1)*zdmfde
        plude(jk) = plude(jk)+plu(jk+1)*zodetr(jk)
        zmfusk = pmfus(jk+1) + zseen - zscde
        zmfuqk = pmfuq(jk+1) + zqeen - zqude
        zmfulk = pmful(jk+1) - plude(jk)
        plu(jk) = zmfulk*(1./MAX(cmfcmin,pmfu(jk)))
        pqu(jk) = zmfuqk*(1./MAX(cmfcmin,pmfu(jk)))
        ptu(jk)=(zmfusk*(1./MAX(cmfcmin,pmfu(jk)))-  &
            pgeoh(jk))*rcpd
        ptu(jk) = MAX(100.,ptu(jk))
        ptu(jk) = MIN(400.,ptu(jk))
        zqold = pqu(jk)
      END IF
!*             DO CORRECTIONS FOR MOIST ASCENT
!*             BY ADJUSTING T,Q AND L IN *CUADJTQ*
!------------------------------------------------
      IK=JK
      ICALL=1
!
      CALL CUADJTQ(KLEV,IK,ZPH,PTU,PQU,LOFLAG,ICALL)
!
      IF(LOFLAG.AND.PQU(JK).NE.ZQOLD) THEN
         KLAB(JK)=2
         PLU(JK)=PLU(JK)+ZQOLD-PQU(JK)
         ZBUO=PTU(JK)*(1.+VTMPC1*PQU(JK)-PLU(JK))-  &
        PTENH(JK)*(1.+VTMPC1*PQENH(JK))
         IF(KLAB(JK+1).EQ.1) ZBUO=ZBUO+ZBUO0
         IF(ZBUO.GT.0..AND.PMFU(JK).GT.0.01*PMFUB.AND. &
                            JK.GE.KCTOP0) THEN
            KCTOP=JK
            LDCUM=.TRUE.
            IF(ZPBASE-PAPH(JK).GE.ZDNOPRC) THEN
               ZPRCON=CPRCON
            ELSE
               ZPRCON=0.
            ENDIF
            ZLNEW=PLU(JK)/(1.+ZPRCON*(PGEOH(JK)-PGEOH(JK+1)))
            PDMFUP(JK)=MAX(0.,(PLU(JK)-ZLNEW)*PMFU(JK))
            PLU(JK)=ZLNEW
         ELSE
            KLAB(JK)=0
            PMFU(JK)=0.
         END IF
      END IF
      IF(LOFLAG) THEN
         PMFUL(JK)=PLU(JK)*PMFU(JK)
         PMFUS(JK)=(CPD*PTU(JK)+PGEOH(JK))*PMFU(JK)
         PMFUQ(JK)=PQU(JK)*PMFU(JK)
      END IF
!
      IF(LMFDUDV) THEN
!
        zdmfen = zdmfen + zoentr(jk)
        zdmfde = zdmfde + zodetr(jk)
           IF(LOFLAG) THEN
              IF(KTYPE.EQ.1.OR.KTYPE.EQ.3) THEN
                 IF(ZDMFEN.LE.1.E-20) THEN
                    ZZ=3.
                 ELSE
                    ZZ=2.
                 ENDIF
              ELSE
                 IF(ZDMFEN.LE.1.0E-20) THEN
                    ZZ=1.
                 ELSE
                    ZZ=0.
                 ENDIF
              END IF
              ZDMFEU=ZDMFEN+ZZ*ZDMFDE
              ZDMFDU=ZDMFDE+ZZ*ZDMFDE
              ZDMFDU=MIN(ZDMFDU,0.75*PMFU(JK+1))
              ZMFUU=ZMFUU+                              &
                       ZDMFEU*PUEN(JK)-ZDMFDU*PUU(JK+1)   
              ZMFUV=ZMFUV+                              &
                       ZDMFEU*PVEN(JK)-ZDMFDU*PVU(JK+1)   
              IF(PMFU(JK).GT.0.) THEN
                 PUU(JK)=ZMFUU*(1./PMFU(JK))
                 PVU(JK)=ZMFUV*(1./PMFU(JK))
              END IF
           END IF
!
        END IF
!
! Compute organized entrainment
! for use at next level
!
       IF (loflag.AND.ktype.EQ.1) THEN
! old scheme
       if(orgen .eq. 1 ) then
        zbuoyz=g*((ptu(jk)-ptenh(jk))/ptenh(jk)+  &
              0.608*(pqu(jk)-pqenh(jk))-plu(jk))
        zbuoyz = MAX(zbuoyz,0.0)
        zdz = (pgeo(jk-1)-pgeo(jk))*zrg
        zdrodz = -LOG(pten(jk-1)/pten(jk))/zdz -  &
                 g/(rd*ptenh(jk))
        zbuoy = zbuoy + zbuoyz*zdz
        zoentr(jk-1) = zbuoyz*0.5/(1.+zbuoy)+zdrodz
        zoentr(jk-1) = MIN(zoentr(jk-1),1.E-3)
        zoentr(jk-1) = MAX(zoentr(jk-1),0.)
       else if(orgen .eq. 2 ) then
! Let's define the fscale
        tt = ptenh(jk-1)
        zqsat = TLUCUA(tt)/paph(jk-1)
        zqsat = zqsat/(1.-VTMPC1*zqsat)
        ttb = ptenh(kcbot)
        zqsatb = TLUCUA(ttb)/paph(kcbot)
        zqsatb = zqsatb/(1.-VTMPC1*zqsatb)
        fscale = (zqsat/zqsatb)**3
! end of defining the fscale
        zoentr(jk-1) = 1.E-3*(1.3-PQEN(jk-1)/PQSEN(jk-1))*fscale
        zoentr(jk-1) = MIN(zoentr(jk-1),1.E-3)
        zoentr(jk-1) = MAX(zoentr(jk-1),0.)
!        write(6,*) "zoentr=",zoentr(jk-1) 
       end if
       END IF
!
  480 CONTINUE
! -----------------------------------------------------------------
!     5.       DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
! -----------------------------------------------------------------
!                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!                         FROM PREVIOUS CALCULATIONS ABOVE)
  500 CONTINUE
      IF(KCTOP.EQ.KLEVM1) LDCUM=.FALSE.
      KCBOT=MAX(KCBOT,KCTOP)
      IS=0
      IF(LDCUM) THEN
         IS=IS+1
      ENDIF
      KCUM=IS
      IF(IS.EQ.0) GO TO 800
      IF(LDCUM) THEN
         JK=KCTOP-1
         ZZDMF=CMFCTOP
         ZDMFDE=(1.-ZZDMF)*PMFU(JK+1)
         PLUDE(JK)=ZDMFDE*PLU(JK+1)
         PMFU(JK)=PMFU(JK+1)-ZDMFDE
         PMFUS(JK)=(CPD*PTU(JK)+PGEOH(JK))*PMFU(JK)
         PMFUQ(JK)=PQU(JK)*PMFU(JK)
         PMFUL(JK)=PLU(JK)*PMFU(JK)
         PLUDE(JK-1)=PMFUL(JK)
         PDMFUP(JK)=0.
      END IF
        IF(LMFDUDV) THEN
           IF(LDCUM) THEN
              JK=KCTOP-1
              PUU(JK)=PUU(JK+1)
              PVU(JK)=PVU(JK+1)
           END IF
        END IF
  800 CONTINUE
      RETURN
      END SUBROUTINE CUASC_NEW
!

!**********************************************
!       SUBROUTINE CUDLFS
!********************************************** 
      SUBROUTINE CUDLFS &
         (KLEV,     KLEVP1,   PTENH,    PQENH,            &
          PUEN,     PVEN,     PGEOH,    PAPH,     PTU,    &
          PQU,      PUU,      PVU,      LDCUM,    KCBOT,  &
          KCTOP,    PMFUB,    PRFL,     PTD,      PQD,    &
          PUD,      PVD,      PMFD,     PMFDS,    PMFDQ,  &
          PDMFDP,   KDTOP,    LDDRAF)
!      THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!      CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!      M.TIEDTKE         E.C.M.W.F.    12/86 MODIF.  12/89
!***PURPOSE.
!   --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!***METHOD.
!  --------
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!***EXTERNALS
!   ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KE,JK,IS,IK,ICALL
      REAL      ZTTEST, ZQTEST, ZBUO, ZMFTOP
      REAL     PTENH(KLEV),       PQENH(KLEV),   &
              PUEN(KLEV),        PVEN(KLEV),     &
              PGEOH(KLEV),       PAPH(KLEVP1),   &
              PTU(KLEV),         PQU(KLEV),      &
              PUU(KLEV),         PVU(KLEV),      &
              PMFUB,            PRFL
      REAL     PTD(KLEV),         PQD(KLEV),     &
              PUD(KLEV),         PVD(KLEV),      &
              PMFD(KLEV),        PMFDS(KLEV),    &
              PMFDQ(KLEV),       PDMFDP(KLEV)    
      REAL     ZTENWB(KLEV),      ZQENWB(KLEV),  &
              ZCOND,            ZPH
      INTEGER  KCBOT,            KCTOP,        &
              KDTOP
      LOGICAL  LDCUM,            LLo2,         &
              LDDRAF
!-----------------------------------------------
!     1.       SET DEFAULT VALUES FOR DOWNDRAFTS
!-----------------------------------------------
  100 CONTINUE
      LDDRAF=.FALSE.
      KDTOP=KLEVP1
      IF(.NOT.LMFDD) GO TO 300
!------------------------------------------------------------
!     2.       DETERMINE LEVEL OF FREE SINKING BY
!              DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!              FOR EVERY POINT AND PROCEED AS FOLLOWS:
!                (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!                (2) DO MIXING WITH CUMULUS CLOUD AIR
!                (3) CHECK FOR NEGATIVE BUOYANCY
!              THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!              OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!              TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!              EVAPORATION OF RAIN AND CLOUD WATER)
!------------------------------------------------------------------
  200 CONTINUE
      KE=KLEV-3
      DO 290 JK=3,KE
!   2.1      CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!            FOR ENVIRONMENTAL AIR IN *CUADJTQ*
! -----------------------------------------------------
  210 CONTINUE
      IS=0
      ZTENWB(JK)=PTENH(JK)
      ZQENWB(JK)=PQENH(JK)
      ZPH=PAPH(JK)
      LLO2=LDCUM.AND.PRFL.GT.0..AND..NOT.LDDRAF.AND. &
              (JK.LT.KCBOT.AND.JK.GT.KCTOP)
      IF(LLO2)THEN
         IS=IS+1
      ENDIF
      IF(IS.EQ.0) GO TO 290
      IK=JK
      ICALL=2
      CALL CUADJTQ(KLEV,IK,ZPH,ZTENWB,ZQENWB,LLO2,ICALL)
!   2.2      DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!            AND CHECK FOR NEGATIVE BUOYANCY.
!            THEN SET VALUES FOR DOWNDRAFT AT LFS.
! -----------------------------------------------------
  220 CONTINUE
      IF(LLO2) THEN
         ZTTEST=0.5*(PTU(JK)+ZTENWB(JK))
         ZQTEST=0.5*(PQU(JK)+ZQENWB(JK))
         ZBUO=ZTTEST*(1.+VTMPC1*ZQTEST)-  &
             PTENH(JK)*(1.+VTMPC1*PQENH(JK))
         ZCOND=PQENH(JK)-ZQENWB(JK)
         ZMFTOP=-CMFDEPS*PMFUB
         IF(ZBUO.LT.0..AND.PRFL.GT.10.*ZMFTOP*ZCOND) THEN
            KDTOP=JK
            LDDRAF=.TRUE.
            PTD(JK)=ZTTEST
            PQD(JK)=ZQTEST
            PMFD(JK)=ZMFTOP
            PMFDS(JK)=PMFD(JK)*(CPD*PTD(JK)+PGEOH(JK))
            PMFDQ(JK)=PMFD(JK)*PQD(JK)
            PDMFDP(JK-1)=-0.5*PMFD(JK)*ZCOND
            PRFL=PRFL+PDMFDP(JK-1)
         END IF
      END IF
         IF(LMFDUDV) THEN
            IF(PMFD(JK).LT.0.) THEN
               PUD(JK)=0.5*(PUU(JK)+PUEN(JK-1))
               PVD(JK)=0.5*(PVU(JK)+PVEN(JK-1))
            END IF
         END IF
  290 CONTINUE
 300  CONTINUE
      RETURN
      END SUBROUTINE CUDLFS
!

!**********************************************
!       SUBROUTINE CUDDRAF
!********************************************** 
      SUBROUTINE CUDDRAF &
         (KLEV,     KLEVP1,   PTENH,    PQENH,           &
          PUEN,     PVEN,     PGEOH,    PAPH,     PRFL,  &
          LDDRAF,   PTD,      PQD,      PUD,      PVD,   &
          PMFD,     PMFDS,    PMFDQ,    PDMFDP)
!     THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!     M.TIEDTKE         E.C.M.W.F.    12/86 MODIF.  12/89
!***PURPOSE.
!   --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!***METHOD.
!  --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!***EXTERNALS
!   ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT
!***REFERENCE
!   ---------
!          (TIEDTKE,1989)
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   JK,IS,ITOPDE, IK, ICALL
      REAL      ZENTR,ZSEEN, ZQEEN, ZSDDE, ZQDDE,ZMFDSK, ZMFDQK
      REAL      ZBUO, ZDMFDP, ZMFDUK, ZMFDVK
      REAL     PTENH(KLEV),       PQENH(KLEV),  &
              PUEN(KLEV),        PVEN(KLEV),    &
              PGEOH(KLEV),       PAPH(KLEVP1) 
      REAL     PTD(KLEV),         PQD(KLEV),    &
              PUD(KLEV),         PVD(KLEV),     &
              PMFD(KLEV),        PMFDS(KLEV),   &
              PMFDQ(KLEV),       PDMFDP(KLEV),  &
              PRFL
      REAL     ZDMFEN,           ZDMFDE,      &
              ZCOND,            ZPH       
      LOGICAL  LDDRAF,           LLO2
!--------------------------------------------------------------
!     1.       CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!                (A) CALCULATING ENTRAINMENT RATES, ASSUMING
!                     LINEAR DECREASE OF MASSFLUX IN PBL
!                 (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!                     AND MOISTENING IS CALCULATED IN *CUADJTQ*
!                 (C) CHECKING FOR NEGATIVE BUOYANCY AND
!                     SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
! ----------------------------------------------------------------
  100 CONTINUE
      DO 180 JK=3,KLEV
      IS=0
      ZPH=PAPH(JK)
      LLO2=LDDRAF.AND.PMFD(JK-1).LT.0.
      IF(LLO2) THEN
         IS=IS+1
      ENDIF
      IF(IS.EQ.0) GO TO 180
      IF(LLO2) THEN
         ZENTR=ENTRDD*PMFD(JK-1)*RD*PTENH(JK-1)/   &
              (G*PAPH(JK-1))*(PAPH(JK)-PAPH(JK-1))
         ZDMFEN=ZENTR
         ZDMFDE=ZENTR
      END IF
      ITOPDE=KLEV-2
         IF(JK.GT.ITOPDE) THEN
            IF(LLO2) THEN
               ZDMFEN=0.
               ZDMFDE=PMFD(ITOPDE)*      &
              (PAPH(JK)-PAPH(JK-1))/     &
              (PAPH(KLEVP1)-PAPH(ITOPDE))
            END IF
         END IF
         IF(LLO2) THEN
            PMFD(JK)=PMFD(JK-1)+ZDMFEN-ZDMFDE
            ZSEEN=(CPD*PTENH(JK-1)+PGEOH(JK-1))*ZDMFEN
            ZQEEN=PQENH(JK-1)*ZDMFEN
            ZSDDE=(CPD*PTD(JK-1)+PGEOH(JK-1))*ZDMFDE
            ZQDDE=PQD(JK-1)*ZDMFDE
            ZMFDSK=PMFDS(JK-1)+ZSEEN-ZSDDE
            ZMFDQK=PMFDQ(JK-1)+ZQEEN-ZQDDE
            PQD(JK)=ZMFDQK*(1./MIN(-CMFCMIN,PMFD(JK)))
            PTD(JK)=(ZMFDSK*(1./MIN(-CMFCMIN,PMFD(JK)))- &
                       PGEOH(JK))*RCPD
            PTD(JK)=MIN(400.,PTD(JK))
            PTD(JK)=MAX(100.,PTD(JK))
            ZCOND=PQD(JK)
         END IF
      IK=JK
      ICALL=2
      CALL CUADJTQ(KLEV,IK,ZPH,PTD,PQD,LLO2,ICALL)
         IF(LLO2) THEN
            ZCOND=ZCOND-PQD(JK)
            ZBUO=PTD(JK)*(1.+VTMPC1*PQD(JK))- &
           PTENH(JK)*(1.+VTMPC1*PQENH(JK))
            IF(ZBUO.GE.0..OR.PRFL.LE.(PMFD(JK)*ZCOND)) THEN
               PMFD(JK)=0.
            ENDIF
            PMFDS(JK)=(CPD*PTD(JK)+PGEOH(JK))*PMFD(JK)
            PMFDQ(JK)=PQD(JK)*PMFD(JK)
            ZDMFDP=-PMFD(JK)*ZCOND
            PDMFDP(JK-1)=ZDMFDP
            PRFL=PRFL+ZDMFDP
         END IF
        IF(LMFDUDV) THEN
             IF(LLO2.AND.PMFD(JK).LT.0.) THEN
                ZMFDUK=PMFD(JK-1)*PUD(JK-1)+   &
               ZDMFEN*PUEN(JK-1)-ZDMFDE*PUD(JK-1)
                ZMFDVK=PMFD(JK-1)*PVD(JK-1)+   &
               ZDMFEN*PVEN(JK-1)-ZDMFDE*PVD(JK-1)
                PUD(JK)=ZMFDUK*(1./MIN(-CMFCMIN,PMFD(JK)))
                PVD(JK)=ZMFDVK*(1./MIN(-CMFCMIN,PMFD(JK)))
             END IF
        END IF
  180 CONTINUE
      RETURN
      END SUBROUTINE CUDDRAF
!

!**********************************************
!       SUBROUTINE CUFLX
!********************************************** 
      SUBROUTINE CUFLX &
         (KLEV,     KLEVP1,   PQEN,    PQSEN,            &
          PTENH,    PQENH,    PAPH,     PGEOH,   KCBOT,    &
          KCTOP,    KDTOP,    KTYPE,    LDDRAF,  LDCUM,  &
          PMFU,     PMFD,     PMFUS,    PMFDS,   PMFUQ,  &
          PMFDQ,    PMFUL,    PLUDE,    PDMFUP,  PDMFDP, &
          PRFL,     PRAIN,    PTEN,     PSFL,    PDPMEL, &
          KTOPM2,   ZTMST,    sig1)
!      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF.  12/89
!***PURPOSE
!   -------
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!***EXTERNALS
!   ---------
!          NONE
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KTOPM2, ITOP, JK, IKB
      REAL      ZTMST, ZCONS1, ZCONS2, ZCUCOV, ZTMELP2
      REAL      ZZP, ZFAC, ZSNMLT, ZRFL, CEVAPCU, ZRNEW
      REAL      ZRMIN, ZRFLN, ZDRFL, ZDPEVAP
      REAL     PQEN(KLEV),        PQSEN(KLEV),  &
              PTENH(KLEV),       PQENH(KLEV),   &
              PAPH(KLEVP1),      PGEOH(KLEV)    
      REAL     PMFU(KLEV),        PMFD(KLEV),   &
              PMFUS(KLEV),       PMFDS(KLEV),   &
              PMFUQ(KLEV),       PMFDQ(KLEV),   &
              PDMFUP(KLEV),      PDMFDP(KLEV),  &
              PMFUL(KLEV),       PLUDE(KLEV),   &
              PRFL,             PRAIN
      REAL     PTEN(KLEV),        PDPMEL(KLEV), &
              PSFL,             ZPSUBCL
      REAL     sig1(KLEV)
      INTEGER  KCBOT,            KCTOP,     &
              KDTOP,            KTYPE
      LOGICAL  LDDRAF,           LDCUM
!*       SPECIFY CONSTANTS
      ZCONS1=CPD/(ALF*G*ZTMST)
      ZCONS2=1./(G*ZTMST)
      ZCUCOV=0.05
      ZTMELP2=TMELT+2.
!*  1.0      DETERMINE FINAL CONVECTIVE FLUXES
!---------------------------------------------
  100 CONTINUE
      ITOP=KLEV
      PRFL=0.
      PSFL=0.
      PRAIN=0.
!     SWITCH OFF SHALLOW CONVECTION
      IF(.NOT.LMFSCV.AND.KTYPE.EQ.2)THEN
        LDCUM=.FALSE.
        LDDRAF=.FALSE.
      ENDIF
      ITOP=MIN(ITOP,KCTOP)
      IF(.NOT.LDCUM.OR.KDTOP.LT.KCTOP) LDDRAF=.FALSE.
      IF(.NOT.LDCUM) KTYPE=0
      KTOPM2=ITOP-2
      DO 120 JK=KTOPM2,KLEV
      IF(LDCUM.AND.JK.GE.KCTOP-1) THEN
         PMFUS(JK)=PMFUS(JK)-PMFU(JK)*  &
                     (CPD*PTENH(JK)+PGEOH(JK))
         PMFUQ(JK)=PMFUQ(JK)-PMFU(JK)*PQENH(JK)
         IF(LDDRAF.AND.JK.GE.KDTOP) THEN
            PMFDS(JK)=PMFDS(JK)-PMFD(JK)*  &
                        (CPD*PTENH(JK)+PGEOH(JK))
            PMFDQ(JK)=PMFDQ(JK)-PMFD(JK)*PQENH(JK)
         ELSE
            PMFD(JK)=0.
            PMFDS(JK)=0.
            PMFDQ(JK)=0.
            PDMFDP(JK-1)=0.
         END IF
      ELSE
         PMFU(JK)=0.
         PMFD(JK)=0.
         PMFUS(JK)=0.
         PMFDS(JK)=0.
         PMFUQ(JK)=0.
         PMFDQ(JK)=0.
         PMFUL(JK)=0.
         PDMFUP(JK-1)=0.
         PDMFDP(JK-1)=0.
         PLUDE(JK-1)=0.
      END IF
  120 CONTINUE
      DO 130 JK=KTOPM2,KLEV
      IF(LDCUM.AND.JK.GT.KCBOT) THEN
         IKB=KCBOT
         ZZP=((PAPH(KLEVP1)-PAPH(JK))/  &
             (PAPH(KLEVP1)-PAPH(IKB)))
         IF(KTYPE.EQ.3) THEN
            ZZP=ZZP**2
         ENDIF
         PMFU(JK)=PMFU(IKB)*ZZP
         PMFUS(JK)=PMFUS(IKB)*ZZP
         PMFUQ(JK)=PMFUQ(IKB)*ZZP
         PMFUL(JK)=PMFUL(IKB)*ZZP
      END IF
!*    2.        CALCULATE RAIN/SNOW FALL RATES
!*              CALCULATE MELTING OF SNOW
!*              CALCULATE EVAPORATION OF PRECIP
!----------------------------------------------
      IF(LDCUM) THEN
         PRAIN=PRAIN+PDMFUP(JK)
         IF(PTEN(JK).GT.TMELT) THEN
            PRFL=PRFL+PDMFUP(JK)+PDMFDP(JK)
            IF(PSFL.GT.0..AND.PTEN(JK).GT.ZTMELP2) THEN
               ZFAC=ZCONS1*(PAPH(JK+1)-PAPH(JK))
               ZSNMLT=MIN(PSFL,ZFAC*(PTEN(JK)-ZTMELP2))
               PDPMEL(JK)=ZSNMLT
               PSFL=PSFL-ZSNMLT
               PRFL=PRFL+ZSNMLT
            END IF
         ELSE
            PSFL=PSFL+PDMFUP(JK)+PDMFDP(JK)
         END IF
      END IF
  130 CONTINUE
        PRFL=MAX(PRFL,0.)
        PSFL=MAX(PSFL,0.)
        ZPSUBCL=PRFL+PSFL
      DO 240 JK=KTOPM2,KLEV
      IF(LDCUM.AND.JK.GE.KCBOT.AND. &
             ZPSUBCL.GT.1.E-20) THEN
          ZRFL=ZPSUBCL
          CEVAPCU=CEVAPCU1*SQRT(CEVAPCU2*SQRT(sig1(JK)))
          ZRNEW=(MAX(0.,SQRT(ZRFL/ZCUCOV)-   &
                  CEVAPCU*(PAPH(JK+1)-PAPH(JK))* &
                MAX(0.,PQSEN(JK)-PQEN(JK))))**2*ZCUCOV
          ZRMIN=ZRFL-ZCUCOV*MAX(0.,0.8*PQSEN(JK)-PQEN(JK)) &
               *ZCONS2*(PAPH(JK+1)-PAPH(JK))
          ZRNEW=MAX(ZRNEW,ZRMIN)
          ZRFLN=MAX(ZRNEW,0.)
          ZDRFL=MIN(0.,ZRFLN-ZRFL)
          PDMFUP(JK)=PDMFUP(JK)+ZDRFL
          ZPSUBCL=ZRFLN
      END IF
  240 CONTINUE
        ZDPEVAP=ZPSUBCL-(PRFL+PSFL)
        PRFL=PRFL+ZDPEVAP*PRFL*  &
                  (1./MAX(1.E-20,PRFL+PSFL))
        PSFL=PSFL+ZDPEVAP*PSFL*  &
                  (1./MAX(1.E-20,PRFL+PSFL))
      RETURN
      END SUBROUTINE CUFLX
!

!**********************************************
!       SUBROUTINE CUDTDQ
!********************************************** 
      SUBROUTINE CUDTDQ &
         (KLEV,     KLEVP1,   KTOPM2,   PAPH,             &
          LDCUM,    PTEN,     PTTE,     PQTE,     PMFUS,  &
          PMFDS,    PMFUQ,    PMFDQ,    PMFUL,    PDMFUP, &
          PDMFDP,   ZTMST,    PDPMEL,   PRAIN,    PRFL,   &
          PSFL,     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT, &
          PRSFC,    PSSFC,    PAPRC,    PAPRSM,   PAPRS,  &
          PQEN,     PQSEN,    PLUDE,    PCTE)
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF.  12/89
!***INTERFACE.
!   ----------
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KTOPM2, JK
      REAL      ZTMST, PSRAIN, PSEVAP, PSHEAT, PSMELT, ZDIAGT, ZDIAGW
      REAL      ZALV, RHK, RHCOE, PLDFD, ZDTDT, ZDQDT
      REAL     PTTE(KLEV),        PQTE(KLEV),  &
              PTEN(KLEV),        PLUDE(KLEV),  &
              PGEO(KLEV),        PAPH(KLEVP1), &
              PAPRC,            PAPRS,       &
              PAPRSM,           PCTE(KLEV),   &
              PRSFC,            PSSFC
      REAL     PMFUS(KLEV),       PMFDS(KLEV), &
              PMFUQ(KLEV),       PMFDQ(KLEV), &
              PMFUL(KLEV),       PQSEN(KLEV), &
              PDMFUP(KLEV),      PDMFDP(KLEV),& 
              PRFL,             PRAIN,      &
              PQEN(KLEV)
      REAL     PDPMEL(KLEV),      PSFL
      REAL     ZSHEAT,           ZMELT
      LOGICAL  LDCUM
!--------------------------------
!*    1.0      SPECIFY PARAMETERS
!--------------------------------
  100 CONTINUE
      ZDIAGT=ZTMST
      ZDIAGW=ZDIAGT/RHOH2O
!--------------------------------------------------
!*    2.0      INCREMENTATION OF T AND Q TENDENCIES
!--------------------------------------------------
  200 CONTINUE
      ZMELT=0.
      ZSHEAT=0.
      DO 250 JK=KTOPM2,KLEV
      IF(JK.LT.KLEV) THEN
         IF(LDCUM) THEN
            IF(PTEN(JK).GT.TMELT) THEN
               ZALV=ALV
            ELSE
               ZALV=ALS
            ENDIF
            RHK=MIN(1.0,PQEN(JK)/PQSEN(JK))
            RHCOE=MAX(0.0,(RHK-RHC)/(RHM-RHC))
            pldfd=MAX(0.0,RHCOE*fdbk*PLUDE(JK))
            ZDTDT=(G/(PAPH(JK+1)-PAPH(JK)))*RCPD*      &
              (PMFUS(JK+1)-PMFUS(JK)+                  &
              PMFDS(JK+1)-PMFDS(JK)-ALF*PDPMEL(JK)  &
              -ZALV*(PMFUL(JK+1)-PMFUL(JK)-pldfd-      &
              (PDMFUP(JK)+PDMFDP(JK))))
            PTTE(JK)=PTTE(JK)+ZDTDT
            ZDQDT=(G/(PAPH(JK+1)-PAPH(JK)))*& 
              (PMFUQ(JK+1)-PMFUQ(JK)+       &
              PMFDQ(JK+1)-PMFDQ(JK)+        &
              PMFUL(JK+1)-PMFUL(JK)-pldfd-  &
              (PDMFUP(JK)+PDMFDP(JK)))
            PQTE(JK)=PQTE(JK)+ZDQDT
            PCTE(JK)=(G/(PAPH(JK+1)-PAPH(JK)))*pldfd
            ZSHEAT=ZSHEAT+ZALV*(PDMFUP(JK)+PDMFDP(JK))
            ZMELT=ZMELT+PDPMEL(JK)
         END IF
      ELSE
         IF(LDCUM) THEN
            IF(PTEN(JK).GT.TMELT) THEN
               ZALV=ALV
            ELSE
               ZALV=ALS
            ENDIF
            RHK=MIN(1.0,PQEN(JK)/PQSEN(JK))
            RHCOE=MAX(0.0,(RHK-RHC)/(RHM-RHC))
            pldfd=MAX(0.0,RHCOE*fdbk*PLUDE(JK))
            ZDTDT=-(G/(PAPH(JK+1)-PAPH(JK)))*RCPD*           &
                (PMFUS(JK)+PMFDS(JK)+ALF*PDPMEL(JK)-ZALV* &
                (PMFUL(JK)+PDMFUP(JK)+PDMFDP(JK)+pldfd))  
            PTTE(JK)=PTTE(JK)+ZDTDT
            ZDQDT=-(G/(PAPH(JK+1)-PAPH(JK)))*                &
                     (PMFUQ(JK)+PMFDQ(JK)+pldfd+             &
                     (PMFUL(JK)+PDMFUP(JK)+PDMFDP(JK)))   
            PQTE(JK)=PQTE(JK)+ZDQDT
            PCTE(JK)=(G/(PAPH(JK+1)-PAPH(JK)))*pldfd
            ZSHEAT=ZSHEAT+ZALV*(PDMFUP(JK)+PDMFDP(JK))
            ZMELT=ZMELT+PDPMEL(JK)
         END IF
      END IF
  250 CONTINUE
!---------------------------------------------------------
!      3.      UPDATE SURFACE FIELDS AND DO GLOBAL BUDGETS
!---------------------------------------------------------
  300 CONTINUE
      PRSFC=PRFL
      PSSFC=PSFL
      PAPRC=PAPRC+ZDIAGW*(PRFL+PSFL)
      PAPRS=PAPRSM+ZDIAGW*PSFL
      PSHEAT=PSHEAT+ZSHEAT
      PSRAIN=PSRAIN+PRAIN
      PSEVAP=PSEVAP-(PRFL+PSFL)
      PSMELT=PSMELT+ZMELT
      PSEVAP=PSEVAP+PSRAIN
      RETURN
      END SUBROUTINE CUDTDQ

!
!**********************************************
!       SUBROUTINE CUDUDV
!********************************************** 
      SUBROUTINE CUDUDV &
         (KLEV,     KLEVP1,   KTOPM2,   KTYPE,            &
          KCBOT,    PAPH,     LDCUM,    PUEN,     PVEN,   &
          PVOM,     PVOL,     PUU,      PUD,      PVU,    &
          PVD,      PMFU,     PMFD,     PSDISS)
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF.  12/89
!***INTERFACE.
!   ----------
!          *CUDUDV* IS CALLED FROM *CUMASTR*
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KTOPM2, JK, IK, IKB
      REAL      PSDISS,ZZP, ZDUDT ,ZDVDT
      REAL     PUEN(KLEV),        PVEN(KLEV),   &
              PVOL(KLEV),        PVOM(KLEV),    &
              PAPH(KLEVP1)
      REAL     PUU(KLEV),         PUD(KLEV),    &
              PVU(KLEV),         PVD(KLEV),     &
              PMFU(KLEV),        PMFD(KLEV)
      REAL     ZMFUU(KLEV),       ZMFDU(KLEV),  &
              ZMFUV(KLEV),       ZMFDV(KLEV),   &
              ZDISS
      INTEGER  KTYPE,            KCBOT
      LOGICAL  LDCUM
!------------------------------------------------------------
!*    1.0      CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
! -----------------------------------------------------------
  100 CONTINUE
      DO 120 JK=KTOPM2,KLEV
      IK=JK-1
      IF(LDCUM) THEN
        ZMFUU(JK)=PMFU(JK)*(PUU(JK)-PUEN(IK))
        ZMFUV(JK)=PMFU(JK)*(PVU(JK)-PVEN(IK))
        ZMFDU(JK)=PMFD(JK)*(PUD(JK)-PUEN(IK))
        ZMFDV(JK)=PMFD(JK)*(PVD(JK)-PVEN(IK))
      END IF
  120 CONTINUE
      DO 140 JK=KTOPM2,KLEV
      IF(LDCUM.AND.JK.GT.KCBOT) THEN
         IKB=KCBOT
         ZZP=((PAPH(KLEVP1)-PAPH(JK))/  &
             (PAPH(KLEVP1)-PAPH(IKB)))
         IF(KTYPE.EQ.3) THEN
            ZZP=ZZP**2
         ENDIF
         ZMFUU(JK)=ZMFUU(IKB)*ZZP
         ZMFUV(JK)=ZMFUV(IKB)*ZZP
         ZMFDU(JK)=ZMFDU(IKB)*ZZP
         ZMFDV(JK)=ZMFDV(IKB)*ZZP
      END IF
  140 CONTINUE
      ZDISS=0.
      DO 190 JK=KTOPM2,KLEV
      IF(JK.LT.KLEV) THEN
            IF(LDCUM) THEN
               ZDUDT=(G/(PAPH(JK+1)-PAPH(JK)))* &
                    (ZMFUU(JK+1)-ZMFUU(JK)+     &
                     ZMFDU(JK+1)-ZMFDU(JK))
               ZDVDT=(G/(PAPH(JK+1)-PAPH(JK)))* &
                    (ZMFUV(JK+1)-ZMFUV(JK)+     &
                     ZMFDV(JK+1)-ZMFDV(JK))
               ZDISS=ZDISS+        &
                        PUEN(JK)*(ZMFUU(JK+1)-ZMFUU(JK)+   &
                                     ZMFDU(JK+1)-ZMFDU(JK))+  &
                        PVEN(JK)*(ZMFUV(JK+1)-ZMFUV(JK)+   &
                                     ZMFDV(JK+1)-ZMFDV(JK))
               PVOM(JK)=PVOM(JK)+ZDUDT
               PVOL(JK)=PVOL(JK)+ZDVDT
            END IF
      ELSE
            IF(LDCUM) THEN
               ZDUDT=-(G/(PAPH(JK+1)-PAPH(JK)))* &
                        (ZMFUU(JK)+ZMFDU(JK))
               ZDVDT=-(G/(PAPH(JK+1)-PAPH(JK)))* &
                        (ZMFUV(JK)+ZMFDV(JK))
               ZDISS=ZDISS-        &
      (PUEN(JK)*(ZMFUU(JK)+ZMFDU(JK))+ &
      PVEN(JK)*(ZMFUV(JK)+ZMFDV(JK)))
               PVOM(JK)=PVOM(JK)+ZDUDT
               PVOL(JK)=PVOL(JK)+ZDVDT
            END IF
       END IF
  190 CONTINUE
      PSDISS=PSDISS+ZDISS
      RETURN
      END SUBROUTINE CUDUDV
!

!#################################################################
!
!                 LEVEL 4 SUBROUTINES
!
!#################################################################
!**************************************************************
!             SUBROUTINE CUBASMC
!**************************************************************
      SUBROUTINE CUBASMC   &
         (KLEV,     KLEVM1,  KK,     PTEN,            &
          PQEN,     PQSEN,    PUEN,    PVEN,   PVERV, &
          PGEO,     PGEOH,    LDCUM,   KTYPE,  KLAB,  &
          PMFU,     PMFUB,    PENTR,   KCBOT,  PTU,   &
          PQU,      PLU,      PUU,     PVU,    PMFUS, &
          PMFUQ,    PMFUL,    PDMFUP,  PMFUU,  PMFUV) 
!      M.TIEDTKE         E.C.M.W.F.     12/89
!***PURPOSE.
!   --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!***METHOD.
!   -------
!          S. TIEDTKE (1989)
!***EXTERNALS
!   ---------
!          NONE
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KLEVM1,KK
      REAL      zzzmb
      REAL     PTEN(KLEV),        PQEN(KLEV),  &
              PUEN(KLEV),        PVEN(KLEV),   &
              PQSEN(KLEV),       PVERV(KLEV),  & 
              PGEO(KLEV),        PGEOH(KLEV)
      REAL     PTU(KLEV),         PQU(KLEV),   &
              PUU(KLEV),         PVU(KLEV),    &
              PLU(KLEV),         PMFU(KLEV),   &
              PMFUB,            PENTR,       &
              PMFUS(KLEV),       PMFUQ(KLEV),  &
              PMFUL(KLEV),       PDMFUP(KLEV), &
              PMFUU,            PMFUV
      INTEGER  KTYPE,            KCBOT,      &
              KLAB(KLEV)
      LOGICAL  LDCUM
!--------------------------------------------------------
!*    1.      CALCULATE ENTRAINMENT AND DETRAINMENT RATES
! -------------------------------------------------------
  100 CONTINUE
          IF( .NOT. LDCUM.AND.KLAB(KK+1).EQ.0.0.AND.  &
             PQEN(KK).GT.0.80*PQSEN(KK)) THEN
            PTU(KK+1)=(CPD*PTEN(KK)+PGEO(KK)-PGEOH(KK+1)) &
                               *RCPD
            PQU(KK+1)=PQEN(KK)
            PLU(KK+1)=0.
            ZZZMB=MAX(CMFCMIN,-PVERV(KK)/G)
            ZZZMB=MIN(ZZZMB,CMFCMAX)
            PMFUB=ZZZMB
            PMFU(KK+1)=PMFUB
            PMFUS(KK+1)=PMFUB*(CPD*PTU(KK+1)+PGEOH(KK+1))
            PMFUQ(KK+1)=PMFUB*PQU(KK+1)
            PMFUL(KK+1)=0.
            PDMFUP(KK+1)=0.
            KCBOT=KK
            KLAB(KK+1)=1
            KTYPE=3
            PENTR=ENTRMID
               IF(LMFDUDV) THEN
                  PUU(KK+1)=PUEN(KK)
                  PVU(KK+1)=PVEN(KK)
                  PMFUU=PMFUB*PUU(KK+1)
                  PMFUV=PMFUB*PVU(KK+1)
               END IF
         END IF
      RETURN
      END SUBROUTINE CUBASMC

!
!**************************************************************
!             SUBROUTINE CUADJTQ
!**************************************************************
      SUBROUTINE CUADJTQ(KLEV,KK,PP,PT,PQ,LDFLAG,KCALL)
!      M.TIEDTKE         E.C.M.W.F.     12/89
!      D.SALMOND         CRAY(UK))      12/8/91
!***PURPOSE.
!   --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G.  CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G.  CUDLFS,CUDDRAF
!***EXTERNALS
!   ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV
      INTEGER   KK, KCALL, ISUM
      REAL      ZQSAT, ZCOR, ZCOND1, TT
      REAL     PT(KLEV),          PQ(KLEV),  &
              ZCOND,            ZQP,       &
              PP
      LOGICAL  LDFLAG
!------------------------------------------------------------------
!     2.      CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!------------------------------------------------------------------
  200 CONTINUE
      IF (KCALL.EQ.1 ) THEN
         ISUM=0
         ZCOND=0.
         IF(LDFLAG) THEN
            ZQP=1./PP
            TT=PT(KK)
            ZQSAT=TLUCUA(TT)*ZQP
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
            ZCOND=MAX(ZCOND,0.)
            PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND
            PQ(KK)=PQ(KK)-ZCOND
            IF(ZCOND.NE.0.0) ISUM=ISUM+1
         END IF
         IF(ISUM.EQ.0) GO TO 230
         IF(LDFLAG.AND.ZCOND.NE.0.) THEN
            TT=PT(KK)
            ZQSAT=TLUCUA(TT)*ZQP
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND1=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
            PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND1
            PQ(KK)=PQ(KK)-ZCOND1
         END IF
  230    CONTINUE
      END IF
      IF(KCALL.EQ.2) THEN
         ISUM=0
         ZCOND=0.
         IF(LDFLAG) THEN
            TT=PT(KK)
            ZQP=1./PP
            ZQSAT=TLUCUA(TT)*ZQP
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
            ZCOND=MIN(ZCOND,0.)
            PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND
            PQ(KK)=PQ(KK)-ZCOND
            IF(ZCOND.NE.0.0) ISUM=ISUM+1
         END IF
         IF(ISUM.EQ.0) GO TO 330
         IF(LDFLAG.AND.ZCOND.NE.0.) THEN
            TT=PT(KK)
            ZQSAT=TLUCUA(TT)*ZQP
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND1=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
            PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND1
            PQ(KK)=PQ(KK)-ZCOND1
         END IF
  330    CONTINUE
      END IF
      IF(KCALL.EQ.0) THEN
         ISUM=0
           TT=PT(KK)
           ZQP=1./PP
           ZQSAT=TLUCUA(TT)*ZQP
           ZQSAT=MIN(0.5,ZQSAT)
           ZCOR=1./(1.-VTMPC1*ZQSAT)
           ZQSAT=ZQSAT*ZCOR
           ZCOND=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
           PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND
           PQ(KK)=PQ(KK)-ZCOND
           IF(ZCOND.NE.0.0) ISUM=ISUM+1
         IF(ISUM.EQ.0) GO TO 430
           TT=PT(KK)
           ZQSAT=TLUCUA(TT)*ZQP
           ZQSAT=MIN(0.5,ZQSAT)
           ZCOR=1./(1.-VTMPC1*ZQSAT)
           ZQSAT=ZQSAT*ZCOR
           ZCOND1=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
           PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND1
           PQ(KK)=PQ(KK)-ZCOND1
  430    CONTINUE
      END IF
      IF(KCALL.EQ.4) THEN
           TT=PT(KK)
           ZQP=1./PP
           ZQSAT=TLUCUA(TT)*ZQP
           ZQSAT=MIN(0.5,ZQSAT)
           ZCOR=1./(1.-VTMPC1*ZQSAT)
           ZQSAT=ZQSAT*ZCOR
           ZCOND=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
           PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND
           PQ(KK)=PQ(KK)-ZCOND
           TT=PT(KK)
           ZQSAT=TLUCUA(TT)*ZQP
           ZQSAT=MIN(0.5,ZQSAT)
           ZCOR=1./(1.-VTMPC1*ZQSAT)
           ZQSAT=ZQSAT*ZCOR
           ZCOND1=(PQ(KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(TT))
           PT(KK)=PT(KK)+TLUCUC(TT)*ZCOND1
           PQ(KK)=PQ(KK)-ZCOND1
      END IF
      RETURN
      END SUBROUTINE CUADJTQ

!
!**********************************************************
!        SUBROUTINE CUENTR_NEW
!**********************************************************
      SUBROUTINE CUENTR_NEW                              &   
         (KLEV,     KLEVP1,   KK,       PTENH,           &
          PAPH,     PAP,      PGEOH,    KLWMIN,   LDCUM, &
          KTYPE,    KCBOT,    KCTOP0,   ZPBASE,   PMFU,  &
          PENTR,    ZDMFEN,   ZDMFDE,   ZODETR,   KHMIN)
!      M.TIEDTKE         E.C.M.W.F.     12/89
!      Y.WANG            IPRC           11/01
!***PURPOSE.
!   --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!***INTERFACE
!   ---------
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!***METHOD.
!  --------
!          S. TIEDTKE (1989), NORDENG(1996)
!***EXTERNALS
!   ---------
!          NONE
! ----------------------------------------------------------------
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER   KLEV, KLEVP1
      INTEGER   KK, IKLWMIN,IKB, IKT, IKH
      REAL      ZRRHO, ZDPRHO, ZPMID, ZENTR, ZZMZK, ZTMZK, ARG, ZORGDE
      REAL     PTENH(KLEV),                           &
              PAP(KLEV),         PAPH(KLEVP1),   &
              PMFU(KLEV),        PGEOH(KLEV),    &
              PENTR,            ZPBASE,        &
              ZDMFEN,           ZDMFDE,        &
              ZODETR(KLEV)
      INTEGER  KLWMIN,           KTYPE,        &
              KCBOT,            KCTOP0,        &
              KHMIN
      LOGICAL  LDCUM,LLO1,LLO2

      real    tt,ttb
      real    zqsat, zqsatb
      real    fscale
!---------------------------------------------------------
!*    1.       CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!---------------------------------------------------------
!*    1.1      SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!----------------------------------------------------------
!*    1.2      SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!-------------------------------------------------------
        zpbase = paph(kcbot)
        zrrho = (rd*ptenh(kk+1))/paph(kk+1)
        zdprho = (paph(kk+1)-paph(kk))*zrg
! old or new choice
        zpmid = 0.5*(zpbase+paph(kctop0))
        zentr = pentr*pmfu(kk+1)*zdprho*zrrho
        llo1 = kk.LT.kcbot.AND.ldcum
! old or new choice
        if(llo1) then
           if(nturben.eq.1) zdmfde = zentr
           if(nturben.eq.2) zdmfde = zentr*1.2
        else
          zdmfde = 0.0
        endif
! old or new choice
        if(nturben .eq. 1) then
          fscale = 1.0
        elseif (nturben .eq. 2) then
! defining the facale
        tt = ptenh(kk+1)
        zqsat = TLUCUA(tt)/paph(kk+1)
        zqsat = zqsat/(1.-VTMPC1*zqsat)
        ttb = ptenh(kcbot)
        zqsatb = TLUCUA(ttb)/zpbase
        zqsatb = zqsatb/(1.-VTMPC1*zqsatb)
        fscale = 4.0*(zqsat/zqsatb)**2
        end if
! end of defining the fscale
        llo2 = llo1.AND.ktype.EQ.2.AND.((zpbase-paph(kk)) &
             .LT.ZDNOPRC.OR.paph(kk).GT.zpmid)
        if(llo2) then
            zdmfen = zentr*fscale
        else
            zdmfen = 0.0
        endif
        iklwmin = MAX(klwmin,kctop0+2)
        llo2 = llo1.AND.ktype.EQ.3.AND.(kk.GE.iklwmin.OR.pap(kk) &
             .GT.zpmid)
        IF (llo2) zdmfen = zentr*fscale
        llo2 = llo1.AND.ktype.EQ.1
! Turbulent entrainment
        IF (llo2) zdmfen = zentr*fscale
! Organized detrainment, detrainment starts at khmin
        ikb = kcbot
        zodetr(kk) = 0.
        IF (llo2.AND.kk.LE.khmin.AND.kk.GE.kctop0) THEN
          ikt = kctop0
          ikh = khmin
          IF (ikh.GT.ikt) THEN
            zzmzk = -(pgeoh(ikh)-pgeoh(kk))*zrg
            ztmzk = -(pgeoh(ikh)-pgeoh(ikt))*zrg
            arg = 3.1415*(zzmzk/ztmzk)*0.5
            zorgde = TAN(arg)*3.1415*0.5/ztmzk
            zdprho = (paph(kk+1)-paph(kk))*(zrg*zrrho)
            zodetr(kk) = MIN(zorgde,1.E-3)*pmfu(kk+1)*zdprho
          END IF
        END IF
!
      RETURN
      END SUBROUTINE CUENTR_NEW

      REAL FUNCTION TLUCUA(TT)
!
!  Set up lookup tables for cloud ascent calculations.
!
      IMPLICIT NONE
      REAL ZCVM3,ZCVM4,TT
!
      IF(TT-TMELT.GT.0.) THEN
         ZCVM3=C3LES
         ZCVM4=C4LES
      ELSE
         ZCVM3=C3IES
         ZCVM4=C4IES
      END IF
      TLUCUA=C2ES*EXP(ZCVM3*(TT-TMELT)*(1./(TT-ZCVM4)))
!
      RETURN
      END FUNCTION TLUCUA
!
      REAL FUNCTION TLUCUB(TT)
!
!  Set up lookup tables for cloud ascent calculations.
!
      IMPLICIT NONE
      REAL Z5ALVCP,Z5ALSCP,ZCVM4,ZCVM5,TT
!
      Z5ALVCP=C5LES*ALV/CPD
      Z5ALSCP=C5IES*ALS/CPD
      IF(TT-TMELT.GT.0.) THEN
         ZCVM4=C4LES
         ZCVM5=Z5ALVCP
      ELSE
         ZCVM4=C4IES
         ZCVM5=Z5ALSCP
      END IF
      TLUCUB=ZCVM5*(1./(TT-ZCVM4))**2
!
      RETURN
      END FUNCTION TLUCUB
!
      REAL FUNCTION TLUCUC(TT)
!
!  Set up lookup tables for cloud ascent calculations.
!
      IMPLICIT NONE
      REAL ZALVDCP,ZALSDCP,TT,ZLDCP
!
      ZALVDCP=ALV/CPD
      ZALSDCP=ALS/CPD
      IF(TT-TMELT.GT.0.) THEN
         ZLDCP=ZALVDCP
      ELSE
         ZLDCP=ZALSDCP
      END IF
      TLUCUC=ZLDCP
!
      RETURN
      END FUNCTION TLUCUC
!

END MODULE module_cu_tiedtke
