MODULE module_cu_emanuel

!**************************************************************************
!****                       SUBROUTINE CONVECT                        *****
!****                          VERSION 4.3c                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************

  use mem_para,    only: myrank
  use misc_coms,   only: io6
  use mem_grid,    only: mza
  use consts_coms, only: r8

  implicit none

  private :: myrank, io6, mza, r8

  real,    parameter, private :: xmbmax = 1.0
  integer, parameter, private :: iupsrc =   0

CONTAINS

SUBROUTINE cuparm_emanuel(iw, dtlong)

  use mem_grid,    only: lpw, zm, zt, dzt, lsw
  use mem_basic,   only: tair, press, rho, rr_v, ue, ve
  use consts_coms, only: t00, grav
  use mem_cuparm , only: thsrc, rtsrc, conprr, cbmf, kcutop, kcubot, kudbot, &
                         qwcon, iactcu, kddbot, kddtop, kddmax, cddf, &
                         umsrc, vmsrc
  use mem_turb,    only: frac_sfc
  use therm_lib,   only: rhovsil

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: dtlong

  real, dimension(mza) :: tc, qc, qsc, u, v, pc, pfc, gz, den, dz
  real, dimension(mza) :: tt, qt, ut, vt, qcldc, mp, qflux

  integer :: k, ka, kc, kp, ks, nd, na, nm, kmax
  real    :: pcprate, wprime, tprime, qprime
  integer :: iflag, kcbase, kctop, kup

  ka = lpw(iw)

  ! Set initial convective top to mza-1

  nd = mza - ka

  ! Redefine top to go no higher than 50mb for convective calculations
  ! to prevent any problems when esat gets near ambient pressure

  do k = ka, nd-1
     if (press(k,iw) < 50.e2) then
        nd = k - ka + 1
        exit
     endif
  enddo

  na = nd + 1
  nm = nd - 1

  ! Temp, water vapor, geopotential, and pressure

  do kc = 1, nd
     k  = kc + ka - 1
     kp = kc + 1

     gz (kc) = grav * (zt(k) - zt(ka)) ! geopotential height relative to 1st level
     tc (kc) = tair(k,iw)
     qc (kc) = max(rr_v(k,iw), 1.e-10)
     pc (kc) = 0.01  *  press(k,iw)
     pfc(kp) = 0.005 * (press(k,iw) + press(k+1,iw))
     den(kc) = rho(k,iw)
     dz (kc)  = dzt(k)
  enddo

  ! Estimate surface pressure

  pfc(1) =  0.01 * (press(ka,iw) + (zt(ka)-zm(ka-1))*rho(ka,iw)*grav)

  do kc = 1, nd
     k  = kc + ka - 1
     u(kc) = ue(k,iw)
     v(kc) = ve(k,iw)
  enddo

  ! Saturated vapor pressure

  do kc = 1, nd
     qsc(kc) = rhovsil(tc(kc)-t00) / den(kc)
     qc (kc) = min(qsc(kc), qc(kc))
  enddo

  call convect43c (     iw,     den,    dz,                                 &
       tc,      qc,     qsc,    u,      v,        pc,    pfc, gz,           &
       nd,      nm,     dtlong, iflag,  tt,       qt,    ut,  vt, mp,       &
       pcprate, wprime, tprime, qprime, cbmf(iw), qcldc, qflux, kup, kcbase, kctop )

  if (iflag == 1 .or. iflag == 4) then

     kcutop(iw) = kctop  + ka - 1
     kcubot(iw) = kcbase + ka - 1
     kudbot(iw) = kup    + ka - 1
     iactcu(iw) = 1

     kmax = maxloc(mp(1:kctop),dim=1)
     if (mp(kmax) > 1.e-6) then

        kddmax(iw) = kmax + ka - 1
        cddf  (iw) = mp(kmax)

        do kc = kctop, kmax+1, -1
           if (mp(kc) > 1.e-7) exit
        enddo
        kddtop(iw) = kc + ka - 1

        do kc = 1, kmax-1
           if (mp(kc) > 1.e-7) exit
        enddo
        kddbot(iw) = max(kc + ka - 2, ka)
     endif

     do kc = 1, kctop
        k  = kc + ka - 1
        thsrc(k,iw) = tt(kc) * den(kc)
        rtsrc(k,iw) = qt(kc) * den(kc)
        qwcon(k,iw) = qcldc(kc)
        umsrc(k,iw) = ut(kc) * den(kc)
        vmsrc(k,iw) = vt(kc) * den(kc)
     enddo

     ! Surface precipitation
     do ks = 1, min(lsw(iw), kctop)
        conprr(iw) = conprr(iw) + frac_sfc(ks,iw) * qflux(kc)
     enddo

  endif

END SUBROUTINE cuparm_emanuel


SUBROUTINE CONVECT43C (  iw,     rho,  dz,                  &
     T,      Q,  QS,     U,      V,    P,      PH, gz,      &
     ND,     NL, DELT,   IFLAG,  FT,   FQ,     FU, FV,  mp, &
     PRECIP, WD, TPRIME, QPRIME, CBMF, QCONDC, qflux, nk, icb, inb )

!-----------------------------------------------------------------------------
!    *** On input:      ***
!
!     T:   Array of absolute temperature (K) of dimension ND, with first
!           index corresponding to lowest model level. Note that this array
!           will be altered by the subroutine if dry convective adjustment
!           occurs and if IPBL is not equal to 0.
!
!     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     QS:  Array of saturation specific humidity of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
!            index corresponding with the lowest model level. Defined at
!            same levels as T. Note that this array will be altered if
!            dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     V:   Same as U but for meridional velocity.
!
!     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
!            where NTRA is the number of different tracers. If no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (ND,1). Tracers are defined at
!            same vertical levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     P:   Array of pressure (mb) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T.
!
!     PH:  Array of pressure (mb) of dimension ND+1, with first index
!            corresponding to lowest level. These pressures are defined at
!            levels intermediate between those of P, T, Q and QS. The first
!            value of PH should be greater than (i.e. at a lower level than)
!            the first value of the array P.
!
!     ND:  The dimension of the arrays T, Q, QS, P, PH, FT and FQ
!
!     NL:  The maximum number of levels to which convection can
!            penetrate, plus 1.
!            NL MUST be less than or equal to ND-1.
!
!     NA:  The parameter na should in general be greater than
!            or equal to ND + 1
!
!     NTRA:The number of different tracers. If no tracer transport
!            is needed, set this equal to 1. (On most compilers, setting
!            NTRA to 0 will bypass tracer calculation, saving some CPU.)
!
!     DELT: The model time step (sec) between calls to CONVECT
!
! -- sb: interface with the cloud parameterization:
!
!     QCONDC: mixing ratio of condensed water within clouds (kg/kg)
!               For use in the Bony-Emanuel cloud parameterization
! sb --
!----------------------------------------------------------------------------
!    ***   On Output:         ***
!
!     IFLAG: An output integer whose value denotes the following:
!
!                VALUE                        INTERPRETATION
!                -----                        --------------
!                  0               No moist convection; atmosphere is not
!                                  unstable, or surface temperature is less
!                                  than 250 K or surface specific humidity
!                                  is non-positive.
!
!                  1               Moist convection occurs.
!
!                  2               No moist convection: lifted condensation
!                                  level is above the 200 mb level.
!
!                  3               No moist convection: cloud base is higher
!                                  then the level NL-1.
!
!                  4               Moist convection occurs, but a CFL condition
!                                  on the subsidence warming is violated. This
!                                  does not cause the scheme to terminate.
!
!     FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
!             grid levels as T, Q, QS and P.
!
!     FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!             defined at same grid levels as T, Q, QS and P.
!
!     FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!             defined at same grid levels as T.
!
!     FV:   Same as FU, but for forcing of meridional velocity.
!
!     FTRA: Array of forcing of tracer content, in tracer mixing ratio per
!             second, defined at same levels as T. Dimensioned (ND,NTRA).
!
!     PRECIP: Scalar convective precipitation rate (mm/s).
!
!     WD:    A convective downdraft velocity scale. For use in surface
!             flux parameterizations. See convect.ps file for details.
!
!     TPRIME: A convective downdraft temperature perturbation scale (K).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!     QPRIME: A convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!              by the calling program between calls to CONVECT.
!
!------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: nd, nl, iw !, ntra
  real,    intent(in) :: t(nd), q(nd), qs(nd), u(nd), v(nd)
  real,    intent(in) :: p(nd), ph(nd+1), gz(nd), delt !, tra(nd,ntra)
  real,    intent(in) :: rho(nd), dz(nd)

  integer, intent(out) :: iflag, icb, inb, nk
  real,    intent(out) :: ft(nd), fq(nd), fu(nd), fv(nd)
  real,    intent(out) :: qcondc(nd), qflux(nd) !, ftra(nd,ntra)
  real,    intent(out) :: wd, tprime, qprime, precip
  real,    intent(inout) :: cbmf
  real,    intent(out) :: mp(nd)

  real :: ad, afac, ahmax, ahmin, alt, altem, am, amde, amp1, anum, asij, &
       awat, b6, bf2, bsum, by, c6, cape, capem, chi, coeff, &
       cpinv, cwat, dbo, dbosum, defrac, dei, delm, delp, &
       delti, denom, dhdp, dpinv, dtma, dtmin, dtpbl, elacrit, epmax, &
       fac, fqold, frac, ftold, fuold, fvold, plcl, qp1, &
       qsm, qstm, qti, rat, revap, rh, scrit, siga, sigt, sjmax, sjmin, smid, &
       smin, stemp, tca, tvaplcl, tvpplcl, wdtrain

  integer :: i, ihmin, j, jtt, k
! integer :: inb1
  integer ::  nent(mza)
  real    ::  uent(mza,mza), vent(mza,mza) !, traent(mza,mza,ntra), tratm(mza)
  real    ::  up(mza), vp(mza) !, trap(mza,ntra)
  real    ::  m(mza), ment(mza,mza), qent(mza,mza), elij(mza,mza)
  real    ::  sij(mza,mza), tvp(mza), tv(mza), water(mza)
  real    ::  qp(mza), ep(mza), wt(mza), evap(mza), clw(mza)
  real    ::  sigp(mza), tp(mza), cpn(mza)
  real    ::  lv(mza), lvcp(mza), h(mza), hp(mza), hm(mza)
  real    ::  qcond(mza), nqcond(mza), wa(mza), ma(mza), ax(mza)

  integer :: nkmax
  real    :: dbmax, deltv
  real(r8):: qsum, tsum, usum, vsum, uav, vav, qav, tav

! -----------------------------------------------------------------------
!
!   ***                     Specify Switches                         ***
!
!   ***   MINORIG: Lowest level from which convection may originate  ***
!   ***     (Should be first model level at which T is defined       ***
!   ***      for models using bulk PBL schemes; otherwise, it should ***
!   ***      be the first model level at which T is defined above    ***
!   ***                      the surface layer)                      ***

  integer, parameter :: MINORIG=3

!------------------------------------------------------------------------------
!
!   ***                    SPECIFY PARAMETERS                        ***
!
!   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!   ***               BETWEEN 0 C AND TLCRIT)                        ***
!   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!   ***                       FORMULATION                            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!   ***                        OF CLOUD                              ***
!   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF RAIN                             ***
!   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF SNOW                             ***
!   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!   ***                         TRANSPORT                            ***
!   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!   ***                   (DAMP MUST BE LESS THAN 1)                 ***

  real, parameter :: &
       ELCRIT=.0011, &
       TLCRIT=-55.0, &
       ENTP=1.0,     & ! 1.5
       SIGD=0.05,    &
       SIGS=0.12,    &
       OMTRAIN=50.0, &
       OMTSNOW=5.5,  &
       COEFFR=1.0,   &
       COEFFS=0.8,   &
       CU=0.7,       &
       BETA=10.0,    &
       DTMAX=0.5,    & ! 0.9
       ALPHA=0.2,    &
       DAMP=0.1,     &
       DELTA=0.01    ! sb (for cloud parameterization)

!   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
!   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
!   ***             THESE SHOULD BE CONSISTENT WITH             ***
!   ***              THOSE USED IN CALLING PROGRAM              ***
!   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***

  real, parameter :: &
       CPD=1005.7,   &
       CPV=1870.0,   &
       CL=2500.0,    &
       RV=461.5,     &
       RD=287.04,    &
       LV0=2.501E6,  &
       G=9.8,        &
       ROWL=1000.0,  &
       CPVMCL=CL-CPV,&
       EPS=RD/RV,    &
       EPSI=1./EPS,  &
       GINV=1.0/G

  DELTI=1.0/DELT

!           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***

  DO I=1,ND
     FT(I)=0.0
     FQ(I)=0.0
     fu(i)=0.0
     fv(i)=0.0

     ! -- sb:
     QCONDC(I)=0.0
     QCOND(I)=0.0
     NQCOND(I)=0.0
     ma(i) = 0.0
     qflux(i) = 0.

!     DO J=1,NTRA
!        FTRA(I,J)=0.0
!     ENDDO
  ENDDO

  PRECIP=0.0
  WD=0.0
  TPRIME=0.0
  QPRIME=0.0
  IFLAG=0

!  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY

!  GZ(1)=0.0
!  CPN(1)=CPD*(1.-Q(1))+Q(1)*CPV
!  H(1)=T(1)*CPN(1)
!  LV(1)=LV0-CPVMCL*(T(1)-273.15)
!  HM(1)=LV(1)*Q(1)
!  TV(1)=T(1)*(1.+Q(1)*EPSI-Q(1))

  DO I=1,NL+1
!    TVX=T(I)*(1.+Q(I)*EPSI-Q(I))
!    TVY=T(I-1)*(1.+Q(I-1)*EPSI-Q(I-1))
!    GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(P(I-1)-P(I))/PH(I)
     CPN(I)=CPD*(1.-Q(I))+CPV*Q(I)
     H(I)=T(I)*CPN(I)+GZ(I)
     LV(I)=LV0-CPVMCL*(T(I)-273.15)
     HM(I)=(CPD*(1.-Q(I))+CL*Q(I))*(T(I)-T(1))+LV(I)*Q(I)+GZ(I)
     TV(I)=T(I)*(1.+Q(I)*EPSI-Q(I))
  enddo

  ! Find level of minimum moist static energy

  AHMIN=1.0E12
  IHMIN=NL-1
  do i = minorig, nl-2
     IF(HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
        AHMIN=HM(I)
        IHMIN=I
     END IF
  ENDDO

!  ***     Find that model level below the level of minimum moist       ***
!  ***  static energy that has the maximum value of moist static energy ***

  AHMAX=0.0
  nk = minorig
  DO I=MINORIG,IHMIN
     IF(HM(I).GT.AHMAX)THEN
        NK=I
        AHMAX=HM(I)
     END IF
  ENDDO

!  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
!  ***                          ARE REASONABLE                         ***
!  ***      Skip convection if HM increases monotonically upward       ***

  IF(T(NK).LT.250.0.OR.Q(NK).LE.0.0.OR.IHMIN.EQ.(NL-1))THEN
     IFLAG=0
     CBMF=0.0
     RETURN
  END IF

  if (nk < ihmin .and. iupsrc == 1) then
     nkmax = nk
     dbmax = -1.e10

     do i = nk, ihmin
        RH=Q(I)/QS(I)
        CHI=T(I)/(1669.0-122.0*RH-T(I))
        PLCL=P(I)*(RH**CHI)
        IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0) cycle

        icb = nl-1
        do k = i+1, nl
           if (p(k) < plcl) then
              icb = k
              exit
           endif
        enddo
        IF (ICB >= NL-1) cycle

        call tlift(p,t,q,qs,gz,icb,i,tvp,tp,clw,nd,nl,1)
        deltv = tvp(icb) - tp(icb)*q(i) - tv(icb)

        if (deltv > dbmax) then
           dbmax = deltv
           nkmax = i
        endif
     enddo

     nk = nkmax
  endif

!   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
!   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***

  RH=Q(NK)/QS(NK)
  CHI=T(NK)/(1669.0-122.0*RH-T(NK))
  PLCL=P(NK)*(RH**CHI)
  IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
     IFLAG=2
     CBMF=0.0
     RETURN
  END IF

!   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***

  ICB=NL-1
  DO I=NK+1,NL-2
     IF(P(I).LT.PLCL)THEN
        ICB=I
        EXIT
     END IF
  ENDDO
  IF(ICB.GE.(NL-1))THEN
     IFLAG=3
     CBMF=0.0
     RETURN
  END IF

!   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***

!   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
!   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
!   ***                   LIQUID WATER CONTENT                             ***

  CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,1)
  DO I=NK,ICB
     TVP(I)=TVP(I)-TP(I)*Q(NK)
  ENDDO

!   ***  If there was no convection at last time step and parcel    ***
!   ***       is stable at ICB then skip rest of calculation        ***

  IF (CBMF < 1.e-6 .AND. TVP(ICB).LE.(TV(ICB)-DTMAX)) THEN
     cbmf  = 0.0
     IFLAG = 0
     RETURN
  END IF

!   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***

  IF(IFLAG.NE.4)IFLAG=1

!   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***

  CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,2)

  if (clw(icb) <= 0.0) then

     do i = icb+1, nk
        if (clw(i) > 1.e-30) then
           icb = i
           exit
        endif
     enddo

     if (i >= nk) then
        iflag = 0
        cbmf  = 0.0
        return
     end if

  endif

!   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
!   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
!   ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***

  DO I=1,NK
     EP(I)=0.0
     SIGP(I)=SIGS
  ENDDO
  DO I=NK+1,NL
     TCA=TP(I)-273.15
     IF(TCA.GE.0.0)THEN
        ELACRIT=ELCRIT
     ELSE
        ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
     END IF
     ELACRIT=MAX(ELACRIT,0.0)
     EPMAX=0.999
     EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
     EP(I)=MAX(EP(I),0.0)
     EP(I)=MIN(EP(I),EPMAX)
     SIGP(I)=SIGS
  ENDDO

!   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
!   ***                    VIRTUAL TEMPERATURE                    ***

  DO I=ICB+1,NL
     TVP(I)=TVP(I)-TP(I)*Q(NK)
  ENDDO
  TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD

!   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***

  DO I=1,NL+1
     HP(I)=H(I)
     NENT(I)=0
     WATER(I)=0.0
     EVAP(I)=0.0
     WT(I)=OMTSNOW
     MP(I)=0.0
     M(I)=0.0
     LVCP(I)=LV(I)/CPN(I)
     DO J=1,NL+1
        QENT(I,J)=Q(J)
        ELIJ(I,J)=0.0
        MENT(I,J)=0.0
        SIJ(I,J)=0.0
        uent(i,j)=u(j)
        vent(i,j)=v(j)
!        DO K=1,NTRA
!           TRAENT(I,J,K)=TRA(J,K)
!        ENDDO
     ENDDO
  ENDDO
  QP(1)=Q(1)
  up(1)=u(1)
  vp(1)=v(1)
!  DO I=1,NTRA
!     TRAP(1,I)=TRA(1,I)
!  ENDDO
  DO I=2,NL+1
     QP(I)=Q(I-1)
     up(i)=u(i-1)
     vp(i)=v(i-1)
!     DO J=1,NTRA
!        TRAP(I,J)=TRA(I-1,J)
!     ENDDO
  ENDDO

!  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
!  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
!  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***

  CAPE=0.0
  INB=ICB+1
! INB1=INB

  DO I=ICB+1,NL-1
     CAPEM=CAPE
     BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)
     CAPE=CAPE+BY
!    IF(BY.GE.0.0)INB1=I+1
     IF(CAPE.GT.0.0)THEN
        INB=I+1
     ELSE
        EXIT
     ENDIF
  ENDDO

! inb1=min(inb1,inb)
  DEFRAC=CAPEM-CAPE
  DEFRAC=MAX(DEFRAC,0.001)
  FRAC=-CAPE/DEFRAC
! FRAC=MIN(FRAC,1.0)
  FRAC=MAX(FRAC,0.0)

  if (frac >= 1.0) then
     inb  = inb - 1
     frac = 00
  endif

  ! require convection to span at least 2 layers
  IF (inb == icb) THEN
     iflag = 0
     cbmf  = 0.0
     RETURN
  END IF

!   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***

  DO I=ICB,INB
     HP(I)=H(NK)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)
  ENDDO

!   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
!   ***                   AT EACH MODEL LEVEL                       ***

  DBOSUM=0.0

!   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
!   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
!	
  TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL) / (CPN(ICB-1)*P(ICB-1))
  TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB)) / (P(ICB)-P(ICB+1))
  DTPBL=0.0
  DO I=NK,ICB-1
     DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))
  ENDDO
  DTPBL=DTPBL/(PH(NK)-PH(ICB))
  DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
  DTMA=DTMIN

!   ***  ADJUST CLOUD BASE MASS FLUX   ***

! CBMFOLD=CBMF
! DELT0=300.0
! DELT0=600.0
! DAMPS=DAMP*DELT/DELT0
! CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA
  CBMF = (1.-DAMP)*CBMF + 0.1*ALPHA*DTMA
  CBMF=MAX(CBMF,0.0)
  CBMF=MIN(CBMF,xmbmax)
!
! cbmf=min(cbmf,0.95*delti*ginv*(ph(minorig)-ph(minorig+1))/(0.01*2.))

!   *** If cloud base mass flux is zero, skip rest of calculation  ***

  IF (CBMF < 1.e-6) THEN
     iflag = 0
     cbmf  = 0.0
     RETURN
  END IF

!   ***   CALCULATE RATES OF MIXING,  M(I)   ***

! M(ICB)=0.0
! DO I=ICB+1,INB
  DO I=ICB,INB
!    K=MIN(I,INB1)
!    DBO=ABS(TV(K)-TVP(K))+ENTP*0.02*(PH(K)-PH(K+1))
     K=I
     DBO=max(TVP(K)-TV(K),0.0)+ENTP*0.02*(PH(K)-PH(K+1))
     DBOSUM=DBOSUM+DBO
     M(I)=CBMF*DBO
  ENDDO
! DO I=ICB+1,INB
  DO I=ICB,INB
     M(I)=M(I)/DBOSUM
  ENDDO

!   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
!   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
!   ***                        FRACTION (SIJ)                          ***

! DO I=ICB+1,INB
  DO I=ICB,INB
     QTI=Q(NK)-EP(I)*CLW(I)
     DO J=ICB,INB
        BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
        ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
        DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
        DEI=DENOM
        IF(ABS(DEI).LT.0.01)DEI=0.01
        SIJ(I,J)=ANUM/DEI
        SIJ(I,I)=1.0
        ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
        ALTEM=ALTEM/BF2
        CWAT=CLW(J)*(1.-EP(J))
        STEMP=SIJ(I,J)
        IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.ALTEM.GT.CWAT).AND.J.GT.I)THEN
           ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
           DENOM=DENOM+LV(J)*(Q(I)-QTI)
           IF(ABS(DENOM).LT.0.01)DENOM=0.01
           SIJ(I,J)=ANUM/DENOM
           ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
           ALTEM=ALTEM-(BF2-1.)*CWAT
        END IF
        IF(SIJ(I,J).GE.0.0.AND.SIJ(I,J).LE.0.9)THEN
           QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI
           uent(i,j)=sij(i,j)*u(i)+(1.-sij(i,j))*u(nk)
           vent(i,j)=sij(i,j)*v(i)+(1.-sij(i,j))*v(nk)
!           DO K=1,NTRA
!              TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1.-SIJ(I,J))*TRA(NK,K)
!           END DO
           ELIJ(I,J)=ALTEM
           ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
           MENT(I,J)=M(I)/(1.-SIJ(I,J))
           NENT(I)=NENT(I)+1
        END IF
        SIJ(I,J)=MAX(0.0,SIJ(I,J))
        SIJ(I,J)=MIN(1.0,SIJ(I,J))
     ENDDO

     ! IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
     ! AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***

     IF(NENT(I).EQ.0)THEN
        MENT(I,I)=M(I)
        QENT(I,I)=Q(NK)-EP(I)*CLW(I)
        uent(i,i)=u(nk)
        vent(i,i)=v(nk)
!        DO J=1,NTRA
!           TRAENT(I,I,J)=TRA(NK,J)
!        END DO
        ELIJ(I,I)=CLW(I)
        SIJ(I,I)=1.0
     END IF
  ENDDO
  SIJ(INB,INB)=1.0

!   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
!   ***              PROBABILITIES OF MIXING                     ***

! DO I=ICB+1,INB
  DO I=ICB,INB
     IF(NENT(I).NE.0)THEN
        QP1=Q(NK)-EP(I)*CLW(I)
        ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
        DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
        IF(ABS(DENOM).LT.0.01)DENOM=0.01
        SCRIT=ANUM/DENOM
        ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
        IF(ALT.LT.0.0)SCRIT=1.0
        SCRIT=MAX(SCRIT,0.0)
        ASIJ=0.0
        SMIN=1.0
        DO J=ICB,INB
           IF(SIJ(I,J).GE.0.0.AND.SIJ(I,J).LE.0.9)THEN
              IF(J.GT.I)THEN
                 SMID=MIN(SIJ(I,J),SCRIT)
                 SJMAX=SMID
                 SJMIN=SMID
                 IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
                    SMIN=SMID
                    SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
                    SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
                    SJMIN=MIN(SJMIN,SCRIT)
                 END IF
              ELSE
                 SJMAX=MAX(SIJ(I,J+1),SCRIT)
                 SMID=MAX(SIJ(I,J),SCRIT)
                 SJMIN=0.0
                 IF(J.GT.1)SJMIN=SIJ(I,J-1)
                 SJMIN=MAX(SJMIN,SCRIT)
              END IF
              DELP=ABS(SJMAX-SMID)
              DELM=ABS(SJMIN-SMID)
              ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))
              MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))
           END IF
        ENDDO
        ASIJ=MAX(1.0E-21,ASIJ)
        ASIJ=1.0/ASIJ
        DO J=ICB,INB
           MENT(I,J)=MENT(I,J)*ASIJ
        ENDDO
        BSUM=0.0
        DO J=ICB,INB
           BSUM=BSUM+MENT(I,J)
        ENDDO
        IF(BSUM.LT.1.0E-18)THEN
           NENT(I)=0
           MENT(I,I)=M(I)
           QENT(I,I)=Q(NK)-EP(I)*CLW(I)
           uent(i,i)=u(nk)
           vent(i,i)=v(nk)
!           DO J=1,NTRA
!              TRAENT(I,I,J)=TRA(NK,J)
!           END DO
           ELIJ(I,I)=CLW(I)
           SIJ(I,I)=1.0
        END IF
     END IF
  ENDDO

!   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
!   ***             DOWNDRAFT CALCULATION                      ***

  IF (EP(INB) >= 0.0001) THEN

     !    ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
     !    ***                AND CONDENSED WATER FLUX                    ***

     JTT=2

     !    ***                    BEGIN DOWNDRAFT LOOP                    ***

     DO I=INB,1,-1

        ! ***              CALCULATE DETRAINED PRECIPITATION             ***

        WDTRAIN=G*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
           DO J=1,I-1
              AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
              AWAT=MAX(0.0,AWAT)
              WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
           ENDDO
        END IF

        ! ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
        ! ***              ESTIMATES OF QP(I)AND QP(I-1)             ***

        ! *** Value of terminal velocity and coefficient of evaporation for snow

        COEFF=COEFFS
        WT(I)=OMTSNOW

        ! ***  Value of terminal velocity and coefficient of evaporation for rain

        IF(T(I).GT.273.0)THEN
           COEFF=COEFFR
           WT(I)=OMTRAIN
        END IF
        QSM=0.5*(Q(I)+QP(I+1))
        AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4+2.0E3*PH(I)*QS(I))
        AFAC=MAX(AFAC,0.0)
        SIGT=SIGP(I)
        SIGT=MAX(0.0,SIGT)
        SIGT=MIN(1.0,SIGT)
        B6=100.*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
        REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
        EVAP(I)=SIGT*AFAC*REVAP
        WATER(I)=REVAP*REVAP

        ! ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
        ! ***              HYDROSTATIC APPROXIMATION                 ***

        IF (I .NE. 1) THEN
           DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))
           DHDP=MAX(DHDP,10.0)
           MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
           MP(I)=MAX(MP(I),0.0)

           ! ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***

           FAC=20.0/(PH(I-1)-PH(I))
           MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)

           ! ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
           ! ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***

           IF(P(I).GT.(0.949*P(1)))THEN
              JTT=MAX(JTT,I)
              MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))
           END IF
        ENDIF

        ! ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***

        IF (I .NE. INB) THEN
           IF(I.EQ.1)THEN
              QSTM=QS(1)
           ELSE
              QSTM=QS(I-1)
           END IF
           IF(MP(I).GT.MP(I+1))THEN
              RAT=MP(I+1)/MP(I)
              QP(I)=QP(I+1)*RAT+Q(I)*(1.0-RAT)+100.*GINV* &
                   SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))
              up(i)=up(i+1)*rat+u(i)*(1.-rat)
              vp(i)=vp(i+1)*rat+v(i)*(1.-rat)
!              DO J=1,NTRA
!               ! TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1.-RAT)
!                 trap(i,j)=trap(i+1,j)*rat+tra(i,j)*(1.-rat)
!              END DO
           ELSE
              IF(MP(I+1).GT.0.0)THEN
                 QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)* &
                      (CL-CPD))+CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))
                 up(i)=up(i+1)
                 vp(i)=vp(i+1)
!                 DO J=1,NTRA
!                    TRAP(I,J)=TRAP(I+1,J)
!                 END DO
              END IF
           END IF
           QP(I)=MIN(QP(I),QSTM)
           QP(I)=MAX(QP(I),0.0)
        ENDIF

     ENDDO

!   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
!   PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)

!   ***  CALCULATE SURFACE PRECIPITATION IN MM/s     ***
    precip=wt(1)*sigd*water(1)/g

    do j = 1, inb
       qflux(j) = wt(j)*sigd*water(j)/g
    enddo

  ENDIF

!   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
!   ***                    WATER VAPOR FLUCTUATIONS                      ***

  WD=BETA*ABS(MP(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))
  QPRIME=0.5*(QP(1)-Q(1))
  TPRIME=LV0*QPRIME/CPD

!   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
!   ***                      AND MIXING RATIO                        ***

  DPINV=0.01/(PH(1)-PH(2))
  AM=0.0
  IF(NK.EQ.1)THEN
     DO K=2,INB
        AM=AM+M(K)
     ENDDO
  END IF
  IF((2.*G*DPINV*AM).GE.DELTI) IFLAG=4
  FT(1)=FT(1)+G*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))
  FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
  FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(T(2)-T(1))*DPINV/CPN(1)
  FQ(1)=FQ(1)+G*MP(2)*(QP(2)-Q(1))*DPINV+SIGD*EVAP(1)
  FQ(1)=FQ(1)+G*AM*(Q(2)-Q(1))*DPINV
  fu(1)=fu(1)+g*dpinv*(mp(2)*(up(2)-u(1))+am*(u(2)-u(1)))
  fv(1)=fv(1)+g*dpinv*(mp(2)*(vp(2)-v(1))+am*(v(2)-v(1)))
!  DO J=1,NTRA
!     FTRA(1,J)=FTRA(1,J)+G*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+AM*(TRA(2,J)-TRA(1,J)))
!  END DO
  AMDE=0.0
  DO J=2,INB
     FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))
     fu(1)=fu(1)+g*dpinv*ment(j,1)*(uent(j,1)-u(1))
     fv(1)=fv(1)+g*dpinv*ment(j,1)*(vent(j,1)-v(1))
!     DO K=1,NTRA
!        FTRA(1,K)=FTRA(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)-TRA(1,K))
!     END DO
  ENDDO

!   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
!   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***

!   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
!   ***                      THROUGH EACH LEVEL                          ***

  DO I=2,INB
     DPINV=0.01/(PH(I)-PH(I+1))
     CPINV=1.0/CPN(I)
     AMP1=0.0
     AD=0.0
     IF(I.GE.NK)THEN
        DO K=I+1,INB+1
           AMP1=AMP1+M(K)
        ENDDO
     END IF
     DO K=1,I
        DO J=I+1,INB+1
           AMP1=AMP1+MENT(K,J)
        ENDDO
     ENDDO
     IF((2.*G*DPINV*AMP1).GE.DELTI) IFLAG=4
     DO K=1,I-1
        DO J=I,INB
           AD=AD+MENT(J,K)
        ENDDO
     ENDDO
     FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))*CPINV) &
          -AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV)) &
          -SIGD*LVCP(I)*EVAP(I)
     FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+ &
          T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV
     FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)* &
          (T(I+1)-T(I))*DPINV*CPINV
     FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))- &
          AD*(Q(I)-Q(I-1)))
     fu(i)=fu(i)+g*dpinv*(amp1*(u(i+1)-u(i))-ad*(u(i)-u(i-1)))
     fv(i)=fv(i)+g*dpinv*(amp1*(v(i+1)-v(i))-ad*(v(i)-v(i-1)))
!     DO K=1,NTRA
!        FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)- &
!              TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
!     END DO
     DO K=1,I-1
        AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
        AWAT=MAX(AWAT,0.0)
        FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
        fu(i)=fu(i)+g*dpinv*ment(k,i)*(uent(k,i)-u(i))
        fv(i)=fv(i)+g*dpinv*ment(k,i)*(vent(k,i)-v(i))
        ! -- sb:
        ! (saturated updrafts resulting from mixing)
         QCOND(I)=QCOND(I)+(ELIJ(K,I)-AWAT)
         NQCOND(I)=NQCOND(I)+1.
         ! sb --
!        DO J=1,NTRA
!           FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-TRA(I,J))
!        END DO
     ENDDO
     DO K=I,INB
        FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
        fu(i)=fu(i)+g*dpinv*ment(k,i)*(uent(k,i)-u(i))
        fv(i)=fv(i)+g*dpinv*ment(k,i)*(vent(k,i)-v(i))
!        DO J=1,NTRA
!           FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-TRA(I,J))
!        END DO
     ENDDO
     FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)* &
          (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
     fu(i)=fu(i)+g*(mp(i+1)*(up(i+1)-u(i))-mp(i)*(up(i)-u(i-1)))*dpinv
     fv(i)=fv(i)+g*(mp(i+1)*(vp(i+1)-v(i))-mp(i)*(vp(i)-v(i-1)))*dpinv

!     DO J=1,NTRA
!        FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))- &
!             MP(I)*(TRAP(I,J)-TRA(I-1,J)))
!     END DO
     ! -- sb:
     ! (saturated downdrafts resulting from mixing)
     DO K=I+1,INB
        QCOND(I)=QCOND(I)+ELIJ(K,I)
        NQCOND(I)=NQCOND(I)+1.
     ENDDO
     ! (particular case: no detraining level is found)
     IF (NENT(I).EQ.0) THEN
        QCOND(I)=QCOND(I)+(1.-EP(I))*CLW(I)
        NQCOND(I)=NQCOND(I)+1.
     ENDIF
     IF (NQCOND(I).GT.1.e-10) THEN
        QCOND(I)=QCOND(I)/NQCOND(I)
     ENDIF
     ! sb --
  ENDDO

!   *** Adjust tendencies at top of convection layer to reflect  ***
!   ***       actual position of the level zero CAPE             ***

  FQOLD=FQ(INB)
  FQ(INB)=FQ(INB)*(1.-FRAC)
  FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PH(INB)-PH(INB+1)) / &
       (PH(INB-1)-PH(INB)))*LV(INB)/LV(INB-1)

  FTOLD=FT(INB)
  FT(INB)=FT(INB)*(1.-FRAC)
  FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PH(INB)-PH(INB+1)) / &
       (PH(INB-1)-PH(INB)))*CPN(INB)/CPN(INB-1)

  fuold = fu(inb)
  fu(inb) = fu(inb)*(1.-frac)
  fu(inb-1) = fu(inb-1)+frac*fuold*((ph(inb)-ph(inb+1)) / &
       (ph(inb-1)-ph(inb)))

  fvold = fv(inb)
  fv(inb) = fv(inb)*(1.-frac)
  fv(inb-1) = fv(inb-1)+frac*fvold*((ph(inb)-ph(inb+1)) / &
       (ph(inb-1)-ph(inb)))

!  DO K=1,NTRA
!     FTRAOLD=FTRA(INB,K)
!     FTRA(INB,K)=FTRA(INB,K)*(1.-FRAC)
!     FTRA(INB-1,K)=FTRA(INB-1,K)+FRAC*FTRAOLD*(PH(INB)-PH(INB+1)) / &
!          (PH(INB-1)-PH(INB))
!  END DO

!   ***   Very slightly adjust tendencies to force exact   ***
!   ***     enthalpy, momentum and tracer conservation     ***

  qav =  precip
  tav = -precip * lv0
  uav = 0.0_r8
  vav = 0.0_r8

  qsum = 0.0_r8
  tsum = 0.0_r8
  usum = 0.0_r8
  vsum = 0.0_r8

  do i=1, inb
     qav = qav + fq(i)*rho(i)*dz(i)
     tav = tav + ft(i)*rho(i)*dz(i)*cpn(i)
     uav = uav + fu(i)*rho(i)*dz(i)
     vav = vav + fv(i)*rho(i)*dz(i)

     qsum = qsum + abs(fq(i))*rho(i)*dz(i)
     tsum = tsum + abs(ft(i))*rho(i)*dz(i)*cpn(i)
     usum = usum + abs(fu(i))*rho(i)*dz(i)
     vsum = vsum + abs(fv(i))*rho(i)*dz(i)
  enddo

  qav = qav / max(qsum, 1.e-20_r8)
  tav = tav / max(tsum, 1.e-20_r8)
  uav = uav / max(usum, 1.e-20_r8)
  vav = vav / max(vsum, 1.e-20_r8)

  do i=1,inb
     fq(i) =  fq(i) - qav * abs(fq(i))
     ft(i) =  ft(i) - tav * abs(ft(i))
     fu(i) = (fu(i) - uav * abs(fu(i))) * (1.-cu)
     fv(i) = (fv(i) - vav * abs(fv(i))) * (1.-cu)
  enddo

!  DO K=1,NTRA
!     TRAAV=0.0
!     DO I=1,INB
!        TRAAV=TRAAV+FTRA(I,K)*(PH(I)-PH(I+1))
!     ENDDO
!     TRAAV=TRAAV/(PH(1)-PH(INB+1))
!     DO I=1,INB
!        FTRA(I,K)=FTRA(I,K)-TRAAV
!     ENDDO
!  ENDDO

  ! IN-CLOUD MIXING RATIO OF CONDENSED WATER :

  ma(inb) = 0.0
  wa(inb) = 1.e-10

  DO I=inb-1,icb,-1
     ma(i) = ma(i+1) + m(i+1)
  ENDDO

  DO I=ICB,INB-1
     AX(I)=0.
     DO J=ICB,I
        AX(I)=AX(I)+RD*(TVP(J)-TV(J))*(PH(J)-PH(J+1))/P(J)
     ENDDO
     IF (AX(I).GT.0.) THEN
        WA(I)=max(SQRT(2.*AX(I)), 1.e-10)
     else
        wa(i)=1.e-10
     endif
  ENDDO

  do i=icb,inb
     siga = ma(i) / ( rho(i) * wa(i) * delta )
     SIGA = MIN(SIGA,0.9)
     SIGA = MAX(SIGA,0.1)
     QCONDC(I) = SIGA * CLW(I)*(1.-EP(I)) + (1.-SIGA) * QCOND(I)
  ENDDO

END SUBROUTINE CONVECT43C


! ---------------------------------------------------------------------------


SUBROUTINE TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ND, NL, KK, ICB, NK
  REAL,    INTENT(IN)  :: GZ(ND), P(ND), T(ND), Q(ND), QS(ND)
  REAL,    INTENT(OUT) :: TVP(ND), TPK(ND), CLW(ND)

  REAL    :: AH0, AHG, ALV, ES, RG, TC, DENOM, QG, CPP, CPINV, TG, S
  INTEGER :: I, J, NSB, NST

!   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***

  REAL, PARAMETER ::  &
       CPD=1005.7,    &
       CPV=1870.0,    &
       CL=2500.0,     &
       RV=461.5,      &
       RD=287.04,     &
       LV0=2.501E6,   &
       CPVMCL=CL-CPV, &
       EPS=RD/RV,     &
       EPSI=1./EPS

!   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***

  AH0=(CPD*(1.-Q(NK))+CL*Q(NK))*T(NK)+Q(NK)*(LV0-CPVMCL* &
       (T(NK)-273.15))+GZ(NK)
  CPP=CPD*(1.-Q(NK))+Q(NK)*CPV
  CPINV=1./CPP

  IF(KK.EQ.1)THEN

!   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***

     DO I=1,ICB-1
        CLW(I)=0.0
     ENDDO
     DO I=NK,ICB-1
        TPK(I)=T(NK)-(GZ(I)-GZ(NK))*CPINV
        TVP(I)=TPK(I)*(1.+Q(NK)*EPSI)
     ENDDO
  END IF

!    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***

  NST=ICB
  NSB=ICB
  IF(KK.EQ.2)THEN
     NST=NL
     NSB=ICB+1
  END IF
  DO I=NSB,NST
     TG=T(I)
     QG=QS(I)
     ALV=LV0-CPVMCL*(T(I)-273.15)
     DO J=1,2
        S=CPD+ALV*ALV*QG/(RV*T(I)*T(I))
        S=1./S
        AHG=CPD*TG+(CL-CPD)*Q(NK)*T(I)+ALV*QG+GZ(I)
        TG=TG+S*(AH0-AHG)
        TG=MAX(TG,35.0)
        TC=TG-273.15
        DENOM=243.5+TC
        IF(TC.GE.0.0)THEN
           ES=6.112*EXP(17.67*TC/DENOM)
        ELSE
           ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
        END IF
        QG=EPS*ES/(P(I)-ES*(1.-EPS))
     ENDDO
     TPK(I)=(AH0-(CL-CPD)*Q(NK)*T(I)-GZ(I)-ALV*QG)/CPD
     CLW(I)=Q(NK)-QG
     CLW(I)=MAX(0.0,CLW(I))
     RG=QG/(1.-Q(NK))
     TVP(I)=TPK(I)*(1.+RG*EPSI)
  ENDDO

  RETURN
END SUBROUTINE TLIFT


END MODULE module_cu_emanuel
