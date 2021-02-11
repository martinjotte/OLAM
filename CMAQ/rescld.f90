
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

subroutine rescld ( mrl )

  use mem_ijtabs, only: jtab_w, jtw_prog

  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iw

  !$omp parallel do private(iw) schedule(guided)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     call rescld1( iw )

  enddo
  !$omp end parallel do

end subroutine rescld



subroutine rescld1 ( iw )

!-----------------------------------------------------------------------
!  FUNCTION: Resolved-scale CLOUD processor Models-3 science process:

!  Revision History:
!      No   Date   Who   What
!      -- -------- ---  -----------------------------------------
!       0 01/15/98 sjr  created program
!       1 03/09/98 sjr  made several revisions: fix to read sub-hourly
!                       rainfall data, reordered some of the code
!       2 12/15/98 David Wong at LM
!           -- changed division of GPKG to multiplication of GPKG reciprocal
!           -- interchanged loops structure in line 317
!       3 03/18/99 David Wong at LM
!           -- replace "* M2PHA * ONE_OVER_GPKG" by "* M2PHA_OVER_GPKG" which
!              is a new constant defined as M2PHA / GPKG
!       4 08/30/99 sjr  revised for new aerosol model (with 2nd moments)
!       5 Dec 00   Jeff move CGRID_MAP into f90 module
!       6 01/04/01 sjr  added QS and QI to total water content calcul.
!       7 Sep 01   Jeff Dyn Alloc - Use HGRD_DEFN
!       8 12/18/03 sjr & jp added QG in the water content calc
!       9 07 Dec 04 J.Young: Vert Dyn Alloc - Use VGRD_DEFN
!      10 31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
!                            domain specifications in one module
!      11 25 Mar 08 sjr fixed bug in the precipitation flux calculation:
!                       layer thickness now included in column integrated
!                       water content and in precipitation flux
!                       calculations (bug reported by Raymond D Wright)
!      12 12 Aug 10 J.Young: replace CGRID mechanism include files with
!                    namelists and merge Shawn Roselle's, Sergey Napelenok's
!                    and Steve Howard's aerosol reengineering
!      13 01 Mar 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
!                    removed deprecated TRIMLEN
!      14 11 May D.Wong: incorporated twoway model implementation
!      15 01 Jul 11 G. Sarwar: calculate zenith angle to determine daytime and
!                    nightime needed for sulfur oxidation via metal catalysis
!      16 02Aug12 S.Roselle:  instrumented to calculate and return
!                             transmissivity for resolved clouds
!      07 Nov 14 J.Bash: Updated call to czangle.F for the ASX_DATA_MOD shared data module.
!-----------------------------------------------------------------------

  USE CGRID_SPCS          ! CGRID mechanism species
  USE UTILIO_DEFN
  USE CONST_DATA
  use aero_data,  only: aer_str, aer_end
  use cgrid_defn, only: cgrid
  use mem_micro,  only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, q6, q7, &
                        pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh
  use micro_coms, only: jnmb, rxmin
  use mem_ijtabs, only: itab_w
  use mem_grid,   only: mza, lpw, dzt
  use misc_coms,  only: dtlm
  use mem_basic,  only: tair, rho, press
  use aq_data,    only: tot_wdep, convf, convb
  use therm_lib,  only: qtc
  use mem_radiate,only: cloud_frac, cosz
  use consts_coms,only: rdry
  use mem_co2,    only: rr_co2, i_co2, co2_sh2ppm, co2_ppm2sh, co2_initppm

  implicit none

  integer, intent(in)  :: iw

  real, parameter :: gpkg = 1.0e+03  ! g/kg
  real, parameter :: kgpg = 1.0e-03  ! kg/g
  real, parameter :: m2pha = 1.0e+04 ! 1 hectare = 1.0e4 m**2
  real, parameter :: m2pha_over_gpkg = m2pha / gpkg
  real, parameter :: cnv4 = 1.e-6 * avo
  real, parameter :: cnv4i = 1.0 / cnv4

  INTEGER       JDATE               ! current model date, coded YYYYDDD
  INTEGER       JTIME               ! current model time, coded HHMMSS

!...........Local Variables:

!      LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! flag for first pass thru
!      LOGICAL, SAVE :: QG_AVAIL = .TRUE.   ! flag for QG available on file
!      LOGICAL, SAVE :: QI_AVAIL = .TRUE.   ! flag for QI available on file
!      LOGICAL, SAVE :: QS_AVAIL = .TRUE.   ! flag for QS available on file


!     CHARACTER( 120 ) :: XMSG           ! Exit status message
!     CHARACTER( 16 ), parameter :: PNAME = 'RESCLD'    ! process name
!     CHARACTER( 16 ) ::  VARNM               ! variable name for IOAPI to get
!     CHARACTER( 16 ), SAVE ::  VNAME_RN

!      INTEGER       COL                 ! column loop counter
!      INTEGER       ROW                 ! row loop counter
!      INTEGER       LAY                 ! layer loop counter
  INTEGER       FINI                ! ending position
!     INTEGER, SAVE :: LOGDEV           ! output log unit number
!      INTEGER       MDATE               ! process date (yyyyddd)
!      INTEGER, SAVE :: MSTEP            ! met file time step (hhmmss)
!      INTEGER       MTIME               ! process time (hhmmss)
!      INTEGER, SAVE :: SDATE            ! met file start date
  INTEGER       SPC                 ! liquid species loop counter
  INTEGER       STRT                ! starting position
!      INTEGER, SAVE :: STIME            ! met file start time
!      INTEGER       TCLD                ! cloud lifetime (sec)
  INTEGER       VAR                 ! variable loop counter
!     INTEGER       ALLOCSTAT           ! memory allocation status

  REAL          AIRM                ! total airmass (mol/m2) in cloudy air
! REAL    ::    EALFAC              ! aerosol scavenging coef
  REAL    ::    EWASH               ! rainout coeficient
  REAL          CTHK1               ! cloud thickness (m)
!     REAL          METSTEP             ! timestep on the met file (hr)
!     SAVE          METSTEP
  REAL          PBARC               ! mean cloud pressure (Pa)
! REAL          PRATE1              ! storm rainfall rate (mm/hr)
! REAL          QCRCOL              ! vert column integrated liquid water content
! REAL          QCRISCOL            ! vert column integrated total water content
! REAL          QRSCOL              ! vert column integrated precip content
  real :: qcicol
  REAL          RAIN                ! non-conv rainfall rate (mm/hr)
  REAL          REMOVAC             ! variable storing H+ deposition
  REAL          TAUCLD              ! cloud lifetime (sec)
  REAL          TBARC               ! mean cloud temp (K)
  REAL          WCBAR               ! liq water content of cloud (kg/m3)
!     REAL          WPBAR               ! precipitation water content (kg/m3)
  REAL          WTBAR               ! total water content of cloud (kg/m3)
  REAL          rhoair

  REAL :: POLC ( NSPCSD+1 )            ! incloud conc (mol/mol)
  REAL :: CEND ( NSPCSD+1 )            ! ending conc (mol/mol)
  REAL :: REMOV( NSPCSD+1 )            ! moles/m2 or mm*mol/lit scavenged
  REAL :: DCONC( NSPCSD+1 )
  REAL :: TRFLX( NSPCSD+1 )

!     REAL          RN   ( NCOLS, NROWS ) ! non-convective rainfall (cm)
!     REAL          DENS ( NCOLS, NROWS, NLAYS )  ! air density (kg/m3)
!     REAL          DZZ  ( NCOLS, NROWS, NLAYS )  ! layer thickness (m)
!     REAL          PRES ( NCOLS, NROWS, NLAYS )  ! air pressure (Pa)
  REAL          QC   ( mza )                  ! cloud water content (kg/kg)
!     REAL          QG   ( NCOLS, NROWS, NLAYS )  ! graupel content (kg/kg)
  REAL          QI   ( mza )                  ! ice content (kg/kg)
  REAL          QR   ( mza )                  ! rain water content (kg/kg)
  REAL          QS   ( mza )                  ! snow content (kg/kg)
!     REAL          QV   ( NCOLS, NROWS, NLAYS )  ! specific humidity (kg/kg)
!     REAL          TA   ( NCOLS, NROWS, NLAYS )  ! air temperature (K)
!     REAL          ZH   ( NCOLS, NROWS, NLAYS )  ! mid-layer height (m)
!     REAL          ZF   ( NCOLS, NROWS, NLAYS )  ! level/layer-face height (m)
  real       :: fracliq( mza )
  real       :: qq, qci
  real       :: fract, fl
  real       :: hplus, hprain
  real       :: frac_act
  real       :: pcpflx( mza ), prate1( mza )

!     INTEGER      GXOFF, GYOFF              ! global origin offset from file
! for INTERPX
!     INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2      !Golam Sarwar, July 1, 2011
!     INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2
!     INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

! Gridded meteorology data:
! Latitude and longitude for zenith angle calculation: Golam Sarwar * July 1, 2011
  LOGICAL       DARK                            ! DARK = TRUE is night,  DARK = FALSE is day

  real    :: tc, pdry, qcmin, ff
  integer :: k, mrlw
  integer :: kctop, kitop, ktop

  real, parameter :: onem = nearest(1.0, -1.0)
  real, parameter :: one3 = 1./3.

  jdate = 0
  jtime = 0

  mrlw = itab_w(iw)%mrlw

  ! Set the cloud lifetime (=adv timestep)
  taucld = dtlm(mrlw)

  tot_wdep(iw,:) = 0.0

  qc = 0.0
  qi = 0.0

  qcmin = max(rxmin(1), rxmin(3), 1.e-9)

  ! cloud water

  if (jnmb(1) >= 1) then
     do k = lpw(iw), mza
        qq = rr_c(k,iw) * rho(k,iw)
        if (qq > qcmin) qc(k) = qq
     enddo
  endif

  do k = mza, lpw(iw), -1
     if (qc(k) > qcmin) exit
  enddo
  kctop = k

  ! pristine ice

  if (jnmb(3) >= 1) then
     do k = lpw(iw), mza
        qq = rr_p(k,iw) * rho(k,iw)
        if (qq > qcmin) qi(k) = qq
     enddo
  endif

  do k = mza, lpw(iw), -1
     if (qi(k) > qcmin) exit
  enddo
  kitop = k

  ! exit if no cloud

  if (kctop < lpw(iw) .and. kitop < lpw(iw)) return

  ktop = max(kctop,kitop)
  qr  = 0.0
  qs  = 0.0

  ! rain

  if (jnmb(2) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_r(k,iw) * rho(k,iw)
        if (qq >= rxmin(2)) qr(k) = qq
     enddo
  endif

  ! snow

  if (jnmb(4) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_s(k,iw) * rho(k,iw)
        if (qq >= rxmin(4)) qs(k) = qq
     enddo
  endif

  ! aggregates

  if (jnmb(5) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_a(k,iw) * rho(k,iw)
        if (qq >= rxmin(5)) qs(k) = qs(k) + qq
     enddo
  endif

  ! graupel

  if (jnmb(6) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_g(k,iw) * rho(k,iw)
        if (qq >= rxmin(6)) then
           call qtc(q6(k,iw)/rr_g(k,iw),tc,fl)
           qr(k) = qr(k) + qq * fl
           qs(k) = qs(k) + qq * (1. - fl)
        endif
     enddo
  endif

  ! hail

  if (jnmb(7) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_h(k,iw) * rho(k,iw)
        if (qq >= rxmin(7)) then
           call qtc(q7(k,iw)/rr_h(k,iw),tc,fl)
           qr(k) = qr(k) + qq * fl
           qs(k) = qs(k) + qq * (1. - fl)
        endif
     enddo
  endif

  ! drizzle

  if (jnmb(8) >= 1) then
     do k = lpw(iw), ktop
        qq = rr_d(k,iw) * rho(k,iw)
        if (qq >= rxmin(8)) qr(k) = qr(k) + qq
     enddo
  endif

  ! Liquid fraction (monotonically increasing for now)

  fracliq(ktop) = (qc(ktop) + qr(ktop)) / (qc(ktop) + qr(ktop) + qi(ktop) + qs(ktop))

  do k = ktop-1, lpw(iw), -1
     qq = qc(k) + qr(k) + qi(k) + qs(k)
     if (qq < rxmin(1)) then
        fracliq(k) = fracliq(k+1)
     else
        fracliq(k) = max(fracliq(k+1), (qc(k) + qr(k)) / qq)
     endif
  enddo

  !...Sum microphysics precipitation rate over species

  rain = 0.

  if (allocated(pcprd)) rain = rain + pcprd(iw)
  if (allocated(pcprr)) rain = rain + pcprr(iw)
  if (allocated(pcprp)) rain = rain + pcprp(iw)
  if (allocated(pcprs)) rain = rain + pcprs(iw)
  if (allocated(pcpra)) rain = rain + pcpra(iw)
  if (allocated(pcprg)) rain = rain + pcprg(iw)
  if (allocated(pcprh)) rain = rain + pcprh(iw)

! Washout rate based on surface rain and integrated cloud water
! until linkage with microphysics

  pcpflx = 0.
  prate1 = 0.

  if (rain >= 1.e-12) then

     qcicol = 0.0
     do k = lpw(iw), ktop
        qcicol = qcicol + dzt(k) * ( qc(k) + qi(k) )
     enddo

     do k = lpw(iw), ktop
        qq = qc(k) + qi(k)
        if (qq >= qcmin) then
           prate1(k) = rain * qq * dzt(k) / qcicol
        endif
     enddo

     do k = ktop, lpw(iw), -1
        pcpflx(k) = pcpflx(k+1) + prate1(k)
     enddo

     if (pcpflx(lpw(iw)) < 1.e-12) then

        rain  = 0.0
        ewash = 1.0
        pcpflx(lpw(iw):ktop) = 0.0
        prate1(lpw(iw):ktop) = 0.0

     else

        ewash = max( exp( -taucld * rain / qcicol ), 0.3)

     endif

  else

     rain  = 0.0
     ewash = 1.0

  endif

!...Calculate the cloud optical depth using a formula derived from
!...  Stephens (1978), JAS(35), pp2111-2132.
!...  only calculate the cloud optical depth when the liquid water
!...  path is >= 10 g/m2

!          LWP = QCICOL * 1000.0  ! converts to g/m2
!          IF ( LWP .GE. 10.0 ) THEN
!             CLOD = 10.0**( 0.2633 + 1.7095 * LOG( LOG10( LWP ) ) )
!          ELSE
!             CLOD = 0.0
!          END IF
!
!C...If no cloud or optical depth < 5, set clear sky values.
!C...  (i.e. don't do anything)
!
!          IF ( CLOD .GE. 5.0 ) THEN
!
!             RESTRANS( COL, ROW ) = ( 5.0 - EXP ( -CLOD ) ) / ( 4.0 + 0.42 * CLOD )
!
!          END IF

!...loop through layers

!        IF ( QCRCOL .GT. 1.e-5 ) THEN

  hprain  = 2.e-6
  trflx   = 0.0
  remov   = 0.0
  removac = 0.0

  DO k = ktop, lpw(iw), -1

     qci = qc(k) + qi(k)

     ! Skip this cell if no clouds and precipitation

     if (qci < qcmin .and. pcpflx(k+1) < 1.e-12) cycle

     ! Get in-cloud pollutant concentrations in moles
     ! per mole air

     do spc = 1, nspcsd
        polc( spc ) = max( cgrid( k, iw, spc ) * convf( spc ),  1.e-36 )
     enddo

     ! Include CO2 for aqueous / wet depostion

     if (i_co2 > 0) then
        polc( nspcsd+1 ) = rr_co2(k,iw) * co2_sh2ppm * 1.e-6 ! convert to mol/molV
     else
        polc( nspcsd+1 ) = co2_initppm * 1.e-6  ! convert to mol/molV
     endif

     cend   = polc
     tbarc  = tair (k,iw)
     pbarc  = press(k,iw)
     cthk1  = dzt  (k)
     rhoair = rho  (k,iw)

     frac_act = 0.0
     prate1   = 0.0

     if (qci >= qcmin) then

        fract = max(cloud_frac(k,iw), 0.2)

        WCBAR = ( QC( K ) + QR( K ) ) / fract

        WTBAR = ( QC( K ) + QR( K ) &
                + QI( K ) + QS( K ) ) / fract

        ! fraction of aerosol activated in mixed clouds, from
        ! Verheggen et al. 2007 JGR
        fl = wcbar / wtbar

        if (fl < .999) then
           frac_act = 0.051 + 0.949 * exp(6.*(fl - 1.0))
        else
           frac_act = 1.0
        endif

        ! limit aerosol activation for small water content?
        if (qci < 1.e-5) then
           frac_act = min(frac_act, sqrt(qci * 1.e5))
        endif

        airm  = rhoair * cthk1 * inv_mwairkg

        ! IN-CLOUD AITKEN SCAVENGING;
        ! SCAVENGED MASS ADDED TO ACCUMULATION MODE, NO DEPOSITION
        call aitken_incloud_scav(cend, wtbar, rhoair, tbarc, pbarc, taucld)

        if (wcbar > 1.e-6) then

           dark = ( cosz(iw) < 0.003 )
           pdry = rhoair * tbarc * rdry

           ! AQUEOUS CHEMISTRY
           call aq_map(k, iw, wcbar, tbarc, pdry, taucld, cend, airm, &
                       frac_act, remov, hplus, removac, ewash, dark )

           ! IN-CLOUD RAINOUT OF GASES NOT INVOLVED IN AQUEOUS (Henry's Law)
           call gas_incloud_rainout_skipaq( ewash, wcbar, tbarc, hplus, &
                                            airm, cend, remov )

        else

           ! Typical H+ value for clean rain
           hplus = 2.e-6

           ! IN-CLOUD RAINOUT OF GASES (Henry's Law)
           call gas_incloud_rainout( ewash, wcbar, tbarc, hplus, &
                                     airm, cend, remov )

        endif

        ! IN-CLOUD RAINOUT / SNOWOUT OF AEROSOLS
        call aero_incloud_rainout( ewash, airm, frac_act, cend, remov )

        ! TAKE INTO ACCOUNT CLOUD FRACTION
        do spc = 1, nspcsd+1
           cend(spc) = polc( spc ) * (1.0 - fract) &
                     + cend( spc ) *        fract
        enddo

        ! CONVERT POLLUTANTS REMOVED INTO A FLUX
        if (ewash < onem) then
           ff = fract / taucld
           do spc = 1, nspcsd+1
              trflx(spc) = trflx(spc) + remov(spc) * ff
           enddo
        endif

     endif  ! qci > 1.e-10)

     ! Washout from precipitation falling from above

     if (pcpflx(k+1) >= 1.e-12) then
        call precip_scav(cend, dconc, remov, tbarc, rhoair, pbarc, cthk1, frac_act, &
                         hprain, pcpflx(k+1), trflx, fracliq(k), qc(k))

        do spc = 1, nspcsd+1
           cend (spc) = cend (spc) + dconc(spc) * taucld
           trflx(spc) = trflx(spc) + remov(spc)
        enddo
     endif

     ! Add rain formed in this cell to precip flux

     if (prate1(k) > 1.e-15) then
        hprain = (hprain * pcpflx(k+1) + hplus * prate1(k)) / pcpflx(k)
     endif

     ! Set cgrid to the ending concentrations

     do spc = 1, nspcsd
        cgrid( k, iw, spc ) = cend( spc ) * convb( spc )
     end do

  enddo

!...convert from moles/m**2 to kg/m**2 and kg/m**2 to kg/hectare

!...  for gases

     STRT = GC_STRT
     FINI = GC_STRT - 1 + N_GC_SPC
     DO VAR = STRT, FINI
        SPC = VAR - STRT + 1
        TRFLX( VAR ) = TRFLX( VAR ) * GC_MOLWT( SPC ) * KGPG
     END DO

!...  for aerosols

     STRT = AER_STR
     FINI = AER_END
     DO VAR = STRT, FINI
        SPC = VAR - STRT + 1
        TRFLX( VAR ) = TRFLX( VAR ) * AE_MOLWT( SPC ) * KGPG
     END DO

!...  for non-reactives

     STRT = NR_STRT
     FINI = NR_STRT - 1 + N_NR_SPC
     DO VAR = STRT, FINI
        SPC = VAR - STRT + 1
        TRFLX( VAR ) = TRFLX( VAR ) * NR_MOLWT( SPC ) * KGPG
     END DO

   END SUBROUTINE RESCLD1
