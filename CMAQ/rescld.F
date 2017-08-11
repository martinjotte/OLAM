
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


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/arc/CCTM/src/cloud/cloud_acm_ae6/rescld.F,v 1.7 2011/10/21 16:10:24 yoj Exp $

C what( 1 ) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      subroutine rescld ( mrl )

      use aq_data,     only: tot_wdep
      use mem_grid,    only: mwa
      use utilio_defn, only: m3exit, xstat1
      use cgrid_spcs,  only: n_gc_wdep, n_ae_wdep, n_nr_wdep
      use mem_ijtabs, only: jtab_w, jtw_prog

      implicit none
      
      integer, intent(in) :: mrl

      logical, save       :: firstime = .true.
      integer             :: allocstat
      integer             :: isize
      integer             :: j, iw
      CHARACTER( 30 )     :: XMSG              ! Exit status message
      real                :: COSZEN( mwa )     ! cosine of zenith angle

      if (firstime) then
         firstime = .false.

         isize = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + 1

         ALLOCATE ( TOT_WDEP( mwa,isize ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating TOT_WDEP'
            CALL M3EXIT( 'RESCLD', 0, 0, XMSG, XSTAT1 )
         ENDIF
         
         tot_wdep = 0.0
      endif

      ! Compute zenith angle; Golam Sarwar, July 1, 2011 
      CALL CZANGLE ( COSZEN )

      !$omp parallel do private(iw)
      do j = 1, jtab_w(jtw_prog)%jend(mrl)
         iw =  jtab_w(jtw_prog)%iw(j)

         call rescld1( iw, coszen )

      enddo
      !$omp end parallel do

      end subroutine rescld

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE RESCLD1 ( iw, coszen )

C-----------------------------------------------------------------------
C  FUNCTION: Resolved-scale CLOUD processor Models-3 science process:
 
C  Revision History:
C      No   Date   Who   What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 03/09/98 sjr  made several revisions: fix to read sub-hourly
C                       rainfall data, reordered some of the code
C       2 12/15/98 David Wong at LM
C           -- changed division of GPKG to multiplication of GPKG reciprocal
C           -- interchanged loops structure in line 317
C       3 03/18/99 David Wong at LM
C           -- replace "* M2PHA * ONE_OVER_GPKG" by "* M2PHA_OVER_GPKG" which
C              is a new constant defined as M2PHA / GPKG
C       4 08/30/99 sjr  revised for new aerosol model (with 2nd moments)
C       5 Dec 00   Jeff move CGRID_MAP into f90 module
C       6 01/04/01 sjr  added QS and QI to total water content calcul.
C       7 Sep 01   Jeff Dyn Alloc - Use HGRD_DEFN
C       8 12/18/03 sjr & jp added QG in the water content calc
C       9 07 Dec 04 J.Young: Vert Dyn Alloc - Use VGRD_DEFN
C      10 31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                            domain specifications in one module
C      11 25 Mar 08 sjr fixed bug in the precipitation flux calculation:
C                       layer thickness now included in column integrated
C                       water content and in precipitation flux 
C                       calculations (bug reported by Raymond D Wright)
C      12 12 Aug 10 J.Young: replace CGRID mechanism include files with
C                    namelists and merge Shawn Roselle's, Sergey Napelenok's
C                    and Steve Howard's aerosol reengineering
C      13 01 Mar 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                    removed deprecated TRIMLEN
C      14 11 May D.Wong: incorporated twoway model implementation
C      15 01 Jul 11 G. Sarwar: calculate zenith angle to determine daytime and  
C                    nightime needed for sulfur oxidation via metal catalysis
C      16 02Aug12 S.Roselle:  instrumented to calculate and return
C                             transmissivity for resolved clouds
C      07 Nov 14 J.Bash: Updated call to czangle.F for the ASX_DATA_MOD shared data module.  
C  Called by:  CLDPROC
 
C  Calls the following subroutines:  SCAVWDEP and AQ_MAP
C-----------------------------------------------------------------------

      USE CGRID_SPCS          ! CGRID mechanism species
      USE UTILIO_DEFN

      use cgrid_defn, only: cgrid
      use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, q6, q7,
     &                       pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh
      use micro_coms, only: jnmb, rxmin
      use mem_ijtabs, only: itab_w
      use mem_grid,   only: mwa, mza, lpw, dzt
      use mem_cuparm, only: kcubot, kcutop, qwcon, conprr, iactcu
      use misc_coms,  only: dtlm
      use mem_basic,  only: tair, rho, press
      use aq_data,    only: tot_wdep

      IMPLICIT NONE

C...........Includes:

      INCLUDE 'CONST.EXT'                ! constants

C...........Arguments:
      integer, intent(in)  :: iw
      REAL,    intent(in)  :: COSZEN( mwa )  ! cosine of zenith angle

C...........Parameters:
      REAL, PARAMETER :: GPKG = 1.0E+03  ! g/kg
      REAL, PARAMETER :: M2PHA = 1.0E+04 ! 1 hectare = 1.0e4 m**2
      REAL, PARAMETER :: M2PHA_OVER_GPKG = M2PHA / GPKG
      REAL, PARAMETER :: CNV4 = 1.E-6 * AVO
      REAL, PARAMETER :: CNV4I = 1.0 / CNV4

C Number of species in CGRID
      INTEGER,              SAVE :: MXSPCS
      INTEGER,              save :: N_SPC_WDEP    ! # of wet deposition species
      INTEGER, ALLOCATABLE, SAVE :: WDEP_MAP( : ) ! wet deposition map to CGRID

      INTEGER       JDATE               ! current model date, coded YYYYDDD
      INTEGER       JTIME               ! current model time, coded HHMMSS

C...........Local Variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! flag for first pass thru
!      LOGICAL, SAVE :: QG_AVAIL = .TRUE.   ! flag for QG available on file
!      LOGICAL, SAVE :: QI_AVAIL = .TRUE.   ! flag for QI available on file
!      LOGICAL, SAVE :: QS_AVAIL = .TRUE.   ! flag for QS available on file


      CHARACTER( 120 ) :: XMSG           ! Exit status message
      CHARACTER( 16 ), parameter :: PNAME = 'RESCLD'    ! process name
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
      INTEGER       ALLOCSTAT           ! memory allocation status

      REAL          AIRM                ! total airmass (mol/m2) in cloudy air
      REAL          ALFAC               ! aerosol scavenging coef
      REAL          CTHK1               ! cloud thickness (m)
!     REAL          METSTEP             ! timestep on the met file (hr)
!     SAVE          METSTEP
      REAL          PBARC               ! mean cloud pressure (Pa)
      REAL          PRATE1              ! storm rainfall rate (mm/hr)
      REAL          QCRCOL              ! vert column integrated liquid water content
      REAL          QCRISCOL            ! vert column integrated total water content
      REAL          QRSCOL              ! vert column integrated precip content
      REAL          RAIN                ! non-conv rainfall rate (mm/hr)
      REAL          REMOVAC             ! variable storing H+ deposition
      REAL          TAUCLD              ! cloud lifetime (sec)
      REAL          TBARC               ! mean cloud temp (K)
      REAL          WCBAR               ! liq water content of cloud (kg/m3)
      REAL          WPBAR               ! precipitation water content (kg/m3)
      REAL          WTBAR               ! total water content of cloud (kg/m3)
      REAL          rhoair

      REAL :: POLC ( N_GC_SPCD + N_AE_SPC + N_NR_SPC ) ! incloud conc (mol/mol)
      REAL :: CEND ( N_GC_SPCD + N_AE_SPC + N_NR_SPC ) ! ending conc (mol/mol)
      REAL :: REMOV( N_GC_SPCD + N_AE_SPC + N_NR_SPC ) ! moles/m2 or mm*mol/lit scavenged

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
      real       :: frac ( mza )                  ! cloud fraction
      real       :: fract

!     INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
!     INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2      !Golam Sarwar, July 1, 2011
!     INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2
!     INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C Gridded meteorology data:
C Latitude and longitude for zenith angle calculation: Golam Sarwar * July 1, 2011 
      REAL          COSZ                            ! local cosine of zenith angle
      LOGICAL       DARK                            ! DARK = TRUE is night,  DARK = FALSE is day

      real    :: tc, fracliq, prate, twc, qcu
      integer :: k, mrlw

      logical, allocatable, save :: isnum_spc(:)
      logical, allocatable, save :: issrf_spc(:)

C...........External Functions:

      INTERFACE
        SUBROUTINE SCAVWDEP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                        CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                        REMOV, REMOVAC, ALFAC, RHOAIR )
           INTEGER, INTENT( IN )  :: JDATE, JTIME
           REAL,    INTENT( IN )  :: WTBAR, WCBAR, TBARC, PBARC,
     &                               CTHK1, AIRM, PRATE1, TAUCLD
           REAL,    INTENT( IN )  :: POLC ( : )
           REAL,    INTENT( OUT ) :: REMOVAC
           REAL,    INTENT( OUT ) :: CEND( : ), REMOV( : )
           REAL,    INTENT( OUT ) :: ALFAC
           REAL,    INTENT( IN )  :: RHOAIR
           
        END SUBROUTINE SCAVWDEP
        SUBROUTINE AQ_MAP( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                     CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                     REMOV, REMOVAC, ALFAC, DARK )
           INTEGER, INTENT( IN )    :: JDATE, JTIME
           REAL,    INTENT( IN )    :: WTBAR, WCBAR, TBARC, PBARC,
     &                                 CTHK1, AIRM, PRATE1, TAUCLD
           REAL,    INTENT( IN )    :: POLC ( : )
           REAL,    INTENT( INOUT ) :: REMOVAC
           REAL,    INTENT( INOUT ) :: CEND( : ), REMOV( : )
           REAL,    INTENT( IN )    :: ALFAC
           LOGICAL, INTENT( IN )    :: DARK                                 !Golam Sarwar, July 1, 2011
        END SUBROUTINE AQ_MAP
      END INTERFACE
C-----------------------------------------------------------------------

C...Initialization

      jdate = 0
      jtime = 0

      !$omp critical
      IF ( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

        N_SPC_WDEP = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP
        ALLOCATE ( WDEP_MAP( N_SPC_WDEP ), STAT = ALLOCSTAT )

        IF ( ALLOCSTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating WDEP_MAP'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        IF ( N_SPC_WDEP .LE. 0 ) THEN
           XMSG = 'No species were specified for wet deposition ' //
     &            'tracking'
           CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
        END IF

        SPC = 0
        STRT = 1
        FINI = N_GC_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = GC_STRT - 1 + GC_WDEP_MAP( SPC )
        END DO

        SPC = 0
        STRT = N_GC_WDEP + 1
        FINI = N_GC_WDEP + N_AE_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = AE_STRT - 1 + AE_WDEP_MAP( SPC )
        END DO

        SPC = 0
        STRT = N_GC_WDEP + N_AE_WDEP + 1
        FINI = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = NR_STRT - 1 + NR_WDEP_MAP( SPC )
        END DO

        MXSPCS = N_GC_SPCD + N_AE_SPC + N_NR_SPC

!        ALLOCATE ( CEND ( MXSPCS ),
!     &             POLC ( MXSPCS ),
!     &             REMOV( MXSPCS ), STAT = ALLOCSTAT )
!        IF ( ALLOCSTAT .NE. 0 ) THEN
!          XMSG = 'Failure allocating CEND, POLC or REMOV'
!          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!        END IF

        allocate(isnum_spc(mxspcs)); isnum_spc(:) = .false.
        allocate(issrf_spc(mxspcs)); issrf_spc(:) = .false.

        SPC = 0
        STRT = AE_STRT
        FINI = AE_STRT - 1 + N_AE_SPC
        DO VAR = STRT, FINI
           SPC = SPC + 1
           isnum_spc( var ) = ( INDEX( AE_SPC( SPC ), 'NUM' ) .GT. 0 )
           issrf_spc( var ) = ( INDEX( AE_SPC( SPC ), 'SRF' ) .GT. 0 )
        ENDDO

      END IF
      !$omp end critical

C...Loop through all grid cells

         mrlw = itab_w(iw)%mrlw

         ! Set the cloud lifetime (=adv timestep)
         taucld = dtlm(mrlw)

         tot_wdep(iw,:) = 0.0

C...Sum microphysics precipitation rate over species and convert from (mm/s) to (mm/hr)

         rain = 0.

         if (allocated(pcprd)) rain = rain + pcprd(iw)
         if (allocated(pcprr)) rain = rain + pcprr(iw)
         if (allocated(pcprp)) rain = rain + pcprp(iw)
         if (allocated(pcprs)) rain = rain + pcprs(iw)
         if (allocated(pcpra)) rain = rain + pcpra(iw)
         if (allocated(pcprg)) rain = rain + pcprg(iw)
         if (allocated(pcprh)) rain = rain + pcprh(iw)

         rain = rain * 3600.

C...calculate the integrated column cloud+rain water content
C... JP & SJR:  treat graupel as rainwater until we have a more
C...            advanced cloud microphysics scheme
C... include layer thickness in the column integration

         qc(:) = 0.0
         qr(:) = 0.0
         qi(:) = 0.0
         qs(:) = 0.0

         ! cloud water

         if (jnmb(1) >= 1) then
            do k = lpw(iw), mza
               if (sh_c(k,iw) >= rxmin(1)) then
                  qc(k) = qc(k) + sh_c(k,iw)
               endif
            enddo
         endif

         ! rain

         if (jnmb(2) >= 1) then
            do k = lpw(iw), mza
               if (sh_r(k,iw) >= rxmin(2)) then
                  qr(k) = qr(k) + sh_r(k,iw)
               endif
            enddo
         endif

         ! pristine ice

         if (jnmb(3) >= 1) then
            do k = lpw(iw), mza
               if (sh_p(k,iw) >= rxmin(3)) then
                  qi(k) = qi(k) + sh_p(k,iw)
               endif
            enddo
         endif

         ! snow

         if (jnmb(4) >= 1) then
            do k = lpw(iw), mza
               if (sh_s(k,iw) >= rxmin(4)) then
                  qs(k) = qs(k) + sh_s(k,iw)
               endif
            enddo
         endif

         ! aggregates

         if (jnmb(5) >= 1) then
            do k = lpw(iw), mza
               if (sh_a(k,iw) >= rxmin(5)) then
                  qs(k) = qs(k) + sh_a(k,iw)
               endif
            enddo
         endif

         ! graupel

         if (jnmb(6) >= 1) then
            do k = lpw(iw), mza
               if (sh_g(k,iw) >= rxmin(6)) then
                  call qtc(q6(k,iw),tc,fracliq)
                  qr(k) = qr(k) + sh_g(k,iw) * fracliq 
                  qs(k) = qs(k) + sh_g(k,iw) * (1. - fracliq)
               endif
            enddo
         endif

         ! hail
         
         if (jnmb(7) >= 1) then
            do k = lpw(iw), mza
               if (sh_h(k,iw) >= rxmin(7)) then
                  call qtc(q7(k,iw),tc,fracliq)
                  qr(k) = qr(k) + sh_h(k,iw) * fracliq 
                  qs(k) = qs(k) + sh_h(k,iw) * (1. - fracliq)
               endif
            enddo
         endif

         ! drizzle

         if (jnmb(8) >= 1) then
            do k = lpw(iw), mza
               if (sh_d(k,iw) >= rxmin(8)) then
                  qr(k) = qr(k) + sh_d(k,iw)
               endif
            enddo
         endif

         ! convective clouds

         if (iactcu(iw)) then

            do k = kcubot(iw), kcutop(iw)
               if (tair(k,iw) < 255.15) then
                  fracliq = 0.0
               elseif (tair(k,iw) < 273.15) then
                  fracliq = (tair(k,iw) - 273.15) / 18.0
               else
                  fracliq = 1.0
               endif

               qcu = max( qwcon(k,iw), 1.e-5 )
               qc(k) = qc(k) + qcu * fracliq
               qi(k) = qi(k) + qcu * (1. - fracliq)
            enddo

            if (conprr(iw) > 1.e-10) then

               prate = conprr(iw) * 3600.0
               rain  = rain + prate

               twc = 0.06 * prate**0.846 / 1000.0 ! tot wat cont kg/m3

               do k = lpw(iw), kcutop(iw)-1
                  if (tair(k,iw) < 255.15) then
                     fracliq = 0.0
                  elseif (tair(k,iw) < 273.15) then
                     fracliq = (tair(k,iw) - 273.15) / 18.0
                  else
                     fracliq = 1.0
                  endif

                  qr(k) = qr(k) + twc / real(rho(k,iw)) * fracliq
                  qs(k) = qs(k) + twc / real(rho(k,iw)) * (1. - fracliq)
               enddo

            endif

         endif

         QCRCOL   = 0.0
         QCRISCOL = 0.0
         QRSCOL   = 0.0

         do k = lpw(iw), mza

            QCRCOL   = QCRCOL   + real( rho( k, iw ) ) * dzt( k )
     &               * ( QC( K ) + QR( K ) )

            QCRISCOL = QCRISCOL + real( rho( k, iw ) ) * dzt( k )
     &               * ( QC( K ) + QR( K ) + QI( K ) + QS( K ) )

            QRSCOL   = QRSCOL   + real( rho( k, iw ) ) * dzt( k )
     &               * ( QR( K ) + QS( K ) )

         end do

C...Calculate the cloud optical depth using a formula derived from
C...  Stephens (1978), JAS(35), pp2111-2132.
C...  only calculate the cloud optical depth when the liquid water
C...  path is >= 10 g/m2

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

C...loop through layers

         IF ( QCRCOL .GT. 1.e-5 ) THEN

            ! Diagnose cloud fraction

            call get_cloud_frac(iw, lpw(iw), frac)

            DO k = lpw(iw), mza

C...Compute cloud quantities

               IF ( ( QC( K ) + QR( K ) ) >= 1.e-5 .and. press(k,iw) > 100.d2) THEN

                  fract = max(frac(k), 0.1)

                  TBARC = tair( k, iw )

                  PBARC = press( k, iw )

                  CTHK1 = dzt( K )

                  rhoair = rho( k, iw )

                  airm = rhoair * cthk1 / mwair * 1000.

                  WCBAR = ( QC( K ) + QR( K ) ) * rhoair / fract

                  WTBAR = ( QC( K ) + QR( K )
     &                  +   QI( K ) + QS( K ) ) * rhoair / fract

C...Compute precipitation flux
C...  include layer thickness in the calculation

!                 IF ( QRSGCOL .GT. 1.e-5 ) THEN
!                    WPBAR = ( QR( K ) + QS( K ) ) * rho( k, iw )
!                    PRATE1 = RAIN * WPBAR * CTHK1 / QRSCOL    ! convert to mm/hr
!                 ELSE
                     PRATE1 = RAIN * WTBAR * CTHK1 / QCRISCOL  ! convert to mm/hr
!                 END IF

C...Finally, get in-cloud pollutant concentrations in moles sp
C...  per mole air

                  DO SPC = 1, NSPCSD
                     if ( isnum_spc(spc) ) then
                        polc ( spc ) = max( cgrid( k, iw, spc ) * cnv4,  1.e-36 )
                     else
                        POLC ( SPC ) = max( CGRID( k, iw, SPC ) * 1.e-6, 1.e-36 )
                     endif
                     CEND ( SPC ) = POLC( SPC )
                     REMOV( SPC ) = 0.0
                  END DO

C...perform scavenging and aqueous chemistry within the cloud
C...  and re-adjust the ending and removed amounts for those species
C...  that were scavenged or that participated in cloud chemistry

                  CALL SCAVWDEP ( JDATE, JTIME, WTBAR, WCBAR, TBARC,
     &                            PBARC, CTHK1, AIRM, PRATE1, TAUCLD,
     &                            POLC, CEND, REMOV, REMOVAC, 
     &                            ALFAC, RHOAIR )

C...if the liquid water content is above the specified threshold
C...  then perform the aqueous chemistry within the cloud and
C...  re-adjust the ending and removed amounts for those species
C...  that participated in cloud chemistry

C...  calculate cosine of zenith angle for the cell and determine day or night; Golam Sarwar 

                  COSZ = COSZEN( iw )

                  IF ( COSZ .LE. 0.003 ) THEN
                     DARK = .TRUE. ! night
                  ELSE
                     DARK = .FALSE. ! day
                  END IF

                  CALL AQ_MAP ( JDATE, JTIME, WTBAR, WCBAR, TBARC,
     &                          PBARC, CTHK1, AIRM, PRATE1, TAUCLD,
     &                          POLC, CEND, REMOV, REMOVAC, ALFAC, 
     &                          DARK )

C...convert from moles/m**2 to kg/m**2 and kg/m**2 to kg/hectare

C...  for gases

                  SPC = 0
                  STRT = GC_STRT
                  FINI = GC_STRT - 1 + N_GC_SPC
                  DO VAR = STRT, FINI
                     SPC = SPC + 1
                     REMOV( VAR ) = REMOV( VAR ) * GC_MOLWT( SPC )
     &                            * M2PHA_OVER_GPKG
                  END DO

C...  for aerosols

                  SPC = 0
                  STRT = AE_STRT
                  FINI = AE_STRT - 1 + N_AE_SPC
                  DO VAR = STRT, FINI
                     SPC = SPC + 1
                     IF ( (.not. isnum_spc( var )) .AND.
     &                    (.not. issrf_spc( var )) ) THEN
                        REMOV( VAR ) = REMOV( VAR ) * AE_MOLWT( SPC )
     &                               * M2PHA_OVER_GPKG
                     ELSE
                        REMOV( VAR ) = REMOV( VAR ) * M2PHA
                     END IF
                  END DO

C...  for non-reactives

                  SPC = 0
                  STRT = NR_STRT
                  FINI = NR_STRT - 1 + N_NR_SPC
                  DO VAR = STRT, FINI
                     SPC = SPC + 1
                     REMOV( VAR ) = REMOV( VAR ) * NR_MOLWT( SPC )
     &                            * M2PHA_OVER_GPKG
                  END DO

C...load deposition amounts into the DEP array

                  DO VAR = 1, N_SPC_WDEP
                     TOT_WDEP( iw, VAR ) = TOT_WDEP( iw, VAR )
     &                                   + REMOV( WDEP_MAP( VAR ) ) * fract
                  END DO

C...load H+ concentration into the deposition array as well

                  TOT_WDEP( iw, N_SPC_WDEP+1 ) = TOT_WDEP( iw, N_SPC_WDEP+1 )
     &                                         + REMOVAC * fract

C...set cgrid to the ending concentrations

                  DO SPC = 1, NSPCSD
                     if ( isnum_spc(spc) ) then
                        CGRID( k, iw, SPC ) = ( POLC( SPC ) * (1.0 - fract)
     &                                        + CEND( SPC ) *        fract  ) * cnv4i
                     else
                        CGRID( k, iw, SPC ) = ( POLC( SPC ) * (1.0 - fract)
     &                                        + CEND( SPC ) *        fract  ) * 1.e6
                     endif
                  END DO

               END IF
            END DO
         END IF

      END SUBROUTINE RESCLD1
