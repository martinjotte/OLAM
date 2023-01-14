
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
C $Header: /project/yoj/arc/CCTM/src/depv/m3dry/DEPV_DEFN.F,v 1.12 2012/01/19 14:22:45 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      MODULE DEPV_DEFN

C Function: deposition velocity interface to the chemistry-transport model

C Revision History:
C   24 Nov 2006 J.Young: initial implementation using Bill Hutzell's and
C                        Tanya Otte's work
C   29 Jan 2010 D. Wong: Eliminate potential race condition at the MYPE = 0
C                        to open the diagnostic file. Eliminate the BARRIER
C                        and the CLOSE3 and OPEN3.
C   21 Jun 2010 J.Young: convert for Namelist redesign
C   27 Jan 2011 J.Young: remove references to AERO3 and AERO4
C   20 Jan 2011 J. Bash: Passed the component fluxes and tstep between vdiff
C                        & m3dry
C   16 Feb 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C    5 May 2011 J.Young: remove writing diagnostic output for initial call
C   22 Sep 2011 D. Wong: - incorporated twoway model implemenation
C                        - used a robust way to compute WTIME
C   11 Oct 2011 J.Young: PUBLIC/PRIVATE sections, re-do WTIME computation
C                        and WTIME computation for MOSAIC
C   13 Dec 2011 J.Young: removed separate SOA variables - no longer used

      IMPLICIT NONE

C model depositon velocites array
!     REAL,    ALLOCATABLE, SAVE :: DEPV( :,:,: )

C model depositon velocites array by land use category
!     REAL,    ALLOCATABLE, SAVE :: DEPVJ( :,:,:,: )

C model depositon velocites array by land use category
!     REAL,    ALLOCATABLE, SAVE :: DEPVJ_FST( :,:,:,: )

C no. of deposition velocity model species
      INTEGER, SAVE :: N_SPC_DEPV

C ddep spc conversion factors
!     REAL,    ALLOCATABLE, SAVE :: DD_CONV( : )

C model-oriented bi-directional production [ppm-m/s]
!     REAL,    ALLOCATABLE, SAVE :: PLDV( :,:,: )
!     REAL,    ALLOCATABLE, SAVE :: CMP ( :,:,: )

C flag for HONO interaction with surfaces within in-lining depv
      LOGICAL, SAVE :: SFC_HONO = .FALSE.

C flag for ammonia bi-directional flux with in-lining depv
      LOGICAL, SAVE :: ABFLUX

C flag for mosaic - luc specific deposition - within in-lining depv
      LOGICAL, SAVE :: MOSAIC = .FALSE.

C flag for mosaic - luc specific deposition - stomatal flux only - within in-lining depv
      LOGICAL, SAVE :: FST = .FALSE.

C gas species (gas-phase and gas non-reactive) deposition velocities
      REAL,    ALLOCATABLE, SAVE :: DEPVEL_GAS( :,: )

C m3dry-oriented bi-directional production [ppm-m/s]
!     REAL,    ALLOCATABLE, SAVE :: PVD( :,:,: )
!     REAL,    ALLOCATABLE, SAVE :: PVD_land( :,: )

C CGRID for m3dry bi-directional flux calculation
!     REAL,    ALLOCATABLE, SAVE :: CGRIDL1( :,:,: )

C flag for diagnostic DEPV file
      LOGICAL, SAVE :: DEPV_DIAG
      INTEGER, SAVE :: SDATE, STIME

      INTEGER, SAVE                      :: N_UNIQUE_GDEPV   ! gas depv
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: UNIQUE_GDEPV( : )
      INTEGER, ALLOCATABLE, SAVE         :: GAS_DEPV_SUR( : )
      LOGICAL, ALLOCATABLE, SAVE         :: GAS_DEPV_FOUND( : )

!     INTEGER, SAVE                      :: N_UNIQUE_ADEPV   ! aero depv
!     CHARACTER( 16 ), ALLOCATABLE, SAVE :: UNIQUE_ADEPV( : )
!     INTEGER, ALLOCATABLE, SAVE         :: AER_DEPV_SUR( : )

      INTEGER, SAVE :: NR_OFFSET            ! Offset to NR start in CGRID
      INTEGER, SAVE :: N_GAS_DEPV

      REAL, ALLOCATABLE :: ADEPV( :,: )    ! Aerosol dep vel from
                                           ! subroutine AERO_DEPV

!     REAL,    ALLOCATABLE, SAVE :: DBUFF( :,: )      ! diagnostic output buffer
!     REAL,    ALLOCATABLE, SAVE :: DBUFF3( :,:,: )   ! diagnostic output buffer with 3 dimensions

      integer :: ie_hono = 0
      integer :: ic_no2  = 0

      INTEGER, PARAMETER :: N_AE_DEP_SPC = 9 ! no. of surrogates for aerosol dry dep velocities

      CHARACTER( 16 ) :: VDAE_NAME( N_AE_DEP_SPC ) ! dep vel surrogate name table
      DATA         VDAE_NAME( 1 ) / 'VNUMATKN' /
      DATA         VDAE_NAME( 2 ) / 'VNUMACC ' /
      DATA         VDAE_NAME( 3 ) / 'VNUMCOR ' /
      DATA         VDAE_NAME( 4 ) / 'VMASSI  ' /
      DATA         VDAE_NAME( 5 ) / 'VMASSJ  ' /
      DATA         VDAE_NAME( 6 ) / 'VMASSC  ' /
      DATA         VDAE_NAME( 7 ) / 'VSRFATKN' /
      DATA         VDAE_NAME( 8 ) / 'VSRFACC ' /
      DATA         VDAE_NAME( 9 ) / 'VSRFCOR ' /

      INTEGER, ALLOCATABLE, SAVE :: AE_DEPV_SUR( : )   ! pointer to surrogate

      CONTAINS

C-----------------------------------------------------------------------

         subroutine DEPV_INIT( )

         USE CGRID_SPCS
         USE DEPVVARS
         USE UTILIO_DEFN
         USE RXNS_DATA
         USE CONST_DATA
         USE AERO_DATA
         use mem_grid,   only: mwa
         use oname_coms, only: nl

         IMPLICIT NONE

C Arguments:
         INTEGER JDATE, JTIME, TSTEP( 3 )
         LOGICAL SUCCESS

C External Functions:
         INTEGER, EXTERNAL :: FINDEX               ! finds index of a number in a list

C Parameters:
         REAL, PARAMETER :: RGAS1 = 1.0 / RGASUNIV ! univ. gas constant reciprocal
         REAL, PARAMETER :: M2PHA = 1.0E+04        ! 1 hectare = 1.0e4 m**2
         REAL, PARAMETER :: CMLMR = 1.0E+06        ! ppmV/Molar Mixing Ratio
         REAL, PARAMETER :: CNVTD = M2PHA / CMLMR / MWAIR ! combined ddep conv factor

C Local variables:
         CHARACTER( 16 ), ALLOCATABLE :: ALL_GAS_DEPV( : )
!        LOGICAL, ALLOCATABLE         :: AE_DEPV_FOUND( : )

C environment variable for in-line DEPV
!        CHARACTER( 16 )  :: CTM_ILDEPV = 'CTM_ILDEPV      '

C environment variable for ammonia bi-directional flux
!        CHARACTER( 16 )  :: CTM_ABFLUX = 'CTM_ABFLUX      '

C environment variable for HONO interaction with surfaces
!        CHARACTER( 16 )  :: CTM_SFC_HONO = 'CTM_SFC_HONO    '

C variables for deposition velocities diagnostic file
!        INTEGER          :: STATUS      ! ENV... status
!        CHARACTER( 80 )  :: VARDESC     ! environment variable description

C environment variable for diagnostic DEPV file
!        CHARACTER( 16 )  :: CTM_DEPV_FILE = 'CTM_DEPV_FILE    '

        CHARACTER( 16 )  :: PNAME = 'DEPV_INIT       '
        CHARACTER( 120 ) :: XMSG = ' '

         INTEGER          LOGDEV, ALLOCSTAT
         INTEGER          N, V, STRT

         INTERFACE
!            SUBROUTINE OPDEPV_DIAG ( JDATE, JTIME, TSTEP,
!     &                               N_GDEPV_NAMES, GDEPV_NAMES,
!     &                               N_ADEPV_NAMES, ADEPV_NAMES )
!!    &                               N_SDEPV_NAMES, SDEPV_NAMES )
!               INTEGER,         INTENT( IN ) :: JDATE, JTIME, TSTEP
!               INTEGER,         INTENT( IN ) :: N_GDEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: GDEPV_NAMES( : )
!               INTEGER,         INTENT( IN ) :: N_ADEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: ADEPV_NAMES( : )
!!              INTEGER,         INTENT( IN ) :: N_SDEPV_NAMES
!!              CHARACTER( 16 ), INTENT( IN ) :: SDEPV_NAMES( : )
!            END SUBROUTINE OPDEPV_DIAG
!            SUBROUTINE OPDEPV_MOS ( JDATE, JTIME, TSTEP,
!     &                              N_GDEPV_NAMES, GDEPV_NAMES,
!     &                              N_ADEPV_NAMES, ADEPV_NAMES )
!!    &                              N_SDEPV_NAMES, SDEPV_NAMES )
!               INTEGER,         INTENT( IN ) :: JDATE, JTIME, TSTEP
!               INTEGER,         INTENT( IN ) :: N_GDEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: GDEPV_NAMES( : )
!               INTEGER,         INTENT( IN ) :: N_ADEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: ADEPV_NAMES( : )
!!              INTEGER,         INTENT( IN ) :: N_SDEPV_NAMES
!!              CHARACTER( 16 ), INTENT( IN ) :: SDEPV_NAMES( : )
!            END SUBROUTINE OPDEPV_MOS
!            SUBROUTINE OPDEPV_FST ( JDATE, JTIME, TSTEP,
!     &                              N_GDEPV_NAMES, GDEPV_NAMES,
!     &                              N_ADEPV_NAMES, ADEPV_NAMES )
!!    &                              N_SDEPV_NAMES, SDEPV_NAMES )
!               INTEGER,         INTENT( IN ) :: JDATE, JTIME, TSTEP
!               INTEGER,         INTENT( IN ) :: N_GDEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: GDEPV_NAMES( : )
!               INTEGER,         INTENT( IN ) :: N_ADEPV_NAMES
!               CHARACTER( 16 ), INTENT( IN ) :: ADEPV_NAMES( : )
!!              INTEGER,         INTENT( IN ) :: N_SDEPV_NAMES
!!              CHARACTER( 16 ), INTENT( IN ) :: SDEPV_NAMES( : )
!            END SUBROUTINE OPDEPV_FST
            SUBROUTINE GAS_DEPV_MAP ( N_DEPV_NAMES, DEPV_NAMES, DEPV_SUR )
               INTEGER,         INTENT( IN )    :: N_DEPV_NAMES
               CHARACTER( 16 ), INTENT( INOUT ) :: DEPV_NAMES( : )
               INTEGER,         INTENT( OUT )   :: DEPV_SUR( : )
            END SUBROUTINE GAS_DEPV_MAP
         END INTERFACE

C-----------------------------------------------------------------------

         if (nl%do_drydep == 0) return

         LOGDEV = INIT3()
         JDATE=0; JTIME=0; TSTEP=0

         SUCCESS = .TRUE.

         SDATE = JDATE; STIME = JTIME

         ie_hono = index1( 'HONO', n_gc_emis, gc_emis )
         if (ie_hono == 0) then
            xmsg = 'Error: HONO is not in the emissions table'
            call m3exit( pname, 0, 0, xmsg, xstat3)
         endif

         ic_no2 = index1( 'NO2', n_gc_spc, gc_spc )
         if (ic_no2 == 0) then
            xmsg = 'Error: NO2 is not in the CGRID'
            call m3exit( pname, 0, 0, xmsg, xstat3)
         endif

         N_GAS_DEPV = N_GC_DEPV  + N_NR_DEPV ! + N_TR_DEPV
         N_SPC_DEPV = N_GAS_DEPV + N_AE_DEPV

!a         ALLOCATE ( DEPV( N_SPC_DEPV+1,NCOLS,NROWS ),
!a     &              PLDV( N_SPC_DEPV,NCOLS,NROWS ), STAT = ALLOCSTAT )
!a         IF ( ALLOCSTAT .NE. 0 ) THEN
!a            XMSG = 'Failure allocating DEPV or PLDV arrays'
!a            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!a            SUCCESS = .FALSE.; RETURN
!a         END IF

C Initialize deposition velocities for nondeposited species to zero:
C setting DEPV( N_SPC_DEPV+1,C,R ) = 0.0 accounts for dry dep. species names
C as a subset of the vert. diffused species list
!a         DEPV = 0.0   ! array assignment
!a         PLDV = 0.0   ! array assignment

!         ALLOCATE ( DD_CONV( N_SPC_DEPV+1 ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating DD_CONV'
!            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!            SUCCESS = .FALSE.; RETURN
!         END IF

C set dry dep conversion factor (done here so that vdiff doesn`t need the include files)
!         S = 0
!         DO V = 1, N_GC_DEPV
!            S = S + 1
!            DD_CONV( S ) = CNVTD * GC_MOLWT( GC_DEPV_MAP( V ) )
!         END DO
!
!         DO V = 1, N_AE_DEPV
!            S = S + 1
!            IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'NUM' ) THEN
!               DD_CONV( S ) = CNVTD * AVO * 1.0E+03    ! --> #/Ha
!            ELSE IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'SRF' ) THEN
!               DD_CONV( S ) = M2PHA * 1.0E+03 / MWAIR  ! --> M**2/Ha
!            ELSE
!               DD_CONV( S ) = CNVTD * AE_MOLWT( AE_DEPV_MAP( V ) )
!            END IF
!         END DO
!
!         DO V = 1, N_NR_DEPV
!            S = S + 1
!            DD_CONV( S ) = CNVTD * NR_MOLWT( NR_DEPV_MAP( V ) )
!         END DO
!
!         DO V = 1, N_TR_DEPV
!            S = S + 1
!            DD_CONV( S ) = CNVTD * TR_MOLWT( TR_DEPV_MAP( V ) )
!         END DO

C Pleim`s ammonia bi-directional flux with in-line deposition velocities calculation?
         ABFLUX = .FALSE.

C Sarwar`s HONO interaction with leaf and building surfaces?
         SFC_HONO = .TRUE.

C Mosaic method to get land-use specific deposition velocities?
         MOSAIC = .FALSE.

C Get diagnostic file flag.
         DEPV_DIAG = .FALSE.    ! default

         IF ( N_GAS_DEPV .GT. 0 ) THEN

            ALLOCATE ( ALL_GAS_DEPV( N_GAS_DEPV ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating ALL_GAS_DEPV'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

            DO V = 1, N_GC_DEPV
               ALL_GAS_DEPV( V ) = GC_DEPV( V )
            END DO
            STRT = N_GC_DEPV
            DO V = 1, N_NR_DEPV
               ALL_GAS_DEPV( V + STRT ) = NR_DEPV( V )
            END DO
            STRT = N_GC_DEPV + N_NR_DEPV
            DO V = 1, N_TR_DEPV
               ALL_GAS_DEPV( V + STRT ) = TR_DEPV( V )
            END DO

C Identify and remove multiple occurrences of gas deposition velocities
            ALLOCATE( GAS_DEPV_FOUND( N_GAS_DEPV ),
     &                GAS_DEPV_SUR  ( N_GAS_DEPV ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating GAS_DEPV_FOUND or GAS_DEPV_SUR'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

            GAS_DEPV_FOUND = .FALSE.
            N_UNIQUE_GDEPV = 0

            LOOP_UNIQUE:
     &      DO N = 1, N_GAS_DEPV
               IF ( GAS_DEPV_FOUND( N ) ) CYCLE LOOP_UNIQUE
               DO V = 1, N_GAS_DEPV
                  IF ( ALL_GAS_DEPV( N ) .EQ. ALL_GAS_DEPV( V ) ) THEN
                     IF ( GAS_DEPV_FOUND( V ) ) CYCLE LOOP_UNIQUE
                     N_UNIQUE_GDEPV = N_UNIQUE_GDEPV + 1
                     GAS_DEPV_FOUND( V ) = .TRUE.
                     CYCLE LOOP_UNIQUE
                  END IF
               END DO
            END DO LOOP_UNIQUE

            ALLOCATE( UNIQUE_GDEPV( N_UNIQUE_GDEPV ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating UNIQUE_GDEPV'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

            V = 1
            DO N = 1, N_GAS_DEPV
               IF ( GAS_DEPV_FOUND( N ) ) THEN
                  UNIQUE_GDEPV( V ) = ALL_GAS_DEPV( N )
!#ifdef Verbose
!                  write( logdev,1019 ) n, all_gas_depv( n ), v, unique_gdepv( v )
!1019              format( 2x, i5, 1x, a16, 2x, i5, 1x, a16 )
!#endif
                  V = V + 1
!#ifdef Verbose
!               else
!                  write( logdev,1019 ) n, all_gas_depv( n ), -1, ' '
!#endif
               END IF
            END DO

C Get pointers to model species (GAS_DEPV_SUR)
            CALL GAS_DEPV_MAP ( N_UNIQUE_GDEPV, UNIQUE_GDEPV, GAS_DEPV_SUR )

C DEPVEL_GAS is used also for non-ILDEPV in the sgtot calculation

            ALLOCATE ( DEPVEL_GAS( N_UNIQUE_GDEPV, MWA ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DEPVEL_GAS arrays'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF
            DEPVEL_GAS = 0.0

!a            IF ( ILDEPV ) THEN
!a               ALLOCATE ( CGRIDL1( N_GAS_DEPV,NCOLS,NROWS ),
!a     &                    PVD    ( N_GAS_DEPV,NCOLS,NROWS ),
!a     &                    CMP    ( LCMP,NCOLS,NROWS),STAT = ALLOCSTAT )

!a            ALLOCATE ( pvd_land( mlandflux, N_GAS_DEPV ), stat = allocstat )
!a            IF ( ALLOCSTAT .NE. 0 ) THEN
!a               XMSG = 'Failure allocating CGRIDL1, CMP or PVD arrays'
!a               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!a               SUCCESS = .FALSE.; RETURN
!a            END IF

!a               CGRIDL1 = 0.0       ! array assignment (PVD initalized in m3dry)
!a            END IF

!            VARDESC = 'Flag for writing the DEPV diagnostic file'
!            DEPV_DIAG = ENVYN( CTM_DEPV_FILE, VARDESC, DEPV_DIAG, STATUS )
!            IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
!            IF ( STATUS .EQ. 0 ) THEN
!               WRITE( LOGDEV,'( 5X, A )' ) VARDESC
!               XMSG = TRIM( PNAME ) // ': writes GAS DEPV values to '
!     &                              // TRIM( CTM_DEPV_FILE )
!               WRITE( LOGDEV,'( 5X, A )' ) XMSG
!            ELSE IF ( STATUS .EQ. 1 ) THEN
!               XMSG = 'Environment variable improperly formatted'
!               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!               SUCCESS = .FALSE.; RETURN
!            ELSE IF ( STATUS .EQ. -1 ) THEN
!               XMSG =
!     &         'Environment variable set, but empty ... Using default:'
!               WRITE( LOGDEV,'( 5X, A, I9)' ) XMSG, JTIME
!            ELSE IF ( STATUS .EQ. -2 ) THEN
!               XMSG = 'Environment variable not set ... Using default:'
!               WRITE( LOGDEV,'( 5X, A, I9 )' ) XMSG, JTIME
!            END IF

         END IF

!         IF ( .NOT. ILDEPV ) THEN
!
!               DEPV_DIAG = .FALSE.   ! make sure to avoid bogus write attempt later
!
!         ELSE   ! ILDEPV
!
!            IF ( DEPV_DIAG ) THEN
!
!               ALLOCATE ( DBUFF( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
!               IF ( ALLOCSTAT .NE. 0 ) THEN
!                  XMSG = 'Failure allocating DBUFF'
!                  CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!                  SUCCESS = .FALSE.; RETURN
!               END IF
!
!               IF ( MYPE .EQ. 0 ) THEN
!                  CALL OPDEPV_DIAG ( JDATE, JTIME, TSTEP( 1 ),
!     &                               N_UNIQUE_GDEPV, UNIQUE_GDEPV,
!     &                               N_UNIQUE_ADEPV, UNIQUE_ADEPV )
!               END IF
!
!            END IF
!
!         END IF

         IF ( N_AE_DEPV .GT. 0 ) THEN

            ALLOCATE ( ADEPV(N_AE_DEP_SPC, MWA), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating ADEPV'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF
            ADEPV = 0.0

            ALLOCATE ( AE_DEPV_SUR( N_AE_DEPV ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DEPV_SUR'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ! Set the dep vel surrogate pointers
            DO V = 1, N_AE_DEPV
               N = INDEX1( AE_DEPV( V ), N_AE_DEP_SPC, VDAE_NAME )
               IF ( N .NE. 0 ) THEN
                  AE_DEPV_SUR( V ) = N
               ELSE
                  XMSG = 'Could not find ' // AE_DEPV( V ) // ' in aerosol' //
     &                   ' surrogate table. >>> Dep vel set to zero <<< '
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!                 DEPV_SUR( V ) = 0
               END IF
            END DO

         END IF

         END SUBROUTINE DEPV_INIT

!-----------------------------------------------------------------------

         subroutine get_depv( )

         use mem_ijtabs, only: jtab_w, jtw_prog
         use mem_grid,   only: lsw
         use oname_coms, only: nl

         implicit none

         integer             :: j, iw

         if (n_spc_depv < 1 .or. nl%do_drydep == 0) return

         !$omp parallel do private(j,iw)
         do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
            call get_depv1(iw, lsw(iw))
         enddo
         !$omp end parallel do

         end subroutine get_depv

!-----------------------------------------------------------------------

         SUBROUTINE GET_DEPV1( iw, nsfc )

         USE CGRID_SPCS          ! CGRID mechanism species
         USE UTILIO_DEFN
         USE CONST_DATA
         use cgrid_defn, only: cgrid
         use misc_coms,  only: dtlm
         use mem_grid,   only: lpw, volt, arw0
         use mem_basic,  only: rho
         use mem_turb,   only: frac_sfc
         use cgrid_defn, only: cgrid, sxfer_gc, sxfer_nr, sxfer_ae

         IMPLICIT NONE

         ! Arguments:
         integer, intent(in) :: iw
         integer, intent(in) :: nsfc

         ! Local variables:
         real    :: fact1(nsfc), fact2(nsfc)
         integer :: v, s, ndx, icg, kw, ks

         REAL :: VDEP_AE_SCALED(N_AE_DEP_SPC)
         REAL :: VDEP_AERO     (N_AE_DEP_SPC,NSFC)

         REAL :: VDEP_GAS_SCALED(N_UNIQUE_GDEPV)
         REAL :: VDEP_GAS       (N_UNIQUE_GDEPV,NSFC)

         do ks = 1, nsfc
            kw = lpw(iw) + ks - 1

            fact1(ks) = real(rho(kw,iw) * dtlm) * frac_sfc(ks,iw) * arw0(iw)
            fact2(ks) = 0.95 * real(volt(kw,iw) * rho(kw,iw))
         enddo


         if (N_AE_DEPV > 0) then

            call aero_depv ( iw, nsfc, vdep_aero )

            ! For diagnostic output
            adepv(:,iw) = 0.0

            do ks = 1, nsfc
               if (frac_sfc(ks,iw) > 1.e-5) then
                  kw = lpw(iw) + ks - 1

                  do s = 1, N_AE_DEP_SPC
                     vdep_ae_scaled(s) = - min(vdep_aero(s,ks) * fact1(ks), fact2(ks))

                     adepv(s,iw) = adepv(s,iw) + vdep_aero(s,ks) * frac_sfc(ks,iw)
                  enddo

                  do v = 1, n_ae_depv
                     NDX = AE_DEPV_SUR(V)
                     ICG = AE_DEPV_MAP(V) + AE_STRT - 1
                     sxfer_ae(ks,iw,v) = vdep_ae_scaled(ndx) * AE_DEPV_FAC(v) * cgrid(kw,iw,icg)
                  enddo

               endif
            enddo

         endif

         if (N_GAS_DEPV > 0) then

            call m3dry( iw, abflux, sfc_hono, nsfc, vdep_gas )

            ! For diagnostic output
            depvel_gas(:,iw) = 0.0

            do ks = 1, nsfc
               if (frac_sfc(ks,iw) > 1.e-5) then
                  kw = lpw(iw) + ks - 1

                  do s = 1, N_UNIQUE_GDEPV
                     vdep_gas_scaled(s) = - min(vdep_gas(s,ks) * fact1(ks), fact2(ks))

                     depvel_gas(s,iw) = depvel_gas(s,iw) + vdep_gas(s,ks) * frac_sfc(ks,iw)
                  enddo

                  do v = 1, n_gc_depv
                     NDX = GAS_DEPV_SUR( v )
                     ICG =  GC_DEPV_MAP( v )
                     sxfer_gc(ks,iw,v) = vdep_gas_scaled(ndx) * GC_DEPV_FAC(v) * cgrid(kw,iw,icg)
                  enddo

                  do v = 1, n_nr_depv
                     NDX = GAS_DEPV_SUR( v + N_GC_DEPV )
                     ICG = NR_STRT - 1 + NR_DEPV_MAP ( v )
                     sxfer_nr(ks,iw,v) = vdep_gas_scaled(ndx) * NR_DEPV_FAC(v) * cgrid(kw,iw,icg)
                  enddo

               endif
            enddo

         endif

         END SUBROUTINE GET_DEPV1

      END MODULE DEPV_DEFN
