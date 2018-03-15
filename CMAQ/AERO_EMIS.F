
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

! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /project/yoj/arc/CCTM/src/aero/aero6/AERO_EMIS.F,v 1.11 2012/01/19 13:20:08 yoj Exp $

! what(1) key, module and SID; SCCS file; date and time of last delta:
! %W% %P% %G% %U%

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE AERO_EMIS

!  Emissions data and code required for the modal aerosol module in CMAQ
!  Based on original codes by Dr. Francis S. Binkowski and J. Young
 
!  Dependent Upon:  NONE
 
!  Revision History:
 
!   30 Aug 01 J.Young:  dyn alloc - Use HGRD_DEFN
!   09 Oct 03 J.Gipson: added MW array for AE emis species to module contents
!   31 Jan 05 J.Young:  dyn alloc - establish both horizontal & vertical
!                       domain specifications in one module, GRID_CONF
!   26 Apr 05 P.Bhave:  removed code supporting the "old type" of emission 
!                        files that had unspeciated PM10 and PM2.5 only
!                       removed need for 'AERO_SPC.EXT' by declaring the 
!                        required variables locally
!   13 Jun 05 P.Bhave:  added vars needed for sea-salt emission processing
!                       inherit N_AE_EMIS,AE_EMIS,AE_EMIS_MAP from AE_EMIS.EXT
!                       moved RHO* parameters from RDEMIS_AE to this module
!                        for use by SSEMIS routine
!   24 Aug 07 J.Young:  Modified to enable in-line plume rise calculation for
!                       3D pt source emissions. Distinguish between PM (primary,
!                       unspeciated, file data) and AE (model speciated). Re-
!                       named RDEMIS_AE to GET_AERO_EMIS.
!   11 Apr 08 J.Kelly:  added code to emit coarse surface area
!    4 Jan 10 J.Young:  restructure; eliminate ref to older AERO versions
!   21 Feb 10 J.Young:  move sea salt emissions to its own module (SSEMIS)
!   23 Apr 10 J.Young:  replace include files with mechanism namelists
!   30 Apr 10 J.Young:  update to use aero_reeng by Steve Howard, Prakash Bhave,
!                       Jeff Young, and Sergey Napelenok
!   23 Jul 10 D.Wong:   remove CLOSE3 and BARRIER
!   24 Feb 11 J.Young:  Reorganized module with initialization and timestepping
!                       procedures
!   25 Feb 11 J.Young:  add windblown dust module
!   25 Mar 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
!   11 May 11 D.Wong: incorporated twoway model implementation
!   18 Aug 11 David Wong: In the merge inline point source PM species calculation,
!                         arrays EMBUFF and PMEMIS_PT have incorrect index values
!   17 Apr 13 J.Young: replace "SPFC ASO4" (found by Havala Pye) with "SPFC_ASO4"
!   07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module. 
!   23 Dec 15 J.Young: turn off default windblown dust option

      USE AERO_DATA,  ONLY: N_MODE
      use cgrid_defn, only: vdemis_gc, vdemis_nr, vdemis_ae
       
      IMPLICIT NONE
      
! aerosol emissions: [ppmv/s] for mass & number spcs, [m2/mol/s] for surface area spcs

      PUBLIC AERO_EMIS_INIT, GET_AERO_EMIS
      PRIVATE

! Emission rate of all aerosol species interpolated to current time
      INTEGER, ALLOCATABLE, SAVE :: EM_MAP  ( :,: )     ! index in namelist table
      REAL,    ALLOCATABLE, SAVE :: EM_FAC  ( :,: )     ! Namelist emissions scaling factor
      INTEGER, ALLOCATABLE, SAVE :: NUM_MAP ( : )       !
      REAL,    ALLOCATABLE, SAVE :: NUM_FAC ( : )       !
      INTEGER, ALLOCATABLE, SAVE :: SRF_MAP ( : )       !
      REAL,    ALLOCATABLE, SAVE :: SRF_FAC ( : )       !
      INTEGER, ALLOCATABLE, SAVE :: SS_MAP  ( :,: )     !
      INTEGER, ALLOCATABLE, SAVE :: DUST_MAP( :,: )     !

! EMBUFF     mass emission rate input global array
! PM_EM      mass emission rate input localized to one grid cell
! EM_RATE    scaled PM_EM
! EM_PART    particulate mass emission rate in one grid cell
! EM_MAP                 mass emission rate index in VDEMIS_AE (CGRID)
! EM_FAC                 mass emission Namelist scaling factor for VDEMIS_AE (CGRID)
! NUM_MAP                number emission rate index in VDEMIS_AE (CGRID)
! NUM_FAC
! SRF_MAP                sfc area emission rate index in VDEMIS_AE (CGRID)
! SRF_FAC
! SS_MAP     sea salt emission rate index in SSOUTM
! DUST_MAP   windblown dust emission rate index in DUSTOUTM

! Factor for converting aerosol emissions from input units ...
      REAL,    SAVE :: CONVEM_PM_MASS ! into [ug/sec]

! *** PVB: temporary hack to speciate PMC
      REAL,    SAVE :: SPFC_ACORS
      REAL,    SAVE :: SPFC_ASO4
      REAL,    SAVE :: SPFC_ANO3
      REAL,    SAVE :: SPFC_ACL
      REAL,    SAVE :: SPFC_AH2O

! Variables for converting emission rates into molar-mixing-ratio units
      REAL,    PARAMETER :: GPKG = 1.0E+03     ! g/kg
      REAL,    PARAMETER :: MGPG = 1.0E+06     ! ug/g

! For handling vapor-phase sulfuric acid emissions
      INTEGER, SAVE :: VSULF               ! index to H2SO4 in VDEMIS array
      INTEGER, SAVE :: VPSO4               ! index to PSO4 in EM_RATE array

! Miscellaneous variables
      LOGICAL, SAVE :: WB_DUST   ! flag for windblown dust option
      INTEGER, SAVE :: LOGDEV

      integer, allocatable, save :: edgar_ae_map(:)

      CONTAINS

!-----------------------------------------------------------------------

         FUNCTION  AERO_EMIS_INIT ( JDATE, JTIME, TSTEP, NEMIS_AE,
     &                              nvars3d_emis, vname3d_emis ) RESULT ( SUCCESS)

!  Revision History:
 
!   30 Aug 01 J.Young:  dynamic allocation - Use INTERPX
!   29 Jul 03 P.Bhave:  added compatibility with emission files that contain 
!                       PM10, PEC, POA, PNO3, PSO4, and PMF, but do not 
!                       contain PMC
!   20 Aug 03 J.Young:  return aero emissions in molar mixing ratio, ppm units
!   09 Oct 03 J.Gipson: added MW array for AE emis species to module contents
!   01 Sep 04 P.Bhave:  changed MW for primary organics from 120 to 220 g/mol,
!                       to match MWPOA in subroutine ORGAER3.
!   31 Jan 05 J.Young:  dyn alloc - removed HGRD_ID, VGRID_ID, and COORD_ID 
!                       include files because those parameters are now 
!                       inherited from the GRID_CONF module
!   26 Apr 05 P.Bhave:  removed code supporting the "old type" of emission 
!                        files that had unspeciated PM10 and PM2.5 only
!                       removed need for 'AERO_CONST.EXT' by declaring the
!                        required variables locally
!                       simplified the CONVM, CONVN, CONVS calculations
!                       updated and enhanced in-line documentation
!   03 May 05 P.Bhave:  fixed bug in the H2SO4 unit conversion, initially
!                        identified by Jinyou Liang of CARB
!   13 Jun 05 P.Bhave:  calculate sea-salt emissions; execute if MECHNAME = AE4
!                        read input fields from new OCEAN_1 file
!                        read extra input fields from MET_CRO_2D and MET_CRO_3D
!                        write diagnostic sea-salt emission file
!                        added TSTEP to call vector for diagnostic output file
!                       inherit MWs from AE_SPC.EXT instead of hardcoding
!                       find pointers to CGRID indices instead of hardcoding
!   08 Mar 07 P.Bhave&   added capability for emission files that contain 
!             S.Roselle:  POC or POA
!   30 Jan 08 P.Bhave:  added compatibility with AE5 mechanisms
!   23 Mar 08 J.Young:  modifications to allow for in-line point source emissions
!   11 Apr 08 J.Kelly:  added code to emit coarse surface area
!   09 Sep 08 P.Bhave:  backward compatibility with AE4 mechanisms
!   20 Feb 10 J.Young:  move ssemis out to its own F90 module
!   24 Feb 11 J.Young:  add windblown dust emissions option
!   25 Mar 11 S.Roselle: Replaced I/O API include files with UTILIO_DEFN
!   07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
!   17 Sep 14 K.Fahey:  Changed geometric mean diameter and geometric
!                       standard deviation of emitted particles according to 
!                       Elleman and Covert (2010)
!   15 Apr 16 J.Young: Use aerosol factors from the AERO_DATA module's named constants;
!                      Moved K.Fahey's mods to geometric mean diameter and standard
!                      deviation to the AERO_DATA module
 
!  References:
!    CRC76,        "CRC Handbook of Chemistry and Physics (76th Ed)",
!                   CRC Press, 1995
!    Elleman & Covert, "Aerosol size distribution modeling with the Community
!                   Multiscale Air Quality modeling system in the Pacific
!                   Northwest: 3. Size distribution of particles emitted
!                   into a mesoscale model", J. Geophys. Res., Vol 115,
!                   No D3, doi:10.1029/2009JD012401, 2010
!    Hobbs, P.V.   "Basic Physical Chemistry for the Atmospheric Sciences",
!                   Cambridge Univ. Press, 206 pp, 1995.
!    Snyder, J.P.  "Map Projections-A Working Manual", U.S. Geological Survey
!                   Paper 1395 U.S.GPO, Washington, DC, 1987.
!    Binkowski & Roselle  Models-3 Community Multiscale Air Quality (CMAQ)
!                   model aerosol component 1: Model Description.  
!                   J. Geophys. Res., Vol 108, No D6, 4183 
!                   doi:10.1029/2001JD001409, 2003

         USE RXNS_DATA           ! chemical mechanism data
         USE CGRID_SPCS          ! CGRID mechanism species
         USE SSEMIS              ! sea-salt emissions
         USE DUST_EMIS           ! wind-blown dust emissions
         USE AERO_DATA
         USE PRECURSOR_DATA
         USE UTILIO_DEFN
         USE CONST_DATA          ! physical and mathematical constants

! Arguments:

         INTEGER, INTENT( IN ) :: JDATE      ! current model date, coded YYYYDDD
         INTEGER, INTENT( IN ) :: JTIME      ! current model time, coded HHMMSS
         INTEGER, INTENT( IN ) :: TSTEP      ! time step vector (HHMMSS)
                                             ! TSTEP(1) = local output step
         INTEGER, INTENT( IN ) :: NEMIS_AE   ! number of aerosol emissions species

         integer,       intent(in) :: nvars3d_emis
         character(24), intent(in) :: vname3d_emis(:)

         LOGICAL SUCCESS

! External Functions:
         INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

! Local Variables:
         LOGICAL :: ZERO_PCSOA

! Miscellaneous variables
         INTEGER STATUS                   ! ENV..., ALLOCATE status
         CHARACTER( 16 ), SAVE :: PNAME = 'AERO_EMIS_INIT  '
         CHARACTER( 16 ) :: VNAME         ! temp var for species names
         INTEGER L, N, S, V, IPOA          ! Loop indices
         INTEGER :: IDX
         CHARACTER( 96 ) :: XMSG = ' '
         INTEGER, EXTERNAL :: SETUP_LOGDEV

! ----------------------------------------------------------------------

         LOGDEV = SETUP_LOGDEV()

! *** Map data modules
         CALL MAP_AERO()
         CALL MAP_PRECURSOR()
!        CALL MAP_PMEMIS()  ! Called in emis_init

! *** set up for sea-salt emission processing
         IF ( .NOT. SSEMIS_INIT( JDATE, JTIME, TSTEP ) ) THEN
            XMSG = 'Failure initializing sea-salt emission processing'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

! *** windblown dust emissions?
         wb_dust = .true.

! *** set up for dust emission processing
         IF ( WB_DUST ) THEN
            IF ( .NOT. DUST_EMIS_INIT( JDATE, JTIME, TSTEP ) ) THEN
               XMSG = 'Failure initializing dust emission processing'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF
         END IF

! *** Find location of sulfuric acid vapor in VDEMIS array
         VNAME = 'SULF'
         IDX = INDEX1 ( VNAME, N_GC_EMIS, GC_EMIS )
         IF ( IDX .NE. 0 ) THEN
            VSULF = IDX  ! index for vapor-phase H2SO4 emissions
         ELSE
            XMSG = 'Could not find ' // VNAME // 'in gas table'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
         END IF

         ALLOCATE ( EM_MAP  ( N_AEROSPC, N_MODE ),
     &              EM_FAC  ( N_AEROSPC, N_MODE ),
     &              NUM_MAP (            N_MODE ),
     &              NUM_FAC (            N_MODE ),
     &              SRF_MAP (            N_MODE ),
     &              SRF_FAC (            N_MODE ),
     &              SS_MAP  ( N_AEROSPC, N_MODE ),
     &              DUST_MAP( N_AEROSPC, N_MODE ), STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = '*** memory allocation failed for '
     &           // 'EM_PART, EM_MAP, EM_FAC, NUM_MAP, NUM_FAC, SRF_MAP, SRF_FAC, SS_MAP or DUST_MAP'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

! *** Map indices of emission species
      write( logdev,* ) '=== aerospc em_map'
         EM_MAP = 0  
         DO S = 1, N_AEROSPC
            IF ( AEROSPC( S )%EMIS .NE. ' ' ) THEN
               ! Verify that aerospc( S )%emis is a pm emis
               DO N = 1, N_MODE
                  IF ( .NOT. AERO_MISSING(S,N) .AND.
     &                 AEROSPC( S )%EMIS_SPLIT( N ) .NE. 0.0 ) THEN
                     EM_MAP( S,N ) = INDEX1( AEROSPC( S )%NAME( N ), N_AE_EMIS, AE_EMIS )
                     write( logdev,* ) ' ', s, n, aerospc( s )%name( n ), em_map( s,n )
!                    if ( em_map( s,n ) .eq. 0 ) then
!                       write( logdev,* ) '*** aerospc em_map: ', s, n, em_map( s,n )
!                    end if
                  END IF
               END DO
            END IF
         END DO

! *** set across-the-board adjustment factors (from Namelist file)
!     Strategy: first check Namelist Emis_Fac for antrhop. emis, then check for non-zero split factor
      write( logdev,* ) 'modified aerospc em_fac'
         EM_FAC = 1.0 ! array assignment
         DO S = 1, N_AEROSPC
            IF ( AEROSPC( S )%EMIS .NE. ' ' ) THEN
               DO N = 1, N_MODE
                  IF ( .NOT. AERO_MISSING(S,N) ) THEN
                     L = INDEX1( AEROSPC( S )%NAME( N ), N_AE_EMIS, AE_EMIS )
                     IF ( L .GT. 0 ) THEN
                        EM_FAC( S,N ) = AE_EMIS_FAC( L )
                        ! Test for POA Species and Aggregate Volatility
                        ! Parameters
                        IF ( AEROSPC( S )%EMIS .EQ. 'POC' ) THEN
                            DO IPOA = 1,N_VBS_BIN
                               IF ( AEROSPC( S )%NAME( N )(2:6) .EQ. POA_NAME( IPOA ) ) THEN
                                  EM_FAC( S,N ) = EM_FAC( S,N ) * POA_OP_VF( IPOA ) * POA_AMF( IPOA )
                               ENDIF
                            END DO
                        ENDIF

                        if ( em_fac( s,n ) .ne. 1.0 ) then
                           write( logdev,* ) ' ', s, n, aerospc( s )%name( n ), em_fac( s,n )
                        end if
                     END IF
                  END IF
               END DO
            ELSE
               DO N = 1, N_MODE
                  IF ( AEROSPC( S )%EMIS_SPLIT( N ) .NE. 0.0 ) THEN
                     L = INDEX1( AEROSPC( S )%NAME( N ), N_AE_EMIS, AE_EMIS )
                     IF ( AE_EMIS( L ) .NE. ' ' ) THEN
                        EM_FAC( S,N ) = AE_EMIS_FAC( L )
                        if ( em_fac( s,n ) .ne. 1.0 ) then
                           write( logdev,* ) ' ', s, n, aerospc( s )%name( n ), em_fac( s,n )
                        end if
                     END IF
                  END IF
               END DO
            END IF
         END DO


         IF ( INDEX( MECHNAME,"NVPOA" ) .EQ. 0 ) THEN
           ! Semivolatile POA  
! *** Let the user set the pcVOC species (which makes pcSOA) emissions
! *** scaling to 0.0
           ! Retrieve the Environment Variable allowing the User to Toggle
           ! the pcVOC emissions
           ZERO_PCSOA = .FALSE.         ! default

           V = INDEX1( 'PCVOC', N_GC_SPC, GC_SPC )
           IF ( V .LT. 1 ) THEN
             WRITE(LOGDEV,99901)'PCVOC','PCVOC',TRIM( MECHNAME )
             SUCCESS = .FALSE.
           END IF

           IF ( ZERO_PCSOA ) THEN
             ! Find the index for the pcVOC species in the emissions
             ! scale-factor vector
             IF ( V .GT. 0 ) THEN
                S = FINDEX( V, N_GC_EMIS, GC_EMIS_MAP )
                IF ( S .GT. 0 )THEN
                   GC_EMIS_FAC( S ) = 0.0
                ELSE
                   WRITE(LOGDEV,99902)'PCVOC'
                END IF
             END IF
           ENDIF

! *** Assign Volaitlity Distribution to Primary Orgnaic Vapor Species
           ! Find the Index for Each Vapor and alter the emissions scaling
           ! to reflect the parameters in AERO_DATA. The GC_EMIS_FAC field
           ! from the GC namelist should account for the conversion from
           ! g/s (POC) to mol/s (Vapors). That factor is then modified
           ! by the volatility distribution itself.
           DO IPOA = 1,N_VBS_BIN
             VNAME = 'V' // POA_NAME( IPOA )
             V = INDEX1( VNAME, N_GC_SPC, GC_SPC )
             S = FINDEX( V, N_GC_EMIS, GC_EMIS_MAP )
             GC_EMIS_FAC( S ) = GC_EMIS_FAC( S ) * POA_OP_VF( IPOA ) * 
     &                          (1.0 - POA_AMF( IPOA ) ) / GC_MOLWT( V )
           END DO
         ENDIF

! *** Map aerosol number emission species and set across-the-board adjustment factors
         NUM_MAP = 0   ! array assignment
         NUM_FAC = 1.0 ! array assignment
         DO N = 1, N_MODE
            NUM_MAP( N ) = INDEX1( AEROMODE( N )%NUM_NAME, N_AE_EMIS, AE_EMIS )
            NUM_FAC( N ) = AE_EMIS_FAC( NUM_MAP( N ) )
#ifdef Verbose1
            write( logdev,* ) 'num_map: ', n, num_map( n ), aeromode( n )%num_name, num_fac( n )
#endif
            IF ( NUM_MAP( N ) .EQ. 0 ) THEN
               XMSG = 'Could not find ' // TRIM( AEROMODE( N )%NUM_NAME )
     &              // ' in AE_EMIS table'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF
         END DO

! *** Map aerosol surface emission species and set across-the-board adjustment factors
         SRF_MAP = 0   ! array assignment
         SRF_FAC = 1.0 ! array assignment
         DO N = 1, N_MODE
            SRF_MAP( N ) = INDEX1( AEROMODE( N )%SRF_NAME, N_AE_EMIS, AE_EMIS )
            SRF_FAC( N ) = AE_EMIS_FAC( SRF_MAP( N ) )
#ifdef Verbose1
            write( logdev,* ) 'srf_map: ', n, srf_map( n ), aeromode( n )%srf_name, srf_fac( n )
#endif
            IF ( SRF_MAP( N ) .EQ. 0 ) THEN
               XMSG = 'Could not find ' // TRIM( AEROMODE( N )%SRF_NAME )
     &              // ' in AE_EMIS table'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF
         END DO

! *** PVB: speciate PMC (from AERO_DATA)
         N = N_MODE
         SPFC_ASO4  = REAL( ACORSEM_ASO4_FAC ) * AEROSPC( ASO4_IDX )%EMIS_SPLIT( N )
         SPFC_ANO3  = REAL( ACORSEM_ANO3_FAC ) * AEROSPC( ANO3_IDX )%EMIS_SPLIT( N )
         SPFC_ACL   = REAL( ACORSEM_ACL_FAC ) * AEROSPC( ACL_IDX )%EMIS_SPLIT( N )
         SPFC_AH2O  = REAL( ACORSEM_AH2O_FAC ) * AEROSPC( AH2O_IDX )%EMIS_SPLIT( N )
         SPFC_ACORS = REAL( ACORS_RENORM ) * AEROSPC( ACORS_IDX )%EMIS_SPLIT( N )

! *** Map sea salt emission species
         SS_MAP = 0   ! array assignment
         DO S = 1, NSSSPC
            IDX = SSSPC_MAP( S )
#ifdef Verbose1
            write( logdev,* ) 'ssspc_map: ', s, idx, ssspc_map( s )
#endif
            DO N = 2, N_MODE   ! add SEASALT indices to EM_MAP
               IF ( .NOT. AERO_MISSING( IDX,N ) ) THEN
                  EM_MAP( IDX,N ) = INDEX1( AEROSPC( IDX )%NAME( N ), N_AE_EMIS, AE_EMIS )
                  SS_MAP( IDX,N ) = S
#ifdef Verbose1
               write( logdev,* ) 'ss_map:    ', s, n, idx, ss_map( idx,n ), aerospc( idx )%name( n )
#endif
               END IF
            END DO
         END DO

         IF ( WB_DUST ) THEN
! *** Map windblown dust emission species
            DUST_MAP = 0   ! array assignment
            DO N = 2, N_MODE
               DO S = 1, NDUST_SPC
                  IF ( DUST_SPC( S )%NAME( N ) .NE. ' ' ) THEN
                     IDX = FINDAERO( DUST_SPC( S )%NAME( N ), .FALSE. )
                     IF ( IDX .NE. 0 ) THEN
                        EM_MAP( IDX,N ) = INDEX1( AEROSPC( IDX )%NAME( N ), N_AE_EMIS, AE_EMIS )
                        DUST_MAP( IDX,N ) = S
                     ELSE
                        WRITE( XMSG, '(A43,A23,A8,A11,A62,A53)' ), 
     &                         'Cannot find aerosol species to accomodate ',
     &                         'dust emission species: ', DUST_SPC( S )%NAME( N ), '. You must ',
     &                         'provide a prognostic variable for this species in both the AE ',
     &                         'namelist and the AERO_DATA table in order to proceed.'
                        CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
                        SUCCESS = .FALSE.
                     END IF
#ifdef Verbose1
                  write( logdev,* ) 'dust_map: ', n, s, dust_map( idx,n ), aerospc( idx )%name( n )
#endif
                  END IF
               END DO
            END DO
         END IF

#ifdef Verbose1
         do s = 1, n_aerospc
            do n = 1, n_mode
               if ( em_map( s,n ) .ne. 0 ) then
                  write( logdev,* ) 'em_map: ', n, s, em_map( s,n ), aerospc( s )%name( n ), em_fac( s,n )
               end if
            end do
         end do
#endif

! *** get index of SO4 emissions
         VPSO4 = 0
         DO S = 1, N_EMIS_PM
            IF ( PMEM_MAP( S ) .EQ. ASO4_IDX ) VPSO4 = S
#ifdef Verbose1
               write( logdev,* ) 'vpso4: ', s, pmem_map( s ), aso4_idx, vpso4
#endif
         END DO

         IF ( VPSO4 .EQ. 0 ) THEN
            XMSG = 'Cannot define index for SO4 emissions: '
     &           // 'aerospc%emis='' '' in AERO_DATA.F'
            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
         END IF

! *** Get scaling factor for converting aerosol emissions from
!     their input units to [ug/s] and then to [ug/m2/s] using layer-
!     specific grid-cell volume.  Calling DESC3( EMIS_1 ) assumes EMIS_1
!     already opened.

         IF ( PMEM_UNITS .EQ. 'G/S' .OR.
     &        PMEM_UNITS .EQ. 'g/s' ) THEN
            CONVEM_PM_MASS = MGPG                  ! (g/s) -> (ug/s)
         ELSE IF ( PMEM_UNITS .EQ. 'KG/HR' .OR.
     &             PMEM_UNITS .EQ. 'kg/hr' ) THEN
            CONVEM_PM_MASS = GPKG * MGPG / 3600.0  ! (kg/hr) -> (ug/s)
         ELSE
            CONVEM_PM_MASS = 1.0
            XMSG = 'Units incorrect on EMIS_1'
!           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
         END IF

! Mapping for edgar emissions

         allocate(edgar_ae_map(n_emis_pm))

         do s = 1, n_emis_pm
            n = index1( pmem_map_name(s), nvars3d_emis, vname3d_emis )
            if (n > 0) then
               edgar_ae_map(s) = n
            else
               xmsg = trim( pmem_map_name(s) ) // ' was expected to be in emissions file.'
               call m3exit( pname, 0, 0, xmsg, xstat2 )
            endif
         enddo

         SUCCESS = .TRUE.; RETURN

99901     FORMAT(/'FATAL: GC namelist used does not include species, ', A, 
     &         / 'Check this namelist contains the same information '
     &         / 'for ', A,' as the repository version of the GC namelist for '
     &         / 'the ', A, ' mechanism.' )
99902     FORMAT(/'WARNING: GC namelist line for ', A, ' does not include an emissions '
     &         / 'surrogate. Simulation cannot evaluate the effect from zeroing out its '
     &         / 'emissions.')
         END FUNCTION  AERO_EMIS_INIT

!-----------------------------------------------------------------------

         subroutine get_aero_emis ( mrl, nvars3d_emis, vname3d_emis, edgar_kmax )

         use ssemis,     only: get_ssemis
         use dust_emis,  only: get_dust_emis
         use mem_ijtabs, only: jtab_w, jtw_prog

         implicit none

         integer,       intent(in) :: mrl
         integer,       intent(in) :: nvars3d_emis
         character(24), intent(in) :: vname3d_emis(:)
         integer,       intent(in) :: edgar_kmax(:)
         integer                   :: j, iw

! *** Get sea salt emissions
         CALL GET_SSEMIS ( )

! *** Get dust emissions
         IF ( WB_DUST ) THEN
            CALL GET_DUST_EMIS ( )
         END IF

         !$omp parallel do private(j,iw)
         do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

            call get_aero_emis1( iw, nvars3d_emis, vname3d_emis, edgar_kmax )

         enddo
         !$omp end parallel do

         end subroutine get_aero_emis

!-----------------------------------------------------------------------

         SUBROUTINE GET_AERO_EMIS1 ( iw, nvars3d_emis, vname3d_emis, edgar_kmax )

!  GET_AERO_EMIS reads PM emissions from gridded input file and
!  optionally from plume rise calculated point source emissions (via VDEMIS
!  array) and converts into molar-mixing-ratio units, for use in the vertical
!  diffusion routines
 
!  Key Subroutines/Functions Called: SSEMIS_INIT, GET_SSEMIS, DUST_INIT, GET_DUST_EMIS

!  Revision History:

!   18 Aug 11 David Wong: In the merge inline point source PM species calculation,
!                         arrays EMBUFF and PMEMIS_PT have incorrect index values
 
         USE CGRID_SPCS
!        USE PT3D_DEFN,    ONLY: PT3DEMIS, PMEMIS_PT  ! from plume rise
!        USE PTMAP,        ONLY: N_SPC_PTPM, PTPM_MAP ! defines emissions file(s) PM mapping
         USE AERO_DATA
         USE PRECURSOR_DATA
         USE SSEMIS
         USE UTILIO_DEFN
         USE DUST_EMIS
         use CONST_DATA,   only: mwairkg, avo
         use mem_grid,     only: lpw, volti
         use mem_ijtabs,   only: itab_w
         use mem_basic,    only: rho
         use mem_sea,      only: itab_ws
         use mem_leaf,     only: itab_wl
         use edgar_42_emis,only: edgar42_emis

! Arguments:
         integer,       intent(in) :: iw
         integer,       intent(in) :: nvars3d_emis
         character(24), intent(in) :: vname3d_emis(:)
         integer,       intent(in) :: edgar_kmax(:)

! Variables for converting emission rates into molar-mixing-ratio units
         REAL, PARAMETER :: RAVO = 1.E6 / AVO   ! reciprocal of Avogadro
         REAL      :: CONVM     ! conversion factor for mass emissions [m3/mol]
         REAL      :: CONVN     ! conversion factor for number emissions [1e6*m3]
         REAL      :: CONVS     ! conversion factor for surface area emissions [m3/mol]
         REAL      :: GSFAC     ! into [ug/m3/sec]

! Third moment emissions rates [m3/m3/s]
         REAL      :: EMISM3( N_MODE )

! Number emissions rates [1/m3/s]
         REAL      :: EM_NUM( N_MODE )

! Surface area emission rates [m2/m3/s]
         REAL      :: EM_SRF( N_MODE )

! Mass emissions rates
         real      :: pm_em  ( n_emis_pm )
         real      :: em_rate( n_emis_pm )
         real      :: em_part( n_aerospc, n_mode )

! Miscellaneous variables
         INTEGER :: N, S
         integer :: k, iws, iwl, jws, jwl
         INTEGER :: IDX
         real    :: emis

         do k = lpw(iw), edgar_kmax(iw)

! *** Store aerosol emission rates for this grid cell in a 1D array
            DO S = 1, N_EMIS_PM
               PM_EM( S ) = edgar42_emis( k,iw,edgar_ae_map(s) )
            END DO

! *** Calculate scaling factor for converting mass emissions into [ug/m3/s]
!     note: RJACM converts grid heights from sigma coordinates to meters
!     Also calculate scaling factors for converting to molar-mixing-ratio units
!           GSFAC = CONVEM_PM( L ) * Met_Data%RJACM( C,R,L )
            GSFAC = CONVEM_PM_MASS * volti(k,iw) ! [ug/g/m^3]
            CONVM = MWAIRKG / rho(k,iw)          ! [m3/mol]
            CONVN = CONVM * RAVO                 ! 10^6 * [m3]
            CONVS = CONVM                        ! [m3/mol]

! *** Calculate speciated mass emission rates for fine aerosol [ug/m3/s]
            DO S = 1, N_EMIS_PM
               EM_RATE( S ) = PM_EM( S ) * GSFAC
            END DO

! *** Extract H2SO4 vapor emission rate from the VDEMIS array, add it to
!     the fine-PM sulfate emissions, and remove it from the gas emissions
            EM_RATE( VPSO4 ) = EM_RATE( VPSO4 )
     &                       + VDEMIS_GC( k,iw,VSULF ) * PRECURSOR_MW( SULF_IDX ) / CONVM
            VDEMIS_GC( k,iw,VSULF ) = 0.0

! *** Calculate emissions rate for third moments [m3/m3/s] of each mode 
!     (excluding sea salt), as in Equation 7a of Binkowski & Roselle (2003).
            EM_PART = 0.0       ! array assignment
            EMISM3  = 0.0       ! array assignment
            DO S = 1, N_EMIS_PM
               IDX = PMEM_MAP( S )
               DO N = 1, N_MODE - 1
                  emis = AEROSPC( IDX )%EMIS_SPLIT( N ) * EM_RATE( S ) * EM_FAC( IDX,N )
                  EM_PART( IDX,N ) = EM_PART( IDX,N ) + emis
                  IF ( .NOT. AEROSPC( IDX )%TRACER ) THEN
                     EMISM3( N ) = EMISM3( N ) + emis * aerospc_m3fac( idx )
                  END IF
               END DO

               ! Apply Emissions to Coarse Aerosol Mode                    
               N = N_MODE
               IF ( IDX .EQ. ASO4_IDX ) THEN
                  EM_PART( IDX,N ) = SPFC_ASO4 * EM_RATE( S ) * EM_FAC( IDX,N )
               ELSE IF ( IDX .EQ. ANO3_IDX ) THEN
                  EM_PART( IDX,N ) = SPFC_ANO3 * EM_RATE( S ) * EM_FAC( IDX,N )
               ELSE IF ( IDX .EQ. ACL_IDX ) THEN
                  EM_PART( IDX,N ) = SPFC_ACL * EM_RATE( S ) * EM_FAC( IDX,N )
               ELSE IF ( IDX .EQ. AH2O_IDX ) THEN
                  EM_PART( IDX,N ) = SPFC_AH2O * EM_RATE( S ) * EM_FAC( IDX,N )
               ELSE IF ( IDX .EQ. ACORS_IDX ) THEN
                  EM_PART( IDX,N ) = SPFC_ACORS * EM_RATE( S ) * EM_FAC( IDX,N )
               ELSE
                  EM_PART( IDX,N ) = AEROSPC( IDX )%EMIS_SPLIT( N )
     &                             * EM_RATE( S ) * EM_FAC( IDX,N )
               END IF
               IF ( .NOT. AEROSPC( IDX )%TRACER ) THEN
                  EMISM3( N ) = EMISM3( N ) + EM_PART( IDX,N ) * aerospc_m3fac( idx )
               END IF
            END DO

! *** Calculate the number emissions rate for each mode [1/m3/s], using 
!     Equation 7b of Binkowski & Roselle (2003).
            DO N = 1, N_MODE
               EM_NUM( N ) = FACNUM( N ) * EMISM3( N )
            END DO

! *** Calculate the surface area emissions rate for the fine modes [m2/m3/s],
!     using Equation 7c of Binkowski & Roselle (2003).  Multiplying by PI 
!     converts 2nd moment to surface area.
            DO N = 1, N_MODE
               EM_SRF( N ) = FACSRF( N ) * EMISM3( N )
            END DO

! *** Convert emission rates into molar-mixing-ratio units, as required by
!     the vertical diffusion routines.  Mass and number emissions are 
!     converted to [ppmv/s].  Surface area emissions are converted to 
!     [m2/mol/s].  Save results in the VDEMIS_AE array.
            DO N = 1, N_MODE
               DO S = 1, N_AEROSPC
                  IF ( EM_MAP( S,N ) .GT. 0 ) THEN
                     vdemis_ae( k,iw,EM_MAP(S,N) ) = EM_PART( S,N )
     &                                             * CONVM * AEROSPC_MWI( S )
                  END IF
               END DO
               VDEMIS_AE( k,iw,NUM_MAP( N ) ) = EM_NUM( N ) * NUM_FAC( N ) * CONVN
               VDEMIS_AE( k,iw,SRF_MAP( N ) ) = EM_SRF( N ) * SRF_FAC( N ) * CONVS
            END DO

         enddo                  ! loop on k

! *** ADD SEA SALT EMSISSIONS

         do jws = 1, itab_w(iw)%nsea
            iws = itab_w(iw)%isea(jws)
            k   = itab_ws(iws)%kw

! *** Get grid cell values for sea salt emission rates
!     Convert units of SSOUTM from [g/m3/s] to [ug/m3/s], and
!     Convert emission rates into molar-mixing-ratio units, as required by
!     the vertical diffusion routines.  Mass and number emissions are 
!     converted to [ppmv/s].  Surface area emissions are converted to 
!     [m2/mol/s].  Save results in the VDEMIS_AE array.

            CONVM = MWAIRKG / rho(k,iw)  ! [m3/mol]
            CONVN = CONVM * RAVO         ! 10^6 * [m3]
            CONVS = CONVM                ! [m3/mol]

            DO N = 1, N_MODE
               DO S = 1, N_AEROSPC
                  IF ( EM_MAP( S,N ) .GT. 0 .and. SS_MAP( S,N ) .NE. 0 ) THEN
                     VDEMIS_AE( k,iw,EM_MAP( S,N ) ) = VDEMIS_AE( k,iw,EM_MAP( S,N ) )
     &                    + SSOUTM( SS_MAP( S,N ),N,iws ) * MGPG
     &                    * CONVM * AEROSPC_MWI( S )
                  END IF
               END DO

               VDEMIS_AE( k,iw,NUM_MAP( N ) ) = VDEMIS_AE( k,iw,NUM_MAP( N ) )
     &                                        + SSOUTN( N,iws ) * NUM_FAC( N ) * CONVN
               VDEMIS_AE( k,iw,SRF_MAP( N ) ) = VDEMIS_AE( k,iw,SRF_MAP( N ) ) 
     &                                        + SSOUTS( N,iws ) * SRF_FAC( N ) * CONVS
            END DO

         enddo

! *** ADD DUST EMSISSIONS

         do jwl = 1, itab_w(iw)%nland
            iwl = itab_w(iw)%iland(jwl)
            k   = itab_wl(iwl)%kw

! *** Get grid cell values for sea salt emission rates
!     Convert units of SSOUTM from [g/m3/s] to [ug/m3/s], and
!     Convert emission rates into molar-mixing-ratio units, as required by
!     the vertical diffusion routines.  Mass and number emissions are 
!     converted to [ppmv/s].  Surface area emissions are converted to 
!     [m2/mol/s].  Save results in the VDEMIS_AE array.

            if (dust_em(iwl) > 1.e-25) then

               CONVM = MWAIRKG / rho(k,iw)  ! [m3/mol]
               CONVN = CONVM * RAVO         ! 10^6 * [m3]
               CONVS = CONVM                ! [m3/mol]

               DO N = 2, N_MODE
                  DO S = 1, N_AEROSPC
                     IF ( EM_MAP( S,N ) .GT. 0 .and. DUST_MAP( S,N ) .NE. 0 ) THEN
                        VDEMIS_AE( k,iw,EM_MAP( S,N ) ) = VDEMIS_AE( k,iw,EM_MAP( S,N ) )
     &                       + DUSTOUTM( DUST_MAP( S,N ),N,iwl ) * MGPG
     &                       * CONVM * AEROSPC_MWI( S )
                     END IF
                  END DO
                  VDEMIS_AE( k,iw,NUM_MAP( N ) ) = VDEMIS_AE( k,iw,NUM_MAP( N ) )
     &                                           + DUSTOUTN( N,iwl ) * NUM_FAC( N ) * CONVN
                  VDEMIS_AE( k,iw,SRF_MAP( N ) ) = VDEMIS_AE( k,iw,SRF_MAP( N ) ) 
     &                                           + DUSTOUTS( N,iwl ) * SRF_FAC( N ) * CONVS
               END DO

            endif
         enddo

         END SUBROUTINE GET_AERO_EMIS1

      END MODULE AERO_EMIS
