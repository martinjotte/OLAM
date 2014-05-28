
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
C $Header: /project/yoj/arc/CCTM/src/emis/emis/EMIS_DEFN.F,v 1.10 2011/10/21 16:10:46 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE EMIS_DEFN

C-----------------------------------------------------------------------
C Function: emissions interface to the chemistry-transport model

C Revision History:
C     28 Jul 2006 J.Young: initial implementation
C     18 Aug 2007 J.Young: move beis part to separate module; add plume rise
C     23 Sep 2009 B.Hutzell: modified algorithm that loads gas emissions from point
C                 sources into VDEMIS array to enable multi-use of an emission species
C     26 Jan 2010 J.Young: fix bug overwriting point source layer 1 NH3; inline rdemis
C     07 Jan 2011 B.Hutzell: updated for namelist definition of model species
C     16 Feb 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                            removed deprecated TRIMLEN
C      6 Apr 2011 J.Young, R.Pinder: add lightning NO emissions capability
C     11 May 2011 D.Wong: incorporated twoway model implementation

C-----------------------------------------------------------------------

!     USE GRID_CONF           ! horizontal & vertical domain specifications
      use max_dims,    only: pathlen
      use consts_coms, only: r8
      use cgrid_defn,  only: vdemis_gc, vdemis_nr, vdemis_ae

      IMPLICIT NONE

!     REAL,    ALLOCATABLE, target :: VDEMIS( :,:,: ) ! total emissions array
!#ifdef Verbose2
!      real,    allocatable, save :: semis( :,:,: )
!#endif

      INTEGER,              SAVE :: N_SPC_EMIS ! global emissions species
      INTEGER,              SAVE :: EMLAYS     ! no. of emis integration layers
                                               ! EMLAYS must be .LT. NLAYS
      INTEGER,              SAVE :: EMLYRS     ! no. of layers on base emis input file -

      PUBLIC N_SPC_EMIS, EMLAYS, EMLYRS, EMIS_INIT, GET_EMIS
      public emis_files, nemisfiles, ctotdate_emis, s1900_emis
      public iemisflg, iupdemis, iemiscyclic, iemisfile
      public nvars3d_emis, vname3d_emis, units3d_emis
      public sea_emisf,  sea_emisp,  sea_emisc
      public land_emisf, land_emisp, land_emisc
!#ifdef Verbose2
!      public semis
!#endif
      PRIVATE
                                           ! a 2d file if inline, otherwise a 3d file
      INTEGER, SAVE :: NEMIS_AE   ! no. of aero emis species

      LOGICAL, SAVE :: EM_TRAC = .FALSE. ! do tracer emissions?
!     REAL,    SAVE :: CNVTE      ! combined conversion factor
      REAL,    SAVE :: DT         ! TSTEP (output) in sec

      REAL,    ALLOCATABLE, SAVE :: BUFF( :,:,: )        ! read buffer

!      INTEGER, SAVE :: STRT_GC, FINI_GC !, STRT_AE, FINI_AE,
!      INTEGER, SAVE :: STRT_NR, FINI_NR !, STRT_TR, FINI_TR
      INTEGER, SAVE :: SDATE, STIME  ! scenario start date/time (beis)
      INTEGER       :: LDATE, LTIME  ! step start date/time (beis)
      INTEGER, SAVE :: NDATE, NTIME  ! step next date/time (beis)

      character(pathlen), allocatable :: emis_files(:)
      integer                         :: nemisfiles, iemisfile
      character(14),      allocatable :: ctotdate_emis(:)
      real(r8),           allocatable :: s1900_emis(:)
      integer                         :: iemisflg, iupdemis, iemiscyclic

      integer :: nvars3d_emis
      character(24), allocatable :: vname3d_emis(:)
      character(24), allocatable :: units3d_emis(:)

      real, allocatable :: sea_emisf(:,:,:),  sea_emisc(:,:,:),  sea_emisp(:,:,:)
      real, allocatable :: land_emisf(:,:,:), land_emisc(:,:,:), land_emisp(:,:,:)

!#ifdef Verbose2
!      character( 16 ), allocatable, save :: emis_spc( : ) ! global emis surrogates
!      integer,         allocatable, save :: emis_map( : )  ! name map to CGRID
!#endif

      CONTAINS

C-----------------------------------------------------------------------
!        FUNCTION EMIS_INIT ( JDATE, JTIME, TSTEP ) RESULT ( SUCCESS )
         subroutine emis_init ()

         USE CGRID_SPCS          ! CGRID mechanism species

!         USE BEIS_DEFN           ! biogenic emissions
!         USE LTNG_DEFN           ! NO emissions from lightning strikes
!         USE PT3D_DEFN           ! plume rise emissions
!         USE UTILIO_DEFN         ! I/O API
         USE AERO_EMIS           ! inherits GRID_CONF
         USE AERO_DATA, ONLY: MAP_PMEMIS  ! access subroutine map_pmemis

         use m3utils,   only: m3exit, xstat1, xstat2, m3warn, index1, upcase
         use misc_coms, only: io6, s1900_sim, idate1, imonth1, iyear1, itime1
         use mem_grid,  only: mwa, mza

         IMPLICIT NONE

C Includes:

         INCLUDE 'CONST.EXT'     ! constants

C Arguments:

          INTEGER :: JDATE, JTIME, TSTEP ! TSTEP is output time step (HHMMSS)
!         LOGICAL :: SUCCESS

C Parameters:

C Local Variables:

         REAL    DX1, DX2        ! CX x1- and x2-cell widths

         CHARACTER( 16 ), SAVE :: CTM_EMLAYS = 'CTM_EMLAYS'    ! env var for no. of
                                                               ! emiss. layers
         CHARACTER( 16 ) :: PNAME = 'EMIS_INIT'
         CHARACTER( 80 ) :: VARDESC   ! env variable description
         CHARACTER( 120 ) :: XMSG = ' '
         INTEGER V, LOGDEV, STATUS
         CHARACTER( 16 ) :: UNITSCK

         REAL, PARAMETER :: CMLMR  = 1.0E+06 ! ppmV/Molar Mixing Ratio
         REAL, parameter :: CONVEM = 1.0E-03 ! conversion for emissions rates to Kg/s
         real, parameter :: CNVTE = CMLMR * CONVEM * MWAIR ! conversion from moles / sec
                                                           ! to rho ppmV Vol / sec

         integer :: n, s
         real :: timefac_emis
         integer, external :: julday

!#ifdef Verbose2
!         integer s
!#endif

C-----------------------------------------------------------------------

!        SUCCESS = .TRUE.

!        LOGDEV = INIT3()
         LOGDEV = io6
         
         jdate = iyear1 * 1000 + julday(imonth1,idate1,iyear1) ! YYYYDDD
         jtime = itime1 * 100                                  ! HHMMSS
         tstep = 0

!         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
!            DX1 = DG2M * XCELL_GD ! in m.
!            DX2 = DG2M * YCELL_GD
!     &          * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( GL_NROWS/2 )))! in m.
!         ELSE
!            DX1 = XCELL_GD        ! in m.
!            DX2 = YCELL_GD        ! in m.
!         END IF

C In-line biogenic emissions?

!         IF ( .NOT. BEIS_INIT( JDATE, JTIME, TSTEP ) ) THEN
!            XMSG = 'Failure initializing biogenics emissions module'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

C Open Emissions files

!         CALL OPEMIS ( JDATE, JTIME, N_GC_EMIS, EM_TRAC, CONVEM, EMLAYS )

         ! HARDWIRE IUPDEMIS=1 for now
         iupdemis = 1

         call emis_database_read(0)
         if (nemisfiles > 1) then
            call emis_database_read(1)
         endif

         DO N = 1, N_GC_EMIS
            V = INDEX1( GC_EMIS( N ), NVARS3D_EMIS, VNAME3D_EMIS )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D_EMIS( V )
               CALL UPCASE ( UNITSCK )
               IF ( UNITSCK .NE. 'MOLES/S'   .AND.
     &              UNITSCK .NE. 'MOLE/S'    .AND.
     &              UNITSCK .NE. 'MOL/S'     .AND.
     &              UNITSCK .NE. 'MOLES/SEC' .AND.
     &              UNITSCK .NE. 'MOLE/SEC'  .AND.
     &              UNITSCK .NE. 'MOL/SEC' ) THEN
                  XMSG = 'GC units incorrect in emissions file'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            ELSE
               XMSG = 'Emissions species '
     &              // TRIM( GC_EMIS( N ) )
     &              // ' not found in emissions file'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF
         END DO

         DO N = 1, N_NR_EMIS
            V = INDEX1( NR_EMIS( N ), NVARS3D_EMIS, VNAME3D_EMIS )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D_EMIS( V )
               CALL UPCASE ( UNITSCK )
               IF ( UNITSCK .NE. 'MOLES/S'   .AND.
     &              UNITSCK .NE. 'MOLE/S'    .AND.
     &              UNITSCK .NE. 'MOL/S'     .AND.
     &              UNITSCK .NE. 'MOLES/SEC' .AND.
     &              UNITSCK .NE. 'MOLE/SEC'  .AND.
     &              UNITSCK .NE. 'MOL/SEC' ) THEN
                  XMSG = 'NR units incorrect in emissions file'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            ELSE
               XMSG = 'Emissions species '
     &              // TRIM( NR_EMIS( N ) )
     &              // ' not found in emissions file'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF
         END DO

! Time interpolation factor for updating EMIS

         timefac_emis = 0.0
         
         if (nemisfiles > 1) then
            timefac_emis = (s1900_sim             - s1900_emis(iemisfile-1))
     &                   / (s1900_emis(iemisfile) - s1900_emis(iemisfile-1))
         endif

         sea_emisc(:,:,:) = sea_emisp(:,:,:)
     &                    + (sea_emisf(:,:,:) - sea_emisp(:,:,:)) * timefac_emis

         land_emisc(:,:,:) = land_emisp(:,:,:)
     &                     + (land_emisf(:,:,:) - land_emisp(:,:,:)) * timefac_emis

C Get number of emissions layers

!         IF ( EMLAYS .LE. 1 ) THEN   ! better be inline plume rise run
!            EMLAYS = NLAYS   ! env var default
!            VARDESC = 'Number of emission layers'
!            EMLAYS = ENVINT( CTM_EMLAYS, VARDESC, EMLAYS, STATUS )
!            IF ( STATUS .NE. 0 ) THEN
!               WRITE( LOGDEV, '(5X, A)' ) VARDESC
!               IF ( STATUS .EQ. 1 ) THEN
!                  XMSG = 'Environment variable improperly formatted'
!                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
!                  SUCCESS = .FALSE.; RETURN
!               ELSE
!                  XMSG = 'Environment variable not set or empty ... Using default:'
!                  WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
!               END IF
!            END IF
!         END IF
!
!         EMLAYS = MIN ( EMLAYS, NLAYS )
!
!         WRITE( LOGDEV,1009 ) EMLAYS, NLAYS
!1009     FORMAT( / 10X, 'Number of Emissions Layers:         ', I3
!     &           / 10X, 'out of total Number of Model Layers:', I3 )

C lightning NO emissions?

!         IF ( .NOT. LTNG_INIT( JDATE, JTIME, TSTEP, EMLAYS ) ) THEN
!            XMSG = 'Failure initializing lightning emissions module'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

C Initialize aerosol emissions

         IF ( N_AE_SPC .GT. 0 ) THEN
            WRITE( LOGDEV,'( /5X, A )' ) 'Aerosol Emissions Processing in '
     &                                // 'Vertical diffusion ...'
            NEMIS_AE = N_AE_EMIS  ! from CGRID_SPCS (originally, AE_EMIS.EXT)
            ! set up PM emitted species mapping
            CALL MAP_PMEMIS ( nvars3d_emis, vname3d_emis, units3d_emis ) 
         ELSE
            NEMIS_AE = 0
         END IF

         PM_EMLYRS = EMLAYS   ! set PM_EMLYRS in AERO_EMIS module

         IF ( .NOT. AERO_EMIS_INIT( JDATE, JTIME, TSTEP, NEMIS_AE ) ) THEN
            XMSG = 'Failure initializing aerosol emissions module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         N_SPC_EMIS = N_GC_EMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

!#ifdef Verbose2
!         allocate ( emis_spc ( n_spc_emis+1 ),
!     &              emis_map ( n_spc_emis+1 ), stat = status )
!         if ( status .ne. 0 ) then
!            xmsg = 'emis_spc or emis_map memory allocation failed'
!            call m3warn ( pname, jdate, jtime, xmsg )
!            success = .false.; return
!         end if
!         emis_spc = ' '   ! array assignment
!         emis_map = 0     ! array assignment
!
!         s = 0
!         do v = 1, n_gc_emis
!            s = s + 1
!            emis_spc( s ) = gc_emis( v )
!            emis_map( s ) = gc_strt - 1 + gc_emis_map( v )
!         end do
!         do v = 1, n_ae_emis
!            s = s + 1
!            emis_spc( s ) = ae_emis( v )
!            emis_map( s ) = ae_strt - 1 + ae_emis_map( v )
!         end do
!         do v = 1, n_nr_emis
!            s = s + 1
!            emis_spc( s ) = nr_emis( v )
!            emis_map( s ) = nr_strt - 1 + nr_emis_map( v )
!         end do
!         do v = 1, n_tr_emis
!            s = s + 1
!            emis_spc( s ) = tr_emis( v )
!            emis_map( s ) = tr_strt - 1 + tr_emis_map( v )
!         end do
!         do v = 1, n_spc_emis   ! debug
!            write( logdev,'(a,2i4,3x,a)' ) ' v, emis_map(v), emis_spc(v): ',
!     &                                       v, emis_map( v ), emis_spc( v )
!         end do
!#endif

C 3D point source emissions (N_SPC_EMIS+1 serves as the null map for PTEM_MAP)
!         IF ( .NOT. PT3D_INIT( N_SPC_EMIS, EMLAYS, JDATE, JTIME, TSTEP ) ) THEN
!            XMSG = 'Failure initializing plume rise emissions module'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

! THIS IS NOW ALLOCATED ALONG WITH CGRID BECAUSE WE NEED EMISSIONS
! IN THE SCALAR TABLE

!         ALLOCATE ( VDEMIS( N_SPC_EMIS+1,EMLAYS,MY_NCOLS,MY_NROWS ),
!        allocate( vdemis( mza, mwa, N_SPC_EMIS+1 ) ) ; vdemis = 0.0
!         allocate( vdemis( mza, mwa, n_gc_emis+n_nr_emis ) ) ; vdemis = 0.0

!#ifdef Verbose2
!     &              semis( emlays,my_ncols,my_nrows ),
!#endif
!     &              STAT = STATUS )
!         IF ( STATUS .NE. 0 ) THEN
!            XMSG = 'VDEMIS memory allocation failed'
!            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
!            SUCCESS = .FALSE.; RETURN
!         END IF
!
!         IF ( PT3DEMIS ) THEN
!            EMLYRS = 1
!         ELSE
!            EMLYRS = EMLAYS
!         END IF
!
         IF ( N_SPC_EMIS .LE. 0 ) THEN
            XMSG = 'There are no emissions: VDEMIS is set to zero' ! below
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
         END IF
!
!         ALLOCATE ( BUFF( NCOLS,NROWS,EMLYRS ), STAT = STATUS )
!         IF ( STATUS .NE. 0 ) THEN
!            XMSG = 'BUFF memory allocation failed'
!            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
!            SUCCESS = .FALSE.; RETURN
!         END IF
!         BUFF = 0.0   ! array assignment

C combined gas emssions conversion factor

!        CNVTE = CMLMR * CONVEM * MWAIR / ( DX1 * DX2 )

!         STRT_GC = 1
!         FINI_GC = N_GC_EMIS
!        STRT_AE = N_GC_EMIS + 1
!        FINI_AE = N_GC_EMIS + NEMIS_AE
!        STRT_NR = N_GC_EMIS + NEMIS_AE + 1
!        FINI_NR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS
!         STRT_NR = N_GC_EMIS + 1
!         FINI_NR = N_GC_EMIS + N_NR_EMIS
!        STRT_TR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS + 1
!        FINI_TR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS + N_TR_EMIS

!        SUCCESS = .TRUE.; RETURN

         END subroutine emis_init

!C-----------------------------------------------------------------------
!         SUBROUTINE GET_EMIS ( JDATE, JTIME, TSTEP, CONVPA,
!     &                         RDX3F, RRHOJ, RJACM )
         
         SUBROUTINE GET_EMIS ( mrl )


!C NEW APPROACH:
!C Apply NML factors to all *input* emissions and inline plume rise, aero emis
!C biog, SeaSalt, Dust, and Lightning NO
!C GET_AERO_EMIS has to apply factors - don`t do it at this level.
!
         USE CGRID_SPCS          ! CGRID mechanism species
!
         USE AERO_EMIS           ! inherits GRID_CONF
!         USE BEIS_DEFN           ! biogenic emissions
!         USE BIOG_EMIS, ONLY: MSPCS
!         USE PT3D_DEFN           ! plume rise emissions
!         USE PTMAP               ! defines emissions file(s) PM mapping
!         USE LTNG_DEFN           ! lightning NO emissions
!         USE UTILIO_DEFN
!
         use m3utils, only: index1
         use consts_coms, only: r8
         use misc_coms,   only: io6, time8, timmax8, initial, s1900_sim, isubdomain
         use mem_sflux,   only: seaflux, landflux, jseaflux, jlandflux
         use mem_grid,    only: mza, volti
         use mem_ijtabs,  only: itabg_w, jtab_w, jtw_prog
         use mem_basic,   only: rho
         use mem_megan,   only: megan_emis, mgn_2_gc_map !, mgn_2_nr_map

         IMPLICIT NONE
!
!C Includes:
!
!        INCLUDE 'FILES_CTM.EXT'  ! file name parameters
         INCLUDE 'CONST.EXT'     ! constants

!C Arguments:

          integer, intent( in ) :: mrl
!         INTEGER, INTENT( IN ) :: JDATE, JTIME  ! date (YYYYDDD), time (HHMMSS)
!         INTEGER, INTENT( IN ) :: TSTEP( 3 )    ! time step vector (HHMMSS)
!         REAL,    INTENT( IN ) :: CONVPA        ! conversion factor to pressure in Pascals
!         REAL,    INTENT( IN ) :: RDX3F( NLAYS )             ! reciprocal layer thickness
!         REAL,    INTENT( IN ) :: RRHOJ( NCOLS,NROWS,NLAYS ) ! reciprocal density X Jacobian
!         REAL,    INTENT( IN ) :: RJACM( NCOLS,NROWS,NLAYS ) ! reciprocal mid-layer Jacobian
!
!C Local Variables:
!
          REAL             CNVTC         ! combined conversion factor
          REAL             CNVTI         ! intermediate combined conv. factor
!         REAL,    SAVE :: CNVTP         ! intermediate combined conv. factor
!         REAL             DELT          ! interpolation factor

         INTEGER          C, R, L, N, S, V
         integer :: j, k, ke, kw, ktop, iw, isf, ilf               ! loop induction variables

!        real, parameter :: convpa = 1.0
         REAL, PARAMETER :: CMLMR  = 1.0E+06 ! ppmV/Molar Mixing Ratio
         REAL, parameter :: CONVEM = 1.0E-03 ! conversion for emissions rates to Kg/s
         real, parameter :: CNVTE  = CMLMR * CONVEM * MWAIR ! conversion from moles / sec
                                                            ! to rho ppmV Vol / sec
         real :: timefac_emis

!         CHARACTER( 16 ) :: PNAME = 'GET_EMIS'
!         CHARACTER( 120 ) :: XMSG = ' '
!#ifdef Verbose2
!         integer, save :: logdev
!         real sumem
!#endif
!         LOGICAL, SAVE :: FIRSTIME = .TRUE.
!         INTEGER          GXOFF, GYOFF         ! global origin offset from file
!         INTEGER, SAVE :: STARTCOL, ENDCOL, STARTROW, ENDROW
!
!C-----------------------------------------------------------------------

         if (iupdemis == 1 .and. nemisfiles > 1) then
            if (s1900_sim > s1900_emis(iemisfile) .and. time8+0.001_r8 < timmax8) then
               call emis_database_read(1)
            endif

            timefac_emis = (s1900_sim             - s1900_emis(iemisfile-1))
     &                   / (s1900_emis(iemisfile) - s1900_emis(iemisfile-1))

            sea_emisc(:,:,:) = sea_emisp(:,:,:)
     &                       + (sea_emisf(:,:,:) - sea_emisp(:,:,:)) * timefac_emis

            land_emisc(:,:,:) = land_emisp(:,:,:)
     &                        + (land_emisf(:,:,:) - land_emisp(:,:,:)) * timefac_emis
         endif

!         IF ( FIRSTIME ) THEN
!            FIRSTIME = .FALSE.
!#ifdef Verbose2
!            logdev = init3()
!            logdev = io6
!#endif
!C Get domain decomp info from the emissions file
!            CALL SUBHFILE ( EMIS_1, GXOFF, GYOFF, STARTCOL, ENDCOL, STARTROW, ENDROW )
!            CNVTP = CONVPA * CNVTE
!         END IF

C read & interpolate emissions (create VDEMIS in the species class order)
         
         do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
            vdemis_gc(:,iw,:) = 0.0
         enddo

         do s = 1, n_gc_emis
            n = INDEX1( GC_EMIS(S), NVARS3D_EMIS, VNAME3D_EMIS )

            if (n /= 0) then

               CNVTI = CNVTE * GC_EMIS_FAC( S )

               do j = 1,jseaflux(1)%jend(mrl)
                  isf = jseaflux(1)%iseaflux(j)
                  kw  = seaflux(isf)%kw
                  iw  = seaflux(isf)%iw ! global index

                  ! If run is parallel, get local rank indices

                  if (isubdomain == 1) then
                     iw = itabg_w(iw)%iw_myrank
                  endif
               
                  ktop = min(kw + emlays - 1, mza-1)
                  do k = kw, ktop
                     ke = k - kw + 1
                     ! convert mol/ sec emissions to ppmV / sec
                     cnvtc = cnvti * volti(k,iw) / rho(k,iw)
                     vdemis_gc(k,iw,s) = vdemis_gc(k,iw,s) + sea_emisc(ke,isf,n) * cnvtc
                  enddo
               enddo

               do j = 1,jlandflux(1)%jend(mrl)
                  ilf = jlandflux(1)%ilandflux(j)
                  kw  = landflux(ilf)%kw
                  iw  = landflux(ilf)%iw ! global index

                  ! If run is parallel, get local rank indices

                  if (isubdomain == 1) then
                     iw = itabg_w(iw)%iw_myrank
                  endif
               
                  ktop = min(kw + emlays - 1, mza-1)
                  do k = kw, ktop
                     ke = k - kw + 1
                     ! convert mol/ sec emissions to ppmV / sec
                     cnvtc = cnvti * volti(k,iw) / rho(k,iw)
                     vdemis_gc(k,iw,s) = vdemis_gc(k,iw,s) + land_emisc(ke,ilf,n) * cnvtc
                  enddo
               enddo

            endif

!#ifdef Verbose2
!            write( logdev,* ) '  GC EMIS: ', s, '  ', gc_emis( s )
!#endif

         END DO   ! s

C beis

!         IF ( BIOGEMIS ) THEN
!
!            CALL GET_BEIS ( JDATE, JTIME, TSTEP )
!
!C tmpbeis312 is based on hourly input - need to interpolate return value
!C merge: VDEMIS_BI in mol/hr. Convert to ppm/sec
!
!#ifdef Verbose2
!            do s = 1, mspcs
!               v = beis_map( s )
!               if ( v .gt. 0 ) then
!                  write( logdev,* ) '  BEIS_MAP: ', s, v, '  ', gc_emis( v ), gc_emis_fac( v )
!               end if
!            end do
!#endif
!            L = 1
!            CNVTI = CONVPA * CNVTE * RDX3F( L )  ! per sec
!            DO R = 1, MY_NROWS
!               DO C = 1, MY_NCOLS
!                  CNVTC = CNVTI * RRHOJ( C,R,L )
!                  N = 0
!                  DO S = 1, MSPCS
!                     V = BEIS_MAP( S )
!                     IF ( V .GT. 0 ) THEN
!                        N = N + 1
!                        VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R )
!     &                                    + CNVTC * GC_EMIS_FAC( V ) * VDEMIS_BI( N,C,R )
!                      END IF
!                  END DO
!               END DO
!            END DO
!
!         END IF

C lightning NO
!         IF ( LTNG_NO ) THEN
!            CALL GET_LTNG ( JDATE, JTIME, TSTEP )
!
!            DO R = 1, MY_NROWS
!               DO C = 1, MY_NCOLS
!                  DO L = 1, EMLAYS
!                     CNVTI = CONVPA * CNVTE * RDX3F( L )  ! per sec
!                     CNVTC = CNVTI * RRHOJ( C,R,L )
!                     V = LTNG_MAP
!                     VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R )
!     &                                 + CNVTC * GC_EMIS_FAC( V ) * VDEMIS_LT( C,R,L )
!                  END DO
!               END DO
!            END DO
!
!         END IF

C plume rise

!         IF ( PT3DEMIS ) THEN
!            CALL GET_PT3D_EMIS ( JDATE, JTIME, TSTEP )
!
!C merge: VDEMIS_PT in mol/sec. Convert to ppm/sec
!C If aerosols, then the PM emissions obtained from GET_PT3D_EMIS need to be
!C accessed by GET_AERO_EMIS to be "speciated" and merged with the AE emissions.
!
!            DO R = 1, MY_NROWS
!               DO C = 1, MY_NCOLS
!                  DO L = 1, EMLAYS
!                     CNVTC = CNVTP * RDX3F( L ) * RRHOJ( C,R,L )
!                     DO S = 1, N_GC_EMIS + N_NR_EMIS + N_TR_EMIS
!                        V = SPC_PTEM_MAP( S ) ! <- pointer to non-PM (non-aerosol) species
!                        N = PTEM_MAP( S )     ! <- pointer to non-PM emissions
!                        IF ( V .LT. 1 ) CYCLE
!                        VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R )
!     &                                    + CNVTC * VDEMIS_PT( C,R,L,N )
!                     END DO
!                  END DO
!               END DO
!            END DO
!
!         END IF
!
!C aerosol emissions - all units conversions done in GET_AERO_EMIS for aerosols
!C GET_AERO_EMIS in f90 module AERO_EMIS
!
         IF ( N_AE_SPC .GT. 0 ) THEN
!           CALL GET_AERO_EMIS ( JDATE, JTIME, TSTEP, RJACM, VDEMIS )
            CALL GET_AERO_EMIS ( mrl, nvars3d_emis, vname3d_emis, sea_emisc, land_emisc )
!
!            DO R = 1, MY_NROWS
!               DO C = 1, MY_NCOLS
!                  DO L = 1, EMLAYS
!                     DO V = STRT_AE, FINI_AE
!                        S = V + 1 - STRT_AE
!                        VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R )
!     &                                    + VDEMIS_AE( S,L,C,R )
!                     END DO
!                  END DO
!               END DO
!            END DO
!
         END IF

C non-reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

         do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
            vdemis_nr(:,iw,:) = 0.0
         enddo

         do s = 1, n_nr_emis
            n = INDEX1( NR_EMIS(S), NVARS3D_EMIS, VNAME3D_EMIS )

            if (n /= 0) then

               CNVTI = CNVTE * NR_EMIS_FAC( S )

               do j = 1,jseaflux(1)%jend(mrl)
                  isf = jseaflux(1)%iseaflux(j)
                  kw  = seaflux(isf)%kw
                  iw  = seaflux(isf)%iw ! global index

                  ! If run is parallel, get local rank indices

                  if (isubdomain == 1) then
                     iw = itabg_w(iw)%iw_myrank
                  endif
               
                  ktop = min(kw + emlays - 1, mza-1)
                  do k = kw, ktop
                     ke = k - kw + 1
                     ! convert mol/ sec emissions to ppmV / sec
                     cnvtc = cnvti * volti(k,iw) / rho(k,iw)
                     vdemis_nr(k,iw,s) = vdemis_nr(k,iw,s) + sea_emisc(ke,isf,n) * cnvtc
                  enddo
               enddo

               do j = 1,jlandflux(1)%jend(mrl)
                  ilf = jlandflux(1)%ilandflux(j)
                  kw  = landflux(ilf)%kw
                  iw  = landflux(ilf)%iw ! global index

                  ! If run is parallel, get local rank indices

                  if (isubdomain == 1) then
                     iw = itabg_w (iw )%iw_myrank
                  endif
               
                  ktop = min(kw + emlays - 1, mza-1)
                  do k = kw, ktop
                     ke = k - kw + 1
                     ! convert mol/ sec emissions to ppmV / sec
                     cnvtc = cnvti * volti(k,iw) / rho(k,iw)
                     vdemis_nr(k,iw,s) = vdemis_nr(k,iw,s) + land_emisc(ke,ilf,n) * cnvtc
                  enddo
               enddo

            endif

!#ifdef Verbose2
!            write( logdev,* ) '  GC EMIS: ', s, '  ', gc_emis( s )
!#endif

         END DO   ! s


!         DO V = STRT_NR, FINI_NR
!            S = V + 1 - STRT_NR
!            BUFF = 0.0   ! array assignment
!C                                       vvvvvvv - emis surrogate names
!            IF ( .NOT. INTERPX( EMIS_1, NR_EMIS( S ), PNAME,
!     &                          STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMLYRS,
!     &                          JDATE, JTIME, BUFF ) ) THEN
!               XMSG = 'Could not read '
!     &              // TRIM( NR_EMIS( S ) ) // ' from ' // EMIS_1
!               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!            CNVTI = CNVTP * NR_EMIS_FAC( S )
!
!            DO R = 1, MY_NROWS
!               DO C = 1, MY_NCOLS
!                  DO L = 1, EMLYRS
!                     CNVTC = CNVTI * RDX3F( L ) * RRHOJ( C,R,L )
!                     VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R )
!     &                                 + CNVTC * BUFF( C,R,L )
!                  END DO
!               END DO
!            END DO
!
!#ifdef Verbose2
!            write( logdev,* ) '  NR EMIS: ', s, '  ', nr_emis( s )
!#endif
!
!         END DO



! MEGAN Emissions

         do s = 1, n_gc_emis
            n = mgn_2_gc_map(s)

            if (n /= 0) then

               do j = 1,jlandflux(1)%jend(mrl)
                  ilf = jlandflux(1)%ilandflux(j)
                  kw  = landflux(ilf)%kw
                  iw  = landflux(ilf)%iw ! global index

                  ! If run is parallel, get local rank indices
                  if (isubdomain == 1) then
                     iw = itabg_w(iw)%iw_myrank
                  endif
               
                  ! convert mol/ sec emissions to ppmV / sec
                  cnvtc = cnvte * volti(kw,iw) / rho(kw,iw)
                  vdemis_gc(kw,iw,s) = vdemis_gc(kw,iw,s) + megan_emis(ilf,n) * cnvtc

               enddo
            endif
         enddo

         ! TODO: soil nox!

!         do s = 1, n_nr_emis
!            n = mgn_2_nr_map(s)
!
!            if (n /= 0) then
!
!               do j = 1,jlandflux(1)%jend(mrl)
!                  ilf = jlandflux(1)%ilandflux(j)
!                  kw  = landflux(ilf)%kw
!                  iw  = landflux(ilf)%iw ! global index
!
!                  ! If run is parallel, get local rank indices
!                  if (isubdomain == 1) then
!                     iw = itabg_w (iw )%iw_myrank
!                  endif
!               
!                  ! convert mol/ sec emissions to ppmV / sec
!                  cnvtc = cnvte * volti(kw,iw) / rho(kw,iw)
!                  write(*,*) ilf, vdemis_nr(kw,iw,s), megan_emis(ilf,n) * cnvtc
!                  vdemis_nr(kw,iw,s) = vdemis_nr(kw,iw,s) + megan_emis(ilf,n) * cnvtc
!
!               enddo
!            endif
!         enddo



!C tracer gases (conversion to ppmv/s)
!
!         IF ( EM_TRAC ) THEN
!
!            DO V = STRT_TR, FINI_TR
!               S = V + 1 - STRT_TR
!               BUFF = 0.0   ! array assignment
!C                                               vvvvvvv - emis surrogate names
!               IF ( .NOT. INTERPX( EMIS_TRAC_1, TR_EMIS( S ), PNAME,
!     &                             STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMLYRS,
!     &                             JDATE, JTIME, BUFF ) ) THEN
!                  XMSG = 'Could not read  '
!     &                 // TRIM( TR_EMIS( S ) )
!     &                 // ' from ' // EMIS_TRAC_1
!                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!               END IF
!               CNVTI = CNVTP * TR_EMIS_FAC( S )
!
!               DO R = 1, MY_NROWS
!                  DO C = 1, MY_NCOLS
!                     DO L = 1, EMLYRS
!                        CNVTC = CNVTI * RDX3F( L ) * RRHOJ( C,R,L )
!                        VDEMIS( V,L,C,R ) = CNVTC * BUFF( C,R,L )
!                     END DO
!                  END DO
!               END DO
!
!#ifdef Verbose2
!               write( logdev,* ) '  TR EMIS: ', s, '  ', tr_emis( s )
!#endif
!
!            END DO
!
!         END IF
!
!C zero out emissions values for species not included in diffused list ...
!C accounts for emissions species names as a subset of the vert. diffused species list
!
!         VDEMIS( N_SPC_EMIS+1,:,:,: ) = 0.0
!
!#ifdef Verbose2
!!        v = index1 ( 'NO', n_gc_emis, gc_emis )
!!        semis = vdemis( v,:,:,: )   ! all layers
!!        write( logdev,* ) '  sumVDEMIS( NO ): ', sum( semis )
!         write( logdev,'(/"@@@0 ", i8, i7.6)') jdate, jtime
!         l = 1
!         do v = strt_gc, fini_gc
!            s = v + 1 - strt_gc
!            sumem = 0.0
!            do r = 1, my_nrows
!               do c = 1, my_ncols
!                  sumem = sumem + vdemis( v,l,c,r )
!               end do
!            end do
!            write( logdev,'("@@@2 ", a, g16.6)') gc_emis( s ), sumem
!         end do
!         do v = strt_nr, fini_nr
!            s = v + 1 - strt_nr
!            sumem = 0.0
!            do r = 1, my_nrows
!               do c = 1, my_ncols
!                  sumem = sumem + vdemis( v,l,c,r )
!               end do
!            end do
!            write( logdev,'("@@@3 ", a, g16.6)') nr_emis( s ), sumem
!         end do
!#endif
!
!         RETURN
!
         END SUBROUTINE GET_EMIS
      END MODULE EMIS_DEFN
