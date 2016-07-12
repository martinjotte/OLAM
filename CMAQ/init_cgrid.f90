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

subroutine init_cgrid ()

  !-----------------------------------------------------------------------
  ! Function:
  !   Initialize simulation time period and time stepping constants for
  !   core model driver
  !   Environment variable can reference a previous CONC file to use as
  !   initial data.
  !   Write initial conc data as step "0" on output conc file
  
  !   revised  6/12/97 by Jerry Gipson: Get ICs by species name, by surrogate
  !            name, or zero
  !   Jeff - Aug 97 - fixed problems, cleaned up
  !   Jeff - Dec 97 - add CMIN
  !   Jeff - Dec 97 - put in aerosol sulfate inititalization
  !   Jeff - Feb 98 - close init cond files after reading
  !   2 October, 1998 by Al Bourgeois at LM: parallel implementation
  !   revised 10/7/99 by Shawn Roselle: added surface area species to
  !            aerosol species types
  !   Jeff - Dec 00 - check if append, split out opconc and load_cgrid
  !                 - move CGRID_MAP into f90 module
  !   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; assumed shape arrays
  !   17 Mar 03 D.Wong: move barrier to avoid race conditions
  !   28 Aug 03 J.Young: following Zion Wang at CERT, only pe 0 closes
  !   30 May 05 J.Young: add call to load RHOJ into CGRID
  !   21 Jun 10 J.Young: convert for Namelist redesign
  !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
  !                        removed deprecated TRIMLEN
  !    2 Sep 11 J.Young: change ICBC_FAC policy to always assigning factor,
  !                      if specified, not just if a surrogate is also specified 
  !   11 Sep 15 B.Murphy: add condition for no surrogate name
  !-----------------------------------------------------------------------

  use utilio_defn
  use misc_coms,  only: imonth1, idate1, iyear1, itime1
  use max_dims,   only: pathlen      
  use cgrid_spcs
  use cgrid_defn, only: cgrid
  use string_lib
  use rxns_data

  IMPLICIT NONE

! Parameters:

  REAL, PARAMETER :: CMIN = 1.0E-30

! Local Variables

  CHARACTER( pathlen )          :: FNAME
  CHARACTER( 256 )              :: LINEIN
  character(  16 ), allocatable :: vname3d( : )

  integer,             external :: julday
  real,             allocatable :: vglevs_in(:)
  real,             allocatable :: vglays_in(:)
  real,             allocatable :: vghgts_in(:)
  real,             allocatable :: inprof( :,: )
  
  logical :: exists
  INTEGER :: LOGDEV       ! FORTRAN unit number for log file
  INTEGER :: JDATE        ! starting date,    format YYYYDDD
  INTEGER :: JTIME        ! starting time,    format HHMMSS
  integer :: nvars3d
  integer :: nlays_in
  integer :: n

!-----------------------------------------------------------------------

  LOGDEV = INIT3()

! indx(:) = 0

  jdate = iyear1 * 1000 + julday(imonth1,idate1,iyear1) ! YYYYDDD
  jtime = itime1 * 100                                  ! HHMMSS

  ! Initialize the CGRID array (now done at allocation in cgrid_defn)

! CGRID = CMIN
  
  FNAME = '../CMAQ/' // to_lower(trim(MECHNAME)) // '/ic_profile.dat'
  
  inquire(file=FNAME, exist=exists)
  if (.not. exists) then
     write(*,'(/,a)') 'init_cgrid: error opening input file ' // trim(fname)
     stop       'init_cgrid: fix location of initial chemical profiles'
  endif

  open( unit=10, file=FNAME, status='old' )
      
  ! Skip the 3-line header
  DO N = 1, 3
     READ( 10, '(A)' ) LINEIN
  END DO

  ! Read the # of layers, # of species, and sigma levels
  READ( 10, '(A)' ) LINEIN
  read( linein, * ) nlays_in, nvars3d
  allocate( vglevs_in( nlays_in + 1 ) )
  allocate( vglays_in( nlays_in ) )
  allocate( vghgts_in( nlays_in ) )
  allocate( vname3d( nvars3d ) )
  allocate( inprof( nlays_in, nvars3d ) )
  read( linein, * ) nlays_in, nvars3d, vglevs_in(1:nlays_in+1)

  ! Consume a date and time line
  READ( 10, * ) 

  ! Get file data
  DO N = 1, NVARS3d
     READ( 10, * ) vname3d(n), inprof(1:nlays_in,n)
  END DO

  close ( 10 )
  
  do n = 1, nlays_in
     vglays_in(n) = 0.5 * (vglevs_in(n) + vglevs_in(n+1))
  enddo

  ! Load CGRID
  
  if ( n_gc_spc .gt. 0 ) then
     call load_cgrid ( 'GC' )
  end if

! Get aerosols IC's

  if ( n_ae_spc .gt. 0 ) then
     call load_cgrid ( 'AE' )
  end if

  ! Get non-reactives IC's

  if ( n_nr_spc .gt. 0 ) then
     call load_cgrid ( 'NR' )
  end if

!!  ! Get tracer IC's
!!  if ( n_tr_spc .gt. 0 ) then
!!     call load_cgrid ( 'TR' )
!!  end if


contains


  subroutine load_cgrid ( spc_cat )

    use utilio_defn, only: index1, m3exit, xstat3
    use mem_grid,    only: mza, mwa, lpw, zm, zt
    use mem_basic,   only: rho, press
    use mem_ijtabs,  only: jtab_w, jtw_init

    implicit none

    CHARACTER( 2 ),  intent(in) :: SPC_CAT

    INCLUDE 'CONST.EXT'       ! constants

    ! minimum aerosol sulfate concentration [ ug/m**3 ]
    REAL, PARAMETER :: AEROCONCMIN = 0.001

    ! The following two factors assume that sulfate density is 1.8e3 [ kg/m**3 ]
    ! and that the geometric mean diameter and geometric standard deviations
    ! for the Aitken mode are 0.01e-6 [ m ] and 1.7 respectively
    ! and are 0.07e-6 and 2.0 respectively for the accumulation mode.
    
    ! factor to calculate aerosol number concentration from aerosol sulfate mass
    ! concentration in the Aitken mode [ ug ].
    REAL, PARAMETER :: NUMFACT_I = 2.988524E11
  
    ! factor to calculate aerosol number concentration from aerosol sulfate mass
    ! concentration in the Accumulation mode [ ug ].
    REAL, PARAMETER :: NUMFACT_J = 3.560191E08

    ! fraction of sulfuric acid vapor taken as aerosol for first time step
    REAL, PARAMETER :: SO4VAPTOAER = 0.999

    ! initial fraction of total aerosol sulfate in the Aitken mode
    REAL, PARAMETER :: IFRACATKN = 0.04

    ! External Functions:

    INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

    ! Local Variables

    REAL    :: MWH2SO4                           ! H2SO4 molec. wt.
    REAL    :: H2SO4CONV                         ! ppm -> ug/m**3
    INTEGER :: LSULF, LO3                        ! Gas chem CGRID index
    INTEGER :: ISO4AJ, ISO4AI, INUMATKN, INUMACC ! CGRID aerosol indices

    INTEGER :: SPC_STRT
    INTEGER :: N_SPCS                     ! no. of species for this call
    INTEGER :: NDX                        ! loop copy of INDX
    INTEGER :: ISUR                       ! surrogate index
    INTEGER :: SPC, V, j, iw, k, ka, kr   ! loop counters
    real    :: zbot, ztop
    real    :: vec (mza)

    character(50), parameter :: fmt9  = '( /, 5X, ''IC/BC Factors used for '', A, /)'
    character(50), parameter :: fmt13 = '( 5X, I3, 2X, A, 1PG13.5 )'
    
    CHARACTER( 16 ) :: PNAME = 'LOAD_CGRID'
    CHARACTER( 16 ) :: VNAME
    CHARACTER( 16 ) :: CONCMIN
    CHARACTER( 96 ) :: XMSG  = ''
    CHARACTER( 24 ) :: ESTR1 = 'No IC for species '
    CHARACTER( 34 ) :: ESTR2 = ''
    CHARACTER( 34 ) :: ESTR3 = ''

    INTEGER :: INDX    ( MAX( N_GC_SPC, N_AE_SPC, N_NR_SPC ) ) ! Variable indices for all IC species
    REAL    :: ICBC_FAC( MAX( N_GC_SPC, N_AE_SPC, N_NR_SPC ) ) ! Factor to be applied to ICs
    
    ! The original policy for using surrogate names is first, check if the Namelist
    ! species is on the I! file; if so ignore any surrogate. If the Namelist species
    ! is not on the IC file, then check if the surrogate name is; if so also use the
    ! scale factor (default = 1.0).
    ! Note: parsing in CGRID_SPCS follows this policy for all the Namelist surrogate
    ! types (EMIS, DEPV, ICBC, and SCAV).
    ! => Change this for ICBC:
    ! First check if there's a surrogate name in the Namelist and use it (and the
    ! corresponding scale factor) if it exists. If it's not on the IC file, which it
    ! wouldn`t be if it were blank, e.g., then look for the Namelist species name. If
    ! that name is found on the IC file, then the default scale factor is applied
    ! (default = 1.0). To use a scale factor other that 1.0, there must be a name in
    ! the surrogate slot; it could be the same as the Namelist main species name.

    indx    (:) = 0
    icbc_fac(:) = 0.0

    WRITE( CONCMIN,'(1PE9.2)' ) CMIN
    ESTR2 = ' in ' // TRIM( FNAME ) // '; Look for '
    ESTR3 = ' in ' // TRIM( FNAME ) // '; set to ' // TRIM( CONCMIN )

    IF ( SPC_CAT .EQ. 'GC' ) THEN

       WRITE( LOGDEV, fmt9 ) 'transported gas-phase species'

       SPC_STRT = GC_STRT
       N_SPCS   = N_GC_SPC

       LO3 = INDEX1( 'O3', N_GC_SPC, GC_SPC )

       DO SPC = 1, N_SPCS

          ! is there a surrogate name?
          ISUR = FINDEX ( SPC, N_GC_ICBC, GC_ICBC_MAP )
          NDX  = 0

          IF ( ISUR .NE. 0 ) THEN

             ! is it on the IC file?
             NDX = INDEX1( GC_ICBC( ISUR ), NVARS3D, VNAME3D )

             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX   ! index in the IC file
                ICBC_FAC( SPC ) = GC_ICBC_FAC( ISUR )
             ELSE
                XMSG = ESTR1 // TRIM( GC_ICBC( ISUR ) ) // ESTR2 &
                     // TRIM( GC_SPC( SPC ) )
                write( logdev, '(a)') xmsg
             END IF

          END IF

          ! Cannot find a surrogate, look for the (main) species name on the IC file
          IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

             NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )

             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX   ! index in the IC file
                ICBC_FAC( SPC ) = 1.0
             ELSE
                XMSG = ESTR1 // TRIM( GC_SPC( SPC ) ) // ESTR3
                write( logdev, '(a)') xmsg
             END IF

          END IF

          IF ( INDX( SPC ) .GT. 0 ) &
               WRITE( LOGDEV, fmt13 ) INDX( SPC ), GC_SPC( SPC ), ICBC_FAC( SPC )

       END DO

     
    ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN
       
       WRITE( LOGDEV, fmt9 ) 'transported aerosol species'
       
       SPC_STRT = AE_STRT
       N_SPCS   = N_AE_SPC

       DO SPC = 1, N_SPCS

          ! is there a surrogate name?
          ISUR = FINDEX ( SPC, N_AE_ICBC, AE_ICBC_MAP )
          NDX  = 0

          IF ( ISUR .NE. 0 ) THEN

             ! is it on the IC file?
             NDX = INDEX1( AE_ICBC( ISUR ), NVARS3D, VNAME3D )

             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX   ! index in the IC file
                ICBC_FAC( SPC ) = AE_ICBC_FAC( ISUR )
             ELSE
                XMSG = ESTR1 // TRIM( AE_ICBC( ISUR ) ) // ESTR2 &
                     // TRIM( AE_SPC( SPC ) )
                write( logdev, '(a)') xmsg
             END IF

          END IF

          ! Cannot find a surrogate, look for the (main) species name on the IC file
          IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

             NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )

             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX
                ICBC_FAC( SPC ) = 1.0
             ELSE
                XMSG = ESTR1 // TRIM( AE_SPC( SPC ) ) // ESTR3
                write( logdev, '(a)') xmsg
             END IF
          END IF

          IF ( INDX( SPC ) .GT. 0 ) &
               WRITE( LOGDEV, fmt13 ) INDX( SPC ), AE_SPC( SPC ), ICBC_FAC( SPC )
        
       END DO
       
    ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN

       WRITE( LOGDEV, fmt9 ) 'transported non-reactive gas species'

       SPC_STRT = NR_STRT
       N_SPCS   = N_NR_SPC

       DO SPC = 1, N_SPCS

          ! is there a surrogate name?
          ISUR = FINDEX ( SPC, N_NR_ICBC, NR_ICBC_MAP )
          NDX  = 0
          
          IF ( ISUR .NE. 0 ) THEN

             ! is it on the IC file?
             NDX = INDEX1( NR_ICBC( ISUR ), NVARS3D, VNAME3D )

             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX   ! index in the IC file
                ICBC_FAC( SPC ) = NR_ICBC_FAC( ISUR )
             ELSE
                XMSG = ESTR1 // TRIM( NR_ICBC( ISUR ) ) // ESTR2 &
                     // TRIM( NR_SPC( SPC ) )
                write( logdev, '(a)') xmsg
             END IF

          END IF

          ! Cannot find a surrogate, look for the (main) species name on the IC file
          IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

             NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )
           
             IF ( NDX .NE. 0 ) THEN
                INDX( SPC ) = NDX
                ICBC_FAC( SPC ) = 1.0
             ELSE
                  XMSG = ESTR1 // TRIM( NR_SPC( SPC ) ) // ESTR3
                write( logdev, '(a)') xmsg
             END IF

          END IF

          IF ( INDX( SPC ) .GT. 0 ) &
               WRITE( LOGDEV, fmt13 ) INDX( SPC ), NR_SPC( SPC ), ICBC_FAC( SPC )
        
       END DO

!!    ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN
!!
!!       WRITE( LOGDEV, fmt9 ) 'transported inert tracer gas species'
!!
!!       SPC_STRT = TR_STRT
!!       N_SPCS   = N_TR_SPC
!!
!!       DO SPC = 1, N_SPCS
!!
!!          ! is there a surrogate name?
!!          ISUR = FINDEX ( SPC, N_TR_ICBC, TR_ICBC_MAP )
!!
!!          IF ( ISUR .NE. 0 ) THEN
!!
!!             ! is it on the IC file?
!!             NDX = INDEX1( TR_ICBC( ISUR ), NVARS3D, VNAME3D )
!!
!!             IF ( NDX .NE. 0 ) THEN
!!                INDX( SPC ) = NDX   ! index in the IC file
!!                ICBC_FAC( SPC ) = TR_ICBC_FAC( ISUR )
!!             ELSE
!!                XMSG = ESTR1 // TRIM( TR_SPC( SPC ) ) // ESTR2
!!                write( logdev, '(a)') xmsg
!!             END IF
!!
!!          ELSE
!!
!!             ! is the (main) species name on the IC file?
!!             NDX = INDEX1( TR_SPC( SPC ), NVARS3D, VNAME3D )
!!             IF ( NDX .NE. 0 ) THEN
!!                INDX( SPC ) = NDX
!!                ICBC_FAC( SPC ) = 1.0
!!             ELSE
!!                XMSG = ESTR1 // TRIM( TR_SPC( SPC ) ) // ESTR2
!!                write( logdev, '(a)') xmsg
!!             END IF
!!        
!!          END IF
!!        
!!          IF ( INDX( SPC ) .GT. 0 ) &
!!               WRITE( LOGDEV, fmt13 ) INDX( SPC ), TR_SPC( SPC ), ICBC_FAC( SPC )
!!
!!       END DO

    ELSE

       XMSG = 'Species categories incorrect for CGRID '
       CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, 0 )

    END IF

    DO SPC = 1, N_SPCS

       V   = SPC_STRT - 1 + SPC
       NDX = INDX( SPC )

       ! Skip ozone, since we use climatological values or reanalysis,
       ! but set to a constant value of 0.025 ppmV in troposphere
       if ( SPC_CAT .EQ. 'GC' ) THEN
          if (V == LO3) then
             
             do j = 1, jtab_w(jtw_init)%jend(1)
                iw   = jtab_w(jtw_init)%iw(j)

                ka = mza
                do k = lpw(iw), mza
                   if (press(k,iw) < 200.0e2) then
                      ka = k
                      exit
                   endif
                enddo

                cgrid(1:ka,iw,v) = 0.025
             enddo

             NDX = 0
          endif
       endif

       IF ( NDX .GT. 0 ) THEN

          do j = 1, jtab_w(jtw_init)%jend(1)
             iw   = jtab_w(jtw_init)%iw(j)
             ka   = lpw(iw)
             ztop = zm(mza)
             zbot = zm(ka-1)

             vghgts_in(:) = ztop - (ztop - zbot) * vglays_in(:) 

             call hintrp_cc ( nlays_in, inprof(1,ndx), vghgts_in, mza, vec, zt )

             cgrid(:,iw,v) = ICBC_FAC( SPC ) * vec(:)

          enddo

       endif
    enddo

    IF ( SPC_CAT .EQ. 'AE' ) THEN

       ! are ASO4J IC`s available on the file?

       VNAME = 'ASO4J'
       NDX = INDEX1( VNAME, NVARS3D, VNAME3D )
     
       IF ( NDX .EQ. 0 ) THEN  ! ASO4J not on file

          ! Set pointers for gas (vapor) phase sulfur species

          NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
          IF ( NDX .NE. 0 ) THEN
             ISO4AJ = AE_STRT - 1 + NDX
          ELSE
             XMSG = 'Could not find ' // VNAME // 'in aerosol table'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
          END IF

          VNAME = 'SULF'
          NDX = INDEX1( VNAME, N_GC_G2AE, GC_G2AE )
          IF ( NDX .NE. 0 ) THEN
             LSULF   = GC_STRT - 1 + GC_G2AE_MAP( NDX )
             MWH2SO4 = GC_MOLWT( GC_G2AE_MAP( NDX ) )
          ELSE
             XMSG = 'Could not find ' // VNAME // 'in gas chem aerosol table'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
          END IF

          VNAME = 'ASO4I'
          NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
          IF ( NDX .NE. 0 ) THEN
             ISO4AI = AE_STRT - 1 + NDX
          ELSE
             XMSG = 'Could not find ' // VNAME // 'in aerosol table'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
          END IF

          VNAME = 'NUMATKN'
          NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
          IF ( NDX .NE. 0 ) THEN
             INUMATKN = AE_STRT - 1 + NDX
          ELSE
             XMSG = 'Could not find ' // VNAME // 'in aerosol table'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
          END IF

          VNAME = 'NUMACC'
          NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
          IF ( NDX .NE. 0 ) THEN
             INUMACC = AE_STRT - 1 + NDX
          ELSE
             XMSG = 'Could not find ' // VNAME // 'in aerosol table'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
          END IF

          ! Partition the aerosol sulfate arrays with a fraction of the initial SO4 

          H2SO4CONV = 1.0E3 * MWH2SO4 / MWAIR * SO4VAPTOAER

          do j = 1, jtab_w(jtw_init)%jend(1) ; iw = jtab_w(jtw_init)%iw(j)
             do k = 1, mza
                kr = max(k, lpw(iw))

                ! total accumulation mode sulfate:

                CGRID( k,iw,ISO4AJ )   = MAX ( AEROCONCMIN,     &
                                         ( 1.0 - IFRACATKN )    &
                                       * H2SO4CONV              &
                                       * rho ( kr, iw )         &
                                       * CGRID( K,IW,LSULF ) )

                ! Accumulation mode number:
    
                CGRID( k,iw,INUMACC )  = NUMFACT_J              &
                                       * CGRID( K,IW,ISO4AJ )

                ! Aitken mode sulfate:
    
                CGRID( k,iw,ISO4AI )   = MAX ( AEROCONCMIN,     &
                                         IFRACATKN              &
                                       * H2SO4CONV              &
                                       * rho ( kr, iw )         &
                                       * CGRID( K,IW,LSULF ) )
    
                ! Aitken mode number:

                CGRID( k,iw,INUMATKN ) = NUMFACT_I              &
                                       * CGRID( K,IW,ISO4AI )
    
                ! correct sulfate vapor concentration for part removed:
    
                CGRID( k,iw,LSULF )    = ( 1.0 - SO4VAPTOAER )  &
                                       * CGRID( K,IW,LSULF )
    
             END DO
          END DO

          write(*,*) 'No IC''s found for aerosol sulfate. '
          write(*,*) 'Gas Chem sulfate used for partitioning.'
            
       END IF  ! NDX .EQ. 0

    END IF  !  SPC_CAT .EQ. 'AE'

  end subroutine load_cgrid

end subroutine init_cgrid
