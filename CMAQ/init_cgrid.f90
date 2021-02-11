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

  use utilio_defn, only: index1, m3exit, xstat3, init3
  use misc_coms,   only: imonth1, idate1, iyear1, itime1
  use max_dims,    only: pathlen
  use cgrid_spcs
  use cgrid_defn,  only: cgrid, cgrid_names
  use string_lib
  use rxns_data
  use const_data
  use mem_grid,    only: mza, lpw
  use mem_basic,   only: rho, press
  use mem_ijtabs,  only: jtab_w, jtw_init
  use aero_data,   only: aer_str, aer_end, aer_num, srf_str, srf_end, &
                         mode_map, aer_trac, aer_wet, aer_m3, num_str, num_end
  IMPLICIT NONE

! Parameters:

  REAL, PARAMETER :: CMIN = 1.0E-30
  REAL, PARAMETER :: PTOP = 50.E2
  REAL, PARAMETER :: ONETHIRD = 1. / 3.
  REAL, PARAMETER :: MAOAVO1000 = 1.0E+03 * MWAIR / AVO

  ! External Functions:

  INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

! Local Variables

  CHARACTER( pathlen )          :: FNAME
  CHARACTER( 256 )              :: LINEIN
  character(  16 ), allocatable :: vname3d( : )

  integer,             external :: julday
  real,             allocatable :: vglevs_in(:)
  real,             allocatable :: vglays_in(:)
  real,             allocatable :: press_in (:)
  real,             allocatable :: densi_in (:)
  real,             allocatable :: vec_in   (:)
  character(  10 ), allocatable :: units_in (:)
  real,             allocatable :: inprof   (:,:)

  integer :: lo3, lch4
  INTEGER :: SPC_STRT
  INTEGER :: N_SPCS                     ! no. of species for this call
  INTEGER :: NDX                        ! loop copy of INDX
  INTEGER :: ISUR                       ! surrogate index
  INTEGER :: SPC, V, j, iw, k, ka       ! loop counters
  integer :: i, len
  real    :: psfc
  real    :: vec(mza), press_out(mza), dens_out(mza)

  character(50), parameter :: fmt9  = '( /, 5X, ''IC/BC Factors used for '', A, /)'
  character(50), parameter :: fmt13 = '( 5X, I3, 2X, A, 1PG13.5 )'
  CHARACTER(18), parameter :: ESTR1 = 'No IC for species '

  CHARACTER( 16 ) :: CONCMIN
  CHARACTER( 34 ) :: ESTR2 = ''
  CHARACTER( 90 ) :: ESTR3 = ''

  INTEGER :: INDX    ( NSPCSD ) ! Variable indices for all IC species
  REAL    :: ICBC_FAC( NSPCSD ) ! Factor to be applied to ICs

  real    :: conc(aer_num), m3(aer_num)
  real    :: m3_dry, m3_wet

  logical :: exists
  INTEGER :: LOGDEV       ! FORTRAN unit number for log file
  INTEGER :: JDATE        ! starting date,    format YYYYDDD
  INTEGER :: JTIME        ! starting time,    format HHMMSS
  integer :: nvars3d
  integer :: nlays_in
  integer :: n, iskip
  logical :: icb6

!-----------------------------------------------------------------------

  LOGDEV = INIT3()

! indx(:) = 0

  jdate = iyear1 * 1000 + julday(imonth1,idate1,iyear1) ! YYYYDDD
  jtime = itime1 * 100                                  ! HHMMSS

  ! Initialize the CGRID array (now done at allocation in cgrid_defn)

! CGRID = CMIN

  icb6 = (index(mechname, 'CB6R3') > 0)

  FNAME = '../CMAQ/' // to_lower(trim(MECHNAME)) // '/ic_profile.dat'

  inquire(file=FNAME, exist=exists)
  if (.not. exists) then
     write(*,'(/,a)') 'init_cgrid: error opening input file ' // trim(fname)
     stop       'init_cgrid: fix location of initial chemical profiles'
  endif

  open( unit=10, file=FNAME, status='old' )

  iskip = 3
  if (icb6) iskip = 11

  ! Skip the header
  DO N = 1, iskip
     READ( 10, '(A)' ) LINEIN
  END DO

  ! Read the # of layers, # of species, and sigma levels
  READ( 10, * ) nlays_in, nvars3d
  allocate( vglevs_in( nlays_in + 1 ) )
  allocate( vglays_in( nlays_in ) )
  allocate( press_in ( nlays_in ) )
  allocate( densi_in ( nlays_in ) )
  allocate( vec_in   ( nlays_in ) )
  allocate( vname3d( nvars3d ) )
  allocate( inprof( nlays_in, nvars3d ) )

  if (icb6) allocate(units_in(nvars3d))

  read( 10, * ) vglevs_in(1:nlays_in+1)

  ! Get file data

  DO N = 1, NVARS3D
     IF (ICB6) THEN
        READ( 10, * ) VNAME3D(N), UNITS_IN(N), INPROF(1:NLAYS_IN,N)
     ELSE
        READ( 10, * ) VNAME3D(N), INPROF(1:NLAYS_IN,N)
     ENDIF
  END DO

  close ( 10 )

  do n = 1, nlays_in
     vglays_in(n) = 0.5 * (vglevs_in(n) + vglevs_in(n+1))
  enddo

  ! Load CGRID

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

  LO3  = INDEX1( 'O3'  , N_GC_SPC, GC_SPC ) + GC_STRT - 1
  LCH4 = INDEX1( 'ECH4', N_GC_SPC, GC_SPC ) + GC_STRT - 1

  WRITE( LOGDEV, fmt9 ) 'initializing CMAQ gas-phase species'

  SPC_STRT = GC_STRT
  N_SPCS   = N_GC_SPC

  DO SPC = 1, N_SPCS
     V = SPC + SPC_STRT - 1

     ! skip dummy RXN counter species
     len = len_trim(GC_SPC(SPC))
     if (len > 2) then
        if (gc_spc(spc)(len-2:len) == 'RXN') cycle
     endif

     ! skip H2SO4
     if (GC_SPC(SPC) == 'SULF') cycle

     ! is there a surrogate name?
     ISUR = FINDEX ( SPC, N_GC_ICBC, GC_ICBC_MAP )
     NDX  = 0

     IF ( ISUR .NE. 0 ) THEN

        ! is it on the IC file?
        NDX = INDEX1( GC_ICBC( ISUR ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX   ! index in the IC file
           ICBC_FAC( V ) = GC_ICBC_FAC( ISUR )
           !  ELSE
           !     XMSG = ESTR1 // TRIM( GC_ICBC( ISUR ) ) // ESTR2 &
           !          // TRIM( GC_SPC( SPC ) )
           !     write( logdev, '(a)') xmsg
        END IF

     END IF

     ! Cannot find a surrogate, look for the (main) species name on the IC file
     IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

        NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX   ! index in the IC file
           ICBC_FAC( V ) = 1.0
        ELSE
           write( logdev, '(a)') ESTR1 // TRIM( GC_SPC( SPC ) ) // trim(ESTR3)
        END IF

     END IF

     IF ( INDX( V ) .GT. 0 ) &
          WRITE( LOGDEV, fmt13 ) INDX( V ), GC_SPC( SPC ), ICBC_FAC( V )

  END DO


  WRITE( LOGDEV, fmt9 ) 'initializing CMAQ aerosol species'

  SPC_STRT = AE_STRT
  N_SPCS   = N_AE_SPC

  DO SPC = 1, N_SPCS
     V = SPC + SPC_STRT - 1

     ! is there a surrogate name?
     ISUR = FINDEX ( SPC, N_AE_ICBC, AE_ICBC_MAP )
     NDX  = 0

     IF ( ISUR .NE. 0 ) THEN

        ! is it on the IC file?
        NDX = INDEX1( AE_ICBC( ISUR ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX   ! index in the IC file
           ICBC_FAC( V ) = AE_ICBC_FAC( ISUR )
           !  ELSE
           !     XMSG = ESTR1 // TRIM( AE_ICBC( ISUR ) ) // ESTR2 &
           !          // TRIM( AE_SPC( SPC ) )
           !     write( logdev, '(a)') xmsg
        END IF

     END IF

     ! Cannot find a surrogate, look for the (main) species name on the IC file
     IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

        NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX
           ICBC_FAC( V ) = 1.0
        ELSE
           write( logdev, '(a)') ESTR1 // TRIM( AE_SPC( SPC ) ) // trim(ESTR3)
        END IF
     END IF

     IF ( INDX( V ) .GT. 0 ) &
          WRITE( LOGDEV, fmt13 ) INDX( V ), AE_SPC( SPC ), ICBC_FAC( V )

  END DO


  WRITE( LOGDEV, fmt9 ) 'initializing CMAQ non-reactive gas species'

  SPC_STRT = NR_STRT
  N_SPCS   = N_NR_SPC

  DO SPC = 1, N_SPCS
     V = SPC + SPC_STRT - 1

     ! is there a surrogate name?
     ISUR = FINDEX ( SPC, N_NR_ICBC, NR_ICBC_MAP )
     NDX  = 0

     IF ( ISUR .NE. 0 ) THEN

        ! is it on the IC file?
        NDX = INDEX1( NR_ICBC( ISUR ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX   ! index in the IC file
           ICBC_FAC( V ) = NR_ICBC_FAC( ISUR )
           ! ELSE
           !    XMSG = ESTR1 // TRIM( NR_ICBC( ISUR ) ) // ESTR2 &
           !         // TRIM( NR_SPC( SPC ) )
           !    write( logdev, '(a)') xmsg
        END IF

     END IF

     ! Cannot find a surrogate, look for the (main) species name on the IC file
     IF ( ISUR .EQ. 0 .OR. NDX .EQ. 0 ) THEN

        NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )

        IF ( NDX .NE. 0 ) THEN
           INDX( V ) = NDX
           ICBC_FAC( V ) = 1.0
        ELSE
           write( logdev, '(a)') ESTR1 // TRIM( NR_SPC( SPC ) ) // trim(ESTR3)
        END IF

     END IF

     IF ( INDX( V ) .GT. 0 ) &
          WRITE( LOGDEV, fmt13 ) INDX( V ), NR_SPC( SPC ), ICBC_FAC( V )

  END DO

  !$omp parallel private(press_in,press_out,densi_in,dens_out,vec_in,vec,conc,m3)
  !$omp do private(iw,ka,psfc,spc,v,ndx,k,i,m3_dry,m3_wet)
  do j = 1, jtab_w(jtw_init)%jend(1)
     iw   = jtab_w(jtw_init)%iw(j)
     ka   = lpw(iw)

     psfc        = press(ka,iw)
     press_in(:) = vglays_in(:) * (psfc-ptop) + ptop

     press_out(1:mza-ka+1) = press(ka:mza,iw)
     dens_out (1:mza-ka+1) = rho(ka:mza,iw)

     ! interpolate model density to the input chemistry levels
     call pintrp_cc(mza-ka+1, dens_out, press_out, nlays_in, densi_in, press_in )

     densi_in = 1.0 / densi_in

     DO V = 1, NSPCSD
        NDX = INDX( V )

        ! do we have this species in the input file?
        IF ( NDX .GT. 0 ) THEN

           if (v >= ae_strt .and. v <= ae_fini) then

              if (v >= num_str .and. v <= num_end) then

                 ! number/m**3 aerosol -> ppmv
                 vec_in = inprof(:,ndx) * densi_in(:) * MAOAVO1000

              elseif (v >= srf_str .and. v <= srf_end) then

                 ! m*2/m**3 aerosol -> m**2/mol air
                 vec_in = inprof(:,ndx) * densi_in(:) * MWAIRKG

              else

                 ! micro-grams/m**3 aerosol -> ppmv
                 vec_in = inprof(:,ndx) * densi_in(:) * ae_molwti(v-aer_str+1) * MWAIRKG

              endif

              call pintrp_cc( nlays_in, vec_in, press_in, mza-ka+1, vec, press_out )

           else

              call pintrp_cc( nlays_in, inprof(:,ndx), press_in, mza-ka+1, vec, press_out )

           endif

        else

           vec = cmin

        endif

        if (V == LO3) then

           do k = ka, mza-1
              if (press(k,iw) < 200.d2) exit
           enddo
           cgrid(ka:k,iw,v) = max( ICBC_FAC( V ) * vec(1:k-ka+1), cmin)

        elseif (V == LCH4) then

!          cgrid(ka:mza,iw,v) = ATM_CH4 + cgrid(1:mza,iw,v)
           cgrid(ka:mza,iw,v) = ATM_CH4

        else

           cgrid(ka:mza,iw,v) = max( ICBC_FAC( V ) * vec(1:mza-ka+1), cmin)

        endif

        cgrid(1:ka-1,iw,v) = cgrid(ka,iw,v)

     enddo

     ! Convert aerosol 2nd moment to wet

     do k = ka, mza
        conc = cgrid(k,iw,aer_str:aer_end)
        m3   = conc * aer_m3 * ae_molwt(1:aer_num) * INV_MWAIRKG * real(rho(k,iw))

        do i = 1, 3
           v = srf_str + i - 1

           m3_dry = sum(m3, mask=(mode_map==i .and. .not. aer_trac .and. .not. aer_wet))
           m3_wet = sum(m3, mask=(mode_map==i .and. .not. aer_trac))
           cgrid(k,iw,v) = cgrid(k,iw,v) * ( (m3_wet / m3_dry)**2 ) ** onethird
        enddo
     enddo

     do v = srf_str, srf_end
        do k = 1, ka-1
           cgrid(k,iw,v) = cgrid(ka,iw,v)
        enddo
     enddo

  enddo
  !$omp end do
  !$omp end parallel

end subroutine init_cgrid
