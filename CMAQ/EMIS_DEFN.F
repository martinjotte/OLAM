
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
! $Header: /project/yoj/arc/CCTM/src/emis/emis/EMIS_DEFN.F,v 1.10 2011/10/21 16:10:46 yoj Exp $

! what(1) key, module and SID; SCCS file; date and time of last delta:
! %W% %P% %G% %U%

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE EMIS_DEFN

!-----------------------------------------------------------------------
! Function: emissions interface to the chemistry-transport model

! Revision History:
!     28 Jul 2006 J.Young: initial implementation
!     18 Aug 2007 J.Young: move beis part to separate module; add plume rise
!     23 Sep 2009 B.Hutzell: modified algorithm that loads gas emissions from point
!                 sources into VDEMIS array to enable multi-use of an emission species
!     26 Jan 2010 J.Young: fix bug overwriting point source layer 1 NH3; inline rdemis
!     07 Jan 2011 B.Hutzell: updated for namelist definition of model species
!     16 Feb 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
!                            removed deprecated TRIMLEN
!      6 Apr 2011 J.Young, R.Pinder: add lightning NO emissions capability
!     11 May 2011 D.Wong: incorporated twoway model implementation

!-----------------------------------------------------------------------

      use max_dims,    only: pathlen
      use consts_coms, only: r8
      use cgrid_defn,  only: vdemis_gc, vdemis_nr, vdemis_ae

      IMPLICIT NONE

      INTEGER, SAVE :: N_SPC_EMIS  ! global emissions species

      PUBLIC N_SPC_EMIS, EMIS_INIT, GET_EMIS
      public emis_files, nemisfiles, ctotdate_emis, s1900_emis
      public iemisflg, iupdemis, iemiscyclic, iemisfile
      public nvars3d_emis, vname3d_emis, units3d_emis

      PRIVATE

      INTEGER, SAVE :: NEMIS_AE   ! no. of aero emis species

      REAL,    SAVE :: DT         ! TSTEP (output) in sec

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

      integer, allocatable :: edgar_kmax(:)

      integer, allocatable, save :: gc_emis_index(:)
      integer, allocatable, save :: nr_emis_index(:)
      logical, allocatable, save :: is_poc       (:)
      integer,              save :: pncom_index

      CONTAINS

!-----------------------------------------------------------------------

      subroutine emis_init ()

         USE CGRID_SPCS          ! CGRID mechanism species
         USE LTNG_DEFN           ! NO emissions from lightning strikes
         USE UTILIO_DEFN         ! I/O API
         USE AERO_EMIS           ! inherits GRID_CONF
         USE AERO_DATA,     ONLY: MAP_PMEMIS  ! access subroutine map_pmemis
         USE CONST_DATA
         use misc_coms,     only: idate1, imonth1, iyear1, itime1
         use mem_grid,      only: mwa, mza, lpw, lsw
         use edgar_42_emis, only: edgar_42_init, vinterp, nsec
         use soil_nox,      only: soil_nox_init
         use mem_ijtabs,    only: jtab_w, jtw_prog
         use oname_coms,    only: nl

         IMPLICIT NONE

! Parameters:

         REAL, PARAMETER :: CMLMR  = 1.0E+06 ! ppmV/Molar Mixing Ratio
         REAL, parameter :: CONVEM = 1.0E-03 ! conversion for emissions rates to Kg/s
         real, parameter :: CNVTE = CMLMR * CONVEM * MWAIR ! conversion from moles / sec
                                                           ! to rho ppmV Vol / sec
! Local Variables:

         INTEGER :: JDATE, JTIME, TSTEP
         CHARACTER( 16 ) :: PNAME = 'EMIS_INIT'
         CHARACTER( 120 ) :: XMSG = ' '
         INTEGER :: V, LOGDEV
         CHARACTER( 16 ) :: UNITSCK

         integer :: n, jw, iw, ka
         integer, external :: julday

         if (nl%do_emis == 0) return

         LOGDEV = INIT3()

         jdate = iyear1 * 1000 + julday(imonth1,idate1,iyear1) ! YYYYDDD
         jtime = itime1 * 100                                  ! HHMMSS
         tstep = 0

         ! Hardwire iupdemis=1 for now
         iupdemis = 1

         call edgar_42_init(nvars3d_emis, vname3d_emis, units3d_emis)

         allocate(edgar_kmax(mwa)) ; edgar_kmax = mza
         do jw = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(jw)
            ka = lpw(iw) + lsw(iw) - 1
            edgar_kmax(iw) = maxval(vinterp(1:nsec,ka)%nlevs) + ka - 1
         enddo

         allocate(is_poc       (n_gc_emis)) ; is_poc        = .false.
         allocate(gc_emis_index(n_gc_emis)) ; gc_emis_index = 0
         allocate(nr_emis_index(n_nr_emis)) ; nr_emis_index = 0

         pncom_index = INDEX1( 'PNCOM', NVARS3D_EMIS, VNAME3D_EMIS )

         DO N = 1, N_GC_EMIS
            V = INDEX1( GC_EMIS( N ), NVARS3D_EMIS, VNAME3D_EMIS )
            gc_emis_index(n) = v
            if (gc_emis(n) == 'POC') is_poc(n) = .true.
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D_EMIS( V )
               CALL UPCASE ( UNITSCK )
               IF (GC_EMIS(n) == 'POC') then
                  if ( unitsck .ne. 'G/S'   .and.
     &                 unitsck .ne. 'G/SEC' ) then
                     XMSG = 'GC units incorrect in emissions file'
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                  endif
               ELSEIF ( UNITSCK .NE. 'MOLES/S'   .AND.
     &                  UNITSCK .NE. 'MOLE/S'    .AND.
     &                  UNITSCK .NE. 'MOL/S'     .AND.
     &                  UNITSCK .NE. 'MOLES/SEC' .AND.
     &                  UNITSCK .NE. 'MOLE/SEC'  .AND.
     &                  UNITSCK .NE. 'MOL/SEC' ) THEN
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
            nr_emis_index(n) = v
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

! lightning NO emissions?

         IF ( .NOT. LTNG_INIT() ) THEN
            XMSG = 'Failure initializing lightning emissions module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

! Initialize aerosol emissions

         IF ( N_AE_SPC .GT. 0 ) THEN
            WRITE( LOGDEV,'( /5X, A )' ) 'Aerosol Emissions Processing in '
     &                                // 'Vertical diffusion ...'
            NEMIS_AE = N_AE_EMIS  ! from CGRID_SPCS (originally, AE_EMIS.EXT)
            ! set up PM emitted species mapping
            CALL MAP_PMEMIS ( nvars3d_emis, vname3d_emis, units3d_emis )
         ELSE
            NEMIS_AE = 0
         END IF

         IF ( .NOT. AERO_EMIS_INIT( JDATE, JTIME, TSTEP, NEMIS_AE,
     &                              nvars3d_emis, vname3d_emis ) ) THEN
            XMSG = 'Failure initializing aerosol emissions module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         N_SPC_EMIS = N_GC_EMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

! 3D point source emissions (N_SPC_EMIS+1 serves as the null map for PTEM_MAP)

!        IF ( .NOT. PT3D_INIT( N_SPC_EMIS, EMLAYS, JDATE, JTIME, TSTEP ) ) THEN
!           XMSG = 'Failure initializing plume rise emissions module'
!           CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!        END IF

! Check the Emissions Species Available on the Input Files Against the
! Surrogates Identified by the User via the Namelists and Stop the model
! or Print Warnings if Mistakes Are Made.
!        CALL EMIS_SPC_CHECK( JDATE, JTIME )

C Warn the User if there are No Emissions Species On Any File
         IF ( N_SPC_EMIS .LE. 0 ) THEN
            XMSG = 'There are no emissions: VDEMIS is set to zero' ! below
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
         END IF

         call soil_nox_init()

      end subroutine emis_init

!-----------------------------------------------------------------------

      subroutine get_emis ( )

         USE CGRID_SPCS          ! CGRID mechanism species
         USE AERO_EMIS           ! inherits GRID_CONF
!        USE PT3D_DEFN           ! plume rise emissions
!        USE PTMAP               ! defines emissions file(s) PM mapping
         USE LTNG_DEFN           ! lightning NO emissions
         USE UTILIO_DEFN
         USE CONST_DATA
         use consts_coms,  only: r8
         use leaf_coms,    only: isfcl
         use mem_grid,     only: mza, volti, lpw
         use mem_ijtabs,   only: jtab_w, jtw_prog, itab_w
         use mem_land,     only: omland
         use mem_sfcg,     only: itab_wsfc
         use mem_basic,    only: rho
         use mem_megan,    only: megan_driver, megan_emis, mgn_2_gc_map
         use edgar_42_emis,only: process_edgar_emis, edgar42_emis
         use soil_nox,     only: get_soil_nox, sn_map, soilnox
         use oname_coms,   only: nl

         IMPLICIT NONE

! Local Variables:

         REAL    :: CNVTC(mza)  ! combined conversion factor
         INTEGER :: N, S
         integer :: j, k, kw, iw, jsfc, iwsfc, jasfc, iland
         real    :: arc, arcoarsfc, arcoariw, arcoarkw   ! Use only one of these

! Parameters

         REAL, PARAMETER :: CMLMR  = 1.0E+06 ! ppmV/Molar Mixing Ratio
         REAL, parameter :: CONVEM = 1.0E-03 ! conversion for emissions rates to Kg/s
         real, parameter :: CNVTE  = CMLMR * CONVEM * MWAIR ! conversion from moles / sec
                                                            ! to rho ppmV Vol / sec
         if (nl%do_emis == 0) return

! Read and/or compute emissions

         call process_edgar_emis()

         if (ltng_no) then
            call get_ltng( )
         endif

         if (isfcl == 1) then
            call get_soil_nox()
            call megan_driver()
         endif

! Apply gaseous emissions (gas chemistry and non-reacting species)

         !$omp parallel private(cnvtc)
         !$omp do private(iw,k,s,n,jsfc,iwsfc,jasfc,kw,iland)
         do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

            ! Zero-out emissions arrays

            vdemis_gc( :,iw,: ) = 0.0
            vdemis_nr( :,iw,: ) = 0.0
            vdemis_ae( :,iw,: ) = 0.0

            do k = lpw(iw), mza
               cnvtc(k) = cnvte * volti(k,iw) / real(rho(k,iw))
            enddo

            ! Apply EDGAR gas chemistry emissions

            do s = 1, n_gc_emis
               n = gc_emis_index(s)

               if (n /= 0) then

                  if (is_poc(s)) then

                     if (pncom_index /= 0) then
                        do k = lpw(iw), edgar_kmax(iw)
                           ! Convert g / sec emissions to ppmV / sec,
                           ! and include non-carbon PNCOM emissions with POCs
                           vdemis_gc(k,iw,s) = (edgar42_emis(k,iw,n) + edgar42_emis(k,iw,pncom_index))
     &                                       * cnvtc(k) * GC_EMIS_FAC( S )
                        enddo
                     else
                        do k = lpw(iw), edgar_kmax(iw)
                           ! Convert g / sec emissions to ppmV / sec
                           vdemis_gc(k,iw,s) = edgar42_emis(k,iw,n) * cnvtc(k) * GC_EMIS_FAC( S )
                        enddo
                     endif

                  else

                     do k = lpw(iw), edgar_kmax(iw)
                        ! convert mol / sec emissions to ppmV / sec
                        vdemis_gc(k,iw,s) = edgar42_emis(k,iw,n) * cnvtc(k)
                     enddo
                  endif

               endif
            enddo

            ! Apply EDGAR non-reactive gas emissions (conversion to ppmv/s)

            do s = 1, n_nr_emis

               n = INDEX1( NR_EMIS(S), NVARS3D_EMIS, VNAME3D_EMIS )
               if (n /= 0) then

                  do k = lpw(iw), edgar_kmax(iw)
                     ! Convert mol / sec emissions to ppmV / sec
                     vdemis_nr(k,iw,s) = edgar42_emis(k,iw,n) * cnvtc(k)
                  enddo

               endif
            enddo

            ! Lightning NOx emissions

            if ( ltng_no ) then

               do k = lpw(iw), mza
                   ! Convert mol/ sec emissions to ppmV / sec
                  vdemis_gc(k,iw,LTNG_MAP) = vdemis_gc(k,iw,LTNG_MAP) + vdemis_lt(k,iw) * cnvtc(k)
               enddo

            endif

            ! Soil NOx emissions

            if (isfcl == 1 .and. SN_MAP /= 0) then

               ! Loop over land cells beneath this atmospheric grid column

               do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
                  iwsfc = itab_w(iw)%iwsfc(jsfc)
                  jasfc = itab_w(iw)%jasfc(jsfc)

                  kw        = itab_wsfc(iwsfc)%kwatm(jasfc)

                  arc       = itab_wsfc(iwsfc)%arc      (jasfc)  ! Use one of these
                  arcoarsfc = itab_wsfc(iwsfc)%arcoarsfc(jasfc)  ! Use one of these
                  arcoariw  = itab_wsfc(iwsfc)%arcoariw (jasfc)  ! Use one of these
                  arcoarkw  = itab_wsfc(iwsfc)%arcoarkw (jasfc)  ! Use one of these

                  iland = iwsfc - omland

!BOB: What are units of [soilnox(iland) * cnvtc(kw)] and vdemis_gc?  In the following,
!     soilnox needs a multiplier of arc, arcoarsfc, arcoarkw, or arcoariw.

                  ! Convert mol / sec emissions to ppmV / sec
                  vdemis_gc(kw,iw,SN_MAP) = vdemis_gc(kw,iw,SN_MAP) + soilnox(iland) * cnvtc(kw)
               enddo
            endif

            ! MEGAN vegetation emissions

            if (isfcl == 1) then

               ! No non-reacting NH3 here?
               do s = 1, n_gc_emis

                  n = mgn_2_gc_map(s)
                  if (n /= 0) then

                     ! Loop over land cells beneath this atmospheric grid column

                     do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
                        iwsfc = itab_w(iw)%iwsfc(jsfc)
                        jasfc = itab_w(iw)%jasfc(jsfc)

                        kw        = itab_wsfc(iwsfc)%kwatm(jasfc)

                        arc       = itab_wsfc(iwsfc)%arc      (jasfc)  ! Use one of these
                        arcoarsfc = itab_wsfc(iwsfc)%arcoarsfc(jasfc)  ! Use one of these
                        arcoariw  = itab_wsfc(iwsfc)%arcoariw (jasfc)  ! Use one of these
                        arcoarkw  = itab_wsfc(iwsfc)%arcoarkw (jasfc)  ! Use one of these

                        iland = iwsfc - omland

!BOB: What are units of [megan_emis(iland,n) * cnvtc(kw)] and vdemis_gc?  In the following,
!     megan_emis needs a multiplier of arc, arcoarsfc, arcoarkw, or arcoariw.

                        ! Convert mol / sec emissions to ppmV / sec
                        vdemis_gc(kw,iw,s) = vdemis_gc(kw,iw,s) + megan_emis(iland,n) * cnvtc(kw)

                     enddo

                  endif
               enddo
            endif

         enddo
         !$omp end do
         !$omp end parallel

! Aerosol emissions - all units conversions done in GET_AERO_EMIS for aerosols
! GET_AERO_EMIS is in f90 module AERO_EMIS

         if ( n_ae_spc .gt. 0 ) then
            CALL GET_AERO_EMIS ( nvars3d_emis, vname3d_emis, edgar_kmax )
         endif

      END SUBROUTINE GET_EMIS

      END MODULE EMIS_DEFN
