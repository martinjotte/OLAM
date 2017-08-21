
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
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Module soa_defn

C  Defines aerosol species arrays and parameters required in SOA processing.

C  Contains:
C     Subroutine map_soa
C     Subroutine extract_soa
C     Subroutine update_soa
C     Function findVapor
C     Function findOrgprod
C     Subroutine orgaer

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C HP  03/11/11 Updated monoterpene SOA alphas and Cstars to Carlton et al. 2010 values
C HP  07/24/11 Changed aromatic SOA alphas for consistency with updated reaction counters
C               BNZ, TOL, XYL numbers now match Ng et al. 2007 Atmos. Chem. Phys.
C 08 Jun 12 J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
C 13 Aug 13 H. Pye: Xylene and toluene low-NOx yields switched. Values now
C               follow experimental data of Ng et al. 2007 ACP as shown in Table 3.
C               Values in Table 6 of Ng et al. (previously used) are incorrect.
C 18 Dec 13 G.Sarwar: added orgprod parent names based on RACM2
C 07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
C 21 Jul 14 B.Hutzell: used ifdef statement to make oligomerization
C                      optional because process represented in chemical
C                      mechanism
C 26 Sep 14 H. Pye: Added isoprene + NO3 SOA (see mech.def, no changes
C                   here). When IEPOX uptake present in gas phase for
C                   cb05e51, replace preivous acid enhanced isoprene SOA
C                   with IEPOX uptake SOA now handled as a heterogeneous
C                   rxn. For saprc07tic_ae6i, perform more detailed
C                   IEPOX and MAE uptake and do not do Carlton et al. 2010
C                   acid enhancement. Note that saprc07tic_ae6i is a research
C                   version and it is unclear how duplicative AISO1+AISO2+their oligomers
C                   are with  IEPOX+MAE uptake and their oligomers.
C                   Both pathways occur with sarpc07tic_ae6i. To turn AISO1+AISO2
C                   oligomers off, set "Decay" in vaporspc to 0.0 for SV_ISO1/2.
C 27 Sep 14 H. Pye: Added alkane and PAH SOA (Pye and Pouliot 2012 ES&T)
C 15 Jul 15 G. Sarwar: updated SOA from alkane, PAH, and isoprene for RACM2
C 03/03/16  D. Luecken: added capability for CB6
C 24 Mar 16 G. Sarwar: updated for CB05EH51
C    May 16 B. Murphy, H. Pye: updated treatment of aerosol moments
C
C-----------------------------------------------------------------------
      Implicit None

      Integer, Parameter :: n_vapor_list = 28   ! # of potential partitioning SVOCs
      Integer, Parameter :: n_orgprod    = 13   ! # of ROG rxn counter species
      Integer, Save      :: n_vapor             ! # of actual simulated SOA vapors

      Type vapor_type
         Character( 16 ) :: name
         Real            :: alpha       ! Mass-based stoichiometric coefficients [ug/m^3]/[ug/m^3]
         Real            :: cstar       ! Effective saturation concentrations [ug/m^3] at 298 K
         Real            :: enth        ! Enthalphy of Vaporization [J/mol]
         Character( 16 ) :: soa_name    ! Species name
         Character(  1 ) :: soa_origin  ! A = Anthropogenic; B = Biogenic
         Real            :: soa_decay   ! Factor used in oligomerization
         Character( 16 ) :: drog_name   ! Precursor name
      End Type vapor_type

      Type( vapor_type ), Allocatable, Save :: vaporspc( : )

      Type( vapor_type ), Parameter :: vaporlist( n_vapor_list ) = (/
C
C                   Name       Alpha     CStar    Enth      SOA          Decay    Drog
C                                                          Name    Origin         Name
C                  -------    -------   -------- -------  ------     --   ----   -------
     & vapor_type('SV_ALK1 ',  0.0334,   0.1472,  53.0E3, 'AALK1J ', 'A', 12.0, 'ALKRXN  '),
     & vapor_type('SV_ALK2 ',  0.2164,  51.8775,  53.0E3, 'AALK2J ', 'A', 12.0, 'ALKRXN  '),
     & vapor_type('SV_XYL1 ',  0.0310,   1.3140,  32.0E3, 'AXYL1J ', 'A',  8.0, 'XYLNRXN '),
     & vapor_type('SV_XYL2 ',  0.0900,  34.4830,  32.0E3, 'AXYL2J ', 'A',  8.0, 'XYLNRXN '),
     & vapor_type('SV_TOL1 ',  0.0580,   2.3260,  18.0E3, 'ATOL1J ', 'A',  7.0, 'TOLNRXN '),
     & vapor_type('SV_TOL2 ',  0.1130,  21.2770,  18.0E3, 'ATOL2J ', 'A',  7.0, 'TOLNRXN '),
     & vapor_type('SV_BNZ1 ',  0.0720,   0.3020,  18.0E3, 'ABNZ1J ', 'A',  6.0, 'BNZNRXN '),
     & vapor_type('SV_BNZ2 ',  0.8880, 111.1100,  18.0E3, 'ABNZ2J ', 'A',  6.0, 'BNZNRXN '),
     & vapor_type('SV_PAH1 ',  0.2100,   1.6598,  18.0E3, 'APAH1J ', 'A', 10.0, 'PAHNRXN '),
     & vapor_type('SV_PAH2 ',  1.0700, 264.6675,  18.0E3, 'APAH2J ', 'A', 10.0, 'PAHNRXN '),
     & vapor_type('SV_TRP1 ',  0.1393,  14.7920,  40.0E3, 'ATRP1J ', 'B', 10.0, 'TRPRXN  '),
     & vapor_type('SV_TRP2 ',  0.4542, 133.7297,  40.0E3, 'ATRP2J ', 'B', 10.0, 'TRPRXN  '),
     & vapor_type('MTNO3   ',     0.0,     12.0,  40.0E3, 'AMTNO3J', 'B', 10.0, '        '),
     & vapor_type('ISOPNN  ',     0.0,     8.9,   40.0E3, 'AISOPNNJ','B', 5.0,  '        '),
     & vapor_type('SV_ISO1 ',  0.2320, 116.0100,  40.0E3, 'AISO1J ', 'B',  5.0, 'ISOPRXN '),
     & vapor_type('SV_ISO2 ',  0.0288,   0.6170,  40.0E3, 'AISO2J ', 'B',  5.0, 'ISOPRXN '),
     & vapor_type('SV_SQT  ',  1.5370,  24.9840,  40.0E3, 'ASQTJ  ', 'B', 15.0, 'SESQRXN '),
     & vapor_type('LV_PCSOG',  1.0   ,   1.0e-5,  40.0E3, 'APCSOJ ', 'A', 10.0, 'PCSOARXN'),
     & vapor_type('VLVPO1  ',  0.0000,    1.e-1,  96.0E3, 'ALVPO1J', 'A',  0.0, '        '),
     & vapor_type('VSVPO1  ',  0.0000,    1.e+0,  85.0E3, 'ASVPO1J', 'A',  0.0, '        '),
     & vapor_type('VSVPO2  ',  0.0000,    1.e+1,  74.0E3, 'ASVPO2J', 'A',  0.0, '        '),
     & vapor_type('VSVPO3  ',  0.0000,    1.e+2,  63.0E3, 'ASVPO3J', 'A',  0.0, '        '),
     & vapor_type('VIVPO1  ',  0.0000,    1.e+3,  52.0E3, 'AIVPO1J', 'A',  0.0, '        '),
     & vapor_type('VLVOO1  ',  0.0000,    1.e-2, 107.0E3, 'ALVOO1J', 'A',  0.0, '        '),
     & vapor_type('VLVOO2  ',  0.0000,    1.e-1,  96.0E3, 'ALVOO2J', 'A',  0.0, '        '),
     & vapor_type('VSVOO1  ',  0.0000,    1.e+0,  85.0E3, 'ASVOO1J', 'A',  0.0, '        '),
     & vapor_type('VSVOO2  ',  0.0000,    1.e+1,  74.0E3, 'ASVOO2J', 'A',  0.0, '        '),
     & vapor_type('VSVOO3  ',  0.0000,    1.e+2,  63.0E3, 'ASVOO3J', 'A',  0.0, '        ')/)

C Required species
      Character( 16 ), Private, Parameter :: req_sviso1  = 'SV_ISO1'
      Character( 16 ), Private, Parameter :: req_sviso2  = 'SV_ISO2'
      Character( 16 ), Private, Parameter :: req_mtno3   = 'MTNO3'
      Character( 16 ), Private, Parameter :: req_isopnn  = 'ISOPNN'
      Character( 16 ), Private, Parameter :: req_aeiso1  = 'AISO1'
      Character( 16 ), Private, Parameter :: req_aeiso2  = 'AISO2'
      Character( 16 ), Private, Parameter :: req_aeiso3  = 'AISO3'
      Character( 16 ), Private, Parameter :: req_aeolga  = 'AOLGA'
      Character( 16 ), Private, Parameter :: req_aeolgb  = 'AOLGB'
      Character( 16 ), Private, Parameter :: req_aeorgc  = 'AORGC'
      Character( 16 ), Private, Parameter :: req_amtno3  = 'AMTNO3'
      Character( 16 ), Private, Parameter :: req_aisopnn = 'AISOPNN'

C Indices of required species
      Integer :: iso1_idx
      Integer :: iso2_idx
      Integer :: mtno3_idx
      Integer :: isopnn_idx
      Integer :: amtno3_idx
      Integer :: aisopnn_idx
      Integer :: aiso1_idx
      Integer :: aiso2_idx
      Integer :: aiso3_idx
      Integer :: aolga_idx
      Integer :: aolgb_idx
      Integer :: aorgc_idx

      Real, Allocatable    :: vapor_mw  ( : )   ! Molecular weights
      Real, Allocatable    :: vapor_mwi ( : )   ! 1 / Molecular weights
      Real, Allocatable    :: vapor_conc( : )   ! vapor concentration
      Integer, Allocatable :: drog_map  ( : )   ! drog mapping assignments

      !$omp threadprivate(vapor_conc)

C orgprod species

      Type orgprod_type
         Character( 16 ) :: name
         Character( 16 ) :: parent
         Character( 16 ) :: aerospc      ! Non-volatile species name
         Real            :: alphaH       ! Stoichiometric coefficients for non-volatile aromatic SOA [ug/m^3]/[ug/m^3]
      End Type orgprod_type

C   This default orgprod table is for CB05, parent names are changed in subroutine
C   extract_soa for the SAPRC99 and SAPRC07T mechanisms
C
      Type( orgprod_type ), Save :: orgprod( n_orgprod ) = (/
C                      Name      Parent     Aerospc  AlphaH
C                     -------    -------    ------   ------
     & orgprod_type( 'ALKRXN  ', 'SOAALK ', '      ', 0.000 ),
     & orgprod_type( 'XYLNRXN ', 'XYL    ', '      ', 0.000 ),
     & orgprod_type( 'XYLHRXN ', 'XYL    ', 'AXYL3J', 0.360 ),
     & orgprod_type( 'TOLNRXN ', 'TOL    ', '      ', 0.000 ),
     & orgprod_type( 'TOLHRXN ', 'TOL    ', 'ATOL3J', 0.300 ),
     & orgprod_type( 'BNZNRXN ', 'BENZENE', '      ', 0.000 ),
     & orgprod_type( 'BNZHRXN ', 'BENZENE', 'ABNZ3J', 0.370 ),
     & orgprod_type( 'PAHNRXN ', 'NAPH   ', '      ', 0.000 ),
     & orgprod_type( 'PAHHRXN ', 'NAPH   ', 'APAH3J', 0.730 ),
     & orgprod_type( 'TRPRXN  ', 'TRPRXN ', '      ', 0.000 ),
     & orgprod_type( 'ISOPRXN ', 'ISOPRXN', '      ', 0.000 ),
     & orgprod_type( 'SESQRXN ', 'SESQRXN', '      ', 0.000 ),
     & orgprod_type( 'PCSOARXN', 'PCVOC  ', '      ', 0.000 ) /)


      Real    :: orgprod_mw     ( n_orgprod )  ! Mol weight of reactive organic gases
      Real    :: orgprod_conc   ( n_orgprod )  ! orgprod concentration
      Integer :: orgprod_aeroMap( n_orgprod )  ! orgprod pointers to aerospc
      Integer, Allocatable :: soa_aeroMap    ( : )    ! soa pointers to aerospc

      !$omp threadprivate(orgprod_conc)

C mapping variables used for loading and unloading to CGRID array
      Integer, Allocatable :: vapor_map( : )     ! pointers of vapor species to CGRID
      Integer              :: orgprod_map( n_orgprod ) ! pointers of orgprod species to CGRID

      Logical, Private, Save :: mapped              = .False.
      Logical, Private, Save :: IEPOX_SOA_mechanism = .False. ! IEPOX uptake in mechanism reactions?
      Logical, Private, Save :: RXNS_eflag          = .False. ! error flag for RXNS modules used
      Logical, Private, Save :: SOA_eflag           = .False. ! error flag for soa_defn

      Character( 16 ), Private, Save :: pname = 'Soa_Defn        '
      Integer, Private, Save :: logdev

      Real, Save :: mwpoa                     ! Molecular weight of POA [ g/mol ]

      Real, Save ::  drog_ppm2ug( n_orgprod ) ! [ ppm per ug/m3 ] for ORGPROD at
                                              ! reference temperature and pressure
      logical, save :: isnvpoa = .false.

      Contains

C-----------------------------------------------------------------------

      Subroutine map_soa( )

      Use rxns_data       ! chemical mechanism data
      Use aero_data       ! aero species
      Use cgrid_spcs      ! CGRID mechanism species
      Use utilio_defn
      Use aeromet_data, only: rgasuniv

C Local Variables:
      Character( 132 ) :: xmsg
      Integer         :: n, g, i, spc
      Integer         :: vapor_to_vaporlist( n_vapor_list )

      Real, Parameter :: rgas1  = 1.0 / rgasuniv ! reciprocal of universal gas constant
      Real, Parameter :: tref   = 298.0          ! reference temperature [ K ]
      Real, Parameter :: trefm1 = 1.0 / tref ! inverse of reference temperature

      Real, Parameter :: convfac_298 = 101325.0 * rgas1 * trefm1 ! P/RT at 1 atm and 298 K [ mole/m**3 ]

      if ( .not. mapped ) then

          logdev = init3()

C Before Anything Else, the VaporSpc Array must be populated with
C the species from VaporList that have been identified in the GC
C namelist or in the NR namelist.
          n_vapor = 0
          Do spc = 1,n_vapor_list
             n = 0
             g = 0
             n = index1( vaporlist( spc )%name, n_nr_n2ae, nr_n2ae )
             g = index1( vaporlist( spc )%name, n_gc_g2ae, gc_g2ae )
             If ( n .Gt. 0 .Or. g .Gt. 0 ) Then
                 n_vapor = n_vapor + 1
                 vapor_to_vaporlist( n_vapor ) = spc
             End If
          End Do

          Allocate ( vaporspc   ( n_vapor ))
          Allocate ( vapor_mw   ( n_vapor ))
          Allocate ( vapor_mwi  ( n_vapor ))
!         Allocate ( vapor_conc ( n_vapor ))
          Allocate ( drog_map   ( n_vapor ))
          Allocate ( soa_aeroMap( n_vapor ))
          Allocate ( vapor_map  ( n_vapor ))

          Do spc = 1,n_vapor
             vaporspc( spc ) = vaporlist( vapor_to_vaporlist( spc ) )
          End Do

C Check GC or NR namelist contains precursor for each vapor used if the vapor has a
C precursor
          Do spc = 1, n_vapor
            If( Len_Trim( vaporspc( spc )%drog_name ) .Lt. 1 ) Cycle
            n = index1( vaporspc( spc )%drog_name, n_gc_spc, gc_spc )
            If ( n .Lt. 1 ) Then
              n = index1( vaporspc( spc )%drog_name, n_nr_spc, nr_spc )
              If ( n .Lt. 1 ) Then
!                 xmsg = 'FATAL: SOA Precursor, ' // Trim( vaporspc( spc )%drog_name )
!     &               // '. is not found among the species in GC or NR namelist'
!     &               // ' used. Required for modeling vapor, '
!     &               // Trim( vaporspc( spc )%name )
                SOA_eflag  = .True.
                write(logdev,99903)Trim( vaporspc( spc )%drog_name ),
     &          Trim( vaporspc( spc )%name )
!                Call m3warn( pname, 0, 0, xmsg )
              End If
            End If
         End Do

C Check that the soa species found in the AE namelist have their vapor
C phase as model species
          Do spc = 1,n_vapor_list
             n = 0
             g = 0
             n = index1( vaporlist( spc )%soa_name, n_ae_spc, ae_spc )
             If ( n .Gt. 0 ) Then ! vapor phase
!                 write(logdev,*)'found soa species:',TRIM( ae_spc(n) ),
!     &           ' Search for vapor species:',Trim( vaporlist( spc )%name )
                 g = index1( vaporlist( spc )%name, n_nr_spc, nr_spc )
                 If( g .Lt. 1 ) Then
                    g = index1( vaporlist( spc )%name, n_gc_spc, gc_spc )
                    If( g .Lt. 1 ) Then
!                      xmsg = 'FATAL: SOA Vapor, ' // Trim( vaporlist( spc )%name )
!     &                     // ', is not found among the species in the NR or GC'
!     &                     // ' namelist used. Required for modeling '
!     &                     // Trim ( vaporlist( spc )%soa_name ) // '.'
                      SOA_eflag  = .True.
                      write(logdev,99902)Trim( vaporlist( spc )%name ), 
     &                Trim ( vaporlist( spc )%soa_name )
!                      Call m3warn( pname, 0, 0, xmsg )
                    End If
                 End If
!                 If( g .Gt. 0 ) Then
!                    write(logdev,*)'found vapor species:', Trim( vaporlist( spc )%name )
!                End If
             End If
         End Do

C Build mapping to CGRID for each vapor species
C Add option to have species in GC.nml (for real semivolatiles)
         vapor_mw  = 0.0
         vapor_mwi = 0.0
         vapor_map = 0
         Do spc = 1, n_vapor
            If( vaporspc( spc )%alpha .Gt. 0.0 ) Then
               n = index1( vaporspc( spc )%name, n_nr_n2ae, nr_n2ae )
               If ( n .Ne. 0 ) Then
                  vapor_map( spc ) = nr_strt - 1 + nr_n2ae_map( n )
                  vapor_mw ( spc ) = nr_molwt ( nr_n2ae_map( n ) )
                  vapor_mwi( spc ) = nr_molwti( nr_n2ae_map( n ) )
               Else
                  xmsg = 'FATAL: ' // Trim( vaporspc( spc )%name )
     &                 // ' is not found in the N2AE values of NR namelist'
     &                 // ' used. '
                  SOA_eflag  = .True.
                  Call m3warn( pname, 0, 0, xmsg )
!                 Call m3exit( pname, 0, 0, xmsg, xstat3 )
               End If
            Else ! vaporspc is a real species in GC (alpha=0 )
               n = index1( vaporspc( spc )%name, n_gc_g2ae, gc_g2ae )
               If ( n .Ne. 0 ) Then
                  vapor_map( spc ) = gc_strt - 1 + gc_g2ae_map( n )
                  vapor_mw ( spc ) = gc_molwt ( gc_g2ae_map( n ) )
                  vapor_mwi( spc ) = gc_molwti( gc_g2ae_map( n ) )
               Else
                  xmsg = 'FATAL: ' // Trim( vaporspc( spc )%name )
     &                 // ' is not found in the G2AE values of GC namelist'
     &                 // ' used. '
                  SOA_eflag  = .True.
                  Call m3warn( pname, 0, 0, xmsg )
!                 Call m3exit( pname, 0, 0, xmsg, xstat3 )
               End If
            End If
         End Do

#ifdef verbose_soa
         write(logdev,'(a16,1x,a4,1x,a7)')"    vaporspc    ","INDX","Mol.Wei"
         Do spc = 1, n_vapor
           write(logdev,'(a16,1x,i4,1x,f7.2)')vaporspc( spc )%name, vapor_map( spc ), 
     &                                        vapor_mw( spc )
         End Do
#endif

C define orgprod parent names based on mechanism dependence
C if they differ from above declarations and
C determine if gas-phase IEPOX uptake occurs for the mechanism
            If ( INDEX( MECHNAME, 'SAPRC07TIC_AE6I' ) .Gt. 0 ) Then
               orgprod( 2 )%parent = 'MXYL'
               orgprod( 3 )%parent = 'MXYL'
               orgprod( 4 )%parent = 'TOLUENE'
               orgprod( 5 )%parent = 'TOLUENE'
               orgprod( 8 )%parent = 'NAPHTHAL'
               orgprod( 9 )%parent = 'NAPHTHAL'
               IEPOX_SOA_mechanism = .True.
            Else If ( INDEX( MECHNAME, 'SAPRC07TB' ) .Gt. 0 ) Then
               orgprod( 2 )%parent     = 'MXYL'
               orgprod( 3 )%parent     = 'MXYL'
               orgprod( 4 )%parent     = 'TOLUENE'
               orgprod( 5 )%parent     = 'TOLUENE'
               orgprod( 8 )%parent     = 'NAPHTHAL'
               orgprod( 9 )%parent     = 'NAPHTHAL'
               IEPOX_SOA_mechanism = .True.
            Else If ( INDEX( MECHNAME, 'CRIe' ) .Gt. 0 ) Then
               orgprod( 2 )%parent     = 'MXYL'
               orgprod( 3 )%parent     = 'MXYL'
               orgprod( 4 )%parent     = 'TOLUENE'
               orgprod( 5 )%parent     = 'TOLUENE'
               orgprod( 8 )%parent     = 'NAPHTHAL'
               orgprod( 9 )%parent     = 'NAPHTHAL'
            Else If ( INDEX( MECHNAME, 'SAPRC07TC' ) .Gt. 0 ) Then
               orgprod( 2 )%parent     = 'MXYL'
               orgprod( 3 )%parent     = 'MXYL'
               orgprod( 4 )%parent     = 'TOLUENE'
               orgprod( 5 )%parent     = 'TOLUENE'
               orgprod( 8 )%parent     = 'NAPHTHAL'
               orgprod( 9 )%parent     = 'NAPHTHAL'
               IEPOX_SOA_mechanism = .True.
            Else If( INDEX(MECHNAME, 'RACM2') .Gt. 0) Then
! Changed by Golam Sarwar, 07/15/2015 for CMAQv51
               orgprod( 2 )%parent     = 'XYM'
               orgprod( 3 )%parent     = 'XYM'
               orgprod( 4 )%parent     = 'TOL'
               orgprod( 5 )%parent     = 'TOL'
               IEPOX_SOA_mechanism = .True.
            Else If( INDEX(MECHNAME, 'CB05TUCL') .Gt. 0) Then
               orgprod( 2 )%parent = 'XYLMN'
               orgprod( 3 )%parent = 'XYLMN'
            Else If( INDEX(MECHNAME, 'CB05E51') .Gt. 0) Then
               orgprod( 2 )%parent = 'XYLMN'
               orgprod( 3 )%parent = 'XYLMN'
               IEPOX_SOA_mechanism = .True.
            Else If( INDEX(MECHNAME, 'CB6R3') .Gt. 0) Then
               orgprod( 2 )%parent = 'XYLMN'
               orgprod( 3 )%parent = 'XYLMN'
               IEPOX_SOA_mechanism = .True.
            Else If( INDEX(MECHNAME, 'CB05EH51') .Gt. 0) Then
               orgprod( 2 )%parent = 'XYLMN'
               orgprod( 3 )%parent = 'XYLMN'
               IEPOX_SOA_mechanism = .True.
            Else 
               xmsg =  'FATAL: The ' // MECHNAME // ' mechanism '
     &             // ' is not supported in this version of the model.'
     &             // ' Change CCTM build script to use a mechanism based on '
     &             // ' mechanisms available in currect CMAQ repository or '
     &             // ' Change MECHNAME to one of the repository mechanisms.'
               RXNS_eflag  = .True.
               Call m3warn( pname, 0, 0, xmsg )
!              Call m3exit( pname, 0, 0, xmsg, xstat3 )
            End If

! SOA_DEFN NO LONGER HAS VERSION CONSTRAINTS ON SPECIES
! This SOA_DEFN works for AE6 and AE6i versions of aero
!            If( INDEX( MECHNAME, '_AE6_' ) .EQ. 0 ) then
!                  xmsg =  'Using AERO6 option and the Mechanism, '
!     &                // Trim( MECHNAME ) // ', is not based on AERO6 '
!                 Call m3exit( pname, 0, 0, xmsg, xstat3 )
!            End If


C Build mapping to CGRID for each ORGPROD species and find MW of parent

         orgprod_map = 0
         Do spc = 1, n_orgprod
            n = index1( orgprod( spc )%name, n_gc_g2ae, gc_g2ae )
            If ( n .Ne. 0 ) Then
               orgprod_map( spc ) = gc_strt - 1 + gc_g2ae_map( n )
               n = index1( orgprod( spc )%parent, n_gc_spc, gc_spc )
               If ( n .Ne. 0 ) Then
                  orgprod_mw( spc ) = gc_molwt( n )
               Else
                  xmsg = 'FATAL: Species ' // orgprod( spc )%parent
     &                 // ' in ORGPROD parent name is not in the G2AE values'
     &                 // ' of the GC namelist for the ' // TRIM ( MECHNAME ) 
     &                 // ' mechanism '
                  SOA_eflag  = .True.
                  Call m3warn( pname, 0, 0, xmsg )
!                 Call m3exit( pname, 0, 0, xmsg, xstat3 )
               End If
            End If
         End Do

#ifdef verbose_soa
         write(logdev,'(2(a16),1x,a4,1x,a7)')"    orgprod    ","    parent    ","INDX","Mol.Wei"
         Do spc = 1, n_orgprod
           If ( orgprod_map( spc ) .Lt. 1 )Cycle
           write(logdev,'(2(a16),1x,i4,1x,f7.2)')orgprod( spc )%name,orgprod( spc )%parent, 
     &                                        orgprod_map( spc ),orgprod_mw( spc ) 
         End Do
#endif

C Find indices of required species
         iso1_idx   = findVapor( req_sviso1 )
         iso2_idx   = findVapor( req_sviso2 )
         aiso1_idx  = findAero( req_aeiso1, .True. )
         aiso2_idx  = findAero( req_aeiso2, .True. )
         aiso3_idx  = findAero( req_aeiso3, .True. )
         aolga_idx  = findAero( req_aeolga, .True. )
         aolgb_idx  = findAero( req_aeolgb, .True. )
         aorgc_idx  = findAero( req_aeorgc, .True. )

#ifdef verbose_soa
      ! Mandatory Species
         Write( logdev,'( /5x, a )' ) 'soa_defn: required species'
         Write( logdev,'( 5x, a, i4 )' ) 'iso1_idx:   ', iso1_idx 
         Write( logdev,'( 5x, a, i4 )' ) 'iso2_idx:   ', iso2_idx 
         Write( logdev,'( 5x, a, i4 )' ) 'aiso1_idx:  ', aiso1_idx
         Write( logdev,'( 5x, a, i4 )' ) 'aiso2_idx:  ', aiso2_idx
         Write( logdev,'( 5x, a, i4 )' ) 'aiso3_idx:  ', aiso3_idx
         Write( logdev,'( 5x, a, i4 )' ) 'aolga_idx:  ', aolga_idx
         Write( logdev,'( 5x, a, i4 )' ) 'aolgb_idx:  ', aolgb_idx
         Write( logdev,'( 5x, a, i4 )' ) 'aorgc_idx:  ', aorgc_idx
#endif
         ! Look for Additional Isoprene Product Species if AE6i is being
         ! used
         If ( INDEX( MECHNAME,"SAPRC07TIC_AE6I" ) .Ne. 0 ) Then
           mtno3_idx  = findVapor( req_mtno3  )
           isopnn_idx = findVapor( req_isopnn )
           amtno3_idx  = findAero( req_amtno3, .True. )
           aisopnn_idx = findAero( req_aisopnn,.True. )
#ifdef verbose_soa
           Write( logdev,'( /5x, a )' )
     &      'SAPRC07TIC_AE6I soa_defn: required species'
           Write( logdev,'( 5x, a, i4 )' ) 'mtno3_idx:   ', mtno3_idx  
           Write( logdev,'( 5x, a, i4 )' ) 'isopnn_idx:  ', isopnn_idx 
           Write( logdev,'( 5x, a, i4 )' ) 'amtno3_idx:  ', amtno3_idx 
           Write( logdev,'( 5x, a, i4 )' ) 'aisopnn_idx: ', aisopnn_idx
#endif
         End If


C Build mapping array of Semivolatile Vapors (vaporspc) to the Aerosol
C Species they correspond to (aerospc is used inside findAero)
         soa_aeroMap = 0
#ifdef verbose_soa
           Write( logdev,'( /5x, a )' )
     &      'soa_defn: vapor species'
#endif     
         Do spc = 1, n_vapor
            If ( vaporspc( spc )%soa_name .Ne. ' ' ) Then
               soa_aeroMap( spc ) = findAero( vaporspc( spc )%soa_name, .True. )
#ifdef verbose_soa
           Write( logdev,'( 5x, a20, i4 )' ) Trim(vaporspc( spc )%soa_name) // ':', 
     &                                       soa_aeroMap( spc )
#endif
               If( soa_aeroMap( spc ) .Lt. 1 ) write(logdev,99904) Trim( vaporspc( spc )%name )
            End If
         End Do

C Build Mapping array of Reactive variables that track production rates
C of condensable gases (orgprod) to any nonvolatile species they yield
C (aerospc is used inside findAero)
         orgprod_aeroMap = 0
#ifdef verbose_soa
           Write( logdev,'( /5x, a )' )
     &      'soa_defn: orgprod_aero species'
#endif     
         Do spc = 1, n_orgprod
            If ( orgprod( spc )%aerospc .Ne. ' ' ) Then
               orgprod_aeroMap( spc ) = findAero( orgprod( spc )%aerospc, .True. )
#ifdef verbose_soa
           Write( logdev,'( 5x, a20, i4 )' ) Trim(orgprod( spc )%aerospc) // ':', 
     &                                       orgprod_aeroMap( spc )
#endif
               If( orgprod_aeroMap( spc ) .Lt. 1 ) write(logdev,99905) Trim( orgprod( spc )%name )
            End If
         End Do

C Build mapping array of sources of Semivolatile Vapors (vaporspc) to the
C Reactive variables that store their production rates (orgprod)
         drog_map = 0
#ifdef verbose_soa
           Write( logdev,'( /5x, a )' )
     &      'soa_defn: drog species'
#endif     
         Do spc = 1, n_vapor
            If ( vaporspc( spc )%drog_name .Ne. ' ' ) Then
               drog_map( spc ) = findOrgprod( vaporspc( spc )%drog_name )
#ifdef verbose_soa
           Write( logdev,'( 5x, a20, i4 )' ) Trim(vaporspc( spc )%drog_name) // ':', 
     &                                       drog_map( spc )
#endif
               If( drog_map( spc ) .Lt. 1 ) write(logdev,99905) Trim( orgprod( spc )%name )
            End If
         End Do

         If( AE_eflag .Or. SOA_eflag .Or. RXNS_eflag )Then
             If( AE_eflag .Or. SOA_eflag )Write(logdev,99901)Trim( mechname )
             xmsg = 'The FATAL errors found in namelist used. Check '
     &          //  'the log of exiting processor if more details are needed.'
             Call m3exit( pname, 0, 0, xmsg, xstat3 )
         End If 

C     Set unit conversion

         Do i = 1, n_orgprod
            drog_ppm2ug( i ) = orgprod_mw( i ) * convfac_298
         End Do

         mwpoa = aerospc_mw( apoc_idx )

         isnvpoa = (INDEX( MECHNAME, 'NVPOA' ) .Ne. 0)

         mapped = .True.

      End If     ! mapping condition

99901 Format(/'FATAL Error(s) found in the namelist(s) used. Check that '
     &     /  'these namelists contain the data as the files in '
     &     /  'the respository version of the mechanism: ', a / )
99902 Format('FATAL: SOA Vapor, ', a, ', is not found among the species in '
     &        / 'the NR or GC namelist used. Required for modeling ', a / )
99903 Format( 'FATAL: SOA Precursor, ', a,  '. is not found among the species',
     &         / 'in GC or NR namelist used. Required for modeling vapor, ', a / )
99904 Format(1X,'aerosol product from namelist species: ', a /)
99905 Format(1X,'aerosol by-product from namelist species : ', a /)
C99906 Format(1X,'precursor from namelist species : ', a /)
  
      End Subroutine Map_soa

C-----------------------------------------------------------------------

      Subroutine extract_soa( conc )

C  Extracts the required soa data from CGRID into the conc array.

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C SH  03/10/11 Renamed met_data to aeromet_data
C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN
C HP  09/27/14 alk_factor removed, updated for alkane/PAH SOA.
C              Conversion of reacted alkane to dodecane equivalent
C              is handled in mech.def. A factor of 0.47 is in use and
C              reflects the fact that alkane SOA precursor
C              emissions are dominated by compounds smaller than dodecane.
C-----------------------------------------------------------------------

      Use aeromet_data          ! air and met variables

      Implicit None

C Arguments:
      Real, Intent( In ) :: conc( : )

C Local Variables:
      Real            :: gasconv
      Real            :: vtmp
      Integer         :: n
      Integer         :: spc

      if (.not. allocated(vapor_conc)) allocate( vapor_conc( n_vapor ) )

C Compute gas conversion constant
      gasconv = airdens * inv_mwair

C Copy grid cell concentrations of vapor species
      vapor_conc = 0.0

      Do spc = 1, n_vapor
         n = vapor_map( spc )
         If ( n .Ne. 0 ) Then
            vtmp = gasconv * vapor_mw( spc )
            vapor_conc( spc ) = Max( conc( n ) * vtmp, min_gasconc )
         End If
      End Do

C Copy grid cell concentrations of ORGPROD species
      orgprod_conc = 0.0
      Do spc = 1, n_orgprod
         n = orgprod_map( spc )
         If ( n .Ne. 0 ) Then
            orgprod_conc( spc ) = Max( conc( n ), min_gasconc )
         End If
      End Do

      End Subroutine extract_soa

C-----------------------------------------------------------------------

      Subroutine update_soa( conc )

C  Populates CGRID from the conc array with updated SOA values.

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C SH  03/10/11 Renamed met_data to aeromet_data
C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN
C-----------------------------------------------------------------------

      Use aeromet_data    ! air and met variables
      Use utilio_defn

      Implicit None

C Arguments:
      Real, Intent( Out ) :: conc( : )

C Local Variables:
      Character( 80 ) :: xmsg
      Real            :: gasconvi
      Real            :: vtmp
      Integer         :: n
      Integer         :: spc

      If ( .Not. mapped ) Then
         xmsg = 'CGRID Species has not been mapped'
         Call m3exit( pname, 0, 0, xmsg, xstat3 )
      End If

C Compute gas conversion constant
      gasconvi = mwairkg / airdens

C Copy vapor_conc back to grid cell concentrations
      Do spc = 1, n_vapor
         n = vapor_map( spc )
         vtmp = vapor_mwi( spc ) * gasconvi
         conc( n ) = Max ( vapor_conc( spc ) * vtmp, min_gasconc )
      End Do

      End Subroutine update_soa

C-----------------------------------------------------------------------

      Function findVapor( vname ) Result ( ndx )

C  Finds the index of 'required' semivolatile species in the vaporspc list

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN
C-----------------------------------------------------------------------

      Use utilio_defn

      Implicit None

C Arguments:
      Character( 16 ) :: vname
      Integer ndx

C Local Variables:
      Character( 80 ) :: xmsg
      Integer         :: spc

      ndx = 0
      Do spc = 1, n_vapor
         If ( vaporspc( spc )%name .Eq. vname ) Then
            ndx = spc
            Return
         End If
      End Do

!      xmsg = 'Required Species ' // Trim( vname ) // ' not found in [vapor names] array'
      xmsg = 'FATAL:' // Trim( vname )
     &    // ' is not found in the G2AE or N2AE values of GC or namelists'
      SOA_eflag  = .True.
      Call m3warn( pname, 0, 0, xmsg )
!     Call m3exit( pname, 0, 0, xmsg, xstat3 )

      End Function findVapor

C-----------------------------------------------------------------------
      Function findOrgprod( vname ) result ( ndx )

C  Finds the index of 'required' organic product species in the orgprod list

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN
C-----------------------------------------------------------------------

      Use utilio_defn        ! i/o api

      Implicit None

C Arguments:
      Character( 16 ) :: vname
      Integer ndx

C Local Variables:
      Character( 80 ) :: xmsg
      Integer         :: spc

      Do spc = 1, n_orgprod
         If ( orgprod( spc )%name .Eq. vname ) Then
            ndx = spc
            Return
         End If
      End Do

!      xmsg = 'Required Species ' // Trim( vname )
!     &     // ' Not found in [orgprod names] array'
      xmsg = 'FATAL:' // Trim( vname )
     &    // ' is not found in the G2AE or N2AE values of GC or namelists'
      SOA_eflag  = .True.
      Call m3warn( pname, 0, 0, xmsg )
!     Call m3exit( pname, 0, 0, xmsg, xstat3 )

      End Function findOrgprod

C-----------------------------------------------------------------------
      Subroutine orgaer( dt, layer )

C Updates CGRID via several pathways for secondary organic aerosol (SOA)
C formation, as recommended by Edney et al. (2007).  These include SOA
C formation from isoprene, monoterpenes, sesquiterpenes, long alkanes, and
C aromatics (incl. benzene).

C Input includes the concentrations of reactive organic gases (ROG)
C that were oxidized during the time step (ORGPROD), the vapor-phase
C concentration of each semi-volatile organic compound, the
C concentration of each SOA species, and the concentration of primary
C organic aerosol (all concentrations are stored in the CBLK array).
C Output includes updated concentrations of SOA species, vapor-phase
C semi-volatile organic compounds, and moments of the accumulation
C mode.  The geometric mean diameter of the accumulation mode is also
C updated.  All SOA formation is restricted to the accumulation mode.

C This code relies on 12 counter species to be incorporated in the
C gas-phase chemical mechanisms to track the amounts of individual
C ROG that reacted during the current time step (i.e., NPREC=10).
C The arrays of length = NPREC include:
C       (1) "long" alkanes  (ALKRXN)
C       (2) low-yield aromatics, high-NOx pathway (XYLNRXN)
C       (3) low-yield aromatics, low-NOx pathway (XYLHRXN)
C       (4) high-yield aromatics, high-NOx pathway (TOLNRXN)
C       (5) high-yield aromatics, low-NOx pathway (TOLHRXN)
C       (6) benzene, high-NOx pathway (BNZNRXN)
C       (7) benzene, low-NOx pathway (BNZHRXN)
C       (8) monoterpenes (TRPRXN)
C       (9) isoprene (ISOPRXN)
C      (10) sesquiterpenes (SESQRXN)
C      (11) PAHs/naphthalene, high-NOx pathway (PAHNRXN)
C      (12) PAHs/naphthalene, low-NOx pathway (PAHHRXN)

C In total, 15 organic species are allowed to partition between the
C vapor and particulate phases (i.e., NCVAP=12). The arrays of
C length = NCVAP include:
C           alkane (2 semi-volatile products)
C           low-yield aromatics, high-NOx pathway (2 products)
C           high-yield aromatics, high-NOx pathway (2 products)
C           benzene, high-NOx pathway (2 products)
C           monoterpenes (2 products)
C           isoprene (2 products)
C           sesquiterpenes (1 product)
C           PAHs/naphthalene, high-NOx pathway (2 products)

C Equilibrium partitioning calculations are based on the absorptive
C partitioning model of Pankow (1994) that was extended by Odum et
C al. (1996).  Saturation vapor pressures (cstar) and mass-based
C stoichiometric yield coefficients (alpha) are obtained from smog-
C chamber studies.  Saturation vapor pressures are modified as a
C function of temperature using eqn 6 of Sheehan & Bowman (2001).

C If the pre-existing organic aerosol concentration is zero,
C gas/particle equilibrium is established only after the organic gas
C concentration reaches the threshold value defined in eqn 9 of
C Schell et al. (2001).  Until this threshold value is reached,
C organic vapors do not partition to the particle phase.  Once the
C organic gas/particle equilibrium has been established, gas and
C particle-phase concentrations of each condensible species are
C calculated iteratively using a globally convergent variation of
C Newton's method (SUBROUTINE NEWT), as described in eqn 8 of Schell
C et al. (2001).

C In addition to the various pathways of semi-volatile SOA formation
C treated in previous versions of the model, four types of non-
C volatile SOA are considered here:
C   (1) aromatic and PAH-derived SOA under low-NOx conditions
C   (2) oligomerization of all particle-phase semi-volatile material
C   (3) SOA formed by in-cloud oxidation  (SUBROUTINE AQCHEM)
C   (4) isoprene IEPOX-derived SOA under acidic conditions (AEROSOL_CHEMISTRY)

C Previous code revision history:
C   Originally coded August 1, 2001 by Dr. Francis S. Binkowski

C   Revised April 4, 2003 by Gerald Gipson to allow for evaporation
C   of organics from aerosols. Now total vapor + aerosol phase is
C   repartitioned at each time step and totorgsw ( Mo ) does not
C   include oldsoa.

C   Revised July 14, 2003 by Dr. Prakash V. Bhave
C   - changed cstar(2,3) from 10.103 & 90.925 to 111.11 & 1000.0
C     because smog chamber data of Kalberer et al. were collected
C     at 298 K (not 310 K, as was previously assumed)
C   - changed mw_vap(9,10) from 184 g/mol to 177 g/mol, to be
C     consistent with mwsoa_b
C   - modified threshold criteria for establishing gas/particle
C     equilibrium by removing the loose criterion involving "mtot"
C   - changed variable names to reflect that the combined vapor +
C     aerosol concentrations are now being repartitioned during
C     each time step (not just the newly formed SVOC's)
C   - added documentation and removed extraneous lines of code

C   Revised December 4, 2003 by Dr. Francis S. Binkowski
C   - output variables ORGRATE and ORGBRATE removed and replaced
C     by SOA_A and SOA_B, the newly equilibrated values of
C     Anthropogenic and Biogenic SOA, respectively.  These are non-
C     negative values.
C   - variable jj also removed

C   Revised January 8, 2004 by Dr. Prakash V. Bhave
C   - removed the output variable YIELD.  It has no physical meaning
C     after the 12/04/2003 revisions.

C   Revised January 12, 2004 by Dr. Chris G. Nolte
C   - for computational efficiency, modified the initial caer guess
C     used as input to NEWT.  If NEWT returns check .eq. true, then
C     NEWT is called again with a guess of caer = 0.5*ctot
C   - removed ITS parameter from NEWT call vector
C   - fixed bug where concentrations less than TOLMIN (i.e., 1.0E-12)
C     were reset to 1.0e-30
C   - removed extraneous code related to "Pandis method" of SVOC
C     partitioning when threshold criterion is not met (i.e.,
C     insufficient organic matter to establish gas/particle
C     equilibrium)  ** results unaffected by this change
C
C   Revised September 7, 2007 by Dr. Sergey L. Napelenok
C   - Replaced old SOA species (SOA_A, SOA_B) with an array of
C     precursor-specific SOA species.  Replaced OLDSOA_A and OLDSOA_B
C     with an array (OLDSOA).  Updated call vector accordingly.
C   - Deleted nole* and nbio* variables (now obsolete)
C   - Increased the dimension of several arrays to accommodate new
C     SOA precursors (benzene, sesquiterpenes) and pathways (low-NOx,
C     acid-catalyzed, oligomers, in-cloud)
C
C   Revised November 29, 2007 by Dr. Prakash V. Bhave
C   - Renamed subroutine from ORGAER3 to ORGAER5
C   - Modified M2 and M3 calculations to account for the updated
C     definition of DRY aerosol (which now includes non-volatile SOA)
C   - Updated Hvap and corresponding cstar values
C   - Added parameters for SOA from isoprene and sesquiterpenes
C   - Updated aromatic SOA scheme to include semi-volatile and non-
C     volatile products that form under high-NOx and low-NOx
C     conditions, respectively
C   - Added oligomerization process
C   - Added enhancement of isoprene SOA under acidic conditions
C
C   Revised June 2, 2008 by Dr. Prakash V. Bhave
C   - Changed h_vap of benzene SOA to match that of toluene SOA, based
C     on consultation with Dr. Ed Edney and Dr. Tad Kleindienst.

C   Revised June 5, 2008 by Drs. Prakash Bhave and Sergey Napelenok
C   - Simplified the code for conserving low-volatility isoprene
C     oxidation products and removed a minor bug in the acid-induced
C     isoprene SOA calculation.
C
C   Revised September 9, 2008 by Dr. Prakash V. Bhave
C   - Increased alpha values for SV_TRP1, SV_TRP2, and SV_SQT by a
C     factor of 1.3 to correct for the implicit assumption of unit
C     density in those SOA yield parameters.
C   - Reduced SOA/SOC ratio of AISO1 and AISO2 from 2.5 to 1.6, and
C     increased SOA/SOC ratio of AISO3 from 2.5 to 2.7.  Accordingly,
C     the molar masses of AISO1 and AISO2 were decreased to 96 g/mol
C     and the molar mass of AISO3 was increased to 162.
C
C   Revised September 26, 2014 by Dr. Havala Pye
C   - Removed previous acid enhanced isoprene SOA. Acid catalyzed
C     isoprene SOA now follows Pye et al. 2013 ES&T uptake of IEPOX.
C     See AEROSOL_CHEMISTRY and mech.def for IEPOX SOA.
C   - Allowed for alternate method to NEWT for solving partitioning
C     equations

C References:
C   1. Edney, E.O., T.E. Kleindienst, M. Lewandowski, and J.H.
C      Offenberg, Updated SOA chemical mechanism for the Community
C      Multi-Scale Air Quality model, EPA 600/X-07/025, U.S. EPA,
C      Research Triangle Park, NC, 2007.

C   2. Pankow, J. F., An absorption model of gas/particle partitioning
C      of organic compounds in the atmosphere, Atmos. Environ., Vol 28,
C      No 2, 185-188, 1994.

C   3. Odum, J. R., T. Hoffmann, F. Bowman, D. Collins, R. C. Flagan,
C      and J. H. Seinfeld, Gas/particle partitioning and secondary
C      organic aerosol yields, Environ. Sci. Technol., Vol 30, No 8,
C      2580-2585, 1996.

C   4. Sheehan, P. E. and F. M. Bowman, Estimated effects of temperature
C      on secondary organic aerosol concentrations, Environ. Sci.
C      Technol., Vol 35, No 11, 2129-2135, 2001.

C   5. Schell, B., I. J. Ackermann, H. Hass, F. S. Binkowski, and
C      A. Abel, Modeling the formation of secondary organic aerosol
C      within a comprehensive air quality modeling system, J. Geophys.
C      Res., Vol 106, No D22, 28275-28293, 2001.

C   6. Strader, R., F. Lurmann, and S. N. Pandis, Evaluation of
C      secondary organic aerosol formation in winter, Atmos. Environ.,
C      Vol 33, 4849-4863, 1999.

C   7. Ng, N. L., J. H. Kroll, A. W. H. Chan, P. S. Chhabra, R. C.
C      Flagan, and J. H. Seinfeld, Secondary organic aerosol formation
C      from m-xylene, toluene, and benzene, Atmos. Chem. Phys., Vol 7,
C      3909-3922, 2007a.

C   8. Griffin, R. J., D. R. Cocker III, R. C. Flagan, and J. H.
C      Seinfeld, Organic aerosol formation from the oxidation of
C      biogenic hydrocarbons, J. Geophys. Res., Vol 104, No D3,
C      3555-3567, 1999.

C   9. Bian, F. and F. M. Bowman, Theoretical method for lumping
C      multicomponent secondary organic aerosol mixtures, Environ.
C      Sci. Technol., Vol 36, No 11, 2491-2497, 2002.

C  10. Offenberg, J. H., T. E. Kleindienst, M. Jaoui, M. Lewandowski,
C      and E. O. Edney, Thermal properties of secondary organic
C      aerosols, Geophys. Res. Lett., Vol 33, L03816, doi:10.1029/
C      2005GL024623, 2006.

C  11. Bahreini, R., M. D. Keywood, N. L. Ng, V. Varutbangkul, S. Gao,
C      R. C. Flagan, J. H. Seinfeld, D. R. Worsnop, and J. L. Jimenez,
C      Measurements of secondary organic aerosol from oxidation of
C      cycloalkenes, terpenes, and m-xylene using an Aerodyne aerosol
C      mass spectrometer, Environ. Sci. Technol., Vol 39, 5674-5688,
C      2005.

C  12. Alfarra, M. R., D. Paulsen, M. Gysel, A. A. Gaforth, J. Dommen,
C      A. S. H. Prevot, D. R. Worsnop, U. Baltensperger, and H. Coe,
C      A mass spectrometric study of secondary organic aerosols formed
C      from the photooxidation of anthropogenic and biogenic precursors
C      in a reaction chamber, Atmos. Chem. Phys., Vol 6, 5279-5293,
C      2006.

C  13. Ng, N. L., P. S. Chhabra, A. W. H. Chan, J. D. Surratt, J. H.
C      Kroll, A. J. Kwan, D. C. McCabe, P. O. Wennberg, A. Sorooshian,
C      S. M. Murphy, N. F. Dalleska, R. C. Flagan, and J. H. Seinfeld,
C      Effect of NOx level on secondary organic aerosol (SOA) formation
C      from the photooxidation of terpenes, Atmos. Chem. Phys., Vol 7,
C      5159-5174, 2007b.

C  14. Kostenidou, E., R. K. Pathak, and S. N. Pandis, An algorithm for
C      the calculation of secondary organic aerosol density combining
C      AMS and SMPS data, Aerosol Sci. Technol., Vol 41, 1002-1010,
C      2007.

C  15. Offenberg, J. H., C. W. Lewis, M. Lewandowski, M. Jaoui, T. E.
C      Kleindienst, and E. O. Edney, Contributions of toluene and
C      alpha-pinene to SOA formed in an irradiated toluene/alpha-pinene/
C      NOx/air mixture: comparison of results using 14C content and SOA
C      organic tracer methods, Environ. Sci. Technol., Vol 41, 3972-
C      3976, 2007.

C  16. Henze, D. K. and J. H. Seinfeld, Global secondary organic aerosol
C      from isoprene oxidation, Geophys. Res. Lett., Vol 33, L09812,
C      doi:10.1029/2006GL025976, 2006.

C  17. Kleindienst, T. E., M. Jaoui, M. Lewandowski, J. H. Offenberg,
C      C. W. Lewis, P. V. Bhave, and E. O. Edney, Estimates of the
C      contributions of biogenic and anthropogenic hydrocarbons to
C      secondary organic aerosol at a southeastern US location, Atmos.
C      Environ., Vol 41, 8288-8300, 2007.

C  18. Kalberer, M., D. Paulsen, M. Sax, M. Steinbacher, J. Dommen,
C      A. S. H. Prevot, R. Fisseha, E. Weingartner, V. Frankevich,
C      R. Zenobi, and U. Baltensperger, Identification of polymers as
C      major components of atmospheric organic aerosols, Science, Vol
C      303, 1659-1662, 2004.

C  19. Turpin, B. J. and H.-J. Lim, Species contributions to PM2.5 mass
C      concentrations: revisiting common assumptions for estimating
C      organic mass, Aero. Sci. Technol., Vol 35, 602-610, 2001.

C  20. Surratt, J. D., M. Lewandowski, J. H. Offenberg, M. Jaoui, T. E.
C      Kleindienst, E. O. Edney, and J. H. Seinfeld, Effect of acidity
C      on secondary organic aerosol formation from isoprene, Environ.
C      Sci. Technol., Vol 41, 5363-5369, 2007.

C  21. Pye et al., Epoxide pathways improve model prediction of isoprene
C      markers and reveal key role of acidity in aerosol formation,
C      Environ. Sci. Technol., 2013.


C Revision History:
C    First orgaer version was coded in April 2010 by Steve Howard with
C    Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C SH  03/10/11 Renamed met_data to aeromet_data
C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN
C HOTP 05Aug15 Made the calculation for total number of organic moles more robust
C BNM 11/09/15 Added Some Comments to the SOA Scheme
C-----------------------------------------------------------------------

C Key Subroutines/Functions called:  newt, soabisection

      Use aero_data
      Use aeromet_data
      Use utilio_defn
      Use rxns_data, only : MECHNAME

      Implicit None

C Arguments:
      Real    :: dt            ! Synchronization time step [ s ]
      Integer :: layer         ! model layer number

C Local variables:
!     Character( 300 ) :: xmsg
      Integer       :: i, im, indx

C Parameters & variables for adjusting cstar to ambient conditions
      Real, Parameter :: tref   = 298.0          ! reference temperature [ K ]
      Real, Parameter :: trefm1 = 1.0 / tref     ! inverse of reference temperature
      Real, Parameter :: prefm1 = 1.0 / stdatmpa ! inverse of reference pressure
      Real, Parameter :: rgas1  = 1.0 / rgasuniv ! reciprocal of universal gas constant
      Real, Parameter :: kolig  = 0.69314718 / 72000.0  ! 20h half-life of oligomerization rate [ 1/s ]
      Real, Parameter :: olgrat = 2.1            ! SOA/SOC ratio for oligomers
      Real, Parameter :: kacid  = 3.61           ! acid-induced enhancement factor
      Real, Parameter :: threshmin = 1.0E-19     ! small positive number
      Real, Parameter :: ctolmin = 1.0E-06
!     Real, Parameter :: convfac_298 = 101325.0 * rgas1 * trefm1  ! P/RT at 1 atm and 298 K [ mole/m**3 ]
      Real, Parameter :: difforg = 9.36e-6  ! Diffusivity of organics [m2 s-1]

      Real convfac
      Real tt1, tt2      ! temperature-related factors
      Real tempcorr      ! temperature correction factor for cstar

C Variables used in oligomerization calculations
#ifdef Oligo_aero
      Real expdt         ! non-dimensional loss coefficient
      Real nsvpa         ! particle-phase anthropogenic SVOC [ umolC/m3 ]
      Real nsvpb         ! particle-phase biogenic SVOC [ umolC/m3 ]
#endif

C Variables used in acid-enhanced isoprene SOA calculations
      Real hplus         ! accumulation-mode H+ conc [ nmol/m3 ]
      Real aiso12        ! particle-phase isoprene SVOC [ ug/m3 ]
      Real vviso         ! vapor-phase isoprene SVOC [ ug/m3 ]
      Real xiso3         ! newly produced AISO3J [ ug/m3 ]
      Real isofrac       ! ratio for depletion of vapor-phase products

C Variables used in equilibrium partitioning calculations
      Real drog( n_orgprod ) ! change in precursor conc [ ug/m3 ]
      Real totrog( n_vapor ) ! drog conc mapped to each SVOC [ ug/m3 ]
      Real(8) GR3( n_vapor,n_mode ) ! 3rd Moment Growth for each mode and compound
      Real GR3FRAC( n_vapor,n_mode ) ! Fraction of each mode growing/shrinking
      Real c0    ( n_vapor ) ! cstar at AIRTEMP [ ug/m3 ]
      Real caer0 ( n_vapor ) ! SVOC conc before current time step [ ug/m3 ]
      Real ctoti ( n_vapor ) ! SVOC conc before current time step [ ug/m3 ]
      Real prod  ( n_vapor ) ! SVOC generated during current step [ ug/m3 ]
      Real ctotf ( n_vapor ) ! SVOC conc after current time step [ ug/m3 ]
      Real caer  ( n_vapor ) ! SVOC conc in aerosol phase [ ug/m3 ]
      Real dcaer             ! Change in SVOC conc after partitioning happens[ ug/m3 ]
      Real cbar_org(n_vapor) ! On-line molecular speed of each organic
      Real dv_org            ! On-line gas-phase diffusivity of each organic
      Real totorgsw          ! POA + non-volatile SOA [ umole/m3 ]
      Real totorg            ! SOA + POA before time step [ umole/m3 ]
      Real threshold         ! criterion for establishing gas/part equil.
!     Real faer              ! fraction of total in aerosol, intermediate value
!     Logical check          ! flag to indicate if NEWT subroutine
                             ! converged to a spurious root
      Real totaer, fact      ! total aerosol-phase mass of each semivolaitle component

      real, parameter :: cbar_fac = sqrt( 8.0 * RGASUNIV * oneovpi * 1.e3 )

C Variables for updating 2nd and 3rd moments
      Real(8) m0_init( 2 )  ! initial 0 moment, wet [ mom/m3 ]
      Real(8) m1wet_init( 2 )  ! initial 1st moment, wet [ mom/m3 ]
      Real(8) m2wet_initD( 2 )  ! initial 2nd moment, wet [ mom/m3 ]
      Real m3wet_init( 2 )  ! initial 3rd moment, wet [ mom/m3 ]
      Real m2wet_init( 2 )  ! initial 2nd moment, wet [ mom/m3 ]

C Added for new SOA bisection (hotp 7/6/11)
      Real               :: upb, orgmoles      ! upper bound, total aerosol moles

C-----------------------------------------------------------------------

C Set temperature factors
      tt1 = tref / airtemp
      tt2 = trefm1 - 1.0 / airtemp
      convfac = tt1 * airpres * prefm1

C Set Mass Transfer Properties (although each condensing organic
C will have its own diffusivity and molecular speed, just assume
C as a first guess that they are all uniform and given by a rep-
C resentative compound.
      Dv_org = difforg * ( STDATMPA / AIRPRES ) * ( AIRTEMP / 273.16 ) ** 1.75
      cbar_org( : ) = cbar_fac * SQRT( AIRTEMP * vapor_mwi(:) )

C Compute 3rd moment, 2nd moment. moment_conc arrays are wet
      Call calcmoments ( .true. )
      m3wet_init( : ) = moment3_conc( 1:2 )
      m2wet_init( : ) = moment2_conc( 1:2 )

      GR3( :,: ) = 1.d-30
      m0_init    ( : ) = moment0_conc( 1:2 )
      m2wet_initD( : ) = m2wet_init  ( 1:2 )

      m1wet_init( : ) = m0_init( : ) * aeromode_diam( 1:2 ) * 
     &                      exp( 0.5 * aeromode_lnsg( 1:2 ) ** 2 )

      Do im = 1,2
        Do i = 1, n_vapor
          If ( aero_missing( soa_aeroMap( i ),im ) ) Cycle
          
          Call HCOND3_l( m1wet_init( im ), m2wet_initD( im ), 
     &                   Dv_org, cbar_org(i), GR3( i,im ) )
        End Do
      End Do

      ! Calculate fraction of mass transfer to/from each mode
      gr3frac(:,1) = gr3(:,1) / ( GR3(:,1) + GR3(:,2) )
      gr3frac(:,2) = 1.0 - gr3frac(:,1)

C Initialize drog from ORGPROD and change units to [ ug / m**3 ]
      drog( : ) = orgprod_conc( : ) * drog_ppm2ug( : ) * convfac

C Assignment of drog to totrog. This moving mass from the reactive
C gas tracer species to the array that will be multiplied by alpha
C to give newly formed semivolatile vapors.
      totrog( : ) = 0.0
      Do i = 1, n_vapor
         If ( drog_map( i ) .Gt. 0 ) Then
            totrog( i ) = drog( drog_map( i ) )
         End If
      End Do

C Transfer non-volatile products directly to aerospc_conc array (2nd mode)
      Do i = 1, n_orgprod
         If ( orgprod_aeroMap( i ) .Gt. 0 ) Then
            aerospc_conc( orgprod_aeroMap( i ),2 ) =
     &      aerospc_conc( orgprod_aeroMap( i ),2 )
     &      + orgprod( i )%alphaH * drog( i )
         End If
      End Do

#ifdef Oligo_aero
C Oligomerization

C Convert semi-volatile SOA to non-volatile oligomers by an exponential decay
C process that conserves carbon mass.

      nsvpa = 0.0
      nsvpb = 0.0

      Do i = 1, n_vapor
         If     ( vaporspc( i )%soa_origin .Eq. 'A' ) then
            nsvpa = nsvpa + aerospc_conc( soa_aeroMap( i ),2 )
     &            * vapor_mwi( i ) * vaporspc( i )%soa_decay
         ElseIf ( vaporspc( i )%soa_origin .Eq. 'B' ) then
            nsvpb = nsvpb + aerospc_conc( soa_aeroMap( i ),2 )
     &            * vapor_mwi( i ) * vaporspc( i )%soa_decay
         Endif
      End Do

      expdt = exp( - kolig * dt )

      Do i = 1, n_vapor
         aerospc_conc( soa_aeroMap( i ),2 ) =
     &   aerospc_conc( soa_aeroMap( i ),2 ) * expdt
      End Do

      aerospc_conc( aolga_idx,2 ) = aerospc_conc( aolga_idx,2 )
     &                            + 12.0 * nsvpa * ( 1.0 - expdt ) * olgrat
      aerospc_conc( aolgb_idx,2 ) = aerospc_conc( aolgb_idx,2 )
     &                            + 12.0 * nsvpb * ( 1.0 - expdt ) * olgrat

#endif
C Compute molar concentrations of non-volatile organic material and total
C organic material; in future code updates, molecular weights of non-volatile
C SOA species should be inherited from the include files or namelists to avoid
C hard-coded values.

      ! First Sum Up POA Species if this the NON-Volatile POA Mechanism
      If ( isnvpoa ) Then
         totorgsw = ( aerospc_conc(   apoc_idx,1 ) + aerospc_conc(   apoc_idx,2 )) * aerospc_mwi( apoc_idx )
     &            + ( aerospc_conc( apncom_idx,1 ) + aerospc_conc( apncom_idx,2 )) * aerospc_mwi( apncom_idx )
      Else
         totorgsw = 0.0
      End If

      ! Sum Up Nonvolatile OA. This is just a couple of SOA products.
      ! They are in the orgprd table.
      Do i = 1, n_aerospc
         If ( aerospc( i )%nonVol_soa ) Then
            totorgsw = totorgsw + aerospc_conc( i,2 ) * aerospc_mwi( i )
         End If
      End Do

      ! Sum Up all other OA species. This is most of them. They are in
      ! the vaporspc table.
      totorg = totorgsw
      Do i = 1, n_vapor
         totorg = totorg + ( aerospc_conc( soa_aeroMap( i ),2 ) + 
     &                       aerospc_conc( soa_aeroMap( i ),1 ) ) * vapor_mwi( i )
      End Do

C Initialize ctoti as sum of vapor-phase and particle-phase SVOCs
C Note: all of these CBLK species are in [ug/m**3]
      Do i = 1, n_vapor
         caer0( i ) = aerospc_conc( soa_aeroMap( i ),2 )  
     &              + aerospc_conc( soa_aeroMap( i ),1 )

         ctoti( i ) = vapor_conc( i ) + caer0( i )
      End Do

C Equilibrium Partitioning Calculations

C Initial guess of caer is computed as follows:
C    From eqn (8) of Schell et al. (2001)
C    caer = ctotf - c0 * (caer/MW) / totorg
C    Assuming totorg doesn't change during the timestep,
C    caer * (1 + c0/MW / totorg) = ctotf
C    caer = ctotf / ( 1 + c0/MW / totorg )

      threshold = 0.0  ! This threshold will be compared to the
                       ! saturation concentration to determine whether
                       ! or not OA partitioning is even likely.
      upb = 0.0  ! Upper bound for bisection method

      Do i = 1, n_vapor
         tempcorr   = tt1 * Exp( vaporspc( i )%enth * rgas1 * tt2 )
         c0( i )    = vaporspc( i )%cstar * tempcorr        ! Satn Conc.
         prod( i )  = vaporspc( i )%alpha * totrog( i )     ! Total Vapor Produced
         ctotf( i ) = ctoti( i ) + prod( i )                ! Vapor + Particle
         threshold  = threshold +  ctotf( i ) / c0( i )
         upb        = upb + ctotf( i ) * vapor_mwi( i )      ! upper: all moles (semi and nonvolatile)
      End Do

C Check If gas/particle equilibrium can be established
      If ( ( threshold .Gt. 1.0 ) .Or. ( totorgsw .Gt. threshmin ) ) Then

C        Calculate new SOA by partitioning.
C        Method uses bisection method to solve for total number of moles
C        of orgaic aerosol. Caer is then calculated based on that number.

C        Solve for total number of organic moles in aerosol (hotp 7/5/11)
!        lowb      = totorgsw   ! lower: moles of nonvolatile
!        upb       = lowb + upb ! upper: all moles (semi and nonvolatile)
!        upb       = upb * ( 1.0 + 1e-7*n_vapor ) ! prevent numerical issues
         If ( upb .lt. threshmin ) Then
            ! no SOA to partition
            orgmoles = totorgsw
         Else
            upb       = totorgsw + upb
            orgmoles  = soabisect( totorgsw, upb, totorgsw, c0, ctotf, vapor_mw )
         End If

C        Determine concentration of each semivoltile (hotp 7/5/11)
C        Constrain caer to values between evapmin and ctotf
         Do i = 1, n_vapor
            caer(i) = max( 1.e-7*caer0(i), min( ctotf(i), ctotf(i) * vapor_mw(i) * orgmoles
     &                                                   / ( c0(i) + vapor_mw(i) * orgmoles ) ) )

            ! Check for overshoot in any of the specific modes
            indx = soa_aeroMap( i )

            ! Calculate total change in OA Species Concentration
            dcaer = caer( i ) - caer0( i )  

            if (dcaer < 0.) then

               ! Sum the total aerosol mass of this compound
               totaer = max( aerospc_conc( indx,1 ) + aerospc_conc( indx,2 ), conmin )
               fact   =  1.0 + max(dcaer, -totaer) / totaer

               ! Evaporate using mode-dependent mass fraction
               if (.not. aero_missing(indx,1)) aerospc_conc( indx,1 ) = aerospc_conc( indx,1 ) * fact
               if (.not. aero_missing(indx,2)) aerospc_conc( indx,2 ) = aerospc_conc( indx,2 ) * fact

            else

               ! Condense using mode-dependent condensation flux fraction
               if (.not. aero_missing(indx,1))
     &              aerospc_conc( indx,1 ) = aerospc_conc( indx,1 ) + dcaer * GR3FRAC( i,1 )
               if (.not. aero_missing(indx,2))
     &              aerospc_conc( indx,2 ) = aerospc_conc( indx,2 ) + dcaer * GR3FRAC( i,2 )

            endif

         End Do
 
      Else   ! threshold not exceeded; no material transferred to aerosol phase

         Do i = 1, n_vapor
            indx = soa_aeroMap( i )
            If ( .Not. aero_missing( indx,1 ) ) aerospc_conc( indx,1 ) = conmin   ! array assignment
            If ( .Not. aero_missing( indx,2 ) ) aerospc_conc( indx,2 ) = conmin   ! array assignment
            caer( i ) = sum( aerospc_conc( indx, 1:2 ) )
         End Do

      End If    ! check on equilibrium threshold

C Update Vapor-Phase Organic Species
      vapor_conc( : ) = max( ctotf( : ) - caer( : ), conmin )


C Acid-induced Isoprene SOA if IEPOX uptake is not considered mechanism reactions
      If ( .Not. IEPOX_SOA_mechanism ) Then

C Particle-phase acidity is parameterized by the air concentration of hydrogen
C ion in the particle phase (hplus).  This is computed by electroneutrality,
C and is bounded so that it cannot fall below zero (acid-neutral) nor above the
C range of experimental conditions (530 nmol/m3) tested by Surratt et al.
C (2007).  In future code updates, molecular weights of inorganic ions should
C be inherited from include or namelist files to avoid hard-coded values.

      hplus = aerospc_conc(aso4_idx,2) * aerospc_chomw(aso4_idx)
     &      + aerospc_conc(ano3_idx,2) * aerospc_chomw(ano3_idx)
     &      + aerospc_conc(acl_idx,2)  * aerospc_chomw(acl_idx)
     &      + aerospc_conc(anh4_idx,2) * aerospc_chomw(anh4_idx)
     &      + aerospc_conc(ana_idx,2)  * aerospc_chomw(ana_idx)
     &      + aerospc_conc(amg_idx,2)  * aerospc_chomw(amg_idx)
     &      + aerospc_conc(ak_idx,2)   * aerospc_chomw(ak_idx)
     &      + aerospc_conc(aca_idx,2)  * aerospc_chomw(aca_idx)

      hplus = Min( .530, Max ( 0.0, hplus ) ) ! restrict extrapolation

C Compute amount of acid-enhanced isoprene SOA that is produced during the current
C time step (xiso3) using the normalized expression derived from Surratt et al. (2007)

      aiso12 = aerospc_conc( aiso1_idx, 2 ) + aerospc_conc( aiso2_idx, 2 )
      xiso3  = kacid * hplus * aiso12 - aerospc_conc( aiso3_idx,2 )

C Update CBLK(VISO3J) while conserving the total concentration of low-volatility
C isoprene-derived products.  Since the acid-enhanced product is assumed to be
C non-volatile, only positive values of xiso3 are considered below.  Also, xiso3
C is not allowed to exceed the vapor-phase concentration of isoprene SVOC.  We assume
C the xiso3 mass is depleted from the 2 vapor-phase isoprene products in quantities
C proportional to their equilibrium concentrations.  The factor of 2.7/1.6 accounts
C for the enhancement in the OM/OC ratio when isoprene SOA is formed in the presence
C of SO2, based on the lab experiments of Kleindienst et al. (2007)

      If ( xiso3 .Gt. 0.0 ) Then
         vviso = vapor_conc( iso1_idx ) + vapor_conc( iso2_idx )
         xiso3  = Min( xiso3, vviso )
         aerospc_conc( aiso3_idx,2 ) = aerospc_conc( aiso3_idx,2 ) + xiso3 * ( 2.7 / 1.6 )
         isofrac = vapor_conc( iso1_idx ) / vviso
         vapor_conc( iso1_idx ) = vapor_conc( iso1_idx ) - isofrac * xiso3
         vapor_conc( iso2_idx ) = vapor_conc( iso2_idx ) - ( 1.0 - isofrac ) * xiso3
      End If

      End If ! IEPOX not present for uptake

C Update 3rd moment, 2nd moment, and Dg in CBLK array by assuming that SOA
C condensation/evaporation does not affect the geometric standard deviation.

      ! Get new third moment of all species, do not get second moment
      ! directly from the calcmoments diagnostic. Instead, update it
      ! manually to account for condensation in ORGAER.
      Call calcmoments( .true. )
      moment2_conc( 1 ) = m2wet_init( 1 ) * ((moment3_conc(1)/m3wet_init(1))**2) ** ( 1.0 / 3.0 )
      moment2_conc( 2 ) = m2wet_init( 2 ) * ((moment3_conc(2)/m3wet_init(2))**2) ** ( 1.0 / 3.0 )

      End Subroutine orgaer

C-----------------------------------------------------------------------

      subroutine hcond3_l( am1, am2, dv, cbar, f )

C  calculates the size-dependent term in the condensational-growth rate
C  expression for the 2nd and 3rd moments of a lognormal aerosol mode
C  using the harmonic mean method.  This code follows Section A2 of
C  Binkowski & Shankar (1995).
 
C  Key Subroutines/Functions called:  none
 
C  Revision History:
C     coded November 7, 2003 by Dr. Francis S. Binkowski
C     Revised November 20, 2003 by F. Binkowski to have am1 and
C     am2 as inputs
 
C  Reference:
C   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
C      model 1. Model description and preliminary results, J. Geophys.
C      Res., Vol 100, No D12, 26101-26209, 1995.

      use aeromet_data, only: pi
      implicit none

C *** Arguments:

      real( 8 ), intent( in ) :: am1   ! first moment of mode   [ m/m**3 ]
      real( 8 ), intent( in ) :: am2   ! second moment of mode  [ m**2/m**3 ]
      real,      intent( in ) :: dv    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      real,      intent( in ) :: cbar  ! kinetic velocity of condensing vapor [ m/s ]

      real( 8 ), intent( out ) :: f    ! size-dependent term in condensational-growth
                                       ! rate: F = 3rd moment [ m**3/m**3 s ]
C *** Local Variables:

      real( 8 ) :: gnc3 ! integrals used to calculate F(2) [m^3 / m^3 s]
      real( 8 ) :: gfm3 !

      real( 8 ), parameter :: twopi = 2.0d0 * pi
      real( 8 ), parameter :: pi4 = 0.25d0 * pi

C *** Implement equation A15 of Binkowski & Shankar (1995) for the
C     2nd and 3rd moments of a lognormal mode of arbitrary size.

      gnc3 = twopi * dv * am1   ! 3rd moment, near-continuum
      gfm3 = pi4 * cbar * am2   ! 3rd moment, free-molecular

C *** Implement equation A13 of Binkowski & Shankar (1995) for a
C     lognormal mode of arbitrary size.  These are the size-dependent
C     terms in the condensational-growth rate expression, given in
C     equation 7a of B&S (1995).

      f = max( gnc3 * gfm3 / ( gnc3 + gfm3 ), 1.e-30 ) ! 3rd moment

      end subroutine hcond3_l

C-----------------------------------------------------------------------

      function soabisect( lowerb, upperb, nonvolmol, cstaratt, totsemivol,
     &                    mlwt) result ( nroot )

!     Determines the root of an equation, nroot, that is located
!     between the lowerb and upperb. The equation of interest is defined
!     in the function soaequation. Before calling soabisect, the threshold 
!     for SOA should have already been checked which should ensure a solution.
!     The recommended lower bound is the number of moles of nonvolatile aerosol
!     and the recommended upper bound is the total number of organic moles
!
!     History
!     Created 7/2011 by HOT Pye

      Use utilio_defn

      implicit none

!     Function arguments
      Real, Intent(IN) :: lowerb, upperb      ! lower bound, upper bound 
      Real, Intent(IN) :: nonvolmol           ! nonvolatile aerosol in moles = POA + non-voltile SOA
      Real, Intent(IN) :: cstaratt(n_vapor)   ! Cstars at T of interest
      Real, Intent(IN) :: totsemivol(n_vapor) ! total semivolatile in ug/m3 (gas+aer+newlyformed)
      Real, Intent(IN) :: mlwt(n_vapor)       ! molecular weight of semivolatiles 

!     Result
      Real             :: nroot               ! solution (final nguess)

!     Parameters for solution convergence
      real, parameter  :: eps     = 1.0e-06   ! solution converged if function value is small
      Real, parameter  :: FRACTOL = 1.0e-06   ! Solution converged if abs(old-new)/new < FRACTOL
      Real, parameter  :: ABSTOL  = 1.0e-08   ! Solution converged if there are less than ABSTOL umol/m3 in aerosol

!     Local variables
      Real             :: xn                  ! current guess for total moles organic aerosol
      Real             :: x1, x2, f1, f2, fn  ! intermediate values
      Integer          :: iter                ! counter to prevent infinite loops

      if (upperb < abstol) then
         nroot = upperb
         return
      endif

!     Initialize bounds
      x1  = lowerb
      x2  = upperb

!     Evaluate at lower bound; check if it's a root
      f1 = soaequation( nonvolmol, cstaratt, totsemivol, mlwt, x1 )
      if (abs(f1) < eps) then
         nroot = x1
         return
      endif

!     Evaluate at upper bound; check if it's a root
      f2 = soaequation( nonvolmol, cstaratt, totsemivol, mlwt, x2 )
      if (abs(f2) < eps) then
         nroot = x2
         return
      endif

!     Iterate until a root is found or until abstol (absolute) amount
      do iter = 1, 15

!        Next guess
         xn = (f2 * x1 - f1 * x2) / (f2 - f1)

!        Function evaluated at new guess
         fn = soaequation( nonvolmol, cstaratt, totsemivol, mlwt, xn )

!        Check and end iteration if we aren't moving much
         if (abs(x2 - xn) < fractol * xn .and. abs(fn) < 1.e-4) then
            nroot = xn
            return
         endif

!       Check and end iteration if we are numerically close to a root
        if (abs(fn) < eps) then
           nroot = xn
           return
        endif

!       Update bounds based on new guess
        if ( fn * f2 < 0.0 ) then

           ! ROOT WAS TRAPPED, SO USE REGULA FALSI
           x1 = x2
           f1 = f2
           x2 = xn
           f2 = fn

        else

           ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
           x2 = xn
           f2 = fn
           f1 = 0.5 * f1

          endIf

          ! End iteration if there is only a tiny amount
          if (max(x1,x2) < abstol) then
             nroot = xn
             return
          endif
 
      enddo

      nroot = xn
!     call m3warn( pname, 0, 0, 'Warning: too many iterations in soa bisection')

      end function soabisect

C-----------------------------------------------------------------------

      Function soaequation( nonvmol, cstar, semivol, mlwt, currentN ) Result ( eqnerror )

C     Evaluates the following function of total aerosol moles, N
C     Function evaluates to zero at equilibrium 
C     
C                      totalsemivol_i        mols nonvolatile POA + SOA
C     f(N) = sum_i (  ----------------   ) + ---------------------------  - 1.0 
C                     cstar_i + mw_i*N                  N
C
C     History
C     Created 7/2011 HOT Pye

      Implicit None

C     Function inputs and output
      Real, INTENT(IN) :: nonvmol, currentN                ! nonvolatile moles, total moles
      Real, INTENT(IN) :: mlwt(n_vapor)                    ! molecular weights of semivolatiles
      Real, INTENT(IN) :: semivol(n_vapor), cstar(n_vapor) ! semivolatile mass, sat conc in ug/m3
      Real             :: eqnerror                         ! result, deviation from 0.0 in equation

C     Local variables
      Integer          :: i                                ! counter

C     Compute function value for current value of N 
      eqnerror = nonvmol / currentN - 1.0
      Do i = 1, n_vapor
         eqnerror = eqnerror + semivol(i) / ( cstar(i) + mlwt(i) * currentN )
      End Do

      End function soaequation

C-----------------------------------------------------------------------

      End Module soa_defn
