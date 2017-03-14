module edgar_42_emis
  implicit none
  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !  Based on code from Chris Nolte, US EPA July 2013
  !  Heavily modified from code by Jia Xing, US EPA, originally
  !  prepared by Wang Litao, Chen Dan and Zhang Qiang
  !  DESE, Tsinghua University
  !  All rights reserved
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  integer, parameter :: nsec    = 17 ! number of edgar emissions sectors
  integer, parameter :: nvars3d = 48
  integer, parameter :: nvoc    = 19 ! number of speciation factors for nmvoc
  integer, parameter :: nevars  =  8 ! number of raw emissions species
  integer, parameter :: npms    = 19 ! number of speciation factors for pm2.5

  real :: molwt(nvars3d)

  ! Gas chemistry species

  integer :: ICO   = 1
  integer :: INO   = 2
  integer :: INO2  = 3
  integer :: INH3  = 4
  integer :: ISO2  = 5
  integer :: ISULF = 6

  ! NMVOC species

  integer :: IALD2    =  7
  integer :: IALDX    =  8
  integer :: IETH     =  9
  integer :: IETHA    = 10
  integer :: IETOH    = 11
  integer :: IFORM    = 12
  integer :: IIOLE    = 13
  integer :: IISOP    = 14
  integer :: IMEOH    = 15    
  integer :: INVOL    = 16
  integer :: IOLE     = 17
  integer :: IPAR     = 18
  integer :: ITERP    = 19
  integer :: ITOL     = 20
  integer :: IUNK     = 21
  integer :: IUNR     = 22
  integer :: IXYL     = 23
  integer :: IBENZENE = 24
  integer :: ISESQ    = 25

! Other gas species (currently not set here)

  integer :: ICL2  = 26
  integer :: IHCL  = 27
  integer :: IHONO = 28

! PM species
  integer :: IPMC    = 29
  integer :: IPEC    = 30
  integer :: IPMFINE = 31
  integer :: IPNO3   = 32
  integer :: IPOC    = 33
  integer :: IPSO4   = 34
  integer :: IPCL    = 35
  integer :: IPNH4   = 36
  integer :: IPNA    = 37
  integer :: IPMG    = 38
  integer :: IPK     = 39
  integer :: IPCA    = 40
  integer :: IPNCOM  = 41
  integer :: IPFE    = 42
  integer :: IPAL    = 43
  integer :: IPSI    = 44
  integer :: IPTI    = 45
  integer :: IPMN    = 46
  integer :: IPH2O   = 47
  integer :: IPMOTHR = 48

  type overlap_vars
     integer              :: ncells = 0
     integer, allocatable :: i(:)
     integer, allocatable :: j(:)
     real,    allocatable :: area(:)
  end type overlap_vars

  type(overlap_vars), allocatable :: emis_w(:)

  integer, parameter :: nx_e42 = 3600
  integer, parameter :: ny_e42 = 1800

!!  integer, parameter :: n_ch4 = 17
!!  integer, parameter :: ch4_sects(n_ch4) = &
!!                 (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 /)
!!
  integer, parameter :: n_co = 11
  integer, parameter :: co_sects(n_co) = &
                                        (/ 1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 16 /)

  integer, parameter :: n_nh3 = 11
  integer, parameter :: nh3_sects(n_nh3) = &
                                       (/ 1, 2, 4, 5, 6, 7, 9, 11, 12, 13, 14 /)

  integer, parameter :: n_nmvoc = 12
  integer, parameter :: nmvoc_sects(n_nmvoc) = &
                                   (/ 1, 2, 3, 4, 5, 6, 7, 10, 13, 14, 15, 16 /)

  integer, parameter :: n_nox = 13
  integer, parameter :: nox_sects(n_nox) = &
                               (/ 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16 /)

  integer, parameter :: n_pm2_5 = 13
  integer, parameter :: pm2_5_sects(n_pm2_5) = &
                               (/ 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16 /)

  integer, parameter :: n_pm10 = 13
  integer, parameter :: pm10_sects(n_pm10) = &
                               (/ 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16 /)

  integer, parameter :: n_so2 = 11
  integer, parameter :: so2_sects(n_so2) = &
                                        (/ 1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 16 /)

  integer, save :: year_stored = -1

  type edgar_vars_and_sects
!!   real :: ch4  (n_ch4)
     real :: co   (n_co)
     real :: nh3  (n_nh3)
     real :: nmvoc(n_nmvoc)
     real :: nox  (n_nox)
     real :: pm2_5(n_pm2_5)
     real :: pm10 (n_pm10)
     real :: so2  (n_so2)
  end type edgar_vars_and_sects

  type(edgar_vars_and_sects), allocatable :: edgar42_vars(:)

  real, allocatable :: edgar42_emis(:,:,:)

  type vert_interp_type
     integer :: nlevs
     real    :: facts(150)
  end type vert_interp_type

  type(vert_interp_type), allocatable :: vinterp(:,:)

  real :: emiscnvt(nvars3d)

  real :: rawmn(nsec, 12)
  real :: rawwk(nsec,  7)
  real :: rawhr(nsec, 24)

  real :: vocspec(nsec, nvoc)
  real :: pmspec (nsec, npms)


contains


  subroutine edgar_42_init(nvars3d_emis, vname3d_emis, units3d_emis)
    use mem_grid,  only: mwa, mza
    use geia_emis, only: geia_init
    implicit none
    
    integer,                    intent(inout) :: nvars3d_emis
    character(24), allocatable, intent(inout) :: vname3d_emis(:)
    character(24), allocatable, intent(inout) :: units3d_emis(:)

    real, parameter :: gpkg = 1.0e3

    integer :: i
    
    allocate(emis_w          (mwa))
    allocate(edgar42_vars    (mwa))
    allocate(edgar42_emis(mza,mwa,nvars3d))
    allocate(vinterp(nsec,mza))

    nvars3d_emis = nvars3d
    allocate(vname3d_emis(nvars3d_emis))
    allocate(units3d_emis(nvars3d_emis))

    call emis_overlap(nx_e42, ny_e42)
    call interp_to_olam()
    call comp_vert_facts()
    call geia_init()

    i = ico
    vname3d_emis(i) = 'CO'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  28.0

    i = ino
    vname3d_emis(i) = 'NO'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  14.0  ! emis units are mass N, so use 14 to convert to moles rather than 30
      
    i = ino2
    vname3d_emis(i) = 'NO2'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  14.0  ! emis units are mass N, so use 14 to convert to moles rather than 46

    i = inh3
    vname3d_emis(i) = 'NH3'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  17.0

    i = iso2
    vname3d_emis(i) = 'SO2'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  32.0  ! emis units are mass S, so use 32 to convert to moles rather than 64

    i = isulf
    vname3d_emis(i) = 'SULF'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  32.0  ! emis units are mass S, so use 32 to convert to moles rather than 64

    i = iald2
    vname3d_emis(I) = 'ALD2'
    units3d_emis(I) = 'MOLES/S'
    molwt       (i) =  44.0

    i = ialdx
    vname3d_emis(i) = 'ALDX'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  44.0

    i = ieth 
    vname3d_emis(i) = 'ETH'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  28.0

    i = ietha
    vname3d_emis(i) = 'ETHA'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  30.0

    i = ietoh
    vname3d_emis(i) = 'ETOH'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  46.0

    i = iform
    vname3d_emis(i) = 'FORM'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  30.0

    i = iiole 
    vname3d_emis(i) = 'IOLE'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  48.0

    i = iisop
    vname3d_emis(i) = 'ISOP'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  68.0

    i = imeoh
    vname3d_emis(i) = 'MEOH'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  32.0

    i = invol
    vname3d_emis(i) = 'NVOL'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  16.0

    i = iole 
    vname3d_emis(i) = 'OLE'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  27.0

    i = ipar 
    vname3d_emis(i) = 'PAR'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  44.0

    i = iterp
    vname3d_emis(i) = 'TERP'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  136.0

    i = itol 
    vname3d_emis(i) = 'TOL'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  92.0

    i = iunk
    vname3d_emis(i) = 'UNK'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  352.0

    i = iunr
    vname3d_emis(i) = 'UNR'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  16.0

    i = ixyl 
    vname3d_emis(i) = 'XYL'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  106.0

    i = ibenzene
    vname3d_emis(i) = 'BENZENE'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  78.0

    i = isesq
    vname3d_emis(i) = 'SESQ'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  1.0   ! CHECK: why is this 1 in Jia's file?  are BVOC emis in gC?

    i = icl2
    vname3d_emis(i) = 'CL2'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  71.0

    i = ihcl
    vname3d_emis(i) = 'HCL'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  36.5

    i = ihono
    vname3d_emis(i) = 'HONO'
    units3d_emis(i) = 'MOLES/S'
    molwt       (i) =  47.0

    i = ipmc
    vname3d_emis(i) = 'PMC'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipec
    vname3d_emis(i) = 'PEC'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipmfine
    vname3d_emis(i) = 'PMFINE'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipno3
    vname3d_emis(i) = 'PNO3'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipoc
    vname3d_emis(i) = 'POC'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipso4
    vname3d_emis(i) = 'PSO4'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipcl
    vname3d_emis(i) = 'PCL'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipnh4
    vname3d_emis(i) = 'PNH4'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipna
    vname3d_emis(i) = 'PNA'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipmg
    vname3d_emis(i) = 'PMG'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipk
    vname3d_emis(i) = 'PK'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipca
    vname3d_emis(i) = 'PCA'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipncom
    vname3d_emis(i) = 'PNCOM'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipfe
    vname3d_emis(i) = 'PFE'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipal
    vname3d_emis(i) = 'PAL'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipsi
    vname3d_emis(i) = 'PSI'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipti
    vname3d_emis(i) = 'PTI'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipmn
    vname3d_emis(i) = 'PMN'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = iph2o
    vname3d_emis(i) = 'PH2O'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    i = ipmothr
    vname3d_emis(i) = 'PMOTHR'
    units3d_emis(i) = 'G/S'
    molwt       (i) =  1.0

    do i = 1, nvars3d
       emiscnvt(i) = gpkg / molwt(i)
    enddo

    ! Monthly distribution factors

    rawmn( 1,:) = (/ 1.20, 1.15, 1.05, 1.00, 0.90, 0.85, 0.80, 0.88, 0.95, 1.00, 1.08, 1.15 /)
    rawmn( 2,:) = (/ 1.10, 1.08, 1.05, 1.00, 0.95, 0.90, 0.93, 0.95, 0.97, 1.00, 1.02, 1.05 /)
    rawmn( 3,:) = (/ 0.88, 0.92, 0.98, 1.03, 1.05, 1.06, 1.01, 1.02, 1.06, 1.05, 1.01, 0.93 /)
    rawmn( 4,:) = (/ 0.88, 0.92, 0.98, 1.03, 1.05, 1.06, 1.01, 1.02, 1.06, 1.05, 1.01, 0.93 /)
    rawmn( 5,:) = (/ 1.70, 1.50, 1.30, 1.00, 0.70, 0.40, 0.20, 0.40, 0.70, 1.05, 1.40, 1.65 /)
    rawmn( 6,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn( 7,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn( 8,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn( 9,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn(10,:) = (/ 0.95, 0.96, 1.02, 1.00, 1.01, 1.03, 1.03, 1.01, 1.04, 1.03, 1.01, 0.91 /)
    rawmn(11,:) = (/ 0.45, 1.30, 2.35, 1.70, 0.85, 0.85, 0.85, 1.00, 1.10, 0.65, 0.45, 0.45 /)
    rawmn(12,:) = (/ 0.45, 1.30, 2.35, 1.70, 0.85, 0.85, 0.85, 1.00, 1.10, 0.65, 0.45, 0.45 /)
    rawmn(13,:) = (/ 0.45, 1.30, 2.35, 1.70, 0.85, 0.85, 0.85, 1.00, 1.10, 0.65, 0.45, 0.45 /)
    rawmn(14,:) = (/ 0.45, 1.30, 2.35, 1.70, 0.85, 0.85, 0.85, 1.00, 1.10, 0.65, 0.45, 0.45 /)
    rawmn(15,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn(16,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawmn(17,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)

    ! Weekly distribution factors per sector

    rawwk( 1,:) = (/ 0.85, 1.06, 1.06, 1.06, 1.06, 1.06, 0.85 /)
    rawwk( 2,:) = (/ 0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80 /)
    rawwk( 3,:) = (/ 0.79, 1.02, 1.06, 1.08, 1.10, 1.14, 0.81 /)
    rawwk( 4,:) = (/ 0.79, 1.02, 1.06, 1.08, 1.10, 1.14, 0.81 /)
    rawwk( 5,:) = (/ 0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80 /)
    rawwk( 6,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk( 7,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk( 8,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk( 9,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(10,:) = (/ 0.50, 1.20, 1.20, 1.20, 1.20, 1.20, 0.50 /)
    rawwk(11,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(12,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(13,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(14,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(15,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(16,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)
    rawwk(17,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /)

    ! Hourly distribution factors per sector

    rawhr( 1,:) = (/ 0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92, 1.08, 1.19, 1.22, &
                     1.21, 1.21, 1.17, 1.15, 1.14, 1.13, 1.10, 1.07, 1.04, 1.02, &
                     1.02, 1.01, 0.96, 0.88 /)
    rawhr( 2,:) = (/ 0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02, 1.09, 1.16, 1.22, &
                     1.28, 1.30, 1.22, 1.24, 1.25, 1.16, 1.08, 1.01, 0.95, 0.90, &
                     0.85, 0.81, 0.78, 0.75 /)
    rawhr( 3,:) = (/ 0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, &
                     1.24, 1.20, 1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, &
                     0.74, 0.62, 0.61, 0.44 /)
    rawhr( 4,:) = (/ 0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, &
                     1.24, 1.20, 1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, &
                     0.74, 0.62, 0.61, 0.44 /)
    rawhr( 5,:) = (/ 0.40, 0.40, 0.40, 0.40, 0.40, 0.50, 1.20, 1.50, 1.60, 1.60, &
                     1.40, 1.20, 1.10, 1.10, 1.00, 1.00, 1.00, 1.10, 1.40, 1.50, &
                     1.40, 1.40, 1.00, 0.40 /)
    rawhr( 6,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr( 7,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr( 8,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr( 9,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr(10,:) = (/ 0.50, 0.35, 0.20, 0.10, 0.10, 0.20, 0.75, 1.25, 1.40, 1.50, &
                     1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.40, 1.25, 1.10, &
                     1.00, 0.90, 0.80, 0.70 /)
    rawhr(11,:) = (/ 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, &
                     1.45, 1.60, 1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, &
                     0.65, 0.60, 0.60, 0.60 /)
    rawhr(12,:) = (/ 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, &
                     1.45, 1.60, 1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, &
                     0.65, 0.60, 0.60, 0.60 /)
    rawhr(13,:) = (/ 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, &
                     1.45, 1.60, 1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, &
                     0.65, 0.60, 0.60, 0.60 /)
    rawhr(14,:) = (/ 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, &
                     1.45, 1.60, 1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, &
                     0.65, 0.60, 0.60, 0.60 /)
    rawhr(15,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr(16,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)
    rawhr(17,:) = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
                     1.00, 1.00, 1.00, 1.00 /)

    vocspec( 1,:) = (/ 0.24724e-5, 0.83117e-6, 0.58162e-3, 0.61433e-2, 0.0       ,  &
                       0.55483e-2, 0.16416e-4, 0.0       , 0.0       , 0.0       ,  &
                       0.15930e-2, 0.30714e-1, 0.0       , 0.12013e-2, 0.0       ,  &
                       0.98334e-2, 0.10760e-2, 0.56493e-3, 0.0                     /)
 
    vocspec( 2,:) = (/ 0.0       , 0.0       , 0.12065e-3, 0.44364e-3, 0.0       ,  &
                       0.20907e-1, 0.89703e-5, 0.0       , 0.0       , 0.0       ,  &
                       0.21292e-3, 0.18985e-1, 0.0       , 0.28711e-3, 0.0       ,  &
                       0.36421e-2, 0.18305e-3, 0.40610e-3, 0.0                     /)
 
    vocspec( 3,:) = (/ 0.37758e-3, 0.18399e-3, 0.20582e-2, 0.72071e-3, 0.11894e-4,  &
                       0.74435e-3, 0.32335e-3, 0.25169e-4, 0.17082e-4, 0.17634e-4,  &
                       0.13015e-2, 0.34423e-1, 0.10655e-4, 0.15034e-2, 0.23314e-6,  &
                       0.47702e-2, 0.12890e-2, 0.50683e-3, 0.17049e-5              /)
 
    vocspec( 4,:) = (/ 0.37758e-3, 0.18399e-3, 0.20582e-2, 0.72071e-3, 0.11894e-4,  &
                       0.74435e-3, 0.32335e-3, 0.25169e-4, 0.17082e-4, 0.17634e-4,  &
                       0.13015e-2, 0.34423e-1, 0.10655e-4, 0.15034e-2, 0.23314e-6,  &
                       0.47702e-2, 0.12890e-2, 0.50683e-3, 0.17049e-5              /)
 
    vocspec( 5,:) = (/ 0.32275e-2, 0.30807e-2, 0.0       , 0.20430e-4, 0.0       ,  &
                       0.62666e-2, 0.47299e-6, 0.0       , 0.0       , 0.11461e-4,  &
                       0.12019e-2, 0.92665e-2, 0.58834e-5, 0.15193e-2, 0.0       ,  &
                       0.46819e-2, 0.43657e-3, 0.16691e-5, 0.94134e-6              /)
 
    vocspec( 6,:) = (/ 0.24724e-5, 0.83117e-6, 0.58162e-3, 0.61433e-2, 0.0       ,  &
                       0.55483e-2, 0.16416e-4, 0.0       , 0.0       , 0.0       ,  &
                       0.15930e-2, 0.30714e-1, 0.0       , 0.12013e-2, 0.0       ,  &
                       0.98334e-2, 0.10760e-2, 0.56493e-3, 0.0                     /)
 
    vocspec( 7,:) = (/ 0.16257e-3, 0.24866e-3, 0.52832e-3, 0.88840e-2, 0.94085e-3,  &
                       0.58609e-2, 0.52268e-4, 0.25074e-4, 0.18672e-3, 0.10465e-3,  &
                       0.10168e-2, 0.27446e-1, 0.29961e-4, 0.41018e-3, 0.0       ,  &
                       0.13330e-1, 0.20880e-3, 0.28897e-3, 0.47938e-5              /)
 
    vocspec( 8,:) = (/ 0.43741e-3, 0.62256e-3, 0.12869e-2, 0.53024e-3, 0.34969e-3,  &
                       0.60290e-3, 0.13877e-3, 0.67461e-4, 0.50238e-3, 0.28157e-3,  &
                       0.19828e-2, 0.29057e-1, 0.58069e-4, 0.73368e-3, 0.0       ,  &
                       0.82535e-2, 0.50463e-3, 0.41431e-3, 0.92910e-5              /)
 
    vocspec( 9,:) = (/ 0.24539e-3, 0.39357e-3, 0.40283e-2, 0.33265e-3, 0.19618e-3,  &
                       0.33864e-3, 0.79623e-4, 0.37847e-4, 0.28184e-3, 0.15796e-3,  &
                       0.20358e-2, 0.31656e-1, 0.54062e-4, 0.71939e-3, 0.0       ,  &
                       0.84159e-2, 0.31972e-3, 0.64305e-3, 0.86499e-5              /)
 
    vocspec(10,:) = (/ 0.19595e-4, 0.12467e-3, 0.57465e-4, 0.68298e-5, 0.11234e-2,  &
                       0.10167e-3, 0.30200e-4, 0.22703e-7, 0.92545e-3, 0.88522e-4,  &
                       0.34678e-3, 0.40462e-1, 0.16402e-3, 0.12774e-2, 0.0       ,  &
                       0.58793e-2, 0.79457e-3, 0.27556e-4, 0.26243e-4              /)
 
    vocspec(11,:) = (/ 0.27735e-2, 0.18805e-2, 0.29772e-2, 0.10280e-3, 0.67794e-4,  &
                       0.62859e-2, 0.23283e-3, 0.90771e-4, 0.97394e-4, 0.54586e-4,  &
                       0.20293e-2, 0.18743e-1, 0.13215e-4, 0.42119e-3, 0.0       ,  &
                       0.37797e-2, 0.17795e-3, 0.27230e-3, 0.21143e-5              /)
 
    vocspec(12,:) = (/ 0.53115e-3, 0.51103e-3, 0.25429e-2, 0.22169e-2, 0.18121e-3,  &
                       0.49473e-2, 0.14082e-3, 0.18312e-4, 0.13677e-3, 0.49815e-4,  &
                       0.14797e-2, 0.28392e-1, 0.44577e-4, 0.79310e-3, 0.31085e-7,  &
                       0.65904e-2, 0.54441e-3, 0.30820e-3, 0.71323e-5              /)
 
    vocspec(13,:) = (/ 0.0       , 0.45794e-3, 0.85461e-2, 0.43683e-2, 0.0       ,  &
                       0.29857e-5, 0.11602e-3, 0.0       , 0.0       , 0.0       ,  &
                       0.21746e-2, 0.28112e-1, 0.15971e-3, 0.39872e-3, 0.0       ,  &
                       0.81625e-2, 0.30230e-3, 0.0       , 0.25553e-4              /)
 
    vocspec(14,:) = (/ 0.0       , 0.45794e-3, 0.85461e-2, 0.43683e-2, 0.0       ,  &
                       0.29857e-5, 0.11602e-3, 0.0       , 0.0       , 0.0       ,  &
                       0.21746e-2, 0.28112e-1, 0.15971e-3, 0.39872e-3, 0.0       ,  &
                       0.81625e-2, 0.30230e-3, 0.0       , 0.25553e-4              /)
 
    vocspec(15,:) = (/ 0.34125e-3, 0.29209e-4, 0.66517e-2, 0.24878e-4, 0.16407e-4,  &
                       0.34901e-3, 0.64866e-3, 0.31652e-5, 0.23571e-4, 0.13211e-4,  &
                       0.30184e-2, 0.44776e-1, 0.27239e-5, 0.34411e-4, 0.0       ,  &
                       0.16987e-2, 0.23228e-4, 0.19431e-4, 0.43583e-6              /)
 
    vocspec(16,:) = (/ 0.0       , 0.0       , 0.12065e-3, 0.44364e-3, 0.0       ,  &
                       0.20907e-1, 0.89703e-5, 0.0       , 0.0       , 0.0       ,  &
                       0.21292e-3, 0.18985e-1, 0.0       , 0.28711e-3, 0.0       ,  &
                       0.36421e-2, 0.18305e-3, 0.40610e-3, 0.0                     /)

    pmspec( 1,:) = (/ 0.5093335E-01, 0.7896125E+00, 0.1130028E-02, 0.4630086E-01, &
                      0.1120233E+00, 0.7211590E-03, 0.3337995E-02, 0.7741310E-04, &
                      0.3904520E-03, 0.4337831E-02, 0.3660728E-01, 0.1848782E-01, &
                      0.2732115E-01, 0.5511979E-01, 0.8349822E-01, 0.4166080E-02, &
                      0.2657730E-03, 0.0000000E+00, 0.5552815E+00  /)
 
    pmspec( 2,:) = (/ 0.9305152E-01, 0.6753209E+00, 0.2679254E-02, 0.8549677E-01, &
                      0.1434515E+00, 0.7028740E-03, 0.2242595E-02, 0.1085700E-06, &
                      0.5534310E-06, 0.3046828E-02, 0.2230610E-01, 0.3394840E-01, &
                      0.1897411E-01, 0.3847033E-01, 0.5918732E-01, 0.2844717E-02, &
                      0.1858590E-03, 0.0000000E+00, 0.4934112E+00  /)
 
    pmspec( 3,:) = (/ 0.6101691E+00, 0.9693286E-01, 0.1155748E-02, 0.2226151E+00, &
                      0.3277256E-02, 0.3359480E-03, 0.1200506E-02, 0.1450480E-03, &
                      0.3722520E-04, 0.4685540E-04, 0.7059460E-03, 0.5582348E-01, &
                      0.4950410E-03, 0.1100320E-03, 0.4725560E-03, 0.1009660E-04, &
                      0.3457080E-05, 0.0000000E+00, 0.3754668E-01 /)
 
    pmspec( 4,:) = (/ 0.6101691E+00, 0.9693286E-01, 0.1155748E-02, 0.2226151E+00, &
                      0.3277256E-02, 0.3359480E-03, 0.1200506E-02, 0.1450480E-03, &
                      0.3722520E-04, 0.4685540E-04, 0.7059460E-03, 0.5582348E-01, &
                      0.4950410E-03, 0.1100320E-03, 0.4725560E-03, 0.1009660E-04, &
                      0.3457080E-05, 0.0000000E+00, 0.3754668E-01 /)
 
    pmspec( 5,:) = (/ 0.5864368E-01, 0.4091476E+00, 0.1990023E-02, 0.5227183E+00, &
                      0.7500431E-02, 0.3031797E-02, 0.1569836E-02, 0.9900380E-03, &
                      0.1294020E-03, 0.9462254E-02, 0.2145300E-03, 0.3636599E+00, &
                      0.1619750E-03, 0.1594740E-03, 0.5064650E-03, 0.6129750E-05, &
                      0.1067850E-05, 0.0000000E+00, 0.2925463E-01 /)
 
    pmspec( 6,:) = (/ 0.5093335E-01, 0.7896125E+00, 0.1130028E-02, 0.4630086E-01, &
                      0.1120233E+00, 0.7211590E-03, 0.3337995E-02, 0.7741310E-04, &
                      0.3904520E-03, 0.4337831E-02, 0.3660728E-01, 0.1848782E-01, &
                      0.2732115E-01, 0.5511979E-01, 0.8349822E-01, 0.4166080E-02, &
                      0.2657730E-03, 0.0000000E+00, 0.5552815E+00 /)
 
    pmspec( 7,:) = (/ 0.4459063E-01, 0.7177606E+00, 0.4337183E-02, 0.1437265E+00, &
                      0.8958507E-01, 0.7313706E-02, 0.6027760E-04, 0.5988345E-02, &
                      0.1362770E-03, 0.6925087E-02, 0.1695260E-02, 0.5749354E-01, &
                      0.3648100E-03, 0.5915320E-03, 0.5089387E-02, 0.2078024E-02, &
                      0.1569250E-03, 0.2083330E-01, 0.6090342E+00 /)
 
    pmspec( 8,:) = (/ 0.5228960E-02, 0.7616066E+00, 0.2224370E-02, 0.4563876E-01, &
                      0.1853013E+00, 0.2118021E+00, 0.0000000E+00, 0.7470303E-01, &
                      0.9194220E-03, 0.1227662E+00, 0.5442225E-02, 0.1824227E-01, &
                      0.7411399E-01, 0.5171750E-03, 0.7125520E-02, 0.0000000E+00, &
                      0.1626227E-02, 0.1055746E-02, 0.2432926E+00 /)
 
    pmspec( 9,:) = (/ 0.1831900E-01, 0.8551678E+00, 0.3498630E-02, 0.9169705E-01, &
                      0.3131752E-01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.3669882E-01, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.7504197E-02, 0.8109648E+00 /)
 
    pmspec(10,:) = (/ 0.1736608E-01, 0.7923931E+00, 0.4582310E-03, 0.1743871E+00, &
                      0.1539545E-01, 0.1765447E-02, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.2125743E-02, 0.2731039E-01, 0.6961630E-01, &
                      0.3602950E-03, 0.2666186E-01, 0.2745451E-01, 0.1837506E+00, &
                      0.0000000E+00, 0.3688202E-02, 0.4496598E+00 /)
 
    pmspec(11,:) = (/ 0.1063594E+00, 0.4947548E+00, 0.3428792E-02, 0.3793308E+00, &
                      0.1612624E-01, 0.8833633E-01, 0.1757649E-01, 0.6425010E-02, &
                      0.8050290E-03, 0.6918893E-01, 0.9166050E-03, 0.2656971E+00, &
                      0.1497949E-02, 0.2462452E-02, 0.6068226E-02, 0.1359610E-03, &
                      0.3155090E-04, 0.9610900E-05, 0.3560352E-01 /)
 
    pmspec(12,:) = (/ 0.1260744E+00, 0.6392629E+00, 0.2052486E-02, 0.1552949E+00, &
                      0.6853534E-01, 0.2105129E-01, 0.2184586E-02, 0.5903437E-02, &
                      0.1897730E-03, 0.1502208E-01, 0.1032118E-01, 0.7554849E-01, &
                      0.1133864E-01, 0.1451952E-01, 0.2217069E-01, 0.1333417E-01, &
                      0.1817300E-03, 0.4846071E-02, 0.4426512E+00 /)
 
    pmspec(13,:) = (/ 0.4410000E-01, 0.8114600E+00, 0.1640000E-02, 0.8770000E-01, &
                      0.5510000E-01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.3510000E-01, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.1320000E-01, 0.7631600E+00 /)
 
    pmspec(14,:) = (/ 0.4410000E-01, 0.8114600E+00, 0.1640000E-02, 0.8770000E-01, &
                      0.5510000E-01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.3510000E-01, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.1320000E-01, 0.7631600E+00 /)
 
    pmspec(15,:) = (/ 0.4410000E-01, 0.8114600E+00, 0.1640000E-02, 0.8770000E-01, &
                      0.5510000E-01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.3510000E-01, &
                      0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                      0.0000000E+00, 0.1320000E-01, 0.7631600E+00 /)
 
    pmspec(16,:) = (/ 0.9305152E-01, 0.6753209E+00, 0.2679254E-02, 0.8549677E-01, &
                      0.1434515E+00, 0.7028740E-03, 0.2242595E-02, 0.1085700E-06, &
                      0.5534310E-06, 0.3046828E-02, 0.2230610E-01, 0.3394840E-01, &
                      0.1897411E-01, 0.3847033E-01, 0.5918732E-01, 0.2844717E-02, &
                      0.1858590E-03, 0.0000000E+00, 0.4934112E+00 /)

  end subroutine edgar_42_init


  subroutine process_edgar_emis()
    use misc_coms,  only: current_time, io6
    use mem_ijtabs, only: jtab_w, jtw_prog, itab_w
    use mem_grid,   only: glonw, lsw, lpw
    use mem_turb,   only: frac_sfc
    use geia_emis,  only: cl_emis, hcl_emis

    implicit none

    integer :: year, month, hour, day
    integer :: itz, dow
    integer :: jw, iw, j, n, ns, k, ks, ka, v, i

    real :: timefac(nsec)
    real :: fact1, fact2, fact3, emisn

    real, parameter :: so2_to_sulf(nsec) = (/ 0.000229, 0.000202, 0.0, 0.0, &
         0.000143, 0.000229, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000202, 0.0 /)

    real, parameter :: gpkg = 1.0e3

    integer, parameter ::  cl_sec = 13
    integer, parameter :: hcl_sec =  1

    integer, external :: day_of_week

    year = current_time%year
    if (year < 1990) year = 1990
    if (year > 2008) year = 2008

    if (year /= year_stored)  call interp_to_olam()

    month = current_time%month
    dow = day_of_week(current_time%month, current_time%date, current_time%year)
    
    do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)

       edgar42_emis(:,iw,:) = 0.0
       
       day  = dow
       itz  = nint( glonw(iw) / 15.0 )
       hour = int( current_time%time / 3600. ) + itz + 1

       if (hour <  1) then
          day  = day - 1
          hour = hour + 24
       elseif (hour > 24) then
          day  = day + 1
          hour = hour - 24
       endif
       
       if (day < 1) day = day + 7
       if (day > 7) day = day - 7

       do n = 1, nsec
          timefac(n) = rawmn(n,month) * rawwk(n,day) * rawhr(n,hour)
       enddo
       
!       if (jw == 1) write(io6,*) "Processing CO..."

       do j = 1, n_co
          n = co_sects(j)

         if (edgar42_vars(iw)%co(j) > 1.e-20) then

             fact2 = emiscnvt(ico) * timefac(n) * edgar42_vars(iw)%co(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   edgar42_emis(k,iw,ico) = edgar42_emis(k,iw,ico) &
                                          + vinterp(n,ka)%facts(ks) * fact3
                enddo
             enddo
         endif
       enddo

!       if (jw == 1) write(io6,*) "Processing NOx..."

       do j = 1, n_nox
          n = nox_sects(j)

          if (edgar42_vars(iw)%nox(j) > 1.e-20) then

             fact2 = emiscnvt(ino) * timefac(n) * edgar42_vars(iw)%nox(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   emisn = vinterp(n,ka)%facts(ks) * fact3
                   edgar42_emis(k,iw,ino ) = edgar42_emis(k,iw,ino ) + 0.9 * emisn
                   edgar42_emis(k,iw,ino2) = edgar42_emis(k,iw,ino2) + 0.1 * emisn
                enddo

             enddo
          endif
       enddo

!       if (jw == 1) write(io6,*) "Processing NH3..."

       do j = 1, n_nh3
          n = nh3_sects(j)

          if (edgar42_vars(iw)%nh3(j) > 1.e-20) then

             fact2 = emiscnvt(inh3) * timefac(n) * edgar42_vars(iw)%nh3(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   edgar42_emis(k,iw,inh3) = edgar42_emis(k,iw,inh3) &
                                           + vinterp(n,ka)%facts(ks) * fact3
                enddo

             enddo
          endif
       enddo

!      if (jw == 1) write(io6,*) "Processing SO2..."

       do j = 1, n_so2
          n = so2_sects(j)

          if (edgar42_vars(iw)%so2(j) > 1.e-20) then

             fact2 = emiscnvt(iso2) * timefac(n) * edgar42_vars(iw)%so2(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   emisn = vinterp(n,ka)%facts(ks) * fact3
                   edgar42_emis(k,iw,iso2 ) = edgar42_emis(k,iw,iso2 ) + emisn
                   edgar42_emis(k,iw,isulf) = edgar42_emis(k,iw,isulf) + emisn * so2_to_sulf(n)
                enddo

             enddo
          endif
       enddo

!      if (jw == 1) write(io6,*) "Processing NMVOC..."

       do j = 1, n_nmvoc
          n = nmvoc_sects(j)

          if (edgar42_vars(iw)%nmvoc(j) > 1.e-20) then

             fact1 = timefac(n) * edgar42_vars(iw)%nmvoc(j) * gpkg

             do v = iald2, isesq
                i = v - iald2 + 1

                if (vocspec(n,i) > 1.e-12) then
                   
                   fact2 = fact1 * vocspec(n,i)

                   do ns = 1, lsw(iw)
                      ka = lpw(iw) + ns - 1

                      fact3 = fact2 * frac_sfc(ns,iw)

                      do ks = 1, vinterp(n,ka)%nlevs
                         k  = ks + ka - 1
                         edgar42_emis(k,iw,v) = edgar42_emis(k,iw,v) &
                                              + vinterp(n,ka)%facts(ks) * fact3
                      enddo
                   enddo
                endif
             enddo
          endif
       enddo

!      if (jw == 1) write(io6,*) "Processing PM10..."

       do j = 1, n_pm10
          n = pm10_sects(j)

          if (edgar42_vars(iw)%pm10(j) > 1.e-20) then

             fact2 = gpkg * timefac(n) * edgar42_vars(iw)%pm10(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   emisn = vinterp(n,ka)%facts(ks) * fact3
                   edgar42_emis(k,iw,ipmc) = edgar42_emis(k,iw,ipmc) + emisn
                enddo

             enddo
          endif
       enddo

!      if (jw == 1) write(io6,*) "Processing PM2.5..."

       do j = 1, n_pm2_5
          n = pm2_5_sects(j)

          if (edgar42_vars(iw)%pm2_5(j) > 1.e-20) then

             fact2 = gpkg * timefac(n) * edgar42_vars(iw)%pm2_5(j)

             do ns = 1, lsw(iw)
                ka = lpw(iw) + ns - 1

                fact3 = fact2 * frac_sfc(ns,iw)

                do ks = 1, vinterp(n,ka)%nlevs
                   k  = ks + ka - 1
                   emisn = vinterp(n,ka)%facts(ks) * fact3
                   edgar42_emis(k,iw,ipmc) = max( edgar42_emis(k,iw,ipmc) - emisn, 0.0)
                   do v = ipec, ipmothr
                      i = v - ipec + 1
                      edgar42_emis(k,iw,v) = edgar42_emis(k,iw,v) + emisn * pmspec(n,i)
                   enddo
                enddo

             enddo
          endif
       enddo

!       if (jw == 1) write(io6,*) "Processing HCL from geia..."

       n = hcl_sec
       if (hcl_emis(iw) > 1.e-20) then

          fact2 = timefac(n) * hcl_emis(iw)

          do ns = 1, lsw(iw)
             ka = lpw(iw) + ns - 1

             fact3 = fact2 * frac_sfc(ns,iw)

             do ks = 1, vinterp(n,ka)%nlevs
                k  = ks + ka - 1
                edgar42_emis(k,iw,ihcl) = edgar42_emis(k,iw,ihcl) &
                                        + vinterp(n,ka)%facts(ks) * fact3
             enddo

          enddo
       endif

!       if (jw == 1) write(io6,*) "Processing particulate CL from geia..."

       n = cl_sec
       if (cl_emis(iw) > 1.e-20) then

          fact2 = timefac(n) * cl_emis(iw)

          do ns = 1, lsw(iw)
             ka = lpw(iw) + ns - 1

             fact3 = fact2 * frac_sfc(ns,iw)

             do ks = 1, vinterp(n,ka)%nlevs
                k  = ks + ka - 1
                edgar42_emis(k,iw,ipcl) = edgar42_emis(k,iw,ipcl) &
                                        + vinterp(n,ka)%facts(ks) * fact3
             enddo

          enddo
       endif

    enddo

  end subroutine process_edgar_emis



  subroutine interp_to_olam()
    use misc_coms,  only: current_time, io6, iparallel
    use cgrid_spcs, only: n_gc_spc, gc_spc
    use mem_ijtabs, only: jtab_w, jtw_prog, itab_w
    use mem_grid,   only: arw0, mwa
    use mem_para,   only: myrank
    use oname_coms, only: nl
    use hdf5_utils

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer :: year, n, j, jw, iw, ier
    real    :: rawdata(nx_e42,ny_e42)
    integer :: ndims, idims(2)
    logical :: exists

    character(200) :: filename
    character( 4 ) :: yyear

    character(50), parameter :: secname(nsec) = (/                             &
                             '1.ENERGY_INDUSTRY_AND_WASTE_INCINERATION      ', &
                             '2.COMBUSTION_IN_MANUFACTURING_INDUSTRY        ', &
                             '3.NON-ROAD_TRANSPORTATION                     ', &
                             '4.ROAD_TRANSPORTATION                         ', &
                             '5.RESIDENTIAL                                 ', &
                             '6.TRANSFORMATION_OIL_PRODUCTION_AND_REFINERIES', &
                             '7.NON-METALLIC_PAPER_CHEMICAL_INDUSTRY        ', &
                             '8.METAL_PROCESSES                             ', &
                             '9.CHEMICAL_INDUSTRY                           ', &
                             '10.SOLVENTS                                   ', &
                             '11.MANURE_MANAGEMENT                          ', &
                             '12.AGRICULTURAL_SOILS                         ', &
                             '13.AGRICULTURAL_WASTE_BURNING                 ', &
                             '14.LARGE_SCALE_BIOMASS_BURNING                ', &
                             '15.SOLID_WASTE_DISPOSAL                       ', &
                             '16.FOSSIL_FUEL_FIRES                          ', &
                             '17.WASTE_WATER                                ' /)

!!    character(40), parameter :: ch4suf(n_ch4) = (/               &
!!                                '_IPCC_1A1_1A2.0.1x0.1.nc4    ', &
!!                                '_IPCC_1B2b.0.1x0.1.nc4       ', &
!!                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4 ', &
!!                                '_IPCC_1A3b.0.1x0.1.nc4       ', &
!!                                '_IPCC_1A4.0.1x0.1.nc4        ', &
!!                                '_IPCC_1B2b.0.1x0.1.nc4       ', &
!!                                '_IPCC_2.0.1x0.1.nc4          ', &
!!                                '_IPCC_1B1.0.1x0.1.nc4        ', &
!!                                '_IPCC_1B2b.0.1x0.1.nc4       ', &
!!                                '_IPCC_4A.0.1x0.1.nc4         ', &
!!                                '_IPCC_4B.0.1x0.1.nc4         ', &
!!                                '_IPCC_4C_4D.0.1x0.1.nc4      ', &
!!                                '_IPCC_4F.0.1x0.1.nc4         ', &
!!                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4', &
!!                                '_IPCC_6A_6C.0.1x0.1.nc4      ', &
!!                                '_IPCC_7A.0.1x0.1.nc4         ', &
!!                                '_IPCC_6B.0.1x0.1.nc4         ' /)

    character(40), parameter :: cosuf(n_co) = (/                   &
                                '_IPCC_1A1a_6.0.1x0.1.nc4       ', &
                                '_IPCC_1A2.0.1x0.1.nc4          ', &
                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4   ', &
                                '_IPCC_1A3b.0.1x0.1.nc4         ', &
                                '_IPCC_1A4.0.1x0.1.nc4          ', &
                                '_IPCC_1B2a_c_1A1b_c.0.1x0.1.nc4', &
                                '_IPCC_2A_2B_2D.0.1x0.1.nc4     ', &
                                '_IPCC_2C.0.1x0.1.nc4           ', &
                                '_IPCC_4F.0.1x0.1.nc4           ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4  ', &
                                '_IPCC_7A.0.1x0.1.nc4           ' /)

    character(40), parameter :: nh3suf(n_nh3) = (/               &
                                '_IPCC_1A1a_6.0.1x0.1.nc4     ', &
                                '_IPCC_1A2.0.1x0.1.nc4        ', &
                                '_IPCC_1A3.0.1x0.1.nc4        ', &
                                '_IPCC_1A4.0.1x0.1.nc4        ', &
                                '_IPCC_1A1b_c.0.1x0.1.nc4     ', &
                                '_IPCC_2A.0.1x0.1.nc4         ', &
                                '_IPCC_2B.0.1x0.1.nc4         ', &
                                '_IPCC_4B.0.1x0.1.nc4         ', &
                                '_IPCC_4C_4D.0.1x0.1.nc4      ', &
                                '_IPCC_4F.0.1x0.1.nc4         ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4' /)

    character(40), parameter :: nmvocsuf(n_nmvoc) = (/                 &
                                '_IPCC_1A1a.0.1x0.1.nc4             ', &
                                '_IPCC_1A2.0.1x0.1.nc4              ', &
                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4       ', &
                                '_IPCC_1A3b.0.1x0.1.nc4             ', &
                                '_IPCC_1A4.0.1x0.1.nc4              ', &
                                '_IPCC_1A1b_c_1B_2C1_2C2.0.1x0.1.nc4', &
                                '_IPCC_2A_B_D_E_F_G.0.1x0.1.nc4     ', &
                                '_IPCC_3.0.1x0.1.nc4                ', &
                                '_IPCC_4F.0.1x0.1.nc4               ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4      ', &
                                '_IPCC_6A_6C.0.1x0.1.nc4            ', &
                                '_IPCC_7A.0.1x0.1.nc4               ' /)

    character(40), parameter :: noxsuf(n_nox) = (/                 &
                                '_IPCC_1A1a.0.1x0.1.nc4         ', &
                                '_IPCC_1A2.0.1x0.1.nc4          ', &
                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4   ', &
                                '_IPCC_1A3b.0.1x0.1.nc4         ', &
                                '_IPCC_1A4.0.1x0.1.nc4          ', &
                                '_IPCC_1B2a_c_1A1b_c.0.1x0.1.nc4', &
                                '_IPCC_2.0.1x0.1.nc4            ', &
                                '_IPCC_4B.0.1x0.1.nc4           ', &
                                '_IPCC_4C_4D.0.1x0.1.nc4        ', &
                                '_IPCC_4F.0.1x0.1.nc4           ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4  ', &
                                '_IPCC_6A_6C.0.1x0.1.nc4        ', &
                                '_IPCC_7A.0.1x0.1.nc4           ' /)

    character(40), parameter :: pm2_5suf(n_pm2_5) = (/ &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4', &
                                '.nc4' /)

    character(40), parameter :: pm10suf(n_pm10) = (/               &
                                '_IPCC_1A1a.0.1x0.1.nc4         ', &
                                '_IPCC_1A2.0.1x0.1.nc4          ', &
                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4   ', &
                                '_IPCC_1A3b.0.1x0.1.nc4         ', &
                                '_IPCC_1A4.0.1x0.1.nc4          ', &
                                '_IPCC_1B2a_c_1A1b_c.0.1x0.1.nc4', &
                                '_IPCC_2_3.0.1x0.1.nc4          ', &
                                '_IPCC_4B.0.1x0.1.nc4           ', &
                                '_IPCC_4C_4D.0.1x0.1.nc4        ', &
                                '_IPCC_4F.0.1x0.1.nc4           ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4  ', &
                                '_IPCC_6A_6C.0.1x0.1.nc4        ', &
                                '_IPCC_7A.0.1x0.1.nc4           ' /)

    character(40), parameter :: so2suf(n_so2) = (/                  &
                                '_IPCC_1A1a_6C.0.1x0.1.nc4       ', &
                                '_IPCC_1A2.0.1x0.1.nc4           ', &
                                '_IPCC_1A3a_c_d_e.0.1x0.1.nc4    ', &
                                '_IPCC_1A3b.0.1x0.1.nc4          ', &
                                '_IPCC_1A4.0.1x0.1.nc4           ', &
                                '_IPCC_1B1_1B2_1A1b_c.0.1x0.1.nc4', &
                                '_IPCC_2B_2D.0.1x0.1.nc4         ', &
                                '_IPCC_2C.0.1x0.1.nc4            ', &
                                '_IPCC_4F.0.1x0.1.nc4            ', &
                                '_IPCC_5A_C_D_F_4E.0.1x0.1.nc4   ', &
                                '_IPCC_7A.0.1x0.1.nc4            ' /)

    year = current_time%year
    if (year < 1990) year = 1990
    if (year > 2008) year = 2008

    if (year == year_stored) return

    do iw = 1, mwa
!!     edgar42_vars(iw)%ch4  (:) = 0.0
       edgar42_vars(iw)%co   (:) = 0.0
       edgar42_vars(iw)%nh3  (:) = 0.0
       edgar42_vars(iw)%nmvoc(:) = 0.0
       edgar42_vars(iw)%nox  (:) = 0.0
       edgar42_vars(iw)%pm2_5(:) = 0.0
       edgar42_vars(iw)%pm10 (:) = 0.0
       edgar42_vars(iw)%so2  (:) = 0.0
    enddo

    year_stored = year
    write(yyear,'(I4)') year
    
!!    ! Average CH4 to olam grid
!!
!!    do j = 1, n_ch4
!!       n = ch4_sects(j)
!!
!!       if (myrank == 0) then
!!
!!          rawdata(:,:) = 0.0
!!
!!          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/ch4/v42_CH4_" // &
!!                     yyear // trim(ch4suf(j))
!!
!!          write(*,'(A)') "reading: ", trim(filename)
!!
!!          inquire(file=filename, exist=exists)
!!
!!          if (.not. exists) then
!!
!!             write(*,'(A)') "Cannot find emissions file ", trim(filename)
!!
!!          else
!!
!!             call shdf5_open(filename, 'R')
!!
!!             ndims = 0
!!             idims = 0
!!             call shdf5_info('emi_ch4', ndims, idims)
!!             
!!             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
!!                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " CH4"
!!             else
!!                call shdf5_irec(ndims, idims, 'emi_ch4', rvar2=rawdata)
!!             endif
!!
!!             call shdf5_close()
!!          endif
!!       endif
!!
!!#ifdef OLAM_MPI
!!       if (iparallel == 1) then
!!          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
!!       endif
!!#endif
!!
!!       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
!!          
!!          ! sum the emission overlaps that contribute to this cell and
!!          ! multiply by area to get emissions in kg/sec
!!
!!          do n = 1, emis_w(iw)%ncells
!!             edgar42_vars(iw)%ch4(j) = edgar42_vars(iw)%ch4(j) + &
!!                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
!!          enddo
!!
!!       enddo
!!    enddo

    ! Average CO to olam grid

    do j = 1, n_co
       n = co_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/co/v42_CO_" // &
                     yyear // trim(cosuf(j))

          write(*,'(A)') "reading: ", trim(filename)

          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_co', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " CO"
             else
                call shdf5_irec(ndims, idims, 'emi_co', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%co(j) = edgar42_vars(iw)%co(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo
       enddo

    enddo

    ! Average NH3 to olam grid

    do j = 1, n_nh3
       n = nh3_sects(j)

       if (myrank == 0) then
          
          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/nh3/v42_NH3_" // &
                     yyear // trim(nh3suf(j))

          write(*,'(A)') "reading: ", trim(filename)

          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_nh3', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " NH3"
             else
                call shdf5_irec(ndims, idims, 'emi_nh3', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%nh3(j) = edgar42_vars(iw)%nh3(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo

       enddo
    enddo

    ! Average NMVOC to olam grid

    do j = 1, n_nmvoc
       n = nmvoc_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/nmvoc/v42_NMVOC_" // &
                     yyear // trim(nmvocsuf(j))

          write(*,'(A)') "reading: ", trim(filename)

          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_nmvoc', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " NMVOC"
             else
                call shdf5_irec(ndims, idims, 'emi_nmvoc', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%nmvoc(j) = edgar42_vars(iw)%nmvoc(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo

       enddo
    enddo

    ! Average NOX to olam grid

    do j = 1, n_nox
       n = nox_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/nox/v42_NOx_" // &
                     yyear // trim(noxsuf(j))

          write(*,'(A)') "reading: ", trim(filename)

          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_nox', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " NOX"
             else
                call shdf5_irec(ndims, idims, 'emi_nox', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%nox(j) = edgar42_vars(iw)%nox(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo
       enddo

    enddo

    ! Average PM2.5 to olam grid

    do j = 1, n_pm2_5
       n = pm2_5_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/pm2.5/v42_PM25_" // &
                     yyear // trim(pm2_5suf(j))

          write(io6,'(A)') "reading: ", trim(filename)
     
          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_pm2.5', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " PM2.5"
             else
                call shdf5_irec(ndims, idims, 'emi_pm2.5', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%pm2_5(j) = edgar42_vars(iw)%pm2_5(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo

       enddo
    enddo

    ! Average PM10 to olam grid

    do j = 1, n_pm10
       n = pm10_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/pm10/v42_PM10_" // &
                     yyear // trim(pm10suf(j))

          write(io6,'(A)') "reading: ", trim(filename)
     
          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_pm10', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " PM10"
             else
                call shdf5_irec(ndims, idims, 'emi_pm10', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif
       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%pm10(j) = edgar42_vars(iw)%pm10(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo
       enddo

    enddo

    ! Average SO2 to olam grid

    do j = 1, n_so2
       n = so2_sects(j)

       if (myrank == 0) then

          rawdata(:,:) = 0.0

          filename = trim(nl%emis_dir) // "/" // trim(secname(n)) // "/so2/v42_SO2_" // &
                     yyear // trim(so2suf(j))

          write(*,'(A)') "reading: ", trim(filename)

          inquire(file=filename, exist=exists)

          if (.not. exists) then

             write(*,'(A)') "Cannot find emissions file ", trim(filename)

          else

             call shdf5_open(filename, 'R')

             ndims = 0
             idims = 0
             call shdf5_info('emi_so2', ndims, idims)
             
             if ( ndims /= 2 .and. idims(1) /= nx_e42 .and. idims(2) /= ny_e42 ) then
                write(*,*) "Cannot find emissions field ", trim(SECNAME(n)), " SO2"
             else
                call shdf5_irec(ndims, idims, 'emi_so2', rvar2=rawdata)
             endif

             call shdf5_close()
          endif
       endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Bcast( rawdata, nx_e42*ny_e42, MPI_REAL, 0, MPI_COMM_WORLD, ier )
       endif
#endif

       do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)
          
          do n = 1, emis_w(iw)%ncells
             edgar42_vars(iw)%so2(j) = edgar42_vars(iw)%so2(j) + &
                  rawdata(emis_w(iw)%i(n),emis_w(iw)%j(n)) * emis_w(iw)%area(n)
          enddo

       enddo
    enddo

  end subroutine interp_to_olam


  subroutine emis_overlap(nlon, nlat)
    use consts_coms, only: erad, pio180, piu180, r8
    use misc_coms,   only: io6
    use mem_ijtabs,  only: jtab_w, itab_w, jtw_prog
    use mem_grid,    only: glatw, glonw, glatm, glonm, arw0, mwa, &
                           xem, yem, zem, xew, yew, zew

    implicit none

    integer, intent(in) :: nlon, nlat
  
    integer, parameter :: maxvert = 7

    integer :: np, i, j, n, ng, nn, ngrp, jw, iw, im
    integer :: js, je, is(2), ie(2), ijsize
    real    :: max180, min180
    real    :: dx, dy, area, x, y, sumarea, areaij
    real    :: lats(4), lons(4)

    real(r8):: xg(4), yg(4), alpha(maxvert)
    real(r8):: xf(maxvert), yf(maxvert)
    real    :: flats(maxvert), flons(maxvert)

    integer :: jtrap
    real :: xtrap (4,30+4+30*4)  ! trapezoid x coordinates
    real :: ytrap (4,30+4+30*4)  ! trapezoid y coordinates
    real :: traparea(30+4+30*4)  ! trapezoid area

    real :: tolerance
    real, parameter :: lat0 = -90.0
    real, parameter :: lon0 =   0.0

    real :: dxe, dye, dze
    real :: coswlon, sinwlon
    real :: coswlat, sinwlat

    integer, allocatable :: ipoints(:), itmp(:)
    integer, allocatable :: jpoints(:), jtmp(:)
    real,    allocatable :: areafrc(:), atmp(:)

    allocate(ipoints(1000))
    allocate(jpoints(1000))
    allocate(areafrc(1000))

    dx = 360.0 / real(nlon)
    dy = 180.0 / real(nlat)

    do jw = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(jw)

       np = itab_w(iw)%npoly
       
       if (mod(jw,1000) == 0) &
            write(io6,*) "Hex cells:", jw, jtab_w(jtw_prog)%jend(1), glatw(iw), glonw(iw)

       ! Skip cells near poles - small emissions and many lat-lon cells

       if (glatw(iw) > -70.0 .and. glatw(iw) < 78.0 ) then

          tolerance = 1.e-4 * arw0(iw)
          sumarea   = 0.0

          do n = 1, np
             im = itab_w(iw)%im(n)
             flats(n) = glatm(im)
             flons(n) = glonm(im)
          enddo

          ! Do any cells straddle 0 degrees longitude?

          if ( any( flons(1:np) < 0.0 .and. flons(1:np) > -30.0) .and. &
               any( flons(1:np) > 0.0 .and. flons(1:np) <  30.0) ) then

             ngrp = 2
             max180 = maxval(flons(1:np))
             min180 = minval(flons(1:np)) + 360.0

          else
             ngrp = 1
          endif

          ! Convert -180 <-> 180 to 0 <-> 360

          do n = 1, np
             if (flons(n) < 0.0) flons(n) = flons(n) + 360.0
          enddo

          ! Calculate which lat-lon emissions cells may overlap with this iw location

          js = 1 + int( (minval(flats(1:np))+90.) / dy)
          je = 1 + int( (maxval(flats(1:np))+90.) / dy)

          if (ngrp == 1) then
             is(1) = 1 + int( minval(flons(1:np)) / dx)
             ie(1) = 1 + int( maxval(flons(1:np)) / dx)
          else
             is(1) = 1
             ie(1) = 1 + int( max180 / dx )

             is(2) = 1 + int( min180 / dx )
             ie(2) = nlon
          endif

          ! Calculate the iw cell vertices on a polar-stereographic tangent plane

          sinwlat = sin(glatw(iw) * pio180)
          coswlat = cos(glatw(iw) * pio180)
          sinwlon = sin(glonw(iw) * pio180)
          coswlon = cos(glonw(iw) * pio180)

          do n = 1, np
             im = itab_w(iw)%im(n)

             dxe = xem(im) - xew(iw)
             dye = yem(im) - yew(iw)
             dze = zem(im) - zew(iw)

             call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x, y)

             xf(n) = x
             yf(n) = y
          enddo

          ng = 0

          ! loop over all emissions cells that may overlap with this olam cell

          do nn = 1, ngrp
             do j = js, je
                do i = is(nn), ie(nn)

                   lons = (/ real(i-1)*dx, real(i-1)*dx,  real(i)*dx, real(i  )*dx /)
                   lats = (/ real(j-1)*dy, real(j  )*dy,  real(j)*dy, real(j-1)*dy /) + lat0
        
                   ! x,y coordinates of emissions cell on tangent plane

                   do n = 1, 4
                      call ll_xy2(lats(n), lons(n), coswlat, sinwlat, coswlon, sinwlon, &
                                  xew(iw), yew(iw), zew(iw), x, y)
                      xg(n) = x
                      yg(n) = y
                   enddo
                 
                   areaij =  pio180 * erad**2 * abs(sin(lats(2)*pio180)-sin(lats(1)*pio180)) * dx

                   alpha(:) = 0.0_r8

                   ! check if emis cell is entirely within the olam cell

                   if (areaij < arw0(iw)) then
                      do n = 1, 4
                         call inout_check(np, xf, yf, xg(n), yg(n), alpha(n))
                         if (alpha(n) < 1.0_r8) exit
                      enddo
                   endif

                   if (all(alpha(1:4) > 1.0_r8)) then

                      ! emis cell is entirely within this olam cell
                      area = areaij

                   else

                      alpha(:) = 0.0_r8

                      ! check if olam cell is entirely within the emis cell

                      if (arw0(iw) < areaij) then
                         do n = 1, np
                            call inout_check(4, xg, yg, xf(n), yf(n), alpha(n))
                            if (alpha(n) < 1.0_r8) exit
                         enddo
                      endif

                      if (all(alpha(1:np) > 1.0_r8)) then

                         ! olam cell is entirely within this emissions cell
                         area = arw0(iw)

                      else

                         ! compute any overlap
                         call polygon_overlap(iw, np, 4, xf, yf, xg, yg, area, &
                              jtrap, xtrap, ytrap, traparea)

                      endif

                   endif
                 
                   if (area > tolerance) then
                      ng = ng + 1
                      sumarea = sumarea + area

                      ijsize = size(ipoints)
                      if (ng > ijsize) then
                         allocate(itmp(ijsize+1000))
                         allocate(jtmp(ijsize+1000))
                         allocate(atmp(ijsize+1000))
                         itmp(1:ijsize) = ipoints(1:ijsize)
                         jtmp(1:ijsize) = jpoints(1:ijsize)
                         atmp(1:ijsize) = areafrc(1:ijsize)
                         call move_alloc(itmp, ipoints)
                         call move_alloc(jtmp, jpoints)
                         call move_alloc(atmp, areafrc)
                      endif
                   
                      ipoints(ng) = i
                      jpoints(ng) = j
                      areafrc(ng) = area
                   endif

                enddo
             enddo
          enddo

       else

          ! We are close to the N/S pole; don't do any emissions
          ng      = 0
          area    = 0.0
          sumarea = 0.0

       endif

       if (ng > 0) then
          emis_w(iw)%ncells = ng

          allocate(emis_w(iw)%i(ng))
          allocate(emis_w(iw)%j(ng))
          allocate(emis_w(iw)%area(ng))
     
          emis_w(iw)%i(1:ng) = ipoints(1:ng)
          emis_w(iw)%j(1:ng) = jpoints(1:ng)
          emis_w(iw)%area(1:ng) = areafrc(1:ng)
       endif

    enddo

  end subroutine emis_overlap


  subroutine comp_vert_facts()
    use mem_grid, only: zm, lpw, mza
    use misc_coms, only: io6
    implicit none

    real,    parameter :: ztop  = 20.e3
    integer, parameter :: nlays = 20

    real, parameter :: zlays(nlays) = (/  19.79,  42.79,  69.03,  99.32, 134.52, &
                                         174.67, 221.44, 274.93, 337.70, 409.93, &
                                         493.47, 590.29, 703.32, 834.82, 988.14, &
                                        1166.97,1377.21,1624.59,1916.68,2264.44 /)
    real :: rawlay(nsec, nlays)
    real :: rawsum(nsec, nlays)
    real :: rawhgt      (nlays)

    real :: zsum(mza)
    real :: zfac(mza)

    integer :: k, n, j, nlvs
    real :: zbot
    
    ! Vertical distribution factors per sector

    rawlay( 1,:) = (/ 0.01, 0.01, 0.01, 0.02, 0.04, 0.04, 0.04, 0.09, 0.09, 0.14, &
                      0.14, 0.11, 0.11, 0.11, 0.02, 0.02, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 2,:) = (/ 0.01, 0.01, 0.01, 0.02, 0.04, 0.04, 0.04, 0.09, 0.09, 0.14, &
                      0.14, 0.11, 0.11, 0.11, 0.02, 0.02, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 3,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 4,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 5,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 6,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 7,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 8,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay( 9,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(10,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(11,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(12,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(13,:) = (/ 0.02, 0.02, 0.03, 0.03, 0.07, 0.08, 0.13, 0.13, 0.14, 0.18, &
                      0.17, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(14,:) = (/ 0.13, 0.12, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, &
                      0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.00 /)

    rawlay(15,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(16,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawlay(17,:) = (/ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)

    rawsum(:,1) = rawlay(:,1)
    do k = 2, 20
       rawsum(:,k) = rawsum(:,k-1) + rawlay(:,k)
    enddo

    do k = 2, mza-1

       zbot = zm(k-1)
       rawhgt(:) = zbot + zlays(:) * (ztop - zbot) / ztop

       do n = 1, nsec
          
          vinterp(n,k)%facts(:) = 0.0
          
          if (rawlay(n,1) > 0.999 .or. zm(k) > 0.5*ztop) then

             ! Put all emissions in bottom layer layer
             vinterp(n,k)%nlevs    = 1
             vinterp(n,k)%facts(1) = 1.0

          else

             ! Distribute emissions in the vertical
             zsum(k-1) = 0.0
             call hintrp_cc(nlays, rawsum(n,:), rawhgt, mza-k, zsum(k), zm(k))
             zfac(1:mza-k) = zsum(k:mza-1) - zsum(k-1:mza-2)

             nlvs = 1
             do j = mza-k, 1, -1
                if (zfac(j) > 1.e-4) then
                   nlvs = min(j,150)
                   exit
                endif
             enddo
             vinterp(n,k)%nlevs         = nlvs
             vinterp(n,k)%facts(1:nlvs) = zfac(1:nlvs)

          endif

       enddo
    enddo

    vinterp(:,1) = vinterp(:,2)

  end subroutine comp_vert_facts

end module edgar_42_emis
