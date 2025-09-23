Module mem_sfcg

  use max_dims,     only: maxnlspoly, maxgrds, maxngrdll, pathlen, maxremote
  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          nest_ud_vars, nest_wd_vars
  implicit none

  private :: maxnlspoly, maxgrds, maxngrdll, pathlen, maxremote
  private :: itab_md_vars, itab_ud_vars, itab_wd_vars
  private :: nest_ud_vars, nest_wd_vars

  character(pathlen) :: sfcgfile

! HEX-GRID INFORMATION FOR LAND/SEA CELLS (ALL GLOBAL INDICES)

  integer :: nmsfc, mmsfc
  integer :: nvsfc, mvsfc
  integer :: nwsfc, mwsfc

  Type itab_msfc_vars
     integer :: imglobe
     integer :: irank   ! rank of parallel process at this MSFC pt
     integer :: ivn(3)  ! array of V neighbors of this M pt
     integer :: iwn(3)  ! array of W neighbors of this M pt
     integer :: imn(3)  ! array of M neighbors of this M pt
  End type itab_msfc_vars

  Type itab_vsfc_vars
     integer :: ivglobe
     integer :: irank   ! rank of parallel process at this VSFC pt
     integer :: imn(2)
     integer :: iwn(2)

     real :: cosv(2)    ! cosine of angle between V and zonal dir (Voronoi)
     real :: sinv(2)    ! sine of angle between V and zonal dir (Voronoi)

     real :: dxps(2)    ! xps (eastward) displacement from neighbor W pts
     real :: dyps(2)    ! yps (northward) displacement from neighbor W pts
  End type itab_vsfc_vars

  Type itab_wsfc_vars
     integer :: iwglobe  ! global sfcg index of this WSFC pt (in parallel run)
     integer :: irank    ! rank of parallel process at this WSFC pt

     integer :: npoly

     integer :: imn(7)
     integer :: ivn(7)
     integer :: iwn(7)

     real :: dirv(7)     ! pos direction of V neighbors

     real :: farm(7)     ! Fraction of arw0 in each M point sector
     real :: farv(7)     ! Fraction of arw0 in each V point sector

     real :: gxps1(7)    ! gradient weight xe component for point 1
     real :: gyps1(7)    ! gradient weight ye component for point 1

     real :: gxps2(7)    ! gradient weight xe component for point 2
     real :: gyps2(7)    ! gradient weight ye component for point 2

     real :: ecvec_vx(7) ! factors converting V to earth cart. velocity
     real :: ecvec_vy(7) ! factors converting V to earth cart. velocity
     real :: ecvec_vz(7) ! factors converting V to earth cart. velocity

     ! The following ITAB_WSFC members are written to the GRIDFILE.  All other
     ! members are written to the SFCGRIDFILE.  Dimension of (8) is estimate
     ! of max possible number of atm cells that could overlap this sfc cell.

     integer :: nwatm          ! number of atm cells coupled to this sfc cell
     integer :: iwatm    (8)   ! local-rank atm index of atm cell coupled to this sfc cell
     integer :: kwatm    (8)   ! vertical index of atm cell coupled to this sfc cell
     real    :: arc      (8)   ! coupling area between this sfc cell and 1 atm cell [m^2]
     real    :: arcoarsfc(8)   ! arc over (divided by) area of sfc cell
     real    :: arcoariw (8)   ! arc over area of atm cell
     real    :: arcoarkw (8)   ! arc over area of atm cell that contacts sfc at kw level
  End type itab_wsfc_vars

  type (itab_msfc_vars), allocatable, target :: itab_msfc(:)
  type (itab_vsfc_vars), allocatable, target :: itab_vsfc(:)
  type (itab_wsfc_vars), allocatable, target :: itab_wsfc(:)

  Type itab_msfc_pd_vars
     integer :: ivn(3)
     integer :: iwn(3)
     integer :: imn(3)
  End type itab_msfc_pd_vars

  Type itab_vsfc_pd_vars
     integer :: imn(2)
     integer :: iwn(2)
  End type itab_vsfc_pd_vars

  Type itab_wsfc_pd_vars
     integer :: npoly = 0

     integer :: imn(7)
     integer :: ivn(7)
     integer :: iwn(7)

     integer :: nwatm    = 0
     integer :: iwatm(8) = 0  ! dimension is estimate of max possible

     integer :: leaf_class = -1
  End type itab_wsfc_pd_vars

  type (itab_msfc_pd_vars), allocatable, target :: itab_msfc_pd(:)
  type (itab_vsfc_pd_vars), allocatable, target :: itab_vsfc_pd(:)
  type (itab_wsfc_pd_vars), allocatable, target :: itab_wsfc_pd(:)

  Type itabg_msfc_vars            ! Global data structure for MSFC pts
     integer :: imsfc_myrank = -1 ! local (parallel subdomain) index of a MSFC pt
     integer :: irank = -1        ! rank of parallel process at a MSFC pt
  End Type itabg_msfc_vars

  Type itabg_vsfc_vars            ! Global data structure for VSFC pts
     integer :: ivsfc_myrank = -1 ! local (parallel subdomain) index of a VSFC pt
     integer :: irank = -1        ! rank of parallel process at a VSFC pt
  End Type itabg_vsfc_vars

  Type itabg_wsfc_vars            ! Global data structure for WSFC pts
     integer :: iwsfc_myrank = -1 ! local (parallel subdomain) WSFC index
     integer :: irank = -1        ! rank of parallel process at a WSFC pt
  End Type itabg_wsfc_vars

  type (itabg_msfc_vars), allocatable, target :: itabg_msfc(:)
  type (itabg_vsfc_vars), allocatable, target :: itabg_vsfc(:)
  type (itabg_wsfc_vars), allocatable, target :: itabg_wsfc(:)

  Type jtab_msfc_vars
     integer, allocatable :: imsfc(:)
     integer              :: jend
  End Type jtab_msfc_vars

  Type jtab_vsfc_vars
     integer, allocatable :: ivsfc(:)
     integer              :: jend
  End Type jtab_vsfc_vars

  Type jtab_wsfc_vars
     integer, allocatable :: iwsfc(:)
     integer              :: jend
  End Type jtab_wsfc_vars

  type (jtab_msfc_vars) :: jtab_msfc_swm
  type (jtab_vsfc_vars) :: jtab_vsfc_swm
  type (jtab_wsfc_vars) :: jtab_wsfc_swm

  Type surface_grid_vars

     ! Surface grid geometry and location

     real, allocatable :: xem  (:) ! earth x coord of sfc M points [m]
     real, allocatable :: yem  (:) ! earth y coord of sfc M points [m]
     real, allocatable :: zem  (:) ! earth z coord of sfc M points [m]
     real, allocatable :: glatm(:) ! latitude of sfc cell M points [deg]
     real, allocatable :: glonm(:) ! longitude of sfc cell M points [deg]
!    real, allocatable :: topm (:) ! topographic height of sfc W points
     real, allocatable :: arm0 (:) ! Area of IM triangle at earth surface [m^2]

     real, allocatable :: xev (:) ! earth x coord of sfc V points [m]
     real, allocatable :: yev (:) ! earth y coord of sfc V points [m]
     real, allocatable :: zev (:) ! earth z coord of sfc V points [m]
     real, allocatable :: dnu (:) ! grid cell distance across U face [m]
     real, allocatable :: dniu(:) ! inverse of dnu [1/m]
     real, allocatable :: dnv (:) ! grid cell distance across V face [m]
     real, allocatable :: dniv(:) ! inverse of dnv [1/m]
     real, allocatable :: unx  (:) ! U face normal unit vector x component [m]
     real, allocatable :: uny  (:) ! U face normal unit vector y component [m]
     real, allocatable :: unz  (:) ! U face normal unit vector z component [m]
     real, allocatable :: vnx  (:) ! V face normal unit vector x component [m]
     real, allocatable :: vny  (:) ! V face normal unit vector y component [m]
     real, allocatable :: vnz  (:) ! V face normal unit vector z component [m]

     real, allocatable :: area    (:) ! cell surface area [m^2]
     real, allocatable :: xew     (:) ! earth x coord of sfc W points [m]
     real, allocatable :: yew     (:) ! earth y coord of sfc W points [m]
     real, allocatable :: zew     (:) ! earth z coord of sfc W points [m]
     real, allocatable :: glatw   (:) ! latitude of sfc cell W points [deg]
     real, allocatable :: glonw   (:) ! longitude of sfc cell W points [deg]
     real, allocatable :: topw    (:) ! topographic height of sfc W points [m]
     real, allocatable :: bathym  (:) ! bathymetric height of sfc W points [m]
     real, allocatable :: wnx     (:) ! norm unit vector x comp of sfc cells [m]
     real, allocatable :: wny     (:) ! norm unit vector y comp of sfc cells [m]
     real, allocatable :: wnz     (:) ! norm unit vector z comp of sfc cells [m]
     real, allocatable :: dzt_bot (:) ! surface similarity grid-height [m]

     real, allocatable :: gxps_coef(:)
     real, allocatable :: gyps_coef(:)

     ! Surface type (land/vegetation, lake, or sea)

     integer, allocatable :: leaf_class (:) ! leaf ("vegetation") class
     integer, allocatable :: ioge       (:) ! integer array for storing database data
     logical, allocatable :: swm_active (:) ! shallow-water/flood model active in these sfcg cells

     ! Atmospheric near-surface properties

     real, allocatable :: vels     (:) ! wind speed [m/s]
     real, allocatable :: prss     (:) ! surface air pressure [Pa]
     real, allocatable :: canexner (:) ! Exner function in canopy [ ]
     real, allocatable :: rhos     (:) ! atmos near-surface density [kg_dry/m^3]
     real, allocatable :: airtemp  (:) ! air temperature [K]
     real, allocatable :: airtheta (:) ! air potential temperature [K]
     real, allocatable :: airrrv   (:) ! air mixing ratio [kg_vap/kg_dryair]
     real, allocatable :: airco2   (:) ! air mixing ratio of CO2 [kg_co2/kg_dryair]

     ! Canopy to atmosphere turbulent flux quantities

     real, allocatable :: ustar   (:) ! friction velocity [m/s]
     real, allocatable :: vkmsfc  (:) ! surface drag coefficient [kg/(m s)]
     real, allocatable :: vkhsfc  (:) ! surface heat and vapor transfer coefficient [kg/(m s)]
     real, allocatable :: sfluxt  (:) ! canopy-to-atm sensible heat flux [W m^-2]
     real, allocatable :: sfluxr  (:) ! canopy-to-atm water vapor flux [kg_vap m^-2 s^-1]
     real, allocatable :: sfluxc  (:)
     real, allocatable :: ggaer   (:) ! surface aerodynamic conductance [m/s]
     real, allocatable :: zeta    (:) ! surface z / M-O length [ ]
     real, allocatable :: wthv    (:) ! surface buoyancy flux [K m/s]

     ! Radiative flux quantities

     real, allocatable :: albedo_beam   (:) ! surface s/w beam albedo [0-1]
     real, allocatable :: albedo_diffuse(:) ! surface s/w diffuse albedo [0-1]
     real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
     real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
     real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
     real, allocatable :: rlong_albedo  (:) ! surface l/w lbedo [0-1]
     real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]

     ! Precipitation fluxes and runoff

     real, allocatable :: pcpg  (:) ! new pcp amount this timestep [kg/m^2]
     real, allocatable :: qpcpg (:) ! new pcp energy this timestep [J/m^2]
     real, allocatable :: dpcpg (:) ! new pcp depth this timestep [m]
     real, allocatable :: runoff(:) ! new runoff mass this timestep [kg/m^2]

     ! Canopy and surface quantities:

     real, allocatable :: can_depth(:) ! "canopy" depth for heat & vap capacity [m]
     real, allocatable :: cantemp  (:) ! "canopy" air temperature [K]
     real, allocatable :: canrrv   (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]
     real, allocatable :: rough    (:) ! roughness height [m]
     real, allocatable :: head1    (:) ! water surface hydraulic head (rel to topo) [m]

     ! Shallow Water Model (SWM) quantities:
     ! (Members of sfcg because V-point variables/indices are required)

     integer, allocatable :: iwdepv(:) ! Advective donor cell IW index for V face

     real, allocatable :: vc(:)     ! Current velocity [m/s]
     real, allocatable :: vmp(:)    ! Past velocity * depth [m^2/s]
     real, allocatable :: vmc(:)    ! Current velocity * depth [m^2/s]
     real, allocatable :: vort(:)   ! Current vorticity [s^-1]

     real, allocatable :: hflux_wat(:) ! Volume flux across V face (VMC * DNU) [m^3/s]
     real, allocatable :: hflux_enr(:) ! Energy flux across V face [J/s]
     real, allocatable :: hflux_vxe(:)
     real, allocatable :: hflux_vye(:)
     real, allocatable :: hflux_vze(:)

  End type surface_grid_vars

  type (surface_grid_vars), target :: sfcg

! TRI-GRID INFORMATION FOR INDEPENDENT REFINING OF SFC GRID

  integer :: nsfcgrids
  integer :: sfcgrid_res_factor
  integer :: nxp_sfc

  integer, target :: nsfcgrdll(maxgrds)
  real   , target :: sfcgrdrad(maxgrds,maxngrdll)
  real   , target :: sfcgrdlat(maxgrds,maxngrdll)
  real   , target :: sfcgrdlon(maxgrds,maxngrdll)

! SWM ZONE INFORMATION

  integer :: nswmzons

  integer, target :: nswmzonll(maxgrds)
  real   , target :: swmzonrad(maxgrds,maxngrdll)
  real   , target :: swmzonlat(maxgrds,maxngrdll)
  real   , target :: swmzonlon(maxgrds,maxngrdll)

Contains

!=========================================================================

  subroutine alloc_sfcgrid1(mmsfc0, mvsfc0, mwsfc0, alloc_xyzew)

  use misc_coms, only: rinit

  implicit none

  integer, intent(in) :: mmsfc0, mvsfc0, mwsfc0

  logical, intent(in), optional :: alloc_xyzew
  logical                       :: do_allocw
  integer                       :: imsfc, ivsfc, iwsfc

  do_allocw = .true.
  if (present(alloc_xyzew)) then
     do_allocw = alloc_xyzew
  endif

  ! Allocate surface grid tables

  allocate( itab_msfc(mmsfc0) )
  allocate( itab_vsfc(mvsfc0) )
  allocate( itab_wsfc(mwsfc0) )

  ! Initialize tables to default values

  !$omp parallel
  !$omp do
  do imsfc = 1, mmsfc0
     itab_msfc(imsfc) = itab_msfc_vars( imglobe =  1, &
                                        irank   = -1, &
                                        ivn     =  1, &
                                        iwn     =  1, &
                                        imn     =  1  )
  enddo
  !$omp end do

  !$omp do
  do ivsfc = 1, mvsfc0
     itab_vsfc(ivsfc) = itab_vsfc_vars( ivglobe =  1,  &
                                        irank   = -1,  &
                                        imn     =  1,  &
                                        iwn     =  1,  &
                                        cosv    =  0., &
                                        sinv    =  0., &
                                        dxps    =  0., &
                                        dyps    =  0.  )
  enddo
  !$omp end do

  !$omp do
  do iwsfc = 1, mwsfc0
     itab_wsfc(iwsfc) = itab_wsfc_vars( iwglobe   =  1,  &
                                        irank     = -1,  &
                                        npoly     =  0,  &
                                        imn       =  1,  &
                                        ivn       =  1,  &
                                        iwn       =  1,  &
                                        dirv      =  0., &
                                        farm      =  0., &
                                        farv      =  0., &
                                        gxps1     =  0., &
                                        gyps1     =  0., &
                                        gxps2     =  0., &
                                        gyps2     =  0., &
                                        ecvec_vx  =  0., &
                                        ecvec_vy  =  0., &
                                        ecvec_vz  =  0., &
                                        nwatm     =  0,  &
                                        iwatm     =  1,  &
                                        kwatm     =  0,  &
                                        arc       =  0., &
                                        arcoarsfc =  0., &
                                        arcoariw  =  0., &
                                        arcoarkw  =  0.  )
  enddo
  !$omp end do
  !$omp end parallel

  ! Allocate sfcg member arrays

  allocate( sfcg%xem  (mmsfc0) )
  allocate( sfcg%yem  (mmsfc0) )
  allocate( sfcg%zem  (mmsfc0) )
  allocate( sfcg%glatm(mmsfc0) )
  allocate( sfcg%glonm(mmsfc0) )
! allocate( sfcg%topm (mmsfc0) )
  allocate( sfcg%arm0 (mmsfc0) )

  allocate( sfcg%xev (mvsfc0) )
  allocate( sfcg%yev (mvsfc0) )
  allocate( sfcg%zev (mvsfc0) )
  allocate( sfcg%dnu (mvsfc0) )
  allocate( sfcg%dniu(mvsfc0) )
  allocate( sfcg%dnv (mvsfc0) )
  allocate( sfcg%dniv(mvsfc0) )
  allocate( sfcg%unx (mvsfc0) )
  allocate( sfcg%uny (mvsfc0) )
  allocate( sfcg%unz (mvsfc0) )
  allocate( sfcg%vnx (mvsfc0) )
  allocate( sfcg%vny (mvsfc0) )
  allocate( sfcg%vnz (mvsfc0) )

  allocate( sfcg%area      (mwsfc0) )
  allocate( sfcg%glatw     (mwsfc0) )
  allocate( sfcg%glonw     (mwsfc0) )
  allocate( sfcg%topw      (mwsfc0) )
  allocate( sfcg%bathym    (mwsfc0) )
  allocate( sfcg%wnx       (mwsfc0) )
  allocate( sfcg%wny       (mwsfc0) )
  allocate( sfcg%wnz       (mwsfc0) )
  allocate( sfcg%dzt_bot   (mwsfc0) )

  allocate( sfcg%gxps_coef (mwsfc0) )
  allocate( sfcg%gyps_coef (mwsfc0) )
  allocate( sfcg%leaf_class(mwsfc0) )
  allocate( sfcg%ioge      (mwsfc0) )
  allocate( sfcg%swm_active(mwsfc0) )

  if (do_allocw) then
     allocate( sfcg%xew(mwsfc0) )
     allocate( sfcg%yew(mwsfc0) )
     allocate( sfcg%zew(mwsfc0) )
  endif

  ! Initialize sfcg arrays to default values

  !$omp parallel
  !$omp do
  do imsfc = 1, mmsfc0
     sfcg%xem  (imsfc) = rinit
     sfcg%yem  (imsfc) = rinit
     sfcg%zem  (imsfc) = rinit
     sfcg%glatm(imsfc) = rinit
     sfcg%glonm(imsfc) = rinit
!    sfcg%topm (imsfc) = rinit
     sfcg%arm0 (imsfc) = 0.
  enddo
  !$omp end do

  !$omp do
  do ivsfc = 1, mvsfc0
     sfcg%xev (ivsfc) = rinit
     sfcg%yev (ivsfc) = rinit
     sfcg%zev (ivsfc) = rinit
     sfcg%dnu (ivsfc) = rinit
     sfcg%dniu(ivsfc) = rinit
     sfcg%dnv (ivsfc) = rinit
     sfcg%dniv(ivsfc) = rinit
     sfcg%unx (ivsfc) = rinit
     sfcg%uny (ivsfc) = rinit
     sfcg%unz (ivsfc) = rinit
     sfcg%vnx (ivsfc) = rinit
     sfcg%vny (ivsfc) = rinit
     sfcg%vnz (ivsfc) = rinit
  enddo
  !$omp end do

  !$omp do
  do iwsfc = 1, mwsfc0
     sfcg%area      (iwsfc) = 0.
     sfcg%glatw     (iwsfc) = rinit
     sfcg%glonw     (iwsfc) = rinit
     sfcg%topw      (iwsfc) = rinit
     sfcg%bathym    (iwsfc) = rinit
     sfcg%wnx       (iwsfc) = rinit
     sfcg%wny       (iwsfc) = rinit
     sfcg%wnz       (iwsfc) = rinit
     sfcg%dzt_bot   (iwsfc) = rinit
     sfcg%gxps_coef (iwsfc) = rinit
     sfcg%gyps_coef (iwsfc) = rinit
     sfcg%leaf_class(iwsfc) = 0
     sfcg%ioge      (iwsfc) = 0
     sfcg%swm_active(iwsfc) = .false.

     if (do_allocw) then
        sfcg%xew(iwsfc) = rinit
        sfcg%yew(iwsfc) = rinit
        sfcg%zew(iwsfc) = rinit
     endif

  enddo
  !$omp end do
  !$omp end parallel

  end subroutine alloc_sfcgrid1

!=========================================================================

  subroutine alloc_sfcgrid2()

     use misc_coms, only: rinit
     use mem_co2,   only: co2flag

     implicit none

     integer :: iwsfc, ivsfc

!    Allocate and initialize surface grid arrays

     allocate (sfcg%vels          (mwsfc))
     allocate (sfcg%prss          (mwsfc))
     allocate (sfcg%canexner      (mwsfc))
     allocate (sfcg%rhos          (mwsfc))
     allocate (sfcg%airtemp       (mwsfc))
     allocate (sfcg%airtheta      (mwsfc))
     allocate (sfcg%airrrv        (mwsfc))

     allocate (sfcg%ustar         (mwsfc))
     allocate (sfcg%vkmsfc        (mwsfc))
     allocate (sfcg%vkhsfc        (mwsfc))
     allocate (sfcg%sfluxt        (mwsfc))
     allocate (sfcg%sfluxr        (mwsfc))

     if (co2flag /= 0) then
        allocate (sfcg%airco2     (mwsfc))
        allocate (sfcg%sfluxc     (mwsfc))
     endif

     allocate (sfcg%ggaer         (mwsfc))
     allocate (sfcg%wthv          (mwsfc))
     allocate (sfcg%albedo_beam   (mwsfc))
     allocate (sfcg%albedo_diffuse(mwsfc))
     allocate (sfcg%rshort        (mwsfc))
     allocate (sfcg%rshort_diffuse(mwsfc))
     allocate (sfcg%rlong         (mwsfc))
     allocate (sfcg%rlong_albedo  (mwsfc))
     allocate (sfcg%rlongup       (mwsfc))
     allocate (sfcg%pcpg          (mwsfc))
     allocate (sfcg%qpcpg         (mwsfc))
     allocate (sfcg%dpcpg         (mwsfc))
     allocate (sfcg%runoff        (mwsfc))
     allocate (sfcg%can_depth     (mwsfc))
     allocate (sfcg%cantemp       (mwsfc))
     allocate (sfcg%canrrv        (mwsfc))
     allocate (sfcg%rough         (mwsfc))
     allocate (sfcg%head1         (mwsfc))

     if (nswmzons > 0) then
        allocate (sfcg%iwdepv     (mvsfc))

        allocate (sfcg%vc         (mvsfc))
        allocate (sfcg%vmp        (mvsfc))
        allocate (sfcg%vmc        (mvsfc))
        allocate (sfcg%vort          (mmsfc)) ; sfcg%vort           = 0.0

        allocate (sfcg%hflux_wat  (mvsfc))
        allocate (sfcg%hflux_enr  (mvsfc))
        allocate (sfcg%hflux_vxe  (mvsfc))
        allocate (sfcg%hflux_vye  (mvsfc))
        allocate (sfcg%hflux_vze  (mvsfc))
     endif

     !$omp parallel do
     do iwsfc = 1, mwsfc

        if ( allocated( sfcg%vels          ) ) sfcg%vels          (iwsfc) = rinit
        if ( allocated( sfcg%prss          ) ) sfcg%prss          (iwsfc) = rinit
        if ( allocated( sfcg%canexner      ) ) sfcg%canexner      (iwsfc) = rinit
        if ( allocated( sfcg%rhos          ) ) sfcg%rhos          (iwsfc) = rinit
        if ( allocated( sfcg%airtemp       ) ) sfcg%airtemp       (iwsfc) = rinit
        if ( allocated( sfcg%airtheta      ) ) sfcg%airtheta      (iwsfc) = rinit
        if ( allocated( sfcg%airrrv        ) ) sfcg%airrrv        (iwsfc) = rinit

        if ( allocated( sfcg%ustar         ) ) sfcg%ustar         (iwsfc) = rinit
        if ( allocated( sfcg%vkmsfc        ) ) sfcg%vkmsfc        (iwsfc) = rinit
        if ( allocated( sfcg%vkhsfc        ) ) sfcg%vkhsfc        (iwsfc) = rinit
        if ( allocated( sfcg%sfluxt        ) ) sfcg%sfluxt        (iwsfc) = rinit
        if ( allocated( sfcg%sfluxr        ) ) sfcg%sfluxr        (iwsfc) = rinit

        if ( allocated( sfcg%airco2        ) ) sfcg%airco2        (iwsfc) = rinit
        if ( allocated( sfcg%sfluxc        ) ) sfcg%sfluxc        (iwsfc) = rinit

        if ( allocated( sfcg%ggaer         ) ) sfcg%ggaer         (iwsfc) = rinit
        if ( allocated( sfcg%wthv          ) ) sfcg%wthv          (iwsfc) = rinit
        if ( allocated( sfcg%albedo_beam   ) ) sfcg%albedo_beam   (iwsfc) = 0.0
        if ( allocated( sfcg%albedo_diffuse) ) sfcg%albedo_diffuse(iwsfc) = 0.0
        if ( allocated( sfcg%rshort        ) ) sfcg%rshort        (iwsfc) = 0.0
        if ( allocated( sfcg%rshort_diffuse) ) sfcg%rshort_diffuse(iwsfc) = 0.0
        if ( allocated( sfcg%rlong         ) ) sfcg%rlong         (iwsfc) = 0.0
        if ( allocated( sfcg%rlong_albedo  ) ) sfcg%rlong_albedo  (iwsfc) = 0.0
        if ( allocated( sfcg%rlongup       ) ) sfcg%rlongup       (iwsfc) = 0.0
        if ( allocated( sfcg%pcpg          ) ) sfcg%pcpg          (iwsfc) = 0.0
        if ( allocated( sfcg%qpcpg         ) ) sfcg%qpcpg         (iwsfc) = 0.0
        if ( allocated( sfcg%dpcpg         ) ) sfcg%dpcpg         (iwsfc) = 0.0
        if ( allocated( sfcg%runoff        ) ) sfcg%runoff        (iwsfc) = 0.0
        if ( allocated( sfcg%can_depth     ) ) sfcg%can_depth     (iwsfc) = rinit
        if ( allocated( sfcg%cantemp       ) ) sfcg%cantemp       (iwsfc) = rinit
        if ( allocated( sfcg%canrrv        ) ) sfcg%canrrv        (iwsfc) = rinit
        if ( allocated( sfcg%rough         ) ) sfcg%rough         (iwsfc) = rinit
        if ( allocated( sfcg%head1         ) ) sfcg%head1         (iwsfc) = 0.0

     enddo
     !$omp end parallel do

     !$omp parallel do
     do ivsfc = 1, mvsfc

        if ( allocated( sfcg%iwdepv    ) ) sfcg%iwdepv    (ivsfc) = 0
        if ( allocated( sfcg%vc        ) ) sfcg%vc        (ivsfc) = 0.0
        if ( allocated( sfcg%vmp       ) ) sfcg%vmp       (ivsfc) = 0.0
        if ( allocated( sfcg%vmc       ) ) sfcg%vmc       (ivsfc) = 0.0
        if ( allocated( sfcg%hflux_wat ) ) sfcg%hflux_wat (ivsfc) = 0.0
        if ( allocated( sfcg%hflux_enr ) ) sfcg%hflux_enr (ivsfc) = 0.0
        if ( allocated( sfcg%hflux_vxe ) ) sfcg%hflux_vxe (ivsfc) = 0.0
        if ( allocated( sfcg%hflux_vye ) ) sfcg%hflux_vye (ivsfc) = 0.0
        if ( allocated( sfcg%hflux_vze ) ) sfcg%hflux_vze (ivsfc) = 0.0

     enddo
     !$omp end parallel do

  end subroutine alloc_sfcgrid2

!=========================================================================

   subroutine filltab_sfcg()

     use var_tables, only: increment_vtable

     implicit none

     if (allocated(sfcg%vels))           call increment_vtable('SFCG%VELS',           'CW', rvar1=sfcg%vels)
     if (allocated(sfcg%prss))           call increment_vtable('SFCG%PRSS',           'CW', rvar1=sfcg%prss)
     if (allocated(sfcg%canexner))       call increment_vtable('SFCG%CANEXNER',       'CW', rvar1=sfcg%canexner)
     if (allocated(sfcg%rhos))           call increment_vtable('SFCG%RHOS',           'CW', rvar1=sfcg%rhos)
     if (allocated(sfcg%airtemp))        call increment_vtable('SFCG%AIRTEMP',        'CW', rvar1=sfcg%airtemp)
     if (allocated(sfcg%airtheta))       call increment_vtable('SFCG%AIRTHETA',       'CW', rvar1=sfcg%airtheta)
     if (allocated(sfcg%airrrv))         call increment_vtable('SFCG%AIRRRV',         'CW', rvar1=sfcg%airrrv)
     if (allocated(sfcg%airco2))         call increment_vtable('SFCG%AIRCO2',         'CW', rvar1=sfcg%airco2)
     if (allocated(sfcg%ustar))          call increment_vtable('SFCG%USTAR',          'CW', rvar1=sfcg%ustar)
     if (allocated(sfcg%vkhsfc))         call increment_vtable('SFCG%VKHSFC',         'CW', rvar1=sfcg%vkhsfc)
     if (allocated(sfcg%sfluxt))         call increment_vtable('SFCG%SFLUXT',         'CW', rvar1=sfcg%sfluxt)
     if (allocated(sfcg%sfluxr))         call increment_vtable('SFCG%SFLUXR',         'CW', rvar1=sfcg%sfluxr)
     if (allocated(sfcg%sfluxc))         call increment_vtable('SFCG%SFLUXC',         'CW', rvar1=sfcg%sfluxc)
     if (allocated(sfcg%ggaer))          call increment_vtable('SFCG%GGAER',          'CW', rvar1=sfcg%ggaer)
     if (allocated(sfcg%zeta))           call increment_vtable('SFCG%ZETA',           'CW', rvar1=sfcg%zeta)
     if (allocated(sfcg%wthv))           call increment_vtable('SFCG%WTHV',           'CW', rvar1=sfcg%wthv)
     if (allocated(sfcg%albedo_beam))    call increment_vtable('SFCG%ALBEDO_BEAM',    'CW', rvar1=sfcg%albedo_beam)
     if (allocated(sfcg%albedo_diffuse)) call increment_vtable('SFCG%ALBEDO_DIFFUSE', 'CW', rvar1=sfcg%albedo_diffuse)
     if (allocated(sfcg%rshort))         call increment_vtable('SFCG%RSHORT',         'CW', rvar1=sfcg%rshort)
     if (allocated(sfcg%rshort_diffuse)) call increment_vtable('SFCG%RSHORT_DIFFUSE', 'CW', rvar1=sfcg%rshort_diffuse)
     if (allocated(sfcg%rlong))          call increment_vtable('SFCG%RLONG',          'CW', rvar1=sfcg%rlong)
     if (allocated(sfcg%rlong_albedo))   call increment_vtable('SFCG%RLONG_ALBEDO',   'CW', rvar1=sfcg%rlong_albedo)
     if (allocated(sfcg%rlongup))        call increment_vtable('SFCG%RLONGUP',        'CW', rvar1=sfcg%rlongup)
     if (allocated(sfcg%pcpg))           call increment_vtable('SFCG%PCPG',           'CW', rvar1=sfcg%pcpg)
     if (allocated(sfcg%qpcpg))          call increment_vtable('SFCG%QPCPG',          'CW', rvar1=sfcg%qpcpg)
     if (allocated(sfcg%dpcpg))          call increment_vtable('SFCG%DPCPG',          'CW', rvar1=sfcg%dpcpg)
     if (allocated(sfcg%runoff))         call increment_vtable('SFCG%RUNOFF',         'CW', rvar1=sfcg%runoff)
!    if (allocated(sfcg%can_depth))      call increment_vtable('SFCG%CAN_DEPTH',      'CW', rvar1=sfcg%can_depth)
     if (allocated(sfcg%cantemp))        call increment_vtable('SFCG%CANTEMP',        'CW', rvar1=sfcg%cantemp)
     if (allocated(sfcg%canrrv))         call increment_vtable('SFCG%CANRRV',         'CW', rvar1=sfcg%canrrv)
     if (allocated(sfcg%rough))          call increment_vtable('SFCG%ROUGH',          'CW', rvar1=sfcg%rough)
     if (allocated(sfcg%head1))          call increment_vtable('SFCG%HEAD1',          'CW', rvar1=sfcg%head1)

     if (allocated(sfcg%vc))             call increment_vtable('SFCG%VC',             'CV', rvar1=sfcg%vc)
     if (allocated(sfcg%vmp))            call increment_vtable('SFCG%VMP',            'CV', rvar1=sfcg%vmp)
     if (allocated(sfcg%vmc))            call increment_vtable('SFCG%VMC',            'CV', rvar1=sfcg%vmc)
     if (allocated(sfcg%vort))           call increment_vtable('SFCG%VORT',           'CM', rvar1=sfcg%vort)

   end subroutine filltab_sfcg

!===============================================================================

  subroutine fill_jtab_sfcg(mwa)

  use mem_ijtabs, only: itab_w
  implicit none

  integer, intent(in) :: mwa

  integer :: iwsfc, ivsfc, imsfc, iw1, iw2, iw3, j
  integer :: iw, ipass, jsfc2, np

  ! JTAB_WSFC_SWM, restricted to SEA points that are SWM-active
  ! (and their immediate neighbors)

  j = 0
  do iwsfc = 2,mwsfc
!    if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) j = j + 1
     if (sfcg%leaf_class(iwsfc) == 0) then

        np = itab_wsfc(iwsfc)%npoly
        if ( sfcg%swm_active(iwsfc) .or. &
             any(sfcg%swm_active( itab_wsfc(iwsfc)%iwn(1:np) )) ) j = j + 1

     endif
  enddo
  jtab_wsfc_swm%jend = j

  allocate( jtab_wsfc_swm%iwsfc(j) )

  !$omp parallel do
  !$ do j = 1, jtab_wsfc_swm%jend
  !$    jtab_wsfc_swm%iwsfc(j) = 0
  !$ enddo
  !$omp end parallel do

  j = 0
  do iwsfc = 2,mwsfc
!    if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
!       j = j + 1
!       jtab_wsfc_swm%iwsfc(j) = iwsfc
!    endif
     if (sfcg%leaf_class(iwsfc) == 0) then

        np = itab_wsfc(iwsfc)%npoly
        if ( sfcg%swm_active(iwsfc) .or. &
             any(sfcg%swm_active( itab_wsfc(iwsfc)%iwn(1:np) )) ) then
           j = j + 1
           jtab_wsfc_swm%iwsfc(j) = iwsfc
        endif

     endif
  enddo

  ! JTAB_VSFC_SWM, includes all VSFC edges adjacent to at least one SWM-active cell

  j = 0
  do ivsfc = 2,mvsfc
     iw1 = itab_vsfc(ivsfc)%iwn(1)
     iw2 = itab_vsfc(ivsfc)%iwn(2)
     if (sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) j = j + 1
  enddo
  jtab_vsfc_swm%jend = j

  allocate( jtab_vsfc_swm%ivsfc(j) )

  !$omp parallel do
  !$ do j = 1, jtab_vsfc_swm%jend
  !$    jtab_vsfc_swm%ivsfc(j) = 0
  !$ enddo
  !$omp end parallel do

  j = 0
  do ivsfc = 2,mvsfc
     iw1 = itab_vsfc(ivsfc)%iwn(1)
     iw2 = itab_vsfc(ivsfc)%iwn(2)
     if (sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) then
        j = j + 1
        jtab_vsfc_swm%ivsfc(j) = ivsfc
     endif
  enddo

  ! JTAB_MSFC_SWM, restricted to M vertices surrounded by 3 SEA cells that are SWM-active

  j = 0
  do imsfc = 2,mmsfc
     iw1 = itab_msfc(imsfc)%iwn(1)
     iw2 = itab_msfc(imsfc)%iwn(2)
     iw3 = itab_msfc(imsfc)%iwn(3)

     if (iw1 > 1 .and. iw2 > 1 .and. iw3 > 1) then
        if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
            (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0) .and. &
            (sfcg%swm_active(iw3) .and. sfcg%leaf_class(iw3) == 0)) j = j + 1
     endif
  enddo
  jtab_msfc_swm%jend = j

  allocate( jtab_msfc_swm%imsfc(j) )

  !$omp parallel do
  !$ do j = 1, jtab_msfc_swm%jend
  !$    jtab_msfc_swm%imsfc(j) = 0
  !$ enddo
  !$omp end parallel do

  j = 0
  do imsfc = 2,mmsfc
     iw1 = itab_msfc(imsfc)%iwn(1)
     iw2 = itab_msfc(imsfc)%iwn(2)
     iw3 = itab_msfc(imsfc)%iwn(3)
     if (iw1 > 1 .and. iw2 > 1 .and. iw3 > 1) then
        if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
            (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0) .and. &
            (sfcg%swm_active(iw3) .and. sfcg%leaf_class(iw3) == 0)) then
           j = j + 1
           jtab_msfc_swm%imsfc(j) = imsfc
        endif
     endif
  enddo

  ! Do two passes through the following code to build lists of attached
  ! land, lake, and sea cells for each IW column

  do ipass = 1, 2

     ! Set land, lake and sea counters to zero for all atmosphere IW columns

     !$omp parallel do
     do iw = 1, mwa
        itab_w(iw)%jsfc2  = 0
        itab_w(iw)%jland2 = 0
        itab_w(iw)%jlake2 = 0
        itab_w(iw)%jsea2  = 0
     enddo
     !$omp end parallel do

     ! Loop over all SFC grid cells and get atmosphere column iw index

     do iwsfc = 2,mwsfc
        do j = 1,itab_wsfc(iwsfc)%nwatm

           iw = itab_wsfc(iwsfc)%iwatm(j) ! local index (or 1 if not on this rank)

           ! Skip this j/iw point if it is not on current node

           if (iw < 2) cycle

           itab_w(iw)%jsfc2 = itab_w(iw)%jsfc2 + 1

           ! Increment land, lake, or sea cell counter for IW column

           if (sfcg%leaf_class(iwsfc) >= 2) then
              itab_w(iw)%jland2 = itab_w(iw)%jland2 + 1
           elseif (sfcg%leaf_class(iwsfc) == 1) then
              itab_w(iw)%jlake2 = itab_w(iw)%jlake2 + 1
           else
              itab_w(iw)%jsea2 = itab_w(iw)%jsea2 + 1
           endif

           ! If second pass, enter SFC cell indices into itab_w(iw) member arrays

           if (ipass == 2) then
              itab_w(iw)%iwsfc(itab_w(iw)%jsfc2) = iwsfc  ! local-rank iwsfc and iw indices
              itab_w(iw)%jasfc(itab_w(iw)%jsfc2) = j
           endif

        enddo
     enddo

     ! If first pass, allocate iwsfc member of itab_w(iw)

     if (ipass == 1) then

        !$omp parallel do private(jsfc2)
        do iw = 2, mwa
           jsfc2 = itab_w(iw)%jsfc2

           allocate(itab_w(iw)%iwsfc(jsfc2))
           allocate(itab_w(iw)%jasfc(jsfc2))

           itab_w(iw)%iwsfc(:) = 0
           itab_w(iw)%jasfc(:) = 0
        enddo
        !$omp end parallel do

     endif

     ! If second pass, convert j-index members of itab_w(iw) if land, lake,
     ! and/or sea points are present

     if (ipass == 2) then

        !$omp parallel do
        do iw = 2, mwa
           itab_w(iw)%jland1 = 1
           itab_w(iw)%jlake1 = itab_w(iw)%jland2 + 1
           itab_w(iw)%jsea1  = itab_w(iw)%jland2 + itab_w(iw)%jlake2 + 1

           if (itab_w(iw)%jsea2  > 0) itab_w(iw)%jsea2  = itab_w(iw)%jland2 + itab_w(iw)%jlake2 + itab_w(iw)%jsea2
           if (itab_w(iw)%jlake2 > 0) itab_w(iw)%jlake2 = itab_w(iw)%jland2 + itab_w(iw)%jlake2

        enddo
        !$omp end parallel do

     endif

  enddo  ! ipass

  end subroutine fill_jtab_sfcg

End Module mem_sfcg
