module umwm_module

  use consts_coms, only: pi2, grav

  implicit none

  real, parameter :: euler    = exp(1.)    ! e
  real, parameter :: eulerinv = 1. / euler ! 1/e
  real, parameter :: rhosw    = 1025.      ! Avg density of seawater at shallow depths [kg/m^3]
  real, parameter :: rhoswi   = 1. / rhosw ! Inverse of rhosw [m^3/kg]

  ! VARIABLES INPUT FROM UMWM NAMELIST (NOW SET AS PARAMETERS HERE)

  !Bob - orig values were: om=37,pm=36,fmax=2.0,fprog=2.0
  !Bob - also, orig lm in stokes = 42; Bob changed lm to 1

  integer, parameter :: om     = 32      ! number of frequency/wavenumber bins
  integer, parameter :: pm     = 32      ! number of direction bins
  logical, parameter :: stokes = .true.  ! Output switch for Stokes drift velocity fields

  real, parameter :: fmin      = 0.035   ! Frequency of lowest frequency bin [Hz]
  real, parameter :: fmax      = 2.0     ! Frequency of highest frequency bin [Hz]
  real, parameter :: fprog     = 2.0     ! Highest prognostic frequency bin [Hz]
  real, parameter :: nu_air    = 1.56E-5 ! Kinematic viscosity of air [m^2/s]
  real, parameter :: nu_water  = 0.90E-6 ! Kinematic viscosity of water [m^2/s]
  real, parameter :: sfct      = 0.07    ! Surface tension [N/m]
  real, parameter :: gustiness = 0.0     ! Random wind gustiness factor (should be between 0 and 0.2)
  real, parameter :: dmin      = 10.0    ! Minimum allowed ocean depth [m]
  real, parameter :: explim    = 0.9     ! Exponent limiter (0.69 ~ 100% growth)
  real, parameter :: sin_fac   = 0.11    ! Input factor from following winds
  real, parameter :: sin_diss1 = 0.10    ! Damping factor from opposing winds
  real, parameter :: sin_diss2 = 0.001   ! Damping factor from swell overrunning wind

! UMWM defaults:
! real, parameter :: sds_fac   = 42.     ! Breaking dissipation factor
! real, parameter :: sds_power = 2.4     ! Saturation spectrum power
! real, parameter :: mss_fac   = 360.    ! Mean-square-slope adjustment to Sds

! Donelan (2012, 2018)
  real, parameter :: sds_fac   = 45.0    ! Breaking dissipation factor
  real, parameter :: sds_power = 2.5     ! Saturation spectrum power
  real, parameter :: mss_fac   = 240.    ! Mean-square-slope adjustment to Sds

  real, parameter :: snl_fac   = 5.0     ! Wave energy downshifting factor
  real, parameter :: sdt_fac   = 0.002   ! Dissipation due to turbulence factor
  real, parameter :: sbf_fac   = 0.003   ! Bottom friction coefficient [m/s]
  real, parameter :: sbp_fac   = 0.003   ! Bottom percolation coefficient [m/s]
  real, parameter :: fice_lth  = 0.05    ! Sea ice fraction - lower threshold for attenuation
  real, parameter :: fice_uth  = 0.60    ! Sea ice fraction - upper threshold for attenuation

  ! Other variables

  integer :: umwmflg       ! if wave model is active, from OLAM namlist
  real    :: dt_olam       ! Time step from OLAM [s]
  real    :: dta           ! Time step for advection and updating [s]
  real    :: dtr           ! Time step for rotation (for now, set equal to dta) [s]

  ! Various constants:

  real, parameter :: dth           = pi2 / real(pm)              ! wave direction increment [rad]
  real, parameter :: oneovdth      = 1. / dth                    ! 1 / dth
  real, parameter :: dthg          = dth * grav                  ! dth * grav
  real, parameter :: dlnf          = log(fmax/fmin) / real(om-1) ! frequency increment
  real, parameter :: fieldscale1   = sin_diss1 / sin_fac         ! sin_diss1 / sin_fac [ ]
  real, parameter :: fieldscale2   = sin_diss2 / sin_diss1       ! sin_diss2 / sin_diss1 [ ]
  real, parameter :: inv_sds_power = 1. / sds_power              ! 1 / sds_power
  real, parameter :: twopisds_fac  = pi2 * sds_fac               ! 2 * pi * sds_fac

  ! diffusion values in 2 frequencies:

  real, parameter :: bfa = exp(-16. * dlnf * dlnf)
  real, parameter :: bfb = exp(-64. * dlnf * dlnf)
  real, parameter :: bf1 = bfa / (bfa + bfb)
  real, parameter :: bf2 = bfb / (bfa + bfb)

  ! The following arrays are included in umwm data structure to prevent any
  ! conflicts or confusion with similar quantities or variable names in OLAM

                                         !I = assigned at initialization only
                                         !P = prognosed or updated as part of prognostic quantitites
                                         !D = diagnostic only for output; not used for prognosing
  Type umwm_vars
     real,    allocatable :: seadep (:) !I (msea) sea depth [m]
     real,    allocatable :: alogzs (:) !I (msea) log( dzt_bot )
     real,    allocatable :: uwind  (:) !P (msea) eastward wind velocity [m/s]
     real,    allocatable :: vwind  (:) !P (msea) northward wind velocity [m/s]
     real,    allocatable :: wspd   (:) !P (msea) wind speed [m/s]
     real,    allocatable :: wspd10m(:) !P (msea) wind speed at 10m height [m/s]
     real,    allocatable :: alogzo (:) !P (msea) log( zo )
     real,    allocatable :: wdir   (:) !P (msea) wind direction [rad, -pi to pi]
     real,    allocatable :: cd     (:) !Psr (msea) drag coefficient of air [ ]
     real,    allocatable :: ustar  (:) !Psr (msea) friction velocity of air [m/s]
     real,    allocatable :: taux   (:) !Psr (msea) stress (momentum flux) [n/m^2]
     real,    allocatable :: tauy   (:) !Psr (msea) stress (momentum flux) [n/m^2]
     logical, allocatable :: iactive(:) !P (msea) if wave model was activated on this cell
  End type umwm_vars

  type (umwm_vars)     :: umwm

  integer, allocatable :: oc      (:) !I (msea) cut-off frequency index (maximum prognostic)
! integer, allocatable :: pl      (:) !I (pm)   anti-clockwise directional index
! integer, allocatable :: pr      (:) !I (pm)   clockwise directional index

  real, allocatable :: freq       (:) !I (om)   Wave frequency of frequency bins ([Hz]
  real, allocatable :: fsdsf      (:) !I (om)   Frequency bin factor for sds forcing

  real, allocatable :: th         (:) !I (pm)   Wave direction [rad; -pi to pi]
  real, allocatable :: cth        (:) !I (pm)   Cosine of wave direction th [ ]
  real, allocatable :: sth        (:) !I (pm)   Sine of wave direction th [ ]
! real, allocatable :: cth2       (:) !I (pm)   cos^2(dth*(p-1))**2 [ ]

! real, allocatable :: fcutoff    (:) !P (msea) cutoff frequency [Hz] = 0.53 * grav / wspd(msea)
! real, allocatable :: shelt      (:) !P (msea) sheltering coefficient [ ]
  real, allocatable :: ucurr      (:) !P (msea) eastward surface ocean current [m/s]
  real, allocatable :: vcurr      (:) !P (msea) northward surface ocean current [m/s]
  real, allocatable :: taux_form  (:) !Psr (msea) form drag from atmosphere, x-component [N/m^2]
  real, allocatable :: tauy_form  (:) !Psr (msea) form drag from atmosphere, y-component [N/m^2]
  real, allocatable :: taux_skin  (:) !Psr (msea) skin-induced stress, x-component [N/m^2]
  real, allocatable :: tauy_skin  (:) !Psr (msea) skin-induced stress, y-component [N/m^2]
! real, allocatable :: tailatmx   (:) !Psr (msea) atmosphere tail stress part, x-component [N/m^2]
! real, allocatable :: tailatmy   (:) !Psr (msea) atmosphere tail stress part, y-component [N/m^2]

! real, allocatable :: tailocnx   (:) !Dsr (msea) ocean tail stress part, x-component [N/m^2]
! real, allocatable :: tailocny   (:) !Dsr (msea) ocean tail stress part, y-component [N/m^2]
! real, allocatable :: taux_diag  (:) !Dsr (msea) diagnostic form drag, x-component [N/m^2]
! real, allocatable :: tauy_diag  (:) !Dsr (msea) diagnostic form drag, y-component [N/m^2]
! real, allocatable :: taux_ocntop(:) !Dsr (msea) breaking-wave stress into ocean top, x-component [N/m^2]
! real, allocatable :: tauy_ocntop(:) !Dsr (msea) breaking-wave stress into ocean top, y-component [N/m^2]
! real, allocatable :: taux_ocnbot(:) !Dsr (msea) stress from waves into ocean bottom, x-component [N/m^2]
! real, allocatable :: tauy_ocnbot(:) !Dsr (msea) stress from waves into ocean bottom, y-component [N/m^2]
! real, allocatable :: taux_snl   (:) !Dsr (msea) momentum loss due to snl, x-component [N/m^2]
! real, allocatable :: tauy_snl   (:) !Dsr (msea) momentum loss due to snl, y-component [N/m^2]
! real, allocatable :: taux1      (:) !Dsr (msea) form drag x-component: wind pushing waves [N/m^2]
! real, allocatable :: tauy1      (:) !Dsr (msea) form drag y-component: wind pushing waves [N/m^2]
! real, allocatable :: taux2      (:) !Dsr (msea) form drag x-component: waves against wind [N/m^2]
! real, allocatable :: tauy2      (:) !Dsr (msea) form drag y-component: waves against wind [N/m^2]
! real, allocatable :: taux3      (:) !Dsr (msea) form drag x-component: waves overrunning wind [N/m^2]
! real, allocatable :: tauy3      (:) !Dsr (msea) form drag y-component: waves overrunning wind [N/m^2]
! real, allocatable :: epsx_atm   (:) !Dsr (msea) wave energy growth flux, x-component [kg/s^3]
! real, allocatable :: epsy_atm   (:) !Dsr (msea) wave energy growth flux, y-component [kg/s^3]
! real, allocatable :: epsx_ocn   (:) !Dsr (msea) wave energy dissipation flux, x-component [kg/s^3]
! real, allocatable :: epsy_ocn   (:) !Dsr (msea) wave energy dissipation flux, y-component [kg/s^3]
  real, allocatable :: dwd        (:) !Ddi (msea) dominant wave direction [rad]
  real, allocatable :: dwl        (:) !Ddi (msea) dominant wavelength [m]
  real, allocatable :: dwp        (:) !Ddi (msea) dominant wave period [s]
  real, allocatable :: dcp0       (:) !Ddi (msea) dominant phase speed, intrinsic [m/s]
  real, allocatable :: dcp        (:) !Ddi (msea) dominant phase speed [m/s]
  real, allocatable :: dcg0       (:) !Ddi (msea) dominant group speed, intrinsic [m/s]
  real, allocatable :: dcg        (:) !Ddi (msea) dominant group speed [m/s]
  real, allocatable :: swh        (:) !Ddi (msea) significant wave height [m]
  real, allocatable :: mss        (:) !Ddi (msea) mean-squared slope [ ]
  real, allocatable :: mwd        (:) !Ddi (msea) mean wave direction [rad]
  real, allocatable :: mwl        (:) !Ddi (msea) mean wavelength [m]
  real, allocatable :: mwp        (:) !Ddi (msea) mean wave period [s]
! real, allocatable :: momx       (:) !Ddi (msea) total wave momentum, x-component [kg m/s]
! real, allocatable :: momy       (:) !Ddi (msea) total wave momentum, y-component [kg m/s]
! real, allocatable :: cgmxx      (:) !Ddi (msea) cg*momentum, xx-component [kg m^2/s^2]
! real, allocatable :: cgmxy      (:) !Ddi (msea) cg*momentum, xy-component [kg m^2/s^2]
! real, allocatable :: cgmyy      (:) !Ddi (msea) cg*momentum, yy-component [kg m^2/s^2]

  integer, allocatable :: opeak   (:)
  integer, allocatable :: ppeak   (:)

  real, allocatable :: cth2pp    (:,:) !I (pm)       utility array used for mss in sds routine:
  real, allocatable :: cthp_dirv (:,:) !I (pm,mvsfc) wave ray directions projected onto sfcg%dirv [ ]
  real, allocatable :: sthp_dirv (:,:) !I (pm,mvsfc) wave ray directions projected onto negative sfcg%diru [ ]
  real, allocatable :: bf1_renorm(:,:) !I (om,msea)  snl downshifting weights, used in snl routine:
  real, allocatable :: bf2_renorm(:,:) !I (om,msea)  snl downshifting weights, used in snl routine:
  real, allocatable :: cp0       (:,:) !I (om,msea)  phase velocity (pi2 * sqrt(seadep(i) / grav)) [m/s]
  real, allocatable :: cg0       (:,:) !I (om,msea)  group velocity [m/s]
  real, allocatable :: cothkd    (:,:) !I (om,msea)  coth(0.2*min(k*seadep,20))
  real, allocatable :: dwn       (:,:) !I (om,msea)  increment between discrete wavenumbers [rad/m]
  real, allocatable :: invcp0    (:,:) !I (om,msea)  1./cp0(o,i)
  real, allocatable :: fkovg     (:,:) !I (om,msea)  f(o)*k(o,i) / grav
  real, allocatable :: k         (:,:) !I (om,msea)  wavenumber [rad/m]
  real, allocatable :: k4        (:,:) !I (om,msea)  k(o,i)**4. [rad^4/m^4]
  real, allocatable :: kdk       (:,:) !I (om,msea)  k(o,i)*dwn(o,i) [rad^2/m^2]
  real, allocatable :: k3dk      (:,:) !I (om,msea)  k(o,i)**3.*dwn(o,i) [rad^4/m^4]
  real, allocatable :: l2        (:,:) !I (om,msea)  half wavelength [m]
  real, allocatable :: logl2overz(:,:) !I (om,msea)  log(min(20,l2) / sfcg%dzt_bot) [ ]
  real, allocatable :: oneoverk4 (:,:) !I (om,msea)  1./k4(o,i) [m^4/rad^4]

! real, allocatable :: sbf       (:,:) !I (om,msea)  bottom friction from roughness + percolation
  real, allocatable :: sdv       (:,:) !I (om,msea)  dissipation due to viscosity
! real, allocatable :: snl_arg   (:,:) !I (om,msea)  exp function argument for snl
! real, allocatable :: sdt       (:,:) !P (om,msea)  dissipation due to turbulence
! real, allocatable :: sice      (:,:) !P (om,msea)  wave attenuation by sea ice function [ ]

  real, allocatable :: evs  (:,:,:) !P (pm,om,msea) wave variance spectrum [m^2/(rad/m)]
  real, allocatable :: evsf (:,:,:) !P (pm,om,msea) wave variance spectrum, forward in time
! real, allocatable :: rotl (:,:,:) !P (pm,om,msea) anti-clockwise rotation, used in refraction
! real, allocatable :: rotr (:,:,:) !P (pm,om,msea) clockwise rotation, used in refraction
! real, allocatable :: sds  (:,:,:) !P (pm,om,msea) wave dissipation sink function
! real, allocatable :: snl  (:,:,:) !P (pm,om,msea) wave downshifting source/sink function
  real, allocatable :: ssin (:,:,:) !P (pm,om,msea) wind input function
! real, allocatable :: dummy(:,:,:) !P (pm,om,msea) dummy array, used in sds, advection and refraction

contains

!===============================================================================

subroutine alloc_umwm()

  use mem_sea,  only: msea
  use mem_sfcg, only: mvsfc

  implicit none

  integer :: i

  ! Allocates UMWM arrays

  allocate(umwm%seadep (msea)) ! depth
  allocate(umwm%alogzs (msea)) ! log( dzt_bot )
  allocate(umwm%alogzo (msea)) ! log( z0 )
  allocate(umwm%uwind  (msea)) ! wind component
  allocate(umwm%vwind  (msea)) ! wind component
  allocate(umwm%wspd   (msea)) ! wind speed
  allocate(umwm%wdir   (msea)) ! wind direction
  allocate(umwm%wspd10m(msea)) ! wind speed at 10m
! allocate(umwm%cd     (msea)) ! air-side drag coefficient
  allocate(umwm%ustar  (msea)) ! air-side friction velocity
  allocate(umwm%taux   (msea)) !
  allocate(umwm%tauy   (msea)) !
  allocate(umwm%iactive(msea)) ! wave model activated the previous timestep

  ! 1-d arrays:
  allocate(cth(pm),sth(pm),th(pm))

  allocate(cthp_dirv(pm,mvsfc))
  allocate(sthp_dirv(pm,mvsfc))

  allocate(oc(msea))

  ! 2-d arrays (remapped):
! allocate(fcutoff(msea))   ! cutoff frequency
  allocate(swh(msea))       ! significant wave height
  allocate(mss(msea))       ! mean-squared slope

! allocate(shelt(msea))     ! sheltering coefficient
! shelt = 0

  ! mean spectrum quantities:
  allocate(mwd(msea)) ! direction
  allocate(mwp(msea)) ! period
  allocate(mwl(msea)) ! wavelength

  allocate(opeak(msea))
  allocate(ppeak(msea))

  ! dominant spectrum quantities:
  allocate(dwd (msea)) !; dwd  = 0 ! direction
  allocate(dwp (msea)) !; dwp  = 0 ! period
  allocate(dwl (msea)) !; dwl  = 0 ! wavelength
  allocate(dcp0(msea)) !; dcp0 = 0 ! intrinsic phase speed
  allocate(dcp (msea)) !; dcp  = 0 ! phase speed
  allocate(dcg0(msea)) !; dcg0 = 0 ! intrinsic group speed
  allocate(dcg (msea)) !; dcg  = 0 ! group speed

! allocate(momx (msea))
! allocate(momy (msea)) ! total wave momentum
! allocate(cgmxx(msea))
! allocate(cgmxy(msea))
! allocate(cgmyy(msea)) ! cg*m

  ! momentum fluxes:
  allocate(taux_form  (msea)) ! taux_form   = 0.
  allocate(tauy_form  (msea)) ! tauy_form   = 0.
  allocate(taux_skin  (msea)) ! taux_skin   = 0.
  allocate(tauy_skin  (msea)) ! tauy_skin   = 0.
! allocate(taux_diag  (msea)) ; taux_diag   = 0.
! allocate(tauy_diag  (msea)) ; tauy_diag   = 0.
! allocate(taux_ocntop(msea)) ; taux_ocntop = 0.
! allocate(tauy_ocntop(msea)) ; tauy_ocntop = 0.
! allocate(taux_ocnbot(msea)) ; taux_ocnbot = 0.
! allocate(tauy_ocnbot(msea)) ; tauy_ocnbot = 0.
! allocate(taux_snl   (msea)) ; taux_snl    = 0.
! allocate(tauy_snl   (msea)) ; tauy_snl    = 0.
! allocate(epsx_atm   (msea)) ; epsx_atm    = 0.
! allocate(epsy_atm   (msea)) ; epsy_atm    = 0.
! allocate(epsx_ocn   (msea)) ; epsx_ocn    = 0.
! allocate(epsy_ocn   (msea)) ; epsy_ocn    = 0.
! allocate(taux1      (msea)) ; taux1       = 0.
! allocate(tauy1      (msea)) ; tauy1       = 0.
! allocate(taux2      (msea)) ; taux2       = 0.
! allocate(tauy2      (msea)) ; tauy2       = 0.
! allocate(taux3      (msea)) ; taux3       = 0.
! allocate(tauy3      (msea)) ; tauy3       = 0.

  allocate(ucurr(msea)) !; ucurr = 0 ! ocean currents
  allocate(vcurr(msea)) !; vcurr = 0 ! ocean currents

  allocate(bf1_renorm(om,msea))
  allocate(bf2_renorm(om,msea))
  allocate(    cothkd(om,msea))
  allocate(       dwn(om,msea))
  allocate(     fkovg(om,msea))
  allocate(    invcp0(om,msea))
  allocate(        k4(om,msea))
  allocate(       kdk(om,msea))
  allocate(      k3dk(om,msea))
  allocate( oneoverk4(om,msea))
! allocate(       sbf(om,msea))
! allocate(       sdt(om,msea))
  allocate(       sdv(om,msea))
! allocate(   snl_arg(om,msea))
  allocate(        l2(om,msea))
  allocate(logl2overz(om,msea))
  allocate(       cp0(om,msea))
  allocate(       cg0(om,msea))
  allocate(         k(om,msea))

  allocate(    freq  (om))
  allocate(    fsdsf (om))

! allocate(      sice(om,msea)) ! wave attenuation by sea ice function

  allocate(  evs(pm,om,msea)) ! wave variance spectrum
  allocate( evsf(pm,om,msea)) ! wave variance spectrum, forward in time
! allocate( rotl(pm,om,msea)) ! anti-clockwise rotation, used in refraction
! allocate( rotr(pm,om,msea)) ! clockwise rotation, used in refraction
! allocate(  sds(pm,om,msea)) ! wave dissipation sink function
! allocate(  snl(pm,om,msea)) ! wave downshifting source/sink function
  allocate( ssin(pm,om,msea)) ! wind input source/sink function
! allocate(dummy(pm,om,msea)) ! dummy array, used in sds, advection and refraction

  ! Initialize:

  !$omp parallel do
  do i = 1, msea
     evs     (:,:,i) = 0.0
     ssin    (:,:,i) = 0.0

     umwm%iactive(i) = .false.

     swh  (i) = 0.0
     mss  (i) = 0.0
     mwd  (i) = 0.0
     mwp  (i) = 0.0
     mwl  (i) = 0.0
     dwd  (i) = 0.0
     dwp  (i) = 0.0
     dwl  (i) = 0.0
     dcp0 (i) = 0.0
     dcp  (i) = 0.0
     dcg0 (i) = 0.0
     dcg  (i) = 0.0
     ucurr(i) = 0.0
     vcurr(i) = 0.0

  enddo
  !$omp end parallel do

end subroutine alloc_umwm

!=========================================================================

subroutine filltab_umwm()

  use var_tables, only: increment_vtable

  implicit none

  if (allocated(evs))          call increment_vtable('EVS',          'SW', rvar3=evs)
  if (allocated(umwm%ustar))   call increment_vtable('UMWM%USTAR',   'SW', rvar1=umwm%ustar)
  if (allocated(umwm%alogzs))  call increment_vtable('UMWM%ALOGZS',  'SW', rvar1=umwm%alogzs)
  if (allocated(umwm%alogzs))  call increment_vtable('UMWM%ALOGZO',  'SW', rvar1=umwm%alogzo)
  if (allocated(umwm%iactive)) call increment_vtable('UMWM%IACTIVE', 'SW', lvar1=umwm%iactive)

end subroutine filltab_umwm

end module umwm_module
