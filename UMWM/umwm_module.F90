module umwm_module

  implicit none

  real, parameter :: euler    = exp(1.)    ! e
  real, parameter :: eulerinv = 1. / euler ! 1/e
  real, parameter :: rhosw    = 1025.      ! Avg density of seawater at shallow depths [kg/m^3]
  real, parameter :: rhoswi   = 1. / rhosw ! Inverse of rhosw [m^3/kg]

  ! VARIABLES INPUT FROM UMWM NAMELIST (NOW FILLED IN SUBROUTINE NMLASSIGN)

  integer :: umwmflg ! main flag to use UMWM (0 = no; 1 = yes)
  integer :: om      ! number of frequency/wavenumber bins
  integer :: pm      ! number of direction bins
  logical :: stokes  ! Output switch for Stokes drift velocity fields

  real :: fmin      ! Frequency of lowest frequency bin [Hz]
  real :: fmax      ! Frequency of highest frequency bin [Hz]
  real :: fprog     ! Highest prognostic frequency bin [Hz]
  real :: nu_air    ! Kinematic viscosity of air [m^2/s]
  real :: nu_water  ! Kinematic viscosity of water [m^2/s]
  real :: sfct      ! Surface tension [N/m]
  real :: gustiness ! Random wind gustiness factor (should be between 0 and 0.2)
  real :: dmin      ! Minimum allowed ocean depth [m]
  real :: explim    ! Exponent limiter (0.69 ~ 100% growth)
  real :: sin_fac   ! Input factor from following winds
  real :: sin_diss1 ! Damping factor from opposing winds
  real :: sin_diss2 ! Damping factor from swell overrunning wind
  real :: sds_fac   ! Breaking dissipation factor
  real :: sds_power ! Saturation spectrum power
  real :: mss_fac   ! Mean-square-slope adjustment to Sds
  real :: snl_fac   ! Wave energy downshifting factor
  real :: sdt_fac   ! Dissipation due to turbulence factor
  real :: sbf_fac   ! Bottom friction coefficient [m/s]
  real :: sbp_fac   ! Bottom percolation coefficient [m/s]
  real :: fice_lth  ! Sea ice fraction - lower threshold for attenuation
  real :: fice_uth  ! Sea ice fraction - upper threshold for attenuation

  ! OTHER VARIABLES

  real :: dt_olam       ! Time step from OLAM [s]
  real :: dta           ! Time step for advection and updating [s]
  real :: dtr           ! Time step for rotation (for now, set equal to dta) [s]
  real :: dth           ! Increment between discrete wave directions [rad]
  real :: oneovdth      ! 1 / dth
  real :: dthg          ! dth * grav
  real :: dlnf          ! (log(fmax)-log(fmin))/real(om-1)
  real :: bf1           ! exp(-16*dlnf*dlnf) / [exp(-16*dlnf*dlnf) + exp(-64*dlnf*dlnf)]
  real :: bf2           ! exp(-64*dlnf*dlnf) / [exp(-16*dlnf*dlnf) + exp(-64*dlnf*dlnf)]
  real :: fieldscale1   ! sin_diss1 / sin_fac [ ]
  real :: fieldscale2   ! sin_diss2 / sin_diss1 [ ]
  real :: inv_sds_power ! 1 / sds_power
  real :: twopisds_fac  ! 2 * pi * sds_fac

  ! The following arrays are included in umwm data structure to prevent any 
  ! conflicts or confusion with similar quantities or variable names in OLAM 

                                         !I = assigned at initialization only
                                         !P = prognosed or updated as part of prognostic quantitites
                                         !D = diagnostic only for output; not used for prognosing
  Type umwm_vars
     real, allocatable :: f          (:) !I (om)   Wave frequency of frequency bins ([Hz]
     real, allocatable :: seadep     (:) !I (msea) sea depth [m]
     real, allocatable :: uwind      (:) !P (msea) eastward wind velocity [m/s]
     real, allocatable :: vwind      (:) !P (msea) northward wind velocity [m/s]
     real, allocatable :: wspd       (:) !P (msea) wind speed [m/s]
     real, allocatable :: wdir       (:) !P (msea) wind direction [rad, -pi to pi]
     real, allocatable :: cd         (:) !Psr (msea) drag coefficient of air [ ]
     real, allocatable :: ustar      (:) !Psr (msea) friction velocity of air [m/s]
     real, allocatable :: taux       (:) !Psr (msea) stress (momentum flux) [n/m^2]
     real, allocatable :: tauy       (:) !Psr (msea) stress (momentum flux) [n/m^2]
  End type

  type (umwm_vars)     :: umwm

  integer, allocatable :: oc      (:) !I (msea) cut-off frequency index (maximum prognostic)
  integer, allocatable :: pl      (:) !I (pm)   anti-clockwise directional index
  integer, allocatable :: pr      (:) !I (pm)   clockwise directional index

  real, allocatable :: th         (:) !I (pm)   Wave direction [rad; -pi to pi]
  real, allocatable :: cth        (:) !I (pm)   Cosine of wave direction th [ ]
  real, allocatable :: sth        (:) !I (pm)   Sine of wave direction th [ ]
  real, allocatable :: cth2       (:) !I (pm)   cos^2(dth*(p-1))**2 [ ]

  real, allocatable :: fcutoff    (:) !P (msea) cutoff frequency [Hz] = 0.53 * grav / wspd(msea)
  real, allocatable :: shelt      (:) !P (msea) sheltering coefficient [ ]
  real, allocatable :: ucurr      (:) !P (msea) eastward surface ocean current [m/s]
  real, allocatable :: vcurr      (:) !P (msea) northward surface ocean current [m/s]
  real, allocatable :: taux_form  (:) !Psr (msea) form drag from atmosphere, x-component [N/m^2]
  real, allocatable :: tauy_form  (:) !Psr (msea) form drag from atmosphere, y-component [N/m^2]
  real, allocatable :: taux_skin  (:) !Psr (msea) skin-induced stress, x-component [N/m^2]
  real, allocatable :: tauy_skin  (:) !Psr (msea) skin-induced stress, y-component [N/m^2]
  real, allocatable :: tailatmx   (:) !Psr (msea) atmosphere tail stress part, x-component [N/m^2]
  real, allocatable :: tailatmy   (:) !Psr (msea) atmosphere tail stress part, y-component [N/m^2]

  real, allocatable :: tailocnx   (:) !Dsr (msea) ocean tail stress part, x-component [N/m^2]
  real, allocatable :: tailocny   (:) !Dsr (msea) ocean tail stress part, y-component [N/m^2]
  real, allocatable :: taux_diag  (:) !Dsr (msea) diagnostic form drag, x-component [N/m^2]
  real, allocatable :: tauy_diag  (:) !Dsr (msea) diagnostic form drag, y-component [N/m^2]
  real, allocatable :: taux_ocntop(:) !Dsr (msea) breaking-wave stress into ocean top, x-component [N/m^2]
  real, allocatable :: tauy_ocntop(:) !Dsr (msea) breaking-wave stress into ocean top, y-component [N/m^2]
  real, allocatable :: taux_ocnbot(:) !Dsr (msea) stress from waves into ocean bottom, x-component [N/m^2]
  real, allocatable :: tauy_ocnbot(:) !Dsr (msea) stress from waves into ocean bottom, y-component [N/m^2]
  real, allocatable :: taux_snl   (:) !Dsr (msea) momentum loss due to snl, x-component [N/m^2]
  real, allocatable :: tauy_snl   (:) !Dsr (msea) momentum loss due to snl, y-component [N/m^2]
  real, allocatable :: taux1      (:) !Dsr (msea) form drag x-component: wind pushing waves [N/m^2]
  real, allocatable :: tauy1      (:) !Dsr (msea) form drag y-component: wind pushing waves [N/m^2]
  real, allocatable :: taux2      (:) !Dsr (msea) form drag x-component: waves against wind [N/m^2]
  real, allocatable :: tauy2      (:) !Dsr (msea) form drag y-component: waves against wind [N/m^2]
  real, allocatable :: taux3      (:) !Dsr (msea) form drag x-component: waves overrunning wind [N/m^2]
  real, allocatable :: tauy3      (:) !Dsr (msea) form drag y-component: waves overrunning wind [N/m^2]
  real, allocatable :: epsx_atm   (:) !Dsr (msea) wave energy growth flux, x-component [kg/s^3]
  real, allocatable :: epsy_atm   (:) !Dsr (msea) wave energy growth flux, y-component [kg/s^3]
  real, allocatable :: epsx_ocn   (:) !Dsr (msea) wave energy dissipation flux, x-component [kg/s^3]
  real, allocatable :: epsy_ocn   (:) !Dsr (msea) wave energy dissipation flux, y-component [kg/s^3]
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
  real, allocatable :: momx       (:) !Ddi (msea) total wave momentum, x-component [kg m/s]
  real, allocatable :: momy       (:) !Ddi (msea) total wave momentum, y-component [kg m/s]
  real, allocatable :: cgmxx      (:) !Ddi (msea) cg*momentum, xx-component [kg m^2/s^2]
  real, allocatable :: cgmxy      (:) !Ddi (msea) cg*momentum, xy-component [kg m^2/s^2]
  real, allocatable :: cgmyy      (:) !Ddi (msea) cg*momentum, yy-component [kg m^2/s^2]

  real, allocatable :: cthp_dirv (:,:) !I (pm,mvsfc) wave ray directions projected onto sfcg%dirv [ ]
  real, allocatable :: sthp_dirv (:,:) !I (pm,mvsfc) wave ray directions projected onto negative sfcg%diru [ ]
  real, allocatable :: cth2pp    (:,:) !I (pm,pm)    utility array used for mss in sds routine:
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
  real, allocatable :: sbf       (:,:) !I (om,msea)  bottom friction from roughness + percolation
  real, allocatable :: sdv       (:,:) !I (om,msea)  dissipation due to viscosity
  real, allocatable :: snl_arg   (:,:) !I (om,msea)  exp function argument for snl

  real, allocatable :: sdt       (:,:) !P (om,msea) dissipation due to turbulence
  real, allocatable :: sice      (:,:) !P (om,msea) wave attenuation by sea ice function [ ]

  real, allocatable :: evs  (:,:,:) !P (om,pm,msea) wave variance spectrum [m^2/(rad/m)]
  real, allocatable :: evsf (:,:,:) !P (om,pm,msea) wave variance spectrum, forward in time
  real, allocatable :: rotl (:,:,:) !P (om,pm,msea) anti-clockwise rotation, used in refraction
  real, allocatable :: rotr (:,:,:) !P (om,pm,msea) clockwise rotation, used in refraction
  real, allocatable :: sds  (:,:,:) !P (om,pm,msea) wave dissipation sink function
  real, allocatable :: snl  (:,:,:) !P (om,pm,msea) wave downshifting source/sink function
  real, allocatable :: ssin (:,:,:) !P (om,pm,msea) wind input function
  real, allocatable :: dummy(:,:,:) !P (om,pm,msea) dummy array, used in sds, advection and refraction

contains

!===============================================================================

subroutine alloc_umwm()

  use mem_sea,  only: msea
  use mem_sfcg, only: mvsfc

  ! Allocates UMWM arrays

  allocate(umwm%f     (om))   !
  allocate(umwm%seadep(msea)) ! depth
  allocate(umwm%uwind (msea)) ! wind component
  allocate(umwm%vwind (msea)) ! wind component
  allocate(umwm%wspd  (msea)) ! wind speed
  allocate(umwm%wdir  (msea)) ! wind direction
  allocate(umwm%cd    (msea)) ! air-side drag coefficient
  allocate(umwm%ustar (msea)) ! air-side friction velocity
  allocate(umwm%taux  (msea)) !
  allocate(umwm%tauy  (msea)) !

  ! 1-d arrays:
  allocate(cth(pm),cth2(pm),pl(pm),pr(pm),sth(pm),th(pm))

  allocate(cthp_dirv(pm,mvsfc))
  allocate(sthp_dirv(pm,mvsfc))

  allocate(oc(msea))

  ! 2-d arrays (remapped):
  allocate(fcutoff(msea))   ! cutoff frequency
  allocate(swh(msea))       ! significant wave height
  allocate(mss(msea))       ! mean-squared slope

  allocate(shelt(msea))     ! sheltering coefficient
  shelt = 0

  ! mean spectrum quantities:
  allocate(mwd(msea)) ! direction
  allocate(mwp(msea)) ! period
  allocate(mwl(msea)) ! wavelength

  ! dominant spectrum quantities:
  allocate(dwd (msea)) ; dwd  = 0 ! direction
  allocate(dwp (msea)) ; dwp  = 0 ! period
  allocate(dwl (msea)) ; dwl  = 0 ! wavelength
  allocate(dcp0(msea)) ; dcp0 = 0 ! intrinsic phase speed
  allocate(dcp (msea)) ; dcp  = 0 ! phase speed
  allocate(dcg0(msea)) ; dcg0 = 0 ! intrinsic group speed
  allocate(dcg (msea)) ; dcg  = 0 ! group speed

  allocate(momx (msea))
  allocate(momy (msea)) ! total wave momentum
  allocate(cgmxx(msea))
  allocate(cgmxy(msea))
  allocate(cgmyy(msea)) ! cg*m

  ! momentum fluxes:
  allocate(taux_form  (msea)) ; taux_form   = 0.
  allocate(tauy_form  (msea)) ; tauy_form   = 0.
  allocate(taux_skin  (msea)) ; taux_skin   = 0.
  allocate(tauy_skin  (msea)) ; tauy_skin   = 0.
  allocate(taux_diag  (msea)) ; taux_diag   = 0.
  allocate(tauy_diag  (msea)) ; tauy_diag   = 0.
  allocate(taux_ocntop(msea)) ; taux_ocntop = 0.
  allocate(tauy_ocntop(msea)) ; tauy_ocntop = 0.
  allocate(taux_ocnbot(msea)) ; taux_ocnbot = 0.
  allocate(tauy_ocnbot(msea)) ; tauy_ocnbot = 0.
  allocate(taux_snl   (msea)) ; taux_snl    = 0.
  allocate(tauy_snl   (msea)) ; tauy_snl    = 0.
  allocate(epsx_atm   (msea)) ; epsx_atm    = 0.
  allocate(epsy_atm   (msea)) ; epsy_atm    = 0.
  allocate(epsx_ocn   (msea)) ; epsx_ocn    = 0.
  allocate(epsy_ocn   (msea)) ; epsy_ocn    = 0.
  allocate(taux1      (msea)) ; taux1       = 0.
  allocate(tauy1      (msea)) ; tauy1       = 0.
  allocate(taux2      (msea)) ; taux2       = 0.
  allocate(tauy2      (msea)) ; tauy2       = 0.
  allocate(taux3      (msea)) ; taux3       = 0.
  allocate(tauy3      (msea)) ; tauy3       = 0.

  allocate(tailatmx(msea),tailatmy(msea))
  allocate(tailocnx(msea),tailocny(msea))

  allocate(ucurr(msea)) ; ucurr = 0 ! ocean currents
  allocate(vcurr(msea)) ; vcurr = 0 ! ocean currents

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
  allocate(       sbf(om,msea))
  allocate(       sdt(om,msea))
  allocate(       sdv(om,msea))
  allocate(   snl_arg(om,msea))
  allocate(        l2(om,msea))
  allocate(logl2overz(om,msea))
  allocate(      sice(om,msea)) ! wave attenuation by sea ice function
  allocate(       cp0(om,msea))
  allocate(       cg0(om,msea))
  allocate(         k(om,msea))

  allocate(  evs(om,pm,msea)) ! wave variance spectrum
  allocate( evsf(om,pm,msea)) ! wave variance spectrum, forward in time
  allocate( rotl(om,pm,msea)) ! anti-clockwise rotation, used in refraction
  allocate( rotr(om,pm,msea)) ! clockwise rotation, used in refraction
  allocate(  sds(om,pm,msea)) ! wave dissipation sink function
  allocate(  snl(om,pm,msea)) ! wave downshifting source/sink function
  allocate( ssin(om,pm,msea)) ! wind input source/sink function
  allocate(dummy(om,pm,msea)) ! dummy array, used in sds, advection and refraction

  ! initialize:
  evs = tiny(evs)

end subroutine alloc_umwm

!=========================================================================

subroutine filltab_umwm()

  use var_tables, only: increment_vtable

  implicit none

  if (allocated(evs))        call increment_vtable('EVS',        'SW', rvar3=evs)
  if (allocated(umwm%ustar)) call increment_vtable('UMWM%USTAR', 'SW', rvar1=umwm%ustar)

end subroutine filltab_umwm

end module umwm_module
