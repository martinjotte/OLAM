subroutine simple_physics (mza, ka, za, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test, iw)
!----------------------------------------------------------------------- 
! 
! Purpose: Simple Physics Package
!
! Author: K. A. Reed (University of Michigan, kareed@umich.edu)
!         version 5 
!         July/8/2012
! 
!  Change log:
!  v2: removal of some NCAR CAM-specific 'use' associations
!  v3: corrected precl(i) computation, the precipitation rate is now computed via a vertical integral, the previous single-level computation in v2 was a bug
!  v3: corrected dtdt(i,1) computation, the term '-(i,1)' was missing the temperature variable: '-t(i,1)'
!  v4: modified and enhanced parameter list to make the routine truly standalone, the number of columns and vertical levels have been added: pcols, pver
!  v4: 'ncol' has been removed, 'pcols' is used instead
!  v5: the sea surface temperature (SST) field Tsurf is now an array, the SST now depends on the latitude
!  v5: addition of the latitude array 'lat' and the flag 'test' in the parameter list
!      if test = 0: constant SST is used, correct setting for the tropical cyclone test case 5-1
!      if test = 1: newly added latitude-dependent SST is used, correct setting for the moist baroclinic wave test with simple-physics (test 4-3)
! 
! Description: Includes large-scale precipitation, surface fluxes and
!              boundary-leyer mixing. The processes are time-split
!              in that order. A partially implicit formulation is
!              used to foster numerical stability.
!              The routine assumes that the model levels are ordered
!              in a top-down approach, e.g. level 1 denotes the uppermost
!              full model level.
!
!              This routine is based on an implementation which was
!              developed for the NCAR Community Atmosphere Model (CAM).
!              Adjustments for other models will be necessary.
!
!              The routine provides both updates of the state variables
!              u, v, T, q (these are local copies of u,v,T,q within this physics
!              routine) and also collects their time tendencies.
!              The latter might be used to couple the physics and dynamics
!              in a process-split way. For a time-split coupling, the final
!              state should be given to the dynamical core for the next time step.
! Test:      0 = Reed and Jablonowski (2011) tropical cyclone test case (test 5-1)
!            1 = Moist baroclinic instability test (test 4-3)
!========================================================================================
!            2 = Moist baroclinic instability test (test 4-2)    CODE CHANGE BY BOB WALKO
!========================================================================================
!
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!-----------------------------------------------------------------------
  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

   implicit none

   integer, parameter :: r8 = selected_real_kind(12)

! Input arguments - MODEL DEPENDENT

   integer, intent(in) :: iw           ! OLAM horizontal index

   integer, intent(in) :: mza  ! number of vertical levels
   integer, intent(in) :: ka   ! lowest predicted level (min value is 2)

   real(r8), intent(in) :: za           ! height of lowest predicted model level
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   real(r8), intent(in) :: lat          ! Latitude 
   integer, intent(in)  :: test         ! Test number

! Input/Output arguments 

   real(r8), intent(inout) :: t(mza)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(mza)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(mza)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(mza)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(in) :: pmid(mza)   ! Pressure is full-model level (Pa)
   real(r8), intent(in) :: pint(mza) ! Pressure at model interfaces (Pa)
   real(r8), intent(in) :: pdel(mza)   ! Layer thickness (Pa)
   real(r8), intent(in) :: rpdel(mza)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(in) :: ps           ! Surface Pressue (Pa) [assumed at z = zm(ka-1)]

! Output arguments 

   real(r8), intent(out) :: precl         ! Precipitation rate (m_water / s)

!---------------------------Local workspace-----------------------------

! Integers for loops

   integer  k                           ! level index

! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Simple Physics Specific Constants 

!++++++++                     
   real(r8) Tsurf                       ! Sea Surface Temperature (constant for tropical cyclone)
!++++++++                                 Tsurf needs to be dependent on latitude for the
                                        ! moist baroclinic wave test 4-3 with simple-physics, adjust
   real(r8) SST_tc                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test

! Physics Tendency Arrays
  real(r8) dtdt(mza)             ! Temperature tendency 
  real(r8) dqdt(mza)             ! Specific humidity tendency
  real(r8) dudt(mza)             ! Zonal wind tendency
  real(r8) dvdt(mza)             ! Meridional wind tendency

! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

! Variables for Boundary Layer Calculation

   real(r8) wind                 ! Magnitude of Wind
   real(r8) Cd                   ! Drag coefficient for momentum
   real(r8) Km(mza)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(mza)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(mza)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(mza)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(mza)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(mza)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(mza)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(mza)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(mza)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(mza)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(mza)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(mza)           ! Matrix Coefficents for PBL Scheme

! Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
! conserve the mass of the dry air

   real(r8) qini(mza)            ! Initial specific humidity

!===============================================================================
!
! Definition of array dimensions - MODEL DEPENDENT
!
!===============================================================================

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Smith and Vogl 2008
      SST_tc   = 302.15_r8      ! Constant Value for SST for tropical cyclone test
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021            ! Maximum specific humidity for baro test

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================

! Calculate hydrostatic height za of the lowest model level

!     dlnpint = log(ps) - log(pint(ka))                 ! ps is identical to pint(ka-1), note: this is the correct sign (corrects typo in JAMES paper)
!     za = rair/gravit*t(ka)*(1._r8+zvir*q(ka))*0.5_r8*dlnpint

! Set Initial Specific Humidity - For dry mass adjustment at the end

     qini(:mza) = q(:mza)

!--------------------------------------------------------------
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics 
!--------------------------------------------------------------

     if (test .eq. 1) then     ! moist baroclinic wave with simple-physics
        Tsurf = (T00 + pi*u0/rair * 1.5_r8 * sin(etav) * (cos(etav))**0.5_r8 *                 &
                ((-2._r8*(sin(lat))**6 * ((cos(lat))**2 + 1._r8/3._r8) + 10._r8/63._r8)* &
                u0 * (cos(etav))**1.5_r8  +                                                    &
                (8._r8/5._r8*(cos(lat))**3 * ((sin(lat))**2 + 2._r8/3._r8) - pi/4._r8)*a*omega*0.5_r8 ))/ &
                (1._r8+zvir*q0*exp(-(lat/latw)**4))

     else
        Tsurf = SST_tc
     end if


!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:mza)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:mza)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:mza)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:mza)  = 0._r8            ! initialize meridional wind tendency with zero
     precl = 0._r8                  ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Precipitation Rate
!
!===============================================================================

! Calculate Tendencies

      do k=ka,mza  ! Loop order does not matter
         qsat = epsilo*e0/pmid(k)*exp(-latvap/rh2o*((1._r8/t(k))-1._r8/T0))  ! saturation specific humidity

if (iw == 46685 .and. k < 15) then
   write(6,'(a,i5,10f12.5)') 'sp1 ',k,qsat,q(k),q(k)-qsat
endif
         if (q(k) > qsat) then                                                 ! saturated?
            tmp  = 1._r8/dtime*(q(k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k)**2)))
            dtdt(k) = dtdt(k)+latvap/cpair*tmp
            dqdt(k) = dqdt(k)-tmp
            precl = precl + tmp*pdel(k)/(gravit*rhow)                    ! precipitation rate, computed via a vertical integral
                                                                                    ! corrected in version 1.3

if (iw == 46685 .and. k < 15) then
   write(6,'(a,i5,10f12.5)') 'sp2 ',k,qsat,q(k),dtdt(k),latvap/cpair*tmp,dqdt(k),dtdt(k)*dtime
endif


         end if
      end do

! Update moisture and temperature fields from Large-Scale Precipitation Scheme

      do k=ka,mza
         t(k) =  t(k) + dtdt(k)*dtime    ! update the state variables T and q
         q(k) =  q(k) + dqdt(k)*dtime
      end do

!=================================================================
if (test == 2) return   ! Code added by Bob Walko for test_case 42
!=================================================================

!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! We are using Simplified Ekman theory to compute the diffusion coefficients 
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================

! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
! they depend on the conditions at the lowest model level and stay constant
! up to the 850 hPa level. Above this level the coefficients are decreased
! and tapered to zero. At the 700 hPa level the strength of the K coefficients
! is about 10% of the maximum strength. 

      wind = sqrt(u(ka)**2+v(ka)**2)    ! wind magnitude at the lowest level
      Ke(ka-1) = C*wind*za
      if( wind .lt. v20) then
         Cd = Cd0+Cd1*wind 
         Km(ka-1) = Cd*wind*za
      else
         Cd = Cm
         Km(ka-1) = Cm*wind*za
      endif

      do k=ka,mza
         if( pint(k) .ge. pbltop) then
            Km(k) = Km(ka-1)                 ! constant Km below 850 hPa level
            Ke(k) = Ke(ka-1)                 ! constant Ke below 850 hPa level
         else
            Km(k) = Km(ka-1)*exp(-(pbltop-pint(k))**2/(pblconst)**2)  ! Km tapered to 0
            Ke(k) = Ke(ka-1)*exp(-(pbltop-pint(k))**2/(pblconst)**2)  ! Ke tapered to 0
         end if 
      end do     

!===============================================================================
! Update the state variables u, v, t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (JAMES, 2012)
!
! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics, adjust
!===============================================================================

     qsats = epsilo*e0/ps*exp(-latvap/rh2o*((1._r8/Tsurf)-1._r8/T0))  ! saturation specific humidity at the surface
     dudt(ka) = dudt(ka) + (u(ka) &
                         /(1._r8+Cd*wind*dtime/za)-u(ka))/dtime
     dvdt(ka) = dvdt(ka) + (v(ka) &
                         /(1._r8+Cd*wind*dtime/za)-v(ka))/dtime
     u(ka)   = u(ka)/(1._r8+Cd*wind*dtime/za)
     v(ka)   = v(ka)/(1._r8+Cd*wind*dtime/za)
     dtdt(ka) = dtdt(ka) +((t(ka)+C*wind*Tsurf*dtime/za) &
                         /(1._r8+C*wind*dtime/za)-t(ka))/dtime 
     t(ka)   = (t(ka)+C*wind*Tsurf*dtime/za) &
                         /(1._r8+C*wind*dtime/za)  
     dqdt(ka) = dqdt(ka) +((q(ka)+C*wind*qsats*dtime/za) &
                         /(1._r8+C*wind*dtime/za)-q(ka))/dtime
     q(ka) = (q(ka)+C*wind*qsats*dtime/za)/(1._r8+C*wind*dtime/za)
!===============================================================================

!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme

      do k=ka+1,mza  ! loop order does not matter
         rho = (pint(k-1)/(rair*(t(k-1)+t(k))/2.0_r8))
         CAm(k)   = rpdel(k)*dtime*gravit*gravit*Km(k-1)*rho*rho   &
                      /(pmid(k-1)-pmid(k))    
         CCm(k-1) = rpdel(k-1)*dtime*gravit*gravit*Km(k-1)*rho*rho &
                      /(pmid(k-1)-pmid(k))
         CA(k)    = rpdel(k)*dtime*gravit*gravit*Ke(k-1)*rho*rho   &
                      /(pmid(k-1)-pmid(k))
         CC(k-1)  = rpdel(k-1)*dtime*gravit*gravit*Ke(k-1)*rho*rho &
                      /(pmid(k-1)-pmid(k))
      end do
      CAm(ka) = 0._r8
      CCm(mza) = 0._r8
      CEm(ka-1) = 0._r8
      CA(ka) = 0._r8
      CC(mza) = 0._r8
      CE(ka-1) = 0._r8
      CFu(ka-1) = 0._r8
      CFv(ka-1) = 0._r8
      CFt(ka-1) = 0._r8
      CFq(ka-1) = 0._r8 
      do k=ka,mza  ! loop order matters
         CE(k)  = CC(k)/(1._r8+CA(k)+CC(k)-CA(k)*CE(k-1)) 
         CEm(k) = CCm(k)/(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k-1))
         CFu(k) = (u(k)+CAm(k)*CFu(k-1)) &
                    /(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k-1))
         CFv(k) = (v(k)+CAm(k)*CFv(k-1)) &
                    /(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k-1))
         CFt(k) = ((p0/pmid(k))**(rair/cpair)*t(k)+CA(k)*CFt(k-1)) &
                    /(1._r8+CA(k)+CC(k)-CA(k)*CE(k-1)) 
         CFq(k) = (q(k)+CA(k)*CFq(k-1)) &
                    /(1._r8+CA(k)+CC(k)-CA(k)*CE(k-1))
      end do

! Calculate the updated temperature, specific humidity and horizontal wind

! First we need to calculate the updates at the top model level

      dudt(mza)  = dudt(mza)+(CFu(mza)-u(mza))/dtime
      dvdt(mza)  = dvdt(mza)+(CFv(mza)-v(mza))/dtime
      u(mza)    = CFu(mza)
      v(mza)    = CFv(mza)
      dtdt(mza)  = dtdt(mza)+(CFt(mza)*(pmid(mza)/p0)**(rair/cpair)-t(mza))/dtime  ! corrected in version 1.3
      t(mza)    = CFt(mza)*(pmid(mza)/p0)**(rair/cpair)
      dqdt(mza)  = dqdt(mza)+(CFq(mza)-q(mza))/dtime
      q(mza)  = CFq(mza)

! Loop over the remaining levels

      do k=mza-1,ka,-1  ! loop order matters
         dudt(k)  = dudt(k)+(CEm(k)*u(k+1)+CFu(k)-u(k))/dtime
         dvdt(k)  = dvdt(k)+(CEm(k)*v(k+1)+CFv(k)-v(k))/dtime
         u(k)    = CEm(k)*u(k+1)+CFu(k) 
         v(k)    = CEm(k)*v(k+1)+CFv(k)
         dtdt(k)  = dtdt(k)+((CE(k)*t(k+1) &
                           *(p0/pmid(k+1))**(rair/cpair)+CFt(k)) &
                           *(pmid(k)/p0)**(rair/cpair)-t(k))/dtime 
         t(k)    = (CE(k)*t(k+1)*(p0/pmid(k+1))**(rair/cpair)+CFt(k)) &
                           *(pmid(k)/p0)**(rair/cpair)
         dqdt(k)  = dqdt(k)+(CE(k)*q(k+1)+CFq(k)-q(k))/dtime
         q(k)  = CE(k)*q(k+1)+CFq(k)

      end do

!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   return
end subroutine simple_physics 

!======================================================================================================

subroutine simple_physics0 (pver, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test, iw, &
                            dtdt, dqdt, dudt, dvdt)
!----------------------------------------------------------------------- 
! 
! Purpose: Simple Physics Package
!
! Author: K. A. Reed (University of Michigan, kareed@umich.edu)
!         version 5 
!         July/8/2012
!
!  Change log:
!  v2: removal of some NCAR CAM-specific 'use' associations
!  v3: corrected precl(i) computation, the precipitation rate is now computed via a vertical integral, the previous single-level computation in v2 was a bug
!  v3: corrected dtdt(i,1) computation, the term '-(i,1)' was missing the temperature variable: '-t(i,1)'
!  v4: modified and enhanced parameter list to make the routine truly standalone, the number of columns and vertical levels have been added: pcols, pver
!  v4: 'ncol' has been removed, 'pcols' is used instead
!  v5: the sea surface temperature (SST) field Tsurf is now an array, the SST now depends on the latitude
!  v5: addition of the latitude array 'lat' and the flag 'test' in the parameter list
!      if test = 0: constant SST is used, correct setting for the tropical cyclone test case 5-1
!      if test = 1: newly added latitude-dependent SST is used, correct setting for the moist baroclinic wave test with simple-physics (test 4-3)
! 
! Description: Includes large-scale precipitation, surface fluxes and
!              boundary-leyer mixing. The processes are time-split
!              in that order. A partially implicit formulation is
!              used to foster numerical stability.
!              The routine assumes that the model levels are ordered
!              in a top-down approach, e.g. level 1 denotes the uppermost
!              full model level.
!
!              This routine is based on an implementation which was
!              developed for the NCAR Community Atmosphere Model (CAM).
!              Adjustments for other models will be necessary.
!
!              The routine provides both updates of the state variables
!              u, v, T, q (these are local copies of u,v,T,q within this physics
!              routine) and also collects their time tendencies.
!              The latter might be used to couple the physics and dynamics
!              in a process-split way. For a time-split coupling, the final
!              state should be given to the dynamical core for the next time step.
! Test:      0 = Reed and Jablonowski (2011) tropical cyclone test case (test 5-1)
!            1 = Moist baroclinic instability test (test 4-3)
!
!
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!-----------------------------------------------------------------------
  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

   implicit none

   integer, parameter :: r8 = selected_real_kind(12)


! Input arguments - MODEL DEPENDENT

   integer, intent(in) :: iw           ! OLAM horizontal index

   integer, intent(in)  :: pver         ! Set number of model levels
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   real(r8), intent(in) :: lat   ! Latitude 
   integer, intent(in)  :: test         ! Test number
   

! Input/Output arguments 

!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere

   real(r8), intent(inout) :: t(pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(inout) :: pmid(pver)   ! Pressure is full-model level (Pa)
   real(r8), intent(inout) :: pint(pver+1) ! Pressure at model interfaces (Pa)
   real(r8), intent(inout) :: pdel(pver)   ! Layer thickness (Pa)
   real(r8), intent(inout) :: rpdel(pver)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(inout) :: ps          ! Surface Pressue (Pa)

! Physics Tendency Arrays
  real(r8), intent(out):: dtdt(pver)             ! Temperature tendency 
  real(r8), intent(out):: dqdt(pver)             ! Specific humidity tendency
  real(r8), intent(out):: dudt(pver)             ! Zonal wind tendency
  real(r8), intent(out):: dvdt(pver)             ! Meridional wind tendency



! Output arguments 

   real(r8), intent(out) :: precl         ! Precipitation rate (m_water / s)


!---------------------------Local workspace-----------------------------


! Integers for loops

   integer  i,k                         ! Longitude, level indices

! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Simple Physics Specific Constants 

!++++++++                     
   real(r8) Tsurf                ! Sea Surface Temperature (constant for tropical cyclone)
!++++++++                                 Tsurf needs to be dependent on latitude for the
                                        ! moist baroclinic wave test 4-3 with simple-physics, adjust

   real(r8) SST_tc                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test

! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

! Variables for Boundary Layer Calculation

   real(r8) wind                 ! Magnitude of Wind
   real(r8) Cd                   ! Drag coefficient for momentum
   real(r8) Km(pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za                   ! Heights at midpoints of first model level
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pver+1)           ! Matrix Coefficents for PBL Scheme


! Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
! conserve the mass of the dry air

   real(r8) qini(pver)            ! Initial specific humidity

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Smith and Vogl 2008
      SST_tc   = 302.15_r8      ! Constant Value for SST for tropical cyclone test
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021            ! Maximum specific humidity for baro test

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================

! Calculate hydrostatic height za of the lowest model level

        dlnpint = log(ps) - log(pint(pver))                 ! ps(i) is identical to pint(i,pver+1), note: this is the correct sign (corrects typo in JAMES paper)
        za = rair/gravit*t(pver)*(1._r8+zvir*q(pver))*0.5_r8*dlnpint

! Set Initial Specific Humidity - For dry mass adjustment at the end

     qini(:pver) = q(:pver)
!--------------------------------------------------------------
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics 
!--------------------------------------------------------------
     if (test .eq. 1) then     ! moist baroclinic wave with simple-physics
           Tsurf = (T00 + pi*u0/rair * 1.5_r8 * sin(etav) * (cos(etav))**0.5_r8 *                 &
                     ((-2._r8*(sin(lat))**6 * ((cos(lat))**2 + 1._r8/3._r8) + 10._r8/63._r8)* &
                     u0 * (cos(etav))**1.5_r8  +                                                    &
                     (8._r8/5._r8*(cos(lat))**3 * ((sin(lat))**2 + 2._r8/3._r8) - pi/4._r8)*a*omega*0.5_r8 ))/ &
                     (1._r8+zvir*q0*exp(-(lat/latw)**4))

     else
           Tsurf = SST_tc
     end if

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pver)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pver)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pver)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pver)  = 0._r8            ! initialize meridional wind tendency with zero
     precl = 0._r8                  ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Precipitation Rate
!
!===============================================================================

! Calculate Tendencies

      do k=1,pver
            qsat = epsilo*e0/pmid(k)*exp(-latvap/rh2o*((1._r8/t(k))-1._r8/T0))  ! saturation specific humidity
            if (q(k) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k)**2)))
               dtdt(k) = dtdt(k)+latvap/cpair*tmp
               dqdt(k) = dqdt(k)-tmp
               precl = precl + tmp*pdel(k)/(gravit*rhow)                    ! precipitation rate, computed via a vertical integral
                                                                                    ! corrected in version 1.3
            end if
      end do

! Update moisture and temperature fields from Large-Scale Precipitation Scheme

      do k=1,pver
            t(k) =  t(k) + dtdt(k)*dtime    ! update the state variables T and q
            q(k) =  q(k) + dqdt(k)*dtime
      end do

!=================================================================
if (test == 2) return   ! Code added by Bob Walko for test_case 42
!=================================================================

!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! We are using Simplified Ekman theory to compute the diffusion coefficients 
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================

! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
! they depend on the conditions at the lowest model level and stay constant
! up to the 850 hPa level. Above this level the coefficients are decreased
! and tapered to zero. At the 700 hPa level the strength of the K coefficients
! is about 10% of the maximum strength. 

        wind = sqrt(u(pver)**2+v(pver)**2)    ! wind magnitude at the lowest level
        Ke(pver+1) = C*wind*za
        if( wind .lt. v20) then
           Cd = Cd0+Cd1*wind 
           Km(pver+1) = Cd*wind*za
        else
           Cd = Cm
           Km(pver+1) = Cm*wind*za
        endif

      do k=1,pver
            if( pint(k) .ge. pbltop) then
               Km(k) = Km(pver+1)                 ! constant Km below 850 hPa level
               Ke(k) = Ke(pver+1)                 ! constant Ke below 850 hPa level
            else
               Km(k) = Km(pver+1)*exp(-(pbltop-pint(k))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(k) = Ke(pver+1)*exp(-(pbltop-pint(k))**2/(pblconst)**2)  ! Ke tapered to 0
            end if 
      end do     


!===============================================================================
! Update the state variables u, v, t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (JAMES, 2012)
!
! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics, adjust
!===============================================================================

        qsats = epsilo*e0/ps*exp(-latvap/rh2o*((1._r8/Tsurf)-1._r8/T0))  ! saturation specific humidity at the surface
        dudt(pver) = dudt(pver) + (u(pver) &
                            /(1._r8+Cd*wind*dtime/za)-u(pver))/dtime
        dvdt(pver) = dvdt(pver) + (v(pver) &
                            /(1._r8+Cd*wind*dtime/za)-v(pver))/dtime
        u(pver)   = u(pver)/(1._r8+Cd*wind*dtime/za)
        v(pver)   = v(pver)/(1._r8+Cd*wind*dtime/za)
        dtdt(pver) = dtdt(pver) +((t(pver)+C*wind*Tsurf*dtime/za) &
                            /(1._r8+C*wind*dtime/za)-t(pver))/dtime 
        t(pver)   = (t(pver)+C*wind*Tsurf*dtime/za) &
                            /(1._r8+C*wind*dtime/za)  
        dqdt(pver) = dqdt(pver) +((q(pver)+C*wind*qsats*dtime/za) &
                            /(1._r8+C*wind*dtime/za)-q(pver))/dtime
        q(pver) = (q(pver)+C*wind*qsats*dtime/za)/(1._r8+C*wind*dtime/za)
!===============================================================================


!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme

      do k=1,pver-1
            rho = (pint(k+1)/(rair*(t(k+1)+t(k))/2.0_r8))
            CAm(k)   = rpdel(k)*dtime*gravit*gravit*Km(k+1)*rho*rho   &
                         /(pmid(k+1)-pmid(k))    
            CCm(k+1) = rpdel(k+1)*dtime*gravit*gravit*Km(k+1)*rho*rho &
                         /(pmid(k+1)-pmid(k))
            CA(k)    = rpdel(k)*dtime*gravit*gravit*Ke(k+1)*rho*rho   &
                         /(pmid(k+1)-pmid(k))
            CC(k+1)  = rpdel(k+1)*dtime*gravit*gravit*Ke(k+1)*rho*rho &
                         /(pmid(k+1)-pmid(k))
      end do
         CAm(pver) = 0._r8
         CCm(1) = 0._r8
         CEm(pver+1) = 0._r8
         CA(pver) = 0._r8
         CC(1) = 0._r8
         CE(pver+1) = 0._r8
         CFu(pver+1) = 0._r8
         CFv(pver+1) = 0._r8
         CFt(pver+1) = 0._r8
         CFq(pver+1) = 0._r8 
         do k=pver,1,-1
            CE(k)  = CC(k)/(1._r8+CA(k)+CC(k)-CA(k)*CE(k+1)) 
            CEm(k) = CCm(k)/(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k+1))
            CFu(k) = (u(k)+CAm(k)*CFu(k+1)) &
                       /(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k+1))
            CFv(k) = (v(k)+CAm(k)*CFv(k+1)) &
                       /(1._r8+CAm(k)+CCm(k)-CAm(k)*CEm(k+1))
            CFt(k) = ((p0/pmid(k))**(rair/cpair)*t(k)+CA(k)*CFt(k+1)) &
                       /(1._r8+CA(k)+CC(k)-CA(k)*CE(k+1)) 
            CFq(k) = (q(k)+CA(k)*CFq(k+1)) &
                       /(1._r8+CA(k)+CC(k)-CA(k)*CE(k+1))
      end do


! Calculate the updated temperature, specific humidity and horizontal wind

! First we need to calculate the updates at the top model level

            dudt(1)  = dudt(1)+(CFu(1)-u(1))/dtime
            dvdt(1)  = dvdt(1)+(CFv(1)-v(1))/dtime
            u(1)    = CFu(1)
            v(1)    = CFv(1)
            dtdt(1)  = dtdt(1)+(CFt(1)*(pmid(1)/p0)**(rair/cpair)-t(1))/dtime  ! corrected in version 1.3
            t(1)    = CFt(1)*(pmid(1)/p0)**(rair/cpair)
            dqdt(1)  = dqdt(1)+(CFq(1)-q(1))/dtime
            q(1)  = CFq(1)

! Loop over the remaining level

         do k=2,pver
            dudt(k)  = dudt(k)+(CEm(k)*u(k-1)+CFu(k)-u(k))/dtime
            dvdt(k)  = dvdt(k)+(CEm(k)*v(k-1)+CFv(k)-v(k))/dtime
            u(k)    = CEm(k)*u(k-1)+CFu(k) 
            v(k)    = CEm(k)*v(k-1)+CFv(k)
            dtdt(k)  = dtdt(k)+((CE(k)*t(k-1) &
                              *(p0/pmid(k-1))**(rair/cpair)+CFt(k)) &
                              *(pmid(k)/p0)**(rair/cpair)-t(k))/dtime 
            t(k)    = (CE(k)*t(k-1)*(p0/pmid(k-1))**(rair/cpair)+CFt(k)) &
                              *(pmid(k)/p0)**(rair/cpair)
            dqdt(k)  = dqdt(k)+(CE(k)*q(k-1)+CFq(k)-q(k))/dtime
            q(k)  = CE(k)*q(k-1)+CFq(k)
         end do

!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   return
end subroutine simple_physics0

