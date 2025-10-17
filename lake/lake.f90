subroutine lakecells(ilake)

  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_lake,    only: lake, omlake
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use consts_coms, only: grav, p00i, rocp, eps_virt, cpi
  use mem_sfcnud,  only: sfcwat_nud, sfctemp_nud, fracliq_nud
  use oname_coms,  only: nl
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in) :: ilake

  ! Local variables

  integer :: iwsfc

  real :: canexneri, cantheta, canthetav
  real :: airthetav, wstars

  real, parameter :: onethird = 1./3.

  iwsfc = ilake + omlake

  ! Skip this cell if running in parallel and cell rank is not MYRANK

  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  ! Update LAKE fields

  if (nl%igw_spinup /= 1) then ! Standard canopy interaction

     ! Evaluate surface layer exchange coefficients vkmsfc and vkhsfc

     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))
     canexneri = 1. / sfcg%canexner(iwsfc)

     cantheta  = sfcg%cantemp(iwsfc) * canexneri
     canthetav = cantheta * (1.0 + eps_virt * sfcg%canrrv(iwsfc))

     wstars = (grav * sfcg%pblh(iwsfc) * max(sfcg%wthv(iwsfc),0.0) / airthetav) ** onethird

     call stars(sfcg%dzt_bot (iwsfc), &
                sfcg%rough   (iwsfc), &
                sfcg%vels    (iwsfc), &
                sfcg%rhos    (iwsfc), &
                wstars              , &
                airthetav           , &
                canthetav           , &
                sfcg%vkmsfc  (iwsfc), &
                sfcg%vkhsfc  (iwsfc), &
                sfcg%ustar   (iwsfc), &
                sfcg%ggaer   (iwsfc)  )

     call lakecell(ilake, iwsfc,             &
                   sfcg%rhos        (iwsfc), &
                   sfcg%ustar       (iwsfc), &
                   sfcg%vkhsfc      (iwsfc), &
                   sfcg%can_depth   (iwsfc), &
                   sfcg%topw        (iwsfc), &
                   sfcg%bathym      (iwsfc), &
                   sfcg%rshort      (iwsfc), &
                   sfcg%rlong       (iwsfc), &
                   sfcg%rlongup     (iwsfc), &
                   sfcg%albedo_beam (iwsfc), &
                   sfcg%pcpg        (iwsfc), &
                   sfcg%qpcpg       (iwsfc), &
                   sfcg%glatw       (iwsfc), &
                   sfcg%glonw       (iwsfc), &
                   sfcg%airtheta    (iwsfc), &
                   sfcg%airrrv      (iwsfc), &
                   sfcg%canexner    (iwsfc), &
                   sfcg%cantemp     (iwsfc), &
                   sfcg%canrrv      (iwsfc), &
                   lake%depth       (ilake), &
                   lake%lake_energy (ilake), &
                   sfcg%head1       (iwsfc), &
                   sfcg%runoff      (iwsfc), &
                   sfcg%sfluxt      (iwsfc), &
                   sfcg%sfluxr      (iwsfc), &
                   sfcg%rough       (iwsfc)  )

! Original calculation of wthv for sfluxt units [kg_dry K m^-2 s^-1]

!     sfcg%wthv(iwsfc) = ( sfcg%sfluxt(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
!          + sfcg%sfluxr(iwsfc) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

! New calculation of wthv for sfluxt units [W m^-2]

     sfcg%wthv(iwsfc) = ( sfcg%sfluxt(iwsfc) * cpi * canexneri * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
          + sfcg%sfluxr(iwsfc) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

  else ! "Fast canopy" nudging

     call lakecell_nud(iwsfc, ilake,            &
                       lake%depth      (ilake), &
                       lake%lake_energy(ilake), &
                       sfcg%topw       (iwsfc), &
                       sfcg%bathym     (iwsfc), &
                       sfcg%head1      (iwsfc), &
                       sfcg%runoff     (iwsfc), &
                       sfcwat_nud      (iwsfc), &
                       sfctemp_nud     (iwsfc), &
                       fracliq_nud     (iwsfc)  )

  endif

end subroutine lakecells

!===============================================================================

subroutine lakecell(ilake, iwsfc, rhos, ustar, vkhsfc, can_depth, topw, &
                    bathym, rshort, rlong, rlongup, albedo_beam, pcpg, qpcpg, &
                    glatw, glonw, airtheta, airrrv, canexner, cantemp, canrrv, depth, &
                    lake_energy, head1, runoff, sfluxt, sfluxr, rough)

  use lake_coms,   only: dt_lake
  use consts_coms, only: cp, grav, t00, cliq1000, alvl, r8
  use therm_lib,   only: rhovsl, qtk
  use sea_swm,     only: depthmin_flux
  use mem_sfcg,    only: sfcg
  use matrix,      only: matrix8_NxN, matrix8_3x3
  use leaf4_canopy,only: sing_print

  implicit none

  integer, intent(in)    :: ilake        ! current lake cell index
  integer, intent(in)    :: iwsfc        ! current sfc grid cell index
  real,    intent(in)    :: rhos         ! air density [kg/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: vkhsfc       ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: topw         ! topographic height of sfc W points [m]
  real,    intent(in)    :: bathym       ! bathymetric height of sfc W points [m]
  real,    intent(in)    :: rshort       ! downward can-top s/w flux [W/m^2]
  real,    intent(in)    :: rlong        ! downward can-top l/w flux [W/m^2]
  real,    intent(in)    :: rlongup      ! upward can-top l/w flux [W/m^2]
  real,    intent(in)    :: albedo_beam  ! surface s/w beam albedo [0-1]
  real,    intent(in)    :: pcpg         ! new pcp amount this timestep [kg/m^2]
  real,    intent(in)    :: qpcpg        ! new pcp energy this timestep [J/m^2]
  real,    intent(in)    :: glatw        ! Latitude of lake cell 'center' [deg]
  real,    intent(in)    :: glonw        ! Longitude of lake cell 'center' [deg]
  real,    intent(in)    :: airtheta     ! atm potential temp [K]
  real,    intent(in)    :: airrrv       ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: canexner     ! canopy Exner function []
  real,    intent(inout) :: cantemp      ! "canopy" air temp [K]
  real,    intent(inout) :: canrrv       ! "canopy" air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(inout) :: depth        ! lake mean depth [m]
  real,    intent(inout) :: lake_energy  ! lake energy lake energy [J/kg]
  real,    intent(inout) :: head1        ! lake water hydraulic head (rel to topo datum) [m]
  real,    intent(inout) :: runoff       ! new runoff mass this timestep [kg/m^2]
  real,    intent(out)   :: sfluxt       ! can_air to atm heat flux [kg_dry K m^-2 s^-1]
  real,    intent(out)   :: sfluxr       ! can_air to atm vapor flux [kg_vap m^-2 s^-1]
  real,    intent(out)   :: rough        ! lake cell roughess height [m]

  ! Local parameters

  real, parameter :: z0fac_water = .016 / grav      ! factor for Charnok roughness height
  real, parameter :: ozo = 1.59e-5                  ! base roughness height in HWRF
  real, parameter :: lake_runoff_time = 3. * 86400. ! time scale for lake runoff [s]
  real, parameter :: lake_head1_thresh = 1.0        ! head1 threshold for lake runoff [m]
  real, parameter :: fcn = 0.75                     ! Crank-Nicolson future time weight

  ! Local variables

  real :: rdi       ! canopy conductance [m/s]
  real :: hxfersc   ! heat xfer from lake surface to can_air this step [J/m^2]
  real :: wxfersc   ! vapor xfer from lake surface to can_air this step [kg_vap/m^2]
  real :: hxferca   ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca   ! vapor xfer from can_air to atm this step [kg/m^2]
  real :: sfc_rhovs ! sat vapor density at lake surface temp [kg_vap/m^3]
  real :: radsfc    ! Radiation absorbed by surface [J/m^2]
  real :: lake_epm2 ! lake energy expressed in units of [J/m^2]
  real :: can_rhov      ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair        ! Canopy air mass [kg/m^2]
  real :: hcapcan       ! Canopy air heat capacity [J/(m^2 K)]
  real :: canairi       ! Inverse of canair
  real :: hcapcani      ! Inverse of hcapcan

  real :: zn1, zn2, zw
  real :: laketemp, fracliq, dheight

  real(r8) :: a5, a6, a9, a10
  real(r8) :: h4, h7, h8
  real(r8) :: y2, y5, y9, y10

  real(r8) :: aa3(3,3), xx3(3), yy3(3)  ! 3x3 matrix equation terms
  real(r8) :: aa4(4,4), xx4(4), yy4(4)  ! 4x4 matrix equation terms

  logical :: sing

  ! Diagnose lake temperature and liquid fraction

  call qtk(lake_energy, laketemp, fracliq)

! Evaluate surface saturation vapor density and mixing ratio of sea surface

  sfc_rhovs    = rhovsl(laketemp - 273.15)

  ! rdi = ustar/5 is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  radsfc = dt_lake * (rshort * (1. - albedo_beam) + rlong - rlongup)

  ! Canopy air quantities

  can_rhov = canrrv * rhos
  canair   = rhos * can_depth
  canairi  = 1. / canair
  hcapcan  = cp * canair
  hcapcani = 1. / hcapcan

  ! Set up and solve a linear system of equations that use trapezoidal-implicit
  ! differencing.  The solution of the system consists of turbulent heat and
  ! water vapor fluxes between canopy air, the lake surface, and the free
  ! atmosphere, and the consequent changes to water and temperature of canopy
  ! air and the lake surface.

  ! It is assumed that the heat capacity of the (upper level of) lake water
  ! is sufficiently high that lake temperature change is negligible in the
  ! context of the matrix solver.

  a5  = dt_lake * rdi   ! sfc sublayer vap xfer coef
  a6  = cp * rhos * a5  ! sfc sublayer heat xfer coef
  a9  = dt_lake * vkhsfc / sfcg%dzt_bot(iwsfc)
  a10 = cp * a9

  h4 = fcn * rhos * canairi  ! = fcn / can_depth
  h7 = fcn * hcapcani
  h8 = fcn * canairi

  y2  = sfc_rhovs - can_rhov
  y5  = laketemp  - cantemp
  y9  = canrrv    - airrrv
  y10 = cantemp   - canexner * airtheta

  aa4(1,1) = 1._r8 + a5 * h4
  aa4(1,2) = 0._r8
  aa4(1,3) =       - a5 * h4
  aa4(1,4) = 0._r8
  yy4(1)   =         a5 * y2   ! WSC row

  aa4(2,1) = 0._r8
  aa4(2,2) = 1._r8 + a6 * h7
  aa4(2,3) = 0._r8
  aa4(2,4) =       - a6 * h7
  yy4(2)   =         a6 * y5   ! HSC row

  aa4(3,1) =       - a9 * h8
  aa4(3,2) = 0._r8
  aa4(3,3) = 1._r8 + a9 * h8
  aa4(3,4) = 0._r8
  yy4(3)   =         a9 * y9   ! WCA row

  aa4(4,1) = 0._r8
  aa4(4,2) =       - a10 * h7
  aa4(4,3) = 0._r8
  aa4(4,4) = 1._r8 + a10 * h7
  yy4(4)   =         a10 * y10 ! HCA row

  call matrix8_NxN(4,aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'lake1',4,aa4,yy4,glatw,glonw)

  if (depth > depthmin_flux .or. xx4(1) < 0.) then

     wxfersc = xx4(1)
     hxfersc = xx4(2)
     wxferca = xx4(3)
     hxferca = xx4(4)

  else

     ! If lake depth is below min depth threshold AND the lake continues to evaporate,
     ! solve a reduced equation set in which evaporation is not permitted.

     aa3(1,1) = 1._r8 + a6 * h7
     aa3(1,2) = 0._r8
     aa3(1,3) =       - a6 * h7
     yy3(1)   =         a6 * y5  ! HSC row

     aa3(2,1) = 0._r8
     aa3(2,2) = 1._r8 + a9 * h8
     aa3(2,3) = 0._r8
     yy3(2)   =         a9 * y9   ! WCA row

     aa3(3,1) =       - a10 * h7
     aa3(3,2) = 0._r8
     aa3(3,3) = 1._r8 + a10 * h7
     yy3(3)   =         a10 * y10 ! HCA row

     call matrix8_3x3(aa3,yy3,xx3,sing); if (sing) call sing_print(iwsfc,'lake2',3,aa3,yy3,glatw,glonw)

     wxfersc = 0.
     hxfersc = xx3(1)
     wxferca = xx3(2)
     hxferca = xx3(3)

  endif

  cantemp = cantemp + (hxfersc - hxferca) * hcapcani
  canrrv  = canrrv  + (wxfersc - wxferca) * canairi

  sfluxt = hxferca / dt_lake
  sfluxr = wxferca / dt_lake

  head1 = depth + bathym - topw

  if (head1 > lake_head1_thresh) then
     runoff = (head1 - lake_head1_thresh) * dt_lake / lake_runoff_time
  else
     runoff = 0.
  endif

  lake_epm2 = lake_energy * depth * 1000.

  ! Update lake water depth, energy, and head1

  dheight = 0.001 * (pcpg - wxfersc) - runoff

  depth = depth + dheight

  lake_epm2 = lake_epm2 + radsfc - hxfersc - wxfersc * alvl + qpcpg

  lake_energy = lake_epm2 / (depth * 1000.) ! water density = 1000 kg/m^3

  head1 = head1 + dheight

  ! Evaluate lake roughness height

  ! Charnok (1955):
  ! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)

  ! Davis et al. (2008) originally used in HWRF
  ! rough = 10. * exp(-10. / ustar ** .333333)
  ! rough = max(.125e-6, min(2.85e-3,rough))

  ! 2012 HWRF scheme; interpolates between the Charnok scheme at low wind
  ! and the Davis et al. curve fit at high wind speeds

  zw    = min( (ustar/1.06)**0.3, 1.0)
  zn1   = 0.011 * ustar * ustar / grav + ozo
  zn2   = 10. * exp(-9.5 * ustar**(-.3333333)) + 1.65e-6 / ustar
  rough = (1.0-zw) * zn1 + zw * zn2
  rough = min( rough, 2.85e-3)

end subroutine lakecell

!===============================================================================

subroutine lakecell_nud(iwsfc, ilake, depth, lake_energy, topw, bathym, head1, &
                    runoff, sfcwat_nud, sfctemp_nud, fracliq_nud)

  use lake_coms,   only: dt_lake
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in)    :: iwsfc       ! current sfc grid cell index
  integer, intent(in)    :: ilake       ! current lake cell index
  real,    intent(inout) :: depth       ! lake mean depth [m]
  real,    intent(inout) :: lake_energy ! lake energy lake energy [J/kg]
  real,    intent(in)    :: topw        ! topographic height of sfc W points [m]
  real,    intent(in)    :: bathym      ! bathymetric height of sfc W points [m]
  real,    intent(inout) :: head1       ! lake water hydraulic head (rel to topo datum) [m]
  real,    intent(inout) :: runoff      ! new runoff mass this timestep [kg/m^2]
  real,    intent(in)    :: sfcwat_nud  !
  real,    intent(in)    :: sfctemp_nud !
  real,    intent(in)    :: fracliq_nud !

  ! Local parameters

  real, parameter :: lake_runoff_time = 3. * 86400. ! time scale for lake runoff [s]
  real, parameter :: lake_head1_thresh = 1.0        ! head1 threshold for lake runoff [m]

  ! Local variables

  real :: lake_epm2 ! lake energy expressed in units of [J/m^2]
  real :: laketemp, fracliq, ediff, eflux

  ! Diagnose lake temperature and liquid fraction

  call qtk(lake_energy, laketemp, fracliq)

  ! Compute surface energy flux based on difference between model and nudging
  ! values of sfcwater_tempk and sfcwater_fracliq

  ediff = sfctemp_nud - laketemp + 80. * (fracliq_nud - fracliq)

  eflux = max(-500.,min(500.,20. * ediff))

  lake_epm2 = lake_energy * depth * 1000.

  ! Update lake water depth, energy, and head1

  depth = max(0., depth + sfcwat_nud * 0.001)

  lake_epm2 = lake_epm2 + eflux * dt_lake

  lake_energy = lake_epm2 / (depth * 1000.) ! water density = 1000 kg/m^3

  head1 = depth + bathym - topw

  ! Compute lake runoff; no change to lake_energy needed since [J/kg] remains the same

  if (head1 > lake_head1_thresh) then
     runoff = (head1 - lake_head1_thresh) * dt_lake / lake_runoff_time
     depth = depth - runoff
     head1 = head1 - runoff
  endif

end subroutine lakecell_nud
