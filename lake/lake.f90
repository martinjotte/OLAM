!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================
subroutine lakecells()

  use mem_ijtabs,  only: itabg_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_lake,    only: lake, mlake, omlake
  use misc_coms,   only: io6, time8, isubdomain
  use mem_para,    only: myrank
  use mem_basic,   only: rho, press, theta, tair, vxe, vye, vze, rr_v
  use mem_micro,   only: rr_c
  use consts_coms, only: grav
  use mem_sfcnud,  only: sfcwat_nud, sfctemp_nud, fracliq_nud
  use oname_coms,  only: nl

  implicit none

! Local variables

  integer :: ilake             ! lake cell loop counter
  integer :: iwsfc

! Loop over ALL LAKE CELLS

  !$omp parallel do private (iwsfc)
  do ilake = 2, mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

! Update LAKE fields

     if (nl%igw_spinup /= 1) then ! Standard canopy interaction

        call lakecell(iwsfc, ilake,             &
                      lake%depth       (ilake), &
                      lake%lake_energy (ilake), &
                      lake%surface_srrv(ilake), &
                      sfcg%rhos        (iwsfc), &
                      sfcg%ustar       (iwsfc), &
                      sfcg%sxfer_t     (iwsfc), &
                      sfcg%sxfer_r     (iwsfc), &
                      sfcg%can_depth   (iwsfc), &
                      sfcg%cantemp     (iwsfc), &
                      sfcg%canrrv      (iwsfc), &
                      sfcg%rough       (iwsfc), &
                      sfcg%head1       (iwsfc), &
                      sfcg%rshort      (iwsfc), &
                      sfcg%rlong       (iwsfc), &
                      sfcg%rlongup     (iwsfc), &
                      sfcg%albedo_beam (iwsfc), &
                      sfcg%pcpg        (iwsfc), &
                      sfcg%qpcpg       (iwsfc), &
                      sfcg%runoff      (iwsfc), &
                      itab_wsfc(iwsfc)%ivoronoi)

     else ! "Fast canopy" nudging

        call lakecell_nud(iwsfc, ilake,            &
                          lake%depth      (ilake), &
                          lake%lake_energy(ilake), &
                          sfcg%head1      (iwsfc), &
                          sfcg%runoff     (iwsfc), &
                          itab_wsfc(iwsfc)%ivoronoi, &
                          sfcwat_nud      (iwsfc), &
                          sfctemp_nud     (iwsfc), &
                          fracliq_nud     (iwsfc)  )

     endif

! Zero out sfcg%SXFER_T(iwsfc) and sfcg%SXFER_R(iwsfc) now that they have 
! been applied to the canopy

     sfcg%sxfer_t(iwsfc) = 0.
     sfcg%sxfer_r(iwsfc) = 0.
   ! sfcg%sxfer_c(iwsfc) = 0.

  enddo
  !$omp end parallel do

end subroutine lakecells

!===============================================================================

subroutine lakecell(iwsfc, ilake, depth, lake_energy, surface_srrv, rhos, ustar, &
                    sxfer_t, sxfer_r, can_depth, cantemp, canrrv, rough, head1, &
                    rshort, rlong, rlongup, albedo_beam, pcpg, qpcpg, runoff, ivoronoi)

  use lake_coms,   only: dt_lake
  use consts_coms, only: cp, grav, t00, cliq1000, alvl
  use misc_coms,   only: io6
  use therm_lib,   only: rhovsil, qtk

  implicit none

  integer, intent(in)    :: iwsfc       ! current sfc grid cell index
  integer, intent(in)    :: ilake       ! current lake cell index
  real,    intent(in)    :: depth       ! lake mean depth [m]
  real,    intent(inout) :: lake_energy ! lake energy lake energy [J/kg]
  real,    intent(out)   :: surface_srrv! lake surface sat mixing ratio [kg_vap/kg_dryair] 
  real,    intent(in)    :: rhos        ! air density [kg/m^3]
  real,    intent(in)    :: ustar       ! friction velocity [m/s]
  real,    intent(in)    :: sxfer_t     ! can_air to atm heat xfer this step [kg_air K/m^2]
  real,    intent(in)    :: sxfer_r     ! can_air to atm vapor xfer this step (kg_vap/m^2]
  real,    intent(in)    :: can_depth   ! "canopy" depth for heat and vap capacity [m]
  real,    intent(inout) :: cantemp     ! "canopy" air temp [K]
  real,    intent(inout) :: canrrv      ! "canopy" air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(out)   :: rough       ! lake cell roughess height [m] 
  real,    intent(inout) :: head1       ! lake water hydraulic head (rel to topo datum) [m]
  real,    intent(in)    :: rshort      ! downward can-top s/w flux [W/m^2]
  real,    intent(in)    :: rlong       ! downward can-top l/w flux [W/m^2]
  real,    intent(in)    :: rlongup     ! upward can-top l/w flux [W/m^2]
  real,    intent(in)    :: albedo_beam ! surface s/w beam albedo [0-1]
  real,    intent(in)    :: pcpg        ! new pcp amount this timestep [kg/m^2]
  real,    intent(in)    :: qpcpg       ! new pcp energy this timestep [J/m^2]
  real,    intent(inout) :: runoff      ! new runoff mass this timestep [kg/m^2]
  integer, intent(in)    :: ivoronoi    ! Voronoi flag (to control water mass eqn.)

  ! Local parameters

  real, parameter :: z0fac_water = .016 / grav      ! factor for Charnok roughness height
  real, parameter :: ozo = 1.59e-5                  ! base roughness height in HWRF
  real, parameter :: lake_runoff_time = 3. * 86400. ! time scale for lake runoff [s]
  real, parameter :: lake_head1_thresh = 1.0        ! head1 threshold for lake runoff [m]

  ! Local variables

  real :: rdi     ! canopy conductance [m/s]
  real :: hxfersc ! heat xfer from lake surface to can_air this step [J/m^2]
  real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxfersc ! vapor xfer from lake surface to can_air this step [kg_vap/m^2]
  real :: radsfc  ! Radiation absorbed by surface [J/m^2]
  real :: del_energy_per_m2 ! change in lake energy this timestep [J/m^2]

  real :: zn1, zn2, zw
  real :: laketemp, fracliq

  ! Diagnose lake temperature and liquid fraction

  call qtk(lake_energy, laketemp, fracliq)

  ! Evaluate surface saturation specific humidity

  surface_srrv = rhovsil(laketemp-t00) / rhos

  ! Update temperature and vapor specific humidity of "canopy" from
  ! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
  ! is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  radsfc = dt_lake * (rshort * (1. - albedo_beam) + rlong - rlongup)

  hxfersc = dt_lake * cp * rhos * rdi * (laketemp - cantemp)
  wxfersc = dt_lake *      rhos * rdi * (surface_srrv - canrrv)

  hxferca = cp * sxfer_t  ! sxfer_t and sxfer_r already incorporate dt_lake

  cantemp = cantemp + (hxfersc - hxferca) / (can_depth * rhos * cp)
  canrrv  = canrrv  + (wxfersc - sxfer_r) / (can_depth * rhos)

  del_energy_per_m2 = radsfc - hxfersc - wxfersc * alvl + qpcpg

  lake_energy = lake_energy + del_energy_per_m2 / (depth * 1000.) ! water density = 1000 kg/m^3

  ! Evaluate lake roughness height

  ! Charnok (1955):
  ! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)

  ! Davis et al. (2008) originally used in HWRF
  ! rough = 10. * exp(-10. / ustar ** .333333)
  ! rough = max(.125e-6, min(2.85e-3,rough))

  ! 2012 HWRF scheme; interpolates between the Charnok scheme at low wind
  ! and the Davis et al. curve fit at high wind speeds

  zw    = min( (ustar/1.06)**0.3, 1.0)
  zn1   = 0.011 * ustar * ustar /grav + ozo
  zn2   = 10. * exp(-9.5 * ustar**(-.3333333)) + 1.65e-6 / ustar
  rough = (1.0-zw) * zn1 + zw * zn2
  rough = min( rough, 2.85e-3)

  ! If ivoronoi = 3, update lake water level (mass)

  if (ivoronoi == 3) then
     if (head1 > lake_head1_thresh) then
        runoff = (head1 - lake_head1_thresh) * dt_lake / lake_runoff_time
     else
        runoff = 0.
     endif

     head1 = head1 + 0.001 * (pcpg - wxfersc) - runoff
  endif

end subroutine lakecell

!===============================================================================

subroutine lakecell_nud(iwsfc, ilake, depth, lake_energy, head1, &
                    runoff, ivoronoi, sfcwat_nud, sfctemp_nud, fracliq_nud)

  use lake_coms,   only: dt_lake
  use misc_coms,   only: io6
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in)    :: iwsfc       ! current sfc grid cell index
  integer, intent(in)    :: ilake       ! current lake cell index
  real,    intent(in)    :: depth       ! lake mean depth [m]
  real,    intent(inout) :: lake_energy ! lake energy lake energy [J/kg]
  real,    intent(inout) :: head1       ! lake water hydraulic head (rel to topo datum) [m]
  real,    intent(inout) :: runoff      ! new runoff mass this timestep [kg/m^2]
  integer, intent(in)    :: ivoronoi    ! Voronoi flag (to control water mass eqn.)
  real,    intent(in)    :: sfcwat_nud  !
  real,    intent(in)    :: sfctemp_nud !
  real,    intent(in)    :: fracliq_nud !

  ! Local parameters

  real, parameter :: lake_runoff_time = 3. * 86400. ! time scale for lake runoff [s]
  real, parameter :: lake_head1_thresh = 1.0        ! head1 threshold for lake runoff [m]

  ! Local variables

  real :: del_energy_per_m2 ! change in lake energy this timestep [J/m^2]
  real :: laketemp, fracliq, ediff, flux

  ! Diagnose lake temperature and liquid fraction

  call qtk(lake_energy, laketemp, fracliq)

  ! Compute surface energy flux based on difference between model and nudging
  ! values of sfcwater_tempk and sfcwater_fracliq

  ediff = sfctemp_nud - laketemp + 80. * (fracliq_nud - fracliq)

  flux = max(-500.,min(500.,20. * ediff))

  del_energy_per_m2 = flux * dt_lake

  lake_energy = lake_energy + del_energy_per_m2 / (depth * 1000.) ! water density = 1000 kg/m^3

  ! If ivoronoi = 3, update lake water level (mass)

  if (ivoronoi == 3) then
     head1 = max(-2.0, head1 + sfcwat_nud * 0.001)

     if (head1 > lake_head1_thresh) then
        runoff = (head1 - lake_head1_thresh) * dt_lake / lake_runoff_time
        head1 = head1 - runoff
     endif
  endif

end subroutine lakecell_nud
