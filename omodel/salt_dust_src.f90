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
subroutine sea_spray(mrl)

  ! Compute sea spray number concentration fluxes, applied as tendencies in
  ! the grid level adjacent to the sea surface.  Based on O’Dowd, C. D.,
  ! Smith, M. H., and Jennings, S. G., (1993): Submicron aerosol, radon and
  ! soot carbon characteristics over the North East Atlantic. J. Geophys.
  ! Res., 98, 1132-1136.  Updated for O'Dowd (1997, 1999).

  use mem_grid,    only: dzt_bot, volti
  use mem_basic,   only: vxe, vye, vze, rho
  use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
  use micro_coms,  only: iccn, igccn
  use ccnbin_coms, only: isalt
  use mem_micro,   only: ccntyp, con_gccn
  use mem_tend,    only: con_gccnt
  use mem_sea,     only: sea, itab_ws

  implicit none

  integer, intent(in) :: mrl

  real, parameter :: timescale_salt = 3600.0 ! relaxation time [s]
  real, parameter :: depthscale_salt = 100.  ! assumed depth scale filled by fluxes [m]
  real, parameter :: dssotss = timescale_salt / depthscale_salt
  integer :: j, iw, nsea, jws, iws, kw
  real :: vels, vels10, film_diag, jet_diag, film_flux, jet_flux

  ! Return if neither ccn nor gccn are prognosed

  if (isalt < 1 .and. igccn < 2) return

  ! Loop over prognosed atmospheric grid columns

  !-----------------------------------------------------------------------------
  !$omp parallel do private(iw,nsea,jws,iws,kw,vels,vels10,film_diag,jet_diag, &
  !$omp                      film_flux,jet_flux)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
  !-----------------------------------------------------------------------------

     ! Check for sea area beneath this atmospheric grid column; cycle if none

     nsea = itab_w(iw)%nsea
     if (nsea == 0) cycle

     ! Loop over sea cells beneath this atmospheric grid column

     do jws = 1,nsea
        iws = itab_w(iw)%isea(jws)
        kw = itab_ws(iws)%kw

        ! Skip if not seawater
        if (sea%leaf_class(iws) /= 0) cycle

        ! Diagnose wind speed in grid cell and at 10 m height, and limit the
        ! latter to a maximum of 20 m/s.

        vels = sqrt(vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2)

        vels10 = min(20., vels * log(10.         / sea%sea_rough(iws)) &
                               / log(dzt_bot(kw) / sea%sea_rough(iws)))

        ! Diagnose equilibrium sea salt concentrations [#/m^3] for 10 m wind speed

        film_diag = 10.**(0.0950 * vels10 + 6.2830) ! Film mode diameter 0.1 micron
        jet_diag  = 10.**(0.0422 * vels10 + 5.7122) ! Jet mode diameter 1.0 micron
        !spm_diag = 10.**(0.0690 * vels10 + 0.1900)

        ! If model surface concentrations of sea salt are less than diagnosed
        ! values, compute surface fluxes [#/(m^2 s)] for film and jet modes
        ! that would bring model concentrations up to diagnosed concentrations
        ! in a time scale of timescale_salt, assuming that the model deficit
        ! applies and must be filled over a depth scale of depthscale_salt.
        ! Fixing the depth scale in this way makes the fluxes independent of
        ! the chosen vertical grid spacing.  

        ! Convert surface fluxes to concentration tendency [#/(m^3 s)] in
        ! atmospheric cells that are adjacent to the sea surface.  Tendency is
        ! obtained by multiplying fluxes by surface area and dividing by grid
        ! cell volume.

        if (isalt > 0) then
           film_flux = dssotss * (film_diag - ccntyp(isalt)%con_ccn(kw,iw) * rho(kw,iw))

           if (film_flux > 0.) then
                ccntyp(isalt)%con_ccnt(kw,iw) &
              = ccntyp(isalt)%con_ccnt(kw,iw) &
              + film_flux * sea%area(iws) * volti(kw,iw) * (1.0 - sea%seaicec(iws))
           endif
        endif

        if (igccn == 2) then
           jet_flux = dssotss * ( jet_diag - con_gccn(kw,iw) * rho(kw,iw))

           if (jet_flux > 0.) then
                con_gccnt(kw,iw) &
              = con_gccnt(kw,iw) &
              + jet_flux * sea%area(iws) * volti(kw,iw) * (1.0 - sea%seaicec(iws))
           endif
        endif

     enddo  ! jws/iws

  enddo  ! iw
  !$omp end parallel do

end subroutine sea_spray

!===============================================================================

subroutine dust_src(mrl)

  use mem_grid,    only: dzt_bot, volti
  use mem_basic,   only: vxe, vye, vze, rho
  use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
  use micro_coms,  only: iccn, igccn
  use ccnbin_coms, only: idust1, idust2
  use mem_micro,   only: ccntyp, con_gccn
  use mem_leaf,    only: land, itab_wl
  use leaf_coms,   only: nzg, slmsts_ch, slmsts_vg
  use mem_tend,    only: con_gccnt

  use nuclei_coms, only: nbin, fact, sp, nodust, amassi, uth, wprime

  implicit none

  integer, intent(in) :: mrl

  integer :: j, iw, nland, jwl, iwl, kw
  integer :: leaf_class, nts, ibin
  real :: vels, vels10, vels10_cm
  real :: gwet, wetfac, flux_m, flux_n1, flux_n2

  ! Return if neither ccn nor gccn are prognosed

  if (idust1 < 1 .and. idust2 < 2) return

  ! Loop over prognosed atmospheric grid columns

  !-----------------------------------------------------------------------------
  !$omp parallel do private(iw,nland,jwl,iwl,leaf_class,nts,gwet,kw,vels, &
  !$omp                     vels10,wetfac,flux_n1,flux_n2,vels10_cm,flux_m)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
  !-----------------------------------------------------------------------------

     ! Check for land area beneath this atmospheric grid column; cycle if none

     nland = itab_w(iw)%nland
     if (nland == 0) cycle

     ! Loop over land cells beneath this atmospheric grid column

     do jwl = 1,nland
        iwl = itab_w(iw)%iland(jwl)
        leaf_class = land%leaf_class(iwl)

        ! If leaf_class is not type that can generate dust, skip this land cell

        if (nodust(leaf_class)) cycle

        ! Diagnose fractional soil wetness, and cycle if >= 0.5

        nts = land%ntext_soil(nzg,iwl)
        if (land%flag_vg(iwl)) then
           gwet = land%soil_water(nzg,iwl) / slmsts_vg(nts)
        else
           gwet = land%soil_water(nzg,iwl) / slmsts_ch(nts)
        endif

        if (gwet >= 0.5) cycle

         kw = itab_wl(iwl)%kw

        ! Diagnose wind speed in grid cell and at 10 m height, and limit the
        ! latter to a maximum of 20 m/s.

        vels = sqrt(vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2)

        vels10 = min(20., vels * log(10.         / land%rough(iwl)) &
                               / log(dzt_bot(kw) / land%rough(iwl)))

        if (gwet * 100. > wprime(nts)) then
           wetfac = sqrt(1. + 1.21 * (gwet * 100. - wprime(nts))**0.68)
        else
           wetfac = 1.
        endif

        flux_n1 = 0.
        flux_n2 = 0.

        do ibin = 1,nbin

           !> Mass flux [g cm^-2 s^-1] from Ginoux et al. 2001 (Eq. 2) 

           vels10_cm = 100. * vels10
           flux_m = fact * sp(ibin) * (vels10_cm**2.) * (vels10_cm - uth(ibin) * wetfac)

           !> Convert mass flux to [#/(m^2 s)]

           if (flux_m > 0.) then
              if (ibin < 5) then
                 flux_n1 = flux_n1 + flux_m * 1.e4 * amassi(ibin)
              else
                 flux_n2 = flux_n2 + flux_m * 1.e4 * amassi(ibin)
              endif
           endif

        enddo

        ! Convert surface fluxes to concentration tendency [#/(m^3 s)] in
        ! atmospheric cells that are adjacent to the land surface.  Tendency is
        ! obtained by multiplying fluxes by surface area and dividing by grid
        ! cell volume.

        if (idust1 > 0 .and. flux_n1 > 0.) then
             ccntyp(idust1)%con_ccnt(kw,iw) &
           = ccntyp(idust1)%con_ccnt(kw,iw) &
           + flux_n1 * land%area(iwl) * volti(kw,iw)
        endif

        if (idust2 == 2 .and. flux_n2 > 0.) then
             ccntyp(idust2)%con_ccnt(kw,iw) &
           = ccntyp(idust2)%con_ccnt(kw,iw) &
           + flux_n2 * land%area(iwl) * volti(kw,iw)
        endif

     enddo  ! jws/iws

  enddo  ! iw
  !$omp end parallel do

end subroutine dust_src

