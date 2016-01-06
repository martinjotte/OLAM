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
Module mem_addgrid

  use consts_coms, only: r8
  use max_dims,    only: maxgrds, maxngrdll

  implicit none

  integer :: nzp_og
  integer :: nxp_og
  integer :: mdomain_og
  integer :: ngrids_og
  integer :: isfcl_og
  integer :: itopoflg_og
  integer :: ndz_og
  integer :: nzg_og
  integer :: nve2_max_og

  real :: deltax_og

  real, allocatable :: hdz_og(:)
  real, allocatable :: dz_og(:)

  integer :: nza_og
  integer :: nma_og
  integer :: nua_og
  integer :: nva_og
  integer :: nwa_og
  integer :: mrls_og
  integer :: nwl_og
  integer :: nws_og

  integer, allocatable :: lpw_og(:)  ! (nwa_og)
  integer, allocatable :: lve2_og(:) ! (nwa_og)

  real, allocatable :: wnx_og(:) ! (nwa_og)
  real, allocatable :: wny_og(:) ! (nwa_og)
  real, allocatable :: wnz_og(:) ! (nwa_og)

  real, allocatable :: xew_og(:) ! (nwa_og)
  real, allocatable :: yew_og(:) ! (nwa_og)
  real, allocatable :: zew_og(:) ! (nwa_og)

  real, allocatable :: xewl_og(:) ! (nwl_og)
  real, allocatable :: yewl_og(:) ! (nwl_og)
  real, allocatable :: zewl_og(:) ! (nwl_og)

  real, allocatable :: xews_og(:) ! (nws_og)
  real, allocatable :: yews_og(:) ! (nws_og)
  real, allocatable :: zews_og(:) ! (nws_og)

  real, allocatable :: vc_og(:,:) ! (nza_og,nwa_og)
  real, allocatable :: wc_og(:,:) ! (nza_og,nwa_og)

  real, allocatable :: vc_accum_og(:,:) ! (nza_og,nwa_og)

  real, allocatable :: vxe_og(:,:) ! (nza_og,nwa_og)
  real, allocatable :: vye_og(:,:) ! (nza_og,nwa_og)
  real, allocatable :: vze_og(:,:) ! (nza_og,nwa_og)

  real, allocatable :: vxe2_og(:,:) ! (nza_og,nwa_og)
  real, allocatable :: vye2_og(:,:) ! (nza_og,nwa_og)
  real, allocatable :: vze2_og(:,:) ! (nza_og,nwa_og)

  real, allocatable :: soil_water_og (:,:) ! (nzg_og,nwl_og)
  real, allocatable :: soil_energy_og(:,:) ! (nzg_og,nwl_og)

  integer, allocatable :: ntext_soil_og(:,:) ! (nzg_og,nwl_og)

  integer, allocatable :: leaf_class_og(:) ! (nwl_og)

  Type itab_wog_vars          ! data structure for OLD grid IW pts (global domain)
     integer :: npoly = 0     ! number of W neighbors of this W pt
     integer :: iv(7)         ! neighbor V pts
     integer :: iw(7)         ! neighbor W pts

     real :: ecvec_vx(7) = 0. ! factors converting V to earth cart. velocity
     real :: ecvec_vy(7) = 0. ! factors converting V to earth cart. velocity
     real :: ecvec_vz(7) = 0. ! factors converting V to earth cart. velocity
  End Type itab_wog_vars

  Type itab_wadd_vars         ! data structure for WADD pts (individual rank)
     integer :: iw_og(3) = 1  ! intrp pts
     real    :: f_og(3) = 0.  ! intrp coeffs
  End Type itab_wadd_vars

  Type itab_wladd_vars        ! data structure for WLADD pts (individual rank)
     integer :: iwl_og = 1  ! indes of nearest land point on OLD grid
  End Type itab_wladd_vars

  Type itab_wsadd_vars      ! data structure for WSADD pts (individual rank)
     integer :: iws_og = 1  ! index of nearest sea point on OLD grid
  End Type itab_wsadd_vars

  type (itab_wog_vars),    allocatable, target :: itab_wog(:)   ! (dimension nwa_og)
  type (itab_wadd_vars),   allocatable, target :: itab_wadd(:)  ! (dimension mwa)
  type (itab_wladd_vars),  allocatable, target :: itab_wladd(:) ! (dimension mwl)
  type (itab_wsadd_vars),  allocatable, target :: itab_wsadd(:) ! (dimension mws)

Contains

!==============================================================================

  subroutine init_addgrid()

  use misc_coms,  only: io6, mdomain, runtype, nxp
  use mem_ijtabs, only: itab_w
  use mem_grid,   only: mwa, xew, yew, zew, glatw, glonw
  use leaf_coms,  only: mwl, isfcl
  use mem_leaf,   only: land, itab_wl
  use sea_coms,   only: mws
  use mem_sea,    only: sea, itab_ws

  implicit none

  real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
  real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
  real :: dist, dist_min, dist_max, xi, yi
  real :: xin(7),yin(7)

  integer :: iw, iw_og, iwn, iwn_og, npoly, j, jmaxneg, jminpos, ngrw
  integer :: iwl, iws, iwl_og, iws_og

  allocate(itab_wadd(mwa))

! Compute a maximum possible distance between adjacent OLD grid IW points.
! It is safe to assume that any NEW grid IW point will be closer than this
! distance to at least one OLD grid IW point.

  dist_max = 7200.e3 / real(nxp)

! Loop over all local-domain IW points in NEW grid

  do iw = 2,mwa

     if (mod(iw,100000) == 0) write(io6,*) 'init_addgrid iw, mwa ',iw,mwa

! Skip points that do not need to be interpolated

     ngrw = 0
     do j = 1,itab_w(iw)%npoly
        iwn = itab_w(iw)%iw(j)

        if (ngrw < itab_w(iwn)%ngr) ngrw = itab_w(iwn)%ngr
     enddo

     if (itab_w(iw)%ngr <= ngrids_og .and. ngrw <= ngrids_og) cycle

! Initialize minimum distance variable

     dist_min = dist_max

! Loop over all global-domain IW points in OLD grid

     do iw_og = 2,nwa_og

! Skip interaction if NEW and OLD grid points are more than dist_min apart
! in any earth-frame direction

        if (abs(xew(iw) - xew_og(iw_og)) > dist_min) cycle
        if (abs(yew(iw) - yew_og(iw_og)) > dist_min) cycle
        if (abs(zew(iw) - zew_og(iw_og)) > dist_min) cycle

! Compute distance between NEW and OLD grid IW points, and reset minimum
! distance and OLD grid point index if new minimum is found

        dist = sqrt((xew(iw) - xew_og(iw_og))**2 &
                  + (yew(iw) - yew_og(iw_og))**2 &
                  + (zew(iw) - zew_og(iw_og))**2)

        if (dist_min > dist) then
           dist_min = dist
           itab_wadd(iw)%iw_og(1) = iw_og
        endif

     enddo

! Now that primary OLD grid point has been found for NEW grid IW point,
! compute its x,y components on a polar stereographic plane tangent
! at NEW grid point (NEW grid W point is at 0,0)

     iw_og = itab_wadd(iw)%iw_og(1)

     call e_ps(xew_og(iw_og),yew_og(iw_og),zew_og(iw_og), &
               glatw(iw),glonw(iw),xi,yi)

! Initialize vecprodz_minpos and vecprodz_maxneg

     vecprodz_maxneg = -1.e15
     vecprodz_minpos =  1.e15

! Loop through nearest polygon neighbors (j, iwn_og) of OLD grid point IW_og

     npoly = itab_wog(iw_og)%npoly

     do j = 1,npoly

! Get nudging point index (iwnudn) for current polygon neighbor of iwnud.

        iwn_og = itab_wog(iw_og)%iw(j)

! Compute x,y components of iwn_og polygon center on a polar stereographic 
! plane tangent at IW point

        call e_ps(xew_og(iwn_og),yew_og(iwn_og),zew_og(iwn_og), &
                  glatw(iw),glonw(iw),xin(j),yin(j))

! Compute z component (in polar stereographic space) of vector product of
! the vector from iwnud to iw (at 0,0) and the vector from iwnud to iwnudn.

        vecprodz = -xi * (yin(j) - yi) + yi * (xin(j) - xi)

! Compute scalar product of the vector from iwnud to iw (at 0,0) and the
! vector from iwnud to iwnudn in polar stereographic space.

        scalprod = -xi * (xin(j) - xi) - yi * (yin(j) - yi)

! Check whether scalar product is positive for current J point.  If so,
! J point is a candidate for the nudging triad for IW.

        if (scalprod > 0.) then

! Identify maximum negative vecprodz among all iwn_og polygon neighbors of
! iw_og.  This iwn_og will be second point of HISTORY grid IW_og triad for IW

           if (vecprodz < 0. .and. vecprodz > vecprodz_maxneg) then
              vecprodz_maxneg = vecprodz
              jmaxneg = j
              itab_wadd(iw)%iw_og(2) = iwn_og
           endif

! Identify minimum positive vecprodz among all iwn_og polygon neighbors of
! iw_og.  This iwn_og will be third point of HISTORY grid IW_og triad for IW

           if (vecprodz >= 0. .and. vecprodz < vecprodz_minpos) then
              vecprodz_minpos = vecprodz
              jminpos = j
              itab_wadd(iw)%iw_og(3) = iwn_og
           endif
        endif

     enddo

! Lastly, fill 3 HISTORY grid weight coefficients for this W point.
! Weights are computed in 2_d polar stereographic space.

! Invert matrix of coordinates

     call matinv3x3(1.,xi,yi, &
                    1.,xin(jmaxneg),yin(jmaxneg), &
                    1.,xin(jminpos),yin(jminpos), &
                    b11,b21,b31,b12,b22,b32,b13,b23,b33)

! Assign coefficients

     itab_wadd(iw)%f_og(1) = b11
     itab_wadd(iw)%f_og(2) = b21
     itab_wadd(iw)%f_og(3) = b31

  enddo

! Check whether LAND/SEA models are used

  if (isfcl == 1) then

! Land and sea cells are active.

  allocate(itab_wladd(mwl))
  allocate(itab_wsadd(mws))

! Loop through all NEW grid land cells, find the nearest land cell on the
! OLD grid, and attach its index to the NEW grid land cell.

     !$omp parallel do private(iw,iwl_og,dist_min,dist)
     do iwl = 2,mwl
        iw = itab_wl(iwl)%iw

        if (mod(iwl,100000) == 0) write(io6,*) 'init_addgrid iwl,mwl ',iwl,mwl

! Initialize minimum distance variable

        dist_min = 1.e12

! Loop over all local-domain IWL points in OLD grid

        do iwl_og = 2,nwl_og

! Skip interaction if NEW and OLD grid points are more than dist_min apart
! in any earth-frame direction

           if (abs(land%xew(iwl) - xewl_og(iwl_og)) > dist_min) cycle
           if (abs(land%yew(iwl) - yewl_og(iwl_og)) > dist_min) cycle
           if (abs(land%zew(iwl) - zewl_og(iwl_og)) > dist_min) cycle

! Compute distance between NEW and OLD grid IWL points, and reset minimum
! distance and HISTORY grid point index if new minimum is found

           dist = sqrt((land%xew(iwl) - xewl_og(iwl_og))**2 &
                     + (land%yew(iwl) - yewl_og(iwl_og))**2 &
                     + (land%zew(iwl) - zewl_og(iwl_og))**2)

           if (dist_min > dist) then
              dist_min = dist
              itab_wladd(iwl)%iwl_og = iwl_og
           endif

        enddo  ! iwl_og

     enddo  ! iwl
     !$omp end parallel do

! Loop through all NEW grid sea cells, find the nearest sea cell on the
! OLD grid, and attach its index to the NEW grid sea cell.

     do iws = 2,mws
        iw = itab_ws(iws)%iw

        if (mod(iws,100000) == 0) write(io6,*) 'init_addgrid iws,mws ',iws,mws

! Initialize minimum distance variable

        dist_min = 1.e12

! Loop over all local-domain IWS points in OLD grid

        do iws_og = 2,nws_og

! Skip interaction if NEW and OLD grid points are more than dist_min apart
! in any earth-frame direction

          if (abs(sea%xew(iws) - xews_og(iws_og)) > dist_min) cycle
          if (abs(sea%yew(iws) - yews_og(iws_og)) > dist_min) cycle
          if (abs(sea%zew(iws) - zews_og(iws_og)) > dist_min) cycle

! Compute distance between NEW and OLD grid IWS points, and reset minimum
! distance and HISTORY grid point index if new minimum is found

          dist = sqrt((sea%xew(iws) - xews_og(iws_og))**2 &
                    + (sea%yew(iws) - yews_og(iws_og))**2 &
                    + (sea%zew(iws) - zews_og(iws_og))**2)

          if (dist_min > dist) then
             dist_min = dist
             itab_wsadd(iws)%iws_og = iws_og
          endif

       enddo  ! iws_og

     enddo  ! iws

  endif  ! isfcl = 1

  end subroutine init_addgrid

!=========================================================================

  subroutine interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            ivara1,ivara2,rvara1,rvara2,dvara1,dvara2, &
                            ivarb1,ivarb2,rvarb1,rvarb2,dvarb1,dvarb2  )

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza, &
                         zt, xew, yew, zew, glatw, glonw
  use mem_ijtabs,  only: itab_w, itab_v
  use mem_leaf,    only: itab_wl, land
  use leaf_coms,   only: nwl, mwl, nzg
  use mem_sea,     only: itab_ws
  use sea_coms,    only: nws, mws

  implicit none

  integer,       intent(in) :: ndims
  integer,       intent(in) :: idims(3), jdims(3)
  character(*),  intent(in) :: varn ! Variable label
  character (2), intent(in) :: stagpt

  integer,  optional, intent(inout) :: ivara1(idims(1)), ivara2(idims(1),idims(2))
  real,     optional, intent(inout) :: rvara1(idims(1)), rvara2(idims(1),idims(2))
  real(r8), optional, intent(inout) :: dvara1(idims(1)), dvara2(idims(1),idims(2))
  integer,  optional, intent(inout) :: ivarb1(jdims(1)), ivarb2(jdims(1),jdims(2))
  real,     optional, intent(inout) :: rvarb1(jdims(1)), rvarb2(jdims(1),jdims(2))
  real(r8), optional, intent(inout) :: dvarb1(jdims(1)), dvarb2(jdims(1),jdims(2))

  integer :: iw, iv, iw1, iw2, iwl, iws, k, j, iwn, ngrw
  integer :: iw_og, iw_og1, iw_og2, iw_og3, iv_og, iwl_og, iws_og, ka_og
  real :: f_og1, f_og2, f_og3, cof

! Check stagpt to identify field type

  if (stagpt == 'AW') then

! Apply lower boundary condition to OLD grid fields in case they have not been
! defined below LPW level

     do iw_og = 2,nwa_og
        ka_og = lpw_og(iw_og)

        if (present(ivara2)) then
           ivara2(1:ka_og-1,iw_og) = ivara2(ka_og,iw_og)
        elseif (present(rvara2)) then
           rvara2(1:ka_og-1,iw_og) = rvara2(ka_og,iw_og)
        elseif (present(dvara2)) then
           if (trim(varn) == 'RHO'        .or. &
               trim(varn) == 'PRESS'      .or. &
               trim(varn) == 'PRESS_ACCUM') then

              ! Extrapolate density and pressure below ground level to give
              ! values roughly in accord with hydrostatic balance

              do k = ka_og-1,2,-1
                 cof = (zt(k) - zt(k+1)) / (zt(k+1) - zt(k+2))
                 dvara2(k,iw_og) = dvara2(k+1,iw_og) * (1. + cof) &
                                 - dvara2(k+2,iw_og) * cof
              enddo
           else
              dvara2(1:ka_og-1,iw_og) = dvara2(ka_og,iw_og)
           endif
        endif
     enddo

! Loop over all IW points on NEW grid

     do iw = 2,mwa

! Check ngr value of current IW point

        ngrw = 0
        do j = 1,itab_w(iw)%npoly
           iwn = itab_w(iw)%iw(j)

           if (ngrw < itab_w(iwn)%ngr) ngrw = itab_w(iwn)%ngr
        enddo

        if (itab_w(iw)%ngr <= ngrids_og .and. ngrw <= ngrids_og) then

! This IW point should be identical between OLD and NEW grids, so set its
! OLD grid index to NEW grid iwglobe value

           iw_og = itab_w(iw)%iwglobe

! Compare IW coordinates of OLD and NEW grids to make sure they are identical

           if (abs(xew(iw) - xew_og(iw_og)) > 1. .or. &
               abs(yew(iw) - yew_og(iw_og)) > 1. .or. &
               abs(zew(iw) - zew_og(iw_og)) > 1.) then

              print*, 'interp_addgrid: copy points not at same location'
              print*, 'iw   ',iw,xew(iw),yew(iw),zew(iw)
              print*, 'iwog ',iw_og,xew_og(iw_og),yew_og(iw_og),zew_og(iw_og)
              stop 'stop interp_addgrid'
           endif

! Copy values from OLD grid to NEW grid

           if (present(ivara1)) then
              ivarb1(iw) = ivara1(iw_og)
           elseif (present(ivara2)) then
              ivarb2(:,iw) = ivara2(:,iw_og)
           elseif (present(rvara1)) then
              rvarb1(iw) = rvara1(iw_og)
           elseif (present(rvara2)) then
              rvarb2(:,iw) = rvara2(:,iw_og)
           elseif (present(dvara1)) then
              dvarb1(iw) = dvara1(iw_og)
           elseif (present(dvara2)) then
              dvarb2(:,iw) = dvara2(:,iw_og)
           endif

        else

! This IW point is located within new mesh refinement so its values must be
! interpolated from a trio of OLD grid points

           iw_og1 = itab_wadd(iw)%iw_og(1); f_og1 = itab_wadd(iw)%f_og(1)
           iw_og2 = itab_wadd(iw)%iw_og(2); f_og2 = itab_wadd(iw)%f_og(2)
           iw_og3 = itab_wadd(iw)%iw_og(3); f_og3 = itab_wadd(iw)%f_og(3)

           if (present(ivara1)) then
              ivarb1(iw) = ivara1(iw_og1)
           elseif (present(ivara2)) then
              ivarb2(:,iw) = ivara2(:,iw_og1)
           elseif (present(rvara1)) then
              rvarb1(iw) = f_og1 * rvara1(iw_og1) &
                         + f_og2 * rvara1(iw_og2) &
                         + f_og3 * rvara1(iw_og3)
           elseif (present(rvara2)) then
              rvarb2(:,iw) = f_og1 * rvara2(:,iw_og1) &
                           + f_og2 * rvara2(:,iw_og2) &
                           + f_og3 * rvara2(:,iw_og3)
           elseif (present(dvara1)) then
              dvarb1(iw) = f_og1 * dvara1(iw_og1) &
                         + f_og2 * dvara1(iw_og2) &
                         + f_og3 * dvara1(iw_og3)
           elseif (present(dvara2)) then
              dvarb2(:,iw) = f_og1 * dvara2(:,iw_og1) &
                           + f_og2 * dvara2(:,iw_og2) &
                           + f_og3 * dvara2(:,iw_og3)
           endif

        endif

     enddo

  elseif (stagpt == 'AV') then

! The following fields are initialized elsewhere and should not be copied
! or interpolated here.

     do iv = 2,mva
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

! Check ngr values of IW points that are adjacent to this IV point

        if (itab_w(iw1)%ngr <= ngrids_og .and. itab_w(iw2)%ngr <= ngrids_og) then

! This IV point should be identical between OLD and NEW grids, so set its OLD
! grid index to NEW grid ivglobe value.  Copy values from OLD grid to NEW grid.

           iv_og = itab_v(iv)%ivglobe

           if (present(rvara2)) then
              rvarb2(:,iv) = rvara2(:,iv_og)
           elseif (present(dvara2)) then
              dvarb2(:,iv) = dvara2(:,iv_og)
           endif

        endif

     enddo

  elseif (stagpt == 'LW') then

     do iwl = 2,mwl
        iw = itab_wl(iwl)%iw

! Copy values to this NEW grid land cell from the nearest land cell of the
! OLD grid.

        iwl_og = itab_wladd(iwl)%iwl_og

        if (present(ivara1)) then
           ivarb1(iwl) = ivara1(iwl_og)
        elseif (present(ivara2)) then
           ivarb2(:,iwl) = ivara2(:,iwl_og)
        elseif (present(rvara1)) then
           rvarb1(iwl) = rvara1(iwl_og)
        elseif (present(rvara2)) then
           rvarb2(:,iwl) = rvara2(:,iwl_og)
        elseif (present(dvara1)) then
           dvarb1(iwl) = dvara1(iwl_og)
        elseif (present(dvara2)) then
           dvarb2(:,iwl) = dvara2(:,iwl_og)
        endif

     enddo

  elseif (stagpt == 'SW') then

     do iws = 2,mws
        iw = itab_ws(iws)%iw

! Copy values to this NEW grid land cell from the nearest land cell of the
! OLD grid.

        iws_og = itab_wsadd(iws)%iws_og

        if (present(ivara1)) then
           ivarb1(iws) = ivara1(iws_og)
        elseif (present(ivara2)) then
           ivarb2(:,iws) = ivara2(:,iws_og)
        elseif (present(rvara1)) then
           rvarb1(iws) = rvara1(iws_og)
        elseif (present(rvara2)) then
           rvarb2(:,iws) = rvara2(:,iws_og)
        elseif (present(dvara1)) then
           dvarb1(iws) = dvara1(iws_og)
        elseif (present(dvara2)) then
           dvarb2(:,iws) = dvara2(:,iws_og)
        endif

     enddo

  endif

  end subroutine interp_addgrid

!===============================================================================

  subroutine vel_t3d_hex_oldgrid()

  use misc_coms, only: io6

  implicit none

  integer :: iw_og, npoly, kb_og, k_og, jv, iv_og, ksw
  real    :: wst

! Horizontal loop over W columns

  do iw_og = 2,nwa_og

     npoly = itab_wog(iw_og)%npoly
     kb_og = lpw_og(iw_og)

! Vertical loop over T levels

     do k_og = kb_og, nza_og

! Diagnose 3D earth-velocity vector at T points; W contribution first

        if (k_og == kb_og) then
           wst = 0.5 * wc_og(k_og,iw_og)
        else
           wst = 0.5 * (wc_og(k_og-1,iw_og) + wc_og(k_og,iw_og))
        endif

        vxe_og(k_og,iw_og) = wst * wnx_og(iw_og)
        vye_og(k_og,iw_og) = wst * wny_og(iw_og)
        vze_og(k_og,iw_og) = wst * wnz_og(iw_og)

     enddo

! Effective contribution from submerged V faces

   do ksw = 1, lve2_og(iw_og)
      k_og = kb_og + ksw - 1

      vxe_og(k_og,iw_og) = vxe_og(k_og,iw_og) + vxe2_og(ksw,iw_og) 
      vye_og(k_og,iw_og) = vye_og(k_og,iw_og) + vye2_og(ksw,iw_og) 
      vze_og(k_og,iw_og) = vze_og(k_og,iw_og) + vze2_og(ksw,iw_og) 
   enddo

! Loop over V neighbors of this W cell

     do jv = 1, npoly

        iv_og = itab_wog(iw_og)%iv(jv)

! Vertical loop over T levels

        do k_og = kb_og, nza_og

! Diagnose 3D earth-velocity vector at T points; VC contribution

           vxe_og(k_og,iw_og) = vxe_og(k_og,iw_og) &
                              + itab_wog(iw_og)%ecvec_vx(jv) * vc_og(k_og,iv_og)
           vye_og(k_og,iw_og) = vye_og(k_og,iw_og) &
                              + itab_wog(iw_og)%ecvec_vy(jv) * vc_og(k_og,iv_og)
           vze_og(k_og,iw_og) = vze_og(k_og,iw_og) &
                              + itab_wog(iw_og)%ecvec_vz(jv) * vc_og(k_og,iv_og)

        enddo
      
     enddo

     vxe_og(1:kb_og-1,iw_og) = vxe_og(kb_og,iw_og)
     vye_og(1:kb_og-1,iw_og) = vye_og(kb_og,iw_og)
     vze_og(1:kb_og-1,iw_og) = vze_og(kb_og,iw_og)

  enddo

  end subroutine vel_t3d_hex_oldgrid


!===============================================================================


  subroutine diagvel_t3d_init_addgrid()

  use mem_basic,  only: wc, vc, vxe2, vye2, vze2
  use mem_ijtabs, only: jtab_w, itab_v, itab_w, jtw_prog
  use mem_grid,   only: lpw, lve2, lpv, wnxo2, wnyo2, wnzo2
  use misc_coms,  only: io6

  implicit none

  integer :: j,iw,npoly,ka,k,jv,iv,ksw,kbv,iw_og, jn, iwn, ngrw

! Horizontal loop over W columns

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     npoly = itab_w(iw)%npoly
     ka = lpw(iw)

     vxe2(:,iw) = 0.
     vye2(:,iw) = 0.
     vze2(:,iw) = 0.

! Check ngr value of current IW point

     ngrw = 0
     do jn = 1,npoly
        iwn = itab_w(iw)%iw(jn)

        if (ngrw < itab_w(iwn)%ngr) ngrw = itab_w(iwn)%ngr
     enddo

     if (itab_w(iw)%ngr <= ngrids_og .and. ngrw <= ngrids_og) then

! This IW point should be identical between OLD and NEW grids, so set its
! OLD grid index to NEW grid iwglobe value

        iw_og = itab_w(iw)%iwglobe

! Compare IW coordinates of OLD and NEW grids to make sure they are identical

        if ( lve2_og(iw_og) /= lve2(iw) ) then

           print*, 'diagvel_t3d_init_addgrid: grid changed'
           print*, 'lve2   ',iw,lve2(iw)
           print*, 'lve2_og',iw_og,lve2_og(iw_og)
           stop 'stop diagvel_t3d_init_addgrid'
        endif

! Copy values from OLD grid to NEW grid

        vxe2(1:lve2(iw),iw) = vxe2_og(1:lve2_og(iw_og),iw_og)
        vye2(1:lve2(iw),iw) = vye2_og(1:lve2_og(iw_og),iw_og)
        vze2(1:lve2(iw),iw) = vze2_og(1:lve2_og(iw_og),iw_og)

     else

! This IW point is located within new mesh refinement, so vxe2, vye2, and vze2
! must be reinitialized from the local VC and WC winds

        vxe2(1,iw) = wnxo2(iw) * wc(ka,iw)
        vye2(1,iw) = wnyo2(iw) * wc(ka,iw)
        vze2(1,iw) = wnzo2(iw) * wc(ka,iw)

! Loop over adjacent V faces

        do jv = 1, npoly

           iv  = itab_w(iw)%iv(jv)
           kbv = lpv(iv)

! Check if any V faces are below ground

           if (ka < kbv) then
              do k = ka, kbv-1
                 ksw = k - ka + 1

! Project INITIAL VC from below-ground V faces back to (vxe2, vye2, vze2)

                 vxe2(ksw,iw) = vxe2(ksw,iw) + itab_w(iw)%ecvec_vx(jv) * vc(kbv,iv)
                 vye2(ksw,iw) = vye2(ksw,iw) + itab_w(iw)%ecvec_vy(jv) * vc(kbv,iv)
                 vze2(ksw,iw) = vze2(ksw,iw) + itab_w(iw)%ecvec_vz(jv) * vc(kbv,iv)
              enddo
           endif

        enddo

     endif

  enddo

  end subroutine diagvel_t3d_init_addgrid

End Module mem_addgrid
