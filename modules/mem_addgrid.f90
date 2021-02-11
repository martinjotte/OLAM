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
  integer :: nland_og
  integer :: nsea_og

  integer, allocatable :: lpw_og(:)  ! (nwa_og)
  integer, allocatable :: lve2_og(:) ! (nwa_og)

  real, allocatable :: wnx_og(:) ! (nwa_og)
  real, allocatable :: wny_og(:) ! (nwa_og)
  real, allocatable :: wnz_og(:) ! (nwa_og)

  real, allocatable :: xew_og(:) ! (nwa_og)
  real, allocatable :: yew_og(:) ! (nwa_og)
  real, allocatable :: zew_og(:) ! (nwa_og)

  real, allocatable :: xewsfc_og(:) ! (nsfc_og)
  real, allocatable :: yewsfc_og(:) ! (nsfc_og)
  real, allocatable :: zewsfc_og(:) ! (nsfc_og)

  real, allocatable :: soil_water_og (:,:) ! (nzg_og,nland_og)
  real, allocatable :: soil_energy_og(:,:) ! (nzg_og,nland_og)

  real, allocatable :: specifheat_drysoil_og(:,:) ! (nzg_og,nland_og)

  integer, allocatable :: leaf_class_og(:) ! (nland_og)

  Type itab_wog_vars          ! data structure for OLD grid IW pts (global domain)
     integer :: npoly = 0     ! number of W neighbors of this W pt
     integer :: iv(7)         ! neighbor V pts
     integer :: iw(7)         ! neighbor W pts
  End Type itab_wog_vars

  Type itab_wadd_vars         ! data structure for WADD pts (individual rank)
     integer :: iw_og(3) = 1  ! intrp pts
     real    :: f_og(3) = 0.  ! intrp coeffs
  End Type itab_wadd_vars

  Type itab_landadd_vars        ! data structure for landADD pts (individual rank)
     integer :: iland_og = 1  ! indes of nearest land point on OLD grid
  End Type itab_landadd_vars

  Type itab_seaadd_vars      ! data structure for seaADD pts (individual rank)
     integer :: isea_og = 1  ! index of nearest sea point on OLD grid
  End Type itab_seaadd_vars

  type (itab_wog_vars),    allocatable, target :: itab_wog(:)   ! (dimension nwa_og)
  type (itab_wadd_vars),   allocatable, target :: itab_wadd(:)  ! (dimension mwa)
  type (itab_landadd_vars),  allocatable, target :: itab_landadd(:) ! (dimension mland)
  type (itab_seaadd_vars),  allocatable, target :: itab_seaadd(:) ! (dimension msea)

Contains

!==============================================================================

  subroutine init_addgrid()

  use misc_coms,  only: io6, mdomain, runtype, nxp
  use mem_ijtabs, only: itab_w
  use mem_sfcg,    only: sfcg
  use mem_grid,   only: mwa, xew, yew, zew, glatw, glonw
  use leaf_coms,  only: isfcl
  use mem_land,   only: land, mland
  use mem_sea,    only: sea, msea, omsea

  implicit none

  real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
  real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
  real :: dist, dist_min, dist_max, xi, yi
  real :: xin(7),yin(7)

  integer :: iw, iw_og, iwn, iwn_og, npoly, j, jmaxneg, jminpos, ngrw
  integer :: iland, isea, iwsfc, iland_og, isea_og, isfc_og

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

  allocate(itab_landadd(mland))
  allocate(itab_seaadd(msea))

! Loop through all NEW grid land cells, find the nearest land cell on the
! OLD grid, and attach its index to the NEW grid land cell.

     !$omp parallel do private(iwsfc,dist_min,iland_og,isfc_og,dist)
     do iland = 2,mland
        iwsfc = iland + 0

        if (mod(iland,100000) == 0) write(io6,*) 'init_addgrid iland,mland ',iland,mland

! Initialize minimum distance variable

        dist_min = 1.e12

! Loop over all local-domain iland points in OLD grid

        do iland_og = 2,nland_og
           isfc_og = iland_og

! Skip interaction if NEW and OLD grid points are more than dist_min apart
! in any earth-frame direction

           if (abs(sfcg%xew(iwsfc) - xewsfc_og(isfc_og)) > dist_min) cycle
           if (abs(sfcg%yew(iwsfc) - yewsfc_og(isfc_og)) > dist_min) cycle
           if (abs(sfcg%zew(iwsfc) - zewsfc_og(isfc_og)) > dist_min) cycle

! Compute distance between NEW and OLD grid iland points, and reset minimum
! distance and HISTORY grid point index if new minimum is found

           dist = sqrt((sfcg%xew(iwsfc) - xewsfc_og(isfc_og))**2 &
                     + (sfcg%yew(iwsfc) - yewsfc_og(isfc_og))**2 &
                     + (sfcg%zew(iwsfc) - zewsfc_og(isfc_og))**2)

           if (dist_min > dist) then
              dist_min = dist
              itab_landadd(iland)%iland_og = iland_og
           endif

        enddo  ! iland_og

     enddo  ! iland
     !$omp end parallel do

! Loop through all NEW grid sea cells, find the nearest sea cell on the
! OLD grid, and attach its index to the NEW grid sea cell.

     do isea = 2,msea
        iwsfc = isea + omsea

        if (mod(isea,100000) == 0) write(io6,*) 'init_addgrid isea,msea ',isea,msea

! Initialize minimum distance variable

        dist_min = 1.e12

! Loop over all local-domain isea points in OLD grid

        do isea_og = 2,nsea_og

! Skip interaction if NEW and OLD grid points are more than dist_min apart
! in any earth-frame direction

          if (abs(sfcg%xew(iwsfc) - xewsfc_og(isfc_og)) > dist_min) cycle
          if (abs(sfcg%yew(iwsfc) - yewsfc_og(isfc_og)) > dist_min) cycle
          if (abs(sfcg%zew(iwsfc) - zewsfc_og(isfc_og)) > dist_min) cycle

! Compute distance between NEW and OLD grid isea points, and reset minimum
! distance and HISTORY grid point index if new minimum is found

          dist = sqrt((sfcg%xew(iwsfc) - xewsfc_og(isfc_og))**2 &
                    + (sfcg%yew(iwsfc) - yewsfc_og(isfc_og))**2 &
                    + (sfcg%zew(iwsfc) - zewsfc_og(isfc_og))**2)

          if (dist_min > dist) then
             dist_min = dist
             itab_seaadd(isea)%isea_og = isea_og
          endif

       enddo  ! isea_og

     enddo  ! isea

  endif  ! isfcl = 1

  end subroutine init_addgrid

!=========================================================================

  subroutine interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            ivara1,ivara2,ivara3,rvara1,rvara2,rvara3,dvara1,dvara2,dvara3, &
                            ivarb1,ivarb2,ivarb3,rvarb1,rvarb2,rvarb3,dvarb1,dvarb2,dvarb3  )

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza, &
                         zt, xew, yew, zew, glatw, glonw
  use mem_ijtabs,  only: itab_w, itab_v
  use mem_land,    only: land, nland, mland, nzg
  use mem_sea,     only: nsea, msea

  implicit none

  integer,       intent(in) :: ndims
  integer,       intent(in) :: idims(3), jdims(3)
  character(*),  intent(in) :: varn ! Variable label
  character (2), intent(in) :: stagpt

  integer,  optional, intent(inout) :: ivara1(idims(1))
  integer,  optional, intent(inout) :: ivara2(idims(1),idims(2))
  integer,  optional, intent(inout) :: ivara3(idims(1),idims(2),idims(3))

  real,     optional, intent(inout) :: rvara1(idims(1))
  real,     optional, intent(inout) :: rvara2(idims(1),idims(2))
  real,     optional, intent(inout) :: rvara3(idims(1),idims(2),idims(3))

  real(r8), optional, intent(inout) :: dvara1(idims(1))
  real(r8), optional, intent(inout) :: dvara2(idims(1),idims(2))
  real(r8), optional, intent(inout) :: dvara3(idims(1),idims(2),idims(3))

  integer,  optional, intent(inout) :: ivarb1(idims(1))
  integer,  optional, intent(inout) :: ivarb2(idims(1),idims(2))
  integer,  optional, intent(inout) :: ivarb3(idims(1),idims(2),idims(3))

  real,     optional, intent(inout) :: rvarb1(idims(1))
  real,     optional, intent(inout) :: rvarb2(idims(1),idims(2))
  real,     optional, intent(inout) :: rvarb3(idims(1),idims(2),idims(3))

  real(r8), optional, intent(inout) :: dvarb1(idims(1))
  real(r8), optional, intent(inout) :: dvarb2(idims(1),idims(2))
  real(r8), optional, intent(inout) :: dvarb3(idims(1),idims(2),idims(3))

  integer :: iw, iv, iw1, iw2, iland, isea, k, j, iwn, ngrw
  integer :: iw_og, iw_og1, iw_og2, iw_og3, iv_og, iland_og, isea_og, ka_og
  real :: f_og1, f_og2, f_og3, cof

! Check stagpt to identify field type

  if (stagpt == 'AW') then

! Apply lower boundary condition to OLD grid fields in case they have not been
! defined below LPW level

     do iw_og = 2,nwa_og
        ka_og = lpw_og(iw_og)

        if     (present(ivara2)) then
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
        elseif (present(ivara3)) then
           do k = 1, ka_og-1
              ivara3(k,:,iw_og) = ivara3(ka_og,:,iw_og)
           enddo
        elseif (present(rvara3)) then
           do k = 1, ka_og-1
              rvara3(k,:,iw_og) = rvara3(ka_og,:,iw_og)
           enddo
        elseif (present(dvara3)) then
           do k = 1, ka_og-1
              dvara3(k,:,iw_og) = dvara3(ka_og,:,iw_og)
           enddo
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
           elseif (present(ivara3)) then
              ivarb3(:,:,iw) = ivara3(:,:,iw_og1)

           elseif (present(rvara1)) then
              rvarb1(iw) = f_og1 * rvara1(iw_og1) &
                         + f_og2 * rvara1(iw_og2) &
                         + f_og3 * rvara1(iw_og3)
           elseif (present(rvara2)) then
              rvarb2(:,iw) = f_og1 * rvara2(:,iw_og1) &
                           + f_og2 * rvara2(:,iw_og2) &
                           + f_og3 * rvara2(:,iw_og3)
           elseif (present(rvara3)) then
              rvarb3(:,:,iw) = f_og1 * rvara3(:,:,iw_og1) &
                             + f_og2 * rvara3(:,:,iw_og2) &
                             + f_og3 * rvara3(:,:,iw_og3)

           elseif (present(dvara1)) then
              dvarb1(iw) = f_og1 * dvara1(iw_og1) &
                         + f_og2 * dvara1(iw_og2) &
                         + f_og3 * dvara1(iw_og3)
           elseif (present(dvara2)) then
              dvarb2(:,iw) = f_og1 * dvara2(:,iw_og1) &
                           + f_og2 * dvara2(:,iw_og2) &
                           + f_og3 * dvara2(:,iw_og3)
           elseif (present(dvara3)) then
              dvarb3(:,:,iw) = f_og1 * dvara3(:,:,iw_og1) &
                             + f_og2 * dvara3(:,:,iw_og2) &
                             + f_og3 * dvara3(:,:,iw_og3)
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

     do iland = 2,mland

! Copy values to this NEW grid land cell from the nearest land cell of the
! OLD grid.

        iland_og = itab_landadd(iland)%iland_og

        if (present(ivara1)) then
           ivarb1(iland) = ivara1(iland_og)
        elseif (present(ivara2)) then
           ivarb2(:,iland) = ivara2(:,iland_og)
        elseif (present(rvara1)) then
           rvarb1(iland) = rvara1(iland_og)
        elseif (present(rvara2)) then
           rvarb2(:,iland) = rvara2(:,iland_og)
        elseif (present(dvara1)) then
           dvarb1(iland) = dvara1(iland_og)
        elseif (present(dvara2)) then
           dvarb2(:,iland) = dvara2(:,iland_og)
        endif

     enddo

  elseif (stagpt == 'SW') then

     do isea = 2,msea

! Copy values to this NEW grid land cell from the nearest land cell of the
! OLD grid.

        isea_og = itab_seaadd(isea)%isea_og

        if (present(ivara1)) then
           ivarb1(isea) = ivara1(isea_og)
        elseif (present(ivara2)) then
           ivarb2(:,isea) = ivara2(:,isea_og)
        elseif (present(rvara1)) then
           rvarb1(isea) = rvara1(isea_og)
        elseif (present(rvara2)) then
           rvarb2(:,isea) = rvara2(:,isea_og)
        elseif (present(dvara1)) then
           dvarb1(isea) = dvara1(isea_og)
        elseif (present(dvara2)) then
           dvarb2(:,isea) = dvara2(:,isea_og)
        endif

     enddo

  endif

  end subroutine interp_addgrid

End Module mem_addgrid
