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

subroutine history_start(action)

  ! This routine initializes the model from the history file

  use misc_coms,  only: io6, hfilin, time8, time_istp8, runtype, &
                        iparallel, isubdomain
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_para,   only: myrank, mgroupsize
  use max_dims,   only: pathlen

  implicit none

  character(*), intent(in) :: action
  logical                  :: exans

  ! List of all allowable runtype, seq-para run, histfile read combinations:

  ! 'HISTORY'     , sequential run, sequential histfile  (isubdomain = 0)
  ! 'HISTORY'     , parallel   run, sequential histfile  (isubdomain = 1)
  ! 'HISTADDGRID' , sequential run, sequential histfile  (isubdomain = 0)
  ! 'HISTADDGRID' , parallel   run, sequential histfile  (isubdomain = 1)
  ! 'PLOTONLY'    , sequential run, sequential histfile  (isubdomain = 0)
  ! 'PLOTONLY'    , parallel   run, sequential histfile  (isubdomain = 1)

  ! Check if history files exist
  
  inquire(file=hfilin, exist=exans)   ! global restart file

  if (exans) then

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening history file '//trim(hfilin)//' for '//trim(action)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(hfilin,'R')

     if (trim(action) == 'COMMIO') then

        ! Read the common variables
        call commio('READ')

     else

        ! Read the model fields
        call shdf5_irec(1, (/1/), 'time8', dvars=time8)

        if (runtype == 'HISTADDGRID') then

! If adding new grids (new mesh refinements), read model fields from old
! history file and map them to new model grid

           print*, 'calling hist_read_addgrid'
           call hist_read_addgrid()

        else

! If NOT adding new grids (new mesh refinements), read model fields from
! standard history file

           print*, 'calling hist_read'
           call hist_read()
        endif

        write(io6,*) 'finished history read'

     endif

     time_istp8 = time8  ! time_istp8 is used for plotting time
     call shdf5_close()

  else

     ! History files do not exist, stop model.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open history file file:'
     write(io6,*) '!!!   '//trim(hfilin)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in history_start'

  endif

end subroutine history_start

!=========================================================================

subroutine hist_read()

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza
  use mem_ijtabs,  only: itab_w, itab_v, itab_m
  use misc_coms,   only: io6, runtype
  use var_tables,  only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils,  only: shdf5_info, shdf5_irec
  use mem_leaf,    only: itab_wl
  use leaf_coms,   only: nwl, mwl, nzg
  use mem_sea,     only: itab_ws
  use sea_coms,    only: nws, mws
  use mem_para,    only: myrank
  use consts_coms, only: r8
  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  implicit none

  integer          :: nv, nvcnt, ndims, idims(3)
  character(32)    :: varn
  character (2)    :: stagpt
  integer, pointer :: ilocal(:)

  nvcnt =  0

  ! Loop through all variables in the model vtables

  do nv = 1, num_var

     ndims = -1
     idims =  0

     ! Skip to next variable if we don't want the current one

     if (.not. vtab_r(nv)%ihist) cycle
     if (vtab_r(nv)%nread) cycle

     varn    = trim(vtab_r(nv)%name)
     stagpt  = vtab_r(nv)%stagpt

     ! We want it...read it if it's in the history file

     call shdf5_info(varn, ndims, idims)

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        write(io6,*)
        write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
        write(io6,*)
        cycle
     endif

     ! Identify the points we want to read from the history file

     if     (stagpt == 'AW' .and. idims(ndims) == nwa) then
        ilocal => itab_w(:)%iwglobe
        idims(ndims) = mwa
     elseif (stagpt == 'AV' .and. idims(ndims) == nva) then
        ilocal => itab_v(:)%ivglobe
        idims(ndims) = mva
     elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
        ilocal => itab_m(:)%imglobe
        idims(ndims) = mma
     elseif (stagpt == 'LW' .and. idims(ndims) == nwl) then
        ilocal => itab_wl(:)%iwglobe
        idims(ndims) = mwl
     elseif (stagpt == 'SW' .and. idims(ndims) == nws) then
        ilocal => itab_ws(:)%iwglobe
        idims(ndims) = mws
     elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
        ilocal => itab_wnud(:)%iwnudglobe
        idims(ndims) = mwnud
     else

        ! TODO: Const values
        ! TODO: Land U and M; Sea U and M? (probably not!)

        stop "invalid array size in history_read"
     endif

     if     (associated(vtab_r(nv)%ivar1_p)) then
        call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%ivar2_p)) then
        call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar2_p, points=ilocal)
     elseif (associated(vtab_r(nv)%rvar1_p)) then
        call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%rvar2_p)) then
        call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar2_p, points=ilocal)
     elseif (associated(vtab_r(nv)%dvar1_p)) then
        call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%dvar2_p)) then
        call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar2_p, points=ilocal)
     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

end subroutine hist_read

!=========================================================================

subroutine hist_read_addgrid()

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza, lpv, &
                         vnx, vny, vnz, lpw
  use mem_ijtabs,  only: itab_w, itab_v, itab_m, jtab_v, jtv_prog
  use misc_coms,   only: io6, runtype, iparallel
  use mem_basic,   only: vc, vp, vmc, vmp, wc, wmc, rho, vxe, vye, vze
  use var_tables,  only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils,  only: shdf5_info, shdf5_irec
  use mem_leaf,    only: itab_wl, land
  use leaf_coms,   only: mwl, nzg, &
                         slzt, slcpd, tai_max, soil_rough, snow_rough
  use leaf4_soil,  only: soil_wat2pot, soil_pot2wat
  use mem_sea,     only: itab_ws
  use sea_coms,    only: nws, mws
  use mem_para,    only: myrank
  use consts_coms, only: r8, cliq1000, cice1000, alli1000
  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  use mem_flux_accum, only: vc_accum

  use mem_addgrid, only: interp_addgrid, vel_t3d_hex_oldgrid, &
                         nza_og, nwa_og, nva_og, nma_og, nwl_og, nws_og, &
                         nzg_og, ngrids_og, lpw_og, itab_wladd, &
                         vc_og, wc_og, vc_accum_og, vxe_og, vye_og, vze_og, &
                         soil_water_og, soil_energy_og, ntext_soil_og, &
                         leaf_class_og, &
                         vxe2_og, vye2_og, vze2_og, nve2_max_og, &
                         diagvel_t3d_init_addgrid

  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

  implicit none

  integer :: j, iw, iw1, iw2, iv, nv, kb, k, mrl
  integer :: nvcnt, ndims, idims(3), jdims(3)
  character(32) :: varn
  character (2) :: stagpt

  integer :: iwl, iwl_og, nts, nts_og, kb_og, iw_og, k_og
  real :: water_frac_ul, water_frac, psi, head, headp
  real :: soil_tempk, soil_tempc, soil_fracliq, tai, tai_og

! Scratch arrays for copying input

  integer,  allocatable :: iscr1(:)
  integer,  allocatable :: iscr2(:,:)
  real,     allocatable :: rscr1(:)
  real,     allocatable :: rscr2(:,:)
  real(r8), allocatable :: dscr1(:)
  real(r8), allocatable :: dscr2(:,:)

! Allocate temporary OLD grid velocity [and momentum] arrays

  allocate (vc_og(nza_og,nva_og))
  allocate (wc_og(nza_og,nwa_og))

  allocate (vc_accum_og(nza_og,nva_og))

  allocate (soil_water_og (nzg_og,nwl_og))
  allocate (soil_energy_og(nzg_og,nwl_og))

  allocate (vxe2_og(nve2_max_og,nwa_og))
  allocate (vye2_og(nve2_max_og,nwa_og))
  allocate (vze2_og(nve2_max_og,nwa_og))

  nvcnt =  0

  ! Loop through all variables in the model vtables

  do nv = 1, num_var

     ndims = -1
     idims =  0

     ! Skip to next variable if we don't want the current one

     if (.not. vtab_r(nv)%ihist) cycle
     if (vtab_r(nv)%nread) cycle

     varn    = trim(vtab_r(nv)%name)
     stagpt  = vtab_r(nv)%stagpt

! In a HISTADDGRID restart, the following fields need not or should not be read
! from the OLD grid history file.  They are diagnostic or auxiliary variables,
! and are initialized or computed elsewhere, no later than on the first
! timestep of the history restart.

     if (trim(varn) == 'THSRC'               .or. &
         trim(varn) == 'RTSRC'               .or. &
         trim(varn) == 'QWCON'               .or. &
         trim(varn) == 'CBMF'                .or. &
         trim(varn) == 'KCUTOP'              .or. &
         trim(varn) == 'KCUBOT'              .or. &
         trim(varn) == 'VXSRC'               .or. &
         trim(varn) == 'VYSRC'               .or. &
         trim(varn) == 'VZSRC'               .or. &
         trim(varn) == 'FTHRD_SW'            .or. &
         trim(varn) == 'FTHRD_LW'            .or. &
         trim(varn) == 'CLOUD_FRAC'          .or. &
         trim(varn) == 'RSHORT'              .or. &
         trim(varn) == 'RLONG'               .or. &
         trim(varn) == 'RLONGUP'             .or. &
         trim(varn) == 'RSHORT_TOP'          .or. &
         trim(varn) == 'RSHORTUP_TOP'        .or. &
         trim(varn) == 'RLONGUP_TOP'         .or. &
         trim(varn) == 'ALBEDT'              .or. &
         trim(varn) == 'COSZ'                .or. &
         trim(varn) == 'HKM'                 .or. &
         trim(varn) == 'SXFER_TK'            .or. &
         trim(varn) == 'SXFER_RK'            .or. &
         trim(varn) == 'VKM_SFC'             .or. &
         trim(varn) == 'SFLUXT'              .or. &
         trim(varn) == 'SFLUXR'              .or. &
         trim(varn) == 'USTAR'               .or. &
         trim(varn) == 'WSTAR'               .or. &
         trim(varn) == 'WTV0'                .or. &
         trim(varn) == 'PBLH'                .or. &
         trim(varn) == 'KPBLH'               .or. &
         trim(varn) == 'FTHPBL'              .or. &
         trim(varn) == 'FQTPBL'              .or. &
         trim(varn) == 'LAND%RHOS'           .or. &
         trim(varn) == 'LAND%VELS'           .or. &
         trim(varn) == 'LAND%PRSS'           .or. &
         trim(varn) == 'LAND%AIRTEMP'        .or. &
         trim(varn) == 'LAND%AIRSHV'         .or. &
         trim(varn) == 'LAND%USTAR'          .or. &
         trim(varn) == 'LAND%VKMSFC'         .or. &
         trim(varn) == 'LAND%SFLUXT'         .or. &
         trim(varn) == 'LAND%SFLUXR'         .or. &
         trim(varn) == 'LAND%SFLUXC'         .or. &
         trim(varn) == 'LAND%SXFER_T'        .or. &
         trim(varn) == 'LAND%SXFER_R'        .or. &
         trim(varn) == 'LAND%SXFER_C'        .or. &
         trim(varn) == 'LAND%SXFER_T'        .or. &
         trim(varn) == 'LAND%SXFER_R'        .or. &
         trim(varn) == 'LAND%SXFER_C'        .or. &
         trim(varn) == 'LAND%ED_GGAER'       .or. &
         trim(varn) == 'LAND%ED_ZETA'        .or. &
         trim(varn) == 'LAND%ED_RIB'         .or. &
         trim(varn) == 'LAND%ALBEDO_BEAM'    .or. &
         trim(varn) == 'LAND%ALBEDO_DIFFUSE' .or. &
         trim(varn) == 'LAND%RSHORT'         .or. &
         trim(varn) == 'LAND%RSHORT_DIFFUSE' .or. &
         trim(varn) == 'LAND%RLONG'          .or. &
         trim(varn) == 'LAND%RLONG_ALBEDO'   .or. &
         trim(varn) == 'LAND%RLONGUP'        .or. &
         trim(varn) == 'LAND%RSHORT_G'       .or. &
         trim(varn) == 'LAND%RSHORT_S'       .or. &
         trim(varn) == 'LAND%RSHORT_V'       .or. &
         trim(varn) == 'LAND%RLONG_G'        .or. &
         trim(varn) == 'LAND%RLONG_S'        .or. &
         trim(varn) == 'LAND%RLONG_V'        .or. &
         trim(varn) == 'LAND%COSZ'           .or. &
         trim(varn) == 'LAND%PCPG'           .or. &
         trim(varn) == 'LAND%QPCPG'          .or. &
         trim(varn) == 'LAND%DPCPG'          .or. &
         trim(varn) == 'LAND%CAN_DEPTH'      .or. &
         trim(varn) == 'LAND%HEAD0'          .or. &
         trim(varn) == 'LAND%HCAPVEG'        .or. &
         trim(varn) == 'LAND%SURFACE_SSH'    .or. &
         trim(varn) == 'LAND%GROUND_SHV'     .or. &
         trim(varn) == 'LAND%ROUGH'          .or. &
         trim(varn) == 'LAND%VEG_FRACAREA'   .or. &
         trim(varn) == 'LAND%VEG_LAI'        .or. &
         trim(varn) == 'LAND%VEG_ROUGH'      .or. &
         trim(varn) == 'LAND%VEG_HEIGHT'     .or. &
         trim(varn) == 'LAND%VEG_ALBEDO'     .or. &
         trim(varn) == 'LAND%VEG_TAI'        .or. &
         trim(varn) == 'LAND%VEG_NDVIC'      .or. &
         trim(varn) == 'LAND%SNOWFAC'        .or. &
         trim(varn) == 'LAND%VF'             .or. &
         trim(varn) == 'SEA%RHOS'            .or. &
         trim(varn) == 'SEA%VELS'            .or. &
         trim(varn) == 'SEA%PRSS'            .or. &
         trim(varn) == 'SEA%AIRTEMP'         .or. &
         trim(varn) == 'SEA%AIRSHV'          .or. &
         trim(varn) == 'SEA%USTAR'           .or. &
         trim(varn) == 'SEA%SEA_USTAR'       .or. &
         trim(varn) == 'SEA%ICE_USTAR'       .or. &
         trim(varn) == 'SEA%VKMSFC'          .or. &
         trim(varn) == 'SEA%SEA_VKMSFC'      .or. &
         trim(varn) == 'SEA%ICE_VKMSFC'      .or. &
         trim(varn) == 'SEA%SFLUXT'          .or. &
         trim(varn) == 'SEA%SEA_SFLUXT'      .or. &
         trim(varn) == 'SEA%ICE_SFLUXT'      .or. &
         trim(varn) == 'SEA%SFLUXR'          .or. &
         trim(varn) == 'SEA%SEA_SFLUXR'      .or. &
         trim(varn) == 'SEA%ICE_SFLUXR'      .or. &
         trim(varn) == 'SEA%SFLUXC'          .or. &
         trim(varn) == 'SEA%SXFER_T'         .or. &
         trim(varn) == 'SEA%SEA_SXFER_T'     .or. &
         trim(varn) == 'SEA%ICE_SXFER_T'     .or. &
         trim(varn) == 'SEA%SXFER_R'         .or. &
         trim(varn) == 'SEA%SEA_SXFER_R'     .or. &
         trim(varn) == 'SEA%ICE_SXFER_R'     .or. &
         trim(varn) == 'SEA%SXFER_C'         .or. &
         trim(varn) == 'SEA%ALBEDO_BEAM'     .or. &
         trim(varn) == 'SEA%ALBEDO_DIFFUSE'  .or. &
         trim(varn) == 'SEA%SEA_ALBEDO'      .or. &
         trim(varn) == 'SEA%ICE_ALBEDO'      .or. &
         trim(varn) == 'SEA%RSHORT'          .or. &
         trim(varn) == 'SEA%RSHORT_DIFFUSE'  .or. &
         trim(varn) == 'SEA%RLONG'           .or. &
         trim(varn) == 'SEA%RLONG_ALBEDO'    .or. &
         trim(varn) == 'SEA%RLONGUP'         .or. &
         trim(varn) == 'SEA%SEA_RLONGUP'     .or. &
         trim(varn) == 'SEA%ICE_RLONGUP'     .or. &
         trim(varn) == 'SEA%ICE_NET_RLONG'   .or. &
         trim(varn) == 'SEA%ICE_NET_RSHORT'  .or. &
         trim(varn) == 'SEA%PCPG'            .or. &
         trim(varn) == 'SEA%QPCPG'           .or. &
         trim(varn) == 'SEA%DPCPG'           .or. &
         trim(varn) == 'SEA%CAN_DEPTH'       .or. &
         trim(varn) == 'SEA%SEATC'           .or. &
         trim(varn) == 'SEA%SURFACE_SSH'     .or. &
         trim(varn) == 'SEA%SEA_SFC_SSH'     .or. &
         trim(varn) == 'SEA%ICE_SFC_SSH')    cycle

     ! We want it...read it if it's in the history file

     call shdf5_info(varn, ndims, idims)

     print*, 'varn ',varn

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        write(io6,*)
        write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
        write(io6,*)
        cycle
     endif

     jdims = idims

     ! Identify the points we want to read from the history file

     if     (stagpt == 'AW' .and. idims(ndims) == nwa_og) then
        jdims(ndims) = mwa
     elseif (stagpt == 'AV' .and. idims(ndims) == nva_og) then
        jdims(ndims) = mva
     elseif (stagpt == 'AM' .and. idims(ndims) == nma_og) then
        jdims(ndims) = mma
     elseif (stagpt == 'LW' .and. idims(ndims) == nwl_og) then
        jdims(ndims) = mwl
     elseif (stagpt == 'SW' .and. idims(ndims) == nws_og) then
        jdims(ndims) = mws
     elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
        jdims(ndims) = mwnud
     else

        ! TODO: Const values
        ! TODO: Land U and M; Sea U and M? (probably not!)

        stop "invalid array size in history_read"
     endif

     if     (associated(vtab_r(nv)%ivar1_p)) then
        allocate(iscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, ivara=iscr1)
        call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            ivara1=iscr1, ivarb1=vtab_r(nv)%ivar1_p)
        deallocate(iscr1)
     elseif (associated(vtab_r(nv)%ivar2_p)) then
        allocate(iscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, ivara=iscr2)
        call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            ivara2=iscr2, ivarb2=vtab_r(nv)%ivar2_p)
        deallocate(iscr2)
     elseif (associated(vtab_r(nv)%rvar1_p)) then
        allocate(rscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, rvara=rscr1)
        call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            rvara1=rscr1, rvarb1=vtab_r(nv)%rvar1_p)
        deallocate(rscr1)
     elseif (associated(vtab_r(nv)%rvar2_p)) then

        allocate(rscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, rvara=rscr2)
        
        if (varn /= 'VXE2' .and. varn /= 'VYE2' .and. varn /= 'VZE2') then
           call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                               rvara2=rscr2, rvarb2=vtab_r(nv)%rvar2_p)
        endif

! Horizontal velocity at IV points that are inside new mesh refinement area
! cannot be copied or interpolated.  Instead, vxe, vye, vze are first 
! interpolated to IW points, and the NEW grid horizontal velocity is then
! obtained by projection onto IV faces.  In preparation, velocity is copied to
! storage arrays here.  Similarly, soil_water and soil_energy require special
! processing, so they are also stored here.

        if (varn == 'VC') then
           vc_og = rscr2
        elseif (varn == 'WC') then
           wc_og = rscr2
        elseif (varn == 'VC_ACCUM') then
           vc_accum_og = rscr2
        elseif (varn == 'LAND%SOIL_WATER') then
           soil_water_og = rscr2
        elseif (varn == 'LAND%SOIL_ENERGY') then
           soil_energy_og = rscr2
        elseif (varn == 'VXE2') then
           vxe2_og = rscr2
        elseif (varn == 'VYE2') then
           vye2_og = rscr2
        elseif (varn == 'VZE2') then
           vze2_og = rscr2
        endif

        deallocate(rscr2)

     elseif (associated(vtab_r(nv)%dvar1_p)) then
        allocate(dscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, dvara=dscr1)
        call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            dvara1=dscr1, dvarb1=vtab_r(nv)%dvar1_p)
        deallocate(dscr1)
     elseif (associated(vtab_r(nv)%dvar2_p)) then
        allocate(dscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, dvara=dscr2)
        call interp_addgrid(ndims, idims, jdims, varn, stagpt, &
                            dvara2=dscr2, dvarb2=vtab_r(nv)%dvar2_p)
        deallocate(dscr2)
     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

! Allocate and diagnose OLD grid vxe_og, vye_og, and vze_og

  allocate (vxe_og(nza_og,nwa_og))
  allocate (vye_og(nza_og,nwa_og))
  allocate (vze_og(nza_og,nwa_og))

  call vel_t3d_hex_oldgrid()

! Interpolate vxe, vye, vze from OLD grid to NEW grid

  ndims = 2
  idims(1) = nza_og
  idims(2) = nwa_og
  jdims(1) = mza
  jdims(2) = mwa

  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vxe_og, rvarb2=vxe)
  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vye_og, rvarb2=vye)
  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vze_og, rvarb2=vze)

! Do MPI communication of vxe, vye, vze between parallel processes

  mrl = 1
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
  endif

! Loop over prognostic IV points

  do j = 1,jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

! Check ngr values of IW points that are adjacent to this IV point

     if (itab_w(iw1)%ngr > ngrids_og .or. itab_w(iw2)%ngr > ngrids_og) then

! This IV point was not filled from OLD history file because at least one
! adjacent IW point is part of NEW mesh refinement region.  Thus, compute
! VC by projecting vxe, vye, vze onto IV face

        kb = lpv(iv)

! Vertical loop over T levels

        do k = kb,mza

! Update VM tendency from turbulent fluxes

           vc(k,iv) = .5 &
                    * (vnx(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                    +  vny(iv) * (vye(k,iw1) + vye(k,iw2)) &
                    +  vnz(iv) * (vze(k,iw1) + vze(k,iw2)))

           vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

           if (allocated(vp)) vp (k,iv) = vc (k,iv)
           vmp(k,iv) = vmc(k,iv)

        enddo

     endif

     vc (1:kb-1,iv) = 0.
     vmc(1:kb-1,iv) = 0.
     vmp(1:kb-1,iv) = 0.
     if (allocated(vp)) vp(1:kb-1,iv) = 0.

  enddo

  ! MPI parallel send/recv of V group

  if (iparallel == 1) then
     if (allocated(vp)) then
        call mpi_send_v(mrl, rvara1=vmc, rvara2=vc, rvara3=vmp, rvara4=vp)
        call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc, rvara3=vmp, rvara4=vp)
     else
        call mpi_send_v(mrl, rvara1=vmc, rvara2=vc, rvara3=vmp)
        call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc, rvara3=vmp)
     endif
  endif

! Now set underground W momentum/velocity to 0

  do iw = 2, mwa
     wmc(1:lpw(iw)-1,iw) = 0.
     wc (1:lpw(iw)-1,iw) = 0.
  enddo

! Reset VXE2, VYE2, and VZE2 in new mesh area where terrain has changed

  call diagvel_t3d_init_addgrid()

! Diagnose OLD grid vxe_accum_og, vye_accum_og, vze_accum_og

  do iw_og = 2,nwa_og
     kb_og = lpw_og(iw_og)

     do k_og = kb_og, nza_og
        vc_og(k_og,iw_og) = vc_accum_og(k_og,iw_og)
     enddo
     vc_og(1:kb_og-1,iw_og) = vc_og(kb_og,iw_og)
  enddo

  wc_og   = 0.
  vxe2_og = 0.
  vye2_og = 0.
  vze2_og = 0.

  call vel_t3d_hex_oldgrid()

! Interpolate vxe_accum, vye_accum, vze_accum from OLD grid to NEW grid

  ndims = 2
  idims(1) = nza_og
  idims(2) = nwa_og
  jdims(1) = mza
  jdims(2) = mwa

  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vxe_og, rvarb2=vxe)
  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vye_og, rvarb2=vye)
  call interp_addgrid(ndims, idims, jdims, varn, 'AW', rvara2=vze_og, rvarb2=vze)

  deallocate (vxe2_og, vye2_og, vze2_og)
  deallocate (vc_og, wc_og, vc_accum_og, vxe_og, vye_og, vze_og)

! Do MPI communication of vxe_accum, vye_accum, vze_accum between parallel processes

  mrl = 1
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
  endif

! Loop over prognostic IV points

  do j = 1,jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

! Check ngr values of IW points that are adjacent to this IV point

     if (itab_w(iw1)%ngr > ngrids_og .or. itab_w(iw2)%ngr > ngrids_og) then

! This IV point was not filled from OLD history file because at least one
! adjacent IW point is part of NEW mesh refinement region.  Thus, compute
! VC_ACCUM by projecting vxe_accum, vye_accum, vze_accum onto IV face.

        kb = lpv(iv)

! Vertical loop over T levels

        do k = kb,mza

! Update VM tendency from turbulent fluxes

           vc_accum(k,iv) = .5 &
                          * (vnx(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                          +  vny(iv) * (vye(k,iw1) + vye(k,iw2)) &
                          +  vnz(iv) * (vze(k,iw1) + vze(k,iw2)))

        enddo

     endif

  enddo


! Loop over all land cells on NEW grid, and on each, find index of identical
! or nearest land cell on OLD grid.

  do iwl = 2,mwl
     iwl_og = itab_wladd(iwl)%iwl_og

! Transfer soil_water and soil_energy from OLD grid to NEW grid under the
! constraint that soil water HEAD (gravity contribution excluded), TEMPERATURE,
! and FRACLIQ are held constant.

     do k = 1,nzg

        nts    = land%ntext_soil(k,iwl)
        nts_og = ntext_soil_og(k,iwl_og)

        ! Get water fraction and head [Assume that old grid value of "flag_vg_og",
        ! which has not been implemented yet, is the same as the new grid value,
        ! flag_vg, for each transfer of soil water from old grid to new grid.
        ! Thus, flag_vg is passed to subroutine soilwatpot in place of flag_vg_og.]

        call soil_wat2pot(iwl_og,nts_og,land%flag_vg(iwl),soil_water_og(k,iwl_og), slzt(k), &
           water_frac_ul, water_frac, psi, head, headp)

        ! Diagnose soil water from given head0, nts, flag_vg, and slzt

        call soil_pot2wat(iwl, nts, land%flag_vg(iwl), head, slzt(k), &
                          water_frac_ul, land%soil_water(k,iwl))

        ! Diagnose temperature and liquid fraction of OLD grid soil water

        call qwtk(soil_energy_og(k,iwl_og),soil_water_og(k,iwl_og)*1.e3,  &
                  slcpd(nts_og),soil_tempk,soil_fracliq)

        ! Compute NEW grid soil energy [J/m^3] given NEW grid soil textural
        ! class and OLD grid temperature, water content, and liquid fraction
        ! (as opposed to ice fraction)

        soil_tempc = soil_tempk - 273.15

        if (soil_tempc > 0.) then

           land%soil_energy(k,iwl)                             &
              = soil_tempc   * slcpd(nts)                        &
              + soil_tempc   * land%soil_water(k,iwl) * cliq1000 &
              + soil_fracliq * land%soil_water(k,iwl) * alli1000
             
        else
      
           land%soil_energy(k,iwl)                             &
              = soil_tempc   * slcpd(nts)                        &
              + soil_tempc   * land%soil_water(k,iwl) * cice1000 &
              + soil_fracliq * land%soil_water(k,iwl) * alli1000
             
        endif

     enddo

! Transfer veg_water from OLD grid to NEW grid, but reduce it proportionally if
! the total area index is lower in the NEW grid than in the OLD grid.

     tai = tai_max(land%leaf_class(iwl))
     tai_og = max(.01,tai_max(leaf_class_og(iwl_og)))

     if (tai < tai_og) then
        land%veg_water(iwl) = land%veg_water(iwl) * tai / tai_og
     endif

! Diagnose land%snowfac and land%rough

     land%snowfac(iwl) = 0.
     do k = 1,land%nlev_sfcwater(iwl)
        land%snowfac(iwl) = land%snowfac(iwl) + land%sfcwater_depth(k,iwl)
     enddo
     land%snowfac(iwl) = land%snowfac(iwl) / max(.001,land%veg_height(iwl))
     if (land%snowfac(iwl) > 0.9) land%snowfac(iwl) = 1.0

     land%rough(iwl) = max(soil_rough,land%veg_rough(iwl)) &
                     * (1. - land%snowfac(iwl)) + snow_rough

  enddo

end subroutine hist_read_addgrid
