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
  use mem_para,   only: mgroupsize
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

           write(io6,*) 'calling hist_read_addgrid'
           call hist_read_addgrid()

        else

! If NOT adding new grids (new mesh refinements), read model fields from
! standard history file

           write(io6,*) 'calling hist_read'
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
  use mem_sfcg,    only: itab_wsfc, nwsfc, mwsfc
  use mem_land,    only: itab_land, nland, mland, nzg
  use mem_lake,    only: itab_lake, nlake, mlake
  use mem_sea,     only: itab_sea, nsea, msea
  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  implicit none

  integer       :: nv, nvcnt, ns, ndims, idims(3)
  character(32) :: varn
  character (2) :: stagpt
  integer       :: ilocal(max(mwa,mva,mma,mwsfc,mland,mlake,msea,mwnud))

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
        ilocal(1:mwa) = itab_w(1:mwa)%iwglobe
        idims(ndims) = mwa
     elseif (stagpt == 'AV' .and. idims(ndims) == nva) then
        ilocal(1:mva) = itab_v(1:mva)%ivglobe
        idims(ndims) = mva
     elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
        ilocal(1:mma) = itab_m(1:mma)%imglobe
        idims(ndims) = mma
     elseif (stagpt == 'CW' .and. idims(ndims) == nwsfc) then
        ilocal(1:mwsfc) = itab_wsfc(1:mwsfc)%iwglobe
        idims(ndims) = mwsfc
     elseif (stagpt == 'LW' .and. idims(ndims) == nland) then
        ilocal(1:mland) = itab_land(1:mland)%iwglobe
        idims(ndims) = mland
     elseif (stagpt == 'RW' .and. idims(ndims) == nlake) then
        ilocal(1:mlake) = itab_lake(1:mlake)%iwglobe
        idims(ndims) = mlake
     elseif (stagpt == 'SW' .and. idims(ndims) == nsea) then
        ilocal(1:msea) = itab_sea(1:msea)%iwglobe
        idims(ndims) = msea
     elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
        ilocal(1:mwnud) = itab_wnud(1:mwnud)%iwnudglobe
        idims(ndims) = mwnud
     else

        ! TODO: Const values
        ! TODO: Land U and M; Sea U and M? (probably not!)

        stop "invalid array size in history_read"
     endif

     ns = idims(ndims)

     if     (associated(vtab_r(nv)%ivar1_p)) then
        call shdf5_irec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%ivar2_p)) then
        call shdf5_irec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%ivar3_p)) then
        call shdf5_irec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, &
                        points=ilocal(1:ns), stagpt=stagpt)

     elseif (associated(vtab_r(nv)%rvar1_p)) then
        call shdf5_irec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%rvar2_p)) then
        call shdf5_irec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%rvar3_p)) then
        call shdf5_irec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, &
                        points=ilocal(1:ns), stagpt=stagpt)

     elseif (associated(vtab_r(nv)%dvar1_p)) then
        call shdf5_irec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%dvar2_p)) then
        call shdf5_irec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     elseif (associated(vtab_r(nv)%dvar3_p)) then
        call shdf5_irec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, &
                        points=ilocal(1:ns), stagpt=stagpt)
     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

end subroutine hist_read

!=========================================================================

subroutine hist_read_addgrid()

  ! This subroutine needs to be revised for OLAM-SOIL version of model.  One
  ! possible decision that can be made in the revision is that only the atmospheric
  ! grid will be refined, while the SURFACE GRID will remain unchanged.  This would 
  ! avoid the (questionable) interpolation/remapping of the soil grid.  The interface
  ! between atmosphere and surface grid cells would, of course, need to be recomputed.

  ! For now, surface-related statements that do not compile have been commented out.

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza, lpv, &
                         vnxo2, vnyo2, vnzo2, lpw
  use mem_ijtabs,  only: itab_w, itab_v, itab_m, jtab_v, jtv_prog
  use misc_coms,   only: io6, runtype, iparallel
  use mem_basic,   only: vc, vmc, wc, wmc, rho, vxe, vye, vze
  use var_tables,  only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils,  only: shdf5_info, shdf5_irec
  use mem_sfcg,    only: sfcg, mwsfc
  use mem_land,    only: land, mland, omland, nzg
  use leaf_coms,   only: tai_max, soil_rough, snow_rough
  use leaf4_soil,  only: soil_wat2pot, soil_pot2wat
  use mem_sea,     only: nsea, msea
  use consts_coms, only: r8, cliq1000, cice1000, alli1000
  use mem_nudge,   only: nwnud, mwnud, itab_wnud
  use therm_lib,   only: qwtk

  use mem_addgrid, only: interp_addgrid, &
                         nza_og, nwa_og, nva_og, nma_og, nland_og, nsea_og, &
                         nzg_og, ngrids_og, lpw_og, itab_landadd, &
                         soil_water_og, soil_energy_og, specifheat_drysoil_og, &
                         leaf_class_og

  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

  implicit none

  integer :: j, iw, iw1, iw2, iv, nv, kb, k, mrl
  integer :: nvcnt, ndims, idims(3)
  character(32) :: varn
  character (2) :: stagpt

  integer :: iland, iland_og, kb_og, iw_og, k_og, iwsfc, iwsfc_og
  real :: soil_watfrac_ul, soil_watfrac, psi, psip
  real :: soil_tempk, soil_tempc, soil_fracliq, tai, tai_og

! Scratch arrays for copying input

  integer,  allocatable :: iscr1(:)
  integer,  allocatable :: iscr2(:,:)
  integer,  allocatable :: iscr3(:,:,:)

  real,     allocatable :: rscr1(:)
  real,     allocatable :: rscr2(:,:)
  real,     allocatable :: rscr3(:,:,:)

  real(r8), allocatable :: dscr1(:)
  real(r8), allocatable :: dscr2(:,:)
  real(r8), allocatable :: dscr3(:,:,:)

  allocate (soil_water_og (nzg_og,nland_og))
  allocate (soil_energy_og(nzg_og,nland_og))

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
         trim(varn) == 'VKM'                 .or. &
         trim(varn) == 'VKH'                 .or. &
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
         trim(varn) == 'SFCG%RHOS'           .or. &
         trim(varn) == 'SFCG%VELS'           .or. &
         trim(varn) == 'SFCG%PRSS'           .or. &
         trim(varn) == 'SFCG%AIRTEMP'        .or. &
         trim(varn) == 'SFCG%AIRSHV'         .or. &
         trim(varn) == 'SFCG%USTAR'          .or. &
         trim(varn) == 'SFCG%VKMSFC'         .or. &
         trim(varn) == 'SFCG%SFLUXT'         .or. &
         trim(varn) == 'SFCG%SFLUXR'         .or. &
         trim(varn) == 'SFCG%SFLUXC'         .or. &
         trim(varn) == 'SFCG%SXFER_T'        .or. &
         trim(varn) == 'SFCG%SXFER_R'        .or. &
         trim(varn) == 'SFCG%SXFER_C'        .or. &
         trim(varn) == 'SFCG%SXFER_T'        .or. &
         trim(varn) == 'SFCG%SXFER_R'        .or. &
         trim(varn) == 'SFCG%SXFER_C'        .or. &
         trim(varn) == 'LAND%ED_GGAER'       .or. &
         trim(varn) == 'LAND%ED_ZETA'        .or. &
         trim(varn) == 'LAND%ED_RIB'         .or. &
         trim(varn) == 'SFCG%ALBEDO_BEAM'    .or. &
         trim(varn) == 'SFCG%ALBEDO_DIFFUSE' .or. &
         trim(varn) == 'SFCG%RSHORT'         .or. &
         trim(varn) == 'SFCG%RSHORT_DIFFUSE' .or. &
         trim(varn) == 'SFCG%RSHORT_CLR'     .or. &
         trim(varn) == 'SFCG%RLONG'          .or. &
         trim(varn) == 'SFCG%RLONG_ALBEDO'   .or. &
         trim(varn) == 'SFCG%RLONGUP'        .or. &
         trim(varn) == 'LAND%RSHORT_G'       .or. &
         trim(varn) == 'LAND%RSHORT_S'       .or. &
         trim(varn) == 'LAND%RSHORT_V'       .or. &
         trim(varn) == 'LAND%RLONG_G'        .or. &
         trim(varn) == 'LAND%RLONG_S'        .or. &
         trim(varn) == 'LAND%RLONG_V'        .or. &
         trim(varn) == 'LAND%COSZ'           .or. &
         trim(varn) == 'SFCG%PCPG'           .or. &
         trim(varn) == 'SFCG%QPCPG'          .or. &
         trim(varn) == 'SFCG%DPCPG'          .or. &
         trim(varn) == 'SFCG%CAN_DEPTH'      .or. &
         trim(varn) == 'LAND%HEAD0'          .or. &
         trim(varn) == 'LAND%HCAPVEG'        .or. &
         trim(varn) == 'LAND%SURFACE_SSH'    .or. &
         trim(varn) == 'LAND%GROUND_SHV'     .or. &
         trim(varn) == 'SFCG%ROUGH'          .or. &
         trim(varn) == 'LAND%VEG_FRACAREA'   .or. &
         trim(varn) == 'LAND%VEG_LAI'        .or. &
         trim(varn) == 'LAND%VEG_ROUGH'      .or. &
         trim(varn) == 'LAND%VEG_HEIGHT'     .or. &
         trim(varn) == 'LAND%VEG_ALBEDO'     .or. &
         trim(varn) == 'LAND%VEG_TAI'        .or. &
         trim(varn) == 'LAND%VEG_NDVIC'      .or. &
         trim(varn) == 'LAND%SNOWFAC'        .or. &
         trim(varn) == 'LAND%VF'             .or. &
         trim(varn) == 'SEA%SEA_USTAR'       .or. &
         trim(varn) == 'SEA%ICE_USTAR'       .or. &
         trim(varn) == 'SEA%SEA_VKMSFC'      .or. &
         trim(varn) == 'SEA%ICE_VKMSFC'      .or. &
         trim(varn) == 'SEA%SEA_SFLUXT'      .or. &
         trim(varn) == 'SEA%ICE_SFLUXT'      .or. &
         trim(varn) == 'SEA%SEA_SFLUXR'      .or. &
         trim(varn) == 'SEA%ICE_SFLUXR'      .or. &
         trim(varn) == 'SEA%SEA_SXFER_T'     .or. &
         trim(varn) == 'SEA%ICE_SXFER_T'     .or. &
         trim(varn) == 'SEA%SEA_SXFER_R'     .or. &
         trim(varn) == 'SEA%ICE_SXFER_R'     .or. &
         trim(varn) == 'SEA%SEA_ALBEDO'      .or. &
         trim(varn) == 'SEA%ICE_ALBEDO'      .or. &
         trim(varn) == 'SEA%SEA_RLONGUP'     .or. &
         trim(varn) == 'SEA%ICE_RLONGUP'     .or. &
         trim(varn) == 'SEA%ICE_NET_RLONG'   .or. &
         trim(varn) == 'SEA%ICE_NET_RSHORT'  .or. &
         trim(varn) == 'SEA%SEATC'           .or. &
         trim(varn) == 'SEA%SURFACE_SSH'     .or. &
         trim(varn) == 'SEA%SEA_SFC_SSH'     .or. &
         trim(varn) == 'SEA%ICE_SFC_SSH')    cycle

     ! We want it...read it if it's in the history file

     call shdf5_info(varn, ndims, idims)

     write(io6,*) 'varn ',varn

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        write(io6,*)
        write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
        write(io6,*)
        cycle
     endif

     if     (associated(vtab_r(nv)%ivar1_p)) then

        allocate(iscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, ivar1=iscr1)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            ivara1=iscr1, ivarb1=vtab_r(nv)%ivar1_p)
        deallocate(iscr1)

     elseif (associated(vtab_r(nv)%ivar2_p)) then

        allocate(iscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, ivar2=iscr2)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            ivara2=iscr2, ivarb2=vtab_r(nv)%ivar2_p)
        deallocate(iscr2)

     elseif (associated(vtab_r(nv)%ivar3_p)) then

        allocate(iscr3(idims(1),idims(2),idims(3)))
        call shdf5_irec(ndims, idims, varn, ivar3=iscr3)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            ivara3=iscr3, ivarb3=vtab_r(nv)%ivar3_p)
        deallocate(iscr3)

     elseif (associated(vtab_r(nv)%rvar1_p)) then

        allocate(rscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, rvar1=rscr1)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            rvara1=rscr1, rvarb1=vtab_r(nv)%rvar1_p)
        deallocate(rscr1)

     elseif (associated(vtab_r(nv)%rvar2_p)) then

        allocate(rscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, rvar2=rscr2)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            rvara2=rscr2, rvarb2=vtab_r(nv)%rvar2_p)

! Horizontal velocity at IV points that are inside new mesh refinement area
! cannot be copied or interpolated.  Instead, vxe, vye, vze are first 
! interpolated to IW points, and the NEW grid horizontal velocity is then
! obtained by projection onto IV faces.  In preparation, velocity is copied to
! storage arrays here.  Similarly, soil_water and soil_energy require special
! processing, so they are also stored here.

        if     (varn == 'LAND%SOIL_WATER') then
           soil_water_og = rscr2
        elseif (varn == 'LAND%SOIL_ENERGY') then
           soil_energy_og = rscr2
        endif

        deallocate(rscr2)

     elseif (associated(vtab_r(nv)%rvar3_p)) then

        allocate(rscr3(idims(1),idims(2),idims(3)))
        call shdf5_irec(ndims, idims, varn, rvar3=rscr3)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            rvara3=rscr3, rvarb3=vtab_r(nv)%rvar3_p)
        deallocate(rscr3)

     elseif (associated(vtab_r(nv)%dvar1_p)) then

        allocate(dscr1(idims(1)))
        call shdf5_irec(ndims, idims, varn, dvar1=dscr1)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            dvara1=dscr1, dvarb1=vtab_r(nv)%dvar1_p)
        deallocate(dscr1)

     elseif (associated(vtab_r(nv)%dvar2_p)) then

        allocate(dscr2(idims(1),idims(2)))
        call shdf5_irec(ndims, idims, varn, dvar2=dscr2)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            dvara2=dscr2, dvarb2=vtab_r(nv)%dvar2_p)
        deallocate(dscr2)

     elseif (associated(vtab_r(nv)%dvar3_p)) then

        allocate(dscr3(idims(1),idims(2),idims(3)))
        call shdf5_irec(ndims, idims, varn, dvar3=dscr3)
        call interp_addgrid(ndims, idims, varn, stagpt, &
                            dvara3=dscr3, dvarb3=vtab_r(nv)%dvar3_p)
        deallocate(dscr3)

     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

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

           vc(k,iv) = vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                    + vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                    + vnzo2(iv) * (vze(k,iw1) + vze(k,iw2))


           vmc(k,iv) = vc(k,iv) * .5 * real(rho(k,iw1) + rho(k,iw2))

        enddo

     endif

     vc (1:kb-1,iv) = 0.
     vmc(1:kb-1,iv) = 0.

  enddo

  ! MPI parallel send/recv of V group

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif

! Now set underground W momentum/velocity to 0

  do iw = 2, mwa
     wmc(1:lpw(iw)-1,iw) = 0.
     wc (1:lpw(iw)-1,iw) = 0.
  enddo

!!! Loop over prognostic IV points
!!
!!  do j = 1,jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)
!!     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!!
!!! Check ngr values of IW points that are adjacent to this IV point
!!
!!     if (itab_w(iw1)%ngr > ngrids_og .or. itab_w(iw2)%ngr > ngrids_og) then
!!
!!! This IV point was not filled from OLD history file because at least one
!!! adjacent IW point is part of NEW mesh refinement region.  Thus, compute
!!! VC_ACCUM by projecting vxe_accum, vye_accum, vze_accum onto IV face.
!!
!!        kb = lpv(iv)
!!
!!! Vertical loop over T levels
!!
!!        do k = kb,mza
!!
!!! Update VM tendency from turbulent fluxes
!!
!!           vc_accum(k,iv) = vnx(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
!!                          + vny(iv) * (vye(k,iw1) + vye(k,iw2)) &
!!                          + vnz(iv) * (vze(k,iw1) + vze(k,iw2))
!!
!!        enddo
!!
!!     endif
!!
!!  enddo

! Loop over all land cells on NEW grid, and on each, find index of identical
! or nearest land cell on OLD grid.

  do iland = 2,mland
     iwsfc = iland + omland
     iland_og = itab_landadd(iland)%iland_og
     iwsfc_og = iland_og

! Transfer soil_water and soil_energy from OLD grid to NEW grid under the
! constraint that soil water potential, TEMPERATURE, and FRACLIQ are held constant.

     do k = 1,nzg

  !! In the following (if we decide to retain the hist_read_addgrid capability),
  !! we will need to consider the new ability to horizontally average soil
  !! physical properties.
 
      ! Get water fraction and potential

  !!      call soil_wat2pot(iland_og,nts_og,soil_water_og(k,iland_og), &
  !!         soil_watfrac_ul, soil_watfrac, psi, psip)

        ! Diagnose soil water from given psi and nts

  !!      call soil_pot2wat(iland, nts, psi, &
  !!                        soil_watfrac_ul, land%soil_water(k,iland))

        ! Diagnose temperature and liquid fraction of OLD grid soil water

 !com       call qwtk(soil_energy_og(k,iland_og),soil_water_og(k,iland_og)*1.e3,  &
 !com                 specifheat_drysoil(k,iland_og),soil_tempk,soil_fracliq)

        ! Compute NEW grid soil energy [J/m^3] given NEW grid soil textural
        ! class and OLD grid temperature, water content, and liquid fraction
        ! (as opposed to ice fraction)

        soil_tempc = soil_tempk - 273.15

        if (soil_tempc > 0.) then

           land%soil_energy(k,iland)                               &
              = soil_tempc   * land%specifheat_drysoil(k,iland)    &
              + soil_tempc   * land%soil_water(k,iland) * cliq1000 &
              + soil_fracliq * land%soil_water(k,iland) * alli1000

        else

           land%soil_energy(k,iland)                               &
              = soil_tempc   * land%specifheat_drysoil(k,iland)    &
              + soil_tempc   * land%soil_water(k,iland) * cice1000 &
              + soil_fracliq * land%soil_water(k,iland) * alli1000

        endif

     enddo

! Transfer veg_water from OLD grid to NEW grid, but reduce it proportionally if
! the total area index is lower in the NEW grid than in the OLD grid.

     tai = tai_max(sfcg%leaf_class(iwsfc))
     tai_og = max(.01,tai_max(leaf_class_og(iwsfc_og)))

     if (tai < tai_og) then
        land%veg_water(iland) = land%veg_water(iland) * tai / tai_og
     endif

! Diagnose land%snowfac and land%rough

     land%snowfac(iland) = 0.
     do k = 1,land%nlev_sfcwater(iland)
        land%snowfac(iland) = land%snowfac(iland) + land%sfcwater_depth(k,iland)
     enddo
     land%snowfac(iland) = land%snowfac(iland) / max(.001,land%veg_height(iland))
     if (land%snowfac(iland) > 0.9) land%snowfac(iland) = 1.0

     sfcg%rough(iwsfc) = max(soil_rough,land%veg_rough(iland)) &
                     * (1. - land%snowfac(iland)) + snow_rough

  enddo

end subroutine hist_read_addgrid
