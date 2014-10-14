!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

Module mem_sflux

  use max_dims,  only: maxremote
  implicit none

  private :: maxremote

! Total number of sea and land flux cells

  integer :: nseaflux,  mseaflux
  integer :: nlandflux, mlandflux

! Land and sea flux variables

  integer, parameter :: maxnpolyf = 30

  Type lflux_vars
     logical :: sendf(maxremote) = .false.
     
     integer :: ifglobe =  1
     integer :: iwrank  = -1
     integer :: iw      =  0
     integer :: kw      =  0
     integer :: iwls    =  0
     integer :: npoly   =  0

     real :: xem(maxnpolyf) = 0.0
     real :: yem(maxnpolyf) = 0.0
     real :: zem(maxnpolyf) = 0.0
     real :: zm (maxnpolyf) = 0.0  ! In case it is needed

     real :: dtf     = 0.0
     real :: area    = 0.0
     real :: xef     = 0.0
     real :: yef     = 0.0
     real :: zef     = 0.0
     real :: glatf   = 0.0
     real :: glonf   = 0.0
     real :: arf_atm = 0.0
     real :: arf_sfc = 0.0
     real :: arf_kw  = 0.0
     real :: sfluxt  = 0.0
     real :: sfluxr  = 0.0
     real :: sfluxc  = 0.0
     real :: ustar   = 0.0
     real :: rhos    = 0.0
     real :: prss    = 0.0
     real :: vels    = 0.0
     real :: airtemp = 0.0
     real :: airshv  = 0.0
     real :: sxfer_t = 0.0
     real :: sxfer_r = 0.0
     real :: sxfer_c = 0.0
     real :: ed_zeta = 0.0
     real :: ed_rib  = 0.0
     real :: ggaer   = 0.0
     real :: pcpg    = 0.0
     real :: qpcpg   = 0.0
     real :: dpcpg   = 0.0
     real :: rlong   = 0.0
     real :: rshort  = 0.0
     real :: rshort_diffuse = 0.0
  End Type lflux_vars

  Type sflux_vars
     logical :: sendf(maxremote) = .false.
     
     integer :: ifglobe =  1
     integer :: iwrank  = -1
     integer :: iw      =  0
     integer :: kw      =  0
     integer :: iwls    =  0
     integer :: npoly   =  0

     real :: xem(maxnpolyf) = 0.0
     real :: yem(maxnpolyf) = 0.0
     real :: zem(maxnpolyf) = 0.0
     real :: zm (maxnpolyf) = 0.0  ! In case it is needed

     real :: dtf     = 0.0
     real :: area    = 0.0
     real :: xef     = 0.0
     real :: yef     = 0.0
     real :: zef     = 0.0
     real :: glatf   = 0.0
     real :: glonf   = 0.0
     real :: arf_atm = 0.0
     real :: arf_sfc = 0.0
     real :: arf_kw  = 0.0
     real :: rhos    = 0.0
     real :: airtemp = 0.0
     real :: airshv  = 0.0
     real :: pcpg    = 0.0
     real :: qpcpg   = 0.0
     real :: dpcpg   = 0.0
     real :: rlong   = 0.0
     real :: rshort  = 0.0
     real :: rshort_diffuse = 0.0

     real ::     sfluxt = 0.0
     real :: sea_sfluxt = 0.0
     real :: ice_sfluxt = 0.0

     real ::     sfluxr = 0.0
     real :: sea_sfluxr = 0.0
     real :: ice_sfluxr = 0.0

     real ::     sfluxc = 0.0

     real ::     ustar  = 0.0
     real :: sea_ustar  = 0.0
     real :: ice_ustar  = 0.0

     real ::     sxfer_t = 0.0
     real :: sea_sxfer_t = 0.0
     real :: ice_sxfer_t = 0.0

     real ::     sxfer_r = 0.0
     real :: sea_sxfer_r = 0.0
     real :: ice_sxfer_r = 0.0

     real ::     sxfer_c = 0.0
     real ::     ed_zeta = 0.0
     real ::      ed_rib = 0.0

     real ::      ggaer  = 0.0
     real ::  sea_ggaer  = 0.0
     real ::  ice_ggaer  = 0.0
  End Type sflux_vars

  type (sflux_vars), allocatable, target :: seaflux(:)
  type (lflux_vars), allocatable, target :: landflux(:)

  type flux_pd_vars
     integer :: iw      =  0
     integer :: iwls    =  0
  end type flux_pd_vars

  type (flux_pd_vars), allocatable, target :: seaflux_pd(:)
  type (flux_pd_vars), allocatable, target :: landflux_pd(:)

!----------------------------------------------------------------------------

  Type seafluxg_vars
     integer :: isf_myrank = -1
     integer :: iwrank = -1
  End Type seafluxg_vars

  Type landfluxg_vars
     integer :: ilf_myrank = -1
     integer :: iwrank = -1
  End Type landfluxg_vars

  type (seafluxg_vars),  allocatable, target :: seafluxg(:)
  type (landfluxg_vars), allocatable, target :: landfluxg(:)

!----------------------------------------------------------------------------

  Type jseaflux_vars
     integer, allocatable :: iseaflux(:)
     integer, allocatable :: jend(:)
  End Type jseaflux_vars

  Type jlandflux_vars
     integer, allocatable :: ilandflux(:)
     integer, allocatable :: jend(:)
  End Type jlandflux_vars

  type (jseaflux_vars)  :: jseaflux(12)  ! 2 + 10 remote nodes
  type (jlandflux_vars) :: jlandflux(12) ! 2 + 10 remote nodes

Contains

!===============================================================================

  subroutine filltab_sflux()

    use var_tables, only: vtab_r, num_var, increment_vtable
    use misc_coms,  only: iparallel, runtype, ipar_out
    implicit none

    if (allocated(seaflux)) then

       if ( (runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE') .or. &
            (iparallel == 1 .and. ipar_out == 0) ) then

          call increment_vtable('SEAFLUX%IFGLOBE', 'SF', noread=.true.)
          vtab_r(num_var)%ivar1_p => seaflux%ifglobe

          call increment_vtable('SEAFLUX%IWRANK',  'SF', noread=.true.)
          vtab_r(num_var)%ivar1_p => seaflux%iwrank

       endif

       call increment_vtable('SEAFLUX%SFLUXT', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sfluxt

       call increment_vtable('SEAFLUX%SEA_SFLUXT', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sea_sfluxt

       call increment_vtable('SEAFLUX%ICE_SFLUXT', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%ice_sfluxt

       call increment_vtable('SEAFLUX%SFLUXR', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sfluxr

       call increment_vtable('SEAFLUX%SEA_SFLUXR', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sea_sfluxr

       call increment_vtable('SEAFLUX%ICE_SFLUXR', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%ice_sfluxr

       call increment_vtable('SEAFLUX%SFLUXC', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sfluxc

    endif

    if (allocated(landflux)) then

       if ( (runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE') .or. &
            (iparallel == 1 .and. ipar_out == 0) ) then

          call increment_vtable('LANDFLUX%IFGLOBE', 'LF', noread=.true.)
          vtab_r(num_var)%ivar1_p => landflux%ifglobe

          call increment_vtable('LANDFLUX%IWRANK',  'LF', noread=.true.)
          vtab_r(num_var)%ivar1_p => landflux%iwrank

       endif

       call increment_vtable('LANDFLUX%SFLUXT', 'LF')
       vtab_r(num_var)%rvar1_p => landflux%sfluxt

       call increment_vtable('LANDFLUX%SFLUXR', 'LF')
       vtab_r(num_var)%rvar1_p => landflux%sfluxr

       call increment_vtable('LANDFLUX%SFLUXC', 'LF')
       vtab_r(num_var)%rvar1_p => landflux%sfluxc

    endif

  end subroutine filltab_sflux

!============================================================================

   subroutine fill_jflux()

   use mem_ijtabs, only: itab_w, mrls, itabg_w
   use misc_coms,  only: io6, dtlm, iparallel, isubdomain
   use leaf_coms,  only: dt_leaf
   use mem_leaf,   only: itab_wl, itabg_wl
   use mem_sea,    only: itab_ws, itabg_ws
   use mem_para,   only: myrank, nsends_wsf, nsends_wlf

   implicit none

   integer :: mrl,mrlw,isf,ilf,iw,iws,iwl,iloop,jsend,jend

! Allocate and initialize JSEAFLUX%JEND

   do iloop = 1,12
      allocate (jseaflux(iloop)%jend(mrls)) ; jseaflux(iloop)%jend = 0
   enddo

! Allocate ISEAFLUX member for first 2 JSEAFLUX data structures

   allocate (jseaflux(1)%iseaflux(mseaflux)) ; jseaflux(1)%iseaflux = 0
   allocate (jseaflux(2)%iseaflux(mseaflux)) ; jseaflux(2)%iseaflux = 0

   do mrl = mrls,1,-1
      do isf = 2,mseaflux
         iw  = seaflux(isf)%iw   ! global index
         iws = seaflux(isf)%iwls ! global index
      
! If run is parallel, convert iw to local domain

         if (isubdomain == 1) then
            iw  = itabg_w (iw )%iw_myrank
            iws = itabg_ws(iws)%iws_myrank
         endif

! In a parallel run on a given rank, ISF ranges through all active flux cells,
! whether the associated atmospheric IW cell, or the sea IWS cell, or both are
! prognosed on the rank.  Even if the IWS cell is not prognosed on that rank,
! it will still exist in memory and will have an index > 1, but the IW cell may
! not exist, in which case IW will have a value of 1.  The following IF 
! statement checks for this contingency.

         if (iw > 1) then
            mrlw = itab_w(iw)%mrlw
         else
            mrlw = 1
         endif

         seaflux(isf)%dtf = min( real(dtlm(mrlw)), dt_leaf )

         if (mrlw == mrl) then
      
            if (iw > 1) then
               if (isubdomain == 0 .or.  &
                  (isubdomain == 1 .and. itab_w(iw)%irank == myrank)) then

                  jseaflux(1)%jend(1:mrl) = jseaflux(1)%jend(1:mrl) + 1
                  jseaflux(1)%iseaflux(jseaflux(1)%jend(1)) = isf
               endif
            endif

            if (isubdomain == 0 .or.  &
               (isubdomain == 1 .and. itab_ws(iws)%irank == myrank)) then

               jseaflux(2)%jend(1:mrl) = jseaflux(2)%jend(1:mrl) + 1
               jseaflux(2)%iseaflux(jseaflux(2)%jend(1)) = isf
            endif

         endif
      enddo
   enddo

! Allocate and initialize JLANDFLUX%JEND

   do iloop = 1,12
      allocate (jlandflux(iloop)%jend(mrls)) ; jlandflux(iloop)%jend = 0
   enddo

! Allocate ILANDFLUX member for first 2 JLANDFLUX structures

   allocate (jlandflux(1)%ilandflux(mlandflux)) ; jlandflux(1)%ilandflux = 0
   allocate (jlandflux(2)%ilandflux(mlandflux)) ; jlandflux(2)%ilandflux = 0

   do mrl = mrls,1,-1
      do ilf = 2,mlandflux
         iw  = landflux(ilf)%iw
         iwl = landflux(ilf)%iwls

! If run is parallel, convert iw to local domain

         if (isubdomain == 1) then
            iw  = itabg_w (iw )%iw_myrank
            iwl = itabg_wl(iwl)%iwl_myrank
         endif

! In a parallel run on a given rank, ILF ranges through all active flux cells,
! whether the associated atmospheric IW cell, or the land IWL cell, or both are
! prognosed on the rank.  Even if the IWL cell is not prognosed on that rank,
! it will still exist in memory and will have an index > 1, but the IW cell may
! not exist, in which case IW will have a value of 1.  The following IF 
! statement checks for this contingency.

         if (iw > 1) then
            mrlw = itab_w(iw)%mrlw
         else
            mrlw = 1
         endif

         landflux(ilf)%dtf = min( real(dtlm(mrlw)), dt_leaf )

         if (mrlw == mrl) then

            if (iw > 1) then
               if (isubdomain == 0 .or.  &
                  (isubdomain == 1 .and. itab_w(iw)%irank == myrank)) then

                  jlandflux(1)%jend(1:mrl) = jlandflux(1)%jend(1:mrl) + 1
                  jlandflux(1)%ilandflux(jlandflux(1)%jend(1)) = ilf
               endif
            endif

            if (isubdomain == 0 .or.  &
               (isubdomain == 1 .and. itab_wl(iwl)%irank == myrank)) then

               jlandflux(2)%jend(1:mrl) = jlandflux(2)%jend(1:mrl) + 1
               jlandflux(2)%ilandflux(jlandflux(2)%jend(1)) = ilf
            endif

         endif
      enddo
   enddo

! Return if run is not parallel

   if (iparallel == 0) return

! Compute and store JSEAFLUX%JEND(1)

   do jsend = 1,nsends_wsf(1)
      jseaflux(2+jsend)%jend(1) = 0
      do isf = 2,mseaflux
         if (seaflux(isf)%sendf(jsend)) then
            jseaflux(2+jsend)%jend(1) = jseaflux(2+jsend)%jend(1) + 1
         endif
      enddo
      jseaflux(2+jsend)%jend(1) = max(1,jseaflux(2+jsend)%jend(1))
   enddo

! Allocate and zero-fill JSEAFLUX%ISEAFLUX

   do jsend = 1,nsends_wsf(1)
      jend = jseaflux(2+jsend)%jend(1)
      allocate (jseaflux(2+jsend)%iseaflux(jend))
                jseaflux(2+jsend)%iseaflux(1:jend) = 0
   enddo

! Initialize JSEAFLUX%JEND counters to zero

   do jsend = 1,nsends_wsf(1)
      jseaflux(2+jsend)%jend(1:mrls) = 0
   enddo

! Compute JSEAFLUX%ISEAFLUX

   do mrl = mrls,1,-1
      do isf = 2,mseaflux

         iw = seaflux(isf)%iw
         if (isubdomain == 1) iw = itabg_w(iw)%iw_myrank

         ! skip this seaflux point if the atm cell does not exist on this node
         if (iw < 2) cycle

         do jsend = 1,nsends_wsf(1)

            if (seaflux(isf)%sendf(jsend) .and. itab_w(iw)%mrlw == mrl) then
               jseaflux(2+jsend)%jend(1:mrl) = jseaflux(2+jsend)%jend(1:mrl) + 1
               jseaflux(2+jsend)%iseaflux(jseaflux(2+jsend)%jend(1)) = isf
            endif

         enddo
      enddo
   enddo

! Compute and store JLANDFLUX%JEND(1)

   do jsend = 1,nsends_wlf(1)
      jlandflux(2+jsend)%jend(1) = 0
      do ilf = 2,mlandflux
         if (landflux(ilf)%sendf(jsend)) then
            jlandflux(2+jsend)%jend(1) = jlandflux(2+jsend)%jend(1) + 1
         endif
      enddo
      jlandflux(2+jsend)%jend(1) = max(1,jlandflux(2+jsend)%jend(1))
   enddo

! Allocate and zero-fill JLANDFLUX%ILANDFLUX

   do jsend = 1,nsends_wlf(1)
      jend = jlandflux(2+jsend)%jend(1)
      allocate (jlandflux(2+jsend)%ilandflux(jend))
                jlandflux(2+jsend)%ilandflux(1:jend) = 0
   enddo

! Initialize JLANDFLUX%JEND counters to zero

   do jsend = 1,nsends_wlf(1)
      jlandflux(2+jsend)%jend(1:mrls) = 0
   enddo

! Compute JLANDFLUX%ILANDFLUX

   do mrl = mrls,1,-1
      do ilf = 2,mlandflux

         iw = landflux(ilf)%iw
         if (isubdomain == 1) iw = itabg_w(iw)%iw_myrank

         ! skip this landflux point if the atm cell does not exist on this node
         if (iw < 2) cycle

         do jsend = 1,nsends_wlf(1)

            if (landflux(ilf)%sendf(jsend) .and. itab_w(iw)%mrlw == mrl) then
               jlandflux(2+jsend)%jend(1:mrl) = jlandflux(2+jsend)%jend(1:mrl) + 1
               jlandflux(2+jsend)%ilandflux(jlandflux(2+jsend)%jend(1)) = ilf
            endif

         enddo
      enddo
   enddo

   return
   end subroutine fill_jflux

!============================================================================

subroutine init_fluxcells()

use mem_grid,   only: nza, nwa, arw0, xem, yem, zem, xew, yew, zew, &
                      zm, topm, topw, arw, volt, nma, dzt, &
                      glatm, glonm, glatw, glonw
use mem_leaf,   only: land, itab_wl
use mem_sea,    only: sea,  itab_ws
use mem_ijtabs, only: itab_w, mrls
use leaf_coms,  only: nwl, ilandgrid
use sea_coms,   only: nws, iseagrid
use misc_coms,  only: io6, rinit
use consts_coms,only: piu180, r8

implicit none

! local variables

integer, parameter :: npmax = 7  ! Heptagons are max polygon for SINGLE GRID LEVEL
                                 ! in atm polygon cell 

integer, parameter :: nqmax = 5  ! Land cells could be up to pentagons because 
                                 ! they are generated by the piecewise-planar
                                 ! topographic surface on the OLAM triangular
                                 ! Delaunay grid and its intersection with
                                 ! model levels

integer, parameter :: nsub = 30

integer :: iws, iwl
integer :: isf, ilf
integer :: iw, im
integer :: j, kpoly, j1, j2, j3

integer :: k, kw0
integer :: mrl
integer :: np, nps, nq
integer :: incr_flux, nflux_est, ifsize
integer :: jtrap, jt, it, ktrap
integer :: iml, jm, jwl, jml
integer :: npoly
integer :: nspoly, ims, jms, nlpoly, jws
integer :: jsub, isub, jc, jpoly

real :: xmean, ymean
real :: xp0, yp0, xq0, yq0
real :: xtrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid x coordinates
real :: ytrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid y coordinates
real :: ztrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid z coordinates
real :: traparea(npmax+nqmax+npmax*nqmax) ! trapezoid area
real :: area
real :: fracx, fracy, xcent, ycent, x1, x2, dx, dy
real :: zmin, zmax, topp, fraczmin, fraczmax, fractopp

real :: xm(0:maxnpolyf), ym(0:maxnpolyf), crossprod, dist, dist12, dist23

real(r8) :: az,bz,cz,alpha, xw, yw
real(r8) :: xp(npmax),yp(npmax)
real(r8) :: xq(nqmax),yq(nqmax),zq(nqmax)

integer :: iter
integer :: istop, i4

real :: raxis

real :: aatmin, aatmax, astmin, astmax, altmin, altmax

! Temporary scratch space during initialization only:

type(lflux_vars), allocatable :: landflux_temp(:)
type(sflux_vars), allocatable :: seaflux_temp (:)

real :: arf_atm_test(nwa)
real :: arf_sea_test(nws)
real :: arf_land_test(nwl)

real :: xeamin(nwa),xelmin(nwl),xesmin(nws)
real :: yeamin(nwa),yelmin(nwl),yesmin(nws)
real :: zeamin(nwa),zelmin(nwl),zesmin(nws)

real :: xeamax(nwa),xelmax(nwl),xesmax(nws)
real :: yeamax(nwa),yelmax(nwl),yesmax(nws)
real :: zeamax(nwa),zelmax(nwl),zesmax(nws)

! Estimate required array sizes and allocate arrays

nflux_est = 50 * nwa
incr_flux = 10 * nwa

allocate (seaflux (nflux_est))
allocate (landflux(nflux_est))

! new - find maximum extents of all atm, land, and sea cells

xeamin(:) = 1.e9
yeamin(:) = 1.e9
zeamin(:) = 1.e9

xelmin(:) = 1.e9
yelmin(:) = 1.e9
zelmin(:) = 1.e9

xesmin(:) = 1.e9
yesmin(:) = 1.e9
zesmin(:) = 1.e9

xeamax(:) = -1.e9
yeamax(:) = -1.e9
zeamax(:) = -1.e9

xelmax(:) = -1.e9
yelmax(:) = -1.e9
zelmax(:) = -1.e9

xesmax(:) = -1.e9
yesmax(:) = -1.e9
zesmax(:) = -1.e9

do iw = 2,nwa
   npoly = itab_w(iw)%npoly
   
   do jm = 1,npoly
      im = itab_w(iw)%im(jm)
      
      if (xeamin(iw) > xem(im)) xeamin(iw) = xem(im)
      if (yeamin(iw) > yem(im)) yeamin(iw) = yem(im)
      if (zeamin(iw) > zem(im)) zeamin(iw) = zem(im)

      if (xeamax(iw) < xem(im)) xeamax(iw) = xem(im)
      if (yeamax(iw) < yem(im)) yeamax(iw) = yem(im)
      if (zeamax(iw) < zem(im)) zeamax(iw) = zem(im)
   enddo
   
enddo

do iwl = 2,nwl
   npoly = itab_wl(iwl)%npoly
   
   do jml = 1,npoly
      iml = itab_wl(iwl)%im(jml)
      
      if (xelmin(iwl) > land%xem(iml)) xelmin(iwl) = land%xem(iml)
      if (yelmin(iwl) > land%yem(iml)) yelmin(iwl) = land%yem(iml)
      if (zelmin(iwl) > land%zem(iml)) zelmin(iwl) = land%zem(iml)

      if (xelmax(iwl) < land%xem(iml)) xelmax(iwl) = land%xem(iml)
      if (yelmax(iwl) < land%yem(iml)) yelmax(iwl) = land%yem(iml)
      if (zelmax(iwl) < land%zem(iml)) zelmax(iwl) = land%zem(iml)
   enddo
   
enddo

do iws = 2,nws
   npoly = itab_ws(iws)%npoly
   
   do jms = 1,npoly
      ims = itab_ws(iws)%im(jms)
      
      if (xesmin(iws) > sea%xem(ims)) xesmin(iws) = sea%xem(ims)
      if (yesmin(iws) > sea%yem(ims)) yesmin(iws) = sea%yem(ims)
      if (zesmin(iws) > sea%zem(ims)) zesmin(iws) = sea%zem(ims)

      if (xesmax(iws) < sea%xem(ims)) xesmax(iws) = sea%xem(ims)
      if (yesmax(iws) < sea%yem(ims)) yesmax(iws) = sea%yem(ims)
      if (zesmax(iws) < sea%zem(ims)) zesmax(iws) = sea%zem(ims)
   enddo
   
enddo

! Initialize sea and land flux cell counters to 1

isf = 1
ilf = 1

landflux(1)%ifglobe = 1
landflux(1)%iw = 1
landflux(1)%iwls = 1

seaflux(1)%ifglobe = 1
seaflux(1)%iw = 1
seaflux(1)%iwls = 1

! Initialize topm and topw to missing values

topm(:) = -1.e6
topw(:) = -1.e6

! Loop over all W points

do iw = 2,nwa

if (mod(iw,1000) == 2) then
   print*, 'init_fluxcells ',iw,nwa
endif

! Loop over all WL points

   do iwl = 2,nwl

! new - skip interaction using new non-overlap check

      if (abs(xew(iw) + land%xew(iwl)) < 12.e6) then  ! Skip if both near X pole
         if (xeamin(iw ) > xelmax(iwl) + 10.) cycle
         if (xelmin(iwl) > xeamax(iw ) + 10.) cycle
      endif

      if (abs(yew(iw) + land%yew(iwl)) < 12.e6) then  ! Skip if both near Y pole
         if (yeamin(iw ) > yelmax(iwl) + 10.) cycle
         if (yelmin(iwl) > yeamax(iw ) + 10.) cycle
      endif

      if (abs(zew(iw) + land%zew(iwl)) < 12.e6) then  ! Skip if both near Z pole
         if (zeamin(iw ) > zelmax(iwl) + 10.) cycle
         if (zelmin(iwl) > zeamax(iw ) + 10.) cycle
      endif

! Evaluate x,y coordinates of current W point on polar stereographic plane
! tangent at IWL

      call e_ps(xew(iw), yew(iw), zew(iw), &
                land%glatw(iwl), land%glonw(iwl), xp0, yp0)

      xw = real(xp0,r8)
      yw = real(yp0,r8)

! Loop over all neighbor M points of this IW

      npoly = itab_w(iw)%npoly

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

! Evaluate x,y coordinates of current M point on polar stereographic plane
! tangent at IWL

         call e_ps(xem(im), yem(im), zem(im), &
                   land%glatw(iwl), land%glonw(iwl), xp0, yp0)

         xp(jm) = real(xp0,r8)
         yp(jm) = real(yp0,r8)
      enddo

! Loop over all neighbor M points of this IWL

      nlpoly = itab_wl(iwl)%npoly

      zmin = 1.e6

      do jml = 1,nlpoly
         iml = itab_wl(iwl)%im(jml)

! Evaluate x,y coordinates of LAND cell M points on polar stereographic plane
! tangent at IWL

         call e_ps(land%xem(iml), land%yem(iml), land%zem(iml), &
                   land%glatw(iwl), land%glonw(iwl), xq0, yq0)

         xq(jml) = real(xq0,r8)
         yq(jml) = real(yq0,r8)

! Store topography height of M point in zq array

         zq(jml) = real(land%zm(iml),r8)

         if (zmin > land%zm(iml)) zmin = land%zm(iml)
      enddo

! Find k level of T cell that intersects topography in this IWL cell
! (In makesfc when topography is within top 1% of a model dzt layer, it is
! moved up to the next zm level, so 0.995 is good threshold to check here.)

      do k = 2,nza-1
         fraczmin = (zmin - zm(k-1)) / (zm(k) - zm(k-1))

         if (fraczmin < 0.995) then
            kw0 = k
            exit
         endif
      enddo

! Evaluate possible overlap of ATM and LAND polygons

      call polygon_overlap(iwl,npoly,nlpoly,xp,yp,xq,yq,area, &
                           jtrap,xtrap,ytrap,traparea)

! Skip iwl if overlap is zero

      if (area < 1.0e-7) cycle

! This IW polygon overlaps with land cell IWL.

! Find fit coefficients for linear elevation surface of land cell IWL
! using 3 vertices of IWL polygon. (It is assumed that IWL polygon is convex, 
! so that selected points are not collinear.)

      j3 = 3
      if (nlpoly > 4) j3 = 4

      call matrix8_3x3(1.0_r8,xq(1),yq(1)   &
                      ,1.0_r8,xq(2),yq(2)   &
                      ,1.0_r8,xq(j3),yq(j3) &
                      ,zq(1),zq(2),zq(j3)   &
                      ,az,bz,cz             )

! Determine TOPW by checking whether IW point is inside (or within 10 m of) 
! current IWL polygon

      call inout_check(nlpoly,xq,yq,xw,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nlpoly,xq,yq,xw-10.0_r8,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nlpoly,xq,yq,xw+10.0_r8,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nlpoly,xq,yq,xw,yw-10.0_r8,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nlpoly,xq,yq,xw,yw+10.0_r8,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

! Determine TOPM values of IW cell by checking whether each IM point is inside
! (or within 10 m of) current IWL polygon

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nlpoly,xq,yq,xp(jm)-10.0_r8,yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nlpoly,xq,yq,xp(jm)+10.0_r8,yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm)-10.0_r8,alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm)+10.0_r8,alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))
         
      enddo

! Skip remaining calculations if overlap is small

      if (area < 1.e-6 * min(arw0(iw),land%area(iwl))) cycle

! Find height of 4 corners of all trapezoids of current IW-IWL overlap,
! and find their overall maximum and minimum values

      zmax = -1.e6
      zmin = 1.e6

      do jt = 1,jtrap
         do jc = 1,4
            ztrap(jc,jt) = real(az) &
                         + real(bz) * xtrap(jc,jt) + real(cz) * ytrap(jc,jt)

            if (zmax < ztrap(jc,jt)) zmax = ztrap(jc,jt)
            if (zmin > ztrap(jc,jt)) zmin = ztrap(jc,jt)
         enddo
      enddo

! Add full-cell contribution to ARW for levels kw0 and above

      arw(kw0:nza-1,iw) = arw(kw0:nza-1,iw) + area

! Add full-cell contribution to VOLT for levels above kw0

      volt(kw0+1:nza-1,iw) = volt(kw0+1:nza-1,iw) + area * dzt(kw0+1:nza-1)

! Make sure that topo height range is (approximately) within height range of
! T(kw0) level.  This is based on new requirement (July 2013) that vertical
! grid be identical between 'MAKESFC' and 'MAKEGRID' runs.

      if (zmin < zm(kw0-1) - 2.0 .or. zmax > zm(kw0) + 1.0) then 
         print*, 'Topography range in flux cell outside DZT range. '
         print*, 'IW, IWL, KW0, zmin, zmax, zm(kw0-1), zm(kw0) '
         print*,  iw, iwl, kw0, zmin, zmax, zm(kw0-1), zm(kw0)

         stop 'stop flux cell topography range '
      endif

! Numerically integrate over all trapezoids of this IW-IWL overlap
! to get VOLT for level kw0 and to get barycenter location of flux cell

      xmean = 0.
      ymean = 0.

! Loop over all trapezoids of current IW-IWL overlap

      do jt = 1,jtrap

! Numerically integrate over trapezoid to get VOLT, xmean, and ymean

         dy = (ytrap(4,jt) - ytrap(1,jt)) / real(nsub)

         do jsub = 1,nsub

            fracy = (real(jsub) - 0.5) / real(nsub)

            ycent = ytrap(1,jt) + (ytrap(4,jt) - ytrap(1,jt)) * fracy

            x1 = xtrap(1,jt) + (xtrap(4,jt) - xtrap(1,jt)) * fracy
            x2 = xtrap(2,jt) + (xtrap(3,jt) - xtrap(2,jt)) * fracy

            dx = (x2 - x1) / real(nsub)

            do isub = 1,nsub

               fracx = (real(isub) - 0.5) / real(nsub)

               xcent = x1 + (x2 - x1) * fracx

               topp = real(az) + real(bz) * xcent + real(cz) * ycent

               fractopp = (topp - zm(kw0-1)) / (zm(kw0) - zm(kw0-1))

! Check that topp is within (zm(k-1),zm(k)) range

               if (fractopp < -.02 .or. fractopp > 1.02) then

                  print*, 'topp outside DZT range. '
                  print*, 'IW, IWL, KW0, jsub, isub '
                  print*,  iw, iwl, kw0, jsub, isub
                  print*, 'zmin, zmax, topp '
                  print*,  zmin, zmax, topp
                  print*, 'zm(kw0-1), zm(kw0), fractopp '
                  print*,  zm(kw0-1), zm(kw0), fractopp

                  stop 'stop topp '
               endif

               volt(kw0,iw) = volt(kw0,iw) &
                            + dx * dy * (zm(kw0) - max(topp,zm(kw0-1)))
               
               xmean = xmean + xcent * dx * dy
               ymean = ymean + ycent * dx * dy

            enddo  ! isub

         enddo  ! jsub

      enddo  ! jt

! Set new landflux cell index, allocate more space if necessary,
! and fill landflux values

      ilf = ilf + 1

      ifsize = size(landflux)
      if (ilf > ifsize) then
         allocate ( landflux_temp(ifsize+incr_flux) )
         landflux_temp(1:ifsize) = landflux(1:ifsize)
         call move_alloc(landflux_temp, landflux)
      endif

! Fill M points for landflux cell from trapezoid corners, skipping repeated pts

      do ktrap = 1,jtrap
         kpoly = ktrap

         xm(kpoly) = xtrap(2,ktrap)
         ym(kpoly) = ytrap(2,ktrap)
      enddo

      kpoly = kpoly + 1

      xm(kpoly) = xtrap(3,jtrap)
      ym(kpoly) = ytrap(3,jtrap)

      dist = sqrt((xtrap(4,jtrap) - xtrap(3,jtrap))**2 &
                + (ytrap(4,jtrap) - ytrap(3,jtrap))**2)

      if (dist > 1.0) then
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(4,jtrap)
         ym(kpoly) = ytrap(4,jtrap)
      endif

      do ktrap = jtrap,2,-1
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(1,ktrap)
         ym(kpoly) = ytrap(1,ktrap)
      enddo

      dist = sqrt((xtrap(2,1) - xtrap(1,1))**2 &
                + (ytrap(2,1) - ytrap(1,1))**2)

      if (dist > 1.0) then
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(1,1)
         ym(kpoly) = ytrap(1,1)
      endif

! Check for collinearity

      jpoly = 0

      do j2 = 1,kpoly
         j1 = j2 - 1
         if (j2 == 1) j1 = kpoly
         j3 = j2 + 1
         if (j2 == kpoly) j3 = 1

         crossprod = xm(j1) * (ym(j2) - ym(j3)) &
                   + xm(j2) * (ym(j3) - ym(j1)) &
                   + xm(j3) * (ym(j1) - ym(j2))

         dist12 = sqrt((xm(j2)-xm(j1))**2 + (ym(j2)-ym(j1))**2)
         dist23 = sqrt((xm(j3)-xm(j2))**2 + (ym(j3)-ym(j2))**2)

! Assume that points are not collinear if either the direction change is
! larger than 1.e-2 radians or middle point is more than 1.e0 m out of line

         if (abs(crossprod) > 1.e-2 * dist12 * dist23 .or. &
             abs(crossprod) > 1.e0  * max(dist12,dist23)) then

            jpoly = jpoly + 1

            call ps_e(landflux(ilf)%xem(jpoly), &
                      landflux(ilf)%yem(jpoly), &
                      landflux(ilf)%zem(jpoly), &
                      land%glatw(iwl),land%glonw(iwl), &
                      xm(j2),ym(j2))

         endif
            
      enddo

      call ps_e(landflux(ilf)%xef, landflux(ilf)%yef, landflux(ilf)%zef, &
                land%glatw(iwl), land%glonw(iwl), xmean/area, ymean/area)

      raxis = sqrt( landflux(ilf)%xef**2 + landflux(ilf)%yef**2 )

      landflux(ilf)%ifglobe = ilf ! full domain index
      landflux(ilf)%iw      = iw  ! full domain index
      landflux(ilf)%kw      = kw0
      landflux(ilf)%iwls    = iwl ! full domain index
      landflux(ilf)%npoly   = jpoly
      landflux(ilf)%area    = area
      landflux(ilf)%arf_atm = area / arw0(iw)
      landflux(ilf)%arf_sfc = area / land%area(iwl)
      landflux(ilf)%glatf   = atan2(landflux(ilf)%zef,raxis) * piu180
      landflux(ilf)%glonf   = atan2(landflux(ilf)%yef,landflux(ilf)%xef) * piu180

   enddo  ! iwl loop

! Loop over all WS points

   do iws = 2,nws

! new - skip interaction using new non-overlap check

      if (abs(xew(iw) + sea%xew(iws)) < 12.e6) then  ! Skip if both near X pole
         if (xeamin(iw ) > xesmax(iws) + 10.) cycle
         if (xesmin(iws) > xeamax(iw ) + 10.) cycle
      endif

      if (abs(yew(iw) + sea%yew(iws)) < 12.e6) then  ! Skip if both near Y pole
         if (yeamin(iw ) > yesmax(iws) + 10.) cycle
         if (yesmin(iws) > yeamax(iw ) + 10.) cycle
      endif

      if (abs(zew(iw) + sea%zew(iws)) < 12.e6) then  ! Skip if both near Z pole
         if (zeamin(iw ) > zesmax(iws) + 10.) cycle
         if (zesmin(iws) > zeamax(iw ) + 10.) cycle
      endif

! Evaluate x,y coordinates of current W point on polar stereographic plane
! tangent at IWS

      call e_ps(xew(iw), yew(iw), zew(iw), &
                sea%glatw(iws), sea%glonw(iws), xp0, yp0)

      xw = real(xp0,r8)
      yw = real(yp0,r8)

! Loop over all neighbor M points of this IW

      npoly = itab_w(iw)%npoly

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

! Evaluate x,y coordinates of current M point on polar stereographic plane
! tangent at IWS

         call e_ps(xem(im), yem(im), zem(im), &
                   sea%glatw(iws), sea%glonw(iws), xp0, yp0)

         xp(jm) = real(xp0,r8)
         yp(jm) = real(yp0,r8)
      enddo

! Loop over all neighbor M points of this IWS

      nspoly = itab_ws(iws)%npoly

      zmin = 1.e6

      do jms = 1,nspoly
         ims = itab_ws(iws)%im(jms)

! Evaluate x,y coordinates of SEA cell M points on polar stereographic plane
! tangent at IWS

         call e_ps(sea%xem(ims), sea%yem(ims), sea%zem(ims), &
                   sea%glatw(iws), sea%glonw(iws), xq0, yq0)

         xq(jms) = real(xq0,r8)
         yq(jms) = real(yq0,r8)

! Store topography height of M point in zq array

         zq(jms) = real(sea%zm(ims),r8)

         if (zmin > sea%zm(ims)) zmin = sea%zm(ims)
      enddo

! Find k level of T cell that intersects topography in this IWS cell
! (In makesfc when topography is within top 1% of a model dzt layer, it is
! moved up to the next zm level, so 0.995 is good threshold to check here.)

      do k = 2,nza-1
         fraczmin = (zmin - zm(k-1)) / (zm(k) - zm(k-1))

         if (fraczmin < 0.995) then
            kw0 = k
            exit
         endif
      enddo

! Evaluate possible overlap of ATM and SEA polygons

      call polygon_overlap(iws,npoly,nspoly,xp,yp,xq,yq,area, &
                           jtrap,xtrap,ytrap,traparea)

! Skip iwl if overlap is zero

      if (area < 1.0e-7) cycle

! This IW polygon overlaps with sea cell IWS.

! Find fit coefficients for linear elevation surface of sea cell IWS
! using 3 vertices of IWS polygon. (It is assumed that IWS polygon is convex, 
! so that selected points are not collinear.)

      j3 = 3
      if (nspoly > 4) j3 = 4

      call matrix8_3x3(1.0_r8,xq(1),yq(1)   &
                      ,1.0_r8,xq(2),yq(2)   &
                      ,1.0_r8,xq(j3),yq(j3) &
                      ,zq(1),zq(2),zq(j3)   &
                      ,az,bz,cz             )

! Determine TOPW by checking whether IW point is inside (or within 10 of)
! current IWS polygon

      call inout_check(nspoly,xq,yq,xw,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nspoly,xq,yq,xw-10.0_r8,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nspoly,xq,yq,xw+10.0_r8,yw,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nspoly,xq,yq,xw,yw-10.0_r8,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

      call inout_check(nspoly,xq,yq,xw,yw+10.0_r8,alpha)
      if (alpha > 1.0_r8) topw(iw) = real(az + bz * xw + cz * yw)

! Determine TOPM values of IW cell by checking whether each IM point is inside
! (or within 10 m of) current IWS polygon

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nspoly,xq,yq,xp(jm)-10.0_r8,yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nspoly,xq,yq,xp(jm)+10.0_r8,yp(jm),alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm)-10.0_r8,alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm)+10.0_r8,alpha)
         if (alpha > 1.0_r8) topm(im) = real(az + bz * xp(jm) + cz * yp(jm))
         
      enddo

! Skip remaining calculations if overlap is small

      if (area < 1.e-6 * min(arw0(iw),sea%area(iws))) cycle

! Find height of 4 corners of all trapezoids of current IW-IWL overlap,
! and find their overall maximum and minimum values

      zmax = -1.e6
      zmin = 1.e6

      do jt = 1,jtrap
         do jc = 1,4
            ztrap(jc,jt) = real(az) &
                         + real(bz) * xtrap(jc,jt) + real(cz) * ytrap(jc,jt)

            if (zmax < ztrap(jc,jt)) zmax = ztrap(jc,jt)
            if (zmin > ztrap(jc,jt)) zmin = ztrap(jc,jt)
         enddo
      enddo

! Add full-cell contribution to ARW for levels kw0 and above

      arw(kw0:nza-1,iw) = arw(kw0:nza-1,iw) + area

! Add full-cell contribution to VOLT for levels above kw0

      volt(kw0+1:nza-1,iw) = volt(kw0+1:nza-1,iw) + area * dzt(kw0+1:nza-1)

! Make sure that topo height range is (approximately) within height range of
! T(kw0) level.  This is based on new requirement (July 2013) that vertical
! grid be identical between 'MAKESFC' and 'MAKEGRID' runs.

      if (zmin < zm(kw0-1) - 2.0 .or. zmax > zm(kw0) + 1.0) then 
         print*, 'Topography range in flux cell outside DZT range. '
         print*, 'IW, IWL, KW0, zmin, zmax, zm(kw0-1), zm(kw0) '
         print*,  iw, iwl, kw0, zmin, zmax, zm(kw0-1), zm(kw0)

         stop 'stop flux cell topography range '
      endif
 
! Numerically integrate over all trapezoids of this IW-IWS overlap
! to get VOLT for level kw0 and to get barycenter location of flux cell

      xmean = 0.
      ymean = 0.

! Loop over all trapezoids of current IW-IWS overlap

      do jt = 1,jtrap

! Numerically integrate over trapezoid to get VOLT, xmean, and ymean

         dy = (ytrap(4,jt) - ytrap(1,jt)) / real(nsub)

         do jsub = 1,nsub

            fracy = (real(jsub) - 0.5) / real(nsub)

            ycent = ytrap(1,jt) + (ytrap(4,jt) - ytrap(1,jt)) * fracy

            x1 = xtrap(1,jt) + (xtrap(4,jt) - xtrap(1,jt)) * fracy
            x2 = xtrap(2,jt) + (xtrap(3,jt) - xtrap(2,jt)) * fracy

            dx = (x2 - x1) / real(nsub)

            do isub = 1,nsub

               fracx = (real(isub) - 0.5) / real(nsub)

               xcent = x1 + (x2 - x1) * fracx

               topp = real(az) + real(bz) * xcent + real(cz) * ycent

               fractopp = (topp - zm(kw0-1)) / (zm(kw0) - zm(kw0-1))

! Check that topp is within (zm(k-1),zm(k)) range

               if (fractopp < -.01 .or. fractopp > 1.01) then

                  print*, 'topp outside DZT range. '
                  print*, 'IW, IWS, KW0, jsub, isub '
                  print*,  iw, iws, kw0, jsub, isub
                  print*, 'zmin, zmax, topp '
                  print*,  zmin, zmax, topp
                  print*, 'zm(kw0-1), zm(kw0) '
                  print*,  zm(kw0-1), zm(kw0)

                  stop 'stop topp '
               endif

               volt(kw0,iw) = volt(kw0,iw) &
                            + dx * dy * (zm(kw0) - max(topp,zm(kw0-1)))
               
               xmean = xmean + xcent * dx * dy
               ymean = ymean + ycent * dx * dy

            enddo  ! isub

         enddo  ! jsub

      enddo  ! jt

! Set new seaflux cell index, allocate more space if necessary,
! and fill seaflux values

      isf = isf + 1

      ifsize = size(seaflux)
      if (isf > ifsize) then
         allocate ( seaflux_temp(ifsize+incr_flux) )
         seaflux_temp(1:ifsize) = seaflux(1:ifsize)
         call move_alloc(seaflux_temp, seaflux)
      endif

! Fill M points for seaflux cell from trapezoid corners, skipping repeated pts

      do ktrap = 1,jtrap
         kpoly = ktrap

         xm(kpoly) = xtrap(2,ktrap)
         ym(kpoly) = ytrap(2,ktrap)
      enddo

      kpoly = kpoly + 1

      xm(kpoly) = xtrap(3,jtrap)
      ym(kpoly) = ytrap(3,jtrap)

      dist = sqrt((xtrap(4,jtrap) - xtrap(3,jtrap))**2 &
                + (ytrap(4,jtrap) - ytrap(3,jtrap))**2)

      if (dist > 1.0) then
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(4,jtrap)
         ym(kpoly) = ytrap(4,jtrap)
      endif

      do ktrap = jtrap,2,-1
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(1,ktrap)
         ym(kpoly) = ytrap(1,ktrap)
      enddo

      dist = sqrt((xtrap(2,1) - xtrap(1,1))**2 &
                + (ytrap(2,1) - ytrap(1,1))**2)

      if (dist > 1.0) then
         kpoly = kpoly + 1

         xm(kpoly) = xtrap(1,1)
         ym(kpoly) = ytrap(1,1)
      endif

! Check for collinearity

      jpoly = 0

      do j2 = 1,kpoly
         j1 = j2 - 1
         if (j2 == 1) j1 = kpoly
         j3 = j2 + 1
         if (j2 == kpoly) j3 = 1

         crossprod = xm(j1) * (ym(j2) - ym(j3)) &
                   + xm(j2) * (ym(j3) - ym(j1)) &
                   + xm(j3) * (ym(j1) - ym(j2))

         dist12 = sqrt((xm(j2)-xm(j1))**2 + (ym(j2)-ym(j1))**2)
         dist23 = sqrt((xm(j3)-xm(j2))**2 + (ym(j3)-ym(j2))**2)

! Assume that points are not collinear if either the direction change is
! larger than 1.e-2 radians or middle point is more than 1.e0 m out of line

         if (abs(crossprod) > 1.e-2 * dist12 * dist23 .or. &
             abs(crossprod) > 1.e0  * max(dist12,dist23)) then

            jpoly = jpoly + 1

            call ps_e(seaflux(isf)%xem(jpoly), &
                      seaflux(isf)%yem(jpoly), &
                      seaflux(isf)%zem(jpoly), &
                      sea%glatw(iws),sea%glonw(iws), &
                      xm(j2),ym(j2))

         endif
            
      enddo

      call ps_e(seaflux(isf)%xef, seaflux(isf)%yef, seaflux(isf)%zef, &
                sea%glatw(iws), sea%glonw(iws), xmean/area, ymean/area)

      raxis = sqrt( seaflux(isf)%xef**2 + seaflux(isf)%yef**2 )

      seaflux(isf)%ifglobe = isf ! full domain index
      seaflux(isf)%iw      = iw  ! full domain index
      seaflux(isf)%kw      = kw0
      seaflux(isf)%iwls    = iws ! full domain index
      seaflux(isf)%npoly   = jpoly
      seaflux(isf)%area    = area
      seaflux(isf)%arf_atm = area / arw0(iw)
      seaflux(isf)%arf_sfc = area / sea%area(iws)
      seaflux(isf)%glatf   = atan2(seaflux(isf)%zef,raxis) * piu180
      seaflux(isf)%glonf   = atan2(seaflux(isf)%yef,seaflux(isf)%xef) * piu180

   enddo ! iws loop

enddo ! iw loop

nseaflux = isf
nlandflux = ilf

mseaflux = isf
mlandflux = ilf

! Check sums of arf values to make sure they are close to 1.0

do iter = 1,2

   arf_atm_test(1:nwa) = 0.
   arf_sea_test(1:nws) = 0.
   arf_land_test(1:nwl) = 0.

   do isf = 2,nseaflux
      iw  = seaflux(isf)%iw
      iws = seaflux(isf)%iwls

      arf_atm_test(iw)  = arf_atm_test(iw)  + seaflux(isf)%arf_atm
      arf_sea_test(iws) = arf_sea_test(iws) + seaflux(isf)%arf_sfc
   enddo

   do ilf = 2,nlandflux
      iw  = landflux(ilf)%iw
      iwl = landflux(ilf)%iwls

      arf_atm_test(iw)   = arf_atm_test(iw)   + landflux(ilf)%arf_atm
      arf_land_test(iwl) = arf_land_test(iwl) + landflux(ilf)%arf_sfc
   enddo

   aatmin =  1.e9
   aatmax = -1.e9

   astmin =  1.e9
   astmax = -1.e9

   altmin =  1.e9
   altmax = -1.e9

   do iw = 2,nwa
      if (aatmin > arf_atm_test(iw)) aatmin = arf_atm_test(iw)
      if (aatmax < arf_atm_test(iw)) aatmax = arf_atm_test(iw)
   enddo

   do iws = 2,nws
      if (astmin > arf_sea_test(iws)) astmin = arf_sea_test(iws)
      if (astmax < arf_sea_test(iws)) astmax = arf_sea_test(iws)
   enddo

   do iwl = 2,nwl
      if (altmin > arf_land_test(iwl)) altmin = arf_land_test(iwl)
      if (altmax < arf_land_test(iwl)) altmax = arf_land_test(iwl)
   enddo

   if (iter == 1) then

      do isf = 2,nseaflux
         iw  = seaflux(isf)%iw
         iws = seaflux(isf)%iwls

         seaflux(isf)%arf_atm = seaflux(isf)%arf_atm / arf_atm_test(iw) 
         seaflux(isf)%arf_sfc = seaflux(isf)%arf_sfc / arf_sea_test(iws)
      enddo

      do ilf = 2,nlandflux
         iw  = landflux(ilf)%iw
         iwl = landflux(ilf)%iwls

         landflux(ilf)%arf_atm = landflux(ilf)%arf_atm  / arf_atm_test(iw) 
         landflux(ilf)%arf_sfc = landflux(ilf)%arf_sfc / arf_land_test(iwl)
      enddo
      
   endif

enddo

! Make sure that all topm and topw values have been assigned

istop = 0

do im = 2,nma
   if (topm(im) < -1.e3) then
      print*, 'In init_fluxcells, topm missing at im = ',im,topm(im),glatm(im),glonm(im)
      istop = 1
   endif
enddo

do iw = 2,nwa
   if (topw(iw) < -1.e3) then
      print*, 'In init_fluxcells, topw missing at iw = ',iw,topw(iw),glatw(iw),glonw(iw)
      istop = 1
   endif
enddo

if (istop > 0) stop 'stop_init_fluxcells_topo'

! If ilandgrid > 1, partition land grid cells from atm grid overlay
! and write supplementary landfile

if (ilandgrid > 1) then
   call partition_land_cells()
   call landfile_write()
endif

! If iseagrid > 1, partition sea grid cells based on atm grid overlay,
! and write supplementary seafile.

if (iseagrid > 1) then
   call partition_sea_cells()

! If iseagrid > 2, combine sea grid cells (and seaflux cells) that share a
! sea-only atm grid column if all do so at the same kw level.

   if (iseagrid > 2) call combine_sea_cells()

   call seafile_write()
endif

return
end subroutine init_fluxcells

!============================================================================

subroutine polygon_overlap(iwls,np,nq,xp,yp,xq,yq,area,ntrap,xtrap,ytrap,traparea)

! Given x,y coordinates of the vertices of polygons p and q, compute the area
! of overlap between the polygons using a sweepline algorithm.

! Method adapted from:
! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

use consts_coms, only: r8

implicit none

integer, intent(in) :: np,nq,iwls ! Number of vertices in p and q
real(r8), intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
real(r8), intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

integer, intent(out) :: ntrap              ! number of trapezoids
real, intent(out) :: xtrap(4,np+nq+np*nq)  ! trapezoid x coordinates
real, intent(out) :: ytrap(4,np+nq+np*nq)  ! trapezoid y coordinates
real, intent(out) :: traparea(np+nq+np*nq) ! trapezoid area

real, intent(out) :: area                  ! area of overlap of p and q

real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

integer :: nev      ! # of events
integer :: nsect    ! # of intersections between strip centerline and p,q edges
real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges 
real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges 
real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges 
real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values

integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,i,ia,ib,is
real(r8) :: p0,q0,dx,dy,dxtrap

real(r8) :: alpha    ! angle swept by ray from point to polygon perimeter during circuit

! Find vertices in p that are not outside q

nev = 0
ntrap = 0
area = 0.

do ip = 1,np
   call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

   if (abs(alpha) > 0.2_r8) then
      nev = nev + 1
      yev(nev) = yp(ip)
   endif
enddo

! Find vertices in q that are not outside p

do iq = 1,nq
   call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

   if (abs(alpha) > 0.2_r8) then
      nev = nev + 1
      yev(nev) = yq(iq)
   endif
enddo

! Find intersecting edges of polygons p and q

do ipa = 1,np
   ipb = ipa + 1
   if (ipa == np) ipb = 1

   do iqa = 1,nq
      iqb = iqa + 1
      if (iqa == nq) iqb = 1

      call intersect(xp(ipa),yp(ipa),xp(ipb),yp(ipb)  &
                    ,xq(iqa),yq(iqa),xq(iqb),yq(iqb),p0,q0,iflag)

      if (iflag == -1) cycle
      if (p0 < -0.000000001_r8 .or. p0 > 1.000000001_r8) cycle
      if (q0 < -0.000000001_r8 .or. q0 > 1.000000001_r8) cycle

! Line segments pa-pb and qa-qb intersect; find y-coordinate of intersection

      nev = nev + 1
      yev(nev) = (1.0_r8 - p0) * yp(ipa) + p0 * yp(ipb)
   enddo
enddo

if (nev == 0) then
!   print*, 'no overlap'
   return
endif

! Sort event list to obtain increasing order of yev

call fltsort(nev,yev)

! Loop over event points

do iev = 1,nev-1

! dy = width of strip between current and next event

   dy = yev(iev+1) - yev(iev)
   
! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)
   
   if (dy < 0.1_r8) cycle

! ymid = y-coordinate of strip centerline.
! Initialize dx (centerline length sum) to 0.

   ymid = yev(iev) + 0.5_r8 * dy
   dx = 0.0_r8

! Find x-coordinate of intersections of strip centerline with edges of p and q

   nsect = 0

   do ia = 1,np
      ib = ia + 1
      if (ia == np) ib = 1

      if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) / (yp(ib) - yp(ia))
      xbot(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) / (yp(ib) - yp(ia))
      xtop(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) / (yp(ib) - yp(ia))

   enddo

   do ia = 1,nq
      ib = ia + 1
      if (ia == nq) ib = 1

      if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) / (yq(ib) - yq(ia))
      xbot(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) / (yq(ib) - yq(ia))
      xtop(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) / (yq(ib) - yq(ia))

   enddo

! Sort xmid values into increasing order

   call fltsort3(nsect,xmid,xbot,xtop)

! see if the segment is inside both polygons

   do is = 1,nsect - 1
      xcent = 0.5_r8 * (xmid(is) + xmid(is+1))

      if (xcent == xmid(is)) cycle          ! if segment length is 0

      call inout_check(np,xp,yp,xcent,ymid,alpha)
      if (abs(alpha) < 1.0_r8) cycle

      call inout_check(nq,xq,yq,xcent,ymid,alpha)
      if (abs(alpha) < 1.0_r8) cycle

      dxtrap = xmid(is+1) - xmid(is)

      dx = dx + dxtrap

! Find x at 4 corners of current trapezoidal region

      ntrap = ntrap + 1
      
      xtrap(1,ntrap) = xbot(is)
      xtrap(2,ntrap) = xbot(is+1)
      xtrap(3,ntrap) = xtop(is+1)
      xtrap(4,ntrap) = xtop(is)
      
      ytrap(1,ntrap) = yev(iev) 
      ytrap(2,ntrap) = yev(iev) 
      ytrap(3,ntrap) = yev(iev+1) 
      ytrap(4,ntrap) = yev(iev+1) 

      traparea(ntrap) = dxtrap * dy
   enddo

   area = area + dx * dy

enddo

return
end subroutine polygon_overlap

!============================================================================

subroutine inout_check(n,x,y,x0,y0,th1)

! Given planar Cartesian coordinates of the vertices of a simple closed polygon,
! x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0, determine
! whether the point is inside, outside, or on the border of the polygon.

! Set:
! th1 = 0 if outside
! th1 = pi if on the border
! th1 = 2 * pi if inside

use consts_coms, only: r8

implicit none

integer, intent(in) :: n  ! number of points in polygon

real(r8), intent(in) :: x(n),y(n)  ! x,y coordinates of polygon points      
real(r8), intent(in) :: x0,y0      ! x,y coordinates of additional point
real(r8), intent(out) :: th1       ! output value

! Local variables

integer :: i
real(r8) :: theta
real(r8) :: xh0, xh1, xh2
real(r8) :: edgemax
real(r8) :: x1(n), y1(n)

real(r8), parameter :: pi1 = 3.1415926535898_r8, pi2 = 2.0_r8 * pi1

do i = 1,n
   x1(i) = x(i) - x0
   y1(i) = y(i) - y0

! Check whether x0,y0 is exactly on a vertex (to within 1.e-6 m)

   if ((abs(x1(i)) < 0.000001_r8) .and. (abs(y1(i)) < 0.000001_r8)) then
      th1 = pi1
      return
   endif
enddo

th1 = 0.0_r8

xh2 = atan2(y1(1),x1(1))
if (xh2 < 0.0_r8) xh2 = xh2 + pi2

xh0 = xh2

do i = 1,n
   
   if (i == n) then
      xh1 = xh0
   else
      xh1 = atan2(y1(i+1),x1(i+1))
      if (xh1 < 0.0_r8) xh1 = xh1 + pi2
   endif

   theta = xh1 - xh2
   if (theta < -pi1) theta = theta + pi2

   if ((abs(abs(theta) - pi1) < 0.000000001_r8)) then
      th1 = pi1
      return
   endif

   if (theta > pi1) theta = theta - pi2
   th1 = th1 + theta
   
   xh2 = xh1
enddo

th1 = abs(th1)

return
end subroutine inout_check

!============================================================================

subroutine intersect(xpa,ypa,xpb,ypb,xqa,yqa,xqb,yqb,p0,q0,iflag)

! Given x,y coordinates of points pa and pb on line p and points qa and qb on
! line q, find parameteric values p0 and q0 where p and q intersect.
! If no intersection, set iflag = -1.

use consts_coms, only: r8

implicit none

real(r8), intent(in) :: xpa,ypa,xpb,ypb  ! x,y coordinates of pa and pb
real(r8), intent(in) :: xqa,yqa,xqb,yqb  ! x,y coordinates of qa and qb

integer, intent(out) :: iflag
real(r8), intent(out) :: p0,q0

real(r8) :: pabxqab  ! cross product of vectors pa-pb and qa-qb 

iflag = -1

pabxqab = (xpb - xpa) * (yqb - yqa) - (ypb - ypa) * (xqb - xqa)

if (pabxqab /= 0.0_r8) then  ! Infinite lines intersect
   iflag = 0
   if (pabxqab < 0.0_r8) iflag = 1

   p0 = ((yqb - yqa) * (xqa - xpa) + (yqa - ypa) * (xqa - xqb)) / pabxqab
   q0 = ((ypb - ypa) * (xqa - xpa) + (yqa - ypa) * (xpa - xpb)) / pabxqab
endif

!print*, 'pabxqab ',pabxqab,p0,q0,iflag

return
end subroutine intersect

!============================================================================

subroutine fltsort(n,f)

! Sort n floating point numbers f into ascending order

use consts_coms, only: r8

implicit none

integer, intent(in) :: n
real(r8), intent(inout) :: f(n)

integer :: i,j

real(r8) :: f0

do i = 1,n-1
   do j = i+1,n
      if (f(j) < f(i)) then
         f0 = f(i)
         f(i) = f(j)
         f(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort
   
!============================================================================

subroutine fltsort2(n,f1,f2)

! Sort n floating point numbers in each of f1 and f2 into ascending order by f1

use consts_coms, only: r8

implicit none

integer, intent(in) :: n
real(r8), intent(inout) :: f1(n),f2(n)

integer :: i,j

real(r8) :: f0

do i = 1,n-1
   do j = i+1,n
      if (f1(j) < f1(i)) then
         f0 = f1(i)
         f1(i) = f1(j)
         f1(j) = f0

         f0 = f2(i)
         f2(i) = f2(j)
         f2(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort2

!============================================================================

subroutine fltsort3(n,f1,f2,f3)

! Sort n floating point numbers in each of f1, f2, and f3 into ascending order by f1

use consts_coms, only: r8

implicit none

integer, intent(in) :: n
real(r8), intent(inout) :: f1(n),f2(n),f3(n)

integer :: i,j

real(r8) :: f0

do i = 1,n-1
   do j = i+1,n
      if (f1(j) < f1(i)) then
         f0 = f1(i)
         f1(i) = f1(j)
         f1(j) = f0

         f0 = f2(i)
         f2(i) = f2(j)
         f2(j) = f0

         f0 = f3(i)
         f3(i) = f3(j)
         f3(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort3

!===============================================================================

subroutine matrix8_3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3,x1,x2,x3)

use consts_coms, only: r8

implicit none

real(r8), intent(in) :: a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3
real(r8), intent(out) :: x1,x2,x3

integer :: i
real(r8), dimension(4,3) :: abr

abr(1,1) = a11
abr(2,1) = a21
abr(3,1) = a31
abr(4,1) = b1

abr(1,2) = a12
abr(2,2) = a22
abr(3,2) = a32
abr(4,2) = b2

abr(1,3) = a13
abr(2,3) = a23
abr(3,3) = a33
abr(4,3) = b3

! Interchange rows if necessary so that first row has 
! largest (magnitude) element of first column

if (abs(abr(1,2)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange8(abr(i,1),abr(i,2))
   enddo
endif

if (abs(abr(1,3)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange8(abr(i,1),abr(i,3))
   enddo
endif

! Add -abr(1,2)/abr(1,1) times first row to second row and
! add -abr(1,3)/abr(1,1) times first row to third row.

do i = 2,4
   abr(i,2) = abr(i,2) - abr(1,2)/abr(1,1)*abr(i,1)
   abr(i,3) = abr(i,3) - abr(1,3)/abr(1,1)*abr(i,1)
enddo

! Interchange rows 2 and 3 if necessary so that second row
! has larger (magnitude) element of second column

if (abs(abr(2,3)) > abs(abr(2,2))) then
   do i = 2,4
      call rchange8(abr(i,2),abr(i,3))
   enddo
endif

! Add -abr(2,3)/abr(2,2) times second row to third row.

do i = 3,4
   abr(i,3) = abr(i,3) - abr(2,3)/abr(2,2)*abr(i,2)
enddo

! Back substitution

x3 = abr(4,3) / abr(3,3)
x2 = (abr(4,2) - abr(3,2) * x3) / abr(2,2)
x1 = (abr(4,1) - abr(2,1) * x2 - abr(3,1) * x3) / abr(1,1)

return
end subroutine matrix8_3x3

!===============================================================================

  subroutine rchange8(r1,r2)

  use consts_coms, only: r8

  implicit none

  real(r8), intent(inout) :: r1,r2
  real(r8) :: c

  c = r1
  r1 = r2
  r2 = c

  return
  end subroutine rchange8

!============================================================================

  subroutine partition_land_cells()

  use mem_leaf,   only: land, land_vars, itab_ml, itab_ul, itab_wl
  use mem_ijtabs, only: itab_w, mrls
  use leaf_coms,  only: nzg, nwl, ilandgrid
  use misc_coms,  only: io6
  use mem_mksfc,  only: itab_wls_vars

  implicit none

  type(itab_wls_vars), allocatable :: ltab_wl(:)

  type(land_vars) :: land_t

  integer :: ilf, iwl

! Deallocate arrays that are not used with ILANDGRID > 1

  deallocate(land%xem)
  deallocate(land%yem)
  deallocate(land%zem)
  deallocate(land%zm)
  deallocate(land%glatm)
  deallocate(land%glonm)

  deallocate(itab_ml)
  deallocate(itab_ul)

! Allocate temporary land arrays to size of landflux arrays

  allocate (ltab_wl(nlandflux))

  allocate (land_t%area      (nlandflux))
  allocate (land_t%ntext_soil(nzg,nlandflux))
  allocate (land_t%leaf_class(nlandflux))
  allocate (land_t%xew       (nlandflux))
  allocate (land_t%yew       (nlandflux))
  allocate (land_t%zew       (nlandflux))
  allocate (land_t%glatw     (nlandflux))
  allocate (land_t%glonw     (nlandflux))
  allocate (land_t%wnx       (nlandflux))
  allocate (land_t%wny       (nlandflux))
  allocate (land_t%wnz       (nlandflux))

! Loop over all landflux cells

  do ilf = 2,nlandflux
     iwl = landflux(ilf)%iwls

! Fill new land values from a combination of landflux and old land values

     ltab_wl(ilf)%npoly = landflux(ilf)%npoly

     land_t%area (ilf) = landflux(ilf)%area
     land_t%xew  (ilf) = landflux(ilf)%xef
     land_t%yew  (ilf) = landflux(ilf)%yef
     land_t%zew  (ilf) = landflux(ilf)%zef
     land_t%glatw(ilf) = landflux(ilf)%glatf
     land_t%glonw(ilf) = landflux(ilf)%glonf

     land_t%ntext_soil(1:nzg,ilf) = land%ntext_soil(1:nzg,iwl)
     land_t%leaf_class(ilf) = land%leaf_class(iwl)

     land_t%wnx(ilf) = land%wnx(iwl)
     land_t%wny(ilf) = land%wny(iwl)
     land_t%wnz(ilf) = land%wnz(iwl)

! Reset landflux land cell index and surface area-fraction

     landflux(ilf)%iwls = ilf
     landflux(ilf)%arf_sfc = 1.0

  enddo  ! ilf

  nwl = nlandflux

! Move temporary arrays to original leaf arrays

  call move_alloc(ltab_wl, itab_wl)

  call move_alloc(land_t%area,       land%area)
  call move_alloc(land_t%ntext_soil, land%ntext_soil)
  call move_alloc(land_t%leaf_class, land%leaf_class)
  call move_alloc(land_t%xew,        land%xew)
  call move_alloc(land_t%yew,        land%yew)
  call move_alloc(land_t%zew,        land%zew)
  call move_alloc(land_t%glatw,      land%glatw)
  call move_alloc(land_t%glonw,      land%glonw)
  call move_alloc(land_t%wnx,        land%wnx)
  call move_alloc(land_t%wny,        land%wny)
  call move_alloc(land_t%wnz,        land%wnz)

  return
  end subroutine partition_land_cells

!============================================================================

  subroutine partition_sea_cells()

  use mem_sea,    only: sea, sea_vars, itab_ms, itab_us, itab_ws
  use mem_ijtabs, only: itab_w, mrls
  use sea_coms,   only: nws, iseagrid
  use misc_coms,  only: io6
  use mem_mksfc,  only: itab_wls_vars

  implicit none

  type(itab_wls_vars), allocatable :: ltab_ws(:)

  type(sea_vars) :: sea_t

  integer :: isf, iws

! Deallocate arrays that are not used with ISEAGRID > 1

  deallocate(sea%xem)
  deallocate(sea%yem)
  deallocate(sea%zem)
  deallocate(sea%zm)
  deallocate(sea%glatm)
  deallocate(sea%glonm)

  deallocate(itab_ms)
  deallocate(itab_us)

! Allocate temporary sea arrays to size of seaflux arrays

  allocate (ltab_ws(nseaflux))

  allocate (sea_t%area      (nseaflux))
  allocate (sea_t%leaf_class(nseaflux))
  allocate (sea_t%xew       (nseaflux))
  allocate (sea_t%yew       (nseaflux))
  allocate (sea_t%zew       (nseaflux))
  allocate (sea_t%glatw     (nseaflux))
  allocate (sea_t%glonw     (nseaflux))

! Loop over all seaflux cells

  do isf = 2,nseaflux
     iws = seaflux(isf)%iwls

! Fill new sea values from a combination of seaflux and old sea values

     ltab_ws(isf)%npoly = seaflux(isf)%npoly

     sea_t%area (isf) = seaflux(isf)%area
     sea_t%xew  (isf) = seaflux(isf)%xef
     sea_t%yew  (isf) = seaflux(isf)%yef
     sea_t%zew  (isf) = seaflux(isf)%zef
     sea_t%glatw(isf) = seaflux(isf)%glatf
     sea_t%glonw(isf) = seaflux(isf)%glonf

     sea_t%leaf_class(isf) = sea%leaf_class(iws)

! Reset seaflux sea cell index and surface area-fraction

     seaflux(isf)%iwls = isf
     seaflux(isf)%arf_sfc = 1.0

  enddo  ! isf

  nws = nseaflux

! Move temporary arrays to original sea arrays

  call move_alloc(ltab_ws, itab_ws)

  call move_alloc(sea_t%area,       sea%area)
  call move_alloc(sea_t%leaf_class, sea%leaf_class)
  call move_alloc(sea_t%xew,        sea%xew)
  call move_alloc(sea_t%yew,        sea%yew)
  call move_alloc(sea_t%zew,        sea%zew)
  call move_alloc(sea_t%glatw,      sea%glatw)
  call move_alloc(sea_t%glonw,      sea%glonw)

  return
  end subroutine partition_sea_cells

!============================================================================

  subroutine combine_sea_cells()

  use mem_sea,    only: sea, sea_vars, itab_ms, itab_us, itab_ws
  use mem_ijtabs, only: itab_w, mrls
  use sea_coms,   only: nws, iseagrid
  use misc_coms,  only: io6
  use mem_mksfc,  only: itab_wls_vars
  use mem_grid,   only: nwa, xem, yem, zem, xew, yew, zew, glatw, glonw, arw0

  implicit none

  type(sflux_vars), allocatable :: seaflux_temp(:)

  type(itab_wls_vars), allocatable :: ltab_ws(:)

  type(sea_vars) :: sea_t

  integer :: isf, iws, ilf, iw, kw, jsf, jpoly, im

  integer, allocatable :: iwflag(:)

! Loop over all landflux cells, and flag the atmosphere grid IW columns that
! are attached to them

  allocate (iwflag(nwa))

  iwflag(:) = 0

  do ilf = 2,nlandflux
     iw = landflux(ilf)%iw

     if (iwflag(iw) == 0) iwflag(iw) = -1
  enddo

! Loop over all seaflux cells, and for those attached to an atmosphere grid
! column that has no landflux cell attached, compare seaflux cell kw levels
! within the same atmosphere grid column.  If any differ, flag atmosphere
! grid column with -1 value.

  do isf = 2,nseaflux
     iw = seaflux(isf)%iw
     kw = seaflux(isf)%kw

     if (iwflag(iw) == 0) then
        iwflag(iw) = kw
     elseif (iwflag(iw) > 0 .and. iwflag(iw) /= kw) then
        iwflag(iw) = -1
     endif
  enddo

! Loop over all seaflux cells and count the number (jsf) to be retained

  jsf = 1

  do isf = 2,nseaflux
     iw  = seaflux(isf)%iw

     if (iwflag(iw) == -1) then

! This seaflux cell shares an atmospheric grid column with one or more landflux
! cells and/or with seaflux cells having a different kw value, so increment jsf
! to retain this seaflux cell as is

        jsf = jsf + 1

     elseif (iwflag(iw) /= -2) then

! This seaflux cell shares an atmospheric grid column only with other seaflux
! cells, and all have the same kw value, so increment jsf only on first
! encounter with atmospheric cell, and reflag atmosphere cell (with -2 value)
! as having been visited

        jsf = jsf + 1
        iwflag(iw) = -2

     endif

  enddo

! Allocate temporary seaflux and sea arrays to new size for permanent arrays

  allocate (seaflux_temp(jsf))

  allocate (ltab_ws(jsf))

  allocate (sea_t%area      (jsf))
  allocate (sea_t%leaf_class(jsf))
  allocate (sea_t%xew       (jsf))
  allocate (sea_t%yew       (jsf))
  allocate (sea_t%zew       (jsf))
  allocate (sea_t%glatw     (jsf))
  allocate (sea_t%glonw     (jsf))

print*, 'csc1 ',nseaflux,jsf,size(sea_t%area)

! Loop again over all seaflux cells, using same logic for incrementing jsf,
! but this time fill temporary seaflux and sea array data values

  jsf = 1

  do isf = 2,nseaflux
     iw  = seaflux(isf)%iw
     iws = seaflux(isf)%iwls

     if (iwflag(iw) == -1) then

! This seaflux cell shares an atmospheric grid column with one or more landflux
! cells and/or with seaflux cells having a different kw value, so increment jsf
! to retain this seaflux cell as is

        jsf = jsf + 1

! Copy seaflux values to jsf index in temporary arrays

        do jpoly = 1,seaflux(isf)%npoly
           seaflux_temp(jsf)%xem(jpoly) = seaflux(isf)%xem(jpoly)
           seaflux_temp(jsf)%yem(jpoly) = seaflux(isf)%yem(jpoly)
           seaflux_temp(jsf)%zem(jpoly) = seaflux(isf)%zem(jpoly)
        enddo

        seaflux_temp(jsf)%xef = seaflux(isf)%xef
        seaflux_temp(jsf)%yef = seaflux(isf)%yef
        seaflux_temp(jsf)%zef = seaflux(isf)%zef

        seaflux_temp(jsf)%ifglobe = jsf
        seaflux_temp(jsf)%iw      = seaflux(isf)%iw
        seaflux_temp(jsf)%kw      = seaflux(isf)%kw
        seaflux_temp(jsf)%iwls    = jsf
        seaflux_temp(jsf)%npoly   = seaflux(isf)%npoly
        seaflux_temp(jsf)%area    = seaflux(isf)%area
        seaflux_temp(jsf)%arf_atm = seaflux(isf)%arf_atm
        seaflux_temp(jsf)%arf_sfc = seaflux(isf)%arf_sfc
        seaflux_temp(jsf)%glatf   = seaflux(isf)%glatf
        seaflux_temp(jsf)%glonf   = seaflux(isf)%glonf

! Copy sea values to jsf index in temporary arrays

        ltab_ws(jsf)%npoly = itab_ws(iws)%npoly

        sea_t%area (jsf) = seaflux(isf)%area
        sea_t%xew  (jsf) = seaflux(isf)%xef
        sea_t%yew  (jsf) = seaflux(isf)%yef
        sea_t%zew  (jsf) = seaflux(isf)%zef
        sea_t%glatw(jsf) = seaflux(isf)%glatf
        sea_t%glonw(jsf) = seaflux(isf)%glonf

        sea_t%leaf_class(jsf) = sea%leaf_class(iws)

     elseif (iwflag(iw) /= -3) then

! This seaflux cell shares an atmospheric grid column only with other seaflux
! cells, and all have the same kw value, so increment jsf only on first
! encounter with atmospheric cell, and reflag atmosphere cell (with -3 value
! this time because -2 was used previously) as having been visited

        jsf = jsf + 1
        iwflag(iw) = -3

! Copy seaflux values to jsf index in temporary arrays

        do jpoly = 1,itab_w(iw)%npoly
           im = itab_w(iw)%im(jpoly)

           seaflux_temp(jsf)%xem(jpoly) = xem(im)
           seaflux_temp(jsf)%yem(jpoly) = yem(im)
           seaflux_temp(jsf)%zem(jpoly) = zem(im)
        enddo

        seaflux_temp(jsf)%xef = xew(iw)
        seaflux_temp(jsf)%yef = yew(iw)
        seaflux_temp(jsf)%zef = zew(iw)

        seaflux_temp(jsf)%ifglobe  = jsf
        seaflux_temp(jsf)%iw       = seaflux(isf)%iw
        seaflux_temp(jsf)%kw       = seaflux(isf)%kw
        seaflux_temp(jsf)%iwls     = jsf
        seaflux_temp(jsf)%npoly    = itab_w(iw)%npoly
        seaflux_temp(jsf)%area     = arw0(iw)
        seaflux_temp(jsf)%arf_atm  = 1.0
        seaflux_temp(jsf)%arf_sfc  = 1.0
        seaflux_temp(jsf)%glatf    = glatw(iw)
        seaflux_temp(jsf)%glonf    = glonw(iw)

! Copy sea values to jsf index in temporary arrays

        ltab_ws(jsf)%npoly = itab_w(iw)%npoly

        sea_t%area (jsf) = arw0(iw)
        sea_t%xew  (jsf) = xew(iw)
        sea_t%yew  (jsf) = yew(iw)
        sea_t%zew  (jsf) = zew(iw)
        sea_t%glatw(jsf) = glatw(iw)
        sea_t%glonw(jsf) = glonw(iw)

        sea_t%leaf_class(jsf) = sea%leaf_class(iws)

     endif

  enddo

  nseaflux = jsf
  nws = jsf

! Move temporary arrays to original sea arrays

  call move_alloc(seaflux_temp, seaflux)

  call move_alloc(ltab_ws, itab_ws)

  call move_alloc(sea_t%area,       sea%area)
  call move_alloc(sea_t%leaf_class, sea%leaf_class)
  call move_alloc(sea_t%xew,        sea%xew)
  call move_alloc(sea_t%yew,        sea%yew)
  call move_alloc(sea_t%zew,        sea%zew)
  call move_alloc(sea_t%glatw,      sea%glatw)
  call move_alloc(sea_t%glonw,      sea%glonw)

print*, 'csc2 ',nseaflux,size(sea_t%area),size(sea%area)

  return
  end subroutine combine_sea_cells

End Module mem_sflux
