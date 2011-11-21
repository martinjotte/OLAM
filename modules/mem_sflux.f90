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

  integer :: nlfpats, mlfpats   ! number of patches for all land flux cells
  integer :: nsfpats, msfpats   ! number of patches for all sea  flux cells

  integer, allocatable :: nlfpatm(:)
  integer, allocatable :: nsfpatm(:)

  real, allocatable :: xemlfpat(:,:)
  real, allocatable :: yemlfpat(:,:)
  real, allocatable :: zemlfpat(:,:)

  real, allocatable :: xemsfpat(:,:)
  real, allocatable :: yemsfpat(:,:)
  real, allocatable :: zemsfpat(:,:)

! Land and sea flux variables

  Type flux_vars
     logical :: sendf(maxremote) = .false.
     
     integer :: ifglobe =  1
     integer :: iwrank  = -1
     integer :: iw      =  0
     integer :: kw      =  0
     integer :: iwls    =  0
     integer :: jpats   =  0  ! number of plot patches in this flux cell
     integer :: ipat    =  1  ! index of first plot patch in this flux cell

     real :: dtf     = 0.0
     real :: area    = 0.0
     real :: xef     = 0.0
     real :: yef     = 0.0
     real :: zef     = 0.0
     real :: arf_atm = 0.0
     real :: arf_sfc = 0.0
     real :: sfluxt  = 0.0
     real :: sfluxr  = 0.0
     real :: ustar   = 0.0
     real :: rhos    = 0.0
     real :: prss    = 0.0
     real :: vels    = 0.0
     real :: airtemp = 0.0
     real :: airshv  = 0.0
     real :: sxfer_t = 0.0
     real :: sxfer_r = 0.0
     real :: pcpg    = 0.0
     real :: qpcpg   = 0.0
     real :: dpcpg   = 0.0
     real :: rlong   = 0.0
     real :: rshort  = 0.0
     real :: rshort_diffuse = 0.0
  End Type flux_vars

  type (flux_vars), allocatable, target :: seaflux(:)
  type (flux_vars), allocatable, target :: landflux(:)

  type flux_pd_vars
     integer :: iw      =  0
     integer :: iwls    =  0
  end type flux_pd_vars

  type (flux_vars), allocatable, target :: seaflux_pd(:)
  type (flux_vars), allocatable, target :: landflux_pd(:)

!----------------------------------------------------------------------------

  Type seafluxg_vars
     integer :: isf_myrank = -1
     integer :: irank = -1
  End Type seafluxg_vars

  Type landfluxg_vars
     integer :: ilf_myrank = -1
     integer :: irank = -1
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
    use misc_coms,  only: iparallel, runtype
    implicit none

    if (allocated(seaflux)) then

       if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

          call increment_vtable('SEAFLUX%IFGLOBE', 'SF', noread=.true.)
          vtab_r(num_var)%ivar1_p => seaflux%ifglobe

          call increment_vtable('SEAFLUX%IWRANK',   'SF', noread=.true.)
          vtab_r(num_var)%ivar1_p => seaflux%iwrank

       endif

       call increment_vtable('SEAFLUX%SFLUXT', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sfluxt

       call increment_vtable('SEAFLUX%SFLUXR', 'SF')
       vtab_r(num_var)%rvar1_p => seaflux%sfluxr

    endif

    if (allocated(landflux)) then

       if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

          call increment_vtable('LANDFLUX%IFGLOBE', 'LF', noread=.true.)
          vtab_r(num_var)%ivar1_p => landflux%ifglobe

          call increment_vtable('LANDFLUX%IWRANK',   'LF', noread=.true.)
          vtab_r(num_var)%ivar1_p => landflux%iwrank

       endif

       call increment_vtable('LANDFLUX%SFLUXT', 'LF')
       vtab_r(num_var)%rvar1_p => landflux%sfluxt

       call increment_vtable('LANDFLUX%SFLUXR', 'LF')
       vtab_r(num_var)%rvar1_p => landflux%sfluxr

    endif

  end subroutine filltab_sflux

!============================================================================

   subroutine fill_jflux()

   use mem_ijtabs, only: itab_w, mrls, itabg_w
   use misc_coms,  only: io6, dtlm, iparallel
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

         if (iparallel == 1) then
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

         seaflux(isf)%dtf = min(dtlm(mrlw),dt_leaf)

         if (mrlw == mrl) then
      
            if (iw > 1) then
               if (iparallel == 0 .or.  &
                  (iparallel == 1 .and. itab_w(iw)%irank == myrank)) then

                  jseaflux(1)%jend(1:mrl) = jseaflux(1)%jend(1:mrl) + 1
                  jseaflux(1)%iseaflux(jseaflux(1)%jend(1)) = isf
               endif
            endif

            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_ws(iws)%irank == myrank)) then

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

         if (iparallel == 1) then
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

         landflux(ilf)%dtf = min(dtlm(mrlw),dt_leaf)

         if (mrlw == mrl) then

            if (iw > 1) then
               if (iparallel == 0 .or.  &
                  (iparallel == 1 .and. itab_w(iw)%irank == myrank)) then

                  jlandflux(1)%jend(1:mrl) = jlandflux(1)%jend(1:mrl) + 1
                  jlandflux(1)%ilandflux(jlandflux(1)%jend(1)) = ilf
               endif
            endif

            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_wl(iwl)%irank == myrank)) then

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
         if (iparallel == 1) iw = itabg_w(iw)%iw_myrank

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
         if (iparallel == 1) iw = itabg_w(iw)%iw_myrank

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
                      glatw, glonw, zm, topm, topw, arw, volt, nma, dzt
use mem_leaf,   only: land, itab_wl
use mem_sea,    only: sea,  itab_ws
use mem_ijtabs, only: itab_w, mrls
use leaf_coms,  only: nwl, dt_leaf
use sea_coms,   only: nws
use misc_coms,  only: io6, rinit

implicit none

! local variables

integer, parameter :: npmax = 7  ! Heptagons are max polygon for SINGLE GRID LEVEL
                                 ! in atm polygon cell 

integer, parameter :: nqmax = 5  ! Land cells could be up to pentagons because 
                                 ! they are generated by the piecewise-planar
                                 ! topographic surface on the OLAM triangular
                                 ! Delaunay grid and its intersection with
                                 ! model levels

integer :: iws, iwl
integer :: isf, ilf
integer :: iw, im
integer :: j, j3

integer :: ks, k, kflag, kk
integer :: mrl
integer :: np, nps, nq
integer :: nseaflux_est, nlandflux_est, nlfpats_est, nsfpats_est
integer :: incr_flux
integer :: incr_pats
integer :: jtrap, jt, it
integer :: iml, jm, jwl, jml
integer :: npoly
integer :: nspoly, ims, jms, nlpoly, jws

real :: xmean,ymean
real :: az,bz,cz,alpha
real :: xp(npmax),yp(npmax)
real :: xq(nqmax),yq(nqmax),zq(nqmax)
real :: area
real :: xtrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid x coordinates
real :: ytrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid y coordinates
real :: ztrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid z coordinates
real :: traparea(npmax+nqmax+npmax*nqmax) ! trapezoid area

integer :: itmax
integer :: iter
real :: trapareamax

real :: aatmin, aatmax, astmin, astmax, altmin, altmax
real :: dists, distn, wts, wtn

! Temporary scratch space during initialization only:

integer, allocatable :: iscrpat(:)
real,    allocatable :: scrpat(:,:)

type(flux_vars), allocatable :: landflux_temp(:)
type(flux_vars), allocatable :: seaflux_temp (:)

real :: arf_atm_test(nwa)
real :: arf_sea_test(nws)
real :: arf_land_test(nwl)

! real :: arw1t(nwa),arw2t(nwa),arw3t(nwa)

real :: xeamin(nwa),xelmin(nwl),xesmin(nws)
real :: yeamin(nwa),yelmin(nwl),yesmin(nws)
real :: zeamin(nwa),zelmin(nwl),zesmin(nws)

real :: xeamax(nwa),xelmax(nwl),xesmax(nws)
real :: yeamax(nwa),yelmax(nwl),yesmax(nws)
real :: zeamax(nwa),zelmax(nwl),zesmax(nws)

! Estimate required array sizes and allocate arrays

incr_flux = 52 * nwa
incr_pats = 55 * nwa

nseaflux_est = incr_flux
nlandflux_est = incr_flux

nlfpats_est = incr_pats
nsfpats_est = incr_pats

allocate (seaflux(nseaflux_est))
allocate (landflux(nlandflux_est))

allocate (nlfpatm (  nlfpats_est)) ; nlfpatm  = 0
allocate (xemlfpat(5,nlfpats_est)) ; xemlfpat = rinit
allocate (yemlfpat(5,nlfpats_est)) ; yemlfpat = rinit
allocate (zemlfpat(5,nlfpats_est)) ; zemlfpat = rinit

allocate (nsfpatm (  nsfpats_est)) ; nsfpatm  = 0
allocate (xemsfpat(5,nsfpats_est)) ; xemsfpat = rinit
allocate (yemsfpat(5,nsfpats_est)) ; yemsfpat = rinit
allocate (zemsfpat(5,nsfpats_est)) ; zemsfpat = rinit

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

nlfpats = 0
nsfpats = 0

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

if (mod(iw,1000) == 0) then
   print*, 'init_fluxcells ',iw,nwa
endif

! Loop over all neighbor M points of this IW

   npoly = itab_w(iw)%npoly

!---------------------------------------------------------------
! Initialize to zero the 3 planar area measures of atm polygon

!   arw1t(iw) = 0.
!   arw2t(iw) = 0.
!   arw3t(iw) = 0.

! end special
!---------------------------------------------------------------

   do jm = 1,npoly
      im = itab_w(iw)%im(jm)

! Evaluate x,y coordinates of current M point on polar stereographic plane
! tangent at IW

      call e_ps(xem(im),yem(im),zem(im),glatw(iw),glonw(iw),xp(jm),yp(jm))

!---------------------------------------------------------------
! special add up area of atm polygon based on xp,yp
!
!   if (jm > 2) then
!
!      area = .5 * abs(xp(1   ) * (yp(jm-1) - yp(jm  ))  &
!                    + xp(jm-1) * (yp(jm  ) - yp(1   ))  &
!                    + xp(jm  ) * (yp(1   ) - yp(jm-1)))
!
!      arw1t(iw) = arw1t(iw) + area
!      
!   endif
!
! end special
!---------------------------------------------------------------

   enddo

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

      nlpoly = itab_wl(iwl)%npoly

! Loop over all neighbor M points of this IWL

      do jml = 1,nlpoly
         iml = itab_wl(iwl)%im(jml)

! Evaluate x,y coordinates of LAND cell M points on polar stereographic plane
! tangent at IW

         call e_ps(land%xem(iml),land%yem(iml),land%zem(iml), &
                   glatw(iw),glonw(iw),xq(jml),yq(jml))

! Store topography height of M point in zq array

         zq(jml) = land%zm(iml)

      enddo

! Evaluate possible overlap of ATM and LAND polygons

      call polygon_overlap(iwl,npoly,nlpoly,xp,yp,xq,yq,area,  &
                           jtrap,xtrap,ytrap,traparea)

!-------------------------------------------------------
! special - add up 2nd measure of ATM cell area
!      arw2t(iw) = arw2t(iw) + area
! end special
!-------------------------------------------------------

! Skip iwl if overlap is zero

      if (area < 1.e-7) cycle

! This IW polygon overlaps with land cell IWL.

! Find fit coefficients for linear elevation surface of land cell IWL
! using 3 vertices of IWL polygon. (It is assumed that IWL polygon is convex, 
! so that selected points are not collinear.)

      j3 = 3
      if (nlpoly > 4) j3 = 4

! Do this set if points 1 and 2 do not coincide

!            if (abs(xtrap(1,jt) - xtrap(2,jt)) > 1.e0) then

      call matrix_3x3(1.,xq(1),yq(1)      &
                     ,1.,xq(2),yq(2)      &
                     ,1.,xq(j3),yq(j3)    &
                     ,zq(1),zq(2),zq(j3)  &
                     ,az,bz,cz            )

! Loop over all trapezoids of current IW-IWL overlap

      do jt = 1,jtrap

! Find height of 4 corners of trapezoid

         ztrap(1:4,jt) = az + bz * xtrap(1:4,jt) + cz * ytrap(1:4,jt)


!-------------------------------------------------------
! special - add up 3rd measure of ATM cell area
!      arw3t(iw) = arw3t(iw) + traparea(jt)
! end special
!-------------------------------------------------------
      
      enddo

! Determine TOPW by checking whether IW point is inside (or within 10 m of) 
! current IWL polygon

      call inout_check(nlpoly,xq,yq,0.,0.,alpha)
      if (alpha > 1.) topw(iw) = az
      call inout_check(nlpoly,xq,yq,-10.,0.,alpha)
      if (alpha > 1.) topw(iw) = az
      call inout_check(nlpoly,xq,yq,10.,0.,alpha)
      if (alpha > 1.) topw(iw) = az
      call inout_check(nlpoly,xq,yq,0.,-10.,alpha)
      if (alpha > 1.) topw(iw) = az
      call inout_check(nlpoly,xq,yq,0.,10.,alpha)
      if (alpha > 1.) topw(iw) = az

! Determine TOPM values of IW cell by checking whether each IM point is inside
! (or within 10 m of) current IWL polygon

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         call inout_check(nlpoly,xq,yq,xp(jm)-10.,yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         call inout_check(nlpoly,xq,yq,xp(jm)+10.,yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm)-10.,alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         call inout_check(nlpoly,xq,yq,xp(jm),yp(jm)+10.,alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         
      enddo

! Skip remaining calculations if overlap is small

      if (area < 1.e-6 * min(arw0(iw),land%area(iwl))) cycle

! Loop over vertical T levels in IW

      do k = 2,nza-1

! Set kflag to zero to indicate that flux cell for (IW,IWL,K) is not yet activated.

         kflag = 0

         xmean = 0.
         ymean = 0.

! Loop over all trapezoids of current IW-IWL overlap

         do jt = 1,jtrap

! Check whether all trapezoid vertex elevations are at or above top of T level

            if (all(ztrap(1:4,jt) >= zm(k))) cycle

! Check whether all trapezoid vertex elevations are below bottom of T level

            if (all(ztrap(1:4,jt) < zm(k-1))) then

               arw (k,iw) = arw (k,iw) + traparea(jt)
               volt(k,iw) = volt(k,iw) + traparea(jt) * dzt(k)

               cycle
            endif

! If we got here, this T level overlaps surface height range of this trapezoid.
! If area is above a small threshold, a flux cell will be activated for the
! current set of (IW,IWL,K).


! Find portion of current trapezoid that overlaps height range of current
! K level, and produce reduced-area polygon that represents this area.

!-----------------------------------------
! First set of 3 trapezoid vertices
!-----------------------------------------

! Skip this set if all 3 vertex elevations are at or above top of T level

            if (ztrap(1,jt) >= zm(k) .and.  &
                ztrap(2,jt) >= zm(k) .and.  &
                ztrap(3,jt) >= zm(k)) go to 11

! Skip this set if all 3 vertex elevations are below bottom of T level

            if (ztrap(1,jt) < zm(k-1) .and.  &
                ztrap(2,jt) < zm(k-1) .and.  &
                ztrap(3,jt) < zm(k-1)) go to 11

! Do this set if points 1 and 2 do not coincide

            if (abs(xtrap(1,jt) - xtrap(2,jt)) > 1.e0) then

               call trap3(iw,k,kflag,npoly,npmax,   &
                          xmean,ymean,xp,yp,              &
                          xtrap(1,jt),xtrap(2,jt),xtrap(3,jt),  &
                          ytrap(1,jt),ytrap(2,jt),ytrap(3,jt),  &
                          ztrap(1,jt),ztrap(2,jt),ztrap(3,jt),  &
                          iwl=iwl,ilf=ilf)

            endif

            11 continue

!-----------------------------------------
! Second set of 3 trapezoid vertices
!-----------------------------------------

! Skip this set if all 3 vertex elevations are at or above top of T level

            if (ztrap(1,jt) >= zm(k) .and.  &
                ztrap(3,jt) >= zm(k) .and.  &
                ztrap(4,jt) >= zm(k)) go to 12

! Skip this set if all 3 vertex elevations are below bottom of T level

            if (ztrap(1,jt) < zm(k-1) .and.  &
                ztrap(3,jt) < zm(k-1) .and.  &
                ztrap(4,jt) < zm(k-1)) go to 12

! Do this set if points 3 and 4 do not coincide

            if (abs(xtrap(3,jt) - xtrap(4,jt)) > 1.e0) then

               call trap3(iw,k,kflag,npoly,npmax,   &
                          xmean,ymean,xp,yp,              &
                          xtrap(3,jt),xtrap(4,jt),xtrap(1,jt),  &
                          ytrap(3,jt),ytrap(4,jt),ytrap(1,jt),  &
                          ztrap(3,jt),ztrap(4,jt),ztrap(1,jt),  &
                          iwl=iwl,ilf=ilf)

            endif

            12 continue

         enddo  ! jt

! Convert (x,y) of center of mass to earth coordinates

         if (kflag == 0) cycle

         call ps_e(landflux(ilf)%xef,         &
                   landflux(ilf)%yef,         &
                   landflux(ilf)%zef,         &
                   glatw(iw),glonw(iw),       &
                   xmean/landflux(ilf)%area,  &
                   ymean/landflux(ilf)%area   )

         landflux(ilf)%arf_atm  = landflux(ilf)%area / arw0(iw)
         landflux(ilf)%arf_sfc = landflux(ilf)%area / land%area(iwl)

      enddo ! k loop

   enddo ! jwl loop

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

      nspoly = itab_ws(iws)%npoly

! Loop over all neighbor M points of this IWS

      do jms = 1,nspoly
         ims = itab_ws(iws)%im(jms)

! Evaluate x,y coordinates of SEA cell M points on polar stereographic plane
! tangent at IW

         call e_ps(sea%xem(ims), sea%yem(ims), sea%zem(ims),  &
                   glatw(iw), glonw(iw), xq(jms), yq(jms))

! Store topography height of M point in zq array

         zq(jms) = sea%zm(ims)

      enddo

! Evaluate possible overlap of ATM and SEA polygons

      call polygon_overlap(iws,npoly,nspoly,xp,yp,xq,yq,area,  &
                           jtrap,xtrap,ytrap,traparea)


!-------------------------------------------------------
! special - add up 2nd measure of ATM cell area
!      arw2t(iw) = arw2t(iw) + area
! end special
!-------------------------------------------------------
      
! Skip iwl if overlap is zero

      if (area < 1.e-7) cycle

! This IW polygon overlaps with sea cell IWS.

! Find fit coefficients for linear elevation surface of sea cell IWS
! using 3 vertices of IWS polygon. (It is assumed that IWS polygon is convex, 
! so that selected points are not collinear.)

      j3 = 3
      if (nspoly > 4) j3 = 4

      call matrix_3x3(1.,xq(1),yq(1)      &
                     ,1.,xq(2),yq(2)      &
                     ,1.,xq(j3),yq(j3)    &
                     ,zq(1),zq(2),zq(j3)  &
                     ,az,bz,cz            )

! Loop over all trapezoids of current IW-IWS overlap

      do jt = 1,jtrap

! Find height of 4 corners of trapezoid

         ztrap(1:4,jt) = az + bz * xtrap(1:4,jt) + cz * ytrap(1:4,jt)

!-------------------------------------------------------
! special - add up 3rd measure of ATM cell area
!      arw3t(iw) = arw3t(iw) + traparea(jt)
! end special
!-------------------------------------------------------
      
      enddo

! Determine TOPW by checking whether IW point is inside (or within 10 of)
! current IWS polygon

      call inout_check(nspoly,xq,yq,0.,0.,alpha)
      if (alpha > 1.) topw(iw) = az

      call inout_check(nspoly,xq,yq,-10.,0.,alpha)
      if (alpha > 1.) topw(iw) = az

      call inout_check(nspoly,xq,yq,10.,0.,alpha)
      if (alpha > 1.) topw(iw) = az

      call inout_check(nspoly,xq,yq,0.,-10.,alpha)
      if (alpha > 1.) topw(iw) = az

      call inout_check(nspoly,xq,yq,0.,10.,alpha)
      if (alpha > 1.) topw(iw) = az

! Determine TOPM values of IW cell by checking whether each IM point is inside
! (or within 10 m of) current IWS polygon

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)

         call inout_check(nspoly,xq,yq,xp(jm)-10.,yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)

         call inout_check(nspoly,xq,yq,xp(jm)+10.,yp(jm),alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm)-10.,alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)

         call inout_check(nspoly,xq,yq,xp(jm),yp(jm)+10.,alpha)
         if (alpha > 1.) topm(im) = az + bz * xp(jm) + cz * yp(jm)
         
      enddo

! Skip remaining calculations if overlap is small

      if (area < 1.e-6 * min(arw0(iw),sea%area(iws))) cycle

! Loop over vertical T levels in IW

      do k = 2,nza-1

! Set kflag to zero to indicate that flux cell for (IW,IWS,K) is not yet activated.

         kflag = 0

         xmean = 0.
         ymean = 0.

! Loop over all trapezoids of current IW-IWS overlap

         do jt = 1,jtrap

! Check whether all trapezoid vertex elevations are at or above top of T level

            if (all(ztrap(1:4,jt) >= zm(k))) cycle

! Check whether all trapezoid vertex elevations are below bottom of T level

            if (all(ztrap(1:4,jt) < zm(k-1))) then

               arw (k,iw) = arw (k,iw) + traparea(jt)
               volt(k,iw) = volt(k,iw) + traparea(jt) * dzt(k)

               cycle
            endif

! If we got here, this T level overlaps surface height range of this trapezoid,
! and a flux cell will be activated for the current set of (IW,IWS,K).

! Find portion of current trapezoid that overlaps height range of current
! K level, and produce reduced-area polygon that represents this area.

!-----------------------------------------
! First set of 3 trapezoid vertices
!-----------------------------------------

! Skip this set if all 3 vertex elevations are at or above top of T level

            if (ztrap(1,jt) >= zm(k) .and.  &
                ztrap(2,jt) >= zm(k) .and.  &
                ztrap(3,jt) >= zm(k)) go to 13

! Skip this set if all 3 vertex elevations are below bottom of T level

            if (ztrap(1,jt) < zm(k-1) .and.  &
                ztrap(2,jt) < zm(k-1) .and.  &
                ztrap(3,jt) < zm(k-1)) go to 13

! Do this set if points 1 and 2 do not coincide

            if (abs(xtrap(1,jt) - xtrap(2,jt)) > 1.e0) then

               call trap3(iw,k,kflag,npoly,npmax,   &
                          xmean,ymean,xp,yp,              &
                          xtrap(1,jt),xtrap(2,jt),xtrap(3,jt),  &
                          ytrap(1,jt),ytrap(2,jt),ytrap(3,jt),  &
                          ztrap(1,jt),ztrap(2,jt),ztrap(3,jt),  &
                          iws=iws,isf=isf)

            endif

            13 continue

!-----------------------------------------
! Second set of 3 trapezoid vertices
!-----------------------------------------

! Skip this set if all 3 vertex elevations are at or above top of T level

            if (ztrap(1,jt) >= zm(k) .and.  &
                ztrap(3,jt) >= zm(k) .and.  &
                ztrap(4,jt) >= zm(k)) go to 14

! Skip this set if all 3 vertex elevations are below bottom of T level

            if (ztrap(1,jt) < zm(k-1) .and.  &
                ztrap(3,jt) < zm(k-1) .and.  &
                ztrap(4,jt) < zm(k-1)) go to 14

! Do this set if points 3 and 4 do not coincide

            if (abs(xtrap(3,jt) - xtrap(4,jt)) > 1.e0) then

               call trap3(iw,k,kflag,npoly,npmax,   &
                          xmean,ymean,xp,yp,              &
                          xtrap(3,jt),xtrap(4,jt),xtrap(1,jt),  &
                          ytrap(3,jt),ytrap(4,jt),ytrap(1,jt),  &
                          ztrap(3,jt),ztrap(4,jt),ztrap(1,jt),  &
                          iws=iws,isf=isf)

            endif

            14 continue

         enddo  ! jt

         if (kflag == 0) cycle

! Convert (x,y) of center of mass to earth coordinates

         call ps_e(seaflux(isf)%xef,         &
                   seaflux(isf)%yef,         &
                   seaflux(isf)%zef,         &
                   glatw(iw),glonw(iw),      &
                   xmean/seaflux(isf)%area,  &
                   ymean/seaflux(isf)%area   )

         seaflux(isf)%arf_atm = seaflux(isf)%area / arw0(iw)
         seaflux(isf)%arf_sfc = seaflux(isf)%area / sea%area(iws)

      enddo ! k loop

   enddo ! jws loop

! If number of flux cells or polygons is getting close to allocated 
! quantity, allocate more space

   if (nseaflux_est - isf < 1000) then 

      nseaflux_est = nseaflux_est + incr_flux
      allocate (seaflux_temp(nseaflux_est))
      seaflux_temp(1:isf) = seaflux(1:isf)
      call move_alloc(seaflux_temp, seaflux)

   endif

   if (nlandflux_est - ilf < 1000) then 

      nlandflux_est = nlandflux_est + incr_flux
      allocate (landflux_temp(nlandflux_est))
      landflux_temp(1:ilf) = landflux(1:ilf)
      call move_alloc(landflux_temp, landflux)

   endif

   if (nsfpats_est - nsfpats < 1000) then 

      nsfpats_est = nsfpats_est + incr_pats

      allocate (iscrpat(nsfpats_est))
      iscrpat(1:nsfpats) = nsfpatm(1:nsfpats)
      call move_alloc(iscrpat, nsfpatm)

      allocate (scrpat(5,nsfpats))

      scrpat(1:5,1:nsfpats) = xemsfpat(1:5,1:nsfpats)
      deallocate (xemsfpat)
      allocate (xemsfpat(5,nsfpats_est))
      xemsfpat(1:5,1:nsfpats) = scrpat(1:5,1:nsfpats)
      
      scrpat(1:5,1:nsfpats) = yemsfpat(1:5,1:nsfpats)
      deallocate (yemsfpat)
      allocate (yemsfpat(5,nsfpats_est))
      yemsfpat(1:5,1:nsfpats) = scrpat(1:5,1:nsfpats)
      
      scrpat(1:5,1:nsfpats) = zemsfpat(1:5,1:nsfpats)
      deallocate (zemsfpat)
      allocate (zemsfpat(5,nsfpats_est))
      zemsfpat(1:5,1:nsfpats) = scrpat(1:5,1:nsfpats)
      
      deallocate (scrpat)
      
   endif      

   if (nlfpats_est - nlfpats < 1000) then 

      nlfpats_est = nlfpats_est + incr_pats

      allocate (iscrpat(nlfpats_est))
      iscrpat(1:nlfpats) = nlfpatm(1:nlfpats)
      call move_alloc(iscrpat, nlfpatm)
      
      allocate (scrpat(5,nlfpats))

      scrpat(1:5,1:nlfpats) = xemlfpat(1:5,1:nlfpats)
      deallocate (xemlfpat)
      allocate (xemlfpat(5,nlfpats_est))
      xemlfpat(1:5,1:nlfpats) = scrpat(1:5,1:nlfpats)
      
      scrpat(1:5,1:nlfpats) = yemlfpat(1:5,1:nlfpats)
      deallocate (yemlfpat)
      allocate (yemlfpat(5,nlfpats_est))
      yemlfpat(1:5,1:nlfpats) = scrpat(1:5,1:nlfpats)
      
      scrpat(1:5,1:nlfpats) = zemlfpat(1:5,1:nlfpats)
      deallocate (zemlfpat)
      allocate (zemlfpat(5,nlfpats_est))
      zemlfpat(1:5,1:nlfpats) = scrpat(1:5,1:nlfpats)
      
      deallocate (scrpat)
      
   endif
      
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

!-----------------------------------------------------------------------
! special - compare areas
!do iw = 2,nwa
!
!   write(6,'(a,i6,4e13.4,10f10.5)') 'arw123 compare ',iw,arw0(iw),arw1t(iw),  &
!       arw2t(iw),arw3t(iw),arw1t(iw)/arw0(iw),  &
!       arw2t(iw)/arw1t(iw),arw3t(iw)/arw1t(iw)
!enddo
!
! end special
!-----------------------------------------------------------------------

! Make sure that all topm and topw values have been assigned

do im = 2,nma
   if (topm(im) < -1.e3) then
      print*, 'In init_fluxcells, topm missing at im = ',im
      stop 'stop_init_fluxcells_topm'
   endif
enddo

do iw = 2,nwa
   if (topw(iw) < -1.e3) then
      print*, 'In init_fluxcells, topw missing at iw = ',iw
      stop 'stop_init_fluxcells_topw'
   endif
enddo

return
end subroutine init_fluxcells

!============================================================================

subroutine trap3(iw,k,kflag,npoly,npmax,   &
                 xmean,ymean,xp,yp,  &
                 xtrap1,xtrap2,xtrap3,     &
                 ytrap1,ytrap2,ytrap3,     &
                 ztrap1,ztrap2,ztrap3,     &
                 iwl,iws,ilf,isf)

use mem_grid,   only: nza,nwa,zm,glatw,glonw,arw,volt
use mem_ijtabs, only: itab_w

implicit none

integer, intent(in) :: iw,k,npoly,npmax
integer, intent(inout) :: kflag

integer, optional, intent(in) :: iwl,iws
integer, optional, intent(inout) :: ilf,isf

real, intent(inout) :: xmean,ymean
real, intent(in) :: xp(npmax),yp(npmax)

real, intent(in) :: xtrap1,xtrap2,xtrap3
real, intent(in) :: ytrap1,ytrap2,ytrap3
real, intent(in) :: ztrap1,ztrap2,ztrap3

integer, parameter :: nsub = 30

integer :: npatm,jpatm
integer :: ipair,isub,jsub

real :: xmpat(5),ympat(5),zmpat(5)

real :: x1,x2,x3,x4,y1,y2,y3,y4,xp0,yp0
real :: del_arw
real :: dt13,dt32,dt14,dt43,dx13,dx32,dx14,dx43,dy13,dy32,dy14,dy43
real :: faci1,faci2,facj1,facj2
real :: top1,top2,top3,top4,topp
real :: flux_area

! Compute intersection of current surface triangle (half-trapezoid) with 
! current atmospheric model level.  This results in one polygonal patch 
! with 3 to 5 vertices.

call sfc_patch(zm(k-1),zm(k),xtrap1,xtrap2,xtrap3,  &
                             ytrap1,ytrap2,ytrap3,  &
                             ztrap1,ztrap2,ztrap3,  &
                             npatm,xmpat,ympat,zmpat) ! zmpat not currently used

! Interpolate 3 topography heights to nsub*nsub sub-triangles

dt13 = ztrap3 - ztrap1
dt32 = ztrap2 - ztrap3

dx13 = xtrap3 - xtrap1
dx32 = xtrap2 - xtrap3

dy13 = ytrap3 - ytrap1
dy32 = ytrap2 - ytrap3

! Planar area of triangle (1/2 times base times height)

del_arw = abs(.5 * (xtrap2 - xtrap1) * (ytrap3 - ytrap1)) / real(nsub*nsub)

! [NEED TO REPLACE ABOVE AREA WITH SPHERICAL AREA ?]

flux_area = 0.

do jsub = 1,nsub

   facj1 = real(jsub-1) / real(nsub)
   facj2 = real(jsub)   / real(nsub)

   do isub = 1,jsub

      faci1 = real(isub-1) / real(nsub)
      faci2 = real(isub)   / real(nsub)

      top1 = ztrap1 + facj1 * dt13 + faci1 * dt32
      top2 = ztrap1 + facj1 * dt13 + faci2 * dt32
      top3 = ztrap1 + facj2 * dt13 + faci1 * dt32
      top4 = ztrap1 + facj2 * dt13 + faci2 * dt32

      x1 = xtrap1 + facj1 * dx13 + faci1 * dx32
      x2 = xtrap1 + facj1 * dx13 + faci2 * dx32
      x3 = xtrap1 + facj2 * dx13 + faci1 * dx32
      x4 = xtrap1 + facj2 * dx13 + faci2 * dx32

      y1 = ytrap1 + facj1 * dy13 + faci1 * dy32
      y2 = ytrap1 + facj1 * dy13 + faci2 * dy32
      y3 = ytrap1 + facj2 * dy13 + faci1 * dy32
      y4 = ytrap1 + facj2 * dy13 + faci2 * dy32

      do ipair = 1,2

         if (ipair == 1) then
            topp = (top1 + top3 + top4) / 3.
            xp0 = (x1 + x3 + x4) / 3.
            yp0 = (y1 + y3 + y4) / 3.
         else
            topp = (top1 + top2 + top4) / 3.
            xp0 = (x1 + x2 + x4) / 3.
            yp0 = (y1 + y2 + y4) / 3.
         endif

         if (topp < zm(k)) then
            arw (k,iw) = arw (k,iw) + del_arw
            volt(k,iw) = volt(k,iw) + del_arw * (zm(k) - max(topp,zm(k-1)))

            if (topp >= zm(k-1)) then
               flux_area = flux_area + del_arw
               
               xmean = xmean + xp0 * del_arw
               ymean = ymean + yp0 * del_arw
            endif
         endif

         if (isub == jsub) exit 

      enddo  ! ipair

   enddo  ! ip

enddo  ! jp

! Check current contribution to flux area

if (flux_area > 1.e2) then

! Convert xpat,ypat from ps to earth coordinates

   if (present(iwl)) then

! If kflag = 0, this is the first plot patch for this flux cell.  Activate flux cell.

      if (kflag == 0) then
         kflag = 1

         ilf = ilf + 1

         landflux(ilf)%ifglobe = ilf ! full domain index
         landflux(ilf)%iw      = iw   ! full domain index
         landflux(ilf)%kw      = k
         landflux(ilf)%iwls    = iwl ! full domain index
         landflux(ilf)%ipat    = nlfpats + 1
         landflux(ilf)%jpats   = 0
      endif

      landflux(ilf)%jpats = landflux(ilf)%jpats + 1
      landflux(ilf)%area  = landflux(ilf)%area + flux_area

      nlfpats = nlfpats + 1

      nlfpatm(nlfpats) = npatm

      do jpatm = 1,npatm
         call ps_e(xemlfpat(jpatm,nlfpats),  &
                   yemlfpat(jpatm,nlfpats),  &
                   zemlfpat(jpatm,nlfpats),  &
                   glatw(iw),glonw(iw),      &
                   xmpat(jpatm),ympat(jpatm) )
      enddo

   elseif (present(iws)) then

! If kflag = 0, this is the first plot patch for this flux cell.  Activate flux cell.

      if (kflag == 0) then
         kflag = 1

         isf = isf + 1

         seaflux(isf)%ifglobe = isf ! full domain index
         seaflux(isf)%iw      = iw  ! full domain index
         seaflux(isf)%kw      = k
         seaflux(isf)%iwls    = iws ! full domain index
         seaflux(isf)%ipat    = nsfpats + 1
         seaflux(isf)%jpats   = 0
      endif

      seaflux(isf)%jpats = seaflux(isf)%jpats + 1
      seaflux(isf)%area  = seaflux(isf)%area + flux_area

      nsfpats = nsfpats + 1

      nsfpatm(nsfpats) = npatm
      
      do jpatm = 1,npatm
         call ps_e(xemsfpat(jpatm,nsfpats),  &
                   yemsfpat(jpatm,nsfpats),  &
                   zemsfpat(jpatm,nsfpats),  &
                   glatw(iw),glonw(iw),      &
                   xmpat(jpatm),ympat(jpatm) )

      enddo

   endif

endif

return
end subroutine trap3

!============================================================================

subroutine sfc_patch(zm1,zm2,x1,x2,x3,y1,y2,y3,z1,z2,z3,np,xp,yp,zp)

implicit none

real, intent(in) :: zm1,zm2  ! Bottom, top of current model T level
real, intent(in) :: z1,z2,z3 ! Topography heights at 3 M pts
real, intent(in) :: x1,x2,x3 ! x ps coords of 3 M pts
real, intent(in) :: y1,y2,y3 ! y ps coords of 3 M pts

integer, intent(out) :: np ! Number of patch pts
real, intent(out) :: xp(5),yp(5),zp(5) ! x,y,z ps coords of patch pts
  
if     (z1 <= z2 .and. z1 <= z3) then
   call sfc_pat(zm1,zm2,z1,z2,z3,x1,x2,x3,y1,y2,y3,np,xp,yp,zp)
elseif (z2 <= z3 .and. z2 <= z1) then
   call sfc_pat(zm1,zm2,z2,z3,z1,x2,x3,x1,y2,y3,y1,np,xp,yp,zp)
else
   call sfc_pat(zm1,zm2,z3,z1,z2,x3,x1,x2,y3,y1,y2,np,xp,yp,zp)
endif

return
end subroutine sfc_patch

!============================================================================

subroutine sfc_pat(zm1,zm2,z1,z2,z3,x1,x2,x3,y1,y2,y3,np,xp,yp,zp)

use misc_coms,  only: io6

implicit none

real, intent(in) :: zm1,zm2  ! Bottom, top of current model T level
real, intent(in) :: z1,z2,z3 ! Topo heights at 3 M pts, cyclic order, z1 lowest
real, intent(in) :: x1,x2,x3 ! x ps coords of 3 M pts
real, intent(in) :: y1,y2,y3 ! y ps coords of 3 M pts

integer, intent(out) :: np ! Number of patch pts
real, intent(out) :: xp(5),yp(5),zp(5) ! x,y,z ps coords of patch pts

if (zm1 <= z1) then

   xp(1) = x1
   yp(1) = y1
   zp(1) = z1

   if (zm2 <= z2) then

      xp(2) = x1 + (x2-x1) * (zm2-z1) / (z2-z1)
      yp(2) = y1 + (y2-y1) * (zm2-z1) / (z2-z1)
      zp(2) = zm2

      if (zm2 <= z3) then

         xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)
         zp(3) = zm2

         np = 3

      else
      
         xp(3) = x3 + (x2-x3) * (zm2-z3) / (z2-z3)
         yp(3) = y3 + (y2-y3) * (zm2-z3) / (z2-z3)
         zp(3) = zm2

         xp(4) = x3
         yp(4) = y3
         zp(4) = z3

         np = 4
         
      endif         
      
   else
      
      xp(2) = x2
      yp(2) = y2
      zp(2) = z2

      if (zm2 <= z3) then

         xp(3) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
         yp(3) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)
         zp(3) = zm2

         xp(4) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)
         zp(4) = zm2

         np = 4

      else
      
         xp(3) = x3
         yp(3) = y3
         zp(3) = z3

         np = 3
         
      endif         
      
   endif

elseif (zm1 <= z2) then

   xp(1) = x1 + (x2-x1) * (zm1-z1) / (z2-z1)
   yp(1) = y1 + (y2-y1) * (zm1-z1) / (z2-z1)
   zp(1) = zm1

   if (zm2 <= z2) then

      xp(2) = x1 + (x2-x1) * (zm2-z1) / (z2-z1)
      yp(2) = y1 + (y2-y1) * (zm2-z1) / (z2-z1)
      zp(2) = zm2

      if (zm2 <= z3) then

         xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)
         zp(3) = zm2

         xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
         zp(4) = zm1

         np = 4
            
      else

         xp(3) = x3 + (x2-x3) * (zm2-z3) / (z2-z3)
         yp(3) = y3 + (y2-y3) * (zm2-z3) / (z2-z3)
         zp(3) = zm2

         if (zm1 <= z3) then

            xp(4) = x3
            yp(4) = y3
            zp(4) = z3

            xp(5) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
            yp(5) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
            zp(5) = zm1

            np = 5
               
         else
            
            xp(4) = x3 + (x2-x3) * (zm1-z3) / (z2-z3)
            yp(4) = y3 + (y2-y3) * (zm1-z3) / (z2-z3)
            zp(4) = zm1

            np = 4
            
         endif
            
      endif
          
   else

      xp(2) = x2
      yp(2) = y2
      zp(2) = z2

      if (zm2 <= z3) then

         xp(3) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
         yp(3) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)
         zp(3) = zm2

         xp(4) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)
         zp(4) = zm2

         xp(5) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
         yp(5) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
         zp(5) = zm1

         np = 5
               
      else
         
         if (zm1 <= z3) then

            xp(3) = x3
            yp(3) = y3
            zp(3) = z3

            xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
            yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
            zp(4) = zm1

            np = 4
               
         else
            
            xp(3) = x3 + (x2-x3) * (zm1-z3) / (z2-z3)
            yp(3) = y3 + (y2-y3) * (zm1-z3) / (z2-z3)
            zp(3) = zm1

            np = 3
            
         endif

      endif
      
   endif

else

   xp(1) = x2 + (x3-x2) * (zm1-z2) / (z3-z2)
   yp(1) = y2 + (y3-y2) * (zm1-z2) / (z3-z2)
   zp(1) = zm1

   if (zm2 <= z3) then

      xp(2) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
      yp(2) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)
      zp(2) = zm2

      xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
      yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)
      zp(3) = zm2

      xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
      yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
      zp(4) = zm1

      np = 4
               
   else   

      xp(2) = x3
      yp(2) = y3
      zp(2) = z3

      xp(3) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
      yp(3) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)
      zp(3) = zm1

      np = 3

   endif
            
endif

return
end subroutine sfc_pat

!============================================================================

subroutine polygon_overlap(iwls,np,nq,xp,yp,xq,yq,area,ntrap,xtrap,ytrap,traparea)

! Given x,y coordinates of the vertices of polygons p and q, compute the area
! of overlap between the polygons using a sweepline algorithm.

! Method adapted from:
! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

implicit none

integer, intent(in) :: np,nq,iwls ! Number of vertices in p and q
real, intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
real, intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

integer, intent(out) :: ntrap              ! number of trapezoids
real, intent(out) :: xtrap(4,np+nq+np*nq)  ! trapezoid x coordinates
real, intent(out) :: ytrap(4,np+nq+np*nq)  ! trapezoid y coordinates
real, intent(out) :: traparea(np+nq+np*nq) ! trapezoid area
real, intent(out) :: area                  ! area of overlap of p and q

real :: yev(np+nq+np*nq)  ! y-coordinates of event

integer :: nev      ! # of events
integer :: nsect    ! # of intersections between strip centerline and p,q edges
real :: ymid        ! y-coord of centerline of strip between consecutive events
real :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges 
real :: x1(np+nq)   ! x-coords where strip .1 line intersects p and q edges 
real :: xcent       ! x-coord of midpoint of segment between xmid values

integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,i,ia,ib,is
real :: p0,q0,dx,dy,y1,dxtrap

real :: alpha    ! angle swept by ray from point to polygon perimeter during circuit

! Find vertices in p that are not outside q

nev = 0
area = 0.
ntrap = 0

do ip = 1,np
   call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

   if (abs(alpha) > .2) then 
      nev = nev + 1
      yev(nev) = yp(ip)
   endif
enddo

! Find vertices in q that are not outside p

do iq = 1,nq
   call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

   if (abs(alpha) > .2) then 
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
      if (p0 < -0.000001 .or. p0 > 1.000001) cycle
      if (q0 < -0.000001 .or. q0 > 1.000001) cycle

! Line segments pa-pb and qa-qb intersect; find y-coordinate of intersection

      nev = nev + 1
      yev(nev) = (1. - p0) * yp(ipa) + p0 * yp(ipb)

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
   
! Reject overlap event if dy is less than 1 meter (threshold ok for single precision)
   
   if (dy < 1.) cycle

! ymid = y-coordinate of strip centerline.
! Initialize dx (centerline length sum) to 0.

   ymid = yev(iev) + .5 * dy
   y1   = yev(iev) + .1 * dy
   dx = 0.

! Find x-coordinate of intersections of strip centerline with edges of p and q

   nsect = 0

   do ia = 1,np
      ib = ia + 1
      if (ia == np) ib = 1

      if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xp(ia)  &
                  + (xp(ib) - xp(ia)) * (ymid - yp(ia)) / (yp(ib) - yp(ia))
      x1(nsect) = xp(ia)  &
                + (xp(ib) - xp(ia)) * (y1 - yp(ia)) / (yp(ib) - yp(ia))
   enddo

   do ia = 1,nq
      ib = ia + 1
      if (ia == nq) ib = 1

      if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xq(ia)  &
                  + (xq(ib) - xq(ia)) * (ymid - yq(ia)) / (yq(ib) - yq(ia))
      x1(nsect) = xq(ia)  &
                + (xq(ib) - xq(ia)) * (y1 - yq(ia)) / (yq(ib) - yq(ia))
   enddo

! Sort xmid values into increasing order

   call fltsort2(nsect,xmid,x1)

! see if the segment is inside both polygons

   do is = 1,nsect - 1
      xcent = .5 * (xmid(is) + xmid(is+1))

      if (xcent == xmid(is)) cycle          ! if segment length is 0

      call inout_check(np,xp,yp,xcent,ymid,alpha)

      if (abs(alpha) < 1.) cycle

      call inout_check(nq,xq,yq,xcent,ymid,alpha)

      if (abs(alpha) < 1.) cycle

      dxtrap = xmid(is+1) - xmid(is)

      dx = dx + dxtrap

! Find x at 4 corners of current trapezoidal region

      ntrap = ntrap + 1
      
      xtrap(1,ntrap) =  1.25 * x1(is)   -  .25 * xmid(is) 
      xtrap(2,ntrap) =  1.25 * x1(is+1) -  .25 * xmid(is+1) 
      xtrap(3,ntrap) = -1.25 * x1(is+1) + 2.25 * xmid(is+1) 
      xtrap(4,ntrap) = -1.25 * x1(is)   + 2.25 * xmid(is)
      
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

implicit none

integer, intent(in) :: n  ! number of points in polygon

real, intent(in) :: x(n),y(n)  ! x,y coordinates of polygon points      
real, intent(in) :: x0,y0      ! x,y coordinates of additional point
real, intent(out) :: th1       ! output value

! Local variables

integer :: i
real :: theta
real :: xh1, xh2
real :: edgemax
real :: x1(n+1), y1(n+1)

real, parameter :: pi1 = 3.1415926536, pi2 = 2. * pi1

do i = 1,n
   x1(i) = x(i) - x0
   y1(i) = y(i) - y0

! Check whether x0,y0 is exactly on a vertex (to within 1.e-3 m)

   if ((abs(x1(i)) < .1) .and. (abs(y1(i)) < 1.e-3)) then
      th1 = pi1
      return
   endif
enddo

x1(n+1) = x1(1)
y1(n+1) = y1(1)
th1 = 0.

xh2 = atan2(y1(1),x1(1))
if (xh2 < 0.) xh2 = xh2 + pi2

do i = 1,n
   xh1 = atan2(y1(i+1),x1(i+1))
   if (xh1 < 0.) xh1 = xh1 + pi2

   theta = xh1 - xh2
   if (theta < -pi1) theta = theta + pi2

   if ((abs(abs(theta) - pi1) < .00001)) then
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

implicit none

real, intent(in) :: xpa,ypa,xpb,ypb  ! x,y coordinates of pa and pb
real, intent(in) :: xqa,yqa,xqb,yqb  ! x,y coordinates of qa and qb

integer, intent(out) :: iflag
real, intent(out) :: p0,q0

real :: pabxqab  ! cross product of vectors pa-pb and qa-qb 

iflag = -1

pabxqab = (xpb - xpa) * (yqb - yqa) - (ypb - ypa) * (xqb - xqa)

if (pabxqab /= 0.) then  ! Infinite lines intersect
   iflag = 0

   p0 = ((yqb - yqa) * (xqa - xpa) + (yqa - ypa) * (xqa - xqb)) / pabxqab
   q0 = ((ypb - ypa) * (xqa - xpa) + (yqa - ypa) * (xpa - xpb)) / pabxqab
endif

!print*, 'pabxqab ',pabxqab,p0,q0,iflag

return
end subroutine intersect

!============================================================================

subroutine fltsort(n,f)

! Sort n floating point numbers f into ascending order

implicit none

integer, intent(in) :: n
real, intent(inout) :: f(n)

integer :: i,j

real :: f0

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

implicit none

integer, intent(in) :: n
real, intent(inout) :: f1(n),f2(n)

integer :: i,j

real :: f0

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

End Module mem_sflux
