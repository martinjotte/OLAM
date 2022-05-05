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

Module mem_ijtabs

  use max_dims, only: maxremote

  implicit none

  private :: maxremote

  integer, parameter :: mloops = 7 ! max # non-para DO loops for M,V,W pts

  ! M, U, V, and W loop indices. The last 4 letters of these indices have the
  ! following meanings:
  !
  ! *_grid is for setting up grid parameters in the ctrlvols subroutines
  ! *_init is for initialization of ATM fields (in ohhi, olhi, fldsisan,
  !        hurricane_init, omic_init)
  ! *_prog is for points where primary quantities such as VMC or THIL are
  !        prognosed
  ! *_wadj is for U, V, or W points that are adjacent to W points where primary
  !        quantities are prognosed
  ! *_wstn is for all U, V, or W points in the stencil of a W point where
  !        primary quantities are prognosed
  ! *_lbcp is for lateral boundary points whose values are copied from an
  !        interior prognostic point
  ! *_vadj is for M points that are adjacent to a prognostic V point, while
  !        jtw_vadj is the representation of jtm_vadj in the pre-hex-grid stage
  !        of model grid initialization)
  ! *_wall is for U or V points along the wall of an x-z channel domain (mdomain = 3)

  integer, parameter :: jtm_grid = 1, jtu_grid = 1, jtv_grid = 1, jtw_grid = 1
  integer, parameter :: jtm_init = 2, jtu_init = 2, jtv_init = 2, jtw_init = 2
  integer, parameter :: jtm_prog = 3, jtu_prog = 3, jtv_prog = 3, jtw_prog = 3
  integer, parameter :: jtm_wadj = 4, jtu_wadj = 4, jtv_wadj = 4, jtw_wadj = 4
  integer, parameter :: jtm_wstn = 5, jtu_wstn = 5, jtv_wstn = 5, jtw_wstn = 5
  integer, parameter :: jtm_lbcp = 6, jtu_lbcp = 6, jtv_lbcp = 6, jtw_lbcp = 6
  integer, parameter :: jtm_vadj = 7,                             jtw_vadj = 7

  integer :: nstp  ! # of finest grid acoustic timesteps in coarse grid dtlong
  integer :: istp  ! Current timestep counter from 1 to nstp
  integer :: mrls  ! Number of active mesh refinement levels (MRLs)

  integer, allocatable :: mrl_begl(:)  ! MRL at beginning of long timestep
  integer, allocatable :: mrl_begr(:)  ! MRL at beginning of RK step
  integer, allocatable :: mrl_begs(:)  ! MRL at beginning of short timestep
  integer, allocatable :: mrl_ends(:)  ! MRL at end of short timestep
  integer, allocatable :: mrl_endr(:)  ! MRL at end of RK step
  integer, allocatable :: mrl_endl(:)  ! MRL at end of long timestep
  real,    allocatable :: dtrk    (:)  ! MRL RK timestep factor

  integer, allocatable :: leafstep(:)  ! flag to run leaf on any sub-timestep

  Type itab_m_vars             ! data structure for M pts (individual rank)
     logical :: loop(mloops) = .false.
     integer :: npoly = 0       ! number of V/W neighbors of this M pt
     integer :: imp = 1         ! M point from which to copy this M pt's values
     integer :: irank = -1      ! rank of parallel process at this M pt
     integer :: imglobe = 1     ! global index of this M pt (in parallel case)
     integer :: mrlm = 0        ! mesh refinement level of this M pt
     integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
     integer :: mrow = 0        ! Full row number outside nest
     integer :: ngr = 0         ! Grid number
     integer :: iv(3) = 1       ! array of V neighbors of this M pt
     integer :: im(3) = 1       ! array of M neighbors of this M pt
     integer :: iw(3) = 1       ! array of W neighbors of this M pt
  End Type itab_m_vars

  Type itab_v_vars             ! data structure for V pts (individual rank)
     logical :: loop(mloops) = .false.
     integer :: ivp = 1       ! V pt from which to copy this V pt's values
     integer :: irank = -1    ! rank of parallel process at this V pt
     integer :: ivglobe = 1   ! global index of this V pt (in parallel case)
     integer :: mrlv = 0      ! mesh refinement level of this V pt
!    integer :: im(6) = 1     ! neighbor M pts of this V pt
     integer :: im(2) = 1     ! neighbor M pts of this V pt
!    integer :: iw(4) = 1     ! neighbor W pts of this V pt
     integer :: iw(2) = 1     ! neighbor W pts of this V pt
!    integer :: iv(4) = 1     ! neighbor V pts

     real :: farw(2) = 0.     ! Interp of ARW to V control volume [VADV + VDIFF]

     real :: cosv(2) = 0.     ! cosine of angle between V and zonal dir (Voronoi)
     real :: sinv(2) = 0.     ! sine of angle between V and zonal dir (Voronoi)

     real :: dxps(2) = 0.     ! xps (eastward) displacement from neighbor W pts
     real :: dyps(2) = 0.     ! yps (northward) displacement from neighbor W pts
  End Type itab_v_vars

  Type itab_w_vars             ! data structure for W pts (individual rank)
     logical :: loop(mloops) = .false.
     integer :: npoly = 0     ! number of M/V neighbors of this W pt
     integer :: iwp = 1       ! W pt from which to copy this W pt's values
     integer :: irank = -1    ! rank of parallel process at this W pt
     integer :: iwglobe = 1   ! global index of this W pt (in parallel run)
     integer :: mrlw = 0      ! mesh refinement level of this W pt
     integer :: mrlw_orig = 0 ! original MRL of this W pt
     integer :: ngr = 0       ! Grid number
     integer :: im(7) = 1     ! neighbor M pts
     integer :: iv(7) = 1     ! neighbor V pts
     integer :: iw(7) = 1     ! neighbor W pts

     real :: dirv(7) = 0.     ! pos direction of V neighbors

     real :: farm(7) = 0.     ! Fraction of arw0 in each M point sector
     real :: farv(7) = 0.     ! Fraction of arw0 in each V point sector

     real :: gxps1(7) = 0.    ! gradient weight xe component for point 1
     real :: gyps1(7) = 0.    ! gradient weight ye component for point 1

     real :: gxps2(7) = 0.    ! gradient weight xe component for point 2
     real :: gyps2(7) = 0.    ! gradient weight ye component for point 2

     real :: unx_w = 0.       ! xe component of eastward unit normal vector
     real :: uny_w = 0.       ! ye component of eastward unit normal vector

     real :: vnx_w = 0.       ! xe component of northward unit normal vector
     real :: vny_w = 0.       ! ye component of northward unit normal vector
     real :: vnz_w = 0.       ! ze component of northward unit normal vector

     real :: ecvec_vx(7) = 0. ! factors converting V to earth cart. velocity
     real :: ecvec_vy(7) = 0. ! factors converting V to earth cart. velocity
     real :: ecvec_vz(7) = 0. ! factors converting V to earth cart. velocity

     integer :: iwnud(3) = 1  ! local nudpoly pts
     real    :: fnud (3) = 0. ! local nudpoly coeffs

     integer :: jsfc2  = 0 ! number of surface cells attached to this W column
     integer :: jland1 = 0 ! beginning land cell counter (if any land cells present)
     integer :: jland2 = 0 ! ending    land cell counter (if any land cells present)
     integer :: jlake1 = 0 ! beginning lake cell counter (if any lake cells present)
     integer :: jlake2 = 0 ! ending    lake cell counter (if any lake cells present)
     integer :: jsea1  = 0 ! beginning sea  cell counter (if any sea  cells present)
     integer :: jsea2  = 0 ! ending    sea  cell counter (if any sea  cells present)

     integer, allocatable :: iwsfc (:) ! local-rank indices of attached surface cells
     integer, allocatable :: jasfc (:) ! atm j index of attached surface cells
  End Type itab_w_vars

  Type itabg_m_vars            ! data structure for M pts (global)
     integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
     integer :: irank     = -1 ! rank of parallel process at this M pt
  End Type itabg_m_vars

  Type itabg_v_vars            ! data structure for V pts (global)
     integer :: iv_myrank = -1 ! local (parallel subdomain) index of this V pt
     integer :: irank     = -1 ! rank of parallel process at this V pt
  End Type itabg_v_vars

  Type itabg_w_vars            ! data structure for W pts (global)
     integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
     integer :: irank     = -1 ! rank of parallel process at this W pt
  End Type itabg_w_vars

  Type jtab_m_vars
     integer, allocatable :: im(:)
     integer, allocatable :: jend(:)
  End Type jtab_m_vars

  Type jtab_v_vars
     integer, allocatable :: iv(:)
     integer, allocatable :: jend(:)
  End Type jtab_v_vars

  Type jtab_w_vars
     integer, allocatable :: iw(:)
     integer, allocatable :: jend(:)
  End Type jtab_w_vars

  Type itab_m_pd_vars      ! data structure for M pts (individual rank) on para_(decomp,init)
     integer :: npoly = 0  ! number of V/W neighbors of this M pt
     integer :: imp   = 1  ! M point from which to copy this M pt's value
     integer :: iv(3) = 1  ! array of V neighbors of this M pt
     integer :: im(3) = 1  ! array of M neighbors of this M pt
     integer :: iw(3) = 1  ! array of W neighbors of this M pt
     integer :: im_myrank_imp = -1 ! local M point that corresponds to global iwp
  End Type itab_m_pd_vars

  Type itab_v_pd_vars      ! data structure for V pts (individual rank) on para_(decomp,init)
     integer :: ivp    = 1 ! V pt from which to copy this V pt's values
!    integer :: im(6)  = 1 ! neighbor M pts of this V pt
     integer :: im(2)  = 1 ! neighbor M pts of this V pt
!    integer :: iv(4)  = 1 ! neighbor V pts
!    integer :: iw(4)  = 1 ! neighbor W pts of this V pt
     integer :: iw(2)  = 1 ! neighbor W pts of this V pt
     integer :: iv_myrank_ivp = -1 ! local V point that corresponds to global ivp
  End Type itab_v_pd_vars

  Type itab_w_pd_vars      ! data structure for W pts (individual rank) on para_(decomp,init)
     integer :: iwp   = 1  ! W pt from which to copy this W pt's values
     integer :: npoly = 0  ! number of M/V neighbors of this W pt
     integer :: im(7) = 1  ! neighbor M pts of this W pt
     integer :: iv(7) = 1  ! neighbor V pts
     integer :: iw(7) = 1  ! neighbor W pts
     integer :: iw_myrank_iwp = -1 ! local W point that corresponds to global iwp
     integer :: iwnud(3) = 1
  End Type itab_w_pd_vars

  type (itab_m_vars),  allocatable, target :: itab_m(:)
  type (itab_v_vars),  allocatable, target :: itab_v(:)
  type (itab_w_vars),  allocatable, target :: itab_w(:)

  type (itab_m_pd_vars), allocatable, target :: itab_m_pd(:)
  type (itab_v_pd_vars), allocatable, target :: itab_v_pd(:)
  type (itab_w_pd_vars), allocatable, target :: itab_w_pd(:)

  type (itabg_m_vars), allocatable, target :: itabg_m(:)
  type (itabg_v_vars), allocatable, target :: itabg_v(:)
  type (itabg_w_vars), allocatable, target :: itabg_w(:)

  type (jtab_m_vars) :: jtab_m(mloops)
  type (jtab_v_vars) :: jtab_v(mloops)
  type (jtab_w_vars) :: jtab_w(mloops)

Contains

!===============================================================================

  subroutine alloc_itabs(mma, mva, mwa, input)

    implicit none

    integer, intent(in) :: mma, mva, mwa, input

    allocate (itab_m(mma))
    allocate (itab_v(mva))
    allocate (itab_w(mwa))

  end subroutine alloc_itabs

!===============================================================================

  subroutine alloc_itabs_pd(nma, nva, nwa)

    implicit none

    integer, intent(in) :: nma, nva, nwa

    allocate (itab_m_pd(nma))
    allocate (itab_v_pd(nva))
    allocate (itab_w_pd(nwa))

    ! Allocate permanent itabg data structures

    allocate (itabg_m(nma))
    allocate (itabg_v(nva))
    allocate (itabg_w(nwa))

  end subroutine alloc_itabs_pd

!===============================================================================

  subroutine fill_jtabs(mma, mva, mwa, input)

    implicit none

    integer, intent(in) :: mma, mva, mwa, input

    integer :: iw, iv, im
    integer :: iloop, j

    do iloop = 1, mloops

       ! Compute JTAB_M%IM

       allocate( jtab_m(iloop)%jend(mrls) )
       jtab_m(iloop)%jend(1:mrls) = count( itab_m(2:mma)%loop(iloop) )

       allocate( jtab_m(iloop)%im( jtab_m(iloop)%jend(1) ) )

       j = 0
       do im = 2, mma
          if (itab_m(im)%loop(iloop)) then
             j = j + 1
             jtab_m(iloop)%im(j) = im
          endif
       enddo

       ! Compute JTAB_M%IV

       allocate( jtab_v(iloop)%jend(mrls) )
       jtab_v(iloop)%jend(1:mrls) = count( itab_v(2:mva)%loop(iloop) )

       allocate( jtab_v(iloop)%iv( jtab_v(iloop)%jend(1) ) )

       j = 0
       do iv = 2, mva
          if (itab_v(iv)%loop(iloop)) then
             j = j + 1
             jtab_v(iloop)%iv(j) = iv
          endif
       enddo

       ! Compute JTAB_M%IW

       allocate( jtab_w(iloop)%jend(mrls) )
       jtab_w(iloop)%jend(1:mrls) = count( itab_w(2:mwa)%loop(iloop) )

       allocate( jtab_w(iloop)%iw( jtab_w(iloop)%jend(1) ) )

       j = 0
       do iw = 2, mwa
          if (itab_w(iw)%loop(iloop)) then
             j = j + 1
             jtab_w(iloop)%iw(j) = iw
          endif
       enddo

    enddo

  end subroutine fill_jtabs

End Module mem_ijtabs

