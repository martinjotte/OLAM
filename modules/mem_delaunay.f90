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

Module mem_delaunay

  use mem_ijtabs, only: mloops

  implicit none

  private :: mloops

  Type itab_md_vars             ! data structure for M pts (individual rank)
                                ! on the Delaunay mesh
     logical :: loop(mloops) = .false. ! flag to perform each DO loop at this M pt
     integer :: npoly = 0       ! number of V/W neighbors of this M pt
     integer :: imp = 1         ! M point from which to copy this M pt's values
     integer :: mrlm = 0        ! mesh refinement level of this M pt
     integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
     integer :: ngr = 0         ! Grid number
     integer :: im(7) = 1       ! array of M neighbors of this M pt
     integer :: iu(7) = 1       ! array of U neighbors of this M pt
     integer :: iw(7) = 1       ! array of W neighbors of this M pt
  End Type itab_md_vars

  Type itab_ud_vars             ! data structure for U pts (individual rank)
                                ! on the Delaunay mesh
     logical :: loop(mloops) = .false. ! flag to perform each DO loop at this M pt
     integer :: iup = 1       ! U pt from which to copy this U pt's values
     integer :: mrlu = 0      ! mesh refinement level of this U pt
     integer :: im(2) = 1     ! neighbor M pts of this U pt
     integer :: iu(12) = 1    ! neighbor U pts
     integer :: iw(6) = 1     ! neighbor W pts
  End Type itab_ud_vars

  Type itab_wd_vars             ! data structure for W pts (individual rank)
                                ! on the Delaunay mesh
     logical :: loop(mloops) = .false. ! flag to perform each DO loop at this W pt
     integer :: npoly = 0     ! number of M/V neighbors of this W pt
     integer :: iwp = 1       ! W pt from which to copy this W pt's values
     integer :: mrlw = 0      ! mesh refinement level of this W pt
     integer :: mrlw_orig = 0 ! original MRL of this W pt
     integer :: mrow = 0      ! Full row number outside nest
     integer :: ngr = 0       ! Grid number
     integer :: im(3) = 1     ! neighbor M pts
     integer :: iu(3) = 1     ! neighbor U pts
     integer :: iw(9) = 1     ! neighbor W pts
  End Type itab_wd_vars

  Type nest_ud_vars        ! temporary U-pt data structure for spawning nested grids
     integer :: im=0, iu=0 ! new M/U pts attached to this U pt
  End Type nest_ud_vars

  Type nest_wd_vars        ! temporary W-pt data structure for spawning nested grids
     integer :: iu(3) = 0  ! new U pts attached to this W pt
     integer :: iw(3) = 0  ! new W pts attached to this W pt
  End Type nest_wd_vars

  type (itab_md_vars), allocatable :: itab_md(:)
  type (itab_ud_vars), allocatable :: itab_ud(:)
  type (itab_wd_vars), allocatable :: itab_wd(:)

  type (itab_md_vars), allocatable :: itab_md_copy(:)
  type (itab_ud_vars), allocatable :: itab_ud_copy(:)
  type (itab_wd_vars), allocatable :: itab_wd_copy(:)

  real, allocatable :: xemd(:), yemd(:), zemd(:)

  real, allocatable :: xemd_copy(:), yemd_copy(:), zemd_copy(:)

  integer :: nmd, nud, nwd

  integer :: nmd_copy = 0
  integer :: nud_copy = 0
  integer :: nwd_copy = 0

  integer, allocatable :: iwdorig(:), iwdorig_temp(:)

Contains

!===============================================================================

  subroutine alloc_itabsd(mma, mua, mwa)

    use misc_coms,  only: rinit
    use mem_ijtabs, only: mloops

    implicit none

    integer, intent(in) :: mma, mua, mwa
    integer :: imd, iud, iwd

    allocate (itab_md(mma))
    allocate (itab_ud(mua))
    allocate (itab_wd(mwa))

    allocate(xemd(mma)) ; xemd = rinit
    allocate(yemd(mma)) ; yemd = rinit
    allocate(zemd(mma)) ; zemd = rinit

    xemd(1) = 0.
    yemd(1) = 0.
    zemd(1) = 0.

  end subroutine alloc_itabsd

!===============================================================================

  subroutine copy_tri_grid()

    implicit none

    ! Save a copy of triangle structure of ATM grid in its current state of
    ! construction for subsequent independent local refinement of SURFACE grid.

    nmd_copy = nmd
    nud_copy = nud
    nwd_copy = nwd

    allocate (xemd_copy(nmd))
    allocate (yemd_copy(nmd))
    allocate (zemd_copy(nmd))

    allocate (itab_md_copy(nmd))
    allocate (itab_ud_copy(nud))
    allocate (itab_wd_copy(nwd))

    xemd_copy = xemd
    yemd_copy = yemd
    zemd_copy = zemd

    itab_md_copy = itab_md
    itab_ud_copy = itab_ud
    itab_wd_copy = itab_wd

  end subroutine copy_tri_grid

!===============================================================================

  subroutine copyback_tri_grid()

    implicit none

    integer :: iw

    ! Save a copy of triangle structure of ATM grid in its current state of
    ! construction for subsequent independent local refinement of SURFACE grid.

    if (nmd_copy == 0 .or. nud_copy == 0 .or. nwd_copy == 0) then
       write(*,*) " Error in copyback_tri_grid:"
       stop ' Delaunay mesh was not previously saved'
    endif

    nmd = nmd_copy
    nud = nud_copy
    nwd = nwd_copy

    call move_alloc(xemd_copy, xemd)
    call move_alloc(yemd_copy, yemd)
    call move_alloc(zemd_copy, zemd)

    call move_alloc(itab_md_copy, itab_md)
    call move_alloc(itab_ud_copy, itab_ud)
    call move_alloc(itab_wd_copy, itab_wd)

    allocate(iwdorig(nwd))
    do iw = 1, nwd
       iwdorig(iw) = iw
    enddo

  end subroutine copyback_tri_grid

!===============================================================================

End Module mem_delaunay
