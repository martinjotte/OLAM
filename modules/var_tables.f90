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

Module var_tables
  use consts_coms, only: r8
  implicit none
  
  private :: r8

  type var_tables_r
     integer,  pointer :: ivar0_p        => null()
     integer,  pointer :: ivar1_p(:)     => null()
     integer,  pointer :: ivar2_p(:,:)   => null()

     integer,  pointer :: ivar3_p(:,:,:) => null()
     real,     pointer :: rvar0_p        => null()
     real,     pointer :: rvar1_p(:)     => null()
     real,     pointer :: rvar2_p(:,:)   => null()
     real,     pointer :: rvar3_p(:,:,:) => null()

     real(r8), pointer :: dvar0_p        => null()
     real(r8), pointer :: dvar1_p(:)     => null()
     real(r8), pointer :: dvar2_p(:,:)   => null()
     real(r8), pointer :: dvar3_p(:,:,:) => null()

     character(32) :: name
     character( 2) :: stagpt

     logical :: ihist = .true.
     logical :: nread = .false.
  end type var_tables_r

  type(var_tables_r), allocatable :: vtab_r(:)
  integer,            allocatable :: nptonv(:)

  integer :: num_var  = 0
  integer :: nvar_par = 0

!-------------------------------------------------------------------

  type scalar_table

     real, pointer :: var_p(:,:)
     real, pointer :: var_t(:,:)

     character (len=32) :: name

  end type scalar_table

  type(scalar_table), allocatable :: scalar_tab(:)

  integer :: num_scalar = 0

!-------------------------------------------------------------------

  type ED_table

     integer,      pointer :: ivar1_p(:)   => null()
     integer,      pointer :: ivar2_p(:,:) => null()
     real,         pointer :: rvar1_p(:)   => null()
     real,         pointer :: rvar2_p(:,:) => null()
     real(kind=8), pointer :: dvar1_p(:)   => null()
     real(kind=8), pointer :: dvar2_p(:,:) => null()

     integer :: ndims
     integer :: idims(3) 

     logical :: mavg = .false. ! write-out on the monthly average?
     logical :: yavg = .false. ! write-out on the yearly average?
     logical :: hist = .false.

     character (len=32) :: name

  end type ED_table

  type(ED_table), allocatable :: vtab_ED(:)

  integer :: num_ED = 0

!-------------------------------------------------------------------

!!  type var_table_par
!!
!!     real, pointer :: rvar2_p(:,:)
!!
!!  end type var_table_par
!!
!!  type(var_table_par), allocatable :: vtab_par(:)
!!
!!  integer :: nvar_par = 0

Contains


!===============================================================================


  subroutine increment_vtable(name, stagpt, hist, noread, mpt1)
    use misc_coms, only: io6
    implicit none

    character(*),      intent(in) :: name, stagpt
    logical, optional, intent(in) :: hist, noread, mpt1

    character(2), parameter :: ptypes(13) = (/  &
         'AU', 'AW', 'AM', 'AN',  &  ! Atmos U, W, M, or nudging array
         'LU', 'LW', 'LM', 'LF',  &  ! Land  U, W, M, or flux array
         'SU', 'SW', 'SM', 'SF',  &  ! Sea   U, W, M, or flux array
         'CN'                     /) ! Contstant (scalar) value

    integer :: ntsize
    integer, parameter :: ialloc = 20 ! Increment to increase tables

    type(var_tables_r),  allocatable :: vtab_copy(:)
    integer,             allocatable :: nptonv_copy(:)

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A NAME

    if (len_trim(name) == 0) then
       write(io6,*) "Error in subroutine vtables:"
       stop         "Vtables called with no variable name."
    endif

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A STAGGER POINT

    if (len_trim(stagpt) == 0) then
       write(io6,*) "Error in subroutine vtables:"
       stop         "vtables called with no stagger point."
    endif
   
    ! INITIAL ALLOCATION OF TABLES IF THIS IS THE FIRST CALL

    if (.not. allocated(vtab_r)) allocate(vtab_r(ialloc))
    if (.not. allocated(nptonv)) allocate(nptonv(ialloc))

    num_var = num_var + 1
    ntsize = size(vtab_r)

    ! INCREASE VTAB SIZE IF NECESSARY
    
    if (num_var > ntsize) then
       allocate (vtab_copy (ntsize+ialloc) )
       vtab_copy(1:ntsize) = vtab_r
       call move_alloc(vtab_copy, vtab_r) 
    endif
    
    vtab_r(num_var)%name   = name
    vtab_r(num_var)%stagpt = stagpt

    if (present(hist))   vtab_r(num_var)%ihist = hist
    if (present(noread)) vtab_r(num_var)%nread = noread

!   Parallel communication table for scalars
!   (Currently only implemented for 2-D real variables)
    
    if (present(mpt1)) then

       nvar_par = nvar_par + 1
       ntsize = size(nptonv)

       ! Increase par table size if necessary

       if (nvar_par > ntsize) then
          allocate( nptonv_copy( ntsize+ialloc ))
          nptonv_copy( 1:ntsize ) = nptonv
          call move_alloc( nptonv_copy, nptonv) 
       endif

       nptonv(nvar_par) = num_var

    endif
  end subroutine increment_vtable


!===============================================================================


  subroutine increment_EDtab(name, mavg, yavg, hist)
    use misc_coms, only: io6
    implicit none

    character(*), intent(in) :: name
    logical, intent(in), optional :: mavg, yavg, hist

    type(ED_table), allocatable :: ED_copy(:)
    integer, parameter :: ialloc = 20 ! Increment to increase tables
    integer :: ntsize

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A NAME

    if (len_trim(name) == 0) then
       write(io6,*) "Error in subroutine vtables_ED:"
       stop         "Vtables_ED called with no variable name."
    endif

    ! INITIAL ALLOCATION OF TABLES IF THIS IS THE FIRST CALL

    if (.not. allocated(vtab_ED))  allocate(vtab_ED(ialloc))

    num_ED = num_ED + 1
    ntsize = size(vtab_ED)

    ! INCREASE VTAB SIZE IF NECESSARY

    if (num_ED > ntsize) then
       allocate(ED_copy(ntsize+ialloc))
       ED_copy(1:ntsize) = vtab_ED
       call move_alloc(ED_copy, vtab_ED) 
    endif

    vtab_ED(num_ED)%name = name

    if (present(hist)) vtab_ED(num_ED)%hist = hist
    if (present(mavg)) vtab_ED(num_ED)%mavg = mavg
    if (present(yavg)) vtab_ED(num_ED)%yavg = yavg

  end subroutine increment_EDtab

!===============================================================================

  subroutine vtables_scalar(varp,vart,name)
    implicit none

    real, target,     intent(in) :: varp(:,:)
    real, target,     intent(in) :: vart(:,:)
    character(len=*), intent(in) :: name

    integer                         :: ntsize
    integer, parameter              :: ialloc = 20
    type(scalar_table), allocatable :: scalar_copy(:)

    ! Initial allocation of tables if this is the first call
    if (.not. allocated(scalar_tab)) allocate(scalar_tab(ialloc))

    num_scalar = num_scalar + 1
    ntsize = size(scalar_tab)

    ! Increase scalar table size if necessary
    if (nvar_par > ntsize) then
       allocate(scalar_copy(ntsize+ialloc))
       scalar_copy(1:ntsize) = scalar_tab
       call move_alloc(scalar_copy, scalar_tab) 
    endif

    scalar_tab(num_scalar)%name = name
    scalar_tab(num_scalar)%var_p => varp
    scalar_tab(num_scalar)%var_t => vart

  end subroutine vtables_scalar

!===============================================================================

  subroutine get_vtab_dims(nv, ndims, idims)
    implicit none

    integer, intent(in)  :: nv
    integer, intent(out) :: ndims
    integer, intent(out) :: idims(3)

    ndims = 0
    idims = 0

    if     (associated(vtab_r(nv)%ivar0_p)) then

       ndims      = 1
       idims(1)   = 1

    elseif (associated(vtab_r(nv)%ivar1_p)) then

       ndims      = 1
       idims(1:1) = shape(vtab_r(nv)%ivar1_p)

    elseif (associated(vtab_r(nv)%ivar2_p)) then

       ndims      = 2
       idims(1:2) = shape(vtab_r(nv)%ivar2_p)

    elseif (associated(vtab_r(nv)%ivar3_p)) then

       ndims      = 3
       idims(1:3) = shape(vtab_r(nv)%ivar3_p)

    elseif (associated(vtab_r(nv)%rvar0_p)) then

       ndims      = 1
       idims(1)   = 1

    elseif (associated(vtab_r(nv)%rvar1_p)) then

       ndims      = 1
       idims(1:1) = shape(vtab_r(nv)%rvar1_p)

    elseif (associated(vtab_r(nv)%rvar2_p)) then

       ndims      = 2
       idims(1:2) = shape(vtab_r(nv)%rvar2_p)

    elseif (associated(vtab_r(nv)%rvar3_p)) then

       ndims      = 3
       idims(1:3) = shape(vtab_r(nv)%rvar3_p)

    elseif (associated(vtab_r(nv)%dvar0_p)) then

       ndims      = 1
       idims(1)   = 1

    elseif (associated(vtab_r(nv)%dvar1_p)) then

       ndims      = 1
       idims(1:1) = shape(vtab_r(nv)%dvar1_p)

    elseif (associated(vtab_r(nv)%dvar2_p)) then

       ndims      = 2
       idims(1:2) = shape(vtab_r(nv)%dvar2_p)

    elseif (associated(vtab_r(nv)%dvar3_p)) then

       ndims      = 3
       idims(1:3) = shape(vtab_r(nv)%dvar3_p)

    endif

  end subroutine get_vtab_dims

End Module var_tables
