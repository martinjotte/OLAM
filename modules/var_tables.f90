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
     logical :: ilite = .false.

  end type var_tables_r

  type(var_tables_r), allocatable :: vtab_r(:)
  integer,            allocatable :: nptonv(:)

  integer :: num_var  = 0
  integer :: nvar_par = 0
  integer :: num_lite = 0

!-------------------------------------------------------------------

  type scalar_table

     real, pointer, contiguous :: var_p(:,:) => null()
     real, pointer, contiguous :: var_t(:,:) => null()
     real, pointer, contiguous :: sxfer(:,:) => null()
     real, pointer, contiguous :: emis (:,:) => null()

     character (len=32) :: name = ""
     logical            :: pdef = .true.

  end type scalar_table

  type(scalar_table), allocatable :: scalar_tab(:)
  integer,            allocatable :: sxfer_map (:)
  integer,            allocatable :: emis_map  (:)
  integer,            allocatable :: cumix_map (:)

  integer :: num_scalar = 0
  integer :: num_sxfer  = 0
  integer :: num_emis   = 0
  integer :: num_cumix  = 0

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

Contains

!===============================================================================

  subroutine increment_vtable(name, stagpt, hist, noread, mpt1, lite, &
                              ivar0, ivar1, ivar2, ivar3,             &
                              rvar0, rvar1, rvar2, rvar3,             &
                              dvar0, dvar1, dvar2, dvar3              )
    use misc_coms, only: io6
    implicit none

    character(*),      intent(in) :: name, stagpt
    logical, optional, intent(in) :: hist, noread, mpt1, lite

    integer,  target, optional, intent(in) :: ivar0
    integer,  target, optional, intent(in) :: ivar1(:)
    integer,  target, optional, intent(in) :: ivar2(:,:)
    integer,  target, optional, intent(in) :: ivar3(:,:,:)

    real,     target, optional, intent(in) :: rvar0
    real,     target, optional, intent(in) :: rvar1(:)
    real,     target, optional, intent(in) :: rvar2(:,:)
    real,     target, optional, intent(in) :: rvar3(:,:,:)

    real(r8), target, optional, intent(in) :: dvar0
    real(r8), target, optional, intent(in) :: dvar1(:)
    real(r8), target, optional, intent(in) :: dvar2(:,:)
    real(r8), target, optional, intent(in) :: dvar3(:,:,:)

    character(2), parameter :: ptypes(11) = (/  &
         'AV', 'AW', 'AM', 'AN',  &  ! Atmos V, W, M, NUDGE array
         'LU', 'LW', 'LM',        &  ! Land  U, W, M array
         'SU', 'SW', 'SM',        &  ! Sea   U, W, M array
         'CN'                     /) ! Contstant (scalar)

    integer :: ntsize, iv
    integer, parameter :: ialloc = 20 ! Increment to increase tables

    type(var_tables_r),  allocatable :: vtab_copy(:)
    integer,             allocatable :: nptonv_copy(:)

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A NAME

    if (len_trim(name) == 0) then
       write(io6,*) "Error in subroutine increment_vtable:"
       stop         "Vtables called with no variable name."
    endif

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A STAGGER POINT

    if (len_trim(stagpt) == 0) then
       write(io6,*) "Error in subroutine increment_vtable:"
       stop         "vtables called with no stagger point."
    endif

    ! ERROR CHECKING: MAKE SURE NAMES ARE UNIQUE
    if (num_var > 0) then
       do iv = 1, num_var
          if (vtab_r(iv)%name == name) then
             write(io6,*) "Error in subroutine increment_vtable:"
             write(io6,*) "Name " // trim(name) // " already used."
             stop
          endif
       enddo
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
    if (present(lite))   vtab_r(num_var)%ilite = lite

     if     (present(ivar0)) then
        vtab_r(num_var)%ivar0_p => ivar0
     elseif (present(ivar1)) then
        vtab_r(num_var)%ivar1_p => ivar1
     elseif (present(ivar2)) then
        vtab_r(num_var)%ivar2_p => ivar2
     elseif (present(ivar3)) then
        vtab_r(num_var)%ivar3_p => ivar3

     elseif (present(rvar0)) then
        vtab_r(num_var)%rvar0_p => rvar0
     elseif (present(rvar1)) then
        vtab_r(num_var)%rvar1_p => rvar1
     elseif (present(rvar2)) then
        vtab_r(num_var)%rvar2_p => rvar2
     elseif (present(rvar3)) then
        vtab_r(num_var)%rvar3_p => rvar3

     elseif (present(dvar0)) then
        vtab_r(num_var)%dvar0_p => dvar0
     elseif (present(dvar1)) then
        vtab_r(num_var)%dvar1_p => dvar1
     elseif (present(dvar2)) then
        vtab_r(num_var)%dvar2_p => dvar2
     elseif (present(dvar3)) then
        vtab_r(num_var)%dvar3_p => dvar3
     endif

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

  subroutine increment_EDtab(name, mavg, yavg, hist, &
                             ivar1, ivar2,           &
                             rvar1, rvar2,           &
                             dvar1, dvar2            )

    use misc_coms, only: io6
    implicit none

    character(*), intent(in) :: name
    logical, intent(in), optional :: mavg, yavg, hist

    type(ED_table), allocatable :: ED_copy(:)
    integer, parameter :: ialloc = 20 ! Increment to increase tables
    integer :: ntsize

    integer,  target, optional, intent(in) :: ivar1(:)
    integer,  target, optional, intent(in) :: ivar2(:,:)

    real,     target, optional, intent(in) :: rvar1(:)
    real,     target, optional, intent(in) :: rvar2(:,:)

    real(r8), target, optional, intent(in) :: dvar1(:)
    real(r8), target, optional, intent(in) :: dvar2(:,:)

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

    if     (present(ivar1)) then
       vtab_r(num_var)%ivar1_p => ivar1
    elseif (present(ivar2)) then
       vtab_r(num_var)%ivar2_p => ivar2

    elseif (present(rvar1)) then
       vtab_r(num_var)%rvar1_p => rvar1
    elseif (present(rvar2)) then
       vtab_r(num_var)%rvar2_p => rvar2

    elseif (present(dvar1)) then
       vtab_r(num_var)%dvar1_p => dvar1
    elseif (present(dvar2)) then
       vtab_r(num_var)%dvar2_p => dvar2
    endif

  end subroutine increment_EDtab

!===============================================================================

  subroutine vtables_scalar(varp,vart,name,sxfer,emis,cu_mix,pos_def)
    use misc_coms, only: io6
    implicit none

    real, target, contiguous,          intent(in) :: varp (:,:)
    real, target, contiguous,          intent(in) :: vart (:,:)
    character(len=*),       intent(in) :: name

    real, target, contiguous, optional, intent(in) :: sxfer(:,:)
    real, target, contiguous, optional, intent(in) :: emis (:,:)
    logical,      optional, intent(in) :: cu_mix
    logical,      optional, intent(in) :: pos_def

    integer                         :: ntsize
    integer,            parameter   :: ialloc = 20
    logical                         :: cumix

    type(scalar_table), allocatable :: scalar_copy(:)
    integer,            allocatable ::  sxfer_copy(:)
    integer,            allocatable ::   emis_copy(:)
    integer,            allocatable ::  cumix_copy(:)

    if (present(cu_mix)) then
       cumix = cu_mix
    else
       cumix = .false.
    endif

    ! Initial allocation of tables if this is the first call

    if (.not. allocated(scalar_tab)) allocate(scalar_tab(ialloc))

    if (present(sxfer)) then
       if (.not. allocated(sxfer_map)) allocate(sxfer_map(ialloc))
    endif

    if (present(emis)) then
       if (.not. allocated(emis_map)) allocate(emis_map(ialloc))
    endif

    if (cumix) then
       if (.not. allocated(cumix_map)) allocate(cumix_map(ialloc))
    endif

    num_scalar = num_scalar + 1
    ntsize = size(scalar_tab)

    ! Increase scalar table size if necessary
    if (num_scalar > ntsize) then
       allocate(scalar_copy(ntsize+ialloc))
       scalar_copy(1:ntsize) = scalar_tab
       call move_alloc(scalar_copy, scalar_tab) 
    endif

    scalar_tab(num_scalar)%name  =  name
    scalar_tab(num_scalar)%var_p => varp
    scalar_tab(num_scalar)%var_t => vart

    if (present(pos_def)) then
       scalar_tab(num_scalar)%pdef = pos_def
    else
       scalar_tab(num_scalar)%pdef = .true.
    endif

    ! If this species has surface transfer, include it in the table
    
    if (present(sxfer)) then

       num_sxfer = num_sxfer + 1
       ntsize = size(sxfer_map)

       ! Increase sxfer table size if necessary
       if (num_sxfer > ntsize) then
          allocate(sxfer_copy(ntsize+ialloc))
          sxfer_copy(1:ntsize) = sxfer_map
          call move_alloc(sxfer_copy, sxfer_map)
       endif

       scalar_tab(num_scalar)%sxfer => sxfer
       sxfer_map(num_sxfer) = num_scalar
       
    endif

    ! If this species has emissions, include it in the table
    
    if (present(emis)) then

       num_emis = num_emis + 1
       ntsize = size(emis_map)

       ! Increase emis table size if necessary
       if (num_emis > ntsize) then
          allocate(emis_copy(ntsize+ialloc))
          emis_copy(1:ntsize) = emis_map
          call move_alloc(emis_copy, emis_map)
       endif

       scalar_tab(num_scalar)%emis => emis
       emis_map(num_emis) = num_scalar
       
       write(io6,*) "Emissions:", scalar_tab(num_scalar)%name, num_emis, num_scalar

    endif

    ! If this species will be mixed by subgrid cumulus, include it in table

    if (cumix) then

       num_cumix = num_cumix + 1
       ntsize = size(cumix_map)

       ! Increase cumix table size if necessary
       if (num_cumix > ntsize) then
          allocate(cumix_copy(ntsize+ialloc))
          cumix_copy(1:ntsize) = cumix_map
          call move_alloc(cumix_copy, cumix_map)
       endif

       cumix_map(num_cumix) = num_scalar

    endif

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
