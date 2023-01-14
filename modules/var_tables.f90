Module var_tables

  use consts_coms, only: r8
  implicit none

  private :: r8

  type var_tables_r

     logical,  pointer             :: lvar0_p        => null()
     logical,  pointer, contiguous :: lvar1_p(:)     => null()
     logical,  pointer, contiguous :: lvar2_p(:,:)   => null()
     logical,  pointer, contiguous :: lvar3_p(:,:,:) => null()

     integer,  pointer             :: ivar0_p        => null()
     integer,  pointer, contiguous :: ivar1_p(:)     => null()
     integer,  pointer, contiguous :: ivar2_p(:,:)   => null()
     integer,  pointer, contiguous :: ivar3_p(:,:,:) => null()

     real,     pointer             :: rvar0_p        => null()
     real,     pointer, contiguous :: rvar1_p(:)     => null()
     real,     pointer, contiguous :: rvar2_p(:,:)   => null()
     real,     pointer, contiguous :: rvar3_p(:,:,:) => null()

     real(r8), pointer             :: dvar0_p        => null()
     real(r8), pointer, contiguous :: dvar1_p(:)     => null()
     real(r8), pointer, contiguous :: dvar2_p(:,:)   => null()
     real(r8), pointer, contiguous :: dvar3_p(:,:,:) => null()

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

  character(2) :: vtypes(16) = &
       [ 'AV', 'AW', 'AM', 'AN', &  ! Atmos V, W, M, NUDGE array
         'CV', 'CW', 'CM',       &  ! "Common" SFC grid V, W, M array
         'LV', 'LW', 'LM',       &  ! Land  V, W, M array
         'RV', 'RW', 'RM',       &  ! Lake  V, W, M array
         'SV', 'SW', 'SM'        ]  ! Sea   V, W, M array

  type scalar_table

     real, pointer, contiguous :: var_p(:,:) => null()
     real, pointer, contiguous :: var_t(:,:) => null()
     real, pointer, contiguous :: sxfer(:,:) => null()
     real, pointer, contiguous :: emis (:,:) => null()

     character (len=32) :: name  = ""
     logical            :: pdef  = .true.
     logical            :: do_sgsmix = .true.
     logical            :: do_cumix  = .false.
     logical            :: do_sxfer  = .false.
     logical            :: do_emis   = .false.

  end type scalar_table

  type(scalar_table), allocatable :: scalar_tab(:)
  integer,            allocatable :: sxfer_map (:)
  integer,            allocatable :: emis_map  (:)
  integer,            allocatable :: pblmix_map(:)
  integer,            allocatable :: cumix_map (:)

  integer :: num_scalar = 0
  integer :: num_sxfer  = 0
  integer :: num_emis   = 0
  integer :: num_pblmix = 0
  integer :: num_cumix  = 0

Contains

!===============================================================================

  subroutine increment_vtable(name, stagpt, hist, noread, mpt1, lite, &
                              lvar0, lvar1, lvar2, lvar3,             &
                              ivar0, ivar1, ivar2, ivar3,             &
                              rvar0, rvar1, rvar2, rvar3,             &
                              dvar0, dvar1, dvar2, dvar3              )
    use misc_coms, only: io6
    implicit none

    character(*),      intent(in) :: name, stagpt
    logical, optional, intent(in) :: hist, noread, mpt1, lite

    logical,  target, optional,             intent(in) :: lvar0
    logical,  target, optional, contiguous, intent(in) :: lvar1(:)
    logical,  target, optional, contiguous, intent(in) :: lvar2(:,:)
    logical,  target, optional, contiguous, intent(in) :: lvar3(:,:,:)

    integer,  target, optional,             intent(in) :: ivar0
    integer,  target, optional, contiguous, intent(in) :: ivar1(:)
    integer,  target, optional, contiguous, intent(in) :: ivar2(:,:)
    integer,  target, optional, contiguous, intent(in) :: ivar3(:,:,:)

    real,     target, optional,             intent(in) :: rvar0
    real,     target, optional, contiguous, intent(in) :: rvar1(:)
    real,     target, optional, contiguous, intent(in) :: rvar2(:,:)
    real,     target, optional, contiguous, intent(in) :: rvar3(:,:,:)

    real(r8), target, optional,             intent(in) :: dvar0
    real(r8), target, optional, contiguous, intent(in) :: dvar1(:)
    real(r8), target, optional, contiguous, intent(in) :: dvar2(:,:)
    real(r8), target, optional, contiguous, intent(in) :: dvar3(:,:,:)

    integer            :: ntsize, iv
    integer, parameter :: ialloc = 20 ! Increment to increase tables

    type(var_tables_r),  allocatable :: vtab_copy(:)
    integer,             allocatable :: nptonv_copy(:)

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A NAME

    if (len_trim(name) == 0) then
       write(io6,*) "Error in subroutine increment_vtable:"
       stop         "Vtables called with no variable name."
    endif

    ! ERROR CHECKING: MAKE SURE VARIABLE HAS A STAGGER POINT

    if (len_trim(stagpt) /= 2) then
       write(io6,*) "Error in subroutine increment_vtable:"
       stop         "vtables called with invalid stagger point."
    endif

    ! ERROR CHECKING: MAKE SURE VARIABLE STAGGER POINT IS VALID

    if ( all( vtypes /= stagpt ) ) then
       write(io6,*) "Error in subroutine increment_vtable:"
       stop         "Vtables called with invalid stagger point."
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

     if     (present(lvar0)) then
        vtab_r(num_var)%lvar0_p => lvar0
     elseif (present(lvar1)) then
        vtab_r(num_var)%lvar1_p => lvar1
     elseif (present(lvar2)) then
        vtab_r(num_var)%lvar2_p => lvar2
     elseif (present(lvar3)) then
        vtab_r(num_var)%lvar3_p => lvar3

     elseif (present(ivar0)) then
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

  subroutine vtables_scalar(varp,vart,name,sxfer,emis,cu_mix,pos_def,pbl_mix)

    use misc_coms, only: io6
    implicit none

    real, target, contiguous,          intent(in) :: varp (:,:)
    real, target, contiguous,          intent(in) :: vart (:,:)
    character(len=*),       intent(in) :: name

    real, target, contiguous, optional, intent(in) :: sxfer(:,:)
    real, target, contiguous, optional, intent(in) :: emis (:,:)
    logical,      optional, intent(in) :: cu_mix
    logical,      optional, intent(in) :: pos_def
    logical,      optional, intent(in) :: pbl_mix

    integer                         :: ntsize
    integer,            parameter   :: ialloc = 20
    logical                         :: cumix
    logical                         :: pblmix

    type(scalar_table), allocatable :: scalar_copy(:)
    integer,            allocatable ::  sxfer_copy(:)
    integer,            allocatable ::   emis_copy(:)
    integer,            allocatable ::  cumix_copy(:)
    integer,            allocatable :: pblmix_copy(:)

    if (present(cu_mix)) then
       cumix = cu_mix
    else
       cumix = .false.
    endif

    if (present(pbl_mix)) then
       pblmix = pbl_mix
    else
       pblmix = .true.
    endif

    if (present(sxfer)) pblmix = .true.
    if (present(emis))  pblmix = .true.

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

    if (pblmix) then
       if (.not. allocated(pblmix_map)) allocate(pblmix_map(ialloc))
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

    scalar_tab(num_scalar)%do_sxfer = .false.

    if (present(sxfer)) then

       scalar_tab(num_scalar)%do_sxfer = .true.

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

    scalar_tab(num_scalar)%do_emis = .false.

    if (present(emis)) then

       scalar_tab(num_scalar)%do_emis = .true.

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

    scalar_tab(num_scalar)%do_cumix = cumix

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

    ! If this species will be mixed by subgrid turbulence, include it in table

    scalar_tab(num_scalar)%do_sgsmix = pblmix

    if (pblmix) then

       num_pblmix = num_pblmix + 1
       ntsize = size(pblmix_map)

       ! Increase pblmix table size if necessary
       if (num_pblmix > ntsize) then
          allocate(pblmix_copy(ntsize+ialloc))
          pblmix_copy(1:ntsize) = pblmix_map
          call move_alloc(pblmix_copy, pblmix_map)
       endif

       pblmix_map(num_pblmix) = num_scalar

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

    if     (associated(vtab_r(nv)%lvar0_p)) then

       ndims      = 1
       idims(1)   = 1

    elseif (associated(vtab_r(nv)%lvar1_p)) then

       ndims      = 1
       idims(1:1) = shape(vtab_r(nv)%lvar1_p)

    elseif (associated(vtab_r(nv)%lvar2_p)) then

       ndims      = 2
       idims(1:2) = shape(vtab_r(nv)%lvar2_p)

    elseif (associated(vtab_r(nv)%lvar3_p)) then

       ndims      = 3
       idims(1:3) = shape(vtab_r(nv)%lvar3_p)

    elseif (associated(vtab_r(nv)%ivar0_p)) then

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
