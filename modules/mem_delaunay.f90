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
     integer :: im_orig = 0     ! im index before spawning nests
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

  integer, allocatable :: impent(:)  ! Scratch array for storing 12 pentagonal IM indices
  integer, allocatable :: impent_copy(:)

  integer, allocatable :: iwsfc_orig(:) ! Used for makesfc
  integer, allocatable :: mrl_wsfc  (:) ! Used for makesfc

Contains

!===============================================================================

  subroutine alloc_itabsd(mma, mua, mwa)

    use misc_coms,  only: rinit
    use mem_ijtabs, only: mloops

    implicit none

    integer, intent(in) :: mma, mua, mwa

    if (.not. allocated(impent)) allocate(impent(12))

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

    allocate (impent_copy(12))

    xemd_copy = xemd
    yemd_copy = yemd
    zemd_copy = zemd

    itab_md_copy = itab_md
    itab_ud_copy = itab_ud
    itab_wd_copy = itab_wd

    impent_copy = impent

  end subroutine copy_tri_grid

!===============================================================================

  subroutine copyback_tri_grid()

    implicit none

    integer :: im

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

    call move_alloc(impent_copy, impent)

    ! reset im_orig so that they refer to the atmospheric grid indices
    do im = 2, nmd
       itab_md(im)%im_orig = im
    enddo

  end subroutine copyback_tri_grid

!===============================================================================

End Module mem_delaunay
