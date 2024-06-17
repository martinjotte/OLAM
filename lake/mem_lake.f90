Module mem_lake

  implicit none

  integer :: nzlake    ! Max number of vertical levels in lake model

  integer :: nlake     ! Total # of lake cell W pts in model domain
  integer :: mlake     ! Total # of lake cell W pts in model parallel sub-domain

  integer :: onlake    ! Offset # of lake cell W pts in model global domain
  integer :: omlake    ! Offset # of lake cell W pts in model parallel sub-domain
                        ! (ilake + omlake = iwsfc)

  real, parameter :: timescale_lake = 2. * 86400. ! Lake cell-to-cell equilibration time scale [s]

! LAKE GRID TABLES

  Type itab_lake_vars
     integer :: iwglobe = 1
  End type

  type (itab_lake_vars), allocatable, target :: itab_lake(:)

! lake MODEL VARIABLES

  Type lake_vars

     real, allocatable :: depth       (:) ! lake water mean depth [m]
     real, allocatable :: lake_energy (:) ! lake energy [J/kg]

  End Type lake_vars

  type (lake_vars), target :: lake

Contains

!=========================================================================

  subroutine alloc_lake(mlake)

     use misc_coms, only: rinit

     implicit none

     integer, intent(in) :: mlake

!    Allocate and initialize lake arrays

     allocate (lake%depth       (mlake)) ; lake%depth        = rinit
     allocate (lake%lake_energy (mlake)) ; lake%lake_energy  = rinit

  end subroutine alloc_lake

!=========================================================================

  subroutine filltab_lake()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(lake%depth))        call increment_vtable('LAKE%DEPTH'       , 'RW', rvar1=lake%depth)
     if (allocated(lake%lake_energy))  call increment_vtable('LAKE%LAKE_ENERGY' , 'RW', rvar1=lake%lake_energy)

  end subroutine filltab_lake

End Module mem_lake
