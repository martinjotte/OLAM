subroutine lake_startup()

use mem_lake,  only: mlake, alloc_lake, filltab_lake

implicit none

! Subroutine LAKE_STARTUP allocates some lake arrays

! THIS SUBROUTINE DOES NOT INITIALIZE canopy temperature and moisture
! values, which depend on atmospheric conditions.

  call alloc_lake(mlake)
  call filltab_lake()

end subroutine lake_startup
