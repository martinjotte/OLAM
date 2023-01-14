Module lake_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

   use max_dims,    only: pathlen
   use consts_coms, only: r8

   implicit none

   ! Don't re-export symbols from other modules
   private :: pathlen, r8

   real :: dt_lake     ! lake timestep [s]

   real,    parameter :: emi    = 1.0    ! emissivity of ice
   real,    parameter :: emw    = 1.0    ! emissivity of water

End Module lake_coms

