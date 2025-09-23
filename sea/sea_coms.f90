Module sea_coms

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

   ! SST VARIABLES
   integer            :: nsstfiles, isstfile, isstflg
   integer            :: iupdsst, isstcyclic
   character(pathlen) :: sst_database

   character(pathlen), allocatable :: fnames_sst  (:)
   character(14),      allocatable :: ctotdate_sst(:)
   real(r8),           allocatable :: s1900_sst   (:)

   ! SEAICE VARIABLES
   integer            :: nseaicefiles, iseaicefile, iseaiceflg
   integer            :: iupdseaice, iseaicecyclic
   character(pathlen) :: seaice_database

   character(pathlen) :: tide_database

   character(pathlen), allocatable :: fnames_seaice  (:)
   character(14),      allocatable :: ctotdate_seaice(:)
   real(r8),           allocatable :: s1900_seaice   (:)

   real :: seatmp     ! default sea sfc temperature [K]
   real :: seaice     ! default sea ice fraction [0-1]
   real :: dt_sea     ! sea timestep [s]

   integer, parameter :: nzi    = 3      ! max number of seaice layers
   real,    parameter :: emi    = 1.0    ! emissivity of ice
   real,    parameter :: emw    = 1.0    ! emissivity of water
   real,    parameter :: t00sea = 271.38 ! Freezing temperature of sea water [K]

End Module sea_coms

