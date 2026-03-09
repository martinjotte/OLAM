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
   integer            :: iupdsst, isstcyclic, pom_idata
   character(pathlen) :: sst_database, pom_database

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

   real               :: seatmp0   ! Default sea sfc temperature [K]
   real               :: salnty0   ! Default sea salinity
   real               :: seaice0   ! Default sea ice fraction [0-1]
   real               :: t00sea0   ! Default freezing temperature of sea water [K]
   real               :: fssat0    ! Default reduction of sea water saturation by salinity
   real               :: dt_sea    ! sea timestep [s]
   real               :: dti_sea   ! 1 / sea timestep [s]

   real,    parameter :: emi = 1.0 ! emissivity of ice
   real,    parameter :: emw = 1.0 ! emissivity of water
   integer, parameter :: nzi = 3   ! max number of seaice layers

End Module sea_coms
