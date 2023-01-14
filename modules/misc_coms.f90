Module misc_coms

  use max_dims,    only: maxsndg, maxgrds, maxngrdll, pathlen
  use consts_coms, only: r8

  implicit none

  private :: maxsndg, maxgrds, maxngrdll, pathlen, r8

  type simtime
     integer  :: year
     integer  :: month
     integer  :: date
     real(r8) :: time
  end type simtime

  type(simtime) :: current_time ! current simulation time

  real :: topodb_cutoff

  character(pathlen) :: tmpdir = "/tmp"
  character(64)      :: expnme
  character(16)      :: runtype
  character(1)       :: timeunit
  character(pathlen) :: gridfile
  character(pathlen) :: hfilin
  character(pathlen) :: hfilepref
  character(pathlen) :: zonclim
  character(pathlen) :: topo_database(2)
  character(pathlen) :: bathym_database

  real     :: rinit  = 0.0
  real(r8) :: rinit8 = 0.0_r8

  logical :: debug_fp  = .false.
  logical :: init_nans = .false.

  integer :: io6
  integer :: initial
  integer :: ngrids
  integer :: ngrids_old
  integer :: nzp
  integer :: ilwrtyp
  integer :: iswrtyp
  integer :: icfrac
  integer :: mdomain
  integer :: nsndg
  integer :: iflag
  integer :: iyear1
  integer :: imonth1
  integer :: idate1
  integer :: ihour1
  integer :: itime1
  integer :: naddsc
  integer :: ipsflg
  integer :: itsflg
  integer :: irtsflg
  integer :: iusflg
  integer :: ioutput
  integer :: iclobber
  integer :: itopoflg
  integer :: ibathflg
  integer :: nzpp
  integer :: nscl
  integer :: nxp
  integer :: iparallel = 0
  integer :: ndz
  integer :: mstp
  integer :: nacoust

  integer :: nrk_wrtv = 3
  integer :: nrk_scal = 2

  integer :: idiffk(maxgrds)
  integer :: nqparm(maxgrds)

  real(r8) :: dtlong
  real(r8) :: frqstate
  real(r8) :: radfrq
  real(r8) :: confrq

  real :: ubmin
  real :: topref
  real :: polelat
  real :: polelon
  real :: deltax
  real :: hdz(10)
  real :: dz(10)
  real :: p_sfc
  real :: cfracrh1
  real :: cfracrh2
  real :: cfraccup

  integer, target :: ngrdll(maxgrds)
  real,    target :: grdrad(maxgrds,maxngrdll)
  real,    target :: grdlat(maxgrds,maxngrdll)
  real,    target :: grdlon(maxgrds,maxngrdll)

  real :: cflxy (maxgrds)
  real :: cflz  (maxgrds)
  real :: csz   (maxgrds)
  real :: csx   (maxgrds)
  real :: akmin (maxgrds)

  real(r8) :: dtlm
  real(r8) :: dtsm

  real, allocatable :: u01d (:)
  real, allocatable :: v01d (:)
  real, allocatable :: pr01d(:)
  real, allocatable :: th01d(:)
  real, allocatable :: dn01d(:)
  real, allocatable :: rt01d(:)

  real :: us  (maxsndg)
  real :: vs  (maxsndg)
  real :: ts  (maxsndg)
  real :: thds(maxsndg)
  real :: ps  (maxsndg)
  real :: hs  (maxsndg)
  real :: rts (maxsndg)

  real(r8) :: time8
  real(r8) :: time8p
  real(r8) :: time_istp8
  real(r8) :: time_istp8p
  real(r8) :: timmax8
  real(r8) :: s1900_init
  real(r8) :: s1900_sim
  real(r8) :: time_prevhist = 0.0_r8
  real(r8) :: time_bias ! A number small compared to the timestep

  integer  :: do_chem = 0
  integer  :: i_o3    = 0
  integer  :: i_ch4   = 0
  integer  :: i_co    = 0

Contains

!===============================================================================

  subroutine alloc_misc(mza)

    implicit none

    integer, intent(in) :: mza

    allocate (u01d (mza))
    allocate (v01d (mza))
    allocate (pr01d(mza))
    allocate (th01d(mza))
    allocate (dn01d(mza))
    allocate (rt01d(mza))

  end subroutine alloc_misc

End Module misc_coms


