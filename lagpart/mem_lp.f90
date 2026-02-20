module mem_lp

  implicit none (external, type)

  integer, parameter :: nent = 1009 ! number of entries in table of normally distributed numbers

  real, allocatable :: relcnt(:)    ! num particles released so far from each source
  real, allocatable :: anorm (:)

  character(len=80) :: lpnmlfile  = "LPIN" ! Namelist file for LP model

  ! Variables read from namelist

  integer :: lpflag = 0  ! Run Lagrangian particles: 0=no, 1=yes
  integer :: lpturb = 0  ! Include turbulence:       0=no, 1=yes
  integer :: lpout  = 0  ! Output particles:         0=no, 1=yes
  integer :: lphist = 0  ! Read particle/source from history file: 0=no, 1=yes

  character(len=80) :: lpfilepref = ' '  ! TODO: write partical/source output files
  character(len=80) :: lphistin   = ' '  ! TODO: read particle/source files for history restart?

  real    :: lpfreq  ! Output file frequency (time interval)
  integer :: nlpsrc  ! Number of sources for this simulation

  ! Source information read from namelist:
  Type oname_lpsrc
     real          :: begtime
     real          :: endtime
     integer       :: numparts
     real          :: zbot
     real          :: ztop
     real          :: xsize
     real          :: ysize
     real          :: rotate
     integer       :: outline
     real          :: centlat
     real          :: centlon
     character(30) :: srcname
     character(30) :: specname
     real          :: wgtmol
     integer       :: ifall
     real          :: szmin
     real          :: szmax
     real          :: szpwr
     character(30) :: units
  End Type oname_lpsrc

  ! Source information stored in memory; includes all fields read in from
  ! namelist and additional derived fields not read from namelist:
  Type, extends(oname_lpsrc) :: lpsrc_vars
     integer       :: iwcent   = 0
     real          :: relperdt = 0.  ! num particles released per second per source
  End Type lpsrc_vars

  type(lpsrc_vars), allocatable :: lpsrc(:)

  ! Particle attribute datatype

  Type part_att
     real :: xep0, yep0, zep0    ! Particle earth coordinates [m]
     real :: zp                  ! Particle height above sea level
     real :: vxepp, vyepp, vzepp ! Particle turbulent earth velocity components [m/s]
     real :: fallvel             ! Particle terminal fall velocity [m/s]
     real :: mass                ! Particle mass
     real :: ppm                 ! Particle concentration
     real :: reltime             ! Particle release time [s]

     integer :: iw, k       ! Particle grid cell indices
     integer :: ncluster    !
     integer :: nsource     !
!    integer :: nspecies    !
  End Type part_att

  type(part_att), allocatable :: atp(:)

  ! Other variables

  integer :: nlp   ! Number of Lagrangian particles released so far
  real    :: dtlp  ! timestep length (s)

  real, allocatable :: vxeh(:,:), vyeh(:,:), vzeh(:,:) ! Half-forward T-cell velocity components

Contains

!==============================================================================

subroutine lp_init()

  use consts_coms, only: pio180, erad
  use mem_grid,    only: mza, mwa, xew, yew, zew, mwa
  use misc_coms,   only: dtlm
  import,          only: dtlp, nlp, nlpsrc, relcnt, vxeh, vyeh, vzeh, anorm, nent, &
                         atp, lpsrc

  implicit none (external, type)

  integer  :: nsrc, iw, nlpall
  real     :: durtot
  real     :: zesrc, xesrc, yesrc, raxis
  real     :: dist, dist0

  external :: normdist

  dtlp = dtlm
  nlp  = 0

  allocate(relcnt(nlpsrc))
  relcnt(:) = 0.

  ! Allocate half-forward T-cell velocity components

  allocate(vxeh(mza,mwa), vyeh(mza,mwa), vzeh(mza,mwa))

  ! Allocate particle data structure to total number of particles

  nlpall = sum( lpsrc(:)%numparts ) !+ nlpsrc

  if (allocated(atp)) deallocate(atp); allocate(atp(nlpall))

  atp(:)%ncluster = 0

  atp(:)%vxepp = 0.
  atp(:)%vyepp = 0.
  atp(:)%vzepp = 0.

  ! fill normally-distributed, random number array

  allocate(anorm(nent))
  call normdist(anorm,nent)

  ! Initialize Lagrangian particle sources and allocate particle attribute arrays

  do nsrc = 1,nlpsrc

     ! Find closest IW point to lpsource center and store in iwcent

     zesrc = erad  * sin(pio180 * lpsrc(nsrc)%centlat)
     raxis = erad  * cos(pio180 * lpsrc(nsrc)%centlat)
     xesrc = raxis * cos(pio180 * lpsrc(nsrc)%centlon)
     yesrc = raxis * sin(pio180 * lpsrc(nsrc)%centlon)

     dist0 = 1.e9
     do iw = 2,mwa
        dist = sqrt((xesrc - xew(iw))**2 + (yesrc - yew(iw))**2 + (zesrc - zew(iw))**2)
        if (dist0 > dist) then
           dist0 = dist
           lpsrc(nsrc)%iwcent = iw
        endif
     enddo

     ! Source duration is not less than one timestep

     if (lpsrc(nsrc)%endtime < lpsrc(nsrc)%begtime + dtlp) then
        lpsrc(nsrc)%endtime = lpsrc(nsrc)%begtime + dtlp
        durtot = dtlp
     else
        durtot = lpsrc(nsrc)%endtime - lpsrc(nsrc)%begtime
     endif

     lpsrc(nsrc)%relperdt = real(lpsrc(nsrc)%numparts) * dtlp / durtot

  enddo

end subroutine lp_init

!==============================================================================

subroutine lp_readnml()

  use misc_coms, only: iparallel, io6
  import,        only: lpflag, lpturb, lpout, lphist, lpfilepref, lphistin, lpfreq, &
                       nlpsrc, lpsrc, oname_lpsrc, lpsrc_vars, lpnmlfile

  implicit none(external, type)

  logical                        :: lexists
  integer                        :: i
  type(oname_lpsrc), allocatable :: lpsrcin(:)

  namelist /LPIN/ lpflag, lpturb, lpout, lphist, lpfilepref, lphistin, lpfreq, nlpsrc
  namelist /SRCIN/ lpsrcin

  ! For now, the LP model only works for a serial/OpenMP run. If this run is parallel,
  ! we will not even attempt to do the particle model. Will re-examine this if the
  ! need arises.

  if (iparallel == 1) return

  ! If LP namelist file does not exist, skip running particle model

  inquire( file=lpnmlfile, exist=lexists )
  if (.not. lexists) return

  ! Namelist file exists; check if we want to run LP model

  open(11, status='OLD', file=lpnmlfile)
  read(11, nml=LPIN)

! TODO: read source and existing particle information from history file?
! if (lpflag == 1 .and. lphist == 1) then
!    close(11)
!    call lp_read_histfile()
! endif

  ! If lpflag/=1 or there are no sources, turn off LP model!

  if (nlpsrc < 1) lpflag = 0

  if (lpflag /= 1) then
     close(11)
     return
  endif

  ! LP model is active; read source specifications

  allocate(lpsrcin(nlpsrc+50))  ! allow more than nlpsrc to be present in namelist
  read(11, nml=SRCIN)
  close(11)

  write(io6,*)
  write(io6,'(A,I0,A)') " LP model is active with ", nlpsrc, " sources:"

  allocate(lpsrc(nlpsrc))
  do i = 1, nlpsrc
     lpsrc(i) = lpsrc_vars( oname_lpsrc=lpsrcin(i) )
     write(io6,*) trim(lpsrc(i)%srcname)//'    ', lpsrc(i)%centlat, lpsrc(i)%centlon
  enddo

end subroutine lp_readnml

end module mem_lp
