Module mem_ijtabs

  use max_dims, only: maxremote

  implicit none

  private :: maxremote

  integer, parameter :: mloops = 7 ! max # non-para DO loops for M,V,W pts

  ! M, U, V, and W loop indices. The last 4 letters of these indices have the
  ! following meanings:
  !
  ! *_grid is for setting up grid parameters in the ctrlvols subroutines
  ! *_init is for initialization of ATM fields (in ohhi, olhi, fldsisan,
  !        hurricane_init, omic_init)
  ! *_prog is for points where primary quantities such as VMC or THIL are
  !        prognosed
  ! *_wadj is for U, V, or W points that are adjacent to W points where primary
  !        quantities are prognosed
  ! *_wstn is for all U, V, or W points in the stencil of a W point where
  !        primary quantities are prognosed
  ! *_lbcp is for lateral boundary points whose values are copied from an
  !        interior prognostic point
  ! *_vadj is for M points that are adjacent to a prognostic V point, while
  !        jtw_vadj is the representation of jtm_vadj in the pre-hex-grid stage
  !        of model grid initialization)
  ! *_wall is for U or V points along the wall of an x-z channel domain (mdomain = 3)

  integer, parameter :: jtm_grid = 1, jtu_grid = 1, jtv_grid = 1, jtw_grid = 1
  integer, parameter :: jtm_init = 2, jtu_init = 2, jtv_init = 2, jtw_init = 2
  integer, parameter :: jtm_prog = 3, jtu_prog = 3, jtv_prog = 3, jtw_prog = 3
  integer, parameter :: jtm_wadj = 4, jtu_wadj = 4, jtv_wadj = 4, jtw_wadj = 4
  integer, parameter :: jtm_wstn = 5, jtu_wstn = 5, jtv_wstn = 5, jtw_wstn = 5
  integer, parameter :: jtm_lbcp = 6, jtu_lbcp = 6, jtv_lbcp = 6, jtw_lbcp = 6
  integer, parameter :: jtm_vadj = 7,                             jtw_vadj = 7

  integer :: nstp  ! # of finest grid acoustic timesteps in coarse grid dtlong
  integer :: istp  ! Current timestep counter from 1 to nstp
  integer :: mrls  ! Number of active mesh refinement levels (MRLs)
  integer :: iip = 1

  integer, allocatable :: mrl_begl(:)  ! MRL at beginning of long timestep
  integer, allocatable :: mrl_begr(:)  ! MRL at beginning of RK step
  integer, allocatable :: mrl_endr(:)  ! MRL at end of RK step
  integer, allocatable :: mrl_endl(:)  ! MRL at end of long timestep
  real,    allocatable :: dtrk    (:)  ! MRL RK timestep factor

  Type itab_m_vars            ! data structure for M pts (individual rank)
     logical :: loop(mloops)
     integer :: npoly         ! number of V/W neighbors of this M pt
     integer :: imp           ! M point from which to copy this M pt's values
     integer :: irank         ! rank of parallel process at this M pt
     integer :: imglobe       ! global index of this M pt (in parallel case)
     integer :: mrlm          ! mesh refinement level of this M pt
     integer :: mrlm_orig     ! original MRL of this M pt (hex only)
     integer :: mrow          ! Full row number outside nest
     integer :: ngr           ! Grid number
     integer :: iv(3)         ! array of V neighbors of this M pt
     integer :: im(3)         ! array of M neighbors of this M pt
     integer :: iw(3)         ! array of W neighbors of this M pt
  End Type itab_m_vars

  Type itab_v_vars            ! data structure for V pts (individual rank)
     logical :: loop(mloops)
     integer :: ivp           ! V pt from which to copy this V pt's values
     integer :: irank         ! rank of parallel process at this V pt
     integer :: ivglobe       ! global index of this V pt (in parallel case)
     integer :: mrlv          ! mesh refinement level of this V pt
     integer :: im(2)         ! neighbor M pts of this V pt
     integer :: iw(2)         ! neighbor W pts of this V pt
     real    :: farw(2)       ! Interp of ARW to V control volume [VADV + VDIFF]
     real    :: cosv(2)       ! cosine of angle between V and zonal dir (Voronoi)
     real    :: sinv(2)       ! sine of angle between V and zonal dir (Voronoi)
     real    :: dxps(2)       ! xps (eastward) displacement from neighbor W pts
     real    :: dyps(2)       ! yps (northward) displacement from neighbor W pts
  End Type itab_v_vars

  Type itab_w_vars            ! data structure for W pts (individual rank)
     logical :: loop(mloops)
     integer :: npoly         ! number of M/V neighbors of this W pt
     integer :: iwp           ! W pt from which to copy this W pt's values
     integer :: irank         ! rank of parallel process at this W pt
     integer :: iwglobe       ! global index of this W pt (in parallel run)
     integer :: mrlw          ! mesh refinement level of this W pt
     integer :: mrlw_orig     ! original MRL of this W pt
     integer :: ngr           ! Grid number
     integer :: im(7)         ! neighbor M pts
     integer :: iv(7)         ! neighbor V pts
     integer :: iw(7)         ! neighbor W pts
     real    :: dirv(7)       ! pos direction of V neighbors
     real    :: farm(7)       ! Fraction of arw0 in each M point sector
     real    :: farv(7)       ! Fraction of arw0 in each V point sector
     real    :: gxps1(7)      ! gradient weight xe component for point 1
     real    :: gyps1(7)      ! gradient weight ye component for point 1
     real    :: gxps2(7)      ! gradient weight xe component for point 2
     real    :: gyps2(7)      ! gradient weight ye component for point 2
     real    :: unx_w         ! xe component of eastward unit normal vector
     real    :: uny_w         ! ye component of eastward unit normal vector
     real    :: vnx_w         ! xe component of northward unit normal vector
     real    :: vny_w         ! ye component of northward unit normal vector
     real    :: vnz_w         ! ze component of northward unit normal vector
     real    :: ecvec_vx(7)   ! factors converting V to earth cart. velocity
     real    :: ecvec_vy(7)   ! factors converting V to earth cart. velocity
     real    :: ecvec_vz(7)   ! factors converting V to earth cart. velocity
     integer :: iwnud(3)      ! local nudpoly pts
     real    :: fnud (3)      ! local nudpoly coeffs
     integer :: jsfc2         ! number of surface cells attached to this W column
     integer :: jland1        ! beginning land cell counter (if any land cells present)
     integer :: jland2        ! ending    land cell counter (if any land cells present)
     integer :: jlake1        ! beginning lake cell counter (if any lake cells present)
     integer :: jlake2        ! ending    lake cell counter (if any lake cells present)
     integer :: jsea1         ! beginning sea  cell counter (if any sea  cells present)
     integer :: jsea2         ! ending    sea  cell counter (if any sea  cells present)

     integer, allocatable :: iwsfc(:) ! local-rank indices of attached surface cells
     integer, allocatable :: jasfc(:) ! atm j index of attached surface cells
  End Type itab_w_vars

  Type itabg_m_vars            ! data structure for M pts (global)
     integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
     integer :: irank     = -1 ! rank of parallel process at this M pt
  End Type itabg_m_vars

  Type itabg_v_vars            ! data structure for V pts (global)
     integer :: iv_myrank = -1 ! local (parallel subdomain) index of this V pt
     integer :: irank     = -1 ! rank of parallel process at this V pt
  End Type itabg_v_vars

  Type itabg_w_vars            ! data structure for W pts (global)
     integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
     integer :: irank     = -1 ! rank of parallel process at this W pt
  End Type itabg_w_vars

  Type jtab_m_vars
     integer, allocatable :: im(:)
     integer              :: jend
  End Type jtab_m_vars

  Type jtab_v_vars
     integer, allocatable :: iv(:)
     integer              :: jend
  End Type jtab_v_vars

  Type jtab_w_vars
     integer, allocatable :: iw(:)
     integer              :: jend
  End Type jtab_w_vars

  Type itab_m_pd_vars      ! data structure for M pts (individual rank) on para_(decomp,init)
     integer :: npoly = 0  ! number of V/W neighbors of this M pt
     integer :: imp   = 1  ! M point from which to copy this M pt's value
     integer :: iv(3) = 1  ! array of V neighbors of this M pt
     integer :: im(3) = 1  ! array of M neighbors of this M pt
     integer :: iw(3) = 1  ! array of W neighbors of this M pt
     integer :: im_myrank_imp = -1 ! local M point that corresponds to global iwp
  End Type itab_m_pd_vars

  Type itab_v_pd_vars      ! data structure for V pts (individual rank) on para_(decomp,init)
     integer :: ivp    = 1 ! V pt from which to copy this V pt's values
     integer :: im(2)  = 1 ! neighbor M pts of this V pt
     integer :: iw(2)  = 1 ! neighbor W pts of this V pt
     integer :: iv_myrank_ivp = -1 ! local V point that corresponds to global ivp
  End Type itab_v_pd_vars

  Type itab_w_pd_vars      ! data structure for W pts (individual rank) on para_(decomp,init)
     integer :: iwp   = 1  ! W pt from which to copy this W pt's values
     integer :: npoly = 0  ! number of M/V neighbors of this W pt
     integer :: im(7) = 1  ! neighbor M pts of this W pt
     integer :: iv(7) = 1  ! neighbor V pts
     integer :: iw(7) = 1  ! neighbor W pts
     integer :: iw_myrank_iwp = -1 ! local W point that corresponds to global iwp
     integer :: iwnud(3) = 1
  End Type itab_w_pd_vars

  type (itab_m_vars),  allocatable, target :: itab_m(:)
  type (itab_v_vars),  allocatable, target :: itab_v(:)
  type (itab_w_vars),  allocatable, target :: itab_w(:)

  type (itab_m_pd_vars), allocatable, target :: itab_m_pd(:)
  type (itab_v_pd_vars), allocatable, target :: itab_v_pd(:)
  type (itab_w_pd_vars), allocatable, target :: itab_w_pd(:)

  type (itabg_m_vars), allocatable, target :: itabg_m(:)
  type (itabg_v_vars), allocatable, target :: itabg_v(:)
  type (itabg_w_vars), allocatable, target :: itabg_w(:)

  type (jtab_m_vars) :: jtab_m(mloops)
  type (jtab_v_vars) :: jtab_v(mloops + 2)
  type (jtab_w_vars) :: jtab_w(mloops)

Contains

!===============================================================================

  subroutine alloc_itabs(mma, mva, mwa)

    implicit none

    integer, intent(in) :: mma, mva, mwa
    integer             :: im,  iv,  iw

    allocate( itab_m(mma) )
    allocate( itab_v(mva) )
    allocate( itab_w(mwa) )

    !$omp parallel
    !$omp do
    do im = 1, mma
       itab_m(im) = itab_m_vars( loop      = .false., &
                                 npoly     =  0,      &
                                 imp       =  1,      &
                                 irank     = -1,      &
                                 imglobe   =  1,      &
                                 mrlm      =  0,      &
                                 mrlm_orig =  0,      &
                                 mrow      =  0,      &
                                 ngr       =  0,      &
                                 iv        =  1,      &
                                 im        =  1,      &
                                 iw        =  1       )
    enddo
    !$omp end do nowait

    !$omp do
    do iv = 1, mva
       itab_v(iv) = itab_v_vars( loop    = .false., &
                                 ivp     =  1,      &
                                 irank   = -1,      &
                                 ivglobe =  1,      &
                                 mrlv    =  0,      &
                                 im      =  1,      &
                                 iw      =  1,      &
                                 farw    =  0.,     &
                                 cosv    =  0.,     &
                                 sinv    =  0.,     &
                                 dxps    =  0.,     &
                                 dyps    =  0.      )
    enddo
    !$omp end do nowait

    !$omp do
    do iw = 1, mwa
       itab_w(iw) = itab_w_vars( loop      = .false., &
                                 npoly     =  0,      &
                                 iwp       =  1,      &
                                 irank     = -1,      &
                                 iwglobe   =  1,      &
                                 mrlw      =  0,      &
                                 mrlw_orig =  0,      &
                                 ngr       =  0,      &
                                 im        =  1,      &
                                 iv        =  1,      &
                                 iw        =  1,      &
                                 dirv      =  0.,     &
                                 farm      =  0.,     &
                                 farv      =  0.,     &
                                 gxps1     =  0.,     &
                                 gyps1     =  0.,     &
                                 gxps2     =  0.,     &
                                 gyps2     =  0.,     &
                                 unx_w     =  0.,     &
                                 uny_w     =  0.,     &
                                 vnx_w     =  0.,     &
                                 vny_w     =  0.,     &
                                 vnz_w     =  0.,     &
                                 ecvec_vx  =  0.,     &
                                 ecvec_vy  =  0.,     &
                                 ecvec_vz  =  0.,     &
                                 iwnud     =  1,      &
                                 fnud      =  0.,     &
                                 jsfc2     =  0,      &
                                 jland1    =  0,      &
                                 jland2    =  0,      &
                                 jlake1    =  0,      &
                                 jlake2    =  0,      &
                                 jsea1     =  0,      &
                                 jsea2     =  0,      &
                                 iwsfc     = null(),  &
                                 jasfc     = null()   )
   enddo
   !$omp end do nowait
   !$omp end parallel

  end subroutine alloc_itabs

!===============================================================================

  subroutine alloc_itabs_pd(nma, nva, nwa)

    implicit none

    integer, intent(in) :: nma, nva, nwa

    allocate( itab_m_pd(nma) )
    allocate( itab_v_pd(nva) )
    allocate( itab_w_pd(nwa) )

    ! Allocate permanent itabg data structures

    allocate( itabg_m(nma) )
    allocate( itabg_v(nva) )
    allocate( itabg_w(nwa) )

  end subroutine alloc_itabs_pd

!===============================================================================

  subroutine fill_jtabs(mma, mva, mwa, myrank, iparallel)

    implicit none

    integer, intent(in) :: mma, mva, mwa
    integer, intent(in) :: myrank, iparallel

    integer :: iw, iv, im, iw1, iw2
    integer :: iloop, j, np, nb

    do iloop = 1, mloops

       ! Compute JTAB_M%IM

       jtab_m(iloop)%jend = count( itab_m(2:mma)%loop(iloop) )

       allocate( jtab_m(iloop)%im( jtab_m(iloop)%jend ) )

       !$omp parallel do
       !$ do j = 1, jtab_m(iloop)%jend
       !$    jtab_m(iloop)%im(j) = 0
       !$ enddo
       !$omp end parallel do

       j = 0
       do im = 2, mma
          if (itab_m(im)%loop(iloop)) then
             j = j + 1
             jtab_m(iloop)%im(j) = im
          endif
       enddo

       ! Compute JTAB_V%IV

       jtab_v(iloop)%jend = count( itab_v(2:mva)%loop(iloop) )

       allocate( jtab_v(iloop)%iv( jtab_v(iloop)%jend ) )

       !$omp parallel do
       !$ do j = 1, jtab_v(iloop)%jend
       !$    jtab_v(iloop)%iv(j) = 0
       !$ enddo
       !$omp end parallel do

       j = 0
       do iv = 2, mva
          if (itab_v(iv)%loop(iloop)) then
             j = j + 1
             jtab_v(iloop)%iv(j) = iv
          endif
       enddo

       ! Compute JTAB_W%IW

       jtab_w(iloop)%jend = count( itab_w(2:mwa)%loop(iloop) )

       allocate( jtab_w(iloop)%iw( jtab_w(iloop)%jend ) )

       !$omp parallel do
       !$ do j = 1, jtab_w(iloop)%jend
       !$    jtab_w(iloop)%iw(j) = 0
       !$ enddo
       !$omp end parallel do

       j = 0
       do iw = 2, mwa
          if (itab_w(iw)%loop(iloop)) then
             j = j + 1
             jtab_w(iloop)%iw(j) = iw
          endif
       enddo

    enddo

    if (iparallel == 1) then

       iip = 2

       np = 0
       do j = 1, jtab_v(jtv_wadj)%jend ; iv = jtab_v(jtv_wadj)%iv(j)
          iw1 = itab_v(iv)%iw(1) ; iw2 = itab_v(iv)%iw(2)
          if (itab_w(iw1)%irank == myrank .and. itab_w(iw2)%irank == myrank) then
             np = np + 1
          endif
       enddo

       jtab_v(mloops+1)%jend = np
       allocate( jtab_v(mloops+1)%iv( jtab_v(mloops+1)%jend ) )

       jtab_v(mloops+2)%jend = jtab_v(jtv_wadj)%jend - np
       allocate( jtab_v(mloops+2)%iv( jtab_v(mloops+2)%jend ) )

       np = 0
       nb = 0
       do j = 1, jtab_v(jtv_wadj)%jend ; iv = jtab_v(jtv_wadj)%iv(j)
          iw1 = itab_v(iv)%iw(1) ; iw2 = itab_v(iv)%iw(2)
          if (itab_w(iw1)%irank == myrank .and. itab_w(iw2)%irank == myrank) then
             np = np + 1
             jtab_v(mloops+1)%iv(np) = iv
          else
             nb = nb + 1
             jtab_v(mloops+2)%iv(nb) = iv
          endif
       enddo

    else

       iip = 1

       jtab_v(mloops+1)%jend = jtab_v(jtv_wadj)%jend
       allocate( jtab_v(mloops+1)%iv( jtab_v(mloops+1)%jend ) )
       jtab_v(mloops+1)%iv(:) = jtab_v(jtv_wadj)%iv(:)

       jtab_v(mloops+2)%jend = 0
       allocate( jtab_v(mloops+2)%iv( jtab_v(mloops+2)%jend ) )

    endif

  end subroutine fill_jtabs

End Module mem_ijtabs

