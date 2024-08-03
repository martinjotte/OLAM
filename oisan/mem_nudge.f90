Module mem_nudge

  use consts_coms, only: r8
  implicit none

  private :: r8

  integer, parameter :: nloops_wnud = 100 ! # WNUD DO loops for para

  Type itab_wnud_vars          ! data structure for WNUD points (individual rank)
     logical :: loop(nloops_wnud) = .false. ! flag to perform DO loop at this WNUD pt
     integer :: npoly = 0      ! number of W neighbors of this WNUD pt
     integer :: irank = -1     ! rank of process at this WNUD pt (for hist write only)
     integer :: iwnud(6) = 1   ! array of WNUD neighbors of this WNUD pt
     integer :: iwnudglobe = 1 ! global index of WNUD point
  End Type itab_wnud_vars

  type (itab_wnud_vars), allocatable, target :: itab_wnud(:)

  Type itabg_wnud_vars            ! data structure for WNUD pts (global)
     integer :: iwnud_myrank = -1 ! local (parallel subdomain) index of this WNUD pt
     integer :: irank = -1        ! rank of process at this WNUD pt
  End Type itabg_wnud_vars

  type (itabg_wnud_vars), allocatable, target :: itabg_wnud(:)

  Type jtab_wnud_vars
     integer, allocatable :: iwnud(:)
     integer              :: jend
  End Type jtab_wnud_vars

  real,    allocatable ::    rho_obsp(:,:)
  real,    allocatable ::  theta_obsp(:,:)
  real,    allocatable ::    rrw_obsp(:,:)
  real,    allocatable :: uzonal_obsp(:,:)
  real,    allocatable :: umerid_obsp(:,:)
  real,    allocatable ::  ozone_obsp(:,:)

  real,    allocatable ::    rho_obsf(:,:)
  real,    allocatable ::  theta_obsf(:,:)
  real,    allocatable ::    rrw_obsf(:,:)
  real,    allocatable :: uzonal_obsf(:,:)
  real,    allocatable :: umerid_obsf(:,:)
  real,    allocatable ::  ozone_obsf(:,:)

  real(r8),allocatable ::   rhot_nud(:,:)

  ! This set is only needed for spectral nudging:
  real,    allocatable ::    rho_obs(:,:)
  real,    allocatable ::  theta_obs(:,:)
  real,    allocatable ::    rrw_obs(:,:)
  real,    allocatable :: uzonal_obs(:,:)
  real,    allocatable :: umerid_obs(:,:)

  ! This set is only needed for spectral nudging:
  real,    allocatable ::    rho_sim(:,:)
  real,    allocatable ::  theta_sim(:,:)
  real,    allocatable ::    rrw_sim(:,:)
  real,    allocatable :: uzonal_sim(:,:)
  real,    allocatable :: umerid_sim(:,:)

  ! This set is only needed for spectral nudging:
  real(r8),allocatable ::   volwnudi(:,:)

  ! This set is only needed for spectral nudging:
  type (jtab_wnud_vars) :: jtab_wnud(nloops_wnud)

  ! This set is only needed for spectral nudging:
  real,    allocatable ::      xewnud(:)
  real,    allocatable ::      yewnud(:)
  real,    allocatable ::      zewnud(:)

  integer :: nudflag
  integer :: nudnxp
  integer :: nnudfiles
  integer :: nwnud = 1
  integer :: mwnud = 1

  real    :: tnudcent

  integer :: o3nudflag = 0
  real    :: tnudi_o3
  real    :: o3nudpress

  real(r8) :: past_mass_fact
  real(r8) :: future_mass_fact

Contains

  subroutine alloc_nudge1(lwnud,ii)

    use misc_coms, only: io6, rinit
    implicit none
    integer, intent(in) :: lwnud, ii

!   Allocate arrays based on options (if necessary)
!   These are only needed for spectral nudging

    write(io6,*) 'allocating nudge1 ', lwnud

    if (ii > 1) then
       allocate (itab_wnud(lwnud))
    endif

    allocate (xewnud(lwnud)) ; xewnud = rinit
    allocate (yewnud(lwnud)) ; yewnud = rinit
    allocate (zewnud(lwnud)) ; zewnud = rinit

  end subroutine alloc_nudge1

!=========================================================================

  subroutine alloc_nudge2(mza,mwa)

    use misc_coms,   only: io6, rinit, rinit8
    implicit none

    integer, intent(in) :: mza, mwa

!   Allocate arrays based on options (if necessary)

    write(io6,*) 'allocating nudge2 ',mza,mwnud

    allocate ( rhot_nud(mza,mwa) )

! The following should only be allocated for spectral nudging

    if ( nudnxp > 0) then

       allocate (   rho_obsp(mza,mwnud))
       allocate ( theta_obsp(mza,mwnud))
       allocate (   rrw_obsp(mza,mwnud))
       allocate (uzonal_obsp(mza,mwnud))
       allocate (umerid_obsp(mza,mwnud))

       allocate (   rho_obsf(mza,mwnud))
       allocate ( theta_obsf(mza,mwnud))
       allocate (   rrw_obsf(mza,mwnud))
       allocate (uzonal_obsf(mza,mwnud))
       allocate (umerid_obsf(mza,mwnud))

       allocate (    rho_obs(mza,mwnud))
       allocate (  theta_obs(mza,mwnud))
       allocate (    rrw_obs(mza,mwnud))
       allocate ( uzonal_obs(mza,mwnud))
       allocate ( umerid_obs(mza,mwnud))

    endif

  end subroutine alloc_nudge2

  !=========================================================================

  subroutine alloc_nudge_o3(mza,mwa)

!   use misc_coms,  only: io6, rinit, i_o3
    implicit none

    integer, intent(in) :: mza, mwa

!   Note:
!   This routine will only be necessary if we do spectral nudging of ozone

!   if (i_o3 /= 0) then
!      write(io6,*) 'allocating nudge_o3 ', mza, mwa
!
!      allocate (ozone_obsp(mza,mwa)) ; ozone_obsp = rinit
!      allocate (ozone_obsf(mza,mwa)) ; ozone_obsf = rinit
!   endif

  end subroutine alloc_nudge_o3

!=========================================================================

  subroutine filltab_nudge()

!   use var_tables, only: increment_vtable
    implicit none

!   NOTE:
!   Obs-nudging arrays do not need to be in the vtables for history
!   read/writes or parallel communication. The gridded observational fields
!   are read in on initialization or history restart

!   character(2) :: stagpt

    ! obs (point-by-point) nudging is done at W points, while spectral
    ! nudging is done on a separate nudging mesh. Currently, parallel
    ! output on the nudging mesh in not implemented

!   if (nudnxp == 0) then
!      stagpt = 'AW'
!   else
!      stagpt = 'AN'
!   endif
!
!   if (allocated(rho_obsp))    call increment_vtable('RHO_OBSP',    stagpt, noread=.true., rvar2=rho_obsp)
!
!   if (allocated(theta_obsp))  call increment_vtable('THETA_OBSP',  stagpt, noread=.true., rvar2=theta_obsp)
!
!   if (allocated(rrw_obsp))    call increment_vtable('RRW_OBSP',    stagpt, noread=.true., rvar2=rrw_obsp)
!
!   if (allocated(uzonal_obsp)) call increment_vtable('UZONAL_OBSP', stagpt, noread=.true., rvar2=uzonal_obsp)
!
!   if (allocated(umerid_obsp)) call increment_vtable('UMERID_OBSP', stagpt, noread=.true., rvar2=umerid_obsp)
!
!   if (allocated(rho_obsf))    call increment_vtable('RHO_OBSF',    stagpt, noread=.true., rvar2=rho_obsf)
!
!   if (allocated(theta_obsf))  call increment_vtable('THETA_OBSF',  stagpt, noread=.true., rvar2=theta_obsf)
!
!   if (allocated(rrw_obsf))    call increment_vtable('RRW_OBSF',    stagpt, noread=.true., rvar2=rrw_obsf)
!
!   if (allocated(uzonal_obsf)) call increment_vtable('UZONAL_OBSF', stagpt, noread=.true., rvar2=uzonal_obsf)
!
!   if (allocated(umerid_obsf)) call increment_vtable('UMERID_OBSF', stagpt, noread=.true., rvar2=umerid_obsf)

  end subroutine filltab_nudge

!=========================================================================

  subroutine filltab_nudge_o3()

!   use var_tables, only: increment_vtable
    implicit none

!   Ozone currently only uses non-spectral nudging, which is only
!   computed at W points.

!   NOTE:
!   Obs-nudging arrays do not need to be in the vtables for history
!   read/writes or parallel communication. The gridded observational fields
!   are read in on initialization or history restart

!   if (allocated(ozone_obsp)) &
!        call increment_vtable('OZONE_OBSP', 'AW', noread=.true., rvar2=ozone_obsp)
!
!   if (allocated(ozone_obsf)) &
!        call increment_vtable('OZONE_OBSF', 'AW', noread=.true., rvar2=ozone_obsf)

  end subroutine filltab_nudge_o3

!===============================================================================

  subroutine fill_jnudge()

    implicit none

    integer :: iwnud, iloop, jend

! Allocate and zero-fill jtab%jend()

    do iloop = 1,nloops_wnud
       jtab_wnud(iloop)%jend = 0
    enddo

! Compute and store jtab%jend

    do iloop = 1,nloops_wnud
       jtab_wnud(iloop)%jend = 0
       do iwnud = 2,mwnud
          if (itab_wnud(iwnud)%loop(iloop)) then
             jtab_wnud(iloop)%jend = jtab_wnud(iloop)%jend + 1
          endif
       enddo
       jtab_wnud(iloop)%jend = max(1,jtab_wnud(iloop)%jend)
    enddo

! Allocate and zero-fill JTAB_WNUD%IWNUD

    do iloop = 1,nloops_wnud
       jend = jtab_wnud(iloop)%jend
       allocate (jtab_wnud(iloop)%iwnud(jend))
       jtab_wnud(iloop)%iwnud(1:jend) = 0
    enddo

! Initialize JTAB%JEND counters to zero

    do iloop = 1,nloops_wnud
       jtab_wnud(iloop)%jend = 0
    enddo

! Compute JTAB_WNUD%IWNUD

    do iwnud = 2,mwnud
       do iloop = 1,nloops_wnud
          if (itab_wnud(iwnud)%loop(iloop)) then
             jtab_wnud(iloop)%jend = jtab_wnud(iloop)%jend + 1
             jtab_wnud(iloop)%iwnud(jtab_wnud(iloop)%jend) = iwnud
          endif
       enddo
    enddo

  end subroutine fill_jnudge

End Module mem_nudge
