Module mem_basic

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable :: vmc  (:,:) ! current V horiz momentum [kg/(m^2 s)]
  real, allocatable :: vc   (:,:) ! current V horiz velocity [m/s]
  real, allocatable :: wmc  (:,:) ! current vert momentum [kg/(m^2 s)]
  real, allocatable :: wc   (:,:) ! current vert velocity [m/s]
  real, allocatable :: rr_w (:,:) ! tot water mixing ratio [kg_wat/kg_dryair]
  real, allocatable :: rr_v (:,:) ! water vapor mixing ratio [kg_vap/kg_dryair]
  real, allocatable :: thil (:,:) ! ice-liquid pot temp [K]
  real, allocatable :: theta(:,:) ! pot temp [K]
  real, allocatable :: tair (:,:) ! temperature [K]

  real, allocatable, target :: vxe  (:,:) ! earth-relative x velocity at T point [m/s]
  real, allocatable, target :: vye  (:,:) ! earth-relative y velocity at T point [m/s]
  real, allocatable         :: vze  (:,:) ! earth-relative z velocity at T point [m/s]

  real, pointer, contiguous :: ue(:,:) ! easterly wind
  real, pointer, contiguous :: ve(:,:) ! northerly wind

  real(r8), allocatable :: press(:,:) ! air pressure [Pa]
  real(r8), allocatable :: rho  (:,:) ! dry air density [kg/m^3]

  ! Half-forrward advecting velocities for scalars averaged over long timestep
  real, allocatable :: wmasc(:,:)
  real, allocatable :: vmasc(:,:)

  real, allocatable :: alpha_press(:,:)
  real, allocatable :: pwfac      (:,:)
  real, allocatable :: pvfac      (:,:)

Contains

!===============================================================================

  subroutine alloc_basic(mza,mva,mwa)

    use misc_coms,   only: rinit, rinit8, mdomain
    use consts_coms, only: pc1, rdry, cpocv, r8
    use mem_grid,    only: zfacit, dniv, dzim

    implicit none

    integer, intent(in) :: mza,mva,mwa
    integer             :: iv, iw
    real,     parameter :: alpha_dry = pc1 * rdry**cpocv

!   Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
!   and initialize allocated arrays to zero

    allocate (vmc  (mza,mva))
    allocate (vc   (mza,mva))
    allocate (pvfac(mza,mva))
    allocate (vmasc(mza,mva))

    !$omp parallel do
    do iv = 1, mva
       vmc  (:,iv) = rinit
       vc   (:,iv) = rinit
       pvfac(:,iv) = zfacit(:) * dniv(iv)
       vmasc(:,iv) = 0.0
    enddo
    !$omp end parallel do

    allocate( rho        (mza,mwa) )
    allocate( press      (mza,mwa) )
    allocate( wmc        (mza,mwa) )
    allocate( wc         (mza,mwa) )
    allocate( thil       (mza,mwa) )
    allocate( theta      (mza,mwa) )
    allocate( tair       (mza,mwa) )
    allocate( rr_w       (mza,mwa) )
    allocate( rr_v       (mza,mwa) )
    allocate( vxe        (mza,mwa) )
    allocate( vye        (mza,mwa) )
    allocate( vze        (mza,mwa) )
    allocate( wmasc      (mza,mwa) )
    allocate( alpha_press(mza,mwa) )
    allocate( pwfac      (mza,mwa) )

    if (mdomain <= 1) then
       allocate( ue(mza,mwa) )
       allocate( ve(mza,mwa) )
    else
       ue => vxe
       ve => vye
    endif

    !$omp parallel do
    do iw = 1, mwa
       rho        (:,iw) = 0.0_r8
       press      (:,iw) = rinit8
       wmc        (:,iw) = rinit
       wc         (:,iw) = 0.0
       thil       (:,iw) = 0.0
       theta      (:,iw) = 0.0
       tair       (:,iw) = 0.0
       rr_w       (:,iw) = 0.0
       rr_v       (:,iw) = 0.0
       vxe        (:,iw) = 0.0
       vye        (:,iw) = 0.0
       vze        (:,iw) = 0.0
       wmasc      (:,iw) = 0.0
       alpha_press(:,iw) = alpha_dry
       pwfac      (:,iw) = dzim(:)

       if (mdomain <= 1) then
          ue      (:,iw) = rinit
          ve      (:,iw) = rinit
       endif

    enddo
    !$omp end parallel do

  end subroutine alloc_basic

!===============================================================================

  subroutine dealloc_basic()

    implicit none

    if (allocated(vmc))   deallocate (vmc)
    if (allocated(vc))    deallocate (vc)

    if (allocated(wmc))   deallocate (wmc)
    if (allocated(wc))    deallocate (wc)
    if (allocated(rho))   deallocate (rho)
    if (allocated(rr_w))  deallocate (rr_w)
    if (allocated(rr_v))  deallocate (rr_v)
    if (allocated(press)) deallocate (press)
    if (allocated(thil))  deallocate (thil)
    if (allocated(theta)) deallocate (theta)
    if (allocated(tair))  deallocate (tair)

    if (allocated(vxe))   deallocate (vxe)
    if (allocated(vye))   deallocate (vye)
    if (allocated(vze))   deallocate (vze)

    if (allocated(wmasc)) deallocate (wmasc)
    if (allocated(vmasc)) deallocate (vmasc)

  end subroutine dealloc_basic

!===============================================================================

  subroutine filltab_basic()

    use var_tables, only: increment_vtable
    use misc_coms,  only: mdomain
    use oname_coms, only: nl

    implicit none

    if (allocated(vmc))   call increment_vtable('VMC',  'AV', rvar2=vmc)

    if (allocated(vc))    call increment_vtable('VC',   'AV', rvar2=vc)

    if (allocated(wmc))   call increment_vtable('WMC',  'AW', rvar2=wmc)

    if (allocated(wc))    call increment_vtable('WC',   'AW', rvar2=wc)

    if (allocated(rr_w))  call increment_vtable('RR_W', 'AW', rvar2=rr_w,  mpt1=.true.)

    if (allocated(rr_v))  call increment_vtable('RR_V', 'AW', rvar2=rr_v,  mpt1=.true.)

    if (allocated(thil))  call increment_vtable('THIL', 'AW', rvar2=thil)

    if (allocated(theta)) call increment_vtable('THETA','AW', rvar2=theta, mpt1=.true.)

    if (allocated(tair))  call increment_vtable('TAIR', 'AW', rvar2=tair,  mpt1=.true.)

    if (allocated(rho))   call increment_vtable('RHO',  'AW', dvar2=rho)

    if (allocated(press)) call increment_vtable('PRESS','AW', dvar2=press)

    if (allocated(vxe))   call increment_vtable('VXE',  'AW', rvar2=vxe)

    if (allocated(vye))   call increment_vtable('VYE',  'AW', rvar2=vye)

    if (allocated(vze))   call increment_vtable('VZE',  'AW', rvar2=vze)

    if (nl%hist_write_ue_ve .and. mdomain <= 1) then

       if (associated(ue)) call increment_vtable('UE', 'AW', rvar2=ue, noread=.true.)
       if (associated(ve)) call increment_vtable('VE', 'AW', rvar2=ve, noread=.true.)

    endif

  end subroutine filltab_basic

End Module mem_basic
