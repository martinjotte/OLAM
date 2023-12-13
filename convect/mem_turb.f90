Module mem_turb

  implicit none

  real,    allocatable :: tkep     (:,:)
  real,    allocatable :: epsp     (:,:)
  real,    allocatable :: vkm      (:,:)
  real,    allocatable :: vkh      (:,:)
  real,    allocatable :: agamma   (:,:)

  real,    allocatable :: akm_sfc  (:,:)
  real,    allocatable :: ustar_k  (:,:)
  real,    allocatable :: wtv0_k   (:,:)

  real,    allocatable :: sfluxt     (:)
  real,    allocatable :: sfluxr     (:)
  real,    allocatable :: ustar      (:)
  real,    allocatable :: wstar      (:)
  real,    allocatable :: moli       (:)
  real,    allocatable :: wtv0       (:)
  real,    allocatable :: pblh       (:)
  real,    allocatable :: vkm_sfc    (:)
  integer, allocatable :: kpblh      (:)

  real,    allocatable :: frac_land  (:)
  real,    allocatable :: frac_sea   (:)
  real,    allocatable :: frac_lake  (:)
  real,    allocatable :: frac_urb   (:)

  real,    allocatable :: frac_sfc (:,:)
  real,    allocatable :: frac_sfck(:,:)
  real,    allocatable :: arw_sfc  (:,:)

  real,    allocatable :: akmodx   (:,:)
  real,    allocatable :: akhodx   (:,:)

  integer, allocatable :: khtop      (:)
  integer, allocatable :: khtopv     (:)

  integer, allocatable :: kmtop      (:)
  integer, allocatable :: kmtopv     (:)

Contains

!===============================================================================

  subroutine alloc_turb(mza, mwa, mva, nsw_max, idiffk, mrls)

    use misc_coms,  only: rinit
    implicit none

    integer, intent(in) :: mza, mwa, mva, nsw_max, mrls, idiffk(mrls)

!   Allocate arrays based on options (if necessary)
!   Initialize arrays to zero

    allocate (akm_sfc  (nsw_max,mwa)) ; akm_sfc   = 0.0
    allocate (ustar_k  (nsw_max,mwa)) ; ustar_k   = 0.0
    allocate (wtv0_k   (nsw_max,mwa)) ; wtv0_k    = 0.0

    allocate (frac_sfc (nsw_max,mwa)) ; frac_sfc  = 0.0
    allocate (frac_sfck(nsw_max,mwa)) ; frac_sfck = 0.0
    allocate (arw_sfc  (nsw_max,mwa)) ; arw_sfc   = 0.0

    allocate (vkm(mza,mwa)) ; vkm = 0.0
    allocate (vkh(mza,mwa)) ; vkh = 0.0

    if (any(idiffk(1:mrls) == 1)) then
       allocate (agamma(mza,mwa)) ; agamma = 0.0
    endif

    allocate (sfluxt    (mwa)) ; sfluxt    = rinit
    allocate (sfluxr    (mwa)) ; sfluxr    = rinit
    allocate (ustar     (mwa)) ; ustar     = rinit
    allocate (wstar     (mwa)) ; wstar     = rinit
    allocate (moli      (mwa)) ; moli      = rinit
    allocate (wtv0      (mwa)) ; wtv0      = rinit
    allocate (pblh      (mwa)) ; pblh      = rinit
    allocate (vkm_sfc   (mwa)) ; vkm_sfc   = rinit
    allocate (kpblh     (mwa)) ; kpblh     = 1
    allocate (frac_urb  (mwa)) ; frac_urb  = 0.0
    allocate (frac_land (mwa)) ; frac_land = 0.0
    allocate (frac_lake (mwa)) ; frac_lake = 0.0
    allocate (frac_sea  (mwa)) ; frac_sea  = 0.0

    allocate (akmodx(mza,mva)) ; akmodx    = 0.0
    allocate (akhodx(mza,mva)) ; akhodx    = 0.0

    allocate (khtop (mwa)) ; khtop = 0
    allocate (kmtop (mwa)) ; kmtop = 0

    allocate (khtopv(mva)) ; khtopv = 0
    allocate (kmtopv(mva)) ; kmtopv = 0

  end subroutine alloc_turb

  !===============================================================================

  subroutine dealloc_turb()

    implicit none

    if (allocated(tkep))      deallocate (tkep)
    if (allocated(epsp))      deallocate (epsp)
    if (allocated(vkm))       deallocate (vkm)
    if (allocated(vkh))       deallocate (vkh)
    if (allocated(akm_sfc))   deallocate (akm_sfc)
    if (allocated(vkm_sfc))   deallocate (vkm_sfc)
    if (allocated(ustar_k))   deallocate (ustar_k)
    if (allocated(wtv0_k))    deallocate (wtv0_k)
    if (allocated(sfluxt))    deallocate (sfluxt)
    if (allocated(sfluxr))    deallocate (sfluxr)
    if (allocated(ustar))     deallocate (ustar)
    if (allocated(wstar))     deallocate (wstar)
    if (allocated(moli))      deallocate (moli)
    if (allocated(wtv0))      deallocate (wtv0)
    if (allocated(pblh))      deallocate (pblh)
    if (allocated(kpblh))     deallocate (kpblh)
    if (allocated(akmodx))    deallocate (akmodx)
    if (allocated(akhodx))    deallocate (akhodx)
    if (allocated(agamma))    deallocate (agamma)

  end subroutine dealloc_turb

!===============================================================================

  subroutine filltab_turb()

    use var_tables, only: increment_vtable
    implicit none

    if (allocated(tkep))     call increment_vtable('TKEP',    'AW', rvar2=tkep, mpt1=.true.)

    if (allocated(epsp))     call increment_vtable('EPSP',    'AW', rvar2=epsp, mpt1=.true.)

    if (allocated(vkm))      call increment_vtable('VKM',     'AW', rvar2=vkm)

    if (allocated(vkh))      call increment_vtable('VKH',     'AW', rvar2=vkh)

!   if (allocated(vkm_sfc))  call increment_vtable('VKM_SFC', 'AW', rvar2=vkm_sfc)

    if (allocated(sfluxt))   call increment_vtable('SFLUXT',  'AW', rvar1=sfluxt)

    if (allocated(sfluxr))   call increment_vtable('SFLUXR',  'AW', rvar1=sfluxr)

    if (allocated(ustar))    call increment_vtable('USTAR',   'AW', rvar1=ustar)

    if (allocated(wstar))    call increment_vtable('WSTAR',   'AW', rvar1=wstar)

    if (allocated(wtv0))     call increment_vtable('WTV0',    'AW', rvar1=wtv0)

    if (allocated(pblh))     call increment_vtable('PBLH',    'AW', rvar1=pblh)

    if (allocated(kpblh))    call increment_vtable('KPBLH',   'AW', ivar1=kpblh)

  end subroutine filltab_turb

End Module mem_turb
