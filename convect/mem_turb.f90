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

  real,    allocatable :: sfluxt     (:) ! canopy-to-atm sensible heat flux [W m^-2]
  real,    allocatable :: sfluxr     (:) ! canopy-to-atm water vapor flux [kg_vap m^-2 s^-1]
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
    integer             :: iw, iv

!   Allocate arrays based on options (if necessary)
!   Initialize arrays to zero

    allocate( akm_sfc  (nsw_max,mwa) )
    allocate( ustar_k  (nsw_max,mwa) )
    allocate( wtv0_k   (nsw_max,mwa) )
    allocate( frac_sfc (nsw_max,mwa) )
    allocate( frac_sfck(nsw_max,mwa) )
    allocate( arw_sfc  (nsw_max,mwa) )

    allocate( vkm(mza,mwa) )
    allocate( vkh(mza,mwa) )

    if (any(idiffk(1:mrls) == 1)) then
       allocate( agamma(mza,mwa) )
    endif

    allocate( sfluxt   (mwa) )
    allocate( sfluxr   (mwa) )
    allocate( ustar    (mwa) )
    allocate( wstar    (mwa) )
    allocate( moli     (mwa) )
    allocate( wtv0     (mwa) )
    allocate( pblh     (mwa) )
    allocate( vkm_sfc  (mwa) )
    allocate( kpblh    (mwa) )
    allocate( frac_urb (mwa) )
    allocate( frac_land(mwa) )
    allocate( frac_lake(mwa) )
    allocate( frac_sea (mwa) )

    allocate( khtop(mwa) )
    allocate( kmtop(mwa) )

    allocate( akmodx(mza,mva) )
    allocate( akhodx(mza,mva) )

    allocate( khtopv(mva) )
    allocate( kmtopv(mva) )

    !$omp parallel
    !$omp do
    do iw = 1, mwa
       if ( allocated( akm_sfc  ) ) akm_sfc  (:,iw) = 0.
       if ( allocated( ustar_k  ) ) ustar_k  (:,iw) = 0.
       if ( allocated( wtv0_k   ) ) wtv0_k   (:,iw) = 0.
       if ( allocated( frac_sfc ) ) frac_sfc (:,iw) = 0.
       if ( allocated( frac_sfck) ) frac_sfck(:,iw) = 0.
       if ( allocated( arw_sfc  ) ) arw_sfc  (:,iw) = 0.
       if ( allocated( vkm      ) ) vkm      (:,iw) = 0.
       if ( allocated( vkh      ) ) vkh      (:,iw) = 0.
       if ( allocated( agamma   ) ) agamma   (:,iw) = 0.
       if ( allocated( sfluxt   ) ) sfluxt     (iw) = rinit
       if ( allocated( sfluxr   ) ) sfluxr     (iw) = rinit
       if ( allocated( ustar    ) ) ustar      (iw) = rinit
       if ( allocated( wstar    ) ) wstar      (iw) = rinit
       if ( allocated( moli     ) ) moli       (iw) = rinit
       if ( allocated( wtv0     ) ) wtv0       (iw) = rinit
       if ( allocated( pblh     ) ) pblh       (iw) = rinit
       if ( allocated( vkm_sfc  ) ) vkm_sfc    (iw) = rinit
       if ( allocated( kpblh    ) ) kpblh      (iw) = 1
       if ( allocated( frac_urb ) ) frac_urb   (iw) = 0.0
       if ( allocated( frac_land) ) frac_land  (iw) = 0.0
       if ( allocated( frac_lake) ) frac_lake  (iw) = 0.0
       if ( allocated( frac_sea ) ) frac_sea   (iw) = 0.0
       if ( allocated( khtop    ) ) khtop      (iw) = 0
       if ( allocated( kmtop    ) ) kmtop      (iw) = 0
    enddo
    !$omp end do nowait

    !$omp do
    do iv = 1, mva
       if ( allocated ( akmodx ) ) akmodx(:,iv) = 0.0
       if ( allocated ( akhodx ) ) akhodx(:,iv) = 0.0
       if ( allocated ( khtopv ) ) khtopv  (iv) = 0
       if ( allocated ( kmtopv ) ) kmtopv  (iv) = 0
    enddo
    !$omp end do
    !$omp end parallel

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
