Module mem_radiate

  integer :: jday ! Julian day

  real :: solfac  ! solar-constant coefficient for variable Earth-Sun dist
  real :: sunx    ! x-component of unit vector pointing to sun [m]
  real :: suny    ! y-component of unit vector pointing to sun [m]
  real :: sunz    ! z-component of unit vector pointing to sun [m]

  ! Note: For better energy conservation, the longwave and shortwave
  ! heating rates have units of energy (divided by the constant Cp)
  ! which equals density times temperature. This is then converted to
  ! a potential temperature tendency when applied each sub-timestep.
  real, allocatable :: fthrd_sw    (:,:)
  real, allocatable :: fthrd_lw    (:,:)

  real, allocatable :: dlong    (:,:)  ! Downward atm longwave flux at W level (W/m^2)
  real, allocatable :: ulong    (:,:)  ! Upward atm longwave flux at W level (W/m^2)
  real, allocatable :: ulong_sfc(:,:)  ! Upward longwave flux at (shaved cell) surface(s) (W/m^2)

  real, allocatable :: Dulong_DT    (:,:) ! Upward atm longwave flux derivative w.r.t. skin temperature (W/m^2/K)
  real, allocatable :: Dulong_sfc_DT(:,:) ! Upward longwave flux derivative at surface w.r.t. skin temperatrue (W/m^2/K)
  real, allocatable :: Tskin_rad    (:,:) ! Skin temperature at time of radiative transfer call

  real, allocatable :: cloud_frac  (:,:)
  real, allocatable :: rshort        (:)
  real, allocatable :: rlong         (:)
  real, allocatable :: rlongup       (:)
  real, allocatable :: rshort_top    (:)
  real, allocatable :: rshortup_top  (:)
  real, allocatable :: rlongup_top   (:)
  real, allocatable :: albedt        (:)
  real, allocatable :: albedt_beam   (:)
  real, allocatable :: albedt_diffuse(:)
  real, allocatable :: cosz          (:) ! current solar zenith angle
  real, allocatable :: cosz_rad      (:) ! solar zenith angle at time of radiation update
  real, allocatable :: rlong_albedo  (:)
  real, allocatable :: rshort_diffuse(:)
  real, allocatable :: par           (:)
  real, allocatable :: par_diffuse   (:)
  real, allocatable :: ppfd          (:)
  real, allocatable :: ppfd_diffuse  (:)
  real, allocatable :: uva           (:)
  real, allocatable :: uvb           (:)
  real, allocatable :: uvc           (:)
  real, allocatable :: pbl_cld_forc  (:)

  ! arrays for transferring radiation quantities
  ! to surface with shaved cells

  real, allocatable :: rshort_ks        (:,:)
  real, allocatable :: rlong_ks         (:,:)
  real, allocatable :: rshort_diffuse_ks(:,:)
! real, allocatable :: par_ks           (:,:)
! real, allocatable :: par_diffuse_ks   (:,:)
  real, allocatable :: ppfd_ks          (:,:)
  real, allocatable :: ppfd_diffuse_ks  (:,:)

  ! clear-sky values

! real, allocatable :: rshort_clr      (:)
! real, allocatable :: rshortup_clr    (:)
! real, allocatable :: rshort_top_clr  (:)
! real, allocatable :: rshortup_top_clr(:)

! real, allocatable :: rlong_clr       (:)
! real, allocatable :: rlongup_clr     (:)
! real, allocatable :: rlongup_top_clr (:)

! RRTMg random cloud overlap seed

  integer, allocatable :: mcica_seed(:,:)

! These are used for adding extra levels at the top with the Mclatchy soundings

  integer, parameter :: maxadd_rad = 10    ! max allowed # of added rad levels
  integer            :: nadd_rad           ! actual # of added radiation levels
  real               :: zmrad      = 30.e3 ! top of radiation grid
  real,    parameter :: cosz_min   = 0.02  ! limit of cos(zenith angle) to call rad

Contains

!===============================================================================

  subroutine alloc_radiate(mza,mwa,nsw_max,ilwrtyp,iswrtyp)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza
    integer, intent(in) :: mwa
    integer, intent(in) :: nsw_max
    integer, intent(in) :: ilwrtyp
    integer, intent(in) :: iswrtyp
    integer             :: iw

    ! Always allocate 3D cloud fraction

    allocate (cloud_frac(mza,mwa))

    ! Allocate arrays based on options (if necessary)

    if (ilwrtyp + iswrtyp > 0)  then

       write(io6,*) 'allocating rad ', mwa, mza

       allocate (fthrd_sw  (mza,mwa))
       allocate (fthrd_lw  (mza,mwa))

       allocate (dlong        (mza,mwa))
       allocate (ulong        (mza,mwa))
       allocate (ulong_sfc(nsw_max,mwa))

       allocate (Dulong_DT        (mza,mwa))
       allocate (Dulong_sfc_DT(nsw_max,mwa))
       allocate (Tskin_rad    (nsw_max,mwa))

       allocate (rshort        (mwa))
       allocate (rlong         (mwa))
       allocate (rlongup       (mwa))
       allocate (rshort_top    (mwa))
       allocate (rshortup_top  (mwa))
       allocate (rlongup_top   (mwa))
       allocate (rshort_diffuse(mwa))

       allocate (par           (mwa))
       allocate (par_diffuse   (mwa))
       allocate (ppfd          (mwa))
       allocate (ppfd_diffuse  (mwa))
       allocate (uva           (mwa))
       allocate (uvb           (mwa))
       allocate (uvc           (mwa))
       allocate (pbl_cld_forc  (mwa))

       allocate (rlong_albedo  (mwa))
       allocate (albedt        (mwa))
       allocate (albedt_beam   (mwa))
       allocate (albedt_diffuse(mwa))
       allocate (cosz          (mwa))
       allocate (cosz_rad      (mwa))

!      allocate (rshort_clr      (mwa))
!      allocate (rshortup_clr    (mwa))
!      allocate (rshort_top_clr  (mwa))
!      allocate (rshortup_top_clr(mwa))
!      allocate (rlong_clr       (mwa))
!      allocate (rlongup_clr     (mwa))
!      allocate (rlongup_top_clr (mwa))

       allocate (mcica_seed(4,mwa))

       allocate (rshort_ks        (nsw_max,mwa))
       allocate (rshort_diffuse_ks(nsw_max,mwa))
       allocate (rlong_ks         (nsw_max,mwa))

!      allocate (par_ks           (nsw_max,mwa))
!      allocate (par_diffuse_ks   (nsw_max,mwa))
       allocate (ppfd_ks          (nsw_max,mwa))
       allocate (ppfd_diffuse_ks  (nsw_max,mwa))

    endif

    !$omp parallel do
    do iw = 1, mwa

       if ( allocated( cloud_frac    ) ) cloud_frac   (:,iw) = 0.0
       if ( allocated( fthrd_sw      ) ) fthrd_sw     (:,iw) = 0.0
       if ( allocated( fthrd_lw      ) ) fthrd_lw     (:,iw) = 0.0

       if ( allocated( dlong         ) ) dlong        (:,iw) = 0.0
       if ( allocated( ulong         ) ) ulong        (:,iw) = 0.0
       if ( allocated( ulong_sfc     ) ) ulong_sfc    (:,iw) = 0.0

       if ( allocated ( Dulong_DT    ) ) Dulong_DT    (:,iw) = 0.0
       if ( allocated ( Dulong_sfc_DT) ) Dulong_sfc_DT(:,iw) = 0.0
       if ( allocated ( Tskin_rad    ) ) Tskin_rad    (:,iw) = 0.0

       if ( allocated( rshort        ) ) rshort        (iw) = 0.0
       if ( allocated( rlong         ) ) rlong         (iw) = 0.0
       if ( allocated( rlongup       ) ) rlongup       (iw) = 0.0
       if ( allocated( rshort_top    ) ) rshort_top    (iw) = 0.0
       if ( allocated( rshortup_top  ) ) rshortup_top  (iw) = 0.0
       if ( allocated( rlongup_top   ) ) rlongup_top   (iw) = 0.0
       if ( allocated( rshort_diffuse) ) rshort_diffuse(iw) = 0.0

       if ( allocated( par           ) ) par           (iw) = 0.0
       if ( allocated( par_diffuse   ) ) par_diffuse   (iw) = 0.0
       if ( allocated( ppfd          ) ) ppfd          (iw) = 0.0
       if ( allocated( ppfd_diffuse  ) ) ppfd_diffuse  (iw) = 0.0
       if ( allocated( uva           ) ) uva           (iw) = 0.0
       if ( allocated( uvb           ) ) uvb           (iw) = 0.0
       if ( allocated( uvc           ) ) uvc           (iw) = 0.0
       if ( allocated( pbl_cld_forc  ) ) pbl_cld_forc  (iw) = 0.0

       if ( allocated( rlong_albedo  ) ) rlong_albedo  (iw) = rinit
       if ( allocated( albedt        ) ) albedt        (iw) = rinit
       if ( allocated( albedt_beam   ) ) albedt_beam   (iw) = rinit
       if ( allocated( albedt_diffuse) ) albedt_diffuse(iw) = rinit
       if ( allocated( cosz          ) ) cosz          (iw) = rinit
       if ( allocated( cosz_rad      ) ) cosz_rad      (iw) = rinit

!      if ( allocated( rshort_clr      ) ) rshort_clr      (iw) = 0.0
!      if ( allocated( rshortup_clr    ) ) rshortup_clr    (iw) = 0.0
!      if ( allocated( rshort_top_clr  ) ) rshort_top_clr  (iw) = 0.0
!      if ( allocated( rshortup_top_clr) ) rshortup_top_clr(iw) = 0.0
!      if ( allocated( rlong_clr       ) ) rlong_clr       (iw) = 0.0
!      if ( allocated( rlongup_clr     ) ) rlongup_clr     (iw) = 0.0
!      if ( allocated( rlongup_top_clr ) ) rlongup_top_clr (iw) = 0.0

       if ( allocated( mcica_seed       ) ) mcica_seed       (:,iw) = 0

       if ( allocated( rshort_ks        ) ) rshort_ks        (:,iw) = 0.0
       if ( allocated( rshort_diffuse_ks) ) rshort_diffuse_ks(:,iw) = 0.0
       if ( allocated( rlong_ks         ) ) rlong_ks         (:,iw) = 0.0

!      if ( allocated( par_ks           ) ) par_ks           (:,iw) = 0.0
!      if ( allocated( par_diffuse_ks   ) ) par_diffuse_ks   (:,iw) = 0.0
       if ( allocated( ppfd_ks          ) ) ppfd_ks          (:,iw) = 0.0
       if ( allocated( ppfd_diffuse_ks  ) ) ppfd_diffuse_ks  (:,iw) = 0.0

    enddo
    !$omp end parallel do

  end subroutine alloc_radiate

!===============================================================================

  subroutine dealloc_radiate()

    implicit none

    if (allocated(fthrd_sw))       deallocate (fthrd_sw)
    if (allocated(fthrd_lw))       deallocate (fthrd_lw)
    if (allocated(cloud_frac))     deallocate (cloud_frac)
    if (allocated(rshort))         deallocate (rshort)
    if (allocated(rlong))          deallocate (rlong)
    if (allocated(rlongup))        deallocate (rlongup)
    if (allocated(rshort_top))     deallocate (rshort_top)
    if (allocated(rshortup_top))   deallocate (rshortup_top)
    if (allocated(rlongup_top))    deallocate (rlongup_top)
    if (allocated(rshort_diffuse)) deallocate (rshort_diffuse)
    if (allocated(rlong_albedo))   deallocate (rlong_albedo)
    if (allocated(albedt))         deallocate (albedt)
    if (allocated(albedt_beam))    deallocate (albedt_beam)
    if (allocated(albedt_diffuse)) deallocate (albedt_diffuse)
    if (allocated(cosz))           deallocate (cosz)
    if (allocated(cosz_rad))       deallocate (cosz_rad)

!   if (allocated(rshort_clr))       deallocate (rshort_clr)
!   if (allocated(rshortup_clr))     deallocate (rshortup_clr)
!   if (allocated(rshort_top_clr))   deallocate (rshort_top_clr)
!   if (allocated(rshortup_top_clr)) deallocate (rshortup_top_clr)

!   if (allocated(rlong_clr))        deallocate (rlong_clr)
!   if (allocated(rlongup_clr))      deallocate (rlongup_clr)
!   if (allocated(rlongup_top_clr))  deallocate (rlongup_top_clr)

    if (allocated(par))              deallocate (par)
    if (allocated(par_diffuse))      deallocate (par_diffuse)
    if (allocated(ppfd))             deallocate (ppfd)
    if (allocated(ppfd_diffuse))     deallocate (ppfd_diffuse)
    if (allocated(uva))              deallocate (uva)
    if (allocated(uvb))              deallocate (uvb)
    if (allocated(uvc))              deallocate (uvc)
    if (allocated(pbl_cld_forc))     deallocate (pbl_cld_forc)

  end subroutine dealloc_radiate

!===============================================================================

  subroutine filltab_radiate()

    use var_tables, only: increment_vtable
    implicit none

    if (allocated(fthrd_sw))         call increment_vtable('FTHRD_SW',        'AW', rvar2=fthrd_sw)
    if (allocated(fthrd_lw))         call increment_vtable('FTHRD_LW',        'AW', rvar2=fthrd_lw)

    ! For incremental updates of longwave between RRTMg calls:
    if (allocated(dlong))            call increment_vtable('DLONG',           'AW', rvar2=dlong)
    if (allocated(ulong))            call increment_vtable('ULONG',           'AW', rvar2=ulong)
    if (allocated(ulong_sfc))        call increment_vtable('ULONG_SFC',       'AW', rvar2=ulong_sfc)
    if (allocated(Dulong_DT))        call increment_vtable('Dulong_DT',       'AW', rvar2=Dulong_DT)
    if (allocated(Dulong_sfc_DT))    call increment_vtable('Dulong_sfc_DT',   'AW', rvar2=Dulong_sfc_DT)
    if (allocated(Tskin_rad))        call increment_vtable('Tskin_rad',       'AW', rvar2=Tskin_rad)

    if (allocated(cloud_frac))       call increment_vtable('CLOUD_FRAC',      'AW', rvar2=cloud_frac)

    if (allocated(rshort))           call increment_vtable('RSHORT',          'AW', rvar1=rshort)

    if (allocated(rlong))            call increment_vtable('RLONG',           'AW', rvar1=rlong)

    if (allocated(rlongup))          call increment_vtable('RLONGUP',         'AW', rvar1=rlongup)

    if (allocated(rshort_top))       call increment_vtable('RSHORT_TOP',      'AW', rvar1=rshort_top)

    if (allocated(rshortup_top))     call increment_vtable('RSHORTUP_TOP',    'AW', rvar1=rshortup_top)

    if (allocated(rlongup_top))      call increment_vtable('RLONGUP_TOP',     'AW', rvar1=rlongup_top)

    if (allocated(albedt))           call increment_vtable('ALBEDT',          'AW', rvar1=albedt)

    if (allocated(cosz_rad))         call increment_vtable('COSZ_RAD',        'AW', rvar1=cosz_rad)

!   if (allocated(rshort_clr))       call increment_vtable('RSHORT_CLR',      'AW', rvar1=rshort_clr)

!   if (allocated(rshortup_clr))     call increment_vtable('RSHORTUP_CLR',    'AW', rvar1=rshortup_clr)

!   if (allocated(rshort_top_clr))   call increment_vtable('RSHORT_TOP_CLR',  'AW', rvar1=rshort_top_clr)

!   if (allocated(rshortup_top_clr)) call increment_vtable('RSHORTUP_TOP_CLR','AW', rvar1=rshortup_top_clr)

!   if (allocated(rlong_clr))        call increment_vtable('RLONG_CLR',       'AW', rvar1=rlong_clr)

!   if (allocated(rlongup_clr))      call increment_vtable('RLONGUP_CLR',     'AW', rvar1=rlongup_clr)

!   if (allocated(rlongup_top_clr))  call increment_vtable('RLONGUP_TOP_CLR', 'AW', rvar1=rlongup_top_clr)

    if (allocated(pbl_cld_forc))     call increment_vtable('PBL_CLD_FORC',    'AW', rvar1=pbl_cld_forc)

    if (allocated(par))              call increment_vtable('PAR_TOT',         'AW', rvar1=par)

    if (allocated(par_diffuse))      call increment_vtable('PAR_DIFFUSE',     'AW', rvar1=par_diffuse)

    if (allocated(ppfd))             call increment_vtable('PPFD_TOT',        'AW', rvar1=ppfd)

    if (allocated(ppfd_diffuse))     call increment_vtable('PPFD_DIFFUSE',    'AW', rvar1=ppfd_diffuse)

    if (allocated(mcica_seed))       call increment_vtable('MCICA_SEED',      'AW', ivar2=mcica_seed)

  end subroutine filltab_radiate

End Module mem_radiate
