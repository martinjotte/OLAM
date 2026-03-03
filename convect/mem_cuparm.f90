Module mem_cuparm

   use consts_coms, only: r8

   real,    allocatable :: thsrc (:,:) ! heat / cp  tend ( rho X temp )
   real,    allocatable :: rtsrc (:,:) ! water mass tend ( rho X rr_w )
   real,    allocatable :: umsrc (:,:) ! rho * ue tend
   real,    allocatable :: vmsrc (:,:) ! rho * ve tend

   real,    allocatable :: cu_pwa(:,:)
   real,    allocatable :: cu_pev(:,:)

   real(r8),allocatable :: aconpr  (:)
   real,    allocatable :: conprr  (:)
   real,    allocatable :: qwcon (:,:) ! convective cloud water

   real,    allocatable :: cbmf    (:) ! updraft mass flux
   real,    allocatable :: cddf    (:) ! downdraft mass flux
   integer, allocatable :: kcutop  (:)
   integer, allocatable :: kcubot  (:)
   integer, allocatable :: kudbot  (:)
   integer, allocatable :: kddtop  (:)
   integer, allocatable :: kddmax  (:)
   integer, allocatable :: kddbot  (:)
   integer, allocatable :: iactcu  (:)
   integer, allocatable :: kstabi  (:)

   private :: r8

Contains

!===============================================================================

  subroutine alloc_cuparm(mza, mwa, mrls, nqparm)

    use consts_coms, only: r8
    use oname_coms,  only: nl

    implicit none

    integer, intent(in) :: mza, mwa, mrls
    integer, intent(in) :: nqparm(mrls)
    integer             :: iw

    allocate( iactcu(mwa) )

    if ( any(nqparm(1:mrls) > 0) ) then

       ! Base tendency arrays for all deep convective schemes

       allocate( thsrc(mza,mwa) )
       allocate( rtsrc(mza,mwa) )
       allocate( aconpr   (mwa) )
       allocate( conprr   (mwa) )
       allocate( qwcon(mza,mwa) )

       if ( any(nqparm(1:mrls) == 1) .or. &
            any(nqparm(1:mrls) == 2) .or. &
            any(nqparm(1:mrls) == 4) .or. &
            any(nqparm(1:mrls) == 5) ) then

          allocate( umsrc(mza,mwa) )
          allocate( vmsrc(mza,mwa) )

       endif

       ! Do we need to save the convective precipitation flux?

       if (nl%do_chem > 0) then
          allocate( cu_pwa(mza,mwa) )
          allocate( cu_pev(mza,mwa) )
       endif

       ! Diagnostic arrays for clouds/radiation/tracer mixing

       allocate( cbmf  (mwa) )
       allocate( cddf  (mwa) )
       allocate( kcutop(mwa) )
       allocate( kcubot(mwa) )
       allocate( kudbot(mwa) )
       allocate( kddtop(mwa) )
       allocate( kddmax(mwa) )
       allocate( kddbot(mwa) )
       allocate( kstabi(mwa) )

    endif

    !$omp parallel do
    do iw = 1, mwa
       if ( allocated( thsrc ) ) thsrc (:,iw) = 0.
       if ( allocated( rtsrc ) ) rtsrc (:,iw) = 0.
       if ( allocated( umsrc ) ) umsrc (:,iw) = 0.
       if ( allocated( vmsrc ) ) vmsrc (:,iw) = 0.
       if ( allocated( cu_pwa) ) cu_pwa(:,iw) = 0.
       if ( allocated( cu_pev) ) cu_pev(:,iw) = 0.
       if ( allocated( aconpr) ) aconpr  (iw) = 0._r8
       if ( allocated( conprr) ) conprr  (iw) = 0.
       if ( allocated( qwcon ) ) qwcon (:,iw) = 0.
       if ( allocated( cbmf  ) ) cbmf    (iw) = 0.
       if ( allocated( cddf  ) ) cddf    (iw) = 0.
       if ( allocated( kcutop) ) kcutop  (iw) = 0
       if ( allocated( kcubot) ) kcubot  (iw) = 0
       if ( allocated( kudbot) ) kudbot  (iw) = 0
       if ( allocated( kddtop) ) kddtop  (iw) = 0
       if ( allocated( kddmax) ) kddmax  (iw) = 0
       if ( allocated( kddbot) ) kddbot  (iw) = 0
       if ( allocated( iactcu) ) iactcu  (iw) = 0
       if ( allocated( kstabi) ) kstabi  (iw) = 0
    enddo
    !$omp end parallel do

  end subroutine alloc_cuparm

!===============================================================================

  subroutine dealloc_cuparm()
    implicit none

    if (allocated(thsrc))    deallocate (thsrc)
    if (allocated(rtsrc))    deallocate (rtsrc)
    if (allocated(umsrc))    deallocate (umsrc)
    if (allocated(vmsrc))    deallocate (vmsrc)
    if (allocated(aconpr))   deallocate (aconpr)
    if (allocated(conprr))   deallocate (conprr)
    if (allocated(qwcon))    deallocate (qwcon)
    if (allocated(cbmf))     deallocate (cbmf)
    if (allocated(cddf))     deallocate (cddf)
    if (allocated(kcutop))   deallocate (kcutop)
    if (allocated(kcubot))   deallocate (kcubot)
    if (allocated(kudbot))   deallocate (kudbot)
    if (allocated(kddtop))   deallocate (kddtop)
    if (allocated(kddmax))   deallocate (kddmax)
    if (allocated(kddbot))   deallocate (kddbot)
    if (allocated(iactcu))   deallocate (iactcu)
    if (allocated(cu_pwa))   deallocate (cu_pwa)
    if (allocated(cu_pev))   deallocate (cu_pev)

  end subroutine dealloc_cuparm

!===============================================================================

  subroutine filltab_cuparm(mrls, nqparm)

    use var_tables, only: increment_vtable

    implicit none

    integer, intent(in) :: mrls
    integer, intent(in) :: nqparm(mrls)

    if (allocated(thsrc))  call increment_vtable('THSRC', 'AW', rvar2=thsrc)

    if (allocated(rtsrc))  call increment_vtable('RTSRC', 'AW', rvar2=rtsrc)

    if (allocated(umsrc))  call increment_vtable('UMSRC', 'AW', rvar2=umsrc)

    if (allocated(vmsrc))  call increment_vtable('VMSRC', 'AW', rvar2=vmsrc)

    if (allocated(aconpr)) call increment_vtable('ACONPR','AW', dvar1=aconpr)

    if (allocated(conprr)) call increment_vtable('CONPRR','AW', rvar1=conprr)

    if (allocated(qwcon))  call increment_vtable('QWCON', 'AW', rvar2=qwcon)

    if (allocated(cbmf))   call increment_vtable('CBMF',  'AW', rvar1=cbmf)

    if (allocated(cddf))   call increment_vtable('CDDF',  'AW', rvar1=cddf)

    if (allocated(kcutop)) call increment_vtable('KCUTOP','AW', ivar1=kcutop)

    if (allocated(kcubot)) call increment_vtable('KCUBOT','AW', ivar1=kcubot)

    if (allocated(kudbot)) call increment_vtable('KUDBOT','AW', ivar1=kudbot)

    if (allocated(kddtop)) call increment_vtable('KDDTOP','AW', ivar1=kddtop)

    if (allocated(kddmax)) call increment_vtable('KDDMAX','AW', ivar1=kddmax)

    if (allocated(kddbot)) call increment_vtable('KDDBOT','AW', ivar1=kddbot)

    if (allocated(kstabi)) call increment_vtable('KSTABI','AW', ivar1=kstabi)

    if (allocated(cu_pwa)) call increment_vtable('CU_PWA','AW', rvar2=cu_pwa)

    if (allocated(cu_pev)) call increment_vtable('CU_PEV','AW', rvar2=cu_pev)

    if ( any(nqparm(1:mrls) > 0) ) then
       if (allocated(iactcu)) call increment_vtable('IACTCU','AW', ivar1=iactcu)
    endif

  end subroutine filltab_cuparm

End Module mem_cuparm
