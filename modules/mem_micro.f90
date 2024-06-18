Module mem_micro

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable :: rr_c(:,:) ! cloud water mixing ratio [kg_cld/kg_dryair]
  real, allocatable :: rr_d(:,:) ! drizzle mixing ratio [kg_driz/kg_dryair]
  real, allocatable :: rr_r(:,:) ! rain mixing ratio [kg_rain/kg_dryair]
  real, allocatable :: rr_p(:,:) ! pristine ice mixing ratio [kg_pris/kg_dryair]
  real, allocatable :: rr_s(:,:) ! snow mixing ratio [kg_snow/kg_dryair]
  real, allocatable :: rr_a(:,:) ! aggregates mixing ratio [kg_agg/kg_dryair]
  real, allocatable :: rr_g(:,:) ! graupel mixing ratio [kg_graup/kg_dryair]
  real, allocatable :: rr_h(:,:) ! hail mixing ratio [kg_hail/kg_dryair]

  real, allocatable :: con_c(:,:) ! cloud drop number conc [#_cld/kg_dryair]
  real, allocatable :: con_d(:,:) ! drizzle number conc [#_driz/kg_dryair]
  real, allocatable :: con_r(:,:) ! rain number conc [#_rain/kg_dryair]
  real, allocatable :: con_p(:,:) ! pristine ice number conc [#_pris/kg_dryair]
  real, allocatable :: con_s(:,:) ! snow number conc [#_snow/kg_dryair]
  real, allocatable :: con_a(:,:) ! aggregates number conc [#_aggr/kg_dryair]
  real, allocatable :: con_g(:,:) ! graupel number conc [#_graup/kg_dryair]
  real, allocatable :: con_h(:,:) ! hail number conc [#_hail/kg_dryair]

  real, allocatable :: q2(:,:) ! rain internal energy [J/kg_dryair]
  real, allocatable :: q6(:,:) ! graupel internal energy [J/kg_dryair]
  real, allocatable :: q7(:,:) ! hail internal energy [J/kg_dryair]

  real(r8), allocatable :: accpd (:) ! sfc drizzle accum [kg/m^2]
  real(r8), allocatable :: accpr (:) ! sfc rain accum [kg/m^2]
  real(r8), allocatable :: accpp (:) ! sfc pristine ice accum [kg/m^2]
  real(r8), allocatable :: accps (:) ! sfc snow accum [kg/m^2]
  real(r8), allocatable :: accpa (:) ! sfc aggregates  accum [kg/m^2]
  real(r8), allocatable :: accpg (:) ! sfc graupel accum [kg/m^2]
  real(r8), allocatable :: accph (:) ! sfc hail accum [kg/m^2]

  real, allocatable :: pcprd (:) ! sfc drizzle pcp rate [kg/(m^2 s)]
  real, allocatable :: pcprr (:) ! sfc rain pcp rate [kg/(m^2 s)]
  real, allocatable :: pcprp (:) ! sfc pristine ice pcp rate [kg/(m^2 s)]
  real, allocatable :: pcprs (:) ! sfc snow pcp rate [kg/(m^2 s)]
  real, allocatable :: pcpra (:) ! sfc aggregates pcp rate [kg/(m^2 s)]
  real, allocatable :: pcprg (:) ! sfc graupel pcp rate [kg/(m^2 s)]
  real, allocatable :: pcprh (:) ! sfc hail pcp rate [kg/(m^2 s)]

  real, allocatable :: cldnum(:) ! cloud drop number conc (climatology) [#_cld/kg_dryair]

  real, allocatable :: con_gccn(:,:)   ! GCCN number conc [#_gccn/kg_dryair]
  real, allocatable :: con_ifn (:,:)   ! IFN  number conc [#_ifn/kg_dryair]

  Type ccntyp_vars
     real, allocatable :: con_ccn (:,:) ! CCN  number conc [#_ccn/kg_dryair]
     real, allocatable :: con_ccnt(:,:) ! CCN  number tendency [#_ccn/(m^3 s)]
  End Type ccntyp_vars

  type (ccntyp_vars), allocatable :: ccntyp(:)

Contains

!===============================================================================

  subroutine alloc_micro(mza,mwa,miclevel,ncat,nccntyp,iccn,igccn,iifn,jnmb)

    use misc_coms, only: rinit
    use oname_coms, only: nl
    use consts_coms, only: r8

    implicit none

    integer, intent(in) :: mza,mwa,miclevel,ncat,nccntyp,iccn,igccn,iifn
    integer, intent(in) :: jnmb(ncat)
    integer             :: ic, iw

    ! Allocate arrays based on options
    ! Initialize arrays to zero

    allocate( cldnum(mwa) )

    if (miclevel >= 2) then
       allocate( rr_c(mza,mwa) )
    endif

    if (miclevel >= 3) then

       if (jnmb(1) == 5) then
          allocate( con_c(mza,mwa) )
       endif

       if (jnmb(3) == 5) then
          allocate( con_p(mza,mwa) )
          allocate( rr_p (mza,mwa) )
          allocate( accpp    (mwa) )
          allocate( pcprp    (mwa) )
       endif

       if (jnmb(8) == 5) then
          allocate( con_d(mza,mwa) )
          allocate( rr_d (mza,mwa) )
          allocate( accpd    (mwa) )
          allocate( pcprd    (mwa) )
       endif

       if (jnmb(2) >= 1)  then
          if (jnmb(2) == 5) then
             allocate( con_r(mza,mwa) )
          endif
          allocate( rr_r(mza,mwa) )
          allocate( q2  (mza,mwa) )
          allocate( accpr   (mwa) )
          allocate( pcprr   (mwa) )
       endif

       if (jnmb(4) >= 1)  then
          if (jnmb(4) == 5) then
             allocate( con_s(mza,mwa) )
          endif
          allocate( rr_s(mza,mwa) )
          allocate( accps   (mwa) )
          allocate( pcprs   (mwa) )
       endif

       if (jnmb(5) >= 1)  then
          if (jnmb(5) == 5) then
             allocate( con_a(mza,mwa) )
          endif
          allocate( rr_a(mza,mwa) )
          allocate( accpa   (mwa) )
          allocate( pcpra   (mwa) )
       endif

       if (jnmb(6) >= 1) then
          if (jnmb(6) == 5) then
             allocate( con_g(mza,mwa) )
          endif
          allocate( rr_g(mza,mwa) )
          allocate( q6  (mza,mwa) )
          allocate( accpg   (mwa) )
          allocate( pcprg   (mwa) )
       endif

       if (jnmb(7) >= 1) then
          if (jnmb(7) == 5) then
             allocate( con_h(mza,mwa) )
          endif
          allocate( rr_h(mza,mwa) )
          allocate( q7  (mza,mwa) )
          allocate( accph   (mwa) )
          allocate( pcprh   (mwa) )
       endif

       allocate( ccntyp(nccntyp) )

       if (iccn >= 2) then
          do ic = 1,nccntyp
             allocate( ccntyp(ic)%con_ccn(mza,mwa) )
          enddo
       endif

       if (igccn == 2) then
          allocate( con_gccn(mza,mwa) )
       endif

       if (iifn == 2) then
          allocate( con_ifn(mza,mwa) )
       endif

    endif

    if (nl%test_case > 13) then ! Encompasses DCMIP & parcel cases (March/2016)

       if (.not. allocated(rr_c))  allocate( rr_c(mza,mwa) )
       if (.not. allocated(rr_r))  allocate( rr_r(mza,mwa) )
       if (.not. allocated(q2))    allocate( q2  (mza,mwa) )
       if (.not. allocated(accpr)) allocate( accpr   (mwa) )
       if (.not. allocated(pcprr)) allocate( pcprr   (mwa) )

    endif

    !$omp parallel do private(ic)
    do iw = 1, mwa

       if ( allocated( rr_c ) ) rr_c(:,iw) = rinit
       if ( allocated( rr_d ) ) rr_d(:,iw) = rinit
       if ( allocated( rr_r ) ) rr_r(:,iw) = rinit
       if ( allocated( rr_p ) ) rr_p(:,iw) = rinit
       if ( allocated( rr_s ) ) rr_s(:,iw) = rinit
       if ( allocated( rr_a ) ) rr_a(:,iw) = rinit
       if ( allocated( rr_g ) ) rr_g(:,iw) = rinit
       if ( allocated( rr_h ) ) rr_h(:,iw) = rinit

       if ( allocated( con_c ) ) con_c(:,iw) = rinit
       if ( allocated( con_d ) ) con_d(:,iw) = rinit
       if ( allocated( con_r ) ) con_r(:,iw) = rinit
       if ( allocated( con_p ) ) con_p(:,iw) = rinit
       if ( allocated( con_s ) ) con_s(:,iw) = rinit
       if ( allocated( con_a ) ) con_a(:,iw) = rinit
       if ( allocated( con_g ) ) con_g(:,iw) = rinit
       if ( allocated( con_h ) ) con_h(:,iw) = rinit

       if ( allocated( q2 ) ) q2(:,iw) = rinit
       if ( allocated( q6 ) ) q6(:,iw) = rinit
       if ( allocated( q7 ) ) q7(:,iw) = rinit

       if ( allocated( accpd ) ) accpd(iw) = 0._r8
       if ( allocated( accpr ) ) accpr(iw) = 0._r8
       if ( allocated( accpp ) ) accpp(iw) = 0._r8
       if ( allocated( accps ) ) accps(iw) = 0._r8
       if ( allocated( accpa ) ) accpa(iw) = 0._r8
       if ( allocated( accpg ) ) accpg(iw) = 0._r8
       if ( allocated( accph ) ) accph(iw) = 0._r8

       if ( allocated( pcprd ) ) pcprd(iw) = 0.
       if ( allocated( pcprr ) ) pcprr(iw) = 0.
       if ( allocated( pcprp ) ) pcprp(iw) = 0.
       if ( allocated( pcprs ) ) pcprs(iw) = 0.
       if ( allocated( pcpra ) ) pcpra(iw) = 0.
       if ( allocated( pcprg ) ) pcprg(iw) = 0.
       if ( allocated( pcprh ) ) pcprh(iw) = 0.

       if ( allocated( cldnum ) ) cldnum(iw) = 0.

       if ( allocated( con_gccn ) ) con_gccn(:,iw) = rinit
       if ( allocated( con_ifn  ) ) con_ifn (:,iw) = rinit

       if ( allocated( ccntyp ) ) then
          do ic = 1, nccntyp
             if ( allocated( ccntyp(ic)%con_ccn ) ) ccntyp(ic)%con_ccn(:,iw) = rinit
          enddo
       endif

    enddo
    !$omp end parallel do

  end subroutine alloc_micro

!===============================================================================

  subroutine dealloc_micro(nccntyp)

    implicit none

    integer, intent(in) :: nccntyp
    integer             :: ic

    if (allocated(rr_c))    deallocate (rr_c)
    if (allocated(rr_d))    deallocate (rr_d)
    if (allocated(rr_r))    deallocate (rr_r)
    if (allocated(rr_p))    deallocate (rr_p)
    if (allocated(rr_s))    deallocate (rr_s)
    if (allocated(rr_a))    deallocate (rr_a)
    if (allocated(rr_g))    deallocate (rr_g)
    if (allocated(rr_h))    deallocate (rr_h)

    if (allocated(con_c))   deallocate (con_c)
    if (allocated(con_d))   deallocate (con_d)
    if (allocated(con_r))   deallocate (con_r)
    if (allocated(con_p))   deallocate (con_p)
    if (allocated(con_s))   deallocate (con_s)
    if (allocated(con_a))   deallocate (con_a)
    if (allocated(con_g))   deallocate (con_g)
    if (allocated(con_h))   deallocate (con_h)

    if (allocated(con_gccn))deallocate (con_gccn)
    if (allocated(con_ifn)) deallocate (con_ifn)

    if (allocated(ccntyp)) then
       do ic = 1,nccntyp
          if (allocated(ccntyp(ic)%con_ccn)) deallocate (ccntyp(ic)%con_ccn)
       enddo
       deallocate(ccntyp)
    endif

    if (allocated(q2))      deallocate (q2)
    if (allocated(q6))      deallocate (q6)
    if (allocated(q7))      deallocate (q7)

    if (allocated(accpd))   deallocate (accpd)
    if (allocated(accpr))   deallocate (accpr)
    if (allocated(accpp))   deallocate (accpp)
    if (allocated(accps))   deallocate (accps)
    if (allocated(accpa))   deallocate (accpa)
    if (allocated(accpg))   deallocate (accpg)
    if (allocated(accph))   deallocate (accph)

    if (allocated(pcprd))   deallocate (pcprd)
    if (allocated(pcprr))   deallocate (pcprr)
    if (allocated(pcprp))   deallocate (pcprp)
    if (allocated(pcprs))   deallocate (pcprs)
    if (allocated(pcpra))   deallocate (pcpra)
    if (allocated(pcprg))   deallocate (pcprg)
    if (allocated(pcprh))   deallocate (pcprh)

  end subroutine dealloc_micro

!===============================================================================

  subroutine filltab_micro(nccntyp)

    use var_tables, only: increment_vtable

    implicit none

    integer, intent(in) :: nccntyp

    integer :: ic
    character(10) :: sname

    if (allocated(rr_c)) call increment_vtable('RR_C', 'AW', mpt1=.true., rvar2=rr_c)

    if (allocated(rr_d)) call increment_vtable('RR_D', 'AW', mpt1=.true., rvar2=rr_d)

    if (allocated(rr_r)) call increment_vtable('RR_R', 'AW', mpt1=.true., rvar2=rr_r)

    if (allocated(rr_p)) call increment_vtable('RR_P', 'AW', mpt1=.true., rvar2=rr_p)

    if (allocated(rr_s)) call increment_vtable('RR_S', 'AW', mpt1=.true., rvar2=rr_s)

    if (allocated(rr_a)) call increment_vtable('RR_A', 'AW', mpt1=.true., rvar2=rr_a)

    if (allocated(rr_g)) call increment_vtable('RR_G', 'AW', mpt1=.true., rvar2=rr_g)

    if (allocated(rr_h)) call increment_vtable('RR_H', 'AW', mpt1=.true., rvar2=rr_h)

    if (allocated(con_c)) call increment_vtable('CON_C', 'AW', mpt1=.true., rvar2=con_c)

    if (allocated(con_d)) call increment_vtable('CON_D', 'AW', mpt1=.true., rvar2=con_d)

    if (allocated(con_r)) call increment_vtable('CON_R', 'AW', mpt1=.true., rvar2=con_r)

    if (allocated(con_p)) call increment_vtable('CON_P', 'AW', mpt1=.true., rvar2=con_p)

    if (allocated(con_s)) call increment_vtable('CON_S', 'AW', mpt1=.true., rvar2=con_s)

    if (allocated(con_a)) call increment_vtable('CON_A', 'AW', mpt1=.true., rvar2=con_a)

    if (allocated(con_g)) call increment_vtable('CON_G', 'AW', mpt1=.true., rvar2=con_g)

    if (allocated(con_h)) call increment_vtable('CON_H', 'AW', mpt1=.true., rvar2=con_h)

    if (allocated(con_gccn)) call increment_vtable('CON_GCCN', 'AW', mpt1=.true., rvar2=con_gccn)

    if (allocated(con_ifn)) call increment_vtable('CON_IFN', 'AW', mpt1=.true., rvar2=con_ifn)

    if (allocated(ccntyp)) then
       do ic = 1,nccntyp
          if (allocated (ccntyp(ic)%con_ccn)) then
             write(sname,'(a7,i3.3)') 'CON_CCN', ic
             call increment_vtable(sname, 'AW', mpt1=.true., rvar2=ccntyp(ic)%con_ccn)
          endif
       enddo
    endif

    if (allocated(q2)) call increment_vtable('Q2', 'AW', mpt1=.true., rvar2=q2)

    if (allocated(q6)) call increment_vtable('Q6', 'AW', mpt1=.true., rvar2=q6)

    if (allocated(q7)) call increment_vtable('Q7', 'AW', mpt1=.true., rvar2=q7)

    if (allocated(accpd)) call increment_vtable('ACCPD', 'AW', dvar1=accpd)

    if (allocated(accpr)) call increment_vtable('ACCPR', 'AW', dvar1=accpr)

    if (allocated(accpp)) call increment_vtable('ACCPP', 'AW', dvar1=accpp)

    if (allocated(accps)) call increment_vtable('ACCPS', 'AW', dvar1=accps)

    if (allocated(accpa)) call increment_vtable('ACCPA', 'AW', dvar1=accpa)

    if (allocated(accpg)) call increment_vtable('ACCPG', 'AW', dvar1=accpg)

    if (allocated(accph)) call increment_vtable('ACCPH', 'AW', dvar1=accph)

    if (allocated(pcprd)) call increment_vtable('PCPRD', 'AW', rvar1=pcprd)

    if (allocated(pcprr)) call increment_vtable('PCPRR', 'AW', rvar1=pcprr)

    if (allocated(pcprp)) call increment_vtable('PCPRP', 'AW', rvar1=pcprp)

    if (allocated(pcprs)) call increment_vtable('PCPRS', 'AW', rvar1=pcprs)

    if (allocated(pcpra)) call increment_vtable('PCPRA', 'AW', rvar1=pcpra)

    if (allocated(pcprg)) call increment_vtable('PCPRG', 'AW', rvar1=pcprg)

    if (allocated(pcprh)) call increment_vtable('PCPRH', 'AW', rvar1=pcprh)

  end subroutine filltab_micro

End Module mem_micro
