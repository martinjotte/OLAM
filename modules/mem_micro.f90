!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================

Module mem_micro

  use consts_coms, only: r8

  real, allocatable, target :: sh_c(:,:) ! cloud water spec dens [kg_cld/kg_air]
  real, allocatable, target :: sh_d(:,:) ! drizzle spec dens [kg_driz/kg_air]
  real, allocatable, target :: sh_r(:,:) ! rain spec dens [kg_rain/kg_air]
  real, allocatable, target :: sh_p(:,:) ! pristine ice spec dens [kg_pris/kg_air]
  real, allocatable, target :: sh_s(:,:) ! snow spec dens [kg_snow/kg_air]
  real, allocatable, target :: sh_a(:,:) ! aggregates spec dens [kg_agg/kg_air]
  real, allocatable, target :: sh_g(:,:) ! graupel spec dens [kg_graup/kg_air]
  real, allocatable, target :: sh_h(:,:) ! hail spec dens [kg_hail/kg_air]

  real, allocatable, target :: con_c(:,:) ! cloud drop number conc [#_cld/kg_air]
  real, allocatable, target :: con_d(:,:) ! drizzle number conc [#_driz/kg_air]
  real, allocatable, target :: con_r(:,:) ! rain number conc [#_rain/kg_air]
  real, allocatable, target :: con_p(:,:) ! pristine ice number conc [#_pris/kg_air]
  real, allocatable, target :: con_s(:,:) ! snow number conc [#_snow/kg_air]
  real, allocatable, target :: con_a(:,:) ! aggregates number conc [#_aggr/kg_air]
  real, allocatable, target :: con_g(:,:) ! graupel number conc [#_graup/kg_air]
  real, allocatable, target :: con_h(:,:) ! hail number conc [#_hail/kg_air]

  real, allocatable, target :: q2(:,:) ! rain internal energy [J/kg]
  real, allocatable, target :: q6(:,:) ! graupel internal energy [J/kg]
  real, allocatable, target :: q7(:,:) ! hail internal energy [J/kg]

  real(r8), allocatable, target :: accpd (:) ! sfc drizzle accum [kg/m^2]
  real(r8), allocatable, target :: accpr (:) ! sfc rain accum [kg/m^2]
  real(r8), allocatable, target :: accpp (:) ! sfc pristine ice accum [kg/m^2]
  real(r8), allocatable, target :: accps (:) ! sfc snow accum [kg/m^2]
  real(r8), allocatable, target :: accpa (:) ! sfc aggregates  accum [kg/m^2]
  real(r8), allocatable, target :: accpg (:) ! sfc graupel accum [kg/m^2]
  real(r8), allocatable, target :: accph (:) ! sfc hail accum [kg/m^2]

  real, allocatable, target :: pcprd (:) ! sfc drizzle pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprr (:) ! sfc rain pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprp (:) ! sfc pristine ice pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprs (:) ! sfc snow pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcpra (:) ! sfc aggregates pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprg (:) ! sfc graupel pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprh (:) ! sfc hail pcp rate [kg/(m^2 s)]

  real, allocatable, target :: cldnum(:) ! cloud drop number conc (climatology) [#_cld/kg_air]

  real, allocatable, target :: con_gccn(:,:)   ! GCCN number conc [#_gccn/kg_air]
  real, allocatable, target :: con_ifn (:,:)   ! IFN  number conc [#_ifn/kg_air]

  Type ccntyp_vars   
     real, allocatable :: con_ccn (:,:) ! CCN  number conc [#_ccn/kg_air]
     real, allocatable :: con_ccnt(:,:) ! CCN  number tendency [#_ccn/(m^3 s)]
  End Type ccntyp_vars
   
  type (ccntyp_vars), allocatable, target :: ccntyp(:)

Contains

!===============================================================================

  subroutine alloc_micro(mza,mwa,miclevel,ncat,nccntyp,iccn,igccn,iifn,jnmb)

    use misc_coms, only: rinit
    use oname_coms, only: nl
    use consts_coms, only: r8

    implicit none

    integer, intent(in) :: mza,mwa,miclevel,ncat,nccntyp,iccn,igccn,iifn
    integer, intent(in) :: jnmb(ncat)
    integer :: ic

    ! Allocate arrays based on options
    ! Initialize arrays to zero

    allocate(cldnum(mwa)); cldnum = 0.

    if (miclevel >= 2) then
       allocate (sh_c(mza,mwa)) ; sh_c = rinit
    endif

    if (miclevel >= 3) then

       if (jnmb(1) == 5) then
          allocate (con_c(mza,mwa)) ; con_c = rinit
       endif

       if (jnmb(3) == 5) then
          allocate (con_p(mza,mwa)) ; con_p = rinit
          allocate (sh_p(mza,mwa)) ; sh_p  = rinit
          allocate (accpp   (mwa)) ; accpp = 0.0_r8
          allocate (pcprp   (mwa)) ; pcprp = rinit
       endif

       if (jnmb(8) == 5) then
          allocate (con_d(mza,mwa)) ; con_d = rinit
          allocate (sh_d(mza,mwa)) ; sh_d  = rinit
          allocate (accpd   (mwa)) ; accpd = 0.0_r8
          allocate (pcprd   (mwa)) ; pcprd = rinit
       endif

       if (jnmb(2) >= 1)  then
          if (jnmb(2) == 5) then
             allocate (con_r(mza,mwa)) ; con_r = rinit
          endif
          allocate (sh_r(mza,mwa)) ; sh_r  = rinit
          allocate (q2  (mza,mwa)) ; q2    = rinit
          allocate (accpr   (mwa)) ; accpr = 0.0_r8
          allocate (pcprr   (mwa)) ; pcprr = rinit
       endif

       if (jnmb(4) >= 1)  then
          if (jnmb(4) == 5) then
             allocate (con_s(mza,mwa)) ; con_s = rinit
          endif
          allocate (sh_s(mza,mwa)) ; sh_s  = rinit
          allocate (accps   (mwa)) ; accps = 0.0_r8
          allocate (pcprs   (mwa)) ; pcprs = rinit
       endif

       if (jnmb(5) >= 1)  then
          if (jnmb(5) == 5) then
             allocate (con_a(mza,mwa)) ; con_a = rinit
          endif
          allocate (sh_a(mza,mwa)) ; sh_a  = rinit
          allocate (accpa   (mwa)) ; accpa = 0.0_r8
          allocate (pcpra   (mwa)) ; pcpra = rinit
       endif

       if (jnmb(6) >= 1) then
          if (jnmb(6) == 5) then
             allocate (con_g(mza,mwa)) ; con_g = rinit
          endif
          allocate (sh_g(mza,mwa)) ; sh_g  = rinit
          allocate (q6  (mza,mwa)) ; q6    = rinit
          allocate (accpg   (mwa)) ; accpg = 0.0_r8
          allocate (pcprg   (mwa)) ; pcprg = rinit
       endif

       if (jnmb(7) >= 1) then
          if (jnmb(7) == 5) then
             allocate (con_h(mza,mwa)) ; con_h = rinit
          endif
          allocate (sh_h(mza,mwa)) ; sh_h  = rinit
          allocate (q7  (mza,mwa)) ; q7    = rinit
          allocate (accph   (mwa)) ; accph = 0.0_r8
          allocate (pcprh   (mwa)) ; pcprh = rinit
       endif

       allocate (ccntyp(nccntyp))

       if (iccn >= 2) then
          do ic = 1,nccntyp
             allocate (ccntyp(ic)%con_ccn(mza,mwa)) ; ccntyp(ic)%con_ccn = rinit
          enddo
       endif

       if (igccn == 2) then
          allocate (con_gccn(mza,mwa)) ; con_gccn = rinit
       endif

       if (iifn == 2) then
          allocate (con_ifn(mza,mwa)) ; con_ifn = rinit
       endif

    endif

    if (nl%test_case > 10) then ! Encompasses DCMIP & parcel cases (March/2016)

       if (.not. allocated(sh_c))   allocate (sh_c(mza,mwa)) ; sh_c  = rinit
       if (.not. allocated(sh_r))   allocate (sh_r(mza,mwa)) ; sh_r  = rinit
       if (.not. allocated(q2))     allocate (q2  (mza,mwa)) ; q2    = rinit
       if (.not. allocated(accpr))  allocate (accpr   (mwa)) ; accpr = 0.0_r8
       if (.not. allocated(pcprr))  allocate (pcprr   (mwa)) ; pcprr = rinit

    endif

  end subroutine alloc_micro

!===============================================================================

  subroutine dealloc_micro(nccntyp)

    implicit none

    integer, intent(in) :: nccntyp
    integer :: ic

    if (allocated(sh_c))    deallocate (sh_c)
    if (allocated(sh_d))    deallocate (sh_d)
    if (allocated(sh_r))    deallocate (sh_r)
    if (allocated(sh_p))    deallocate (sh_p)
    if (allocated(sh_s))    deallocate (sh_s)
    if (allocated(sh_a))    deallocate (sh_a)
    if (allocated(sh_g))    deallocate (sh_g)
    if (allocated(sh_h))    deallocate (sh_h)

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

    if (allocated(sh_c)) call increment_vtable('SH_C', 'AW', mpt1=.true., rvar2=sh_c)

    if (allocated(sh_d)) call increment_vtable('SH_D', 'AW', mpt1=.true., rvar2=sh_d)

    if (allocated(sh_r)) call increment_vtable('SH_R', 'AW', mpt1=.true., rvar2=sh_r)

    if (allocated(sh_p)) call increment_vtable('SH_P', 'AW', mpt1=.true., rvar2=sh_p)

    if (allocated(sh_s)) call increment_vtable('SH_S', 'AW', mpt1=.true., rvar2=sh_s)

    if (allocated(sh_a)) call increment_vtable('SH_A', 'AW', mpt1=.true., rvar2=sh_a)

    if (allocated(sh_g)) call increment_vtable('SH_G', 'AW', mpt1=.true., rvar2=sh_g)

    if (allocated(sh_h)) call increment_vtable('SH_H', 'AW', mpt1=.true., rvar2=sh_h)

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

    do ic = 1,nccntyp
       if (allocated (ccntyp(ic)%con_ccn)) then
          write(sname,'(a7,i3.3)') 'CON_CCN', ic
          call increment_vtable(sname, 'AW', mpt1=.true., rvar2=ccntyp(ic)%con_ccn)
       endif
    enddo

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
