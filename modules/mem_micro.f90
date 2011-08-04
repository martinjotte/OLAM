!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

Module mem_micro

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

  real, allocatable, target :: con_ccn (:,:) ! CCN number conc [#_ccn/kg_air]
  real, allocatable, target :: con_gccn(:,:) ! GCCN number conc [#_gccn/kg_air]
  real, allocatable, target :: con_ifn (:,:) ! IFN number conc [#_ifn/kg_air]

  real, allocatable, target :: q2(:,:) ! rain internal energy [J/kg]
  real, allocatable, target :: q6(:,:) ! graupel internal energy [J/kg]
  real, allocatable, target :: q7(:,:) ! hail internal energy [J/kg]

  real, allocatable, target :: accpd (:) ! sfc drizzle accum [kg/m^2]
  real, allocatable, target :: accpr (:) ! sfc rain accum [kg/m^2]
  real, allocatable, target :: accpp (:) ! sfc pristine ice accum [kg/m^2]
  real, allocatable, target :: accps (:) ! sfc snow accum [kg/m^2]
  real, allocatable, target :: accpa (:) ! sfc aggregates  accum [kg/m^2]
  real, allocatable, target :: accpg (:) ! sfc graupel accum [kg/m^2]
  real, allocatable, target :: accph (:) ! sfc hail accum [kg/m^2]

  real, allocatable, target :: pcprd (:) ! sfc drizzle pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprr (:) ! sfc rain pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprp (:) ! sfc pristine ice pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprs (:) ! sfc snow pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcpra (:) ! sfc aggregates pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprg (:) ! sfc graupel pcp rate [kg/(m^2 s)]
  real, allocatable, target :: pcprh (:) ! sfc hail pcp rate [kg/(m^2 s)]

  real, allocatable, target :: pcpgr (:) ! sfc total pcp RATE [kg/(m^2 s)]
  real, allocatable, target :: qpcpgr(:) ! sfc total pcp energy FLUX [J/(m^2 s)]
  real, allocatable, target :: dpcpgr(:) ! sfc total pcp depth accum RATE [m/s]

Contains

!===============================================================================

  subroutine alloc_micro(mza,mwa,level,ncat, &
       icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail)

    use misc_coms, only: rinit
    implicit none

    integer, intent(in) :: mza,mwa,level,ncat
    integer, intent(in) :: icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail

    ! Allocate arrays based on options
    ! Initialize arrays to zero

    if (level >= 2) then
       allocate (sh_c(mza,mwa)) ; sh_c = rinit
    endif

    if (level >= 3) then

       if (idriz >= 1)  then
          allocate (sh_d(mza,mwa)) ; sh_d  = rinit
          allocate (accpd   (mwa)) ; accpd = rinit
          allocate (pcprd   (mwa)) ; pcprd = rinit
       endif

       if (irain >= 1)  then
          allocate (sh_r(mza,mwa)) ; sh_r  = rinit
          allocate (q2  (mza,mwa)) ; q2    = rinit
          allocate (accpr   (mwa)) ; accpr = rinit
          allocate (pcprr   (mwa)) ; pcprr = rinit
       endif

       if (ipris >= 1)  then
          allocate (sh_p(mza,mwa)) ; sh_p  = rinit
          allocate (accpp   (mwa)) ; accpp = rinit
          allocate (pcprp   (mwa)) ; pcprp = rinit
       endif

       if (isnow >= 1)  then
          allocate (sh_s(mza,mwa)) ; sh_s  = rinit
          allocate (accps   (mwa)) ; accps = rinit
          allocate (pcprs   (mwa)) ; pcprs = rinit
       endif

       if (iaggr >= 1)  then
          allocate (sh_a(mza,mwa)) ; sh_a  = rinit
          allocate (accpa   (mwa)) ; accpa = rinit
          allocate (pcpra   (mwa)) ; pcpra = rinit
       endif

       if (igraup >= 1) then
          allocate (sh_g(mza,mwa)) ; sh_g  = rinit
          allocate (q6  (mza,mwa)) ; q6    = rinit
          allocate (accpg   (mwa)) ; accpg = rinit
          allocate (pcprg   (mwa)) ; pcprg = rinit
       endif

       if (ihail >= 1) then
          allocate (sh_h(mza,mwa)) ; sh_h  = rinit
          allocate (q7  (mza,mwa)) ; q7    = rinit
          allocate (accph   (mwa)) ; accph = rinit
          allocate (pcprh   (mwa)) ; pcprh = rinit
       endif

       if (icloud >= 5) then
          allocate (con_c(mza,mwa)) ; con_c = rinit
       endif

       if (idriz >= 5) then
          allocate (con_d(mza,mwa)) ; con_d = rinit
       endif

       if (irain == 5) then
          allocate (con_r(mza,mwa)) ; con_r = rinit
       endif

       if (ipris >= 5) then
          allocate (con_p(mza,mwa)) ; con_p = rinit
       endif

       if (isnow == 5) then
          allocate (con_s(mza,mwa)) ; con_s = rinit
       endif

       if (iaggr == 5) then
          allocate (con_a(mza,mwa)) ; con_a = rinit
       endif

       if (igraup == 5) then
          allocate (con_g(mza,mwa)) ; con_g = rinit
       endif

       if (ihail == 5) then
          allocate (con_h(mza,mwa)) ; con_h = rinit
       endif

       if (icloud == 7) then
          allocate (con_ccn(mza,mwa)) ; con_ccn = rinit
       endif

       if (idriz == 7) then
          allocate (con_gccn(mza,mwa)) ; con_gccn = rinit
       endif

       if (ipris == 7) then
          allocate (con_ifn(mza,mwa)) ; con_ifn = rinit
       endif

       allocate (pcpgr (mwa)) ; pcpgr  = rinit
       allocate (qpcpgr(mwa)) ; qpcpgr = rinit
       allocate (dpcpgr(mwa)) ; dpcpgr = rinit

    endif

  end subroutine alloc_micro

!===============================================================================

  subroutine dealloc_micro()
    implicit none

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

    if (allocated(con_ccn)) deallocate (con_ccn)
    if (allocated(con_gccn))deallocate (con_gccn)
    if (allocated(con_ifn)) deallocate (con_ifn)

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

    if (allocated(pcpgr))   deallocate (pcpgr)
    if (allocated(qpcpgr))  deallocate (qpcpgr)
    if (allocated(dpcpgr))  deallocate (dpcpgr)

  end subroutine dealloc_micro

!===============================================================================

  subroutine filltab_micro()

    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(sh_c)) then
       call increment_vtable('SH_C', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_c
    endif

    if (allocated(sh_d)) then
       call increment_vtable('SH_D', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_d
    endif

    if (allocated(sh_r)) then
       call increment_vtable('SH_R', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_r
    endif

    if (allocated(sh_p)) then
       call increment_vtable('SH_P', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_p
    endif

    if (allocated(sh_s)) then
       call increment_vtable('SH_S', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_s
    endif

    if (allocated(sh_a)) then
       call increment_vtable('SH_A', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_a
    endif

    if (allocated(sh_g)) then
       call increment_vtable('SH_G', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_g
    endif

    if (allocated(sh_h)) then
       call increment_vtable('SH_H', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_h
    endif

    if (allocated(con_c)) then
       call increment_vtable('CON_C', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_c
    endif

    if (allocated(con_d)) then
       call increment_vtable('CON_D', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_d
    endif

    if (allocated(con_r)) then
       call increment_vtable('CON_R', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_r
    endif

    if (allocated(con_p)) then
       call increment_vtable('CON_P', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_p
    endif

    if (allocated(con_s)) then
       call increment_vtable('CON_S', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_s
    endif

    if (allocated(con_a)) then
       call increment_vtable('CON_A', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_a
    endif

    if (allocated(con_g)) then
       call increment_vtable('CON_G', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_g
    endif

    if (allocated(con_h)) then
       call increment_vtable('CON_H', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_h

    endif

    if (allocated(con_ccn)) then
       call increment_vtable('CON_CCN', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_ccn
    endif

    if (allocated(con_gccn)) then
       call increment_vtable('CON_GCCN', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_gccn
    endif

    if (allocated(con_ifn)) then
       call increment_vtable('CON_IFN', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => con_ifn
    endif

    if (allocated(q2)) then
       call increment_vtable('Q2', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => q2
    endif

    if (allocated(q6)) then
       call increment_vtable('Q6', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => q6
    endif

    if (allocated(q7)) then
       call increment_vtable('Q7', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => q7
    endif

    if (allocated(accpd)) then
       call increment_vtable('ACCPD', 'AW')
       vtab_r(num_var)%rvar1_p => accpd
    endif

    if (allocated(accpr)) then
       call increment_vtable('ACCPR', 'AW')
       vtab_r(num_var)%rvar1_p => accpr
    endif

    if (allocated(accpp)) then
       call increment_vtable('ACCPP', 'AW')
       vtab_r(num_var)%rvar1_p => accpp
    endif

    if (allocated(accps)) then
       call increment_vtable('ACCPS', 'AW')
       vtab_r(num_var)%rvar1_p => accps
    endif

    if (allocated(accpa)) then
       call increment_vtable('ACCPA', 'AW')
       vtab_r(num_var)%rvar1_p => accpa
    endif

    if (allocated(accpg)) then
       call increment_vtable('ACCPG', 'AW')
       vtab_r(num_var)%rvar1_p => accpg
    endif

    if (allocated(accph)) then
       call increment_vtable('ACCPH', 'AW')
       vtab_r(num_var)%rvar1_p => accph
    endif

    if (allocated(pcprd)) then
       call increment_vtable('PCPRD', 'AW')
       vtab_r(num_var)%rvar1_p => pcprd
    endif

    if (allocated(pcprr)) then
       call increment_vtable('PCPRR', 'AW')
       vtab_r(num_var)%rvar1_p => pcprr
    endif

    if (allocated(pcprp)) then
       call increment_vtable('PCPRP', 'AW')
       vtab_r(num_var)%rvar1_p => pcprp
    endif

    if (allocated(pcprs)) then
       call increment_vtable('PCPRS', 'AW')
       vtab_r(num_var)%rvar1_p => pcprs
    endif

    if (allocated(pcpra)) then
       call increment_vtable('PCPRA', 'AW')
       vtab_r(num_var)%rvar1_p => pcpra
    endif

    if (allocated(pcprg)) then
       call increment_vtable('PCPRG', 'AW')
       vtab_r(num_var)%rvar1_p => pcprg
    endif

    if (allocated(pcprh)) then
       call increment_vtable('PCPRH', 'AW')
       vtab_r(num_var)%rvar1_p => pcprh
    endif

    if (allocated(pcpgr)) then
       call increment_vtable('PCPGR',  'AW')
       vtab_r(num_var)%rvar1_p => pcpgr
    endif

    if (allocated(qpcpgr)) then
       call increment_vtable('QPCPGR', 'AW')
       vtab_r(num_var)%rvar1_p => qpcpgr
    endif

    if (allocated(dpcpgr)) then
       call increment_vtable('DPCPGR', 'AW')
       vtab_r(num_var)%rvar1_p => dpcpgr
    endif

  end subroutine filltab_micro

End Module mem_micro
