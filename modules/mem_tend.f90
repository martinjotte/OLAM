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
Module mem_tend

   real, allocatable :: vmxet(:,:) ! Earth-cartesian x momentum tend [kg/(m^2 s^2)]
   real, allocatable :: vmyet(:,:) ! Earth-cartesian y momentum tend [kg/(m^2 s^2)]
   real, allocatable :: vmzet(:,:) ! Earth-cartesian z momentum tend [kg/(m^2 s^2)]

   real, target, allocatable :: thilt   (:,:) ! (rho * thil) tend [kg_air K / (m^3 s)]
   real, target, allocatable :: sh_wt   (:,:) ! water mass tend [kg_wat/(m^3 s)]

   real, target, allocatable :: sh_ct   (:,:) ! cloud mass tend [kg_cld/(m^3 s)]
   real, target, allocatable :: sh_rt   (:,:) ! rain mass tend [kg_rain/(m^3 s)]
   real, target, allocatable :: sh_pt   (:,:) ! pristine ice mass tend [kg_pris/(m^3 s)]
   real, target, allocatable :: sh_st   (:,:) ! snow mass tend [kg_snow/(m^3 s)]
   real, target, allocatable :: sh_at   (:,:) ! aggregates mass tend [kg_agg/(m^3 s)]
   real, target, allocatable :: sh_gt   (:,:) ! graupel mass tend [kg_grp/(m^3 s)]
   real, target, allocatable :: sh_ht   (:,:) ! hail mass tend [kg_hail/(m^3 s)]
   real, target, allocatable :: sh_dt   (:,:) ! drizzle mass tend [kg_drz/(m^3 s)]

   real, target, allocatable :: con_ct  (:,:) ! cloud number tend [#/(m^3 s)]]
   real, target, allocatable :: con_rt  (:,:) ! rain number tend [#/(m^3 s)]]
   real, target, allocatable :: con_pt  (:,:) ! pristine ice number tend [#/(m^3 s)]]
   real, target, allocatable :: con_st  (:,:) ! snow number tend [#/(m^3 s)]]
   real, target, allocatable :: con_at  (:,:) ! aggregates number tend [#/(m^3 s)]]
   real, target, allocatable :: con_gt  (:,:) ! graupel number tend [#/(m^3 s)]]
   real, target, allocatable :: con_ht  (:,:) ! hail number tend [#/(m^3 s)]]
   real, target, allocatable :: con_dt  (:,:) ! drizzle number tend [#/(m^3 s)]]

   real, target, allocatable :: q2t     (:,:) ! rain internal energy tend [J/(kg s)]
   real, target, allocatable :: q6t     (:,:) ! graupel internal energy tend [J/(kg s)]
   real, target, allocatable :: q7t     (:,:) ! hail internal energy tend [J/(kg s)]

   real, target, allocatable :: con_gccnt(:,:) ! GCCN number tendency [#_gccn/(m^3 s)]
   real, target, allocatable :: con_ifnt (:,:) ! IFN  number tendency [#_ifn/(m^3 s)]

   real, target, allocatable :: tket    (:,:) ! subgrid-scale turb KE tend [m^2/s^3]
   real, target, allocatable :: epst    (:,:) ! subgrid dissipation rate tend [m^2/s^4]

   integer :: num_omic = 0
   
Contains

!===============================================================================

   subroutine alloc_tend(lza,lva,lwa,naddsc,nccntyp)

   use mem_turb,   only: tkep, epsp
   use mem_basic,  only: vmc, wmc, thil, sh_w, vxe, vye, vze
   use mem_addsc,  only: addsc
   use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,        &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h,&
                         ccntyp, con_ifn, con_gccn, q2, q6, q7
   use misc_coms,  only: io6
   use oname_coms, only: nl
   
   implicit none

   integer, intent(in) :: lza,lva,lwa,naddsc,nccntyp

   integer :: iaddsc, ic

   logical :: qxtrans = .true. ! Deprecated namelist flag for optional Q2/Q6/Q7 transport

   write(io6,*) 'enter alloc_tend'

! Find the maximum number of grid points needed for any grid.

   if (allocated(vxe))     allocate (vmxet(lza,lwa))
   if (allocated(vye))     allocate (vmyet(lza,lwa))
   if (allocated(vze))     allocate (vmzet(lza,lwa))

   if (allocated(thil))    allocate (thilt(lza,lwa)) ; thilt = 0.
   if (allocated(sh_w))    allocate (sh_wt(lza,lwa)) ; sh_wt = 0.

   if (allocated(sh_c))    allocate (sh_ct(lza,lwa))
   if (allocated(sh_r))    allocate (sh_rt(lza,lwa))
   if (allocated(sh_p))    allocate (sh_pt(lza,lwa))
   if (allocated(sh_s))    allocate (sh_st(lza,lwa))
   if (allocated(sh_a))    allocate (sh_at(lza,lwa))
   if (allocated(sh_g))    allocate (sh_gt(lza,lwa))
   if (allocated(sh_h))    allocate (sh_ht(lza,lwa))
   if (allocated(sh_d))    allocate (sh_dt(lza,lwa))

   if (allocated(con_c))   allocate (con_ct(lza,lwa))
   if (allocated(con_r))   allocate (con_rt(lza,lwa))
   if (allocated(con_p))   allocate (con_pt(lza,lwa))
   if (allocated(con_s))   allocate (con_st(lza,lwa))
   if (allocated(con_a))   allocate (con_at(lza,lwa))
   if (allocated(con_g))   allocate (con_gt(lza,lwa))
   if (allocated(con_h))   allocate (con_ht(lza,lwa))
   if (allocated(con_d))   allocate (con_dt(lza,lwa))

   if (allocated(con_gccn))allocate (con_gccnt(lza,lwa))
   if (allocated(con_ifn)) allocate (con_ifnt (lza,lwa))

   do ic = 1,nccntyp
      if       (allocated(ccntyp(ic)%con_ccn) .and.  &
         (.not. allocated(ccntyp(ic)%con_ccnt)))      &
                allocate (ccntyp(ic)%con_ccnt(lza,lwa))
   enddo

   if (qxtrans) then
      if (allocated(q2))   allocate (q2t(lza,lwa))
      if (allocated(q6))   allocate (q6t(lza,lwa))
      if (allocated(q7))   allocate (q7t(lza,lwa))
   endif

   if (allocated(tkep))    allocate (tket(lza,lwa))
   if (allocated(epsp))    allocate (epst(lza,lwa))

   do iaddsc = 1,naddsc
      if       (allocated(addsc(iaddsc)%sclp) .and.  &
         (.not. allocated(addsc(iaddsc)%sclt)))      &
                allocate (addsc(iaddsc)%sclt(lza,lwa))
   enddo

   end subroutine alloc_tend

!===============================================================================
               
   subroutine dealloc_tend(naddsc,nccntyp)
   
   use mem_addsc, only: addsc
   use mem_micro, only: ccntyp

   implicit none

   integer, intent(in) :: naddsc, nccntyp
   
   integer :: iaddsc, ic

! Deallocate all tendency arrays
 
   if (allocated(vmxet))    deallocate (vmxet)
   if (allocated(vmyet))    deallocate (vmyet)
   if (allocated(vmzet))    deallocate (vmzet)

   if (allocated(thilt))    deallocate (thilt)
   if (allocated(sh_wt))    deallocate (sh_wt)

   if (allocated(sh_ct))    deallocate (sh_ct)
   if (allocated(sh_rt))    deallocate (sh_rt)
   if (allocated(sh_pt))    deallocate (sh_pt)
   if (allocated(sh_st))    deallocate (sh_st)
   if (allocated(sh_at))    deallocate (sh_at)
   if (allocated(sh_gt))    deallocate (sh_gt)
   if (allocated(sh_ht))    deallocate (sh_ht)
   if (allocated(sh_dt))    deallocate (sh_dt)

   if (allocated(con_ct))   deallocate (con_ct)
   if (allocated(con_rt))   deallocate (con_rt)
   if (allocated(con_pt))   deallocate (con_pt)
   if (allocated(con_st))   deallocate (con_st)
   if (allocated(con_at))   deallocate (con_at)
   if (allocated(con_gt))   deallocate (con_gt)
   if (allocated(con_ht))   deallocate (con_ht)
   if (allocated(con_dt))   deallocate (con_dt)

   if (allocated(con_gccnt))deallocate (con_gccnt)
   if (allocated(con_ifnt)) deallocate (con_ifnt)

   do ic = 1,nccntyp
      if (allocated(ccntyp(ic)%con_ccn)) deallocate (ccntyp(ic)%con_ccn)
   enddo
        
   if (allocated(q2t))      deallocate (q2t)
   if (allocated(q6t))      deallocate (q6t)
   if (allocated(q7t))      deallocate (q7t)

   if (allocated(tket))     deallocate (tket)
   if (allocated(epst))     deallocate (epst)

   do iaddsc = 1,naddsc
      if (allocated(addsc(iaddsc)%sclt)) deallocate (addsc(iaddsc)%sclt)
   enddo

   end subroutine dealloc_tend

!===============================================================================
               
   subroutine filltab_tend(naddsc,nccntyp)

   use mem_turb,   only: tkep, epsp, sxfer_tk, sxfer_rk
   use mem_basic,  only: thil, sh_w
   use mem_addsc,  only: addsc
   use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,        &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h,&
                         ccntyp, con_ifn, con_gccn, q2, q6, q7
   use var_tables, only: vtables_scalar, num_var
   use misc_coms,  only: do_chem
   use cgrid_defn, only: cgrid_scalar_tabs

   implicit none

   integer, intent(in) :: naddsc, nccntyp
   
   integer :: iaddsc, ic
   character (len=10) :: sname

! Fill pointers to scalar arrays into scalar tables

   if (allocated(thilt))    call vtables_scalar (thil,thilt, 'THIL', sxfer=sxfer_tk)
   if (allocated(sh_wt))    call vtables_scalar (sh_w,sh_wt, 'SH_W', sxfer=sxfer_rk)

   if (allocated(sh_ct))    call vtables_scalar (sh_c, sh_ct, 'SH_C')
   if (allocated(sh_rt))    call vtables_scalar (sh_r, sh_rt, 'SH_R')
   if (allocated(sh_pt))    call vtables_scalar (sh_p, sh_pt, 'SH_P')
   if (allocated(sh_st))    call vtables_scalar (sh_s, sh_st, 'SH_S')
   if (allocated(sh_at))    call vtables_scalar (sh_a, sh_at, 'SH_A')
   if (allocated(sh_gt))    call vtables_scalar (sh_g, sh_gt, 'SH_G')
   if (allocated(sh_ht))    call vtables_scalar (sh_h, sh_ht, 'SH_H')
   if (allocated(sh_dt))    call vtables_scalar (sh_d, sh_dt, 'SH_D')

   if (allocated(con_ct))   call vtables_scalar (con_c, con_ct, 'CON_C')
   if (allocated(con_rt))   call vtables_scalar (con_r, con_rt, 'CON_R')
   if (allocated(con_pt))   call vtables_scalar (con_p, con_pt, 'CON_P')
   if (allocated(con_st))   call vtables_scalar (con_s, con_st, 'CON_S')
   if (allocated(con_at))   call vtables_scalar (con_a, con_at, 'CON_A')
   if (allocated(con_gt))   call vtables_scalar (con_g, con_gt, 'CON_G')
   if (allocated(con_ht))   call vtables_scalar (con_h, con_ht, 'CON_H')
   if (allocated(con_dt))   call vtables_scalar (con_d, con_dt, 'CON_D')

   if (allocated(con_gccnt))call vtables_scalar (con_gccn, con_gccnt, 'CON_GCCN', cu_mix=.true.)
   if (allocated(con_ifnt)) call vtables_scalar (con_ifn,  con_ifnt,  'CON_IFN',  cu_mix=.true.)

   do ic = 1,nccntyp
      write(sname,'(a7,i3.3)') 'CON_CCN',ic
      
      if (allocated(ccntyp(ic)%con_ccn)) then
         call vtables_scalar (ccntyp(ic)%con_ccn, ccntyp(ic)%con_ccnt, sname, cu_mix=.true.)
      endif
   enddo

   if (allocated(q2t))      call vtables_scalar (q2, q2t, 'Q2')
   if (allocated(q6t))      call vtables_scalar (q6, q6t, 'Q6')
   if (allocated(q7t))      call vtables_scalar (q7, q7t, 'Q7')

   num_omic = num_var

   if (allocated(tket))     call vtables_scalar (tkep, tket, 'TKEP')
   if (allocated(epst))     call vtables_scalar (epsp, epst, 'EPSP')

   if (do_chem == 1) call cgrid_scalar_tabs()

   do iaddsc = 1,naddsc
      write(sname,'(a4,i3.3)') 'SCLP',iaddsc
      
      if (allocated(addsc(iaddsc)%sclt)) then
         call vtables_scalar (addsc(iaddsc)%sclp, addsc(iaddsc)%sclt, trim(sname), cu_mix=.true.)
      endif
   enddo

   end subroutine filltab_tend

End Module mem_tend
