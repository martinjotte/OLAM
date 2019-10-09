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
   real, allocatable :: vmt  (:,:) ! V momentum tend [kg/(m^2 s^2)]

   real, allocatable :: thilt   (:,:) ! (rho * thil) tend [kg_air K / (m^3 s)]
   real, allocatable :: rr_wt   (:,:) ! water mass tend [kg_wat/(m^3 s)]

   real, allocatable :: rr_ct   (:,:) ! cloud mass tend [kg_cld/(m^3 s)]
   real, allocatable :: rr_rt   (:,:) ! rain mass tend [kg_rain/(m^3 s)]
   real, allocatable :: rr_pt   (:,:) ! pristine ice mass tend [kg_pris/(m^3 s)]
   real, allocatable :: rr_st   (:,:) ! snow mass tend [kg_snow/(m^3 s)]
   real, allocatable :: rr_at   (:,:) ! aggregates mass tend [kg_agg/(m^3 s)]
   real, allocatable :: rr_gt   (:,:) ! graupel mass tend [kg_grp/(m^3 s)]
   real, allocatable :: rr_ht   (:,:) ! hail mass tend [kg_hail/(m^3 s)]
   real, allocatable :: rr_dt   (:,:) ! drizzle mass tend [kg_drz/(m^3 s)]

   real, allocatable :: con_ct  (:,:) ! cloud number tend [#/(m^3 s)]]
   real, allocatable :: con_rt  (:,:) ! rain number tend [#/(m^3 s)]]
   real, allocatable :: con_pt  (:,:) ! pristine ice number tend [#/(m^3 s)]]
   real, allocatable :: con_st  (:,:) ! snow number tend [#/(m^3 s)]]
   real, allocatable :: con_at  (:,:) ! aggregates number tend [#/(m^3 s)]]
   real, allocatable :: con_gt  (:,:) ! graupel number tend [#/(m^3 s)]]
   real, allocatable :: con_ht  (:,:) ! hail number tend [#/(m^3 s)]]
   real, allocatable :: con_dt  (:,:) ! drizzle number tend [#/(m^3 s)]]

   real, allocatable :: q2t     (:,:) ! rain internal energy tend [J/(kg s)]
   real, allocatable :: q6t     (:,:) ! graupel internal energy tend [J/(kg s)]
   real, allocatable :: q7t     (:,:) ! hail internal energy tend [J/(kg s)]

   real, allocatable :: con_gccnt(:,:) ! GCCN number tendency [#_gccn/(m^3 s)]
   real, allocatable :: con_ifnt (:,:) ! IFN  number tendency [#_ifn/(m^3 s)]

   real, allocatable :: tket    (:,:) ! subgrid-scale turb KE tend [m^2/s^3]
   real, allocatable :: epst    (:,:) ! subgrid dissipation rate tend [m^2/s^4]

   real, allocatable :: rr_co2t (:,:) ! CO2 mass tend [kg_CO2/(m^3 s)]

   integer :: num_omic = 0

Contains

!===============================================================================

   subroutine alloc_tend(lza,lva,lwa,naddsc,nccntyp)

   use mem_turb,   only: tkep, epsp
   use mem_basic,  only: thil, rr_w, vxe, vye, vze, vc
   use mem_addsc,  only: addsc
   use mem_micro,  only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h,        &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h,&
                         ccntyp, con_ifn, con_gccn, q2, q6, q7
   use micro_coms, only: miclevel
   use misc_coms,  only: io6
   use mem_co2,    only: rr_co2

   implicit none

   integer, intent(in) :: lza,lva,lwa,naddsc,nccntyp

   integer :: iaddsc, ic

   write(io6,*) 'enter alloc_tend'

! Find the maximum number of grid points needed for any grid.

   if (allocated(vxe))     allocate (vmxet(lza,lwa))
   if (allocated(vye))     allocate (vmyet(lza,lwa))
   if (allocated(vze))     allocate (vmzet(lza,lwa))

   if (allocated(vc))      allocate (vmt(lza,lva))

   if (allocated(thil))    allocate (thilt(lza,lwa)) ; thilt = 0.
   if (allocated(rr_w))    allocate (rr_wt(lza,lwa)) ; rr_wt = 0.

   if (miclevel > 2) then
      if (allocated(rr_c)) allocate (rr_ct(lza,lwa))
   endif

   if (allocated(rr_r))    allocate (rr_rt(lza,lwa))
   if (allocated(rr_p))    allocate (rr_pt(lza,lwa))
   if (allocated(rr_s))    allocate (rr_st(lza,lwa))
   if (allocated(rr_a))    allocate (rr_at(lza,lwa))
   if (allocated(rr_g))    allocate (rr_gt(lza,lwa))
   if (allocated(rr_h))    allocate (rr_ht(lza,lwa))
   if (allocated(rr_d))    allocate (rr_dt(lza,lwa))

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

   if (allocated(ccntyp)) then
      do ic = 1,nccntyp
         if       (allocated(ccntyp(ic)%con_ccn) .and.   &
            (.not. allocated(ccntyp(ic)%con_ccnt)))      &
                   allocate (ccntyp(ic)%con_ccnt(lza,lwa))
      enddo
   endif

   if (allocated(q2))   allocate (q2t(lza,lwa))
   if (allocated(q6))   allocate (q6t(lza,lwa))
   if (allocated(q7))   allocate (q7t(lza,lwa))

   if (allocated(tkep)) allocate (tket(lza,lwa))
   if (allocated(epsp)) allocate (epst(lza,lwa))

   do iaddsc = 1,naddsc
      if       (allocated(addsc(iaddsc)%sclp) .and.  &
         (.not. allocated(addsc(iaddsc)%sclt)))      &
                allocate (addsc(iaddsc)%sclt(lza,lwa))
   enddo

   if (allocated(rr_co2)) allocate (rr_co2t(lza,lwa))

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

   if (allocated(vmt))      deallocate (vmt)

   if (allocated(thilt))    deallocate (thilt)
   if (allocated(rr_wt))    deallocate (rr_wt)

   if (allocated(rr_ct))    deallocate (rr_ct)
   if (allocated(rr_rt))    deallocate (rr_rt)
   if (allocated(rr_pt))    deallocate (rr_pt)
   if (allocated(rr_st))    deallocate (rr_st)
   if (allocated(rr_at))    deallocate (rr_at)
   if (allocated(rr_gt))    deallocate (rr_gt)
   if (allocated(rr_ht))    deallocate (rr_ht)
   if (allocated(rr_dt))    deallocate (rr_dt)

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

   if (allocated(ccntyp)) then
      do ic = 1,nccntyp
         if (allocated(ccntyp(ic)%con_ccn)) deallocate (ccntyp(ic)%con_ccn)
      enddo
   endif

   if (allocated(q2t))      deallocate (q2t)
   if (allocated(q6t))      deallocate (q6t)
   if (allocated(q7t))      deallocate (q7t)

   if (allocated(tket))     deallocate (tket)
   if (allocated(epst))     deallocate (epst)

   do iaddsc = 1,naddsc
      if (allocated(addsc(iaddsc)%sclt)) deallocate (addsc(iaddsc)%sclt)
   enddo

   if (allocated(rr_co2t)) deallocate(rr_co2t)

   end subroutine dealloc_tend

!===============================================================================

   subroutine filltab_tend(naddsc,nccntyp)

   use mem_turb,   only: tkep, epsp, sxfer_rk
   use mem_basic,  only: rr_w
   use mem_addsc,  only: addsc
   use mem_micro,  only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h,        &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h,&
                         ccntyp, con_ifn, con_gccn, q2, q6, q7
   use var_tables, only: vtables_scalar, num_scalar, scalar_tab
   use misc_coms,  only: do_chem, i_o3, i_co, i_ch4
   use cgrid_defn, only: cgrid_scalar_tabs
   use mem_co2,    only: rr_co2, i_co2

   implicit none

   integer, intent(in) :: naddsc, nccntyp

   integer :: iaddsc, ic, n
   character (len=10) :: sname

! Fill pointers to scalar arrays into scalar tables

   if (allocated(rr_wt))    call vtables_scalar (rr_w,rr_wt, 'RR_W', sxfer=sxfer_rk)

   if (allocated(rr_ct))    call vtables_scalar (rr_c, rr_ct, 'RR_C')
   if (allocated(rr_rt))    call vtables_scalar (rr_r, rr_rt, 'RR_R')
   if (allocated(rr_pt))    call vtables_scalar (rr_p, rr_pt, 'RR_P')
   if (allocated(rr_st))    call vtables_scalar (rr_s, rr_st, 'RR_S')
   if (allocated(rr_at))    call vtables_scalar (rr_a, rr_at, 'RR_A')
   if (allocated(rr_gt))    call vtables_scalar (rr_g, rr_gt, 'RR_G')
   if (allocated(rr_ht))    call vtables_scalar (rr_h, rr_ht, 'RR_H')
   if (allocated(rr_dt))    call vtables_scalar (rr_d, rr_dt, 'RR_D')

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

   if (allocated(ccntyp)) then
      do ic = 1,nccntyp
         write(sname,'(a7,i3.3)') 'CON_CCN',ic
         if (allocated(ccntyp(ic)%con_ccn)) then
            call vtables_scalar (ccntyp(ic)%con_ccn, ccntyp(ic)%con_ccnt, sname, cu_mix=.true.)
         endif
      enddo
   endif

   if (allocated(q2t)) call vtables_scalar (q2, q2t, 'Q2', pos_def=.false.)
   if (allocated(q6t)) call vtables_scalar (q6, q6t, 'Q6', pos_def=.false.)
   if (allocated(q7t)) call vtables_scalar (q7, q7t, 'Q7', pos_def=.false.)

   num_omic = num_scalar

   if (allocated(tket)) call vtables_scalar (tkep, tket, 'TKEP')
   if (allocated(epst)) call vtables_scalar (epsp, epst, 'EPSP')

   if (do_chem == 1) call cgrid_scalar_tabs()

   do iaddsc = 1,naddsc
      write(sname,'(a4,i3.3)') 'SCLP',iaddsc

      if (allocated(addsc(iaddsc)%sclt)) then
         call vtables_scalar (addsc(iaddsc)%sclp, addsc(iaddsc)%sclt, trim(sname), cu_mix=.true.)
      endif
   enddo

   if (allocated(rr_co2t)) then
      call vtables_scalar (rr_co2, rr_co2t, 'RR_CO2', cu_mix=.true.)
      i_co2 = num_scalar  ! save index of co2 in scalar table
   endif

    ! If any prognostic scalars is ozone, save its scalar index for nudging and radiation

    i_o3 = 0
    do n = 1, num_scalar
       if ( (scalar_tab(n)%name == 'O3'   ) .or. &
            (scalar_tab(n)%name == 'o3'   ) .or. &
            (scalar_tab(n)%name == 'OZONE') .or. &
            (scalar_tab(n)%name == 'ozone') ) then
          i_o3 = n
          exit
       endif
    enddo

    ! If any prognostic scalars is methane, save its index for radiation

    i_ch4 = 0
    do n = 1, num_scalar
       if ( (scalar_tab(n)%name == 'CH4'    ) .or. &
            (scalar_tab(n)%name == 'ch4'    ) .or. &
            (scalar_tab(n)%name == 'METHANE') .or. &
            (scalar_tab(n)%name == 'methane') ) then
          i_ch4 = n
          exit
       endif
    enddo

    ! If any prognostic scalars is carbon monoxide, save its index for radiation

    i_co = 0
    do n = 1, num_scalar
       if ( (scalar_tab(n)%name == 'CO'    ) .or. &
            (scalar_tab(n)%name == 'co'    ) ) then
          i_co = n
          exit
       endif
    enddo

 end subroutine filltab_tend

End Module mem_tend
