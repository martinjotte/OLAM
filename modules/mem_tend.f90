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
Module mem_tend

   real, allocatable :: umt     (:,:) ! U-mom density tend [kg/(m^2 s^2)]
   real, allocatable :: vmt     (:,:) ! V-mom density tend [kg/(m^2 s^2)]
   real, allocatable :: wmt     (:,:) ! W-mom density tend [kg/(m^2 s^2)]
   real, allocatable :: sh_wt   (:,:) ! water mass tend [kg_wat/(m^3 s)]
   real, target, allocatable :: thilt(:,:) ! (rho * thil) tend [kg_air K / (m^3 s)]
   real, allocatable :: sh_ct   (:,:) ! cloud mass tend [kg_cld/(m^3 s)]
   real, allocatable :: sh_rt   (:,:) ! rain mass tend [kg_rain/(m^3 s)]
   real, allocatable :: sh_pt   (:,:) ! pristine ice mass tend [kg_pris/(m^3 s)]
   real, allocatable :: sh_st   (:,:) ! snow mass tend [kg_snow/(m^3 s)]
   real, allocatable :: sh_at   (:,:) ! aggregates mass tend [kg_agg/(m^3 s)]
   real, allocatable :: sh_gt   (:,:) ! graupel mass tend [kg_grp/(m^3 s)]
   real, allocatable :: sh_ht   (:,:) ! hail mass tend [kg_hail/(m^3 s)]
   real, allocatable :: con_ct  (:,:) ! cloud number tend [#/(m^3 s)]]
   real, allocatable :: con_rt  (:,:) ! rain number tend [#/(m^3 s)]]
   real, allocatable :: con_pt  (:,:) ! pristine ice number tend [#/(m^3 s)]]
   real, allocatable :: con_st  (:,:) ! snow number tend [#/(m^3 s)]]
   real, allocatable :: con_at  (:,:) ! aggregates number tend [#/(m^3 s)]]
   real, allocatable :: con_gt  (:,:) ! graupel number tend [#/(m^3 s)]]
   real, allocatable :: con_ht  (:,:) ! hail number tend [#/(m^3 s)]]
   real, allocatable :: con_ccnt(:,:) ! ccn number tend [#/(m^3 s)]
   real, allocatable :: con_ifnt(:,:) ! ifn number tend [#/(m^3 s)]
   real, allocatable :: tket    (:,:) ! subgrid-scale turb KE [m^2/s^2]
   real, allocatable :: epst    (:,:) ! 
   
Contains

!===============================================================================

   subroutine alloc_tend(lza,lua,lva,lwa,naddsc)

   use mem_turb,   only: tkep, epsp
   use mem_basic,  only: umc, vmc, wmc, thil, sh_w
   use mem_addsc,  only: addsc
   use mem_micro,  only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,  &
                         con_c, con_r, con_p, con_s, con_a, con_g, con_h,  &
                         con_ccn, con_ifn
   use misc_coms,  only: io6
   
   implicit none

   integer, intent(in) :: lza,lua,lva,lwa,naddsc

   integer :: iaddsc

write(io6,*) 'enter alloc_tend'

! Find the maximum number of grid points needed for any grid.

   if (allocated(umc))     allocate (umt(lza,lua))
   if (allocated(vmc))     allocate (vmt(lza,lva))
   if (allocated(wmc))     allocate (wmt(lza,lwa))

   if (allocated(thil))    allocate (thilt(lza,lwa))
   if (allocated(sh_w))    allocate (sh_wt(lza,lwa))
   if (allocated(sh_c))    allocate (sh_ct(lza,lwa))
   if (allocated(sh_r))    allocate (sh_rt(lza,lwa))
   if (allocated(sh_p))    allocate (sh_pt(lza,lwa))
   if (allocated(sh_s))    allocate (sh_st(lza,lwa))
   if (allocated(sh_a))    allocate (sh_at(lza,lwa))
   if (allocated(sh_g))    allocate (sh_gt(lza,lwa))
   if (allocated(sh_h))    allocate (sh_ht(lza,lwa))
   if (allocated(con_c))   allocate (con_ct(lza,lwa))
   if (allocated(con_r))   allocate (con_rt(lza,lwa))
   if (allocated(con_p))   allocate (con_pt(lza,lwa))
   if (allocated(con_s))   allocate (con_st(lza,lwa))
   if (allocated(con_a))   allocate (con_at(lza,lwa))
   if (allocated(con_g))   allocate (con_gt(lza,lwa))
   if (allocated(con_h))   allocate (con_ht(lza,lwa))
   if (allocated(con_ccn)) allocate (con_ccnt(lza,lwa))
   if (allocated(con_ifn)) allocate (con_ifnt(lza,lwa))
   if (allocated(tkep))    allocate (tket(lza,lwa))
   if (allocated(epsp))    allocate (epst(lza,lwa))

   do iaddsc = 1,naddsc

      if       (allocated(addsc(iaddsc)%sclp) .and.  &
         (.not. allocated(addsc(iaddsc)%sclt)))      &
                allocate (addsc(iaddsc)%sclt(lza,lwa))
   enddo

   return
   end subroutine alloc_tend

!===============================================================================
               
   subroutine dealloc_tend(naddsc)
   
   use mem_addsc, only: addsc

   implicit none

   integer, intent(in) :: naddsc   
   
   integer :: iaddsc

! Deallocate all tendency arrays
 
   if (allocated(umt))      deallocate (umt)
   if (allocated(vmt))      deallocate (vmt)
   if (allocated(wmt))      deallocate (wmt)
   if (allocated(thilt))    deallocate (thilt)
   if (allocated(sh_wt))    deallocate (sh_wt)
   if (allocated(sh_ct))    deallocate (sh_ct)
   if (allocated(sh_rt))    deallocate (sh_rt)
   if (allocated(sh_pt))    deallocate (sh_pt)
   if (allocated(sh_st))    deallocate (sh_st)
   if (allocated(sh_at))    deallocate (sh_at)
   if (allocated(sh_gt))    deallocate (sh_gt)
   if (allocated(sh_ht))    deallocate (sh_ht)
   if (allocated(con_ct))   deallocate (con_ct)
   if (allocated(con_rt))   deallocate (con_rt)
   if (allocated(con_pt))   deallocate (con_pt)
   if (allocated(con_st))   deallocate (con_st)
   if (allocated(con_at))   deallocate (con_at)
   if (allocated(con_gt))   deallocate (con_gt)
   if (allocated(con_ht))   deallocate (con_ht)
   if (allocated(con_ccnt)) deallocate (con_ccnt)
   if (allocated(con_ifnt)) deallocate (con_ifnt)

   if (allocated(tket))     deallocate (tket)
   if (allocated(epst))     deallocate (epst)

   do iaddsc = 1,naddsc
      if (allocated(addsc(iaddsc)%sclt)) deallocate (addsc(iaddsc)%sclt)
   enddo
        
   return
   end subroutine dealloc_tend

!===============================================================================
               
   subroutine filltab_tend(naddsc)

   use mem_turb,   only: tkep, epsp
   use mem_basic,  only: thil, sh_w
   use mem_addsc,  only: addsc
   use mem_micro,  only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,  &
                         con_c, con_r, con_p, con_s, con_a, con_g, con_h,  &
                         con_ccn, con_ifn
   use var_tables, only: vtables_scalar
   use misc_coms,  only: io6

   implicit none

   integer, intent(in) :: naddsc
   
   integer :: iaddsc
   character (len=7) :: sname

! Fill pointers to scalar arrays into scalar tables

   if (allocated(sh_wt))    call vtables_scalar (sh_w,sh_wt,'SH_W')
   if (allocated(sh_ct))    call vtables_scalar (sh_c,sh_ct,'SH_C')
   if (allocated(sh_rt))    call vtables_scalar (sh_r,sh_rt,'SH_R')
   if (allocated(sh_pt))    call vtables_scalar (sh_p,sh_pt,'SH_P')
   if (allocated(sh_st))    call vtables_scalar (sh_s,sh_st,'SH_S')
   if (allocated(sh_at))    call vtables_scalar (sh_a,sh_at,'SH_A')
   if (allocated(sh_gt))    call vtables_scalar (sh_g,sh_gt,'SH_G')
   if (allocated(sh_ht))    call vtables_scalar (sh_h,sh_ht,'SH_H')
   if (allocated(con_ct))   call vtables_scalar (con_c,con_ct,'CON_C')
   if (allocated(con_rt))   call vtables_scalar (con_r,con_rt,'CON_R')
   if (allocated(con_pt))   call vtables_scalar (con_p,con_pt,'CON_P')
   if (allocated(con_st))   call vtables_scalar (con_s,con_st,'CON_S')
   if (allocated(con_at))   call vtables_scalar (con_a,con_at,'CON_A')
   if (allocated(con_gt))   call vtables_scalar (con_g,con_gt,'CON_G')
   if (allocated(con_ht))   call vtables_scalar (con_h,con_ht,'CON_H')
   if (allocated(con_ccnt)) call vtables_scalar (con_ccn,con_ccnt,'CON_CCN')
   if (allocated(con_ifnt)) call vtables_scalar (con_ifn,con_ifnt,'CON_IFN')

   if (allocated(tket))     call vtables_scalar (tkep,tket,'TKEP')
   if (allocated(epst))     call vtables_scalar (epsp,epst,'EPSP')
        
   do iaddsc = 1,naddsc
      write(sname,'(a4,i3.3)') 'SCLP',iaddsc
      
      if (allocated(addsc(iaddsc)%sclt))  &
         call vtables_scalar (addsc(iaddsc)%sclp,addsc(iaddsc)%sclt,sname)
   enddo

   return
   end subroutine filltab_tend

End Module mem_tend
