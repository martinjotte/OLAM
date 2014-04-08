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
subroutine olam_mem_alloc()

  use mem_basic,   only: alloc_basic, filltab_basic
  use mem_cuparm,  only: alloc_cuparm, filltab_cuparm
  use mem_micro,   only: alloc_micro, filltab_micro
  use mem_radiate, only: alloc_radiate, filltab_radiate
  use mem_addsc,   only: alloc_addsc, filltab_addsc
  use mem_tend,    only: alloc_tend, filltab_tend
  use mem_turb,    only: alloc_turb, filltab_turb
  use mem_sflux,   only: mseaflux, mlandflux, filltab_sflux
  use mem_grid,    only: mza, nsw_max, mma, mua, mva, mwa
  use mem_nudge,   only: nudflag, nudnxp, mwnud, alloc_nudge2, filltab_nudge
  use mem_ijtabs,  only: mrls, filltab_itabs
  use oname_coms,  only: nl
  use mem_thuburn, only: alloc_thuburn

  use misc_coms,   only: io6, naddsc, initial, idiffk, ilwrtyp, iswrtyp,  &
                         nqparm, dtsm, meshtype

  use micro_coms,  only: level, ncat, jnmb, &
                         icloud, idriz, irain, ipris, isnow, iaggr, igraup, ihail
                       
  use leaf_coms,   only: mwl
  use sea_coms,    only: mws
  use mem_sflux,   only: mseaflux, mlandflux
  use mem_plot,    only: alloc_plot

  use mem_timeavg, only: alloc_timeavg, filltab_timeavg

  use mem_average_vars, only: alloc_average_vars

  implicit none 

  integer :: ng,nv,ntpts

! Allocate basic memory and fill variable tables

  call filltab_itabs()  ! Already allocated

  call alloc_basic(meshtype,mza,mua,mva,mwa)
  call filltab_basic()

  call alloc_cuparm(mza, mwa, mrls, nqparm)
  call filltab_cuparm() 

  call alloc_micro(mza,mwa,level,ncat, &
     icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail)
  call filltab_micro()

  call alloc_radiate(mza,mwa,ilwrtyp,iswrtyp)
  call filltab_radiate()

  call alloc_turb(mza,mwa,nsw_max,idiffk(1),mrls)
  call filltab_turb()

  call alloc_timeavg(mza,mwa)
  call filltab_timeavg()

! Allocate field average arrays

  call alloc_average_vars(mwa,mwl,mws,mlandflux,mseaflux)

  call filltab_sflux()  ! Already allocated

  if (initial == 2 .and. nudflag > 0) then

! If doing point-by-point (non-spectral) nudging, define mwnud here

      if (nudnxp == 0) mwnud = mwa

      call alloc_nudge2(mza)
      call filltab_nudge()
  endif

! Allocate any added Scalar types 

  write(io6,*) 'start addsc alloc'

  if (naddsc > 0) then
     call alloc_addsc(mza,mwa,naddsc)
     call filltab_addsc(naddsc)
  endif

! Allocate tendency data type and fill scalar table.  These subroutine
! calls must occur after allocating all prognostic variables.

  write(io6,*) 'start tendency alloc'

  call alloc_tend(mza,mua,mva,mwa,naddsc)
  call filltab_tend(naddsc) 

! Extra memory for Thuburn's monotonic advection

  call alloc_thuburn(meshtype, nl%iscal_monot, mza, mva, mwa)

! Memory for storing past values for plotting

  call alloc_plot()

  write(io6,*) 'end alloc'

  return
end subroutine olam_mem_alloc

!===============================================================================

!subroutine dealloc_all()

!use mem_ijtabs
!use mem_basic
!use mem_cuparm
!use mem_grid
!use mem_leaf
!use mem_micro
!use mem_radiate
!use mem_addsc
!use mem_tend
!use mem_turb
!use mem_nudge

!use var_tables

!use misc_coms

! deallocate all model memory.  Used on dynamic balance

!call dealloc_basic
!call dealloc_cuparm
!call dealloc_micro
!call dealloc_radiate
!call dealloc_turb

!call dealloc_tend(naddsc)
!call dealloc_addsc(naddsc)

!return
!end
