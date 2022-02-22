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
subroutine olam_mem_alloc()

  use mem_basic,   only: alloc_basic, filltab_basic
  use mem_cuparm,  only: alloc_cuparm, filltab_cuparm
  use mem_micro,   only: alloc_micro, filltab_micro
  use mem_co2,     only: alloc_co2, filltab_co2
  use mem_radiate, only: alloc_radiate, filltab_radiate
  use mem_addsc,   only: alloc_addsc, filltab_addsc
  use mem_tend,    only: alloc_tend, filltab_tend
  use mem_turb,    only: alloc_turb, filltab_turb
  use mem_grid,    only: alloc_grid_other, mza, nsw_max, mva, mwa
  use mem_nudge,   only: nudflag, nudnxp, mwnud, alloc_nudge2, filltab_nudge, &
                         o3nudflag, alloc_nudge_o3, filltab_nudge_o3
  use mem_ijtabs,  only: mrls
  use oname_coms,  only: nl
  use wrtv_rk,     only: alloc_wrtv_rk
  use wrtv_orig,   only: alloc_wrtv_orig

  use misc_coms,   only: io6, naddsc, initial, idiffk, ilwrtyp, iswrtyp,  &
                         nqparm, do_chem, nrk_wrtv

  use micro_coms,  only: miclevel, ncat, jnmb, iccn, igccn, iifn
  use ccnbin_coms, only: nccntyp

  use leaf_coms,   only: isfcl
  use mem_sfcg,    only: mwsfc
  use mem_land,    only: mland
  use mem_sea,     only: msea

  use mem_flux_accum, only: alloc_flux_accum, filltab_flux_accum

  use mem_average_vars, only: alloc_average_vars

  use cgrid_defn,  only: alloc_cgrid, filltab_cgrid

  use mem_megan,   only: alloc_megan, filltab_megan
  use soil_nox,    only: alloc_soil_nox, filltab_soil_nox

  implicit none

! Allocate basic memory and fill variable tables

  call alloc_grid_other()

  call alloc_basic(mza,mva,mwa)
  call filltab_basic()

  call alloc_cuparm(mza, mwa, mrls, nqparm)
  call filltab_cuparm()

  call alloc_micro(mza,mwa,miclevel,ncat,nccntyp,iccn,igccn,iifn,jnmb)
  call filltab_micro(nccntyp)

  call alloc_co2(mza,mwa)
  call filltab_co2()

  call alloc_radiate(mza,mwa,nsw_max,ilwrtyp,iswrtyp)
  call filltab_radiate()

  call alloc_turb(mza,mwa,mva,nsw_max,idiffk,mrls)
  call filltab_turb()

  call alloc_flux_accum(mza,mva,mwa,mwsfc,mland)
  call filltab_flux_accum()

  if (do_chem == 1) then
     call alloc_cgrid(mza,mwa)
     call filltab_cgrid()

     if (isfcl > 0) then
        call alloc_megan(mland,mwa)
        call filltab_megan()

        call alloc_soil_nox()
        call filltab_soil_nox()
     endif
  endif

! Allocate field average arrays

  call alloc_average_vars(mwa,mwsfc,mland)

  if (initial == 2 .and. nudflag > 0) then

! If doing point-by-point (non-spectral) nudging, define mwnud here

      if (nudnxp == 0) mwnud = mwa

      call alloc_nudge2(mza,mwa)
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

  call alloc_tend(mza,mva,mwa,naddsc,nccntyp)
  call filltab_tend(naddsc,nccntyp)

! Extra memory for prog_wrtv

  if (nrk_wrtv == 1) then
     call alloc_wrtv_orig()
  else
     call alloc_wrtv_rk()
  endif

! Extra memory if nudging ozone. Must be called after scalar tables are set up

  if (initial == 2 .and. o3nudflag > 0) then

     ! point-by-point (non-spectral) nudging of ozone

      call alloc_nudge_o3(mza,mwa)
      call filltab_nudge_o3()
  endif

  write(io6,*) 'end alloc'

end subroutine olam_mem_alloc

!===============================================================================

!subroutine dealloc_all()

!use mem_ijtabs
!use mem_basic
!use mem_cuparm
!use mem_grid
!use mem_land
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
