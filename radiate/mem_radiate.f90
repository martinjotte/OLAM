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

Module mem_radiate

  integer :: jday ! Julian day

  real :: solfac  ! solar-constant coefficient for variable Earth-Sun dist
  real :: sunx    ! x-component of unit vector pointing to sun [m]
  real :: suny    ! y-component of unit vector pointing to sun [m]
  real :: sunz    ! z-component of unit vector pointing to sun [m]

  real, allocatable, target :: fthrd_sw    (:,:)
  real, allocatable, target :: fthrd_lw    (:,:)
  real, allocatable, target :: cloud_frac  (:,:)

  real, allocatable, target :: rshort        (:)
  real, allocatable, target :: rlong         (:)
  real, allocatable, target :: rlongup       (:)
  real, allocatable, target :: rshort_top    (:)
  real, allocatable, target :: rshortup_top  (:)
  real, allocatable, target :: rlongup_top   (:)
  real, allocatable, target :: albedt        (:)
  real, allocatable         :: albedt_beam   (:)
  real, allocatable         :: albedt_diffuse(:)
  real, allocatable, target :: cosz          (:)
  real, allocatable         :: rlong_albedo  (:)
  real, allocatable         :: rshort_diffuse(:)
  real, allocatable         :: par           (:)
  real, allocatable         :: par_diffuse   (:)
  real, allocatable         :: ppfd          (:)
  real, allocatable         :: ppfd_diffuse  (:)
  real, allocatable         :: uva           (:)
  real, allocatable         :: uvb           (:)
  real, allocatable         :: uvc           (:)
  real, allocatable         :: pbl_cld_forc  (:)

  ! arrays for transferring radiation quantities
  ! to surface with shaved cells

  real, allocatable         :: rshort_ks        (:,:)
  real, allocatable         :: rlong_ks         (:,:)
  real, allocatable         :: rshort_diffuse_ks(:,:)
! real, allocatable         :: par_ks           (:,:)
! real, allocatable         :: par_diffuse_ks   (:,:)
  real, allocatable, target :: ppfd_ks          (:,:)
  real, allocatable, target :: ppfd_diffuse_ks  (:,:)

  ! clear-sky values

  real, allocatable, target :: rshort_clr      (:)
  real, allocatable, target :: rshortup_clr    (:)
  real, allocatable, target :: rshort_top_clr  (:)
  real, allocatable, target :: rshortup_top_clr(:)

  real, allocatable, target :: rlong_clr       (:)
  real, allocatable, target :: rlongup_clr     (:)
  real, allocatable, target :: rlongup_top_clr (:)

! RRTMg random cloud overlap seed

  integer, allocatable, target :: mcica_seed(:,:)

! These are used for adding extra levels at the top with the Mclatchy soundings

  integer, parameter :: maxadd_rad = 10 ! max allowed # of added rad levels
  integer            :: nadd_rad        ! actual # of added radiation levels
  real               :: zmrad = 30.e3   ! top of radiation grid

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

!   Allocate arrays based on options (if necessary)

    if (ilwrtyp + iswrtyp > 0)  then

       write(io6,*) 'allocating rad ', mwa, mza

       allocate (fthrd_sw  (mza,mwa)) ; fthrd_sw       = 0.0
       allocate (fthrd_lw  (mza,mwa)) ; fthrd_lw       = 0.0
       allocate (rshort        (mwa)) ; rshort         = 0.0
       allocate (rlong         (mwa)) ; rlong          = 0.0
       allocate (rlongup       (mwa)) ; rlongup        = 0.0
       allocate (rshort_top    (mwa)) ; rshort_top     = 0.0
       allocate (rshortup_top  (mwa)) ; rshortup_top   = 0.0
       allocate (rlongup_top   (mwa)) ; rlongup_top    = 0.0
       allocate (rshort_diffuse(mwa)) ; rshort_diffuse = 0.0

       allocate (par           (mwa)) ; par            = 0.0
       allocate (par_diffuse   (mwa)) ; par_diffuse    = 0.0
       allocate (ppfd          (mwa)) ; ppfd           = 0.0
       allocate (ppfd_diffuse  (mwa)) ; ppfd_diffuse   = 0.0
       allocate (uva           (mwa)) ; uva            = 0.0
       allocate (uvb           (mwa)) ; uvb            = 0.0
       allocate (uvc           (mwa)) ; uvc            = 0.0
       allocate (pbl_cld_forc  (mwa)) ; pbl_cld_forc   = 0.0

       allocate (rlong_albedo  (mwa)) ; rlong_albedo   = rinit
       allocate (albedt        (mwa)) ; albedt         = rinit
       allocate (albedt_beam   (mwa)) ; albedt_beam    = rinit
       allocate (albedt_diffuse(mwa)) ; albedt_diffuse = rinit
       allocate (cosz          (mwa)) ; cosz           = rinit

       allocate (cloud_frac(mza,mwa)) ; cloud_frac  = 0.0

       allocate (rshort_clr      (mwa)) ; rshort_clr       = 0.0
       allocate (rshortup_clr    (mwa)) ; rshortup_clr     = 0.0
       allocate (rshort_top_clr  (mwa)) ; rshort_top_clr   = 0.0
       allocate (rshortup_top_clr(mwa)) ; rshortup_top_clr = 0.0

       allocate (rlong_clr      (mwa)) ; rlong_clr       = 0.0
       allocate (rlongup_clr    (mwa)) ; rlongup_clr     = 0.0
       allocate (rlongup_top_clr(mwa)) ; rlongup_top_clr = 0.0

       allocate (mcica_seed   (4,mwa)) ; mcica_seed = 0

       allocate (rshort_ks        (nsw_max,mwa)) ; rshort_ks         = 0.0
       allocate (rshort_diffuse_ks(nsw_max,mwa)) ; rshort_diffuse_ks = 0.0
       allocate (rlong_ks         (nsw_max,mwa)) ; rlong_ks          = 0.0

!      allocate (par_ks           (nsw_max,mwa)) ; par_ks            = 0.0
!      allocate (par_diffuse_ks   (nsw_max,mwa)) ; par_diffuse_ks    = 0.0
       allocate (ppfd_ks          (nsw_max,mwa)) ; ppfd_ks           = 0.0
       allocate (ppfd_diffuse_ks  (nsw_max,mwa)) ; ppfd_diffuse_ks   = 0.0

    endif

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
    if (allocated(cloud_frac))     deallocate (cloud_frac)

    if (allocated(rshort_clr))       deallocate (rshort_clr)
    if (allocated(rshortup_clr))     deallocate (rshortup_clr)
    if (allocated(rshort_top_clr))   deallocate (rshort_top_clr)
    if (allocated(rshortup_top_clr)) deallocate (rshortup_top_clr)

    if (allocated(rlong_clr))        deallocate (rlong_clr)
    if (allocated(rlongup_clr))      deallocate (rlongup_clr)
    if (allocated(rlongup_top_clr))  deallocate (rlongup_top_clr)

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

    if (allocated(cloud_frac))       call increment_vtable('CLOUD_FRAC',      'AW', rvar2=cloud_frac)

    if (allocated(rshort))           call increment_vtable('RSHORT',          'AW', rvar1=rshort)

    if (allocated(rlong))            call increment_vtable('RLONG',           'AW', rvar1=rlong)

    if (allocated(rlongup))          call increment_vtable('RLONGUP',         'AW', rvar1=rlongup)

    if (allocated(rshort_top))       call increment_vtable('RSHORT_TOP',      'AW', rvar1=rshort_top)

    if (allocated(rshortup_top))     call increment_vtable('RSHORTUP_TOP',    'AW', rvar1=rshortup_top)

    if (allocated(rlongup_top))      call increment_vtable('RLONGUP_TOP',     'AW', rvar1=rlongup_top)

    if (allocated(albedt))           call increment_vtable('ALBEDT',          'AW', rvar1=albedt)

    if (allocated(cosz))             call increment_vtable('COSZ',            'AW', rvar1=cosz)

    if (allocated(rshort_clr))       call increment_vtable('RSHORT_CLR',      'AW', rvar1=rshort_clr)

    if (allocated(rshortup_clr))     call increment_vtable('RSHORTUP_CLR',    'AW', rvar1=rshortup_clr)

    if (allocated(rshort_top_clr))   call increment_vtable('RSHORT_TOP_CLR',  'AW', rvar1=rshort_top_clr)

    if (allocated(rshortup_top_clr)) call increment_vtable('RSHORTUP_TOP_CLR','AW', rvar1=rshortup_top_clr)

    if (allocated(rlong_clr))        call increment_vtable('RLONG_CLR',       'AW', rvar1=rlong_clr)

    if (allocated(rlongup_clr))      call increment_vtable('RLONGUP_CLR',     'AW', rvar1=rlongup_clr)

    if (allocated(rlongup_top_clr))  call increment_vtable('RLONGUP_TOP_CLR', 'AW', rvar1=rlongup_top_clr)

    if (allocated(pbl_cld_forc))     call increment_vtable('PBL_CLD_FORC',    'AW', rvar1=pbl_cld_forc)

    if (allocated(par))              call increment_vtable('PAR_TOT',         'AW', rvar1=par)

    if (allocated(par_diffuse))      call increment_vtable('PAR_DIFFUSE',     'AW', rvar1=par_diffuse)

    if (allocated(ppfd))             call increment_vtable('PPFD_TOT',        'AW', rvar1=ppfd)

    if (allocated(ppfd_diffuse))     call increment_vtable('PPFD_DIFFUSE',    'AW', rvar1=ppfd_diffuse)

    if (allocated(mcica_seed))       call increment_vtable('MCICA_SEED',      'AW', ivar2=mcica_seed)

  end subroutine filltab_radiate

End Module mem_radiate
