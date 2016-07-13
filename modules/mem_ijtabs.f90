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
Module mem_ijtabs

   use max_dims, only: maxgrds, maxremote
   implicit none
   
   private :: maxgrds, maxremote

   integer, parameter :: mloops = 7 ! max # non-para DO loops for M,V,W pts

   integer, parameter :: nloops_m = mloops + maxremote ! # M DO loops incl para
   integer, parameter :: nloops_v = mloops + maxremote ! # V DO loops incl para
   integer, parameter :: nloops_w = mloops + maxremote ! # W DO loops incl para

   ! M, U, V, and W loop indices. The last 4 letters of these indices have the
   ! following meanings:
   !
   ! *_grid is for setting up grid parameters in the ctrlvols subroutines
   ! *_init is for initialization of ATM fields (in ohhi, olhi, fldsisan, 
   !        hurricane_init, omic_init)
   ! *_prog is for points where primary quantities such as VMC or THIL are
   !        prognosed
   ! *_wadj is for U, V, or W points that are adjacent to W points where primary
   !        quantities are prognosed
   ! *_wstn is for all U, V, or W points in the stencil of a W point where
   !        primary quantities are prognosed
   ! *_lbcp is for lateral boundary points whose values are copied from an
   !        interior prognostic point
   ! *_vadj is for M points that are adjacent to a prognostic V point, while
   !        jtw_vadj is the representation of jtm_vadj in the pre-hex-grid stage
   !        of model grid initialization)
   ! *_wall is for U or V points along the wall of an x-z channel domain (mdomain = 3)

   integer, parameter :: jtm_grid = 1, jtu_grid = 1, jtv_grid = 1, jtw_grid = 1
   integer, parameter :: jtm_init = 2, jtu_init = 2, jtv_init = 2, jtw_init = 2
   integer, parameter :: jtm_prog = 3, jtu_prog = 3, jtv_prog = 3, jtw_prog = 3
   integer, parameter :: jtm_wadj = 4, jtu_wadj = 4, jtv_wadj = 4, jtw_wadj = 4
   integer, parameter :: jtm_wstn = 5, jtu_wstn = 5, jtv_wstn = 5, jtw_wstn = 5
   integer, parameter :: jtm_lbcp = 6, jtu_lbcp = 6, jtv_lbcp = 6, jtw_lbcp = 6
   integer, parameter :: jtm_vadj = 7, jtu_wall = 7, jtv_wall = 7, jtw_vadj = 7

   integer :: nstp  ! # of finest grid acoustic timesteps in coarse grid dtlong
   integer :: istp  ! Current timestep counter from 1 to nstp
   integer :: mrls  ! Number of active mesh refinement levels (MRLs)

   integer, allocatable :: mrl_begl(:)  ! MRL at beginning of long timestep
   integer, allocatable :: mrl_begr(:)  ! MRL at beginning of RK step
   integer, allocatable :: mrl_begs(:)  ! MRL at beginning of short timestep
   integer, allocatable :: mrl_ends(:)  ! MRL at end of short timestep
   integer, allocatable :: mrl_endr(:)  ! MRL at end of RK step
   integer, allocatable :: mrl_endl(:)  ! MRL at end of long timestep
   real,    allocatable :: dtrk    (:)  ! MRL RK timestep factor

   integer, allocatable :: leafstep(:)  ! flag to run leaf on any sub-timestep

   Type itab_m_vars             ! data structure for M pts (individual rank)
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this M pt

      integer :: npoly = 0       ! number of V/W neighbors of this M pt
      integer :: imp = 1         ! M point from which to copy this M pt's values
      integer :: imglobe = 1     ! global index of this M pt (in parallel case)
      integer :: mrlm = 0        ! mesh refinement level of this M pt
      integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
      integer :: mrow = 0        ! Full row number outside nest
      integer :: ngr = 0         ! Grid number
      integer :: iv(3) = 1       ! array of V neighbors of this M pt
      integer :: iw(3) = 1       ! array of W neighbors of this M pt
   End Type itab_m_vars

   Type itab_v_vars             ! data structure for V pts (individual rank)
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this V pt

      integer :: ivp = 1       ! V pt from which to copy this V pt's values
      integer :: irank = -1    ! rank of parallel process at this V pt
      integer :: ivglobe = 1   ! global index of this V pt (in parallel case)
      integer :: mrlv = 0      ! mesh refinement level of this V pt
      integer :: im(6) = 1     ! neighbor M pts of this V pt
      integer :: iw(4) = 1     ! neighbor W pts of this V pt
      integer :: iv(4) = 1     ! neighbor V pts

      real :: farw(2) = 0.     ! Interp of ARW to V control volume [VADV + VDIFF]

      real :: cosv(2) = 0.     ! cosine of angle between V and zonal dir (Voronoi)
      real :: sinv(2) = 0.     ! sine of angle between V and zonal dir (Voronoi)

      real :: dxps(2) = 0.     ! xps (eastward) displacement from neighbor W pts
      real :: dyps(2) = 0.     ! yps (northward) displacement from neighbor W pts
   End Type itab_v_vars

   Type itab_w_vars             ! data structure for W pts (individual rank)
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this W pt

      integer :: npoly = 0     ! number of M/V neighbors of this W pt
      integer :: iwp = 1       ! W pt from which to copy this W pt's values
      integer :: irank = -1    ! rank of parallel process at this W pt
      integer :: iwglobe = 1   ! global index of this W pt (in parallel run)
      integer :: mrlw = 0      ! mesh refinement level of this W pt
      integer :: mrlw_orig = 0 ! original MRL of this W pt
      integer :: ngr = 0       ! Grid number
      integer :: im(7) = 1     ! neighbor M pts
      integer :: iv(7) = 1     ! neighbor V pts
      integer :: iw(7) = 1     ! neighbor W pts

      real :: dirv(7) = 0.     ! pos direction of V neighbors

      real :: farm(7) = 0.     ! Fraction of arw0 in each M point sector 
      real :: farv(7) = 0.     ! Fraction of arw0 in each V point sector 

      real :: gxps1(7) = 0.    ! gradient weight xe component for point 1
      real :: gyps1(7) = 0.    ! gradient weight ye component for point 1

      real :: gxps2(7) = 0.    ! gradient weight xe component for point 2
      real :: gyps2(7) = 0.    ! gradient weight ye component for point 2

      real :: unx_w = 0.       ! xe component of eastward unit normal vector
      real :: uny_w = 0.       ! ye component of eastward unit normal vector

      real :: vnx_w = 0.       ! xe component of northward unit normal vector
      real :: vny_w = 0.       ! ye component of northward unit normal vector
      real :: vnz_w = 0.       ! ze component of northward unit normal vector

      real :: ecvec_vx(7) = 0. ! factors converting V to earth cart. velocity
      real :: ecvec_vy(7) = 0. ! factors converting V to earth cart. velocity
      real :: ecvec_vz(7) = 0. ! factors converting V to earth cart. velocity

      integer :: iwnud(3) = 1  ! local nudpoly pts
      real    :: fnud(3) = 0.  ! local nudpoly coeffs

      integer :: nland         ! number of land cells attached to this W column
      integer :: nsea          ! number of sea cells attached to this W column

      integer, allocatable :: iland(:) ! indices of attached land cells
      integer, allocatable :: isea(:)  ! indices of attached sea cells
   End Type itab_w_vars

   Type itab_md_vars             ! data structure for M pts (individual rank)
                                 ! on the Delaunay mesh
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this M pt

      integer :: npoly = 0       ! number of V/W neighbors of this M pt
      integer :: imp = 1         ! M point from which to copy this M pt's values
      integer :: mrlm = 0        ! mesh refinement level of this M pt
      integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
      integer :: ngr = 0         ! Grid number
      integer :: im(7) = 1       ! array of M neighbors of this M pt
      integer :: iu(7) = 1       ! array of U neighbors of this M pt
      integer :: iw(7) = 1       ! array of W neighbors of this M pt
   End Type itab_md_vars

   Type itab_ud_vars             ! data structure for U pts (individual rank)
                                 ! on the Delaunay mesh
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this M pt

      integer :: iup = 1       ! U pt from which to copy this U pt's values
      integer :: mrlu = 0      ! mesh refinement level of this U pt
      integer :: im(2) = 1     ! neighbor M pts of this U pt
      integer :: iu(12) = 1    ! neighbor U pts
      integer :: iw(6) = 1     ! neighbor W pts
   End Type itab_ud_vars

   Type itab_wd_vars             ! data structure for W pts (individual rank)
                                 ! on the Delaunay mesh
      logical, allocatable :: loop(:) ! flag to perform each DO loop at this W pt

      integer :: npoly = 0     ! number of M/V neighbors of this W pt
      integer :: iwp = 1       ! W pt from which to copy this W pt's values
      integer :: mrlw = 0      ! mesh refinement level of this W pt
      integer :: mrlw_orig = 0 ! original MRL of this W pt
      integer :: mrow = 0      ! Full row number outside nest
      integer :: ngr = 0       ! Grid number
      integer :: im(3) = 1     ! neighbor M pts
      integer :: iu(3) = 1     ! neighbor U pts
      integer :: iw(9) = 1     ! neighbor W pts
   End Type itab_wd_vars

   Type itabg_m_vars            ! data structure for M pts (global)
      integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
      integer :: irank = -1     ! rank of parallel process at this M pt
      integer :: im_myrank_imp = -1 ! local M point that corresponds to global imp
   End Type itabg_m_vars

   Type itabg_v_vars            ! data structure for V pts (global)
      integer :: iv_myrank = -1 ! local (parallel subdomain) index of this V pt
      integer :: irank = -1     ! rank of parallel process at this V pt
      integer :: iv_myrank_ivp = -1 ! local V point that corresponds to global ivp
   End Type itabg_v_vars

   Type itabg_w_vars            ! data structure for W pts (global)
      integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
      integer :: irank = -1     ! rank of parallel process at this W pt
      integer :: iw_myrank_iwp = -1 ! local W point that corresponds to global iwp
   End Type itabg_w_vars

   Type nest_ud_vars        ! temporary U-pt data structure for spawning nested grids
      integer :: im=0, iu=0 ! new M/U pts attached to this U pt 
   End Type nest_ud_vars
   
   Type nest_wd_vars        ! temporary W-pt data structure for spawning nested grids
      integer :: iu(3) = 0  ! new U pts attached to this W pt
      integer :: iw(3) = 0  ! new W pts attached to this W pt
   End Type nest_wd_vars

   Type jtab_m_vars
      integer, allocatable :: im(:)
      integer, allocatable :: jend(:)
   End Type jtab_m_vars

   Type jtab_v_vars
      integer, allocatable :: iv(:)
      integer, allocatable :: jend(:)
   End Type jtab_v_vars

   Type jtab_w_vars
      integer, allocatable :: iw(:)
      integer, allocatable :: jend(:)
   End Type jtab_w_vars

   Type itab_m_pd_vars      ! data structure for M pts (individual rank) on para_(decomp,init)
      integer :: npoly = 0  ! number of V/W neighbors of this M pt
      integer :: imp   = 1  ! M point from which to copy this M pt's value
      integer :: iv(3) = 1  ! array of V neighbors of this M pt
      integer :: iw(3) = 1  ! array of W neighbors of this M pt
   End Type itab_m_pd_vars

   Type itab_v_pd_vars      ! data structure for V pts (individual rank) on para_(decomp,init)
      integer :: ivp    = 1 ! V pt from which to copy this V pt's values
      integer :: im(6)  = 1 ! neighbor M pts of this V pt
      integer :: iv(4)  = 1 ! neighbor V pts
      integer :: iw(4)  = 1 ! neighbor W pts of this V pt
   End Type itab_v_pd_vars

   Type itab_w_pd_vars      ! data structure for W pts (individual rank) on para_(decomp,init)
      integer :: iwp   = 1  ! W pt from which to copy this W pt's values
      integer :: npoly = 0  ! number of M/V neighbors of this W pt
      integer :: im(7) = 1  ! neighbor M pts of this W pt
      integer :: iv(7) = 1  ! neighbor V pts
      integer :: iw(7) = 1  ! neighbor W pts
      integer :: iwnud(3) = 1  ! local nudpoly pts
   End Type itab_w_pd_vars

   type (itab_md_vars), allocatable :: itab_md(:)
   type (itab_ud_vars), allocatable :: itab_ud(:)
   type (itab_wd_vars), allocatable :: itab_wd(:)

   type (itab_m_vars),  allocatable, target :: itab_m(:)
   type (itab_v_vars),  allocatable, target :: itab_v(:)
   type (itab_w_vars),  allocatable, target :: itab_w(:)

   type (itab_m_pd_vars), allocatable, target :: itab_m_pd(:)
   type (itab_v_pd_vars), allocatable, target :: itab_v_pd(:)
   type (itab_w_pd_vars), allocatable, target :: itab_w_pd(:)

   type (itab_md_vars), allocatable :: ltab_md(:)
   type (itab_ud_vars), allocatable :: ltab_ud(:)
   type (itab_wd_vars), allocatable :: ltab_wd(:)

   type (itabg_m_vars), allocatable, target :: itabg_m(:)
   type (itabg_v_vars), allocatable, target :: itabg_v(:)
   type (itabg_w_vars), allocatable, target :: itabg_w(:)

   type (nest_ud_vars), allocatable :: nest_ud(:)
   type (nest_wd_vars), allocatable :: nest_wd(:)
   
   type (jtab_m_vars) :: jtab_m(nloops_m)
   type (jtab_v_vars) :: jtab_v(nloops_v)
   type (jtab_w_vars) :: jtab_w(nloops_w)

Contains

!===============================================================================

   subroutine alloc_itabsd(mma, mua, mwa)

   implicit none

   integer, intent(in) :: mma, mua, mwa
   integer :: imd, iud, iwd

   allocate (itab_md(mma))
   allocate (itab_ud(mua))
   allocate (itab_wd(mwa))

   do imd = 1,mma
      allocate(itab_md(imd)%loop(mloops))
      itab_md(imd)%loop(1:mloops) = .false.
   enddo

   do iud = 1,mua
      allocate(itab_ud(iud)%loop(mloops))
      itab_ud(iud)%loop(1:mloops) = .false.
   enddo

   do iwd = 1,mwa
      allocate(itab_wd(iwd)%loop(mloops))
      itab_wd(iwd)%loop(1:mloops) = .false.
   enddo

   end subroutine alloc_itabsd

!===============================================================================

   subroutine alloc_itabs(mma, mva, mwa, input)

   implicit none

   integer, intent(in) :: mma, mva, mwa, input
   integer :: im, iv, iw

   allocate (itab_m(mma))
   allocate (itab_v(mva))
   allocate (itab_w(mwa))

   if (input == 0) then

      do im = 1,mma
         allocate(itab_m(im)%loop(mloops))
         itab_m(im)%loop(1:mloops) = .false.
      enddo

      do iv = 1,mva
         allocate(itab_v(iv)%loop(mloops))
         itab_v(iv)%loop(1:mloops) = .false.
      enddo

      do iw = 1,mwa
         allocate(itab_w(iw)%loop(mloops))
         itab_w(iw)%loop(1:mloops) = .false.
      enddo

   else

      do im = 1,mma
         allocate(itab_m(im)%loop(nloops_m))
         itab_m(im)%loop(1:nloops_m) = .false.
      enddo

      do iv = 1,mva
         allocate(itab_v(iv)%loop(nloops_v))
         itab_v(iv)%loop(1:nloops_v) = .false.
      enddo

      do iw = 1,mwa
         allocate(itab_w(iw)%loop(nloops_w))
         itab_w(iw)%loop(1:nloops_w) = .false.
      enddo

   endif

   end subroutine alloc_itabs

!===============================================================================

   subroutine alloc_itabs_pd(mma, mva, mwa)

   implicit none

   integer, intent(in) :: mma, mva, mwa

   allocate (itab_m_pd(mma))
   allocate (itab_v_pd(mva))
   allocate (itab_w_pd(mwa))

   end subroutine alloc_itabs_pd

!===============================================================================

   subroutine fill_jtabs(mma, mva, mwa, input)

   implicit none

   integer, intent(in) :: mma, mva, mwa, input

   integer :: iw, iv, im, mrl
   integer :: iloop, jend
   integer :: nlm, nlv, nlw

   if (input == 0) then
      nlm = mloops
      nlv = mloops
      nlw = mloops
   else
      nlm = nloops_m
      nlv = nloops_v
      nlw = nloops_w
   endif

! Allocate and zero-fill jtab%jend()

   do iloop = 1,nlm
      allocate (jtab_m(iloop)%jend(mrls))
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo
   
   if (allocated(itab_v)) then 
      do iloop = 1,nlv
         allocate (jtab_v(iloop)%jend(mrls))
         jtab_v(iloop)%jend(1:mrls) = 0
      enddo
   endif
   
   do iloop = 1,nlw
      allocate (jtab_w(iloop)%jend(mrls))
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! Compute and store jtab%jend(1)

   do iloop = 1,nlm
      jtab_m(iloop)%jend(1) = 0
      do im = 2,mma
         if (itab_m(im)%loop(iloop)) then
            jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
         endif
      enddo
      jtab_m(iloop)%jend(1) = max(1,jtab_m(iloop)%jend(1))
   enddo

   do iloop = 1,nlv
      jtab_v(iloop)%jend(1) = 0
      do iv = 2,mva
         if (itab_v(iv)%loop(iloop)) then
            jtab_v(iloop)%jend(1) = jtab_v(iloop)%jend(1) + 1
         endif
      enddo
      jtab_v(iloop)%jend(1) = max(1,jtab_v(iloop)%jend(1))
   enddo
   
   do iloop = 1,nlw
      jtab_w(iloop)%jend(1) = 0
      do iw = 2,mwa
         if (itab_w(iw)%loop(iloop)) then
            jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
         endif
      enddo
      jtab_w(iloop)%jend(1) = max(1,jtab_w(iloop)%jend(1))
   enddo

! Allocate and zero-fill JTAB_M%IM, JTAB_V%IV, JTAB_W%IW

   do iloop = 1,nlm
      jend = jtab_m(iloop)%jend(1)
      allocate (jtab_m(iloop)%im(jend))
      jtab_m(iloop)%im(1:jend) = 0
   enddo

   do iloop = 1,nlv
      jend = jtab_v(iloop)%jend(1)
      allocate (jtab_v(iloop)%iv(jend))
      jtab_v(iloop)%iv(1:jend) = 0
   enddo
   
   do iloop = 1,nlw
      jend = jtab_w(iloop)%jend(1)
      allocate (jtab_w(iloop)%iw(jend))
      jtab_w(iloop)%iw(1:jend) = 0
   enddo

! Initialize JTAB%JEND counters to zero

   do iloop = 1,nlm
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo

   do iloop = 1,nlv
      jtab_v(iloop)%jend(1:mrls) = 0
   enddo
   
   do iloop = 1,nlw
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_M%IM
! ///////////////////////////////////////////////////////////////////////////

   do mrl = mrls,1,-1
      do im = 2,mma
         do iloop = 1,nlm
            if (itab_m(im)%loop(iloop) .and. itab_m(im)%mrlm == mrl) then
               jtab_m(iloop)%jend(1:mrl) = jtab_m(iloop)%jend(1:mrl) + 1
               jtab_m(iloop)%im(jtab_m(iloop)%jend(1)) = im
            endif
         enddo
      enddo
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_V%IV
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   do mrl = mrls,1,-1
      do iv = 2,mva
         do iloop = 1,nlv
            if (itab_v(iv)%loop(iloop) .and. itab_v(iv)%mrlv == mrl) then
               jtab_v(iloop)%jend(1:mrl) = jtab_v(iloop)%jend(1:mrl) + 1
               jtab_v(iloop)%iv(jtab_v(iloop)%jend(1)) = iv
            endif
         enddo
      enddo
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_W%IW
! ///////////////////////////////////////////////////////////////////////////

   do mrl = mrls,1,-1
      do iw = 2,mwa
         do iloop = 1,nlw
            if (itab_w(iw)%loop(iloop) .and. itab_w(iw)%mrlw == mrl) then
               jtab_w(iloop)%jend(1:mrl) = jtab_w(iloop)%jend(1:mrl) + 1
               jtab_w(iloop)%iw(jtab_w(iloop)%jend(1)) = iw
            endif
         enddo
      enddo
   enddo

   end subroutine fill_jtabs

End Module mem_ijtabs

