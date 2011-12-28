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
Module mem_ijtabs

   use max_dims, only: maxgrds, maxremote
   implicit none
   
   private :: maxgrds, maxremote

   integer, parameter :: mloops_m =  4 ! max # non-para DO loops for M pts
   integer, parameter :: mloops_u = 25 ! max # non-para DO loops for U pts
   integer, parameter :: mloops_v = 25 ! max # non-para DO loops for V pts
   integer, parameter :: mloops_w = 35 ! max # non-para DO loops for W pts

   integer, parameter :: nloops_m = mloops_m             ! (no para loops for M)
   integer, parameter :: nloops_u = mloops_u + maxremote ! # U DO loops incl para
   integer, parameter :: nloops_v = mloops_v + maxremote ! # V DO loops incl para
   integer, parameter :: nloops_w = mloops_w + maxremote ! # W DO loops incl para
   
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
      logical :: loop(nloops_m) = .false. ! flag to perform each DO loop at this M pt

      integer :: npoly = 0       ! number of V/W neighbors of this M pt
      integer :: itopm = 1       ! M point from which to copy this M pt's topo
      integer :: imglobe = 1     ! global index of this M pt (in parallel case)
      integer :: mrlm = 0        ! mesh refinement level of this M pt
      integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
      integer :: mrow=0, mrowh=0 ! Full and half row number outside nest
      integer :: iu(7) = 1       ! array of U neighbors of this M pt (Delaunay)
      integer :: iv(7) = 1       ! array of V neighbors of this M pt (Voronoi)
      integer :: iw(7) = 1       ! array of W neighbors of this M pt (Del or Vor)
      real    :: fmw(7) = 0.     ! Interp coefs of W values to M point
   End Type itab_m_vars

   Type itab_u_vars             ! data structure for U pts (individual rank)
      logical :: loop(nloops_u) = .false. ! flag to perform each DO loop at this M pt

      integer :: iup=1         ! U pt from which to copy this U pt's values
      integer :: irank = -1    ! rank of parallel process at this U pt
      integer :: iuglobe = 1   ! global index of this U pt (in parallel case)
      integer :: mrlu = 0      ! mesh refinement level of this U pt
      integer :: im(2) = 1     ! neighbor M pts of this U pt
      integer :: iu(12) = 1    ! neighbor U pts
      integer :: iw(6) = 1     ! neighbor W pts
      real    :: diru(4) = 0.  ! pos direction of U neighbors

! Some of the following will change when V is routinely diagnosed on TRI grid

      real :: fuu(12) = 0.   ! proj coefs of U neighbors
      real :: fuw(6) = 0.    ! proj coefs of W neighbors
      real :: tuu(4) = 0.    ! proj coefs of U neighbors
      real :: guw(4)  = 0.   ! gradient coefs of outer W neighbors of U
      real :: gcf36 = 0., gcf45 = 0. ! UMC grad coefs for intrp to U edge center
      real :: pgc12=0., pgc45=0., pgc63=0., pgc12b=0.    ! PGF proj coefs
      real :: pgc45b=0., pgc12c=0., pgc63c=0., pgc12d=0. ! PGF proj coefs
      real :: vxu_u(4) = 0.  ! proj coefs
      real :: vyu_u(4) = 0.  ! proj coefs
      real :: vxw_u(2) = 0.  ! proj coefs
      real :: vyw_u(2) = 0.  ! proj coefs
      real :: crossmm = 0.
      real :: crossww = 0.
   End Type itab_u_vars

   Type itab_v_vars             ! data structure for V pts (individual rank)
      logical :: loop(nloops_v) = .false. ! flag to perform each DO loop at this V pt

      integer :: ivp=1         ! V pt from which to copy this V pt's values
      integer :: irank = -1    ! rank of parallel process at this V pt
      integer :: ivglobe = 1   ! global index of this V pt (in parallel case)
      integer :: mrlv=0        ! mesh refinement level of this V pt
      integer :: im(6)=1       ! neighbor M pts of this V pt
      integer :: iv(16)=1      ! neighbor V pts
      integer :: iw(4)=1       ! neighbor W pts of this V pt

      real :: fvv(12)=0.       ! Proj of V5-12 onto V    [HADV + HDIFF]
      real :: fvw (4)=0.       ! Proj of W1-4 onto V     [HADV + HDIFF]
      real :: fuv(16)=0.       ! Interp of V1-16 onto U  [CORF]
      real :: farw(2)=0.       ! Interp of ARW to V control volume [VADV + VDIFF]

      real :: cosv(2)=0.       ! cosine of angle between V and zonal dir (Voronoi)
      real :: sinv(2)=0.       ! sine of angle between V and zonal dir (Voronoi)

      real :: dxps(2)=0.       ! xps (eastward) displacement from neighbor W pts
      real :: dyps(2)=0.       ! yps (northward) displacement from neighbor W pts
   End Type itab_v_vars

   Type itab_w_vars             ! data structure for W pts (individual rank)
      logical :: loop(nloops_w) = .false. ! flag to perform each DO loop at this W pt

      integer :: npoly = 0   ! number of M/V neighbors of this W pt
      integer :: iwp=1       ! W pt from which to copy this W pt's values
      integer :: irank = -1  ! rank of parallel process at this W pt
      integer :: iwglobe = 1 ! global index of this W pt (in parallel run)
      integer :: mrlw=0      ! mesh refinement level of this W pt
      integer :: mrlw_orig=0 ! mesh refinement level of this W pt
      integer :: mrow=0      ! Full row number outside nest (Delaunay)
      integer :: mrowh=0     ! Half row number outside nest (Delaunay)
      integer :: im(7)=1     ! neighbor M pts of this W pt
      integer :: iu(9)=1     ! neighbor U pts (9 Delaunay, 7 Voronoi)
      integer :: iv(7)=1     ! neighbor V pts (Voronoi)
      integer :: iw(9)=1     ! neighbor W pts (9 Delaunay, 7 Voronoi)

      real :: diru(3)=0.     ! pos direction of U neighbors (Delaunay)
      real :: dirv(7)=0.     ! pos direction of V neighbors (Voronoi)
      real :: fwv (7)=0.     ! Proj of V1-7 onto W [HADV + HDIFF]
      real :: fww (7)=0.     ! Proj of W1-7 onto W [HADV + HDIFF]
      
      real :: farm(7)=0.     ! Fraction of arw0 in each M point sector 
      real :: farv(7)=0.     ! Fraction of arw0 in each V point sector 

      real :: gxps1(7)=0.    ! gradient weight xe component for point 1
      real :: gyps1(7)=0.    ! gradient weight ye component for point 1

      real :: gxps2(7)=0.    ! gradient weight xe component for point 2
      real :: gyps2(7)=0.    ! gradient weight ye component for point 2

      real :: unx_w=0.       ! xe component of eastward unit normal vector
      real :: uny_w=0.       ! ye component of eastward unit normal vector

      real :: vnx_w=0.       ! xe component of northward unit normal vector
      real :: vny_w=0.       ! ye component of northward unit normal vector
      real :: vnz_w=0.       ! ze component of northward unit normal vector

! Delaunay section - will change some

      real :: fwu(9) = 0.  ! proj coefs of U neighbors
!      real :: fww(3) = 0.  ! proj coefs of W neighbors
      real :: vxu(3) = 0.  ! proj coefs of U neighbors
      real :: vyu(3) = 0.  ! proj coefs of U neighbors
      real :: vzu(3) = 0.  ! proj coefs of U neighbors
      real :: vxw = 0.     ! proj coefs of W neighbors
      real :: vyw = 0.     ! proj coefs of W neighbors
      real :: vzw = 0.     ! proj coefs of W neighbors
      real :: vxu_w(3) = 0.  ! proj coefs of U neighbors
      real :: vyu_w(3) = 0.  ! proj coefs of U neighbors

      integer :: iwnud(3) = 1  ! local nudpoly pts
      real    :: fnud(3) = 0. ! local nudpoly coeffs
   End Type itab_w_vars

   Type itabg_m_vars            ! data structure for M pts (global)
      integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
      integer :: irank = -1     ! rank of parallel process at this M pt
   End Type itabg_m_vars

   Type itabg_u_vars            ! data structure for U pts (global)
      integer :: iu_myrank = -1 ! local (parallel subdomain) index of this U pt
      integer :: irank = -1     ! rank of parallel process at this U pt
   End Type itabg_u_vars

   Type itabg_v_vars            ! data structure for V pts (global)
      integer :: iv_myrank = -1 ! local (parallel subdomain) index of this V pt
      integer :: irank = -1     ! rank of parallel process at this V pt
   End Type itabg_v_vars

   Type itabg_w_vars            ! data structure for W pts (global)
      integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
      integer :: irank = -1     ! rank of parallel process at this W pt
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

   Type jtab_u_vars
      integer, allocatable :: iu(:)
      integer, allocatable :: jend(:)
   End Type jtab_u_vars

   Type jtab_v_vars
      integer, allocatable :: iv(:)
      integer, allocatable :: jend(:)
   End Type jtab_v_vars

   Type jtab_w_vars
      integer, allocatable :: iw(:)
      integer, allocatable :: jend(:)
   End Type jtab_w_vars

   Type itab_m_pd_vars         ! data structure for M pts (individual rank) on para_(decomp,init)
      integer :: npoly = 0     ! number of V/W neighbors of this M pt
      integer :: itopm = 1     ! M point from which to copy this M pt's topo
      integer :: iu(7) = 1     ! array of U neighbors of this M pt (Delaunay)
      integer :: iv(7) = 1     ! array of V neighbors of this M pt (Voronoi)
      integer :: iw(7) = 1     ! array of W neighbors of this M pt (Del or Vor)
   End Type itab_m_pd_vars

   Type itab_u_pd_vars         ! data structure for U pts (individual rank) on para_(decomp,init)
      integer :: iup    = 1    ! U pt from which to copy this U pt's values
      integer :: im(2)  = 1    ! neighbor M pts of this U pt
      integer :: iu(12) = 1    ! neighbor U pts
      integer :: iw(6)  = 1    ! neighbor W pts
   End Type itab_u_pd_vars

   Type itab_v_pd_vars         ! data structure for V pts (individual rank) on para_(decomp,init)
      integer :: ivp    = 1    ! V pt from which to copy this V pt's values
      integer :: im(6)  = 1    ! neighbor M pts of this V pt
      integer :: iv(16) = 1    ! neighbor V pts
      integer :: iw(4)  = 1    ! neighbor W pts of this V pt
   End Type itab_v_pd_vars

   Type itab_w_pd_vars         ! data structure for W pts (individual rank) on para_(decomp,init)
      integer :: iwp   = 1     ! W pt from which to copy this W pt's values
      integer :: npoly = 0     ! number of M/V neighbors of this W pt
      integer :: im(7) = 1     ! neighbor M pts of this W pt
      integer :: iu(9) = 1     ! neighbor U pts (9 Delaunay, 7 Voronoi)
      integer :: iv(7) = 1     ! neighbor V pts (Voronoi)
      integer :: iw(9) = 1     ! neighbor W pts (9 Delaunay, 7 Voronoi)
   End Type itab_w_pd_vars

   type (itab_m_vars), allocatable :: itab_md(:)
   type (itab_u_vars), allocatable :: itab_ud(:)
   type (itab_w_vars), allocatable :: itab_wd(:)

   type (itab_m_vars),  allocatable, target :: itab_m(:)
   type (itab_u_vars),  allocatable, target :: itab_u(:)
   type (itab_v_vars),  allocatable, target :: itab_v(:)
   type (itab_w_vars),  allocatable, target :: itab_w(:)

   type (itab_m_pd_vars), allocatable, target :: itab_m_pd(:)
   type (itab_u_pd_vars), allocatable, target :: itab_u_pd(:)
   type (itab_v_pd_vars), allocatable, target :: itab_v_pd(:)
   type (itab_w_pd_vars), allocatable, target :: itab_w_pd(:)

   type (itab_m_vars), allocatable :: ltab_md(:)
   type (itab_u_vars), allocatable :: ltab_ud(:)
   type (itab_w_vars), allocatable :: ltab_wd(:)

   type (itabg_m_vars), allocatable, target :: itabg_m(:)
   type (itabg_u_vars), allocatable, target :: itabg_u(:)
   type (itabg_v_vars), allocatable, target :: itabg_v(:)
   type (itabg_w_vars), allocatable, target :: itabg_w(:)

   type (nest_ud_vars), allocatable :: nest_ud(:)
   type (nest_wd_vars), allocatable :: nest_wd(:)
   
   type (jtab_m_vars) :: jtab_m(nloops_m)
   type (jtab_u_vars) :: jtab_u(nloops_u)
   type (jtab_v_vars) :: jtab_v(nloops_v)
   type (jtab_w_vars) :: jtab_w(nloops_w)

Contains

!===============================================================================

   subroutine alloc_itabsd(mma,mua,mwa)

   implicit none

   integer, intent(in) :: mma,mua,mwa

   allocate (itab_md(mma))
   allocate (itab_ud(mua))
   allocate (itab_wd(mwa))

   return
   end subroutine alloc_itabsd

!===============================================================================

   subroutine alloc_itabs(meshtype,mma,mua,mva,mwa)

   implicit none

  integer, intent(in) :: meshtype,mma,mua,mva,mwa

   if (meshtype == 1) then
      allocate (itab_u(mua))
   elseif (meshtype == 2) then
      allocate (itab_v(mva))
   endif
   
   allocate (itab_m(mma))
   allocate (itab_w(mwa))

   return
   end subroutine alloc_itabs

!===============================================================================

   subroutine alloc_itabs_pd(meshtype,mma,mua,mva,mwa)

   implicit none

  integer, intent(in) :: meshtype,mma,mua,mva,mwa

  allocate (itab_m_pd(mma))
   if (meshtype == 1) then
      allocate (itab_u_pd(mua))
   elseif (meshtype == 2) then
      allocate (itab_v_pd(mva))
   endif
   
   allocate (itab_w_pd(mwa))

   return
   end subroutine alloc_itabs_pd

!===============================================================================

   subroutine filltab_itabs()

     use var_tables, only: vtab_r, num_var, increment_vtable
     use misc_coms,  only: iparallel, runtype, meshtype

     implicit none

!  THESE ONLY NEED TO BE WRITTEN TO HISTORY FILE FOR PARALLEL RUNS,
!  AND READ FOR PLOTONLY OR PARCOMBINE RUNS.

   if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

      if (allocated(itab_u)) then

         call increment_vtable('IUGLOBE', 'AU', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_u%iuglobe
            
         call increment_vtable('IRANKU', 'AU', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_u%irank
            
      endif

      if (allocated(itab_v)) then
            
         call increment_vtable('IVGLOBE', 'AU', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_v%ivglobe
            
         call increment_vtable('IRANKV', 'AU', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_v%irank

      endif

      if (allocated(itab_m)) then

         call increment_vtable('IMGLOBE', 'AM', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_m%imglobe

      endif

      if (allocated(itab_w)) then

         call increment_vtable('IWGLOBE', 'AW', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_w%iwglobe

         call increment_vtable('IRANKW', 'AW', noread=.true.)
         vtab_r(num_var)%ivar1_p => itab_w%irank

      endif

   endif

   end subroutine filltab_itabs

!===============================================================================

   subroutine fill_jtabs(mma,mua,mva,mwa)

   use misc_coms,  only: io6, nqparm

   implicit none

   integer, intent(in) :: mma,mua,mva,mwa

   integer :: iw,iu,iv,im,k,nl,mrl
   integer :: iloop,iw1,iw2,mrl0,jend

! Allocate and zero-fill jtab%jend()

   do iloop = 1,nloops_m
      allocate (jtab_m(iloop)%jend(mrls))
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo
   
   if (allocated(itab_u)) then 
      do iloop = 1,nloops_u
         allocate (jtab_u(iloop)%jend(mrls))
         jtab_u(iloop)%jend(1:mrls) = 0
      enddo
   endif

   if (allocated(itab_v)) then 
      do iloop = 1,nloops_v
         allocate (jtab_v(iloop)%jend(mrls))
         jtab_v(iloop)%jend(1:mrls) = 0
      enddo
   endif
   
   do iloop = 1,nloops_w
      allocate (jtab_w(iloop)%jend(mrls))
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! Compute and store jtab%jend(1)

   do iloop = 1,nloops_m
      jtab_m(iloop)%jend(1) = 0
      do im = 2,mma
         if (itab_m(im)%loop(iloop)) then
            jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
         endif
      enddo
      jtab_m(iloop)%jend(1) = max(1,jtab_m(iloop)%jend(1))
   enddo

   if (allocated(itab_u)) then 
      do iloop = 1,nloops_u
         jtab_u(iloop)%jend(1) = 0
         do iu = 2,mua
            if (itab_u(iu)%loop(iloop)) then
               jtab_u(iloop)%jend(1) = jtab_u(iloop)%jend(1) + 1
            endif
         enddo
         jtab_u(iloop)%jend(1) = max(1,jtab_u(iloop)%jend(1))
      enddo
   endif

   if (allocated(itab_v)) then 
      do iloop = 1,nloops_v
         jtab_v(iloop)%jend(1) = 0
         do iv = 2,mva
            if (itab_v(iv)%loop(iloop)) then
               jtab_v(iloop)%jend(1) = jtab_v(iloop)%jend(1) + 1
            endif
         enddo
         jtab_v(iloop)%jend(1) = max(1,jtab_v(iloop)%jend(1))
      enddo
   endif
   
   do iloop = 1,nloops_w
      jtab_w(iloop)%jend(1) = 0
      do iw = 2,mwa
         if (itab_w(iw)%loop(iloop)) then
            jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
         endif
      enddo
      jtab_w(iloop)%jend(1) = max(1,jtab_w(iloop)%jend(1))
   enddo

! Allocate and zero-fill JTAB_M%IM, JTAB_V%IV, JTAB_W%IW

   do iloop = 1,nloops_m
      jend = jtab_m(iloop)%jend(1)
      allocate (jtab_m(iloop)%im(jend))
      jtab_m(iloop)%im(1:jend) = 0
   enddo

   if (allocated(itab_u)) then 
      do iloop = 1,nloops_u
         jend = jtab_u(iloop)%jend(1)
         allocate (jtab_u(iloop)%iu(jend))
         jtab_u(iloop)%iu(1:jend) = 0
      enddo
   endif

   if (allocated(itab_v)) then 
      do iloop = 1,nloops_v
         jend = jtab_v(iloop)%jend(1)
         allocate (jtab_v(iloop)%iv(jend))
         jtab_v(iloop)%iv(1:jend) = 0
      enddo
   endif
   
   do iloop = 1,nloops_w
      jend = jtab_w(iloop)%jend(1)
      allocate (jtab_w(iloop)%iw(jend))
      jtab_w(iloop)%iw(1:jend) = 0
   enddo

! Initialize JTAB%JEND counters to zero

   do iloop = 1,nloops_m
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo

   if (allocated(itab_u)) then 
      do iloop = 1,nloops_u
         jtab_u(iloop)%jend(1:mrls) = 0
      enddo
   endif

   if (allocated(itab_v)) then 
      do iloop = 1,nloops_v
         jtab_v(iloop)%jend(1:mrls) = 0
      enddo
   endif
   
   do iloop = 1,nloops_w
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_M%IM
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   do im = 2,mma
      do iloop = 1,2
         if (itab_m(im)%loop(iloop)) then
            jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
            jtab_m(iloop)%im(jtab_m(iloop)%jend(1)) = im
         endif
      enddo
   enddo
   
! MRL-dependent loops

   do mrl = mrls,1,-1
      do im = 2,mma
         do iloop = 3,nloops_m
            if (itab_m(im)%loop(iloop) .and. itab_m(im)%mrlm == mrl) then
               jtab_m(iloop)%jend(1:mrl) = jtab_m(iloop)%jend(1:mrl) + 1
               jtab_m(iloop)%im(jtab_m(iloop)%jend(1)) = im
            endif
         enddo
      enddo
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_U%IU
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   if (allocated(itab_u)) then 
      do iu = 2,mua
         do iloop = 1,10
            if (itab_u(iu)%loop(iloop)) then
               jtab_u(iloop)%jend(1) = jtab_u(iloop)%jend(1) + 1
               jtab_u(iloop)%iu(jtab_u(iloop)%jend(1)) = iu
            endif
         enddo      
      enddo

! MRL-dependent loops

      do mrl = mrls,1,-1
         do iu = 2,mua
            do iloop = 11,nloops_u
               if (itab_u(iu)%loop(iloop) .and. itab_u(iu)%mrlu == mrl) then
                  jtab_u(iloop)%jend(1:mrl) = jtab_u(iloop)%jend(1:mrl) + 1
                  jtab_u(iloop)%iu(jtab_u(iloop)%jend(1)) = iu
               endif
            enddo
         enddo
      enddo
   endif

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_V%IV
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   if (allocated(itab_v)) then 
      do iv = 2,mva
         do iloop = 1,10
            if (itab_v(iv)%loop(iloop)) then
               jtab_v(iloop)%jend(1) = jtab_v(iloop)%jend(1) + 1
               jtab_v(iloop)%iv(jtab_v(iloop)%jend(1)) = iv
            endif
         enddo      
      enddo

! MRL-dependent loops

      do mrl = mrls,1,-1
         do iv = 2,mva
            do iloop = 11,nloops_v
               if (itab_v(iv)%loop(iloop) .and. itab_v(iv)%mrlv == mrl) then
                  jtab_v(iloop)%jend(1:mrl) = jtab_v(iloop)%jend(1:mrl) + 1
                  jtab_v(iloop)%iv(jtab_v(iloop)%jend(1)) = iv
               endif
            enddo
         enddo
      enddo
   endif
   
! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_W%IW
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   do iw = 2,mwa
      do iloop = 1,10
         if (itab_w(iw)%loop(iloop)) then
            jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
            jtab_w(iloop)%iw(jtab_w(iloop)%jend(1)) = iw
         endif
      enddo      
   enddo

! MRL-dependent loops

   do mrl = mrls,1,-1
      do iw = 2,mwa
         do iloop = 11,nloops_w

            if (itab_w(iw)%loop(iloop) .and. itab_w(iw)%mrlw == mrl) then
               jtab_w(iloop)%jend(1:mrl) = jtab_w(iloop)%jend(1:mrl) + 1
               jtab_w(iloop)%iw(jtab_w(iloop)%jend(1)) = iw
            endif

         enddo
      enddo
   enddo

! NOTE:  For cumulus parameterization (W loop 15), MRL dependence of the
! parameterization computations is selected in the parameterization itself.
! Here, loop 15 MRL dependence is as for other loops in order to copy
! THSRC and RTSRC to tendency arrays at the proper time for each MRL.

   end subroutine fill_jtabs

End Module mem_ijtabs

