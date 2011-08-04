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

Module mem_nudge

  integer, allocatable         ::   itab_nudp(:,:)

  real,    allocatable, target ::      xenudp(:)
  real,    allocatable, target ::      yenudp(:)
  real,    allocatable, target ::      zenudp(:)

  real,    allocatable, target ::    rho_obsp(:,:)
  real,    allocatable, target ::  theta_obsp(:,:)
  real,    allocatable, target ::    shw_obsp(:,:)
  real,    allocatable, target :: uzonal_obsp(:,:)
  real,    allocatable, target :: umerid_obsp(:,:)

  real,    allocatable, target ::    rho_obsf(:,:)
  real,    allocatable, target ::  theta_obsf(:,:)
  real,    allocatable, target ::    shw_obsf(:,:)
  real,    allocatable, target :: uzonal_obsf(:,:)
  real,    allocatable, target :: umerid_obsf(:,:)

  real,    allocatable         ::    rho_obs(:,:)
  real,    allocatable         ::  theta_obs(:,:)
  real,    allocatable         ::    shw_obs(:,:)
  real,    allocatable         :: uzonal_obs(:,:)
  real,    allocatable         :: umerid_obs(:,:)

  real,    allocatable         ::    rho_sim(:,:)
  real,    allocatable         ::  theta_sim(:,:)
  real,    allocatable         ::    shw_sim(:,:)
  real,    allocatable         :: uzonal_sim(:,:)
  real,    allocatable         :: umerid_sim(:,:)

  integer :: nudflag
  integer :: nudnxp
  integer :: nnudfiles
  integer :: nnudfl
  integer :: nnudp
  integer :: mnudp

  real    :: tnudcent
  real    :: wt_nudge_uv
  real    :: wt_nudge_th
  real    :: wt_nudge_pi
  real    :: wt_nudge_rt
  real    :: wt_nudge_grid

Contains

  subroutine alloc_nudge1(lnudp)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: lnudp

!   Allocate arrays based on options (if necessary)

    write(io6,*) 'allocating nudge1 ', lnudp

    allocate (itab_nudp(lnudp,6)) ; itab_nudp = 0

    allocate (xenudp(lnudp)) ; xenudp = rinit
    allocate (yenudp(lnudp)) ; yenudp = rinit
    allocate (zenudp(lnudp)) ; zenudp = rinit

  end subroutine alloc_nudge1

!=========================================================================

  subroutine alloc_nudge2(mza)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza

!   Allocate arrays based on options (if necessary)

    write(io6,*) 'allocating nudge2 ',mza,mnudp

    allocate (   rho_obsp(mza,mnudp)) ; rho_obsp    = rinit
    allocate ( theta_obsp(mza,mnudp)) ; theta_obsp  = rinit
    allocate (   shw_obsp(mza,mnudp)) ; shw_obsp    = rinit
    allocate (uzonal_obsp(mza,mnudp)) ; uzonal_obsp = rinit
    allocate (umerid_obsp(mza,mnudp)) ; umerid_obsp = rinit

    allocate (   rho_obsf(mza,mnudp)) ; rho_obsf    = rinit
    allocate ( theta_obsf(mza,mnudp)) ; theta_obsf  = rinit
    allocate (   shw_obsf(mza,mnudp)) ; shw_obsf    = rinit
    allocate (uzonal_obsf(mza,mnudp)) ; uzonal_obsf = rinit
    allocate (umerid_obsf(mza,mnudp)) ; umerid_obsf = rinit

    allocate (    rho_obs(mza,mnudp)) ; rho_obs     = rinit
    allocate (  theta_obs(mza,mnudp)) ; theta_obs   = rinit
    allocate (    shw_obs(mza,mnudp)) ; shw_obs     = rinit
    allocate ( uzonal_obs(mza,mnudp)) ; uzonal_obs  = rinit
    allocate ( umerid_obs(mza,mnudp)) ; umerid_obs  = rinit

    allocate (    rho_sim(mza,mnudp)) ; rho_sim     = rinit
    allocate (  theta_sim(mza,mnudp)) ; theta_sim   = rinit
    allocate (    shw_sim(mza,mnudp)) ; shw_sim     = rinit
    allocate ( uzonal_sim(mza,mnudp)) ; uzonal_sim  = rinit
    allocate ( umerid_sim(mza,mnudp)) ; umerid_sim  = rinit

  end subroutine alloc_nudge2

  !=========================================================================

  subroutine filltab_nudge()

    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(rho_obsp)) then
         call increment_vtable('RHO_OBSP', 'AN')
         vtab_r(num_var)%rvar2_p => rho_obsp
      endif

    if (allocated(theta_obsp)) then
         call increment_vtable('THETA_OBSP', 'AN')
         vtab_r(num_var)%rvar2_p => theta_obsp
      endif

    if (allocated(shw_obsp)) then
         call increment_vtable('SHW_OBSP', 'AN')
         vtab_r(num_var)%rvar2_p => shw_obsp
      endif

    if (allocated(uzonal_obsp)) then
         call increment_vtable('UZONAL_OBSP', 'AN')
         vtab_r(num_var)%rvar2_p => uzonal_obsp
      endif

    if (allocated(umerid_obsp)) then
         call increment_vtable('UMERID_OBSP', 'AN')
         vtab_r(num_var)%rvar2_p => umerid_obsp
      endif

    if (allocated(rho_obsf)) then
         call increment_vtable('RHO_OBSF', 'AN')
         vtab_r(num_var)%rvar2_p => rho_obsf
      endif

    if (allocated(theta_obsf)) then
         call increment_vtable('THETA_OBSF', 'AN')
         vtab_r(num_var)%rvar2_p => theta_obsf
      endif

    if (allocated(shw_obsf)) then
         call increment_vtable('SHW_OBSF', 'AN')
         vtab_r(num_var)%rvar2_p => shw_obsf
      endif

    if (allocated(uzonal_obsf)) then
         call increment_vtable('UZONAL_OBSF', 'AN')
         vtab_r(num_var)%rvar2_p => uzonal_obsf
      endif

    if (allocated(umerid_obsf)) then
         call increment_vtable('UMERID_OBSF', 'AN')
         vtab_r(num_var)%rvar2_p => umerid_obsf
      endif

  end subroutine filltab_nudge

End Module mem_nudge
