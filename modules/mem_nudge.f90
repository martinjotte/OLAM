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

  Type itab_wnudge          ! data structure for nudging points (individual rank)
    integer :: npoly = 0    ! number of W neighbors of this WNUD pt
    integer :: iwnud(6) = 1 ! array of WNUD neighbors of this WNUD pt
  End Type itab_wnudge

  type (itab_wnudge), allocatable :: itab_wnud(:)

  real,    allocatable, target ::      xewnud(:)
  real,    allocatable, target ::      yewnud(:)
  real,    allocatable, target ::      zewnud(:)

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
  integer :: nwnud, mwnud

  real    :: tnudcent

Contains

  subroutine alloc_nudge1(lwnud)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: lwnud

!   Allocate arrays based on options (if necessary)

    write(io6,*) 'allocating nudge1 ', lwnud

    allocate (itab_wnud(lwnud))

    allocate (xewnud(lwnud)) ; xewnud = rinit
    allocate (yewnud(lwnud)) ; yewnud = rinit
    allocate (zewnud(lwnud)) ; zewnud = rinit

  end subroutine alloc_nudge1

!=========================================================================

  subroutine alloc_nudge2(mza)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza

!   Allocate arrays based on options (if necessary)

    write(io6,*) 'allocating nudge2 ',mza,mwnud

    allocate (   rho_obsp(mza,mwnud)) ; rho_obsp    = rinit
    allocate ( theta_obsp(mza,mwnud)) ; theta_obsp  = rinit
    allocate (   shw_obsp(mza,mwnud)) ; shw_obsp    = rinit
    allocate (uzonal_obsp(mza,mwnud)) ; uzonal_obsp = rinit
    allocate (umerid_obsp(mza,mwnud)) ; umerid_obsp = rinit

    allocate (   rho_obsf(mza,mwnud)) ; rho_obsf    = rinit
    allocate ( theta_obsf(mza,mwnud)) ; theta_obsf  = rinit
    allocate (   shw_obsf(mza,mwnud)) ; shw_obsf    = rinit
    allocate (uzonal_obsf(mza,mwnud)) ; uzonal_obsf = rinit
    allocate (umerid_obsf(mza,mwnud)) ; umerid_obsf = rinit

    allocate (    rho_obs(mza,mwnud)) ; rho_obs     = rinit
    allocate (  theta_obs(mza,mwnud)) ; theta_obs   = rinit
    allocate (    shw_obs(mza,mwnud)) ; shw_obs     = rinit
    allocate ( uzonal_obs(mza,mwnud)) ; uzonal_obs  = rinit
    allocate ( umerid_obs(mza,mwnud)) ; umerid_obs  = rinit

    allocate (    rho_sim(mza,mwnud)) ; rho_sim     = rinit
    allocate (  theta_sim(mza,mwnud)) ; theta_sim   = rinit
    allocate (    shw_sim(mza,mwnud)) ; shw_sim     = rinit
    allocate ( uzonal_sim(mza,mwnud)) ; uzonal_sim  = rinit
    allocate ( umerid_sim(mza,mwnud)) ; umerid_sim  = rinit

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
