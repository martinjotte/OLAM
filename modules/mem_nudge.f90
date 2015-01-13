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

Module mem_nudge

  implicit none
   
  integer, parameter :: nloops_wnud = 100 ! # WNUD DO loops for para

  Type itab_wnud_vars          ! data structure for WNUD points (individual rank)
    logical :: loop(nloops_wnud) = .false. ! flag to perform DO loop at this WNUD pt
    integer :: npoly = 0      ! number of W neighbors of this WNUD pt
    integer :: irank = -1     ! rank of process at this WNUD pt (for hist write only)
    integer :: iwnud(6) = 1   ! array of WNUD neighbors of this WNUD pt
    integer :: iwnudglobe = 1 ! global index of WNUD point
  End Type itab_wnud_vars

  type (itab_wnud_vars), allocatable, target :: itab_wnud(:)

  Type itabg_wnud_vars            ! data structure for WNUD pts (global)
     integer :: iwnud_myrank = -1 ! local (parallel subdomain) index of this WNUD pt
     integer :: irank = -1        ! rank of process at this WNUD pt
  End Type itabg_wnud_vars

  type (itabg_wnud_vars), allocatable, target :: itabg_wnud(:)

  Type jtab_wnud_vars
     integer, allocatable :: iwnud(:)
     integer              :: jend
  End Type jtab_wnud_vars

  type (jtab_wnud_vars) :: jtab_wnud(nloops_wnud)

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
  integer :: nwnud = 1
  integer :: mwnud = 1

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

    if (nudnxp /= 0) then
       allocate (    rho_sim(mza,mwnud)) ; rho_sim     = rinit
       allocate (  theta_sim(mza,mwnud)) ; theta_sim   = rinit
       allocate (    shw_sim(mza,mwnud)) ; shw_sim     = rinit
       allocate ( uzonal_sim(mza,mwnud)) ; uzonal_sim  = rinit
       allocate ( umerid_sim(mza,mwnud)) ; umerid_sim  = rinit
    endif

  end subroutine alloc_nudge2

  !=========================================================================

  subroutine filltab_nudge()

    use var_tables, only: increment_vtable
    implicit none

    character(2) :: stagpt
    
    ! obs (point-by-point) nudging is done at W points, while spectral
    ! nudging is done on a separate nudging mesh. Currently, parallel
    ! output on the nudging mesh in not implemented

    if (nudnxp == 0) then
       stagpt = 'AW'
    else
       stagpt = 'AN'
    endif

    if (allocated(rho_obsp))    call increment_vtable('RHO_OBSP',    stagpt, noread=.true., rvar2=rho_obsp)

    if (allocated(theta_obsp))  call increment_vtable('THETA_OBSP',  stagpt, noread=.true., rvar2=theta_obsp)

    if (allocated(shw_obsp))    call increment_vtable('SHW_OBSP',    stagpt, noread=.true., rvar2=shw_obsp)

    if (allocated(uzonal_obsp)) call increment_vtable('UZONAL_OBSP', stagpt, noread=.true., rvar2=uzonal_obsp)

    if (allocated(umerid_obsp)) call increment_vtable('UMERID_OBSP', stagpt, noread=.true., rvar2=umerid_obsp)

    if (allocated(rho_obsf))    call increment_vtable('RHO_OBSF',    stagpt, noread=.true., rvar2=rho_obsf)

    if (allocated(theta_obsf))  call increment_vtable('THETA_OBSF',  stagpt, noread=.true., rvar2=theta_obsf)

    if (allocated(shw_obsf))    call increment_vtable('SHW_OBSF',    stagpt, noread=.true., rvar2=shw_obsf)

    if (allocated(uzonal_obsf)) call increment_vtable('UZONAL_OBSF', stagpt, noread=.true., rvar2=uzonal_obsf)

    if (allocated(umerid_obsf)) call increment_vtable('UMERID_OBSF', stagpt, noread=.true., rvar2=umerid_obsf)

  end subroutine filltab_nudge

!===============================================================================

   subroutine fill_jnudge()

   use misc_coms,  only: io6

   implicit none

   integer :: iwnud, iloop, jend

! Allocate and zero-fill jtab%jend()

   do iloop = 1,nloops_wnud
      jtab_wnud(iloop)%jend = 0
   enddo

! Compute and store jtab%jend(1)

   do iloop = 1,nloops_wnud
      jtab_wnud(iloop)%jend = 0
      do iwnud = 2,mwnud
         if (itab_wnud(iwnud)%loop(iloop)) then
            jtab_wnud(iloop)%jend = jtab_wnud(iloop)%jend + 1
         endif
      enddo
      jtab_wnud(iloop)%jend = max(1,jtab_wnud(iloop)%jend)
   enddo

! Allocate and zero-fill JTAB_WNUD%IWNUD

   do iloop = 1,nloops_wnud
      jend = jtab_wnud(iloop)%jend
      allocate (jtab_wnud(iloop)%iwnud(jend))
      jtab_wnud(iloop)%iwnud(1:jend) = 0
   enddo

! Initialize JTAB%JEND counters to zero

   do iloop = 1,nloops_wnud
      jtab_wnud(iloop)%jend = 0
   enddo

! Compute JTAB_WNUD%IWNUD

   do iwnud = 2,mwnud
      do iloop = 1,nloops_wnud
         if (itab_wnud(iwnud)%loop(iloop)) then
            jtab_wnud(iloop)%jend = jtab_wnud(iloop)%jend + 1
            jtab_wnud(iloop)%iwnud(jtab_wnud(iloop)%jend) = iwnud
         endif
      enddo
   enddo

   end subroutine fill_jnudge

End Module mem_nudge
