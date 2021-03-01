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

!=========================================================================

subroutine init_spec_nudge()

  use mem_ijtabs,   only: jtab_w, itab_w, jtw_init
  use mem_grid,     only: mza, volt
  use olam_mpi_atm, only: mpi_send_wnud, mpi_recv_wnud
  use misc_coms,    only: iparallel
  use consts_coms,  only: r8
  use misc_coms,    only: rinit
  use mem_nudge,    only: mwnud, nudflag, nudnxp, volwnudi, &
                          rho_sim, theta_sim, rrw_sim, uzonal_sim, umerid_sim

  implicit none

  integer  :: j, iw, iwnud, k
  real(r8) :: volni(mza,mwnud)

  if (nudflag < 1) return
  if (nudnxp  < 1) return

  !$omp parallel

  !$omp sections
  allocate (    rho_sim(mza,mwnud)) ; rho_sim    = rinit
  !$omp section
  allocate (  theta_sim(mza,mwnud)) ; theta_sim  = rinit
  !$omp section
  allocate (    rrw_sim(mza,mwnud)) ; rrw_sim    = rinit
  !$omp section
  allocate ( uzonal_sim(mza,mwnud)) ; uzonal_sim = rinit
  !$omp section
  allocate ( umerid_sim(mza,mwnud)) ; umerid_sim = rinit
  !$omp section
  allocate (   volwnudi(mza,mwnud)) ; volwnudi   = 0.0_r8
  !$omp section
  volni = 0._r8
  !$omp end sections nowait

  !$omp do private(iw,iwnud,k) reduction(+:volni)
  do j = 1, jtab_w(jtw_init)%jend(1)

     iw    = jtab_w(jtw_init)%iw(j)
     iwnud = itab_w(iw)%iwnud(1)

     do k = 2, mza
        volni(k,iwnud) = volni(k,iwnud) + volt(k,iw)
     enddo

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! MPI SEND/RECV of nudging arrays

  if (iparallel == 1) call mpi_send_wnud(dvara1=volni)
  if (iparallel == 1) call mpi_recv_wnud(dvara1=volni)

  !$omp parallel do private(k)
  do iwnud = 2, mwnud
     do k = 2, mza
        volwnudi(k,iwnud) = 1._r8 / max(1._r8, volni(k,iwnud))
     enddo
  enddo
  !$omp end parallel do

end subroutine init_spec_nudge

!===============================================================================

subroutine nudge_prep_spec(iaction, o_rho, o_theta, o_rrw, o_uzonal, o_umerid)

  use mem_nudge,   only: mwnud, volwnudi, &
                         rho_obsp, theta_obsp, rrw_obsp, &
                         rho_obsf, theta_obsf, rrw_obsf, &
                         uzonal_obsp, umerid_obsp, &
                         uzonal_obsf, umerid_obsf

  use mem_grid,    only: mza, mwa, volt
  use misc_coms,   only: iparallel
  use mem_ijtabs,  only: jtab_w, itab_w, jtw_init
  use olam_mpi_atm,only: mpi_send_wnud, mpi_recv_wnud
  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: iaction

  real(r8), intent(in)  :: o_rho   (mza,mwa)
  real,     intent(in)  :: o_theta (mza,mwa)
  real,     intent(in)  :: o_rrw   (mza,mwa)
  real,     intent(in)  :: o_uzonal(mza,mwa)
  real,     intent(in)  :: o_umerid(mza,mwa)

  integer               :: j,iw,k,iwnud,iwnud1

  real(r8), allocatable :: drho   (:,:)
  real(r8), allocatable :: dtheta (:,:)
  real(r8), allocatable :: drrw   (:,:)
  real(r8), allocatable :: duzonal(:,:)
  real(r8), allocatable :: dumerid(:,:)

  !$omp parallel

  ! Allocate and zero out nudging polygon arrays

  !$omp sections
  allocate(drho   (mza,mwnud)) ; drho    = 0._r8
  !$omp section
  allocate(dtheta (mza,mwnud)) ; dtheta  = 0._r8
  !$omp section
  allocate(drrw   (mza,mwnud)) ; drrw    = 0._r8
  !$omp section
  allocate(duzonal(mza,mwnud)) ; duzonal = 0._r8
  !$omp section
  allocate(dumerid(mza,mwnud)) ; dumerid = 0._r8
  !$omp end sections nowait

  ! Swap future data time into past data time if necessary.

  if (iaction == 1) then

     !$omp sections
     rho_obsp = rho_obsf
     !$omp section
     theta_obsp = theta_obsf
     !$omp section
     rrw_obsp = rrw_obsf
     !$omp section
     uzonal_obsp = uzonal_obsf
     !$omp section
     umerid_obsp = umerid_obsf
     !$omp end sections

  endif

  ! If doing spectral nudging, sum data to nudging polygon arrays

  !----------------------------------------------------------------------
  !$omp do private(iw,iwnud1,k) reduction(+:drho)  &
  !$omp    reduction(+:dtheta)  reduction(+:drrw)  &
  !$omp    reduction(+:duzonal) reduction(+:dumerid)
  do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
     iwnud1 = itab_w(iw)%iwnud(1)
  !---------------------------------------------------------------------

     do k = 2, mza
        drho   (k,iwnud1) =    drho(k,iwnud1) + o_rho   (k,iw) * volt(k,iw)
        dtheta (k,iwnud1) =  dtheta(k,iwnud1) + o_theta (k,iw) * volt(k,iw)
        drrw   (k,iwnud1) =    drrw(k,iwnud1) + o_rrw   (k,iw) * volt(k,iw)
        duzonal(k,iwnud1) = duzonal(k,iwnud1) + o_uzonal(k,iw) * volt(k,iw)
        dumerid(k,iwnud1) = dumerid(k,iwnud1) + o_umerid(k,iw) * volt(k,iw)
     enddo

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! MPI SEND/RECV of nudging arrays

  if (iparallel == 1) then

     call mpi_send_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrw, dvara4=duzonal, dvara5=dumerid)

     call mpi_recv_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrw, dvara4=duzonal, dvara5=dumerid)
  endif

  ! Normalize nudging point sums to get average values
  ! Horizontal loop over nudging polygons

  !$omp parallel
  !$omp do private(k)
  do iwnud = 2, mwnud
     do k = 2, mza
        rho_obsf   (k,iwnud) =    drho(k,iwnud) * volwnudi(k,iwnud)
        theta_obsf (k,iwnud) =  dtheta(k,iwnud) * volwnudi(k,iwnud)
        rrw_obsf   (k,iwnud) =    drrw(k,iwnud) * volwnudi(k,iwnud)
        uzonal_obsf(k,iwnud) = duzonal(k,iwnud) * volwnudi(k,iwnud)
        umerid_obsf(k,iwnud) = dumerid(k,iwnud) * volwnudi(k,iwnud)
     enddo
  enddo
  !$omp end do

  !$omp sections
  deallocate(drho)
  !$omp section
  deallocate(dtheta)
  !$omp section
  deallocate(drrw)
  !$omp section
  deallocate(duzonal)
  !$omp section
  deallocate(dumerid)
  !$omp end sections nowait
  !$omp end parallel

end subroutine nudge_prep_spec

!==============================================================================

subroutine spec_nudge(mrl)

  use mem_nudge, only:   tnudcent,      mwnud,    volwnudi,   rhot_nud,  &
                          rho_sim,    rho_obs,    rho_obsp,    rho_obsf, &
                        theta_sim,  theta_obs,  theta_obsp,  theta_obsf, &
                          rrw_sim,    rrw_obs,    rrw_obsp,    rrw_obsf, &
                       uzonal_sim, uzonal_obs, uzonal_obsp, uzonal_obsf, &
                       umerid_sim, umerid_obs, umerid_obsp, umerid_obsf

  use mem_basic,   only: rho, theta, rr_w, vxe, vye, vze
  use mem_grid,    only: vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns, &
                         mza, lpw, volt
  use misc_coms,   only: s1900_sim, iparallel
  use mem_ijtabs,  only: jtab_w, itab_w, jtv_prog, jtw_prog
  use consts_coms, only: r8
  use mem_tend,    only: thilt, rr_wt, vmxet, vmyet, vmzet
  use isan_coms,   only: ifgfile, s1900_fg
  use olam_mpi_atm,only: mpi_send_wnud, mpi_recv_wnud

  implicit none

  ! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data
  ! using polygon filtering

  integer, intent(in) :: mrl

  integer :: iwnud,k,j,iw,iwnud1,iwnud2,iwnud3,kb

  real :: umzonalt, ummeridt
  real :: uzonal, umerid
  real :: tp, tf, tnudi, rho4
  real :: fnud1, fnud2, fnud3

  real(r8), allocatable :: drho   (:,:)
  real(r8), allocatable :: dtheta (:,:)
  real(r8), allocatable :: drrw   (:,:)
  real(r8), allocatable :: duzonal(:,:)
  real(r8), allocatable :: dumerid(:,:)

!----------------------------------------------------------------------
! EXAMPLE - DEFINE OPTIONAL SPATIAL NUDGING MASK
!
! integer, save :: icall = 0
! real, save, allocatable :: wtnud(:)
! real :: xw,yw,dist
!
! if (icall /= 1) then
!    icall = 1
!
!    allocate (wtnud(mwa))
!
!! Horizontal loop over T points
!
!    do iw = 2,mwa
!
!! Default: Uniform nudging weight = 1
!
!       wtnud(iw) = 1.
!
!! Sample code for modifying nudging weight
!
!! Transform current IW point to polar stereographic coordinates using specified
!! pole point location (pole point lat/lon = 4th & 5th arguments of e_ps)
!
!       call e_ps(xew(iw),yew(iw),zew(iw),37.,-117.,xw,yw)
!
!       dist = sqrt(xw ** 2 + yw ** 2)
!
!       if (dist > 5000.e3) then
!          wtnud(iw) = 1.
!       elseif (dist < 3600.e3) then
!          wtnud(iw) = 0.
!       else
!          wtnud(iw) = ((dist - 3600.e3) / 1400.e3) ** 2
!       endif
!
!    enddo
!
! endif
!----------------------------------------------------------------------

  ! Check whether it is time to nudge

  if (mrl < 1) return

  ! Time interpolation coefficients

  tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
            / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

  tp = 1. - tf

  !$omp parallel

  ! If doing spectral nudging, zero out nudging polygon arrays and volume counter
  ! prior to summing

  ! Allocate and zero out nudging polygon arrays

  !$omp sections
  allocate(drho   (mza,mwnud)) ; drho    = 0._r8
  !$omp section
  allocate(dtheta (mza,mwnud)) ; dtheta  = 0._r8
  !$omp section
  allocate(drrw   (mza,mwnud)) ; drrw    = 0._r8
  !$omp section
  allocate(duzonal(mza,mwnud)) ; duzonal = 0._r8
  !$omp section
  allocate(dumerid(mza,mwnud)) ; dumerid = 0._r8
  !$omp end sections

  ! Horizontal loop over W columns
  !----------------------------------------------------------------------
  !$omp do private(iw,iwnud1,kb,k,uzonal,umerid)   &
  !$omp    reduction(+:drho)                       &
  !$omp    reduction(+:dtheta)  reduction(+:drrw)  &
  !$omp    reduction(+:duzonal) reduction(+:dumerid)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
     iwnud1 = itab_w(iw)%iwnud(1)
  !---------------------------------------------------------------------

     ! Sum model fields to nudging polygon arrays
     do k = kb, mza

        uzonal = vxe(k,iw) * vxn_ew(iw) + vye(k,iw) * vyn_ew(iw)
        umerid = vxe(k,iw) * vxn_ns(iw) + vye(k,iw) * vyn_ns(iw) &
               + vze(k,iw) * vzn_ns(iw)

        drho   (k,iwnud1) =    drho(k,iwnud1) + rho  (k,iw) * volt(k,iw)
        dtheta (k,iwnud1) =  dtheta(k,iwnud1) + theta(k,iw) * volt(k,iw)
        drrw   (k,iwnud1) =    drrw(k,iwnud1) + rr_w (k,iw) * volt(k,iw)
        duzonal(k,iwnud1) = duzonal(k,iwnud1) + uzonal      * volt(k,iw)
        dumerid(k,iwnud1) = dumerid(k,iwnud1) + umerid      * volt(k,iw)
     enddo

  enddo
  !$omp end do
  !$omp end parallel

! MPI SEND/RECV of nudging arrays

  if (iparallel == 1) then

     call mpi_send_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrw, dvara4=duzonal, dvara5=dumerid)

     call mpi_recv_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrw, dvara4=duzonal, dvara5=dumerid)
  endif

  ! Horizontal loop over nudging polygons

  !$omp parallel
  !$omp do private(k)
  do iwnud = 2,mwnud

     ! If doing spectral nudging, normalize nudging point sums to get average values
     ! and interpolate observational fields in time

     ! Vertical loop over nudging polygons
     do k = 2,mza

           rho_sim(k,iwnud) =    drho(k,iwnud) * volwnudi(k,iwnud)
         theta_sim(k,iwnud) =  dtheta(k,iwnud) * volwnudi(k,iwnud)
           rrw_sim(k,iwnud) =    drrw(k,iwnud) * volwnudi(k,iwnud)
        uzonal_sim(k,iwnud) = duzonal(k,iwnud) * volwnudi(k,iwnud)
        umerid_sim(k,iwnud) = dumerid(k,iwnud) * volwnudi(k,iwnud)

           rho_obs(k,iwnud) = tp *    rho_obsp(k,iwnud) + tf *    rho_obsf(k,iwnud)
         theta_obs(k,iwnud) = tp *  theta_obsp(k,iwnud) + tf *  theta_obsf(k,iwnud)
           rrw_obs(k,iwnud) = tp *    rrw_obsp(k,iwnud) + tf *    rrw_obsf(k,iwnud)
        uzonal_obs(k,iwnud) = tp * uzonal_obsp(k,iwnud) + tf * uzonal_obsf(k,iwnud)
        umerid_obs(k,iwnud) = tp * umerid_obsp(k,iwnud) + tf * umerid_obsf(k,iwnud)

     enddo

  enddo
  !$omp end do

  ! Loop over all W columns, find 3 neighboring polygon points for each,
  ! and interpolate (obs - model) differences at each polygon point to the W point

  !----------------------------------------------------------------------
  !$omp do private (iw,iwnud1,iwnud2,iwnud3,fnud1,fnud2,fnud3,tnudi, &
  !$omp             rho4,umzonalt,ummeridt)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
     iwnud1 = itab_w(iw)%iwnud(1);  fnud1 = itab_w(iw)%fnud(1)
     iwnud2 = itab_w(iw)%iwnud(2);  fnud2 = itab_w(iw)%fnud(2)
     iwnud3 = itab_w(iw)%iwnud(3);  fnud3 = itab_w(iw)%fnud(3)
  !---------------------------------------------------------------------

     ! Inverse of nudging time scale
     ! Use spatial nudging mask defined above in this subroutine

   ! tnudi = wtnud(iw) / tnudcent
     tnudi = 1.0 / tnudcent

     do k = lpw(iw), mza

        rho4 = real(rho(k,iw))

        rhot_nud(k,iw) = tnudi * &
                       ( fnud1 * (rho_obs(k,iwnud1) - rho_sim(k,iwnud1)) &
                       + fnud2 * (rho_obs(k,iwnud2) - rho_sim(k,iwnud2)) &
                       + fnud3 * (rho_obs(k,iwnud3) - rho_sim(k,iwnud3)) )

        thilt(k,iw)    = thilt(k,iw) + tnudi * &
                       ( fnud1 * ( rho_obs(k,iwnud1) * theta_obs(k,iwnud1)   &
                                 - rho_sim(k,iwnud1) * theta_sim(k,iwnud1) ) &
                       + fnud2 * ( rho_obs(k,iwnud2) * theta_obs(k,iwnud2)   &
                                 - rho_sim(k,iwnud2) * theta_sim(k,iwnud2) ) &
                       + fnud3 * ( rho_obs(k,iwnud3) * theta_obs(k,iwnud3)   &
                                 - rho_sim(k,iwnud3) * theta_sim(k,iwnud3) ) )

        rr_wt(k,iw)    = tnudi * rho4 * &
                       ( fnud1 * (rrw_obs(k,iwnud1) - rrw_sim(k,iwnud1)) &
                       + fnud2 * (rrw_obs(k,iwnud2) - rrw_sim(k,iwnud2)) &
                       + fnud3 * (rrw_obs(k,iwnud3) - rrw_sim(k,iwnud3)) )

        umzonalt       = tnudi * rho4 * &
                       ( fnud1 * (uzonal_obs(k,iwnud1) - uzonal_sim(k,iwnud1)) &
                       + fnud2 * (uzonal_obs(k,iwnud2) - uzonal_sim(k,iwnud2)) &
                       + fnud3 * (uzonal_obs(k,iwnud3) - uzonal_sim(k,iwnud3)) )

        ummeridt       = tnudi  * rho4 * &
                       ( fnud1 * (umerid_obs(k,iwnud1) - umerid_sim(k,iwnud1)) &
                       + fnud2 * (umerid_obs(k,iwnud2) - umerid_sim(k,iwnud2)) &
                       + fnud3 * (umerid_obs(k,iwnud3) - umerid_sim(k,iwnud3)) )

        vmxet(k,iw) = vmxet(k,iw) + vxn_ew(iw) * umzonalt + vxn_ns(iw) * ummeridt &
                    + vxe(k,iw) * rhot_nud(k,iw)

        vmyet(k,iw) = vmyet(k,iw) + vyn_ew(iw) * umzonalt + vyn_ns(iw) * ummeridt &
                    + vye(k,iw) * rhot_nud(k,iw)

        vmzet(k,iw) = vmzet(k,iw) + vzn_ns(iw) * ummeridt &
                    + vze(k,iw) * rhot_nud(k,iw)

     enddo

  enddo
  !$omp end do

  !$omp sections
  deallocate(drho)
  !$omp section
  deallocate(dtheta)
  !$omp section
  deallocate(drrw)
  !$omp section
  deallocate(duzonal)
  !$omp section
  deallocate(dumerid)
  !$omp end sections
  !$omp end parallel

end subroutine spec_nudge
