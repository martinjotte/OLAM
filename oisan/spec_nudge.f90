subroutine init_spec_nudge()

  use mem_ijtabs,   only: jtab_w, itab_w, jtw_init
  use mem_grid,     only: mza, volt
  use olam_mpi_nud, only: mpi_send_wnud, mpi_recv_wnud
  use misc_coms,    only: iparallel
  use consts_coms,  only: r8
  use misc_coms,    only: rinit
  use mem_nudge,    only: mwnud, nudflag, nudnxp, volwnudi, &
                          rho_sim, theta_sim, rrv_sim, uzonal_sim, umerid_sim

  implicit none

  integer  :: j, iw, iwnud, k
  real(r8) :: volni(mza,mwnud)

  if (nudflag < 1) return
  if (nudnxp  < 1) return

  allocate (    rho_sim(mza,mwnud)) ; rho_sim    = rinit
  allocate (  theta_sim(mza,mwnud)) ; theta_sim  = rinit
  allocate (    rrv_sim(mza,mwnud)) ; rrv_sim    = rinit
  allocate ( uzonal_sim(mza,mwnud)) ; uzonal_sim = rinit
  allocate ( umerid_sim(mza,mwnud)) ; umerid_sim = rinit
  allocate (   volwnudi(mza,mwnud)) ; volwnudi   = 0.0_r8

  volni = 0._r8

  do j = 1, jtab_w(jtw_init)%jend
     iw    = jtab_w(jtw_init)%iw(j)
     iwnud = itab_w(iw)%iwnud(1)

     do k = 2, mza
        volni(k,iwnud) = volni(k,iwnud) + volt(k,iw)
     enddo
  enddo

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

subroutine nudge_prep_spec ( iaction )

  use mem_nudge,   only: mwnud, volwnudi, &
                         rho_obsp, theta_obsp, rrv_obsp, &
                         rho_obsf, theta_obsf, rrv_obsf, &
                         uzonal_obsp, umerid_obsp, &
                         uzonal_obsf, umerid_obsf

  use isan_coms,   only: o_rho, o_theta, o_rrw, o_uzonal, o_umerid
  use mem_grid,    only: mza, mwa, volt
  use misc_coms,   only: iparallel
  use mem_ijtabs,  only: jtab_w, itab_w, jtw_init
  use olam_mpi_nud,only: mpi_send_wnud, mpi_recv_wnud
  use consts_coms, only: r8

  implicit none

  integer,   intent(in) :: iaction
  integer               :: j,iw,k,iwnud,iwnud1

  real(r8), allocatable :: drho   (:,:)
  real(r8), allocatable :: dtheta (:,:)
  real(r8), allocatable :: drrv   (:,:)
  real(r8), allocatable :: duzonal(:,:)
  real(r8), allocatable :: dumerid(:,:)

  ! Allocate and zero out nudging polygon arrays

  allocate(drho   (mza,mwnud)) ; drho    = 0._r8
  allocate(dtheta (mza,mwnud)) ; dtheta  = 0._r8
  allocate(drrv   (mza,mwnud)) ; drrv    = 0._r8
  allocate(duzonal(mza,mwnud)) ; duzonal = 0._r8
  allocate(dumerid(mza,mwnud)) ; dumerid = 0._r8

  ! Swap future data time into past data time if necessary.

  if (iaction == 1) then

     rho_obsp    = rho_obsf
     theta_obsp  = theta_obsf
     rrv_obsp    = rrv_obsf
     uzonal_obsp = uzonal_obsf
     umerid_obsp = umerid_obsf

  endif

  ! If doing spectral nudging, sum data to nudging polygon arrays

  !----------------------------------------------------------------------
  do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
     iwnud1 = itab_w(iw)%iwnud(1)
  !---------------------------------------------------------------------

     do k = 2, mza
        drho   (k,iwnud1) =    drho(k,iwnud1) + o_rho   (k,iw) * volt(k,iw)
        dtheta (k,iwnud1) =  dtheta(k,iwnud1) + o_theta (k,iw) * volt(k,iw)
        drrv   (k,iwnud1) =    drrv(k,iwnud1) + o_rrw   (k,iw) * volt(k,iw)
        duzonal(k,iwnud1) = duzonal(k,iwnud1) + o_uzonal(k,iw) * volt(k,iw)
        dumerid(k,iwnud1) = dumerid(k,iwnud1) + o_umerid(k,iw) * volt(k,iw)
     enddo

  enddo

  ! MPI SEND/RECV of nudging arrays

  if (iparallel == 1) then

     call mpi_send_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrv, dvara4=duzonal, dvara5=dumerid)

     call mpi_recv_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrv, dvara4=duzonal, dvara5=dumerid)
  endif

  ! Normalize nudging point sums to get average values
  ! Horizontal loop over nudging polygons

  !$omp parallel do private(k)
  do iwnud = 2, mwnud
     do k = 2, mza
        rho_obsf   (k,iwnud) =    drho(k,iwnud) * volwnudi(k,iwnud)
        theta_obsf (k,iwnud) =  dtheta(k,iwnud) * volwnudi(k,iwnud)
        rrv_obsf   (k,iwnud) =    drrv(k,iwnud) * volwnudi(k,iwnud)
        uzonal_obsf(k,iwnud) = duzonal(k,iwnud) * volwnudi(k,iwnud)
        umerid_obsf(k,iwnud) = dumerid(k,iwnud) * volwnudi(k,iwnud)
     enddo
  enddo
  !$omp end parallel do

  deallocate(drho)
  deallocate(dtheta)
  deallocate(drrv)
  deallocate(duzonal)
  deallocate(dumerid)

end subroutine nudge_prep_spec

!==============================================================================

subroutine spec_nudge()

  use mem_nudge, only:   tnudcent,      mwnud,    volwnudi,   rhot_nud,  &
                          rho_sim,    rho_obs,    rho_obsp,    rho_obsf, &
                        theta_sim,  theta_obs,  theta_obsp,  theta_obsf, &
                          rrv_sim,    rrv_obs,    rrv_obsp,    rrv_obsf, &
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
  use olam_mpi_nud,only: mpi_send_wnud, mpi_recv_wnud

  implicit none

  ! Nudge selected model fields (rho, thil, rr_w, vmc) to observed data
  ! using polygon filtering

  integer :: iwnud,k,j,iw,iwnud1,iwnud2,iwnud3

  real :: umzonalt, ummeridt
  real :: uzonal, umerid
  real :: tp, tf, tnudi, rho4
  real :: fnud1, fnud2, fnud3

  real(r8), allocatable :: drho   (:,:)
  real(r8), allocatable :: dtheta (:,:)
  real(r8), allocatable :: drrv   (:,:)
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

  ! Time interpolation coefficients

  tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
            / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

  tp = 1. - tf

  ! If doing spectral nudging, zero out nudging polygon arrays and volume counter
  ! prior to summing

  ! Allocate and zero out nudging polygon arrays

  allocate(drho   (mza,mwnud)) ; drho    = 0._r8
  allocate(dtheta (mza,mwnud)) ; dtheta  = 0._r8
  allocate(drrv   (mza,mwnud)) ; drrv    = 0._r8
  allocate(duzonal(mza,mwnud)) ; duzonal = 0._r8
  allocate(dumerid(mza,mwnud)) ; dumerid = 0._r8

  ! Horizontal loop over W columns
  !----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
     iwnud1 = itab_w(iw)%iwnud(1)
  !---------------------------------------------------------------------

     ! Sum model fields to nudging polygon arrays
     do k = 2, mza

        uzonal = vxe(k,iw) * vxn_ew(iw) + vye(k,iw) * vyn_ew(iw)
        umerid = vxe(k,iw) * vxn_ns(iw) + vye(k,iw) * vyn_ns(iw) &
               + vze(k,iw) * vzn_ns(iw)

        drho   (k,iwnud1) =    drho(k,iwnud1) + rho  (k,iw) * volt(k,iw)
        dtheta (k,iwnud1) =  dtheta(k,iwnud1) + theta(k,iw) * volt(k,iw)
        drrv   (k,iwnud1) =    drrv(k,iwnud1) + rr_w (k,iw) * volt(k,iw)
        duzonal(k,iwnud1) = duzonal(k,iwnud1) + uzonal      * volt(k,iw)
        dumerid(k,iwnud1) = dumerid(k,iwnud1) + umerid      * volt(k,iw)
     enddo

  enddo

! MPI SEND/RECV of nudging arrays

  if (iparallel == 1) then

     call mpi_send_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrv, dvara4=duzonal, dvara5=dumerid)

     call mpi_recv_wnud(dvara1=drho, dvara2=dtheta,  &
                        dvara3=drrv, dvara4=duzonal, dvara5=dumerid)
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
           rrv_sim(k,iwnud) =    drrv(k,iwnud) * volwnudi(k,iwnud)
        uzonal_sim(k,iwnud) = duzonal(k,iwnud) * volwnudi(k,iwnud)
        umerid_sim(k,iwnud) = dumerid(k,iwnud) * volwnudi(k,iwnud)

           rho_obs(k,iwnud) = tp *    rho_obsp(k,iwnud) + tf *    rho_obsf(k,iwnud)
         theta_obs(k,iwnud) = tp *  theta_obsp(k,iwnud) + tf *  theta_obsf(k,iwnud)
           rrv_obs(k,iwnud) = tp *    rrv_obsp(k,iwnud) + tf *    rrv_obsf(k,iwnud)
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
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
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
                       ( fnud1 * (rrv_obs(k,iwnud1) - rrv_sim(k,iwnud1)) &
                       + fnud2 * (rrv_obs(k,iwnud2) - rrv_sim(k,iwnud2)) &
                       + fnud3 * (rrv_obs(k,iwnud3) - rrv_sim(k,iwnud3)) )

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
  !$omp end parallel

  deallocate(drho)
  deallocate(dtheta)
  deallocate(drrv)
  deallocate(duzonal)
  deallocate(dumerid)

end subroutine spec_nudge
