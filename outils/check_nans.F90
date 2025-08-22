module check_nan

  interface mass_sum_from_rho
     module procedure         &
          mass_sum_from_rho4, &
          mass_sum_from_rho8
  end interface mass_sum_from_rho

  interface mass_sum_from_mixrat
     module procedure             &
          mass_sum_from_mixrat44, &
          mass_sum_from_mixrat84
  end interface mass_sum_from_mixrat

Contains

!===============================================================================

subroutine check_nans(icall,rvara1,rvara2)

  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

  use mem_basic,  only: rr_w,rho,thil,rr_v,wc,wmc,press,vmc,vc,vxe,vye,vze
  use mem_grid,   only: mza,mwa,lpw,mva,lpv,glatw,glonw,glatv,glonv
  use mem_tend,   only: thilt,rr_wt,vmxet,vmyet,vmzet
  use misc_coms,  only: io6
  use mem_para,   only: myrank

  implicit none

  integer :: i,k,istop
  integer, intent(in) :: icall

  real, optional, intent(in) :: rvara1 (mza,mwa)
  real, optional, intent(in) :: rvara2 (mza,mwa)

  istop = 0

  do i = 2,mwa
     do k = lpw(i),mza
        if ( ieee_is_nan(rr_w (k,i)) .or.  &
             ieee_is_nan(rho  (k,i)) .or.  &
             ieee_is_nan(thil (k,i)) .or.  &
             ieee_is_nan(press(k,i)) .or.  &
             ieee_is_nan(wc   (k,i)) .or.  &
             ieee_is_nan(wmc  (k,i)) .or.  &
             ieee_is_nan(thilt(k,i)) .or.  &
             ieee_is_nan(rr_wt(k,i)) .or.  &
             ieee_is_nan(vmxet(k,i)) .or.  &
             ieee_is_nan(vmyet(k,i)) .or.  &
             ieee_is_nan(vmzet(k,i)) .or.  &
             ieee_is_nan(vxe  (k,i)) .or.  &
             ieee_is_nan(vye  (k,i)) .or.  &
             ieee_is_nan(vze  (k,i)) .or.  &
             ieee_is_nan(rr_v (k,i)) ) then

           write(*,*) 'Node ', myrank
           write(*,*) 'NaN at: k, i, icall ', k, i, icall
           write(io6,*) ''
           write(io6,*) 'check_nans: k, i, icall, lat, lon ',k,i,icall,glatw(i),glonw(i)
           write(io6,*) 'rr_w,rho ',rr_w(k,i),rho(k,i)
           write(io6,*) 'thil,press ',thil(k,i),press(k,i)
           write(io6,*) 'wc, wmc, rr_v ',wc(k,i),wmc(k,i),rr_v(k,i)
           write(io6,*) 'vxe, vye, vze ',vxe(k,i),vye(k,i),vze(k,i)
           write(io6,*) 'thilt, rr_wt ',thilt(k,i),rr_wt(k,i)
           write(io6,*) 'vmxet, vmyet, vmzet ',vmxet(k,i),vmyet(k,i),vmzet(k,i)
           istop = 1
        endif

     enddo
  enddo

  if (present(rvara1)) then
     do i = 2,mwa
        do k = lpw(i),mza
           if ( ieee_is_nan(rvara1 (k,i))) then
              write(*,*) 'Node ', myrank
              write(*,*) 'NaN at: k, i, icall ', k, i, icall
              write(io6,*) ''
              write(io6,*) 'check_nans: k, i, icall, lat, lon ',k,i,icall,glatw(i),glonw(i)
              write(io6,*) 'rvara1 ',rvara1(k,i)
              istop = 1
           endif

        enddo
     enddo
  endif

  if (present(rvara2)) then
     do i = 2,mwa
        do k = lpw(i),mza
           if ( ieee_is_nan(rvara2 (k,i))) then
              write(*,*) 'Node ', myrank
              write(*,*) 'NaN at: k, i, icall ',k,i,icall
              write(io6,*) ''
              write(io6,*) 'check_nans: k, i, icall, lat, lon  ',k,i,icall,glatw(i),glonw(i)
              write(io6,*) 'rvara2 ',rvara2(k,i)
              istop = 1
           endif

        enddo
     enddo
  endif

  do i = 2,mva
     do k = lpv(i), mza
        if (ieee_is_nan(vmc(k,i)) .or.  &
             ieee_is_nan(vc (k,i))) then
           write(io6,*) ''
           write(*,*) 'Node ', myrank
           write(*,*) 'NaN at ', k, i, icall
           write(io6,*) 'check_nans: k, i, icall, lat, lon ',k,i,icall,glatv(i),glonv(i)
           write(io6,*) 'vmc, vc: ', vmc(k,i), vc(k,i)
           istop = 1
        endif
     enddo
  enddo

  if (istop > 0) stop 'stop: check_nans '

end subroutine check_nans

!===============================================================================

subroutine check_pos(icall)

  use mem_ijtabs, only: jtab_w, jtw_prog, istp, mrl_endl
  use mem_grid,   only: lpw, mza
  use var_tables, only: num_scalar, scalar_tab
  use oname_coms, only: nl
  use mem_para,   only: myrank

  implicit none

  integer, intent(in) :: icall

  integer :: j, iw, k, n

  if (nl%iscal_monot < 1) return

  ! Return if this is not the end of the long timestep

  if (mrl_endl(istp) == 0) return

  do n = 1, num_scalar

     if (.not. scalar_tab(n)%pdef) cycle

     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        do k = lpw(iw), mza
           if ( scalar_tab(n)%var_p(k,iw) < 0.0 .or. &
                scalar_tab(n)%var_p(k,iw) /= scalar_tab(n)%var_p(k,iw) ) then
              write(*,*) trim(scalar_tab(n)%name), " equals ",  &
                   scalar_tab(n)%var_p(k,iw), " at ", k, iw, icall, myrank
           endif
        enddo
     enddo
  enddo

end subroutine check_pos

!===============================================================================

subroutine tracer_maxmin()

  use mem_ijtabs, only: jtab_w, jtw_prog
  use mem_grid,   only: lpw, mza
  use var_tables, only: num_scalar, scalar_tab
  use mem_para,   only: myrank
  use misc_coms,  only: iparallel

#ifdef OLAM_MPI
  use mpi_f08
#endif

  implicit none

  real    :: smax(2)
  integer :: n, j, iw

  do n = 1, num_scalar

     smax(:) = -huge(1.)

     !$omp parallel do private(iw) reduction(max:smax)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        smax(1) = max(smax(1),  maxval( scalar_tab(n)%var_p(lpw(iw):mza,iw) ) )
        smax(2) = max(smax(2), -minval( scalar_tab(n)%var_p(lpw(iw):mza,iw) ) )
     enddo
     !$omp end parallel do

#ifdef OLAM_MPI
     if (iparallel == 1) then
        if (myrank==0) call MPI_Reduce(MPI_IN_PLACE, smax, 2, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD)
        if (myrank/=0) call MPI_Reduce(smax,         0,    2, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD)
     endif
#endif

     if (myrank == 0) write(*,*) "Max, Min scp: ", smax(1), -smax(2), n, trim(scalar_tab(n)%name)

  enddo

end subroutine tracer_maxmin

!===============================================================================

subroutine compute_mass_sums()

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, volt
  use misc_coms,   only: naddsc, iparallel
  use mem_addsc,   only: addsc
  use consts_coms, only: r8
  use mem_basic,   only: rho, rr_w
  use mem_para,    only: myrank, mgroupsize
  use mem_co2,     only: i_co2, rr_co2

#ifdef OLAM_MPI
  use mpi_f08
#endif

  implicit none

  real(r8) :: tot_mass_sum
  real(r8) :: dry_mass_sum
  real(r8) :: wat_mass_sum
  real(r8) :: co2_mass_sum
  real(r8) :: scp_mass_sum

  integer  :: j, iw, k

#ifdef OLAM_MPI
  real(r8) :: tmasses(4)
  integer  :: nv
#endif

  logical,  save :: firstime = .true.
  real(r8), save :: tot_mass_sum0
  real(r8), save :: dry_mass_sum0
  real(r8), save :: wat_mass_sum0
  real(r8), save :: co2_mass_sum0
  real(r8), save :: scp_mass_sum0

  dry_mass_sum = 0.0_r8
  wat_mass_sum = 0.0_r8
  co2_mass_sum = 0.0_r8
  scp_mass_sum = 0.0_r8

  !$omp parallel do private(iw,k) reduction(+:dry_mass_sum,wat_mass_sum,co2_mass_sum,scp_mass_sum)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw),mza

        dry_mass_sum = dry_mass_sum + rho(k,iw) * volt(k,iw)
        wat_mass_sum = wat_mass_sum + rho(k,iw) * volt(k,iw) * rr_w(k,iw)

        if (i_co2 > 0) then
           co2_mass_sum = co2_mass_sum + rho(k,iw) * volt(k,iw) * rr_co2(k,iw)
        endif

        if (naddsc > 0) then
           scp_mass_sum = scp_mass_sum + rho(k,iw) * volt(k,iw) * addsc(1)%sclp(k,iw)
        endif

     enddo
  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     tmasses(1) = dry_mass_sum
     tmasses(2) = wat_mass_sum
     nv = 2

     if (i_co2 > 0) then
        nv = nv + 1
        tmasses(nv) = co2_mass_sum
     endif

     if (naddsc > 0) then
        nv = nv + 1
        tmasses(nv) = scp_mass_sum
     endif

     if (myrank==0) call MPI_Reduce( MPI_IN_PLACE, tmasses, nv, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
     if (myrank/=0) call MPI_Reduce( tmasses,      0,       nv, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )

     if (myrank == 0) then
        dry_mass_sum = tmasses(1)
        wat_mass_sum = tmasses(2)
        if (i_co2  > 0) co2_mass_sum = tmasses(3)
        if (naddsc > 0) scp_mass_sum = tmasses(nv)
     endif

  endif
#endif

  if (myrank == 0) then

     tot_mass_sum = dry_mass_sum + wat_mass_sum

     if (firstime) then
        firstime = .false.

        dry_mass_sum0 = dry_mass_sum
        tot_mass_sum0 = tot_mass_sum
        wat_mass_sum0 = wat_mass_sum
        if (i_co2  > 0) co2_mass_sum0 = co2_mass_sum
        if (naddsc > 0) scp_mass_sum0 = scp_mass_sum
     endif

     write(*,'(a,3(g20.13,1x))') ' mass1: tot,wet,dry ', tot_mass_sum, wat_mass_sum, dry_mass_sum
     write(*,'(a,3(g20.13,1x))') ' percent change:    ', &
          (tot_mass_sum - tot_mass_sum0) / tot_mass_sum0 * 100._r8, &
          (wat_mass_sum - wat_mass_sum0) / wat_mass_sum0 * 100._r8, &
          (dry_mass_sum - dry_mass_sum0) / dry_mass_sum0 * 100._r8

     if (i_co2 > 0) then
        write(*,'(a,3(g20.13,1x))') ' co21 mass: ', co2_mass_sum
        write(*,'(a,3(g20.13,1x))') ' % change:   ', &
             (co2_mass_sum - co2_mass_sum0) / co2_mass_sum0 * 100._r8
     endif

     if (naddsc > 0) then
        write(*,'(a,3(g20.13,1x))') ' scp1 mass: ', scp_mass_sum
        write(*,'(a,3(g20.13,1x))') ' % change:   ', &
             (scp_mass_sum - scp_mass_sum0) / scp_mass_sum0 * 100._r8
     endif

  endif

end subroutine compute_mass_sums

!===============================================================================

subroutine mass_sum_from_rho4(mass_sum, rho4, allnodes, mask)

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, mwa, volt
  use consts_coms, only: r8

#ifdef OLAM_MPI
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use mpi_f08
#endif

  implicit none

  real(r8),          intent(out) :: mass_sum
  real,              intent(in)  :: rho4(mza,mwa)
  logical, optional, intent(in)  :: allnodes
  logical, optional, intent(in)  :: mask(mwa)

  integer :: j, iw, k
  logical :: iall

  iall = .false.
  if (present(allnodes)) iall = allnodes

  mass_sum = 0.0_r8

  !$omp parallel do private(iw,k) reduction(+:mass_sum)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Option to maskout certain columns
     if (present(mask)) then
        if (mask(iw)) cycle
     endif

     do k = lpw(iw),mza
        mass_sum = mass_sum + rho4(k,iw) * volt(k,iw)
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (iall) then
        call MPI_Allreduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD )
     else
        if (myrank==0) call MPI_Reduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
        if (myrank/=0) call MPI_Reduce( mass_sum,     0,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
     endif
  endif
#endif

end subroutine mass_sum_from_rho4

!===============================================================================

subroutine mass_sum_from_rho8(mass_sum, rho8, allnodes, mask)

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, mwa, volt
  use consts_coms, only: r8

#ifdef OLAM_MPI
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use mpi_f08
#endif

  implicit none

  real(r8),          intent(out) :: mass_sum
  real(r8),          intent(in)  :: rho8(mza,mwa)
  logical, optional, intent(in)  :: allnodes
  logical, optional, intent(in)  :: mask(mwa)

  integer :: j, iw, k
  logical :: iall

  iall = .false.
  if (present(allnodes)) iall = allnodes

  mass_sum = 0.0_r8

  !$omp parallel do private(iw,k) reduction(+:mass_sum)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Option to maskout certain columns
     if (present(mask)) then
        if (mask(iw)) cycle
     endif

     do k = lpw(iw),mza
        mass_sum = mass_sum + rho8(k,iw) * volt(k,iw)
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (iall) then
        call MPI_Allreduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD )
     else
        if (myrank==0) call MPI_Reduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
        if (myrank/=0) call MPI_Reduce( mass_sum,     0,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
     endif
  endif
#endif

end subroutine mass_sum_from_rho8

!===============================================================================

subroutine mass_sum_from_mixrat44(mass_sum, rho4, rr4, allnodes)

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, mwa, volt
  use consts_coms, only: r8

#ifdef OLAM_MPI
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use mpi_f08
#endif

  implicit none

  real(r8),          intent(out) :: mass_sum
  real,              intent(in)  :: rho4(mza,mwa)
  real,              intent(in)  :: rr4 (mza,mwa)
  logical, optional, intent(in)  :: allnodes

  integer :: j, iw, k
  logical :: iall

  iall = .false.
  if (present(allnodes)) iall = allnodes

  mass_sum = 0.0_r8

  !$omp parallel do private(iw,k) reduction(+:mass_sum)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw),mza
        mass_sum = mass_sum + rho4(k,iw) * rr4(k,iw) * volt(k,iw)
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (iall) then
        call MPI_Allreduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD )
     else
        if (myrank==0) call MPI_Reduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
        if (myrank/=0) call MPI_Reduce( mass_sum,     0,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
     endif
  endif
#endif

end subroutine mass_sum_from_mixrat44

!===============================================================================

subroutine mass_sum_from_mixrat84(mass_sum, rho8, rr4, allnodes)

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, mwa, volt
  use consts_coms, only: r8

#ifdef OLAM_MPI
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use mpi_f08
#endif

  implicit none

  real(r8),          intent(out) :: mass_sum
  real(r8),          intent(in)  :: rho8(mza,mwa)
  real,              intent(in)  :: rr4 (mza,mwa)
  logical, optional, intent(in)  :: allnodes

  integer :: j, iw, k
  logical :: iall

  iall = .false.
  if (present(allnodes)) iall = allnodes

  mass_sum = 0.0_r8

  !$omp parallel do private(iw,k) reduction(+:mass_sum)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw),mza
        mass_sum = mass_sum + rho8(k,iw) * rr4(k,iw) * volt(k,iw)
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (iall) then
        call MPI_Allreduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD )
     else
        if (myrank==0) call MPI_Reduce( MPI_IN_PLACE, mass_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
        if (myrank/=0) call MPI_Reduce( mass_sum,     0,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD )
     endif
  endif
#endif

end subroutine mass_sum_from_mixrat84

!===============================================================================

end module check_nan
