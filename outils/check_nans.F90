module check_nan

Contains

!===============================================================================

subroutine check_nans(icall,rvara1,rvara2)

  use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,vmc,vc,vxe,vye,vze
  use mem_micro,  only: sh_c,sh_r
  use mem_grid,   only: mza,mwa,lpw,volti,mva,lpv,glatw,glonw,glatv,glonv
  use mem_tend,   only: thilt,sh_wt,vmxet,vmyet,vmzet
  use misc_coms,  only: io6, iparallel
  use mem_ijtabs, only: itab_v, itab_w
  use mem_para,   only: myrank

#ifdef IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif

  implicit none

  integer :: i,k,istop
  integer, intent(in) :: icall

  real, optional, intent(in) :: rvara1 (mza,mwa)
  real, optional, intent(in) :: rvara2 (mza,mwa)

#ifdef IEEE_ARITHMETIC

  istop = 0

  do i = 2,mwa
     do k = lpw(i),mza
        if ( ieee_is_nan(sh_w (k,i)) .or.  &
             ieee_is_nan(rho  (k,i)) .or.  &
             ieee_is_nan(thil (k,i)) .or.  &
             ieee_is_nan(press(k,i)) .or.  &
             ieee_is_nan(wc   (k,i)) .or.  &
             ieee_is_nan(wmc  (k,i)) .or.  &
             ieee_is_nan(thilt(k,i)) .or.  &
             ieee_is_nan(sh_wt(k,i)) .or.  &
             ieee_is_nan(vmxet(k,i)) .or.  &
             ieee_is_nan(vmyet(k,i)) .or.  &
             ieee_is_nan(vmzet(k,i)) .or.  &
             ieee_is_nan(vxe  (k,i)) .or.  &
             ieee_is_nan(vye  (k,i)) .or.  &
             ieee_is_nan(vze  (k,i)) .or.  &
             ieee_is_nan(sh_v (k,i)) ) then
           
           write(*,*) 'Node ', myrank
           write(*,*) 'NaN at: k, i, icall ', k, i, icall
           write(io6,*) ''
           write(io6,*) 'check_nans: k, i, icall, lat, lon ',k,i,icall,glatw(i),glonw(i)
           write(io6,*) 'sh_w,rho ',sh_w(k,i),rho(k,i)
           write(io6,*) 'thil,press ',thil(k,i),press(k,i)
           write(io6,*) 'wc, wmc, sh_v ',wc(k,i),wmc(k,i),sh_v(k,i)
           write(io6,*) 'vxe, vye, vze ',vxe(k,i),vye(k,i),vze(k,i)
           write(io6,*) 'thilt, sh_wt ',thilt(k,i),sh_wt(k,i)
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

#endif
end subroutine check_nans



subroutine check_pos(icall)

  use mem_ijtabs, only: jtab_w, jtw_prog, itab_w, istp, mrl_endl
  use mem_grid,   only: lpw, mza
  use var_tables, only: num_scalar, scalar_tab
  use oname_coms, only: nl
  use mem_para,   only: myrank

  implicit none

  integer, intent(in) :: icall
  
  integer :: j, iw, k, mrl, n

  if (nl%iscal_monot < 1) return

  ! Return if this is not the end of the long timestep on any MRL
  mrl = mrl_endl(istp)
  if (mrl == 0) return

  do n = 1, num_scalar

     if (.not. scalar_tab(n)%pdef) cycle

     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
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




subroutine compute_mass_sums()

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw, mza, volt
  use mem_para,    only: myrank, mgroupsize
  use misc_coms,   only: naddsc, iparallel
  use mem_addsc,   only: addsc
  use consts_coms, only: r8
  use mem_basic,   only: rho, sh_w
#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real(r8) :: tot_mass_sum
  real(r8) :: dry_mass_sum
  real(r8) :: wat_mass_sum
  real(r8) :: sclp_mass_sum

  real(r8) :: tmasses(mgroupsize)
  integer  :: j, iw, k, ier

  logical,  save :: firstime = .true.
  real(r8), save :: tot_mass_sum0
  real(r8), save :: dry_mass_sum0
  real(r8), save :: wat_mass_sum0
  real(r8), save :: sclp_mass_sum0

  tot_mass_sum  = 0.0_r8
  wat_mass_sum  = 0.0_r8
  sclp_mass_sum = 0.0_r8

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw),mza

        tot_mass_sum = tot_mass_sum + rho(k,iw) * volt(k,iw)
        wat_mass_sum = wat_mass_sum + rho(k,iw) * volt(k,iw) * sh_w(k,iw)
        
        if (naddsc > 0) then
           sclp_mass_sum = sclp_mass_sum + rho(k,iw) * volt(k,iw) * addsc(1)%sclp(k,iw)
        endif

     enddo
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
 
     call MPI_Gather( tot_mass_sum, 1, MPI_DOUBLE, &
                      tmasses,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ier )
     
     if (myrank == 0) tot_mass_sum = sum(tmasses)

     call MPI_Gather( wat_mass_sum, 1, MPI_DOUBLE, &
                      tmasses,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ier )
     
     if (myrank == 0) wat_mass_sum = sum(tmasses)

     if (naddsc > 0) then

        call MPI_Gather( sclp_mass_sum, 1, MPI_DOUBLE, &
                         tmasses,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ier )
     
        if (myrank == 0) sclp_mass_sum = sum(tmasses)
     endif

  endif
#endif

  if (myrank == 0) then

     dry_mass_sum = tot_mass_sum - wat_mass_sum

     if (firstime) then
        dry_mass_sum0 = dry_mass_sum
        tot_mass_sum0 = tot_mass_sum
        wat_mass_sum0 = wat_mass_sum
        if (naddsc > 0) sclp_mass_sum0 = sclp_mass_sum
     endif

     write(*,'(a,3(g20.13,1x))') ' mass1: tot,wet,dry ', tot_mass_sum, wat_mass_sum, dry_mass_sum
     write(*,'(a,3(g20.13,1x))') ' percent change:    ', &
          (tot_mass_sum - tot_mass_sum0) / tot_mass_sum0 * 100._r8, &
          (wat_mass_sum - wat_mass_sum0) / wat_mass_sum0 * 100._r8, &
          (dry_mass_sum - dry_mass_sum0) / dry_mass_sum0 * 100._r8

     if (naddsc > 0) then
        write(*,'(a,3(g20.13,1x))') ' sclp1 mass: ', sclp_mass_sum
        write(*,'(a,3(g20.13,1x))') ' % change:   ', &
             (sclp_mass_sum - sclp_mass_sum0) / sclp_mass_sum0 * 100._r8
     endif

  endif

  firstime = .false.

end subroutine compute_mass_sums

end module check_nan
