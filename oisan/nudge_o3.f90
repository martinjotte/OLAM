subroutine nudge_prep_o3(iaction, o_ozone)

  use mem_nudge,  only: ozone_obsp, ozone_obsf
  use mem_grid,   only: mza, mwa
  use mem_ijtabs, only: jtab_w, jtw_init
  use misc_coms,  only: i_o3

  implicit none

  integer, intent(in) :: iaction
  real,    intent(in) :: o_ozone(mza,mwa)
  integer             :: j, iw, k

  ! Do nothing if we do not prognose ozone

  if (i_o3 == 0) return

! Swap future data time into past data time if necessary.

  if (iaction == 1) then

     !$omp parallel do private(iw,k)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
        do k = 2, mza
           ozone_obsp(k,iw) = ozone_obsf(k,iw)
        enddo
     enddo
     !omp end parallel do

  endif

! Fill future observational arrays

  !$omp parallel do private(iw,k)
  do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
     do k = 2, mza
        ozone_obsf(k,iw) = o_ozone(k,iw)
     enddo
  enddo
  !omp end parallel do

end subroutine nudge_prep_o3

!===========================================================================

subroutine obs_nudge_o3()

  use mem_nudge,  only: ozone_obsp, ozone_obsf, o3nudflag, tnudi_o3, o3nudpress
  use mem_basic,  only: rho, press
  use mem_grid,   only: mza, lpw
  use mem_ijtabs, only: jtab_w, jtw_prog
  use isan_coms,  only: ifgfile, s1900_fg
  use misc_coms,  only: s1900_sim, i_o3
  use var_tables, only: scalar_tab

  implicit none

  integer :: k, j, iw, kbot
  real    :: tp, tf

  ! Skip if we don't nudge

  if (o3nudflag == 0) return

  ! Do nothing if we do not prognose ozone

  if (i_o3 == 0) return

  ! Time interpolation coefficients

  tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
            / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

  tp = 1. - tf

  ! Compute ozone nudging tendency using point-by-point (non-spectral) information

  !$omp parallel do private(iw,k,kbot)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Find level above which we nudge ozone

     do k = mza, lpw(iw), -1
        if (press(k,iw) > o3nudpress) exit
     enddo
     kbot = k+1

     !dir$ ivdep
     do k = kbot, mza

        scalar_tab(i_o3)%var_t(k,iw) = scalar_tab(i_o3)%var_t(k,iw) + tnudi_o3 * real(rho(k,iw)) * &
             ( (tp * ozone_obsp(k,iw) + tf * ozone_obsf(k,iw)) - scalar_tab(i_o3)%var_p(k,iw) )

     enddo

  enddo
  !$omp end parallel do

end subroutine obs_nudge_o3
