subroutine thermo()

  use mem_ijtabs, only: jtab_w, istp, mrl_endl, jtw_prog
  use micro_coms, only: miclevel

  implicit none

  integer iw,j

  if (mrl_endl(istp) > 0) then

     ! Horizontal loop over W/T points
     !$omp parallel do private (iw)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (miclevel <= 1) then
           call drythrm(iw)
        elseif (miclevel == 2) then
           call satadjst(iw)
        elseif (miclevel == 3) then
           call wetthrm3(iw)
        else
           stop 'Thermo option not supported...MICLEVEL out of bounds'
        endif

     enddo
     !$omp end parallel do

  endif

end subroutine thermo

!===============================================================================

subroutine drythrm(iw)

! This routine calculates theta and rv for the case where no condensate is
! allowed.

  use mem_basic,  only: theta, thil, tair, rr_v, rr_w, press
  use micro_coms, only: miclevel
  use mem_grid,   only: lpw, mza
  use consts_coms,only: p00i, rocp

  implicit none

  integer, intent(in) :: iw

  integer k

  do k = lpw(iw),mza
     theta(k,iw) = thil(k,iw)
     tair (k,iw) = thil(k,iw) * (real(press(k,iw)) * p00i) ** rocp
  enddo

  if (miclevel == 1) then
     do k = lpw(iw),mza
        rr_v(k,iw) = rr_w(k,iw)
     enddo
  endif

end subroutine drythrm

!===============================================================================

subroutine satadjst(iw)

! This routine diagnoses theta, rr_v, and rr_c using a saturation adjustment
! for the case when rr_c is the only allowed condensate

  use mem_basic,   only: thil, theta, tair, rho, rr_w, rr_v, press
  use mem_micro,   only: rr_c
  use consts_coms, only: p00i, cp, alvlocp, rocp
  use mem_grid,    only: lpw, mza
  use therm_lib,   only: rhovsl
  implicit none

  integer, intent(in) :: iw

  integer :: iterate, k
  real    :: temp(mza), t_il(mza), exner(mza), cond(mza), rhoi(mza), rhovs

  do k = lpw(iw), mza
     exner(k) = (real(press(k,iw)) * p00i) ** rocp  ! Defined WITHOUT CP factor
     t_il (k) = thil(k,iw) * exner(k)
     temp (k) = t_il(k)
     rhoi (k) = 1.0 / real(rho(k,iw))
  enddo

  do iterate = 1, 20
     do k = lpw(iw), mza
        rhovs   = rhovsl(temp(k)-273.15)
        cond(k) = max(0., rr_w(k,iw) - rhovs * rhoi(k))
        temp(k) = 0.7 * temp(k)  &
                + 0.3 * t_il(k) * (1. + alvlocp * cond(k) / max(temp(k),253.))
     enddo
  enddo

  do k = lpw(iw), mza
     rr_c (k,iw) = cond(k)
     rr_v (k,iw) = rr_w(k,iw) - cond(k)
     theta(k,iw) = temp(k) / exner(k)
     tair (k,iw) = temp(k)
  enddo

end subroutine satadjst

!===============================================================================

subroutine wetthrm3(iw)

! This routine calculates theta and rr_v for "miclevel 3 microphysics"
! given prognosed theta_il, cloud, drizzle, rain, pristine ice, snow,
! aggregates, graupel, hail, q6, and q7.

  use mem_basic,   only: press, theta, thil, tair, rr_v, rr_w
  use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, q6, q7
  use micro_coms,  only: jnmb, rxmin
  use consts_coms, only: p00i, rocp, alvl, alvi, cpi4, cp253i
  use mem_grid,    only: mza, lpw
  use therm_lib,   only: qtc

  implicit none

  integer, intent(in) :: iw

  integer :: k
  real :: tcoal,fracliq,tairstr

  real :: exner   (mza)  ! automatic array
  real :: til     (mza)  ! automatic array
  real :: totliq  (mza)  ! automatic array
  real :: totice  (mza)  ! automatic array
  real :: qhydm   (mza)  ! automatic array

  do k = lpw(iw),mza
     exner(k) = (real(press(k,iw)) * p00i) ** rocp  ! exner WITHOUT CP factor
     til(k) = thil(k,iw) * exner(k)          ! ice-liquid temperature T_il
     totliq(k) = 0.                          ! total liquid spec. density
     totice(k) = 0.                          ! total ice mixing ratio
  enddo

  if (jnmb(1) >= 1) then
     do k = lpw(iw),mza
        totliq(k) = totliq(k) + rr_c(k,iw)
     enddo
  endif

  if (jnmb(2) >= 1) then
     do k = lpw(iw),mza
        totliq(k) = totliq(k) + rr_r(k,iw)
     enddo
  endif

  if (jnmb(3) >= 1) then
     do k = lpw(iw),mza
        totice(k) = totice(k) + rr_p(k,iw)
     enddo
  endif

  if (jnmb(4) >= 1) then
     do k = lpw(iw),mza
        totice(k) = totice(k) + rr_s(k,iw)
     enddo
  endif

  if (jnmb(5) >= 1) then
     do k = lpw(iw),mza
        totice(k) = totice(k) + rr_a(k,iw)
     enddo
  endif

  if (jnmb(6) >= 1) then
     do k = lpw(iw),mza
        if (rr_g(k,iw) > rxmin(6)) then
           call qtc(q6(k,iw)/rr_g(k,iw),tcoal,fracliq)
           totliq(k) = totliq(k) + rr_g(k,iw) * fracliq
           totice(k) = totice(k) + rr_g(k,iw) * (1. - fracliq)
        endif
     enddo
  endif

  if (jnmb(7) >= 1) then
     do k = lpw(iw),mza
        if (rr_h(k,iw) > rxmin(7)) then
           call qtc(q7(k,iw)/rr_g(k,iw),tcoal,fracliq)
           totliq(k) = totliq(k) + rr_h(k,iw) * fracliq
           totice(k) = totice(k) + rr_h(k,iw) * (1. - fracliq)
        endif
     enddo
  endif

  if (jnmb(8) >= 1) then
     do k = lpw(iw),mza
        totliq(k) = totliq(k) + rr_d(k,iw)
     enddo
  endif

  do k = lpw(iw),mza
     qhydm(k) = alvl * totliq(k) + alvi * totice(k)
     rr_v(k,iw) = rr_w(k,iw) - totliq(k) - totice(k)
  enddo

  do k = lpw(iw),mza
     tairstr = .5 * (til(k) + sqrt( til(k) * (til(k) + cpi4 * qhydm(k)) ) )
     if (tairstr < 253.0) tairstr = til(k) * (1. + qhydm(k) * cp253i)

     theta(k,iw) = tairstr / exner(k)
     tair (k,iw) = tairstr
  enddo

end subroutine wetthrm3
