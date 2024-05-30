module umwm_physics

  use umwm_module

  implicit none

contains

!===============================================================================

subroutine source()

  use misc_coms,   only: time_istp8
  use consts_coms, only: r8
  use mem_sea,     only: msea

  integer :: i,o,p
  real :: physics_time_step

  !$omp parallel do private(o,p)
  do i = 2,msea
     do p = 1,pm
        do o = 1,oc(i)
           ! calculate the exponential argument:
           evsf(o,p,i) = ssin(o,p,i) - sds(o,p,i) * snl_arg(o,i) &
                       - sbf(o,i) - sdt(o,i) - sdv(o,i) + sice(o,i)
        enddo
        evsf(oc(i)+1:om,p,i) = 0.
     enddo
  enddo
  !$omp end parallel do

  ! Compute maximum physics time step for diagnostics.
  physics_time_step = explim / maxval(abs(evsf))

  ! Set timestep to value from OLAM

  dta = dt_olam

  ! if first time step, set time step to zero and integrate only the diagnostic part
  if (time_istp8 < 1.e-3_r8) dta = 0.

  if (dta > physics_time_step) then
     write(6,'(a,2f10.1)') 'dta exceeds physics_time_step ', dta, physics_time_step
     write(6,'(a)')        'Reducing dta to physics time step in subroutine source'
     dta = physics_time_step
  endif

  !$omp parallel do private(o,p)
  do i = 2,msea
     do p = 1,pm
        ! integrate source terms for the prognostic range (o <= ol)
        do o = 1,oc(i)
           evsf(o,p,i) = evs(o,p,i) * exp(dta * (ssin(o,p,i) - sds(o,p,i) - sbf(o,i) &
                                               - sdt(o,i) - sdv(o,i) + sice(o,i))) &
                       + dta * snl(o,p,i)
        enddo

        ! integrate source terms for the diagnostic range (o > ol)
        do o = oc(i) + 1, om
           if (ssin(o,p,i) - sdt(o,i) - sdv(o,i) + sice(o,i) >= 0) then
              evsf(o,p,i) = oneoverk4(o,i) * ((ssin(o,p,i) - sdt(o,i) - sdv(o,i) + sice(o,i)) &
                          / (twopisds_fac * umwm%f(o) * dummy(o,p,i) * cothkd(o,i)))**inv_sds_power
           endif
        enddo
     enddo
     evs(:,:,i) = 0.5 * (evs(:,:,i) + evsf(:,:,i))
  enddo
  !$omp end parallel do

end subroutine source

!===============================================================================

subroutine umwm_diag()

  use mem_sea,     only: msea, omsea
  use mem_sfcg,    only: sfcg
  use consts_coms, only: pi2

  integer :: o, p, i, opeak, ppeak

  real :: mag, xcomp, ycomp
  real :: m0(msea)
  real :: m2(msea)
  real :: evskdk(om,pm)

  real :: evskdkovcp, evskdkmax

  !$omp parallel do private(evskdkmax, p, o, evskdk, opeak, ppeak, evskdkovcp, xcomp, ycomp, mag)
  do i = 2,msea

     m0   (i) = 0.
     m2   (i) = 0.
     momx (i) = 0.
     momy (i) = 0.
     cgmxx(i) = 0.
     cgmxy(i) = 0.
     cgmyy(i) = 0.

     evskdkmax = -1.e6
     do p = 1,pm
        do o = om,1,-1
           evskdk(o,p) = evs(o,p,i) * kdk(o,i)
           m0(i) = m0(i) + evskdk(o,p)
           m2(i) = m2(i) + umwm%f(o)**2 * evskdk(o,p)

           if (evskdkmax < evskdk(o,p)) then
              evskdkmax = evskdk(o,p)
              opeak = o
              ppeak = p
           endif
        enddo
     enddo
     mwp(i) = sqrt(m0(i) / (m2(i) + tiny(m2))) ! mean wave period

     ! total wave momentum:
     do p = 1,pm
        do o = 1,oc(i)
           evskdkovcp = evskdk(o,p) * invcp0(o,i)

           momx(i) = momx(i) + evskdkovcp * cth(p)
           momy(i) = momy(i) + evskdkovcp * sth(p)

           cgmxx(i) = cgmxx(i) + cg0(o,i) * evskdkovcp * cth(p)**2
           cgmxy(i) = cgmxy(i) + cg0(o,i) * evskdkovcp * cth(p) * sth(p)
           cgmyy(i) = cgmyy(i) + cg0(o,i) * evskdkovcp * sth(p)**2
        enddo
     enddo

     momx(i)  = momx (i) * rhosw * dthg
     momy(i)  = momy (i) * rhosw * dthg
     cgmxx(i) = cgmxx(i) * rhosw * dthg
     cgmxy(i) = cgmxy(i) * rhosw * dthg
     cgmyy(i) = cgmyy(i) * rhosw * dthg

     ! significant wave height:
     swh(i) = 0.
     do p = 1,pm
        do o = 1,om
           swh(i) = swh(i) + evskdk(o,p)
        enddo
     enddo
     swh(i) = 4. * sqrt(swh(i) * dth)

     ! mean wave direction:
     xcomp = 0.
     ycomp = 0.
     do p = 1,pm
        mag   = sum(evskdk(:,p))
        xcomp = xcomp + mag * cth(p)
        ycomp = ycomp + mag * sth(p)
     enddo
     mwd(i) = atan2(ycomp,xcomp)

     ! wavenumber spectrum moments:
     m0(i) = 0.
     m2(i) = 0.
     do p = 1,pm
        do o = 1,om
           m0(i) = m0(i) + evskdk(o,p)
           m2(i) = m2(i) + evs(o,p,i) * k3dk(o,i)
        enddo
     enddo

     mss(i) = m2(i) * dth ! mean-squared slope

     mwl(i) = pi2 * sqrt(m0(i) / (m2(i) + tiny(m2(i)))) ! mean wavelength

     dwd(i) = th(ppeak)                       ! dominant wave direction
     dwl(i) = pi2 / k(opeak,i)                ! dominant wave length
     dwp(i) = 1. / umwm%f(opeak)              ! dominant wave period

     dcp0(i) = cp0(opeak,i)                   ! dominant phase speed, intrinsic
     dcg0(i) = cg0(opeak,i)                   ! dominant group speed, intrinsic

     dcp(i) = dcp0(i) + ucurr(i) * cth(ppeak) + vcurr(i) * sth(ppeak) ! dominant phase speed
     dcg(i) = dcg0(i) + ucurr(i) * cth(ppeak) + vcurr(i) * sth(ppeak) ! dominant group speed

if (i + omsea == 91457) then 
  print*, ' '
  write(6,'(190a)') '   dcp    dcg    dcp0   dcg0    dwp   dwl    dwd    mwl    mwd    swh    mwp     mss       momx      momy umwm%taux  umwm%tauy  umwm%ustar  sfcg%ustar umwm%wspd   umwm%wdir'
  write(6,'(11f7.2,f8.3,9f10.2)') dcp(i), dcg(i), dcp0(i), dcg0(i), dwp(i), dwl(i), dwd(i), &
                                  mwl(i), mwd(i), swh(i), mwp(i), mss(i), momx(i), momy(i), &
                                 umwm%taux(i), umwm%tauy(i), umwm%ustar(i), sfcg%ustar(i+omsea), umwm%wspd(i), umwm%wdir(i)
  print*, ' '
endif

  enddo
  !$omp end parallel do

end subroutine umwm_diag

!======================================================================

pure function meanwaveperiod(i) result(mwp)

  ! given a spatial grid index i, returns mean wave period at that
  ! location.

  use umwm_module, only: evs, umwm, kdk, om, pm

  integer, intent(in) :: i

  integer :: o, p
  real    :: m0, m2, mwp

  m0 = 0
  m2 = 0
  do p = 1,pm
     do o = 1,om
        m0 = m0 + evs(o,p,i) * kdk(o,i)
        m2 = m2 + umwm%f(o)**2 * evs(o,p,i) * kdk(o,i)
     enddo
  enddo
  mwp = sqrt(m0 / m2)

end function meanwaveperiod

end module umwm_physics
