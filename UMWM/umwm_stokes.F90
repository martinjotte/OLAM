module umwm_stokes

  implicit none

  real, allocatable  :: ustoksfc(:) ! stokes drift x-component [m/s]
  real, allocatable  :: vstoksfc(:) ! stokes drift y-component [m/s]
  real, allocatable  :: util(:,:,:) ! utility array
!!real, allocatable  :: dstok   (:) ! unused

!!Bob: The following depths copied from STOKES namelist in main.nml
!! integer, parameter :: lm = 42     ! number of values in depths array
!! real :: depths(lm) = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, &
!!                        1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, &
!!                        4.5, 5.0,  6.,  7.,  8.,  9., 10., 12., 14., 16., &
!!                        18., 20., 24., 28., 32., 40., 50., 60., 70., 80., &
!!                        90., 100. ]

!Bob: To speed up the code, compute Stokes drift velocity only for 1 layer, which is at the surface.
!     Only this value is used in computing the relative velocity between wind and water.

  integer, parameter :: lm = 1     ! number of values in depths array
  real :: depths(lm) = [ 0.1 ]

contains

!===============================================================================

subroutine stokes_init()

  ! Computes wave-induced Stokes drift

  use umwm_module, only: umwm, freq, k, dwn, dth, om
  use mem_sea,     only: msea
  use consts_coms, only: pi4, pi2
  use misc_coms,   only: io6

  implicit none

  integer :: i, l, o
  real    :: arg, kd

  ! allocate stokes velocities and utility array

  allocate(ustoksfc  (msea))
  allocate(vstoksfc  (msea))
  allocate(util(om,lm,msea))
!!allocate(dstok     (msea))

  !$omp parallel do private(l,o,arg,kd)
  do i = 2, msea

     ustoksfc(i) = 0.
     vstoksfc(i) = 0.

     ! compute exponent
     do l = 1, lm

        if (depths(l) > umwm%seadep(i)) then

           util(:,l,i) = 0.0

        else

           do o = 1, om
              arg = 2. * k(o,i) * (depths(l) + umwm%seadep(i))
              kd  = k(o,i) * umwm%seadep(i)

              if ( arg > 50. .or. kd > 50.) then
                 ! hyperbolic trig. functions would overflow;
                 ! use deep water approximation instead
                 util(o,l,i) = pi4 * freq(o) * k(o,i)**2 &
                             * exp(2. * k(o,i) * depths(l)) * dwn(o,i) * dth
              else
                 ! first order approximation for arbitrary depth
                 util(o,l,i) = pi2 * freq(o) * k(o,i)**2 &
                             * cosh(arg) / sinh(kd)**2 * dwn(o,i) * dth
              endif
           enddo

        endif
     enddo

  enddo
  !$omp end parallel do

  write(io6,*)
  write(io6,*) 'umwm: stokes_drift: initialized'

end subroutine stokes_init

!===============================================================================

subroutine stokes_drift(i)

  ! Computes wave-induced Stokes drift

  use umwm_module, only: evs, om, pm, cth, sth, umwm

  implicit none

  integer, intent(in) :: i
  integer             :: o
  real                :: evxsum(om), evysum(om)
! real                :: ustokmag(lm), ustok_efolding

  ! Stokes velocities at sea surface

  if (umwm%iactive(i)) then

     do o = 1, om
        evxsum(o) = sum( evs(1:pm,o,i) * cth(1:pm) )
        evysum(o) = sum( evs(1:pm,o,i) * sth(1:pm) )
     enddo

     ustoksfc(i) = sum( evxsum(1:om) * util(1:om,1,i) )
     vstoksfc(i) = sum( evysum(1:om) * util(1:om,1,i) )

  else

     ustoksfc(i) = 0.
     vstoksfc(i) = 0.

  endif

  ! Stokes e-folding depth not used

! dstok(i) = 0.
!
! do l = 1, lm
!    ustokmag(l) = sqrt( ustok(l,i)**2 + vstok(l,i)**2)
! enddo
!
! if (ustokmag(1,i) >= 1.e-6) then
!    ustok_efolding = ustokmag(1,i) * eulerinv
!
!    depth_loop: do l = 2, lm
!      if (ustokmag(l,i) < ustok_efolding) then
!          dstok(i) = (abs(ustokmag(l-1,i) - ustok_efolding) * depths(l)  &
!                   +  abs(ustokmag(l  ,i) - ustok_efolding) * depths(l-1)) &
!                   /     (ustokmag(l-1,i) - ustokmag(l,i))
!          exit depth_loop
!       endif
!    enddo depth_loop
!    dstok(i) = -dstok(i)
! enddo

end subroutine stokes_drift

end module umwm_stokes
