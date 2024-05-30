module umwm_stokes

  implicit none

  real, allocatable  :: dstok     (:) ! stokes drift e-folding depth [m]
  real, allocatable  :: ustok   (:,:) ! stokes drift x-component [m/s]
  real, allocatable  :: vstok   (:,:) ! stokes drift y-component [m/s]
  real, allocatable  :: ustokmag(:,:) ! stokes drift velocity magnitude [m/s]
  real, allocatable  :: util  (:,:,:)
  real, allocatable  :: arg   (:,:,:)

!Bob: The following depths copied from STOKES namelist in main.nml
!  integer, parameter :: lm = 42     ! number of values in depths array
!  real :: depths(lm) = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, &
!                         1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, &
!                         4.5, 5.0,  6.,  7.,  8.,  9., 10., 12., 14., 16., &
!                         18., 20., 24., 28., 32., 40., 50., 60., 70., 80., &
!                         90., 100. ]
!Bob: To speed up the code, compute Stokes drift velocity only for 1 layer, which is at the surface.
!     Only this value is used in computing the relative velocity between wind and water.
  integer, parameter :: lm = 1     ! number of values in depths array
  real :: depths(lm) = [ 0.1 ]

contains

!===============================================================================

subroutine stokes_drift(option)

  ! Computes wave-induced Stokes drift

  use umwm_module, only: eulerinv, umwm, evs, k, dwn, dth, om, pm, cth, sth
  use mem_sea, only: msea
  use consts_coms, only: pi2

  character(4), intent(in), optional :: option
  integer :: i, l, o, p
  real :: ustok_efolding
  real, allocatable :: kd(:,:)

  if (present(option)) then
     if (option == 'init') then

        ! allocate stokes velocities and utility array
        allocate(ustok   (lm,msea)); ustok(:,:) = 0.
        allocate(vstok   (lm,msea)); vstok(:,:) = 0.
        allocate(ustokmag(lm,msea))
        allocate(dstok      (msea))
        allocate(util (om,lm,msea))
        allocate(arg  (om,lm,msea))
        allocate(kd   (om,   msea))

        do concurrent (o=1:om, i=2:msea) 
           kd(o,i) = k(o,i) * umwm%seadep(i)
        end do

        ! compute exponent
        do l = 1, lm
           do i = 2,msea
              do o = 1, om

                 arg(o,l,i) = 2 * k(o,i) * (depths(l) + umwm%seadep(i))

                 if(abs(arg(o,l,i)) > 50 .or. kd(o,i) > 50) then
                    ! hyperbolic trig. functions would overflow;
                    ! use deep water approximation instead
                    util(o,l,i) = pi2 * umwm%f(o) * 2 * k(o,i)**2 &
                                * exp(2 * k(o,i) * depths(l)) * dwn(o,i) * dth
                 else
                    ! first order approximation for arbitrary depth
                    util(o,l,i) = pi2 * umwm%f(o) * k(o,i)**2 &
                                * cosh(2 * k(o,i) * (depths(l) + umwm%seadep(i)))&
                                / sinh(kd(o,i))**2 * dwn(o,i) * dth
                 endif

              enddo

              if (abs(depths(l)) > umwm%seadep(i)) util(:,l,i) = 0

           enddo
        enddo

        deallocate(arg, kd)

        write(*,'(a)') 'umwm: stokes_drift: initialized'

     endif ! option == 'init'
  endif ! present(option)

  ! stokes velocities
  !$omp parallel do private(l,p,o,ustok_efolding)
  do i = 2,msea
     dstok(i) = 0.

     do l = 1, lm
        ustok(l,i) = 0.
        vstok(l,i) = 0.
        do p = 1, pm
           do o = 1, om
              ustok(l,i) = ustok(l,i) + util(o,l,i) * evs(o,p,i) * cth(p)
              vstok(l,i) = vstok(l,i) + util(o,l,i) * evs(o,p,i) * sth(p)
           enddo
        enddo

        ! Stokes drift magnitude
        ustokmag(l,i) = sqrt(ustok(l,i)**2 + vstok(l,i)**2)
     enddo

     ! Stokes e-folding depth
     if (ustokmag(1,i) < 1.e-6) cycle

     ustok_efolding = ustokmag(1,i) * eulerinv

     depth_loop: do l = 2,lm
        if (ustokmag(l,i) < ustok_efolding) then
           dstok(i) = (abs(ustokmag(l-1,i) - ustok_efolding) * depths(l)  &
                    +  abs(ustokmag(l  ,i) - ustok_efolding) * depths(l-1)) &
                    /     (ustokmag(l-1,i) - ustokmag(l,i))

           exit depth_loop
        endif
     enddo depth_loop
     dstok(i) = -dstok(i)
  enddo
  !$omp end parallel do

end subroutine stokes_drift

end module umwm_stokes
