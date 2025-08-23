module umwm_init

  implicit none

  private
  public :: nmlassign, init

contains

!===============================================================================

subroutine nmlassign()

  use umwm_module, only: pm, fmin, fmax, fprog, gustiness, dmin
  use misc_coms,   only: io6

  implicit none

  logical :: namelistok

  ! check namelist values:

  namelistok = .true.

  if (mod(pm,4) /= 0) then
     write(io6,*) 'ERROR (umwm namelist) - pm must be divisible by 4'
     namelistok = .false.
  elseif (mod(pm,8) /= 0) then
     write(io6,*) 'WARNING (umwm namelist) - pm should be divisible by 8 for optimal propagation properties'
  endif

  if (fmin <= 0 .or. fmax <= 0 .or. fprog <= 0) then
     write(io6,*) 'ERROR (umwm namelist) - fmin, fmax and fprog must be > 0'
     namelistok = .false.
  endif

  if (fmin >= fmax) then
     write(io6,*) 'ERROR (umwm namelist) - fmin must be < fmax'
     namelistok = .false.
  endif

  if (fprog > fmax) then
     write(io6,*) 'ERROR (umwm namelist) - fprog must be <= fmax; using highest allowed value'
     namelistok = .false.
  endif

  if (gustiness < 0) then
     write(io6,*) 'ERROR (umwm namelist) - gustiness must be positive'
     namelistok = .false.
  elseif (gustiness > 0.2) then
     write(io6,*) 'WARNING (umwm namelist) - gustiness > 0.2; proceed with caution'
  endif

  if (dmin <= 0) then
     write(io6,*) 'ERROR (umwm namelist) - dmin must be > 0'
     namelistok = .false.
  endif

  if (.not. namelistok) stop 'STOP umwm namelist'

end subroutine nmlassign

!===============================================================================

subroutine init()

  ! Initialize model variables such as frequencies, direction angles,
  ! phase speed and group velocity, wave numbers, etc.

  use consts_coms,    only: piu180, erad
  use mem_sfcg,       only: sfcg, mvsfc, itab_vsfc
  use mem_sea,        only: msea, omsea
  use umwm_oforcing,  only: oforcing
  use umwm_module,    only: om, pm, umwm, freq, fsdsf, fmin, dlnf, th, cth, sth, cp0, k, &
                            cg0, cthp_dirv, sthp_dirv, cth2pp, dth, dmin, twopisds_fac
  use umwm_advection, only: init_refraction, init_propagation
  use misc_coms,      only: io6

  implicit none

  integer :: i, o, p, pp, ind, iw1, iw2, iwsfc
  real    :: raxis, dx, dy, thdirv
  real    :: cth2(pm)

  ! Set frequency bins:
  do o = 1, om
     freq (o) = exp(log(fmin) + real(o-1) * dlnf)
     fsdsf(o) = twopisds_fac * freq(o)
  enddo

  ! Calculate wave ray directions
  do p = 1, pm
     th (p) = (real(p) - 0.5 * real(pm + 1)) * dth  ! angles {-175,-165,...,165,175 deg} but in radians
     cth(p) = cos(th(p))  ! cosines
     sth(p) = sin(th(p))  ! sines
  enddo

  !$omp parallel
  !$omp do private(iwsfc)
  do i = 2, msea
     iwsfc = i + omsea
     umwm%seadep(i) = max(dmin, sfcg%topw(iwsfc) - sfcg%bathym(iwsfc))
     umwm%alogzs(i) = log( sfcg%dzt_bot(iwsfc) )
  enddo
  !$omp end do nowait

  !$omp do private(iw1,iw2,raxis,dx,dy,thdirv,p)
  do i = 2, mvsfc
     iw1 = itab_vsfc(i)%iwn(1)
     iw2 = itab_vsfc(i)%iwn(2)

     ! Do only for sfcg%v faces adjacent to a sea cell

     if (sfcg%leaf_class(iw1) == 0 .or. sfcg%leaf_class(iw2) == 0) then

        ! Find azimuth angle of sfcg%v face normal vector (West = +/-PI azimuth)

        raxis = sqrt(sfcg%xev(i)**2 + sfcg%yev(i)**2)  ! dist from earth axis

        if (raxis > 1.e3) then
           dx = (sfcg%vny(i) * sfcg%xev(i) - sfcg%vnx(i) * sfcg%yev(i)) / raxis
           dy =  sfcg%vnz(i) * raxis / erad &
              - (sfcg%vnx(i) * sfcg%xev(i) + sfcg%vny(i) * sfcg%yev(i)) * sfcg%zev(i) / (raxis * erad)
        else
           dx = 0.
           dy = 0.
        endif

        thdirv = atan2(dy,dx)

        ! Compute projection of th(p) onto thddirv (need to consider grid curvature?)

        do p = 1,pm
           cthp_dirv(p,i) = cos(th(p) - thdirv)
           sthp_dirv(p,i) = sin(th(p) - thdirv)
        enddo

     endif
  enddo
  !$omp end do
  !$omp end parallel

! ! "left" and "right" directional indices for refraction:
! do p = 1, pm
!    pl(p) = p + 1
!    pr(p) = p - 1
! enddo
! pl(pm) = 1
! pr(1)  = pm

  do p = 1, pm
     cth2(p) = cos(dth * real(p-1))**2 * dth
  enddo

  allocate(cth2pp(pm,pm))

  do pp = 1, pm
     do p = 1, pm
        ind = pp - p + 1
        if (ind <= 0) ind = pm + ind
        cth2pp(p,pp) = cth2(ind)
     enddo
  enddo

! call oforcing()

  ! compute wave numbers, phase speeds, and group velocities:
  call dispersion(1.e-2)

  call init_propagation()
  call init_refraction()

! do i = 2, msea
!
!    ! initialize drag coefficient (Large and Pond, 1981):
!    umwm%cd(i) = max(1.2e3, .49e-3 + .065e-3 * umwm%wspd(i))
!
!    ! initialize friction velocity:
!    umwm%ustar(i) = sqrt(umwm%cd(i)) * umwm%wspd(i)
!
! enddo

  write(io6,*)
  write(io6,fmt=102)
  write(io6,*)'initialization summary:'
  write(io6,fmt=103)
  write(io6,*)'bin  f[hz]    t[s]    min(k)   max(k)   min(c)   max(c)   min(cg)   max(cg)'
  write(io6,fmt=103)
  do o=1,om
    write(io6,fmt=101)o,freq(o),1./freq(o),                              &
                    minval(k  (o,2:msea),dim=1),maxval(k  (o,2:msea),dim=1), &
                    minval(cp0(o,2:msea),dim=1),maxval(cp0(o,2:msea),dim=1), &
                    minval(cg0(o,2:msea),dim=1),maxval(cg0(o,2:msea),dim=1)
  enddo
  write(io6,fmt=102)
  write(io6,*)

  101 format(1x,i2,8(2x,f7.4))
  102 format('!',77('='),'!')
  103 format('!',77('-'),'!')

!Bob: print out initialized (om,msea) arrays for a single isea point

! i = 91457 - omsea
! print*, ' '
! write(6,*) '               o bf1_renorm bf2_renorm cp0     cg0     cothkd     dwn    fkovg     l2   logl2overz     k       k4          kdk         k3dk        sbf'
! print*, ' '
! do o = 1, om
!    write(6,'(a,i6,11f9.3,8e12.3)') 'init(o,i) ', o, &
!                      bf1_renorm(o,i), bf2_renorm(o,i),   cp0(o,i), cg0(o,i), &
!                          cothkd(o,i),        dwn(o,i), fkovg(o,i),  l2(o,i), &
!                      logl2overz(o,i),          k(o,i),    k4(o,i),  &
!                             kdk(o,i),       k3dk(o,i),   sbf(o,i)
! enddo
! print*, ' '

end subroutine init

!===============================================================================

subroutine dispersion(tol)

  use consts_coms, only: pi2, grav, gravi
  use mem_sea,     only: msea, omsea
  use mem_sfcg,    only: sfcg
  use umwm_module, only: om, k, cp0, cg0, dwn, l2, logl2overz, k4, oneoverk4, kdk, k3dk, freq, &
                         fkovg, cothkd, invcp0, bf1_renorm, bf2_renorm, umwm, dmin, sfct, sdv, &
                         rhosw, dlnf, sbf_fac, sbp_fac, nu_water, snl_fac, bf1, bf2
  use misc_coms,   only: io6

  implicit none

  ! Iteratively solve the dispersion relation by iteration,
  ! and compute phase and group velocities in absolute reference frame (w/ currents)

  real, intent(in) :: tol

  integer :: counter, i, o, iwsfc
  real    :: dk, b, t, kd
  real    :: f_nd(om)

  write(io6,*) 'UMWM dispersion: solving for dispersion relationship'

  ! set seadep(1) to the minimum water depth. It will be used for propagation
  ! and refraction when a sea cell borders land
  umwm%seadep(1) = dmin

  !$omp parallel do private(iwsfc,o,f_nd,b,counter,dk,t,kd)
  do i = 1, msea

     ! non-dimesionalize frequencies, and use deep water limit as initial guess:
     do o = 1, om
        cp0 (o,i) = pi2 * sqrt(umwm%seadep(i) / grav)
        f_nd(o)   = cp0(o,i) * freq(o)
        k   (o,i) = f_nd(o) * f_nd(o)
     enddo

     ! non-dimesionalize surface tension:
     b = sfct / (rhosw * grav * umwm%seadep(i)**2)

     do o = 1, om

        counter = 1
        dk = 2. * tol

        do while (abs(dk) > tol) ! newton-raphson iteration loop

           t = tanh( k(o,i) )

           dk = -(f_nd(o) * f_nd(o) - k(o,i) * t        &
                * (1. + b * k(o,i) * k(o,i)))           &
                / (3. * b * k(o,i) * k(o,i) * t + t     &
                + k(o,i) * (1. + b * k(o,i) * k(o,i)) * (1. - t * t))
           k(o,i) = k(o,i) - dk

           if (counter == 1000) exit ! escape if stuck
           counter = counter + 1

        enddo

        k(o,i) = abs(k(o,i)) / umwm%seadep(i) ! freq(k)=freq(-k), so k>0 == k<0 roots

     enddo

     iwsfc = i + omsea

     do o = 1, om

        ! phase/group velocity
        kd       = min( k(o,i) * umwm%seadep(i), 20. )  ! avoid overflow
        cp0(o,i) = pi2 * freq(o) / k(o,i)
        cg0(o,i) = cp0(o,i) * ( 0.5 + k(o,i) * umwm%seadep(i) / sinh(2. * kd) &
                              + sfct * k(o,i)**2 / (rhosw * grav + sfct * k(o,i)**2) )

        ! compute some frequently used arrays
        dwn       (o,i) = pi2 * dlnf * freq(o) / abs(cg0(o,i))           ! dk
        l2        (o,i) = 0.5 * abs(cp0(o,i)) / freq(o)                  ! lambda/2 (half wavelength)
        logl2overz(o,i) = log(min(20.,l2(o,i)) / sfcg%dzt_bot(iwsfc))
        k4        (o,i) = k(o,i)**4                                      ! k^4
        oneoverk4 (o,i) = 1. / k4(o,i)                                   ! k^-4
        kdk       (o,i) = k(o,i) * dwn(o,i)                              ! k*dk
        k3dk      (o,i) = k(o,i)**3 * dwn(o,i)                           ! k*k*k*dk
        fkovg     (o,i) = freq(o) * k(o,i) / grav                        ! freq*k/g
        cothkd    (o,i) = cosh(0.2 * kd) / sinh(0.2 * kd)                ! coth(0.2*kd)
        invcp0    (o,i) = 1. / cp0(o,i)                                  ! 1/cp
        sdv       (o,i) = sbf_fac * k(o,i) / sinh(2. * kd) &             ! bottom friction
                        + sbp_fac * k(o,i) / cosh(kd)**2 &               ! bottom percolation
                        + 4. * nu_water * k(o,i)**2                      ! viscosity
     enddo

     ! compute renormalization factors for snl:

     do o = om-1, om
        bf1_renorm(o,i) = 0.
        bf2_renorm(o,i) = 0.
     enddo

     do o = 1, om-2
!       bf1_renorm(o,i) = snl_fac * bf1 * kdk(o+1,i) / kdk(o,i)
!       bf2_renorm(o,i) = snl_fac * bf2 * kdk(o+2,i) / kdk(o,i)
!       snl_arg   (o,i) = 1. - (bf1_renorm(o,i) + bf2_renorm(o,i))

        bf1_renorm(o,i) = bf1 * kdk(o+1,i) / kdk(o,i)
        bf2_renorm(o,i) = bf2 * kdk(o+2,i) / kdk(o,i)
!       snl_arg   (o,i) = 1. - snl_fac * (bf1_renorm(o,i) + bf2_renorm(o,i))
     enddo

  enddo
  !$omp end parallel do

  write(io6,*) 'UMWM dispersion: dispersion relationship done;'
  write(io6,*)

end subroutine dispersion

end module umwm_init
