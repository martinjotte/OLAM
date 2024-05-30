module umwm_init

use umwm_module

implicit none

contains

!===============================================================================

subroutine nmlassign()

  use misc_coms, only: io6

  logical :: namelistok

  ! local variables for reading from namelist

  ! namelist /domain/

  !Bob - orig values were: om=37,pm=36,fmax=2.0,fprog=2.0
  !Bob - also, orig lm in stokes = 42; Bob changed lm to 1

  om        = 25      ! Number of frequency bins
  pm        = 32      ! Number of directions
  fmin      = 0.0313  ! Lowest frequency bin [Hz]
  fmax      = 0.5     ! Highest frequency bin [Hz]
  fprog     = 0.5     ! Highest prognostic frequency bin [Hz]

  ! namelist /physics/

  nu_air    = 1.56E-5 ! Kinematic viscosity of air [m^2/s]
  nu_water  = 0.90E-6 ! Kinematic viscosity of water [m^2/s]
  sfct      = 0.07    ! Surface tension [N/m]
  gustiness = 0.0     ! Random wind gustiness factor (should be between 0 and 0.2)
  dmin      = 10.     ! Ocean depth lower limit [m]
  explim    = 0.9     ! Exponent limiter (0.69 ~ 100% growth)
  sin_fac   = 0.11    ! Input factor from following winds
  sin_diss1 = 0.10    ! Damping factor from opposing winds
  sin_diss2 = 0.001   ! Damping factor from swell overrunning wind
  sds_fac   = 42.     ! Breaking dissipation factor
  sds_power = 2.4     ! Saturation spectrum power
  mss_fac   = 360     ! Mean-square-slope adjustment to Sds
  snl_fac   = 5.0     ! Wave energy downshifting factor
  sdt_fac   = 0.002   ! Dissipation due to turbulence factor
  sbf_fac   = 0.003   ! Bottom friction coefficient [m/s]
  sbp_fac   = 0.003   ! Bottom percolation coefficient [m/s]

  ! namelist /forcing_constant/

  fice_lth  = 0.30    ! Sea ice fraction - lower threshold for attenuation
  fice_uth  = 0.75    ! Sea ice fraction - upper threshold for attenuation

  ! namelist /output/

  stokes    = .true.  ! Compute Stokes drift velocity fields

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

  use consts_coms,   only: piu180, erad, pi2, pio4, grav
  use mem_sfcg,      only: sfcg, mvsfc, itab_vsfc
  use mem_sea,       only: sea, msea, omsea
  use umwm_oforcing, only: oforcing

  integer :: i, n, o, p, pp, ind, iw1, iw2, iwsfc
  real    :: raxis, dx, dy, thdirv, bfa, bfb

  ! set frequency increment:
  dlnf = (log(fmax) - log(fmin)) / real(om-1)

  ! set frequency bins:
  do o = 1,om
    umwm%f(o) = exp(log(fmin) + (o-1) * dlnf)
  enddo

  ! define various constants:
  dth           = pi2 / real(pm)
  dthg          = dth * grav
  oneovdth      = 1. / dth
  twopisds_fac  = pi2 * sds_fac
  fieldscale1   = sin_diss1 / sin_fac
  fieldscale2   = sin_diss2 / sin_diss1
  inv_sds_power = 1. / sds_power

  ! compute diffusion values in 2 frequencies:
  bfa = exp(-16. * dlnf * dlnf)
  bfb = exp(-64. * dlnf * dlnf)
  bf1 = bfa / (bfa + bfb)
  bf2 = bfb / (bfa + bfb)

  ! calculate wave ray directions
  do p = 1,pm
     th(p) = (p - 0.5 * (pm + 1)) * dth ! angles  {-175,-165,...,165,175 deg} but in radians
    cth(p) = cos(th(p)) ! cosines
    sth(p) = sin(th(p)) ! sines
  enddo

  do i = 2,msea
    iwsfc = i + omsea
    umwm%seadep(i) = max(dmin, sfcg%topw(iwsfc) - sfcg%bathym(iwsfc))
  enddo

  do i = 2,mvsfc
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

  ! "left" and "right" directional indices for refraction:
  do p = 1,pm
    pl(p) = p + 1
    pr(p) = p - 1
  enddo
  pl(pm) = 1
  pr(1)  = pm

  do p = 1,pm
    cth2(p) = cos(dth * (p-1))**2 ! same as cos(th)**2?
  enddo

  allocate(cth2pp(pm,pm))

  do p=1,pm
    do pp=1,pm
      ind = pp - p + 1
      if (ind <= 0) ind = pm + ind
      cth2pp(pp,p) = cth2(ind) * dth
    enddo
  enddo

  call oforcing()

  ! compute wave numbers, phase speeds, and group velocities:
  call dispersion(1e-2)

  ! initialize drag coefficient (Large and Pond, 1981):
  umwm%cd = 1.2e-3
  do concurrent (i=2:msea, umwm%wspd(i) > 11)
    umwm%cd(i) = (0.49 + 0.065 * umwm%wspd(i)) * 1e-3
  enddo

  ! initialize friction velocity:
  do i = 2,msea
    umwm%ustar(i) = sqrt(umwm%cd(i)) * umwm%wspd(i)
  enddo

  write(*,fmt=102)
  write(*,*)'initialization summary:'
  write(*,fmt=103)
  write(*,*)'bin  f[hz]    t[s]    min(k)   max(k)   min(c)   max(c)   min(cg)   max(cg)'
  write(*,fmt=103)
  do o=1,om
    write(*,fmt=101)o,umwm%f(o),1./umwm%f(o),                                          &
                    minval(k  (o,2:msea),dim=1),maxval(k  (o,2:msea),dim=1), &
                    minval(cp0(o,2:msea),dim=1),maxval(cp0(o,2:msea),dim=1), &
                    minval(cg0(o,2:msea),dim=1),maxval(cg0(o,2:msea),dim=1)
  enddo
  101 format(1x,i2,8(2x,f7.4))
  write(*,fmt=102)

  102 format('!',77('='),'!')
  103 format('!',77('-'),'!')

!Bob: print out initialized (om,msea) arrays for a single isea point

i = 91457 - omsea
print*, ' '
write(6,'(180a)') '               o bf1_renorm bf2_renorm cp0     cg0     cothkd     dwn    fkovg     l2   logl2overz  snl_arg     k       k4          kdk         k3dk        sbf         sdv'
print*, ' '
do o = 1,om
   write(6,'(a,i6,11f9.3,9e12.3)') 'init(o,i) ', o, &
                     bf1_renorm(o,i), bf2_renorm(o,i),   cp0(o,i), cg0(o,i), &
                         cothkd(o,i),        dwn(o,i), fkovg(o,i),  l2(o,i), &
                     logl2overz(o,i),    snl_arg(o,i),     k(o,i),  k4(o,i), &
                            kdk(o,i),       k3dk(o,i),   sbf(o,i), sdv(o,i)
enddo
print*, ' '

end subroutine init

!===============================================================================

subroutine dispersion(tol)

  use consts_coms, only: pi2, grav
  use mem_sea,     only: msea, omsea
  use mem_sfcg,    only: sfcg

  ! Iteratively solve the dispersion relation by iteration,
  ! and compute phase and group velocities in absolute reference frame (w/ currents)

  integer :: src, dest

  integer :: counter, i, o, iwsfc
  real,intent(in) :: tol
  real :: dk

  real :: dom (om)
  real :: b   (msea)
  real :: f_nd(om,msea)
  real :: kd  (om,msea)
  real :: t   (om,msea)

  write(*,'(a)')'umwm: dispersion: solving for dispersion relationship;'

  ! non-dimesionalize frequencies, and use deep water limit
  ! as initial guess:
  do concurrent (o=1:om, i=2:msea)
     cp0(o,i)  = pi2 * sqrt(umwm%seadep(i) / grav)
     f_nd(o,i) = cp0(o,i) * umwm%f(o)
     k(o,i)    = f_nd(o,i) * f_nd(o,i)
  enddo

  ! non-dimesionalize surface tension:
  b(2:msea) = sfct / (rhosw * grav * umwm%seadep(2:msea)**2)

  do i = 2,msea
     do o = 1,om

        counter = 1
        dk = 2. * tol

        do while (abs(dk) > tol) ! newton-raphson iteration loop

           t(o,i) = tanh(k(o,i))

           dk = -(f_nd(o,i) * f_nd(o,i) - k(o,i) * t(o,i)        &
                * (1. + b(i) * k(o,i) * k(o,i)))                 &
                / (3. * b(i) * k(o,i) * k(o,i) * t(o,i) + t(o,i) &
                + k(o,i) * (1. + b(i) * k(o,i) * k(o,i)) * (1. - t(o,i) * t(o,i)))
           k(o,i) = k(o,i) - dk

           if (counter == 1000) exit ! escape if stuck
           counter = counter + 1

        enddo

        k(o,i) = abs(k(o,i)) / umwm%seadep(i) ! umwm%f(k)=umwm%f(-k), so k>0 == k<0 roots

     enddo
  enddo

  write(*,'(a)')'umwm: dispersion: dispersion relationship done;'


  do concurrent (o=1:om, i=2:msea)
     kd(o,i) = k(o,i) * umwm%seadep(i)
  enddo

  ! limit kd to avoid floating overflow in transcendental functions:
  where (kd > 20.) kd = 20.

  ! phase speed and group velocity:
  cp0 = tiny(cp0)
  cg0 = tiny(cg0)

  do concurrent (o=1:om, i=2:msea)

     cp0(o,i) = pi2 * umwm%f(o) / k(o,i)
     cg0(o,i) = cp0(o,i) * (0.5 + k(o,i) * umwm%seadep(i) / sinh(2. * kd(o,i)) &
                         + sfct * k(o,i) * k(o,i) / (rhosw * grav + sfct * k(o,i) * k(o,i)))
  enddo

print*, 'grav,sfct,rhosw ', grav, sfct, rhosw

print*, ' '
write(6,'(140a)') '           o    cp0     cg0    l2     umwm%f    kd     cothkd   dwn       k        kdk         k4         sbf         sdv        fkovg'
print*, ' '

  ! compute some frequently used arrays:
  do i = 2,msea
     iwsfc = i + omsea
     do o = 1,om

        dwn       (o,i) = pi2 * dlnf * umwm%f(o) / abs(cg0(o,i))             ! dk
        l2        (o,i) = 0.5 * abs(cp0(o,i)) / umwm%f(o)                    ! lambda/2 (half wavelength)
        logl2overz(o,i) = log(min(20.,l2(o,i)) / sfcg%dzt_bot(iwsfc))
        k4        (o,i) = k(o,i)**4                                          ! k^4
        oneoverk4 (o,i) = 1. / k4(o,i)                                       ! k^-4
        kdk       (o,i) = k(o,i) * dwn(o,i)                                  ! k*dk
        k3dk      (o,i) = k(o,i)**3 * dwn(o,i)                               ! k*k*k*dk
        fkovg     (o,i) = umwm%f(o) * k(o,i) / grav                          ! umwm%f*k/g
        cothkd    (o,i) = cosh(0.2 * kd(o,i)) / sinh(0.2 * kd(o,i))          ! coth(0.2*kd)
        invcp0    (o,i) = 1. / cp0(o,i)                                      ! 1/cp
        sbf       (o,i) = sbf_fac * k(o,i) / (sinh(2. * kd(o,i))) &          ! bottom friction
                        + sbp_fac * k(o,i) / (cosh(kd(o,i)) * cosh(kd(o,i))) ! bottom percolation
        sdv       (o,i) = 4. * nu_water * k(o,i)**2                          ! viscosity

if (i == 1000) then
   write(6,'(a,i6,6f8.2,2f8.3,9e12.3)') 'init7 ',o,cp0(o,i),cg0(o,i),l2(o,i),umwm%f(o),kd(o,i), &
        cothkd(o,i),dwn(o,i),k(o,i),kdk(o,i),k4(o,i),sbf(o,i),sdv(o,i),fkovg(o,i)
endif

     enddo
  enddo

  ! compute renormalization factors for snl:
  bf1_renorm = 0.
  bf2_renorm = 0.
  snl_arg    = 0.

  do i = 2,msea
     do o = 1,om-2
        bf1_renorm(o,i) = snl_fac * bf1 * kdk(o+1,i) / kdk(o,i)
        bf2_renorm(o,i) = snl_fac * bf2 * kdk(o+2,i) / kdk(o,i)
        snl_arg   (o,i) = 1. - (bf1_renorm(o,i) + bf2_renorm(o,i))
     enddo
  enddo

  ! handle land points for cp and cg (needed for advection/refraction)
  do o = 1,om
     cp0(o,1) = minval(cp0(o,2:msea))
     cg0(o,1) = minval(cg0(o,2:msea))
  enddo

end subroutine dispersion

end module umwm_init
