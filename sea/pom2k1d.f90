Module pom2k1d

  implicit none

  integer :: nzpom

  real :: dtl2, pom_dztop, pom_depth

  integer, parameter :: ntp = 2  ! Jerlov water type (1=i,2=ia,3=ib,4=ii,5=iii) used in proft

  real, parameter :: cbcmin = .0025e0 ! Minimum bottom friction coefficient
  real, parameter :: cbcmax = 1.e0    ! Maximum bottom friction coefficient
  real, parameter :: umol   = 2.e-5   ! Background vertical diffusivity
  real, parameter :: rhoref = 1025.e0 ! Reference density (seawater) [kg/m^3]
  real, parameter :: small  = 1.e-9   ! Small value []
  real, parameter :: z0b    = .01e0   ! Bottom roughness [m]
  real, parameter :: smoth  = 0.10e0  ! Constant in temporal filter (to prevent splitting)
  real, parameter :: grav   = 9.806e0 ! gravity constant [m/s^2]
  real, parameter :: vonk   = 0.4     ! von Karman constant []

  real, allocatable :: y(:), yy(:), dy(:), dyy(:), cbc(:) ! cbc = bottom friction coefficient

  Type pom_vars
     real, allocatable :: potmp  (:,:) ! current potential temp     [deg C]
     real, allocatable :: salin  (:,:) ! current salinity           [psu = g/kg]
     real, allocatable :: q2     (:,:) ! current tke*2              [m^2/s^2]
     real, allocatable :: q2l    (:,:) ! current q2*tlen            [m^3/s^2]
     real, allocatable :: u      (:,:) ! current eastward velocity  [m/s]
     real, allocatable :: v      (:,:) ! current northward velocity [m/s]

     real, allocatable :: potmpb (:,:) ! past potential temp     [deg C]
     real, allocatable :: salinb (:,:) ! past salinity           [psu or g/kg]
     real, allocatable :: q2b    (:,:) ! past tke*2              [m^2/s^2]
     real, allocatable :: q2lb   (:,:) ! past q2*tlen            [m^3/s^2]
     real, allocatable :: ub     (:,:) ! past eastward velocity  [m/s]
     real, allocatable :: vb     (:,:) ! past northward velocity [m/s]

     real, allocatable :: km     (:,:) ! vertical kinematic viscosity [m^2/s]
     real, allocatable :: kh     (:,:) ! vertical diffusivity         [m^2/s]
     real, allocatable :: kq     (:,:) ! vertical q diffusivity       [m^2/s]

     real, allocatable :: wubot    (:) ! bottom U momentum flux      [m^2/s^2]
     real, allocatable :: wvbot    (:) ! bottom V momentum flux      [m^2/s^2]
     real, allocatable :: cor      (:) ! Coriolis parameter          [1/s]

     integer, allocatable :: kba   (:) ! number of levels in each column []
  End type

  type (pom_vars) :: pom

Contains

!===============================================================================

  subroutine alloc_pomgrid(mpom)

  implicit none

  integer, intent(in) :: mpom

  allocate (y  (nzpom)); y   = 0.
  allocate (yy (nzpom)); yy  = 0.
  allocate (dy (nzpom)); dy  = 0.
  allocate (dyy(nzpom)); dyy = 0.

  allocate (pom%kba (mpom)); pom%kba  = 0

  end subroutine alloc_pomgrid

!===============================================================================

  subroutine pom_levels()

  implicit none

  integer :: iter, k
  real :: srati, thick

  ! Iteration to find POM1D vertical stretch rate

  write(*,'(/,a,/)') 'Iterations for POM1D vertical stretch ratio'

  srati = 0.5
  do iter = 1,20
     srati = 1. / (pom_depth * (1. - srati) / (srati * pom_dztop) + 1.)**(1./real(nzpom-1))
     print*, 'iter, stretch ratio ',iter,1./srati
  enddo

  ! Compute POM1D grid levels

  thick = pom_dztop

  y(1) = 0.

  do k = 2, nzpom
     y (k)   = y(k-1) - thick
     yy(k-1) = 0.5 * (y(k) + y(k-1))

     thick = thick / srati
  enddo

  yy(nzpom) = 2.e0 * yy(nzpom-1) - yy(nzpom-2)

  do k = 1, nzpom-1
     dy (k) = y (k) - y (k+1)
     dyy(k) = yy(k) - yy(k+1)
  enddo

  dy (nzpom) = 0.
  dyy(nzpom) = 0.

  write(6,1)
  1 format(/3x,'k',7x,'y',9x,'yy',9x,'dy',9x,'dyy',/)

  do k = 1,nzpom
     write(6,2) k,y(k),yy(k),dy(k),dyy(k)
     2 format((' ',i3,4f11.2))
  enddo

  write(6,3)
  3 format(//)

  end subroutine pom_levels

!===============================================================================

  subroutine alloc_pom(mpom)

  implicit none

  integer, intent(in) :: mpom

  allocate (pom%potmp (nzpom,mpom)); pom%potmp  = 0.
  allocate (pom%salin (nzpom,mpom)); pom%salin  = 0.
  allocate (pom%q2    (nzpom,mpom)); pom%q2     = 0.
  allocate (pom%q2l   (nzpom,mpom)); pom%q2l    = 0.
  allocate (pom%u     (nzpom,mpom)); pom%u      = 0.
  allocate (pom%v     (nzpom,mpom)); pom%v      = 0.

  allocate (pom%potmpb(nzpom,mpom)); pom%potmpb = 0.
  allocate (pom%salinb(nzpom,mpom)); pom%salinb = 0.
  allocate (pom%q2b   (nzpom,mpom)); pom%q2b    = 0.
  allocate (pom%q2lb  (nzpom,mpom)); pom%q2lb   = 0.
  allocate (pom%ub    (nzpom,mpom)); pom%ub     = 0.
  allocate (pom%vb    (nzpom,mpom)); pom%vb     = 0.

  allocate (pom%km    (nzpom,mpom)); pom%km     = 0.
  allocate (pom%kh    (nzpom,mpom)); pom%kh     = 0.
  allocate (pom%kq    (nzpom,mpom)); pom%kq     = 0.

  allocate (pom%wubot       (mpom)); pom%wubot  = 0.
  allocate (pom%wvbot       (mpom)); pom%wvbot  = 0.
  allocate (pom%cor         (mpom)); pom%cor    = 0.

  allocate (cbc(nzpom)); cbc = 0.

  end subroutine alloc_pom

!===============================================================================

  subroutine filltab_pom()

  use var_tables, only: increment_vtable

  implicit none

  if (allocated(pom%potmp))  call increment_vtable('POM%POTMP' ,'SW', rvar2=pom%potmp)
  if (allocated(pom%salin))  call increment_vtable('POM%SALIN' ,'SW', rvar2=pom%salin)
  if (allocated(pom%q2))     call increment_vtable('POM%Q2'    ,'SW', rvar2=pom%q2)
  if (allocated(pom%q2l))    call increment_vtable('POM%Q2L'   ,'SW', rvar2=pom%q2l)
  if (allocated(pom%u))      call increment_vtable('POM%U'     ,'SW', rvar2=pom%u)
  if (allocated(pom%v))      call increment_vtable('POM%V'     ,'SW', rvar2=pom%v)

  if (allocated(pom%potmpb)) call increment_vtable('POM%POTMPB','SW', rvar2=pom%potmpb)
  if (allocated(pom%salinb)) call increment_vtable('POM%SALINB','SW', rvar2=pom%salinb)
  if (allocated(pom%q2b))    call increment_vtable('POM%Q2B'   ,'SW', rvar2=pom%q2b)
  if (allocated(pom%q2lb))   call increment_vtable('POM%Q2LB'  ,'SW', rvar2=pom%q2lb)
  if (allocated(pom%ub))     call increment_vtable('POM%UB'    ,'SW', rvar2=pom%ub)
  if (allocated(pom%vb))     call increment_vtable('POM%VB'    ,'SW', rvar2=pom%vb)

  if (allocated(pom%km))     call increment_vtable('POM%KM'    ,'SW', rvar2=pom%km)
  if (allocated(pom%kh))     call increment_vtable('POM%KH'    ,'SW', rvar2=pom%kh)
  if (allocated(pom%kq))     call increment_vtable('POM%KQ'    ,'SW', rvar2=pom%kq)

  if (allocated(pom%wubot))  call increment_vtable('POM%WUBOT' ,'SW', rvar1=pom%wubot)
  if (allocated(pom%wvbot))  call increment_vtable('POM%WVBOT' ,'SW', rvar1=pom%wvbot)

  end subroutine filltab_pom

!===============================================================================

  subroutine pom_column(ipom, kb, wusurf, wvsurf, wtsurf, wssurf, swrad)

  implicit none

  integer, intent(in) :: ipom, kb
  integer :: k

  real, intent(in) :: wusurf  ! surface U wind stress  [m^2/s^2]
  real, intent(in) :: wvsurf  ! surface V wind stress  [m^2/s^2]
  real, intent(in) :: wtsurf  ! surface sens heat flux [m/s K]  (includes net LW rad)
  real, intent(in) :: wssurf  ! surface salinity flux  [m/s psu]
  real, intent(in) :: swrad   ! net SW rad flux        [W/m^2]

  real :: potmpf(nzpom) ! future potential temperature [deg C]
  real :: salinf(nzpom) ! future salinity              [psu or g/kg]
  real :: q2f   (nzpom) ! future tke*2                 [m^2/s^2]
  real :: q2lf  (nzpom) ! future q2*tlen               [m^3/s^2]
  real :: uf    (nzpom) ! future eastward velocity     [m/s]
  real :: vf    (nzpom) ! future northward velocity    [m/s]
  real :: rho   (nzpom) ! water density                [kg/m^3]

  ! RETURN if kb < 4 because (1) not all lines of code can be executed with
  ! fewer levels and (2) the ocean depth where kb < 4 is shallow enough that
  ! the layers are likely well mixed anyway.

  if (kb < 4) RETURN

  call dens(kb, pom%potmpb(:,ipom), pom%salinb(:,ipom), rho)

  do k = 1,kb
     q2f (k) = pom%q2b (k,ipom) ! Remnant of subroutine advq called here
     q2lf(k) = pom%q2lb(k,ipom) ! Remnant of subroutine advq called here
  enddo

  call profq(kb, ipom, wusurf, wvsurf, rho, q2f, q2lf)

  ! Step q2 and q2l forward in time

  do k = 1,kb
     pom%q2 (k,ipom) = pom%q2 (k,ipom) + .5e0 * smoth &
                     * (q2f (k) + pom%q2b (k,ipom) - 2.e0 * pom%q2 (k,ipom))
     pom%q2l(k,ipom) = pom%q2l(k,ipom) + .5e0 * smoth &
                     * (q2lf(k) + pom%q2lb(k,ipom) - 2.e0 * pom%q2l(k,ipom))

     pom%q2b (k,ipom) = pom%q2 (k,ipom)
     pom%q2lb(k,ipom) = pom%q2l(k,ipom)
     pom%q2  (k,ipom) = q2f    (k)
     pom%q2l (k,ipom) = q2lf   (k)
  enddo

  pom%potmp (kb,ipom) = pom%potmp (kb-1,ipom) ! Remnant of subroutine advct1 called here
  pom%salin (kb,ipom) = pom%salin (kb-1,ipom) ! Remnant of subroutine advct1 called here
  pom%potmpb(kb,ipom) = pom%potmpb(kb-1,ipom) ! Remnant of subroutine advct1 called here
  pom%salinb(kb,ipom) = pom%salinb(kb-1,ipom) ! Remnant of subroutine advct1 called here

  do k = 1,kb-1
     potmpf(k) = pom%potmpb(k,ipom) ! Remnant of subroutine advct1 called here
     salinf(k) = pom%salinb(k,ipom) ! Remnant of subroutine advct1 called here
  enddo

  call proft(kb, ipom, potmpf, wtsurf, swrad, 2)
  call proft(kb, ipom, salinf, wssurf, swrad, 1)

  ! Step potmp and salin forward in time

  do k = 1,kb
     pom%potmp(k,ipom) = pom%potmp(k,ipom) + .5e0 * smoth &
                       * (potmpf(k) + pom%potmpb(k,ipom) - 2.e0 * pom%potmp(k,ipom))
     pom%salin(k,ipom) = pom%salin(k,ipom) + .5e0 * smoth &
                       * (salinf(k) + pom%salinb(k,ipom) - 2.e0 * pom%salin(k,ipom))

     pom%potmpb(k,ipom) = pom%potmp(k,ipom)
     pom%salinb(k,ipom) = pom%salin(k,ipom)
     pom%potmp (k,ipom) = potmpf   (k)
     pom%salin (k,ipom) = salinf   (k)
  enddo

  call dens(kb, pom%potmp(:,ipom), pom%salin(:,ipom), rho)

  ! Coriolis force (remnant of subroutines advu and advv)

  do k = 1,kb-1
     uf(k) = pom%ub(k,ipom) + dtl2 * pom%cor(ipom) * pom%v(k,ipom)
     vf(k) = pom%vb(k,ipom) - dtl2 * pom%cor(ipom) * pom%u(k,ipom)
  enddo

  call profuv(kb, ipom, uf, vf, wusurf, wvsurf)

  ! Step u and v forward in time

  do k = 1,kb-1
     pom%u(k,ipom) = pom%u(k,ipom) + .5e0 * smoth &
                   * (uf(k) + pom%ub(k,ipom) - 2.e0 * pom%u(k,ipom))

     pom%v(k,ipom) = pom%v(k,ipom) + .5e0 * smoth &
                   * (vf(k) + pom%vb(k,ipom) - 2.e0 * pom%v(k,ipom))
  enddo

  do k = 1,kb
     pom%ub(k,ipom) = pom%u(k,ipom)
     pom%vb(k,ipom) = pom%v(k,ipom)
     pom%u (k,ipom) = uf   (k)
     pom%v (k,ipom) = vf   (k)
  enddo

  end subroutine pom_column

!===============================================================================

  subroutine dens(kb, potmp, salin, rhoo)

  ! Calculates (density-1000.)/rhoref.
  !    (Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611.)

  ! If using 32 bit precision, it is recommended that cr, p, rhor, sr, tr,
  ! tr2, tr3, and tr4 be made double precision, and the "e"s in the
  ! constants be changed to "d"s.

  implicit none

  integer, intent(in) :: kb

  real, intent(in)  :: potmp(kb), salin(kb)
  real, intent(out) :: rhoo(kb)

  integer :: k
  real :: cr, p, rhor, sr, tr, tr2, tr3, tr4, tr5

  do k = 1,kb-1

     tr = potmp(k)
     sr = salin(k)
     tr2 = tr * tr
     tr3 = tr2 * tr
     tr4 = tr3 * tr
     tr5 = tr4 * tr

     ! Approximate pressure in units of bars:

     p = grav * rhoref * (-yy(k)) * 1.e-5

     rhor = -0.157406e0 + 6.793952e-2 * tr - 9.095290e-3 * tr2              &
          +  1.001685e-4 * tr3 - 1.120083e-6 * tr4 + 6.536332e-9 * tr5

     rhor = rhor + (0.824493e0 - 4.0899e-3 * tr + 7.6438e-5 * tr2           &
          - 8.2467e-7 * tr3 + 5.3875e-9 * tr4) * sr                         &
          + (-5.72466e-3 + 1.0227e-4 * tr - 1.6546e-6 * tr2) * abs(sr)**1.5 &
          + 4.8314e-4 * sr * sr

     cr = 1449.1e0 + .0821e0 * p + 4.55e0 * tr - .045e0 * tr2 &
        + 1.34e0 * (sr - 35.e0)

     rhor = rhor + 1.e5 * p / (cr * cr) * (1.e0 - 2.e0 * p / (cr * cr))

     rhoo(k) = rhor / rhoref

  enddo

  end subroutine dens

!===============================================================================

  subroutine profq(kb, ipom, wusurf, wvsurf, rho, q2f, q2lf)

  ! Updated: Sep. 24, 2003
  ! Solves for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
  ! length scale), km (vertical kinematic viscosity), and kh (vertical
  ! kimatic diffusivity), using a simplified version of the level 2 1/2 model
  ! of Mellor and Yamada (1982).  In this version, the Craig-Banner sub-model
  ! whereby breaking wave tke is injected into the surface is included.
  ! However, we use an analytical solution to the near surface tke equation to
  ! solve for q2 at the surface giving the same result as C-B diffusion. The
  ! new scheme is simpler and more robust than the latter scheme.

  ! Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence in the
  !     ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, 1994.
  ! Ezer, T., On the seasonal mixed-layer simulated by a basin-scale ocean
  !     ocean model and the Mellor-Yamada turbulence scheme,
  !     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.
  ! Mellor, G.L. and T. Yamada, Development of a turbulence closure model for
  !     geophysical fluid fluid problems, Rev. Geophys. Space Phys., 20,
  !     851-875, 1982.
  ! Mellor, G. L., One-dimensional, ocean surface layer modeling, a problem
  !     and a solution. J. Phys. Oceanogr., 31(3), 790-809, 2001.
  ! Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface thermal
  !     response, J. Phys. Oceanogr., 2003.
  ! Stacey, M. W., Simulations of the wind-forced near-surface circulation in
  !     Knight Inlet: a parameterization of the roughness length.
  !     J. Phys. Oceanogr., 29, 1363-1367, 1999.

  implicit none

  integer, intent(in) :: kb, ipom
  real, intent(in)    :: wusurf, wvsurf
  real, intent(in)    :: rho(kb)
  real, intent(inout) :: q2f(kb), q2lf(kb)

  real :: a(kb), c(kb), ee(kb), gg(kb), cc(kb), boygr(kb), l(kb), gh(kb)
  real :: prod(kb), stf(kb), dtef(kb), sh(kb), sm(kb), kn(kb)
  real :: coef1, coef2, coef3, coef4, coef5
  real :: ghc
  real :: tp, sp, p
  real :: utau2, l0

  integer k,ki

  real, parameter :: a1=0.92e0, b1=16.6e0, a2=0.74e0, b2=10.1e0, c1=0.08e0
  real, parameter :: e1=1.8e0, e2=1.33e0
  real, parameter :: sef=1.e0, const1=(16.6e0**(2.e0/3.e0))
  real, parameter :: cbcnst=100., surfl=2.e5, shiw=0.0

  do k = 2,kb-1
     a(k) = -dtl2 * ((pom%kq(k+1,ipom) + pom%kq(k,ipom)) * .5e0 + umol) / (dyy(k-1) * dy(k  ))
     c(k) = -dtl2 * ((pom%kq(k-1,ipom) + pom%kq(k,ipom)) * .5e0 + umol) / (dyy(k-1) * dy(k-1))
  enddo

  ! The following section solves the equation:  dt2*(kq*q2')' - q2*(2.*dt2*dtef+1.) = -q2b

  ! Surface and bottom boundary conditions:

  utau2 = sqrt(wusurf**2 + wvsurf**2)
  ! Wave breaking energy- a variant of Craig & Banner (1994); see Mellor and Blumberg, 2003.
  ee(1) = 0.e0
  gg(1) = (15.8 * cbcnst)**(2./3.) * utau2
  ! Surface length scale following Stacey (1999).
  l0 = surfl * utau2 / grav

  q2f(kb) = sqrt(pom%wubot(ipom)**2 + pom%wvbot(ipom)**2) * const1

  ! Calculate speed of sound squared:

  do k = 1,kb-1
     tp = pom%potmp(k,ipom)
     sp = pom%salin(k,ipom)
     p = grav * rhoref * (-yy(k)) * 1.e-4  ! pressure in units of decibars:
     cc(k) = 1449.1e0 + .00821e0 * p + 4.55e0 * tp - .045e0 * tp**2 + 1.34e0 * (sp - 35.0e0)
     cc(k) = cc(k) / sqrt((1.e0 - .01642e0 * p / cc(k)) * (1.e0 - 0.40e0 * p / cc(k)**2))
  enddo

  ! Calculate buoyancy gradient:

  do k = 2,kb-1
     pom%q2b (k,ipom) = abs(pom%q2b (k,ipom))
     pom%q2lb(k,ipom) = abs(pom%q2lb(k,ipom))
     boygr(k) = grav * (rho(k-1) - rho(k)) / dyy(k-1) &
                   + grav**2 * 2.e0 / (cc(k-1)**2 + cc(k)**2)
  enddo

  do k = 2,kb-1
     l(k) = abs(pom%q2lb(k,ipom) / pom%q2b(k,ipom))
     if (y(k) > 0.5 * y(kb)) l(k) = max(l(k), vonk * l0)
     gh(k) = (l(k)**2) * boygr(k) / pom%q2b(k,ipom)
     gh(k) = min(gh(k), .028e0)
  enddo

  l(1) = vonk * l0
  l(kb)  = 0.e0
  gh(1)  = 0.e0
  gh(kb) = 0.e0

  ! Calculate production of turbulent kinetic energy:

  do k = 2,kb-1
     prod(k) = pom%km(k,ipom) * sef * ((pom%u(k,ipom) - pom%u(k-1,ipom))**2  &
                                     + (pom%v(k,ipom) - pom%v(k-1,ipom))**2) &
             / dyy(k-1)**2 &
     ! Add shear due to internal wave field
             - shiw * pom%km(k,ipom) * boygr(k)
     prod(k) = prod(k) + pom%kh(k,ipom) * boygr(k)
  enddo

  ! NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000),
  ! depends on ghc the critical number (empirical -6 to -2) to increase mixing.

  ghc = -6.0e0
  do k = 1,kb
     stf(k) = 1.e0
     ! It is unclear yet if diss. corr. is needed when surf. waves are included.
     !           if(gh(i,j,k).lt.0.e0)
     !    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
     !           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
     dtef(k) = sqrt(abs(pom%q2b(k,ipom))) * stf(k) / (b1 * l(k) + small)
  enddo

  do k = 2,kb-1
     gg(k) = 1.e0 / (a(k) + c(k) * (1.e0 - ee(k-1)) - (2.e0 * dtl2 * dtef(k) + 1.e0))
     ee(k) = a(k) * gg(k)
     gg(k) = (-2.e0 * dtl2 * prod(k) + c(k) * gg(k-1) - q2f(k)) * gg(k)
  enddo

  do k = 1,kb-1
     ki = kb - k
     q2f(ki) = ee(ki) * q2f(ki+1) + gg(ki)
  enddo

  ! The following section solves the equation:  dt2(kq*q2l')' - q2l*(dt2*dtef+1.) = -q2lb

  ee(2)  = 0.e0
  gg(2)  = 0.e0
  q2lf(kb) = 0.e0

  do k = 2,kb-1
     dtef(k) = dtef(k) * (1.e0 + e2 * ((1.e0 / abs(y(k) - y(1)) &
             + 1.e0 / abs(y(k) - y(kb))) * l(k) / vonk)**2)
     gg(k) = 1.e0 / (a(k) + c(k) * (1.e0 - ee(k-1)) - (dtl2 * dtef(k) + 1.e0))
     ee(k) = a(k) * gg(k)
     gg(k) = (dtl2 * (-prod(k) * l(k) * e1) + c(k) * gg(k-1) - q2lf(k)) * gg(k)
  enddo

  do k = 1,kb-2
     ki = kb - k
     q2lf(ki) = ee(ki) * q2lf(ki+1) + gg(ki)
  enddo

  do k = 2,kb-1
     if (q2f(k) <= small .or. q2lf(k) <= small) then
        q2f(k) = small
        q2lf(k) = 0.1 * y(kb) * small
     endif
  enddo

  ! The following section solves for km and kh:

  coef4 = 18.e0 * a1 * a1 + 9.e0 * a1 * a2
  coef5 = 9.e0 * a1 * a2

  ! Note that sm and sh limit to infinity when gh approaches 0.0288:

  do k = 1,kb
     coef1 = a2 * (1.e0 - 6.e0 * a1 / b1 * stf(k))
     coef2 = 3.e0 * a2 * b2 / stf(k) + 18.e0 * a1 * a2
     coef3 = a1 * (1.e0 - 3.e0 * c1 - 6.e0 * a1 / b1 * stf(k))
     sh(k) = coef1 / (1.e0 - coef2 * gh(k))
     sm(k) = coef3 + sh(k) * coef4 * gh(k)
     sm(k) = sm(k) / (1.e0 - coef5 * gh(k))
  enddo

  do k = 1,kb
     kn(k) = l(k) * sqrt(abs(pom%q2(k,ipom)))
     pom%kq(k,ipom) = (kn(k) * .41e0 * sh(k) + pom%kq(k,ipom)) * .5e0
     pom%km(k,ipom) = (kn(k)         * sm(k) + pom%km(k,ipom)) * .5e0
     pom%kh(k,ipom) = (kn(k)         * sh(k) + pom%kh(k,ipom)) * .5e0
  enddo

  end subroutine profq

!===============================================================================

  subroutine proft(kb, ipom, fut, wfsurf, swrad, nbc)

  ! Solves for vertical diffusion of x-momentum using method described by
  ! Richmeyer and Morton.  Irradiance parameters (Jerlov water type)
  ! are from Paulson and Simpson.

  ! Richtmeyer R.D., and K.W. Morton, 1967. Difference Methods for
  ! Initial-Value Problems, 2nd edition, Interscience, New York, 198-201.

  ! Paulson, C. A., and J. Simpson, 1977: Irradiance measurements in the
  ! upper ocean, J. Phys. Oceanogr., 7, 952-956.

  ! NOTES: (1) wfsurf and swrad are negative values when water column is
  !            warming or salt is being added.
  !        (2) nbc may only be 1 and 3 for salinity.

  implicit none

  integer, intent(in) :: kb, ipom, nbc  ! nbc=1 for salin; nbc=2 for potmp

  real, intent(inout) :: fut(kb) ! Future value of pot temp or salinity
  real, intent(in)    :: wfsurf  ! Surface flux of pot temp or salinity
  real, intent(in)    :: swrad   ! Surface SW rad flux

  real :: a(kb), c(kb), rad(kb), ee(kb), gg(kb)
  integer :: k, ki

  ! Jerlov water type:           i       ia      ib      ii     iii

  real, parameter :: r  (5) = [.58e0,  .62e0,  .67e0,  .77e0,  .78e0]
  real, parameter :: ad1(5) = [.35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0]
  real, parameter :: ad2(5) = [23.e0,  20.e0,  17.e0,  14.e0,  7.9e0]

  ! The following section solves the equation:  dt2*(kh*f')'-f=-fb

  do k = 2,kb-1
     a(k-1) = -dtl2 * (pom%kh(k,ipom) + umol) / (dy(k-1) * dyy(k-1))
     c(k)   = -dtl2 * (pom%kh(k,ipom) + umol) / (dy(k)   * dyy(k-1))
  enddo

  ! Calculate penetrative radiation. At the bottom any unattenuated
  ! radiation is deposited in the bottom layer:

  do k = 1,kb
     rad(k) = 0.e0
  enddo

  if (nbc == 2) then
     do k = 1,kb-1
        rad(k) = swrad * (r(ntp) * exp(y(k) / ad1(ntp)) &
               + (1.e0 - r(ntp)) * exp(y(k) / ad2(ntp)))
     enddo
  endif

  ee(1) = a(1) / (a(1) - 1.e0)
  gg(1) = dtl2 * (wfsurf + rad(1) - rad(2)) / dy(1) - fut(1)
  gg(1) = gg(1) / (a(1) - 1.e0)

  do k = 2,kb-2
     gg(k) = 1.e0 / (a(k) + c(k) * (1.e0 - ee(k-1)) - 1.e0)
     ee(k) = a(k) * gg(k)
     gg(k) = (c(k) * gg(k-1) - fut(k) + dtl2 * (rad(k) - rad(k+1)) / dy(k)) * gg(k)
  enddo

  ! Bottom adiabatic boundary condition:

  fut(kb-1) = (c(kb-1) * gg(kb-2) - fut(kb-1) + dtl2 * (rad(kb-1) - rad(kb)) / dy(kb-1)) &
            / (c(kb-1) * (1.e0 - ee(kb-2)) - 1.e0)

  do k = 2,kb-1
     ki = kb - k
     fut(ki) = (ee(ki) * fut(ki+1) + gg(ki))
  enddo

  end subroutine proft

!===============================================================================

  subroutine profuv(kb, ipom, uf, vf, wusurf, wvsurf)

  ! Solves for vertical diffusion of u- and v-momentum using method described by:
  ! Richtmeyer R.D., and K.W. Morton, 1967. Difference  Methods for
  ! Initial-Value Problems, 2nd edition, Interscience, New York, 198-201.

  ! NOTE that wusurf and wvsurf have the opposite sign to the wind velocity.

  implicit none

  integer, intent(in) :: kb, ipom

  real, intent(inout) :: uf(kb), vf(kb)
  real, intent(in)    :: wusurf, wvsurf

  real :: aa, bb, ff
  real :: a(kb), c(kb), ee(kb), gg(kb), hh(kb)

  integer :: k, ki

  ! The following section solves the equation:  dt2*(km*u')'-u=-ub

  do k = 2,kb-1
     a(k-1) = -dtl2 * (pom%km(k,ipom) + umol) / (dy(k-1) * dyy(k-1))
     c(k)   = -dtl2 * (pom%km(k,ipom) + umol) / (dy(k)   * dyy(k-1))
  enddo

  ee(1) = a(1) / (a(1) - 1.e0)
  gg(1) = (-dtl2 * wusurf / (-dy(1)) - uf(1)) / (a(1) - 1.e0)
  hh(1) = (-dtl2 * wvsurf / (-dy(1)) - vf(1)) / (a(1) - 1.e0)

  do k = 2,kb-2
     ff = 1.e0 / (a(k) + c(k) * (1.e0 - ee(k-1)) - 1.e0)
     ee(k) = a(k) * ff
     gg(k) = (c(k) * gg(k-1) - uf(k)) * ff
     hh(k) = (c(k) * hh(k-1) - vf(k)) * ff
  enddo

  aa = cbc(kb-1) * sqrt(pom%ub(kb-1,ipom)**2 + pom%vb(kb-1,ipom)**2)
  bb = 1.0e0 / (aa * dtl2 / (-dy(kb-1)) - 1.e0 - (ee(kb-2) - 1.e0) * c(kb-1))

  uf(kb-1) = (c(kb-1) * gg(kb-2) - uf(kb-1)) * bb
  vf(kb-1) = (c(kb-1) * hh(kb-2) - vf(kb-1)) * bb

  do k = 2,kb-1
     ki = kb - k
     uf(ki) = (ee(ki) * uf(ki+1) + gg(ki))
     vf(ki) = (ee(ki) * vf(ki+1) + hh(ki))
  enddo

  pom%wubot(ipom) = -aa * uf(kb-1)
  pom%wvbot(ipom) = -aa * vf(kb-1)

  end subroutine profuv

End module pom2k1d

