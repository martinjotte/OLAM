subroutine acmcld_uvmix( iw, dtl )

!-----------------------------------------------------------------------
!  FUNCTION:  Subroutine to compute convective mixing in the CBL
!             according to the Asymmetrical Convective Model (ACM).
!             Ref: Pleim snd Chang (1992)
!
!  SUMMARY:
!   ACM is based on the Blackadar non-local convective model which is
!   used in HIRPBL where upward mixing similar to Blackadar but
!   downward mixing is to the next lower level representing more
!   realistic gradual subsidence.
!-----------------------------------------------------------------------

  use mem_grid,    only: mza, zm, zt, xew, yew, zew, arw0, volti, lpw, dzt
  use mem_cuparm,  only: kcutop, kcubot, vxsrc, vysrc, vzsrc, cbmf, thsrc
  use mem_basic,   only: vxe, vye, vze, rho, theta
  use consts_coms, only: eradi
  use mem_ijtabs,  only: itab_w
  use tridiag,     only: acm_matrix

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc

  real :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real :: di(mza,2), ui(mza,2)
  real :: aflux(mza,2)
  real :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real :: clddepth
  real :: raxis, raxisi
  real :: u(mza), v(mza)
  real :: ut(mza), vt(mza)
  real :: dtom(mza)
  real :: uvtr, dtheta, dztop

  real, parameter :: ef = 0.333 ! mixing efficiency for momentum
  real, parameter :: fl = 0.5   ! extra local mixing

  ! First make sure that there was convection

  if ( kcubot(iw) >= lpw(iw) .and. kcutop(iw) > kcubot(iw) .and. &
         cbmf(iw) > 1.e-6 ) then

     kb = kcubot(iw)
     kt = kcutop(iw)

     ! For some schemes, mixing extends 1 level above cloud
     if (thsrc(kt+1,iw) < -1.e-10) then
        kt = kt + 1
     endif

  else

     ! No significant convection
     return

  endif

  nlev = kt - kb + 1
  dtheta = theta(kt,iw) - theta(kt-1,iw)

! Estimate overshoot/entrainment at top

  if (dtheta < 1.e-5) then
     dztop = dzt(kt)
     clddepth = zm(kt) - zm(kb)
  elseif (thsrc(kt,iw) > 0.0) then
     dztop = zt(kt) - zm(kt-1)
     clddepth = zt(kt) - zm(kb)
  else
     dztop = -1.5e4 * dzt(kt) * thsrc(kt,iw) / dtheta
     dztop = min(dztop, dzt(kt))
     dztop = max(dztop, 0.0)
     clddepth = dztop + zm(kt-1) - zm(kb)
  endif

! Compute ACM missing rates

  massflx_acm(kb-1) = 0.0
  massflx_acm(kb)   = ef * cbmf(iw) * arw0(iw)
  massflx_acm(kt)   = 0.0

  ! Loop over W (flux) levels
  do k = kb+1, kt-1

     massflx_acm(k) = massflx_acm(kb) * rho(k,iw) / rho(kb,iw) * &
                  (clddepth - (zm(k) - zm(kb))) / clddepth
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do k = kb, kt-1
     massflx_loc(k) = massflx_acm(kb) * rho(k,iw) / rho(kb,iw) - massflx_acm(k)
     massflx_loc(k) = min(massflx_loc(k), fl*massflx_acm(k))
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo
     
! Compute zonal and meridional winds

  raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
  raxisi = 1.0 / max(raxis, 1.e-12)

  if (raxis > 1.e3) then
     do k = kb, kt
        u(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
        v(k) = vze(k,iw) * raxis * eradi  &
                - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) &
                * zew(iw) * raxisi * eradi
     enddo
  else
     do k = kb, kt
        u(k) = vxe(k,iw)
        v(k) = vye(k,iw)
     enddo
  endif

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Loop over T levels
  do k = kb, kt
     kc = k - kb + 1
     dtom(k) = dtl * volti(k,iw) / rho(k,iw)

     di(kc,1) = u(k)
     di(kc,2) = v(k)

     bi(kc)   = 1.0 + dtom(k) * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =     - dtom(k) *  massflx_tot(k)
     ai(kc,1) =       dtom(k) * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =     - dtom(k) *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + dtom(kb) * massflx_acm(kb)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, 2, 1 )

! Now, soln contains future(t+1) values. Compute fluxes
  
  ! Loop over W (flux) levels
  do k = kb, kt-1
     kc = k - kb + 1

     aflux(k,1) = massflx_acm(k) * ui(1,1) + massflx_loc(k) * ui(kc,1) &
                - massflx_tot(k) * ui(kc+1,1)
     aflux(k,2) = massflx_acm(k) * ui(1,2) + massflx_loc(k) * ui(kc,2) &
                - massflx_tot(k) * ui(kc+1,2)
  enddo

  aflux(kb-1,1) = 0.0
  aflux(kb-1,2) = 0.0

  aflux(kt,1) = 0.0
  aflux(kt,2) = 0.0

! Compute tendencies due to internal convective fluxes

  ! Loop over T levels
  do k = kb, kt
     ut(k) = (aflux(k-1,1) - aflux(k,1)) * volti(k,iw) / rho(k,iw)
     vt(k) = (aflux(k-1,2) - aflux(k,2)) * volti(k,iw) / rho(k,iw)
  enddo

! Convert U,V tendencies to earth cartesian grid

  if (raxis > 1.e3) then

     do k = kb, kt
        uvtr = -vt(k) * zew(iw) * eradi
        vxsrc(k,iw) = (-ut(k) * yew(iw) + uvtr * xew(iw)) * raxisi
        vysrc(k,iw) = ( ut(k) * xew(iw) + uvtr * yew(iw)) * raxisi
        vzsrc(k,iw) =   vt(k) * raxis * eradi 
     enddo

  else

     do k = kb, kt
        vxsrc(k,iw) = ut(k)
        vysrc(k,iw) = vt(k)
        vzsrc(k,iw) = 0.0
     enddo

  endif

end subroutine acmcld_uvmix



subroutine acmcld_tracermix( iw, dtl )

!-----------------------------------------------------------------------
!  FUNCTION:  Subroutine to compute convective mixing in the CBL
!             according to the Asymmetrical Convective Model (ACM).
!             Ref: Pleim snd Chang (1992)
!
!  SUMMARY:
!   ACM is based on the Blackadar non-local convective model which is
!   used in HIRPBL where upward mixing similar to Blackadar but
!   downward mixing is to the next lower level representing more
!   realistic gradual subsidence.
!-----------------------------------------------------------------------

  use mem_grid,    only: mza, zm, zt, xew, yew, zew, arw0, volti, lpw, dzt
  use mem_cuparm,  only: kcutop, kcubot, vxsrc, vysrc, vzsrc, cbmf, thsrc
  use mem_basic,   only: vxe, vye, vze, rho, theta
  use consts_coms, only: eradi
  use mem_ijtabs,  only: itab_w
  use tridiag,     only: acm_matrix
  use var_tables,  only: num_cumix, cumix_map, scalar_tab
  use oname_coms,  only: nl

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc, ns, nc

  real :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real :: di(mza,num_cumix), ui(mza,num_cumix)
  real :: aflux(mza)
  real :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real :: clddepth
  real :: dtom(mza)
  real :: dtheta, dztop

  real, parameter :: ef = 1.0  ! mixing efficiency for tracers
  real, parameter :: fl = 0.5  ! extra local mixing

  ! Do we have any tracers to mix?

  if (num_cumix < 1) return

  ! First make sure that there was convection

  if ( kcubot(iw) >= lpw(iw) .and. kcutop(iw) > kcubot(iw) .and. &
         cbmf(iw) > 1.e-6 ) then

     kb = max( kcubot(iw) - 1, lpw(iw) )
     kt = kcutop(iw)

     ! For some schemes, mixing extends 1 level above cloud
     if (thsrc(kt+1,iw) < -1.e-10) then
        kt = kt + 1
     endif

  else

     ! No significant convection
     return

  endif

  nlev = kt - kb + 1

  ! Assume cloud extends to middle of top layer  
  clddepth = zt(kt) - zm(kb)

  ! Compute ACM missing rates

  massflx_acm(kb-1) = 0.0
  massflx_acm(kb)   = ef * cbmf(iw) * arw0(iw)
  massflx_acm(kt)   = 0.0

  ! Loop over W (flux) levels
  do k = kb+1, kt-1
     massflx_acm(k) = massflx_acm(kb) * rho(k,iw) / rho(kb,iw) * &
                  (clddepth - (zm(k) - zm(kb))) / clddepth
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do k = kb, kt-1
     massflx_loc(k) = massflx_acm(kb) * rho(k,iw) / rho(kb,iw) - massflx_acm(k)
     massflx_loc(k) = min(massflx_loc(k), fl*massflx_acm(k))
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo
     
! Define arrays A,B,E which make up matrix and D which is RHS

  do nc = 1, num_cumix
     ns = cumix_map(nc)
     if (nl%split_scalars > 0) then
        do k = kb, kt
           kc = k - kb + 1
           di(kc,nc) = scalar_tab(ns)%var_p(k,iw) + dtl * scalar_tab(ns)%var_t(k,iw) / rho(k,iw)
        enddo
     else
        do k = kb, kt
           kc = k - kb + 1
           di(kc,nc) = scalar_tab(ns)%var_p(k,iw)
        enddo
     endif
  enddo

  ! Loop over T levels
  do k = kb, kt
     kc = k - kb + 1
     dtom(k) = dtl * volti(k,iw) / rho(k,iw)

     bi(kc)   = 1.0 + dtom(k) * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =     - dtom(k) *  massflx_tot(k)
     ai(kc,1) =       dtom(k) * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =     - dtom(k) *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + dtom(kb) * massflx_acm(kb)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, num_cumix, 1 )

! Now, soln contains future(t+1) values. Compute fluxes

  do nc = 1, num_cumix
     ns = cumix_map(nc)

     ! Loop over W (flux) levels
     do k = kb, kt-1
        kc = k - kb + 1

        aflux(k) = massflx_acm(k) * ui(1,nc) + massflx_loc(k) * ui(kc,nc) &
                 - massflx_tot(k) * ui(kc+1,nc)
     enddo

     aflux(kb-1) = 0.0
     aflux(kt)   = 0.0

! Compute tendencies due to internal convective fluxes

     ! Loop over T levels
     do k = kb, kt
        scalar_tab(ns)%var_t(k,iw) = scalar_tab(ns)%var_t(k,iw) &
                                   + (aflux(k-1) - aflux(k)) * volti(k,iw)
     enddo
  enddo

end subroutine acmcld_tracermix
