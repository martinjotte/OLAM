subroutine acmcld_uvmix( iw, dtl )

  use mem_cuparm, only: iactcu, cbmf, cddf

  implicit none

  ! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

  ! Mixing due to convective updrafts

  if (iactcu(iw) > 0 .and. cbmf(iw) > 1.e-6 ) then
     
     call acmcld_uvmix_up( iw, dtl )

  endif

  ! Mixing due to convective downdrafts

  if (iactcu(iw) > 0 .and. cddf(iw) > 1.e-6 ) then
     
     call acmcld_uvmix_dd( iw, dtl )

  endif

end subroutine acmcld_uvmix



subroutine acmcld_uvmix_up( iw, dtl )

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

  use mem_grid,    only: mza, zm, zt, xew, yew, zew, arw0, volt, volti, lpw
  use mem_cuparm,  only: kcutop, kcubot, vxsrc, vysrc, vzsrc, cbmf
  use mem_basic,   only: vxe, vye, vze, rho, theta
  use consts_coms, only: eradi
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
  real :: dtom
  real :: uvtr

  real, parameter :: ef = 0.333 ! mixing efficiency for momentum
  real, parameter :: fl = 0.5   ! extra local mixing

  ! First make sure that there was convection

  if ( kcubot(iw) >= lpw(iw) .and. kcutop(iw) >= kcubot(iw) ) then

     kb = max( kcubot(iw) - 1, lpw(iw) )
     kt = kcutop(iw)

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

     massflx_acm(k) = massflx_acm(kb) * real(rho(k,iw)) / real(rho(kb,iw)) &
                    * (clddepth - (zm(k) - zm(kb))) /  clddepth
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do k = kb, kt-1
     massflx_loc(k) = massflx_acm(kb) * real(rho(k,iw)) / real(rho(kb,iw)) - massflx_acm(k)
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
     dtom = dtl / real(volt(k,iw) * rho(k,iw))

     di(kc,1) = u(k)
     di(kc,2) = v(k)

     bi(kc)   = 1.0 + dtom * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =     - dtom *  massflx_tot(k)
     ai(kc,1) =       dtom * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =     - dtom *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, 2, 1 )

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
     ut(k) = (aflux(k-1,1) - aflux(k,1)) * real(volti(k,iw))
     vt(k) = (aflux(k-1,2) - aflux(k,2)) * real(volti(k,iw))
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
     enddo

  endif

end subroutine acmcld_uvmix_up



subroutine acmcld_uvmix_dd( iw, dtl )

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

  use mem_grid,    only: mza, zm, zt, xew, yew, zew, arw, volt, volti, lpw
  use mem_cuparm,  only: kddtop, vxsrc, vysrc, vzsrc, cddf
  use mem_basic,   only: vxe, vye, vze, rho, theta
  use consts_coms, only: eradi
  use tridiag,     only: acm_matrix

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc

  real :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real :: di(mza,2), ui(mza,2)
  real :: aflux(0:mza,2)
  real :: massflx_acm(0:mza), massflx_loc(0:mza), massflx_tot(0:mza)
  real :: clddepth
  real :: raxis, raxisi
  real :: u(mza), v(mza)
  real :: ut(mza), vt(mza)
  real :: dtom
  real :: uvtr

  real, parameter :: ef = 0.333 ! mixing efficiency for momentum
  real, parameter :: fl = 0.5   ! extra local mixing

  ! First make sure that there was convection

  if ( kddtop(iw) > lpw(iw) ) then

     nlev = kddtop(iw) - lpw(iw) + 1
     kb   = 1
     kt   = nlev

  else

     ! No significant convection
     return

  endif

  ! Assume cloud extends to bottom of top layer  
  clddepth = zm( kddtop(iw) - 1) - zm( lpw(iw) - 1 )

! Compute ACM missing rates

  massflx_acm(kb-1) = 0.0
  massflx_acm(kb)   = ef * cddf(iw) * arw(kddtop(iw),iw)
  massflx_acm(kt)   = 0.0

  ! Loop over W (flux) levels
  do kc = kb+1, kt-1
     k  = kddtop(iw) - kc + 1

     massflx_acm(kc) = massflx_acm(kb) * real(rho(k,iw)) * arw(k,iw) &
                     * (clddepth - (zm(kddtop(iw)-1)-zm(k-1))) &
                     / ( real(rho(kddtop(iw),iw)) * arw(kddtop(iw),iw) * clddepth )
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do kc = kb, kt-1
     k  = kddtop(iw) - kc + 1
     massflx_loc(kc) = massflx_acm(kb) * real(rho(k,iw)) / real(rho(kddtop(iw),iw)) - massflx_acm(kc)
     massflx_loc(kc) = min(massflx_loc(kc), fl*massflx_acm(kc))
     massflx_tot(kc) = massflx_acm(kc) + massflx_loc(kc)
  enddo

! Compute zonal and meridional winds

  raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
  raxisi = 1.0 / max(raxis, 1.e-12)

  if (raxis > 1.e3) then
     do k = lpw(iw), kddtop(iw)
        u(k) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
        v(k) = vze(k,iw) * raxis * eradi  &
                - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) &
                * zew(iw) * raxisi * eradi
     enddo
  else
     do k = lpw(iw), kddtop(iw)
        u(k) = vxe(k,iw)
        v(k) = vye(k,iw)
     enddo
  endif

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Loop over T levels
  do kc = kb, kt
     k  = kddtop(iw) - kc + 1
     dtom = dtl / real(volt(k,iw) * rho(k,iw))

     di(kc,1) = u(k)
     di(kc,2) = v(k)

     bi(kc)   = 1.0 + dtom * (massflx_tot(kc-1) + massflx_loc(kc))
     ei(kc)   =     - dtom *  massflx_tot(kc)
     ai(kc,1) =       dtom * (massflx_acm(kc) - massflx_acm(kc-1))
     ci(kc)   =     - dtom *  massflx_loc(kc-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, 2, 1 )

! Now, soln contains future(t+1) values. Compute fluxes
  
  ! Loop over W (flux) levels
  do kc = kb, kt-1

     aflux(kc,1) = massflx_acm(kc) * ui(1,1) + massflx_loc(kc) * ui(kc,1) &
                 - massflx_tot(kc) * ui(kc+1,1)
     aflux(kc,2) = massflx_acm(kc) * ui(1,2) + massflx_loc(kc) * ui(kc,2) &
                 - massflx_tot(kc) * ui(kc+1,2)
  enddo

  aflux(kb-1,1) = 0.0
  aflux(kb-1,2) = 0.0

  aflux(kt,1) = 0.0
  aflux(kt,2) = 0.0

! Compute tendencies due to internal convective fluxes

  ! Loop over T levels
  do kc = kb, kt
     k  = kddtop(iw) - kc + 1
     ut(k) = (aflux(kc-1,1) - aflux(kc,1)) * real(volti(k,iw))
     vt(k) = (aflux(kc-1,2) - aflux(kc,2)) * real(volti(k,iw))
  enddo

! Convert U,V tendencies to earth cartesian grid

  if (raxis > 1.e3) then

     do k = lpw(iw), kddtop(iw)
        uvtr = -vt(k) * zew(iw) * eradi
        vxsrc(k,iw) = vxsrc(k,iw) + (-ut(k) * yew(iw) + uvtr * xew(iw)) * raxisi
        vysrc(k,iw) = vysrc(k,iw) + ( ut(k) * xew(iw) + uvtr * yew(iw)) * raxisi
        vzsrc(k,iw) = vzsrc(k,iw) +   vt(k) * raxis * eradi 
     enddo

  else

     do k = lpw(iw), kddtop(iw)
        vxsrc(k,iw) = vxsrc(k,iw) + ut(k)
        vysrc(k,iw) = vysrc(k,iw) + vt(k)
     enddo

  endif

end subroutine acmcld_uvmix_dd



subroutine acmcld_tracermix( iw, dtl )

  use mem_cuparm, only: iactcu, cbmf, cddf
  use var_tables, only: num_cumix
  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

  ! Do we have any tracers to mix?

  if (num_cumix < 1) return

  ! Mixing due to convective updrafts

  if (iactcu(iw) > 0 .and. cbmf(iw) > 1.e-6 ) then
     
     call acmcld_tracermix_up( iw, dtl )

  endif

  ! Mixing due to convective downdrafts

  if (iactcu(iw) > 0 .and. cddf(iw) > 1.e-6 ) then
     
     call acmcld_tracermix_dd( iw, dtl )

  endif

end subroutine acmcld_tracermix



subroutine acmcld_tracermix_up( iw, dtl )

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

  use mem_grid,    only: mza, zm, zt, arw0, volt, volti, lpw
  use mem_cuparm,  only: kcutop, kcubot, cbmf
  use mem_basic,   only: rho
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
  real :: dtom

  real, parameter :: ef = 1.0  ! mixing efficiency for tracers
  real, parameter :: fl = 0.5  ! extra local mixing

  ! First make sure that there was convection

  if ( kcubot(iw) >= lpw(iw) .and. kcutop(iw) >= kcubot(iw) ) then

     kb = max( kcubot(iw) - 1, lpw(iw) )
     kt = kcutop(iw)

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
     massflx_acm(k) = massflx_acm(kb) * real(rho(k,iw)) / real(rho(kb,iw)) &
                    * (clddepth - (zm(k) - zm(kb))) / clddepth
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do k = kb, kt-1
     massflx_loc(k) = massflx_acm(kb) * real(rho(k,iw)) / real(rho(kb,iw)) - massflx_acm(k)
     massflx_loc(k) = min(massflx_loc(k), fl*massflx_acm(k))
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo
     
! Define arrays A,B,E which make up matrix and D which is RHS

  do nc = 1, num_cumix
     ns = cumix_map(nc)
     do k = kb, kt
        kc = k - kb + 1
        di(kc,nc) = scalar_tab(ns)%var_p(k,iw) + dtl * scalar_tab(ns)%var_t(k,iw) / real(rho(k,iw))
     enddo
  enddo

  ! Loop over T levels
  do k = kb, kt
     kc = k - kb + 1
     dtom = dtl / real(volt(k,iw) * rho(k,iw))

     bi(kc)   = 1.0 + dtom * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =     - dtom *  massflx_tot(k)
     ai(kc,1) =       dtom * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =     - dtom *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_cumix, 1 )

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
                                   + (aflux(k-1) - aflux(k)) * real(volti(k,iw))
     enddo
  enddo

end subroutine acmcld_tracermix_up




subroutine acmcld_tracermix_dd( iw, dtl )

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

  use mem_grid,    only: mza, zm, zt, arw, volt, volti, lpw
  use mem_cuparm,  only: kddtop, cddf
  use mem_basic,   only: rho
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
  real :: aflux(0:mza)
  real :: massflx_acm(0:mza), massflx_loc(0:mza), massflx_tot(0:mza)
  real :: clddepth
  real :: dtom

  real, parameter :: ef = 1.0  ! mixing efficiency for tracers
  real, parameter :: fl = 0.5  ! extra local mixing

  ! First make sure that there was convection

  if ( kddtop(iw) > lpw(iw) ) then

     nlev = kddtop(iw) - lpw(iw) + 1
     kb   = 1
     kt   = nlev

  else

     ! No significant convection
     return

  endif

  ! Assume cloud extends to bottom of top layer  
  clddepth = zm( kddtop(iw) - 1) - zm( lpw(iw) - 1 )

  ! Compute ACM missing rates

  massflx_acm(kb-1) = 0.0
  massflx_acm(kb)   = ef * cddf(iw) * arw(kddtop(iw),iw)
  massflx_acm(kt)   = 0.0

  ! Loop over W (flux) levels
  do kc = kb+1, kt-1
     k  = kddtop(iw) - kc + 1

     massflx_acm(kc) = massflx_acm(kb) * real(rho(k,iw)) * arw(k,iw) &
                     * (clddepth - (zm(kddtop(iw)-1)-zm(k-1))) &
                     / ( real(rho(kddtop(iw),iw)) * arw(kddtop(iw),iw) * clddepth )
  enddo

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

  do kc = kb, kt-1
     k  = kddtop(iw) - kc + 1
     massflx_loc(kc) = massflx_acm(kb) * rho(k,iw) / rho(kddtop(iw),iw) - massflx_acm(kc)
     massflx_loc(kc) = min(massflx_loc(kc), fl*massflx_acm(kc))
     massflx_tot(kc) = massflx_acm(kc) + massflx_loc(kc)
  enddo

! Define arrays A,B,E which make up matrix and D which is RHS

  do nc = 1, num_cumix
     ns = cumix_map(nc)
     do kc = kb, kt
        k  = kddtop(iw) - kc + 1
        di(kc,nc) = scalar_tab(ns)%var_p(k,iw) + dtl * scalar_tab(ns)%var_t(k,iw) / real(rho(k,iw))
     enddo
  enddo

  ! Loop over T levels
  do kc = kb, kt
     k  = kddtop(iw) - kc + 1
     dtom = dtl / real(volt(k,iw) * rho(k,iw))

     bi(kc)   = 1.0 + dtom * (massflx_tot(kc-1) + massflx_loc(kc))
     ei(kc)   =     - dtom *  massflx_tot(kc)
     ai(kc,1) =       dtom * (massflx_acm(kc) - massflx_acm(kc-1))
     ci(kc)   =     - dtom *  massflx_loc(kc-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_cumix, 1 )

! Now, soln contains future(t+1) values. Compute fluxes

  do nc = 1, num_cumix
     ns = cumix_map(nc)

     ! Loop over W (flux) levels
     do kc = kb, kt-1

        aflux(kc) = massflx_acm(kc) * ui(1,nc) + massflx_loc(kc) * ui(kc,nc) &
                  - massflx_tot(kc) * ui(kc+1,nc)
     enddo

     aflux(kb-1) = 0.0
     aflux(kt)   = 0.0

! Compute tendencies due to internal convective fluxes

     ! Loop over T levels
     do kc = kb, kt
        k  = kddtop(iw) - kc + 1
        scalar_tab(ns)%var_t(k,iw) = scalar_tab(ns)%var_t(k,iw) &
                                   + (aflux(kc-1) - aflux(kc)) * real(volti(k,iw))
     enddo
  enddo

end subroutine acmcld_tracermix_dd
