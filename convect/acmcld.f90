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

  use mem_grid,    only: mza, zm, zt, arw, volti, lpw, &
                         vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
  use mem_cuparm,  only: kcutop, kcubot, vxsrc, vysrc, vzsrc, cbmf
  use mem_basic,   only: ue, ve, rho
  use consts_coms, only: r8
  use tridiag,     only: acm_matrix

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc

  real(r8) :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,2), ui(mza,2)

  real :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real :: clddepth, massflx0, massflxk
  real :: ut, vt, dtom, dti

  real, parameter :: ef = 0.333 ! mixing efficiency for momentum
  real, parameter :: fl = 0.5   ! extra local mixing

  ! First make sure that there was convection

  if ( kcubot(iw) < lpw(iw) .or. kcutop(iw) < kcubot(iw) &
       .or. cbmf(iw) < 1.e-10 ) return

  kb = max( kcubot(iw) - 1, lpw(iw) )
  kt = kcutop(iw)
  nlev = kt - kb + 1
  dti = 1. / dtl

  ! Assume cloud extends to middle of top layer

  clddepth = zt(kt) - zm(kb)

  ! Compute ACM missing rates

  massflx0 = ef * cbmf(iw) / real( rho(kb,iw) )

  ! Loop over W (flux) levels
  do k = kb, kt-1
     massflxk       = massflx0 * real(rho(k,iw)) * arw(k,iw)
     massflx_acm(k) = massflxk * (zt(kt) - zm(k)) / clddepth
     massflx_loc(k) = min(massflxk - massflx_acm(k), fl * massflx_acm(k))
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo

  massflx_acm(kb-1) = 0.0
  massflx_acm(kt)   = 0.0

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Loop over T levels
  do k = kb, kt
     kc = k - kb + 1
     dtom = dtl * volti(k,iw) / real(rho(k,iw))

     di(kc,1) = ue(k,iw)
     di(kc,2) = ve(k,iw)

     bi(kc)   = 1._r8 + dtom * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =       - dtom *  massflx_tot(k)
     ai(kc,1) =         dtom * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =       - dtom *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, 2, 1 )

  ! Convert U,V tendencies to earth cartesian grid

  do k = kb, kt
     kc = k - kb + 1

     ut = real( ui(kc,1) - di(kc,1) ) * dti
     vt = real( ui(kc,2) - di(kc,2) ) * dti

     vxsrc(k,iw) = ut * vxn_ew(iw) + vt * vxn_ns(iw)
     vysrc(k,iw) = ut * vyn_ew(iw) + vt * vyn_ns(iw)
     vzsrc(k,iw) =                   vt * vzn_ns(iw)
  enddo

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

  use mem_grid,    only: mza, zm, arw, volti, lpw, &
                         vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
  use mem_cuparm,  only: kddtop, vxsrc, vysrc, vzsrc, cddf
  use mem_basic,   only: ue, ve, rho
  use consts_coms, only: r8
  use tridiag,     only: acm_matrix

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc

  real(r8) :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,2), ui(mza,2)

  real :: massflx_acm(0:mza), massflx_loc(0:mza), massflx_tot(0:mza)
  real :: clddepth, massflx0, massflxk
  real :: ut, vt, dtom, dti

  real, parameter :: ef = 0.333 ! mixing efficiency for momentum
  real, parameter :: fl = 0.5   ! extra local mixing

  ! First make sure that there was convection

  if ( kddtop(iw) < lpw(iw) .or. cddf(iw) < 1.e-10) return

  nlev = kddtop(iw) - lpw(iw) + 1
  kb   = 1
  kt   = nlev
  dti  = 1.0 / dtl

  ! Assume cloud extends to bottom of top layer

  clddepth = zm(kddtop(iw) - 1) - zm(lpw(iw) - 1)

  ! Compute ACM missing rates

  massflx0 = ef * cddf(iw) / real( rho(kddtop(iw),iw) )

  ! Loop over W (flux) levels
  do kc = kb, kt-1
     k  = kddtop(iw) - kc + 1

     massflxk        = massflx0 * real(rho(k,iw)) * arw(k,iw)
     massflx_acm(kc) = massflxk * (zm(k-1) - zm(lpw(iw)-1)) / clddepth
     massflx_loc(kc) = min(massflxk - massflx_acm(kc), fl * massflx_acm(kc))
     massflx_tot(kc) = massflx_acm(kc) + massflx_loc(kc)
  enddo

  massflx_acm(kb-1) = 0.0
  massflx_acm(kt)   = 0.0

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Loop over T levels
  do kc = kb, kt
     k  = kddtop(iw) - kc + 1
     dtom = dtl * volti(k,iw) / real(rho(k,iw))

     di(kc,1) = ue(k,iw)
     di(kc,2) = ve(k,iw)

     bi(kc)   = 1._r8 + dtom * (massflx_tot(kc-1) + massflx_loc(kc))
     ei(kc)   =       - dtom *  massflx_tot(kc)
     ai(kc,1) =         dtom * (massflx_acm(kc) - massflx_acm(kc-1))
     ci(kc)   =       - dtom *  massflx_loc(kc-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, 2, 1 )

! Convert U,V tendencies to earth cartesian grid

  do k = lpw(iw), kddtop(iw)
     kc = kddtop(iw) - k + 1

     ut = real( ui(kc,1) - di(kc,1) ) * dti
     vt = real( ui(kc,2) - di(kc,2) ) * dti

     vxsrc(k,iw) = vxsrc(k,iw) + ut * vxn_ew(iw) + vt * vxn_ns(iw)
     vysrc(k,iw) = vysrc(k,iw) + ut * vyn_ew(iw) + vt * vyn_ns(iw)
     vzsrc(k,iw) = vzsrc(k,iw)                   + vt * vzn_ns(iw)
  enddo

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

  use mem_grid,    only: mza, zm, zt, arw, volti, lpw
  use mem_cuparm,  only: kcutop, kcubot, cbmf
  use mem_basic,   only: rho
  use consts_coms, only: r8
  use tridiag,     only: acm_matrix
  use var_tables,  only: num_cumix, cumix_map, scalar_tab

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc, ns, nc

  real(r8) :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,num_cumix), ui(mza,num_cumix)

  real :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real :: clddepth, massflx0, massflxk
  real :: dtom, dti

  real, parameter :: ef = 1.0  ! mixing efficiency for tracers
  real, parameter :: fl = 0.5  ! extra local mixing

  ! First make sure that there was convection

  if ( kcubot(iw) < lpw(iw) .or. kcutop(iw) < kcubot(iw) &
       .or. cbmf(iw) < 1.e-10 ) return

  kb = max( kcubot(iw) - 1, lpw(iw) )
  kt = kcutop(iw)
  nlev = kt - kb + 1
  dti = 1.0 / dtl

  ! Assume cloud extends to middle of top layer  

  clddepth = zt(kt) - zm(kb)

  ! Compute ACM missing rates

  massflx0 = ef * cbmf(iw) / real( rho(kb,iw) )

  ! Loop over W (flux) levels
  do k = kb, kt-1
     massflxk       = massflx0 * real(rho(k,iw)) * arw(k,iw)
     massflx_acm(k) = massflxk * (zt(kt) - zm(k)) / clddepth
     massflx_loc(k) = min(massflxk - massflx_acm(k), fl * massflx_acm(k))
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo

  massflx_acm(kb-1) = 0.0
  massflx_acm(kt)   = 0.0

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

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
     dtom = dtl * volti(k,iw) / real(rho(k,iw))

     bi(kc)   = 1._r8 + dtom * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =       - dtom *  massflx_tot(k)
     ai(kc,1) =         dtom * (massflx_acm(k) - massflx_acm(k-1))
     ci(kc)   =       - dtom *  massflx_loc(k-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_cumix, 1 )

! Now, soln contains future(t+1) values. Compute fluxes

  do nc = 1, num_cumix
     ns = cumix_map(nc)

     ! Loop over T levels
     do k = kb, kt
        kc = k - kb + 1
        scalar_tab(ns)%var_t(k,iw) = scalar_tab(ns)%var_t(k,iw) &
                                   + dti * real( ui(kc,nc) - di(kc,nc) )
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

  use mem_grid,    only: mza, zm, arw, volti, lpw
  use mem_cuparm,  only: kddtop, cddf
  use mem_basic,   only: rho
  use consts_coms, only: r8
  use tridiag,     only: acm_matrix
  use var_tables,  only: num_cumix, cumix_map, scalar_tab

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, nlev, kc, ns, nc

  real(r8) :: ai(mza,1), bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,num_cumix), ui(mza,num_cumix)

  real :: massflx_acm(0:mza), massflx_loc(0:mza), massflx_tot(0:mza)
  real :: clddepth, massflx0, massflxk
  real :: dtom, dti

  real, parameter :: ef = 1.0  ! mixing efficiency for tracers
  real, parameter :: fl = 0.5  ! extra local mixing

  ! First make sure that there was convection

  if ( kddtop(iw) < lpw(iw) .or. cddf(iw) < 1.e-10 ) return

  nlev = kddtop(iw) - lpw(iw) + 1
  kb   = 1
  kt   = nlev
  dti  = 1.0 / dtl

  ! Assume downdraft extends to surface

  clddepth = zm(kddtop(iw) - 1) - zm(lpw(iw) - 1)

  ! Compute ACM missing rates

  massflx0 = ef * cddf(iw) / real( rho(kddtop(iw),iw) )

  ! Loop over W (flux) levels
  do kc = kb, kt-1
     k  = kddtop(iw) - kc + 1

     massflxk        = massflx0 * real(rho(k,iw)) * arw(k,iw)
     massflx_acm(kc) = massflxk * (zm(k-1) - zm(lpw(iw)-1)) / clddepth
     massflx_loc(kc) = min(massflxk - massflx_acm(kc), fl * massflx_acm(kc))
     massflx_tot(kc) = massflx_acm(kc) + massflx_loc(kc)
  enddo

  massflx_acm(kb-1) = 0.0
  massflx_acm(kt)   = 0.0

  massflx_loc(kb-1) = 0.0
  massflx_loc(kt)   = 0.0

  massflx_tot(kb-1) = 0.0
  massflx_tot(kt)   = 0.0

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
     dtom = dtl * volti(k,iw) / real(rho(k,iw))

     bi(kc)   = 1._r8 + dtom * (massflx_tot(kc-1) + massflx_loc(kc))
     ei(kc)   =       - dtom *  massflx_tot(kc)
     ai(kc,1) =         dtom * (massflx_acm(kc) - massflx_acm(kc-1))
     ci(kc)   =       - dtom *  massflx_loc(kc-1)
  enddo

  bi(1) = bi(1) + ai(1,1)
  ci(2) = ci(2) + ai(2,1)

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_cumix, 1 )

! Now, soln contains future(t+1) values. Compute tendencies

  do nc = 1, num_cumix
     ns = cumix_map(nc)

     ! Loop over T levels
     do k = lpw(iw), kddtop(iw)
        kc = kddtop(iw) - k + 1
        scalar_tab(ns)%var_t(k,iw) = scalar_tab(ns)%var_t(k,iw) &
                                   + dti * real( ui(kc,nc) - di(kc,nc) )
     enddo
  enddo

end subroutine acmcld_tracermix_dd
