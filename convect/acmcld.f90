subroutine acmcld_mix( iw, dtl )

  use mem_cuparm,  only: iactcu, cbmf, cddf
  implicit none

  ! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

  ! Mixing due to convective updrafts

  if (iactcu(iw) > 0 .and. cbmf(iw) > 1.e-8) call acmcldmix_up( iw, dtl )

  ! Mixing due to convective downdrafts

  if (iactcu(iw) > 0 .and. cddf(iw) > 1.e-8 ) call acmcldmix_dd( iw, dtl )

end subroutine acmcld_mix



subroutine acmcldmix_up( iw, dtl )

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
                          vxn_ew, vxn_ns, vyn_ew, vyn_ns, vzn_ns
  use mem_cuparm,  only: kcutop, kcubot, kudbot, cbmf
  use mem_basic,   only: rho, ue, ve
  use tridiag,     only: acm_matrix
  use var_tables,  only: scalar_tab, num_cumix, cumix_map
  use mem_tend,    only: vmxet, vmyet, vmzet
  use consts_coms, only: r8
  use oname_coms,  only: nl

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, ka, nlev, kc, ns, n, kk
  integer :: num_mix, nsrc

  real(r8) :: bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,num_cumix+2), ui(mza,num_cumix+2)
  real(r8) :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real(r8) :: massflx

  real :: dtom(mza)
  real :: massflx0, massflxk, detrain
  real :: ut, vt, dti, frac, fsum, fact

  real, parameter :: efv = 1.0  ! mixing efficiency for momentum
  real, parameter :: fl  = 0.7  ! extra local mixing

  real,     allocatable :: fracsrc(:)
  real,     allocatable :: depthi (:)
  real(r8), allocatable :: ai   (:,:)
  real(r8), allocatable :: temp (:,:)

  ! First make sure there is convection

  if ( kcubot(iw) < lpw(iw) .or. kcutop(iw) < kcubot(iw) ) return

  kt = kcutop(iw)
  kb = min( kcubot(iw)+1, kcutop(iw)-1 )
  ka = max(lpw(iw), min( kudbot(iw), kcubot(iw)-1 ) )

  dti = 1.0 / dtl

  nlev = kt - ka + 1
  nsrc = kb - ka + 1

  if (nl%conv_uv_mix == 1) then
     num_mix = num_cumix + 2
  else
     num_mix = num_cumix
  endif

  allocate(fracsrc (nsrc))
  allocate(depthi  (nsrc))
  allocate(ai  (mza,nsrc))
  allocate(temp(mza,nsrc))

  detrain = 0.08 * (zm(kcutop(iw)) - zm(kcubot(iw)-1))

  ! Compute ACM missing rates

  massflx0 = cbmf(iw) / real( rho(kb,iw) )

  if (nsrc > 1) then

     fsum = 0.0
     do kk = 1, nsrc
        k  = kk + ka - 1
        depthi(kk) = 1.0 / (zm(kt) - zm(kk+ka-1) + detrain)
        fsum = fsum + arw(k,iw) * depthi(kk)
     enddo

     frac = arw(kb,iw) / ( (zm(kt) - zm(kb) + detrain) * fsum )

  else

     depthi(1) = 1.0 / (zm(kt) - zm(ka) + detrain)
     frac      = 1.0

  endif

  temp        = 0.0_r8
  massflx_loc = 0.0_r8

  massflxk = frac * massflx0

  ! assumes arw always increases or stays the same going up
  do kk = 1, nsrc
     fact = massflxk * arw(kk+ka-1,iw) * depthi(kk)
     do k = kk + ka - 1, kt-1
        temp(k,kk) = rho(k,iw) * fact * (zm(kt) - zm(k) + detrain)
     enddo
  enddo

  do k = ka, kt-1
     massflx_acm(k) = sum(temp(k,1:nsrc))
  enddo

  do k = kb, kt-1
     massflx        = rho(k,iw) * massflx0 * arw(k,iw)
     massflx_loc(k) = min(massflx - massflx_acm(k), fl * massflx_acm(k))
  enddo

  do k = ka, kt-1
     massflx_loc(k) = max( massflx_loc(k), 0.05_r8 * massflx_acm(k) )
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo

  massflx_acm(ka-1) = 0.0_r8
  massflx_acm(kt)   = 0.0_r8

  massflx_tot(ka-1) = 0.0_r8
  massflx_tot(kt)   = 0.0_r8

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Load scalars into RHS array

  do k = ka, kt
     dtom (k) = dtl * volti(k,iw) / real(rho(k,iw))
  enddo

  !dir$ nounroll_and_jam
  do ns = 1, num_cumix
     n = cumix_map(ns)
     do k = ka, kt
        kc = k - ka + 1
        di(kc,ns) = scalar_tab(n)%var_p(k,iw) &
                  + scalar_tab(n)%var_t(k,iw) * dtl / real(rho(k,iw))
     enddo
  enddo

  ! Load earth-Cartesian velocities into RHS array

  do k = ka, kt
     kc = k - ka + 1
     di(kc,num_cumix+1) = ue(k,iw)
     di(kc,num_cumix+2) = ve(k,iw)
  enddo

  ! Loop over T levels
  do k = ka, kt
     kc = k - ka + 1
     bi(kc)   = 1._r8 + dtom(k) * (massflx_tot(k-1) + massflx_loc(k))
     ei(kc)   =       - dtom(k) *  massflx_tot(k)
     ci(kc)   =       - dtom(k) *  massflx_loc(k-1)
  enddo

  do kk = 1, nsrc
     do k = kk + ka, kt
        kc = k - ka + 1
        ai(kc,kk) = dtom(k) * (temp(k,kk) - temp(k-1,kk))
     enddo
  enddo

  do kc = 1, nsrc
     k  = kc + ka - 1
     bi(kc)   = bi(kc)   + dtom(k) * temp(k,kc)
     ci(kc+1) = ci(kc+1) + ai(kc+1,kc)
  enddo

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_mix, nsrc )

  ! Compute tendencies due to convective fluxes

  !dir$ nounroll_and_jam
  do ns = 1, num_cumix
     n = cumix_map(ns)
     do k = ka, kt
        kc = k - ka + 1
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + dti * real( rho(k,iw) * (ui(kc,ns) - di(kc,ns)) )
     enddo
  enddo

  if (nl%conv_uv_mix == 1) then
     do k = ka, kt
        kc = k - ka + 1

        ut = dti * real( rho(k,iw) * (ui(kc,num_cumix+1) - di(kc,num_cumix+1)) )
        vt = dti * real( rho(k,iw) * (ui(kc,num_cumix+2) - di(kc,num_cumix+2)) )

        vmxet(k,iw) = vmxet(k,iw) + ut * vxn_ew(iw) + vt * vxn_ns(iw)

        vmyet(k,iw) = vmyet(k,iw) + ut * vyn_ew(iw) + vt * vyn_ns(iw)

        vmzet(k,iw) = vmzet(k,iw)                   + vt * vzn_ns(iw)
     enddo
  endif

end subroutine acmcldmix_up




subroutine acmcldmix_dd( iw, dtl )

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
                         vxn_ew, vxn_ns, vyn_ew, vyn_ns, vzn_ns
  use mem_cuparm,  only: kddtop, kddmax, kddbot, cddf
  use mem_basic,   only: rho, ue, ve
  use tridiag,     only: acm_matrix
  use var_tables,  only: scalar_tab, num_cumix, cumix_map
  use mem_tend,    only: vmxet, vmyet, vmzet
  use consts_coms, only: r8
  use oname_coms,  only: nl

  implicit none

! Arguments

  integer, intent(in) :: iw
  real,    intent(in) :: dtl

! Local variables

  integer :: k, kt, kb, ka, kk, nlev, kc, ns, n
  integer :: num_mix, nsrc

  real(r8) :: bi(mza), ci(mza), ei(mza)
  real(r8) :: di(mza,num_cumix+2), ui(mza,num_cumix+2)
  real(r8) :: massflx_acm(mza), massflx_loc(mza), massflx_tot(mza)
  real(r8) :: massflx

  real :: dtom(mza)
  real :: massflx0, massflxk
  real :: ut, vt, dti, fsum, frac, fact

  real, parameter :: efv = 1.0  ! mixing efficiency for momentum
  real, parameter :: fl  = 0.5  ! extra local mixing

  real,     allocatable :: fracsrc(:)
  real,     allocatable :: depthi(:)
  real(r8), allocatable :: ai(:,:)
  real(r8), allocatable :: temp(:,:)

  ! First make sure that there was convection

  if ( kddtop(iw) < lpw(iw) ) return

  kt = min(kddtop(iw) + 1, mza)
  kb = kddmax(iw)
  ka = kddbot(iw)

  dti = 1.0 / dtl

  nlev = kt - ka + 1
  nsrc = kt - kb + 1

  if (nl%conv_uv_mix == 1) then
     num_mix = num_cumix + 2
  else
     num_mix = num_cumix
  endif

  allocate(fracsrc (nsrc))
  allocate(depthi  (nsrc))
  allocate(ai  (mza,nsrc))
  allocate(temp(mza,nsrc))

  ! Compute ACM missing rates

  massflx0 = cddf(iw) / real( rho(kb,iw) )

  if (nsrc > 1) then

     fsum = 0.0
     do kk = 1, nsrc
        k  = kk + ka - 1
        depthi(kk) = 1.0 / (zm(kt-kk) - zm(ka-1))
        fsum = fsum + arw(k,iw) * depthi(kk)
     enddo

     frac = arw(kb,iw) / ( (zm(kb-1) - zm(ka-1)) * fsum )

  else

     depthi(1) = 1.0 / (zm(kt-1) - zm(ka-1))
     frac      = 1.0

  endif

  temp        = 0.0_r8
  massflx_loc = 0.0_r8

  massflxk = frac * massflx0

  ! assumes arw always decreases or stays the same going down
  do kk = 1, nsrc
     fact = massflxk * depthi(kk)
     do k = kt-kk, ka, -1
        temp(k,kk) = rho(k,iw) * fact * arw(k,iw) * (zm(k) - zm(ka-1))
     enddo
  enddo

  do k = ka, kt-1
     massflx_acm(k) = sum(temp(k,1:nsrc))
  enddo

  do k = kb-1, ka, -1
     massflx        = rho(k,iw) * massflx0 * arw(k,iw)
     massflx_loc(k) = min(massflx - massflx_acm(k), fl * massflx_acm(k))
  enddo

  do k = kt-1, ka, -1
     massflx_loc(k) = max( massflx_loc(k), 0.05_r8 * massflx_acm(k) )
     massflx_tot(k) = massflx_acm(k) + massflx_loc(k)
  enddo

  massflx_acm(kt)   = 0._r8
  massflx_acm(ka-1) = 0._r8

  massflx_tot(kt)   = 0._r8
  massflx_tot(ka-1) = 0._r8

! Define arrays A,B,E which make up matrix and D which is RHS

  ! Load scalars into RHS array

  do k = ka, kt
     dtom(k) = dtl * volti(k,iw) / real(rho(k,iw))
  enddo

  !dir$ nounroll_and_jam
  do ns = 1, num_cumix
     n = cumix_map(ns)
     do kc = 1, nlev
        k  = kt - kc + 1
        di(kc,ns) = scalar_tab(n)%var_p(k,iw) + dtl * scalar_tab(n)%var_t(k,iw) / real(rho(k,iw))
     enddo
  enddo

  ! Load THIL and earth-Cartesian velocities into RHS array

  do kc = 1, nlev
     k  = kt - kc + 1
     di(kc,num_cumix+1) = ue(k,iw)
     di(kc,num_cumix+2) = ve(k,iw)
  enddo

  ! Loop over T levels
  do kc = 1, nlev
     k  = kt - kc + 1
     bi(kc)   = 1._r8 + dtom(k) * (massflx_tot(k) + massflx_loc(k-1))
     ei(kc)   =       - dtom(k) *  massflx_tot(k-1)
     ci(kc)   =       - dtom(k) *  massflx_loc(k)
  enddo

  do kk = 1, nsrc
     do kc = kk+1, nlev
        k  = kt - kc + 1
        ai(kc,kk) = dtom(k) * (temp(k-1,kk) - temp(k,kk))
     enddo
  enddo

  do kc = 1, nsrc
     k  = kt - kc + 1
     bi(kc)   = bi(kc)   + dtom(k) * temp(k-1,kc)
     ci(kc+1) = ci(kc+1) + ai(kc+1,kc)
  enddo

  call acm_matrix( ai, bi, ci, di, ei, ui, nlev, mza, num_mix, nsrc )

 ! Now, soln contains future(t+1) values. Compute fluxes

  !dir$ nounroll_and_jam
  do ns = 1, num_cumix
     n = cumix_map(ns)
     do k = ka, kt
        kc = kt - k + 1
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + dti * real( rho(k,iw) * (ui(kc,ns) - di(kc,ns)) )
     enddo
  enddo

  if (nl%conv_uv_mix == 1) then
     do k = ka, kt
        kc = kt - k + 1

        ut = dti * real( rho(k,iw) * (ui(kc,num_cumix+1) - di(kc,num_cumix+1)) )
        vt = dti * real( rho(k,iw) * (ui(kc,num_cumix+2) - di(kc,num_cumix+2)) )

        vmxet(k,iw) = vmxet(k,iw) + ut * vxn_ew(iw) + vt * vxn_ns(iw)

        vmyet(k,iw) = vmyet(k,iw) + ut * vyn_ew(iw) + vt * vyn_ns(iw)

        vmzet(k,iw) = vmzet(k,iw)                   + vt * vzn_ns(iw)
     enddo
  endif

end subroutine acmcldmix_dd
