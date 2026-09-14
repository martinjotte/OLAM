module pbl_drivers

  ! Implictness weightings for tracers and momentum
  real, parameter, private :: fcnt = 1.55
  real, parameter, private :: fcnm = 1.55

  ! This forward-weights the implicit fluxes to damp 2-DT and 2-DZ
  ! oscillations that result when the eddy diffusivities are functions
  ! of T, Q, and shear through the Richardson number, which
  ! makes the diffusion equation nonlinear. See Hammarstrand, MWR, 1997.
  ! ECMWF uses a implicitness weighting of 1.5 to avoid this instability.

contains

!===============================================================================

subroutine pbl_driver()

  use mem_grid,       only: mza, lpw
  use misc_coms,      only: idiffk, iparallel
  use mem_turb,       only: vkm, vkh, kmtop, khtop, kpblh, akm_sfc, ustar, pblh
  use consts_coms,    only: grav, vonk, eps_virt, alvlocp, r8
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog
  use module_bl_acm2, only: acm2_pblhgt, acm2_eddyx
  use smagorinsky,    only: turb_k
  use oname_coms,     only: nl
  use olam_mpi_atm,   only: mpi_send_w, mpi_recv_w
  use obnd,           only: lbcopy_w

  implicit none

  integer :: j, k, ka, iw, mrlw

! Loop over all W/T points where PBL parameterization may be done

!----------------------------------------------------------------------
  !$omp parallel do private(iw,mrlw,ka,k)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

     ! MRL for current IW column

     mrlw = itab_w(iw)%mrlw
     ka   = lpw(iw)

     ! Diagnose PBL height regardless of scheme

     call acm2_pblhgt(iw)

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 1) then

        ! non-local convective tranport scheme

        call acm2_eddyx(iw)

     else if (idiffk(mrlw) == 2 .or. idiffk(mrlw) == 3) then

        ! Smagorinsky scheme

        call turb_k(iw, mrlw)

     else

        vkm(:,iw) = 0.0
        vkh(:,iw) = 0.0

     endif

     ! For DCMIP 2015 baroclinic and tropical cyclone tests, return here so that
     ! vertical mixing is not done by standard OLAM code; instead it will
     ! (optionally) be done in dcmip_physics routine.

     if (nl%test_case == 110 .or. &
         nl%test_case == 111 .or. &
         nl%test_case == 112 .or. &
         nl%test_case == 113 .or. &
         nl%test_case == 114 .or. &
         nl%test_case == 121 .or. &
         nl%test_case == 122) cycle

     if (idiffk(mrlw) == 0) then

        kmtop(iw) = 1
        khtop(iw) = 1

     else

        do k = mza-1, ka, -1
           if (vkm(k,iw) > 1.e-20) exit
        enddo
        if (k < ka) then
           kmtop(iw) = 1
        else
           kmtop(iw) = k
        endif

        do k = mza-1, ka, -1
           if (vkh(k,iw) > 1.e-20) exit
        enddo
        if (k < ka) then
           khtop(iw) = 1
        else
           khtop(iw) = k
        endif

     endif

     ! Compute tendency of scalars due to vertical turbulent transport

     if (khtop(iw) >= ka) then
        call solve_eddy_diff_scalars(iw)
     endif

  enddo
  !$omp end parallel do

  ! MPI send/recv of computed K's

  if (iparallel == 1) then

     call mpi_send_w(rvara1=vkm, rvara2=vkh, swvar1=akm_sfc, &
                     r1dvara1=ustar, r1dvara2=pblh, &
                     i1dvara1=khtop, i1dvara2=kmtop, i1dvara3=kpblh)

     call mpi_recv_w(rvara1=vkm, rvara2=vkh, swvar1=akm_sfc, &
                     r1dvara1=ustar, r1dvara2=pblh, &
                     i1dvara1=khtop, i1dvara2=kmtop, i1dvara3=kpblh)
  endif

 ! Lateral boundary copy of computed K's

  call lbcopy_w(a1=vkm, a2=vkh, s1=akm_sfc, v1=ustar, v2=pblh, &
                iv1=khtop, iv2=kmtop, iv3=kpblh)

end subroutine pbl_driver

!===============================================================================

subroutine comp_horiz_k()

  use mem_ijtabs, only: jtab_v, jtv_wadj
  implicit none

  integer             :: j, iv

  !$omp parallel do private(iv)
  do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

     call comp_horiz_k_column(iv)

  enddo
  !$omp end parallel do

end subroutine comp_horiz_k

!===============================================================================

subroutine comp_horiz_k_column(iv)

  use mem_grid,    only: arw0, lpv, mza, arv, volt, dniv
  use misc_coms,   only: dtlm, akmin
  use mem_ijtabs,  only: itab_v, itab_w
  use mem_turb,    only: vkm, vkh, akmodx, akhodx, kmtop, khtop, kmtopv, khtopv
  use mem_basic,   only: rho

  implicit none

  integer, intent(in) :: iv
  integer             :: iw1, iw2, k, mrl
  real                :: bkmin, dtl
  real                :: hkm, hkh, stab1, stab2
  real                :: fact1, fact2
  real                :: hk0(mza), akstab(mza)

  ! Compute ARV * K / DX, and ensure that horizontal
  ! diffusion is stable over the long timestep

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)
  mrl = itab_v(iv)%mrlv

  dtl = dtlm

  fact1 = 0.95 / (dtl * itab_w(iw1)%npoly )
  fact2 = 0.95 / (dtl * itab_w(iw2)%npoly )

  if (akmin(mrl) > 1.e-7) then

     bkmin = akmin(mrl) * .075 * min(arw0(iw1),arw0(iw2)) ** .66666666
     khtopv(iv) = mza
     kmtopv(iv) = mza

     do k = lpv(iv), mza
        hk0(k) = 0.5 * real(rho(k,iw1) + rho(k,iw2)) * bkmin
     enddo

  else

     khtopv(iv) = max(khtop(iw1), khtop(iw2))
     if (khtopv(iv) >= lpv(iv)) then
        khtopv(iv) = khtopv(iv) + 1
     else
        khtopv(iv) = 1
     endif

     kmtopv(iv) = max(kmtop(iw1), kmtop(iw2))
     if (kmtopv(iv) >= lpv(iv)) then
        kmtopv(iv) = kmtopv(iv) + 1
     else
        kmtopv(iv) = 1
     endif

  endif

  do k = lpv(iv), max(kmtopv(iv), khtopv(iv))
     stab1 = rho(k,iw1) * real(volt(k,iw1)) * fact1
     stab2 = rho(k,iw2) * real(volt(k,iw2)) * fact2
     akstab(k) = min(stab1,stab2)
  enddo

  do k = lpv(iv), kmtopv(iv)
     hkm = 0.25 * (vkm(k,iw1) + vkm(k,iw2) + vkm(k-1,iw1) + vkm(k-1,iw2))

     if (akmin(mrl) > 1.e-7) then
        hkm = max(hkm, hk0(k))
     endif

     akmodx(k,iv) = min( dniv(iv) * arv(k,iv) * hkm, akstab(k) )
  enddo

  do k = lpv(iv), khtopv(iv)
     hkh = 0.25 * (vkh(k,iw1) + vkh(k,iw2) + vkh(k-1,iw1) + vkh(k-1,iw2))

     if (akmin(mrl) > 1.e-7) then
        hkh = max(hkh, hk0(k))
     endif

     akhodx(k,iv) = min( dniv(iv) * arv(k,iv) * hkh, akstab(k) )
  enddo

end subroutine comp_horiz_k_column

!===============================================================================

subroutine pbl_init()

  use mem_grid,     only: lsw, lpw, mwa, arw, zfacim2, zfacm2, arw0
  use mem_ijtabs,   only: jtab_w, jtw_prog, itab_w
  use mem_turb,     only: frac_urb, frac_land, frac_sea, frac_lake, arw_sfc,  &
                          frac_sfc, frac_sfck, frac_norm, ustar, wstar, wtv0, &
                          moli, pblh
  use mem_sfcg,      only: sfcg, itab_wsfc
  use misc_coms,     only: runtype, iparallel, mdomain
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: isfcl
  use consts_coms,   only: eps_virt, alvlocp
  use obnd,          only: lbcopy_w
  use olam_mpi_atm,  only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j, iw, iwsfc, k, ka, ks, kw, jsfc, jasfc

! Populate urban and land, sea, and lake fraction arrays

  if (isfcl > 0) then

     ! Loop over all primary ATM grid cells

     !$omp parallel do private(iw,jsfc,iwsfc,jasfc,ks)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        frac_urb  (iw) = 0.
        frac_land (iw) = 0.
        frac_sea  (iw) = 0.
        frac_lake (iw) = 0.
        arw_sfc (:,iw) = 0.
        frac_sfc(:,iw) = 0.

        ! Loop over surface cells that couple to current ATM cell

        do jsfc = 1,itab_w(iw)%jsfc2
           iwsfc = itab_w(iw)%iwsfc(jsfc)
           jasfc = itab_w(iw)%jasfc(jsfc)
           ks    = itab_wsfc(iwsfc)%kwatm(jasfc) - lpw(iw) + 1

           ! Store the total surface coupling area at each height, and the
           ! total surface coupling area divided by the atm cell area at each height
           arw_sfc (ks,iw) = arw_sfc (ks,iw) + itab_wsfc(iwsfc)%arc     (jasfc)
           frac_sfc(ks,iw) = frac_sfc(ks,iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)

           if (sfcg%leaf_class(iwsfc) == 0) then
              frac_sea(iw) = frac_sea(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)
           elseif (sfcg%leaf_class(iwsfc) == 1) then
              frac_lake(iw) = frac_lake(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)
           else
              frac_land(iw) = frac_land(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)

              if (any(sfcg%leaf_class(iwsfc) == (/ 19, 21 /))) then
                 frac_urb(iw) = frac_urb(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)
              endif
           endif

        enddo
     enddo

     if (iparallel == 1) then
        call mpi_send_w(r1dvara1=frac_urb, r1dvara2=frac_land, &
                        r1dvara3=frac_sea, r1dvara4=frac_lake, &
                        swvar1=arw_sfc,    swvar2=frac_sfc)
        call mpi_recv_w(r1dvara1=frac_urb, r1dvara2=frac_land, &
                        r1dvara3=frac_sea, r1dvara4=frac_lake, &
                        swvar1=arw_sfc,    swvar2=frac_sfc)
     endif

     ! Currently the surface model only uses a global domain
     ! call lbcopy_w(v1=frac_urb, v2=frac_land, v3=frac_sea, v4=frac_lake, &
     !               s1=arw_sfc,  s2=frac_sfc)

     !$omp parallel do private(ks,k)
     do iw = 2, mwa

        frac_sfck(:,iw) = 0.
        frac_norm(:,iw) = 0.

        if (lsw(iw) == 1) then

           frac_sfck(1,iw) = 1.0
           frac_norm(1,iw) = 1.0

        else

           ! Store the fraction of the current layer that intersects with the surface

           frac_sfck(1,iw) = 1.0

           do ks = 2, lsw(iw)
              k  = lpw(iw) + ks - 1
              frac_sfck(ks,iw) = arw_sfc(ks,iw) / (arw(k,iw) * zfacim2(k) * zfacm2(k-1))
           enddo

           ! Frac_norm is used to normalize the surface overlaps below lsw

           do ks = 1, lsw(iw)
              frac_norm(ks,iw) = 1.0 / sum( frac_sfc(1:ks,iw) )
           enddo

        endif

     enddo
     !$omp end parallel do

  else  ! isfcl == 0

     !$omp parallel do private(ka,ks,kw)
     do iw = 2, mwa
        ka = lpw(iw)

        frac_urb  (iw) = 0.
        frac_land (iw) = 1.
        frac_sea  (iw) = 0.
        frac_lake (iw) = 0.

        if ( lsw(iw) == 1 ) then

           frac_sfc (1,iw) = 1.
           frac_sfck(1,iw) = 1.
           frac_norm(1,iw) = 1.

           if (mdomain > 1) then
              arw_sfc(1,iw) = arw(ka,iw) * zfacim2(ka) * zfacm2(ka-1)
           else
              arw_sfc(1,iw) = arw(ka,iw)
           endif

        else

           do ks = 1, lsw(iw)
              kw = ka + ks - 1

              if (mdomain > 1) then
                 arw_sfc  (ks,iw) = arw(kw,iw) * zfacim2(kw) * zfacm2(kw-1) - arw(kw-1,iw)
                 frac_sfc (ks,iw) = arw_sfc(ks,iw) * zfacim2(kw-1) / arw0(iw)
                 frac_sfck(ks,iw) = arw_sfc(ks,iw) / (arw(kw,iw) * zfacim2(kw) * zfacm2(kw-1))
              else
                 arw_sfc  (ks,iw) = arw(kw,iw) - arw(kw-1,iw)
                 frac_sfc (ks,iw) = arw_sfc(ks,iw) / arw0(iw)
                 frac_sfck(ks,iw) = arw_sfc(ks,iw) / arw(kw,iw)
              endif

           enddo

           do ks = 1, lsw(iw)
              frac_norm(ks,iw) = 1.0 / sum( frac_sfc(1:ks,iw) )
           enddo

        endif

     enddo
     !$omp end parallel do

  endif  ! isfcl

  if (runtype /= 'INITIAL') return

! Initialize PBL height and some PBL quantities

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ustar(iw) = 0.2
     wstar(iw) = 0.0
     wtv0 (iw) = 0.0
     moli (iw) = 0.0

     call acm2_pblhgt( iw )

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_w(r1dvara1=pblh)
     call mpi_recv_w(r1dvara1=pblh)
  endif

  call lbcopy_w(v1=pblh)

end subroutine pbl_init

!===============================================================================

subroutine solve_eddy_diff_scalars( iw )

  use mem_grid,   only: mza, arw, volti, lpw, lsw, dzim, volt
  use mem_turb,   only: kpblh, vkh, sfluxr, gam_rkh, khtop
  use mem_basic,  only: rho
  use mem_ijtabs, only: itab_w
  use misc_coms,  only: dtlm
  use tridiag,    only: tridv
  use var_tables, only: num_pblmix, pblmix_map, scalar_tab
  use oname_coms, only: nl
  use consts_coms,only: r8

  use supercell_testm, only: rr_w_init

  implicit none

  integer, intent(in) :: iw

  real    :: dtomass(mza), dtorho(mza), akodz(mza), agamma(mza)
  integer :: ka, k, n, ns, kmax, ks
  real    :: vctr5(mza), vctr6(mza), vctr7(mza)
  real    :: soln(mza,num_pblmix), rhs(mza,num_pblmix), varp(mza)
  real    :: dtl, dtli, wc0
  real(r8):: flux(mza)
  logical :: do_nonloc

  ! For DCMIP 2015 baroclinic and tropical cyclone tests, return here so that
  ! vertical mixing is not done by standard OLAM code; instead it will
  ! (optionally) be done in dcmip_physics routine.

  if (nl%test_case == 110 .or. &
      nl%test_case == 111 .or. &
      nl%test_case == 112 .or. &
      nl%test_case == 113 .or. &
      nl%test_case == 114 .or. &
      nl%test_case == 121 .or. &
      nl%test_case == 122) return

  if (khtop(iw) < lpw(iw) .or. num_pblmix == 0) return

  ka   = lpw(iw)
  dtl  = dtlm
  dtli = 1.0 / dtl
  kmax = khtop(iw) + 1

  ! Vertical loop over T levels
  do k = ka, kmax
     dtorho (k) = dtl / real(rho(k,iw))
     dtomass(k) = fcnt * dtorho(k) * volti(k,iw)
  enddo

  ! Vertical loop over W levels
  do k = ka, khtop(iw)
     akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
     vctr5(k) = -akodz(k) * dtomass(k)
     vctr7(k) = -akodz(k) * dtomass(k+1)
     vctr6(k) = 1. - vctr5(k) - vctr7(k)
  enddo

  ! Non local PBL fluxes

  do_nonloc = .false.

  if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
     if ( gam_rkh(ka,iw) > 1.e-7 ) then

        do_nonloc = .true.

        do k = ka, kpblh(iw)-1
           agamma(k) = arw(k,iw) * gam_rkh(k,iw)
        enddo

     endif
  endif

  ! Load scalars into RHS arrays, including any nonlocal terms

  do ns = 1, num_pblmix
     n = pblmix_map(ns)

     ! Scalar variables with long-timestep forcing included

     do k = ka, kmax
        varp(k) = scalar_tab(n)%var_p(k,iw) &
                + dtorho(k) * scalar_tab(n)%var_t(k,iw)
     enddo

     ! For supercell test case, subtract initial field from rr_w

     if (nl%test_case == 131 .and. ns == 1) then
        do k = ka, kmax
           varp(k) = varp(k) - rr_w_init(k,iw)
        enddo
     endif

     ! SGS flux

     do k = ka, khtop(iw)
        rhs(k,ns) = akodz(k) * (varp(k) - varp(k+1))
     enddo

     ! Non local PBL fluxes

     if (do_nonloc) then

        if (ns == 1) then

           ! water vapor
           wc0 = sfluxr(iw) / real(rho(ka,iw))

        else

           ! estimate surface flux from near-surface tendency

           ks = lsw(iw)+1

           wc0 = sum( real(volt(ka:ka+ks,iw)) * scalar_tab(n)%var_t(ka:ka+ks,iw) &
                    / real(rho(ka:ka+ks,iw)) ) / arw(ka+ks,iw)

        endif

        if (abs(wc0) > 1.e-25) then

           do k = ka, kpblh(iw)-1
              rhs(k,ns) = rhs(k,ns) + wc0 * agamma(k)
           enddo

        endif

     endif

  enddo

  ! Solve tri-diagonal matrix equation for scalars

  call tridv(vctr5, vctr6, vctr7, rhs, soln, ka, khtop(iw), mza, num_pblmix)

  ! Set bottom and top vertical internal turbulent fluxes to zero
  flux(ka-1) = 0._r8
  flux(kmax) = 0._r8

  do ns = 1, num_pblmix
     n = pblmix_map(ns)

     ! Soln contains future(t+1) fluxes. Compute scalar tendencies

     do k = ka, kmax-1
        flux(k) = real( soln(k,ns), r8 )
     enddo

     do k = ka, kmax
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + volti(k,iw) * real( flux(k-1) - flux(k) )
     enddo

  enddo

end subroutine solve_eddy_diff_scalars

!===============================================================================

subroutine solve_eddy_diff_heat(iw, thilt)

  use mem_turb,        only: vkh, kpblh, sfluxt, gam_rkh, khtop
  use mem_basic,       only: rho, thil, theta, tair
  use misc_coms,       only: dtsm
  use consts_coms,     only: cp, r8
  use mem_grid,        only: mza, lpw, volti, dzim, arw
  use oname_coms,      only: nl
  use tridiag,         only: tridiffo
  use mem_ijtabs,      only: itab_w
  use supercell_testm, only: thil_init

  implicit none

  integer, intent(in )   :: iw
  real,    intent(inout) :: thilt(mza)  ! note: thilt assumed scaled by volume

  integer :: ka, k, kmax
  real    :: dts, dtom, akodz, wt0
  real    :: dtomass(mza), varp(mza)
  real    :: vctr5(mza), vctr6(mza), vctr7(mza), rhs(mza), soln(mza)
  real(r8):: flux(mza)

  ! For DCMIP 2015 baroclinic and tropical cyclone tests, return here so that
  ! vertical mixing is not done by standard OLAM code; instead it will
  ! (optionally) be done in dcmip_physics routine.

  if (nl%test_case == 110 .or. &
      nl%test_case == 111 .or. &
      nl%test_case == 112 .or. &
      nl%test_case == 113 .or. &
      nl%test_case == 114 .or. &
      nl%test_case == 121 .or. &
      nl%test_case == 122) return

  if (khtop(iw) < lpw(iw)) return

  ka   = lpw(iw)
  dts  = dtsm
  kmax = khtop(iw) + 1

  ! Vertical loop over T levels
  do k = ka, kmax
     dtom       = dts * volti(k,iw) / real(rho(k,iw))
     dtomass(k) = fcnt * dtom
     varp   (k) = thil(k,iw) + dtom * thilt(k)
  enddo

  ! For supercell test case, subtract initial field from thil

  if (nl%test_case == 131) then
     do k = ka, kmax
        varp(k) = varp(k) - thil_init(k,iw)
     enddo
  endif

  ! Vertical loop over W levels
  do k = ka, khtop(iw)
     akodz    = arw(k,iw) * vkh(k,iw) * dzim(k)
     rhs  (k) =  akodz * (varp(k) - varp(k+1))
     vctr5(k) = -akodz * dtomass(k)
     vctr7(k) = -akodz * dtomass(k+1)
     vctr6(k) = 1. - vctr5(k) - vctr7(k)
  enddo

  ! Non local PBL term

  if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
     if ( gam_rkh(ka,iw) > 1.e-7 ) then

        ! New calculation of wt0 for sfluxt in [W m^-2]
        wt0 = sfluxt(iw) * theta(ka,iw) / ( real(rho(ka,iw)) * cp * tair(ka,iw) )

        do k = ka, kpblh(iw)-1
           rhs(k) = rhs(k) + wt0 * arw(k,iw) * gam_rkh(k,iw)
        enddo

     endif
  endif

  ! Solve tri-diagonal matrix equation

  call tridiffo(mza, ka, khtop(iw), vctr5, vctr6, vctr7, rhs, soln)

  ! Set bottom and top vertical internal turbulent fluxes to zero

  flux(ka-1) = 0._r8
  flux(kmax) = 0._r8

  do k = ka, kmax-1
     flux(k) = real( soln(k), r8 )
  enddo

  ! Compute temperature tendency (weighted by volume)

  do k = ka, kmax
     thilt(k) = thilt(k) + real( flux(k-1) - flux(k) )
  enddo

end subroutine solve_eddy_diff_heat

!===============================================================================

subroutine solve_eddy_diff_momentum(iw, vmxet, vmyet, vmzet)

  use mem_grid,    only: mza, lpw, lsw, dzit_bot, volti, arw, dzim, arw0i, zfacim2
  use mem_basic,   only: rho, vxe, vye, vze
  use mem_turb,    only: akm_sfc, vkm, kmtop, gam_rkm, kpblh
  use misc_coms,   only: dtsm
  use tridiag,     only: tridv
  use oname_coms,  only: nl
  use mem_ijtabs,  only: itab_w

  implicit none

  integer, intent(in )   :: iw
  real,    intent(inout) :: vmxet(mza)  ! note: thilt assumed scaled by volume
  real,    intent(inout) :: vmyet(mza)  ! note: thilt assumed scaled by volume
  real,    intent(inout) :: vmzet(mza)  ! note: thilt assumed scaled by volume

  integer :: ka, k, kmax, ks
  real    :: dts, dtom, akodz, fract
  real    :: dtomass(mza)
  real    :: vctr5(mza), vctr6(mza), vctr7(mza), rhs(mza,3), soln(mza,3)
  real    :: xflxsfc(mza), yflxsfc(mza), zflxsfc(mza)
  real    :: varx(mza), vary(mza), varz(mza), cddtomass(mza)
  real    :: wux0, wuy0, wuz0, fact, agamma

  real, parameter :: vxe_sfc = 0.0  ! Placeholder for ocean current
  real, parameter :: vye_sfc = 0.0  ! Placeholder for ocean current
  real, parameter :: vze_sfc = 0.0  ! Placeholder for ocean current

  ! For DCMIP 2015 baroclinic and tropical cyclone tests, return here so that
  ! vertical mixing is not done by standard OLAM code; instead it will
  ! (optionally) be done in dcmip_physics routine.

  ! Note: this if statement needs to be checked since surface drag is now
  ! included in the implicit solver
  if (nl%test_case == 110 .or. &
      nl%test_case == 111 .or. &
      nl%test_case == 112 .or. &
      nl%test_case == 113 .or. &
      nl%test_case == 114 .or. &
      nl%test_case == 121 .or. &
      nl%test_case == 122) return

  ka   = lpw(iw)
  dts  = dtsm
  kmax = max(kmtop(iw) + 1, ka+lsw(iw)-1)

  ! Vertical loop over T levels
  do k = ka, kmax
     dtom       = dts * volti(k,iw) / real(rho(k,iw))
     dtomass(k) = fcnm * dtom
     varx   (k) = vxe(k,iw) + dtom * vmxet(k)
     vary   (k) = vye(k,iw) + dtom * vmyet(k)
     varz   (k) = vze(k,iw) + dtom * vmzet(k)
  enddo

  ! Surface momentum drag at ka-1
  k = ka-1
  akodz      = akm_sfc(1,iw) * dzit_bot(ka)
  rhs  (k,1) =  akodz * (vxe_sfc - varx(ka))
  rhs  (k,2) =  akodz * (vye_sfc - vary(ka))
  rhs  (k,3) =  akodz * (vze_sfc - varz(ka))
  vctr5(k)   =  0.0
  vctr7(k)   = -akodz * dtomass(ka)
  vctr6(k)   = 1. - vctr7(k)

  ! Special for shaved cells
  do ks = 2, lsw(iw)
     k  = ka + ks - 1

     akodz         = akm_sfc(ks,iw) * dzit_bot(k)
     cddtomass(ks) = dtomass(k) * akodz
     fract         = akodz / ( 1. + cddtomass(ks) )

     xflxsfc  (ks) = -varx(k) * fract
     yflxsfc  (ks) = -vary(k) * fract
     zflxsfc  (ks) = -varz(k) * fract

     varx      (k) = varx(k) + dtomass(k) * xflxsfc(ks)
     vary      (k) = vary(k) + dtomass(k) * yflxsfc(ks)
     varz      (k) = varz(k) + dtomass(k) * zflxsfc(ks)
     dtomass   (k) = dtomass(k) * (1. - dtomass(k) * fract)
  enddo

  ! Vertical loop over W levels
  do k = ka, kmax-1
     akodz      = arw(k,iw) * vkm(k,iw) * dzim(k)
     rhs  (k,1) =  akodz * (varx(k) - varx(k+1))
     rhs  (k,2) =  akodz * (vary(k) - vary(k+1))
     rhs  (k,3) =  akodz * (varz(k) - varz(k+1))
     vctr5(k)   = -akodz * dtomass(k)
     vctr7(k)   = -akodz * dtomass(k+1)
     vctr6(k)   = 1. - vctr5(k) - vctr7(k)
  enddo

  ! Non local PBL term

  if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
     if ( gam_rkm(ka,iw) > 1.e-7 ) then

        wux0 = 0.
        wuy0 = 0.
        wuz0 = 0.

        do ks = 1, lsw(iw)
           k  = ks + ka - 1

           fact = akm_sfc(ks,iw) * dzit_bot(k) * zfacim2(k-1) / real(rho(k,iw)) * arw0i(iw)
           wux0 = wux0 - fact * vxe(k,iw)
           wuy0 = wuy0 - fact * vye(k,iw)
           wuz0 = wuz0 - fact * vze(k,iw)
        enddo

        do k = ka, kpblh(iw)-1
           agamma   = arw(k,iw) * gam_rkm(k,iw)
           rhs(k,1) = rhs(k,1) + wux0 * agamma
           rhs(k,2) = rhs(k,2) + wuy0 * agamma
           rhs(k,3) = rhs(k,3) + wuz0 * agamma
        enddo

     endif
  endif

  ! Solve tri-diagonal matrix equation

  call tridv( vctr5, vctr6, vctr7, rhs, soln, ka-1, kmax-1, mza, 3)

  soln(kmax,:) = 0.

  ! Compute momentum tendency (weighted by volume)

  do k = ka, kmax
     vmxet(k) = vmxet(k) + soln(k-1,1) - soln(k,1)
     vmyet(k) = vmyet(k) + soln(k-1,2) - soln(k,2)
     vmzet(k) = vmzet(k) + soln(k-1,3) - soln(k,3)
  enddo

  ! Add surface drag for shaved cells

  do ks = 2, lsw(iw)
     k  = ka + ks - 1

     fract       = cddtomass(ks) / (1. + cddtomass(ks))
     xflxsfc(ks) = xflxsfc(ks) + fract * (soln(k,1) - soln(k-1,1))
     yflxsfc(ks) = yflxsfc(ks) + fract * (soln(k,2) - soln(k-1,2))
     zflxsfc(ks) = zflxsfc(ks) + fract * (soln(k,3) - soln(k-1,3))

     ! momentum tendency weighted by volume
     vmxet(k) = vmxet(k) + xflxsfc(ks)
     vmyet(k) = vmyet(k) + yflxsfc(ks)
     vmzet(k) = vmzet(k) + zflxsfc(ks)
  enddo

end subroutine solve_eddy_diff_momentum

!===============================================================================

end module pbl_drivers
