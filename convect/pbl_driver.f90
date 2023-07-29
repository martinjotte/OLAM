module pbl_drivers

  real, parameter, private :: forw_imp = 1.5

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

     ! Apply explicit surface fluxes and emissions for scalars

     call apply_surface_fluxes( iw )

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

  use mem_grid,      only: lsw, lpw, mwa, arw, zfacim2
  use mem_ijtabs,    only: jtab_w, jtw_prog, itab_w
  use mem_turb,      only: frac_urb, frac_land, frac_sea, frac_lake, &
                           frac_sfc, frac_sfck, ustar, wstar, wtv0, moli
  use mem_sfcg,       only: sfcg, itab_wsfc
  use misc_coms,     only: runtype, iparallel
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: isfcl
  use consts_coms,   only: eps_virt, alvlocp
  use obnd,          only: lbcopy_w
  use olam_mpi_atm,  only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j, iw, iwsfc, k, km, ks, jsfc, jasfc

! Populate urban and land, sea, and lake fraction arrays

  frac_urb (:) = 0.0
  frac_land(:) = 0.0
  frac_sea (:) = 0.0
  frac_lake(:) = 0.0

  if (isfcl > 0) then

     ! Loop over all ATM grid cells

     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Loop over surface cells that couple to current ATM cell

        do jsfc = 1,itab_w(iw)%jsfc2
           iwsfc = itab_w(iw)%iwsfc(jsfc)
           jasfc = itab_w(iw)%jasfc(jsfc)

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
                        r1dvara3=frac_sea, r1dvara4=frac_lake  )
        call mpi_recv_w(r1dvara1=frac_urb, r1dvara2=frac_land, &
                        r1dvara3=frac_sea, r1dvara4=frac_lake  )
     endif

     call lbcopy_w(v1=frac_urb, v2=frac_land, v3=frac_sea, v4=frac_lake)

  endif

  !$omp parallel do private(jsfc,iwsfc,jasfc,ks,k,km)
  do iw = 2, mwa

     frac_sfc (:,iw) = 0.
     frac_sfck(:,iw) = 0.

     if (lsw(iw) == 1) then

        frac_sfc (1,iw) = 1.0
        frac_sfck(1,iw) = 1.0

     else

        ! Store the fraction of the total surface that intersects with each layer (a+ - a-) / a0

        do jsfc = 1, itab_w(iw)%jsfc2
           iwsfc = itab_w(iw)%iwsfc(jsfc)
           jasfc = itab_w(iw)%jasfc(jsfc)
           ks    =  itab_wsfc(iwsfc)%kwatm(jasfc) - lpw(iw) + 1

           frac_sfc(ks,iw) = frac_sfc(ks,iw) + itab_wsfc(iwsfc)%arcoariw(jasfc)
        enddo

        ! Store the fraction of the current layer that intersects with the surface (a+ - a-) / a+

        frac_sfck(1,iw) = 1.0

        do ks = 2, lsw(iw)
           k  = lpw(iw) + ks - 1
           km = k - 1
           frac_sfck(ks,iw) = &
                max(arw(k,iw) * zfacim2(k) - arw(km,iw) * zfacim2(km),0.) / (arw(k,iw) * zfacim2(k))
        enddo

     endif

  enddo
  !$omp end parallel do

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

end subroutine pbl_init

!===============================================================================

subroutine apply_surface_fluxes( iw )

  use mem_grid,   only: mza, lpw, lsw, volti
  use misc_coms,  only: dtlm
  use mem_basic,  only: rho
  use var_tables, only: sxfer_map, num_sxfer, emis_map, num_emis, scalar_tab
  use mem_tend,   only: thilt
  use mem_turb,   only: sxfer_tk

  implicit none

  integer, intent(in) :: iw
  integer             :: ns, n, ks, k
  real                :: dtli

  dtli = 1.0 / real(dtlm)

  ! Apply surface heat flux directly to tendency arrays

  !dir$ ivdep
  do ks = 1, lsw(iw)
     k = lpw(iw) + ks - 1
     thilt(k,iw) = thilt(k,iw) + dtli * volti(k,iw) * sxfer_tk(ks,iw)
  enddo

  ! Apply surface moisture and scalar fluxes directly to tendency arrays

  do ns = 1, num_sxfer
     n = sxfer_map(ns)

     !dir$ ivdep
     do ks = 1, lsw(iw)
        k = lpw(iw) + ks - 1
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + dtli * volti(k,iw) * scalar_tab(n)%sxfer(ks,iw)
     enddo
  enddo

  ! Apply emissions directly to the tendency arrays
  ! (emis units are concentration / sec )

  do ns = 1, num_emis
     n = emis_map(ns)

     !dir$ ivdep
     do k = lpw(iw), mza-1
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + real(rho(k,iw)) * scalar_tab(n)%emis(k,iw)
     enddo
  enddo

end subroutine apply_surface_fluxes

!===============================================================================

subroutine solve_eddy_diff_scalars( iw )

  use mem_grid,   only: mza, arw, volti, lpw, lsw, dzim, arw0i, volt
  use mem_turb,   only: kpblh, vkh, sfluxr, agamma, khtop
  use mem_basic,  only: rho
  use mem_ijtabs, only: itab_w
  use misc_coms,  only: dtlm
  use tridiag,    only: tridv
  use var_tables, only: num_pblmix, pblmix_map, scalar_tab
  use oname_coms, only: nl

  use supercell_testm, only: rr_w_init

  implicit none

  integer, intent(in) :: iw

  real    :: dtomass(mza), dtorho(mza), akodz(mza)
  integer :: ka, k, n, ns, kmax
  real    :: vctr5(mza), vctr6(mza), vctr7(mza)
  real    :: soln(mza,num_pblmix), rhs(mza,num_pblmix), varp(mza)
  real    :: dtl, dtli, wc0

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
     dtomass(k) = forw_imp * dtorho(k) * volti(k,iw)
  enddo

  ! Vertical loop over W levels
  do k = ka, khtop(iw)
     akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
     vctr5(k) = -akodz(k) * dtomass(k)
     vctr7(k) = -akodz(k) * dtomass(k+1)
     vctr6(k) = 1. - vctr5(k) - vctr7(k)
  enddo

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

     if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
        if (agamma(ka,iw) / arw(ka,iw) > 1.e-6) then

           wc0 = 0.0

           if (ns == 1) then

              ! water vapor
              wc0 = sfluxr(iw) / real(rho(ka,iw))

           elseif (scalar_tab(n)%do_sxfer) then

              ! other scalars with surface flux
              wc0 = sum( scalar_tab(n)%sxfer(1:lsw(iw),iw) / real(rho(ka:ka+lsw(iw)-1,iw)) ) &
                    * dtli * arw0i(iw)

           endif

           if (scalar_tab(n)%do_emis) then

              ! other scalars with near-surface emissions
              wc0 = wc0 + sum( real(volt(ka:ka+2,iw)) * scalar_tab(n)%emis(ka:ka+2,iw) ) &
                        / arw(ka+2,iw)

           endif

           if (abs(wc0) > 1.e-25) then

              do k = ka, kpblh(iw)
                 rhs(k,ns) = rhs(k,ns) + wc0 * agamma(k,iw)
              enddo

           endif

        endif
     endif

  enddo

  ! Solve tri-diagonal matrix equation for scalars

  call tridv(vctr5, vctr6, vctr7, rhs, soln, ka, khtop(iw), mza, num_pblmix)

  do ns = 1, num_pblmix
     n = pblmix_map(ns)

     ! Set bottom and top vertical internal turbulent fluxes to zero

     soln(ka-1,ns) = 0.
     soln(kmax,ns) = 0.

     ! Soln contains future(t+1) fluxes. Compute scalar tendencies

     do k = ka, kmax
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + volti(k,iw) * (soln(k-1,ns) - soln(k,ns))
     enddo

  enddo

end subroutine solve_eddy_diff_scalars

!===============================================================================

subroutine solve_eddy_diff_heat(iw, thilt)

  use mem_turb,        only: vkh, kpblh, sfluxt, agamma, khtop
  use mem_basic,       only: rho, thil
  use misc_coms,       only: dtsm
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
     dtomass(k) = forw_imp * dtom
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
     if (agamma(ka,iw) > 1.e-7) then

        wt0 = sfluxt(iw) / real(rho(ka,iw))
        do k = ka, kpblh(iw)
           rhs(k) = rhs(k) + wt0 * agamma(k,iw)
        enddo

     endif
  endif

  ! Solve tri-diagonal matrix equation

  call tridiffo(mza, ka, khtop(iw), vctr5, vctr6, vctr7, rhs, soln)

  ! Set bottom and top vertical internal turbulent fluxes to zero

  soln(ka-1) = 0.0
  soln(kmax) = 0.0

  ! Compute temperature tendency (weighted by volumne)

  do k = ka, kmax
     thilt(k) = thilt(k) + soln(k-1) - soln(k)
  enddo

end subroutine solve_eddy_diff_heat

!===============================================================================

subroutine solve_eddy_diff_vc( iv, vmt )

  use mem_grid,    only: mza, lpv, lpw, lsw, dzit_bot, volvi, arw, volt, dzim, &
                         unx, uny, unz
  use mem_basic,   only: rho, vc, vxe, vye, vze
  use mem_turb,    only: akm_sfc, vkm, ustar, kmtop
  use misc_coms,   only: dtsm
  use mem_ijtabs,  only: itab_v
  use tridiag,     only: tridiffo
  use consts_coms, only: vonk
  use oname_coms,  only: nl

  implicit none

  integer, intent(in)    :: iv
  real,    intent(inout) :: vmt(mza)

  integer :: k, ks, iw1, iw2, kb, ka, kc, kd, kmax
  real    :: mass, aksodz, dts, uc, wspdv2, wspdw2

  real :: akodz(mza), dtomass(mza), soln(mza), vctr2(mza)
  real :: vctr3(mza), vctr5(mza), vctr6(mza), vctr7(mza)

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

  dts = dtsm
  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  ka = min(lpw(iw1),lpw(iw2))
  kb = lpv(iv)
  kc = max(lpw(iw1),lpw(iw2))
  kd = max(lpw(iw1)+lsw(iw1), lpw(iw2)+lsw(iw2)) - 1

  kmax = max(kmtop(iw1)+1, kmtop(iw2)+1, kb, kd)

! do k = ka, mza-1
  do k = ka, kmax-1
     akodz(k) = (arw(k,iw1) * vkm(k,iw1) + arw(k,iw2) * vkm(k,iw2)) * dzim(k)
  enddo

! special for underground cells: set minimum diffusivity

  do k = ka, kb-1
     akodz(k) = max(akodz(k), vonk * ( arw(k,iw1) * ustar(iw1) * real(rho(k,iw1)) &
                                     + arw(k,iw2) * ustar(iw2) * real(rho(k,iw2)) ))
  enddo

  akodz(ka-1) = 0.0
! akodz(mza)  = 0.0
  akodz(kmax) = 0.0

  ! Vertical loop over T levels

! do k = ka, mza
  do k = ka, kmax
     mass       = volt(k,iw1) * rho(k,iw1) + volt(k,iw2) * rho(k,iw2)
     dtomass(k) = forw_imp * dts / mass
     vctr5  (k) = -dtomass(k) * akodz(k-1)
     vctr7  (k) = -dtomass(k) * akodz(k)
     vctr6  (k) = 1. - vctr5(k) - vctr7(k)
  enddo

  ! Surface drag over multiple levels

  vctr3(ka:kmax) = 0.

  do k = lpw(iw1), lpw(iw1) + lsw(iw1) - 1
     ks = k - lpw(iw1) + 1

     uc = unx(iv) * vxe(k,iw1) &
        + uny(iv) * vye(k,iw1) &
        + unz(iv) * vze(k,iw1)

     wspdv2 = vc(k,iv)**2 + uc**2 + 1.e-8
     wspdw2 = vxe(k,iw1)**2 + vye(k,iw1)**2 + vze(k,iw1)**2 + 1.e-8

     aksodz   = akm_sfc(ks,iw1) * dzit_bot(k) * sqrt(max(1., wspdv2 / wspdw2))
     vctr3(k) = vctr3(k) + aksodz
     vctr6(k) = vctr6(k) + dtomass(k) * aksodz
  enddo

  do k = lpw(iw2), lpw(iw2) + lsw(iw2) - 1
     ks = k - lpw(iw2) + 1

     uc = unx(iv) * vxe(k,iw2) &
        + uny(iv) * vye(k,iw2) &
        + unz(iv) * vze(k,iw2)

     wspdv2 = vc(k,iv)**2 + uc**2 + 1.e-8
     wspdw2 = vxe(k,iw2)**2 + vye(k,iw2)**2 + vze(k,iw2)**2 + 1.e-8

     aksodz   = akm_sfc(ks,iw2) * dzit_bot(k) * sqrt(max(1., wspdv2 / wspdw2))
     vctr3(k) = vctr3(k) + aksodz
     vctr6(k) = vctr6(k) + dtomass(k) * aksodz
  enddo

  ! Solve tri-diagonal matrix equation

  if (ka <= kmax) then
     call tridiffo(mza, ka, kmax, vctr5, vctr6, vctr7, vc(:,iv), soln)
  endif

  ! Now, soln contains (t+1) values.
  ! Compute internal vertical turbulent fluxes

! do k = ka, mza-1
  do k = ka, kmax-1
     vctr2(k) = akodz(k) * (soln(k) - soln(k+1))
  enddo

  ! Set bottom and top internal fluxes to zero

  vctr2(ka-1) = 0.
! vctr2(mza ) = 0.
  vctr2(kmax) = 0.

  ! Include surface drag

  do k = ka, kd
     vmt(k) = vmt(k) + volvi(k,iv) * &
            (vctr2(k-1) - vctr2(k) - vctr3(k) * soln(k))
  enddo

  ! Above surface

! do k = kd+1, mza
  do k = kd+1, kmax
     vmt(k) = vmt(k) + volvi(k,iv) * (vctr2(k-1) - vctr2(k))
  enddo

end subroutine solve_eddy_diff_vc


!===============================================================================

subroutine solve_eddy_diff_vxe( iw, vmxet, vmyet, vmzet )

  use mem_grid,    only: mza, lpw, lsw, dzit_bot, volti, arw, dzim, dzit_bot
  use mem_basic,   only: rho, vxe, vye, vze
  use mem_turb,    only: akm_sfc, vkm, kmtop
  use misc_coms,   only: dtsm
  use tridiag,     only: tridv
  use oname_coms,  only: nl

  implicit none

  integer, intent(in)    :: iw
  real,    intent(inout) :: vmxet(mza), vmyet(mza), vmzet(mza)

  integer :: k, ka, ks, kmax
  real    :: dtom, dts

  real :: akodz(mza), dtomass(mza), soln(mza,3), rhs(mza,3), vctr2(mza,3)
  real :: vctr3(mza), vctr5(mza), vctr6(mza), vctr7(mza)

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

  ka  = lpw(iw)
  dts = dtsm

  kmax = max(kmtop(iw)+1, ka + lsw(iw) - 1)

  akodz(ka-1)             = 0.0
  akodz(kmtop(iw)+1:kmax) = 0.0
! akodz(mza)  = 0.0

  ! Vertical loop over W levels

! do k = ka, mza-1
  do k = ka, kmtop(iw)
     akodz(k) = arw(k,iw) * vkm(k,iw) * dzim(k)
  enddo

  ! Vertical loop over T levels

! do k = ka, mza
  do k = ka, kmax
     dtom         = dts * volti(k,iw) / rho(k,iw)
     dtomass(k)   = forw_imp * dtom
     vctr5  (k)   = -dtomass(k) * akodz(k-1)
     vctr7  (k)   = -dtomass(k) * akodz(k)
     vctr6  (k)   = 1. - vctr5(k) - vctr7(k)
     rhs    (k,1) = vxe(k,iw)
     rhs    (k,2) = vye(k,iw)
     rhs    (k,3) = vze(k,iw)
  enddo

  ! Surface drag over multiple levels

  do k = ka, ka + lsw(iw) - 1
     ks = k - ka + 1

     vctr3(k) = akm_sfc(ks,iw) * dzit_bot(k)
     vctr6(k) = vctr6(k) + dtomass(k) * vctr3(k)
  enddo

  ! Solve tri-diagonal matrix equation

  if (ka <= kmax) then
     call tridv(vctr5, vctr6, vctr7, rhs, soln, ka, kmax, mza, 3)
  endif

  ! Now, soln contains (t+1) values.
  ! Compute internal vertical turbulent fluxes

! do k = ka, mza-1
  do k = ka, kmax-1
     vctr2(k,1) = akodz(k) * (soln(k,1) - soln(k+1,1))
     vctr2(k,2) = akodz(k) * (soln(k,2) - soln(k+1,2))
     vctr2(k,3) = akodz(k) * (soln(k,3) - soln(k+1,3))
  enddo

  ! Set bottom and top internal fluxes to zero

  vctr2(ka-1,:) = 0.
  vctr2(kmax,:) = 0.
! vctr2(mza ,:) = 0.

  ! Surface drag and turbulent fluxes at surface (scaled by volt)

  do k = ka, ka + lsw(iw) - 1
     vmxet(k) = vmxet(k) + (vctr2(k-1,1) - vctr2(k,1) - vctr3(k) * soln(k,1))

     vmyet(k) = vmyet(k) + (vctr2(k-1,2) - vctr2(k,2) - vctr3(k) * soln(k,2))

     vmzet(k) = vmzet(k) + (vctr2(k-1,3) - vctr2(k,3) - vctr3(k) * soln(k,3))
  enddo

  ! Just turbulent fluxes above surface (scaled by volt)

! do k = kd+1, mza
  do k = ka + lsw(iw), kmax
     vmxet(k) = vmxet(k) + vctr2(k-1,1) - vctr2(k,1)

     vmyet(k) = vmyet(k) + vctr2(k-1,2) - vctr2(k,2)

     vmzet(k) = vmzet(k) + vctr2(k-1,3) - vctr2(k,3)
  enddo

end subroutine solve_eddy_diff_vxe


end module pbl_drivers
