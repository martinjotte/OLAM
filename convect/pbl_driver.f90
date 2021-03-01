!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================

module pbl_drivers

  real, parameter, private :: forw_imp = 1.5

contains

!===============================================================================

subroutine pbl_driver(mrl)

  use mem_grid,       only: mza, lpw
  use misc_coms,      only: idiffk, iparallel
  use mem_turb,       only: vkm, vkh, kmtop, khtop, akm_sfc, ustar
  use consts_coms,    only: grav, vonk, eps_virt, alvlocp, r8
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog
  use module_bl_acm2, only: acm2_pblhgt, acm2_eddyx
  use smagorinsky,    only: turb_k
  use oname_coms,     only: nl
  use olam_mpi_atm,   only: mpi_send_w, mpi_recv_w
  use obnd,           only: lbcopy_w

  implicit none

  integer, intent(in) :: mrl

  integer :: j, k, ka, iw, mrlw

! Loop over all W/T points where PBL parameterization may be done

!----------------------------------------------------------------------
  !$omp parallel do private(iw,mrlw,ka,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
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
     ! (except heat/humidity/momentum which are done implicitly)

     call apply_surface_fluxes( iw )

     ! Compute tendency of scalars due to vertical turbulent transport

     call solve_eddy_diff_scalars(iw)

  enddo
  !$omp end parallel do

  ! MPI send/recv of computed K's

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vkm, rvara2=vkh, svara1=akm_sfc, r1dvara1=ustar, &
                          i1dvara1=khtop, i1dvara2=kmtop)
     call mpi_recv_w(mrl, rvara1=vkm, rvara2=vkh, svara1=akm_sfc, r1dvara1=ustar, &
                          i1dvara1=khtop, i1dvara2=kmtop)
  endif

 ! Lateral boundary copy of computed K's

  call lbcopy_w(mrl, a1=vkm, a2=vkh, s1=akm_sfc, v1=ustar, &
                     iv1=khtop, iv2=kmtop)

end subroutine pbl_driver

!===============================================================================

subroutine comp_horiz_k(mrl)

  use mem_ijtabs, only: jtab_v, jtv_wadj
  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iv

  !$omp parallel do private(iv)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

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

  dtl = dtlm(mrl)

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

  use mem_grid,      only: lsw, lpw, mwa, arw, arw0, zfacim2
  use mem_ijtabs,    only: jtab_w, jtw_prog, itabg_w, itab_w
  use mem_turb,      only: frac_urb, frac_land, frac_sea, frac_lake, &
                           frac_sfc, frac_sfck, ustar, wstar, wtv0, moli
  use mem_sfcg,       only: sfcg, itab_wsfc
  use misc_coms,     only: isubdomain, runtype, iparallel
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: isfcl
  use mem_land,      only: mland
  use mem_sea,       only: msea
  use consts_coms,   only: eps_virt, alvlocp
  use obnd,          only: lbcopy_w1d
  use olam_mpi_atm,  only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j, iw, iwsfc, k, km, ks, nland, nsea, jsfc, jasfc

! Populate urban and land, sea, and lake fraction arrays

  frac_urb (:) = 0.0
  frac_land(:) = 0.0
  frac_sea (:) = 0.0
  frac_lake(:) = 0.0

  if (isfcl > 0) then

     ! Loop over all ATM grid cells

     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

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
        call mpi_send_w(1, r1dvara1=frac_urb, r1dvara2=frac_land, &
                           r1dvara3=frac_sea, r1dvara4=frac_lake  )
        call mpi_recv_w(1, r1dvara1=frac_urb, r1dvara2=frac_land, &
                           r1dvara3=frac_sea, r1dvara4=frac_lake  )
     endif

     call lbcopy_w1d(1, a1=frac_urb, a2=frac_land, a3=frac_sea, a4=frac_lake)

  endif

  !$omp parallel do private(ks,k,km)
  do iw = 2, mwa

     ! Store the fraction of the total surface that intersects with each layer

     if (lsw(iw) == 1) then
        frac_sfc(1,iw) = 1.0
     else
        do ks = 1, lsw(iw)
           k  = lpw(iw) + ks - 1
           km = k - 1
           frac_sfc(ks,iw) = &
                max(arw(k,iw) * zfacim2(k) - arw(km,iw) * zfacim2(km),0.) / arw0(iw)
        enddo
     endif

     ! Store the fraction of the current cell that intersects with the surface

     frac_sfck(1,iw) = 1.0

     do ks = 2, lsw(iw)
       k  = lpw(iw) + ks - 1
        km = k - 1
        frac_sfck(ks,iw) = &
             max(arw(k,iw) * zfacim2(k) - arw(km,iw) * zfacim2(km),0.) / (arw(k,iw) * zfacim2(k))
     enddo

  enddo
  !$omp end parallel do

  if (runtype /= 'INITIAL') return

! Initialize PBL height and some PBL quantities

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

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
  use mem_ijtabs, only: itab_w
  use var_tables, only: sxfer_map, num_sxfer, emis_map, num_emis, scalar_tab

  implicit none

  integer, intent(in) :: iw
  integer             :: ns, n, ks, k
  real                :: dtl, dtli

  dtl  = dtlm(itab_w(iw)%mrlw)
  dtli = 1.0 / dtl

  ! Apply surface moisture and scalar fluxes directly to tendency arrays

  if (num_sxfer > 0) then

     do ns = 1, num_sxfer
        n = sxfer_map(ns)

        do ks = 1, lsw(iw)
           k = lpw(iw) + ks - 1

           scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                     + dtli * volti(k,iw) * scalar_tab(n)%sxfer(ks,iw)
        enddo
     enddo

  endif

  ! Apply emissions directly to the tendency arrays
  ! (emis units are concentration / sec )

  if (num_emis > 0) then

     do ns = 1, num_emis
        n = emis_map(ns)

        !dir$ ivdep
        do k = lpw(iw), mza-1
           scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                     + real(rho(k,iw)) * scalar_tab(n)%emis(k,iw)
        enddo
     enddo

  endif

end subroutine apply_surface_fluxes

!===============================================================================

subroutine solve_eddy_diff_scalars( iw )

  use mem_grid,   only: mza, arw, volti, lpw, lsw, dzim, arw0i, volt
  use mem_turb,   only: kpblh, vkh, sfluxr, sxfer_rk, akhs_dzi, agamma, khtop
  use mem_basic,  only: rho, rr_v
  use mem_ijtabs, only: itab_w
  use misc_coms,  only: dtlm
  use tridiag,    only: tridv, tridiffo
  use var_tables, only: num_pblmix, pblmix_map, scalar_tab
  !use mem_leaf,   only: itab_wl, land
  !use mem_sea,    only: itab_ws, sea
  use leaf_coms,  only: isfcl
  use mem_micro,  only: rr_c, rr_p
  use mem_tend,   only: rr_wt
  use oname_coms, only: nl

  use supercell_testm, only: rr_w_init

  implicit none

  integer, intent(in) :: iw

  real    :: akodz(mza), dtomass(mza), dtorho(mza), dgam(mza)
  integer :: ka, k, ks, jwl, iwl, jws, iws, kbot, n, ns
  real    :: vctr2(mza), vctr5(mza), vctr6(mza), vctr7(mza)
  real    :: soln(mza,num_pblmix), rhs(mza,num_pblmix)
  real    :: dtl, dtli, s0(num_pblmix)
  integer :: kmax, kmax_scal

  ka   = lpw(iw)
  dtl  = dtlm(1)
  dtli = 1.0 / dtl
  s0   = 0.

  kmax = max(khtop(iw)+1, ka+lsw(iw)-1)

  if (khtop(iw) < ka) then
     kmax_scal = ka-1
  else
     kmax_scal = khtop(iw)+1
  endif

  akodz(ka-1) = 0.
! akodz(mza)  = 0.
  akodz(khtop(iw)+1:kmax) = 0.

  ! Vertical loop over W levels

! do k = ka, mza-1
  do k = ka, khtop(iw)
     akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
  enddo

  ! Vertical loop over T levels

! do k = ka, mza
  do k = ka, kmax
     dtorho (k) = dtl / rho(k,iw)
     dtomass(k) = forw_imp * dtorho(k) * volti(k,iw)
     vctr5  (k) = -dtomass(k) * akodz(k-1)
     vctr7  (k) = -dtomass(k) * akodz(k)
     vctr6  (k) = 1. - vctr5(k) - vctr7(k)
  enddo

  ! Water mixing ratio

! do k = ka, mza
  do k = ka, kmax
     rhs(k,1) = rr_v(k,iw) + dtorho(k) * rr_wt(k,iw)
     if (allocated(rr_c)) rhs(k,1) = rhs(k,1) + rr_c(k,iw)
     if (allocated(rr_p)) rhs(k,1) = rhs(k,1) + rr_p(k,iw)

     ! For supercell test case, subtract initial field
     if (nl%test_case == 131) rhs(k,1) = rhs(k,1) - rr_w_init(k,iw)
  enddo

  ! Scalar variables with long-timestep forcing included

  do ns = 2, num_pblmix
     n = pblmix_map(ns)
!    do k = ka, mza
     do k = ka, kmax_scal
        rhs(k,ns) = scalar_tab(n)%var_p(k,iw)  &
                 + dtorho(k) * scalar_tab(n)%var_t(k,iw)
     enddo
  enddo

  ! Non local PBL terms

  if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
     if (agamma(ka,iw) > 1.e-7) then

        do k = ka, kpblh(iw)
           dgam(k) = dtomass(k) * (agamma(k-1,iw) - agamma(k,iw) )
        enddo

        ! water vapor (assumed to always be first scalar)
        s0(1) = sfluxr(iw)

        do ns = 2, num_pblmix
           n = pblmix_map(ns)

           ! other scalars with surface flux
           if (scalar_tab(n)%do_sxfer) then
              s0(ns) = sum( scalar_tab(n)%sxfer(1:lsw(iw),iw) / rho(ka:ka+lsw(iw)-1,iw) ) &
                     * dtli * arw0i(iw)
           endif

           ! other scalars with near-surface emissions
           if (scalar_tab(n)%do_emis) then
              s0(ns) = s0(n) + sum( real(volt(ka:ka+2,iw)) * scalar_tab(n)%emis(ka:ka+2,iw) ) &
                             / arw(ka+2,iw)
           endif
        enddo

        do ns = 1, num_pblmix
           if (abs(s0(ns)) > 1.e-7) then
              do k = ka, kpblh(iw)
                 rhs(k,ns) = rhs(k,ns) + s0(ns) * dgam(k)
              enddo
           endif
        enddo

     endif
  endif

  ! Solve tri-diagonal matrix equation for scalars

! if (ka <= mza) then
  if (ka <= kmax_scal .and. num_pblmix > 1) then
     call tridv(vctr5, vctr6, vctr7, rhs(:,2:), soln(:,2:), ka, kmax_scal, mza, num_pblmix-1)
  endif

  ! Special for water mixing ratio: implicit surface terms

  do k = ka, ka + lsw(iw) - 1
     ks = k - ka + 1
     vctr6(k) = vctr6(k) + dtomass(k) * akhs_dzi(ks,iw)
     rhs(k,1) = rhs(k,1) + dtomass(k) * sxfer_rk(ks,iw)
  enddo

  ! Solve tri-diagonal matrix equation for water mixing ratio

  if (ka <= kmax) then
     call tridiffo(mza, ka, kmax, vctr5, vctr6, vctr7, rhs(:,1), soln(:,1))
  endif

  vctr2(ka-1) = 0.
! vctr2(mza ) = 0.
  vctr2(khtop(iw)+1:kmax) = 0.

  do ns = 1, num_pblmix
     n = pblmix_map(ns)

     ! Now, soln contains future(t+1) values.
     ! Compute internal vertical turbulent fluxes

     if (abs(s0(ns)) > 1.e-7) then
        do k = ka, kpblh(iw)-1
           vctr2(k) = akodz(k) * (soln(k,ns) - soln(k+1,ns)) + s0(ns) * agamma(k,iw)
        enddo
        kbot = kpblh(iw)
     else
        kbot = ka
     endif

!    do k = kbot, mza-1
     do k = kbot, khtop(iw)
        vctr2(k) = akodz(k) * (soln(k,ns) - soln(k+1,ns))
     enddo

     ! Water vapor with implicitly-computed surface fluxes

     kbot = ka

     if (n == 1) then
        do k = ka, ka + lsw(iw) - 1
           ks = k - ka + 1

           scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                + volti(k,iw) * ( vctr2(k-1) - vctr2(k) &
                                - akhs_dzi(ks,iw) * soln(k,1) + sxfer_rk(ks,iw) )
        enddo

        kbot = ka + lsw(iw)
     endif

     ! Scalar tendencies from turbulent fluxes

!    do k = kbot, mza
     do k = kbot, kmax_scal
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + volti(k,iw) * (vctr2(k-1) - vctr2(k))
     enddo

  enddo

! Bob thinks the following section is not needed since fluxes are applied to land/lake/sea 
! in subroutine surface_fluxes.  Thus, commenting out the following...

  ! Compute surface fluxes of humidity for leaf given the implicit solution

!  if (isfcl > 0) then

!     do jwl = 1, itab_w(iw)%nland
!        iwl = itab_w(iw)%iland(jwl)
!        k   = itab_wl(iwl)%kw
!        ks  = k - ka + 1

!        land%sxfer_r(iwl) = land%ggaer(iwl) * (land%canshv(iwl) - soln(k,1)) * dtl * rho(k,iw)
!     enddo

!     do jws = 1, itab_w(iw)%nsea
!        iws = itab_w(iw)%isea(jws)
!        k   = itab_ws(iws)%kw
!        ks  = k - ka + 1

!        sea%sea_sxfer_r(iws) = sea%sea_ggaer(iws) * (sea%sea_canshv(iws) - soln(k,1)) * dtl * rho(k,iw)

!        if (sea%nlev_seaice(iws) > 0) then
!           sea%ice_sxfer_r(iws) = sea%ice_ggaer(iws) * (sea%ice_canshv(iws) - soln(k,1)) * dtl * rho(k,iw)
!        endif
!     enddo

!  endif

end subroutine solve_eddy_diff_scalars

!===============================================================================

subroutine solve_eddy_diff_heat(iw, thilt_short)

  use mem_turb,    only: vkh, sxfer_tk, akhs_dzi, kpblh, sfluxt, &
                         agamma, khtop
  use mem_basic,   only: rho, theta, tair
  use mem_micro,   only: rr_c, rr_p
  use misc_coms,   only: dtsm
  use mem_grid,    only: mza, lpw, lsw, volti, dzim, arw, nsw_max
  use oname_coms,  only: nl
  use tridiag,     only: tridiffo
  use mem_ijtabs,  only: itab_w
  use consts_coms, only: t00, p00, rocp, alvlocp, alviocp
  !use mem_leaf,    only: itab_wl, land
  !use mem_sea,     only: itab_ws, sea

  use supercell_testm, only: thil_init

  implicit none

  integer, intent(in )   :: iw
  real,    intent(inout) :: thilt_short(mza)  ! note: thilt assumed scaled by volume

  integer :: ka, k, ks, jwl, iwl, jws, iws, kbot, kmax
  real    :: dts, dtom
  real    :: akodz(mza), dtomass(mza), soln(mza), rhs(mza)
  real    :: vctr2(mza), vctr5(mza), vctr6(mza), vctr7(mza)
  real    :: exner(nsw_max)
  logical :: nonlocal

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
  dts = dtsm(1)

  kmax = max(khtop(iw)+1, ka+lsw(iw)-1)

  nonlocal = .false.

  akodz(ka-1)             = 0.
  akodz(khtop(iw)+1:kmax) = 0.
! akodz(mza)  = 0.

  ! Vertical loop over W levels

! do k = ka, mza-1
  do k = ka, khtop(iw)
     akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
  enddo

  ! Vertical loop over T levels

! do k = ka, mza
  do k = ka, kmax
     dtom       = dts * volti(k,iw) / rho(k,iw)
     dtomass(k) = forw_imp * dtom
     rhs    (k) = theta(k,iw) + dtom * thilt_short(k)
     vctr5  (k) = -dtomass(k) * akodz(k-1)
     vctr7  (k) = -dtomass(k) * akodz(k)
     vctr6  (k) = 1. - vctr5(k) - vctr7(k)

!    if (allocated(rr_c)) then
!       rhs(k) = rhs(k) - theta(k,iw) * alvlocp * rr_c(k,iw) / max(tair(k,iw),253.)
!    endif

!    if (allocated(rr_p)) then
!       rhs(k) = rhs(k) - theta(k,iw) * alviocp * rr_p(k,iw) / max(tair(k,iw),253.)
!    endif
  enddo

  ! Distribution of surface fluxes over multiple levels

  do k = ka, ka + lsw(iw) - 1
     ks = k - ka + 1

!    thl_cor(ks) = rhs(k) - varp(k)

     vctr6(k) = vctr6(k) + dtomass(k) * akhs_dzi(ks,iw)
     rhs  (k) = rhs  (k) + dtomass(k) * sxfer_tk(ks,iw)
!    rhs  (k) = rhs  (k) + dtomass(k) * (sxfer_tk(ks,iw) + akhs_dzi(ks,iw) * thl_cor(ks))
  enddo

  ! For supercell test case, for thil, subtract initial field

  if (nl%test_case == 131) then
     do k = ka, kmax
!    do k = ka, mza
        rhs(k) = rhs(k) - thil_init(k,iw)
     enddo
  endif

  ! Non local PBL term

  if (nl%idiffk(itab_w(iw)%mrlw) == 1) then
     if (agamma(ka,iw) > 1.e-7) then
        nonlocal = .true.

        do k = ka, kpblh(iw)
           rhs(k) = rhs(k) + dtomass(k) * sfluxt(iw) &
                           * ( agamma(k-1,iw) - agamma(k,iw) )
        enddo

     endif
  endif

  ! Solve tri-diagonal matrix equation

! if (ka <= mza) then
  if (ka <= kmax) then
     call tridiffo(mza, ka, kmax, vctr5, vctr6, vctr7, rhs, soln)
  endif

  ! Now, soln contains future(t+1) values.
  ! Compute internal vertical turbulent fluxes

  if (nonlocal) then
     do k = ka, kpblh(iw)-1
        vctr2(k) = akodz(k) * (soln(k) - soln(k+1)) + sfluxt(iw) * agamma(k,iw)
     enddo
     kbot = kpblh(iw)
  else
     kbot = ka
  endif

! do k = kbot, mza-1
  do k = kbot, khtop(iw)
     vctr2(k) = akodz(k) * (soln(k) - soln(k+1))
  enddo

  ! Set bottom and top internal fluxes to zero

  vctr2(ka-1) = 0.
! vctr2(mza ) = 0.
  vctr2(khtop(iw)+1:kmax) = 0.

  ! Include surface fluxes at ground

  do k = ka, ka + lsw(iw) - 1
     ks = k - ka + 1

!    airtheta(ks) = soln(k) - thl_cor(ks)

     thilt_short(k) = thilt_short(k) + &
          (vctr2(k-1) - vctr2(k) - akhs_dzi(ks,iw) * soln(k) + sxfer_tk(ks,iw))
!         (vctr2(k-1) - vctr2(k) - akhs_dzi(ks,iw) * airtheta(ks) + sxfer_tk(ks,iw))

     exner(ks) = dts * rho(k,iw) * tair(k,iw) / theta(k,iw)
  enddo

  ! Above surface

! do k = ka + lsw(iw), mza
  do k = ka + lsw(iw), kmax
     thilt_short(k) = thilt_short(k) + vctr2(k-1) - vctr2(k)
  enddo

! Bob thinks the following section is not needed since fluxes are applied to land/lake/sea 
! in subroutine surface_fluxes.  Thus, commenting out the following...

  ! Canopy heat fluxes

!  do jwl = 1, itab_w(iw)%nland
!     iwl = itab_w(iw)%iland(jwl)
!     k   = itab_wl(iwl)%kw
!     ks  = k - ka + 1

!     land%sxfer_t(iwl) = land%sxfer_t(iwl) &
!          + land%ggaer(iwl) * (land%cantheta(iwl) - soln(k)) * exner(ks)
!!         + land%ggaer(iwl) * (land%cantheta(iwl) - airtheta(ks)) * exner(ks)
!  enddo

!  do jws = 1, itab_w(iw)%nsea
!     iws = itab_w(iw)%isea(jws)
!     k   = itab_ws(iws)%kw
!     ks  = k - ka + 1

!     sea%sea_sxfer_t(iws) = sea%sea_sxfer_t(iws) &
!          + sea%sea_ggaer(iws) * (sea%sea_cantheta(iws) - soln(k)) * exner(ks)
!!         + sea%sea_ggaer(iws) * (sea%sea_cantheta(iws) - airtheta(ks)) * exner(ks)

!     if (sea%nlev_seaice(iws) > 0) then
!        sea%ice_sxfer_t(iws) = sea%ice_sxfer_t(iws) &
!             + sea%ice_ggaer(iws) * (sea%ice_cantheta(iws) - soln(k)) * exner(ks)
!!            + sea%ice_ggaer(iws) * (sea%ice_cantheta(iws) - airtheta(ks)) * exner(ks)
!     endif
!  enddo

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

  dts = dtsm(1)
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
  dts = dtsm(1)

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
