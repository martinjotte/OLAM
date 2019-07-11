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

subroutine pbl_driver(mrl,rhot)

  use mem_grid,       only: mza, mwa, lpw, lsw, volti, xew, yew, zew
  use misc_coms,      only: idiffk, dtlm, mdomain
  use mem_tend,       only: thilt, sh_wt
  use mem_basic,      only: vxe, vye, vze, rho
  use mem_turb,       only: vkm, vkh, sxfer_rk, ue, ve, fqtpbl, fthpbl
  use consts_coms,    only: grav, vonk, eps_virt, alvlocp, r8, eradi
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog, mrls
  use module_bl_acm2, only: acm2_pblhgt, acm2_driver
  use smagorinsky,    only: turb_k

  implicit none

  integer, intent(in)    :: mrl
  real,    intent(inout) :: rhot(mza,mwa)

  integer :: j, k, ka, iw, mrlw, ks
  real    :: dtli, raxis, raxisi

  ! Compute zonal and meridional wind components

  if (any(idiffk(1:mrls) == 3)) then

     !$omp parallel do private(iw,raxis,raxisi,k)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
        raxisi = 1.0 / max(raxis, 1.e-12)

        if (mdomain <= 1 .and. raxis > 1.e3) then

           do k = lpw(iw), mza
              ue(k,iw) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
              ve(k,iw) = vze(k,iw) * raxis * eradi  &
                       - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) &
                         * zew(iw) * raxisi * eradi
           enddo

        else

           do k = lpw(iw), mza
              ue(k,iw) = vxe(k,iw)
              ve(k,iw) = vye(k,iw)
           enddo

        endif

     enddo
     !$omp end parallel do
  endif

! Loop over all W/T points where PBL parameterization may be done

!----------------------------------------------------------------------
  !$omp parallel do private(iw,mrlw,ka,k,ks,dtli)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

     ! MRL for current IW column

     mrlw = itab_w(iw)%mrlw
     ka   = lpw(iw)

     ! Save input tendencies of thil and qt for later determining PBL tendencies

     do k = ka, mza
        fthpbl(k,iw) = thilt(k,iw)
        fqtpbl(k,iw) = sh_wt(k,iw)
     enddo

     ! Diagnose PBL height regardless of scheme

     call acm2_pblhgt(iw)

     ! Apply temperature and scalar surface fluxes and emissions

     call apply_surface_fluxes( iw )

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 1) then

        ! ACM2 non-local convective tranport scheme

        call acm2_driver( iw )

     else if (idiffk(mrlw) == 2 .or. idiffk(mrlw) == 3) then

        ! Smagorinsky scheme

        call turb_k(iw, mrlw)

     else

        ! If no PBL diffusion, just apply momentum fluxes

        vkm(:,iw) = 0.0
        vkh(:,iw) = 0.0
        call apply_momentum_fluxes( iw )

     endif

     vkm(lpw(iw)-1,iw) = vkm(lpw(iw),iw)
     vkh(lpw(iw)-1,iw) = vkh(lpw(iw),iw)

     ! Add surface vapor flux to total density tendency

     dtli = 1.0 / real(dtlm(mrlw))

     do ks = 1, lsw(iw)
        k = lpw(iw) + ks - 1
        rhot(k,iw) = rhot(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)
     enddo

     ! Save PBL tendencies of qt needed by some convective schemes

     do k = lpw(iw), mza
        fqtpbl(k,iw) = (sh_wt(k,iw) - fqtpbl(k,iw)) / real(rho(k,iw))
        fthpbl(k,iw) = (thilt(k,iw) - fthpbl(k,iw)) / real(rho(k,iw))
     enddo

  enddo
  !$omp end parallel do

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



subroutine comp_horiz_k_column(iv)

  use consts_coms, only: r8
  use mem_grid,    only: arw0, lpv, mza, arv, volt, dniv, zm, zt
  use misc_coms,   only: dtlm, akmin
  use mem_ijtabs,  only: itab_v, itab_w
  use mem_turb,    only: vkm, vkh, akmodx, akhodx
  use mem_basic,   only: rho
  use mem_cuparm,  only: kcutop, kcubot, cbmf, iactcu, conprr
  use oname_coms,  only: nl

  implicit none

  integer, intent(in) :: iv
  integer             :: iw1, iw2, k, km, mrl
  real                :: bkmin, dens
  real                :: hkm, hkh, tempm, temph, stab1, stab2
  real(r8)            :: fact1, fact2
  real                :: hcm, hkc(mza)
  real, parameter     :: cbmf0 = 0.066

  ! Compute ARV * K / DX, and ensure that horizontal
  ! diffusion is stable over the long timestep

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)
  mrl = itab_v(iv)%mrlv

  fact1 = 0.95_r8 / ( dtlm(itab_w(iw1)%mrlw) * itab_w(iw1)%npoly )
  fact2 = 0.95_r8 / ( dtlm(itab_w(iw2)%mrlw) * itab_w(iw2)%npoly )

  bkmin = akmin(mrl) * .075 * min(arw0(iw1),arw0(iw2)) ** .66666666

  hkc(:) = 0.0

  ! Extra mixing around deep convection???

  if (nl%conv_akmin(mrl) > 1.e-7) then

     if (iactcu(iw1) > 0 .and. conprr(iw1) > 1.e-8) then
        hcm = 0.0375 * nl%conv_akmin(mrl) * arw0(iw1) ** .66666666 &
                     * cbmf(iw1) / cbmf0

        do k = lpv(iv), min(kcutop(iw1) + 1, mza)
           hkc(k) = hkc(k) + hcm
        enddo
     endif

     if (iactcu(iw2) > 0 .and. conprr(iw2) > 1.e-8) then
        hcm = 0.0375 * nl%conv_akmin(mrl) * arw0(iw2) ** .66666666 &
                     * cbmf(iw2) / cbmf0

        do k = lpv(iv), min(kcutop(iw2) + 1, mza)
           hkc(k) = hkc(k) + hcm
        enddo
     endif

  endif

  do k = lpv(iv), mza-1
    dens = 0.5 * real(rho(k,iw1) + rho(k,iw2))

    hkm = max(0.25 * (vkm(k,iw1) + vkm(k,iw2) + vkm(k-1,iw1) + vkm(k-1,iw2)), &
              hkc(k) * dens, bkmin * dens)
    hkh = max(0.25 * (vkh(k,iw1) + vkh(k,iw2) + vkh(k-1,iw1) + vkh(k-1,iw2)), &
              hkc(k) * dens, bkmin * dens)

    tempm = dniv(iv) * arv(k,iv) * hkm
    temph = dniv(iv) * arv(k,iv) * hkh

    stab1 = rho(k,iw1) * volt(k,iw1) * fact1
    stab2 = rho(k,iw2) * volt(k,iw2) * fact2

    akmodx(k,iv) = min( tempm, stab1, stab2 )
    akhodx(k,iv) = min( temph, stab1, stab2 )
 enddo

 akmodx(mza,iv) = akmodx(mza-1,iv)
 akhodx(mza,iv) = akhodx(mza-1,iv)

end subroutine comp_horiz_k_column

!===============================================================================

subroutine pbl_init()

  use mem_grid,      only: lsw, lpw, mwa, arw, arw0, zfacim2
  use mem_ijtabs,    only: jtab_w, jtw_prog, itabg_w
  use mem_turb,      only: frac_urb, frac_land, frac_sea, frac_lake, frac_sfc, &
                           ustar, wstar, wtv0, fthpbl, fqtpbl, moli, frac_sfck
  use mem_leaf,      only: land, itab_wl
  use mem_sea,       only: sea, itab_ws
  use misc_coms,     only: isubdomain, runtype, iparallel
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: mwl, isfcl
  use sea_coms,      only: mws
  use consts_coms,   only: eps_virt, alvlocp
  use obnd,          only: lbcopy_w1d
  use olam_mpi_atm,  only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j, iw, iwl, iws, k, km, ks

  if (allocated(fthpbl)) fthpbl(:,:) = 0.0
  if (allocated(fqtpbl)) fqtpbl(:,:) = 0.0

! Populate urban and land, sea, and lake fraction arrays

  frac_urb (:) = 0.0
  frac_land(:) = 0.0
  frac_sea (:) = 0.0
  frac_lake(:) = 0.0

  if (isfcl > 0) then

     do iwl = 2, mwl
        iw  = itab_wl(iwl)%iw   ! global index
        if (isubdomain == 1) then
           iw = itabg_w(iw)%iw_myrank  ! local index
        endif

        if (any(land%leaf_class(iwl) == (/ 19, 21 /))) then
           frac_urb(iw) = frac_urb(iw) + itab_wl(iwl)%arf_iw
        endif

        frac_land(iw) = frac_land(iw) + itab_wl(iwl)%arf_iw
     enddo

     do iws = 2, mws
        iw  = itab_ws(iws)%iw   ! global index
        if (isubdomain == 1) then
           iw  = itabg_w(iw)%iw_myrank  ! local index
        endif

        if (sea%leaf_class(iws) == 0) then
           frac_sea (iw) = frac_sea (iw) + itab_ws(iws)%arf_iw
        else
           frac_lake(iw) = frac_lake(iw) + itab_ws(iws)%arf_iw
        endif
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

  use mem_grid,    only: mza, lpw, lsw, volti
  use consts_coms, only: r8
  use misc_coms,   only: dtlm
  use mem_basic,   only: rho
  use mem_ijtabs,  only: itab_w
  use var_tables,  only: sxfer_map, num_sxfer, emis_map, num_emis, scalar_tab
  use mem_turb,    only: sxfer_tk
  use mem_tend,    only: thilt

  implicit none

  integer, intent(in) :: iw
  integer             :: ns, n, ks, k
  real(r8)            :: dtli

  dtli = 1.0_r8 / dtlm(itab_w(iw)%mrlw)

  ! Apply surface heat flux directly to thil tendency

  do ks = 1, lsw(iw)
     k = lpw(iw) + ks - 1
     thilt(k,iw) = thilt(k,iw) + dtli * volti(k,iw) * sxfer_tk(ks,iw)
  enddo

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

        do k = lpw(iw), mza-1
           scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                     + real(rho(k,iw)) * scalar_tab(n)%emis(k,iw)
        enddo
     enddo

  endif

end subroutine apply_surface_fluxes

!===============================================================================

subroutine apply_momentum_fluxes( iw )

  use mem_grid,    only: lpw, lsw, volti, arw, dzim
  use mem_tend,    only: vmxet, vmyet, vmzet
  use mem_basic,   only: vxe, vye, vze
  use mem_turb,    only: vkm_sfc

  implicit none

  integer, intent(in) :: iw
  integer             :: ks, k
  real                :: fact

  do k = lpw(iw), lpw(iw) + lsw(iw) - 1
     ks = k - lpw(iw) + 1

     fact = 2.0 * vkm_sfc(ks,iw) * dzim(k-1) * (arw(k,iw) - arw(k-1,iw)) * real(volti(k,iw))

     vmxet(k,iw) = vmxet(k,iw) - fact * vxe(k,iw)
     vmyet(k,iw) = vmyet(k,iw) - fact * vye(k,iw)
     vmzet(k,iw) = vmzet(k,iw) - fact * vze(k,iw)
  enddo

end subroutine apply_momentum_fluxes

!===============================================================================
