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

  use mem_grid,       only: mza, mwa, lpw, lsw, volti
  use misc_coms,      only: idiffk, dtlm
  use mem_tend,       only: thilt, sh_wt
  use mem_basic,      only: vxe, vye, vze, thil, theta, tair, sh_w, sh_v, rho
  use mem_turb,       only: vkm, vkh, sxfer_rk, ustar, wstar, wtv0, &
                            frac_sfc, pblh, kpblh, fqtpbl, fthpbl, moli
  use consts_coms,    only: grav, vonk, eps_virt, alvlocp, r8
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog
  use mem_radiate,    only: pbl_cld_forc
  use module_bl_acm2, only: acm2_pblhgt, acm2_eddyx, acm2_scalars, acm2_momentum
  use smagorinsky,    only: turb_k
  use mem_micro,      only: sh_c

  implicit none

  integer, intent(in)    :: mrl
  real,    intent(inout) :: rhot(mza,mwa)

  integer :: j, k, ka, iw, mrlw, ks

  real    :: qc(mza)
  real    :: thlv(mza)
  real    :: dtli

! Loop over all W/T points where PBL parameterization may be done

!----------------------------------------------------------------------
  !$omp parallel do private(iw,mrlw,ka,k,qc,thlv)
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

     ! Define a cloud-liquid-water virtual potential temperature for diagnosing buoyancy

     if (allocated(sh_c)) then
        do k = ka, mza
           qc  (k) = max( sh_w(k,iw) - sh_v(k,iw), 0.0 )
           thlv(k) = theta(k,iw) * (1.0 + eps_virt * sh_v(k,iw) - qc(k)) &
                   / ( 1.0 + alvlocp * sh_c(k,iw) / max(tair(k,iw), 253.0) )
        enddo
     else
        do k = ka, mza
           qc  (k) = 0.0
           thlv(k) = theta(k,iw) * (1.0 + eps_virt * sh_v(k,iw))
        enddo
     endif

     ! Diagnose PBL height regardless of scheme

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), ka, mza-1, lsw(iw),    &
                       frac_sfc(:,iw), thlv, vxe(:,iw), vye(:,iw), vze(:,iw), &
                       kpblh(iw), pblh(iw) )

     ! Apply temperature and scalar surface fluxes and emissions

     call apply_surface_fluxes( iw )

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 1) then

        ! ACM2 non-local convective tranport scheme

        call acm2_eddyx( iw, moli(iw), ustar(iw), wstar(iw), pblh(iw),        &
                         kpblh(iw), ka, mza, vkm(:,iw), vkh(:,iw),            &
                         vxe(:,iw), vye(:,iw), vze(:,iw),                     &
                         sh_w(:,iw), sh_v(:,iw), qc, thil(:,iw), theta(:,iw), &
                         thlv, tair(:,iw), pbl_cld_forc(iw), rho(:,iw)        )

        call acm2_scalars ( iw, moli(iw), pblh(iw), kpblh(iw), vkh(:,iw) )
        call acm2_momentum( iw, vkm(:,iw) )

     else if (idiffk(mrlw) == 2 .or. idiffk(mrlw) == 3) then

        ! Smagorinsky scheme

        call turb_k(iw, mrlw, thlv, vkh(:,iw), vkm(:,iw))

     else

        ! If no PBL diffusion, just apply momentum fluxes

        vkm(:,iw) = 0.0
        vkh(:,iw) = 0.0
        call apply_momentum_fluxes( iw )

     endif

     ! Add surface vapor flux to total density tendency

     dtli = 1.0 / dtlm(mrlw)

     do ks = 1, lsw(iw)
        k = lpw(iw) + ks - 1
        rhot(k,iw) = rhot(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)
     enddo

     ! Save PBL tendencies of qt needed by some convective schemes

     do k = lpw(iw), mza
        fqtpbl(k,iw) = (sh_wt(k,iw) - fqtpbl(k,iw)) / rho(k,iw)
        fthpbl(k,iw) = (thilt(k,iw) - fthpbl(k,iw)) / rho(k,iw)
     enddo

  enddo
  !$omp end parallel do

end subroutine pbl_driver


!===============================================================================


subroutine comp_horiz_k(mrl)

  use consts_coms, only: r8
  use mem_grid,    only: arw0, lpv, mza, arv, volt, dniv, zm, zt
  use misc_coms,   only: dtlm, akmin
  use mem_ijtabs,  only: jtab_v, jtv_wadj, itab_v, itab_w
  use mem_turb,    only: vkm, vkh, akmodx, akhodx
  use mem_basic,   only: rho
  use mem_cuparm,  only: kcutop, kcubot, cbmf, iactcu

  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iv, iw1, iw2, k, km
  real                :: bkmin, dens
  real                :: hkm, hkh, tempm, temph, stab1, stab2
  real(r8)            :: fact1, fact2
  real                :: hcm, hch, hkc1(mza), hkc2(mza), eta

  ! Loop over V columns to compute ARV * K / DX, and make sure
  ! horizontal diffusion is stable over the long timestep

  !$omp parallel do private(iv,iw1,iw2,fact1,fact2,bkmin,k,km,&
  !$                        dens,hkm,hkh,tempm,temph,stab1,stab2,&
  !$                        hcm,hch,hkc1,hkc2,eta)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     fact1 = 0.95_r8 / ( dtlm(itab_w(iw1)%mrlw) * itab_w(iw1)%npoly )
     fact2 = 0.95_r8 / ( dtlm(itab_w(iw2)%mrlw) * itab_w(iw2)%npoly )

     bkmin = akmin(itab_v(iv)%mrlv) * .075 * min(arw0(iw1),arw0(iw2)) ** .66666666

!   ! Extra mixing around convection???
!    hkc1(:) = 0.0
!    hkc2(:) = 0.0
!
!    if (iactcu(iw1)) then
!       do k = kcubot(iw1), kcutop(iw1)
!          eta = ( zt(k) - zm(kcubot(iw1)-1) ) / ( zm(kcutop(iw1)) - zm(kcubot(iw1)-1) )
!          hkc1(k) = cbmf(iw1) * 2200.0 * eta * (1.0-eta)**2 * rho(k,iw1)
!       enddo
!    endif
!
!    if (iactcu(iw2)) then
!       do k = kcubot(iw2), kcutop(iw2)
!          eta = ( zt(k) - zm(kcubot(iw2)-1) ) / ( zm(kcutop(iw2)) - zm(kcubot(iw2)-1) )
!          hkc2(k) = cbmf(iw2) * 2200.0 * eta * (1.0-eta)**2 * rho(k,iw1)
!       enddo
!    endif

     do k = lpv(iv), mza

        if (k == lpv(iv)) then
           km = lpv(iv)
        else
           km = k - 1
        endif

        dens = 0.5 * (rho(k,iw1) + rho(k,iw2))

!       ! Extra mixing around convection???
!       hch = 0.5 * ( hkc1(k)+hkc2(k) )
!       hcm = 0.5 * hch
!       hkm = max(0.25 * (vkm(k,iw1) + vkm(k,iw2) + vkm(km,iw1) + vkm(km,iw2)) + hcm, bkmin * dens)
!       hkh = max(0.25 * (vkh(k,iw1) + vkh(k,iw2) + vkh(km,iw1) + vkh(km,iw2)) + hch, bkmin * dens)

        hkm = max(0.25 * (vkm(k,iw1) + vkm(k,iw2) + vkm(km,iw1) + vkm(km,iw2)), bkmin * dens)
        hkh = max(0.25 * (vkh(k,iw1) + vkh(k,iw2) + vkh(km,iw1) + vkh(km,iw2)), bkmin * dens)

        tempm = dniv(iv) * arv(k,iv) * hkm
        temph = dniv(iv) * arv(k,iv) * hkh

        stab1 = rho(k,iw1) * volt(k,iw1) * fact1
        stab2 = rho(k,iw2) * volt(k,iw2) * fact2

        akmodx(k,iv) = min( tempm, stab1, stab2 )
        akhodx(k,iv) = min( temph, stab1, stab2 )
     enddo

  enddo
  !$omp end parallel do

end subroutine comp_horiz_k

!===============================================================================

subroutine pbl_init()

  use mem_grid,      only: lsw, lpw, mza, mwa, arw, arw0, zfacim2
  use mem_ijtabs,    only: jtab_w, jtw_prog, itabg_w
  use mem_turb,      only: frac_urb, frac_land, frac_sea, frac_lake, frac_sfc, &
                           ustar, wstar, wtv0, pblh, kpblh, fthpbl, fqtpbl, moli
  use mem_leaf,      only: land, itab_wl
  use mem_sea,       only: sea, itab_ws
  use mem_basic,     only: vxe, vye, vze, theta, sh_v
  use misc_coms,     only: isubdomain, runtype, iparallel
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: mwl, isfcl
  use sea_coms,      only: mws
  use consts_coms,   only: eps_virt, alvlocp
  use obnd,          only: lbcopy_w1d
  use olam_mpi_atm,  only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j, iw, iwl, iws, k, km, ks
  real    :: thlv(mza)

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

! Store the fraction of the total surface that intersects with each layer

  do iw = 2, mwa
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
  enddo

  if (runtype /= 'INITIAL') return

! Initialize PBL height and some PBL quantities
    
!$omp parallel do private(iw,k,thlv)
  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     ustar(iw) = 0.2
     wstar(iw) = 0.0
     wtv0 (iw) = 0.0
     moli (iw) = 0.0

     do k = lpw(iw), mza
        thlv(k) = theta(k,iw) * (1.0 + eps_virt * sh_v(k,iw))
     enddo

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), lpw(iw), mza-1, lsw(iw), &
                       frac_sfc(:,iw), thlv, vxe(:,iw), vye(:,iw), vze(:,iw),  &
                       kpblh(iw), pblh(iw)                                      )
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
                                     + rho(k,iw) * scalar_tab(n)%emis(k,iw)
        enddo
     enddo
     
  endif

end subroutine apply_surface_fluxes

!===============================================================================

subroutine apply_momentum_fluxes( iw )

  use mem_grid,    only: lpw, lsw, volti, arw, dzim
  use consts_coms, only: r8
  use misc_coms,   only: dtlm
  use mem_tend,    only: vmxet, vmyet, vmzet
  use mem_basic,   only: vxe, vye, vze
  use mem_ijtabs,  only: itab_w
  use mem_turb,    only: vkm_sfc, sxfer_tk

  implicit none

  integer, intent(in) :: iw
  integer             :: ks, k
  real(r8)            :: dtli
  real                :: fact

  dtli = 1.0_r8 / dtlm(itab_w(iw)%mrlw)

  do k = lpw(iw), lpw(iw) + lsw(iw) - 1
     ks = k - lpw(iw) + 1

     fact = 2.0 * vkm_sfc(ks,iw) * dzim(k-1) * (arw(k,iw) - arw(k-1,iw)) * volti(k,iw)

     vmxet(k,iw) = vmxet(k,iw) - fact * vxe(k,iw)
     vmyet(k,iw) = vmyet(k,iw) - fact * vye(k,iw)
     vmzet(k,iw) = vmzet(k,iw) - fact * vze(k,iw)
  enddo

end subroutine apply_momentum_fluxes
