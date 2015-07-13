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

subroutine pbl_driver(mrl)

  use mem_grid,       only: mwa, mza, lpw, lsw,lpv, arv, volt, dniv
  use misc_coms,      only: io6, idiffk, dtlm, iparallel
  use mem_tend,       only: thilt, sh_wt
  use mem_basic,      only: vxe, vye, vze, thil, theta, tair, sh_w, sh_v, rho
  use mem_turb,       only: hkm, hkh, sxfer_tk, sxfer_rk, ustar, wstar, wtv0, &
                            frac_sfc, pblh, kpblh, fthpbl, fqtpbl, akmodx, akhodx
  use consts_coms,    only: grav, vonk, eps_virt, alvlocp, r8
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog, jtab_v, jtv_wadj, itab_v
  use mem_radiate,    only: fthrd_sw, fthrd_lw
  use module_bl_acm2, only: acm2_pblhgt, acm2_eddyx, acm2
  use smagorinsky,    only: turb_k
  use mem_micro,      only: sh_c
  use mem_grid,       only: zm
  use olam_mpi_atm,   only: mpi_send_w, mpi_recv_w  
  use obnd,           only: lbcopy_w

  implicit none

  integer, intent(in) :: mrl
  integer :: j, k, ka, iw, mrlw, ks, iw1, iw2, iv

  real     :: qc    (mza)
  real     :: thlv  (mza)
  real     :: vkh   (mza)
  real     :: vkm   (mza)
  real     :: fthrd (mza)
  real     :: moli
  real(r8) :: fact1, fact2
  real     :: tempm, temph, stab1, stab2

! Loop over all W/T points where PBL parameterization may be done

!----------------------------------------------------------------------
  !$omp parallel do private(iw,mrlw,ka,k,qc,thlv,moli,vkh,vkm,fthrd) 
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

     ! MRL for current IW column

     mrlw = itab_w(iw)%mrlw
     ka   = lpw(iw)

     ! Save input tendencies of theta and qt for later determining PBL tendencies

     do k = ka, mza
        fthpbl(k,iw) = thilt(k,iw)
        fqtpbl(k,iw) = sh_wt(k,iw)
     enddo

     ! Define a cloud-water virtual potential temperature for diagnosing buoyancy

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

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), ka, mza-1, lsw(iw),     &
                       frac_sfc(:,iw), thlv, vxe(:,iw), vye(:,iw), vze(:,iw), &
                       kpblh(iw), pblh(iw) )

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 1) then

        ! ACM2 non-local convective tranport scheme

        moli = - grav * vonk * wtv0(iw) / ustar(iw)**3 / thlv(ka)

        fthrd(ka:mza) = fthrd_sw(ka:mza,iw) + fthrd_lw(ka:mza,iw)

        call acm2_eddyx( iw, moli, ustar(iw), pblh(iw), kpblh(iw), ka, mza, &
                         vkh, vxe(:,iw), vye(:,iw), vze(:,iw), sh_w(:,iw),    &
                         sh_v(:,iw), qc, thil(:,iw), theta(:,iw), thlv,       &
                         tair(:,iw), fthrd, rho(:,iw)                         )

        call acm2( iw, moli, ustar(iw), pblh(iw), kpblh(iw), vkh )

        ! temporarily get horizontal diffusion coefficient from vkh

        hkm(ka,iw) = vkh(ka)
        hkh(ka,iw) = hkm(ka,iw)
        do k = ka+1, mza
           hkm(k,iw) = 0.5 * (vkh(k-1) + vkh(k))
           hkh(k,iw) = hkm(k,iw)
        enddo

     else if (idiffk(mrlw) == 2 .or. idiffk(mrlw) == 3) then

        ! Smagorinsky scheme

        call turb_k(iw, mrlw, thlv, vkh, vkm)

     endif

  enddo
  !$omp end parallel do

  ! Parallel send of horizontal K's

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=hkm, rvara2=hkh)
  endif

!----------------------------------------------------------------------
  !$omp parallel do private(iw,ks,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

     ! Zero out sxfer arrays now that they have been transferred to the atm
     
     do ks = 1, lsw(iw)
        sxfer_tk(ks,iw) = 0.0
        sxfer_rk(ks,iw) = 0.0
     enddo

     ! Save PBL tendencies of theta and qt needed by some convective schemes

     do k = lpw(iw), mza
        fthpbl(k,iw) = (thilt(k,iw) - fthpbl(k,iw)) / rho(k,iw)
        fqtpbl(k,iw) = (sh_wt(k,iw) - fqtpbl(k,iw)) / rho(k,iw)
     enddo

  enddo
  !$omp end parallel do

! MPI Recv and LBC copy of K's

  if (iparallel == 1) then
     call mpi_recv_w(mrl, rvara1=hkm, rvara2=hkh)
  endif
  call lbcopy_w(mrl, a1=hkm, a2=hkh)

! Loop over V columns to compute ARV * K / DX, and make sure
! horizontal diffusion is stable over the long timestep

!----------------------------------------------------------------------
  !$omp parallel do private(iv,iw1,iw2,fact1,fact2,k,tempm,temph,stab1,stab2)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)
!----------------------------------------------------------------------

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     fact1 = 1.0_r8 / ( dtlm(itab_w(iw1)%mrlw) * itab_w(iw1)%npoly )
     fact2 = 1.0_r8 / ( dtlm(itab_w(iw2)%mrlw) * itab_w(iw2)%npoly )

     do k = lpv(iv), mza
        tempm = 0.5 * dniv(iv) * arv(k,iv) * (hkm(k,iw1) + hkm(k,iw2))
        temph = 0.5 * dniv(iv) * arv(k,iv) * (hkh(k,iw1) + hkh(k,iw2))
        stab1 = rho(k,iw1) * volt(k,iw1) * fact1
        stab2 = rho(k,iw2) * volt(k,iw2) * fact2
        akmodx(k,iv) = min( tempm, stab1, stab2 )
        akhodx(k,iv) = min( temph, stab1, stab2 )
     enddo

  enddo
  !$omp end parallel do

end subroutine pbl_driver

!===============================================================================

subroutine pbl_init()

  use mem_grid,      only: lsw, lpw, mza, mwa, arw, arw0, zfacm
  use mem_ijtabs,    only: itab_w, jtab_w, jtw_prog, itabg_w
  use mem_turb,      only: frac_urb, frac_land, frac_sea, frac_lake, frac_sfc, &
                           ustar, wstar, wtv0, pblh, kpblh, fthpbl, fqtpbl
  use mem_leaf,      only: land, itab_wl
  use mem_sea,       only: sea, itab_ws
  use mem_basic,     only: vxe, vye, vze, theta, tair, sh_v
  use misc_coms,     only: io6, isubdomain, runtype, iparallel
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
                (arw(k,iw) / zfacm(k)**2 - arw(km,iw) / zfacm(km)**2) / arw0(iw)
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

     do k = lpw(iw), mza
        thlv(k) = theta(k,iw) * (1.0 + eps_virt * sh_v(k,iw))
     enddo

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), lpw(iw), mza-1, lsw(iw), &
                       frac_sfc(:,iw), thlv, vxe(:,iw), vye(:,iw), vze(:,iw),  &
                       kpblh(iw), pblh(iw)                                      )
  enddo
!$omp end parallel do

end subroutine pbl_init
