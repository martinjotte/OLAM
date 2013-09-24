!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

subroutine pbl_driver(rhot, mrl)

  use mem_grid,       only: mwa, mza, lpw, lsw
  use misc_coms,      only: io6, idiffk, do_chem
  use mem_tend,       only: thilt, sh_wt
  use mem_basic,      only: vxe, vye, vze, thil, theta, tair, sh_w, sh_v, rho
  use mem_turb,       only: hkm, sxfer_tk, sxfer_rk, ustar, wstar, wtv0, &
                            frac_urb, frac_land, frac_sfc, pblh, kpblh
  use consts_coms,    only: grav, vonk, eps_virt
  use mem_ijtabs,     only: jtab_w, itab_w, jtw_prog
  use mem_radiate,    only: fthrd
  use module_bl_acm2, only: acm2_pblhgt, acm2_eddyx, acm2
  use smagorinsky,    only: turb_k
  use emis_defn,      only: get_emis
  use depv_defn,      only: get_depv
  use mem_megan,      only: megan_driver

!$use omp_lib

  implicit none

  real, intent(inout) :: rhot(mza,mwa)
  integer, intent(in) :: mrl
  integer :: j, k, ka, iw, mrlw, ks

  real    :: qc    (mza)
  real    :: thetav(mza)
  real    :: thilv (mza)
  real    :: vkh   (mza)
  real    :: vkm   (mza)
  real    :: moli

! Get chemical emisions for this timestep

  if (do_chem) then
     call megan_driver(mrl)
     call get_emis(mrl)
     call get_depv(mrl)
  endif

! Loop over all W/T points where PBL parameterization may be done

  call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,mrlw,ka,k,ks,qc,thetav,thilv,moli,vkh,vkm) 
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
     call qsub('W',iw)

     ! MRL for current IW column

     mrlw = itab_w(iw)%mrlw
     ka   = lpw(iw)

     do k = ka-1, mza-1
        hkm(k,iw) = 0.
!       vkm(k,iw) = 0.
!       vkh(k,iw) = 0.
     enddo

     ! Diagnose PBL height regardless of scheme

     do k = ka, mza-1
        qc    (k) = max( sh_w(k,iw) - sh_v(k,iw), 0.0 )
        thetav(k) = theta(k,iw) * (1.0 + eps_virt * sh_v(k,iw) - qc(k))
        thilv (k) = thil (k,iw) * (1.0 + eps_virt * sh_v(k,iw) - qc(k))
     enddo

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), ka, mza-2, lsw(iw),     &
                       frac_sfc(:,iw), thilv, vxe(:,iw), vye(:,iw), vze(:,iw), &
                       kpblh(iw), pblh(iw) )

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 1) then

        ! ACM2 non-local convective tranport scheme

        moli = - grav * vonk * wtv0(iw) / ustar(iw)**3 / thetav(ka)

        call acm2_eddyx( iw, moli, ustar(iw), pblh(iw), kpblh(iw), ka, mza-1, &
                         vkh, vxe(:,iw), vye(:,iw), vze(:,iw), sh_w(:,iw),    &
                         sh_v(:,iw), qc, thil(:,iw), theta(:,iw), thetav,     &
                         tair(:,iw), fthrd(:,iw), rho(:,iw)                   )

        call acm2( iw, rhot, moli, ustar(iw), pblh(iw), kpblh(iw), thetav, vkh )

        ! get horizontal diffusion coefficient from vkh
        
        hkm(ka,iw) = vkh(ka)
        do k = ka, mza-1
           hkm(k,iw) = 0.5 * (vkh(k-1) + vkh(k))
        enddo

     else if (idiffk(mrlw) == 2 .or. idiffk(mrlw) == 3) then
   
        ! Smagorinsky scheme

        call turb_k(iw, mrlw, rhot, thetav, vkh, vkm)

     endif

     ! Zero out sxfer arrays now that they have been transferred to the atm
     
     do ks = 1, lsw(iw)
        sxfer_tk(ks,iw) = 0.0
        sxfer_rk(ks,iw) = 0.0
     enddo
   
  enddo
!$omp end parallel do

end subroutine pbl_driver

!===============================================================================

subroutine pbl_init()

  use mem_sflux,     only: landflux, jlandflux
  use mem_grid,      only: lsw, lpw, mza, mwa, arw
  use mem_ijtabs,    only: itab_w, jtab_w, jtw_prog, itabg_w
  use mem_turb,      only: frac_urb, frac_land, frac_sfc, ustar, wstar, wtv0, &
                           pblh, kpblh
  use mem_leaf,      only: land, itabg_wl
  use mem_basic,     only: vxe, vye, vze, thil, sh_v
  use misc_coms,     only: io6, iparallel
  use module_bl_acm2,only: acm2_pblhgt
  use leaf_coms,     only: isfcl
  use consts_coms,   only: eps_virt
  
  implicit none

  integer :: j, ilf, iw, iwl, k
  real    :: thilv(mza)

! Populate urban and land fraction arrays if used
  
  if (allocated(frac_urb )) frac_urb (:) = 0.0
  if (allocated(frac_land)) frac_land(:) = 0.0

  if (isfcl > 0) then
     do j = 1,jlandflux(1)%jend(1)
        ilf = jlandflux(1)%ilandflux(j)

        iw  = landflux(ilf)%iw   ! global index
        iwl = landflux(ilf)%iwls ! global index

        if (iparallel == 1) then
           iw  = itabg_w (iw )%iw_myrank  ! local index
           iwl = itabg_wl(iwl)%iwl_myrank ! local index
        endif
   
        if (allocated(frac_urb)) then
           if (any(land%leaf_class(iwl) == (/ 19, 21 /))) then
              frac_urb(iw) = frac_urb(iw) + landflux(ilf)%arf_atm
           endif
        endif
       
        if (allocated(frac_land)) then
           frac_land(iw) = frac_land(iw) + landflux(ilf)%arf_atm
        endif
     enddo
  endif

! Store the fraction of the total surface that intersects with each layer

  do iw = 2, mwa
     if (lsw(iw) == 1) then
        frac_sfc(1,iw) = 1.0
     else
        do k = 1, lsw(iw)
           frac_sfc(k,iw) = (arw(lpw(iw)+k-1,iw) - arw(lpw(iw)+k-2,iw)) &
                          / arw(lpw(iw)+lsw(iw)-1,iw)
        enddo
     endif
  enddo

! Initialize PBL height and some PBL quantities
    
  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     ustar(iw) = 0.2
     wstar(iw) = 0.0
     wtv0 (iw) = 0.0

     do k = lpw(iw), mza-1
        thilv (k) = thil(k,iw) * (1.0 + eps_virt * sh_v(k,iw))
     enddo

     call acm2_pblhgt( ustar(iw), wstar(iw), wtv0(iw), lpw(iw), mza-2, lsw(iw), &
                       frac_sfc(:,iw), thilv, vxe(:,iw), vye(:,iw), vze(:,iw),  &
                       kpblh(iw), pblh(iw)                                      )
  enddo

end subroutine pbl_init








