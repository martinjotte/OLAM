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
  
  use mem_grid,   only: mwa, mza, lpw, lsw
  use misc_coms,  only: io6, idiffk
  use mem_tend,   only: thilt, sh_wt
  use mem_basic,  only: wc, rho
  use mem_turb,   only: hkm, vkm, vkh, sxfer_tk, sxfer_rk
  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog

  use smagorinsky

!$use omp_lib

  implicit none

  real, intent(inout) :: rhot(mza,mwa)
  integer, intent(in) :: mrl
  integer :: j, k, ka, iw, mrlw, ks

  
! Loop over all W/T points where PBL parameterization may be done

  call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,mrlw,ka,k,ks) 
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

     ! Select PBL scheme based on MRL of current IW column

     if (idiffk(mrlw) == 2) then
   
        ! Smagorinsky scheme

        call turb_k(iw, mrlw, rhot)

     endif

     ! Zero out sxfer arrays now that they have been transferred to the atm
     
     do ks = 1, lsw(iw)
        sxfer_tk(ks,iw) = 0.0
        sxfer_rk(ks,iw) = 0.0
     enddo
   
  enddo

end subroutine pbl_driver
