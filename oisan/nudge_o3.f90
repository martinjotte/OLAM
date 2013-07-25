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

subroutine nudge_prep_o3(iaction, o_ozone)

  use mem_nudge,  only: ozone_obsp, ozone_obsf
  use mem_grid,   only: mza, mwa
  use mem_ijtabs, only: jtab_w, jtw_init
  use misc_coms,  only: do_chem

  implicit none

  integer, intent(in) :: iaction
  real,    intent(in) :: o_ozone(mza,mwa)
  integer             :: j, iw, k

  if (.not. do_chem) return

! Swap future data time into past data time if necessary.

  if (iaction == 1) then

     !$omp parallel do private(iw,k)
     do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
        do k = 2, mza-1
           ozone_obsp(k,iw) = ozone_obsf(k,iw)
        enddo
     enddo
     !omp end parallel do

  endif

! Fill future observational arrays

  !$omp parallel do private(iw,k)
  do j = 1, jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
     do k = 2, mza-1
        ozone_obsf(k,iw) = o_ozone(k,iw)
     enddo
  enddo
  !omp end parallel do

end subroutine nudge_prep_o3

!===========================================================================

subroutine obs_nudge_o3()

  use mem_nudge,  only: ozone_obsp, ozone_obsf, o3nudflag, tnudi_o3, o3nudpress
  use mem_basic,  only: rho, press
  use mem_grid,   only: mza, lpw
  use mem_ijtabs, only: istp, jtab_w, mrl_begl, jtw_prog
  use mem_tend,   only: umt, vmt, thilt, sh_wt
  use isan_coms,  only: ifgfile, s1900_fg
  use cgrid_defn, only: ns_o3
  use misc_coms,  only: io6, time8, s1900_sim, do_chem
  use var_tables, only: scalar_tab

  implicit none

  integer :: k, j, iw, mrl, npoly, kbot
  real    :: tp, tf, tnudi

  ! Nudge ozone to observed data

  if (.not. do_chem) return

  ! Check whether it is time to nudge

  mrl = mrl_begl(istp)
  if (mrl < 1) return

  ! Time interpolation coefficients

  tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
            / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

  tp = 1. - tf

  ! Compute ozone tendency using point-by-point (non-spectral) information

  !$omp parallel do private(iw,k,kbot)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Find level above which we nudge ozone

     do k = mza-1, lpw(iw), -1
        if (press(k,iw) > o3nudpress) exit
     enddo

     kbot = max( k+1, lpw(iw) )

     do k = kbot, mza-1

        scalar_tab(ns_o3)%var_t(k,iw) = scalar_tab(ns_o3)%var_t(k,iw) + tnudi_o3 * rho(k,iw) * &
             ( (tp * ozone_obsp(k,iw) + tf * ozone_obsf(k,iw)) - scalar_tab(ns_o3)%var_p(k,iw) )

     enddo

  enddo
  !$omp end parallel do

end subroutine obs_nudge_o3
