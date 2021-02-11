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
subroutine lake_init_atm()

  use mem_lake,    only: lake, mlake, omlake
  use lake_coms,   only: dt_lake
  use mem_basic,   only: rho, press, vxe, vye, vze, tair, rr_v, theta
  use mem_micro,   only: rr_c
  use misc_coms,   only: io6, time8, s1900_sim, iparallel, isubdomain, runtype
  use mem_ijtabs,  only: itabg_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_para,    only: myrank
  use consts_coms, only: t00, cliq, alli, p00i, grav, rocp
  use therm_lib,   only: rhovsl, rhovsil

  implicit none

  integer :: iw
  integer :: kw
  integer :: ilake, iwsfc, j

  real :: vels, psfc, laketemp

  ! Initialize lake quantities that do not depend on atmospheric conditions

  !$omp parallel do private (iwsfc)
  do ilake = 2,mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

  ! Initialize lake depth and canopy depth

     lake%depth(ilake) = 20. ! SET TO ROUGHLY THE SUMMERTIME THERMOCLINE DEPTH; EVENTUALLY, USE MULTILEVEL LAKE MODEL
     sfcg%can_depth(iwsfc) = 20. * max(1.,.025 * dt_lake)

  enddo  ! ilake
  !$omp end parallel do

  ! End of initialization that does not depend on atmospheric conditions

  if (runtype /= "INITIAL") return

  ! Sea quantities that are initialized only on 'INITIAL' run

  !$omp parallel do private (iwsfc,laketemp)
  do ilake = 2,mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply initial atmospheric properties to lake "canopy"

     sfcg%cantemp(iwsfc) = sfcg%airtheta(iwsfc) * (sfcg%prss(iwsfc) * p00i)**rocp
     sfcg%canrrv (iwsfc) = sfcg%airrrv(iwsfc)
     sfcg%ustar  (iwsfc) = 0.1
     sfcg%ggaer  (iwsfc) = 0.0
     sfcg%wthv   (iwsfc) = 0.0

     ! Lake water quantities (eventually enable lake_energy to be initialized from spin-up simulation)

     laketemp                 = sfcg%cantemp(iwsfc)
     lake%lake_energy (ilake) = (laketemp - t00) * cliq + alli
     lake%surface_srrv(ilake) = rhovsl(laketemp - t00) / sfcg%rhos(iwsfc)
     sfcg%rough       (iwsfc) = .001
     sfcg%ustar       (iwsfc) = 0.1
     sfcg%ggaer       (iwsfc) = 0.0
     sfcg%wthv        (iwsfc) = 0.0

  enddo
  !$omp end parallel do

end subroutine lake_init_atm
