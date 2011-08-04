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
Module lphys_interface

Interface

subroutine lphysiol_full(T_L,  &
     e_A,  &
     C_A,  &
     PAR,  &
     rb,  &
     adens,  &
     A_open,  &
     A_cl,  &
     rsw_open,  &
     rsw_cl,  &
     pft,  &
     prss,  &
     leaf_resp,  &
     green_leaf_factor,  &
     leaf_aging_factor, &
     old_st_data)

  use c34constants
  use pft_coms, only: D0, cuticular_cond, dark_respiration_factor, stomatal_slope,  &
       quantum_efficiency, photosyn_pathway, Vm0, Vm_low_temp
  use ed_structure_defs
  use ed_options, only: istoma_scheme

  implicit none

  real, intent(in) :: T_L
  real, intent(in) :: e_A
  real, intent(in) :: C_A
  real, intent(in) :: PAR
  real, intent(in) :: rb
  real, intent(in) :: adens
  real, intent(out) :: A_open
  real, intent(out) :: A_cl
  real, intent(out) :: rsw_open
  real, intent(out) :: rsw_cl
  integer, intent(in) :: pft
  real, intent(in) :: prss
  real, intent(in) :: green_leaf_factor
  real, intent(in) :: leaf_aging_factor
  real, intent(out) :: leaf_resp
  type(stoma_data), intent(inout) :: old_st_data

end subroutine lphysiol_full

subroutine c3solver(met,apar,gsdata,sol,ilimit)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(solution), intent(inout) :: sol
  integer, intent(out) :: ilimit

end subroutine c3solver

subroutine c4solver(met,apar,gsdata,sol,ilimit)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(solution), intent(inout) :: sol
  integer, intent(out) :: ilimit

end subroutine c4solver

subroutine setapar_c3(gsdata,met,apar,i)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  integer, intent(in) :: i
end subroutine setapar_c3

subroutine setapar_c4(gsdata,met,apar,i)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  integer, intent(in) :: i
end subroutine setapar_c4

real function aflux_c3(apar,gsw,ci)
  use c34constants
  implicit none
  type(glim), intent(in) :: apar
  real, intent(in) :: gsw
  real, intent(in) :: ci

end function aflux_c3

real function aflux_c4(apar,met,gsw)
  use c34constants
  implicit none
  type(glim), intent(in) :: apar
  real, intent(in) :: gsw
  type(metdat), intent(in) :: met

end function aflux_c4

real function csc_c4(apar,met,gsw)
  use c34constants
  implicit none
  type(metdat), intent(in) :: met
  type(glim), intent(in) :: apar
  real, intent(in) :: gsw

end function csc_c4

real function csc_c3(met,a)
  use c34constants
  implicit none
  type(metdat), intent(in) :: met
  real, intent(in) :: a

end function csc_c3

real function residual_c3(gsdata,met,apar,x)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  real, intent(in) :: x
end function residual_c3


real function residual_c4(gsdata,met,apar,x)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  real, intent(in) :: x
end function residual_c4

subroutine zbrak_c3(gsdata,met,apar,x1,x2,n,xb1,xb2,nb)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  real, intent(in) :: x1
  real, intent(in) :: x2
  integer, intent(in) :: n
  integer, intent(in) :: nb
  real, dimension(nb), intent(out) :: xb1,xb2
end subroutine zbrak_c3

subroutine zbrak_c4(gsdata,met,apar,x1,x2,n,xb1,xb2,nb)
  use c34constants
  implicit none

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  real, intent(in) :: x1
  real, intent(in) :: x2
  integer, intent(in) :: n
  integer, intent(in) :: nb
  real, dimension(nb), intent(out) :: xb1,xb2
end subroutine zbrak_c4

real function zbrent_c3(gsdata,met,apar,x1,x2,tol,success_flag)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  real, intent(in) :: x1
  real, intent(in) :: x2
  real, intent(in) :: tol
  integer, intent(inout) :: success_flag
end function zbrent_c3

real function zbrent_c4(gsdata,met,apar,x1,x2,tol,success_flag)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  real, intent(in) :: x1
  real, intent(in) :: x2
  real, intent(in) :: tol
  integer, intent(inout) :: success_flag
end function zbrent_c4

subroutine solve_closed_case_c3(gsdata,met,apar,sol,ilimit)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(solution), intent(inout) :: sol
  integer, intent(in) :: ilimit

end subroutine solve_closed_case_c3

subroutine solve_closed_case_c4(gsdata,met,apar,sol,ilimit)
  use c34constants
  implicit none

  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(solution), intent(inout) :: sol
  integer, intent(in) :: ilimit

end subroutine solve_closed_case_c4

subroutine solve_open_case_c3(gsdata,met,apar,sol,ilimit,success_flag)
  use c34constants
  implicit none
  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(solution), intent(inout) :: sol
  integer, intent(in) :: ilimit
  integer, intent(inout) :: success_flag
end subroutine solve_open_case_c3

subroutine solve_open_case_c4(gsdata,met,apar,sol,ilimit,success_flag)
  use c34constants
  implicit none
  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(solution), intent(inout) :: sol
  integer, intent(in) :: ilimit
  integer, intent(inout) :: success_flag
end subroutine solve_open_case_c4

subroutine closed2open(sol,ilimit)
  use c34constants
  implicit none
  type(solution), intent(inout) :: sol
  integer, intent(in) :: ilimit
end subroutine closed2open

real function quad4ci(gsdata,met,apar,x,success)
  use c34constants
  implicit none
  type(glim), intent(in) :: apar
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  real, intent(in) :: x
  integer, intent(out) :: success
end function quad4ci

subroutine prep_lphys_solution(photosyn_pathway, Vm0, met, Vm_low_temp,  &
     leaf_aging_factor, green_leaf_factor, leaf_resp, gsdata, apar)
  use c34constants
  implicit none
  integer, intent(in) :: photosyn_pathway
  real, intent(in) :: Vm0
  real, intent(in) :: Vm_low_temp
  type(metdat), intent(in) :: met
  real, intent(in) :: leaf_aging_factor
  real, intent(in) :: green_leaf_factor
  real, intent(out) :: leaf_resp
  type(farqdata), intent(in) :: gsdata
  type(glim), intent(inout) :: apar
end subroutine prep_lphys_solution

subroutine exact_lphys_solution(photosyn_pathway, met, apar, gsdata, sol,   &
     ilimit)
  use c34constants
  implicit none
  integer, intent(in) :: photosyn_pathway
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(solution), intent(inout) :: sol
  integer, intent(out) :: ilimit
end subroutine exact_lphys_solution

subroutine store_exact_lphys_solution(old_st_data, met, prss,   &
     leaf_aging_factor, green_leaf_factor, sol, ilimit, gsdata, apar,  &
     photosyn_pathway, Vm0, Vm_low_temp)
  use c34constants
  implicit none
  integer, intent(in) :: photosyn_pathway
  type(stoma_data), intent(inout) :: old_st_data
  type(metdat), intent(inout) :: met
  real, intent(in) :: prss
  real, intent(in) :: leaf_aging_factor
  real, intent(in) :: green_leaf_factor
  type(solution), intent(in) :: sol
  integer, intent(in) :: ilimit
  type(farqdata), intent(in) :: gsdata
  type(glim), intent(inout) :: apar
  real, intent(in) :: Vm0
  real, intent(in) :: Vm_low_temp
end subroutine store_exact_lphys_solution

subroutine fill_lphys_sol_exact(A_open, rsw_open, A_cl, rsw_cl, sol, adens)
  use c34constants
  implicit none
  real, intent(out) :: A_open
  real, intent(out) :: A_cl
  real, intent(out) :: rsw_open
  real, intent(out) :: rsw_cl
  type(solution), intent(in) :: sol
  real, intent(in) :: adens
end subroutine fill_lphys_sol_exact

subroutine fill_lphys_sol_approx(gsdata, met, apar, old_st_data, sol,   &
     A_cl, rsw_cl, adens, rsw_open, A_open, photosyn_pathway, prss)
  use c34constants
  implicit none
  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(stoma_data), intent(in) :: old_st_data
  type(solution), intent(inout) :: sol
  real, intent(out) :: A_cl
  real, intent(out) :: A_open
  real, intent(out) :: rsw_cl 
  real, intent(out) :: rsw_open
  real, intent(in) :: adens
  integer, intent(in) :: photosyn_pathway
  real, intent(in) :: prss
end subroutine fill_lphys_sol_approx

subroutine apply_disturbances(cs)
  use ed_structure_defs
  implicit none
  type(site), target :: cs
end subroutine apply_disturbances

subroutine insert_survivors(np, cp, q, area_fac)
  use ed_structure_defs
  implicit none
  type(patch), target :: np
  type(patch)         :: cp
  integer, intent(in) :: q
  real,    intent(in) :: area_fac
end subroutine insert_survivors

subroutine plant_patch(np, pft, density, height_factor)
  use ed_structure_defs
  implicit none
  type(patch), target :: np
  integer, intent(in) :: pft
  real,    intent(in) :: density
  real,    intent(in) :: height_factor
end subroutine plant_patch

subroutine split_cohorts(cp)
  use ed_structure_defs
  implicit none
  type(patch),  target  :: cp
  type(cohort), pointer :: cc
  type(cohort), pointer :: copyc
end subroutine split_cohorts

subroutine insert_cohort(pcc,cp)
  use ed_structure_defs
  implicit none
  type(patch)           :: cp
  type(cohort), target  :: pcc
end subroutine insert_cohort

subroutine fuse_2_patches(p1, p2)
  use ed_structure_defs
  implicit none
  type(patch), target :: p1
  type(patch), target :: p2
end subroutine fuse_2_patches

subroutine apply_forestry(cs, year)
  use ed_structure_defs
  implicit none
  type(site), target  :: cs
  integer, intent(in) :: year
end subroutine apply_forestry

subroutine reproduction(cs, month)
  use ed_structure_defs
  implicit none
  type(site), target  :: cs
  integer, intent(in) :: month
end subroutine reproduction

End Interface

End Module lphys_interface

