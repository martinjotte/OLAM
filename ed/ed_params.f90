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
subroutine load_ed_ecosystem_params()

use pft_coms, only: include_pft, n_pft

implicit none

call initialize_canopy_radiation_params()

! PLANT FUNCTIONAL TYPES (PFTs):
! 1 - C4 grass
! 2 - early tropical
! 3 - mid tropical
! 4 - late tropical
! 5 - C3 grass
! 6 - northern pines
! 7 - southern pines
! 8 - late conifers
! 9 - early temperate deciduous
! 10 - mid temperate deciduous
! 11 - late temperate deciduous

! Include all plant functional types
include_pft(1:4) = 1
include_pft(5:11) = 0

call initialize_pft_photo_params()
call initialize_pft_resp_params()
call initialize_pft_mort_params()
call initialize_pft_alloc_params()
call initialize_pft_nitro_params()
call initialize_pft_leaf_params()
call initialize_pft_repro_params()

call initialize_pft_derived_params()

return
end subroutine load_ed_ecosystem_params

!=================================================================
subroutine initialize_canopy_radiation_params()

use canopy_radiation_coms, only: leaf_reflect_nir,leaf_trans_nir,  &
     leaf_scatter_nir,leaf_reflect_vis_temperate,leaf_trans_vis_temperate, &
     leaf_scatter_vis,leaf_reflect_vis_tropics, leaf_trans_vis_tropics,  &
     diffuse_backscatter_vis, diffuse_backscatter_nir, emis_v

use pft_coms, only: n_pft, phenology

implicit none

integer :: ipft
real :: leaf_scatter_vis_temperate
real :: leaf_scatter_vis_tropics
real :: diffuse_backscatter_vis_temperate
real :: diffuse_backscatter_vis_tropics

leaf_scatter_nir = leaf_reflect_nir + leaf_trans_nir

leaf_scatter_vis_temperate = leaf_reflect_vis_temperate +   &
     leaf_trans_vis_temperate

leaf_scatter_vis_tropics = leaf_reflect_vis_tropics +   &
     leaf_trans_vis_tropics

diffuse_backscatter_vis_temperate = (2.0 * leaf_reflect_vis_temperate -   &
     leaf_trans_vis_temperate) / ( 3.0 * leaf_scatter_vis_temperate )

diffuse_backscatter_vis_tropics = (2.0 * leaf_reflect_vis_tropics -   &
     leaf_trans_vis_tropics) / ( 3.0 * leaf_scatter_vis_tropics )

diffuse_backscatter_nir = (2.0 * leaf_reflect_nir -   &
     leaf_trans_nir) / ( 3.0 * leaf_scatter_nir )

leaf_scatter_vis(1:4) = leaf_scatter_vis_tropics
leaf_scatter_vis(5:11) = leaf_scatter_vis_temperate

diffuse_backscatter_vis(1:4) = diffuse_backscatter_vis_tropics
diffuse_backscatter_vis(5:11) = diffuse_backscatter_vis_temperate

emis_v(1) = 0.96d0
emis_v(2:4) = 0.95d0
emis_v(5) = 0.96d0
emis_v(6:8) = 0.97d0
emis_v(9:11) = 0.95d0

return
end subroutine initialize_canopy_radiation_params

!=====================================================================
subroutine initialize_pft_photo_params()

use pft_coms, only: D0, Vm_low_temp, Vm0, stomatal_slope, cuticular_cond, &
     quantum_efficiency, photosyn_pathway

implicit none

D0 = 0.01 ! same for all PFTs

Vm_low_temp(1:4) = 5.0 ! tropical PFTs
Vm_low_temp(5:11) = 6.2 ! temperate PFTs

Vm0(1) = 12.5
Vm0(2) = 18.8
Vm0(3) = 12.5
Vm0(4) = 6.25
Vm0(5) = 18.3
Vm0(6) = 24.8
Vm0(7) = 24.8
Vm0(8) = 9.91
Vm0(9) = 29.4
Vm0(10) = 25.1
Vm0(11) = 10.0

stomatal_slope(1) = 10.0
stomatal_slope(2:4) = 8.0
stomatal_slope(5:11) = 6.19

cuticular_cond(1:5) = 10000.0
cuticular_cond(6:8) = 1000.0
cuticular_cond(9:11) = 20000.0

quantum_efficiency(1) = 0.06
quantum_efficiency(2:11) = 0.08

photosyn_pathway(1) = 4
photosyn_pathway(2:11) = 3

return
end subroutine initialize_pft_photo_params

!=====================================================================
subroutine initialize_pft_resp_params()

use pft_coms, only: growth_resp_factor, leaf_turnover_rate,   &
     dark_respiration_factor, storage_turnover_rate,  &
     root_respiration_factor

implicit none

growth_resp_factor(1:5) = 0.333
growth_resp_factor(6:8) = 0.5325
growth_resp_factor(9:11) = 0.0

leaf_turnover_rate(1) = 2.0
leaf_turnover_rate(2) = 1.0
leaf_turnover_rate(3) = 0.5
leaf_turnover_rate(4) = 0.333
leaf_turnover_rate(5) = 2.0
leaf_turnover_rate(6:8) = 0.333
leaf_turnover_rate(9:11) = 1.0

dark_respiration_factor(1) = 0.04
dark_respiration_factor(2:4) = 0.02
dark_respiration_factor(5) = 0.04
dark_respiration_factor(6:11) = 0.02

storage_turnover_rate(1:8) = 0.0
storage_turnover_rate(9:11) = 0.428

root_respiration_factor = 0.528

return
end subroutine initialize_pft_resp_params

!=====================================================================
subroutine initialize_pft_mort_params()

use pft_coms, only: mort1, mort2, mort3, seedling_mortality, treefall_s_gtht,  &
     treefall_s_ltht, plant_min_temp

implicit none

mort1(1:4) = 10.0
mort1(5:11) = 1.0

mort2 = 20.0

mort3(1) = 0.037
mort3(2) = 0.037
mort3(3) = 0.019
mort3(4) = 0.0
mort3(5) = 0.066
mort3(6) = 0.0034
mort3(7) = 0.0043
mort3(8) = 0.0024
mort3(9) = 0.0061
mort3(10) = 0.0038
mort3(11) = 0.0043

seedling_mortality = 0.95

treefall_s_gtht = 0.0

treefall_s_ltht(1) = 0.25
treefall_s_ltht(2:4) = 0.1
treefall_s_ltht(5) = 0.25
treefall_s_ltht(6:11) = 0.1

plant_min_temp(1:4) = 0.0
plant_min_temp(5:6) = -80.0
plant_min_temp(7) = -10.0
plant_min_temp(8) = -60.0
plant_min_temp(9) = -80.0
plant_min_temp(10:11) = -20.0

return
end subroutine initialize_pft_mort_params

!=====================================================================
subroutine initialize_pft_alloc_params()

use pft_coms, only: rho, SLA, q, qsw, hgt_min, b1Ht, b2Ht, b1Bs,  &
     b2Bs, b1Bl, b2Bl, C2B, leaf_turnover_rate

implicit none

rho(1:2) = 0.53
rho(3) = 0.71
rho(4) = 0.9
rho(5:11) = 0.0

SLA(1:4) = 10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(1:4)))) * C2B * 0.1
SLA(5) = 22.0
SLA(6) = 6.0
SLA(7) = 9.0
SLA(8) = 10.0
SLA(9) = 30.0
SLA(10) = 24.2
SLA(11) = 60.0

q(1:5) = 1.0
q(6:8) = 0.417
q(9:11) = 1.33

qsw = SLA / 3900.0

hgt_min = 1.5

b1Ht(1:4) = 0.0
b1Ht(5) = 0.4778
b1Ht(6) = 27.14
b1Ht(7) = 27.14
b1Ht(8) = 22.79
b1Ht(9) = 22.6799
b1Ht(10) = 25.18
b1Ht(11) = 23.3874

b2Ht(1:4) = 0.0
b2Ht(5) = -0.75
b2Ht(6) = -0.03884
b2Ht(7) = -0.03884
b2Ht(8) = -0.04445 
b2Ht(9) = -0.06534
b2Ht(10) = -0.04964
b2Ht(11) = -0.05404

b1Bl(1:4) = 0.0
b1Bl(5) = 0.08
b1Bl(6) = 0.024
b1Bl(7) = 0.024
b1Bl(8) = 0.0454
b1Bl(9) = 0.0129
b1Bl(10) = 0.048
b1Bl(11) = 0.017

b2Bl(1:4) = 0.0
b2Bl(5) = 1.0
b2Bl(6) = 1.899
b2Bl(7) = 1.899
b2Bl(8) = 1.6829
b2Bl(9) = 1.7477
b2Bl(10) = 1.455
b2Bl(11) = 1.731

b1Bs(1:4) = 0.0 
b1Bs(5) = 1.0e-5
b1Bs(6) = 0.147
b1Bs(7) = 0.147
b1Bs(8) = 0.1617
b1Bs(9) = 0.02648
b1Bs(10) = 0.1617
b1Bs(11) = 0.235

b2Bs(1:4) = 0.0
b2Bs(5) = 1.0
b2Bs(6) = 2.238
b2Bs(7) = 2.238
b2Bs(8) = 2.1536
b2Bs(9) = 2.95954
b2Bs(10) = 2.4572
b2Bs(11) = 2.2518

return
end subroutine initialize_pft_alloc_params

!=====================================================================
subroutine initialize_pft_nitro_params()

use pft_coms, only: c2n_leaf, Vm0, SLA, water_conductance

implicit none

c2n_leaf = 1000.0 / ((0.11289 + 0.12947 * Vm0) * SLA)

water_conductance(1:4) = 0.001904
water_conductance(5:11) = 0.004639

return
end subroutine initialize_pft_nitro_params

!=====================================================================
subroutine initialize_pft_leaf_params()

use pft_coms, only: phenology, clumping_factor, leaf_width

implicit none

phenology(1:5) = 1
phenology(6:8) = 0
phenology(9:11) = 2

clumping_factor(1) = 1.0
clumping_factor(2:4) = 0.735
clumping_factor(5) = 0.84
clumping_factor(6:8) = 0.735
clumping_factor(9:11) = 0.84

leaf_width(1:4) = 0.2
leaf_width(5:11) = 0.05

return
end subroutine initialize_pft_leaf_params

!=====================================================================
subroutine initialize_pft_repro_params()

use pft_coms, only: r_fract, seed_rain, nonlocal_dispersal, repro_min_h

implicit none

r_fract(1) = 1.0
r_fract(2:4) = 0.3
r_fract(5) = 1.0
r_fract(6:11) = 0.3

seed_rain = 0.01

nonlocal_dispersal(1:5) = 1.0
nonlocal_dispersal(6:7) = 0.766
nonlocal_dispersal(8) = 0.001
nonlocal_dispersal(9) = 1.0
nonlocal_dispersal(10) = 0.325
nonlocal_dispersal(11) = 0.074

repro_min_h(1:5) = 0.0
repro_min_h(6:11) = 5.0

return
end subroutine initialize_pft_repro_params

!=====================================================================
subroutine initialize_pft_derived_params()

use pft_coms, only: root_turnover_rate, c2n_leaf, max_dbh, b1Ht, b2Ht,  &
     n_pft, hgt_min, q, qsw, c2n_recruit, c2n_stem
use decomposition_coms, only: f_labile

implicit none

integer :: ipft
real :: dbh,h2dbh,bleaf,dbh2bl,bdead,dbh2bd,balive

root_turnover_rate(1) = 2.0
root_turnover_rate(2) = 1.0
root_turnover_rate(3) = 0.5
root_turnover_rate(4:5) = 0.333
root_turnover_rate(6:11) = 3.5169 * c2n_leaf(10) / c2n_leaf(6:11)  ! leaves and fine roots have the same c2n.

max_dbh(1) = 0.498
max_dbh(2:4) = 68.31
max_dbh(5:11) = log(1.0-(0.9*(b1Ht(5:11)+1.3)-1.3)/b1Ht(5:11))/b2Ht(5:11)

do ipft = 1,n_pft
   dbh = h2dbh(hgt_min(ipft),ipft)
   bleaf = dbh2bl(dbh,ipft)
   bdead = dbh2bd(dbh,hgt_min(ipft),ipft)
   balive = bleaf * (1.0 + q(ipft) + qsw(ipft) * hgt_min(ipft))
   c2n_recruit(ipft) = (balive + bdead) / (balive * (f_labile(ipft) /   &
        c2n_leaf(ipft) + (1.0 - f_labile(ipft)) / c2n_stem) +   &
        bdead/c2n_stem)
enddo

end subroutine initialize_pft_derived_params

