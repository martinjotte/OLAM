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
Module canopy_radiation_coms

use mem_ed, only: n_pft

real, parameter :: mubar     = 1.0  ! Leaf angle distribution parameter (dimensionless).  Let mu' be the cosine of leaf angle and G(mu') be the distribution of mu'.  Then, mubar = (integral from 0 to 1) (d mu'   mu' / G(mu')).  See, for example, Dickinson 1983.

real, parameter :: Watts2Ein = 4.6e-6 ! Converts PAR radiation from  watts to Einsteins (units are Ein/watt of PAR).

real, parameter :: visible_fraction = 0.45 ! fraction of solar radiation in the PAR band.  Used when you don't know the direct/diffuse breakdown.

real :: visible_fraction_dir = 0.43 ! fraction of direct solar radiation in the PAR band.

real :: visible_fraction_dif = 0.52 ! fraction of diffuse solar radiation in the PAR band.  Used when you don't know the direct/diffuse breakdown.

! The following two values are based on Dickinson (1983)
real, parameter :: leaf_reflect_nir = 0.577  !  fraction of scattered NIR that is reflected.  

real, parameter :: leaf_trans_nir = 0.248  !  fraction of scattered NIR that is transmitted.

real :: leaf_scatter_nir ! sum of reflected plus scattered NIR

! The following two values are from Baldocchi et al.
real, parameter :: leaf_reflect_vis_temperate = 0.11  !  fraction of scattered PAR that is reflected.  Value obtained for temperate trees.

real, parameter :: leaf_trans_vis_temperate = 0.16  !  fraction of scattered PAR that is transmitted.  Value obtained for temperate trees.

! The following two values are from Poorter et al.
real, parameter :: leaf_reflect_vis_tropics = 0.062  !  fraction of scattered PAR that is reflected.  Value obtained for tropical trees.

real, parameter :: leaf_trans_vis_tropics = 0.028  !  fraction of scattered PAR that is transmitted.  Value obtained for tropical trees.

real, dimension(n_pft) :: leaf_scatter_vis  ! sum of reflected plus scattered PAR

real, dimension(n_pft) :: diffuse_backscatter_vis !  Backscatter parameter for diffuse PAR

real :: diffuse_backscatter_nir !  Backscatter parameter for diffuse NIR

real, parameter :: lai_min = 1.0e-5  ! Minimum LAI used in the radiative transfer scheme.

real(kind=8), dimension(n_pft) :: emis_v  ! Emissivity of the vegetation

End Module canopy_radiation_coms
