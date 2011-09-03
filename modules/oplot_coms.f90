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
Module oplot_coms

use max_dims, only: maxnplt, maxpltfiles, pathlen
implicit none

private :: maxnplt, maxpltfiles, pathlen

Type oplot_vars

! Variables copied or computed from $MODEL_PLOT namelist

   integer :: nplt       ! Number of fields to plot
   integer :: nplt_files ! Number of files to plot from
   integer :: plttype    ! Plot output type (ncgm, ps, or pdf)
   integer :: pltorient  ! Landscape or portrait
   integer :: vec_maxmrl ! Highest mesh level to plot wind vectors

   real :: frqplt    ! Time interval between plots
   real :: dtvec     ! Scaling time (seconds) for plotted velocity vectors
   real :: headspeed ! Scaling speed [m/s] for velocity vector arrow head length
   real :: stemlength ! Map distance [m] of windbarb stem for size scaling
   real :: xmin      ! Frame window minimum x coord in x/y/z or lat/lon/z space
   real :: xmax      ! Frame window maximum x coord in x/y/z or lat/lon/z space
   real :: ymin      ! Frame window minimum y coord in x/y/z or lat/lon/z space
   real :: ymax      ! Frame window maximum y coord in x/y/z or lat/lon/z space
   real :: plat3     ! Plot projection pole latitude
   real :: plon3     ! Plot projection pole longitude
   real :: conelat   ! Plot cone center latitude
   real :: conelon   !  Plot cone center longitude
   real :: coneang   ! Plot cone angle
   real :: viewazim  ! Plot view azimuth for vertical cross sections

   character(pathlen) :: pltname
   character(pathlen) :: plt_files(maxpltfiles)
   character(20)      :: fldname(maxnplt) ! Name of plotted field

   integer :: icolortab(maxnplt) ! Color table number for plotting given field

   character(len=1), dimension(maxnplt) ::  &
       projectn              & ! Plot projection & cross section ['P','E','O','C','X','Y']
      ,contrtyp = 'N'        & ! Contour type ['T','F','L','O']
      ,prtval = 'N'          & ! Flag to print value ['P','N']
      ,pltindx = 'N'         & ! Print ITABLE index ['I','J','N']
      ,vectbarb = 'N'        & ! Plot vectors or windbarbs ['B','U','V','v']
      ,pltgrid = 'N'         & ! Plot grid lines ['G','N']
      ,pltgrid_landsea = 'N' & ! Plot land/sea grid lines ['g','N']
      ,pltdualgrid = 'N'     & ! Plot grid lines ['D','N']
      ,pltborder = 'N'       & ! Plot border, border + ticks + labels ['b','t','N']
      ,maptyp = 'N'          & ! Type of map plot ['M','m','N']
      ,pltcone = 'N'         & ! Flag to plot cone circle ['C','N']
      ,windowin = 'N'        & ! Flag to window in ['W','N']
      ,frameoff = 'N'        & ! Flag to suppress frame call ['f','N']
      ,panel = 'N'           & ! Panel number if reduced-size plot ['1','2','3','4']
      ,colorbar = 'N'        & ! Print colorbar ['c','N']
      ,labelbar = 'N'        & ! Print title, title + info block ['t','i','N']
      ,pltlev = 'N'          & ! Plot on const press level or near sfc ['p','s','N']
      ,ext = 'N'               ! Flag indicating "external" plot field ['e','N']

   real :: slabloc(maxnplt) ! Z-coord of plot slab in x/y/z or lat/lon/z space

! Plot layout coordinates

   real :: h1  ! Frame window minimum x coord in window space (range 0 to 1)
   real :: h2  ! Frame window maximum x coord in window space (range 0 to 1)
   real :: v1  ! Frame window minimum y coord in window space (range 0 to 1)
   real :: v2  ! Frame window maximum y coord in window space (range 0 to 1)
   real :: hp1 ! Panel window minimum x coord in window space (range 0 to 1)
   real :: hp2 ! Panel window maximum x coord in window space (range 0 to 1)
   real :: vp1 ! Panel window minimum y coord in window space (range 0 to 1)
   real :: vp2 ! Panel window maximum y coord in window space (range 0 to 1)
   real :: fx1 ! plot frame left side x coord
   real :: fx2 ! plot frame right side x coord
   real :: fy1 ! plot frame bottom y coord
   real :: fy2 ! plot frame top y coord

   real :: cbx1 ! colorbar left side x coord
   real :: cbx2 ! colorbar right side x coord
   real :: cblx ! colorbar label left end x coord

   real :: fnamey ! field name y coord

   real :: xlaby  ! x-axis label y coord
   real :: xtlaby ! x-axis tick label y coord

   real :: ylabx  ! y-axis label x coord
   real :: ytlabx ! y-axis tick label right end x coord

   real :: timex  ! elapsed time (sec) left end x coord
   real :: timsy  ! elapsed time (sec) y coord
   real :: timdy  ! elapsed time (day) y coord

   real :: slabx  ! slab left end x coord
   real :: slaby  ! slab y coord
   real :: sminy  ! slab min y coord
   real :: smaxy  ! slab max y coord

! Variables set in default subroutine

   integer :: loopplot  ! Flag to plot DO loop indices [1 = YES, 0 = NO]
   integer :: icigrnd   ! "Underground" color index
   integer :: iplotback ! Flag indicating whether plotback has been called since
                        !    last frame call [1 = YES, 0 = 'NO']

   real :: psiz   ! Plot letter/number size
   real :: vsprd  ! Plot vertical distance to neighbor indices

! Variables set in oplot_lib subroutine

   integer :: ifill  ! Flag set to 1 if contrtyp = 'FILL'; 0 otherwise

   real :: fldval_min  ! Minimum value in slab plot
   real :: fldval_max  ! Maximum value in slab plot
   real :: fldvalv_min ! Minimum value in vector slab plot
   real :: fldvalv_max ! Maximum value in vector slab plot

   character(len=40) :: units  ! Units of plotted field
   character(len=40) :: stagpt ! Location in staggered grid of value to be plotted
   character(len=40) :: dimens ! Dimensionality of array to be plotted
   character(len=40) :: label  ! Label of field to be plotted

! Variables used to save the current graphics attributes

   integer :: i_att(14)

   real :: r_att(7)

   real, allocatable :: extfld(:,:)

   character(len=20) :: extfldname = 'N' ! Name of external plotted field

End Type

type (oplot_vars), save :: op

! Miscellaneous

   real :: xepc(2),yepc(2),zepc(2)  ! Earth coords of 2 pts of intersection of
                                    ! plot cone, earth surface, and sides of
                                    ! triangular grid column

End Module
