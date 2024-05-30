Module oplot_coms

  use max_dims,    only: maxnplt, maxpltfiles, pathlen
  use consts_coms, only: r8
  implicit none

  private :: maxnplt, maxpltfiles, pathlen

! Variables copied or computed from $MODEL_PLOT namelist

  Type oplot_vars

     integer :: nplt       ! Number of fields to plot
     integer :: nplt_files ! Number of files to plot from
     integer :: plttype    ! Plot output type (ncgm, ps, or pdf)
     integer :: pltorient  ! Landscape or portrait
     integer :: vec_maxmrl ! Highest mesh level to plot wind vectors

     integer :: mapcolor = 10  ! Color of map outlines
     integer :: llcolor  = 10  ! Color of lat/lon lines

     real(r8) :: frqplt    ! Time interval between plots

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
     real :: conelon   ! Plot cone center longitude
     real :: coneang   ! Plot cone angle
     real :: viewazim  ! Plot view azimuth for vertical cross sections

     real :: pxe
     real :: pye
     real :: pze
     real :: sinplat
     real :: cosplat
     real :: sinplon
     real :: cosplon

     character(pathlen) :: pltname
     character(pathlen) :: plt_files(maxpltfiles)

     character(30) :: fldname(maxnplt) ! Name of plotted field
     character(10) :: prtval_size      ! Printed char size ['small', 'medium', large]

     integer :: icolortab(maxnplt) ! Color table number for plotting given field

     character(len=1), dimension(maxnplt) ::  &
          projectn    = ' ' & ! Plot projection & cross section ['L','P','G','O','Z','C','V']
         ,contrtyp    = 'N' & ! Contour type ['T','F','L','O']
         ,prtval      = 'N' & ! Flag to print value ['P']
         ,pltindx1    = 'N' & ! Print ATM grid index ['I','J']
         ,pltindx2    = 'N' & ! Print SFC grid index ['i','j']
         ,vectbarb    = 'N' & ! Plot vectors or windbarbs for winds ['B','V','w']
         ,vectsea     = 'N' & ! Plot vectors for ocean currents or wave propagation ['Y','y','z']
         ,pltgrid     = 'N' & ! Plot grid lines ['G']
         ,pltgrid_sfc = 'N' & ! Plot surface grid lines ['g']
         ,pltdualgrid = 'N' & ! Plot grid lines ['D']
         ,pltborder   = 'N' & ! Plot border, +ticks+outer_axlabs, +inner_axlabs ['b','t','a']
         ,labelbar    = 'N' & ! Print field name, field name + info block ['n','o']
         ,colorbar    = 'N' & ! Print colorbar, no colorbar but reserve space ['c','r']
         ,maptyp      = 'N' & ! Flag to plot land, country, state outlines ['m']
         ,pltll       = 'N' & ! Flag to plot lat/lon lines ['l']
         ,pltcone     = 'N' & ! Flag to plot cone circle ['C']
         ,pltlev      = 'N' & ! Plot on const press level or near sfc ['p','s']
         ,windowin    = 'N' & ! Flag to window in ['W']
         ,ext         = 'N' & ! Flag indicating "external" plot field ['e']
         ,frameoff    = 'N' & ! Flag to suppress frame call ['f','h']
         ,panel       = 'N' & ! Panel number if reduced-size plot ['1' to '9']
         ,noundrg     = 'N'   ! Suppress masking underground cells ['u']

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

     logical :: has_med_res = .false. ! Indicates the neweer medium-resolution maps
                                      ! are available in the installed NCAR Graphics

     logical :: has_high_res = .false. ! Indicates the newest high-resolution maps
                                       ! are available in the installed NCAR Graphics
  End Type oplot_vars

  type (oplot_vars) :: op

! Miscellaneous

  integer :: iplt_file

  real :: xepc(2), yepc(2), zepc(2)  ! Earth coords of 2 pts of intersection of
                                     ! plot cone, earth surface, and sides of
                                     ! triangular grid column
End Module oplot_coms
