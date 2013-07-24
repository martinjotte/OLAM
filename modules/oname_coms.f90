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
Module oname_coms

   use max_dims,    only: nzgmax, maxsndg, maxgrds, maxisdirs, maxnplt,  &
                          maxpltfiles, maxngrdll, pathlen
   use consts_coms, only: r8

   Type oname_plot
      ! Derived type to hold the components of the plot specification fields
      character(20) :: fldname    = ''
      character(1)  :: projectn   = ''
      integer       :: icolortab  = 0
      character(20) :: pltspec2   = 'N'
      real          :: plotcoord1 = 0.0
      real          :: plotcoord2 = 0.0
      real          :: slabloc    = 0.0
      real          :: plotwid    = 0.0
      real          :: viewazim   = 0.0
   End Type oname_plot

   Type oname_vars

!!    RUNTYPE/NAME

      character(64) :: expnme  = 'OLAM test'
      character(16) :: runtype = ''

!!    SIMULATION ENDING TIME

      character(1) :: timeunit = ''
      real(r8)     :: timmax   = 0.0_r8
      
!!    START OF SIMULATION OR ISAN PROCESSING

      integer :: itime1  = 0
      integer :: idate1  = 0
      integer :: imonth1 = 0
      integer :: iyear1  = 0

!!    GRID SPECIFICATIONS

      integer :: mdomain  = 0
      integer :: meshtype = 0
      integer :: nzp      = 0
      integer :: nxp      = 0

      real(r8) :: dtlong = 0.0_r8
      real     :: deltax = 0.0
      
      integer :: ndz = 0
      
      real :: hdz(10) = 0.
      real :: dz (10) = 0.

      real :: zz(maxsndg) = 0.0

!!    NESTED GRID DEFINITION

      integer :: ngrids   = 1
      integer :: nconcave = 1

      integer :: ngrdll(maxgrds) = 0
      real    :: grdrad(maxgrds) = 0.0

      real    :: grdlat(maxgrds,maxngrdll) = 0.0
      real    :: grdlon(maxgrds,maxngrdll) = 0.0

!!    TIMESTEP RATIOS

      integer :: ndtrat (maxgrds) = 1
      integer :: nacoust(maxgrds) = 1
      
!!    VARIABLE INITIALIZATION INPUT

      integer            :: initial  = 0
      character(pathlen) :: zonclim  = '../etc/ZONAVG_CLIMATE'

!!    NUDGING PARAMETERS

      integer :: nudflag  = 0
      integer :: nudnxp   = 0
      real    :: tnudcent = 0.0

      integer :: o3nudflag  = 0
      real    :: o3nudpress = 150.0
      real    :: o3tnudcent = 86400.0
      
!!    GRID, HISTORY FILES

      integer :: ioutput      = 1
      integer :: ioutput_mavg = 1
      integer :: ioutput_davg = 1
      integer :: iclobber  = 0
      integer :: icompress = 0
      integer :: ipar_out  = 0
      integer :: iquiet    = 0
      real(r8):: frqstate  = 3600.0_r8

      character(pathlen) :: gridfile  = 'sfcfile/gridfile_0'
      character(pathlen) :: hfilin    = ''
      character(pathlen) :: hfilepref = 'hist/'

!!    TOPOGRAPHY INITIALIZATION

      integer            :: itopoflg      = 2
      character(pathlen) :: topo_database = ''

!!    MODEL/NUMERICAL OPTIONS

      integer :: naddsc      = 0
      integer :: icorflg     = 1
      integer :: ithil_monot = 0
      integer :: iwind_monot = 0
      integer :: iscal_monot = 0
      logical :: debug_fp    = .false.
      logical :: init_nans   = .false.

!!    RALEIGH FRICTION PARAMETERS

      real :: rayf_zmin    = 30000.0
      real :: rayf_distim  = 0.0
      real :: rayf_expon   = 1.0

      real :: rayfw_zmin   = 30000.0
      real :: rayfw_distim = 0.0
      real :: rayfw_expon  = 1.0
      
      real :: rayfdiv_zmin   = 30000.0
      real :: rayfdiv_distim = 0.0
      real :: rayfdiv_expon  = 1.0

!!    RADIATION PARAMETERIZATION PARAMETERS

      integer :: iswrtyp = 0
      integer :: ilwrtyp = 0
      real(r8):: radfrq  = 1800.0_r8
      
!!    CUMULUS PARAMETERIZATION PARAMETERS

      integer :: nqparm   (maxgrds) = 0
      integer :: nqparm_sh(maxgrds) = 0

      real(r8) :: confrq = 1800.0_r8
      real     :: wcldbs = 0.01

!!    EDDY DIFFUSION PARAMETERS

      integer :: idiffk(maxgrds) = 2

      real :: zkhkm(maxgrds) = 3.0
      real :: xkhkm(maxgrds) = 3.0
      real :: csx  (maxgrds) = 0.2
      real :: csz  (maxgrds) = 0.2
      real :: akmin(maxgrds) = 0.0
                        
!!    MICROPHYSICS PARAMETERS

      integer :: level  = 1
      integer :: icloud = 4
      integer :: idriz  = 4
      integer :: irain  = 2
      integer :: ipris  = 5
      integer :: isnow  = 2
      integer :: iaggr  = 2
      integer :: igraup = 2
      integer :: ihail  = 2
      integer :: iccnlev = 0
      
      real :: cparm  = 300.e6 ! [#/kg]
      real :: dparm  = 1.e1   ! [#/kg]
      real :: rparm  = .001   ! [m]
      real :: pparm  = 1.e5   ! [#/kg]
      real :: sparm  = .003   ! [m]
      real :: aparm  = .003   ! [m]
      real :: gparm  = .003   ! [m]
      real :: hparm  = .01    ! [m]
      real :: cnparm = .04e-6 ! [m]
      real :: gnparm = 3.0e-6 ! [m]
      
!!    SOUNDING SPECIFICATION

      integer :: nsndg   = 0
      integer :: ipsflg  = 1
      integer :: itsflg  = 0
      integer :: irtsflg = 3
      integer :: iusflg  = 0

      real :: hs    = 0.0
      real :: p_sfc = 1000.0

      real :: sounding(5,maxsndg) = 0.0

!!    LEAF/SEA VARIABLES

      integer :: isfcl    =  1
      integer :: nzg      = 11
      integer :: nzs      =  1

      integer :: ivegflg    = 0
      integer :: isoilflg   = 0
      integer :: ndviflg    = 0
      integer :: isstflg    = 0
      integer :: iseaiceflg = 0

      integer :: isoilstateinit = 0
      integer :: isoildepthflg  = 0
      integer :: iwatertabflg   = 0

      integer :: iupdndvi   = 0
      integer :: iupdsst    = 0
      integer :: iupdseaice = 0
      integer :: nvgcon     = 8
      integer :: nslcon     = 6

      real    :: seatmp = 280.0
      real    :: seaice =   0.0

      real :: slz(nzgmax) = (/ -1.00, -.85, -.70, -.60, -.50, -.40, &
           -.30, -.20, -.15, -.10, -.05, (0.0, i=12,nzgmax) /)

      real :: slmstr(nzgmax) = (/ .35, .35, .35, .35, .35, .35, .35, &
           .35, .35, .35, .35, (0.0, i=12,nzgmax) /)

      character(pathlen) :: landusefile = 'sfcfiles/landh'
      character(pathlen) :: seafile     = 'sfcfiles/seah'

      character(pathlen) :: veg_database    = ''
      character(pathlen) :: soil_database   = ''
      character(pathlen) :: ndvi_database   = ''
      character(pathlen) :: sst_database    = ''
      character(pathlen) :: seaice_database = ''

      character(pathlen) :: soilstate_db    = ''
      character(pathlen) :: soildepth_db    = ''
      character(pathlen) :: watertab_db     = ''

!!    ED MODEL VARIABLES

      character(pathlen) :: ed2_namelist = ''
      integer            :: ed2_active   = 0

!!    ISENTROPIC CONTROL

      character(pathlen) :: iapr(maxisdirs) = ''

!!    MODEL_PLOT VARIABLES

      integer :: nplt        = 0
      integer :: nplt_files  = 0
      integer :: nmavg_files = 0
      integer :: ndavg_files = 0
      real(r8):: frqplt     = 3600.0_r8
      real    :: dtvec      = 1200.0
      real    :: headspeed  = 3.0
      real    :: stemlength = 3.0e3
      integer :: plttype    = 0
      integer :: pltorient  = 0
      integer :: vec_maxmrl = maxgrds
      real    :: zplot_min = -1.0
      real    :: zplot_max = -1.0

      character(pathlen) :: pltname     = 'gmeta'
      character(10)      :: prtval_size = 'medium'

!!    THE LIST OF FILES TO PLOT FROM

      character(pathlen) ::  plt_files(maxpltfiles) = ''
      character(pathlen) :: mavg_files(maxpltfiles) = ''
      character(pathlen) :: davg_files(maxpltfiles) = ''

!!    THE ARRAYS OF FIELDS TO PLOT

      type(oname_plot) :: plotspecs(maxnplt)

!!    NCAR DCMIP 2012 PARAMETERS

      integer :: test_case

   End Type oname_vars

   type (oname_vars), save :: nl

   character(40) :: cmdlne_runtype = ''

End Module oname_coms
