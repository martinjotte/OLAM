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
subroutine emis_database_read(iaction)

!!use mem_sea,     only: sea, itab_ws
!!
!!use sea_coms,    only: mms, mws, iupdsst, isstcyclic, nsstfiles,  &
!!                       fnames_sst, ctotdate_sst, s1900_sst,       &
!!                       isstfile, sst_database, isstflg
!!
use misc_coms,   only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                       time8, runtype, s1900_init, s1900_sim
!!
use consts_coms, only: erad, pio180
!!use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
!!use max_dims,    only: maxsstfiles
use mem_sflux, only: mseaflux, mlandflux, seaflux, landflux

use emis_defn, only: emis_files, nemisfiles, ctotdate_emis, s1900_emis,  &
                     iemisflg, iupdemis, iemiscyclic, iemisfile, emlays, &
                     nvars3d_emis, vname3d_emis, units3d_emis,           &
                     sea_emisf,  sea_emisp,  sea_emisc,                  &
                     land_emisf, land_emisp, land_emisc

use netcdf

implicit none

integer, intent(in) :: iaction

integer :: iemisy,iemism,iemisd,iemish
integer :: oemisy,oemism,oemisd,oemish
integer :: iyears, imonths, idates, ihours

!real, allocatable :: dato(:,:,:,:), data(:,:,:)
real, allocatable :: dato(:,:,:), data(:,:,:)

integer :: nio, njo, nis, njs, nie, nje
integer :: i, j, k, ii, jj
integer :: nx, ny
integer :: io1, io2, jo1, jo2
integer :: nf, n, nv
integer :: isf, ilf
integer :: ntimes, jtime
integer :: slen
integer :: ndims, idims(2)
integer :: indx

real :: wio1, wio2, wjo1, wjo2
real :: xoffpix, yoffpix
real :: glat, glon, arf
real :: rio, rjo
real :: xperdeg, yperdeg
real :: rlat1, rlat2, dlon, dlat
!!real, allocatable :: areai(:)

character(len=128) :: flnm
character(len=10)  :: sdate

character( 3) :: juldstring
character( 4) :: yearstring
character(20) :: header = 'emis2/v8_EM1L_'

logical :: nocall, exists
integer :: rcode, ncid

integer :: nlays_emis, nrows_emis, ncols_emis, nvars_emis
real    :: xorig, yorig
character(nf90_max_name) :: name, units

integer, external :: julday

! This subroutine is simpler than topm_database because it assumes that 
! each emis_database file covers the entire geographic area of the model.
! If this ever changes, this subroutine must be modified.

! Check type of call to emis_database_read

if (iaction == 0) then

   nemisfiles = 2
   allocate(s1900_emis(2))
   allocate(emis_files(2))

! Convert current model time from s1900 to years, months, dates, hours

   call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

   write(juldstring,'(I3)') julday( imonths, idates, iyears )
   write(yearstring,'(I4)') iyears

   indx = ihours / 4 + 1

   iemisy = iyears
   iemism = imonths
   iemisd = idates
   iemish = (indx - 1) * 4
   call date_abs_secs2(iemisy,iemism,iemisd,iemish*100,s1900_emis(1))

   iemisfile = 1
   flnm = trim(header) // yearstring // juldstring !! // '.nc'

   emis_files(1) = flnm
   emis_files(2) = flnm
   s1900_emis(2) = s1900_emis(1)

! HARD-WIRE EMISSIONS DATASET FOR NOW!!

!!   nemisfiles = 2  ! one per month
!!
!!   allocate(emis_files   (nemisfiles+2))
!!   allocate(ctotdate_emis(nemisfiles+2))

!  allocate(s1900_emis   (nemisfiles+2))

!!   allocate(s1900_emis(2))

!!   emis_files( 1) = 'emis2/global_emis_cb05_01.nc4'
!!   emis_files( 2) = 'emis2/global_emis_cb05_02.nc4'
!!   emis_files( 3) = 'emis2/global_emis_cb05_03.nc4'
!!   emis_files( 4) = 'emis2/global_emis_cb05_04.nc4'
!!   emis_files( 5) = 'emis2/global_emis_cb05_05.nc4'
!!   emis_files( 6) = 'emis2/global_emis_cb05_06.nc4'
!!   emis_files( 7) = 'emis2/global_emis_cb05_07.nc4'
!!   emis_files( 8) = 'emis2/global_emis_cb05_08.nc4'
!!   emis_files( 9) = 'emis2/global_emis_cb05_09.nc4'
!!   emis_files(10) = 'emis2/global_emis_cb05_10.nc4'
!!   emis_files(11) = 'emis2/global_emis_cb05_11.nc4'
!!   emis_files(12) = 'emis2/global_emis_cb05_12.nc4'

!!! Initialize sst cyclic flag to zero
!!
!!   iemiscyclic = 0
!!   nemisfiles  = 0
!!
!!   ! new schema to avoid call systems
!!   IF (isstflg == 0) THEN
!!      flnm = TRIM(sst_database)
!!      nocall = .TRUE.
!!      isstflg = 1
!!   ELSE
!!      flnm = TRIM(sst_database)//'??????????.h5'
!!      nocall = .FALSE.
!!   ENDIF
!!
!!   write(io6,*) 'Checking for sst database files'
!!   CALL OLAM_filelist(fnames_sst, maxsstfiles, flnm, nemisfiles, nocall)
!!
!!   if (nemisfiles < 1) then
!!      write(io6,*) 'SST database files '//flnm//' were not found.'
!!      write(io6,*) 'Stopping run.'
!!      stop 'stop: no sea file'
!!   endif
   
!!   ntimes = nemisfiles
!!
!!   do jtime = 1, ntimes
!!
!!      ! Assume the emissions file names aways end with MM.ncf, and
!!      ! use the file name to infer the file date
!!
!!      flnm  = emis_files(jtime)
!!      slen  = len_trim(flnm)
!!      sdate = flnm(slen-5:slen-4)
!!
!!      read(sdate,'(i2)') iemism
!!
!!      ! Hardwired for cyclic years and setting date to the middle of the month
!!
!!      iemisy =  0
!!      iemisd = 15
!!      iemish =  0
!!      
!!      ! If file year is read as zero, emis data is expected to be cyclic
!!      ! over 1 year. Increment iemiscyclic to indicate this and use current
!!      ! simulation year for emis database file times.
!!
!!      if (iemisy == 0) then
!!         iemiscyclic = iemiscyclic + 1
!!         iemisy = iyears
!!      endif
!!
!!      call date_make_big (iemisy,iemism,iemisd,iemish*100,ctotdate_emis(jtime))
!!      call date_abs_secs2(iemisy,iemism,iemisd,iemish*100,s1900_emis(jtime))
!!   enddo
!!
!!   ! Make sure files are sorted by date
!!   call dintsort28(nemisfiles,ctotdate_emis,emis_files,s1900_emis)
!!
!!! If emis cyclic flag > 0, check its value against ntimes and stop if they 
!!! are unequal.  If they are equal, reset emis cyclic flag to 1 and augment
!!! emis file arrays by 1 at each end.
!!
!!   if (iemiscyclic > 0) then
!!      if (iemiscyclic /= ntimes) then
!!         write(io6,'(/,a)') 'Some but not all emis database files do not have'
!!         write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
!!         stop 'stop_emis_inv'
!!      endif
!!      iemiscyclic = 1
!!      nemisfiles = ntimes + 2
!!
!!! Shift emis data file names and times by one array element
!!
!!      do jtime = ntimes,1,-1
!!         emis_files  (jtime+1) = emis_files  (jtime)
!!         ctotdate_emis(jtime+1) = ctotdate_emis(jtime)
!!         s1900_emis   (jtime+1) = s1900_emis   (jtime)
!!      enddo      
!!
!!! Add new emis member at beginning of time sequence
!!
!!      emis_files(1) = emis_files(ntimes+1)
!!
!!      call date_unmake_big(iemisy,iemism,iemisd,iemish,ctotdate_emis(ntimes+1))
!!      call date_make_big(iemisy-1,iemism,iemisd,iemish,ctotdate_emis(1))
!!      call date_abs_secs2(iemisy-1,iemism,iemisd,iemish,s1900_emis(1))
!!
!!! Add new emis member at end of time sequence
!!
!!      emis_files(ntimes+2) = emis_files(2)
!!
!!      call date_unmake_big(iemisy,iemism,iemisd,iemish,ctotdate_emis(2))
!!      call date_make_big(iemisy+1,iemism,iemisd,iemish,ctotdate_emis(ntimes+2))
!!      call date_abs_secs2(iemisy+1,iemism,iemisd,iemish,s1900_emis(ntimes+2))
!!
!!   endif
!!
!!! Loop over number of EMIS_DATABASE file times and search for the one that
!!! corresponds to current or most recent model time.
!!
!!   iemisfile = 0
!!   do nf = 1, nemisfiles
!!
!!      write(io6,*) 'nemisf0 ',nf,s1900_emis(nf),' ',s1900_sim
!!
!!      if (s1900_emis(nf) <= s1900_sim) then
!!         iemisfile = nf
!!      endif
!!   enddo
!!
!!   if (iemisfile < 1) then
!!      write(io6,*) ' '
!!      write(io6,*) 'Unable to find previous or current emis file for current'
!!      write(io6,*) 'model time.  Stopping model.'
!!      stop         'stop: no current emis file'
!!   endif

elseif (iaction == 1) then

   ! Processing next emis file (only called with iaction = 1 if iupdemis = 1)

   ! next emissions data will be 4 hours later

   call date_secs_ymdt(s1900_emis(1), oemisy, oemism, oemisd, oemish)

   call date_add_to8(oemisy, oemism, oemisd, oemish*100, 4.0_8, 'h', &
                     iemisy, iemism, iemisd, iemish)

   iemish = iemish / 100
   indx   = iemish / 4 + 1

   s1900_emis(1) = s1900_emis(2)
   emis_files(1) = emis_files(2)

   call date_abs_secs2(iemisy,iemism,iemisd,iemish*100,s1900_emis(2))

   write(juldstring,'(I3)') julday( iemism, iemisd, iemisy )
   write(yearstring,'(I4)') iemisy

   iemisfile = 2
   flnm = trim(header) // yearstring // juldstring !! // '.nc'
   emis_files(2) = flnm

!!  iemisfile = iemisfile + 1
!!   
!!   if (iemisfile > nemisfiles) then
!!      if (iemiscyclic == 0) then
!!         write(io6,*) ' '
!!         write(io6,*) 'No future emis file is available'
!!         write(io6,*) 'Stopping model '
!!         stop         'stop: no future emis file'
!!      else
!!         iemisfile = 3
!!         do jtime = 1, nemisfiles
!!            call date_unmake_big(iemisy,iemism,iemisd,iemish,ctotdate_emis(jtime))
!!            call date_make_big(iemisy+1,iemism,iemisd,iemish,ctotdate_emis(jtime))
!!            call date_abs_secs2(iemisy+1,iemism,iemisd,iemish,s1900_emis(jtime))
!!         enddo
!!      endif
!!   endif

    sea_emisp =  sea_emisf
   land_emisp = land_emisf

endif

! Open and read emis_database file

write(io6,*) 'emis_database_read ', iemisfile, indx, trim(emis_files(iemisfile))

inquire(file=emis_files(iemisfile), exist=exists)
if (.not. exists) then
   write(*,'(a)') ' emis_database_read: error opening input file ' // &
                  trim(emis_files(iemisfile))
   stop       'emis_database_read: fix location of emissions files'
endif

rcode = nf90_open(emis_files(iemisfile), NF90_NOWRITE, ncid )

if (rcode /= NF90_NOERR) then
   write(io6,*) nf90_strerror(rcode)
   write(io6,*) 'emis_database_read: error opening emissions file ' // &
                 trim(emis_files(iemisfile))
   stop
endif

rcode = nf90_get_att(ncid, NF90_GLOBAL, "NLAYS", nlays_emis)
rcode = nf90_get_att(ncid, NF90_GLOBAL, "NROWS", nrows_emis)
rcode = nf90_get_att(ncid, NF90_GLOBAL, "NCOLS", ncols_emis)
rcode = nf90_get_att(ncid, NF90_GLOBAL, "NVARS", nvars_emis)

rcode = nf90_get_att(ncid, NF90_GLOBAL, "XORIG", xorig)
rcode = nf90_get_att(ncid, NF90_GLOBAL, "YORIG", yorig)

xoffpix = 0.5
nio = ncols_emis + 2
nx  = ncols_emis
nis = 2
nie = ncols_emis + 1

yoffpix = 0.5
njo = nrows_emis + 2
ny  = nrows_emis
njs = 2
nje = nrows_emis + 1

xperdeg = real(nx) / 360.0
yperdeg = real(ny) / 180.0

dlon = 360.0 / real(nx)
dlat = 180.0 / real(ny)

!!allocate(areai(nrows_emis))
!!
!!do j = 1, ny
!!   rlat1 = pio180 * (yorig + dlat * real(j-1))
!!   rlat2 = pio180 * (yorig + dlat * real(j))
!!   areai(j) = 1.0 / (pio180 * erad**2 * abs(sin(rlat2)-sin(rlat1)) * dlon)
!!enddo

allocate(data(ncols_emis, nrows_emis, nlays_emis))
allocate(dato(nio, njo, nlays_emis))
!allocate(dato(nio, njo, nlays_emis, nvars_emis))

if (iaction == 0) then
   emlays = nlays_emis
   nvars3d_emis = nvars_emis

   allocate(vname3d_emis(nvars3d_emis))
   allocate(units3d_emis(nvars3d_emis))

   allocate( sea_emisf(emlays,mseaflux,nvars3d_emis))
   allocate(land_emisf(emlays,mlandflux,nvars3d_emis))

   allocate( sea_emisp(emlays,mseaflux,nvars3d_emis))
   allocate(land_emisp(emlays,mlandflux,nvars3d_emis))

   allocate( sea_emisc(emlays,mseaflux,nvars3d_emis))
   allocate(land_emisc(emlays,mlandflux,nvars3d_emis))

    sea_emisp(:,1,:) = 0.0
   land_emisp(:,1,:) = 0.0

    sea_emisf(:,1,:) = 0.0
   land_emisf(:,1,:) = 0.0

    sea_emisc(:,1,:) = 0.0
   land_emisc(:,1,:) = 0.0
endif

do n = 2, nvars3d_emis+1
   nv = n - 1

   rcode = nf90_inquire_variable(ncid, n, name)

   if (iaction == 0) then
      rcode = nf90_get_att(ncid, n, "units", units)
      vname3d_emis(nv) = name
      if (units == 'moles/m2-s') then
         units = 'moles/s'
      endif
      units3d_emis(nv) = units
   else
      if (name /= vname3d_emis(nv)) then
         write(io6,*) "Error in emis_database read:"
         write(io6,*) "Emissions file format has changed."
         stop
      endif
   endif

   write(io6,*) "Reading " // trim(name) // " from emissions file."

   rcode = nf90_get_var( ncid, n, data, start=(/1,1,1,indx/), &
                         count=(/ncols_emis,nrows_emis,nlays_emis,1/) )

   ! new emissions are per m^s

   do k = 1, nlays_emis
      do j = 1, ny
         do i = 1, nx
            ii = i + 1
            jj = j + 1
!           dato(ii,jj,k,nv) = data(i,j,k) * areai(j)
            dato(ii,jj,k) = data(i,j,k)
         enddo
      enddo
   enddo

! If data is offset 1/2 delta longitude from the lower left point,
! add extra cyclic columns around the border for easier interpolation

   if (mod(ncols_emis,2) == 0) then
      dato(  1,njs:nje,:) = dato(nie,njs:nje,:)
      dato(nio,njs:nje,:) = dato(nis,njs:nje,:)
   endif

! If data is offset 1/2 delta latitude from the lower left point,
! add an extra row at top and bottom for easier interpolation

   if (mod(nrows_emis,2) == 0) then
      dato(:,  1,:) = dato(:,  2,:)
      dato(:,njo,:) = dato(:,nje,:)
   endif

! Fill emis arrays

   do isf = 2, mseaflux

      glat = seaflux(isf)%glatf
      glon = seaflux(isf)%glonf
      arf  = seaflux(isf)%area

!     glon = max(-179.999,min(179.999,glon))

      ! Convert OLAM's -180,180 longitude range to 0,360 

      if (glon < 0.0) glon = 360.0 + glon 
      glon = max(0.001,min(359.999,glon))

!     rio = 1. + (glon + 180.) * xperdeg + xoffpix
      rio = 1. + (glon      ) * xperdeg + xoffpix
      rjo = 1. + (glat + 90.) * yperdeg + yoffpix

      io1 = int(rio)
      jo1 = int(rjo)
         
      wio2 = rio - real(io1)
      wjo2 = rjo - real(jo1)
           
      wio1 = 1. - wio2
      wjo1 = 1. - wjo2

      io2 = min(nio, io1 + 1)
      jo2 = min(njo, jo1 + 1)

      do k = 1, emlays
         sea_emisf(k,isf,nv) = &
              ( wio1 * (wjo1 * dato(io1,jo1,k) + wjo2 * dato(io1,jo2,k)) &
              + wio2 * (wjo1 * dato(io2,jo1,k) + wjo2 * dato(io2,jo2,k)) ) * arf
      enddo
   enddo

   do ilf = 2, mlandflux
   
      glat = landflux(ilf)%glatf
      glon = landflux(ilf)%glonf
      arf  = landflux(ilf)%area

      ! Convert OLAM's -180,180 longitude range to 0,360 

      if (glon < 0.0) glon = 360.0 + glon 
      glon = max(0.001,min(359.999,glon))

!     rio = 1. + (glon + 180.) * xperdeg + xoffpix
      rio = 1. + (glon       ) * xperdeg + xoffpix
      rjo = 1. + (glat +  90.) * yperdeg + yoffpix

      io1 = int(rio)
      jo1 = int(rjo)
         
      wio2 = rio - real(io1)
      wjo2 = rjo - real(jo1)
           
      wio1 = 1. - wio2
      wjo1 = 1. - wjo2

      io2 = min(nio, io1 + 1)
      jo2 = min(njo, jo1 + 1)

      do k = 1, emlays
         land_emisf(k,ilf,nv) = &
              ( wio1 * (wjo1 * dato(io1,jo1,k) + wjo2 * dato(io1,jo2,k)) &
              + wio2 * (wjo1 * dato(io2,jo1,k) + wjo2 * dato(io2,jo2,k)) ) * arf
      enddo

   enddo

enddo

rcode = nf90_close(ncid)

deallocate(data)
deallocate(dato)

if (iaction == 0) then
    sea_emisp =  sea_emisf
   land_emisp = land_emisf
endif

return
end subroutine emis_database_read
