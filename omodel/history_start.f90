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
subroutine history_start(action)

! This routine initializes the model from the history file

use misc_coms,  only: io6, hfilin, time8, time_istp8, runtype, iparallel
use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank, mgroupsize
use max_dims,   only: pathlen

implicit none

character(*), intent(in) :: action

logical                  :: exans, exanz, exanp
integer                  :: nfl, iter
character(pathlen+10)    :: hfilinp
character(10)            :: number

! Check if history files exist

exanp = .false.
nfl = len_trim(hfilin)
hfilinp = hfilin(1:nfl-3)//'_r0.h5'

inquire(file=hfilin,  exist=exans)   ! serial restart file
inquire(file=hfilinp, exist=exanz)   ! rank 0 restart file

if (runtype == 'PARCOMBINE' .and. .not. exanz) then
   write(io6,*)
   write(io6,*) "No parallel history file exists."
   write(io6,*) "Looking for ", trim(hfilinp)
   write(io6,*) "Stopping PARCOMBINE run."
   stop
endif

! PARCOMBINE RUN SHOULD NOT READ SERIAL HISTORY FILES

if (runtype == 'PARCOMBINE') exans = .false.

if (iparallel == 1 .and. exanz) then
   do iter = 0, mgroupsize-1
      write(number,'(i6)') iter
      hfilinp = hfilin(1:nfl-3)//'_r'//trim(adjustl(number))//'.h5'
      inquire(file=hfilinp, exist=exanp)
      if (.not. exanp) exit
   enddo
endif
      
if ((exans .and. iparallel == 0) .or. (exanp .and. iparallel == 1)) then

! Serial run reading its sequential history file, or a
! parallel run reading its history file

   if (iparallel == 1) then
      write(number,'(i6)') myrank
      hfilinp = hfilin(1:nfl-3)//'_r'//trim(adjustl(number))//'.h5'
   else
      hfilinp = hfilin
   endif
    
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(io6,*) 'Opening history file '//trim(hfilinp)//' for '//trim(action)
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   call shdf5_open(trim(hfilinp),'R')

   if (trim(action) == 'COMMIO') then

      ! Read the common variables
      call commio('READ')
   else

      ! Read the model fields
      call shdf5_irec(1, (/1/), 'time8', dvars=time8)
      call hist_read()
      write(io6,*) 'back from hist_read'

   endif

   time_istp8 = time8  ! time_istp8 is used for plotting time
   call shdf5_close()

elseif (exanz .and. iparallel == 0) then

! Serial run reading parallel history files

   if (trim(action) == 'COMMIO') then

      ! Commons only need to be read from one parallel file
      hfilinp = hfilin(1:nfl-3)//'_r0.h5'
      write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(io6,*) 'Opening history file '//trim(hfilinp)//' for '//trim(action)
      write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      call shdf5_open(trim(hfilinp),'R')
      call commio('READ')
      call shdf5_close()

   else

      iter = -1
      do
         iter = iter + 1
         write(number,'(i6)') iter
         hfilinp = hfilin(1:nfl-3)//'_r'//trim(adjustl(number))//'.h5'
      
         ! check if history file exists from parallel run, and exit loop
         ! when file does not exist
         inquire(file=hfilinp, exist=exanp)
         if (.not. exanp) exit

         ! Parallel history file exists.  Open, read, and close file.
         write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(io6,*) 'Opening history file '//trim(hfilinp)//' for '//trim(action)
         write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         call shdf5_open(trim(hfilinp),'R')
         call shdf5_irec(1, (/1/), 'time8', dvars=time8)
         call hist_read_p(iter)
         time_istp8 = time8  ! time_istp8 is used for plotting time
         call shdf5_close()
      enddo

   endif

elseif (exans .and. iparallel == 1) then

! Parallel run reading serial history file

   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(io6,*) 'Opening history file '//trim(hfilin)//' for '//trim(action)
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   call shdf5_open(trim(hfilin),'R')

   if (trim(action) == 'COMMIO') then
      call commio('READ')
   else
      call shdf5_irec(1, (/1/), 'time8', dvars=time8)
      call hist_read_s()
   endif
   time_istp8 = time8  ! time_istp8 is used for plotting time
   call shdf5_close()

else

   ! History files do not exist, stop model.
   
   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(io6,*) '!!!   Trying to open history file file:'
   if (iparallel == 1) then
   write(io6,*) '!!!   '//trim(hfilinp)
   else
   write(io6,*) '!!!   '//trim(hfilin)
   endif
   write(io6,*) '!!!   but it does not exist. The run is ended.'
   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'in history_start'
   
endif

return
end subroutine history_start

!=========================================================================

subroutine hist_read()

use misc_coms,  only: io6, hfilin, runtype
use var_tables, only: num_var, vtab_r, get_vtab_dims
use hdf5_utils, only: shdf5_open, shdf5_info, shdf5_irec, shdf5_close

implicit none

integer           :: nv, nvcnt, ndims, ndims_m, idims(3), idims_m(3)
character(len=32) :: varn

nvcnt = 0
ndims = 0
idims = 0

! Loop through all variables in the model vtables

do nv = 1,num_var

! Skip to next variable if we don't want the current one

   if (.not. vtab_r(nv)%ihist) cycle
   if (vtab_r(nv)%nread .and. runtype == 'HISTORY') cycle
   
   varn = trim(vtab_r(nv)%name)
   call get_vtab_dims(nv, ndims_m, idims_m)

! We want it...read it if it's in the history file

   call shdf5_info(varn, ndims, idims)

! Skip to next variable if the current one is not in the history file

   if (ndims == 0 .and. idims(1) == 0) then
      write(io6,*) 'Variable '//trim(varn)//' is not in the history file.'
      cycle
   endif

! Report an error if the variable dimensions are different

   if (ndims /= ndims_m .or. any(idims(1:ndims) /= idims_m(1:ndims))) then
      write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(io6,*) '!!!   Trying to read variable '//trim(varn)
      write(io6,*) '!!!   from the history file '//trim(hfilin)
      write(io6,*) '!!!   but the array size is different. '
      write(io6,*) '!!!   ndims, ndims_m = ',ndims,ndims_m
      write(io6,*) '!!!   idims  (1:ndims) = ',idims  (1:ndims)
      write(io6,*) '!!!   idims_m(1:ndims) = ',idims_m(1:ndims)
      write(io6,*) '!!!   The run is ended'
      write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop    'stop in hist_read'
   endif

   if     (associated(vtab_r(nv)%ivar1_p)) then
      call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar1_p)
   elseif (associated(vtab_r(nv)%ivar2_p)) then
      call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar2_p)
   elseif (associated(vtab_r(nv)%ivar3_p)) then
      call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar3_p)

   elseif (associated(vtab_r(nv)%rvar1_p)) then
      call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar1_p)
   elseif (associated(vtab_r(nv)%rvar2_p)) then
      call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar2_p)
   elseif (associated(vtab_r(nv)%rvar3_p)) then
      call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar3_p)

   elseif (associated(vtab_r(nv)%dvar1_p)) then
      call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar1_p)
   elseif (associated(vtab_r(nv)%dvar2_p)) then
      call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar2_p)
   elseif (associated(vtab_r(nv)%dvar3_p)) then
      call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar3_p)
   endif
   
   nvcnt = nvcnt + 1
   write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
      'Read: ', nvcnt, trim(varn), idims(1:ndims)
enddo

return
end subroutine hist_read

!=========================================================================

subroutine hist_read_p(myrank)

! Serial run reading parallel history files.
! Also used for parcombine and plotonly runs

use misc_coms,   only: io6, runtype, meshtype
use var_tables,  only: num_var, vtab_r, get_vtab_dims
use hdf5_utils,  only: shdf5_info, shdf5_irec
use mem_leaf,    only: land
use mem_sea,     only: sea
use consts_coms, only: r8
implicit none

integer, intent(in) :: myrank
integer :: nv, nvcnt, ndims, ndims_m, idims(3), idims_m(3)
integer :: i,k,ig
character(len=32) :: varn

! GLOBAL INDICES AND RANKS:
integer :: idgm(3),  idgu(3),  idgv(3), idgw(3)
integer, target, allocatable :: imglobe(:), irankm(:)
integer, target, allocatable :: iuglobe(:), iranku(:)
integer, target, allocatable :: ivglobe(:), irankv(:)
integer, target, allocatable :: iwglobe(:), irankw(:)

! LAND INDICES AND RANKS:
integer :: idlm(3),  idlu(3),  idlw(3)
integer, target, allocatable :: ilmglobe(:), ilrankm(:)
integer, target, allocatable :: iluglobe(:), ilranku(:)
integer, target, allocatable :: ilwglobe(:), ilrankw(:)

! SEA INDICES AND RANKS:
integer :: idsm(3),  idsu(3),  idsw(3)
integer, target, allocatable :: ismglobe(:), isrankm(:)
integer, target, allocatable :: isuglobe(:), isranku(:)
integer, target, allocatable :: iswglobe(:), isrankw(:)

! LANDFLUX INDICES AND RANKS:
integer :: idlf(3)
integer, target, allocatable :: ilfglobe(:), ilfrank(:)

! SEAFLUX INDICES AND RANKS:
integer :: idsf(3)
integer, target, allocatable :: isfglobe(:), isfrank(:)

! NUDGING INDICES
integer :: idn(3)
integer, target, allocatable :: inglobe(:), inrank(:)

character(len=2)   :: stagpt
integer, pointer :: iglobe(:), irank(:)

integer, allocatable :: iscr1(:)
integer, allocatable :: iscr2(:,:)
integer, allocatable :: iscr3(:,:,:)

real, allocatable :: scr1(:)
real, allocatable :: scr2(:,:)
real, allocatable :: scr3(:,:,:)

real(kind=r8), allocatable :: dscr1(:)
real(kind=r8), allocatable :: dscr2(:,:)
real(kind=r8), allocatable :: dscr3(:,:,:)

! Get dimensions, allocate, and read indices and processor ranks:

call shdf5_info('IMGLOBE',ndims,idgm)
if (idgm(1) > 0) then
   allocate (imglobe(idgm(1)))
   call shdf5_irec(ndims, idgm, 'IMGLOBE', ivara=imglobe)
   ! we don't have a IRANKM variable, so create a temporary one
   allocate (irankm(idgm(1)))
   irankm(:) = myrank
endif

if (meshtype == 1) then

   call shdf5_info('IUGLOBE',ndims,idgu)
   if (idgu(1) > 0) then
      allocate (iuglobe(idgu(1)))
      call shdf5_irec(ndims, idgu, 'IUGLOBE', ivara=iuglobe)
   endif

   idgv = 0
else

   call shdf5_info('IVGLOBE',ndims,idgv)
   if (idgv(1) > 0) then
      allocate (ivglobe(idgv(1)))
      call shdf5_irec(ndims, idgv, 'IVGLOBE', ivara=ivglobe)
   endif

   idgu = 0
endif

call shdf5_info('IWGLOBE',ndims,idgw)
if (idgw(1) > 0) then
   allocate (iwglobe(idgw(1)))
   call shdf5_irec(ndims, idgw, 'IWGLOBE', ivara=iwglobe)
endif

call shdf5_info('IMGLOBE_L',ndims,idlm)
if (idlm(1) > 0) then
   allocate (ilmglobe(idlm(1)))
   call shdf5_irec(ndims, idlm, 'IMGLOBE_L', ivara=ilmglobe)
   ! we don't have a IRANKM_L variable, so create a temporary one
   allocate (ilrankm(idlm(1)))
   ilrankm(:) = myrank
endif

call shdf5_info('IUGLOBE_L',ndims,idlu)
if (idlu(1) > 0) then
   allocate (iluglobe(idlu(1)))
   call shdf5_irec(ndims, idlu, 'IUGLOBE_L', ivara=iluglobe)
endif

call shdf5_info('IWGLOBE_L',ndims,idlw)
if (idlw(1) > 0) then
   allocate (ilwglobe(idlw(1)))
   call shdf5_irec(ndims, idlw, 'IWGLOBE_L', ivara=ilwglobe)
endif

call shdf5_info('IMGLOBE_S',ndims,idsm)
if (idsm(1) > 0) then
   allocate (ismglobe(idsm(1)))
   call shdf5_irec(ndims, idsm, 'IMGLOBE_S', ivara=ismglobe)
   ! we don't have a IRANKM_S variable, so create a temporary one
   allocate (isrankm(idsm(1)))
   isrankm(:) = myrank
endif

call shdf5_info('IUGLOBE_S',ndims,idsu)
if (idsu(1) > 0) then
   allocate (isuglobe(idsu(1)))
   call shdf5_irec(ndims, idsu, 'IUGLOBE_S', ivara=isuglobe)
endif

call shdf5_info('IWGLOBE_S',ndims,idsw)
if (idsw(1) > 0) then
   allocate (iswglobe(idsw(1)))
   call shdf5_irec(ndims, idsw, 'IWGLOBE_S', ivara=iswglobe)
endif

call shdf5_info('LANDFLUX%IFGLOBE',ndims,idlf)
if (idlf(1) > 0) then
   allocate (ilfglobe(idlf(1)))
   call shdf5_irec(ndims, idlf, 'LANDFLUX%IFGLOBE', ivara=ilfglobe)
endif

call shdf5_info('SEAFLUX%IFGLOBE',ndims,idsf)
if (idsf(1) > 0) then
   allocate (isfglobe(idsf(1)))
   call shdf5_irec(ndims, idsf, 'SEAFLUX%IFGLOBE', ivara=isfglobe)
endif

if (meshtype == 1) then

   call shdf5_info('IRANKU',ndims,idims)
   if (idims(1) > 0) then
      allocate (iranku(idims(1)))
      call shdf5_irec(ndims, idims, 'IRANKU', ivara=iranku)
   endif

else

   call shdf5_info('IRANKV',ndims,idims)
   if (idims(1) > 0) then
      allocate (irankv(idims(1)))
      call shdf5_irec(ndims, idims, 'IRANKV', ivara=irankv)
   endif

endif

call shdf5_info('IRANKW',ndims,idims)
if (idims(1) > 0) then
   allocate (irankw(idims(1)))
   call shdf5_irec(ndims, idims, 'IRANKW', ivara=irankw)
endif

call shdf5_info('IRANKU_L',ndims,idims)
if (idims(1) > 0) then
   allocate (ilranku(idims(1)))
   call shdf5_irec(ndims, idims, 'IRANKU_L', ivara=ilranku)
endif

call shdf5_info('IRANKW_L',ndims,idims)
if (idims(1) > 0) then
   allocate (ilrankw(idims(1)))
   call shdf5_irec(ndims, idims, 'IRANKW_L', ivara=ilrankw)
endif

call shdf5_info('IRANKU_S',ndims,idims)
if (idims(1) > 0) then
   allocate (isranku(idims(1)))
   call shdf5_irec(ndims, idims, 'IRANKU_S', ivara=isranku)
endif

call shdf5_info('IRANKW_S',ndims,idims)
if (idims(1) > 0) then
   allocate (isrankw(idims(1)))
   call shdf5_irec(ndims, idims, 'IRANKW_S', ivara=isrankw)
endif

call shdf5_info('LANDFLUX%IWRANK',ndims,idims)
if (idims(1) > 0) then
   allocate (ilfrank(idims(1)))
   call shdf5_irec(ndims, idims, 'LANDFLUX%IWRANK', ivara=ilfrank)
endif

call shdf5_info('SEAFLUX%IWRANK',ndims,idims)
if (idims(1) > 0) then
   allocate (isfrank(idims(1)))
   call shdf5_irec(ndims, idims, 'SEAFLUX%IWRANK', ivara=isfrank)
endif

! Loop through all variables in the model vtables

nvcnt = 0
do nv = 1,num_var

! Skip to next variable if we don't want the current one

   if (.not. vtab_r(nv)%ihist) cycle
   if (vtab_r(nv)%nread .and. runtype == 'HISTORY') cycle

   varn    = trim(vtab_r(nv)%name)
   stagpt  = vtab_r(nv)%stagpt
   call get_vtab_dims(nv, ndims_m, idims_m)

! We want it...read it if it's in the history file

   call shdf5_info(varn,ndims,idims)

! Skip to next variable if the current one is not in the history file

   if (ndims == 0 .and. idims(1) == 0) then
      write(io6,*) 'Variable '//trim(varn)//' is not in the history file.'
      cycle
   endif

   if (    stagpt == 'AW' .and. idims(ndims) == idgw(1)) then
      iglobe => iwglobe
      irank  => irankw
   elseif (meshtype == 1 .and. stagpt == 'AU' .and. idims(ndims) == idgu(1)) then
      iglobe => iuglobe
      irank  => iranku
   elseif (meshtype == 2 .and. stagpt == 'AU' .and. idims(ndims) == idgv(1)) then
      iglobe => ivglobe
      irank  => irankv
   elseif (stagpt == 'AM' .and. idims(ndims) == idgm(1)) then
      iglobe => imglobe
      irank  => irankm
   elseif (stagpt == 'LW' .and. idims(ndims) == idlw(1)) then
      iglobe => ilwglobe
      irank  => ilrankw
   elseif (stagpt == 'LU' .and. idims(ndims) == idlu(1)) then
      iglobe => iluglobe
      irank  => ilranku
   elseif (stagpt == 'LM' .and. idims(ndims) == idlm(1)) then
      iglobe => ilmglobe
      irank  => ilrankm
   elseif (stagpt == 'SW' .and. idims(ndims) == idsw(1)) then
      iglobe => iswglobe
      irank  => isrankw
   elseif (stagpt == 'SU' .and. idims(ndims) == idsu(1)) then
      iglobe => isuglobe
      irank  => isranku
   elseif (stagpt == 'SM' .and. idims(ndims) == idsm(1)) then
      iglobe => ismglobe
      irank  => isrankm
   elseif (stagpt == 'LF' .and. idims(ndims) == idlf(1)) then
      iglobe => ilfglobe
      irank  => ilfrank
   elseif (stagpt == 'SF' .and. idims(ndims) == idsf(1)) then
      iglobe => isfglobe
      irank  => isfrank
!  elseif (stagpt == 'AN') then
!
!     TODO: NUDGING DOES NOT YET WORK IN PARALLEL.
!     ASSUME FOR NOW THAT EACH PROCESSOR WILL HAVE THE ENTIRE
!     COARSE GLOBAL OBSERVATIONS???
!
!     allocate(inglobe(idims(ndims)))
!     allocate(inrank(idims(ndims)))
!     do i=1,idims(ndims)
!        inglobe(i) = i
!        inrank(i) = 0
!     enddo
!     iglobe => inglobe
!     irank  => inrank
!  elseif (stagpt == 'CN') then
!
!     TODO: READING CONSTANT VALUES
!     DOES NOT WORK YET IN PARALLEL
!
   else
      write(io6,*)
      write(io6,*) 'Error for variable '//trim(varn)//':'
      stop ' Invalid array size in history_start_p.'
   endif

   if (associated(vtab_r(nv)%ivar1_p)) then
      allocate (iscr1(idims(1)))
      call shdf5_irec(ndims, idims, trim(varn), ivara=iscr1)

      do i = 1,idims(1)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         vtab_r(nv)%ivar1_p(ig) = iscr1(i)
      enddo
      deallocate (iscr1)

   elseif (associated(vtab_r(nv)%ivar2_p)) then

      allocate (iscr2(idims(1),idims(2)))
      call shdf5_irec(ndims, idims, trim(varn), ivara=iscr2)

      do i = 1,idims(2)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         do k = 1,idims(1)
            vtab_r(nv)%ivar2_p(k,ig) = iscr2(k,i)
         enddo
      enddo
      deallocate (iscr2)

   elseif (associated(vtab_r(nv)%rvar1_p)) then

      allocate (scr1(idims(1)))
      call shdf5_irec(ndims, idims, trim(varn), rvara=scr1)

      do i = 1,idims(1)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         vtab_r(nv)%rvar1_p(ig) = scr1(i)
      enddo
      deallocate (scr1)

   elseif (associated(vtab_r(nv)%rvar2_p)) then

      allocate (scr2(idims(1),idims(2)))
      call shdf5_irec(ndims, idims, trim(varn), rvara=scr2)

      do i = 1,idims(2)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         do k = 1,idims(1)
            vtab_r(nv)%rvar2_p(k,ig) = scr2(k,i)
         enddo
      enddo
      deallocate (scr2)

   elseif (associated(vtab_r(nv)%dvar1_p)) then

      allocate (dscr1(idims(1)))
      call shdf5_irec(ndims, idims, trim(varn), dvara=dscr1)

      do i = 1,idims(1)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         vtab_r(nv)%dvar1_p(ig) = dscr1(i)
      enddo
      deallocate (dscr1)

   elseif (associated(vtab_r(nv)%dvar2_p)) then

      allocate (dscr2(idims(1),idims(2)))
      call shdf5_irec(ndims, idims, trim(varn), dvara=dscr2)

      do i = 1,idims(2)
         ig = iglobe(i)
         if (irank(i) /= myrank) cycle
         do k = 1,idims(1)
            vtab_r(nv)%dvar2_p(k,ig) = dscr2(k,i)
         enddo
      enddo
      deallocate (dscr2)

   endif
   
   nvcnt = nvcnt + 1
   write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
      'Read: ', nvcnt, trim(varn), idims(1:ndims)

enddo

if (allocated(imglobe))  deallocate (imglobe)
if (allocated(iuglobe))  deallocate (iuglobe)
if (allocated(iwglobe))  deallocate (iwglobe)
if (allocated(ilmglobe)) deallocate (ilmglobe)
if (allocated(iluglobe)) deallocate (iluglobe)
if (allocated(ilwglobe)) deallocate (ilwglobe)
if (allocated(ismglobe)) deallocate (ismglobe)
if (allocated(isuglobe)) deallocate (isuglobe)
if (allocated(iswglobe)) deallocate (iswglobe)
if (allocated(ilfglobe)) deallocate (ilfglobe)
if (allocated(isfglobe)) deallocate (isfglobe)
if (allocated(inglobe))  deallocate (inglobe)

if (allocated(irankm))  deallocate (irankm)
if (allocated(iranku))  deallocate (iranku)
if (allocated(irankw))  deallocate (irankw)
if (allocated(ilrankm)) deallocate (ilrankm)
if (allocated(ilranku)) deallocate (ilranku)
if (allocated(ilrankw)) deallocate (ilrankw)
if (allocated(isrankm)) deallocate (isrankm)
if (allocated(isranku)) deallocate (isranku)
if (allocated(isrankw)) deallocate (isrankw)
if (allocated(ilfrank)) deallocate (ilfrank)
if (allocated(isfrank)) deallocate (isfrank)
if (allocated(inrank))  deallocate (inrank)

return
end subroutine hist_read_p

!=========================================================================

subroutine hist_read_s()

! Parallel run reading serial history file

use mem_grid,    only: nwa, nua, nva, nma, mwa, mua, mva, mma, mza
use mem_ijtabs,  only: itabg_w, itab_w, itab_u, itab_v, itab_m
use misc_coms,   only: io6, meshtype, runtype
use var_tables,  only: num_var, vtab_r, get_vtab_dims
use hdf5_utils,  only: shdf5_info, shdf5_irec
use mem_sflux,   only: nlandflux, nseaflux, mlandflux, mseaflux, &
                       landflux,  seaflux
use mem_leaf,    only: itab_wl
use leaf_coms,   only: nwl, mwl, nzg
use mem_sea,     only: itab_ws
use sea_coms,    only: nws, mws
use mem_para,    only: myrank
use consts_coms, only: r8
implicit none

integer          :: nv, nvcnt, ndims, idims(3)
integer          :: i, k, il, iw, isf, ilf
character(32)    :: varn
character (2)    :: stagpt
integer, pointer :: ilocal(:)

nvcnt = 0
ndims = 0
idims = 0

! Loop through all variables in the model vtables

do nv = 1,num_var

! Skip to next variable if we don't want the current one

   if (.not. vtab_r(nv)%ihist) cycle
   if (vtab_r(nv)%nread .and. runtype == 'HISTORY') cycle

   varn    = trim(vtab_r(nv)%name)
   stagpt  = vtab_r(nv)%stagpt

! We want it...read it if it's in the history file

   call shdf5_info(varn, ndims, idims)

! Skip to next variable if the current one is not in the history file

   if (ndims == 0 .and. idims(1) == 0) then
      write(io6,*) 'Variable '//trim(varn)//' is not in the history file.'
      cycle
   endif

   if     (stagpt == 'AW' .and. idims(ndims) == nwa) then
      ilocal => itab_w(:)%iwglobe
      idims(ndims) = mwa
   elseif (stagpt == 'AU' .and. idims(ndims) == nua .and. meshtype == 1) then
      ilocal => itab_u(:)%iuglobe
      idims(ndims) = mua
   elseif (stagpt == 'AU' .and. idims(ndims) == nva .and. meshtype == 2) then
      ilocal => itab_v(:)%ivglobe
      idims(ndims) = mva
   elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
      ilocal => itab_m(:)%imglobe
      idims(ndims) = mma
   elseif (stagpt == 'LW' .and. idims(ndims) == nwl) then
      ilocal => itab_wl(:)%iwglobe
      idims(ndims) = mwl
   elseif (stagpt == 'SW' .and. idims(ndims) == nws) then
      ilocal => itab_ws(:)%iwglobe
      idims(ndims) = mws
   elseif (stagpt == 'LF' .and. idims(ndims) == nlandflux) then
      ilocal => landflux(:)%ifglobe
      idims(ndims) = mlandflux
   elseif (stagpt == 'SF' .and. idims(ndims) == nseaflux) then
      ilocal => seaflux(:)%ifglobe
      idims(ndims) = mseaflux
   else

      ! TODO: Nudging values, const values
      ! TODO: Land U and M; Sea U and M?

      stop "invalid array size in history_start_s"
   endif

   if     (associated(vtab_r(nv)%ivar1_p)) then
      call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar1_p, points=ilocal)
   elseif (associated(vtab_r(nv)%ivar2_p)) then
      call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar2_p, points=ilocal)
   elseif (associated(vtab_r(nv)%rvar1_p)) then
      call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar1_p, points=ilocal)
   elseif (associated(vtab_r(nv)%rvar2_p)) then
      call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar2_p, points=ilocal)
   elseif (associated(vtab_r(nv)%dvar1_p)) then
      call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar1_p, points=ilocal)
  elseif (associated(vtab_r(nv)%dvar2_p)) then
      call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar2_p, points=ilocal)
   endif

   ! THE FOLLOWING SECTION IS TO ZERO OUT BOUNDARY CELLS FOR FIELDS
   ! THAT ARE TYPICALLY ONLY COMPUTED ON "PRIMARY" POINTS. WHILE NOT
   ! STRICTLY NECESSARY, IT WILL ENSURE THAT THE OUTPUT IS THE SAME
   ! FOR A HISTORY RUN READING EITHER A PARALLEL OR SEQUENTIAL OUTPUT

   if (any( varn == (/'FTHRD   ', 'FTHRD_LW'/) )) then
      do i = 2,mwa
         if (itab_w(i)%irank /= myrank) then
            vtab_r(nv)%rvar2_p(1:mza,i) = 0.0
         endif
      enddo
   endif
   
   if (any( varn == (/'RSHORT ', 'RLONG  ', 'RLONGUP', 'ALBEDT ',  &
                      'COSZ   ', 'SFLUX_T', 'SFLUX_R', 'USTAR  '/) )) then
      do i = 2,mwa
         if (itab_w(i)%irank /= myrank) then
            vtab_r(nv)%rvar1_p(i) = 0.0
         endif
      enddo
   endif

   if (any( varn == (/'LANDFLUX%SFLUXR', 'LANDFLUX%SFLUXT'/) )) then
      do i = 2,mlandflux
         iw = landflux(i)%iw
         if (itabg_w(iw)%irank /= myrank) then
            vtab_r(nv)%rvar1_p(i) = 0.0
         endif
      enddo
   endif

   if (any( varn == (/'SEAFLUX%SFLUXR', 'SEAFLUX%SFLUXT'/) )) then
      do i = 2,mseaflux
         iw = seaflux(i)%iw
         if (itabg_w(iw)%irank /= myrank) then
            vtab_r(nv)%rvar1_p(i) = 0.0
         endif
      enddo
   endif

   if (any( varn == (/'SEA%SEATC         ', 'SEA%SEAICEC       ',  &
                      'SEA%CAN_DEPTH     ', 'SEA%RLONG         ',  &
                      'SEA%SURFACE_SSH   ', 'SEA%RSHORT_DIFFUSE',  &
                      'SEA%RSHORT        '/) )) then
      do i = 2,mws
         if (itab_ws(i)%irank /= myrank) then
            vtab_r(nv)%rvar1_p(i) = 0.0
         endif
      enddo
   endif

   if (any( varn == (/'LAND%CAN_DEPTH   ', 'LAND%GROUND_SHV  ',  &
                      'LAND%RLONG       ', 'LAND%RLONG_G     ',  &
                      'LAND%STOM_RESIST ', 'LAND%SURFACE_SSH ',  &
                      'LAND%VEG_ALBEDO  ', 'LAND%VEG_FRACAREA',  &
                      'LAND%VEG_HEIGHT  ', 'LAND%VEG_LAI     ',  &
                      'LAND%VEG_TAI     ', 'LAND%VEG_TEMP    ',  &
                      'LAND%HCAPVEG     ', 'LAND%RLONG_V     ',  &
                      'LAND%VEG_ROUGH   ', 'LAND%VF          '/) )) then
      do i = 2,mwl
         if (itab_wl(i)%irank /= myrank) then
            vtab_r(nv)%rvar1_p(i) = 0.0
         endif
      enddo
   endif

   if (any( varn == (/'LAND%SOIL_ENERGY', 'LAND%SOIL_WATER '/) )) then
      do i = 2,mwl
         if (itab_wl(i)%irank /= myrank) then
            vtab_r(nv)%rvar2_p(1:nzg,i) = 0.0
         endif
      enddo
   endif

   nvcnt = nvcnt + 1
   write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
      'Read: ', nvcnt, trim(varn), idims(1:ndims)

enddo

return
end subroutine hist_read_s
