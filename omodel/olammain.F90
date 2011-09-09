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
program main

use misc_coms,    only: io6, iparallel
use mem_para,     only: myrank, mgroupsize
use olam_mpi_atm, only: olam_mpi_init, olam_mpi_finalize
use misc_coms,    only: tmpdir
use max_dims,     only: pathlen

#ifdef OLAM_HDF5_FORTRAN
use hdf5
#endif

implicit none

character(pathlen) :: name_name = 'OLAMIN'
character(pathlen) :: cargv
character(30)      :: io6file
integer :: numarg
integer :: i,n
integer :: bad = 0
integer :: hdferr

! Determine if this run is parallel, and determine myrank and mgroupsize

call olam_mpi_init()

iparallel = 0
if (mgroupsize > 1) iparallel = 1

! If run is sequential, default choice is to set io6 to standard output unit 6.

io6 = 6

if (iparallel == 1 .and. myrank > 0) then

! If run is parallel, default choice is to attach output unit io6 to separate files

   io6 = 20
   
   write (io6file,'(i10)') myrank
   io6file = 'o.io6_r'//trim(adjustl(io6file))

! Output file is replaced if it exists
   open(unit=io6, file=io6file, status='replace', form='formatted')

#ifdef __PGI
   ! Prevent Portland Group compiler from buffering io6
   call setvbuf3f(io6,1,1024)
#endif

endif

write(io6,'(/,a,i6)') ' myrank     = ',myrank
write(io6,'(  a,i6)') ' mgroupsize = ',mgroupsize
write(io6,'(  a,i6)') ' iparallel  = ',iparallel

! initialize HDF5 library
#ifdef OLAM_HDF5_FORTRAN
call h5open_f(hdferr)
#endif

numarg = command_argument_count()
write(io6,*) 'numarg:', numarg

! Parse the command line arguments

i = 1
do while (i <= numarg)
   call get_command_argument(i,cargv)
   ! write(io6,*) 'args: ',i,cargv

   if (cargv(1:1) == '-') then
      if (cargv(2:2) == 'f') then
         call get_command_argument(i+1,name_name)
         if (len_trim(name_name) < 1) bad = bad + 1
         i = i + 2
      else
         ! write(io6,*) 'OLAM unknown option: ', cargv
         i = i + 1
      endif
   else
      ! write(io6,*) 'OLAM unknown option: ', cargv
      i = i + 1
   endif
enddo

if (bad > 0) then
   write(io6,*) 'OLAM usage: ''exec name'' '
   write(io6,*) '  [-f ''Namelist file''] '
   stop 'bad command line arguments'
endif

write(io6,*) 'OLAM input namelist file: ',trim(name_name)

! Get unix tmp directory

call get_environment_variable("TMPDIR", cargv, length=i)
if ( i > 0 ) then
   if (cargv(i:i) == '/') then
      tmpdir = cargv(1:i-1)
   else
      tmpdir = cargv
   endif
else
   tmpdir = "/tmp"
endif

write(io6,*) "Using ", trim(tmpdir), " as the temporary directory."

! Initialize, execute, and end olam run

call olam_run(name_name)

! stop HDF5 library
#ifdef OLAM_HDF5_FORTRAN
call h5close_f(hdferr)
#endif

! If this run is parallel, finalize MPI and close io6 file

call olam_mpi_finalize()
if (iparallel == 1) then
   close(io6)
endif

stop 'olam_end'
end program main
