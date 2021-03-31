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
program main

use misc_coms,    only: io6, iparallel
use mem_para,     only: myrank, mgroupsize
use olam_mpi_atm, only: olam_mpi_init, olam_mpi_finalize, olam_mpi_barrier
use misc_coms,    only: tmpdir
use max_dims,     only: pathlen
use hdf5,         only: h5open_f, h5close_f
use oname_coms,   only: cmdlne_runtype, cmdlne_fields, numcf, nl
use sys_utils,    only: set_environment_variable

implicit none

character(pathlen) :: name_name = 'OLAMIN'
character(pathlen) :: cargv
character(30)      :: io6file
character(10)      :: dowrite
integer :: numarg
integer :: i, n
integer :: bad = 0
integer :: hdferr

! If run is sequential, default choice is to set io6 to standard output unit 6.

io6 = 6

! Determine if this run is parallel, and determine myrank and mgroupsize

call olam_mpi_init()

iparallel = 0
if (mgroupsize > 1) iparallel = 1

! For parallel HDF5
! call set_environment_variable("HDF5_USE_FILE_LOCKING", "FALSE")

! Initialize HDF5 library

call h5open_f(hdferr)

! Parse the command line arguments

numarg = command_argument_count()
if (myrank == 0) write(*,*) 'numarg:', numarg

allocate(cmdlne_fields(numarg/2 + 2))

i = 1
do while (i <= numarg)
   call get_command_argument(i,cargv)
   ! write(io6,*) 'args: ',i,cargv

   if (cargv(1:1) == '-') then

      if (cargv(2:2) == 'f') then

         call get_command_argument(i+1,name_name)
         if (len_trim(name_name) < 1) bad = bad + 1
         i = i + 2

      elseif (cargv(2:2) == 'r') then

         call get_command_argument(i+1,cmdlne_runtype)
         if (len_trim(cmdlne_runtype) < 1) bad = bad + 1
         i = i + 2

      elseif (cargv(2:2) == 'z') then

         if (numcf < size(cmdlne_fields)) then

            numcf = numcf + 1

            call get_command_argument(i+1,cmdlne_fields(numcf))
            if ( len_trim( cmdlne_fields(numcf) ) < 1  .or. &
                 scan( cmdlne_fields(numcf), '=') == 0 ) then
               bad = bad + 1
               numcf = numcf - 1
            endif

         else

            if (myrank == 0) write(*,*) "OLAM too many 'z' arguments"

         endif
         i = i + 2

      else
         i = i + 1
      endif
   else
      i = i + 1
   endif
enddo

if (bad > 0) then
   if (myrank == 0) write(*,*) 'OLAM usage: ''exec name'' '
   if (myrank == 0) write(*,*) '  [-f ''Namelist file''] '
   if (iparallel == 1) call olam_mpi_barrier()
   if (myrank == 0) then
      stop 'Bad command line arguments'
   else
      stop
   endif
endif

! Read Fortran namelist

if (myrank == 0) write(*,*) 'OLAM input namelist file: ',trim(name_name)
if (myrank == 0) write(*,'(/,a)') 'olammain calling namelist read'
call read_nl(name_name)

! If run is parallel, default choice is to attach output unit io6 to separate files

if (iparallel == 1 .and. myrank > 0) then

   io6 = 20

   if (nl%save_node_logs) then

      write (io6file,'(i10)') myrank
      io6file = 'o.io6_r'//trim(adjustl(io6file))

      ! Output file is replaced if it exists
      open(unit=io6, file=io6file, status='replace', form='formatted', action='write')

#ifdef __PGI
      ! Prevent Portland Group compiler from buffering io6
      call setvbuf3f(io6,1,1024)
#endif

   else

      io6file = '/dev/null'

      ! Work around an issue with IBM compiler and MPI
#if defined(__xlC__) || defined(__ibmxl__)
      inquire(file=io6file, write=dowrite, number=i)
      if (dowrite == 'YES' .and. i > 0) close(i)
#endif

      open(unit=io6, file=io6file, action='write')

   endif

endif

write(io6,'(/,a,i6)') ' myrank     = ',myrank
write(io6,'(  a,i6)') ' mgroupsize = ',mgroupsize
write(io6,'(  a,i6)') ' iparallel  = ',iparallel

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

call h5close_f(hdferr)

! If this run is parallel, finalize MPI and close io6 file

call olam_mpi_finalize()

write(io6,'(A)')
write(io6,'(A)') "olam_end"

if (iparallel == 1) then
   close(io6)
endif

end program main
