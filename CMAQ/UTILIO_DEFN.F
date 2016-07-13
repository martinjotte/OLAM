      module utilio_defn
      
      INTEGER, PARAMETER :: FSREAD3 =   1 !  OPEN3() flag:  "old read-only" file
      INTEGER, PARAMETER :: FSRDWR3 =   2 !  "old read-write" file
      INTEGER, PARAMETER :: FSNEW3  =   3 !  "new (read-write)" file

      integer, parameter :: xstat1 = 1
      integer, parameter :: xstat2 = 2
      integer, parameter :: xstat3 = 3

!     integer, private   :: ncid

      contains

      INTEGER FUNCTION INDEX1 (NAME, N, NLIST)

C***********************************************************************
C Version "@(#)$Header$"
C EDSS/Models-3 I/O API.
C Copyright (C) 1992-2002 MCNC and Carlie J. Coats, Jr., and
C (C) 2003-2010 Baron Advanced Meteorological Systems, LLC.
C Distributed under the GNU LESSER GENERAL PUBLIC LICENSE version 2.1
C See file "LGPL.txt" for conditions of use.
C.........................................................................
C  INDEX1    subroutine body starts at line 53
C  INDEXINT1 subroutine body starts at line 99
C
C  FUNCTION:
C
C       Search for character-string or integer key NAME or IKEY in list NLIST
C       and return the subscript (1...N) at which it is found, or return 0
C       when not found in NLIST
C
C  PRECONDITIONS REQUIRED:
C       none
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       none
C
C  REVISION HISTORY:
C    INDEX1:
C       5/1988   Modified for ROMNET
C       9/1994   Modified for Models-3 by CJC
C    INDEXINT1:
C       Prototype 11/2004 by CJC:  MODULE M3UTILIO for I/O API v3
C       Modified   3/2006 by CJC:  moved INDEXINT1() to file "index1.f"
C
C       Modified  03/2010 by CJC: F9x changes for I/O API v3.1
C***********************************************************************

      IMPLICIT NONE

C.......   Arguments and their descriptions:

      CHARACTER*(*), INTENT(IN   ) :: NAME        !  Character string being searched for
      INTEGER      , INTENT(IN   ) :: N           !  Length of array to be searched
      CHARACTER*(*), INTENT(IN   ) :: NLIST(*)    !  array to be searched

C.......   Local variable:

      INTEGER       I   !  loop counter

C.....................................................................
C.......   begin body of INDEX1()

      DO I = 1, N
         IF ( NAME .EQ. NLIST( I ) ) THEN ! Found NAME in NLIST
            INDEX1 = I
            RETURN
         ENDIF
      END DO

      INDEX1 = 0                !  not found
      RETURN
      
      END FUNCTION INDEX1


      SUBROUTINE M3EXIT( CALLER, JDATE, JTIME, MSGTXT, EXITSTAT )

C***********************************************************************
C Version "@(#)$Header$"
C EDSS/Models-3 I/O API.
C Copyright (C) 1992-2002 MCNC and Carlie J. Coats, Jr.,
C (C) 2003-2010 by Baron Advanced Meteorological Systems.
C Distributed under the GNU LESSER GENERAL PUBLIC LICENSE version 2.1
C See file "LGPL.txt" for conditions of use.
C.........................................................................
C  subroutine body starts at line  76
C
C  FUNCTION:  Generate simple error messages for Models-3 core;
C             terminate program execution via CALL EXIT( EXITSTAT )
C
C  PRECONDITIONS REQUIRED:  JDATE:JTIME represented as YYYYDDD:HHMMSS
C
C  SUBROUTINES AND FUNCTIONS CALLED:  DT2STR, INIT3, SHUT3
C
C  REVISION  HISTORY:	
C	prototype 5/1992 by CJC
C       Revised   8/1996 to close currently-open POSIX-OK Fortran units.
C       Modified  1/1997 by CJC to trim trailing blanks from MSGTXT
C       Modified  2/1997 by CJC:  conditional definition of EXIT under AIX 
C       Modified 10/1998 by CJC:  Factor all output through m3msg2(), to
C       get around SGI OpenMP bug.
C       Caution:  messages may not be sequenced correctly when called by
C       multiple threads simultaneously
C       Modified  5/1999 by ALT for coupling-mode operation
C       Modified  4/2004 by CJC:  factor through M3PARAG() -- fixes
C       multi-thread sequencing problem.
C       Modified 03/2010 by CJC: F9x changes for I/O API v3.1
C***********************************************************************

!#ifdef   _AIX
!#define  EXIT exit_
!#endif

      IMPLICIT NONE

C...........   INCLUDES:

!       INCLUDE 'IODECL3.EXT'

C...........   ARGUMENTS and their descriptions:

      CHARACTER*(*), INTENT(IN   ) :: CALLER !  name of the caller
      INTEGER      , INTENT(IN   ) :: JDATE, JTIME !  model date&time for the error
      CHARACTER*(*), INTENT(IN   ) :: MSGTXT !  error message
      INTEGER      , INTENT(IN   ) :: EXITSTAT !  exit status for program

                 
C...........   EXTERNAL FUNCTIONS and their descriptions:

!       CHARACTER*24, EXTERNAL :: DT2STR


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

!        CHARACTER*256   MESG( 6 )
!        INTEGER         LINE, IDEV, ISTAT
!        LOGICAL         OFLAG


C***********************************************************************
C   begin body of subroutine  M3EXIT

!#ifdef IOAPICPL
!        IF ( EXITSTAT .NE. 0 ) THEN
!            CALL SETSYNCHRO3V(1)
!        END IF
!#endif

!        LINE = 1
!        MESG( LINE ) = ' '
!
!        IF ( .NOT. SHUT3() ) THEN
!            LINE = LINE + 1
!            MESG( LINE ) ='Could not shut down I/O API files correctly'
!        END IF
!
!        LINE = LINE + 1
!        IF ( EXITSTAT .NE. 0 ) THEN
!            MESG( LINE ) = '*** ERROR ABORT in subroutine ' // CALLER
!        ELSE    !  errstat = 0:  successful completion
!            MESG( LINE ) = '--->> Normal Completion of program '//CALLER
!        END IF
!
!        LINE = LINE + 1
!        MESG( LINE ) = MSGTXT
!
!        IF ( JDATE .GT. 0  .OR.  JTIME .GT. 0 ) THEN
!        LINE = LINE + 1
!            WRITE( MESG( LINE ), '( A, 2X, A, 2X, A, I7, A, I6.6, A )' ) 
!     &      'Date and time', DT2STR( JDATE, JTIME ),  
!     &      '(', JDATE, ':', JTIME, ')'
!        END IF
!
!        LINE = LINE + 1
!        MESG( LINE ) = ' '
!        CALL M3PARAG( LINE, MESG )
!
!        DO  IDEV = 10, 99
!            INQUIRE( UNIT = IDEV, OPENED = OFLAG, IOSTAT = ISTAT )
!            !!  paranoia to deal with certain rare situations on SGI
!            IF ( ISTAT .EQ. 0 ) THEN
!                IF ( OFLAG ) CLOSE( UNIT = IDEV, IOSTAT = ISTAT )
!                IF ( ISTAT .NE. 0 ) THEN
!                    WRITE( MESG, '( A, 2X, I4, A, I10 )' )
!     &                  'Error closing unit', IDEV,
!     &                  ' status =', ISTAT
!                    CALL M3MSG2( MESG )
!                END IF
!            END IF
!        END DO

!       CALL EXIT( EXITSTAT )

      write(*,'(A)'      ) " Error in " // trim(caller)

      if (jdate /= 0 .and. jtime /= 0)
     &     write(*,'(2(A,I0))') " At date ", JDATE, " and time ", JTIME

      if (len_trim(msgtxt) > 0) then
         write(*,'(1x,A)') trim(msgtxt)
      endif

      if     (exitstat == 1) then
         stop 1
      elseif (exitstat == 2) then
         stop 2
      elseif (exitstat == 3) then
         stop 3
      else
         stop
      endif

      END SUBROUTINE M3EXIT



!      logical function open3( fname, fstat, pname )
!      
!      use netcdf
!      implicit none
!
!      CHARACTER(*), INTENT( IN ) :: FNAME  !  logical name of file to be opened
!      INTEGER     , INTENT( IN ) :: FSTAT  !  read-only, read-write, new, or unknown
!      CHARACTER(*), INTENT( IN ) :: PNAME  !  name of calling program
!
!      logical                    :: exists
!      integer                    :: mode
!      integer                    :: rcode
!
!      open3 = .true.
!
!      if (fstat == fsread3 .or. fstat == fsrdwr3) then
!         
!         inquire(file=fname, exist=exists)
!         if (.not. exists) then
!            write(*,*) 'open3: error opening input file ' // trim(fname)
!            write(*,*) '   called from routine ' // pname
!            open3 = .false.
!            return
!         endif
!
!         if (fstat == fsread3) then
!            mode = nf90_nowrite
!         else
!            mode = nf90_write
!         endif
!
!         rcode = nf90_open( fname, mode, ncid )
!
!         if ( rcode /= nf90_noerr ) then
!            write(*,*) 'open3: error opening input file ' // trim(fname)
!            write(*,*) '   called from routine ' // pname
!            open3 = .false.
!            return
!         endif
!
!      elseif (fstat == fsnew3 ) then
!
!         mode = nf90_clobber
!
!         rcode = nf90_create( fname, mode, ncid )
!
!         if ( rcode /= nf90_noerr ) then
!            write(*,*) 'open3: error creating file ' // trim(fname)
!            write(*,*) '   called from routine ' // pname
!            open3 = .false.
!            return
!         endif
!
!      else
!
!         write(*,*) 'open3: wrong file status argument'
!         write(*,*) '   called from routine ' // pname
!         open3 = .false.
!         return
!
!      endif
!
!      end function open3



        SUBROUTINE  UPCASE ( BUFFER )

C***********************************************************************
C Version "@(#)$Header$"
C EDSS/Models-3 I/O API.
C Copyright (C) 1992-2002 MCNC and Carlie J. Coats, Jr.,
C (C) 2003-2010 by Baron Advanced Meteorological Systems.
C Distributed under the GNU LESSER GENERAL PUBLIC LICENSE version 2.1
C See file "LGPL.txt" for conditions of use.
C.........................................................................
C  subroutine body starts at line  49
C
C  FUNCTION:  upcase the text in BUFFER
C
C  PRECONDITIONS REQUIRED:  text is ASCII
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 1/91 by CJC
C
C       Modified 03/2010 by CJC: F90 changes for I/O API v3.1
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*), INTENT(INOUT) :: BUFFER


C...........   PARAMETER:  ASCII for 'a', 'z', 'A'

        INTEGER, PARAMETER :: IA    =  97
        INTEGER, PARAMETER :: IZ    = 122
        INTEGER, PARAMETER :: AADIF =  32


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER       I, L
        INTEGER       C


C***********************************************************************
C   begin body of subroutine  UPCASE

        L  =  LEN ( BUFFER )
        DO  111  I = 1 , L
            C = ICHAR ( BUFFER ( I:I ) )
            IF ( C .GE. IA  .AND.  C .LE. IZ ) THEN
                BUFFER ( I:I ) = CHAR ( C - AADIF )
            END IF
111     CONTINUE        !  end loop on I

        RETURN

        END SUBROUTINE  UPCASE


      SUBROUTINE M3WARN( CALLER, JDATE, JTIME, MSGTXT )
      use misc_coms, only: io6
      IMPLICIT NONE

C...........   INCLUDES:

!       INCLUDE 'IODECL3.EXT'

C...........   ARGUMENTS and their descriptions:

      CHARACTER*(*), INTENT(IN   ) :: CALLER !  name of the caller
      INTEGER      , INTENT(IN   ) :: JDATE, JTIME !  model date&time for the error
      CHARACTER*(*), INTENT(IN   ) :: MSGTXT !  warning message
                 
C...........   EXTERNAL FUNCTIONS and their descriptions:

!       CHARACTER*24, EXTERNAL :: DT2STR


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

!        CHARACTER*256   MESG( 6 )
!        INTEGER         LINE, IDEV, ISTAT
!        LOGICAL         OFLAG


C***********************************************************************
C   begin body of subroutine  M3EXIT

!#ifdef IOAPICPL
!        IF ( EXITSTAT .NE. 0 ) THEN
!            CALL SETSYNCHRO3V(1)
!        END IF
!#endif

!        LINE = 1
!        MESG( LINE ) = ' '
!
!        IF ( .NOT. SHUT3() ) THEN
!            LINE = LINE + 1
!            MESG( LINE ) ='Could not shut down I/O API files correctly'
!        END IF
!
!        LINE = LINE + 1
!        IF ( EXITSTAT .NE. 0 ) THEN
!            MESG( LINE ) = '*** ERROR ABORT in subroutine ' // CALLER
!        ELSE    !  errstat = 0:  successful completion
!            MESG( LINE ) = '--->> Normal Completion of program '//CALLER
!        END IF
!
!        LINE = LINE + 1
!        MESG( LINE ) = MSGTXT
!
!        IF ( JDATE .GT. 0  .OR.  JTIME .GT. 0 ) THEN
!        LINE = LINE + 1
!            WRITE( MESG( LINE ), '( A, 2X, A, 2X, A, I7, A, I6.6, A )' ) 
!     &      'Date and time', DT2STR( JDATE, JTIME ),  
!     &      '(', JDATE, ':', JTIME, ')'
!        END IF
!
!        LINE = LINE + 1
!        MESG( LINE ) = ' '
!        CALL M3PARAG( LINE, MESG )
!
!        DO  IDEV = 10, 99
!            INQUIRE( UNIT = IDEV, OPENED = OFLAG, IOSTAT = ISTAT )
!            !!  paranoia to deal with certain rare situations on SGI
!            IF ( ISTAT .EQ. 0 ) THEN
!                IF ( OFLAG ) CLOSE( UNIT = IDEV, IOSTAT = ISTAT )
!                IF ( ISTAT .NE. 0 ) THEN
!                    WRITE( MESG, '( A, 2X, I4, A, I10 )' )
!     &                  'Error closing unit', IDEV,
!     &                  ' status =', ISTAT
!                    CALL M3MSG2( MESG )
!                END IF
!            END IF
!        END DO

!       CALL EXIT( EXITSTAT )

      write(io6,'(A)'      ) " Warning in " // trim(caller)

      if (jdate /= 0 .and. jtime /= 0)
     &     write(io6,'(2(A,I0))') " At date ", JDATE, " and time ", JTIME

      if (len_trim(msgtxt) > 0) then
         write(io6,'(1x,A)') trim(msgtxt)
      endif

      END SUBROUTINE M3WARN



      SUBROUTINE M3MESG( MSGTXT )
      use misc_coms, only: io6
      IMPLICIT NONE

C...........   INCLUDES:

!       INCLUDE 'IODECL3.EXT'

C...........   ARGUMENTS and their descriptions:

      CHARACTER*(*), INTENT(IN   ) :: MSGTXT !  warning message

      if (len_trim(msgtxt) > 0) then
         write(io6,'(1x,A)') trim(msgtxt)
      endif

      END SUBROUTINE M3MESG

      

      INTEGER FUNCTION INIT3 ()
      use misc_coms, only: io6

      implicit none

      init3 = io6
      END FUNCTION INIT3

      end module utilio_defn


      
      


