module sys_utils

contains

!===============================================================================

subroutine set_environment_variable(name, value, status)

  use,intrinsic :: iso_c_binding
  implicit none

  ! call setenv(3c) to set environment variable"

  character(*),      intent(in)  :: name
  character(*),      intent(in)  :: value
  integer, optional, intent(out) :: status
  integer                        :: loc_err

  interface
     integer(kind=c_int) function c_setenv(name,value,flag) bind(C,NAME="setenv")
       import :: c_int, c_char
       character(kind=c_char), intent(in)        :: name (*)
       character(kind=c_char), intent(in)        :: value(*)
       integer  (kind=c_int ), intent(in), value :: flag
     end function c_setenv
  end interface

  loc_err = c_setenv(trim(name)//c_null_char, trim(value)//c_null_char, 1)
  if (present(STATUS)) STATUS = loc_err

end subroutine set_environment_variable

!===============================================================================

subroutine get_command_as_string(command, str)

  use, intrinsic :: iso_c_binding
  implicit none

  interface

     type(c_ptr) function popen(command, mode) bind(C,name='popen')
       import :: c_char, c_ptr
       character(kind=c_char), dimension(*) :: command
       character(kind=c_char), dimension(*) :: mode
     end function popen

     type(c_ptr) function fgets(s, siz, stream) bind(C,name='fgets')
       import :: c_char, c_ptr, c_int
       character(kind=c_char), dimension(*) :: s
       integer  (kind=c_int),  value        :: siz
       type(c_ptr),            value        :: stream
     end function fgets

     integer(c_int) function pclose(stream) bind(C,name='pclose')
       import :: c_ptr, c_int
       type(c_ptr), value :: stream
     end function pclose

  end interface

  character(*), intent(in)    :: command
  character(*), intent(inout) :: str

  integer            :: s, i, j
  integer, parameter :: buffer_length = 12

  type(c_ptr) :: h
  integer     :: istat
  character(kind=c_char, len=len(str)+1) :: line

  str = ' '
  s = len_trim(command)
  if (s < 1) return

  h = c_null_ptr
  h = popen(command(1:s)//c_null_char, 'r'//c_null_char)

  if (c_associated(h)) then
     i = 0
     do while (c_associated(fgets(line,len(str)+1,h)) .and. i < len(str))
        j = max(index(line,c_null_char) - 1, 0)
        str = str(1:i) // line(1:j)
        i = i + j
     end do
     istat = pclose(h)
  end if

end subroutine get_command_as_string

!===============================================================================

subroutine sleep(seconds)

  use, intrinsic :: iso_c_binding
  implicit none

  interface
     integer(c_int) function c_sleep(seconds) bind(C,name='sleep')
       import :: c_ptr, c_int
       integer(c_int), value :: seconds
     end function c_sleep
  end interface

  integer, intent(in) :: seconds
  integer             :: i

  i = c_sleep(seconds)

end subroutine sleep

!===============================================================================

end module sys_utils
