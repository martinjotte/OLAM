module sys_utils


contains


subroutine set_environment_variable(NAME, VALUE, STATUS)

  use,intrinsic :: iso_c_binding
  implicit none

  ! call setenv(3c) to set environment variable"

   character(len=*)               :: NAME
   character(len=*)               :: VALUE
   integer, optional, intent(out) :: STATUS
   integer                        :: loc_err

   interface
      integer(kind=c_int) function c_setenv(c_name,c_VALUE) bind(C,NAME="setenv")
        use,intrinsic :: iso_c_binding
        character(kind=c_char)   :: c_name(*)
        character(kind=c_char)   :: c_VALUE(*)
      end function c_setenv
   end interface

   loc_err = c_setenv(str2arr(trim(NAME)),str2arr(VALUE))
   if (present(STATUS)) STATUS = loc_err

end subroutine set_environment_variable



pure function arr2str(array) result (string)

  implicit none

  ! copies null-terminated char array to string"

  character(1), contiguous, intent(in)  :: array(:)
  character(size(array))                :: string
  integer                               :: i

  string=' '
  do i = 1, size(array)
     if (array(i).eq.char(0)) then
        exit
     else
        string(i:i) = array(i)
     endif
  enddo

end function arr2str



pure function str2arr(string) result (array)

  use,intrinsic :: iso_c_binding
  implicit none

  ! copies string to null terminated char array"

  character(len=*), intent(in)  :: string
  character(len=1, kind=c_char) :: array(len(string)+1)
  integer                       :: i

  do i = 1, len_trim(string)
     array(i) = string(i:i)
  enddo
  array(size(array)) = c_null_char

end function str2arr



end module sys_utils
