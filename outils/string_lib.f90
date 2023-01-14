module string_lib

contains

!===============================================================================

function to_upper(string) result(upper)

  implicit none

  character(len=*), intent(in) :: string
  character(len=len(string))   :: upper

  integer :: j

  do j = 1, len(string)
     if (string(j:j) >= "a" .and. string(j:j) <= "z") then
        upper(j:j) = achar(iachar(string(j:j)) - 32)
     else
        upper(j:j) = string(j:j)
     endif
  enddo

end function to_upper

!===============================================================================

function to_lower(string) result(lower)

  implicit none

  character(len=*), intent(in) :: string
  character(len=len(string))   :: lower

  integer :: j

  do j = 1, len(string)
     if (string(j:j) >= "A" .and. string(j:j) <= "Z") then
        lower(j:j) = achar(iachar(string(j:j)) + 32)
     else
        lower(j:j) = string(j:j)
     endif
  enddo

end function to_lower

!===============================================================================

subroutine lowercase(string)

  implicit none

  character(len=*), intent(inout) :: string

  integer :: j

  do j = 1, len_trim(string)
     if (string(j:j) >= "A" .and. string(j:j) <= "Z") then
        string(j:j) = achar(iachar(string(j:j)) + 32)
     endif
  enddo

end subroutine lowercase

!===============================================================================

subroutine strip_char(string,zz)

  implicit none

  character(*), intent(inout) :: string
  character(1), intent(in)    :: zz
  integer                     :: nn, ll, i

  nn = 1
  ll = len_trim(string)

  do while(nn <= ll)

     do i = nn, ll
        if (string(i:i) == zz) then
           string(i:ll-1) = string(i+1:ll)
           string(ll:ll) = ' '
           ll = len_trim(string)
           exit
        endif
        nn = nn + 1
     enddo

  enddo

end subroutine strip_char

!===============================================================================

subroutine strip_chars(string,zz)

  implicit none

  character(*), intent(inout) :: string
  character(1), intent(in)    :: zz(:)
  integer                     :: nn, ll, i

  nn = 1
  ll = len_trim(string)

  do while(nn <= ll)

     do i = nn, ll
        if (any(string(i:i) == zz)) then
           string(i:ll-1) = string(i+1:ll)
           string(ll:ll) = ' '
           ll = len_trim(string)
           exit
        endif
        nn = nn + 1
     enddo

  enddo

end subroutine strip_chars

!===============================================================================

end module string_lib
