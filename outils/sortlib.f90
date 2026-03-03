module sortlib

  ! This module contains a number of routines for in-place sorting of real,
  ! integer, and double precision arrays into ascending order using an insertion
  ! sort algorithm and/or a shell sort. Insertion sort is extremely fast for
  ! small sizes but its speed falls off with order(N^2), whereas shell sort is
  ! much more efficient for large sorts.
  !
  ! The default sort routine here uses an insertion sort for small arrays with
  ! size < 95 and a shell sort for larger sizes, but there are routines to do
  ! either a shell sort or insertion short.
  !
  ! The generic interfaces sort, shell_sort, and insertion_sort are overloaded
  ! to handle real, double precision, or integer arrays and to sort additional
  ! real, double precision, or integer arrays based on the first input array.
  !
  ! A useful option is to define an integer array with the original indices:
  ! iarray(1:n) = [ (i,i=1,n) ]
  ! and then calling sort_fltint(array,iarray) will return in iarray the new
  ! ordering so that existing arrays can be converted to match the new sorted
  ! order of array.

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  private :: r8

  interface sort
     module procedure  &
          sort_flt,    &
          sort_fltint, &
          sort_flt2,   &
          sort_flt3,   &
          sort_dbl,    &
          sort_dblint, &
          sort_dbl2,   &
          sort_dbl3,   &
          sort_int,    &
          sort_int2,   &
          sort_int3
  end interface sort

  interface sort_rev
     module procedure      &
          sort_flt_rev,    &
          sort_fltint_rev, &
          sort_flt2_rev,   &
          sort_flt3_rev,   &
          sort_dbl_rev,    &
          sort_dblint_rev, &
          sort_dbl2_rev,   &
          sort_dbl3_rev,   &
          sort_int_rev,    &
          sort_int2_rev,   &
          sort_int3_rev
  end interface sort_rev

  interface insertion_sort
     module procedure            &
          insertion_sort_flt,    &
          insertion_sort_fltint, &
          insertion_sort_flt2,   &
          insertion_sort_flt3,   &
          insertion_sort_dbl,    &
          insertion_sort_dblint, &
          insertion_sort_dbl2,   &
          insertion_sort_dbl3,   &
          insertion_sort_int,    &
          insertion_sort_int2,   &
          insertion_sort_int3
  end interface insertion_sort

  interface insertion_sort_rev
     module procedure                &
          insertion_sort_flt_rev,    &
          insertion_sort_fltint_rev, &
          insertion_sort_flt2_rev,   &
          insertion_sort_flt3_rev,   &
          insertion_sort_dbl_rev,    &
          insertion_sort_dblint_rev, &
          insertion_sort_dbl2_rev,   &
          insertion_sort_dbl3_rev,   &
          insertion_sort_int_rev,    &
          insertion_sort_int2_rev,   &
          insertion_sort_int3_rev
  end interface insertion_sort_rev

  interface shell_sort
     module procedure        &
          shell_sort_flt,    &
          shell_sort_fltint, &
          shell_sort_flt2,   &
          shell_sort_flt3,   &
          shell_sort_dbl,    &
          shell_sort_dblint, &
          shell_sort_dbl2,   &
          shell_sort_dbl3,   &
          shell_sort_int,    &
          shell_sort_int2,   &
          shell_sort_int3
  end interface shell_sort

  interface shell_sort_rev
     module procedure            &
          shell_sort_flt_rev,    &
          shell_sort_fltint_rev, &
          shell_sort_flt2_rev,   &
          shell_sort_flt3_rev,   &
          shell_sort_dbl_rev,    &
          shell_sort_dblint_rev, &
          shell_sort_dbl2_rev,   &
          shell_sort_dbl3_rev,   &
          shell_sort_int_rev,    &
          shell_sort_int2_rev,   &
          shell_sort_int3_rev
  end interface shell_sort_rev


contains



subroutine insertion_sort_flt(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array(:)
  integer                         :: n, i, j
  real                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) <= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_flt



subroutine insertion_sort_fltint(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an insertion sort algorithm

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)
  integer                            :: n, i, j
  real                               :: temp
  integer                            :: itmp

  n = size(array)
  if (size(indic) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp = array(i)
     itmp = indic(i)
     do j = i-1, 1, -1
        if (array(j) <= temp) exit
        array(j+1) = array(j)
        indic(j+1) = indic(j)
     enddo
     array(j+1) = temp
     indic(j+1) = itmp
  enddo

end subroutine insertion_sort_fltint



subroutine insertion_sort_flt2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  integer                         :: n, i, j
  real                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_flt2



subroutine insertion_sort_flt3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)
  integer                         :: n, i, j
  real                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_flt3



subroutine insertion_sort_dbl(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer                             :: n, i, j
  real(r8)                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) <= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_dbl



subroutine insertion_sort_dblint(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)
  integer                             :: n, i, j
  real(r8)                            :: temp
  integer                             :: itmp

  n = size(array)
  if (size(indic) < n) stop 'invalid array size in dblsort'

  do i = 2, n
     temp = array(i)
     itmp = indic(i)
     do j = i-1, 1, -1
        if (array(j) <= temp) exit
        array(j+1) = array(j)
        indic(j+1) = indic(j)
     enddo
     array(j+1) = temp
     indic(j+1) = itmp
  enddo

end subroutine insertion_sort_dblint



subroutine insertion_sort_dbl2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  integer                             :: n, i, j
  real(r8)                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in dblsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_dbl2



subroutine insertion_sort_dbl3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)
  integer                             :: n, i, j
  real(r8)                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_dbl3



subroutine insertion_sort_int(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array(:)
  integer                            :: n, i, j
  integer                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) <= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_int



subroutine insertion_sort_int2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer                            :: n, i, j
  integer                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in intsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_int2



subroutine insertion_sort_int3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)
  integer                            :: n, i, j
  integer                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) <= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_int3



subroutine shell_sort_flt(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array(:)
  integer                         :: i, interval, last

  last = size(array)
  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt(array(i:last:interval))
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt



subroutine shell_sort_fltint(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an shell sort algorithm

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)
  integer                            :: i, interval, last
  real                               :: temp
  integer                            :: itmp

  last = size(array)
  if (size(indic) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_fltint( array(i:last:interval), &
                                    indic(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_fltint



subroutine shell_sort_flt2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  integer                         :: i, interval, last
  real                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt2( array1(i:last:interval), &
                                  array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt2



subroutine shell_sort_flt3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)
  integer                         :: i, interval, last
  real                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in fltsort'
  if (size(array3) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt3( array1(i:last:interval), &
                                  array2(i:last:interval), &
                                  array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt3



subroutine shell_sort_dbl(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp

  last = size(array)

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl( array(i:last:interval) )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl



subroutine shell_sort_dblint(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp
  integer                             :: itmp

  last = size(array)
  if (size(indic) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dblint( array(i:last:interval), &
                                    indic(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dblint



subroutine shell_sort_dbl2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl2( array1(i:last:interval), &
                                  array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl2



subroutine shell_sort_dbl3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in dblsort'
  if (size(array3) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl3( array1(i:last:interval), &
                                  array2(i:last:interval), &
                                  array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl3



subroutine shell_sort_int(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array(:)
  integer                            :: i, interval, last
  integer                            :: temp

  last = size(array)

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int( array(i:last:interval) )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int



subroutine shell_sort_int2(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer                            :: i, interval, last
  integer                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in intsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int2( array1(i:last:interval), &
                                  array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int2



subroutine shell_sort_int3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)
  integer                            :: i, interval, last
  integer                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in intsort'
  if (size(array3) < last) stop 'invalid array size in intsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int3( array1(i:last:interval), &
                                  array2(i:last:interval), &
                                  array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int3



subroutine sort_flt(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  real, contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_flt(array)
  else
     call shell_sort_flt(array)
  endif

end subroutine sort_flt



subroutine sort_fltint(array,indic)

  ! Sort array and indic into ascending order by array using insertion sort
  ! for a small array and a shell sort for a larger array

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)

  if (size(array) < 95) then
     call insertion_sort_fltint(array,indic)
  else
     call shell_sort_fltint(array,indic)
  endif

end subroutine sort_fltint



subroutine sort_flt2(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using an insertion
  ! sort algorithm for a small array and a shell sort for a larger array

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_flt2(array1,array2)
  else
     call shell_sort_flt2(array1,array2)
  endif

end subroutine sort_flt2



subroutine sort_flt3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm for a small array1 and a
  ! shell sort for a larger array1

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_flt3(array1,array2,array3)
  else
     call shell_sort_flt3(array1,array2,array3)
  endif

end subroutine sort_flt3



subroutine sort_dbl(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_dbl(array)
  else
     call shell_sort_dbl(array)
  endif

end subroutine sort_dbl



subroutine sort_dblint(array,indic)

  ! Sort array and indic into ascending order by array using insertion sort
  ! for a small array and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)

  if (size(array) < 95) then
     call insertion_sort_dblint(array,indic)
  else
     call shell_sort_dblint(array,indic)
  endif

end subroutine sort_dblint



subroutine sort_dbl2(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using an insertion
  ! sort algorithm for a small array and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_dbl2(array1,array2)
  else
     call shell_sort_dbl2(array1,array2)
  endif

end subroutine sort_dbl2



subroutine sort_dbl3(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm for a small array1 and a
  ! shell sort for a larger array1

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_dbl3(array1,array2,array3)
  else
     call shell_sort_dbl3(array1,array2,array3)
  endif

end subroutine sort_dbl3



subroutine sort_int(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  integer, contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_int(array)
  else
     call shell_sort_int(array)
  endif

end subroutine sort_int



subroutine sort_int2(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using insertion
  ! sort for a small array1 and a shell sort for larger array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_int2(array1,array2)
  else
     call shell_sort_int2(array1,array2)
  endif

end subroutine sort_int2



subroutine sort_int3(array1,array2,array3)

  ! Sort array1, array2, and array3 into ascending order by array1 using
  ! insertion sort for a small array1 and a shell sort for larger array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_int3(array1,array2,array3)
  else
     call shell_sort_int3(array1,array2,array3)
  endif

end subroutine sort_int3



subroutine insertion_sort_flt_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array(:)
  integer                         :: n, i, j
  real                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) >= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_flt_rev



subroutine insertion_sort_fltint_rev(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an insertion sort algorithm

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)
  integer                            :: n, i, j
  real                               :: temp
  integer                            :: itmp

  n = size(array)
  if (size(indic) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp = array(i)
     itmp = indic(i)
     do j = i-1, 1, -1
        if (array(j) >= temp) exit
        array(j+1) = array(j)
        indic(j+1) = indic(j)
     enddo
     array(j+1) = temp
     indic(j+1) = itmp
  enddo

end subroutine insertion_sort_fltint_rev



subroutine insertion_sort_flt2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  integer                         :: n, i, j
  real                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_flt2_rev



subroutine insertion_sort_flt3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)
  integer                         :: n, i, j
  real                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_flt3_rev



subroutine insertion_sort_dbl_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer                             :: n, i, j
  real(r8)                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) >= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_dbl_rev



subroutine insertion_sort_dblint_rev(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)
  integer                             :: n, i, j
  real(r8)                            :: temp
  integer                             :: itmp

  n = size(array)
  if (size(indic) < n) stop 'invalid array size in dblsort'

  do i = 2, n
     temp = array(i)
     itmp = indic(i)
     do j = i-1, 1, -1
        if (array(j) >= temp) exit
        array(j+1) = array(j)
        indic(j+1) = indic(j)
     enddo
     array(j+1) = temp
     indic(j+1) = itmp
  enddo

end subroutine insertion_sort_dblint_rev



subroutine insertion_sort_dbl2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  integer                             :: n, i, j
  real(r8)                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in dblsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_dbl2_rev



subroutine insertion_sort_dbl3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)
  integer                             :: n, i, j
  real(r8)                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_dbl3_rev



subroutine insertion_sort_int_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an insertion sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array(:)
  integer                            :: n, i, j
  integer                            :: temp

  n = size(array)

  do i = 2, n
     temp = array(i)
     do j = i-1, 1, -1
        if (array(j) >= temp) exit
        array(j+1) = array(j)
     enddo
     array(j+1) = temp
  enddo

end subroutine insertion_sort_int_rev



subroutine insertion_sort_int2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an insertion sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer                            :: n, i, j
  integer                            :: temp1, temp2

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in intsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
  enddo

end subroutine insertion_sort_int2_rev



subroutine insertion_sort_int3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)
  integer                            :: n, i, j
  integer                            :: temp1, temp2, temp3

  n = size(array1)
  if (size(array2) < n) stop 'invalid array size in fltsort'
  if (size(array3) < n) stop 'invalid array size in fltsort'

  do i = 2, n
     temp1 = array1(i)
     temp2 = array2(i)
     temp3 = array3(i)
     do j = i-1, 1, -1
        if (array1(j) >= temp1) exit
        array1(j+1) = array1(j)
        array2(j+1) = array2(j)
        array3(j+1) = array3(j)
     enddo
     array1(j+1) = temp1
     array2(j+1) = temp2
     array3(j+1) = temp3
  enddo

end subroutine insertion_sort_int3_rev



subroutine shell_sort_flt_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array(:)
  integer                         :: i, interval, last

  last = size(array)
  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt_rev(array(i:last:interval))
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt_rev



subroutine shell_sort_fltint_rev(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an shell sort algorithm

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)
  integer                            :: i, interval, last
  real                               :: temp
  integer                            :: itmp

  last = size(array)
  if (size(indic) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_fltint_rev( array(i:last:interval), &
                                        indic(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_fltint_rev



subroutine shell_sort_flt2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  integer                         :: i, interval, last
  real                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt2_rev( array1(i:last:interval), &
                                      array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt2_rev



subroutine shell_sort_flt3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an shell sort algorithm

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)
  integer                         :: i, interval, last
  real                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in fltsort'
  if (size(array3) < last) stop 'invalid array size in fltsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_flt3_rev( array1(i:last:interval), &
                                      array2(i:last:interval), &
                                      array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_flt3_rev



subroutine shell_sort_dbl_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp

  last = size(array)

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl_rev( array(i:last:interval) )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl_rev



subroutine shell_sort_dblint_rev(array,indic)

  ! Sort n floating point numbers of array and n integers of array2 into
  ! ascending order by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp
  integer                             :: itmp

  last = size(array)
  if (size(indic) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dblint_rev( array(i:last:interval), &
                                        indic(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dblint_rev



subroutine shell_sort_dbl2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl2_rev( array1(i:last:interval), &
                                      array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl2_rev



subroutine shell_sort_dbl3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an shell sort algorithm

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)
  integer                             :: i, interval, last
  real(r8)                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in dblsort'
  if (size(array3) < last) stop 'invalid array size in dblsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_dbl3_rev( array1(i:last:interval), &
                                      array2(i:last:interval), &
                                      array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_dbl3_rev



subroutine shell_sort_int_rev(array)

  ! Sort n floating point numbers of array into ascending order,
  ! using an shell sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array(:)
  integer                            :: i, interval, last
  integer                            :: temp

  last = size(array)

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int_rev( array(i:last:interval) )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int_rev



subroutine shell_sort_int2_rev(array1,array2)

  ! Sort n floating point numbers of array1 and array2 into ascending order
  ! by array1 using an shell sort algorithm

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer                            :: i, interval, last
  integer                            :: temp1, temp2

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in intsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int2_rev( array1(i:last:interval), &
                                      array2(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int2_rev



subroutine shell_sort_int3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)
  integer                            :: i, interval, last
  integer                            :: temp1, temp2, temp3

  last = size(array1)
  if (size(array2) < last) stop 'invalid array size in intsort'
  if (size(array3) < last) stop 'invalid array size in intsort'

  interval = 1

  do while ( interval < last/3 )
     interval = 3*interval + 1
  enddo

  do while ( interval > 0 )
     do i = 1, interval
        call insertion_sort_int3_rev( array1(i:last:interval), &
                                      array2(i:last:interval), &
                                      array3(i:last:interval)  )
     enddo
     interval = (interval-1)/3
  enddo

end subroutine shell_sort_int3_rev



subroutine sort_flt_rev(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  real, contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_flt_rev(array)
  else
     call shell_sort_flt_rev(array)
  endif

end subroutine sort_flt_rev



subroutine sort_fltint_rev(array,indic)

  ! Sort array and indic into ascending order by array using insertion sort
  ! for a small array and a shell sort for a larger array

  implicit none

  real,    contiguous, intent(inout) :: array(:)
  integer, contiguous, intent(inout) :: indic(:)

  if (size(array) < 95) then
     call insertion_sort_fltint_rev(array,indic)
  else
     call shell_sort_fltint_rev(array,indic)
  endif

end subroutine sort_fltint_rev



subroutine sort_flt2_rev(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using an insertion
  ! sort algorithm for a small array and a shell sort for a larger array

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_flt2_rev(array1,array2)
  else
     call shell_sort_flt2_rev(array1,array2)
  endif

end subroutine sort_flt2_rev



subroutine sort_flt3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm for a small array1 and a
  ! shell sort for a larger array1

  implicit none

  real, contiguous, intent(inout) :: array1(:)
  real, contiguous, intent(inout) :: array2(:)
  real, contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_flt3_rev(array1,array2,array3)
  else
     call shell_sort_flt3_rev(array1,array2,array3)
  endif

end subroutine sort_flt3_rev



subroutine sort_dbl_rev(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_dbl_rev(array)
  else
     call shell_sort_dbl_rev(array)
  endif

end subroutine sort_dbl_rev



subroutine sort_dblint_rev(array,indic)

  ! Sort array and indic into ascending order by array using insertion sort
  ! for a small array and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array(:)
  integer,  contiguous, intent(inout) :: indic(:)

  if (size(array) < 95) then
     call insertion_sort_dblint_rev(array,indic)
  else
     call shell_sort_dblint_rev(array,indic)
  endif

end subroutine sort_dblint_rev



subroutine sort_dbl2_rev(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using an insertion
  ! sort algorithm for a small array and a shell sort for a larger array

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_dbl2_rev(array1,array2)
  else
     call shell_sort_dbl2_rev(array1,array2)
  endif

end subroutine sort_dbl2_rev



subroutine sort_dbl3_rev(array1,array2,array3)

  ! Sort n floating point numbers of array1, array2, and array3 into ascending
  ! order by array1 using an insertion sort algorithm for a small array1 and a
  ! shell sort for a larger array1

  implicit none

  real(r8), contiguous, intent(inout) :: array1(:)
  real(r8), contiguous, intent(inout) :: array2(:)
  real(r8), contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_dbl3_rev(array1,array2,array3)
  else
     call shell_sort_dbl3_rev(array1,array2,array3)
  endif

end subroutine sort_dbl3_rev



subroutine sort_int_rev(array)

  ! Sort array into ascending order using insertion sort for a small array
  ! and a shell sort for a larger array

  implicit none

  integer, contiguous, intent(inout) :: array(:)

  if (size(array) < 95) then
     call insertion_sort_int_rev(array)
  else
     call shell_sort_int_rev(array)
  endif

end subroutine sort_int_rev



subroutine sort_int2_rev(array1,array2)

  ! Sort array1 and array2 into ascending order by array1 using insertion
  ! sort for a small array1 and a shell sort for larger array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)

  if (size(array1) < 95) then
     call insertion_sort_int2_rev(array1,array2)
  else
     call shell_sort_int2_rev(array1,array2)
  endif

end subroutine sort_int2_rev



subroutine sort_int3_rev(array1,array2,array3)

  ! Sort array1, array2, and array3 into ascending order by array1 using
  ! insertion sort for a small array1 and a shell sort for larger array1

  implicit none

  integer, contiguous, intent(inout) :: array1(:)
  integer, contiguous, intent(inout) :: array2(:)
  integer, contiguous, intent(inout) :: array3(:)

  if (size(array1) < 95) then
     call insertion_sort_int3_rev(array1,array2,array3)
  else
     call shell_sort_int3_rev(array1,array2,array3)
  endif

end subroutine sort_int3_rev



end module sortlib
