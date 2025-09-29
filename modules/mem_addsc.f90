Module mem_addsc

  implicit none

  ! Added scalar variables and tendencies

  Type addsc_vars
     real, allocatable :: sclp(:,:) ! somethings per kg_air
     real, allocatable :: sclt(:,:) ! somethings per m^3_air per sec
  End Type addsc_vars

  type (addsc_vars), allocatable :: addsc(:)

Contains

!===============================================================================

  subroutine alloc_addsc(mza,mwa,naddsc)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza,mwa,naddsc
    integer             :: iaddsc, iw

    ! Allocate arrays based on options (if necessary).

    allocate (addsc(naddsc))

    do iaddsc = 1,naddsc

       write(io6,*) 'alloc_addsc ', iaddsc, mza, mwa, naddsc
       allocate( addsc(iaddsc)%sclp(mza,mwa) )

       !$omp parallel do
       do iw = 1, mwa
          addsc(iaddsc)%sclp(:,iw) = rinit
       enddo
       !$omp end parallel do

    enddo

  end subroutine alloc_addsc

!===============================================================================

  subroutine dealloc_addsc(naddsc)

    implicit none

    integer, intent(in) :: naddsc
    integer             :: iaddsc

    !  Deallocate arrays

    if (allocated(addsc)) then

       do iaddsc = 1,naddsc
          if ( allocated( addsc(iaddsc)%sclp ) ) deallocate( addsc(iaddsc)%sclp )
       enddo

       deallocate(addsc)

    endif

  end subroutine dealloc_addsc

!===============================================================================

  subroutine filltab_addsc(naddsc)

    use var_tables, only: increment_vtable
    implicit none

    integer, intent(in) :: naddsc

    integer       :: iaddsc
    character(14) :: sname

    do iaddsc = 1,naddsc

       if (allocated (addsc(iaddsc)%sclp)) then

          write(sname,'(a4,i0.3)') 'SCLP', iaddsc
          call increment_vtable(sname, 'AW', mpt1=.true., rvar2=addsc(iaddsc)%sclp)

       endif

    enddo

  end subroutine filltab_addsc

End Module mem_addsc
