module cmaq_cloud

  integer, parameter :: ncats = 9

  character(len=10)  :: cldtypes(ncats) = (/ 'cloud    ', &
                                             'drizzle  ', &
                                             'rain     ', &
                                             'hail     ', &
                                             'aggregate', &
                                             'hollowcol', &
                                             'solidcol ', &
                                             'hexplate ', &
                                             'rosette  '  /)
  type cloud_vars
     character(len=10) :: name = ''
     integer           :: num  = 0
     real              :: start
     real              :: end
     real              :: delr
     real, allocatable :: extsw(:,:)
     real, allocatable :: ssasw(:,:)
     real, allocatable :: asysw(:,:)
     real, allocatable :: abslw(:,:)
  end type cloud_vars

  type (cloud_vars) :: cloud_props(ncats)

contains

  subroutine cmaq_cld_optics_init()

    use hdf5_utils, only: shdf5_open, shdf5_close
    use oname_coms, only: nl
    use max_dims,   only: pathlen

    implicit none

    integer :: i
    logical :: exists

    character(20),     parameter :: omic_cmaq_file = "omic2cmaq.h5"
    character(pathlen)           :: inputfile

    inputfile = trim(nl%rrtmg_datadir) // "/" // trim(omic_cmaq_file)

    inquire(file=inputfile, exist=exists)
    if (.not. exists) then
       write(*,*) "rrtmg_cmaq_cloud_init: Error opening data file " // trim(inputfile)
       stop       "RRTMg shortwave datafile cannot be found"
    endif

    call shdf5_open(inputfile, 'R', trypario=.true.)

    do i = 1, ncats

       cloud_props(i)%name = cldtypes(i)

       call read_cmaq_file(cloud_props(i)%num,   cloud_props(i)%start, cloud_props(i)%end,   &
                           cloud_props(i)%delr,  cloud_props(i)%extsw, cloud_props(i)%ssasw, &
                           cloud_props(i)%asysw, cloud_props(i)%name                         )
    enddo

    call shdf5_close()

  contains

    subroutine read_cmaq_file(num, start, end, delr, extsw, ssasw, asysw, name)
      use csqy_data,  only: nwl
      use hdf5_utils, only: shdf5_irec
      implicit none

      integer,           intent(out)   :: num
      real,              intent(out)   :: start, end, delr
      real, allocatable, intent(inout) :: extsw(:,:), ssasw(:,:), asysw(:,:)
      character(len=*),  intent(in)    :: name
      real, allocatable                :: t1(:,:)

      call shdf5_irec( 1, (/1/), trim(name) // '_numr'  , ivars=num   )
      call shdf5_irec( 1, (/1/), trim(name) // '_rstart', rvars=start )
      call shdf5_irec( 1, (/1/), trim(name) // '_rend'  , rvars=end   )

      delr = (end - start) / real(num - 1)

      allocate(extsw(nwl, num))
      allocate(ssasw(nwl, num))
      allocate(asysw(nwl, num))

      allocate(t1(num, nwl))

      call shdf5_irec( 2, (/num, nwl/), trim(name) // '_ext', rvar2=t1 )
      extsw = transpose(t1)

      call shdf5_irec( 2, (/num, nwl/), trim(name) // '_ssa', rvar2=t1 )
      ssasw = transpose(t1)

      call shdf5_irec( 2, (/num, nwl/), trim(name) // '_asy', rvar2=t1 )
      asysw = transpose(t1)

      deallocate(t1)

    end subroutine read_cmaq_file

  end subroutine cmaq_cld_optics_init

end module cmaq_cloud
