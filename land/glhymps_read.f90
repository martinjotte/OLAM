subroutine glhymps_read()

  use mem_land,    only: land, nland, onland
  use mem_sfcg,    only: sfcg
  use consts_coms, only: piu180
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use leaf_coms,   only: glhymps_database
  use misc_coms,   only: io6

  implicit none

  real, allocatable :: dato(:,:)

  integer :: nio, njo
  integer :: io, jo
  integer :: idataset, iland, iwsfc
  integer :: ndims, idims(2)

  real :: xoffpix, yoffpix
  real :: glat, glon
  real :: xperdeg, yperdeg

  character(pathlen) :: fname, fullname

  logical :: exists

  ! This subroutine is simpler than land_database_read because it assumes that
  ! each glhymps database file covers the entire geographic area of the model.
  ! If this ever changes, this subroutine must be modified.

  ! Loop over the three glhymps datasets

  do idataset = 1,3

     if     (idataset == 1) then
        fname = 'permeability_no_permafrost_5km.h5'
     elseif (idataset == 2) then
        fname = 'permeability_permafrost_5km.h5'
        cycle  !! skip reading for now; not used
     elseif (idataset == 3) then
        fname = 'porosity_5km.h5'
     endif

     ! Make full filename

     fullname = trim(glhymps_database)//trim(fname)

     ! Open the sst filelist and read through the files

     write(io6,*) 'glhymps_read reading ', trim(fullname)

     inquire(file=fullname, exist=exists)
     if (.not. exists) then
        write(io6,*) 'glhymps: error opening glhymps file ' // trim(fullname)
        stop
     endif

     ! Open and read glhymps_database file

     call shdf5_open(fullname,'R',trypario=.true.)
     call shdf5_info('Band1',ndims,idims)

     nio = idims(1)
     njo = idims(2)
     allocate(dato(nio,njo))

     call shdf5_irec(ndims,idims,'Band1',rvar2=dato)
     call shdf5_close()

     ! glhymps data files do NOT have an added perimeter, and data values are
     ! offset by 1/2 pixel from integer lat-lon values.

     xoffpix = 0.5
     yoffpix = 0.5

     xperdeg = real(nio) / 360.0
     yperdeg = real(njo) / 180.0

     ! Fill olam-soil array from nearest neighbor in dataset field

     do iland = 2, nland
        iwsfc = iland + onland

        glat = sfcg%glatw(iwsfc)
        glon = sfcg%glonw(iwsfc)

        glon = max(-179.999,min(179.999,glon))

        io = nint((glon + 180.) * xperdeg + xoffpix)
        jo = nint((glat +  90.) * yperdeg + yoffpix)

        ! Exclude "undefined" and "missing data values.
        ! For idataset = 1 and 2, dato = log10(permeability). Convert this to ksat;
        ! ksat is 1.e7 * permeability

        if     (idataset == 1) then
           if (dato(io,jo) >= -20. .and. dato(io,jo) < -10.) then
              land%glhymps_ksat(iland) = 10.0**(dato(io,jo) + 7.0)
           else
              land%glhymps_ksat(iland) = 10.0**(-20.0 + 7.0) ! Require a nonzero minimum value
           endif
!!      elseif (idataset == 2) then
!!         if (dato(io,jo) >= -20. .and. dato(io,jo) < -10.) then
!!            land%glhymps_ksat_pfr(iland) = 10.0**(dato(io,jo) + 7.0)
!!         else
!!            land%glhymps_ksat_pfr(iland) = 10.0**(-20.0 + 7.0) ! Require a nonzero minimum value
!!         endif
        elseif (idataset == 3) then
           if (dato(io,jo) >= 0.01 .and. dato(io,jo) < 0.5) then
              land%glhymps_poros(iland) = dato(io,jo)
           else
              land%glhymps_poros(iland) = 0.01 ! Require a nonzero minimum value
           endif
        endif
     enddo

     deallocate(dato)

  enddo ! idataset

end subroutine glhymps_read
