Module ll_bins

  implicit none

  integer, parameter :: nlat(4) = (/ 18,   90,  450, 1800 /) ! Number of bins spanning pole-to-pole
                                                             ! distance for each bin set

  integer, parameter :: maxb(4) = (/ 200,  200,  200, 10000/)  ! Max allowed population in a bin
                                                               ! for each bin set

  real, parameter :: delat(4) = (/ real(nlat(:)) / 180. /) ! # bins per deg of lat for each bin set

  Type bin_vars
     integer              :: nw
     integer, allocatable :: iw(:)
  End type bin_vars

  Type blat_vars
     integer                     :: nlon    ! # bins per lon circle for each bin set and lat band
     real                        :: delon   ! # bins per deg of lon for each bin set and lat band
     type(bin_vars), allocatable :: bins(:)
  End type

  Type binset_vars
     type(blat_vars), allocatable :: blat(:)
  End type

  type(binset_vars), allocatable :: bset(:)

  ! itab_w0 stores a copy of the OLD ATM grid structure for a HISTREGRID runtype,
  ! and it stores the normal ATM grid structure for computing overlay with SFCGRID.

  Type itab_w0_vars
     integer :: npoly = 0 ! number of W neighbors of this W pt
     integer :: iw(7)     ! neighbor W pts
  End Type itab_w0_vars

  type(itab_w0_vars), allocatable :: itab_w0(:)

Contains

!=========================================================================

  subroutine latlon_bins(nwa, glatw, glonw)

  implicit none

  integer, intent(in) :: nwa
  real,    intent(in) :: glatw(nwa), glonw(nwa)

  integer :: iset, jsup, i, j, ipass, nlato3, nlon
  integer :: jw, iw, jwn, iwn, npoly, niwtemp

  real :: hlatw, hlonw

  integer, allocatable :: iwtemp(:), iwtemp2(:)

  allocate(bset(4))

  do iset = 1,4
     allocate(bset(iset)%blat(nlat(iset)))
     nlato3 = nlat(iset) / 3

     do j = 1, nlat(iset)
        jsup = nlat(iset) + 1 - j

        nlon = 6 * min(j, jsup, nlato3)

        bset(iset)%blat(j)%nlon  = nlon
        bset(iset)%blat(j)%delon = real(nlon) / 360.

        allocate(bset(iset)%blat(j)%bins(nlon)); bset(iset)%blat(j)%bins(:)%nw = 0
     enddo
  enddo

  ! Loop through all ATM iw cells and count them into 4 sets of global lat-lon
  ! bins spanning a range of spatial resolution

  do iw = 2, nwa
     hlonw = max(-179.9999,min(179.9999,glonw(iw)))
     hlatw = max( -89.9999,min( 89.9999,glatw(iw)))

     do iset = 1,4
        j = int(delat(iset) * (hlatw +  90.)) + 1
        i = int(bset(iset)%blat(j)%delon * (hlonw + 180.)) + 1

        bset(iset)%blat(j)%bins(i)%nw = bset(iset)%blat(j)%bins(i)%nw + 1
     enddo
  enddo

  ! Loop through all bins of all sets and allocate IW(:) array if NW is
  ! greater than 0 and less than maxb

  do iset = 1,4
     do j = 1, nlat(iset)
        do i = 1, bset(iset)%blat(j)%nlon
           if (bset(iset)%blat(j)%bins(i)%nw > 0 .and. bset(iset)%blat(j)%bins(i)%nw < maxb(iset)) then
              allocate(bset(iset)%blat(j)%bins(i)%iw( bset(iset)%blat(j)%bins(i)%nw ))
           endif
           bset(iset)%blat(j)%bins(i)%nw = 0
        enddo
     enddo
  enddo

  ! Loop through all ATM iw cells and sort them into 4 sets of global lat-lon
  ! bins if bin is allocated

  do iw = 2, nwa
     hlonw = max(-179.9999,min(179.9999,glonw(iw)))
     hlatw = max( -89.9999,min( 89.9999,glatw(iw)))

     do iset = 1,4
        j = int(delat(iset) * (hlatw +  90.)) + 1
        i = int(bset(iset)%blat(j)%delon * (hlonw + 180.)) + 1

        if (allocated(bset(iset)%blat(j)%bins(i)%iw)) then
           bset(iset)%blat(j)%bins(i)%nw = bset(iset)%blat(j)%bins(i)%nw + 1
           bset(iset)%blat(j)%bins(i)%iw(  bset(iset)%blat(j)%bins(i)%nw ) = iw
        endif
     enddo
  enddo

  ! Loop TWICE through all allocated bins and augment their IW membership with
  ! nearest IW neighbors that were not already in bin.  This adds two perimeter
  ! rows of IW cells to each bin.

  do ipass = 1,2
     do iset = 1,4
        do j = 1, nlat(iset)
           do i = 1, bset(iset)%blat(j)%nlon

              if (allocated(bset(iset)%blat(j)%bins(i)%iw)) then

                 allocate (iwtemp(maxb(iset)+1000))
                 niwtemp = 0

                 do jw = 1, bset(iset)%blat(j)%bins(i)%nw
                    iw = bset(iset)%blat(j)%bins(i)%iw(jw)

                    npoly = itab_w0(iw)%npoly

                    do jwn = 1, npoly
                       iwn = itab_w0(iw)%iw(jwn)

                       if (any(iwn == bset(iset)%blat(j)%bins(i)%iw(:))) cycle
                       if (niwtemp > 0 .and. any(iwn == iwtemp(1:niwtemp))) cycle

                       niwtemp = niwtemp + 1
                       iwtemp(niwtemp) = iwn

                    enddo
                 enddo

                 allocate (iwtemp2(bset(iset)%blat(j)%bins(i)%nw + niwtemp))

                 iwtemp2(1:bset(iset)%blat(j)%bins(i)%nw) &
                         = bset(iset)%blat(j)%bins(i)%iw( 1:bset(iset)%blat(j)%bins(i)%nw )

                 iwtemp2(bset(iset)%blat(j)%bins(i)%nw+1:bset(iset)%blat(j)%bins(i)%nw + niwtemp) &
                         = iwtemp(1:niwtemp)

                 deallocate (iwtemp)

                 call move_alloc(iwtemp2, bset(iset)%blat(j)%bins(i)%iw)

                 bset(iset)%blat(j)%bins(i)%nw = bset(iset)%blat(j)%bins(i)%nw + niwtemp

              endif
           enddo    ! i
        enddo    ! j
     enddo    ! iset
  enddo    ! ipass

  end subroutine latlon_bins

end module ll_bins

