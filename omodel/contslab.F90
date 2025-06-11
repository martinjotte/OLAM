subroutine contslab_horiz_mp(iplt)

  use oplot_coms,   only: op
  use mem_grid,     only: mma, mwa, xem, yem, zem, xew, yew, zew, &
                          glonm, glatm, glonw, glatw
  use mem_ijtabs,   only: itab_w, jtab_w, jtw_prog, jtab_m, jtm_prog
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank, mgroupsize, nbytes_int, nbytes_real
  use olam_mpi_atm, only: mpi_recv_m, mpi_send_m

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: npoly,j,iw,jw,im,jm
  integer :: iflag180
  integer :: ipwx1,ipwx2,ipwy1,ipwy2

  real :: hpt,vpt
  real :: hcpn(7),vcpn(7),fldvals(7)

  integer :: ktf(mwa),km(mma)
  real :: wtbot(mma),wttop(mma)

  integer :: notavail(mma)
  real    :: opltvals(mma)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 21 * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mma) / 5. )
  else
     inc = mma
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! For parallel runs, we need to compute plot values only on primary points
  ! and then communicate to the border cells

  !$omp parallel do private(im)
  do jm = 1, jtab_m(jtm_prog)%jend
     im = jtab_m(jtm_prog)%im(jm)

     call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                    opltvals(im),notavail(im))

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_m(r1dvara1=opltvals, i1dvara1=notavail)
     call mpi_recv_m(r1dvara1=opltvals, i1dvara1=notavail)
  endif

  ! Loop over W points for contouring M points

  wloop: do jw = 1, jtab_w(jtw_prog)%jend
            iw = jtab_w(jtw_prog)%iw(jw)

     ! Skip this W point if it is underground

     if (ktf(iw) /= 0) cycle

     ! Get plot coordinates of current W point.

     call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),glonw(iw),glatw(iw),hpt,vpt)

     npoly = itab_w(iw)%npoly

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     ! Loop over all M points that surround current W point

     do j = 1, npoly

        ! Current M point index

        im = itab_w(iw)%im(j)

        ! Skip current M point if index < 2

        if (im < 2) cycle wloop

        ! Get plot coordinates of current M point

        call oplot_transform(iplt,xem(im),yem(im),zem(im),glonm(im),glatm(im),hcpn(j),vcpn(j))

        ! Skip this W point if current M point is far outside plot window
        ! (which means that orthographic projection returned large value
        ! that indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle wloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any M point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

        fldvals(j) = opltvals(im)

        if (notavail(im) > 0) cycle wloop
     enddo

     ! If any window flag is zero, all M points for this W point are outside
     ! the same window boundary, so skip this W point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle wloop

     ! Contour plot cell of 2-D or 3-D field

     if (myrank == 0) then

        call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
     ! at other end

     if (iflag180 /= 0) then

        do j = 1,npoly
           hcpn(j) = hcpn(j) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

        else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif

  enddo wloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, npoly,   1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, npoly, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

! NOTE: For Hex grid, may in the future want underground to always be
!       triangle areas.  In that case, first do underground plot and then
!       call contpolyg selectively for trios of points in remaining
!       (open) sectors, using U/V point at face of underground block.

end subroutine contslab_horiz_mp

!===============================================================================

subroutine contslab_topmw(iplt)

! Special routine to plot topm and topw values together

  use oplot_coms, only: op
  use mem_grid,   only: mma, xem, yem, zem, xew, yew, zew, topm, topw, &
                        glonm, glatm, glonw, glatw
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use misc_coms,  only: iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: npoly,j,iw,jw,im,jn
  integer :: iflag180
  integer :: ipwx1,ipwx2,ipwy1,ipwy2

  real :: hpt,vpt
  real :: hcpn(7),vcpn(7),fldvals(7)
  real :: hcpn3(3),vcpn3(3),fldvals3(3)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  nu   = 0
  ipos = 0

  base = 9 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mma) / 5. )
  else
     inc = mma
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over W points

  wloop: do jw = 1, jtab_w(jtw_prog)%jend
            iw = jtab_w(jtw_prog)%iw(jw)

     ! Get plot coordinates of current W point.

     call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),glonw(iw),glatw(iw),hpt,vpt)

     npoly = itab_w(iw)%npoly

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     ! Loop over all M points that surround current W point

     do j = 1, npoly

        ! Current M point index

        im = itab_w(iw)%im(j)

        ! Skip current M point if index < 2

        if (im < 2) cycle wloop

        ! Get plot coordinates of current M point

        call oplot_transform(iplt,xem(im),yem(im),zem(im),glonm(im),glatm(im),hcpn(j),vcpn(j))

        ! Skip this W point if current M point is far outside plot window
        ! (which means that orthographic projection returned large value
        ! that indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle wloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any M point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

     enddo

     ! If any window flag is zero, all M points for this W point are outside
     ! the same window boundary, so skip this W point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle wloop

     ! Loop over all M points that surround current W point and fill field values

     do j = 1, npoly
        im = itab_w(iw)%im(j)
        fldvals(j) = topm(im)
     enddo

     do j = 1, npoly
        jn = j + 1
        if (jn > npoly) jn = 1

        hcpn3(1) = hpt
        vcpn3(1) = vpt
        fldvals3(1) = topw(iw)

        hcpn3(2) = hcpn(j)
        vcpn3(2) = vcpn(j)
        fldvals3(2) = fldvals(j)

        hcpn3(3) = hcpn(jn)
        vcpn3(3) = vcpn(jn)
        fldvals3(3) = fldvals(jn)

        if (myrank == 0) then

           call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(hcpn3,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vcpn3,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldvals3, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
        ! at other end

        if (iflag180 /= 0) then

           hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(hcpn3,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn3,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals3, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif

     enddo ! j/m loop

  enddo wloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn3,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn3,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals3, npoly, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

end subroutine contslab_topmw

!===============================================================================

subroutine contslab_horiz_vn(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mma, mva, mwa, xem, yem, zem, &
                        xev, yev, zev, xew, yew, zew, &
                        glonm, glatm, glonv, glatv, glonw, glatw
  use mem_ijtabs, only: itab_m, itab_w, itab_v, jtab_m, jtm_prog, &
                        jtab_w, jtw_prog
  use misc_coms,  only: iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: npoly,j,jm,jn,jnn,im,notavail,iw,iv,iw1,iw2,iv1,iv2,jw
  integer :: iflag180
  integer :: ipwx1,ipwx2,ipwy1,ipwy2

  real :: hpt,vpt
  real :: hcpn(7),vcpn(7),fldvals(7)
  real :: hcpn3(3),vcpn3(3),fldvals3(3)
  real :: avail, avg

  integer :: ktf(mwa),kv(mva)
  real :: wtbot(mva),wttop(mva)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 9 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mma) / 5. )
  else
     inc = mma
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

!------------------------------------------------------
! FIRST LOOP is over M points for contouring V points
!------------------------------------------------------

  mloop: do jm = 1, jtab_m(jtm_prog)%jend
            im = jtab_m(jtm_prog)%im(jm)

     ! Get plot coordinates of current M point.

     call oplot_transform(iplt,xem(im),yem(im),zem(im),glonm(im),glatm(im),hpt,vpt)

     npoly = itab_m(im)%npoly

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     avail = 0.
     avg = 0.

     ! Loop over all U or V points that surround current M point

     do j = 1, npoly

        ! Current U/V point index

        iv  = itab_m(im)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! Skip current V cell if index < 2

        if (iv < 2) cycle mloop

        ! Get plot coordinates of current V point

        call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),glonv(iv),glatv(iv),hcpn(j),vcpn(j))

        ! Skip this M point if current V point is far outside plot window
        ! (which means that orthographic projection returned large value that
        ! indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle  mloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any V point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

        if (ktf(iw1) /= 0 .and. ktf(iw2) /= 0) cycle

        call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                       fldvals(j),notavail)

        if (ktf(iw1) == 0) then
           avail = avail + .5
           avg = avg + .5 * fldvals(j)
        endif

        if (ktf(iw2) == 0) then
           avail = avail + .5
           avg = avg + .5 * fldvals(j)
        endif

     enddo

     if (avail < .1) cycle mloop

     avg = avg / avail

     ! If any window flag is zero, all V points for this M point are outside
     ! the same window boundary, so skip this M point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle mloop

     ! Loop over all V points that surround current M point and contour plot
     ! each available sector

     do j = 1, npoly

        jn = j + 1
        if (jn > npoly) jn = 1
        jnn = jn + 1
        if (jnn > npoly) jnn = 1

        iv1 = itab_m(im)%iv(j)
        iv2 = itab_m(im)%iv(jn)

        ! Specific way to get IW since ordering of W and U/V neighbors of M
        ! is not identical for both grid systems

        iw = itab_m(im)%iw(jnn)

        if (ktf(iw) == 0) then

           hcpn3(1) = hpt
           vcpn3(1) = vpt
           fldvals3(1) = avg

           hcpn3(2) = hcpn(j)
           vcpn3(2) = vcpn(j)
           fldvals3(2) = fldvals(j)

           hcpn3(3) = hcpn(jn)
           vcpn3(3) = vcpn(jn)
           fldvals3(3) = fldvals(jn)

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(hcpn3,    npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn3,    npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals3, npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

           ! If lat/lon plot and this polygon crosses +/- 180 degrees longitude,
           ! plot again at other end of plot

           if (iflag180 /= 0) then

              hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180

              if (myrank == 0) then

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

              else

#ifdef OLAM_MPI
                 nu = nu + 1
                 if (buffsize < ipos + base) then
                    allocate( bcopy (buffsize + inc * base) )
                    bcopy(1:buffsize) = buffer
                    call move_alloc(bcopy, buffer)
                    buffsize = size(buffer)
                 endif

                 call MPI_Pack(hcpn3,    npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vcpn3,    npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(fldvals3, npoly, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

              endif

           endif

        endif

     enddo

  enddo mloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn3,    npoly, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn3,    npoly, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals3, npoly, MPI_REAL, MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  nu   = 0
  ipos = 0

  base = 21 * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mwa) / 5. )
  else
     inc = mwa
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

!------------------------------------------------------
! SECOND LOOP is over W points for contouring V points
!------------------------------------------------------

  wloop: do jw = 1, jtab_w(jtw_prog)%jend
            iw = jtab_w(jtw_prog)%iw(jw)

     if (ktf(iw) /= 0) cycle

     ! Get plot coordinates of current W point.

     call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),glonw(iw),glatw(iw),hpt,vpt)

     npoly = itab_w(iw)%npoly

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     ! Loop over all V points that surround current W point

     do j = 1, npoly

        ! Current V point index

        iv = itab_w(iw)%iv(j)

        ! Skip current V cell if index < 2

        if (iv < 2) cycle wloop

        ! Get plot coordinates of current V point

        call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),glonv(iv),glatv(iv),hcpn(j),vcpn(j))

        ! Skip this M point if current V point is far outside plot window
        ! (which means that orthographic projection returned large value that
        ! indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle wloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any M point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

     enddo

     ! If any window flag is zero, all M points for this W point are outside
     ! the same window boundary, so skip this W point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle wloop

     ! Loop over all V points that surround current W point and fill field values

     do j = 1,npoly
        iv = itab_w(iw)%iv(j)

        call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                       fldvals(j),notavail)
        if (notavail > 0) cycle
     enddo

     ! Contour plot cell of 2-D or 3-D field

     if (myrank == 0) then

        call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot,
     ! re-plot at other end

     if (iflag180 /= 0) then

        do j = 1, npoly
           hcpn(j) = hcpn(j) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif

  enddo wloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, npoly,   1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, npoly, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

end subroutine contslab_horiz_vn

!===============================================================================

subroutine contslab_horiz_tw(iplt)

  use oplot_coms,   only: op
  use mem_grid,     only: mma, mwa, xem, yem, zem, &
                          xev, yev, zev, xew, yew, zew, &
                          glonm, glatm, glonv, glatv, glonw, glatw
  use mem_ijtabs,   only: itab_m, jtab_m, jtm_prog, jtab_w, jtw_prog
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank, mgroupsize, nbytes_real
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: npoly,j,jm,jn,jnn,im,iv,jw,iw,iw1,iw2
  integer :: iflag180,avail
  integer :: ipwx1,ipwx2,ipwy1,ipwy2

  real :: hpt,vpt
  real :: hcpn (3),vcpn (3),fldvals (3)
  real :: hcpn3(3),vcpn3(3),fldvals3(3)
  real :: avg

  integer :: ktf(mwa), kw(mwa)
  real    :: wtbot(mwa), wttop(mwa)
  real    :: opltvals(mwa)
  integer :: notavail(mwa)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 9 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mma) / 5. )
  else
     inc = mma
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! For parallel runs, we need to compute plot values only on primary points
  ! and then communicate to the border cells

  !$omp parallel do private(iw)
  do jw = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(jw)

     if (ktf(iw) /= 0) then
        notavail(iw) = 3
     else
        call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                       opltvals(iw),notavail(iw))
     endif

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_w(r1dvara1=opltvals, i1dvara1=notavail)
     call mpi_recv_w(r1dvara1=opltvals, i1dvara1=notavail)
  endif

  ! Loop over primary M points for contouring W points

  mloop: do jm = 1, jtab_m(jtm_prog)%jend
            im = jtab_m(jtm_prog)%im(jm)

     ! Get plot coordinates of current M point.

     call oplot_transform(iplt,xem(im),yem(im),zem(im),glonm(im),glatm(im),hpt,vpt)

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     avail = 0
     avg   = 0.

     npoly = itab_m(im)%npoly

     ! Loop over all W points that surround current M point

     do j = 1, npoly

        ! Current W point index

        iw = itab_m(im)%iw(j)

        ! Get plot coordinates of current W point

        call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),glonw(iw),glatw(iw),hcpn(j),vcpn(j))

        ! Skip this M point if current W point is far outside plot window
        ! (which means that orthographic projection returned large value that
        ! indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle mloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any W point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

        if (notavail(iw) > 0) cycle

        fldvals(j) = opltvals(iw)

        avail = avail + 1
        avg = avg + fldvals(j)

     enddo

     if (avail == 0) cycle mloop

     avg = avg / real(avail)

     ! If any window flag is zero, all W points for this M point are outside
     ! the same window boundary, so skip this W point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle mloop

     ! If all W points around this M point are available, plot them together

     if (avail == 3) then

        if (myrank == 0) then

           call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(hcpn,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vcpn,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldvals, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
        ! again at other end of plot

        if (iflag180 /= 0) then

           hcpn(1:3) = hcpn(1:3) + 360. * iflag180

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(hcpn,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
           endif

        endif

     ! If only some W points are available, plot them by sectors

     else

        ! Loop over all W points that surround current M point and contour plot each
        ! available sector

        do j = 1, npoly
           jn = j + 1
           if (jn > npoly) jn = 1
           jnn = jn + 1
           if (jnn > npoly) jnn = 1

           iw1 = itab_m(im)%iw(j)
           iw2 = itab_m(im)%iw(jn)

           hcpn3(1) = hpt
           vcpn3(1) = vpt
           fldvals3(1) = avg

           if (notavail(iw1) == 0 .and. notavail(iw2) == 0) then

              hcpn3(2) = hcpn(j)
              vcpn3(2) = vcpn(j)
              fldvals3(2) = fldvals(j)

              hcpn3(3) = hcpn(jn)
              vcpn3(3) = vcpn(jn)
              fldvals3(3) = fldvals(jn)

           elseif (notavail(iw1) == 0) then

              ! Specific way to get IV since ordering of W and U/V neighbors of M is not
              ! identical for both grid systems

              iv = itab_m(im)%iv(jnn)
              call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),glonv(iv),glatv(iv),hcpn3(3),vcpn3(3))

              fldvals3(3) = fldvals(j)

              hcpn3(2) = hcpn(j)
              vcpn3(2) = vcpn(j)
              fldvals3(2) = fldvals(j)

           elseif (notavail(iw2) == 0) then

              ! Specific way to get IV since ordering of W and U/V neighbors of M is not
              ! identical for both grid systems

              iv = itab_m(im)%iv(jnn)
              call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),glonv(iv),glatv(iv),hcpn3(2),vcpn3(2))

              fldvals3(2) = fldvals(jn)

              hcpn3(3) = hcpn(jn)
              vcpn3(3) = vcpn(jn)
              fldvals3(3) = fldvals(jn)

           else

              cycle

           endif

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(hcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals3, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

           ! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
           ! again at other end of plot

           if (iflag180 /= 0) then

              hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180

              if (myrank == 0) then

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

              else

#ifdef OLAM_MPI
                 nu = nu + 1
                 if (buffsize < ipos + base) then
                    allocate( bcopy (buffsize + inc * base) )
                    bcopy(1:buffsize) = buffer
                    call move_alloc(bcopy, buffer)
                    buffsize = size(buffer)
                 endif

                 call MPI_Pack(hcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(fldvals3, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
              endif

           endif
        enddo

     endif

  enddo mloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    3, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    3, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, 3, MPI_REAL, MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

end subroutine contslab_horiz_tw

!===============================================================================

subroutine contslab_horiz_sfc(iplt)

  use oplot_coms,   only: op
  use mem_sfcg,     only: mwsfc, mmsfc, sfcg, itab_msfc, itab_wsfc
  use mem_land,     only: nzg
  use leaf_coms,    only: nzs
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank, mgroupsize, nbytes_real
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: k,j,imsfc,iwn,jn,jnn,iv,iw1,iw2
  integer :: iflag180
  integer :: ipwx1,ipwx2,ipwy1,ipwy2

  real    :: opltvals(mwsfc)
  integer :: notavail(mwsfc)

  real :: hpt,vpt
  real :: hcpn (3),vcpn (3),fldvals (3)
  real :: hcpn3(3),vcpn3(3),fldvals3(3)
  real :: avg

  integer :: avail

  real, parameter :: wtbot = 1., wttop = 0.

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer, parameter :: npoly = 3

  nu   = 0
  ipos = 0

  base = 9 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mmsfc) / 5. )
  else
     inc = mmsfc
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Find K level to plot if field is 3d

  if (op%dimens == '3G') then
     k = min(nzg,max(1,nint(op%slabloc(iplt))))
  elseif (op%dimens == '3S') then
     k = min(nzs,max(1,nint(op%slabloc(iplt))))
  else
     k = 1
  endif

  ! Compute plot values only on primary points for a parallel run

  !$omp parallel do
  do iwn = 2, mwsfc

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (iparallel == 1 .and. itab_wsfc(iwn)%irank /= myrank) cycle

     ! Skip this cell if wrong surface type

     if (op%stagpt == 'L' .and. sfcg%leaf_class(iwn) <  2) then
        notavail(iwn) = 3
        cycle
     endif

     if (op%stagpt == 'R' .and. sfcg%leaf_class(iwn) /= 1) then
        notavail(iwn) = 3
        cycle
     endif

     if (op%stagpt == 'S' .and. sfcg%leaf_class(iwn) /= 0) then
        notavail(iwn) = 3
        cycle
     endif

     call oplot_lib( k, iwn, 'VALUE', op%fldname(iplt), wtbot, wttop, &
                     opltvals(iwn), notavail(iwn) )
  enddo
  !$omp end paralell do

  ! Communicate plot values to border cells

  if (iparallel == 1) then
     call mpi_send_wsfc('plot_avg', rvar=opltvals, ivar=notavail)
     call mpi_recv_wsfc('plot_avg', rvar=opltvals, ivar=notavail)
  endif

  ! Loop over M points for contouring W points

  mloop: do imsfc = 2, mmsfc

     ! Skip border cells if running in parallel

     if (iparallel == 1 .and. itab_msfc(imsfc)%irank /= myrank) cycle

     ! Get plot coordinates of current MSFC point.

     call oplot_transform(iplt,sfcg%xem(imsfc),sfcg%yem(imsfc),sfcg%zem(imsfc), &
                          sfcg%glonm(imsfc),sfcg%glatm(imsfc),hpt,vpt)

     ! Initialize iflag180 and plot window flags to zero

     iflag180 = 0

     ipwx1 = 0
     ipwx2 = 0
     ipwy1 = 0
     ipwy2 = 0

     avail = 0
     avg   = 0.

     ! Loop over all WSFC points that surround current MSFC point

     do j = 1, npoly

        ! Current W point index

        iwn = itab_msfc(imsfc)%iwn(j)

        ! Get plot coordinates of current W point

        call oplot_transform(iplt,sfcg%xew(iwn),sfcg%yew(iwn),sfcg%zew(iwn), &
                             sfcg%glonw(iwn),sfcg%glatw(iwn),hcpn(j),vcpn(j))

        ! Skip this MSFC point if current IWN point is far outside plot window
        ! (which means that orthographic projection returned large value that
        ! indicates that point is on other side of Earth)

        if (abs(hcpn(j)) > 1.e11) cycle mloop

        ! Avoid wrap-around and set iflag180

        if (op%projectn(iplt)== 'L') then
           call ll_unwrap(hpt,hcpn(j))
           if (hcpn(j) < -180.001) iflag180 =  1
           if (hcpn(j) >  180.001) iflag180 = -1
        endif

        ! Set plot window flag to 1 if any W point is on window side of
        ! respective boundary

        if (hcpn(j) >= op%xmin) ipwx1 = 1
        if (hcpn(j) <= op%xmax) ipwx2 = 1
        if (vcpn(j) >= op%ymin) ipwy1 = 1
        if (vcpn(j) <= op%ymax) ipwy2 = 1

        ! Skip land/lake/sea cells

        if ( notavail(iwn) > 0) cycle

        fldvals(j) = opltvals(iwn)

        avail = avail + 1
        avg = avg + fldvals(j)

     enddo

     if (avail == 0) cycle

     avg = avg / real(avail)

     ! If any window flag is zero, all W points for this M point are outside
     ! the same window boundary, so skip this W point

     if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) cycle mloop

     ! If all W points around this M point are available, plot them together

     if (avail == 3) then

        if (myrank == 0) then

           call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

!          call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
        ! again at other end of plot

        if (iflag180 /= 0) then

           hcpn(1:npoly) = hcpn(1:npoly) + 360. * iflag180

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

!             call MPI_Pack(npoly,   1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(hcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn,    npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals, npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif

     else

        ! Loop over all W points that surround current M point and contour plot each
        ! available sector

        do j = 1, npoly
           jn = j + 1
           if (jn > npoly) jn = 1
           jnn = jn + 1
           if (jnn > npoly) jnn = 1

           iw1 = itab_msfc(imsfc)%iwn(j)
           iw2 = itab_msfc(imsfc)%iwn(jn)

           hcpn3(1) = hpt
           vcpn3(1) = vpt
           fldvals3(1) = avg

           if (notavail(iw1) == 0 .and. notavail(iw2) == 0) then

              hcpn3(2) = hcpn(j)
              vcpn3(2) = vcpn(j)
              fldvals3(2) = fldvals(j)

              hcpn3(3) = hcpn(jn)
              vcpn3(3) = vcpn(jn)
              fldvals3(3) = fldvals(jn)

           elseif (notavail(iw1) == 0) then

              ! Specific way to get IV since ordering of W and U/V neighbors of M is not
              ! identical for both grid systems

              iv = itab_msfc(imsfc)%ivn(jnn)
              call oplot_transform_xyz(iplt,sfcg%xev(iv),sfcg%yev(iv),sfcg%zev(iv),hcpn3(3),vcpn3(3))

              fldvals3(3) = fldvals(j)

              hcpn3(2) = hcpn(j)
              vcpn3(2) = vcpn(j)
              fldvals3(2) = fldvals(j)

           elseif (notavail(iw2) == 0) then

              ! Specific way to get IV since ordering of W and U/V neighbors of M is not
              ! identical for both grid systems

              iv = itab_msfc(imsfc)%ivn(jnn)
              call oplot_transform_xyz(iplt,sfcg%xev(iv),sfcg%yev(iv),sfcg%zev(iv),hcpn3(2),vcpn3(2))

              fldvals3(2) = fldvals(jn)

              hcpn3(3) = hcpn(jn)
              vcpn3(3) = vcpn(jn)
              fldvals3(3) = fldvals(jn)

           else

              cycle

           endif

           if (myrank == 0) then

              call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(hcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvals3, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

           ! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
           ! again at other end of plot

           if (iflag180 /= 0) then

              hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180

              if (myrank == 0) then

                 call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

              else

#ifdef OLAM_MPI
                 nu = nu + 1
                 if (buffsize < ipos + base) then
                    allocate( bcopy (buffsize + inc * base) )
                    bcopy(1:buffsize) = buffer
                    call move_alloc(bcopy, buffer)
                    buffsize = size(buffer)
                 endif

                 call MPI_Pack(hcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vcpn3,    3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(fldvals3, 3, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
              endif

           endif
        enddo

     endif

  enddo mloop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

!                call MPI_Unpack(buffer, buffsize, ipos, npoly,   1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, npoly, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

end subroutine contslab_horiz_sfc

!===============================================================================

subroutine contslab_vert_v(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mza, lpw, zt
  use misc_coms,  only: iparallel
  use mem_ijtabs, only: jtab_w, jtw_prog
  use consts_coms,only: pio180
  use mem_para,   only: myrank, mgroupsize, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: k,jw,iw,iv1,iv2,iok,notavail
  real :: hpt
  real :: hcpn(4),vcpn(4),fldvals(4)
  real :: topo1,topo2

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  real, parameter :: wtbot = 1., wttop = 0.

  nu   = 0
  ipos = 0

  base = 12 * nbytes_real
  inc = ceiling( real(mwa) / 10. ) * mza

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  do jw = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(jw)

     ! Get horizontal plot coordinates for IW point

     if (op%projectn(iplt) == 'C') then
        call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,hcpn)
     elseif (op%projectn(iplt)(1:1) == 'V') then
        call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,hcpn) ! Need to fix for hex??
     endif

     ! Skip current IW point if it does not intersect plot cone

     if (iok /= 1) cycle

     ! Avoid wrap_around

     if (op%projectn(iplt) == 'C') then
        call ll_unwrap(hcpn(1),hcpn(2))
     endif

     hpt = .5 * (hcpn(1) + hcpn(2))

     ! Skip current IV point if entire cell is outside plot window.

     if ( (hcpn(1) < op%xmin .and. hcpn(2) < op%xmin) .or.  &
          (hcpn(1) > op%xmax .and. hcpn(2) > op%xmax) ) cycle

     hcpn(3) = hcpn(2)
     hcpn(4) = hcpn(1)

     do k = lpw(iw),mza-1   ! Loop is over W levels

        ! Skip this K point if entire cell is above or below plot window

        if (zt(k+1) < op%ymin .or. zt(k) > op%ymax) cycle

        ! Get T-cell vertical coordinates

        vcpn(1) = zt(k)
        vcpn(2) = vcpn(1)
        vcpn(3) = zt(k+1)
        vcpn(4) = vcpn(3)

        ! Fill field values of 4 T points around current M point

        call oplot_lib(k,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                       fldvals(1),notavail)
        if (notavail > 0) cycle
        call oplot_lib(k,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                       fldvals(2),notavail)
        if (notavail > 0) cycle
        call oplot_lib(k+1,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                       fldvals(3),notavail)
        if (notavail > 0) cycle
        call oplot_lib(k+1,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                       fldvals(4),notavail)
        if (notavail > 0) cycle

        if (myrank == 0) then

           ! Contour plot cell around current M point
           call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(hcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldvals, 4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
        endif

     enddo

  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, 4, MPI_REAL, MPI_COMM_WORLD, ier)

                 ! Contour plot cell around current M point
                 call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  ! Now plot underground T cells with underground color
  call plot_underground_w(iplt,(/0/))

end subroutine contslab_vert_v

!===============================================================================

subroutine contslab_vert_t(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mza, lpw, zm, zt, &
                        xem, yem, zem, xev, yev, zev, xew, yew, zew
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use misc_coms,  only: mdomain, iparallel
  use consts_coms,only: erad, pio180
  use mem_para,   only: myrank, mgroupsize, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: ilev
  integer :: k,iw,iw1,iw2,iw3,iv2,im1,im2,jv1,jv2,jv3,jw
  integer :: npoly,iok,notavail,itri
  integer :: ka, ka1, ka2, ka3
  real :: radcone
  real :: hcpn(4), vcpn(4), fldvals(4)
  real :: valw(mza), valw1(mza), valw2(mza), valw3(mza)
  real :: val1(mza), val2(mza), val3(mza)
  real :: wta1, wta2, wta3, wtb1, wtb2, wtb3

  real, parameter :: wtbot = 1., wttop = 0.

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  nu   = 0
  ipos = 0

  base = 12 * nbytes_real
  inc = ceiling( real(mwa) / 10. ) * mza

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop is over W for contouring W points
  ! Limit to primary W point or else we will go out of bounds for a parallel run

  do jw = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(jw)

     npoly = itab_w(iw)%npoly

     ! Loop over U/V neighbors of W, and for each, get indices for
     ! M endpoints and opposite W neighbor

     do jv2 = 1,npoly

        jv1 = jv2 - 1
        if (jv2 == 1) jv1 = npoly

        jv3 = jv2 + 1
        if (jv2 == npoly) jv3 = 1

        iv2 = itab_w(iw)%iv(jv2)

        im1 = itab_w(iw)%im(jv1) ! check these
        im2 = itab_w(iw)%im(jv2) ! check these

        iw1 = itab_w(iw)%iw(jv1)
        iw2 = itab_w(iw)%iw(jv2)
        iw3 = itab_w(iw)%iw(jv3)

        ka  = lpw(iw)
        ka1 = lpw(iw1)
        ka2 = lpw(iw2)
        ka3 = lpw(iw3)

        ! Loop over both triangles in current JV sector

        do itri = 1,2

           ! Check for intersection of cone surface and current triangle

           if (itri == 1) then
              call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
                   xem(im1),yem(im1),zem(im1),xev(iv2),yev(iv2),zev(iv2), &
                   wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
           else
              call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
                   xev(iv2),yev(iv2),zev(iv2),xem(im2),yem(im2),zem(im2), &
                   wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
           endif

           ! Jump to end of loop if cone surface does not intersect triangle

           if (iok /= 1) cycle

           ! Skip triangle if it crosses +/- 180 degree point along cone circle
           ! (In future, maybe re-program to truncate cells)

           if (mdomain < 2 .and. op%projectn(iplt) == 'C') then
              radcone = erad * sin(op%coneang * pio180)

              if (abs(hcpn(2) - hcpn(1)) > 3. * radcone) cycle
           endif

           ! Skip this triangle if entire cell is outside plot window

           if ( (hcpn(1) < op%xmin .and. hcpn(2) < op%xmin) .or.  &
                (hcpn(1) > op%xmax .and. hcpn(2) > op%xmax) ) cycle

           ! Fill arrays with field points

           do k = ka,mza
              call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw(k),notavail)
           enddo

           do k = ka2,mza
              call oplot_lib(k,iw2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw2(k),notavail)
           enddo

           if (itri == 1) then
              do k = ka1,mza
                 call oplot_lib(k,iw1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                valw1(k),notavail)
              enddo
           else
              do k = ka3,mza
                 call oplot_lib(k,iw3,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                valw3(k),notavail)
              enddo
           endif

           do k = ka-1,mza
              val1(k) = valw(max(k,ka))

              if (itri == 1) then
                 if (k >= ka .and. k >= ka1 .and. k >= ka2) then
                    val2(k) = (valw(k) + valw1(k) + valw2(k)) / 3.
                 elseif (k >= ka .and. k >= ka1) then
                    val2(k) = (valw(k) + valw1(k)) / 2.
                 elseif (k >= ka .and. k >= ka2) then
                    val2(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= ka1 .and. k >= ka2) then
                    val2(k) = (valw1(k) + valw2(k)) / 2.
                 elseif (k >= ka) then
                    val2(k) = valw(k)
                 elseif (k >= ka1) then
                    val2(k) = valw1(k)
                 elseif (k >= ka2) then
                    val2(k) = valw2(k)
                 elseif (ka == ka1 .and. ka == ka2) then
                    val2(k) = (valw(k+1) + valw1(k+1) + valw2(k+1)) / 3.
                 elseif (ka == ka1) then
                    val2(k) = (valw(k+1) + valw1(k+1)) / 2.
                 elseif (ka == ka2) then
                    val2(k) = (valw(k+1) + valw2(k+1)) / 2.
                 else
                    val2(k) = valw(k+1)
                 endif

                 if (k >= ka .and. k >= ka2) then
                    val3(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= ka) then
                    val3(k) = valw(k)
                 elseif (k >= ka2) then
                    val3(k) = valw2(k)
                 elseif (ka == ka2) then
                    val3(k) = (valw(k+1) + valw2(k+1)) / 2.
                 else
                    val3(k) = valw(k+1)
                 endif

              else

                 if (k >= ka .and. k >= ka2) then
                    val2(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= ka) then
                    val2(k) = valw(k)
                 elseif (k >= ka2) then
                    val2(k) = valw2(k)
                 elseif (ka == ka2) then
                    val2(k) = (valw(k+1) + valw2(k+1)) / 2.
                 else
                    val2(k) = valw(k+1)
                 endif

                 if (k >= ka .and. k >= ka2 .and. k >= ka3) then
                    val3(k) = (valw(k) + valw2(k) + valw3(k)) / 3.
                 elseif (k >= ka .and. k >= ka2) then
                    val3(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= ka .and. k >= ka3) then
                    val3(k) = (valw(k) + valw3(k)) / 2.
                 elseif (k >= ka2 .and. k >= ka3) then
                    val3(k) = (valw2(k) + valw3(k)) / 2.
                 elseif (k >= ka) then
                    val3(k) = valw(k)
                 elseif (k >= ka2) then
                    val3(k) = valw2(k)
                 elseif (k >= ka3) then
                    val3(k) = valw3(k)
                 elseif (ka == ka2 .and. ka == ka3) then
                    val3(k) = (valw(k+1) + valw2(k+1) + valw3(k+1)) / 3.
                 elseif (ka == ka2) then
                    val3(k) = (valw(k+1) + valw2(k+1)) / 2.
                 elseif (ka == ka3) then
                    val3(k) = (valw(k+1) + valw3(k+1)) / 2.
                 else
                    val3(k) = valw(k+1)
                 endif

              endif

           enddo

           do k = ka,mza   ! Loop is over T levels

              do ilev = 1,2

              ! Define vertical plot coordinates depending on whether top or bottom half
              ! (ilev = 1 or 2) of T cell is being contoured

                 if (ilev == 1) then

                    vcpn(1:2) = zm(k-1)
                    vcpn(3:4) = zt(k)

                    fldvals(1) = wta1 * val1(k-1) + wta2 * val2(k-1) + wta3 * val3(k-1)
                    fldvals(2) = wtb1 * val1(k-1) + wtb2 * val2(k-1) + wtb3 * val3(k-1)
                    fldvals(3) = wtb1 * val1(k  ) + wtb2 * val2(k  ) + wtb3 * val3(k  )
                    fldvals(4) = wta1 * val1(k  ) + wta2 * val2(k  ) + wta3 * val3(k  )

                    fldvals(1) = 0.5 * (fldvals(1) + fldvals(4))
                    fldvals(2) = 0.5 * (fldvals(2) + fldvals(3))

                 else ! ilev = 2

                    vcpn(1:2) = zt(k)
                    vcpn(3:4) = zm(k)

                    fldvals(1) = wta1 * val1(k) + wta2 * val2(k) + wta3 * val3(k)
                    fldvals(2) = wtb1 * val1(k) + wtb2 * val2(k) + wtb3 * val3(k)

                    if (k < mza) then
                       fldvals(3) = wtb1 * val1(k+1) + wtb2 * val2(k+1) + wtb3 * val3(k+1)
                       fldvals(4) = wta1 * val1(k+1) + wta2 * val2(k+1) + wta3 * val3(k+1)

                       fldvals(3) = 0.5 * (fldvals(3) + fldvals(2))
                       fldvals(4) = 0.5 * (fldvals(4) + fldvals(1))
                    else
                       fldvals(3) = fldvals(2)
                       fldvals(4) = fldvals(1)
                    endif

                 endif

              ! Skip this triangle if entire cell is above or below plot window

                 if ( all(vcpn(1:4) < op%ymin) .or. &
                      all(vcpn(1:4) > op%ymax) ) cycle

                 if (myrank == 0) then

                    ! Contour plot cell around current M point
                    call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)

                 else

#ifdef OLAM_MPI
                    nu = nu + 1
                    if (buffsize < ipos + base) then
                       allocate( bcopy (buffsize + inc * base) )
                       bcopy(1:buffsize) = buffer
                       call move_alloc(bcopy, buffer)
                       buffsize = size(buffer)
                    endif

                    call MPI_Pack(hcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(vcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(fldvals, 4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

                 endif

              enddo ! ilev

           enddo ! K

        enddo ! ITRI

     enddo ! JV

  enddo ! IW

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, 4, MPI_REAL, MPI_COMM_WORLD, ier)

                 ! Contour plot cell around current M point
                 call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  ! Now plot underground T cells with underground color
  call plot_underground_w(iplt,(/0/))

end subroutine contslab_vert_t

!===============================================================================

subroutine contslab_vert_w(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mza, lpw, zm, &
                        xem, yem, zem, xev, yev, zev, xew, yew, zew
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use misc_coms,  only: mdomain, iparallel
  use consts_coms,only: erad, pio180
  use mem_para,   only: myrank, mgroupsize, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: k,iw,iw1,iw2,iw3,iv2,im1,im2,jv1,jv2,jv3,jw
  integer :: npoly,iok,notavail,itri
  integer :: ka, ka1, ka2, ka3, kaw, kaw1, kaw2, kaw3
  real :: radcone
  real :: hcpn(4), vcpn(4), fldvals(4)
  real :: valw(mza), valw1(mza), valw2(mza), valw3(mza)
  real :: val1(mza), val2(mza), val3(mza)
  real :: wta1, wta2, wta3, wtb1, wtb2, wtb3

  real, parameter :: wtbot = 1., wttop = 0.

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  nu   = 0
  ipos = 0

  base = 12 * nbytes_real
  inc = ceiling( real(mwa) / 10. ) * mza

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop is over W for contouring W points
  ! Limit to primary W point or else we will go out of bounds for a parallel run

  do jw = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(jw)

     npoly = itab_w(iw)%npoly

     ! Loop over U/V neighbors of W, and for each, get indices for
     ! M endpoints and opposite W neighbor

     do jv2 = 1,npoly

        jv1 = jv2 - 1
        if (jv2 == 1) jv1 = npoly

        jv3 = jv2 + 1
        if (jv2 == npoly) jv3 = 1

        iv2 = itab_w(iw)%iv(jv2)

        im1 = itab_w(iw)%im(jv1) ! check these
        im2 = itab_w(iw)%im(jv2) ! check these

        iw1 = itab_w(iw)%iw(jv1)
        iw2 = itab_w(iw)%iw(jv2)
        iw3 = itab_w(iw)%iw(jv3)

        ka  = lpw(iw)
        ka1 = lpw(iw1)
        ka2 = lpw(iw2)
        ka3 = lpw(iw3)

        kaw  = lpw(iw)  - 1
        kaw1 = lpw(iw1) - 1
        kaw2 = lpw(iw2) - 1
        kaw3 = lpw(iw3) - 1

        ! Loop over both triangles in current JV sector

        do itri = 1,2

           ! Check for intersection of cone surface and current triangle

           if (itri == 1) then
              call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
                   xem(im1),yem(im1),zem(im1),xev(iv2),yev(iv2),zev(iv2), &
                   wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
           else
              call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
                   xev(iv2),yev(iv2),zev(iv2),xem(im2),yem(im2),zem(im2), &
                   wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
           endif

           ! Jump to end of loop if cone surface does not intersect triangle

           if (iok /= 1) cycle

           ! Skip triangle if it crosses +/- 180 degree point along cone circle
           ! (In future, maybe re-program to truncate cells)

           if (mdomain < 2 .and. op%projectn(iplt) == 'C') then
              radcone = erad * sin(op%coneang * pio180)

              if (abs(hcpn(2) - hcpn(1)) > 3. * radcone) cycle
           endif

           ! Skip this triangle if entire cell is outside plot window

           if ( (hcpn(1) < op%xmin .and. hcpn(2) < op%xmin) .or.  &
                (hcpn(1) > op%xmax .and. hcpn(2) > op%xmax) ) cycle

           ! Fill arrays with field points

           do k = kaw2,mza
              call oplot_lib(k,iw2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw2(k),notavail)
           enddo

           if (itri == 1) then
              do k = kaw1,mza
                 call oplot_lib(k,iw1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                valw1(k),notavail)
              enddo
           else
              do k = kaw3,mza
                 call oplot_lib(k,iw3,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                valw3(k),notavail)
              enddo
           endif

           do k = kaw,mza
              call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw(k),notavail)

              val1(k) = valw(k)

              if (itri == 1) then
                 if (k >= kaw1 .and. k >= kaw2) then
                    val2(k) = (valw(k) + valw1(k) + valw2(k)) / 3.
                    val3(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= kaw1) then
                    val2(k) = (valw(k) + valw1(k)) / 2.
                    val3(k) =  valw(k)
                 elseif (k >= kaw2) then
                    val2(k) = (valw(k) + valw2(k)) / 2.
                    val3(k) = (valw(k) + valw2(k)) / 2.
                 else
                    val2(k) = valw(k)
                    val3(k) = valw(k)
                 endif
              else
                 if (k >= kaw2 .and. k >= kaw3) then
                    val2(k) = (valw(k) + valw2(k)) / 2.
                    val3(k) = (valw(k) + valw2(k) + valw3(k)) / 3.
                 elseif (k >= kaw2) then
                    val2(k) = (valw(k) + valw2(k)) / 2.
                    val3(k) = (valw(k) + valw2(k)) / 2.
                 elseif (k >= kaw3) then
                    val2(k) =  valw(k)
                    val3(k) = (valw(k) + valw3(k)) / 2.
                 else
                    val2(k) = valw(k)
                    val3(k) = valw(k)
                 endif
              endif
           enddo

           do k = ka,mza   ! Loop is over T levels

              ! Define vertical plot coordinates

              vcpn(1:2) = zm(k-1)
              vcpn(3:4) = zm(k)

              ! Skip this triangle if entire cell is above or below plot window

              if ( all(vcpn(1:4) < op%ymin) .or. &
                   all(vcpn(1:4) > op%ymax) ) cycle

              ! W values interpolated to A and B points for plotting

              fldvals(1) = wta1 * valw(k-1) + wta2 * val2(k-1) + wta3 * val3(k-1)
              fldvals(2) = wtb1 * valw(k-1) + wtb2 * val2(k-1) + wtb3 * val3(k-1)
              fldvals(3) = wtb1 * valw(k  ) + wtb2 * val2(k  ) + wtb3 * val3(k  )
              fldvals(4) = wta1 * valw(k  ) + wta2 * val2(k  ) + wta3 * val3(k  )

              if (myrank == 0) then

                 ! Contour plot cell around current M point
                 call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)

              else

#ifdef OLAM_MPI
                 nu = nu + 1
                 if (buffsize < ipos + base) then
                    allocate( bcopy (buffsize + inc * base) )
                    bcopy(1:buffsize) = buffer
                    call move_alloc(bcopy, buffer)
                    buffsize = size(buffer)
                 endif

                 call MPI_Pack(hcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vcpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(fldvals, 4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

              endif

           enddo ! K

        enddo ! ITRI

     enddo ! JV

  enddo ! IW

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, hcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vcpn,    4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldvals, 4, MPI_REAL, MPI_COMM_WORLD, ier)

                 ! Contour plot cell around current M point
                 call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  ! Now plot underground T cells with underground color
  call plot_underground_w(iplt,(/0/))

end subroutine contslab_vert_w
