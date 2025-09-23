subroutine contslab_grid_tw(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa
  use mem_ijtabs,   only: jtab_w, jtw_prog, itabg_m, itab_m
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll
  use plotcolors,   only: clrtab

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: ktf(mwa), kw(mwa)
  real :: wtbot(mwa), wttop(mwa)

  real :: fldval(mwa)
  integer :: notavail(mwa)

  real, allocatable :: xf(:), xg(:)
  real, allocatable :: yf(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  real,    allocatable :: fieldxy(:,:)
  logical              :: plot_underg

  real :: xf3(3), yf3(3), zf3(3)
  real :: xf4(4), yf4(4), zf4(4)

  logical :: lf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     real, allocatable :: rvals(:)
  End type receive_bufs

  real,               allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer, allocatable :: icolor(:)

  integer :: j, iw, nx, ny, ix, iy, npnts, nrecvs, ix0, ixp, ilev, itab
  integer :: i, irecv, ierr, im, itmp, ii, irank, n, nm, nu, jm, jp
  real :: dx, dy, wtu, wta, val, fmax, fmin

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7
  real,    parameter :: rmiss = -.999e30
  real,    parameter :: rundg = -.999e29

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  nx = nx+1
  ny = ny+1

  allocate(yf(ny))

  do iy = 1, ny
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  allocate(xf(nx))

  do ix = 1, nx
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        allocate(xg(nx))
        do ix = 1, nx
           xg(ix) = xf(ix) + op%plon3
        enddo

        call get_weights_lonlat(xg,yf,wts,nx,ny)
        deallocate(xg)

     endif

     call sort_plot_wts_parallel(nx*ny,wts,op_grid(iplt))

  endif  ! initialized

#ifdef OLAM_MPI
  if (iparallel == 1 .and. myrank == 0) then

     nrecvs = op_grid(iplt)%nrecvs
     if (nrecvs > 0) then

        allocate(recv_buf(nrecvs))
        allocate(ireqs   (nrecvs)) ; ireqs = MPI_REQUEST_NULL

        do i = 1, nrecvs
           npnts = op_grid(iplt)%iremote(i)%npnts
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%rvals( npnts ) )
           call MPI_Irecv( recv_buf(i)%rvals, npnts, MPI_REAL, irank, itag, &
                           MPI_COMM_WORLD, ireqs(i), ierr )
        enddo

     endif
  endif
#endif

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

  ! Loop over W points

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     if ( ktf(iw) /= 0) then
        notavail(iw) = 3
     else
        call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                       fldval(iw),notavail(iw))
     endif

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_w(r1dvara1=fldval, i1dvara1=notavail)
     call mpi_recv_w(r1dvara1=fldval, i1dvara1=notavail)
  endif

  if (myrank == 0) allocate(fieldxy(nx,ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts))

     do i = 1, npnts

        if (op_grid(iplt)%wgts(i)%imglobe == 1) then

           ! We are beyond the limit of the projection
           val = rmiss

        else

           im = itabg_m( op_grid(iplt)%wgts(i)%imglobe )%im_myrank

           wtu = 0.0
           if (plot_underg) then
              if ( any( ktf( itab_m(im)%iw(1:3) ) == 1 ) ) then
                 do j = 1, 3
                    iw = itab_m(im)%iw(j)
                    if (ktf(iw) == 1) wtu = wtu + op_grid(iplt)%wgts(i)%wt(j)
                 enddo
              endif
           endif

           ! Point is mostly underground
           if (wtu > 0.6) then

              val = rundg

           else

              wta = 0.0
              val = 0.0
              do j = 1, 3
                 iw = itab_m(im)%iw(j)
                 if ( (plot_underg .and. ktf(iw) == 1) .or. (notavail(iw) > 0) ) cycle
                 val = val + fldval(iw) * op_grid(iplt)%wgts(i)%wt(j)
                 wta = wta + op_grid(iplt)%wgts(i)%wt(j)
              enddo

              if (wta > 0.4) then

                 ! Point is available and above ground
                 val = val / wta

                 op%fldval_min = min(val,op%fldval_min)
                 op%fldval_max = max(val,op%fldval_max)

              else

                 ! Value is unavailable or below ground
                 val = rmiss

              endif

           endif
        endif

        if (iparallel == 1) then

           if (myrank == 0) then
              ii = op_grid(iplt)%iremote(0)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = val
           else
              send_buf(i) = val
           endif

        else

           iy = (i-1) / nx + 1
           ix = mod(i-1,nx) + 1
           fieldxy(ix,iy) = val

        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        call MPI_Send(send_buf, npnts, MPI_REAL, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i  = 1, npnts
              ii = op_grid(iplt)%iremote(irecv)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = recv_buf(irecv)%rvals(i)
           enddo

        enddo
     endif

  endif
#endif

  if (myrank == 0) then

     op%fldval_min = minval(fieldxy, mask = fieldxy > -.99e29)
     op%fldval_max = maxval(fieldxy, mask = fieldxy > -.99e29)

     allocate(icolor(nx))

     itab = op%icolortab(iplt)

     do iy = 1, ny-1

        icolor(:) = -1
        do ix = 1, nx-1

           zf4 = [ fieldxy(ix:ix+1,iy:iy+1) ]
           nm = count( zf4 == rmiss)
           nu = count( zf4 == rundg)

           if (nu >= 2) then

              icolor(ix) = op%icigrnd

           elseif (nm >= 2) then

              icolor(ix) = icback

           elseif ( nm == 0 .and. nu == 0 ) then

              fmax = maxval( zf4 )
              fmin = minval( zf4 )
              call getclev(itab,fmin,ilev)

              if ( clrtab(itab)%vals(ilev) <= fmin .or. clrtab(itab)%vals(ilev) >= fmax ) then
                 icolor(ix) = clrtab(itab)%ipal(ilev)
              endif

           endif

        enddo

        ix0 = 1
        do ix = 1, nx-1
           ixp = min(ix+1,nx-1)

           if (ix == nx-1 .or. icolor(ix) < 0 .or. icolor(ixp) < 0 .or. icolor(ixp) /= icolor(ix)) then

              if (icolor(ix) > 0) then

                 if (op%ifill > 0 .and. icolor(ix) /= icback) then
                    xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]
                    call o_sfsgfa(xf4,yf4,4,icolor(ix))
                 endif

              else

                 zf4 = [ fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy) ]
                 lf4 = ( zf4 /= rmiss .and. zf4 /= rundg )
                 n   = count( lf4 )

                 if (n == 4) then

                    ! contour top left triangle if not missing or underground

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), &
                         xf(ix), xf(ix),   xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy+1) )

                    ! contour bot right triangle if not missing or underground

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy), &
                         xf(ix), xf(ix+1), xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy) )

                 elseif (n == 3) then

                    ! contour the triangle defined by the 3 points that are available

                    xf4 = [ xf(ix), xf(ix)  , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]

                    xf3 = pack( xf4, mask=lf4 )
                    yf3 = pack( yf4, mask=lf4 )
                    zf3 = pack( zf4, mask=lf4 )

                    call cont3( op%icolortab(iplt), op%ifill, &
                                zf3(1), zf3(2), zf3(3), &
                                xf3(1), xf3(2), xf3(3), &
                                yf3(1), yf3(2), yf3(3)  )

                    ! plots the other half of this cell if it is underground

                    do j = 1, 4
                       if ( zf4(j) == rundg ) then
                          jm = modulo(j-2,4) + 1
                          jp = modulo(j  ,4) + 1

                          xf3 = [ xf4(jm), xf4(j), xf4(jp) ]
                          yf3 = [ yf4(jm), yf4(j), yf4(jp) ]
                          call o_sfsgfa(xf3,yf3,3,op%icigrnd)

                          exit
                       endif
                    enddo

                 endif

              endif

              ix0 = ixp
           endif

        enddo

     enddo

  endif

end subroutine contslab_grid_tw

!===============================================================================

subroutine contslab_grid_mp(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa, mma
  use mem_ijtabs,   only: jtab_m, jtm_prog, itabg_m, itab_m
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
! use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll
  use plotcolors,   only: clrtab

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: ktf(mwa), km(mma)
  real :: wtbot(mma), wttop(mma)

  real :: fldval(mma)
  integer :: notavail(mma)

  real, allocatable :: xf(:), xg(:)
  real, allocatable :: yf(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  real,    allocatable :: fieldxy(:,:)
  logical              :: plot_underg

  real :: xf3(3), yf3(3), zf3(3)
  real :: xf4(4), yf4(4), zf4(4)
  logical :: lf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     real, allocatable :: rvals(:)
  End type receive_bufs

  real,               allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer, allocatable :: icolor(:)

  integer :: j, iw, nx, ny, ix, iy, npnts, nrecvs, ix0, ixp, ilev, itab
  integer :: i, irecv, ierr, im, itmp, ii, irank, n, nm, nu, jm, jp
  real :: dx, dy, wtu, val, fmax, fmin

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7
  real,    parameter :: rmiss = -.999e30
  real,    parameter :: rundg = -.999e29

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  nx = nx+1
  ny = ny+1

  allocate(yf(ny))

  do iy = 1, ny
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  allocate(xf(nx))

  do ix = 1, nx
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        allocate(xg(nx))
        do ix = 1, nx
           xg(ix) = xf(ix) + op%plon3
        enddo

        call get_weights_lonlat(xg,yf,wts,nx,ny)
        deallocate(xg)

     endif

     call sort_plot_wts_parallel(nx*ny,wts,op_grid(iplt))

  endif  ! initialized

#ifdef OLAM_MPI
  if (iparallel == 1 .and. myrank == 0) then

     nrecvs = op_grid(iplt)%nrecvs
     if (nrecvs > 0) then

        allocate(recv_buf(nrecvs))
        allocate(ireqs   (nrecvs)) ; ireqs = MPI_REQUEST_NULL

        do i = 1, nrecvs
           npnts = op_grid(iplt)%iremote(i)%npnts
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%rvals( npnts ) )
           call MPI_Irecv( recv_buf(i)%rvals, npnts, MPI_REAL, irank, itag, &
                           MPI_COMM_WORLD, ireqs(i), ierr )
        enddo

     endif
  endif
#endif

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

  ! Loop over W points

  !$omp parallel do private(im)
  do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

     call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                    fldval(im),notavail(im))

  enddo
  !$omp end parallel do

! Not sure this is necessary
! if (iparallel == 1) then
!    call mpi_send_m(r1dvara1=fldval, i1dvara1=notavail)
!    call mpi_recv_m(r1dvara1=fldval, i1dvara1=notavail)
! endif

  if (myrank == 0) allocate(fieldxy(nx,ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts))

     do i = 1, npnts

        if (op_grid(iplt)%wgts(i)%imglobe == 1) then

           ! We are beyond the limit of the projection
           val = rmiss

        else

           im = itabg_m( op_grid(iplt)%wgts(i)%imglobe )%im_myrank

           wtu = 0.0
           if (plot_underg) then
              if ( any( ktf( itab_m(im)%iw(1:3) ) == 1 ) ) then
                 do j = 1, 3
                    iw = itab_m(im)%iw(j)
                    if (ktf(iw) == 1) wtu = wtu + op_grid(iplt)%wgts(i)%wt(j)
                 enddo
              endif
           endif

           if (wtu > 0.6) then

              ! Point is mostly underground
              val = rundg

           elseif (notavail(im) > 0) then

              ! Point is unavailable
              val = rmiss

           else

              ! Point is available and above ground
              val = fldval(im)

              op%fldval_min = min(val,op%fldval_min)
              op%fldval_max = max(val,op%fldval_max)
           endif

        endif

        if (iparallel == 1) then

           if (myrank == 0) then
              ii = op_grid(iplt)%iremote(0)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = val
           else
              send_buf(i) = val
           endif

        else

           iy = (i-1) / nx + 1
           ix = mod(i-1,nx) + 1
           fieldxy(ix,iy) = val

        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        call MPI_Send(send_buf, npnts, MPI_REAL, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i  = 1, npnts
              ii = op_grid(iplt)%iremote(irecv)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = recv_buf(irecv)%rvals(i)
           enddo

        enddo
     endif

  endif
#endif

  if (myrank == 0) then

     op%fldval_min = minval(fieldxy, mask = fieldxy > -.99e29)
     op%fldval_max = maxval(fieldxy, mask = fieldxy > -.99e29)

     allocate(icolor(nx))

     itab = op%icolortab(iplt)

     do iy = 1, ny-1

        icolor(:) = -1
        do ix = 1, nx-1

           zf4 = [ fieldxy(ix:ix+1,iy:iy+1) ]
           nm = count( zf4 == rmiss)
           nu = count( zf4 == rundg)

           if (nu >= 2) then

              icolor(ix) = op%icigrnd

           elseif (nm >= 2) then

              icolor(ix) = icback

           elseif ( nm == 0 .and. nu == 0 ) then

              fmax = maxval( zf4 )
              fmin = minval( zf4 )
              call getclev(itab,fmin,ilev)

              if ( clrtab(itab)%vals(ilev) <= fmin .or. clrtab(itab)%vals(ilev) >= fmax ) then
                 icolor(ix) = clrtab(itab)%ipal(ilev)
              endif

           endif

        enddo

        ix0 = 1
        do ix = 1, nx-1
           ixp = min(ix+1,nx-1)

           if (ix == nx-1 .or. icolor(ix) < 0 .or. icolor(ixp) < 0 .or. icolor(ixp) /= icolor(ix)) then

              if (icolor(ix) > 0) then

                 if (op%ifill > 0 .and. icolor(ix) /= icback) then
                    xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]
                    call o_sfsgfa(xf4,yf4,4,icolor(ix))
                 endif

              else

                 zf4 = [ fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy) ]
                 lf4 = ( zf4 /= rmiss .and. zf4 /= rundg )
                 n   = count( lf4 )

                 if (n == 4) then

                    ! contour top left triangle if not missing or underground

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), &
                         xf(ix), xf(ix),   xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy+1) )

                    ! contour bot right triangle if not missing or underground

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy), &
                         xf(ix), xf(ix+1), xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy) )

                 elseif (n == 3) then

                    ! contour the triangle defined by the 3 points that are available

                    xf4 = [ xf(ix), xf(ix)  , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]

                    xf3 = pack( xf4, mask=lf4 )
                    yf3 = pack( yf4, mask=lf4 )
                    zf3 = pack( zf4, mask=lf4 )

                    call cont3( op%icolortab(iplt), op%ifill, &
                                zf3(1), zf3(2), zf3(3), &
                                xf3(1), xf3(2), xf3(3), &
                                yf3(1), yf3(2), yf3(3)  )

                    ! plots the other half of this cell if it is underground

                    do j = 1, 4
                       if ( zf4(j) == rundg ) then
                          jm = modulo(j-2,4) + 1
                          jp = modulo(j  ,4) + 1

                          xf3 = [ xf4(jm), xf4(j), xf4(jp) ]
                          yf3 = [ yf4(jm), yf4(j), yf4(jp) ]
                          call o_sfsgfa(xf3,yf3,3,op%icigrnd)

                          exit
                       endif
                    enddo

                 endif

              endif

              ix0 = ixp
           endif

        enddo

     enddo

  endif

end subroutine contslab_grid_mp

!===============================================================================

subroutine contslab_grid_wsfc(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa
  use mem_ijtabs,   only: jtab_w, jtw_prog, itabg_m, itab_m, itab_w
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll
  use plotcolors,   only: clrtab
  use mem_sfcg,     only: mwsfc, itab_wsfc, sfcg
  use mem_land,     only: nzg

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  real :: sfcval(mwsfc)
  real :: fldval(mwa)
  integer :: notavail(mwa)

  real, allocatable :: xf(:), xg(:)
  real, allocatable :: yf(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  real,    allocatable :: fieldxy(:,:)
  logical              :: plot_underg

  real :: xf4(4), yf4(4), zf4(4)
  real :: xf3(3), yf3(3), zf3(3)
  logical :: lf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     real, allocatable :: rvals(:)
  End type receive_bufs

  real,               allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer, allocatable :: icolor(:)

  integer :: j, iw, nx, ny, ix, iy, npnts, nrecvs, ix0, ixp, ilev, itab, n
  integer :: i, irecv, ierr, im, itmp, ii, irank, ks, iwsfc, notavailsfc
  real :: dx, dy, wta, val, fmax, fmin

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7
  real,    parameter :: rmiss = -.999e30
  real,    parameter :: wtbot = 1., wttop = 0.

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  nx = nx+1
  ny = ny+1

  allocate(yf(ny))

  do iy = 1, ny
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  allocate(xf(nx))

  do ix = 1, nx
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xf, yf, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        allocate(xg(nx))
        do ix = 1, nx
           xg(ix) = xf(ix) + op%plon3
        enddo

        call get_weights_lonlat(xg,yf,wts,nx,ny)
        deallocate(xg)

     endif

     call sort_plot_wts_parallel(nx*ny,wts,op_grid(iplt))

  endif  ! initialized

#ifdef OLAM_MPI
  if (iparallel == 1 .and. myrank == 0) then

     nrecvs = op_grid(iplt)%nrecvs
     if (nrecvs > 0) then

        allocate(recv_buf(nrecvs))
        allocate(ireqs   (nrecvs)) ; ireqs = MPI_REQUEST_NULL

        do i = 1, nrecvs
           npnts = op_grid(iplt)%iremote(i)%npnts
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%rvals( npnts ) )
           call MPI_Irecv( recv_buf(i)%rvals, npnts, MPI_REAL, irank, itag, &
                           MPI_COMM_WORLD, ireqs(i), ierr )
        enddo

     endif
  endif
#endif

  ! Find K level to plot if field is 3d

  if (trim(op%dimens) == '3G') then
     ks = min(nzg,max(1,nint(op%slabloc(iplt))))
  else
     ks = 1
  endif

  ! Loop over WSFC points

  !$omp parallel do private(iwsfc)
  do iwsfc = 2, mwsfc
     notavail(iwsfc) = 3
     fldval  (iwsfc) = -999.

     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     if (op%stagpt == 'L' .and. sfcg%leaf_class(iwsfc) <  2) cycle
     if (op%stagpt == 'R' .and. sfcg%leaf_class(iwsfc) /= 1) cycle
     if (op%stagpt == 'S' .and. sfcg%leaf_class(iwsfc) /= 0) cycle

     call oplot_lib(ks,iwsfc,'VALUE',op%fldname(iplt),wtbot,wttop, &
                    sfcval(iwsfc),notavailsfc)

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_wsfc('plot_avg', rvar=sfcval)
     call mpi_recv_wsfc('plot_avg', rvar=sfcval)
  endif

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     if (op%stagpt == 'L') then

        if (itab_w(iw)%jland2 == 0) then
           notavail(iw) = 3
        else
           notavail(iw) = 0
           call landcell_pltaverage(iw, sfcval, fldval(iw))
        endif

     elseif (op%stagpt == 'R') then

        if (itab_w(iw)%jlake2 == 0) then
           notavail(iw) = 3
        else
           notavail(iw) = 0
           call lakecell_pltaverage(iw, sfcval, fldval(iw))
        endif

     elseif (op%stagpt == 'S') then

        if (itab_w(iw)%jsea2 == 0) then
           notavail(iw) = 3
        else
           notavail(iw) = 0
           call seacell_pltaverage(iw, sfcval, fldval(iw))
        endif

     else

        notavail(iw) = 0
        call sfccell_pltaverage(iw, sfcval, fldval(iw))

     endif

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_w(r1dvara1=fldval, i1dvara1=notavail)
     call mpi_recv_w(r1dvara1=fldval, i1dvara1=notavail)
  endif

  if (myrank == 0) allocate(fieldxy(nx,ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts))

     do i = 1, npnts

        if (op_grid(iplt)%wgts(i)%imglobe == 1) then

           ! We are beyond the limit of the projection
           val = rmiss

        else

           im = itabg_m( op_grid(iplt)%wgts(i)%imglobe )%im_myrank

           wta = 0.0
           val = 0.0
           do j = 1, 3
              iw = itab_m(im)%iw(j)
              if ( notavail(iw) > 0 ) cycle
              val = val + fldval(iw) * op_grid(iplt)%wgts(i)%wt(j)
              wta = wta + op_grid(iplt)%wgts(i)%wt(j)
           enddo

           if (wta > 0.2) then

              ! Point is available and above ground
              val = val / wta

              op%fldval_min = min(val,op%fldval_min)
              op%fldval_max = max(val,op%fldval_max)

           else

              ! Value is unavailable or below ground
              val = rmiss

           endif
        endif

        if (iparallel == 1) then

           if (myrank == 0) then
              ii = op_grid(iplt)%iremote(0)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = val
           else
              send_buf(i) = val
           endif

        else

           iy = (i-1) / nx + 1
           ix = mod(i-1,nx) + 1
           fieldxy(ix,iy) = val

        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        call MPI_Send(send_buf, npnts, MPI_REAL, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i  = 1, npnts
              ii = op_grid(iplt)%iremote(irecv)%ipnts(i)
              iy = (ii-1) / nx + 1
              ix = mod(ii-1,nx) + 1
              fieldxy(ix,iy) = recv_buf(irecv)%rvals(i)
           enddo

        enddo
     endif

  endif
#endif

  if (myrank == 0) then

     op%fldval_min = minval(fieldxy, mask = fieldxy > -.99e29)
     op%fldval_max = maxval(fieldxy, mask = fieldxy > -.99e29)

     allocate(icolor(nx))

     itab = op%icolortab(iplt)

     do iy = 1, ny-1

        icolor(:) = -1
        do ix = 1, nx-1

           zf4 = [ fieldxy(ix:ix+1,iy:iy+1) ]
           n = count( zf4 == rmiss)

           if ( n >= 2 ) then

              icolor(ix) = icback

           elseif ( n == 0 ) then

              fmax = maxval( zf4 )
              fmin = minval( zf4 )
              call getclev(itab,fmin,ilev)

              if ( clrtab(itab)%vals(ilev) <= fmin .or. clrtab(itab)%vals(ilev) >= fmax ) then
                 icolor(ix) = clrtab(itab)%ipal(ilev)
              endif

           endif

        enddo

        ix0 = 1
        do ix = 1, nx-1
           ixp = min(ix+1,nx-1)

           if (ix == nx-1 .or. icolor(ix) < 0 .or. icolor(ixp) < 0 .or. icolor(ixp) /= icolor(ix)) then

              if (icolor(ix) > 0) then

                 if (op%ifill > 0 .and. icolor(ix) /= icback) then
                    xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]
                    call o_sfsgfa(xf4,yf4,4,icolor(ix))
                 endif

              else

                 zf4 = [ fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy) ]
                 lf4 = ( zf4 /= rmiss )
                 n   = count( lf4 )

                 if (n == 4) then

                    ! contour top left triangle if not missing

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix,iy+1), fieldxy(ix+1,iy+1), &
                         xf(ix), xf(ix),   xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy+1) )

                    ! contour bot right triangle if not missing

                    call cont3( op%icolortab(iplt), op%ifill, &
                         fieldxy(ix,iy), fieldxy(ix+1,iy+1), fieldxy(ix+1,iy), &
                         xf(ix), xf(ix+1), xf(ix+1), &
                         yf(iy), yf(iy+1), yf(iy) )

                 elseif (n == 3) then

                    ! contour the triangle defined by the 3 points that are available

                    xf4 = [ xf(ix), xf(ix)  , xf(ix+1), xf(ix+1) ]
                    yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]

                    xf3 = pack( xf4, mask=lf4 )
                    yf3 = pack( yf4, mask=lf4 )
                    zf3 = pack( zf4, mask=lf4 )

                    call cont3( op%icolortab(iplt), op%ifill, &
                                zf3(1), zf3(2), zf3(3), &
                                xf3(1), xf3(2), xf3(3), &
                                yf3(1), yf3(2), yf3(3)  )

                 endif

              endif

              ix0 = ixp
           endif

        enddo

     enddo

  endif

end subroutine contslab_grid_wsfc

!===============================================================================

subroutine getclev(itab,fldval,ilev)

  use plotcolors, only: clrtab
  implicit none

  integer, intent(in)  :: itab
  real,    intent(in)  :: fldval
  integer, intent(out) :: ilev

  ilev = 1
  do while (fldval > clrtab(itab)%vals(ilev) .and. &
            ilev   < clrtab(itab)%nvals            )
     ilev = ilev + 1
  enddo

end subroutine getclev

!===============================================================================
