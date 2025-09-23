subroutine tileslab_grid_tw(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa
  use mem_ijtabs,   only: jtab_w, jtw_prog, itabg_m, itab_m
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll, get_weights_lonlat

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: ktf(mwa), kw(mwa)
  real :: wtbot(mwa), wttop(mwa)

  real :: fldval(mwa)
  integer :: notavail(mwa)

  real, allocatable :: xf(:), xc(:)
  real, allocatable :: yf(:), yc(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  integer, allocatable :: icolors(:)
  logical              :: plot_underg

  real :: xf4(4), yf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     integer, allocatable :: ivals(:)
  End type receive_bufs

  integer,            allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer :: j, iw, nx, ny, ix, iy, nrecvs, irank, ixp, ix0, ii, iip
  integer :: i, irecv, ierr, im, icolor, itmp, npnts
  real :: dx, dy, wtu, wta, val

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  allocate(xf(nx+1))
  allocate(xc(nx))

  allocate(yf(ny+1))
  allocate(yc(ny))

  do iy = 1, ny+1
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  do iy = 1, ny
     yc(iy) = real(iy-1) * dy + (op%ymin + 0.5 * dy)
  enddo

  do ix = 1, nx+1
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  do ix = 1, nx
     xc(ix) = real(ix-1) * dx + (op%xmin + 0.5 * dx)
  enddo

!!  if (op%projectn(iplt) == 'L') then
!!     if (op%windowin(iplt) /= 'W') then
!!        xf(:) = xf(:) - op%plon3
!!!       xc(:) = xc(:) - op%plon3
!!     endif
!!  endif

!!  write(*,*) op%plon3, op%plat3
!!  write(*,*) op%xmax, op%xmin
!!  write(*,*) op%ymax, op%ymin
!!
!!  write(*,*) xf
!!  write(*,*)
!!  write(*,*) yf
!!  stop

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        do ix = 1, nx
           xc(ix) = xc(ix) + op%plon3
        enddo

        call get_weights_lonlat(xc,yc,wts,nx,ny)

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
           npnts = op_grid(iplt)%iremote(i)%npnts + 2
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%ivals( npnts ) )

           call MPI_Irecv( recv_buf(i)%ivals, npnts, MPI_INT, irank, itag, &
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

     if ( ktf(iw) /= 0 ) then
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

  if (myrank == 0) allocate(icolors(nx*ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts+2))

     do ii = 1, npnts

        if (op_grid(iplt)%wgts(ii)%imglobe == 1) then

           ! We are beyond the limit of the projection
           icolor = icback

        else

           im = itabg_m( op_grid(iplt)%wgts(ii)%imglobe )%im_myrank

           wtu = 0.0
           if (plot_underg) then
              if ( any( ktf( itab_m(im)%iw(1:3) ) == 1 ) ) then
                 do i = 1, 3
                    iw = itab_m(im)%iw(i)
                    if (ktf(iw) == 1) wtu = wtu + op_grid(iplt)%wgts(ii)%wt(i)
                 enddo
              endif
           endif

           if (wtu > 0.6) then

              icolor = op%icigrnd

           else

              wta = 0.0
              val = 0.0
              do i = 1, 3
                 iw = itab_m(im)%iw(i)
                 if ( (plot_underg .and. ktf(iw) == 1) .or. (notavail(iw) > 0) ) cycle
                 val = val + fldval(iw) * op_grid(iplt)%wgts(ii)%wt(i)
                 wta = wta + op_grid(iplt)%wgts(ii)%wt(i)
              enddo

              if (wta > 0.4) then

                 val = val / wta
                 call getcolor(iplt,val,icolor)

                 op%fldval_min = min(val,op%fldval_min)
                 op%fldval_max = max(val,op%fldval_max)

              else

                 icolor = 7

              endif

           endif
        endif

        if (iparallel == 1) then
           if (myrank == 0) then
              icolors( op_grid(iplt)%iremote(0)%ipnts(ii) ) = icolor
           else
              send_buf(ii) = icolor
           endif
        else
           icolors(ii) = icolor
        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        send_buf(npnts+1:npnts+2) = transfer( [op%fldval_min,op%fldval_max], send_buf(npnts+1:npnts+2) )
        call MPI_Send(send_buf, npnts+2, MPI_INT, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i = 1, npnts
              icolors( op_grid(iplt)%iremote(irecv)%ipnts(i) ) = recv_buf(irecv)%ivals(i)
           enddo

           op%fldval_min = min( op%fldval_min, transfer( recv_buf(irecv)%ivals(npnts+1), op%fldval_min) )
           op%fldval_max = max( op%fldval_max, transfer( recv_buf(irecv)%ivals(npnts+2), op%fldval_max) )
        enddo
     endif

  endif
#endif

  if (myrank == 0) then

!    do ii = 1, nx*ny
!       iy = (ii-1) / nx + 1
!       ix = mod(ii-1,nx) + 1
!       xf4 = [ xf(ix), xf(ix),   xf(ix+1), xf(ix+1) ]
!       yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]
!       call o_sfsgfa(xf4,yf4,4,icolors(ii))
!    enddog

     do iy = 1, ny
        ix0 = 1
        do ix = 1, nx
           ii  = (iy-1)*nx + ix
           ixp = min(ix+1,nx)
           iip = (iy-1)*nx + ixp

           if (ix == nx .or. icolors(iip) /= icolors(ii)) then

              if (icolors(ii) /= icback) then
                 xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                 yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]

                 call o_sfsgfa(xf4,yf4,4,icolors(ii))
              endif
              ix0 = ixp

           endif

        enddo
     enddo

  endif

end subroutine tileslab_grid_tw

!===============================================================================

subroutine tileslab_grid_mp(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa, mma
  use mem_ijtabs,   only: jtab_m, jtm_prog, itabg_m, itab_m
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
! use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll, get_weights_lonlat

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: ktf(mwa), km(mma)
  real :: wtbot(mma), wttop(mma)

  real :: fldval(mma)
  integer :: notavail(mma)

  real, allocatable :: xf(:), xc(:)
  real, allocatable :: yf(:), yc(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  integer, allocatable :: icolors(:)
  logical              :: plot_underg

  real :: xf4(4), yf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     integer, allocatable :: ivals(:)
  End type receive_bufs

  integer,            allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer :: j, im, nx, ny, ix, iy, nrecvs, irank, ixp, ix0, ii, iip
  integer :: i, iw, irecv, ierr, icolor, itmp, npnts
  real :: dx, dy, wtu

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  allocate(xf(nx+1))
  allocate(xc(nx))

  allocate(yf(ny+1))
  allocate(yc(ny))

  do iy = 1, ny+1
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  do iy = 1, ny
     yc(iy) = real(iy-1) * dy + (op%ymin + 0.5 * dy)
  enddo

  do ix = 1, nx+1
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  do ix = 1, nx
     xc(ix) = real(ix-1) * dx + (op%xmin + 0.5 * dx)
  enddo

!!  if (op%projectn(iplt) == 'L') then
!!     if (op%windowin(iplt) /= 'W') then
!!        xf(:) = xf(:) - op%plon3
!!!       xc(:) = xc(:) - op%plon3
!!     endif
!!  endif

!!  write(*,*) op%plon3, op%plat3
!!  write(*,*) op%xmax, op%xmin
!!  write(*,*) op%ymax, op%ymin
!!
!!  write(*,*) xf
!!  write(*,*)
!!  write(*,*) yf
!!  stop

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        do ix = 1, nx
           xc(ix) = xc(ix) + op%plon3
        enddo

        call get_weights_lonlat(xc,yc,wts,nx,ny)

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
           npnts = op_grid(iplt)%iremote(i)%npnts + 2
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%ivals( npnts ) )

           call MPI_Irecv( recv_buf(i)%ivals, npnts, MPI_INT, irank, itag, &
                           MPI_COMM_WORLD, ireqs(i), ierr )
        enddo

     endif
  endif
#endif

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

  ! Loop over M points

  !$omp parallel do private(im)
  do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

     call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                    fldval(im),notavail(im))

  enddo
  !$omp end parallel do

! Not sure this is needed
! if (iparallel == 1) then
!    call mpi_send_w(r1dvara1=fldval, i1dvara1=notavail)
!    call mpi_recv_w(r1dvara1=fldval, i1dvara1=notavail)
! endif

  if (myrank == 0) allocate(icolors(nx*ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts+2))

     do ii = 1, npnts

        if (op_grid(iplt)%wgts(ii)%imglobe == 1) then

           ! We are beyond the limit of the projection
           icolor = icback

        else

           im = itabg_m( op_grid(iplt)%wgts(ii)%imglobe )%im_myrank

           wtu = 0.0
           if (plot_underg) then
              if ( any( ktf( itab_m(im)%iw(1:3) ) == 1 ) ) then
                 do i = 1, 3
                    iw = itab_m(im)%iw(i)
                    if (ktf(iw) == 1) wtu = wtu + op_grid(iplt)%wgts(ii)%wt(i)
                 enddo
              endif
           endif

           if (wtu > 0.6) then

              ! Cell is majority underground
              icolor = op%icigrnd

           elseif (notavail(im) > 0) then

              icolor = 7

           else

              call getcolor(iplt,fldval(im),icolor)

              op%fldval_min = min(fldval(im),op%fldval_min)
              op%fldval_max = max(fldval(im),op%fldval_max)

           endif
        endif

        if (iparallel == 1) then
           if (myrank == 0) then
              icolors( op_grid(iplt)%iremote(0)%ipnts(ii) ) = icolor
           else
              send_buf(ii) = icolor
           endif
        else
           icolors(ii) = icolor
        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        send_buf(npnts+1:npnts+2) = transfer( [op%fldval_min,op%fldval_max], send_buf(npnts+1:npnts+2) )
        call MPI_Send(send_buf, npnts+2, MPI_INT, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i = 1, npnts
              icolors( op_grid(iplt)%iremote(irecv)%ipnts(i) ) = recv_buf(irecv)%ivals(i)
           enddo

           op%fldval_min = min( op%fldval_min, transfer( recv_buf(irecv)%ivals(npnts+1), op%fldval_min) )
           op%fldval_max = max( op%fldval_max, transfer( recv_buf(irecv)%ivals(npnts+2), op%fldval_max) )
        enddo
     endif

  endif
#endif

  if (myrank == 0) then

!    do ii = 1, nx*ny
!       iy = (ii-1) / nx + 1
!       ix = mod(ii-1,nx) + 1
!       xf4 = [ xf(ix), xf(ix),   xf(ix+1), xf(ix+1) ]
!       yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]
!       call o_sfsgfa(xf4,yf4,4,icolors(ii))
!    enddog

     do iy = 1, ny
        ix0 = 1
        do ix = 1, nx
           ii  = (iy-1)*nx + ix
           ixp = min(ix+1,nx)
           iip = (iy-1)*nx + ixp

           if (ix == nx .or. icolors(iip) /= icolors(ii)) then

              if (icolors(ii) /= icback) then
                 xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                 yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]

                 call o_sfsgfa(xf4,yf4,4,icolors(ii))
              endif
              ix0 = ixp

           endif

        enddo
     enddo

  endif

end subroutine tileslab_grid_mp

!===============================================================================

subroutine tileslab_grid_wsfc(iplt)

  use oplot_coms,   only: op, op_grid, sort_plot_wts_parallel
  use mem_grid,     only: mwa
  use mem_ijtabs,   only: jtab_w, jtw_prog, itabg_m, itab_m, itab_w
  use misc_coms,    only: iparallel
  use mem_para,     only: myrank
  use olam_mpi_sfc, only: mpi_send_wsfc, mpi_recv_wsfc
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use minterp_lib,  only: minterp_wghts, get_weights_lonlat
  use map_proj,     only: or_ll, ps_ll, gn_ll, get_weights_lonlat
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

  real, allocatable :: xf(:), xc(:)
  real, allocatable :: yf(:), yc(:)

  real,    allocatable :: olat(:,:), olon(:,:)
  logical, allocatable :: iskp(:,:)
  integer, allocatable :: icolors(:)
  logical              :: plot_underg

  real :: xf4(4), yf4(4)

  integer, allocatable :: ireqs(:)

  Type receive_bufs
     integer, allocatable :: ivals(:)
  End type receive_bufs

  integer,            allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer :: j, iw, nx, ny, ix, iy, nrecvs, irank, ixp, ix0, ii, iip, iwsfc
  integer :: i, irecv, ierr, im, icolor, itmp, npnts, ks, notavailsfc
  real :: dx, dy, wta, val

  integer, parameter :: itag = 999

  type(minterp_wghts), allocatable :: wts(:)

  integer, parameter :: icback = 7
  real, parameter :: wtbot = 1., wttop = 0.

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_grid
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  allocate(xf(nx+1))
  allocate(xc(nx))

  allocate(yf(ny+1))
  allocate(yc(ny))

  do iy = 1, ny+1
     yf(iy) = real(iy-1) * dy + op%ymin
  enddo

  do iy = 1, ny
     yc(iy) = real(iy-1) * dy + (op%ymin + 0.5 * dy)
  enddo

  do ix = 1, nx+1
     xf(ix) = real(ix-1) * dx + op%xmin
  enddo

  do ix = 1, nx
     xc(ix) = real(ix-1) * dx + (op%xmin + 0.5 * dx)
  enddo

!!  if (op%projectn(iplt) == 'L') then
!!     if (op%windowin(iplt) /= 'W') then
!!        xf(:) = xf(:) - op%plon3
!!!       xc(:) = xc(:) - op%plon3
!!     endif
!!  endif

!!  write(*,*) op%plon3, op%plat3
!!  write(*,*) op%xmax, op%xmin
!!  write(*,*) op%ymax, op%ymin
!!
!!  write(*,*) xf
!!  write(*,*)
!!  write(*,*) yf
!!  stop

  if (.not. op_grid(iplt)%initialized) then

     op_grid(iplt)%initialized = .true.

     allocate(wts(nx*ny))

     if (op%projectn(iplt) == 'O') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))
        allocate(iskp(nx,ny))

        call or_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny, iskp)
        call get_weights_lonlat(olon,olat,wts,nx,ny,iskp)

        deallocate(olat,olon,iskp)

     elseif (op%projectn(iplt) == 'P') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call ps_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'G') then

        allocate(olat(nx,ny))
        allocate(olon(nx,ny))

        call gn_ll(olat, olon, op%plat3, op%plon3, xc, yc, nx, ny)
        call get_weights_lonlat(olon,olat,wts,nx,ny)

        deallocate(olat,olon)

     elseif (op%projectn(iplt) == 'L') then

        do ix = 1, nx
           xc(ix) = xc(ix) + op%plon3
        enddo

        call get_weights_lonlat(xc,yc,wts,nx,ny)

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
           npnts = op_grid(iplt)%iremote(i)%npnts + 2
           irank = op_grid(iplt)%iremote(i)%irank
           allocate( recv_buf(i)%ivals( npnts ) )

           call MPI_Irecv( recv_buf(i)%ivals, npnts, MPI_INT, irank, itag, &
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

  if (myrank == 0) allocate(icolors(nx*ny))

  npnts = op_grid(iplt)%npnts

  if (npnts > 0) then

     if (myrank /= 0) allocate(send_buf(npnts+2))

     do ii = 1, npnts

        if (op_grid(iplt)%wgts(ii)%imglobe == 1) then

           ! We are beyond the limit of the projection
           icolor = icback

        else

           im = itabg_m( op_grid(iplt)%wgts(ii)%imglobe )%im_myrank

           wta = 0.0
           val = 0.0
           do i = 1, 3
              iw = itab_m(im)%iw(i)
              if ( notavail(iw) > 0 ) cycle
              val = val + fldval(iw) * op_grid(iplt)%wgts(ii)%wt(i)
              wta = wta + op_grid(iplt)%wgts(ii)%wt(i)
           enddo

           if (wta >= 0.2) then

              val = val / wta
              call getcolor(iplt,val,icolor)

              op%fldval_min = min(val,op%fldval_min)
              op%fldval_max = max(val,op%fldval_max)

           else

              icolor = 7

           endif

        endif

        if (iparallel == 1) then
           if (myrank == 0) then
              icolors( op_grid(iplt)%iremote(0)%ipnts(ii) ) = icolor
           else
              send_buf(ii) = icolor
           endif
        else
           icolors(ii) = icolor
        endif

     enddo
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     if (myrank > 0 .and. npnts > 0) then
        send_buf(npnts+1:npnts+2) = transfer( [op%fldval_min,op%fldval_max], send_buf(npnts+1:npnts+2) )
        call MPI_Send(send_buf, npnts+2, MPI_INT, 0, itag, MPI_COMM_WORLD, ierr)
     endif

     if (myrank == 0) then
        do itmp = 1, nrecvs
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)

           npnts = op_grid(iplt)%iremote(irecv)%npnts
           do i = 1, npnts
              icolors( op_grid(iplt)%iremote(irecv)%ipnts(i) ) = recv_buf(irecv)%ivals(i)
           enddo

           op%fldval_min = min( op%fldval_min, transfer( recv_buf(irecv)%ivals(npnts+1), op%fldval_min) )
           op%fldval_max = max( op%fldval_max, transfer( recv_buf(irecv)%ivals(npnts+2), op%fldval_max) )
        enddo
     endif

  endif
#endif

  if (myrank == 0) then

!    do ii = 1, nx*ny
!       iy = (ii-1) / nx + 1
!       ix = mod(ii-1,nx) + 1
!       xf4 = [ xf(ix), xf(ix),   xf(ix+1), xf(ix+1) ]
!       yf4 = [ yf(iy), yf(iy+1), yf(iy+1), yf(iy)   ]
!       call o_sfsgfa(xf4,yf4,4,icolors(ii))
!    enddog

     do iy = 1, ny
        ix0 = 1
        do ix = 1, nx
           ii  = (iy-1)*nx + ix
           ixp = min(ix+1,nx)
           iip = (iy-1)*nx + ixp

           if (ix == nx .or. icolors(iip) /= icolors(ii)) then

              if (icolors(ii) /= icback) then
                 xf4 = [ xf(ix0), xf(ix0) , xf(ix+1), xf(ix+1) ]
                 yf4 = [ yf(iy),  yf(iy+1), yf(iy+1), yf(iy)   ]

                 call o_sfsgfa(xf4,yf4,4,icolors(ii))
              endif
              ix0 = ixp

           endif

        enddo
     enddo

  endif

end subroutine tileslab_grid_wsfc

!===============================================================================

subroutine getcolor(iplt,fldval,icolor)

  use oplot_coms, only: op
  use plotcolors, only: clrtab

  implicit none

  integer, intent(in)  :: iplt
  real,    intent(in)  :: fldval
  integer, intent(out) :: icolor

  real    :: fldval1
  integer :: itab, ival

  itab = op%icolortab(iplt)

  ! Cyclic treatment of color palette (used with integer-type data)

  fldval1 = fldval

  if (clrtab(itab)%ifmt(1) == 20) &
!      fldval1 = mod(fldval-1.,real(clrtab(itab)%nvals-2)) - 1. ! table 124 has 2 neg vals
       fldval1 = mod(fldval,   real(clrtab(itab)%nvals-2)) - 1. ! table 124 has 2 neg vals

  ! Case for color table 130

  if (clrtab(itab)%ifmt(1) == 30) &
       fldval1 = mod(fldval,100.)

  ! Extract contour color from color table

  ival = 1
  do while (fldval1 > clrtab(itab)%vals(ival) .and. &
            ival < clrtab(itab)%nvals               )
     ival = ival + 1
  enddo
  icolor = clrtab(itab)%ipal(ival)

end subroutine getcolor

!===============================================================================

subroutine sfccell_pltaverage(iw, field, rvalue)

  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: mwsfc, itab_wsfc

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: field(mwsfc)
  real,    intent(out) :: rvalue

  integer :: j, iwsfc, jasfc

  rvalue = 0.

  do j = 1, itab_w(iw)%jsfc2
     iwsfc = itab_w(iw)%iwsfc(j)
     jasfc = itab_w(iw)%jasfc(j)

     rvalue = rvalue + field(iwsfc) * itab_wsfc(iwsfc)%arcoariw(jasfc)
  enddo

end subroutine sfccell_pltaverage

!===============================================================================

subroutine landcell_pltaverage(iw, field, rvalue)

  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: sfcg, mwsfc

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: field(mwsfc)
  real,    intent(out) :: rvalue

  integer :: j, iwsfc
  real    :: area

  rvalue = 0.
  area   = 0.

  do j = itab_w(iw)%jland1, itab_w(iw)%jland2
     iwsfc = itab_w(iw)%iwsfc(j)

     rvalue = rvalue + field(iwsfc) * sfcg%area(iwsfc)
     area   = area   +                sfcg%area(iwsfc)
  enddo

  rvalue = rvalue / area

end subroutine landcell_pltaverage

!===============================================================================

subroutine seacell_pltaverage(iw, field, rvalue)

  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: sfcg, mwsfc

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: field(mwsfc)
  real,    intent(out) :: rvalue

  integer :: j, iwsfc
  real    :: area

  rvalue = 0.
  area   = 0.

  do j = itab_w(iw)%jsea1, itab_w(iw)%jsea2
     iwsfc = itab_w(iw)%iwsfc(j)

     rvalue = rvalue + field(iwsfc) * sfcg%area(iwsfc)
     area   = area   +                sfcg%area(iwsfc)
  enddo

  rvalue = rvalue / area

end subroutine seacell_pltaverage

!===============================================================================

subroutine lakecell_pltaverage(iw, field, rvalue)

  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: sfcg, mwsfc

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: field(mwsfc)
  real,    intent(out) :: rvalue

  integer :: j, iwsfc
  real    :: area

  rvalue = 0.
  area   = 0.

  do j = itab_w(iw)%jlake1, itab_w(iw)%jlake2
     iwsfc = itab_w(iw)%iwsfc(j)

     rvalue = rvalue + field(iwsfc) * sfcg%area(iwsfc)
     area   = area   +                sfcg%area(iwsfc)
  enddo

  rvalue = rvalue / area

end subroutine lakecell_pltaverage

!===============================================================================
