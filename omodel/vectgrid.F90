subroutine vectgrid_horiz_w(iplt)

  use oplot_coms,  only: op, op_vect, sort_plot_wts_parallel
  use mem_grid,    only: mza, mwa, xew, yew, zew, wnx, wny, wnz
  use mem_ijtabs,  only: jtab_w, jtw_prog, itabg_m, itab_m
  use mem_basic,   only: vxe, vye, vze
  use consts_coms, only: eradi
  use misc_coms,   only: mdomain, iparallel
  use mem_para,    only: myrank, mgroupsize, nbytes_int, nbytes_real
  use minterp_lib, only: minterp_wghts, get_weights_lonlat
  use map_proj,    only: or_ll, ps_ll, gn_ll, or_de, ps_de, gn_de, ll_ec

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: jw,iw

  real :: pointx,pointy,tailx,taily
  real :: headlen,head1x,head1y,head2x,head2y
  real :: tailxe,tailye,tailze,stemlen
  real :: stemx,stemy,stemz,snx,sny,snz,rnx,rny,rnz
  real :: head1xe,head1ye,head1ze,head2xe,head2ye,head2ze
  real :: vx,vy,vz,vert_vel

  integer :: ktf(mwa),kw(mwa)
  real :: wtbot(mwa),wttop(mwa)

  ! Wind barbs

  real, parameter :: pf = .3  ! full tick length as a fraction of shaft length
  real, parameter :: pi = .15 ! distance between ticks as a fraction of shaft length
  real, parameter :: ba = 50. ! value for which triangle is drawn
  real, parameter :: bb = 10. ! value for which tick is drawn - half tick for half of bb

  integer, parameter :: maxtris =  4 ! max wind 245 m/s
  integer, parameter :: maxbarb =  4
  integer, parameter :: maxhalf =  1
  integer            :: ntris, nbarb, nhalf

  real :: xtris(3,maxtris), ytris(3,maxtris)
  real :: xbarb(2,maxbarb), ybarb(2,maxbarb)
  real :: xhalf(2,maxhalf), yhalf(2,maxhalf)

  real :: speed, pc, xt, xea, yea, zea, xeb, yeb, zeb, xec, yec, zec
  real :: xa, ya, xb, yb, xd, yd

  real,    allocatable :: olat(:,:), olon(:,:)
  real,    allocatable :: oxe(:,:), oye(:,:), oze(:,:)
  logical, allocatable :: iskp(:,:)
  logical              :: plot_underg

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, j, n, is, kp
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  integer :: nx, ny, ix, iy, nrecvs, irank, ixp, ix0, ii, iip
  integer :: i, irecv, ierr, im, icolor, itmp, ic, npnts
  real :: dx, dy, wtu, wta, val, rr

  real :: vxw, vyw, vzw, wt, xeg, yeg, zeg

  integer :: ncount

  real, allocatable :: xc(:), yc(:), xg(:)
  type(minterp_wghts), allocatable :: wts(:)

  Type receive_bufs
     real, pointer :: rvals(:)
  End type receive_bufs

  real,       target, allocatable :: send_buf(:)
  type(receive_bufs), allocatable :: recv_buf(:)

  integer, allocatable :: ireqs(:)

  logical :: iskip

  ! Find cell K indices on the given plot surface

  op%stagpt = 'T'
  op%dimens = '3'

  plot_underg = ( op%noundrg(iplt) /= 'u' )

  nx = op%nx_vect
  dx = (op%xmax - op%xmin) / real(nx)

  ny = nint( (op%ymax - op%ymin) / dx )
  dy = (op%ymax - op%ymin) / real(ny)

  allocate(xc(nx))
  allocate(yc(ny))

  do iy = 1, ny
     yc(iy) = real(iy-1) * dy + (op%ymin + 0.5 * dy)
  enddo

  do ix = 1, nx
     xc(ix) = real(ix-1) * dx + (op%xmin + 0.5 * dx)
  enddo

  if (.not. op_vect(iplt)%initialized) then

     op_vect(iplt)%initialized = .true.

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

        allocate(xg(nx))
        do ix = 1, nx
           xg(ix) = xc(ix) + op%plon3
        enddo

        call get_weights_lonlat(xg,yc,wts,nx,ny)
        deallocate(xg)

     endif

     call sort_plot_wts_parallel(nx*ny,wts,op_vect(iplt))

  endif  ! initialized

!!  if (op%vectskip /= 1 .and. op_vect(iplt)%npnts > 0) then
!!     allocate(dopoint( op_vect(iplt)%npnts )) ; dopoint = .false.
!!
!!     ncount = 0
!!     do i  = 1, op_vect(iplt)%npnt
!!        ii = op_vect(iplt)%ipnts(i)
!!        iy = (ii-1) / nx + 1
!!        ix = mod(ii-1,nx) + 1
!!
!!        if (contplot .and. (ix == 1 .or. ix == nx)) cycle
!!        if (contplot .and. (iy == 1 .or. iy == ny)) cycle
!!
!!        if (mod(ix,op%vectskip) == 0 .and. mod(iy,op%vectskip) == 0) then
!!           ncount = ncount + 1
!!           dopoint(ii) = .true.
!!        endif
!!     enddo
!!  endif

  allocate( send_buf( op_vect(iplt)%npnts * 7 + 1 ) )

  if (myrank == 0) then
     allocate(recv_buf(0:op_vect(iplt)%nrecvs))
     recv_buf(0)%rvals => send_buf(:)
  endif

  nrecvs = 0

#ifdef OLAM_MPI
  if (iparallel == 1 .and. myrank == 0) then
     nrecvs = op_vect(iplt)%nrecvs
     if (nrecvs > 0) then

        allocate(ireqs(nrecvs)) ; ireqs = MPI_REQUEST_NULL

        do i = 1, nrecvs
           npnts = op_vect(iplt)%iremote(i)%npnts * 7 + 1
           irank = op_vect(iplt)%iremote(i)%irank
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

!!  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
!!
!!     call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
!!                    fldval(iw),notavail(iw))
!!
!!  enddo
!!
!!  if (iparallel == 1) then
!!     call mpi_send_w(r1dvara1=fldval, i1dvara1=notavail)
!!     call mpi_recv_w(r1dvara1=fldval, i1dvara1=notavail)
!!  endif


  ncount = 0
  ipos   = 0

  ! Earth coordinates of left and right head tips
  headlen = op%headspeed * op%dtvec

  do i  = 1, op_vect(iplt)%npnts
     ii = op_vect(iplt)%ipnts(i)
     iy = (ii-1) / nx + 1
     ix = mod(ii-1,nx) + 1

     ! We are beyond the limit of the projection
     if (op_vect(iplt)%wgts(i)%imglobe == 1) cycle

     ! Local im triangle corresponding to this point
     im = itabg_m( op_vect(iplt)%wgts(i)%imglobe )%im_myrank

     wt = 0.0
     vx = 0.0
     vy = 0.0
     vz = 0.0

     do n = 1, 3
        iw = itab_m(im)%iw(n)
        if (plot_underg .and. ktf(iw) /= 0) cycle
        kp = min(kw(iw)+1,mza)

        vxw = wtbot(iw) * vxe(kw(iw),iw) + wttop(iw) * vxe(kp,iw)
        vyw = wtbot(iw) * vye(kw(iw),iw) + wttop(iw) * vye(kp,iw)
        vzw = wtbot(iw) * vze(kw(iw),iw) + wttop(iw) * vze(kp,iw)

        ! Projection onto local vertical vector
        vert_vel = vx * wnx(iw) + vy * wny(iw) + vz * wnz(iw)

        ! Subtract local vertical wind vector to give horizontal wind
        vxw = vxw - vert_vel * wnx(iw)
        vyw = vyw - vert_vel * wny(iw)
        vzw = vzw - vert_vel * wnz(iw)

        vx = vx + vxw * op_vect(iplt)%wgts(i)%wt(n)
        vy = vy + vyw * op_vect(iplt)%wgts(i)%wt(n)
        vz = vz + vzw * op_vect(iplt)%wgts(i)%wt(n)
        wt = wt +       op_vect(iplt)%wgts(i)%wt(n)
     enddo

     ! skip if majority of cell is underground or unavailable
     if (wt < 0.4) cycle

     vx = vx / wt
     vy = vy / wt
     vz = vz / wt

     speed = sqrt(vx**2 + vy**2 + vz**2)

     if (speed < 1.e-5) cycle

     pointx = xc(ix)
     pointy = yc(iy)

     iskip = .false.

     if (op%projectn(iplt) == 'O') then

        call or_de( xeg, yeg, zeg, op%cosplat, op%sinplat, &
                    op%cosplon, op%sinplon, pointx, pointy, iskip )
        xeg = xeg + op%pxe
        yeg = yeg + op%pye
        zeg = zeg + op%pze

     elseif (op%projectn(iplt) == 'P') then

        call ps_de( xeg, yeg, zeg, op%cosplat, op%sinplat, &
                    op%cosplon, op%sinplon, pointx, pointy )
        xeg = xeg + op%pxe
        yeg = yeg + op%pye
        zeg = zeg + op%pze

     elseif (op%projectn(iplt) == 'G') then

        call gn_de( xeg, yeg, zeg, op%cosplat, op%sinplat, &
                    op%cosplon, op%sinplon, pointx, pointy )
        xeg = xeg + op%pxe
        yeg = yeg + op%pye
        zeg = zeg + op%pze

     elseif (op%projectn(iplt) == 'L') then

        call ll_ec( pointx+op%plon3, pointy, xeg, yeg, zeg )

     endif

     if (iskip) cycle

     ncount = ncount + 1
     rr = transfer(ii,rr)

     ! Copy to send buffer
     send_buf( (ncount-1)*7+2:(ncount-1)*7+8 ) = [rr,vx,vy,vz,xeg,yeg,zeg]

  enddo

  send_buf(1) = transfer(ncount,send_buf(1))

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (myrank > 0 .and. op_vect(iplt)%npnts > 0) then
        call MPI_Send(send_buf, ncount*7+1, MPI_REAL, 0, itag, MPI_COMM_WORLD, ierr)
     endif
  endif
#endif

  if (myrank == 0) then

     do itmp = 0, nrecvs
        irecv = itmp

#ifdef OLAM_MPI
        if (itmp > 0) then
           call MPI_Waitany(nrecvs, ireqs, irecv, MPI_STATUS_IGNORE, ierr)
        endif
#endif

        npnts = transfer(recv_buf(irecv)%rvals(1), npnts)

        do i = 1, npnts

           ii = transfer( recv_buf(irecv)%rvals( (i-1)*7+2 ), ii)
           iy = (ii-1) / nx + 1
           ix = mod(ii-1,nx) + 1

           pointx = xc(ix)
           pointy = yc(iy)

           vx  = recv_buf(irecv)%rvals( (i-1)*7+3 )
           vy  = recv_buf(irecv)%rvals( (i-1)*7+4 )
           vz  = recv_buf(irecv)%rvals( (i-1)*7+5 )
           xeg = recv_buf(irecv)%rvals( (i-1)*7+6 )
           yeg = recv_buf(irecv)%rvals( (i-1)*7+7 )
           zeg = recv_buf(irecv)%rvals( (i-1)*7+8 )

           speed = sqrt( vx*vx + vy*vy + vz*vz )

           ! 3D vector displacement (in time interval op%dtvec)...

           if (op%vectbarb(iplt) == 'w') then

              ! Plot total horizontal vector at W point

              stemx = vx * op%dtvec
              stemy = vy * op%dtvec
              stemz = vz * op%dtvec

           else    ! case for op%vectbarb(iplt) == 'B'

              ! Plot wind barb at W point

              stemx = vx * op%stemlength / speed
              stemy = vy * op%stemlength / speed
              stemz = vz * op%stemlength / speed

           endif

           ! Vector length and unit components

           stemlen = sqrt(stemx**2 + stemy**2 + stemz**2)

           snx = stemx / stemlen
           sny = stemy / stemlen
           snz = stemz / stemlen

           ! "Right" unit components

           if (mdomain <= 1) then  ! Spherical geometry case
              rnx = (sny * zeg - snz * yeg) * eradi
              rny = (snz * xeg - snx * zeg) * eradi
              rnz = (snx * yeg - sny * xeg) * eradi
           else                    ! Cartesian case
              rnx =  sny
              rny = -snx
              rnz =  0.
           endif

           ! Earth coordinates of tail

           tailxe = xeg - stemx
           tailye = yeg - stemy
           tailze = zeg - stemz

           if (op%vectbarb(iplt) /= 'B') then

              head1xe = xeg + rnx * .42 * headlen - snx * .91 * headlen
              head1ye = yeg + rny * .42 * headlen - sny * .91 * headlen
              head1ze = zeg + rnz * .42 * headlen - snz * .91 * headlen

              head2xe = xeg - rnx * .42 * headlen - snx * .91 * headlen
              head2ye = yeg - rny * .42 * headlen - sny * .91 * headlen
              head2ze = zeg - rnz * .42 * headlen - snz * .91 * headlen

              ! Transform other tail and coordinates

              call oplot_transform_xyz(iplt,tailxe, tailye, tailze, tailx, taily)
              call oplot_transform_xyz(iplt,head1xe,head1ye,head1ze,head1x,head1y)
              call oplot_transform_xyz(iplt,head2xe,head2ye,head2ze,head2x,head2y)

              ! Avoid wrap-around

              if (op%projectn(iplt) == 'L') then
                 call ll_unwrap(pointx,tailx)
                 call ll_unwrap(pointx,head1x)
                 call ll_unwrap(pointx,head2x)
              endif

              ! Draw vector

              call o_frstpt(tailx,taily)
              call o_vector(pointx,pointy)
              call o_frstpt(head1x,head1y)
              call o_vector(pointx,pointy)
              call o_vector(head2x,head2y)

           else    ! case for op%vectbarb(iplt) == 'B'

              ! Transform stem tail

              call oplot_transform_xyz(iplt,tailxe,tailye,tailze,tailx,taily)

              ! Draw stem

              call o_frstpt(tailx,taily)
              call o_vector(pointx,pointy)

              pc = 1.
              xt = speed + .25 * bb

              ! Draw triangles (if any)

              ntris = 0
              do while (xt >= ba)

                 xea = xeg + (tailxe - xeg) * pc
                 yea = yeg + (tailye - yeg) * pc
                 zea = zeg + (tailze - zeg) * pc

                 xeb = xea - rnx * pf * op%stemlength
                 yeb = yea - rny * pf * op%stemlength
                 zeb = zea - rnz * pf * op%stemlength

                 pc = pc - pi

                 xec = xeg + (tailxe - xeg) * pc
                 yec = yeg + (tailye - yeg) * pc
                 zec = zeg + (tailze - zeg) * pc

                 ! Transform triangle coordinates

                 call oplot_transform_xyz(iplt,xea,yea,zea,xa,ya)
                 call oplot_transform_xyz(iplt,xeb,yeb,zeb,xb,yb)
                 call oplot_transform_xyz(iplt,xec,yec,zec,xd,yd)

                 ! Avoid wrap-around

                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)
                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xd)

                 call o_frstpt(xa,ya)
                 call o_vector(xb,yb)
                 call o_vector(xd,yd)

                 pc = pc - pi
                 xt = xt - ba

              enddo

              ! Draw barbs (if any)

              nbarb = 0
              do while (xt >= bb)

                 xea = xeg + (tailxe - xeg) * pc
                 yea = yeg + (tailye - yeg) * pc
                 zea = zeg + (tailze - zeg) * pc

                 xeb = xea - rnx * pf * op%stemlength
                 yeb = yea - rny * pf * op%stemlength
                 zeb = zea - rnz * pf * op%stemlength

                 ! Transform triangle coordinates

                 call oplot_transform_xyz(iplt,xea,yea,zea,xa,ya)
                 call oplot_transform_xyz(iplt,xeb,yeb,zeb,xb,yb)

                 ! Avoid wrap-around

                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

                 call o_frstpt(xa,ya)
                 call o_vector(xb,yb)

                 pc = pc - pi
                 xt = xt - bb

              enddo

              ! Draw half barb (if any)

              nhalf = 0
              if (xt >= .5*bb) then

                 xea = xeg + (tailxe - xeg) * pc
                 yea = yeg + (tailye - yeg) * pc
                 zea = zeg + (tailze - zeg) * pc

                 xeb = xea - rnx * .5 * pf * op%stemlength
                 yeb = yea - rny * .5 * pf * op%stemlength
                 zeb = zea - rnz * .5 * pf * op%stemlength

                 ! Transform triangle coordinates

                 call oplot_transform_xyz(iplt,xea,yea,zea,xa,ya)
                 call oplot_transform_xyz(iplt,xeb,yeb,zeb,xb,yb)

                 ! Avoid wrap-around

                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
                 if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

                 call o_frstpt(xa,ya)
                 call o_vector(xb,yb)

              endif

           endif

        enddo

     enddo

  endif

end subroutine vectgrid_horiz_w
