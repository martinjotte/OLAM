subroutine expand_delaunay_mesh(ires_factor, iatm)

  use mem_delaunay, only: nmd, nud, nwd
  implicit none

  integer, intent(in) :: ires_factor
  logical, intent(in) :: iatm

  integer :: nn, n2, n3, i

  n2 = 0
  n3 = 0

  nn = ires_factor

  if (nn > 1) then

     do while( mod(nn,2) == 0 )
        n2 = n2 + 1
        nn = nn / 2
     enddo

     do while( mod(nn,3) == 0 )
        n3 = n3 + 1
        nn = nn / 3
     enddo

     if (nn /= 1) then
        stop "Error: sfcgrid_res_factor must be 1 or have factors of 2 and 3."
     endif

     if (iatm) then
        write(*,*) "Before mesh expansion:"
     else
        write(*,*) "Before mesh expansion for sfc grid:"
     endif
     write(*,'(A,I9)') " nmd = ", nmd
     write(*,'(A,I9)') " nud = ", nud
     write(*,'(A,I9)') " nwd = ", nwd

     do i = 1, n3
        call expand_global3(iatm)
     enddo

     do i = 1, n2
        call expand_global2(iatm)
     enddo

  endif

end subroutine expand_delaunay_mesh



subroutine expand_global2(iatmgrid)

  ! This subroutine adds nested grid regions at the beginning of a simulation.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          nest_ud_vars, nest_wd_vars, alloc_itabsd, &
                          nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd, iwdorig, iwdorig_temp

  use mem_ijtabs,   only: mloops, &
                          jtm_grid, jtu_grid, jtw_grid, &
                          jtm_init, jtu_init, jtw_init, &
                          jtm_prog, jtu_prog, jtw_prog, &
                          jtm_wadj, jtu_wadj, jtw_wadj, &
                          jtm_wstn, jtu_wstn, jtw_wstn, &
                          jtm_vadj,           jtw_vadj

  use misc_coms,    only: mdomain
  use consts_coms,  only: pio180, erad, pi1, pi2

  implicit none

  logical, intent(in) :: iatmgrid

  type (itab_md_vars), allocatable :: ltab_md(:)
  type (itab_ud_vars), allocatable :: ltab_ud(:)
  type (itab_wd_vars), allocatable :: ltab_wd(:)

  type (nest_ud_vars), allocatable :: nest_ud(:)
  type (nest_wd_vars), allocatable :: nest_wd(:)

  integer :: iu,iw,im,iw1,iw2,im1,im2
  integer :: iu1,iu2,iu3,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1
  integer :: iu4,iu5,iu6,iw3,ngr,mrlo,mrloo,mrow

  integer :: nmd0,nud0,nwd0,iskip

  real :: expansion

  real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

  real :: xp1, xp2, xq1, xq2
  real :: yp1, yp2, yq1, yq2

  ! Make duplicate of current grid dimensions

  nmd0 = nmd
  nud0 = nud
  nwd0 = nwd

  ! Allocate temporary tables

  allocate (nest_ud(nud), nest_wd(nwd))  ! Nest relations

  ! Copy ITAB information and M-point coordinates to temporary tables

  call move_alloc(itab_md, ltab_md)
  call move_alloc(itab_ud, ltab_ud)
  call move_alloc(itab_wd, ltab_wd)

  call move_alloc(xemd, xem_temp)
  call move_alloc(yemd, yem_temp)
  call move_alloc(zemd, zem_temp)

  if (.not. iatmgrid) call move_alloc(iwdorig, iwdorig_temp)

  do iw = 2, nwd
     nest_wd(iw)%iu(1) = nud0 + 1
     nest_wd(iw)%iu(2) = nud0 + 2
     nest_wd(iw)%iu(3) = nud0 + 3

     nest_wd(iw)%iw(1) = nwd0 + 1
     nest_wd(iw)%iw(2) = nwd0 + 2
     nest_wd(iw)%iw(3) = nwd0 + 3

     nud0 = nud0 + 3
     nwd0 = nwd0 + 3
  enddo

  do iu = 2,nud
     nest_ud(iu)%im = nmd0 + 1
     nest_ud(iu)%iu = nud0 + 1

     nmd0 = nmd0 + 1
     nud0 = nud0 + 1
  enddo

  ! Allocate main tables to expanded size
  ! Initialize all neighbor indices to zero

  call alloc_itabsd(nmd0,nud0,nwd0)

  if (.not. iatmgrid) allocate(iwdorig(nwd0))

  ! Memory copy to main tables

  if (.not. iatmgrid) then
     iwdorig(1:nwd) = iwdorig_temp
  endif

  xemd(1:nmd) = xem_temp
  yemd(1:nmd) = yem_temp
  zemd(1:nmd) = zem_temp

  do im = 1,nmd
     if (iatmgrid) itab_md(im)%loop(1:mloops) = ltab_md(im)%loop(1:mloops)
     itab_md(im)%imp       = ltab_md(im)%imp
     itab_md(im)%mrlm      = ltab_md(im)%mrlm
     itab_md(im)%mrlm_orig = ltab_md(im)%mrlm_orig
     itab_md(im)%ngr       = ltab_md(im)%ngr
  enddo

  itab_ud(1:nud) = ltab_ud
  itab_wd(1:nwd) = ltab_wd

  !$omp parallel
  !$omp do private(im,im1,im2,expansion,mrlo,mrloo,ngr)
  do iu = 2,nud
     im  = nest_ud(iu)%im
     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     ! Average coordinates to new M points

     xemd(im) = .5 * (xemd(im1) + xemd(im2))
     yemd(im) = .5 * (yemd(im1) + yemd(im2))
     zemd(im) = .5 * (zemd(im1) + zemd(im2))

     ! If mdomain <= 1, push M point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt( xemd(im) ** 2 &
                               + yemd(im) ** 2 &
                               + zemd(im) ** 2 )

        xemd(im) = xemd(im) * expansion
        yemd(im) = yemd(im) * expansion
        zemd(im) = zemd(im) * expansion

     endif

     mrlo  = itab_md(im1)%mrlm
     mrloo = itab_md(im1)%mrlm_orig
     ngr   = itab_md(im1)%ngr

     itab_md(im)%mrlm      = mrlo
     itab_md(im)%mrlm_orig = mrloo
     itab_md(im)%ngr       = ngr
  enddo
  !$omp end do

  ! Contruct tables for new fully subdivided triangles

  !$omp do private(iu1o,iu2o,iu3o,mrlo,mrloo,ngr,mrow,iu1o_iw1,iu2o_iw1, &
  !$omp            iu3o_iw1,iu1,iu2,iu3,iu4,iu5,iu6,iw1,iw2,iw3)
  do iw = 2,nwd

     ! Mapping of original undivided triangles

     iu1o  = ltab_wd(iw)%iu(1)
     iu2o  = ltab_wd(iw)%iu(2)
     iu3o  = ltab_wd(iw)%iu(3)
     mrlo  = ltab_wd(iw)%mrlw
     mrloo = ltab_wd(iw)%mrlw_orig
     ngr   = ltab_wd(iw)%ngr
     mrow  = ltab_wd(iw)%mrow

     iu1o_iw1 = ltab_ud(iu1o)%iw(1)
     iu2o_iw1 = ltab_ud(iu2o)%iw(1)
     iu3o_iw1 = ltab_ud(iu3o)%iw(1)

     ! Mapping of new divided triangles

     iu1 = nest_wd(iw)%iu(1)
     iu2 = nest_wd(iw)%iu(2)
     iu3 = nest_wd(iw)%iu(3)

     iu4 = nest_ud(iu1o)%iu
     iu5 = nest_ud(iu2o)%iu
     iu6 = nest_ud(iu3o)%iu

     iw1 = nest_wd(iw)%iw(1)
     iw2 = nest_wd(iw)%iw(2)
     iw3 = nest_wd(iw)%iw(3)

     ! Fill tables with new values

     itab_wd(iw)%iu(1) = iu1
     itab_wd(iw)%iu(2) = iu2
     itab_wd(iw)%iu(3) = iu3

     itab_wd(iw1)%iu(1) = iu1
     itab_wd(iw2)%iu(1) = iu2
     itab_wd(iw3)%iu(1) = iu3

     itab_wd(iw1)%mrlw = mrlo
     itab_wd(iw2)%mrlw = mrlo
     itab_wd(iw3)%mrlw = mrlo

     itab_wd(iw1)%mrlw_orig = mrloo
     itab_wd(iw2)%mrlw_orig = mrloo
     itab_wd(iw3)%mrlw_orig = mrloo

     itab_wd(iw1)%ngr = ngr
     itab_wd(iw2)%ngr = ngr
     itab_wd(iw3)%ngr = ngr

     itab_wd(iw1)%mrow = mrow
     itab_wd(iw2)%mrow = mrow
     itab_wd(iw3)%mrow = mrow

     itab_ud(iu1o)%im(2) = nest_ud(iu1o)%im
     itab_ud(iu4)%im(1)  = nest_ud(iu1o)%im
     itab_ud(iu4)%im(2)  = ltab_ud(iu1o)%im(2)

     itab_ud(iu2o)%im(2) = nest_ud(iu2o)%im
     itab_ud(iu5)%im(1)  = nest_ud(iu2o)%im
     itab_ud(iu5)%im(2)  = ltab_ud(iu2o)%im(2)

     itab_ud(iu3o)%im(2) = nest_ud(iu3o)%im
     itab_ud(iu6)%im(1)  = nest_ud(iu3o)%im
     itab_ud(iu6)%im(2)  = ltab_ud(iu3o)%im(2)

     if (.not. iatmgrid) then
        iwdorig(iw1) = iwdorig_temp(iw)
        iwdorig(iw2) = iwdorig_temp(iw)
        iwdorig(iw3) = iwdorig_temp(iw)
     endif

     if (iw == iu1o_iw1) then

        itab_wd(iw3)%iu(2) = iu1o
        itab_wd(iw2)%iu(3) = iu4

        itab_ud(iu1)%im(1) = nest_ud(iu2o)%im
        itab_ud(iu1)%im(2) = nest_ud(iu3o)%im
        itab_ud(iu1)%iw(1) = iw1
        itab_ud(iu1)%iw(2) = iw

        itab_ud(iu1o)%iw(1) = iw3
        itab_ud(iu4)%iw(1) = iw2

     else

        itab_wd(iw3)%iu(2) = iu4
        itab_wd(iw2)%iu(3) = iu1o

        itab_ud(iu1)%im(1) = nest_ud(iu3o)%im
        itab_ud(iu1)%im(2) = nest_ud(iu2o)%im
        itab_ud(iu1)%iw(1) = iw
        itab_ud(iu1)%iw(2) = iw1

        itab_ud(iu1o)%iw(2) = iw2
        itab_ud(iu4)%iw(2) = iw3

     endif

     if (iw == iu2o_iw1) then

        itab_wd(iw1)%iu(2) = iu2o
        itab_wd(iw3)%iu(3) = iu5

        itab_ud(iu2)%im(1) = nest_ud(iu3o)%im
        itab_ud(iu2)%im(2) = nest_ud(iu1o)%im
        itab_ud(iu2)%iw(1) = iw2
        itab_ud(iu2)%iw(2) = iw

        itab_ud(iu2o)%iw(1) = iw1
        itab_ud(iu5)%iw(1) = iw3

     else

        itab_wd(iw1)%iu(2) = iu5
        itab_wd(iw3)%iu(3) = iu2o

        itab_ud(iu2)%im(1) = nest_ud(iu1o)%im
        itab_ud(iu2)%im(2) = nest_ud(iu3o)%im
        itab_ud(iu2)%iw(1) = iw
        itab_ud(iu2)%iw(2) = iw2

        itab_ud(iu2o)%iw(2) = iw3
        itab_ud(iu5)%iw(2) = iw1

     endif

     if (iw == iu3o_iw1) then

        itab_wd(iw2)%iu(2) = iu3o
        itab_wd(iw1)%iu(3) = iu6

        itab_ud(iu3)%im(1) = nest_ud(iu1o)%im
        itab_ud(iu3)%im(2) = nest_ud(iu2o)%im
        itab_ud(iu3)%iw(1) = iw3
        itab_ud(iu3)%iw(2) = iw

        itab_ud(iu3o)%iw(1) = iw2
        itab_ud(iu6)%iw(1) = iw1

     else

        itab_wd(iw2)%iu(2) = iu6
        itab_wd(iw1)%iu(3) = iu3o

        itab_ud(iu3)%im(1) = nest_ud(iu2o)%im
        itab_ud(iu3)%im(2) = nest_ud(iu1o)%im
        itab_ud(iu3)%iw(1) = iw
        itab_ud(iu3)%iw(2) = iw3

        itab_ud(iu3o)%iw(2) = iw1
        itab_ud(iu6)%iw(2) = iw2

     endif

  enddo    ! end of iw loop
  !$omp end do
  !$omp end parallel

  ! Fill itabs loop tables for newly spawned points

  if (iatmgrid) then

     do im = nmd+1,nmd0
        itab_md(im)%imp = im
        call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_vadj, jtm_wstn)
     enddo

     do iu = nud+1,nud0
        itab_ud(iu)%iup = iu
        call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     enddo

     do iw = nwd+1,nwd0
        itab_wd(iw)%iwp = iw
        call wdloopf('f',iw, jtw_grid, jtw_prog, jtw_wadj, jtw_vadj, 0, 0)
     enddo

  endif

  ! Copy new counter values

  nmd = nmd0
  nud = nud0
  nwd = nwd0

  deallocate (ltab_md,ltab_ud,ltab_wd)
  deallocate (nest_ud,nest_wd)
  deallocate (xem_temp,yem_temp,zem_temp)

  if (.not. iatmgrid) deallocate (iwdorig_temp)

  call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

  ! Plot grid lines

  if (.false.) then

     call o_reopnwk()
     call plotback()

     call oplot_set(1)

     do iu = 2,nud
        im1 = itab_ud(iu)%im(1)
        im2 = itab_ud(iu)%im(2)

        call oplot_transform(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
        call oplot_transform(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        call o_frstpt (xq1,yq1)
        call o_vector (xq2,yq2)
     enddo

     call o_frame()
     call o_clswk()

  endif

  write(*,*)
  write(*,*) "After doubling mesh:"
  write(*,'(a,i9)')   ' nmd = ',nmd
  write(*,'(a,i9)')   ' nud = ',nud
  write(*,'(a,i9)')   ' nwd = ',nwd

end subroutine expand_global2


subroutine expand_global3(iatmgrid)

  ! This subroutine adds nested grid regions at the beginning of a simulation.

  use mem_delaunay, only: nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd, alloc_itabsd, &
                          itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          iwdorig, iwdorig_temp

  use mem_ijtabs,   only: mloops, &
                          jtm_grid, jtu_grid, jtw_grid, &
                          jtm_init, jtu_init, jtw_init, &
                          jtm_prog, jtu_prog, jtw_prog, &
                          jtm_wadj, jtu_wadj, jtw_wadj, &
                          jtm_wstn, jtu_wstn, jtw_wstn, &
                          jtm_vadj,           jtw_vadj

  use misc_coms,    only: mdomain
  use consts_coms,  only: pio180, erad, pi1, pi2

  implicit none

  logical, intent(in) :: iatmgrid

  type (itab_md_vars), allocatable :: ltab_md(:)
  type (itab_ud_vars), allocatable :: ltab_ud(:)
  type (itab_wd_vars), allocatable :: ltab_wd(:)

  Type nest_ud3_vars       ! temporary U-pt data structure for spawning nested grids
     integer :: im(2) = 0  ! new M pts attached to this U pt
     integer :: iu(2) = 0  ! new U pts attached to this U pt
  End Type nest_ud3_vars

  Type nest_wd3_vars       ! temporary W-pt data structure for spawning nested grids
     integer :: iu(9) = 0  ! new U pts attached to this W pt
     integer :: iw(8) = 0  ! new W pts attached to this W pt
     integer :: im    = 0  ! new M pnt attached to this W pt
  End Type nest_wd3_vars

  type (nest_ud3_vars), allocatable :: nest_ud(:)
  type (nest_wd3_vars), allocatable :: nest_wd(:)

  integer :: iu, iw, im, im1, im2, im3, iu1o, iu2o, iu3o
  integer :: iu1o_iw1, iu2o_iw1, iu3o_iw1
  integer :: j, ngr, mrlo, mrloo, mrow
  integer :: iun(15)
  integer :: iwn( 8)

  integer :: nmd0,nud0,nwd0,iskip

  real :: expansion, wgt1, wgt2

  real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

  real :: xp1, xp2, xq1, xq2
  real :: yp1, yp2, yq1, yq2

  ! Make duplicate of current grid dimensions

  nmd0 = nmd
  nud0 = nud
  nwd0 = nwd

  write(*,'(/,a)') 'Expanding global mesh'

  ! Allocate temporary tables

  allocate (nest_ud(nud), nest_wd(nwd))  ! Nest relations

  ! Copy ITAB information and M-point coordinates to temporary tables

  call move_alloc(itab_md, ltab_md)
  call move_alloc(itab_ud, ltab_ud)
  call move_alloc(itab_wd, ltab_wd)

  call move_alloc(xemd, xem_temp)
  call move_alloc(yemd, yem_temp)
  call move_alloc(zemd, zem_temp)

  if (.not. iatmgrid) call move_alloc(iwdorig, iwdorig_temp)

  do iw = 2, nwd
     nest_wd(iw)%iu(1) = nud0 + 1
     nest_wd(iw)%iu(2) = nud0 + 2
     nest_wd(iw)%iu(3) = nud0 + 3
     nest_wd(iw)%iu(4) = nud0 + 4
     nest_wd(iw)%iu(5) = nud0 + 5
     nest_wd(iw)%iu(6) = nud0 + 6
     nest_wd(iw)%iu(7) = nud0 + 7
     nest_wd(iw)%iu(8) = nud0 + 8
     nest_wd(iw)%iu(9) = nud0 + 9

     nest_wd(iw)%iw(1) = nwd0 + 1
     nest_wd(iw)%iw(2) = nwd0 + 2
     nest_wd(iw)%iw(3) = nwd0 + 3
     nest_wd(iw)%iw(4) = nwd0 + 4
     nest_wd(iw)%iw(5) = nwd0 + 5
     nest_wd(iw)%iw(6) = nwd0 + 6
     nest_wd(iw)%iw(7) = nwd0 + 7
     nest_wd(iw)%iw(8) = nwd0 + 8

     nest_wd(iw)%im    = nmd0 + 1

     nwd0 = nwd0 + 8
     nud0 = nud0 + 9
     nmd0 = nmd0 + 1
  enddo

  do iu = 2,nud
     nest_ud(iu)%im(1) = nmd0 + 1
     nest_ud(iu)%im(2) = nmd0 + 2

     nest_ud(iu)%iu(1) = nud0 + 1
     nest_ud(iu)%iu(2) = nud0 + 2

     nmd0 = nmd0 + 2
     nud0 = nud0 + 2
  enddo

  ! Allocate main tables to expanded size
  ! Initialize all neighbor indices to zero

  call alloc_itabsd(nmd0,nud0,nwd0)

  if (.not. iatmgrid) allocate(iwdorig(nwd0))

  ! Memory copy to main tables

  if (.not. iatmgrid) then
     iwdorig(1:nwd) = iwdorig_temp
  endif

  xemd(1:nmd) = xem_temp
  yemd(1:nmd) = yem_temp
  zemd(1:nmd) = zem_temp

  do im = 1,nmd
     if (iatmgrid) itab_md(im)%loop(1:mloops) = ltab_md(im)%loop(1:mloops)
     itab_md(im)%imp       = ltab_md(im)%imp
     itab_md(im)%mrlm      = ltab_md(im)%mrlm
     itab_md(im)%mrlm_orig = ltab_md(im)%mrlm_orig
     itab_md(im)%ngr       = ltab_md(im)%ngr
  enddo

  itab_ud(1:nud) = ltab_ud(1:nud)
  itab_wd(1:nwd) = ltab_wd(1:nwd)

  do iu = 2,nud
     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     do j = 1, 2
        im  = nest_ud(iu)%im(j)

        wgt1 = real(j) / 3.
        wgt2 = 1.0 - wgt1

        ! Average coordinates to new M points

        xemd(im) = wgt2 * xemd(im1) + wgt1 * xemd(im2)
        yemd(im) = wgt2 * yemd(im1) + wgt1 * yemd(im2)
        zemd(im) = wgt2 * zemd(im1) + wgt1 * zemd(im2)

        ! If mdomain <= 1, push M point coordinates out to earth radius

        if (mdomain <= 1) then

           expansion = erad / sqrt( xemd(im) ** 2 &
                                  + yemd(im) ** 2 &
                                  + zemd(im) ** 2 )

           xemd(im) = xemd(im) * expansion
           yemd(im) = yemd(im) * expansion
           zemd(im) = zemd(im) * expansion

        endif

        mrlo  = itab_md(im1)%mrlm
        mrloo = itab_md(im1)%mrlm_orig
        ngr   = itab_md(im1)%ngr

        itab_md(im)%mrlm      = mrlo
        itab_md(im)%mrlm_orig = mrloo
        itab_md(im)%ngr       = ngr

     enddo
  enddo

  ! Contruct tables for new fully subdivided triangles

  do iw = 2,nwd

     im  = nest_wd(iw)%im
     im1 = ltab_wd(iw)%im(1)
     im2 = ltab_wd(iw)%im(2)
     im3 = ltab_wd(iw)%im(3)

     xemd(im) = (xemd(im1) + xemd(im2) + xemd(im3)) / 3.
     yemd(im) = (yemd(im1) + yemd(im2) + yemd(im3)) / 3.
     zemd(im) = (zemd(im1) + zemd(im2) + zemd(im3)) / 3.

     ! If mdomain <= 1, push M point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt( xemd(im) ** 2 &
                               + yemd(im) ** 2 &
                               + zemd(im) ** 2 )

        xemd(im) = xemd(im) * expansion
        yemd(im) = yemd(im) * expansion
        zemd(im) = zemd(im) * expansion

     endif

     ! Mapping of original undivided triangles

     iu1o  = ltab_wd(iw)%iu(1)
     iu2o  = ltab_wd(iw)%iu(2)
     iu3o  = ltab_wd(iw)%iu(3)
     mrlo  = ltab_wd(iw)%mrlw
     mrloo = ltab_wd(iw)%mrlw_orig
     ngr   = ltab_wd(iw)%ngr
     mrow  = ltab_wd(iw)%mrow


     iu1o_iw1 = ltab_ud(iu1o)%iw(1)
     iu2o_iw1 = ltab_ud(iu2o)%iw(1)
     iu3o_iw1 = ltab_ud(iu3o)%iw(1)

     ! Mapping of new divided triangles

     iun(1:9) = nest_wd(iw)%iu(1:9)

     iun(10:11) = nest_ud(iu1o)%iu(1:2)
     iun(12:13) = nest_ud(iu2o)%iu(1:2)
     iun(14:15) = nest_ud(iu3o)%iu(1:2)

     iwn(1:8) = nest_wd(iw)%iw(1:8)

     ! Fill tables with new values

     itab_wd(iw)%iu(1) = iun(1)
     itab_wd(iw)%iu(2) = iun(2)
     itab_wd(iw)%iu(3) = iun(3)

     itab_wd(iwn(1))%iu(1) = iun(1)

     itab_wd(iwn(2))%iu(1) = iun(2)
     itab_wd(iwn(2))%iu(2) = iun(14)
     itab_wd(iwn(2))%iu(3) = iun(4)

     itab_wd(iwn(3))%iu(1) = iun(3)
     itab_wd(iwn(3))%iu(2) = iun(7)
     itab_wd(iwn(3))%iu(3) = iun(12)

     itab_wd(iwn(4))%iu(1) = iun(5)

     itab_wd(iwn(5))%iu(1) = iun(4)
     itab_wd(iwn(5))%iu(2) = iun(5)
     itab_wd(iwn(5))%iu(3) = iun(6)

     itab_wd(iwn(6))%iu(1) = iun(10)
     itab_wd(iwn(6))%iu(2) = iun(8)
     itab_wd(iwn(6))%iu(3) = iun(6)

     itab_wd(iwn(7))%iu(1) = iun(7)
     itab_wd(iwn(7))%iu(2) = iun(8)
     itab_wd(iwn(7))%iu(3) = iun(9)

     itab_wd(iwn(8))%iu(1) = iun(9)

     do j = 1, 8
        itab_wd(iwn(j))%mrlw      = mrlo
        itab_wd(iwn(j))%mrlw_orig = mrloo
        itab_wd(iwn(j))%ngr       = ngr
        itab_wd(iwn(j))%mrow      = mrow

        if (.not. iatmgrid) then
           iwdorig(iwn(j)) = iwdorig_temp(iw)
        endif
     enddo

     itab_ud(iu1o   )%im(2) = nest_ud(iu1o)%im(1)
     itab_ud(iun(10))%im(1) = nest_ud(iu1o)%im(1)
     itab_ud(iun(10))%im(2) = nest_ud(iu1o)%im(2)
     itab_ud(iun(11))%im(1) = nest_ud(iu1o)%im(2)
     itab_ud(iun(11))%im(2) = ltab_ud(iu1o)%im(2)

     itab_ud(iu2o   )%im(2) = nest_ud(iu2o)%im(1)
     itab_ud(iun(12))%im(1) = nest_ud(iu2o)%im(1)
     itab_ud(iun(12))%im(2) = nest_ud(iu2o)%im(2)
     itab_ud(iun(13))%im(1) = nest_ud(iu2o)%im(2)
     itab_ud(iun(13))%im(2) = ltab_ud(iu2o)%im(2)

     itab_ud(iu3o   )%im(2) = nest_ud(iu3o)%im(1)
     itab_ud(iun(14))%im(1) = nest_ud(iu3o)%im(1)
     itab_ud(iun(14))%im(2) = nest_ud(iu3o)%im(2)
     itab_ud(iun(15))%im(1) = nest_ud(iu3o)%im(2)
     itab_ud(iun(15))%im(2) = ltab_ud(iu3o)%im(2)

     if (iw == iu1o_iw1) then

        itab_wd(iwn(8))%iu(2) = iu1o
        itab_wd(iwn(4))%iu(3) = iun(11)

        itab_ud(iun(7))%im(2) = nest_wd(iw)%im
        itab_ud(iun(4))%im(1) = nest_wd(iw)%im

        if (iw == iu2o_iw1) then
           itab_ud(iun(1))%im(1) = nest_ud(iu2o)%im(1)
           itab_ud(iun(7))%im(1) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iun(1))%im(1) = nest_ud(iu2o)%im(2)
           itab_ud(iun(7))%im(1) = nest_ud(iu2o)%im(1)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iun(1))%im(2) = nest_ud(iu3o)%im(2)
           itab_ud(iun(4))%im(2) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iun(1))%im(2) = nest_ud(iu3o)%im(1)
           itab_ud(iun(4))%im(2) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iun(1))%iw(1) = iwn(1)
        itab_ud(iun(1))%iw(2) = iw

        itab_ud(iun(7))%iw(1) = iwn(3)
        itab_ud(iun(7))%iw(2) = iwn(7)

        itab_ud(iun(4))%iw(1) = iwn(2)
        itab_ud(iun(4))%iw(2) = iwn(5)

        itab_ud(iu1o   )%iw(1) = iwn(8)
        itab_ud(iun(10))%iw(1) = iwn(6)
        itab_ud(iun(11))%iw(1) = iwn(4)

     else

        itab_wd(iwn(8))%iu(2) = iun(11)
        itab_wd(iwn(4))%iu(3) = iu1o

        itab_ud(iun(7))%im(1) = nest_wd(iw)%im
        itab_ud(iun(4))%im(2) = nest_wd(iw)%im

        if (iw == iu2o_iw1) then
           itab_ud(iun(1))%im(2) = nest_ud(iu2o)%im(1)
           itab_ud(iun(7))%im(2) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iun(1))%im(2) = nest_ud(iu2o)%im(2)
           itab_ud(iun(7))%im(2) = nest_ud(iu2o)%im(1)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iun(1))%im(1) = nest_ud(iu3o)%im(2)
           itab_ud(iun(4))%im(1) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iun(1))%im(1) = nest_ud(iu3o)%im(1)
           itab_ud(iun(4))%im(1) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iun(1))%iw(1) = iw
        itab_ud(iun(1))%iw(2) = iwn(1)

        itab_ud(iun(7))%iw(1) = iwn(7)
        itab_ud(iun(7))%iw(2) = iwn(3)

        itab_ud(iun(4))%iw(1) = iwn(5)
        itab_ud(iun(4))%iw(2) = iwn(2)

        itab_ud(iu1o   )%iw(2) = iwn(4)
        itab_ud(iun(10))%iw(2) = iwn(6)
        itab_ud(iun(11))%iw(2) = iwn(8)

     endif

     if (iw == iu2o_iw1) then

        itab_wd(iwn(1))%iu(2) = iu2o
        itab_wd(iwn(8))%iu(3) = iun(13)

        itab_ud(iun(2))%im(2) = nest_wd(iw)%im
        itab_ud(iun(8))%im(1) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iun(5))%im(2) = nest_ud(iu1o)%im(2)
           itab_ud(iun(8))%im(2) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iun(5))%im(2) = nest_ud(iu1o)%im(1)
           itab_ud(iun(8))%im(2) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iun(2))%im(1) = nest_ud(iu3o)%im(2)
           itab_ud(iun(5))%im(1) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iun(2))%im(1) = nest_ud(iu3o)%im(1)
           itab_ud(iun(5))%im(1) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iun(2))%iw(1) = iwn(2)
        itab_ud(iun(2))%iw(2) = iw

        itab_ud(iun(8))%iw(1) = iwn(6)
        itab_ud(iun(8))%iw(2) = iwn(7)

        itab_ud(iun(5))%iw(1) = iwn(4)
        itab_ud(iun(5))%iw(2) = iwn(5)

        itab_ud(iu2o   )%iw(1) = iwn(1)
        itab_ud(iun(12))%iw(1) = iwn(3)
        itab_ud(iun(13))%iw(1) = iwn(8)

     else

        itab_wd(iwn(1))%iu(2) = iun(13)
        itab_wd(iwn(8))%iu(3) = iu2o

        itab_ud(iun(2))%im(1) = nest_wd(iw)%im
        itab_ud(iun(8))%im(2) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iun(5))%im(1) = nest_ud(iu1o)%im(2)
           itab_ud(iun(8))%im(1) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iun(5))%im(1) = nest_ud(iu1o)%im(1)
           itab_ud(iun(8))%im(1) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iun(2))%im(2) = nest_ud(iu3o)%im(2)
           itab_ud(iun(5))%im(2) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iun(2))%im(2) = nest_ud(iu3o)%im(1)
           itab_ud(iun(5))%im(2) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iun(2))%iw(1) = iw
        itab_ud(iun(2))%iw(2) = iwn(2)

        itab_ud(iun(8))%iw(1) = iwn(7)
        itab_ud(iun(8))%iw(2) = iwn(6)

        itab_ud(iun(5))%iw(1) = iwn(5)
        itab_ud(iun(5))%iw(2) = iwn(4)

        itab_ud(iu2o   )%iw(2) = iwn(8)
        itab_ud(iun(12))%iw(2) = iwn(3)
        itab_ud(iun(13))%iw(2) = iwn(1)
     endif

     if (iw == iu3o_iw1) then

        itab_wd(iwn(4))%iu(2) = iu3o
        itab_wd(iwn(1))%iu(3) = iun(15)

        itab_ud(iun(3))%im(1) = nest_wd(iw)%im
        itab_ud(iun(6))%im(2) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iun(6))%im(1) = nest_ud(iu1o)%im(2)
           itab_ud(iun(9))%im(1) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iun(6))%im(1) = nest_ud(iu1o)%im(1)
           itab_ud(iun(9))%im(1) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu2o_iw1) then
           itab_ud(iun(3))%im(2) = nest_ud(iu2o)%im(1)
           itab_ud(iun(9))%im(2) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iun(3))%im(2) = nest_ud(iu2o)%im(2)
           itab_ud(iun(9))%im(2) = nest_ud(iu2o)%im(1)
        endif

        itab_ud(iun(3))%iw(1) = iwn(3)
        itab_ud(iun(3))%iw(2) = iw

        itab_ud(iun(6))%iw(1) = iwn(6)
        itab_ud(iun(6))%iw(2) = iwn(5)

        itab_ud(iun(9))%iw(1) = iwn(8)
        itab_ud(iun(9))%iw(2) = iwn(7)

        itab_ud(iu3o   )%iw(1) = iwn(4)
        itab_ud(iun(14))%iw(1) = iwn(2)
        itab_ud(iun(15))%iw(1) = iwn(1)

     else

        itab_wd(iwn(4))%iu(2) = iun(15)
        itab_wd(iwn(1))%iu(3) = iu3o

        itab_ud(iun(3))%im(2) = nest_wd(iw)%im
        itab_ud(iun(6))%im(1) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iun(6))%im(2) = nest_ud(iu1o)%im(2)
           itab_ud(iun(9))%im(2) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iun(6))%im(2) = nest_ud(iu1o)%im(1)
           itab_ud(iun(9))%im(2) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu2o_iw1) then
           itab_ud(iun(3))%im(1) = nest_ud(iu2o)%im(1)
           itab_ud(iun(9))%im(1) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iun(3))%im(1) = nest_ud(iu2o)%im(2)
           itab_ud(iun(9))%im(1) = nest_ud(iu2o)%im(1)
        endif

        itab_ud(iun(3))%iw(1) = iw
        itab_ud(iun(3))%iw(2) = iwn(3)

        itab_ud(iun(6))%iw(1) = iwn(5)
        itab_ud(iun(6))%iw(2) = iwn(6)

        itab_ud(iun(9))%iw(1) = iwn(7)
        itab_ud(iun(9))%iw(2) = iwn(8)

        itab_ud(iu3o   )%iw(2) = iwn(1)
        itab_ud(iun(14))%iw(2) = iwn(2)
        itab_ud(iun(15))%iw(2) = iwn(4)

     endif

  enddo    ! end of iw loop

  ! Fill itabs loop tables for newly spawned points

  if (iatmgrid) then

     do im = nmd+1,nmd0
        itab_md(im)%imp = im
        call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_wstn, 0)
     enddo

     do iu = nud+1,nud0
        itab_ud(iu)%iup = iu
        call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     enddo

     do iw = nwd+1,nwd0
        itab_wd(iw)%iwp = iw
        call wdloopf('f',iw, jtw_grid, jtw_vadj, 0, 0, 0, 0)
     enddo

  endif

  ! Copy new counter values

  nmd = nmd0
  nud = nud0
  nwd = nwd0

  deallocate (ltab_md,ltab_ud,ltab_wd)
  deallocate (nest_ud,nest_wd)
  deallocate (xem_temp,yem_temp,zem_temp)

  if (.not. iatmgrid) deallocate (iwdorig_temp)

  call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

  ! Plot grid lines

  if (.false.) then

     call o_reopnwk()
     call plotback()

     call oplot_set(1)

     do iu = 2,nud
        im1 = itab_ud(iu)%im(1)
        im2 = itab_ud(iu)%im(2)

        call oplot_transform(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
        call oplot_transform(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        call o_frstpt (xq1,yq1)
        call o_vector (xq2,yq2)
     enddo

     call o_frame()
     call o_clswk()

  endif

  write(*,*)
  write(*,*) "After tripling mesh"
  write(*,'(a,i9)')   ' nmd = ',nmd
  write(*,'(a,i9)')   ' nud = ',nud
  write(*,'(a,i9)')   ' nwd = ',nwd

end subroutine expand_global3
