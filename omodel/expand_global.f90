subroutine expand_delaunay_mesh(ires_factor, iatm)

  use mem_delaunay, only: nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd
  use misc_coms,    only: runtype, nxp

  implicit none

  integer, intent(in) :: ires_factor
  logical, intent(in) :: iatm

  integer :: nn, n2, n3, i, nxp00

  n2 = 0
  n3 = 0

  nn = ires_factor

  if (iatm) then
     nxp00 = nxp / ires_factor
  else
     nxp00 = nxp
  endif

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
        write(*,*)
        write(*,*) "Before mesh expansion:"
     else
        write(*,*)
        write(*,*) "Before mesh expansion for sfc grid:"
     endif

     write(*,'(A,I9)') " nmd = ", nmd
     write(*,'(A,I9)') " nud = ", nud
     write(*,'(A,I9)') " nwd = ", nwd

     do i = 1, n3
        call expand_global3(iatm)
        nxp00 = nxp00 * 3
     enddo

     do i = 1, n2
        call expand_global2(iatm)
        nxp00 = nxp00 * 2
     enddo

  endif

end subroutine expand_delaunay_mesh



subroutine expand_global2(iatmgrid)

  ! This subroutine adds nested grid regions at the beginning of a simulation.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          nest_ud_vars, nest_wd_vars, alloc_itabsd, &
                          nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd, impent

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

  integer :: iu,iw,im,iw1,iw2,im1,im2,iun,imn,iwn,j
  integer :: iu1,iu2,iu3,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1
  integer :: iu4,iu5,iu6,iw3,ngr,mrlo,mrloo,mrow,iu1n,iu2n,iu3n

  integer :: nmd0,nud0,nwd0,iskip

  real :: expansion

  real,    allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)
  logical, allocatable :: iwdiv(:), iudiv(:)
  logical, allocatable :: ismnest(:), isunest(:), iswnest(:)

  integer :: imnew(nmd), imnext
  integer :: iwnew(nwd), iwnext
  integer :: iunew(nud), iunext

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

  allocate(iwdiv(nwd))
  iwdiv = .false.

  iunext = 2
  iunew(1) = 1

  do iu = 2, nud
     iunew(iu) = iunext

     iunext = iunext + 1
     nest_ud(iu)%iu = iunext

     do j = 1, 2
        iw = ltab_ud(iu)%iw(j)
        if (.not. iwdiv(iw)) then
           iwdiv(iw) = .true.

           iunext = iunext + 1
           nest_wd(iw)%iu(1) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(2) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(3) = iunext
        endif
     enddo

     iunext = iunext + 1
  enddo

  nud0 = iunext - 1
  deallocate(iwdiv)

  iwnext = 2
  iwnew(1) = 1

  do iw = 2, nwd
     iwnew(iw) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(1) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(2) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(3) = iwnext

     iwnext = iwnext + 1
  enddo

  nwd0 = iwnext - 1

  allocate(iudiv(nud))
  iudiv = .false.

  imnext = 2
  imnew(1) = 1

  do im = 2, nmd
     imnew(im) = imnext

     do j = 1, ltab_md(im)%npoly
        iu = ltab_md(im)%iu(j)
        if (.not. iudiv(iu)) then
           iudiv(iu) = .true.

           imnext = imnext + 1
           nest_ud(iu)%im = imnext
        endif
     enddo

     imnext = imnext + 1
  enddo

  nmd0 = imnext - 1
  deallocate(iudiv)

  do im = 1, 12
     impent(im) = imnew(impent(im))
  enddo

  ! Allocate main tables to expanded size
  ! Initialize all neighbor indices to zero

  call alloc_itabsd(nmd0,nud0,nwd0)

  if (iatmgrid) then
     allocate(ismnest(nmd0)) ; ismnest = .true.
     allocate(isunest(nud0)) ; isunest = .true.
     allocate(iswnest(nwd0)) ; iswnest = .true.
  endif

  ! Memory copy to main tables

  do im = 2, nmd
     imn = imnew(im)

     xemd(imn) = xem_temp(im)
     yemd(imn) = yem_temp(im)
     zemd(imn) = zem_temp(im)

     itab_md(imn)%mrlm      = ltab_md(im)%mrlm
     itab_md(imn)%mrlm_orig = ltab_md(im)%mrlm_orig
     itab_md(imn)%ngr       = ltab_md(im)%ngr
     itab_md(imn)%npoly     = ltab_md(im)%npoly

     if (iatmgrid) then
        itab_md(imn)%loop(1:mloops) = ltab_md(im)%loop(1:mloops)
        itab_md(imn)%imp     = imnew( ltab_md(im)%imp )
        ismnest(imn)         = .false.
        itab_md(imn)%im_orig = im
     else
        itab_md(imn)%im_orig = ltab_md(im)%im_orig
     endif

     do j = 1, 7
        itab_md(imn)%im(j) = imnew( ltab_md(im)%im(j) )
        itab_md(imn)%iu(j) = iunew( ltab_md(im)%iu(j) )
        itab_md(imn)%iw(j) = iwnew( ltab_md(im)%iw(j) )
     enddo
  enddo

  deallocate (xem_temp,yem_temp,zem_temp)

  do iu = 2, nud
     iun = iunew(iu)

     itab_ud(iun)%mrlu = ltab_ud(iu)%mrlu

     if (iatmgrid) then
        itab_ud(iun)%loop(1:mloops) = ltab_ud(iu)%loop(1:mloops)
        itab_ud(iun)%iup     = iunew( ltab_ud(iu)%iup )
        isunest(iun)         = .false.
     endif

     do j = 1, 2
        itab_ud(iun)%im(j) = imnew( ltab_ud(iu)%im(j) )
     enddo
     do j = 1, 6
        itab_ud(iun)%iw(j) = iwnew( ltab_ud(iu)%iw(j) )
     enddo
     do j = 1, 12
        itab_ud(iun)%iu(j) = iunew( ltab_ud(iu)%iu(j) )
     enddo

     im  = nest_ud(iu )%im
     im1 = itab_ud(iun)%im(1)
     im2 = itab_ud(iun)%im(2)

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

  ! Contruct tables for new fully subdivided triangles

  do iw = 2, nwd
     iwn = iwnew(iw)

     itab_wd(iwn)%npoly     = ltab_wd(iw)%npoly
     itab_wd(iwn)%mrow      = ltab_wd(iw)%mrow
     itab_wd(iwn)%ngr       = ltab_wd(iw)%ngr
     itab_wd(iwn)%mrlw      = ltab_wd(iw)%mrlw
     itab_wd(iwn)%mrlw_orig = ltab_wd(iw)%mrlw_orig

     do j = 1, 3
        itab_wd(iwn)%im(j) = imnew( ltab_wd(iw)%im(j) )
        itab_wd(iwn)%iu(j) = iunew( ltab_wd(iw)%iu(j) )
     enddo
     do j = 1, 9
        itab_wd(iwn)%iw(j) = iwnew( ltab_wd(iw)%iw(j) )
     enddo

     if (iatmgrid) then
        itab_wd(iwn)%loop(1:mloops) = ltab_wd(iw)%loop(1:mloops)
        itab_wd(iwn)%iwp     = iwnew( ltab_wd(iw)%iwp )
        iswnest(iwn)         = .false.
     endif

     ! Mapping of original undivided triangles

     iu1o  = ltab_wd(iw)%iu(1)
     iu2o  = ltab_wd(iw)%iu(2)
     iu3o  = ltab_wd(iw)%iu(3)

     iu1n  = iunew(iu1o)
     iu2n  = iunew(iu2o)
     iu3n  = iunew(iu3o)

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

     itab_wd(iwn)%iu(1) = iu1
     itab_wd(iwn)%iu(2) = iu2
     itab_wd(iwn)%iu(3) = iu3

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

     itab_ud(iu1n)%im(2) = nest_ud(iu1o)%im
     itab_ud(iu4 )%im(1) = nest_ud(iu1o)%im
     itab_ud(iu4 )%im(2) = imnew( ltab_ud(iu1o)%im(2) )

     itab_ud(iu2n)%im(2) = nest_ud(iu2o)%im
     itab_ud(iu5 )%im(1) = nest_ud(iu2o)%im
     itab_ud(iu5 )%im(2) = imnew( ltab_ud(iu2o)%im(2) )

     itab_ud(iu3n)%im(2) = nest_ud(iu3o)%im
     itab_ud(iu6 )%im(1) = nest_ud(iu3o)%im
     itab_ud(iu6 )%im(2) = imnew( ltab_ud(iu3o)%im(2) )

     if (iw == iu1o_iw1) then

        itab_wd(iw3)%iu(2) = iu1n
        itab_wd(iw2)%iu(3) = iu4

        itab_ud(iu1)%im(1) = nest_ud(iu2o)%im
        itab_ud(iu1)%im(2) = nest_ud(iu3o)%im
        itab_ud(iu1)%iw(1) = iw1
        itab_ud(iu1)%iw(2) = iwn

        itab_ud(iu1n)%iw(1) = iw3
        itab_ud(iu4 )%iw(1) = iw2

     else

        itab_wd(iw3)%iu(2) = iu4
        itab_wd(iw2)%iu(3) = iu1n

        itab_ud(iu1)%im(1) = nest_ud(iu3o)%im
        itab_ud(iu1)%im(2) = nest_ud(iu2o)%im
        itab_ud(iu1)%iw(1) = iwn
        itab_ud(iu1)%iw(2) = iw1

        itab_ud(iu1n)%iw(2) = iw2
        itab_ud(iu4 )%iw(2) = iw3

     endif

     if (iw == iu2o_iw1) then

        itab_wd(iw1)%iu(2) = iu2n
        itab_wd(iw3)%iu(3) = iu5

        itab_ud(iu2)%im(1) = nest_ud(iu3o)%im
        itab_ud(iu2)%im(2) = nest_ud(iu1o)%im
        itab_ud(iu2)%iw(1) = iw2
        itab_ud(iu2)%iw(2) = iwn

        itab_ud(iu2n)%iw(1) = iw1
        itab_ud(iu5 )%iw(1) = iw3

     else

        itab_wd(iw1)%iu(2) = iu5
        itab_wd(iw3)%iu(3) = iu2n

        itab_ud(iu2)%im(1) = nest_ud(iu1o)%im
        itab_ud(iu2)%im(2) = nest_ud(iu3o)%im
        itab_ud(iu2)%iw(1) = iwn
        itab_ud(iu2)%iw(2) = iw2

        itab_ud(iu2n)%iw(2) = iw3
        itab_ud(iu5 )%iw(2) = iw1

     endif

     if (iw == iu3o_iw1) then

        itab_wd(iw2)%iu(2) = iu3n
        itab_wd(iw1)%iu(3) = iu6

        itab_ud(iu3)%im(1) = nest_ud(iu1o)%im
        itab_ud(iu3)%im(2) = nest_ud(iu2o)%im
        itab_ud(iu3)%iw(1) = iw3
        itab_ud(iu3)%iw(2) = iwn

        itab_ud(iu3n)%iw(1) = iw2
        itab_ud(iu6 )%iw(1) = iw1

     else

        itab_wd(iw2)%iu(2) = iu6
        itab_wd(iw1)%iu(3) = iu3n

        itab_ud(iu3)%im(1) = nest_ud(iu2o)%im
        itab_ud(iu3)%im(2) = nest_ud(iu1o)%im
        itab_ud(iu3)%iw(1) = iwn
        itab_ud(iu3)%iw(2) = iw3

        itab_ud(iu3n)%iw(2) = iw1
        itab_ud(iu6 )%iw(2) = iw2

     endif

  enddo    ! end of iw loop

  ! Fill itabs loop tables for newly spawned points

  if (iatmgrid) then

     do im = 2, nmd0
        if (ismnest(im)) then
           itab_md(im)%imp = im
           call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_vadj, jtm_wstn)
        endif
     enddo
     deallocate (ismnest)

     do iu = 2, nud0
        if (isunest(iu)) then
           itab_ud(iu)%iup = iu
           call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
        endif
     enddo
     deallocate (isunest)

     do iw = 2, nwd0
        if (iswnest(iw)) then
           itab_wd(iw)%iwp = iw
           call wdloopf('f',iw, jtw_grid, jtw_prog, jtw_wadj, jtw_vadj, 0, 0)
        endif
     enddo
     deallocate (iswnest)

  endif

  ! Copy new counter values

  nmd = nmd0
  nud = nud0
  nwd = nwd0

  deallocate (ltab_md,ltab_ud,ltab_wd)
  deallocate (nest_ud,nest_wd)

  call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

  ! Plot grid lines

  if (.false.) then

     call o_reopnwk()
     call plotback()

     call oplot_set(1)

     do iu = 2,nud
        im1 = itab_ud(iu)%im(1)
        im2 = itab_ud(iu)%im(2)

        call oplot_transform_xyz(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
        call oplot_transform_xyz(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

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

  use mem_delaunay, only: nmd, nud, nwd, xemd, yemd, zemd, impent, &
                          itab_md, itab_ud, itab_wd, alloc_itabsd, &
                          itab_md_vars, itab_ud_vars, itab_wd_vars

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
  integer :: iu1o_iw1, iu2o_iw1, iu3o_iw1, iu1n, iu2n, iu3n
  integer :: j, ngr, mrlo, mrloo, mrow, imn, iun, iwn
  integer :: iux(15)
  integer :: iwx( 8)

  integer :: nmd0,nud0,nwd0,iskip

  real :: expansion, wgt1, wgt2

  real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)
  logical, allocatable :: iwdiv(:), iudiv(:)
  logical, allocatable :: ismnest(:), isunest(:), iswnest(:)

  integer :: imnew(nmd), imnext
  integer :: iwnew(nwd), iwnext
  integer :: iunew(nud), iunext

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

  allocate(iwdiv(nwd))
  iwdiv = .false.

  iunext = 2
  iunew(1) = 1

  do iu = 2, nud
     iunew(iu) = iunext

     iunext = iunext + 1
     nest_ud(iu)%iu(1) = iunext

     iunext = iunext + 1
     nest_ud(iu)%iu(2) = iunext

     do j = 1, 2
        iw = ltab_ud(iu)%iw(j)
        if (.not. iwdiv(iw)) then
           iwdiv(iw) = .true.

           iunext = iunext + 1
           nest_wd(iw)%iu(1) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(2) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(3) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(4) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(5) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(6) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(7) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(8) = iunext

           iunext = iunext + 1
           nest_wd(iw)%iu(9) = iunext
        endif
     enddo

     iunext = iunext + 1
  enddo

  nud0 = iunext - 1
  deallocate(iwdiv)

  iwnext = 2
  iwnew(1) = 1

  do iw = 2, nwd
     iwnew(iw) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(1) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(2) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(3) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(4) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(5) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(6) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(7) = iwnext

     iwnext = iwnext + 1
     nest_wd(iw)%iw(8) = iwnext

     iwnext = iwnext + 1
  enddo

  nwd0 = iwnext - 1

  allocate(iudiv(nud))
  iudiv = .false.

  allocate(iwdiv(nwd))
  iwdiv = .false.

  imnext = 2
  imnew(1) = 1

  do im = 2, nmd
     imnew(im) = imnext

     do j = 1, ltab_md(im)%npoly

        iu = ltab_md(im)%iu(j)
        if (.not. iudiv(iu)) then
           iudiv(iu) = .true.

           imnext = imnext + 1
           nest_ud(iu)%im(1) = imnext

           imnext = imnext + 1
           nest_ud(iu)%im(2) = imnext
        endif

        iw = ltab_md(im)%iw(j)
        if (.not. iwdiv(iw)) then
           iwdiv(iw) = .true.

           imnext = imnext + 1
           nest_wd(iw)%im = imnext
        endif

     enddo

     imnext = imnext + 1
  enddo

  nmd0 = imnext - 1
  deallocate(iwdiv)
  deallocate(iudiv)

  do im = 1, 12
     impent(im) = imnew(impent(im))
  enddo

  ! Allocate main tables to expanded size
  ! Initialize all neighbor indices to zero

  call alloc_itabsd(nmd0,nud0,nwd0)

  if (iatmgrid) then
     allocate(ismnest(nmd0)) ; ismnest = .true.
     allocate(isunest(nud0)) ; isunest = .true.
     allocate(iswnest(nwd0)) ; iswnest = .true.
  endif

  ! Memory copy to main tables

  do im = 2, nmd
     imn = imnew(im)

     xemd(imn) = xem_temp(im)
     yemd(imn) = yem_temp(im)
     zemd(imn) = zem_temp(im)

     itab_md(imn)%mrlm      = ltab_md(im)%mrlm
     itab_md(imn)%mrlm_orig = ltab_md(im)%mrlm_orig
     itab_md(imn)%ngr       = ltab_md(im)%ngr
     itab_md(imn)%npoly     = ltab_md(im)%npoly

     if (iatmgrid) then
        itab_md(imn)%loop(1:mloops) = ltab_md(im)%loop(1:mloops)
        itab_md(imn)%imp     = imnew( ltab_md(im)%imp )
        ismnest(imn)         = .false.
     endif

     do j = 1, 7
        itab_md(imn)%im(j) = imnew( ltab_md(im)%im(j) )
        itab_md(imn)%iu(j) = iunew( ltab_md(im)%iu(j) )
        itab_md(imn)%iw(j) = iwnew( ltab_md(im)%iw(j) )
     enddo
  enddo

  deallocate (xem_temp,yem_temp,zem_temp)

  do iu = 2, nud
     iun = iunew(iu)

     itab_ud(iun)%mrlu = ltab_ud(iu)%mrlu

     if (iatmgrid) then
        itab_ud(iun)%loop(1:mloops) = ltab_ud(iu)%loop(1:mloops)
        itab_ud(iun)%iup     = iunew( ltab_ud(iu)%iup )
        isunest(iun)         = .false.
     endif

     do j = 1, 2
        itab_ud(iun)%im(j) = imnew( ltab_ud(iu)%im(j) )
     enddo
     do j = 1, 6
        itab_ud(iun)%iw(j) = iwnew( ltab_ud(iu)%iw(j) )
     enddo
     do j = 1, 12
        itab_ud(iun)%iu(j) = iunew( ltab_ud(iu)%iu(j) )
     enddo

     im1 = itab_ud(iun)%im(1)
     im2 = itab_ud(iun)%im(2)

     do j = 1, 2
        im = nest_ud(iu)%im(j)

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

  do iw = 2, nwd
     iwn = iwnew(iw)

     itab_wd(iwn)%npoly     = ltab_wd(iw)%npoly
     itab_wd(iwn)%mrow      = ltab_wd(iw)%mrow
     itab_wd(iwn)%ngr       = ltab_wd(iw)%ngr
     itab_wd(iwn)%mrlw      = ltab_wd(iw)%mrlw
     itab_wd(iwn)%mrlw_orig = ltab_wd(iw)%mrlw_orig

     do j = 1, 3
        itab_wd(iwn)%im(j) = imnew( ltab_wd(iw)%im(j) )
        itab_wd(iwn)%iu(j) = iunew( ltab_wd(iw)%iu(j) )
     enddo
     do j = 1, 9
        itab_wd(iwn)%iw(j) = iwnew( ltab_wd(iw)%iw(j) )
     enddo

     if (iatmgrid) then
        itab_wd(iwn)%loop(1:mloops) = ltab_wd(iw)%loop(1:mloops)
        itab_wd(iwn)%iwp     = iwnew( ltab_wd(iw)%iwp )
        iswnest(iwn)         = .false.
     endif

     ! New M point

     im  = nest_wd(iw)%im
     im1 = itab_wd(iwn)%im(1)
     im2 = itab_wd(iwn)%im(2)
     im3 = itab_wd(iwn)%im(3)

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

     mrlo  = itab_md(im1)%mrlm
     mrloo = itab_md(im1)%mrlm_orig
     ngr   = itab_md(im1)%ngr

     itab_md(im)%mrlm      = mrlo
     itab_md(im)%mrlm_orig = mrloo
     itab_md(im)%ngr       = ngr

     ! Mapping of original undivided triangles

     iu1o  = ltab_wd(iw)%iu(1)
     iu2o  = ltab_wd(iw)%iu(2)
     iu3o  = ltab_wd(iw)%iu(3)

     iu1n  = iunew(iu1o)
     iu2n  = iunew(iu2o)
     iu3n  = iunew(iu3o)

     mrlo  = ltab_wd(iw)%mrlw
     mrloo = ltab_wd(iw)%mrlw_orig
     ngr   = ltab_wd(iw)%ngr
     mrow  = ltab_wd(iw)%mrow

     iu1o_iw1 = ltab_ud(iu1o)%iw(1)
     iu2o_iw1 = ltab_ud(iu2o)%iw(1)
     iu3o_iw1 = ltab_ud(iu3o)%iw(1)

     ! Mapping of new divided triangles

     iux(1:9) = nest_wd(iw)%iu(1:9)

     iux(10:11) = nest_ud(iu1o)%iu(1:2)
     iux(12:13) = nest_ud(iu2o)%iu(1:2)
     iux(14:15) = nest_ud(iu3o)%iu(1:2)

     iwx(1:8) = nest_wd(iw)%iw(1:8)

     ! Fill tables with new values

     itab_wd(iwn)%iu(1) = iux(1)
     itab_wd(iwn)%iu(2) = iux(2)
     itab_wd(iwn)%iu(3) = iux(3)

     itab_wd(iwx(1))%iu(1) = iux(1)

     itab_wd(iwx(2))%iu(1) = iux(2)
     itab_wd(iwx(2))%iu(2) = iux(14)
     itab_wd(iwx(2))%iu(3) = iux(4)

     itab_wd(iwx(3))%iu(1) = iux(3)
     itab_wd(iwx(3))%iu(2) = iux(7)
     itab_wd(iwx(3))%iu(3) = iux(12)

     itab_wd(iwx(4))%iu(1) = iux(5)

     itab_wd(iwx(5))%iu(1) = iux(4)
     itab_wd(iwx(5))%iu(2) = iux(5)
     itab_wd(iwx(5))%iu(3) = iux(6)

     itab_wd(iwx(6))%iu(1) = iux(10)
     itab_wd(iwx(6))%iu(2) = iux(8)
     itab_wd(iwx(6))%iu(3) = iux(6)

     itab_wd(iwx(7))%iu(1) = iux(7)
     itab_wd(iwx(7))%iu(2) = iux(8)
     itab_wd(iwx(7))%iu(3) = iux(9)

     itab_wd(iwx(8))%iu(1) = iux(9)

     do j = 1, 8
        itab_wd(iwx(j))%mrlw      = mrlo
        itab_wd(iwx(j))%mrlw_orig = mrloo
        itab_wd(iwx(j))%ngr       = ngr
        itab_wd(iwx(j))%mrow      = mrow
     enddo

     itab_ud(iu1n   )%im(2) = nest_ud(iu1o)%im(1)
     itab_ud(iux(10))%im(1) = nest_ud(iu1o)%im(1)
     itab_ud(iux(10))%im(2) = nest_ud(iu1o)%im(2)
     itab_ud(iux(11))%im(1) = nest_ud(iu1o)%im(2)
     itab_ud(iux(11))%im(2) = imnew( ltab_ud(iu1o)%im(2) )

     itab_ud(iu2n   )%im(2) = nest_ud(iu2o)%im(1)
     itab_ud(iux(12))%im(1) = nest_ud(iu2o)%im(1)
     itab_ud(iux(12))%im(2) = nest_ud(iu2o)%im(2)
     itab_ud(iux(13))%im(1) = nest_ud(iu2o)%im(2)
     itab_ud(iux(13))%im(2) = imnew( ltab_ud(iu2o)%im(2) )

     itab_ud(iu3n   )%im(2) = nest_ud(iu3o)%im(1)
     itab_ud(iux(14))%im(1) = nest_ud(iu3o)%im(1)
     itab_ud(iux(14))%im(2) = nest_ud(iu3o)%im(2)
     itab_ud(iux(15))%im(1) = nest_ud(iu3o)%im(2)
     itab_ud(iux(15))%im(2) = imnew( ltab_ud(iu3o)%im(2) )

     if (iw == iu1o_iw1) then

        itab_wd(iwx(8))%iu(2) = iu1n
        itab_wd(iwx(4))%iu(3) = iux(11)

        itab_ud(iux(7))%im(2) = nest_wd(iw)%im
        itab_ud(iux(4))%im(1) = nest_wd(iw)%im

        if (iw == iu2o_iw1) then
           itab_ud(iux(1))%im(1) = nest_ud(iu2o)%im(1)
           itab_ud(iux(7))%im(1) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iux(1))%im(1) = nest_ud(iu2o)%im(2)
           itab_ud(iux(7))%im(1) = nest_ud(iu2o)%im(1)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iux(1))%im(2) = nest_ud(iu3o)%im(2)
           itab_ud(iux(4))%im(2) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iux(1))%im(2) = nest_ud(iu3o)%im(1)
           itab_ud(iux(4))%im(2) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iux(1))%iw(1) = iwx(1)
        itab_ud(iux(1))%iw(2) = iwn

        itab_ud(iux(7))%iw(1) = iwx(3)
        itab_ud(iux(7))%iw(2) = iwx(7)

        itab_ud(iux(4))%iw(1) = iwx(2)
        itab_ud(iux(4))%iw(2) = iwx(5)

        itab_ud(iu1n   )%iw(1) = iwx(8)
        itab_ud(iux(10))%iw(1) = iwx(6)
        itab_ud(iux(11))%iw(1) = iwx(4)

     else

        itab_wd(iwx(8))%iu(2) = iux(11)
        itab_wd(iwx(4))%iu(3) = iu1n

        itab_ud(iux(7))%im(1) = nest_wd(iw)%im
        itab_ud(iux(4))%im(2) = nest_wd(iw)%im

        if (iw == iu2o_iw1) then
           itab_ud(iux(1))%im(2) = nest_ud(iu2o)%im(1)
           itab_ud(iux(7))%im(2) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iux(1))%im(2) = nest_ud(iu2o)%im(2)
           itab_ud(iux(7))%im(2) = nest_ud(iu2o)%im(1)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iux(1))%im(1) = nest_ud(iu3o)%im(2)
           itab_ud(iux(4))%im(1) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iux(1))%im(1) = nest_ud(iu3o)%im(1)
           itab_ud(iux(4))%im(1) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iux(1))%iw(1) = iwn
        itab_ud(iux(1))%iw(2) = iwx(1)

        itab_ud(iux(7))%iw(1) = iwx(7)
        itab_ud(iux(7))%iw(2) = iwx(3)

        itab_ud(iux(4))%iw(1) = iwx(5)
        itab_ud(iux(4))%iw(2) = iwx(2)

        itab_ud(iu1n   )%iw(2) = iwx(4)
        itab_ud(iux(10))%iw(2) = iwx(6)
        itab_ud(iux(11))%iw(2) = iwx(8)

     endif

     if (iw == iu2o_iw1) then

        itab_wd(iwx(1))%iu(2) = iu2n
        itab_wd(iwx(8))%iu(3) = iux(13)

        itab_ud(iux(2))%im(2) = nest_wd(iw)%im
        itab_ud(iux(8))%im(1) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iux(5))%im(2) = nest_ud(iu1o)%im(2)
           itab_ud(iux(8))%im(2) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iux(5))%im(2) = nest_ud(iu1o)%im(1)
           itab_ud(iux(8))%im(2) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iux(2))%im(1) = nest_ud(iu3o)%im(2)
           itab_ud(iux(5))%im(1) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iux(2))%im(1) = nest_ud(iu3o)%im(1)
           itab_ud(iux(5))%im(1) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iux(2))%iw(1) = iwx(2)
        itab_ud(iux(2))%iw(2) = iwn

        itab_ud(iux(8))%iw(1) = iwx(6)
        itab_ud(iux(8))%iw(2) = iwx(7)

        itab_ud(iux(5))%iw(1) = iwx(4)
        itab_ud(iux(5))%iw(2) = iwx(5)

        itab_ud(iu2n   )%iw(1) = iwx(1)
        itab_ud(iux(12))%iw(1) = iwx(3)
        itab_ud(iux(13))%iw(1) = iwx(8)

     else

        itab_wd(iwx(1))%iu(2) = iux(13)
        itab_wd(iwx(8))%iu(3) = iu2n

        itab_ud(iux(2))%im(1) = nest_wd(iw)%im
        itab_ud(iux(8))%im(2) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iux(5))%im(1) = nest_ud(iu1o)%im(2)
           itab_ud(iux(8))%im(1) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iux(5))%im(1) = nest_ud(iu1o)%im(1)
           itab_ud(iux(8))%im(1) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu3o_iw1) then
           itab_ud(iux(2))%im(2) = nest_ud(iu3o)%im(2)
           itab_ud(iux(5))%im(2) = nest_ud(iu3o)%im(1)
        else
           itab_ud(iux(2))%im(2) = nest_ud(iu3o)%im(1)
           itab_ud(iux(5))%im(2) = nest_ud(iu3o)%im(2)
        endif

        itab_ud(iux(2))%iw(1) = iwn
        itab_ud(iux(2))%iw(2) = iwx(2)

        itab_ud(iux(8))%iw(1) = iwx(7)
        itab_ud(iux(8))%iw(2) = iwx(6)

        itab_ud(iux(5))%iw(1) = iwx(5)
        itab_ud(iux(5))%iw(2) = iwx(4)

        itab_ud(iu2n   )%iw(2) = iwx(8)
        itab_ud(iux(12))%iw(2) = iwx(3)
        itab_ud(iux(13))%iw(2) = iwx(1)
     endif

     if (iw == iu3o_iw1) then

        itab_wd(iwx(4))%iu(2) = iu3n
        itab_wd(iwx(1))%iu(3) = iux(15)

        itab_ud(iux(3))%im(1) = nest_wd(iw)%im
        itab_ud(iux(6))%im(2) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iux(6))%im(1) = nest_ud(iu1o)%im(2)
           itab_ud(iux(9))%im(1) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iux(6))%im(1) = nest_ud(iu1o)%im(1)
           itab_ud(iux(9))%im(1) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu2o_iw1) then
           itab_ud(iux(3))%im(2) = nest_ud(iu2o)%im(1)
           itab_ud(iux(9))%im(2) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iux(3))%im(2) = nest_ud(iu2o)%im(2)
           itab_ud(iux(9))%im(2) = nest_ud(iu2o)%im(1)
        endif

        itab_ud(iux(3))%iw(1) = iwx(3)
        itab_ud(iux(3))%iw(2) = iwn

        itab_ud(iux(6))%iw(1) = iwx(6)
        itab_ud(iux(6))%iw(2) = iwx(5)

        itab_ud(iux(9))%iw(1) = iwx(8)
        itab_ud(iux(9))%iw(2) = iwx(7)

        itab_ud(iu3n   )%iw(1) = iwx(4)
        itab_ud(iux(14))%iw(1) = iwx(2)
        itab_ud(iux(15))%iw(1) = iwx(1)

     else

        itab_wd(iwx(4))%iu(2) = iux(15)
        itab_wd(iwx(1))%iu(3) = iu3n

        itab_ud(iux(3))%im(2) = nest_wd(iw)%im
        itab_ud(iux(6))%im(1) = nest_wd(iw)%im

        if (iw == iu1o_iw1) then
           itab_ud(iux(6))%im(2) = nest_ud(iu1o)%im(2)
           itab_ud(iux(9))%im(2) = nest_ud(iu1o)%im(1)
        else
           itab_ud(iux(6))%im(2) = nest_ud(iu1o)%im(1)
           itab_ud(iux(9))%im(2) = nest_ud(iu1o)%im(2)
        endif

        if (iw == iu2o_iw1) then
           itab_ud(iux(3))%im(1) = nest_ud(iu2o)%im(1)
           itab_ud(iux(9))%im(1) = nest_ud(iu2o)%im(2)
        else
           itab_ud(iux(3))%im(1) = nest_ud(iu2o)%im(2)
           itab_ud(iux(9))%im(1) = nest_ud(iu2o)%im(1)
        endif

        itab_ud(iux(3))%iw(1) = iwn
        itab_ud(iux(3))%iw(2) = iwx(3)

        itab_ud(iux(6))%iw(1) = iwx(5)
        itab_ud(iux(6))%iw(2) = iwx(6)

        itab_ud(iux(9))%iw(1) = iwx(7)
        itab_ud(iux(9))%iw(2) = iwx(8)

        itab_ud(iu3n   )%iw(2) = iwx(1)
        itab_ud(iux(14))%iw(2) = iwx(2)
        itab_ud(iux(15))%iw(2) = iwx(4)

     endif

  enddo    ! end of iw loop

  ! Fill itabs loop tables for newly spawned points

  if (iatmgrid) then

     do im = 2, nmd0
        if (ismnest(im)) then
           itab_md(im)%imp = im
           call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_vadj, jtm_wstn)
        endif
     enddo
     deallocate (ismnest)

     do iu = 2, nud0
        if (isunest(iu)) then
           itab_ud(iu)%iup = iu
           call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
        endif
     enddo
     deallocate (isunest)

     do iw = 2, nwd0
        if (iswnest(iw)) then
           itab_wd(iw)%iwp = iw
           call wdloopf('f',iw, jtw_grid, jtw_prog, jtw_wadj, jtw_vadj, 0, 0)
        endif
     enddo
     deallocate (iswnest)

  endif

  ! Copy new counter values

  nmd = nmd0
  nud = nud0
  nwd = nwd0

  deallocate (ltab_md,ltab_ud,ltab_wd)
  deallocate (nest_ud,nest_wd)

  call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

  ! Plot grid lines

  if (.false.) then

     call o_reopnwk()
     call plotback()

     call oplot_set(1)

     do iu = 2,nud
        im1 = itab_ud(iu)%im(1)
        im2 = itab_ud(iu)%im(2)

        call oplot_transform_xyz(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
        call oplot_transform_xyz(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

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
