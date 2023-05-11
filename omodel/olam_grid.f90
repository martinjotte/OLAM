subroutine gridinit()

  use misc_coms,   only: io6, mdomain, ngrids, nxp, alloc_misc, runtype
  use consts_coms, only: r8
  use mem_ijtabs,  only: mrls, fill_jtabs, itab_w

  use mem_delaunay,only: itab_md, itab_ud, itab_wd, xemd, yemd, zemd, &
                         copy_tri_grid, copyback_tri_grid, nmd, nud, nwd

  use mem_grid,    only: nza, nma, nua, nva, nwa, mma, mva, mwa, zm, &
                         xew, yew, zew, glatw, glonw, arw0, &
                         alloc_grid1, alloc_grid2, alloc_gridz_other

  use mem_nudge,   only: nudflag, nudnxp, nwnud, itab_wnud, alloc_nudge1, &
                         xewnud, yewnud, zewnud

  use mem_sfcg,    only: sfcg, sfcgrid_res_factor, itab_wsfc, &
                         nmsfc, nvsfc, nwsfc, mmsfc, mvsfc, mwsfc

  implicit none

  integer :: npoly
  integer :: j, jmaxneg, jminpos
  integer :: imd,imd1,imd2,iud,iwnud,iwnudn,iw
  integer :: kw_sea, kw_land, iwsfc
  integer :: nn, nxp00

  real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
  real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
  real :: dist_wnud,dist,dist_min,xi,yi
  real :: xin(6),yin(6)

  real(r8), allocatable :: tot_area(:)
  real(r8)              :: area_tot

  ! Vertical ATM grid coordinate setup

  write(io6,'(/,a,/)') 'gridinit calling gridset2'
  call gridset2()

  write(io6,'(a,f8.1)') ' Model top height = ',zm(nza)

  write(io6,'(/,a)') 'gridinit calling alloc_gridz_other'
  call alloc_gridz_other()

  ! Horizontal grid setup

  if (mdomain == 0) then

 ! If using nudging on global domain with independent nudging grid,
 ! generate nudging grid here

     if (nudflag > 0 .and. nudnxp > 0 .and. runtype /= 'MAKEGRID_PLOT') then

        write(io6,'(/,a)') 'gridinit calling icosahedron for nudging grid'

        call icosahedron(nudnxp)  ! global spherical domain; calls 2 allocs

 ! Set dimensions of nudging grid arrays, which are the Voronoi dual-grid
 ! counterpart to the Delaunay triangle mesh just generated in icosahedron
 ! and dimensioned by nwa and nma

        nwnud = nmd

 ! Allocate nudging grid structure arrays, copy temporary grid structure to
 ! these arrays, and deallocate temporary arrays

        call alloc_nudge1(nwnud,2)

 ! Loop over nudging grid W points

        do iwnud = 2,nwnud
           imd = iwnud

 ! Set nudging grid W point npoly value and earth coordinates identical to
 ! triangle mesh MD points

           itab_wnud(iwnud)%npoly = itab_md(imd)%npoly

           xewnud(iwnud) = xemd(imd)
           yewnud(iwnud) = yemd(imd)
           zewnud(iwnud) = zemd(imd)

 ! Loop over IUD neighbors of IMD, and get IMD endpoint indices of each

           do j = 1,itab_md(imd)%npoly
              iud = itab_md(imd)%iu(j)

              imd1 = itab_ud(iud)%im(1)
              imd2 = itab_ud(iud)%im(2)

 ! Assign W neighbors of W points on nudging grid based on IMD

              if (imd1 == imd) then
                 itab_wnud(iwnud)%iwnud(j) = imd2
              else
                 itab_wnud(iwnud)%iwnud(j) = imd1
              endif

           enddo
        enddo

 ! Delaunay triangle mesh that was set up for nuding mesh is no longer required

        deallocate (itab_md, itab_ud, itab_wd)
        deallocate (xemd, yemd, zemd)

        write(io6,'(/,a)') 'gridinit after nudging grid construction'
        write(io6,'(a,i8)')    ' nwnud = ',nwnud

     endif

! Now generate global atmospheric grid

     write(io6,'(/,a)') 'gridinit calling icosahedron'

     mrls = 1  ! default value

     nxp00 = nxp
     nn = 1

     do while( mod(nxp00,3)==0 .and. nxp00/3 >= 21 )
        nxp00 = nxp00 / 3
        nn = nn * 3
     enddo

     do while( mod(nxp00,2)==0 .and. nxp00/2 >= 21 )
        nxp00 = nxp00 / 2
        nn = nn * 2
     enddo

     call icosahedron(nxp00)  ! global spherical domain; calls 2 allocs

     if (nn > 1) then
        call expand_delaunay_mesh(nn, .true.)
     endif

     write(io6,'(/,a)') 'gridinit after icosahedron'
     write(io6,'(a,i0)')    ' nmd = ',nmd
     write(io6,'(a,i0)')    ' nud = ',nud
     write(io6,'(a,i0)')    ' nwd = ',nwd

  elseif (mdomain == 4) then

     write(io6,'(/,a)') 'gridinit calling cartesian_3d'
     call cart4_hex()    ! 3D cartesian channel domain; calls 2 allocs

     write(io6,'(/,a)') 'gridinit after cartesian'
     write(io6,'(a,i0)')    ' nma = ',nma
     write(io6,'(a,i0)')    ' nua = ',nua
     write(io6,'(a,i0)')    ' nwa = ',nwa

  elseif (mdomain == 5) then

     write(io6,'(/,a)') 'gridinit calling cart_hex'
     call cart_hex()        ! 3D hexagonal domain with hexagon cells and
                            ! periodic lateral boundaries; calls 2 allocs

     write(io6,'(/,a)') 'gridinit after cartesian'
     write(io6,'(a,i0)')    ' nmd = ',nmd
     write(io6,'(a,i0)')    ' nud = ',nud
     write(io6,'(a,i0)')    ' nwd = ',nwd

  endif

  if (ngrids > 1 .and. (mdomain == 0 .or. mdomain == 5)) then

     write(io6,'(/,a,i0)') 'gridinit calling spawn_nest for ngrids = ',ngrids

     call spawn_nest(.true.)  ! calls 2 allocs

     !write(io6,'(/,a)') 'gridinit after spawn_nest'
     !write(io6,'(a,i8)')   ' nmd = ',nma
     !write(io6,'(a,i8)')   ' nud = ',nua
     !write(io6,'(a,i8)')   ' nwd = ',nwa

  endif

  if (runtype == 'MAKEGRID_PLOT') then
     nma = nwd
     nua = nud
     nva = nud
     nwa = nmd
  endif

  if (mdomain /= 4 .and. runtype /= 'MAKEREGRID') then
     ! Store a temporaty copy of the full Delaunay mesh
     ! to be used later to construct the surface grid
     call copy_tri_grid()
  endif

  if (runtype /= 'MAKEGRID_PLOT') then

     if (mdomain /= 4) then

        call voronoi()
        call pcvt()

     endif

     write(io6,'(/,a)') 'gridinit after voronoi'
     write(io6,'(a,i8)')   ' nma = ',nma
     write(io6,'(a,i8)')   ' nua = ',nua
     write(io6,'(a,i8)')   ' nwa = ',nwa

     ! Allocate remaining GRID FOOTPRINT arrays for full domain

     write(io6,'(/,a)') 'gridinit calling alloc_grid1 for full domain'

     call alloc_grid1(nma, nva, nwa)

     ! Initialize long and short timesteps, and compute the timestep schedule for all operations

     write(io6,'(/,a)') 'gridinit calling modsched'

     call modsched()

     write(io6,'(/,a)') 'gridinit calling fill_jtabs'

     call fill_jtabs(nma,nva,nwa)

     ! Fill remaining GRID FOOTPRINT geometry for full domain

     write(io6,'(/,a)') 'gridinit calling grid_geometry'

     call grid_geometry_hex()

  endif

  ! Generate SFCGRID only if runtype = 'MAKEGRID' or 'MAKEGRID_PLOT'

  if (runtype == 'MAKEGRID' .or. runtype == 'MAKEGRID_PLOT') then

     if (mdomain /= 4) then
        call copyback_tri_grid()

        if (sfcgrid_res_factor > 1) then
           call expand_delaunay_mesh(sfcgrid_res_factor, .false.)
        endif
     endif

     write(io6,'(/,a)') 'gridinit calling makesfc3'

     if (mdomain <= 1) call makesfc3()

     if (runtype == 'MAKEGRID_PLOT') then
        nmsfc = nwd
        nvsfc = nud
        nwsfc = nmd
        return
     endif

  endif

  ! Allocate remaining unstructured grid geometry arrays

  write(io6,'(/,a)') 'gridinit calling alloc_grid2'

  call alloc_grid2(nma, nva, nwa)

  ! read selected SFCGRID values for MAKEREGRID run
  if (runtype == 'MAKEREGRID') call sfcgfile_read_makeregrid()

  ! Compute overlay of ATM and SFC grids
  call sfc_atm_hex_overlay()

  ! Vertical index for sea (ocean) cells

  kw_sea = 2
  do while(zm(kw_sea) < 0.1) ! .1-meter threshold
     kw_sea = kw_sea + 1
  enddo

  do iwsfc = 2, nwsfc
     if (sfcg%leaf_class(iwsfc) == 0) then

        itab_wsfc(iwsfc)%kwatm = kw_sea

     else

        kw_land = 2
        do while(zm(kw_land) < .1 + sfcg%topw(iwsfc)) ! .1-meter threshold
           kw_land = kw_land + 1
        enddo

        itab_wsfc(iwsfc)%kwatm = kw_land

     endif
  enddo

  allocate(tot_area(nwa))
  tot_area(:) = 0.0_r8

  ! Loop over all surface cells

  do iwsfc = 2,nwsfc

     ! Loop over ATM coupling areas of current surface cell

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw = itab_wsfc(iwsfc)%iwatm(j)
        tot_area(iw) = tot_area(iw) + itab_wsfc(iwsfc)%arc(j)
     enddo
  enddo

  ! Scale surface cell and ATM coupling areas slightly so that coupling areas sum
  ! almost exactly to arw0 within an ATM column

  do iwsfc = 2,nwsfc
     area_tot = 0.

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw  = itab_wsfc(iwsfc)%iwatm(j)

        itab_wsfc(iwsfc)%arc(j) = itab_wsfc(iwsfc)%arc(j) * arw0(iw) / tot_area(iw)

        area_tot = area_tot + itab_wsfc(iwsfc)%arc(j)
     enddo

     sfcg%area(iwsfc) = area_tot
  enddo

  deallocate(tot_area)

  ! Set up control volumes in atmospheric grid

  write(io6,'(/,a)') 'gridinit calling ctrlvols'

  if (mdomain <= 1) then
     call ctrlvols_hex()
  else
     call ctrlvols_hex_nosfcg()
  endif

 ! If using nudging on global grid with independent nudging grid,
 ! compute nudging indices and coefficients of W points

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

 ! Compute a mean distance between adjacent nudging points.  It is safe to assume
 ! that any ATM grid IW point will be closer than this distance to at least one
 ! nuding point

     dist_wnud = 7200.e3 / real(nudnxp)

 ! Loop over all W points

     do iw = 2,nwa

 ! Initialize minimum distance variable to large number

        dist_min = 1.e12

 ! Loop over all nudging points

        do iwnud = 2,nwnud

 ! Skip interaction if ATM grid point and nudging point are more than
 ! dist_wnud apart in any earth-frame direction

           if (abs(xew(iw) - xewnud(iwnud)) > dist_wnud) cycle
           if (abs(yew(iw) - yewnud(iwnud)) > dist_wnud) cycle
           if (abs(zew(iw) - zewnud(iwnud)) > dist_wnud) cycle

 ! Compute distance between current ATM IW grid point and nudging point,
 ! and reset minimum distance and nudging point index if new minimum is found

           dist = sqrt((xew(iw) - xewnud(iwnud))**2 &
                     + (yew(iw) - yewnud(iwnud))**2 &
                     + (zew(iw) - zewnud(iwnud))**2)

           if (dist_min > dist) then
              dist_min = dist
              itab_w(iw)%iwnud(1) = iwnud
           endif

        enddo

 ! Now that primary nudging point has been found for current IW point,
 ! compute its x,y components on a polar stereographic plane tangent
 ! at W point (W point is at 0,0)

        iwnud = itab_w(iw)%iwnud(1)

        call e_ps(xewnud(iwnud),yewnud(iwnud),zewnud(iwnud), &
                  glatw(iw),glonw(iw),xi,yi)

 ! Initialize vecprodz_minpos and vecprodz_maxneg

        vecprodz_maxneg = -1.e15
        vecprodz_minpos =  1.e15

 ! Loop through nearest polygon neighbors (j, iwnudn) of nudging point IWNUD

        npoly = itab_wnud(iwnud)%npoly

        do j = 1,npoly

 ! Get nudging point index (iwnudn) for current polygon neighbor of iwnud.

           iwnudn = itab_wnud(iwnud)%iwnud(j)

 ! Compute x,y components of iwnudn polygon center on a polar stereographic
 ! plane tangent at IW point

           call e_ps(xewnud(iwnudn),yewnud(iwnudn),zewnud(iwnudn), &
                     glatw(iw),glonw(iw),xin(j),yin(j))

 ! Compute z component (in polar stereographic space) of vector product of
 ! the vector from iwnud to iw (at 0,0) and the vector from iwnud to iwnudn.

           vecprodz = -xi * (yin(j) - yi) + yi * (xin(j) - xi)

 ! Compute scalar product of the vector from iwnud to iw (at 0,0) and the
 ! vector from iwnud to iwnudn in polar stereographic space.

           scalprod = -xi * (xin(j) - xi) - yi * (yin(j) - yi)

 ! Check whether scalar product is positive for current J point.  If so,
 ! J point is a candidate for the nudging triad for IW.

           if (scalprod > 0.) then

 ! Identify maximum negative vecprodz among all iwnudn polygon neighbors of
 ! iwnud.  This iwnudn will be second point of nudging triad for IW

              if (vecprodz < 0. .and. vecprodz > vecprodz_maxneg) then
                 vecprodz_maxneg = vecprodz
                 jmaxneg = j
                 itab_w(iw)%iwnud(2) = iwnudn
              endif

 ! Identify minimum positive vecprodz among all iwnudn polygon neighbors of
 ! iwnud.  This iwnudn will be third point of nudging triad for IW

              if (vecprodz >= 0. .and. vecprodz < vecprodz_minpos) then
                 vecprodz_minpos = vecprodz
                 jminpos = j
                itab_w(iw)%iwnud(3) = iwnudn
               endif
           endif

        enddo

 ! Lastly, fill 3 nudging weight coefficients for this W point.
 ! Weights are computed in 2_d polar stereographic space.

 ! Invert matrix of coordinates

        call matinv3x3(1.,xi,yi, &
                       1.,xin(jmaxneg),yin(jmaxneg), &
                       1.,xin(jminpos),yin(jminpos), &
                       b11,b21,b31,b12,b22,b32,b13,b23,b33)

 ! Assign coefficients

        itab_w(iw)%fnud(1) = b11
        itab_w(iw)%fnud(2) = b21
        itab_w(iw)%fnud(3) = b31

     enddo

  endif

  ! Write GRIDFILE and SFCGRILE

  write(io6,'(/,a)') 'gridinit calling gridfile_write'
  call gridfile_write()

  if (mdomain <= 1 .and. runtype /= 'MAKEREGRID') then
     write(io6,'(/,a)') 'calling sfcgfile_write'
     call sfcgfile_write()
  endif

  ! Do the following in case plotting will be done

  mma = nma
  mva = nva
  mwa = nwa

  mmsfc = nmsfc
  mvsfc = nvsfc
  mwsfc = nwsfc

  write(io6,'(/,a)') 'gridinit completed'

end subroutine gridinit

!===============================================================================

subroutine gridset2()

  use mem_grid,    only: nza, mza, &
                         zm, zt, dzm, dzt, dzim, dzit, &
                         zfacm, zfact, zfacim, zfacit, zfacm2, zfacim2, &
                         alloc_gridz
  use misc_coms,   only: io6, nzp, ndz, hdz, dz, mdomain
  use consts_coms, only: erad
  use oname_coms,  only: nl

  implicit none

  integer :: idz, kvec, nseries, iseries, k, kd

  real :: zend, dzend, dzbeg, ztarg, dzr, rdzr, zshift

  real, allocatable :: zmvec(:),ztvec(:)

  ! Allocate zmvec and ztvec arrays to large size.

  allocate (zmvec(0:500),ztvec(0:500))

  ! If NDZ = 1, then use ZZ values from namelist

  if (ndz == 1) then

     nza = nzp
     mza = nzp

     zmvec(0) = 2. * nl%zz(1) - nl%zz(2)
     zmvec(1:nza) = nl%zz(1:nza)
     zmvec(nza+1) = 2. * zmvec(nza) - zmvec(nza-1)

     ! compute zt values by geometric interpolation.

     ztvec(1) = .5 * (zmvec(0) + zmvec(1))
     do k = 2,nza
        dzr = sqrt(sqrt((zmvec(k+1) - zmvec(k)) / (zmvec(k-1) - zmvec(k-2))))
        ztvec(k) = zmvec(k-1) + (zmvec(k) - zmvec(k-1)) / (1. + dzr)
     enddo
     ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

  else

  ! IF NDZ > 1, this subroutine will compute number of vertical levels.
  ! Fill starting values of zmvec

     zmvec(0) = hdz(1) - dz(1)
     zmvec(1) = hdz(1)
     zmvec(2) = hdz(1) + dz(1)

     ztvec(1) = hdz(1) - dz(1) * .5
     ztvec(2) = hdz(1) + dz(1) * .5

     kvec = 2

     ! Loop through (hdz,dz) values

     do idz = 2,ndz

        zend = hdz(idz)
        dzend = dz(idz)

        dzbeg = zmvec(kvec) - zmvec(kvec-1)

        ! Target height for next series

        ztarg = zend + .4 * dzend

        ! Stretch ratio based on geometric series sum

        dzr = (ztarg - zmvec(kvec)) / (ztarg - zmvec(kvec-1) - dzend)

        ! NSERIES from stretch ratio

        if (abs(dzr - 1.) > 1.e-3) then
           nseries = max(1,nint(log(dzend/dzbeg) / log(dzr)))
        else
           nseries = max(1,nint((ztarg - zmvec(kvec)) / dzbeg))
        endif

        ! Stretch ratio based on N

        dzr = (dzend / dzbeg) ** (1./real(nseries))
        rdzr = sqrt(dzr)

        ! Loop over series and compute new levels

        do iseries = 1,nseries

           if (kvec >= nzp - 1) GO TO 5

           kvec = kvec + 1
           dzbeg = dzbeg * dzr

           zmvec(kvec) = zmvec(kvec-1) + dzbeg
           ztvec(kvec) = zmvec(kvec-1) + (zmvec(kvec) - zmvec(kvec-1)) / (1. + rdzr)

        enddo

     enddo

     5 continue

     ! If hdz(1) < 0 in order to accomodate geographic areas that are below
     ! sea level, determine the height of the lowest zmvec level that is above
     ! sea level.  Reduce all zmvec and ztvec values by the amount of that height
     ! so that one zmvec value exactly coincides with sea level.  Then, determine
     ! how many zmvec values are below the value of hdz(1).  Discard all such
     ! levels and copy the remainder to arrays zm and zt.

     kd = 0

     if (hdz(1) < 0) then
        k = 2
        do while (zmvec(k) < 0.)
           k = k + 1
        enddo

        zshift = zmvec(k)

        zmvec(:) = zmvec(:) - zshift
        ztvec(:) = ztvec(:) - zshift

        kd = 0
        do while (zmvec(kd+2) < hdz(1))
           kd = kd + 1
        enddo
     endif

     ! Fill NZA and MZA values

     nza = kvec - kd
     mza = nza

     if (nza < 2) then
        write(io6,'(a)')      'OLAM requires that NZA >= 2.'
        write(io6,'(a,i5,a)') 'NZA = ',nza,' in subroutine gridset2.  Stopping model.'
        stop 'stop nza in gridset2'
     endif

     ! Fill top 2 ZMVEC and ZTVEC values

     zmvec(kvec)   = zmvec(kvec-1) * 2. - zmvec(kvec-2)
     zmvec(kvec+1) = zmvec(kvec)   * 2. - zmvec(kvec-1)

     ztvec(kvec)   = .5 * (zmvec(kvec-1) + zmvec(kvec))
     ztvec(kvec+1) = .5 * (zmvec(kvec) + zmvec(kvec+1))

  endif

  ! Allocate main grid arrays

  call alloc_gridz()

  ! Other vertical coordinate values

  do k = 1,nza
     zm(k) = zmvec(k+kd)
     zt(k) = ztvec(k+kd)

     dzm(k) = ztvec(k+kd+1) - ztvec(k+kd)
     dzt(k) = zmvec(k+kd) - zmvec(k+kd-1)

     dzim(k) = 1. / dzm(k)
     dzit(k) = 1. / dzt(k)

     if (mdomain < 2) then
        zfacm(k) = (erad + zm(k)) / erad
        zfact(k) = (erad + zt(k)) / erad

        !-------------------------------------------
        ! DCMIP modifications [may not be required in some of these cases]
        if (nl%test_case ==  11 .or. &
            nl%test_case ==  12 .or. &
            nl%test_case ==  13 .or. &
            nl%test_case ==  20 .or. &
            nl%test_case ==  31 .or. &
            nl%test_case ==  42 .or. &
            nl%test_case ==  43 .or. &
            nl%test_case ==  51 .or. &
          ! nl%test_case == 131 .or. & ! First using deep atmos for this case; try shallow later
            nl%test_case ==  52) then

           zfacm(k) = 1.
           zfact(k) = 1.
       endif
       !-------------------------------------------

     else
        zfacm(k) = 1.
        zfact(k) = 1.
     endif

     zfacm2(k) = zfacm(k)**2

     zfacim (k) = 1. / zfacm (k)
     zfacit (k) = 1. / zfact (k)
     zfacim2(k) = 1. / zfacm2(k)
  enddo

  deallocate (zmvec,ztvec)

end subroutine gridset2

!===============================================================================

subroutine grav_init()

  use mem_grid,    only: zfacim, zfacit, gravm, gravt
  use consts_coms, only: grav
  use oname_coms,  only: nl

  implicit none

  ! Currently, variable gravity only for DCMIP 2016 baroclinic wave test case;
  ! later adopt for general case

  if (nl%test_case == 110 .or. &
      nl%test_case == 111 .or. &
      nl%test_case == 112 .or. &
      nl%test_case == 113 .or. &
      nl%test_case == 114) then

     gravm(:) = grav * zfacim(:)**2
     gravt(:) = grav * zfacit(:)**2
  else
     gravm(:) = grav
     gravt(:) = grav
  endif

end subroutine grav_init

!===============================================================================

subroutine gridset_print()

use misc_coms, only: io6, runtype, mdomain
use mem_grid,  only: nza, zm, zt, dzt, nma, nva, nwa
use mem_sfcg,  only: nwsfc
use mem_land,  only: nland
use mem_lake,  only: nlake
use mem_sea,   only: nsea

implicit none

integer :: k

! Print grid structure

write(io6,'(/,a)') '================================================='
write(io6,'(a)'  ) '         OLAM VERTICAL GRID STRUCTURE'
write(io6,'(a)'  ) '================================================='
write(io6,'(a)'  ) '   k    zm(m)     k     zt(m)    dzm(m)    dzrat'
write(io6,'(a,/)') '================================================='

do k = nza,1,-1
   if (k == nza) then
      write(io6,11) k, zm(k), dzt(k)/dzt(k)
   elseif (k == 1) then
      write(io6,11) k, zm(k), dzt(k+1)/dzt(k)
   else
      write(io6,12) k, zm(k), dzt(k+1)/dzt(k)
   endif
   if (k > 1) write(io6,13) k,zt(k),dzt(k)
enddo

11 format (i4,f10.2,1x,3('========='),f6.3)
12 format (i4,f10.2,1x,3('---------'),f6.3)
13 format (15x,i4,2f10.2)

write(io6,'(/,a)' ) 'Global grid indices:'
write(io6,'(a,i9)') '  nwa = ',nwa
write(io6,'(a,i9)') '  nva = ',nva
write(io6,'(a,i9)') '  nma = ',nma
write(io6,'(a,i9)') '  nza = ',nza

if (mdomain <= 1) then

   if (runtype == 'MAKEGRID_PLOT') then

      write(io6,'(a,i9)') 'nwsfc = ',nwsfc

   else

      write(io6,'(/,a)' ) 'Global sfcg/land/lake/sea indices:'
      write(io6,'(a,i9)') '  nwsfc = ',nwsfc
      write(io6,'(a,i9)') '  nland = ',nland
      write(io6,'(a,i9)') '  nlake = ',nlake
      write(io6,'(a,i9)') '  nsea  = ',nsea
   endif

endif

end subroutine gridset_print

!===============================================================================

subroutine gridfile_write()

  use max_dims,   only: maxngrdll, pathlen
  use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                        iclobber, itopoflg, deltax, ndz, hdz, dz, &
                        ngrdll, grdrad, grdlat, grdlon
  use mem_ijtabs, only: mloops, mrls, itab_m, itab_v, itab_w
  use mem_grid,   only: nza, nma, nua, nva, nwa, nsw_max, &
                        zm, zt, dzm, dzt, dzim, dzit, &
                        zfacm, zfact, zfacim, zfacit, zfacm2, zfacim2, &
                        lpm, lpv, lpw, lsw, lve2, nve2_max, &
                        topm, topw, xem, yem, zem, &
                        xev, yev, zev, xew, yew, zew, &
                        unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                        dnu, dniu, dnv, dniv, arw0, arm0, &
                        glatw, glonw, glatm, glonm, glatv, glonv, &
                        arv, arw, volt
  use mem_sfcg,   only: sfcg, nwsfc, itab_wsfc
  use leaf_coms,  only: isfcl

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

  use mem_nudge,  only: nudflag, nudnxp, nwnud, itab_wnud, &
                        xewnud, yewnud, zewnud

  implicit none

  ! This routine writes the grid variables to the gridfile.

  integer :: im, iv, iw, iwnud

  integer :: ndims, idims(2)

  ! Scratch arrays for copying output

  integer, allocatable :: iscr1(:)
  real,    allocatable :: rscr1(:)
  integer, allocatable :: iscr2(:,:)
  real,    allocatable :: rscr2(:,:)
  logical, allocatable :: lscr2(:,:)

  character(pathlen) :: flnm

  ! Open gridfile

  flnm = trim(gridfile)//'.h5'

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'grid_write: opening file:', flnm
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'W',iclobber)

  ! Write the gridfile information that exists in namelist

  ndims    = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_orec(ndims, idims, 'NZP'     , ivars=nzp)
  call shdf5_orec(ndims, idims, 'NXP'     , ivars=nxp)
  call shdf5_orec(ndims, idims, 'MDOMAIN' , ivars=mdomain)
  call shdf5_orec(ndims, idims, 'NGRIDS'  , ivars=ngrids)
  call shdf5_orec(ndims, idims, 'ISFCL'   , ivars=isfcl)
  call shdf5_orec(ndims, idims, 'ITOPOFLG', ivars=itopoflg)
  call shdf5_orec(ndims, idims, 'DELTAX'  , rvars=deltax)
  call shdf5_orec(ndims, idims, 'NDZ'     , ivars=ndz)
  call shdf5_orec(ndims, idims, 'NUDFLAG' , ivars=nudflag)
  call shdf5_orec(ndims, idims, 'NUDNXP'  , ivars=nudnxp)

  if (ndz > 1) then
     ndims    = 1
     idims(1) = ndz

     call shdf5_orec(ndims, idims, 'HDZ'  , rvar1=hdz)
     call shdf5_orec(ndims, idims, 'DZ'   , rvar1=dz)
  endif

  ndims    = 1
  idims(1) = ngrids

  call shdf5_orec(ndims, idims, 'NGRDLL' , ivar1=ngrdll)

  ndims    = 2
  idims(1) = ngrids
  idims(2) = maxngrdll

  call shdf5_orec(ndims, idims, 'GRDRAD' ,rvar2=grdrad(1:ngrids,:))
  call shdf5_orec(ndims, idims, 'GRDLAT', rvar2=grdlat(1:ngrids,:))
  call shdf5_orec(ndims, idims, 'GRDLON', rvar2=grdlon(1:ngrids,:))

  ! Write the grid dimensions

  ndims    = 1
  idims(1) = 1

  call shdf5_orec(ndims, idims, 'NZA'     , ivars=nza)
  call shdf5_orec(ndims, idims, 'NMA'     , ivars=nma)
  call shdf5_orec(ndims, idims, 'NUA'     , ivars=nua)
  call shdf5_orec(ndims, idims, 'NVA'     , ivars=nva)
  call shdf5_orec(ndims, idims, 'NWA'     , ivars=nwa)
  call shdf5_orec(ndims, idims, 'NSW_MAX' , ivars=nsw_max)
  call shdf5_orec(ndims, idims, 'NVE2_MAX', ivars=nve2_max)
  call shdf5_orec(ndims, idims, 'MRLS'    , ivars=mrls)
  call shdf5_orec(ndims, idims, 'NWNUD'   , ivars=nwnud)

  ! Write grid structure variables

  ndims    = 1
  idims(1) = nza

  call shdf5_orec(ndims, idims, 'ZM'    , rvar1=zm)
  call shdf5_orec(ndims, idims, 'ZT'    , rvar1=zt)
  call shdf5_orec(ndims, idims, 'DZM'   , rvar1=dzm)
  call shdf5_orec(ndims, idims, 'DZT'   , rvar1=dzt)
  call shdf5_orec(ndims, idims, 'DZIM'  , rvar1=dzim)
  call shdf5_orec(ndims, idims, 'DZIT'  , rvar1=dzit)
  call shdf5_orec(ndims, idims, 'ZFACM' , rvar1=zfacm)
  call shdf5_orec(ndims, idims, 'ZFACT' , rvar1=zfact)
  call shdf5_orec(ndims, idims, 'ZFACIM', rvar1=zfacim)
  call shdf5_orec(ndims, idims, 'ZFACIT', rvar1=zfacit)
  call shdf5_orec(ndims, idims, 'ZFACM2' , rvar1=zfacm2)
  call shdf5_orec(ndims, idims, 'ZFACIM2', rvar1=zfacim2)

  ndims    = 1
  idims(1) = nma

  call shdf5_orec(ndims, idims, 'LPM'  , ivar1=lpm)
  call shdf5_orec(ndims, idims, 'TOPM' , rvar1=topm)
  call shdf5_orec(ndims, idims, 'XEM'  , rvar1=xem)
  call shdf5_orec(ndims, idims, 'YEM'  , rvar1=yem)
  call shdf5_orec(ndims, idims, 'ZEM'  , rvar1=zem)
  call shdf5_orec(ndims, idims, 'ARM0' , rvar1=arm0)
  call shdf5_orec(ndims, idims, 'GLATM', rvar1=glatm)
  call shdf5_orec(ndims, idims, 'GLONM', rvar1=glonm)

  ndims    = 1
  idims(1) = nva

  call shdf5_orec(ndims, idims, 'LPV'  , ivar1=lpv)
  call shdf5_orec(ndims, idims, 'XEV'  , rvar1=xev)
  call shdf5_orec(ndims, idims, 'YEV'  , rvar1=yev)
  call shdf5_orec(ndims, idims, 'ZEV'  , rvar1=zev)
  call shdf5_orec(ndims, idims, 'GLATV', rvar1=glatv)
  call shdf5_orec(ndims, idims, 'GLONV', rvar1=glonv)

  call shdf5_orec(ndims, idims, 'DNU' , rvar1=dnu)
  call shdf5_orec(ndims, idims, 'DNV' , rvar1=dnv)
  call shdf5_orec(ndims, idims, 'DNIU', rvar1=dniu)
  call shdf5_orec(ndims, idims, 'DNIV', rvar1=dniv)
  call shdf5_orec(ndims, idims, 'UNX' , rvar1=unx)
  call shdf5_orec(ndims, idims, 'UNY' , rvar1=uny)
  call shdf5_orec(ndims, idims, 'UNZ' , rvar1=unz)
  call shdf5_orec(ndims, idims, 'VNX' , rvar1=vnx)
  call shdf5_orec(ndims, idims, 'VNY' , rvar1=vny)
  call shdf5_orec(ndims, idims, 'VNZ' , rvar1=vnz)

  ndims    = 1
  idims(1) = nwa

  call shdf5_orec(ndims, idims, 'LPW'  , ivar1=lpw)
  call shdf5_orec(ndims, idims, 'LSW'  , ivar1=lsw)
  call shdf5_orec(ndims, idims, 'LVE2' , ivar1=lve2)
  call shdf5_orec(ndims, idims, 'XEW'  , rvar1=xew)
  call shdf5_orec(ndims, idims, 'YEW'  , rvar1=yew)
  call shdf5_orec(ndims, idims, 'ZEW'  , rvar1=zew)
  call shdf5_orec(ndims, idims, 'TOPW' , rvar1=topw)
  call shdf5_orec(ndims, idims, 'ARW0' , rvar1=arw0)
  call shdf5_orec(ndims, idims, 'GLATW', rvar1=glatw)
  call shdf5_orec(ndims, idims, 'GLONW', rvar1=glonw)
  call shdf5_orec(ndims, idims, 'WNX'  , rvar1=wnx)
  call shdf5_orec(ndims, idims, 'WNY'  , rvar1=wny)
  call shdf5_orec(ndims, idims, 'WNZ'  , rvar1=wnz)

  ndims    = 2
  idims(1) = nza
  idims(2) = nva

  call shdf5_orec(ndims, idims, 'ARV'  , rvar2=arv)

  ndims    = 2
  idims(1) = nza
  idims(2) = nwa

  call shdf5_orec(ndims, idims, 'ARW'  , rvar2=arw)
  call shdf5_orec(ndims, idims, 'VOLT' , dvar2=volt)

  ! Write ITAB_M SCALARS

  ndims    = 1
  idims(1) = nma

  allocate (iscr1(nma)) ; iscr1 = 0
  do im = 1,nma
     iscr1(im) = itab_m(im)%npoly
  enddo
  call shdf5_orec(ndims,idims,'itab_m%npoly'    ,ivar1=iscr1)

  do im = 1,nma
     iscr1(im) = itab_m(im)%imp
  enddo
  call shdf5_orec(ndims,idims,'itab_m%imp'      ,ivar1=iscr1)

  do im = 1,nma
     iscr1(im) = itab_m(im)%mrlm
  enddo
  call shdf5_orec(ndims,idims,'itab_m%mrlm'     ,ivar1=iscr1)

  do im = 1,nma
     iscr1(im) = itab_m(im)%mrlm_orig
  enddo
  call shdf5_orec(ndims,idims,'itab_m%mrlm_orig',ivar1=iscr1)

  do im = 1,nma
     iscr1(im) = itab_m(im)%mrow
  enddo
  call shdf5_orec(ndims,idims,'itab_m%mrow'     ,ivar1=iscr1)

  do im = 1,nma
     iscr1(im) = itab_m(im)%ngr
  enddo
  call shdf5_orec(ndims,idims,'itab_m%ngr'      ,ivar1=iscr1)
  deallocate (iscr1)

  ! Write ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = nma

  allocate (lscr2(mloops,nma)) ; lscr2 = .false.
  do im = 1,nma
     lscr2(1:mloops,im) = itab_m(im)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%loop',lvar2=lscr2)
  deallocate (lscr2)

  ndims    = 2
  idims(1) = 3
  idims(2) = nma

  allocate (iscr2(3,nma)) ; iscr2 = 0
  do im = 1,nma
     iscr2(1:3,im) = itab_m(im)%iv(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%iv',ivar2=iscr2)

  do im = 1,nma
     iscr2(1:3,im) = itab_m(im)%iw(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%iw',ivar2=iscr2)

  do im = 1,nma
     iscr2(1:3,im) = itab_m(im)%im(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%im',ivar2=iscr2)
  deallocate (iscr2)

  ! Write ITAB_V SCALARS

  ndims    = 1
  idims(1) = nva
  idims(2) = 1

  allocate (iscr1(nva)) ; iscr1 = 0
  do iv = 1,nva
     iscr1(iv) = itab_v(iv)%ivp
  enddo
  call shdf5_orec(ndims,idims,'itab_v%ivp'    ,ivar1=iscr1)

  do iv = 1,nva
     iscr1(iv) = itab_v(iv)%mrlv
  enddo
  call shdf5_orec(ndims,idims,'itab_v%mrlv'   ,ivar1=iscr1)
  deallocate (iscr1)

  ! Write ITAB_V ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = nva

  allocate (lscr2(mloops,nva))
  do iv = 1,nva
     lscr2(1:mloops,iv) = itab_v(iv)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%loop',lvar2=lscr2)
  deallocate(lscr2)

  idims(1) = 2
  allocate (rscr2(2,nva))

  do iv = 1,nva
     rscr2(1:2,iv) = itab_v(iv)%farw(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%farw',rvar2=rscr2)

  do iv = 1,nva
     rscr2(1:2,iv) = itab_v(iv)%cosv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%cosv',rvar2=rscr2)

  do iv = 1,nva
     rscr2(1:2,iv) = itab_v(iv)%sinv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%sinv',rvar2=rscr2)

  do iv = 1,nva
     rscr2(1:2,iv) = itab_v(iv)%dxps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%dxps',rvar2=rscr2)

  do iv = 1,nva
     rscr2(1:2,iv) = itab_v(iv)%dyps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%dyps',rvar2=rscr2)

  deallocate(rscr2)

! idims(1) = 4
! allocate (iscr2(4,nva))
!
! do iv = 1,nva
!    iscr2(1:4,iv) = itab_v(iv)%iw(1:4)
! enddo
! call shdf5_orec(ndims,idims,'itab_v%iw',ivar2=iscr2)
!
! do iv = 1,nva
!    iscr2(1:4,iv) = itab_v(iv)%iv(1:4)
! enddo
! call shdf5_orec(ndims,idims,'itab_v%iv',ivar2=iscr2)
!
! deallocate(iscr2)
!
! idims(1) = 6
! allocate (iscr2(6,nva))
! do iv = 1,nva
!    iscr2(1:6,iv) = itab_v(iv)%im(1:6)
! enddo
! call shdf5_orec(ndims,idims,'itab_v%im',ivar2=iscr2)
! deallocate(iscr2)

  idims(1) = 2
  allocate (iscr2(2,nva))

  do iv = 1,nva
     iscr2(1:2,iv) = itab_v(iv)%iw(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%iw',ivar2=iscr2)

  do iv = 1,nva
     iscr2(1:2,iv) = itab_v(iv)%im(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%im',ivar2=iscr2)

  deallocate(iscr2)

  ! Write ITAB_W SCALARS

  ndims    = 1
  idims(1) = nwa

  allocate (iscr1(nwa)) ; iscr1 = 0
  do iw = 1,nwa
     iscr1(iw) = itab_w(iw)%npoly
  enddo
  call shdf5_orec(ndims,idims,'itab_w%npoly'    ,ivar1=iscr1)

  do iw = 1,nwa
     iscr1(iw) = itab_w(iw)%iwp
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iwp'      ,ivar1=iscr1)

  do iw = 1,nwa
     iscr1(iw) = itab_w(iw)%mrlw
  enddo
  call shdf5_orec(ndims,idims,'itab_w%mrlw'     ,ivar1=iscr1)

  do iw = 1,nwa
     iscr1(iw) = itab_w(iw)%mrlw_orig
  enddo
  call shdf5_orec(ndims,idims,'itab_w%mrlw_orig',ivar1=iscr1)

  do iw = 1,nwa
     iscr1(iw) = itab_w(iw)%ngr
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ngr'      ,ivar1=iscr1)
  deallocate(iscr1)

  allocate (rscr1(nwa)) ; rscr1 = 0
  do iw = 1,nwa
     rscr1(iw) = itab_w(iw)%unx_w
  enddo
  call shdf5_orec(ndims,idims,'itab_w%unx_w' ,rvar1=rscr1)

  do iw = 1,nwa
     rscr1(iw) = itab_w(iw)%uny_w
  enddo
  call shdf5_orec(ndims,idims,'itab_w%uny_w' ,rvar1=rscr1)

  do iw = 1,nwa
     rscr1(iw) = itab_w(iw)%vnx_w
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vnx_w' ,rvar1=rscr1)

  do iw = 1,nwa
     rscr1(iw) = itab_w(iw)%vny_w
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vny_w' ,rvar1=rscr1)

  do iw = 1,nwa
     rscr1(iw) = itab_w(iw)%vnz_w
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vnz_w' ,rvar1=rscr1)
  deallocate(rscr1)

  ! Write ITAB_W ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = nwa

  allocate (lscr2(mloops,nwa))
  do iw = 1,nwa
     lscr2(1:mloops,iw) = itab_w(iw)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%loop',lvar2=lscr2)
  deallocate(lscr2)

  idims(1) = 3

  allocate (iscr2(3,nwa))
  do iw = 1,nwa
     iscr2(1:3,iw) = itab_w(iw)%iwnud(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iwnud',ivar2=iscr2)
  deallocate(iscr2)

  allocate (rscr2(3,nwa))
  do iw = 1,nwa
     rscr2(1:3,iw) = itab_w(iw)%fnud(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%fnud',rvar2=rscr2)
  deallocate(rscr2)

  idims(1) = 7
  allocate (iscr2(7,nwa))

  do iw = 1,nwa
     iscr2(1:7,iw) = itab_w(iw)%im(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%im',ivar2=iscr2)

  do iw = 1,nwa
     iscr2(1:7,iw) = itab_w(iw)%iv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iv',ivar2=iscr2)

  do iw = 1,nwa
     iscr2(1:7,iw) = itab_w(iw)%iw(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iw',ivar2=iscr2)

  deallocate(iscr2)
  allocate (rscr2(7,nwa))

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%dirv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%dirv',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%farm(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%farm',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%farv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%farv',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%gxps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gxps1',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%gyps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gyps1',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%gxps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gxps2',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%gyps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gyps2',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%ecvec_vx(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vx',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%ecvec_vy(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vy',rvar2=rscr2)

  do iw = 1,nwa
     rscr2(1:7,iw) = itab_w(iw)%ecvec_vz(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vz',rvar2=rscr2)

  deallocate(rscr2)

  ! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     ndims    = 1
     idims(1) = nwnud

     call shdf5_orec(ndims, idims, 'XEWNUD'  , rvar1=xewnud)
     call shdf5_orec(ndims, idims, 'YEWNUD'  , rvar1=yewnud)
     call shdf5_orec(ndims, idims, 'ZEWNUD'  , rvar1=zewnud)

     allocate (iscr1(nwnud)) ; iscr1 = 0
     do iwnud = 1,nwnud
        iscr1(iwnud) = itab_w(iwnud)%npoly
     enddo
     call shdf5_orec(ndims,idims,'itab_wnud%npoly',ivar1=iscr1)
     deallocate(iscr1)

     allocate (iscr2(6,nwnud))

     do iwnud = 1,nwnud
        iscr2(1:6,iwnud) = itab_wnud(iwnud)%iwnud(1:6)
     enddo

     ndims    = 2
     idims(1) = 6
     idims(2) = nwnud

     call shdf5_orec(ndims,idims,'itab_wnud%iwnud' ,ivar2=iscr2)

     deallocate(iscr2)

  endif

  ! If SFCGRID is active, write ATM-SFCGRID coupling information

  if (mdomain <= 1) then

     ! Write the grid dimensions

     ndims = 1
     idims(1) = 1

     call shdf5_orec(ndims, idims, 'nwsfc'  , ivars=nwsfc)

     ndims    = 1
     idims(1) = nwsfc

     call shdf5_orec(ndims, idims, 'dzt_bot'   , rvar1=sfcg%dzt_bot)

     allocate (iscr1(nwsfc))
     do iw = 1,nwsfc
        iscr1(iw) = itab_wsfc(iw)%nwatm
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%nwatm'    ,ivar1=iscr1)
     deallocate(iscr1)

     ndims = 2
     idims(1) = 8
     idims(2) = nwsfc

     allocate (iscr2(8,nwsfc))
     do iw = 1,nwsfc
        iscr2(1:8,iw) = itab_wsfc(iw)%iwatm(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%iwatm',ivar2=iscr2)

     do iw = 1,nwsfc
        iscr2(1:8,iw) = itab_wsfc(iw)%kwatm(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%kwatm',ivar2=iscr2)
     deallocate(iscr2)

     allocate (rscr2(8,nwsfc))
     do iw = 1,nwsfc
        rscr2(1:8,iw) = itab_wsfc(iw)%arc(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%arc',rvar2=rscr2)

     do iw = 1,nwsfc
        rscr2(1:8,iw) = itab_wsfc(iw)%arcoarsfc(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%arcoarsfc',rvar2=rscr2)

     do iw = 1,nwsfc
        rscr2(1:8,iw) = itab_wsfc(iw)%arcoariw(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%arcoariw',rvar2=rscr2)

     do iw = 1,nwsfc
        rscr2(1:8,iw) = itab_wsfc(iw)%arcoarkw(1:8)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%arcoarkw',rvar2=rscr2)
     deallocate(rscr2)

  endif

  ! Close GRIDFILE

  call shdf5_close()

end subroutine gridfile_write

!===============================================================================

! This subroutine reads all the scalars values, and only a few arrays needed
! by para_decomp (arrays ended with _pd)

subroutine gridfile_read_pd()

  use max_dims,   only: maxngrdll, pathlen
  use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                        itopoflg, deltax, ndz, hdz, dz, &
                        ngrdll, grdrad, grdlat, grdlon
  use mem_ijtabs, only: mloops, mrls,  &
                        itab_m_pd, itab_v_pd, itab_w_pd, alloc_itabs_pd
  use mem_grid,   only: nza, nma, nua, nva, nwa, &
                        mza, mma, mua, mva, mwa, nsw_max, nve2_max, &
                        xew, yew, zew, arw0, alloc_xyzem, alloc_xyzew
  use mem_sfcg,   only: nwsfc, itab_wsfc_pd
  use leaf_coms,  only: isfcl

  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close

  use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, xewnud, yewnud, zewnud, &
                        alloc_nudge1

  ! This subroutine checks for the existence of a gridfile, and if it exists,
  ! also checks for agreement of grid configuration between the file and the
  ! current model run.  If the file does not exist or does not match grid
  ! configuration, the run is stopped.

  implicit none

  integer :: im, iv, iw, idz, iwsfc

  integer :: ierr

  integer :: ngr, i
  integer :: ndims, idims(2)

  integer :: ngrids0, mdomain0, nxp0, nzp0, itopoflg0, isfcl0, ndz0
  integer :: nudflag0, nudnxp0

  real    :: deltax0

  logical :: exans

  integer, allocatable :: ngrdll0(:)
  real,    allocatable :: grdrad0(:,:)
  real,    allocatable :: grdlat0(:,:)
  real,    allocatable :: grdlon0(:,:)
  real,    allocatable :: hdz0(:)
  real,    allocatable :: dz0(:)

  ! Scratch arrays for copying input

  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  character(pathlen) :: flnm

  ! Check if gridfile exists

  flnm = trim(gridfile)//'.h5'

  inquire(file=flnm, exist=exans)

  if (.not. exans) then

  ! Gridfile does not exist.

     write(io6,*)
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!  Gridfile does not exist:'
     write(io6,*) '!!!  '//gridfile
     write(io6,*) '!!!  Stopping run'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

     stop 'stop - no gridfile'

  endif

  ! Gridfile exists; open it

  write(io6,*)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'Opening grid file ', trim(flnm)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R',trypario=.true.)

  ! Read the grid information that exists in namelist

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp0)
  call shdf5_irec(ndims, idims, 'NXP'     , ivars=nxp0)
  call shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain0)
  call shdf5_irec(ndims, idims, 'NGRIDS'  , ivars=ngrids0)
  call shdf5_irec(ndims, idims, 'ISFCL'   , ivars=isfcl0)
  call shdf5_irec(ndims, idims, 'ITOPOFLG', ivars=itopoflg0)
  call shdf5_irec(ndims, idims, 'DELTAX'  , rvars=deltax0)

  call shdf5_irec(ndims, idims, 'NDZ'     , ivars=ndz0)

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
     call shdf5_irec(ndims, idims, 'NUDFLAG' , ivars=nudflag0)
     call shdf5_irec(ndims, idims, 'NUDNXP'  , ivars=nudnxp0)
  endif

  allocate( ngrdll0 (ngrids0) )
  allocate( grdrad0 (ngrids0, maxngrdll) )
  allocate( grdlat0 (ngrids0, maxngrdll) )
  allocate( grdlon0 (ngrids0, maxngrdll) )

  if (ndz0 > 1) then
     allocate( hdz0(ndz0))
     allocate( dz0 (ndz0))

     idims(1) = ndz0
     call shdf5_irec(ndims, idims, 'HDZ' , rvar1=hdz0)
     call shdf5_irec(ndims, idims, 'DZ'  , rvar1=dz0)
  endif

  ndims    = 1
  idims(1) = ngrids0

  call shdf5_irec(ndims, idims, 'NGRDLL' , ivar1=ngrdll0)

  ndims    = 2
  idims(1) = ngrids0
  idims(2) = maxngrdll

  call shdf5_irec(ndims, idims, 'GRDRAD', rvar2=grdrad0)
  call shdf5_irec(ndims, idims, 'GRDLAT', rvar2=grdlat0)
  call shdf5_irec(ndims, idims, 'GRDLON', rvar2=grdlon0)

  ! Check equality between grid file information and namelist variables

  ierr = 0

  if (nzp0      /= nzp     ) ierr = 1
  if (nxp0      /= nxp     ) ierr = 1
  if (mdomain0  /= mdomain ) ierr = 1
  if (ngrids0   /= ngrids  ) ierr = 1
  if (isfcl0    /= isfcl   ) ierr = 1
  if (itopoflg0 /= itopoflg) ierr = 1
  if (ndz0      /= ndz     ) ierr = 1

  if (abs(deltax0 - deltax) > 1.e-3) ierr = 1

  do ngr = 2, min(ngrids0,ngrids)
     if (abs(ngrdll0 (ngr) - ngrdll (ngr)) > 1.e1 ) ierr = 1

     do i = 1,ngrdll0(ngr)
        if (abs(grdrad0(ngr,i) - grdrad(ngr,i)) > 1.e1 ) ierr = 1
        if (abs(grdlat0(ngr,i) - grdlat(ngr,i)) > 1.e-3) ierr = 1
        if (abs(grdlon0(ngr,i) - grdlon(ngr,i)) > 1.e-3) ierr = 1
     enddo
  enddo

  if (ndz0 > 1 .and. ndz > 1) then
     do idz = 1, min(ndz0,ndz)
        if (abs(hdz0(idz) - hdz(idz)) > 1.e1 ) ierr = 1
        if (abs(dz0 (idz) - dz (idz)) > 1.e1 ) ierr = 1
     enddo
  endif

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
     if (nudflag0 /= nudflag) ierr = 1
     if (nudnxp0  /= nudnxp ) ierr = 1
  endif

  if (ierr == 1) then

     write(io6,*) 'GRIDFILE mismatch with OLAMIN namelist: Stopping model run'
     write(io6,*) 'Values: gridfile, namelist'
     write(io6,*) '-----------------------------------------------'
     write(io6,*)              'nzp:      ',nzp0     ,nzp
     write(io6,*)              'nxp:      ',nxp0     ,nxp
     write(io6,*)              'mdomain:  ',mdomain0 ,mdomain
     write(io6,*)              'ngrids:   ',ngrids0  ,ngrids
     write(io6,*)              'isfcl:    ',isfcl0   ,isfcl
     write(io6,*)              'itopoflg: ',itopoflg0,itopoflg
     write(io6,*)              'deltax:   ',deltax0  ,deltax
     write(io6,*)              'ndz:      ',ndz0     ,ndz
     write(io6,*) ' '

     if (ndz0 > 1 .and. ndz > 1) then
        write(io6, '(a,20f10.1)') 'hdz0:     ',hdz0 (1:ndz)
        write(io6, '(a,20f10.1)') 'hdz:      ',hdz  (1:ndz)
        write(io6,*) ' '
        write(io6, '(a,20f10.1)') 'dz0:      ',dz0 (1:ndz)
        write(io6, '(a,20f10.1)') 'dz:       ',dz  (1:ndz)
     endif

     write(io6,*) ' '
     write(io6, '(a,20i12)')   'ngrdll0:  ',ngrdll0 (1:ngrids)
     write(io6, '(a,20i12)')   'ngrdll:   ',ngrdll  (1:ngrids)
     write(io6,*) ' '

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        write(io6,*)              'nudflag:  ',nudflag0 ,nudflag
        write(io6,*)              'nudnxp:   ',nudnxp0  ,nudnxp
     endif

     do ngr = 2, min(ngrids0,ngrids)
        write(io6, '(a,i5)') 'ngr: ',ngr
        write(io6,*) ' '
        write(io6, '(a,20f12.1)') 'grdrad0:  ',grdrad0(ngr,1:ngrdll(ngr))
        write(io6, '(a,20f12.1)') 'grdrad:   ',grdrad (ngr,1:ngrdll(ngr))
        write(io6,*) ' '
        write(io6, '(a,20f10.3)') 'grdlat0: ',grdlat0(ngr,1:ngrdll(ngr))
        write(io6, '(a,20f10.3)') 'grdlat:  ',grdlat (ngr,1:ngrdll(ngr))
        write(io6,*) ' '
        write(io6, '(a,20f10.3)') 'grdlon0: ',grdlon0(ngr,1:ngrdll(ngr))
        write(io6, '(a,20f10.3)') 'grdlon:  ',grdlon (ngr,1:ngrdll(ngr))
        write(io6,*) ' '
     enddo

     write(io6,*) '-----------------------------------------------'

     stop 'stop - gridfile mismatch'

  endif

  deallocate (ngrdll0, grdrad0, grdlat0, grdlon0)

  ! Read the grid dimensions

  ndims    = 1
  idims(1) = 1

  call shdf5_irec(ndims, idims, 'NZA'     , ivars=nza)
  call shdf5_irec(ndims, idims, 'NMA'     , ivars=nma)
  call shdf5_irec(ndims, idims, 'NUA'     , ivars=nua)
  call shdf5_irec(ndims, idims, 'NVA'     , ivars=nva)
  call shdf5_irec(ndims, idims, 'NWA'     , ivars=nwa)
  call shdf5_irec(ndims, idims, 'NSW_MAX' , ivars=nsw_max)
  call shdf5_irec(ndims, idims, 'NVE2_MAX', ivars=nve2_max)
  call shdf5_irec(ndims, idims, 'MRLS'    , ivars=mrls)
  call shdf5_irec(ndims, idims, 'NWNUD'   , ivars=nwnud)

  ! Copy grid dimensions

  mza = nza
  mma = nma
  mva = nva
  mua = nua
  mwa = nwa
  mwnud = nwnud

  ! Allocate and read grid structure variables

  call alloc_itabs_pd(nma,nva,nwa)

  call alloc_xyzew(nwa)

  ndims    = 1
  idims(1) = nwa

  call shdf5_irec(ndims, idims, 'XEW', rvar1=xew)
  call shdf5_irec(ndims, idims, 'YEW', rvar1=yew)
  call shdf5_irec(ndims, idims, 'ZEW', rvar1=zew)

  call shdf5_irec(ndims, idims, 'ARW0', rvar1=arw0)

  ! Read ITAB_M SCALARS

  ndims    = 1
  idims(1) = nma

  allocate (iscr1(nma))
  call shdf5_irec(ndims,idims,'itab_m%npoly',ivar1=iscr1)
  do im = 1,nma
     itab_m_pd(im)%npoly = iscr1(im)
  enddo

  call shdf5_irec(ndims,idims,'itab_m%imp'  ,ivar1=iscr1)
  do im = 1,nma
     itab_m_pd(im)%imp = iscr1(im)
  enddo
  deallocate(iscr1)

  ! Read ITAB_M ARRAYS

  ndims    = 2
  idims(2) = nma
  idims(1) = 3

  allocate (iscr2(3,nma))

  call shdf5_irec(ndims,idims,'itab_m%iv',ivar2=iscr2)
  do im = 1,nma
     itab_m_pd(im)%iv(1:3) = iscr2(1:3,im)
  enddo

  call shdf5_irec(ndims,idims,'itab_m%iw',ivar2=iscr2)
  do im = 1,nma
     itab_m_pd(im)%iw(1:3) = iscr2(1:3,im)
  enddo

  call shdf5_irec(ndims,idims,'itab_m%im',ivar2=iscr2)
  do im = 1,nma
     itab_m_pd(im)%im(1:3) = iscr2(1:3,im)
  enddo
  deallocate (iscr2)

  ! Read ITAB_V SCALARS

  ndims    = 1
  idims(1) = nva

  allocate (iscr1(nva))
  call shdf5_irec(ndims,idims,'itab_v%ivp',ivar1=iscr1)
  do iv = 1,nva
     itab_v_pd(iv)%ivp = iscr1(iv)
  enddo
  deallocate(iscr1)

  ! Read ITAB_V ARRAYS

  ndims    = 2
  idims(2) = nva

! idims(1) = 4
! allocate (iscr2(4,nva))
!
! call shdf5_irec(ndims,idims,'itab_v%iw'  ,ivar2=iscr2)
! do iv = 1,nva
!    itab_v_pd(iv)%iw(1:4) = iscr2(1:4,iv)
! enddo
!
! call shdf5_irec(ndims,idims,'itab_v%iv',ivar2=iscr2)
! do iv = 1,nva
!    itab_v_pd(iv)%iv(1:4) = iscr2(1:4,iv)
! enddo
!
! deallocate (iscr2)
!
! idims(1) = 6
! allocate (iscr2(6,nva))
!
! call shdf5_irec(ndims,idims,'itab_v%im',ivar2=iscr2)
! do iv = 1,nva
!    itab_v_pd(iv)%im(1:6) = iscr2(1:6,iv)
! enddo
!
! deallocate (iscr2)

  idims(1) = 2
  allocate (iscr2(2,nva))

  call shdf5_irec(ndims,idims,'itab_v%iw',ivar2=iscr2)
  do iv = 1,nva
     itab_v_pd(iv)%iw(1:2) = iscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_v%im',ivar2=iscr2)
  do iv = 1,nva
     itab_v_pd(iv)%im(1:2) = iscr2(1:2,iv)
  enddo

  deallocate (iscr2)

  ! Read ITAB_W SCALARS

  ndims    = 1
  idims(1) = nwa
  idims(2) = 1

  allocate (iscr1(nwa))
  call shdf5_irec(ndims,idims,'itab_w%npoly'   ,ivar1=iscr1)
  do iw = 1,nwa
     itab_w_pd(iw)%npoly = iscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%iwp'     ,ivar1=iscr1)
  do iw = 1,nwa
     itab_w_pd(iw)%iwp = iscr1(iw)
  enddo
  deallocate(iscr1)

  ! Read ITAB_W ARRAYS

  ndims    = 2
  idims(2) = nwa

  idims(1) = 7
  allocate (iscr2(7,nwa))

  call shdf5_irec(ndims,idims,'itab_w%im',ivar2=iscr2)
  do iw = 1,nwa
     itab_w_pd(iw)%im(1:7) = iscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%iw',ivar2=iscr2)
  do iw = 1,nwa
     itab_w_pd(iw)%iw(1:7) = iscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%iv',ivar2=iscr2)
  do iw = 1,nwa
     itab_w_pd(iw)%iv(1:7) = iscr2(1:7,iw)
  enddo

  deallocate (iscr2)

  ! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     ndims    = 2
     idims(1) = 3
     idims(2) = nwa

     allocate (iscr2(3,mwa))
     call shdf5_irec(ndims,idims,'itab_w%iwnud',ivar2=iscr2)
     do iw = 1,nwa
        itab_w_pd(iw)%iwnud(1:3) = iscr2(1:3,iw)
     enddo
     deallocate (iscr2)

     call alloc_nudge1(nwnud,1)

     ndims    = 1
     idims(1) = nwnud

     call shdf5_irec(ndims, idims, 'XEWNUD', rvar1=xewnud)
     call shdf5_irec(ndims, idims, 'YEWNUD', rvar1=yewnud)
     call shdf5_irec(ndims, idims, 'ZEWNUD', rvar1=zewnud)

  endif

  ! If SFCGRID is active, read ATM-SFCGRID coupling information

  if (mdomain <= 1) then

     ! Read the grid dimensions

     ndims = 1
     idims(1) = 1

     call shdf5_irec(ndims, idims, 'nwsfc'  , ivars=nwsfc)

   !!  allocate (itab_wsfc_pd(nwsfc))

     ! Read ATM-SFCGRID coupling information

     idims(1) = nwsfc

     allocate (iscr1(nwsfc))
     call shdf5_irec(ndims,idims,'itab_wsfc%nwatm'     ,ivar1=iscr1)
     do iwsfc = 1,nwsfc
        itab_wsfc_pd(iwsfc)%nwatm = iscr1(iwsfc)
     enddo
     deallocate(iscr1)

     ndims    = 2
     idims(1) = 8
     idims(2) = nwsfc

     allocate (iscr2(8,nwsfc))
     call shdf5_irec(ndims,idims,'itab_wsfc%iwatm',ivar2=iscr2)
     do iwsfc = 1,nwsfc
        itab_wsfc_pd(iwsfc)%iwatm(1:8) = iscr2(1:8,iwsfc)
     enddo
     deallocate(iscr2)

  endif

  ! Close the GRIDFILE

  call shdf5_close()

  write(io6,*)
  write(io6,*) 'end of gridfile_read_pd '

end subroutine gridfile_read_pd

!===============================================================================

subroutine gridfile_read()

  use max_dims,   only: pathlen
  use misc_coms,  only: io6, gridfile, mdomain
  use mem_ijtabs, only: mloops,itab_m, itab_v, itab_w
  use mem_grid,   only: nza, mma, mva, mwa, &
                        zm, zt, dzm, dzt, dzim, dzit, &
                        zfacm, zfact, zfacim, zfacit, zfacm2, zfacim2, &
                        lpm, lpv, lpw, lsw, lve2, topm, topw, &
                        xem, yem, zem, xev, yev, zev, &
                        unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                        dnu, dniu, dnv, dniv, arm0, &
                        glatw, glonw, glatm, glonm, glatv, glonv, &
                        arv, arw, volt
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_nudge,  only: nudflag, nudnxp

  ! This subroutine checks for the existence of a gridfile, and if it exists,
  ! also checks for agreement of grid configuration between the file and the
  ! current model run.  If the file does not exist or does not match grid
  ! configuration, the run is stopped.

  implicit none

  integer :: im, iv, iw
  integer :: ndims, idims(2)
  logical :: exans

  character(2) :: type

  ! Scratch arrays for copying input

  integer, allocatable :: iscr1(:)
  real,    allocatable :: rscr1(:)
  real,    allocatable :: rscr2(:,:)
  logical, allocatable :: lscr2(:,:)

  ! Pointers to the global index of the local point

  integer :: lgma(mma)
  integer :: lgva(mva)
  integer :: lgwa(mwa)

  character(pathlen) :: flnm

  lgma   = itab_m(1:mma)%imglobe
  lgva   = itab_v(1:mva)%ivglobe
  lgwa   = itab_w(1:mwa)%iwglobe

  ! Check if gridfile exists

  flnm = trim(gridfile)//'.h5'

  inquire(file=flnm, exist=exans)

  if (.not. exans) then

     ! Gridfile does not exist.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!  Gridfile does not exist:'
     write(io6,*) '!!!  '//flnm
     write(io6,*) '!!!  Stopping run'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

     stop 'stop - no gridfile'

  endif

  ! Gridfile exists; open it

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'Opening grid file ', trim(flnm)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R',trypario=.true.)

  ndims    = 1
  idims(1) = nza

  call shdf5_irec(ndims, idims, 'ZM'    , rvar1=zm)
  call shdf5_irec(ndims, idims, 'ZT'    , rvar1=zt)
  call shdf5_irec(ndims, idims, 'DZM'   , rvar1=dzm)
  call shdf5_irec(ndims, idims, 'DZT'   , rvar1=dzt)
  call shdf5_irec(ndims, idims, 'DZIM'  , rvar1=dzim)
  call shdf5_irec(ndims, idims, 'DZIT'  , rvar1=dzit)
  call shdf5_irec(ndims, idims, 'ZFACM' , rvar1=zfacm)
  call shdf5_irec(ndims, idims, 'ZFACT' , rvar1=zfact)
  call shdf5_irec(ndims, idims, 'ZFACIM', rvar1=zfacim)
  call shdf5_irec(ndims, idims, 'ZFACIT', rvar1=zfacit)
  call shdf5_irec(ndims, idims, 'ZFACM2' , rvar1=zfacm2)
  call shdf5_irec(ndims, idims, 'ZFACIM2', rvar1=zfacim2)

  ndims    = 1
  idims(1) = mma
  type     = 'AM'

  call shdf5_irec(ndims, idims, 'LPM'  , ivar1=lpm  , points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'TOPM' , rvar1=topm , points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'ARM0' , rvar1=arm0 , points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLATM', rvar1=glatm, points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLONM', rvar1=glonm, points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'XEM'  , rvar1=xem  , points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'YEM'  , rvar1=yem  , points=lgma, stagpt=type)
  call shdf5_irec(ndims, idims, 'ZEM'  , rvar1=zem  , points=lgma, stagpt=type)

  ndims    = 1
  idims(1) = mva
  type     = 'AV'

  call shdf5_irec(ndims, idims, 'LPV'  , ivar1=lpv  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'XEV'  , rvar1=xev  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'YEV'  , rvar1=yev  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'ZEV'  , rvar1=zev  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLATV', rvar1=glatv, points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLONV', rvar1=glonv, points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'DNU'  , rvar1=dnu  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'DNV'  , rvar1=dnv  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'DNIU' , rvar1=dniu , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'DNIV' , rvar1=dniv , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'UNX'  , rvar1=unx  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'UNY'  , rvar1=uny  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'UNZ'  , rvar1=unz  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'VNX'  , rvar1=vnx  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'VNY'  , rvar1=vny  , points=lgva, stagpt=type)
  call shdf5_irec(ndims, idims, 'VNZ'  , rvar1=vnz  , points=lgva, stagpt=type)

  ndims    = 1
  idims(1) = mwa
  type     = 'AW'

  call shdf5_irec(ndims, idims, 'LPW'  , ivar1=lpw  , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'LSW'  , ivar1=lsw  , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'LVE2' , ivar1=lve2 , points=lgwa, stagpt=type)
! call shdf5_irec(ndims, idims, 'XEW'  , rvar1=xew  , points=lgwa, stagpt=type)
! call shdf5_irec(ndims, idims, 'YEW'  , rvar1=yew  , points=lgwa, stagpt=type)
! call shdf5_irec(ndims, idims, 'ZEW'  , rvar1=zew  , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'WNX'  , rvar1=wnx  , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'WNY'  , rvar1=wny  , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'WNZ'  , rvar1=wnz  , points=lgwa, stagpt=type)
! call shdf5_irec(ndims, idims, 'ARW0' , rvar1=arw0 , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'TOPW' , rvar1=topw , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLATW', rvar1=glatw, points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'GLONW', rvar1=glonw, points=lgwa, stagpt=type)

  ndims    = 2
  idims(1) = nza
  idims(2) = mva
  type     = 'AV'

  call shdf5_irec(ndims, idims, 'ARV', rvar2=arv, points=lgva, stagpt=type)

  ndims    = 2
  idims(1) = nza
  idims(2) = mwa
  type     = 'AW'

  call shdf5_irec(ndims, idims, 'ARW' , rvar2=arw , points=lgwa, stagpt=type)
  call shdf5_irec(ndims, idims, 'VOLT', dvar2=volt, points=lgwa, stagpt=type)

  ! Read ITAB_M SCALARS

  ndims    = 1
  idims(1) = mma
  type     = 'AM'

  allocate (iscr1(mma))

  call shdf5_irec(ndims, idims, 'itab_m%mrlm'     ,ivar1=iscr1, points=lgma, stagpt=type)
  do im = 1,mma
     itab_m(im)%mrlm = iscr1(im)
  enddo

  call shdf5_irec(ndims, idims, 'itab_m%mrlm_orig',ivar1=iscr1, points=lgma, stagpt=type)
  do im = 1,mma
     itab_m(im)%mrlm_orig = iscr1(im)
  enddo

  call shdf5_irec(ndims, idims, 'itab_m%mrow'     ,ivar1=iscr1, points=lgma, stagpt=type)
  do im = 1,mma
     itab_m(im)%mrow = iscr1(im)
  enddo

  call shdf5_irec(ndims, idims, 'itab_m%ngr'      ,ivar1=iscr1, points=lgma, stagpt=type)
  do im = 1,mma
     itab_m(im)%ngr = iscr1(im)
  enddo
  deallocate(iscr1)

  ! Read ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = mma
  type     = 'AM'

  allocate (lscr2(mloops,mma))
  call shdf5_irec(ndims, idims, 'itab_m%loop', lvar2=lscr2, points=lgma, stagpt=type)
  do im = 1,mma
     itab_m(im)%loop(1:mloops) = lscr2(1:mloops,im)
  enddo
  deallocate (lscr2)

  ! Read ITAB_V SCALARS

  ndims    = 1
  idims(1) = mva
  type     = 'AV'

  allocate (iscr1(mva))

  call shdf5_irec(ndims,idims,'itab_v%mrlv'  , ivar1=iscr1, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%mrlv = iscr1(iv)
  enddo
  deallocate (iscr1)

  ! Read ITAB_V ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = mva
  type     = 'AV'

  allocate (lscr2(mloops,mva))
  call shdf5_irec(ndims,idims,'itab_v%loop',lvar2=lscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%loop(1:mloops) = lscr2(1:mloops,iv)
  enddo
  deallocate (lscr2)

  idims(1) = 2
  allocate (rscr2(2,mva))

  call shdf5_irec(ndims,idims,'itab_v%farw',rvar2=rscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%farw(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_v%cosv',rvar2=rscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%cosv(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_v%sinv',rvar2=rscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%sinv(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_v%dxps',rvar2=rscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%dxps(1:2) = rscr2(1:2,iv)
  enddo

  call shdf5_irec(ndims,idims,'itab_v%dyps',rvar2=rscr2, points=lgva, stagpt=type)
  do iv = 1,mva
     itab_v(iv)%dyps(1:2) = rscr2(1:2,iv)
  enddo

  deallocate (rscr2)

  ! Read ITAB_W SCALARS

  ndims    = 1
  idims(1) = mwa
  type     = 'AW'

  allocate (iscr1(mwa))

  call shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivar1=iscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%mrlw = iscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivar1=iscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%mrlw_orig = iscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ngr'      ,ivar1=iscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%ngr = iscr1(iw)
  enddo
  deallocate (iscr1)

  allocate (rscr1(mwa))
  call shdf5_irec(ndims,idims,'itab_w%unx_w' ,rvar1=rscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%unx_w = rscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%uny_w' ,rvar1=rscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%uny_w = rscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%vnx_w' ,rvar1=rscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%vnx_w = rscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%vny_w' ,rvar1=rscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%vny_w = rscr1(iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%vnz_w' ,rvar1=rscr1, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%vnz_w = rscr1(iw)
  enddo
  deallocate (rscr1)

  ! Read ITAB_W ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = mwa
  type     = 'AW'

  allocate (lscr2(mloops,mwa))
  call shdf5_irec(ndims,idims,'itab_w%loop',lvar2=lscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%loop(1:mloops) = lscr2(1:mloops,iw)
  enddo
  deallocate (lscr2)

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     idims(1) = 3
     allocate (rscr2(3,mwa))
     call shdf5_irec(ndims,idims,'itab_w%fnud',rvar2=rscr2, points=lgwa, stagpt=type)
     do iw = 1,mwa
        itab_w(iw)%fnud(1:3) = rscr2(1:3,iw)
     enddo
     deallocate (rscr2)

  endif

  idims(1) = 7
  allocate (rscr2(7,mwa))

  call shdf5_irec(ndims,idims,'itab_w%dirv',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%dirv(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%farm',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%farm(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%farv',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%farv(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%gxps1',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%gxps1(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%gyps1',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%gyps1(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%gxps2',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%gxps2(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%gyps2',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%gyps2(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vx',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vx(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vy',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vy(1:7) = rscr2(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vz',rvar2=rscr2, points=lgwa, stagpt=type)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vz(1:7) = rscr2(1:7,iw)
  enddo

  deallocate (rscr2)

  ! Close the GRIDFILE

  call shdf5_close()

  write(io6,*) 'end of gridfile_read '

end subroutine gridfile_read

!===============================================================================

subroutine gridfile_read_sfc()

  use max_dims,   only: pathlen
  use misc_coms,  only: io6, gridfile, mdomain
  use mem_sfcg,   only: sfcg, mwsfc, itab_wsfc
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close

  ! This subroutine checks for the existence of a gridfile, and if it exists,
  ! also checks for agreement of grid configuration between the file and the
  ! current model run.  If the file does not exist or does not match grid
  ! configuration, the run is stopped.

  implicit none

  integer      :: iw, ndims, idims(2)
  logical      :: exans
  character(2) :: type

  ! Scratch arrays for copying input

  integer, allocatable :: iscr2(:,:)
  real,    allocatable :: rscr2(:,:)

  ! Pointers to the global index of the local point

  integer :: lgwsfc(mwsfc)

  character(pathlen) :: flnm

  ! If SFCGRID is active, read ATM-SFCGRID coupling information

  if (mdomain > 1) return

  ! Check if gridfile exists

  flnm = trim(gridfile)//'.h5'

  inquire(file=flnm, exist=exans)

  if (.not. exans) then

     ! Gridfile does not exist.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!  Gridfile does not exist:'
     write(io6,*) '!!!  '//flnm
     write(io6,*) '!!!  Stopping run'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

     stop 'stop - no gridfile'

  endif

  ! Gridfile exists; open it

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'Opening grid file ', trim(flnm)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R',trypario=.true.)

  lgwsfc = itab_wsfc(1:mwsfc)%iwglobe

  ndims    = 1
  idims(1) = mwsfc
  type     = 'CW'

  call shdf5_irec(ndims, idims, 'dzt_bot' , rvar1=sfcg%dzt_bot, points=lgwsfc, stagpt=type)

  ndims    = 2
  idims(1) = 8
  idims(2) = mwsfc
  type     = 'CW'

  allocate (iscr2(8,mwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%kwatm',ivar2=iscr2, points=lgwsfc, stagpt=type)
  do iw = 1,mwsfc
     itab_wsfc(iw)%kwatm(1:8) = iscr2(1:8,iw)
  enddo
  deallocate(iscr2)

  allocate (rscr2(8,mwsfc))
  call shdf5_irec(ndims,idims,'itab_wsfc%arc',rvar2=rscr2, points=lgwsfc, stagpt=type)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arc(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoarsfc',rvar2=rscr2, points=lgwsfc, stagpt=type)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoarsfc(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoariw',rvar2=rscr2, points=lgwsfc, stagpt=type)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoariw(1:8) = rscr2(1:8,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_wsfc%arcoarkw',rvar2=rscr2, points=lgwsfc, stagpt=type)
  do iw = 1,mwsfc
     itab_wsfc(iw)%arcoarkw(1:8) = rscr2(1:8,iw)
  enddo
  deallocate(rscr2)

  ! Close the GRIDFILE

  call shdf5_close()

  write(io6,*) 'end of gridfile_read_sfc'

end subroutine gridfile_read_sfc

!===============================================================================

subroutine gridfile_read_nudge()

  use max_dims,   only: pathlen
  use misc_coms,  only: io6, gridfile, mdomain
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_nudge,  only: nudflag, nudnxp, mwnud, itab_wnud

  ! This subroutine checks for the existence of a gridfile, and if it exists,
  ! also checks for agreement of grid configuration between the file and the
  ! current model run.  If the file does not exist or does not match grid
  ! configuration, the run is stopped.

  implicit none

  integer :: iwnud
  integer :: ndims, idims(2)
  logical :: exans

  character(2) :: type

  ! Scratch arrays for copying input

  integer, allocatable :: iscr1 (:)
  integer, allocatable :: iscr2 (:,:)
  integer, allocatable :: lgwnud(:)

  character(pathlen) :: flnm

  ! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     ! Check if gridfile exists

     flnm = trim(gridfile)//'.h5'

     inquire(file=flnm, exist=exans)

     if (.not. exans) then

        ! Gridfile does not exist.

        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(io6,*) '!!!  Gridfile does not exist:'
        write(io6,*) '!!!  '//flnm
        write(io6,*) '!!!  Stopping run'
        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

        stop 'stop - no gridfile'

     endif

     ! Gridfile exists; open it

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening nuding grid file ', trim(flnm)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(flnm,'R',trypario=.true.)

     allocate(lgwnud(mwnud))

     ndims    = 1
     idims(1) = mwnud
     type     = 'AN'

     lgwnud = itab_wnud(1:mwnud)%iwnudglobe

     allocate (iscr1(mwnud))
     call shdf5_irec(ndims, idims, 'itab_wnud%npoly', ivar1=iscr1, points=lgwnud, stagpt=type)

     do iwnud = 1,mwnud
        itab_wnud(iwnud)%npoly = iscr1(iwnud)
     enddo
     deallocate(iscr1)

     ndims    = 2
     idims(1) = 6
     idims(2) = mwnud
     type     = 'AN'

     allocate (iscr2(6,mwnud))
     call shdf5_irec(ndims,idims,'itab_wnud%iwnud',ivar2=iscr2, points=lgwnud, stagpt=type)

     do iwnud = 1,mwnud
        itab_wnud(iwnud)%iwnud(1:6) = iscr2(1:6,iwnud)
     enddo
     deallocate(iscr2)

  endif

  ! Close the GRIDFILE

  call shdf5_close()

  write(io6,*) 'end of gridfile_read_nudge'

end subroutine gridfile_read_nudge

!===============================================================================

! This subroutine reads a few quantities from the OLD gridfile for a HISTREGRID
! start, i.e., with a modified horizontal structure of the ATM grid.

subroutine gridfile_read_oldgrid()

  use max_dims,   only: maxngrdll, pathlen
  use misc_coms,  only: io6, gridfile, mdomain, nzp, ndz, hdz, dz
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_regrid, only: nzp_og, mdomain_og, ndz_og, hdz_og, dz_og, &
                        nza_og, nwa_og, xew_og, yew_og, zew_og, &
                        wnx_og, wny_og, wnz_og, glatw_og, glonw_og, &
                        lpw_og
  use ll_bins,    only: itab_w0

  implicit none

  integer :: idz, iw_og
  integer :: ierr
  integer :: ndims, idims(2)
  logical :: exans

  ! Scratch array for copying input

  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  character(pathlen) :: flnm

  ! Check if OLD gridfile exists

  flnm = trim(gridfile)//'-OG'//'.h5'

  inquire(file=flnm, exist=exans)

  if (.not. exans) then

  ! OLD gridfile does not exist.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!  OLD GRIDFILE does not exist:'
     write(io6,*) '!!!  '//flnm
     write(io6,*) '!!!  Stopping run'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

     stop 'stop - no old gridfile'

  endif

  ! Gridfile exists; open it

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'Opening OLD GRIDFILE ', flnm
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R',trypario=.true.)

  ! Read the grid information that exists in namelist

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp_og)
  call shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain_og)
  call shdf5_irec(ndims, idims, 'NDZ'     , ivars=ndz_og)

  allocate( hdz_og(ndz_og))
  allocate( dz_og (ndz_og))

  if (ndz_og > 1) then
     idims(1) = ndz_og

     call shdf5_irec(ndims, idims, 'HDZ' , rvar1=hdz_og)
     call shdf5_irec(ndims, idims, 'DZ'  , rvar1=dz_og)
  endif

  ! Check equality between grid file information and namelist variables

  ierr = 0

  if (nzp_og      /= nzp     ) ierr = 1
  if (mdomain_og  /= mdomain ) ierr = 1
  if (ndz_og      /= ndz     ) ierr = 1

  do idz = 1, min(ndz_og,ndz)
     if (abs(hdz_og(idz) - hdz(idz)) > 1.e1 ) ierr = 1
     if (abs(dz_og (idz) - dz (idz)) > 1.e1 ) ierr = 1
  enddo

  if (ierr == 1) then

     write(io6,*) 'OLDGRIDFILE mismatch with OLAMIN namelist: Stopping model run'
     write(io6,*) 'Values: oldgridfile, namelist'
     write(io6,*) '-----------------------------------------------'
     write(io6,*)              'nzp:      ',nzp_og     ,nzp
     write(io6,*)              'mdomain:  ',mdomain_og ,mdomain
     write(io6,*)              'ndz:      ',ndz_og     ,ndz
     write(io6,*) ' '
     write(io6, '(a,20f10.1)') 'hdz_og:   ',hdz_og(1:ndz_og)
     write(io6, '(a,20f10.1)') 'hdz:      ',hdz   (1:ndz_og)
     write(io6,*) ' '
     write(io6, '(a,20f10.1)') 'dz_og:    ',dz_og(1:ndz_og)
     write(io6, '(a,20f10.1)') 'dz:       ',dz   (1:ndz_og)
     write(io6,*) ' '
     write(io6,*) '-----------------------------------------------'

     stop 'stop - gridfile mismatch'

  endif

  ! Read the grid dimensions

  ndims    = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'NZA'     , ivars=nza_og)
  call shdf5_irec(ndims, idims, 'NWA'     , ivars=nwa_og)

  allocate (lpw_og(nwa_og))

  allocate (xew_og(nwa_og))
  allocate (yew_og(nwa_og))
  allocate (zew_og(nwa_og))

  allocate (wnx_og(nwa_og))
  allocate (wny_og(nwa_og))
  allocate (wnz_og(nwa_og))

  allocate (glatw_og(nwa_og))
  allocate (glonw_og(nwa_og))

  allocate (itab_w0(nwa_og))

  idims(1) = nwa_og

  call shdf5_irec(ndims, idims, 'LPW',  ivar1=lpw_og)

  call shdf5_irec(ndims, idims, 'XEW', rvar1=xew_og)
  call shdf5_irec(ndims, idims, 'YEW', rvar1=yew_og)
  call shdf5_irec(ndims, idims, 'ZEW', rvar1=zew_og)

  call shdf5_irec(ndims, idims, 'WNX', rvar1=wnx_og)
  call shdf5_irec(ndims, idims, 'WNY', rvar1=wny_og)
  call shdf5_irec(ndims, idims, 'WNZ', rvar1=wnz_og)

  call shdf5_irec(ndims, idims, 'GLATW', rvar1=glatw_og)
  call shdf5_irec(ndims, idims, 'GLONW', rvar1=glonw_og)

  allocate (iscr1(nwa_og))
  call shdf5_irec(ndims,idims,'itab_w%npoly',ivar1=iscr1)
  do iw_og = 1,nwa_og
     itab_w0(iw_og)%npoly = iscr1(iw_og)
  enddo
  deallocate(iscr1)

  ndims    = 2
  idims(1) = 7
  idims(2) = nwa_og

  allocate (iscr2(7,nwa_og))

  call shdf5_irec(ndims,idims,'itab_w%iw',ivar2=iscr2)
  do iw_og = 1,nwa_og
     itab_w0(iw_og)%iw(1:7) = iscr2(1:7,iw_og)
  enddo

  deallocate (iscr2)

  call shdf5_close()

  write(io6,*) 'end of gridfile_read_oldgrid '

end subroutine gridfile_read_oldgrid
