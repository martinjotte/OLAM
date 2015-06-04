!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================
subroutine gridinit()

use misc_coms,   only: io6, mdomain, ngrids, initial, nxp, nzp, &
                       timmax8, alloc_misc, &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim

use leaf_coms,   only: nwl
use sea_coms,    only: nws

use mem_ijtabs,  only: istp, mrls, fill_jtabs, itab_md, itab_ud, itab_wd, itab_w
use oplot_coms,  only: op
use mem_grid,    only: nza, nma, nua, nva, nwa, &
                       zm, zt, dzt, xem, yem, zem, xew, yew, zew, glatw, glonw, &
                       alloc_grid1, alloc_grid2
use mem_nudge,   only: nudflag, nudnxp, nwnud, itab_wnud, alloc_nudge1, &
                       xewnud, yewnud, zewnud

implicit none

integer :: npoly
integer :: j, jmaxneg, jminpos
integer :: imd,imd1,imd2,iud,iwd,iwnud,iwnudn,iw

real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
real :: dist_wnud,dist,dist_min,xi,yi
real :: xin(6),yin(6)

! Vertical grid coordinate setup

  write(io6,'(/,a)') 'gridinit calling gridset2'
  call gridset2()

  write(io6,'(a,f8.1)') ' Model top height = ',zm(nza)

! Horizontal grid setup

  if (mdomain == 0) then

! If using nudging on global domain with independent nudging grid,
! generate nudging grid here

     if (nudflag > 0 .and. nudnxp > 0) then

        write(io6,'(/,a)') 'gridinit calling icosahedron for nudging grid'

        call icosahedron(nudnxp)  ! global spherical domain; calls 2 allocs

! Set dimensions of nudging grid arrays, which are the Voronoi dual-grid
! counterpart to the Delaunay triangle mesh just generated in icosahedron
! and dimensioned by nwa and nma

        nwnud = nma

! Allocate nudging grid structure arrays, copy temporary grid structure to
! these arrays, and deallocate temporary arrays

        call alloc_nudge1(nwnud)

! Loop over nudging grid W points

        do iwnud = 2,nwnud
           imd = iwnud

! Set nudging grid W point npoly value and earth coordinates identical to
! triangle mesh MD points

           itab_wnud(iwnud)%npoly = itab_md(imd)%npoly

           xewnud(iwnud) = xem(imd)
           yewnud(iwnud) = yem(imd)
           zewnud(iwnud) = zem(imd)

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
        deallocate (xem, yem, zem)

        write(io6,'(/,a)') 'gridinit after nudging grid construction'
        write(io6,'(a,i8)')    ' nwnud = ',nwnud

     endif

! Now generate global atmospheric grid

     write(io6,'(/,a)') 'gridinit calling icosahedron'
     call icosahedron(nxp)  ! global spherical domain; calls 2 allocs

  elseif (mdomain == 2 .or. mdomain == 3) then

     write(io6,'(/,a)') 'gridinit calling cartesian_2d'
     call cartesian_2d()    ! 2D cartesian channel domain; calls 2 allocs

  elseif (mdomain == 4) then

     write(io6,'(/,a)') 'gridinit calling cartesian_3d'
     call cart4_hex()    ! 3D cartesian channel domain; calls 2 allocs

  elseif (mdomain == 5) then

     write(io6,'(/,a)') 'gridinit calling cart_hex'
     call cart_hex()        ! 3D hexagonal domain with hexagon cells and
                            ! periodic lateral boundaries; calls 2 allocs

  endif

  write(io6,'(/,a)') 'gridinit after icosahedron or cartesian'
  write(io6,'(a,i8)')    ' nma = ',nma
  write(io6,'(a,i8)')    ' nua = ',nua
  write(io6,'(a,i8)')    ' nwa = ',nwa

  if (ngrids > 1 .and. mdomain /= 2 .and. mdomain /= 3) then

     write(io6,'(/,a,i5)') 'gridinit calling spawn_nest; ngrids = ',ngrids

     call spawn_nest()  ! calls 2 allocs

     write(io6,'(/,a)') 'gridinit after spawn_nest'
     write(io6,'(a,i8)')   ' nma = ',nma
     write(io6,'(a,i8)')   ' nua = ',nua
     write(io6,'(a,i8)')   ' nwa = ',nwa

  endif

  if (mdomain /= 4) then

     call voronoi()
     call pcvt()

  endif

! Allocate remaining GRID FOOTPRINT arrays for full domain

  write(io6,'(/,a)') 'gridinit calling alloc_grid1 for full domain'

  call alloc_grid1(nma, nva, nwa)

! Initialize dtlm, dtsm, ndtrat, and nacoust,
! and compute the timestep schedule for all grid operations.

  write(io6,'(/,a)') 'gridinit calling modsched'

  call modsched()

  write(io6,'(/,a)') 'gridinit calling fill_jtabs'

  call fill_jtabs(nma,nva,nwa,0)

! Fill remaining GRID FOOTPRINT geometry for full domain

  write(io6,'(/,a)') 'gridinit calling grid_geometry'

  call grid_geometry_hex()

! Set up topographic surface and its intersection of the 3D atmospheric grid.

  write(io6,'(/,a)') 'gridinit calling makesfc2'

  call makesfc2()

! Allocate remaining unstructured grid geometry arrays

  write(io6,'(/,a)') 'gridinit calling alloc_grid2'

  call alloc_grid2(nma, nva, nwa)

! Set up control volumes in atmospheric grid

  write(io6,'(/,a)') 'gridinit calling ctrlvols'

  call ctrlvols_hex()

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

! Write GRIDFILE, LANDFILE, and SEAFILE

  write(io6,'(/,a)') 'gridinit calling gridfile_write'
  call gridfile_write()

  write(io6,'(/,a)') 'calling landfile_write'
  call landfile_write()

  write(io6,'(/,a)') 'calling seafile_write'
  call seafile_write()

  write(io6,'(/,a)') 'gridinit completed'

end subroutine gridinit

!===============================================================================

subroutine gridset2()

use mem_grid,    only: nza, mza, &
                       zm, zt, dzm, dzt, dzim, dzit, &
                       zfacm, zfact, zfacim, zfacit, &
                       alloc_gridz
use misc_coms,   only: io6, nzp, ndz, hdz, dz, mdomain
use consts_coms, only: erad
use oname_coms,  only: nl

implicit none

integer :: idz, kvec, nseries, iseries, k

real :: zend, dzend, dzbeg, ztarg, dzr, rdzr

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

! Fill NZA and MZA values

   nza = kvec
   mza = nza

   if (nza < 2) then
      write(io6,'(a)')      'OLAM requires that NZA >= 2.'
      write(io6,'(a,i5,a)') 'NZA = ',nza,' in subroutine gridset2.  Stopping model.'
      stop 'stop nza in gridset2'
   endif

! Fill top 2 ZMVEC and ZTVEC values

   zmvec(nza)   = zmvec(nza-1) * 2. - zmvec(nza-2)
   zmvec(nza+1) = zmvec(nza)   * 2. - zmvec(nza-1)

   ztvec(nza)   = .5 * (zmvec(nza-1) + zmvec(nza))
   ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

endif

! Allocate main grid arrays

call alloc_gridz()

! Other vertical coordinate values

do k = 1,nza
   zm(k) = zmvec(k)
   zt(k) = ztvec(k)

   dzm(k) = ztvec(k+1) - ztvec(k)
   dzt(k) = zmvec(k) - zmvec(k-1)

   dzim(k) = 1. / dzm(k)
   dzit(k) = 1. / dzt(k)

   if (mdomain < 2) then
      zfacm(k) = (erad + zm(k)) / erad
      zfact(k) = (erad + zt(k)) / erad
   else
      zfacm(k) = 1.
      zfact(k) = 1.
   endif

   zfacim(k) = 1. / zfacm(k)
   zfacit(k) = 1. / zfact(k)
enddo

deallocate (zmvec,ztvec)

end subroutine gridset2

!===============================================================================

subroutine gridset_print()

use misc_coms, only: io6
use mem_grid,  only: nza, zm, zt, dzt, nma, nua, nva, nwa
use leaf_coms, only: nwl, isfcl
use sea_coms,  only: nws

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
write(io6,'(a,i8)') '  nma = ',nma
write(io6,'(a,i8)') '  nua = ',nua
write(io6,'(a,i8)') '  nva = ',nva
write(io6,'(a,i8)') '  nwa = ',nwa
write(io6,'(a,i8)') '  nza = ',nza

if (isfcl == 1) then
   write(io6,'(/,a)' ) 'Global land/sea indices:'
   write(io6,'(a,i8)') '  nwl       = ',nwl
   write(io6,'(a,i8)') '  nws       = ',nws
endif

end subroutine gridset_print

!===============================================================================

subroutine gridfile_write()

  use max_dims,   only: maxngrdll, pathlen
  use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
       iclobber, itopoflg, deltax, ndz, hdz, dz, ngrdll, grdrad, grdlat, grdlon
  use mem_ijtabs, only: mloops, mrls, &
       itab_m, itab_v, itab_w
  use mem_grid,   only: nza, nma, nua, nva, nwa, nsw_max, &
       zm, zt, dzm, dzt, dzim, dzit, &
       zfacm, zfact, zfacim, zfacit, &
       lpm, lpv, lpw, lsw, &
       topm, topw, xem, yem, zem, &
       xev, yev, zev, xew, yew, zew, &
       unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
       dnu, dniu, dnv, dniv, arw0, arm0, &
       glatw, glonw, glatm, glonm, glatv, glonv, &
       arv, arw, volt, volti
  use leaf_coms,  only: isfcl

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

  use mem_nudge,  only: nudflag, nudnxp, nwnud, itab_wnud, &
                        xewnud, yewnud, zewnud

  implicit none

  ! This routine writes the grid variables to the gridfile.

  integer :: im, iv, iw, iwnud

  integer :: ndims, idims(2)

  ! Scratch arrays for copying output

  logical, allocatable :: lscr(:,:)
  integer, allocatable :: iscr(:,:)

  real, allocatable :: rscr(:,:)

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

     call shdf5_orec(ndims, idims, 'HDZ'  , rvara=hdz)
     call shdf5_orec(ndims, idims, 'DZ'   , rvara=dz)
  endif

  ndims    = 1
  idims(1) = ngrids

  call shdf5_orec(ndims, idims, 'NGRDLL' , ivara=ngrdll)
  call shdf5_orec(ndims, idims, 'GRDRAD' , rvara=grdrad)

  ndims    = 2
  idims(1) = ngrids
  idims(2) = maxngrdll

  call shdf5_orec(ndims, idims, 'GRDLAT', rvara=grdlat(1:ngrids,:))
  call shdf5_orec(ndims, idims, 'GRDLON', rvara=grdlon(1:ngrids,:))

  ! Write the grid dimensions

  ndims    = 1
  idims(1) = 1

  call shdf5_orec(ndims, idims, 'NZA'     , ivars=nza)
  call shdf5_orec(ndims, idims, 'NMA'     , ivars=nma)
  call shdf5_orec(ndims, idims, 'NUA'     , ivars=nua)
  call shdf5_orec(ndims, idims, 'NVA'     , ivars=nva)
  call shdf5_orec(ndims, idims, 'NWA'     , ivars=nwa)
  call shdf5_orec(ndims, idims, 'NSW_MAX' , ivars=nsw_max)
  call shdf5_orec(ndims, idims, 'MRLS'    , ivars=mrls)
  call shdf5_orec(ndims, idims, 'NWNUD'   , ivars=nwnud)

  ! Write grid structure variables

  ndims    = 1
  idims(1) = nza

  call shdf5_orec(ndims, idims, 'ZM'    , rvara=zm)
  call shdf5_orec(ndims, idims, 'ZT'    , rvara=zt)
  call shdf5_orec(ndims, idims, 'DZM'   , rvara=dzm)
  call shdf5_orec(ndims, idims, 'DZT'   , rvara=dzt)
  call shdf5_orec(ndims, idims, 'DZIM'  , rvara=dzim)
  call shdf5_orec(ndims, idims, 'DZIT'  , rvara=dzit)
  call shdf5_orec(ndims, idims, 'ZFACM' , rvara=zfacm)
  call shdf5_orec(ndims, idims, 'ZFACT' , rvara=zfact)
  call shdf5_orec(ndims, idims, 'ZFACIM', rvara=zfacim)
  call shdf5_orec(ndims, idims, 'ZFACIT', rvara=zfacit)

  ndims    = 1
  idims(1) = nma

  call shdf5_orec(ndims, idims, 'LPM'  , ivara=lpm)
  call shdf5_orec(ndims, idims, 'TOPM' , rvara=topm)
  call shdf5_orec(ndims, idims, 'XEM'  , rvara=xem)
  call shdf5_orec(ndims, idims, 'YEM'  , rvara=yem)
  call shdf5_orec(ndims, idims, 'ZEM'  , rvara=zem)
  call shdf5_orec(ndims, idims, 'ARM0' , rvara=arm0)
  call shdf5_orec(ndims, idims, 'GLATM', rvara=glatm)
  call shdf5_orec(ndims, idims, 'GLONM', rvara=glonm)

  ndims    = 1
  idims(1) = nva

  call shdf5_orec(ndims, idims, 'LPV'  , ivara=lpv)
  call shdf5_orec(ndims, idims, 'XEV'  , rvara=xev)
  call shdf5_orec(ndims, idims, 'YEV'  , rvara=yev)
  call shdf5_orec(ndims, idims, 'ZEV'  , rvara=zev)
  call shdf5_orec(ndims, idims, 'GLATV', rvara=glatv)
  call shdf5_orec(ndims, idims, 'GLONV', rvara=glonv)

  call shdf5_orec(ndims, idims, 'DNU' , rvara=dnu)
  call shdf5_orec(ndims, idims, 'DNV' , rvara=dnv)
  call shdf5_orec(ndims, idims, 'DNIU', rvara=dniu)
  call shdf5_orec(ndims, idims, 'DNIV', rvara=dniv)
  call shdf5_orec(ndims, idims, 'UNX' , rvara=unx)
  call shdf5_orec(ndims, idims, 'UNY' , rvara=uny)
  call shdf5_orec(ndims, idims, 'UNZ' , rvara=unz)
  call shdf5_orec(ndims, idims, 'VNX' , rvara=vnx)
  call shdf5_orec(ndims, idims, 'VNY' , rvara=vny)
  call shdf5_orec(ndims, idims, 'VNZ' , rvara=vnz)

  ndims    = 1
  idims(1) = nwa

  call shdf5_orec(ndims, idims, 'LPW'  , ivara=lpw)
  call shdf5_orec(ndims, idims, 'LSW'  , ivara=lsw)
  call shdf5_orec(ndims, idims, 'XEW'  , rvara=xew)
  call shdf5_orec(ndims, idims, 'YEW'  , rvara=yew)
  call shdf5_orec(ndims, idims, 'ZEW'  , rvara=zew)
  call shdf5_orec(ndims, idims, 'TOPW' , rvara=topw)
  call shdf5_orec(ndims, idims, 'ARW0' , rvara=arw0)
  call shdf5_orec(ndims, idims, 'GLATW', rvara=glatw)
  call shdf5_orec(ndims, idims, 'GLONW', rvara=glonw)
  call shdf5_orec(ndims, idims, 'WNX'  , rvara=wnx)
  call shdf5_orec(ndims, idims, 'WNY'  , rvara=wny)
  call shdf5_orec(ndims, idims, 'WNZ'  , rvara=wnz)

  ndims    = 2
  idims(1) = nza
  idims(2) = nva

  call shdf5_orec(ndims, idims, 'ARV'  , rvara=arv)

  ndims    = 2
  idims(1) = nza
  idims(2) = nwa

  call shdf5_orec(ndims, idims, 'ARW'  , rvara=arw)
  call shdf5_orec(ndims, idims, 'VOLT' , dvara=volt)
  call shdf5_orec(ndims, idims, 'VOLTI', dvara=volti)

  ! Write ITAB_M SCALARS

  ndims    = 1
  idims(1) = nma

  call shdf5_orec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
  call shdf5_orec(ndims,idims,'itab_m%imp'      ,ivara=itab_m(:)%imp)
  call shdf5_orec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
  call shdf5_orec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
  call shdf5_orec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
  call shdf5_orec(ndims,idims,'itab_m%ngr'      ,ivara=itab_m(:)%ngr)

  ! Write ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = nma

  allocate (lscr(mloops,nma)) ; lscr = .false.
  do im = 1,nma
     lscr(1:mloops,im) = itab_m(im)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%loop',lvara=lscr)
  deallocate (lscr)

  ndims    = 2
  idims(1) = 7
  idims(2) = nma

  allocate (iscr(7,nma)) ; iscr = 0
  do im = 1,nma
     iscr(1:7,im) = itab_m(im)%iv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%iv',ivara=iscr)
  deallocate (iscr)

  allocate (iscr(7,nma)) ; iscr = 0
  do im = 1,nma
     iscr(1:7,im) = itab_m(im)%iw(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%iw',ivara=iscr)
  deallocate (iscr)

  ! Write ITAB_V SCALARS

  ndims    = 1
  idims(1) = nva
  idims(2) = 1

  call shdf5_orec(ndims,idims,'itab_v%ivp'    ,ivara=itab_v(:)%ivp)
  call shdf5_orec(ndims,idims,'itab_v%mrlv'   ,ivara=itab_v(:)%mrlv)

  ! Write ITAB_V ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = nva

  allocate (lscr(mloops,nva))
  do iv = 1,nva
     lscr(1:mloops,iv) = itab_v(iv)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%loop',lvara=lscr)
  deallocate(lscr)

  idims(1) = 2

  allocate (rscr(2,nva))
  do iv = 1,nva
     rscr(1:2,iv) = itab_v(iv)%farw(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%farw',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(2,nva))
  do iv = 1,nva
     rscr(1:2,iv) = itab_v(iv)%cosv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%cosv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(2,nva))
  do iv = 1,nva
     rscr(1:2,iv) = itab_v(iv)%sinv(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%sinv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(2,nva))
  do iv = 1,nva
     rscr(1:2,iv) = itab_v(iv)%dxps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%dxps',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(2,nva))
  do iv = 1,nva
     rscr(1:2,iv) = itab_v(iv)%dyps(1:2)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%dyps',rvara=rscr)
  deallocate(rscr)

  idims(1) = 4

  allocate (iscr(4,nva))
  do iv = 1,nva
     iscr(1:4,iv) = itab_v(iv)%iw(1:4)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%iw',ivara=iscr)
  deallocate(iscr)

  idims(1) = 6

  allocate (iscr(6,nva))
  do iv = 1,nva
     iscr(1:6,iv) = itab_v(iv)%im(1:6)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%im',ivara=iscr)
  deallocate(iscr)

  idims(1) = 4

  allocate (iscr(4,nva))
  do iv = 1,nva
     iscr(1:4,iv) = itab_v(iv)%iv(1:4)
  enddo
  call shdf5_orec(ndims,idims,'itab_v%iv',ivara=iscr)
  deallocate(iscr)

  ! Write ITAB_W SCALARS

  ndims    = 1
  idims(1) = nwa

  call shdf5_orec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly)
  call shdf5_orec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)
  call shdf5_orec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
  call shdf5_orec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
  call shdf5_orec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
  call shdf5_orec(ndims,idims,'itab_w%ngr'      ,ivara=itab_w(:)%ngr)

  call shdf5_orec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w)
  call shdf5_orec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w)

  call shdf5_orec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w)
  call shdf5_orec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w)
  call shdf5_orec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w)

  ! Write ITAB_W ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = nwa

  allocate (lscr(mloops,nwa))
  do iw = 1,nwa
     lscr(1:mloops,iw) = itab_w(iw)%loop(1:mloops)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%loop',lvara=lscr)
  deallocate(lscr)

  idims(1) = 3

  allocate (iscr(3,nwa))
  do iw = 1,nwa
     iscr(1:3,iw) = itab_w(iw)%iwnud(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iwnud',ivara=iscr)
  deallocate(iscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%fnud(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%fnud',rvara=rscr)
  deallocate(rscr)

  idims(1) = 7

  allocate (iscr(7,nwa))
  do iw = 1,nwa
     iscr(1:7,iw) = itab_w(iw)%im(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%im',ivara=iscr)
  deallocate(iscr)

  allocate (iscr(7,nwa))
  do iw = 1,nwa
     iscr(1:7,iw) = itab_w(iw)%iv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iv',ivara=iscr)
  deallocate(iscr)

  allocate (iscr(7,nwa))
  do iw = 1,nwa
     iscr(1:7,iw) = itab_w(iw)%iw(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iw',ivara=iscr)
  deallocate(iscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%dirv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%dirv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%farm(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%farm',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%farv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%farv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gxps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gxps1',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gyps1(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gyps1',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gxps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gxps2',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gyps2(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%gyps2',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))

  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%ecvec_vx(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vx',rvara=rscr)

  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%ecvec_vy(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vy',rvara=rscr)

  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%ecvec_vz(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%ecvec_vz',rvara=rscr)

  deallocate(rscr)

  ! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     ndims    = 1
     idims(1) = nwnud

     call shdf5_orec(ndims, idims, 'XEWNUD'  , rvara=xewnud)
     call shdf5_orec(ndims, idims, 'YEWNUD'  , rvara=yewnud)
     call shdf5_orec(ndims, idims, 'ZEWNUD'  , rvara=zewnud)

     call shdf5_orec(ndims,idims,'itab_wnud%npoly',ivara=itab_wnud(:)%npoly)

     allocate (iscr(6,nwnud))

     do iwnud = 1,nwnud
        iscr(1:6,iwnud) = itab_wnud(iwnud)%iwnud(1:6)
     enddo

     ndims    = 2
     idims(1) = 6
     idims(2) = nwnud

     call shdf5_orec(ndims,idims,'itab_wnud%iwnud' ,ivara=iscr)

     deallocate(iscr)

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
                      mza, mma, mua, mva, mwa, nsw_max, &
                      xem, yem, zem, xew, yew, zew, &
                      alloc_xyzem, alloc_xyzew
use leaf_coms,  only: isfcl

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank

use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, xewnud, yewnud, zewnud, &
                      alloc_nudge1

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iv, iw, iwnud, idz

integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, nxp0, nzp0, itopoflg0, isfcl0, ndz0

real    :: deltax0

logical :: exans

integer, allocatable :: ngrdll0(:)
real,    allocatable :: grdrad0(:)
real,    allocatable :: grdlat0(:,:)
real,    allocatable :: grdlon0(:,:)
real,    allocatable :: hdz0(:)
real,    allocatable :: dz0(:)

! Scratch arrays for copying input

logical, allocatable :: lscr(:,:)
integer, allocatable :: iscr(:,:)
real,    allocatable :: rscr(:,:)

  character(pathlen) :: flnm

! Check if gridfile exists

  flnm = trim(gridfile)//'.h5'

  inquire(file=flnm, exist=exans)

  if (.not. exans) then

! Gridfile does not exist.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!  Gridfile does not exist:'
     write(io6,*) '!!!  '//gridfile
     write(io6,*) '!!!  Stopping run'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   
     stop 'stop - no gridfile'
   
  endif

! Gridfile exists; open it

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'Opening grid file ', trim(flnm)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R')

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

  call shdf5_irec(ndims, idims, 'NUDFLAG' , ivars=nudflag)
  call shdf5_irec(ndims, idims, 'NUDNXP'  , ivars=nudnxp)

  allocate( ngrdll0 (ngrids0) )
  allocate( grdrad0 (ngrids0) )
  allocate( grdlat0 (ngrids0, maxngrdll) )
  allocate( grdlon0 (ngrids0, maxngrdll) )
  allocate( hdz0(ndz0))
  allocate( dz0 (ndz0))

  if (ndz0 > 1) then
     idims(1) = ndz0

     call shdf5_irec(ndims, idims, 'HDZ' , rvara=hdz0)
     call shdf5_irec(ndims, idims, 'DZ'  , rvara=dz0)
  endif

  ndims    = 1
  idims(1) = ngrids0

  call shdf5_irec(ndims, idims, 'NGRDLL' , ivara=ngrdll0)
  call shdf5_irec(ndims, idims, 'GRDRAD' , rvara=grdrad0)

  ndims    = 2
  idims(1) = ngrids0
  idims(2) = maxngrdll

  call shdf5_irec(ndims, idims, 'GRDLAT', rvara=grdlat0)
  call shdf5_irec(ndims, idims, 'GRDLON', rvara=grdlon0)

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
     if (abs(grdrad0 (ngr) - grdrad (ngr)) > 1.e1 ) ierr = 1

     do i = 1,ngrdll0(ngr)
        if (abs(grdlat0(ngr,i) - grdlat(ngr,i)) > 1.e-3) ierr = 1
        if (abs(grdlon0(ngr,i) - grdlon(ngr,i)) > 1.e-3) ierr = 1
     enddo
  enddo

  do idz = 1, min(ndz0,ndz)
     if (abs(hdz0(idz) - hdz(idz)) > 1.e1 ) ierr = 1
     if (abs(dz0 (idz) - dz (idz)) > 1.e1 ) ierr = 1
  enddo

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
     write(io6, '(a,20f10.1)') 'hdz0:     ',hdz0 (1:ndz)
     write(io6, '(a,20f10.1)') 'hdz:      ',hdz  (1:ndz)
     write(io6,*) ' '
     write(io6, '(a,20f10.1)') 'dz0:      ',dz0 (1:ndz)
     write(io6, '(a,20f10.1)') 'dz:       ',dz  (1:ndz)
     write(io6,*) ' '
     write(io6, '(a,20i12)')   'ngrdll0:  ',ngrdll0 (1:ngrids)
     write(io6, '(a,20i12)')   'ngrdll:   ',ngrdll  (1:ngrids)
     write(io6,*) ' '
     write(io6, '(a,20f12.1)') 'grdrad0:  ',grdrad0 (1:ngrids)
     write(io6, '(a,20f12.1)') 'grdrad:   ',grdrad  (1:ngrids)
     write(io6,*) ' '

     do ngr = 2, min(ngrids0,ngrids)
        write(io6, '(a,i5)') 'ngr: ',ngr
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

  call shdf5_irec(ndims, idims, 'NZA'    , ivars=nza)
  call shdf5_irec(ndims, idims, 'NMA'    , ivars=nma)
  call shdf5_irec(ndims, idims, 'NUA'    , ivars=nua)
  call shdf5_irec(ndims, idims, 'NVA'    , ivars=nva)
  call shdf5_irec(ndims, idims, 'NWA'    , ivars=nwa)
  call shdf5_irec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
  call shdf5_irec(ndims, idims, 'MRLS'   , ivars=mrls)
  call shdf5_irec(ndims, idims, 'NWNUD'  , ivars=nwnud)

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

  call shdf5_irec(ndims, idims, 'XEW', rvara=xew)
  call shdf5_irec(ndims, idims, 'YEW', rvara=yew)
  call shdf5_irec(ndims, idims, 'ZEW', rvara=zew)

! Read ITAB_M SCALARS

  ndims    = 1
  idims(1) = nma

  call shdf5_irec(ndims,idims,'itab_m%npoly',ivara=itab_m_pd(:)%npoly)
  call shdf5_irec(ndims,idims,'itab_m%imp'  ,ivara=itab_m_pd(:)%imp)

! Read ITAB_M ARRAYS

  ndims    = 2
  idims(2) = nma
  idims(1) = 7

  allocate (iscr(7,nma))
  call shdf5_irec(ndims,idims,'itab_m%iv',ivara=iscr)
  do im = 1,nma
     itab_m_pd(im)%iv(1:7) = iscr(1:7,im)
  enddo
  deallocate (iscr)

  allocate (iscr(7,nma))
  call shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr)
  do im = 1,nma
     itab_m_pd(im)%iw(1:7) = iscr(1:7,im)
  enddo
  deallocate (iscr)

! Read ITAB_V SCALARS

  ndims    = 1
  idims(1) = nva
      
  call shdf5_irec(ndims,idims,'itab_v%ivp',ivara=itab_v_pd(:)%ivp)

! Read ITAB_V ARRAYS

  ndims    = 2
  idims(2) = nva

  idims(1) = 4

  allocate (iscr( 4,nva))
  call shdf5_irec(ndims,idims,'itab_v%iw'  ,ivara=iscr)
  do iv = 1,nva
     itab_v_pd(iv)%iw(1:4) = iscr(1:4,iv)
  enddo
  deallocate (iscr)

  idims(1) = 4

  allocate (iscr(4,nva))
  call shdf5_irec(ndims,idims,'itab_v%iv',ivara=iscr)
  do iv = 1,nva
     itab_v_pd(iv)%iv(1:4) = iscr(1:4,iv)
  enddo
  deallocate (iscr)

  idims(1) = 6

  allocate (iscr(6,nva))
  call shdf5_irec(ndims,idims,'itab_v%im',ivara=iscr)
  do iv = 1,nva
     itab_v_pd(iv)%im(1:6) = iscr(1:6,iv)
  enddo
  deallocate (iscr)

! Read ITAB_W SCALARS

  ndims    = 1
  idims(1) = nwa
  idims(2) = 1

  call shdf5_irec(ndims,idims,'itab_w%npoly'   ,ivara=itab_w_pd(:)%npoly)
  call shdf5_irec(ndims,idims,'itab_w%iwp'     ,ivara=itab_w_pd(:)%iwp)

! Read ITAB_W ARRAYS

  ndims    = 2
  idims(2) = nwa

  idims(1) = 3

  allocate (iscr(3,nwa))
  call shdf5_irec(ndims,idims,'itab_w%iwnud',ivara=iscr)
  do iw = 1,nwa
     itab_w_pd(iw)%iwnud(1:3) = iscr(1:3,iw)
  enddo
  deallocate (iscr)

  idims(1) = 7

  allocate (iscr(7,nwa))
  call shdf5_irec(ndims,idims,'itab_w%im',ivara=iscr)
  do iw = 1,nwa
     itab_w_pd(iw)%im(1:7) = iscr(1:7,iw)
  enddo
  deallocate(iscr)

  idims(1) = 7

  allocate (iscr(7,nwa))
  call shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr)
  do iw = 1,nwa
     itab_w_pd(iw)%iw(1:7) = iscr(1:7,iw)
  enddo
  deallocate (iscr)

  idims(1) = 7

  allocate (iscr(7,nwa))
  call shdf5_irec(ndims,idims,'itab_w%iv',ivara=iscr)
  do iw = 1,nwa
     itab_w_pd(iw)%iv(1:7) = iscr(1:7,iw)
  enddo
  deallocate (iscr)

! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     call alloc_nudge1(nwnud)

     ndims    = 1
     idims(1) = nwnud

     call shdf5_irec(ndims, idims, 'XEWNUD'  , rvara=xewnud)
     call shdf5_irec(ndims, idims, 'YEWNUD'  , rvara=yewnud)
     call shdf5_irec(ndims, idims, 'ZEWNUD'  , rvara=zewnud)

  endif

! Close the GRIDFILE

  call shdf5_close()

  write(io6,*) 'end of gridfile_read '

end subroutine gridfile_read_pd

!===============================================================================

subroutine gridfile_read()

use max_dims,   only: maxngrdll, pathlen
use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                      itopoflg, deltax, ndz, hdz, dz, &
                      ngrdll, grdrad, grdlat, grdlon
use mem_ijtabs, only: mloops, mrls, &
                      itab_m, itab_v, itab_w
use mem_grid,   only: nza, &
                      mza, mma, mua, mva, mwa, nsw_max, &
                      zm, zt, dzm, dzt, dzim, dzit, &
                      zfacm, zfact, zfacim, zfacit, &
                      lpm, lpv, lpw, lsw, &
                      topm, topw, &
                      xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                      dnu, dniu, dnv, dniv, arw0, arm0, &
                      glatw, glonw, glatm, glonm, glatv, glonv, &
                      arv, arw, volt, volti
use leaf_coms,  only: isfcl

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank

use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, &
                      xewnud, yewnud, zewnud

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iv, iw, iwnud

integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, nxp0, nzp0, itopoflg0, isfcl0

real    :: deltax0, deltaz0, dzrat0, dzmax0, zbase0

logical :: exans

! Scratch arrays for copying input

logical, allocatable :: lscr(:,:)
integer, allocatable :: iscr(:,:)
real,    allocatable :: rscr(:,:)

! Pointers to the global index of the local point

integer, pointer :: lgma(:)
integer, pointer :: lgua(:)
integer, pointer :: lgva(:)
integer, pointer :: lgwa(:)

integer, pointer :: lgwnud(:)

  character(pathlen) :: flnm

  lgma => itab_m%imglobe
  lgva => itab_v%ivglobe
  lgua => null()
  lgwa => itab_w%iwglobe
  lgwnud => itab_wnud%iwnudglobe

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

  call shdf5_open(flnm,'R')

  ndims    = 1
  idims(1) = nza

  call shdf5_irec(ndims, idims, 'ZM'    , rvara=zm)
  call shdf5_irec(ndims, idims, 'ZT'    , rvara=zt)
  call shdf5_irec(ndims, idims, 'DZM'   , rvara=dzm)
  call shdf5_irec(ndims, idims, 'DZT'   , rvara=dzt)
  call shdf5_irec(ndims, idims, 'DZIM'  , rvara=dzim)
  call shdf5_irec(ndims, idims, 'DZIT'  , rvara=dzit)
  call shdf5_irec(ndims, idims, 'ZFACM' , rvara=zfacm)
  call shdf5_irec(ndims, idims, 'ZFACT' , rvara=zfact)
  call shdf5_irec(ndims, idims, 'ZFACIM', rvara=zfacim)
  call shdf5_irec(ndims, idims, 'ZFACIT', rvara=zfacit)

  ndims    = 1
  idims(1) = mma

  call shdf5_irec(ndims, idims, 'LPM'  , ivara=lpm, points=lgma)
  call shdf5_irec(ndims, idims, 'TOPM' , rvara=topm, points=lgma)
  call shdf5_irec(ndims, idims, 'ARM0' , rvara=arm0, points=lgma)
  call shdf5_irec(ndims, idims, 'GLATM', rvara=glatm, points=lgma)
  call shdf5_irec(ndims, idims, 'GLONM', rvara=glonm, points=lgma)
  call shdf5_irec(ndims, idims, 'XEM'  , rvara=xem, points=lgma)
  call shdf5_irec(ndims, idims, 'YEM'  , rvara=yem, points=lgma)
  call shdf5_irec(ndims, idims, 'ZEM'  , rvara=zem, points=lgma)
   
  ndims    = 1
  idims(1) = mva

  call shdf5_irec(ndims, idims, 'LPV' , ivara=lpv, points=lgva)
  call shdf5_irec(ndims, idims, 'XEV' , rvara=xev, points=lgva)
  call shdf5_irec(ndims, idims, 'YEV' , rvara=yev, points=lgva)
  call shdf5_irec(ndims, idims, 'ZEV' , rvara=zev, points=lgva)
  call shdf5_irec(ndims, idims, 'GLATV', rvara=glatv, points=lgva)
  call shdf5_irec(ndims, idims, 'GLONV', rvara=glonv, points=lgva)
  call shdf5_irec(ndims, idims, 'DNU' , rvara=dnu, points=lgva)
  call shdf5_irec(ndims, idims, 'DNV' , rvara=dnv, points=lgva)
  call shdf5_irec(ndims, idims, 'DNIU', rvara=dniu, points=lgva)
  call shdf5_irec(ndims, idims, 'DNIV', rvara=dniv, points=lgva)
  call shdf5_irec(ndims, idims, 'UNX' , rvara=unx, points=lgva)
  call shdf5_irec(ndims, idims, 'UNY' , rvara=uny, points=lgva)
  call shdf5_irec(ndims, idims, 'UNZ' , rvara=unz, points=lgva)
  call shdf5_irec(ndims, idims, 'VNX' , rvara=vnx, points=lgva)
  call shdf5_irec(ndims, idims, 'VNY' , rvara=vny, points=lgva)
  call shdf5_irec(ndims, idims, 'VNZ' , rvara=vnz, points=lgva)

  ndims    = 1
  idims(1) = mwa

  call shdf5_irec(ndims, idims, 'LPW'  , ivara=lpw, points=lgwa)
  call shdf5_irec(ndims, idims, 'LSW'  , ivara=lsw, points=lgwa)
  call shdf5_irec(ndims, idims, 'XEW'  , rvara=xew, points=lgwa)
  call shdf5_irec(ndims, idims, 'YEW'  , rvara=yew, points=lgwa)
  call shdf5_irec(ndims, idims, 'ZEW'  , rvara=zew, points=lgwa)
  call shdf5_irec(ndims, idims, 'WNX'  , rvara=wnx, points=lgwa)
  call shdf5_irec(ndims, idims, 'WNY'  , rvara=wny, points=lgwa)
  call shdf5_irec(ndims, idims, 'WNZ'  , rvara=wnz, points=lgwa)
  call shdf5_irec(ndims, idims, 'ARW0' , rvara=arw0, points=lgwa)
  call shdf5_irec(ndims, idims, 'TOPW' , rvara=topw, points=lgwa)
  call shdf5_irec(ndims, idims, 'GLATW', rvara=glatw, points=lgwa)
  call shdf5_irec(ndims, idims, 'GLONW', rvara=glonw, points=lgwa)

  ndims    = 2
  idims(1) = nza
  idims(2) = mva

  call shdf5_irec(ndims, idims, 'ARV'  , rvara=arv,   points=lgva)
   
  ndims    = 2
  idims(1) = nza
  idims(2) = mwa

  call shdf5_irec(ndims, idims, 'ARW'  , rvara=arw, points=lgwa)
  call shdf5_irec(ndims, idims, 'VOLT' , dvara=volt, points=lgwa)
  call shdf5_irec(ndims, idims, 'VOLTI', dvara=volti, points=lgwa)
   
! Read ITAB_M SCALARS

  ndims    = 1
  idims(1) = mma

  call shdf5_irec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly, points=lgma)
  call shdf5_irec(ndims,idims,'itab_m%imp'      ,ivara=itab_m(:)%imp, points=lgma)
  call shdf5_irec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm, points=lgma)
  call shdf5_irec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig, points=lgma)
  call shdf5_irec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow, points=lgma)
  call shdf5_irec(ndims,idims,'itab_m%ngr'      ,ivara=itab_m(:)%ngr, points=lgma)

! Read ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops
  idims(2) = mma

  allocate (lscr(mloops,mma))
  call shdf5_irec(ndims,idims,'itab_m%loop',lvara=lscr, points=lgma)
  do im = 1,mma
     itab_m(im)%loop(1:mloops) = lscr(1:mloops,im)
  enddo
  deallocate (lscr)

  idims(1) = 7

  allocate (iscr(7,mma))
  call shdf5_irec(ndims,idims,'itab_m%iv',ivara=iscr, points=lgma)
  do im = 1,mma
     itab_m(im)%iv(1:7) = iscr(1:7,im)
  enddo
  deallocate (iscr)

  allocate (iscr(7,mma))
  call shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr, points=lgma)
  do im = 1,mma
     itab_m(im)%iw(1:7) = iscr(1:7,im)
  enddo
  deallocate (iscr)

! Read ITAB_V SCALARS

  ndims    = 1
  idims(1) = mva

  call shdf5_irec(ndims,idims,'itab_v%ivp'      ,ivara=itab_v(:)%ivp, points=lgva)
  call shdf5_irec(ndims,idims,'itab_v%mrlv'     ,ivara=itab_v(:)%mrlv, points=lgva)

! Read ITAB_V ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = mva

  allocate (lscr(mloops,mva))
  call shdf5_irec(ndims,idims,'itab_v%loop',lvara=lscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%loop(1:mloops) = lscr(1:mloops,iv)
  enddo
  deallocate (lscr)

  idims(1) = 2

  allocate (rscr(2,mva))
  call shdf5_irec(ndims,idims,'itab_v%farw',rvara=rscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%farw(1:2) = rscr(1:2,iv)
  enddo
  deallocate (rscr)

  allocate (rscr(2,mva))
  call shdf5_irec(ndims,idims,'itab_v%cosv',rvara=rscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%cosv(1:2) = rscr(1:2,iv)
  enddo
  deallocate (rscr)

  allocate (rscr(2,mva))
  call shdf5_irec(ndims,idims,'itab_v%sinv',rvara=rscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%sinv(1:2) = rscr(1:2,iv)
  enddo
  deallocate (rscr)

  allocate (rscr(2,mva))
  call shdf5_irec(ndims,idims,'itab_v%dxps',rvara=rscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%dxps(1:2) = rscr(1:2,iv)
  enddo
  deallocate (rscr)

  allocate (rscr(2,mva))
  call shdf5_irec(ndims,idims,'itab_v%dyps',rvara=rscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%dyps(1:2) = rscr(1:2,iv)
  enddo
  deallocate (rscr)

  idims(1) = 4

  allocate (iscr(4,mva))
  call shdf5_irec(ndims,idims,'itab_v%iw',ivara=iscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%iw(1:4) = iscr(1:4,iv)
  enddo
  deallocate (iscr)

  idims(1) = 6

  allocate (iscr(6,mva))
  call shdf5_irec(ndims,idims,'itab_v%im',ivara=iscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%im(1:6) = iscr(1:6,iv)
  enddo
  deallocate (iscr)

  idims(1) = 4

  allocate (iscr(4,mva))
  call shdf5_irec(ndims,idims,'itab_v%iv',ivara=iscr, points=lgva)
  do iv = 1,mva
     itab_v(iv)%iv(1:4) = iscr(1:4,iv)
  enddo
  deallocate (iscr)

! Read ITAB_W SCALARS

  ndims    = 1
  idims(1) = mwa

  call shdf5_irec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%ngr'      ,ivara=itab_w(:)%ngr, points=lgwa)

  call shdf5_irec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w, points=lgwa)

  call shdf5_irec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w, points=lgwa)
  call shdf5_irec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w, points=lgwa)

! Read ITAB_W ARRAYS

  ndims    = 2
  idims(1) = mloops
  idims(2) = mwa

  allocate (lscr(mloops,mwa))
  call shdf5_irec(ndims,idims,'itab_w%loop',lvara=lscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%loop(1:mloops) = lscr(1:mloops,iw)
  enddo
  deallocate (lscr)

  idims(1) = 3

  allocate (iscr(3,mwa))
  call shdf5_irec(ndims,idims,'itab_w%iwnud',ivara=iscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%iwnud(1:3) = iscr(1:3,iw)
  enddo
  deallocate (iscr)

  allocate (rscr(3,mwa))
  call shdf5_irec(ndims,idims,'itab_w%fnud',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%fnud(1:3) = rscr(1:3,iw)
  enddo
  deallocate (rscr)

  idims(1) = 7

  allocate (iscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%im',ivara=iscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%im(1:7) = iscr(1:7,iw)
  enddo
  deallocate (iscr)

  allocate (iscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%iv',ivara=iscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%iv(1:7) = iscr(1:7,iw)
  enddo
  deallocate (iscr)

  allocate (iscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%iw(1:7) = iscr(1:7,iw)
  enddo
  deallocate (iscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%dirv',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%dirv(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%farm',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%farm(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%farv',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%farv(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%gxps1',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%gxps1(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%gyps1',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%gyps1(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%gxps2',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%gxps2(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)

  allocate (rscr(7,mwa))
  call shdf5_irec(ndims,idims,'itab_w%gyps2',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%gyps2(1:7) = rscr(1:7,iw)
  enddo
  deallocate (rscr)
   
  allocate (rscr(7,mwa))

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vx',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vx(1:7) = rscr(1:7,iw)
  enddo
      
  call shdf5_irec(ndims,idims,'itab_w%ecvec_vy',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vy(1:7) = rscr(1:7,iw)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vz',rvara=rscr, points=lgwa)
  do iw = 1,mwa
     itab_w(iw)%ecvec_vz(1:7) = rscr(1:7,iw)
  enddo

  deallocate (rscr)

! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     ndims    = 1
     idims(1) = mwnud

     call shdf5_irec(ndims, idims, 'XEWNUD'  , rvara=xewnud, points=lgwnud)
     call shdf5_irec(ndims, idims, 'YEWNUD'  , rvara=yewnud, points=lgwnud)
     call shdf5_irec(ndims, idims, 'ZEWNUD'  , rvara=zewnud, points=lgwnud)

     call shdf5_irec(ndims, idims, 'itab_wnud%npoly' ,ivara=itab_wnud(:)%npoly, points=lgwnud)

     allocate (iscr(6,mwnud))

     ndims    = 2
     idims(1) = 6
     idims(2) = mwnud

     call shdf5_irec(ndims,idims,'itab_wnud%iwnud',ivara=iscr, points=lgwnud)

     do iwnud = 1,mwnud
        itab_wnud(iwnud)%iwnud(1:6) = iscr(1:6,iwnud)
     enddo

     deallocate(iscr)

  endif

! Close the GRIDFILE

  call shdf5_close()

  write(io6,*) 'end of gridfile_read '

end subroutine gridfile_read

!===============================================================================

! This subroutine reads a few quantities from the OLD gridfile for a history
! start on which new grids (local mesh refinements) are being added.

subroutine gridfile_read_oldgrid()

use max_dims,   only: maxngrdll, pathlen
use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                      itopoflg, deltax, ndz, hdz, dz, &
                      ngrdll, grdrad, grdlat, grdlon
use mem_grid,   only: xew, yew, zew
use leaf_coms,  only: isfcl

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_addgrid,only: nzp_og, nxp_og, mdomain_og, ngrids_og, isfcl_og, &
                      itopoflg_og, ndz_og, deltax_og, hdz_og, dz_og, &
                      nza_og, nma_og, nua_og, nva_og, nwa_og, mrls_og, &
                      xew_og, yew_og, zew_og, wnx_og, wny_og, wnz_og, &
                      itab_wog, lpw_og

implicit none

integer :: iw, idz, iw_og
integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

logical :: exans

  integer, allocatable :: ngrdll_og(:)
  real,    allocatable :: grdrad_og(:)
  real,    allocatable :: grdlat_og(:,:)
  real,    allocatable :: grdlon_og(:,:)

! Scratch array for copying input

  integer, allocatable :: iscr(:,:)
  real,    allocatable :: rscr(:,:)

  character(pathlen) :: flnm

! Check if gridfile exists

  flnm = trim(gridfile)//'-OG'//'.h5'

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
  write(io6,*) 'Opening grid file ', flnm
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(flnm,'R')

! Read the grid information that exists in namelist

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp_og)
  call shdf5_irec(ndims, idims, 'NXP'     , ivars=nxp_og)
  call shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain_og)
  call shdf5_irec(ndims, idims, 'NGRIDS'  , ivars=ngrids_og)
  call shdf5_irec(ndims, idims, 'ISFCL'   , ivars=isfcl_og)
  call shdf5_irec(ndims, idims, 'ITOPOFLG', ivars=itopoflg_og)
  call shdf5_irec(ndims, idims, 'NDZ'     , ivars=ndz_og)
  call shdf5_irec(ndims, idims, 'DELTAX'  , rvars=deltax_og)

  allocate( hdz_og(ndz_og))
  allocate( dz_og (ndz_og))

  allocate( ngrdll_og (ngrids_og) )
  allocate( grdrad_og (ngrids_og) )
  allocate( grdlat_og (ngrids_og, maxngrdll) )
  allocate( grdlon_og (ngrids_og, maxngrdll) )

  if (ndz_og > 1) then
     idims(1) = ndz_og

     call shdf5_irec(ndims, idims, 'HDZ' , rvara=hdz_og)
     call shdf5_irec(ndims, idims, 'DZ'  , rvara=dz_og)
  endif

  ndims    = 1
  idims(1) = ngrids_og

  call shdf5_irec(ndims, idims, 'NGRDLL' , ivara=ngrdll_og)
  call shdf5_irec(ndims, idims, 'GRDRAD' , rvara=grdrad_og)

  ndims    = 2
  idims(1) = ngrids_og
  idims(2) = maxngrdll

  call shdf5_irec(ndims, idims, 'GRDLAT', rvara=grdlat_og)
  call shdf5_irec(ndims, idims, 'GRDLON', rvara=grdlon_og)

! Check equality between grid file information and namelist variables

  ierr = 0

  if (nzp_og      /= nzp     ) ierr = 1 
  if (nxp_og      /= nxp     ) ierr = 1 
  if (mdomain_og  /= mdomain ) ierr = 1 
  if (isfcl_og    /= isfcl   ) ierr = 1 
  if (itopoflg_og /= itopoflg) ierr = 1 
  if (ndz_og      /= ndz     ) ierr = 1 

  if (abs(deltax_og - deltax) > 1.e-3) ierr = 1 

  do ngr = 2, min(ngrids_og,ngrids)
     if (abs(ngrdll_og (ngr) - ngrdll (ngr)) > 1.e1 ) ierr = 1
     if (abs(grdrad_og (ngr) - grdrad (ngr)) > 1.e1 ) ierr = 1

     do i = 1,ngrdll_og(ngr)
        if (abs(grdlat_og(ngr,i) - grdlat(ngr,i)) > 1.e-3) ierr = 1
        if (abs(grdlon_og(ngr,i) - grdlon(ngr,i)) > 1.e-3) ierr = 1
     enddo
  enddo

  do idz = 1, min(ndz_og,ndz)
     if (abs(hdz_og(idz) - hdz(idz)) > 1.e1 ) ierr = 1
     if (abs(dz_og (idz) - dz (idz)) > 1.e1 ) ierr = 1
  enddo

  if (ierr == 1) then

     write(io6,*) 'OLDGRIDFILE mismatch with OLAMIN namelist: Stopping model run'
     write(io6,*) 'Values: oldgridfile, namelist'
     write(io6,*) '-----------------------------------------------'
     write(io6,*)              'nzp:      ',nzp_og     ,nzp
     write(io6,*)              'nxp:      ',nxp_og     ,nxp
     write(io6,*)              'mdomain:  ',mdomain_og ,mdomain
     write(io6,*)              'isfcl:    ',isfcl_og   ,isfcl
     write(io6,*)              'itopoflg: ',itopoflg_og,itopoflg
     write(io6,*)              'deltax:   ',deltax_og  ,deltax
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

  call shdf5_irec(ndims, idims, 'NZA'    , ivars=nza_og)
  call shdf5_irec(ndims, idims, 'NMA'    , ivars=nma_og)
  call shdf5_irec(ndims, idims, 'NUA'    , ivars=nua_og)
  call shdf5_irec(ndims, idims, 'NVA'    , ivars=nva_og)
  call shdf5_irec(ndims, idims, 'NWA'    , ivars=nwa_og)
  call shdf5_irec(ndims, idims, 'MRLS'   , ivars=mrls_og)

  allocate (lpw_og(nwa_og))

  allocate (xew_og(nwa_og))
  allocate (yew_og(nwa_og))
  allocate (zew_og(nwa_og))

  allocate (wnx_og(nwa_og))
  allocate (wny_og(nwa_og))
  allocate (wnz_og(nwa_og))

  allocate (itab_wog(nwa_og))

  idims(1) = nwa_og

  call shdf5_irec(ndims, idims, 'LPW', ivara=lpw_og)

  call shdf5_irec(ndims, idims, 'XEW', rvara=xew_og)
  call shdf5_irec(ndims, idims, 'YEW', rvara=yew_og)
  call shdf5_irec(ndims, idims, 'ZEW', rvara=zew_og)

  call shdf5_irec(ndims, idims, 'WNX', rvara=wnx_og)
  call shdf5_irec(ndims, idims, 'WNY', rvara=wny_og)
  call shdf5_irec(ndims, idims, 'WNZ', rvara=wnz_og)

  call shdf5_irec(ndims,idims,'itab_w%npoly',ivara=itab_wog(:)%npoly)

  ndims    = 2
  idims(1) = 7
  idims(2) = nwa_og

  allocate (iscr(7,nwa_og))

  call shdf5_irec(ndims,idims,'itab_w%iv',ivara=iscr)
  do iw_og = 1,nwa_og
     itab_wog(iw_og)%iv(1:7) = iscr(1:7,iw_og)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr)
  do iw_og = 1,nwa_og
     itab_wog(iw_og)%iw(1:7) = iscr(1:7,iw_og)
  enddo

  deallocate (iscr)

  allocate (rscr(7,nwa_og))

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vx',rvara=rscr)
  do iw_og = 1,nwa_og
     itab_wog(iw_og)%ecvec_vx(1:7) = rscr(1:7,iw_og)
  enddo
      
  call shdf5_irec(ndims,idims,'itab_w%ecvec_vy',rvara=rscr)
  do iw_og = 1,nwa_og
     itab_wog(iw_og)%ecvec_vy(1:7) = rscr(1:7,iw_og)
  enddo

  call shdf5_irec(ndims,idims,'itab_w%ecvec_vz',rvara=rscr)
  do iw_og = 1,nwa_og
     itab_wog(iw_og)%ecvec_vz(1:7) = rscr(1:7,iw_og)
  enddo

  deallocate (rscr)

  call shdf5_close()

  write(io6,*) 'end of gridfile_read_oldgrid '

end subroutine gridfile_read_oldgrid
