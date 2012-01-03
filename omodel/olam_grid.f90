!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida.
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine gridinit()

use misc_coms,   only: io6, runtype, mdomain, ngrids, initial, nxp, nzp, &
                       timmax8, alloc_misc, iparallel, meshtype, &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim

use leaf_coms,   only: nzg, nzs, isfcl, nwl
use sea_coms,    only: nws

use mem_ijtabs,  only: istp, mrls, fill_jtabs, itab_md, itab_ud, itab_wd, itab_w
use oplot_coms,  only: op
use mem_grid,    only: nza, nma, nua, nva, nwa, &
                       zm, zt, dzt, xem, yem, zem, xew, yew, zew, glatw, glonw, &
                       alloc_grid1, alloc_grid2
use mem_sflux,   only: nlandflux, nseaflux
use mem_nudge,   only: nudflag, nudnxp, nwnud, itab_wnud, alloc_nudge1, &
                       xewnud, yewnud, zewnud

implicit none

real, allocatable :: quarter_kite(:,:)

integer :: npoly
integer :: j, jmaxneg, jminpos
integer :: imd,imd1,imd2,iud,iwd,iwnud,iwnudn,iw

real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
real :: dist_wnud,dist,dist_min,xi,yi
real :: xin(6),yin(6)

! Read LAND and SEA files

if (runtype /= 'MAKESFC' .and. isfcl == 1) then

   write(io6,'(/,a)') 'gridinit calling landfile_read'
   call landfile_read()

   write(io6,'(/,a)') 'gridinit calling seafile_read'
   call seafile_read()

endif

! Generate OLAM grid structure for 'MAKESFC' or 'MAKEGRID' runtype

if (runtype == 'MAKESFC' .or. runtype == 'MAKEGRID') then

! Vertical grid coordinate setup

   write(io6,'(/,a)') 'gridinit calling gridset2'
   call gridset2()

   write(io6,'(a,f8.1)') ' Model top height = ',zm(nza-1)

! Print out vertical grid structure 

   call gridset_print()

! Horizontal grid setup

   if (mdomain == 0) then

! If doing 'MAKEGRID' run and using nudging on global domain, generate
! nudging grid here

      if (runtype == 'MAKEGRID' .and. nudflag > 0) then

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
      call cartesian_3d()    ! 3D cartesian channel domain; calls 2 allocs

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

! 'MAKESFC' run uses OLAM triangular atmos grid configuration
! to generate land/sea cells

   if (runtype == 'MAKESFC') then
      write(io6,'(/,a)') 'gridinit calling makesfc'
      call makesfc()
      return
   endif

! CHECK WHETHER THIS IS TRIANGLE OR HEXAGON MESH

   if (meshtype == 1) then
      call delaunay()

      write(io6,'(/,a)') 'gridinit after delaunay'
      write(io6,'(a,i8)')    ' nma = ',nma
      write(io6,'(a,i8)')    ' nua = ',nua
      write(io6,'(a,i8)')    ' nva = ',nva
      write(io6,'(a,i8)')    ' nwa = ',nwa
   else
      call voronoi()
      call pcvt()

      write(io6,'(/,a)') 'gridinit after voronoi and pcvt'
      write(io6,'(a,i8)')    ' nma = ',nma
      write(io6,'(a,i8)')    ' nua = ',nua
      write(io6,'(a,i8)')    ' nva = ',nva
      write(io6,'(a,i8)')    ' nwa = ',nwa
   endif

! Allocate remaining GRID FOOTPRINT arrays for full domain

   write(io6,'(/,a)') 'gridinit calling alloc_grid1 for full domain'

   call alloc_grid1(meshtype, nma, nua, nva, nwa)

! Fill remaining GRID FOOTPRINT geometry for full domain

   write(io6,'(/,a)') 'gridinit calling grid_geometry'

   if (meshtype == 1) then
      call grid_geometry_tri()
   else
      allocate (quarter_kite(2,nva))
      call grid_geometry_hex(quarter_kite)
   endif

! Initialize dtlm, dtsm, ndtrat, and nacoust,
! and compute the timestep schedule for all grid operations.

   write(io6,'(/,a)') 'gridinit calling modsched'

   call modsched()

   write(io6,'(/,a)') 'gridinit calling fill_jtabs'

   call fill_jtabs(nma,nua,nva,nwa)

! Allocate remaining unstructured grid geometry arrays

   write(io6,'(/,a)') 'gridinit calling alloc_grid2'

   call alloc_grid2(meshtype, nma, nua, nva, nwa)

! Set up control volumes

   write(io6,'(/,a)') 'gridinit calling ctrlvols'

   if (meshtype == 1) then
      call ctrlvols_tri()
   else
      call ctrlvols_hex(quarter_kite)
      deallocate (quarter_kite)
   endif

! If doing 'MAKEGRID' run and using nudging on global grid, compute
! nudging indices and coefficients of W points

   if (runtype == 'MAKEGRID' .and. mdomain == 0 .and. nudflag > 0) then

! Compute a mean distance between adjacent nudging points.  It is safe to assume
! that any ATM grid IW point will be closer than this distance to at least one
! nuding point

      dist_wnud = 7100.e3 / real(nudnxp)

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
! plane tangent at W point

!write(6,'(a,3i7,8e13.3)') 'og1 ',iw,j,iwnudn, &
!          xewnud(iwnudn),yewnud(iwnudn),zewnud(iwnudn), &
!          glatw(iw),glonw(iw)

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

! Write GRIDFILE

   write(io6,'(/,a)') 'gridinit calling gridfile_write'
   call gridfile_write()

   if (isfcl == 1) then
      write(io6,'(/,a)') 'gridinit before return to olam_run'
      write(io6,'(a,i8)')   ' nwl       = ',nwl
      write(io6,'(a,i8)')   ' nws       = ',nws
      write(io6,'(a,i8)')   ' nlandflux = ',nlandflux
      write(io6,'(a,i8)')   ' nseaflux  = ',nseaflux
   endif

endif

return
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

      ztarg = zend + .5 * dzend

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

   nza = kvec + 1
   mza = nza

   if (nza < 3) then
      write(io6,'(a)')      'OLAM requires that NZA >= 3.'
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

return
end subroutine gridset2

!===============================================================================

subroutine gridset_print()

use misc_coms, only: io6
use mem_grid,  only: nza, zm, zt, dzt

implicit none

integer :: k

! Print vertical grid structure 

write(io6,'(/,a)') '================================================='
write(io6,'(a)'  ) '         OLAM VERTICAL GRID STRUCTURE'
write(io6,'(a)'  ) '================================================='
write(io6,'(a)'  ) '   k    zm(m)     k     zt(m)    dzm(m)    dzrat'
write(io6,'(a,/)') '================================================='

do k = nza-1,1,-1
   if (k == nza-1 .or. k == 1) then
      write(io6,11) k, zm(k), dzt(k+1)/dzt(k)
   else
      write(io6,12) k, zm(k), dzt(k+1)/dzt(k)
   endif
   if (k > 1) write(io6,13) k,zt(k),dzt(k)
enddo
   
11 format (i4,f10.2,1x,3('========='),f6.3)
12 format (i4,f10.2,1x,3('---------'),f6.3)
13 format (15x,i4,2f10.2)

end subroutine gridset_print

!===============================================================================

subroutine topo_init(nqa,topq,glatq,glonq,xeq,yeq,zeq)

use misc_coms,   only: io6, deltax
use consts_coms, only: pi1, pio180

implicit none

integer, intent(in) :: nqa
real, intent(in) :: glatq(nqa),glonq(nqa),xeq(nqa),yeq(nqa),zeq(nqa)

real, intent(out) :: topq(nqa)

integer :: iq

real :: hfwid
real :: hgt
real :: hfwid2

real :: r, r0

! Fill the TOPQ array with a default value of 0 or modify it as desired.
! If itopoflg is set to 1, these values will be overridden in the call to
! topo_database, which inputs a standard OLAM topography dataset.

hfwid = 10000.

r0 = pi1 / 9.

! dudhia expts
! hfwid = 5. * .866 * deltax

! hgt = 405.
! hgt = 1012.
! end dudhia expts

hfwid2 = hfwid**2

do iq = 2,nqa
   topq(iq) = 0.

!   topq(iq) = 200. * mod(iq,4)

! SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   topq(iq) = max(0.,hgt * hfwid2 / (hfwid2 + xeq(iq)**2) - 1.)
!   write(io6,*) 'topq ',iq,xeq(iq),topq(iq)
! TOPQ = 0 AT LARGE DISTANCE FROM HILL
! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SPECIAL WM5 EXPT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find lat/lon of current M point

!   r = sqrt((glonq(iq) * pio180 + 0.5 * pi1) ** 2 &
!          + (glatq(iq) * pio180 - pi1 / 6.) ** 2)

!   topq(iq) = max(0., 2000. * (1. - r / r0))

!   print*, 'topoinit ',iq,r,r0,topq(iq)

! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine topo_init

!===============================================================================

subroutine gridfile_write()

  use max_dims,   only: maxngrdll
  use misc_coms,  only: io6, ngrids, gridfile, mdomain, meshtype, nzp, nxp, &
       iclobber, itopoflg, deltax, ndz, hdz, dz, ngrdll, grdrad, grdlat, grdlon
  use mem_ijtabs, only: mloops_m, mloops_u, mloops_v, mloops_w, mrls, &
       itab_m, itab_u, itab_v, itab_w
  use mem_grid,   only: nza, nma, nua, nva, nwa, nsw_max, &
       zm, zt, dzm, dzt, dzim, dzit, &
       zfacm, zfact, zfacim, zfacit, &
       lpm, lpu, lcu, lpv, lcv, lpw, lsw, &
       topm, topw, xem, yem, zem, xeu, yeu, zeu, &
       xev, yev, zev, xew, yew, zew, &
       unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
       dnu, dniu, dnv, dniv, arw0, arm0, &
       glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
       aru, arv, volui, volvi, arw, volwi, volt, volti
  use leaf_coms,  only: isfcl
  use mem_sflux,  only: nseaflux, nlandflux, seaflux, landflux, &
       nsfpats, nlfpats, nsfpatm, nlfpatm, &
       xemsfpat, yemsfpat, zemsfpat, &
       xemlfpat, yemlfpat, zemlfpat

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

  use mem_nudge,  only: nudflag, nudnxp, nwnud, itab_wnud, &
                        xewnud, yewnud, zewnud

 implicit none

  ! This routine writes the grid variables to the grid file.

  integer :: im, iu, iv, iw, iwnud

  integer :: ndims, idims(2)

  ! Scratch arrays for copying output

  logical, allocatable :: lscr(:,:)
  integer, allocatable :: iscr(:,:)

  real, allocatable :: rscr(:,:)

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'grid_write: opening file:', trim(gridfile)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(gridfile,'W',iclobber)

  ! Write the gridfile information that exists in namelist

  ndims    = 1
  idims(1) = 1
  idims(2) = 1

  call shdf5_orec(ndims, idims, 'NZP'     , ivars=nzp)
  call shdf5_orec(ndims, idims, 'NXP'     , ivars=nxp)
  call shdf5_orec(ndims, idims, 'MDOMAIN' , ivars=mdomain)
  call shdf5_orec(ndims, idims, 'MESHTYPE', ivars=meshtype)
  call shdf5_orec(ndims, idims, 'NGRIDS'  , ivars=ngrids)
  call shdf5_orec(ndims, idims, 'ISFCL'   , ivars=isfcl)
  call shdf5_orec(ndims, idims, 'ITOPOFLG', ivars=itopoflg)
  call shdf5_orec(ndims, idims, 'DELTAX'  , rvars=deltax)
  call shdf5_orec(ndims, idims, 'NDZ'     , ivars=ndz)

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

  call shdf5_orec(ndims, idims, 'NZA'    , ivars=nza)
  call shdf5_orec(ndims, idims, 'NMA'    , ivars=nma)
  call shdf5_orec(ndims, idims, 'NUA'    , ivars=nua)
  call shdf5_orec(ndims, idims, 'NVA'    , ivars=nva)
  call shdf5_orec(ndims, idims, 'NWA'    , ivars=nwa)
  call shdf5_orec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
  call shdf5_orec(ndims, idims, 'MRLS'   , ivars=mrls)

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

  if (meshtype == 1) then

     ndims    = 1
     idims(1) = nua

     call shdf5_orec(ndims, idims, 'LPU'  , ivara=lpu)
     call shdf5_orec(ndims, idims, 'LCU'  , ivara=lcu)
     call shdf5_orec(ndims, idims, 'XEU'  , rvara=xeu)
     call shdf5_orec(ndims, idims, 'YEU'  , rvara=yeu)
     call shdf5_orec(ndims, idims, 'ZEU'  , rvara=zeu)
     call shdf5_orec(ndims, idims, 'GLATU', rvara=glatu)
     call shdf5_orec(ndims, idims, 'GLONU', rvara=glonu)

  else

     ndims    = 1
     idims(1) = nva

     call shdf5_orec(ndims, idims, 'LPV'  , ivara=lpv)
     call shdf5_orec(ndims, idims, 'LCV'  , ivara=lcv)
     call shdf5_orec(ndims, idims, 'XEV'  , rvara=xev)
     call shdf5_orec(ndims, idims, 'YEV'  , rvara=yev)
     call shdf5_orec(ndims, idims, 'ZEV'  , rvara=zev)
     call shdf5_orec(ndims, idims, 'GLATV', rvara=glatv)
     call shdf5_orec(ndims, idims, 'GLONV', rvara=glonv)

  endif

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

  if (meshtype == 1) then

     ndims    = 2
     idims(1) = nza
     idims(2) = nua

     call shdf5_orec(ndims, idims, 'ARU'  , rvara=aru)
     call shdf5_orec(ndims, idims, 'VOLUI', rvara=volui)

  else

     ndims    = 2
     idims(1) = nza
     idims(2) = nva

     call shdf5_orec(ndims, idims, 'ARV'  , rvara=arv)
     call shdf5_orec(ndims, idims, 'ARU'  , rvara=aru)
     call shdf5_orec(ndims, idims, 'VOLVI', rvara=volvi)

  endif

  ndims    = 2
  idims(1) = nza
  idims(2) = nwa

  call shdf5_orec(ndims, idims, 'ARW'  , rvara=arw)
  call shdf5_orec(ndims, idims, 'VOLWI', rvara=volwi)
  call shdf5_orec(ndims, idims, 'VOLT' , dvara=volt)
  call shdf5_orec(ndims, idims, 'VOLTI', dvara=volti)

  ! Write ITAB_M SCALARS

  ndims    = 1
  idims(1) = nma

  call shdf5_orec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
  call shdf5_orec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm)
  call shdf5_orec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe)
  call shdf5_orec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
  call shdf5_orec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
  call shdf5_orec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
  call shdf5_orec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh)

  ! Write ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops_m
  idims(2) = nma

  allocate (lscr(mloops_m,nma)) ; lscr = .false.
  do im = 1,nma
     lscr(1:mloops_m,im) = itab_m(im)%loop(1:mloops_m)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%loop',lvara=lscr)
  deallocate (lscr)

  ndims    = 2
  idims(1) = 7
  idims(2) = nma

  allocate (iscr(7,nma)) ; iscr = 0
  if (meshtype == 1) then
     do im = 1,nma
        iscr(1:7,im) = itab_m(im)%iu(1:7)
     enddo
     call shdf5_orec(ndims,idims,'itab_m%iu',ivara=iscr)
  else
     do im = 1,nma
        iscr(1:7,im) = itab_m(im)%iv(1:7)
     enddo
     call shdf5_orec(ndims,idims,'itab_m%iv',ivara=iscr)
  endif
  deallocate (iscr)

  allocate (iscr(7,nma)) ; iscr = 0
  do im = 1,nma
     iscr(1:7,im) = itab_m(im)%iw(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%iw',ivara=iscr)
  deallocate (iscr)

  allocate (rscr(7,nma)) ; rscr = 0.0
  do im = 1,nma
     rscr(1:7,im) = itab_m(im)%fmw(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_m%fmw',rvara=rscr)
  deallocate (rscr)


  if (meshtype == 1) then

     ! Write ITAB_U SCALARS

     ndims    = 1
     idims(1) = nua

     call shdf5_orec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)
     call shdf5_orec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)
     call shdf5_orec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)

     call shdf5_orec(ndims,idims,'itab_u%gcf36'  ,rvara=itab_u(:)%gcf36)
     call shdf5_orec(ndims,idims,'itab_u%gcf45'  ,rvara=itab_u(:)%gcf45)

     call shdf5_orec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
     call shdf5_orec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
     call shdf5_orec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
     call shdf5_orec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
     call shdf5_orec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
     call shdf5_orec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
     call shdf5_orec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
     call shdf5_orec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)
     call shdf5_orec(ndims,idims,'itab_u%crossmm',rvara=itab_u(:)%crossmm)
     call shdf5_orec(ndims,idims,'itab_u%crossww',rvara=itab_u(:)%crossww)

     ! Write ITAB_U ARRAYS

     ndims    = 2
     idims(1) = mloops_u
     idims(2) = nua

     allocate (lscr(mloops_u,nua))
     do iu = 1,nua
        lscr(1:mloops_u,iu) = itab_u(iu)%loop(1:mloops_u)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%loop',lvara=lscr)
     deallocate (lscr)

     idims(1) = 2

     allocate (iscr(2,nua))
     do iu = 1,nua
        iscr(1:2,iu) = itab_u(iu)%im(1:2)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%im',ivara=iscr)
     deallocate (iscr)

     allocate (rscr(2,nua))
     do iu = 1,nua
        rscr(1:2,iu) = itab_u(iu)%vxw_u(1:2)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%vxw_u',rvara=rscr)
     deallocate (rscr)

     allocate (rscr(2,nua))
     do iu = 1,nua
        rscr(1:2,iu) = itab_u(iu)%vyw_u(1:2)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%vyw_u',rvara=rscr)
     deallocate (rscr)

     idims(1) = 4

     allocate (rscr(4,nua))
     do iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%diru(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%diru' ,rvara=rscr)
     deallocate (rscr)

     allocate (rscr(4,nua))
     do iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%tuu(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%tuu',rvara=rscr)
     deallocate (rscr)

     allocate (rscr(4,nua))
     do iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%guw(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%guw',rvara=rscr)
     deallocate (rscr)

     allocate (rscr(4,nua))
     do iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%vxu_u(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%vxu_u',rvara=rscr)
     deallocate (rscr)

     allocate (rscr(4,nua))
     do iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%vyu_u(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%vyu_u',rvara=rscr)
     deallocate (rscr)

     idims(1) = 6

     allocate (rscr(6,nua))
     do iu = 1,nua
        rscr(1:6,iu) = itab_u(iu)%fuw(1:6)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%fuw',rvara=rscr)
     deallocate (rscr)

     allocate (iscr(6,nua))
     do iu = 1,nua
        iscr(1:6,iu) = itab_u(iu)%iw(1:6)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%iw' ,ivara=iscr)
     deallocate (iscr)

     idims(1) = 12

     allocate (iscr(12,nua))
     do iu = 1,nua
        iscr(1:12,iu) = itab_u(iu)%iu(1:12)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%iu',ivara=iscr)
     deallocate (iscr)

     allocate (rscr(12,nua))
     do iu = 1,nua
        rscr(1:12,iu) = itab_u(iu)%fuu(1:12)
     enddo
     call shdf5_orec(ndims,idims,'itab_u%fuu',rvara=rscr)
     deallocate (rscr)

  endif

  if (meshtype == 2) then

     ! Write ITAB_V SCALARS

     ndims    = 1
     idims(1) = nva
     idims(2) = 1

     call shdf5_orec(ndims,idims,'itab_v%ivp'    ,ivara=itab_v(:)%ivp)
     call shdf5_orec(ndims,idims,'itab_v%mrlv'   ,ivara=itab_v(:)%mrlv)
     call shdf5_orec(ndims,idims,'itab_v%ivglobe',ivara=itab_v(:)%ivglobe)

     ! Write ITAB_V ARRAYS

     ndims    = 2
     idims(1) = mloops_v
     idims(2) = nva

     allocate (lscr(mloops_v,nva))
     do iv = 1,nva
        lscr(1:mloops_v,iv) = itab_v(iv)%loop(1:mloops_v)
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

     allocate (rscr(4,nva))
     do iv = 1,nva
        rscr(1:4,iv) = itab_v(iv)%fvw(1:4)
     enddo
     call shdf5_orec(ndims,idims,'itab_v%fvw',rvara=rscr)
     deallocate(rscr)

     idims(1) = 6

     allocate (iscr(6,nva))
     do iv = 1,nva
        iscr(1:6,iv) = itab_v(iv)%im(1:6)
     enddo
     call shdf5_orec(ndims,idims,'itab_v%im',ivara=iscr)
     deallocate(iscr)

     idims(1) = 12

     allocate (rscr(12,nva))
     do iv = 1,nva
        rscr(1:12,iv) = itab_v(iv)%fvv(1:12)
     enddo
     call shdf5_orec(ndims,idims,'itab_v%fvv',rvara=rscr)
     deallocate(rscr)

     idims(1) = 16

     allocate (iscr(16,nva))
     do iv = 1,nva
        iscr(1:16,iv) = itab_v(iv)%iv(1:16)
     enddo
     call shdf5_orec(ndims,idims,'itab_v%iv',ivara=iscr)
     deallocate(iscr)

     allocate (rscr(16,nva))
     do iv = 1,nva
        rscr(1:16,iv) = itab_v(iv)%fuv(1:16)
     enddo
     call shdf5_orec(ndims,idims,'itab_v%fuv',rvara=rscr)
     deallocate(rscr)

  endif

  ! Write ITAB_W SCALARS

  ndims    = 1
  idims(1) = nwa

  call shdf5_orec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly)
  call shdf5_orec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)
  call shdf5_orec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
  call shdf5_orec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
  call shdf5_orec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
  call shdf5_orec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
  call shdf5_orec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

  call shdf5_orec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)
  call shdf5_orec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)
  call shdf5_orec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

  call shdf5_orec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w)
  call shdf5_orec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w)

  call shdf5_orec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w)
  call shdf5_orec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w)
  call shdf5_orec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w)

  ! Write ITAB_W ARRAYS

  ndims = 2
  idims(1) = mloops_w
  idims(2) = nwa

  allocate (lscr(mloops_w,nwa))
  do iw = 1,nwa
     lscr(1:mloops_w,iw) = itab_w(iw)%loop(1:mloops_w)
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

  allocate (rscr (3,nwa))
  do iw = 1,nwa
     rscr (1:3,iw) = itab_w(iw)%diru (1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%diru',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vxu(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vxu',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vyu(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vyu',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vzu(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vzu',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vxu_w(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vxu_w',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(3,nwa))
  do iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vyu_w(1:3)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%vyu_w',rvara=rscr)
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

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%dirv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%dirv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%fwv(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%fwv',rvara=rscr)
  deallocate(rscr)

  allocate (rscr(7,nwa))
  do iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%fww(1:7)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%fww',rvara=rscr)
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

  idims(1) = 9

  allocate (iscr(9,nwa))
  do iw = 1,nwa
     iscr(1:9,iw) = itab_w(iw)%iu(1:9)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iu',ivara=iscr)
  deallocate(iscr)

  allocate (iscr(9,nwa))
  do iw = 1,nwa
     iscr(1:9,iw) = itab_w(iw)%iw(1:9)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%iw',ivara=iscr)
  deallocate(iscr)

  allocate (rscr(9,nwa))
  do iw = 1,nwa
     rscr(1:9,iw) = itab_w(iw)%fwu(1:9)
  enddo
  call shdf5_orec(ndims,idims,'itab_w%fwu',rvara=rscr)
  deallocate(rscr)

  ! Check whether LAND/SEA models are used

  if (isfcl == 1) then

     ! Write SEAFLUX VALUES

     ndims    = 1
     idims(1) = 1

     call shdf5_orec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
     call shdf5_orec(ndims, idims, 'NSFPATS',ivars=nsfpats)

     ndims    = 1
     idims(1) = nseaflux

     call shdf5_orec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%ifglobe)
     call shdf5_orec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
     call shdf5_orec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
     call shdf5_orec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iwls)
     call shdf5_orec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
     call shdf5_orec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
     call shdf5_orec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
     call shdf5_orec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
     call shdf5_orec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
     call shdf5_orec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
     call shdf5_orec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
     call shdf5_orec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sfc)

     if (nsfpats > 0) then

        ndims    = 1
        idims(1) = nsfpats

        call shdf5_orec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

        ndims    = 2
        idims(1) = 5
        idims(2) = nsfpats

        call shdf5_orec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
        call shdf5_orec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
        call shdf5_orec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

     endif

     ! Write LANDFLUX VALUES

     ndims    = 1
     idims(1) = 1

     call shdf5_orec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
     call shdf5_orec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

     ndims    = 1
     idims(1) = nlandflux

     call shdf5_orec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ifglobe)
     call shdf5_orec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
     call shdf5_orec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
     call shdf5_orec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwls)
     call shdf5_orec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
     call shdf5_orec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
     call shdf5_orec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
     call shdf5_orec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
     call shdf5_orec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
     call shdf5_orec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
     call shdf5_orec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
     call shdf5_orec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_sfc)

     if (nlfpats > 0) then

        ndims    = 1
        idims(1) = nlfpats

        call shdf5_orec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

        ndims    = 2
        idims(1) = 5
        idims(2) = nlfpats

        call shdf5_orec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
        call shdf5_orec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
        call shdf5_orec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

     endif

  endif

  ! Check whether NUDGING arrays are used

  if (mdomain == 0 .and. nudflag > 0) then

     ndims    = 1
     idims(1) = 1

     call shdf5_orec(ndims, idims, 'NUDNXP' , ivars=nudnxp)
     call shdf5_orec(ndims, idims, 'NWNUD'  , ivars=nwnud)

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

  return
end subroutine gridfile_write

!===============================================================================

! This subroutine reads all the scalars values, and only a few arrays needed
! by para_decomp (arrays ended with _pd)

subroutine gridfile_read_pd()

use max_dims,   only: maxngrdll
use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                      itopoflg, deltax, ndz, hdz, dz, &
                      ngrdll, grdrad, grdlat, grdlon, meshtype
use mem_ijtabs, only: mloops_m, mloops_u, mloops_v, mloops_w, mrls,  &
                      itab_m_pd, itab_u_pd, itab_v_pd, itab_w_pd, alloc_itabs_pd
use mem_grid,   only: nza, nma, nua, nva, nwa, &
                      mza, mma, mua, mva, mwa, nsw_max, &
                      xem, yem, zem, &
                      alloc_xyzem
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, mseaflux, mlandflux,  &
                      nsfpats, nlfpats, msfpats, mlfpats,  &
                      seaflux_pd, landflux_pd

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank

use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, &
                       xewnud, yewnud, zewnud, alloc_nudge1

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iu, iv, iw, iwnud, idz

integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, meshtype0, nxp0, nzp0, itopoflg0, isfcl0, ndz0

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

real, allocatable :: rscr(:,:)

! Check if grid file exists

inquire(file=gridfile, exist=exans)

if (exans) then

! Grid file exists.  Open, read, and close file.

   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(io6,*) 'Opening grid file ', trim(gridfile)
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   call shdf5_open(trim(gridfile),'R')

! Read the grid information that exists in namelist

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp0)
   call shdf5_irec(ndims, idims, 'NXP'     , ivars=nxp0)
   call shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain0)
   call shdf5_irec(ndims, idims, 'MESHTYPE', ivars=meshtype0)
   call shdf5_irec(ndims, idims, 'NGRIDS'  , ivars=ngrids0)
   call shdf5_irec(ndims, idims, 'ISFCL'   , ivars=isfcl0)
   call shdf5_irec(ndims, idims, 'ITOPOFLG', ivars=itopoflg0)
   call shdf5_irec(ndims, idims, 'DELTAX'  , rvars=deltax0)

   call shdf5_irec(ndims, idims, 'NDZ'     , ivars=ndz0)

   allocate( ngrdll0 (ngrids0) )
   allocate( grdrad0 (ngrids0) )
   allocate( grdlat0 (ngrids0, maxngrdll) )
   allocate( grdlon0 (ngrids0, maxngrdll) )
   allocate( hdz0(ndz0))
   allocate( dz0 (ndz0))

   if (ndz > 1) then
      idims(1) = ndz

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
   if (meshtype0 /= meshtype) ierr = 1 
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
      write(io6,*)              'meshtype: ',meshtype0,meshtype
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

! Copy grid dimensions

   mza = nza
   mma = nma
   mva = nva
   mua = nua
   mwa = nwa

! Allocate and read grid structure variables

   call alloc_itabs_pd(meshtype,nma,nua,nva,nwa)
   call alloc_xyzem(nma)

   ndims    = 1
   idims(1) = nma

   call shdf5_irec(ndims, idims, 'XEM', rvara=xem)
   call shdf5_irec(ndims, idims, 'YEM', rvara=yem)
   call shdf5_irec(ndims, idims, 'ZEM', rvara=zem)

! Read ITAB_M SCALARS

   ndims    = 1
   idims(1) = nma

   call shdf5_irec(ndims,idims,'itab_m%npoly',ivara=itab_m_pd(:)%npoly)
   call shdf5_irec(ndims,idims,'itab_m%itopm',ivara=itab_m_pd(:)%itopm)

! Read ITAB_M ARRAYS

   ndims    = 2
   idims(2) = nma
   idims(1) = 7

   allocate (iscr(7,nma))
   if (meshtype == 1) then
      call shdf5_irec(ndims,idims,'itab_m%iu',ivara=iscr)
      do im = 1,nma
         itab_m_pd(im)%iu(1:7) = iscr(1:7,im)
      enddo
   else
      call shdf5_irec(ndims,idims,'itab_m%iv',ivara=iscr)
      do im = 1,nma
         itab_m_pd(im)%iv(1:7) = iscr(1:7,im)
      enddo
   endif
   deallocate (iscr)

   allocate (iscr(7,nma))
   call shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr)
   do im = 1,nma
      itab_m_pd(im)%iw(1:7) = iscr(1:7,im)
   enddo
   deallocate (iscr)

   if (meshtype == 1) then

! Read ITAB_U SCALARS

      ndims    = 1
      idims(1) = nua
      idims(2) = 1
      
      call shdf5_irec(ndims,idims,'itab_u%iup',ivara=itab_u_pd(:)%iup)

! Read ITAB_U ARRAYS

      ndims    = 2
      idims(2) = nua

      idims(1) = 6

      allocate (iscr( 6,nua))

      call shdf5_irec(ndims,idims,'itab_u%iw',ivara=iscr)

      do iu = 1,nua
         itab_u_pd(iu)%iw(1: 6) = iscr(1:6,iu)
      enddo

      deallocate (iscr)

      idims(1) = 12

      allocate (iscr(12,nua))
      call shdf5_irec(ndims,idims,'itab_u%iu',ivara=iscr)
      do iu = 1,nua
         itab_u_pd(iu)%iu(1:12) = iscr(1:12,iu)
      enddo
      deallocate (iscr)

      idims(1) = 2

      allocate (iscr(2,nua))
      call shdf5_irec(ndims,idims,'itab_u%im',ivara=iscr)
      do iu = 1,nua
         itab_u_pd(iu)%im(1:2) = iscr(1:2,iu)
      enddo
      deallocate (iscr)

   elseif (meshtype == 2) then

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

      idims(1) = 16

      allocate (iscr(16,nva))
      call shdf5_irec(ndims,idims,'itab_v%iv',ivara=iscr)
      do iv = 1,nva
         itab_v_pd(iv)%iv(1:16) = iscr(1:16,iv)
      enddo
      deallocate (iscr)

      idims(1) = 6

      allocate (iscr(6,nva))
      call shdf5_irec(ndims,idims,'itab_v%im',ivara=iscr)
      do iv = 1,nva
         itab_v_pd(iv)%im(1:6) = iscr(1:6,iv)
      enddo
      deallocate (iscr)

   endif

! Read ITAB_W SCALARS

   ndims    = 1
   idims(1) = nwa
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_w%npoly',ivara=itab_w_pd(:)%npoly)
   call shdf5_irec(ndims,idims,'itab_w%iwp'  ,ivara=itab_w_pd(:)%iwp)

! Read ITAB_W ARRAYS

   ndims    = 2
   idims(2) = nwa

   idims(1) = 7

   allocate (iscr(7,nwa))
   call shdf5_irec(ndims,idims,'itab_w%im',ivara=iscr)
   do iw = 1,nwa
      itab_w_pd(iw)%im(1:7) = iscr(1:7,iw)
   enddo
   deallocate(iscr)

   idims(1) = 9

   allocate (iscr(9,nwa))
   call shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr)
   do iw = 1,nwa
      itab_w_pd(iw)%iw(1:9) = iscr(1:9,iw)
   enddo
   deallocate (iscr)

   if (meshtype == 1) then

      idims(1) = 9

      allocate (iscr(9,nwa))
      call shdf5_irec(ndims,idims,'itab_w%iu',ivara=iscr)
      do iw = 1,nwa
         itab_w_pd(iw)%iu(1:9) = iscr(1:9,iw)
      enddo
      deallocate (iscr)

   else

      idims(1) = 7

      allocate (iscr(7,nwa))
      call shdf5_irec(ndims,idims,'itab_w%iv',ivara=iscr)
      do iw = 1,nwa
         itab_w_pd(iw)%iv(1:7) = iscr(1:7,iw)
      enddo
      deallocate (iscr)

   endif

! Check whether LAND/SEA models are used

   if (isfcl == 1) then

! Read SEAFLUX VALUES

      ndims    = 1
      idims(1) = 1

      call shdf5_irec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
      call shdf5_irec(ndims, idims, 'NSFPATS' ,ivars=nsfpats)

      mseaflux = nseaflux
      msfpats = nsfpats

      allocate (seaflux_pd(nseaflux))

      ndims    = 1
      idims(1) = nseaflux

      call shdf5_irec(ndims,idims,'seaflux%iw'      ,ivara=seaflux_pd(:)%iw)
      call shdf5_irec(ndims,idims,'seaflux%iws'     ,ivara=seaflux_pd(:)%iwls)

! Read LANDFLUX VALUES

      ndims    = 1
      idims(1) = 1

      call shdf5_irec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
      call shdf5_irec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

      mlandflux = nlandflux
      mlfpats = nlfpats

      allocate (landflux_pd(nlandflux))

      ndims    = 1
      idims(1) = nlandflux

      call shdf5_irec(ndims,idims,'landflux%iw'      ,ivara=landflux_pd(:)%iw)
      call shdf5_irec(ndims,idims,'landflux%iwl'     ,ivara=landflux_pd(:)%iwls)

   endif

! Close the GRIDFILE

   call shdf5_close()

else

! Grid file does not exist.

   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(io6,*) '!!!  Gridfile does not exist:'
   write(io6,*) '!!!  '//trim(gridfile)
   write(io6,*) '!!!  Stopping run'
   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   
   stop 'stop - no gridfile'
   
endif

write(io6,*) 'end of gridfile_read '

return
end subroutine gridfile_read_pd

!===============================================================================

subroutine gridfile_read()

use max_dims,   only: maxngrdll
use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp, &
                      itopoflg, deltax, ndz, hdz, dz, &
                      ngrdll, grdrad, grdlat, grdlon, meshtype
use mem_ijtabs, only: mloops_m, mloops_u, mloops_v, mloops_w, mrls, &
                      itab_m, itab_u, itab_v, itab_w
use mem_grid,   only: nza, &
                      mza, mma, mua, mva, mwa, nsw_max, &
                      zm, zt, dzm, dzt, dzim, dzit, &
                      zfacm, zfact, zfacim, zfacit, &
                      lpm, lpu, lcu, lpv, lcv, lpw, lsw, &
                      topm, topw, xeu, yeu, zeu, &
                      xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                      dnu, dniu, dnv, dniv, arw0, arm0, &
                      glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
                      aru, arv, volui, volvi, arw, volwi, volt, volti
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, mseaflux, mlandflux, &
                      nsfpats, nlfpats, msfpats, mlfpats, nsfpatm, nlfpatm, &
                      seaflux, landflux, &
                      xemsfpat, yemsfpat, zemsfpat, &
                      xemlfpat, yemlfpat, zemlfpat

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank

use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, &
                       xewnud, yewnud, zewnud, alloc_nudge1

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iu, iv, iw, iwnud

integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, meshtype0, nxp0, nzp0, itopoflg0, isfcl0

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

lgma => itab_m%imglobe
lgua => itab_u%iuglobe
lgva => itab_v%ivglobe
lgwa => itab_w%iwglobe

! Check if grid file exists

inquire(file=gridfile, exist=exans)

if (exans) then

! Grid file exists.  Open, read, and close file.

   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(io6,*) 'Opening grid file ', trim(gridfile)
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   call shdf5_open(trim(gridfile),'R')

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
   
   if (meshtype == 1) then

      ndims    = 1
      idims(1) = mua

      call shdf5_irec(ndims, idims, 'LPU'  , ivara=lpu, points=lgua)
      call shdf5_irec(ndims, idims, 'LCU'  , ivara=lcu, points=lgua)
      call shdf5_irec(ndims, idims, 'XEU'  , rvara=xeu, points=lgua)
      call shdf5_irec(ndims, idims, 'YEU'  , rvara=yeu, points=lgua)
      call shdf5_irec(ndims, idims, 'ZEU'  , rvara=zeu, points=lgua)
      call shdf5_irec(ndims, idims, 'GLATU', rvara=glatu, points=lgua)
      call shdf5_irec(ndims, idims, 'GLONU', rvara=glonu, points=lgua)
      call shdf5_irec(ndims, idims, 'DNU' , rvara=dnu, points=lgua)
      call shdf5_irec(ndims, idims, 'DNV' , rvara=dnv, points=lgua)
      call shdf5_irec(ndims, idims, 'DNIU', rvara=dniu, points=lgua)
      call shdf5_irec(ndims, idims, 'DNIV', rvara=dniv, points=lgua)
      call shdf5_irec(ndims, idims, 'UNX' , rvara=unx, points=lgua)
      call shdf5_irec(ndims, idims, 'UNY' , rvara=uny, points=lgua)
      call shdf5_irec(ndims, idims, 'UNZ' , rvara=unz, points=lgua)
      call shdf5_irec(ndims, idims, 'VNX' , rvara=vnx, points=lgua)
      call shdf5_irec(ndims, idims, 'VNY' , rvara=vny, points=lgua)
      call shdf5_irec(ndims, idims, 'VNZ' , rvara=vnz, points=lgua)

   else

      ndims    = 1
      idims(1) = mva

      call shdf5_irec(ndims, idims, 'LPV' , ivara=lpv, points=lgva)
      call shdf5_irec(ndims, idims, 'LCV' , ivara=lcv, points=lgva)
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

   endif


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

   if (meshtype == 1) then

      ndims    = 2
      idims(1) = nza
      idims(2) = mua

      call shdf5_irec(ndims, idims, 'VOLUI', rvara=volui, points=lgua)
      call shdf5_irec(ndims, idims, 'ARU'  , rvara=aru,   points=lgua)
   else

      ndims    = 2
      idims(1) = nza
      idims(2) = mva

      call shdf5_irec(ndims, idims, 'ARV'  , rvara=arv,   points=lgva)
      call shdf5_irec(ndims, idims, 'ARU'  , rvara=aru,   points=lgva)
      call shdf5_irec(ndims, idims, 'VOLVI', rvara=volvi, points=lgva)
   
   endif

   ndims    = 2
   idims(1) = nza
   idims(2) = mwa

   call shdf5_irec(ndims, idims, 'ARW'  , rvara=arw, points=lgwa)
   call shdf5_irec(ndims, idims, 'VOLWI', rvara=volwi, points=lgwa)
   call shdf5_irec(ndims, idims, 'VOLT' , dvara=volt, points=lgwa)
   call shdf5_irec(ndims, idims, 'VOLTI', dvara=volti, points=lgwa)
   
! Read ITAB_M SCALARS

   ndims    = 1
   idims(1) = mma

   call shdf5_irec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow, points=lgma)
   call shdf5_irec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh, points=lgma)

! Read ITAB_M ARRAYS

   ndims = 2
   idims(1) = mloops_m
   idims(2) = mma

   allocate (lscr(mloops_m,mma))
   call shdf5_irec(ndims,idims,'itab_m%loop',lvara=lscr, points=lgma)
   do im = 1,mma
      itab_m(im)%loop(1:mloops_m) = lscr(1:mloops_m,im)
   enddo
   deallocate (lscr)

   idims(1) = 7

   allocate (iscr(7,mma))
   if (meshtype == 1) then
      call shdf5_irec(ndims,idims,'itab_m%iu',ivara=iscr, points=lgma)
      do im = 1,mma
         itab_m(im)%iu(1:7) = iscr(1:7,im)
      enddo
   else
      call shdf5_irec(ndims,idims,'itab_m%iv',ivara=iscr, points=lgma)
      do im = 1,mma
         itab_m(im)%iv(1:7) = iscr(1:7,im)
      enddo
   endif
   deallocate (iscr)

   allocate (iscr(7,mma))
   call shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr, points=lgma)
   do im = 1,mma
      itab_m(im)%iw(1:7) = iscr(1:7,im)
   enddo
   deallocate (iscr)

   allocate (rscr(7,mma))
   call shdf5_irec(ndims,idims,'itab_m%fmw',rvara=rscr, points=lgma)
   do im = 1,mma
      itab_m(im)%fmw(1:7) = rscr(1:7,im)
   enddo
   deallocate (rscr)

   if (meshtype == 1) then

! Read ITAB_U SCALARS

      ndims = 1
      idims(1) = mua
      idims(2) = 1

      call shdf5_irec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe, points=lgua)

      call shdf5_irec(ndims,idims,'itab_u%gcf36'  ,rvara=itab_u(:)%gcf36, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%gcf45'  ,rvara=itab_u(:)%gcf45, points=lgua)

      call shdf5_irec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%crossmm',rvara=itab_u(:)%crossmm, points=lgua)
      call shdf5_irec(ndims,idims,'itab_u%crossww',rvara=itab_u(:)%crossww, points=lgua)

! Read ITAB_U ARRAYS

      ndims = 2
      idims(1) = mloops_u
      idims(2) = mua

      allocate (lscr(mloops_u,mua))
      call shdf5_irec(ndims,idims,'itab_u%loop',lvara=lscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%loop(1:mloops_u) = lscr(1:mloops_u,iu)
      enddo
      deallocate (lscr)

      idims(1) = 2

      allocate (iscr(2,mua))
      call shdf5_irec(ndims,idims,'itab_u%im',ivara=iscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%im(1:2) = iscr(1:2,iu)
      enddo
      deallocate (iscr)

      allocate (rscr(2,mua))
      call shdf5_irec(ndims,idims,'itab_u%vxw_u',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%vxw_u(1:2) = rscr(1:2,iu)
      enddo
      deallocate (rscr)

      allocate (rscr(2,mua))
      call shdf5_irec(ndims,idims,'itab_u%vyw_u',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%vyw_u(1:2) = rscr(1:2,iu)
      enddo
      deallocate (rscr)

      idims(1) = 4

      allocate (rscr(4,mua))
      call shdf5_irec(ndims,idims,'itab_u%diru',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%diru(1:4) = rscr(1:4,iu)
      enddo
      deallocate (rscr)

      allocate (rscr(4,mua))
      call shdf5_irec(ndims,idims,'itab_u%tuu',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%tuu(1:4) = rscr(1:4,iu)
      enddo
      deallocate (rscr)

      allocate (rscr(4,mua))
      call shdf5_irec(ndims,idims,'itab_u%guw',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%guw(1:4) = rscr(1:4,iu)
      enddo
      deallocate (rscr)

      allocate (rscr(4,mua))
      call shdf5_irec(ndims,idims,'itab_u%vxu_u',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%vxu_u(1:4) = rscr(1:4,iu)
      enddo
      deallocate (rscr)

      allocate (rscr(4,mua))
      call shdf5_irec(ndims,idims,'itab_u%vyu_u',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%vyu_u(1:4) = rscr(1:4,iu)
      enddo
      deallocate (rscr)

      idims(1) = 6

      allocate (iscr(6,mua))
      call shdf5_irec(ndims,idims,'itab_u%iw',ivara=iscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%iw(1:6) = iscr(1:6,iu)
      enddo
      deallocate (iscr)

      allocate (rscr(6,mua))
      call shdf5_irec(ndims,idims,'itab_u%fuw',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%fuw(1:6) = rscr(1:6,iu)
      enddo
      deallocate (rscr)

      idims(1) = 12

      allocate (iscr(12,mua))
      call shdf5_irec(ndims,idims,'itab_u%iu',ivara=iscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%iu(1:12) = iscr(1:12,iu)
      enddo
      deallocate (iscr)

      allocate (rscr(12,mua))
      call shdf5_irec(ndims,idims,'itab_u%fuu',rvara=rscr, points=lgua)
      do iu = 1,mua
         itab_u(iu)%fuu(1:12) = rscr(1:12,iu)
      enddo
      deallocate (rscr)

   endif

   if (meshtype == 2) then

! Read ITAB_V SCALARS

      ndims    = 1
      idims(1) = mva

      call shdf5_irec(ndims,idims,'itab_v%ivp'      ,ivara=itab_v(:)%ivp, points=lgva)
      call shdf5_irec(ndims,idims,'itab_v%mrlv'     ,ivara=itab_v(:)%mrlv, points=lgva)
      call shdf5_irec(ndims,idims,'itab_v%ivglobe'  ,ivara=itab_v(:)%ivglobe, points=lgva)

! Read ITAB_V ARRAYS

      ndims    = 2
      idims(1) = mloops_v
      idims(2) = mva

      allocate (lscr(mloops_v,mva))
      call shdf5_irec(ndims,idims,'itab_v%loop',lvara=lscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%loop(1:mloops_v) = lscr(1:mloops_v,iv)
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

      allocate (rscr(4,mva))
      call shdf5_irec(ndims,idims,'itab_v%fvw',rvara=rscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%fvw (1:4) = rscr(1:4,iv)
      enddo
      deallocate (rscr)

      idims(1) = 6

      allocate (iscr(6,mva))
      call shdf5_irec(ndims,idims,'itab_v%im',ivara=iscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%im(1:6) = iscr(1:6,iv)
      enddo
      deallocate (iscr)

      idims(1) = 12

      allocate (rscr(12,mva))
      call shdf5_irec(ndims,idims,'itab_v%fvv',rvara=rscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%fvv(1:12) = rscr(1:12,iv)
      enddo
      deallocate (rscr)

      idims(1) = 16

      allocate (iscr(16,mva))
      call shdf5_irec(ndims,idims,'itab_v%iv',ivara=iscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%iv(1:16) = iscr(1:16,iv)
      enddo
      deallocate (iscr)

      allocate (rscr(16,mva))
      call shdf5_irec(ndims,idims,'itab_v%fuv',rvara=rscr, points=lgva)
      do iv = 1,mva
         itab_v(iv)%fuv(1:16) = rscr(1:16,iv)
      enddo
      deallocate (rscr)

   endif

! Read ITAB_W SCALARS

   ndims    = 1
   idims(1) = mwa

   call shdf5_irec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh, points=lgwa)

   call shdf5_irec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw, points=lgwa)

   call shdf5_irec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w, points=lgwa)

   call shdf5_irec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w, points=lgwa)
   call shdf5_irec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w, points=lgwa)

! Read ITAB_W ARRAYS

   ndims    = 2
   idims(1) = mloops_w
   idims(2) = mwa

   allocate (lscr(mloops_w,mwa))
   call shdf5_irec(ndims,idims,'itab_w%loop',lvara=lscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%loop(1:mloops_w) = lscr(1:mloops_w,iw)
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

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%diru',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%diru(1:3) = rscr(1:3,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%vxu',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%vxu(1:3) = rscr(1:3,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%vyu',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%vyu(1:3) = rscr(1:3,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%vzu',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%vzu(1:3) = rscr(1:3,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%vxu_w',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%vxu_w(1:3) = rscr(1:3,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(3,mwa))
   call shdf5_irec(ndims,idims,'itab_w%vyu_w',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%vyu_w(1:3) = rscr(1:3,iw)
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

   allocate (rscr(7,mwa))
   call shdf5_irec(ndims,idims,'itab_w%dirv',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%dirv(1:7) = rscr(1:7,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(7,mwa))
   call shdf5_irec(ndims,idims,'itab_w%fwv',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%fwv(1:7) = rscr(1:7,iw)
   enddo
   deallocate (rscr)

   allocate (rscr(7,mwa))
   call shdf5_irec(ndims,idims,'itab_w%fww',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%fww(1:7) = rscr(1:7,iw)
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

   idims(1) = 9

   allocate (iscr(9,mwa))
   call shdf5_irec(ndims,idims,'itab_w%iu',ivara=iscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%iu(1:9) = iscr(1:9,iw)
   enddo
   deallocate (iscr)

   allocate (iscr(9,mwa))
   call shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%iw(1:9) = iscr(1:9,iw)
   enddo
   deallocate (iscr)

   allocate (rscr(9,mwa))
   call shdf5_irec(ndims,idims,'itab_w%fwu',rvara=rscr, points=lgwa)
   do iw = 1,mwa
      itab_w(iw)%fwu(1:9) = rscr(1:9,iw)
   enddo
   deallocate (rscr)

! Check whether LAND/SEA models are used

   if (isfcl == 1) then

! Read SEAFLUX VALUES

      ndims    = 1
      idims(1) = 1

      call shdf5_irec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
      call shdf5_irec(ndims, idims, 'NSFPATS' ,ivars=nsfpats)

      mseaflux = nseaflux
      msfpats = nsfpats

      allocate (seaflux(nseaflux))

      allocate (nsfpatm(nsfpats))

      allocate (xemsfpat(5,nsfpats))
      allocate (yemsfpat(5,nsfpats))
      allocate (zemsfpat(5,nsfpats))

      ndims    = 1
      idims(1) = nseaflux

      call shdf5_irec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%ifglobe)
      call shdf5_irec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
      call shdf5_irec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
      call shdf5_irec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iwls)
      call shdf5_irec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
      call shdf5_irec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
      call shdf5_irec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
      call shdf5_irec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
      call shdf5_irec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
      call shdf5_irec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
      call shdf5_irec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sfc)

      if (nsfpats > 0) then

         ndims    = 1
        idims(1) = nsfpats

         call shdf5_irec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

         ndims    = 2
         idims(1) = 5
         idims(2) = nsfpats

         call shdf5_irec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
         call shdf5_irec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
         call shdf5_irec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

      endif

! Read LANDFLUX VALUES

      ndims    = 1
      idims(1) = 1

      call shdf5_irec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
      call shdf5_irec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

      mlandflux = nlandflux
      mlfpats = nlfpats

      allocate (landflux(nlandflux))

      allocate (nlfpatm(nlfpats))

      allocate (xemlfpat(5,nlfpats))
      allocate (yemlfpat(5,nlfpats))
      allocate (zemlfpat(5,nlfpats))

      ndims    = 1
      idims(1) = nlandflux

      call shdf5_irec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ifglobe)
      call shdf5_irec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
      call shdf5_irec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
      call shdf5_irec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwls)
      call shdf5_irec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
      call shdf5_irec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
      call shdf5_irec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
      call shdf5_irec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
      call shdf5_irec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
      call shdf5_irec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
      call shdf5_irec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_sfc)

      if (nlfpats > 0) then

         ndims    = 1
         idims(1) = nlfpats

         call shdf5_irec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

         ndims    = 2
         idims(1) = 5
         idims(2) = nlfpats

         call shdf5_irec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
         call shdf5_irec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
         call shdf5_irec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

      endif

   endif

! Check whether NUDGING arrays are used

   if (mdomain == 0 .and. nudflag > 0) then

      ndims    = 1
      idims(1) = 1

      call shdf5_irec(ndims, idims, 'NUDNXP' , ivars=nudnxp)
      call shdf5_irec(ndims, idims, 'NWNUD'  , ivars=nwnud)

      mwnud = nwnud

      ndims    = 1
      idims(1) = nwnud

      call alloc_nudge1(nwnud)

      call shdf5_irec(ndims, idims, 'XEWNUD'  , rvara=xewnud)
      call shdf5_irec(ndims, idims, 'YEWNUD'  , rvara=yewnud)
      call shdf5_irec(ndims, idims, 'ZEWNUD'  , rvara=zewnud)

      call shdf5_irec(ndims,idims,'itab_wnud%npoly' ,ivara=itab_wnud(:)%npoly)

      allocate (iscr(6,nwnud))

      ndims    = 2
      idims(1) = 6
      idims(2) = nwnud

      call shdf5_irec(ndims,idims,'itab_wnud%iwnud',ivara=iscr)

      do iwnud = 1,nwnud
         itab_wnud(iwnud)%iwnud(1:6) = iscr(1:6,iwnud)
      enddo

      deallocate(iscr)

   endif

! Close the GRIDFILE

   call shdf5_close()

else

! Grid file does not exist.

   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(io6,*) '!!!  Gridfile does not exist:'
   write(io6,*) '!!!  '//trim(gridfile)
   write(io6,*) '!!!  Stopping run'
   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   
   stop 'stop - no gridfile'
   
endif

write(io6,*) 'end of gridfile_read '

return
end subroutine gridfile_read



