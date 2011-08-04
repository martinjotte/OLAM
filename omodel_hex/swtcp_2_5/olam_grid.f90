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

use misc_coms,   only: io6, runtype, mdomain, ngrids, initial, &
                       nxp, nzp, timmax8, alloc_misc, iparallel, meshtype, &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim

use leaf_coms,   only: nzg, nzs, isfcl, nwl, nslcon
use sea_coms,    only: nws

use mem_ijtabs,  only: istp, mrls, fill_jtabs, itab_md, itab_ud, itab_wd
use oplot_coms,  only: op
use mem_grid,    only: nza, nma, nua, nva, nwa, zm, zt, xem, yem, zem, &
                       alloc_grid1, alloc_grid2
use mem_sflux,   only: nlandflux, nseaflux
use mem_nudge,   only: nudnxp, nnudp, xenudp, yenudp, zenudp, alloc_nudge1

use mem_wm5_refsoln_cubic

implicit none

real, allocatable :: quarter_kite(:,:)

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

      write(io6,'(/,a)') 'gridinit calling gridset'
      call gridset()

      write(io6,'(a,f8.1)')    ' Model top height = ',zm(nza-1)

! Horizontal grid setup

   if (mdomain == 0) then
   
! If using nudging on global grid, generate nudging grid here

      if (nudnxp > 0) then   
         write(io6,'(/,a)') 'gridinit calling icosahedron for nudging grid'
         call icosahedron(nudnxp)  ! global spherical domain; calls 2 allocs

! Allocate nudging grid structure arrays, copy temporary grid structure to 
! these arrays, and deallocate temporary arrays 

         nnudp = nma
         call alloc_nudge1(nnudp)
         
         xenudp(:) = xem(:)
         yenudp(:) = yem(:)
         zenudp(:) = zem(:)

         deallocate (itab_md, itab_ud, itab_wd)
         deallocate (xem, yem, zem)      
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

   call alloc_grid1(meshtype)

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

   call alloc_grid2(meshtype)

! Set up control volumes

   write(io6,'(/,a)') 'gridinit calling ctrlvols'

   if (meshtype == 1) then
      call ctrlvols_tri()
   else
      call ctrlvols_hex(quarter_kite)
      deallocate (quarter_kite)
   endif

! Write GRIDFILE

   write(io6,'(/,a)') 'gridinit calling gridfile_write'
   call gridfile_write()

else

! Read atmos grid for INITIAL/HISTORY/PLOTONLY/PARCOMBINE run

   write(io6,'(/,a)') 'gridinit calling gridfile_read'
   call gridfile_read()

!------------------------------------------------------------
! SPECIAL - Call fill_wm5 to read in and initialize reference
! solution for time = 0 and time = 15d.

   if (nslcon > 2) then
      call fill_wm5()
   endif
!------------------------------------------------------------

endif

if (isfcl == 1) then
   write(io6,'(/,a)') 'gridinit before return to olam_run'
   write(io6,'(a,i8)')   ' nwl       = ',nwl
   write(io6,'(a,i8)')   ' nws       = ',nws
   write(io6,'(a,i8)')   ' nlandflux = ',nlandflux
   write(io6,'(a,i8)')   ' nseaflux  = ',nseaflux
endif

return
end subroutine gridinit

!===============================================================================

subroutine gridset()

use mem_grid,    only: nza, mza, &
                       zm, zt, dzm, dzt, dzim, dzit, &
                       zfacm, zfact, zfacim, zfacit, &
                       alloc_gridz
use misc_coms,   only: io6, nzp, deltaz, dzrat, dzmax, ztop, zbase, mdomain
use consts_coms, only: erad
use oname_coms,  only: nl

implicit none

integer :: ifm,icm,k,iinc,icnt,if1,jinc  &
   ,jcnt,jf,kinc,kcnt,kf,nrat,i,j,kcy,kcw,kk,ng,npg
integer :: nidiag,ijcorner,numd
real :: centx1,centy1,centx,centy,dzr,dsum,dzrcm,dzrfm,tsum,dzrati
real :: dxmax
real, allocatable, dimension(:) :: zmvec,ztvec

nza = nzp
mza = nzp

call alloc_gridz()
allocate (zmvec(-1:mza+1),ztvec(-1:mza+1))

! calculate zm

if ( deltaz < spacing(0.) ) then
   zmvec(1:nzp) = nl%zz(1:nzp)
   zmvec(nzp+1) = 2. * zmvec(nzp) - zmvec(nzp-1)
else
   zmvec(1) = zbase
   zmvec(2) = zbase + deltaz
   do k = 3,nzp+1
      zmvec(k) = zmvec(k-1)  &
               + min(dzrat * (zmvec(k-1) - zmvec(k-2)),max(deltaz,dzmax))
   enddo
endif
dzrati = (zmvec(2) - zmvec(1)) / (zmvec(3) - zmvec(2))
zmvec(0) = zmvec(1) - (zmvec(2) - zmvec(1)) * dzrati
zmvec(-1) = zmvec(0) - (zmvec(1) - zmvec(0)) * dzrati

ztop = zmvec(nzp-1)

! compute zt values by geometric interpolation.

do k = 1,nza
   dzrfm = sqrt(sqrt((zmvec(k+1) - zmvec(k)) / (zmvec(k-1) - zmvec(k-2))))
   ztvec(k) = zmvec(k-1) + (zmvec(k) - zmvec(k-1)) / (1. + dzrfm)
enddo
ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

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
end subroutine gridset

!===============================================================================

subroutine topo_init(nqa,topq,glatq,glonq,xeq,yeq,zeq)

use misc_coms,   only: io6, deltax
use consts_coms, only: pi1, pio180

! The following NSLCON is used as a namelist flag to indicate shallow water 
! test case 1, 2, or 5

use leaf_coms, only: nslcon

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

!   r = sqrt((glonq(iq) * pio180 + 0.5 * pi1) ** 2  &
!          + (glatq(iq) * pio180 - pi1 / 6.) ** 2)

!   topq(iq) = max(0., 2000. * (1. - r / r0))
   
if (nslcon > 2) then

   r = sqrt((glonq(iq) * pio180 + 0.5 * pi1) ** 2  &
          + (glatq(iq) * pio180 - pi1 / 6.) ** 2)

   topq(iq) = max(0., 2000. * (1. - r / r0))
   
endif

!   print*, 'topoinit ',iq,r,r0,topq(iq)
   

! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine topo_init

!===============================================================================

subroutine gridfile_write()

use max_dims,   only: maxngrdll
use misc_coms,  only: io6, ngrids, gridfile, mdomain, meshtype, nzp, nxp,  &
                      iclobber, itopoflg,  &
                      deltax, deltaz, dzmax, dzrat, zbase,  &
                      ngrdll, grdrad, grdlat, grdlon, meshtype
use mem_ijtabs, only: mloops_m, mloops_u, mloops_v, mloops_w, mrls,  &
                      itab_m, itab_u, itab_v, itab_w
use mem_grid,   only: nza, nma, nua, nva, nwa, nsw_max,  &
                      zm, zt, dzm, dzt, dzim, dzit,  &
                      zfacm, zfact, zfacim, zfacit, &
                      lpm, lpu, lcu, lpv, lcv, lpw, lsw,  &
                      topm, topw, xem, yem, zem, xeu, yeu, zeu,  &
                      xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz,  &
                      dnu, dniu, dnv, dniv, arw0, arm0,  &
                      glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
                      aru, arv, volui, volvi, arw, volwi, volt, volti
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, seaflux, landflux,  &
                      nsfpats, nlfpats, nsfpatm, nlfpatm,  &
                      xemsfpat, yemsfpat, zemsfpat,  &
                      xemlfpat, yemlfpat, zemlfpat

use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

implicit none

! This routine writes the grid variables to the grid file.

integer :: im, iu, iv, iw

integer :: ndims, idims(2)

! Scratch arrays for copying output

logical, allocatable :: lscr(:,:)
integer, allocatable :: iscr1(:,:),iscr2(:,:),iscr3(:,:),iscr4(:,:)

real, allocatable :: rscr1(:,:),rscr2(:,:),rscr3(:,:),rscr4(:,:),rscr5(:,:), &
                     rscr6(:,:),rscr7(:,:),rscr8(:,:),rscr9(:,:),rscr10(:,:), &
                     rscr11(:,:),rscr12(:,:),rscr13(:,:),rscr14(:,:), &
                     rscr15(:,:),rscr16(:,:)

write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(io6,*) 'grid_write: opening file:', trim(gridfile)
write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

call shdf5_open(gridfile,'W',iclobber)

! Write the gridfile information that exists in namelist

ndims = 1
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
call shdf5_orec(ndims, idims, 'DELTAZ'  , rvars=deltaz)
call shdf5_orec(ndims, idims, 'DZRAT'   , rvars=dzrat)
call shdf5_orec(ndims, idims, 'DZMAX'   , rvars=dzmax)
call shdf5_orec(ndims, idims, 'ZBASE'   , rvars=zbase)

idims(1) = ngrids

call shdf5_orec(ndims, idims, 'NGRDLL' , ivara=ngrdll)
call shdf5_orec(ndims, idims, 'GRDRAD' , rvara=grdrad)

ndims = 2
idims(1) = maxngrdll
idims(2) = ngrids

call shdf5_orec(ndims, idims, 'GRDLAT', rvara=grdlat)
call shdf5_orec(ndims, idims, 'GRDLON', rvara=grdlon)

! Write the grid dimensions

ndims = 1
idims(1) = 1

call shdf5_orec(ndims, idims, 'NZA'    , ivars=nza)
call shdf5_orec(ndims, idims, 'NMA'    , ivars=nma)
call shdf5_orec(ndims, idims, 'NUA'    , ivars=nua)
call shdf5_orec(ndims, idims, 'NVA'    , ivars=nva)
call shdf5_orec(ndims, idims, 'NWA'    , ivars=nwa)
call shdf5_orec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
call shdf5_orec(ndims, idims, 'MRLS'   , ivars=mrls)

! Write grid structure variables

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

   idims(1) = nua

   call shdf5_orec(ndims, idims, 'LPU'  , ivara=lpu)
   call shdf5_orec(ndims, idims, 'LCU'  , ivara=lcu)
   call shdf5_orec(ndims, idims, 'XEU'  , rvara=xeu)
   call shdf5_orec(ndims, idims, 'YEU'  , rvara=yeu)
   call shdf5_orec(ndims, idims, 'ZEU'  , rvara=zeu)
   call shdf5_orec(ndims, idims, 'GLATU', rvara=glatu)
   call shdf5_orec(ndims, idims, 'GLONU', rvara=glonu)

else

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

ndims = 2
idims(1) = nza

if (meshtype == 1) then

   idims(2) = nua

   call shdf5_orec(ndims, idims, 'VOLUI', rvara=volui)

else

   idims(2) = nva

   call shdf5_orec(ndims, idims, 'ARV'  , rvara=arv)
   call shdf5_orec(ndims, idims, 'VOLVI', rvara=volvi)

endif

call shdf5_orec(ndims, idims, 'ARU'  , rvara=aru)

idims(2) = nwa

call shdf5_orec(ndims, idims, 'ARW'  , rvara=arw)
call shdf5_orec(ndims, idims, 'VOLWI', rvara=volwi)
call shdf5_orec(ndims, idims, 'VOLT' , dvara=volt)
call shdf5_orec(ndims, idims, 'VOLTI', dvara=volti)

! Write ITAB_M SCALARS

ndims = 1
idims(1) = nma
idims(2) = 1

call shdf5_orec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
call shdf5_orec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm)
call shdf5_orec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe)
call shdf5_orec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
call shdf5_orec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
call shdf5_orec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
call shdf5_orec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh)

! Write ITAB_M ARRAYS

allocate (lscr(mloops_m,nma))
allocate (iscr1(7,nma))
allocate (iscr2(7,nma))
allocate (iscr3(7,nma))

allocate (rscr1(7,nma))

do im = 1,nma
   lscr(1:mloops_m,im) = itab_m(im)%loop(1:mloops_m)
   iscr1(1:7,im) = itab_m(im)%iu(1:7)
   iscr2(1:7,im) = itab_m(im)%iv(1:7)
   iscr3(1:7,im) = itab_m(im)%iw(1:7)

   rscr1(1:7,im) = itab_m(im)%fmw(1:7)
enddo

ndims = 2
idims(1) = mloops_m
idims(2) = nma

call shdf5_orec(ndims,idims,'itab_m%loop',lvara=lscr)

idims(1) = 7

if (meshtype == 1) then
   call shdf5_orec(ndims,idims,'itab_m%iu',ivara=iscr1)
else
   call shdf5_orec(ndims,idims,'itab_m%iv' ,ivara=iscr2)
endif

call shdf5_orec(ndims,idims,'itab_m%iw',ivara=iscr3)
call shdf5_orec(ndims,idims,'itab_m%fmw',rvara=rscr1)

deallocate (lscr,iscr1,iscr2,iscr3,rscr1)

if (meshtype == 1) then

! Write ITAB_U SCALARS

   ndims = 1
   idims(1) = nua
   idims(2) = 1

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

   allocate (lscr(mloops_u,nua))

   allocate (iscr1( 2,nua))
   allocate (iscr2(12,nua))
   allocate (iscr3( 6,nua))

   allocate (rscr1 ( 4,nua))
   allocate (rscr2 (12,nua))
   allocate (rscr3 ( 6,nua))
   allocate (rscr4 ( 4,nua))
   allocate (rscr5 ( 4,nua))
   allocate (rscr6 ( 4,nua))
   allocate (rscr7 ( 2,nua))
   allocate (rscr8 ( 2,nua))
   allocate (rscr9 ( 4,nua))

   do iu = 1,nua
      lscr(1:mloops_u,iu) = itab_u(iu)%loop(1:mloops_u)

      iscr1(1: 2,iu) = itab_u(iu)%im(1: 2)
      iscr2(1:12,iu) = itab_u(iu)%iu(1:12)
      iscr3(1: 6,iu) = itab_u(iu)%iw(1: 6)

      rscr1 (1: 4,iu) = itab_u(iu)%diru (1: 4)
      rscr2 (1:12,iu) = itab_u(iu)%fuu  (1:12)
      rscr3 (1: 6,iu) = itab_u(iu)%fuw  (1: 6)
      rscr4 (1: 4,iu) = itab_u(iu)%tuu  (1: 4)
      rscr5 (1: 4,iu) = itab_u(iu)%vxu_u(1: 4)
      rscr6 (1: 4,iu) = itab_u(iu)%vyu_u(1: 4)
      rscr7 (1: 2,iu) = itab_u(iu)%vxw_u(1: 2)
      rscr8 (1: 2,iu) = itab_u(iu)%vyw_u(1: 2)
      rscr9 (1: 4,iu) = itab_u(iu)%guw  (1: 4)
   enddo

   ndims = 2
   idims(1) = mloops_u
   idims(2) = nua

   call shdf5_orec(ndims,idims,'itab_u%loop',lvara=lscr)

   idims(1) = 2

   call shdf5_orec(ndims,idims,'itab_u%im'   ,ivara=iscr1)
   call shdf5_orec(ndims,idims,'itab_u%vxw_u',rvara=rscr7)
   call shdf5_orec(ndims,idims,'itab_u%vyw_u',rvara=rscr8)

   idims(1) = 4

   call shdf5_orec(ndims,idims,'itab_u%diru' ,rvara=rscr1)
   call shdf5_orec(ndims,idims,'itab_u%tuu'  ,rvara=rscr4)
   call shdf5_orec(ndims,idims,'itab_u%guw'  ,rvara=rscr9)
   call shdf5_orec(ndims,idims,'itab_u%vxu_u',rvara=rscr5)
   call shdf5_orec(ndims,idims,'itab_u%vyu_u',rvara=rscr6)

   idims(1) = 6

   call shdf5_orec(ndims,idims,'itab_u%fuw',rvara=rscr3)
   call shdf5_orec(ndims,idims,'itab_u%iw' ,ivara=iscr3)

   idims(1) = 12

   call shdf5_orec(ndims,idims,'itab_u%iu' ,ivara=iscr2)
   call shdf5_orec(ndims,idims,'itab_u%fuu',rvara=rscr2)

   deallocate (lscr,iscr1,iscr2,iscr3,  &
               rscr1,rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9)

endif

if (meshtype == 2) then

! Write ITAB_V SCALARS

   ndims = 1
   idims(1) = nva
   idims(2) = 1

   call shdf5_orec(ndims,idims,'itab_v%ivp'    ,ivara=itab_v(:)%ivp)
   call shdf5_orec(ndims,idims,'itab_v%mrlv'   ,ivara=itab_v(:)%mrlv)
   call shdf5_orec(ndims,idims,'itab_v%ivglobe',ivara=itab_v(:)%ivglobe)

! Write ITAB_V ARRAYS

   allocate (lscr(mloops_v,nva))

   allocate (iscr1( 6,nva))
   allocate (iscr2(16,nva))
   allocate (iscr3( 4,nva))

   allocate (rscr2(12,nva))
   allocate (rscr3( 4,nva))
   allocate (rscr4(16,nva))
   allocate (rscr5( 2,nva))
   allocate (rscr6( 2,nva))
   allocate (rscr7( 2,nva))
   allocate (rscr8( 2,nva))
   allocate (rscr9( 2,nva))

   do iv = 1,nva
      lscr(1:mloops_v,iv) = itab_v(iv)%loop(1:mloops_v)

      iscr1(1: 6,iv) = itab_v(iv)%im(1: 6)
      iscr2(1:16,iv) = itab_v(iv)%iv(1:16)
      iscr3(1: 4,iv) = itab_v(iv)%iw(1: 4)

      rscr2(1:12,iv) = itab_v(iv)%fvv (1:12)
      rscr3(1: 4,iv) = itab_v(iv)%fvw (1: 4)
      rscr4(1:16,iv) = itab_v(iv)%fuv (1:16)
      rscr5(1: 2,iv) = itab_v(iv)%farw(1: 2)
      rscr6(1: 2,iv) = itab_v(iv)%cosv(1: 2)
      rscr7(1: 2,iv) = itab_v(iv)%sinv(1: 2)
      rscr8(1: 2,iv) = itab_v(iv)%dxps(1: 2)
      rscr9(1: 2,iv) = itab_v(iv)%dyps(1: 2)
   enddo

   ndims = 2
   idims(1) = mloops_v
   idims(2) = nva

   call shdf5_orec(ndims,idims,'itab_v%loop',lvara=lscr)

   idims(1) = 2

   call shdf5_orec(ndims,idims,'itab_v%farw',rvara=rscr5)
   call shdf5_orec(ndims,idims,'itab_v%cosv',rvara=rscr6)
   call shdf5_orec(ndims,idims,'itab_v%sinv',rvara=rscr7)
   call shdf5_orec(ndims,idims,'itab_v%dxps',rvara=rscr8)
   call shdf5_orec(ndims,idims,'itab_v%dyps',rvara=rscr9)

   idims(1) = 4

   call shdf5_orec(ndims,idims,'itab_v%iw'  ,ivara=iscr3)
   call shdf5_orec(ndims,idims,'itab_v%fvw' ,rvara=rscr3)

   idims(1) = 6

   call shdf5_orec(ndims,idims,'itab_v%im'  ,ivara=iscr1)

   idims(1) = 12

   call shdf5_orec(ndims,idims,'itab_v%fvv' ,rvara=rscr2)

   idims(1) = 16

   call shdf5_orec(ndims,idims,'itab_v%iv'  ,ivara=iscr2)
   call shdf5_orec(ndims,idims,'itab_v%fuv' ,rvara=rscr4)

   deallocate (lscr,iscr1,iscr2,iscr3, &
               rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9)

endif

! Write ITAB_W SCALARS

ndims = 1
idims(1) = nwa
idims(2) = 1

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

allocate (lscr(mloops_w,nwa))

allocate (iscr1(7,nwa))
allocate (iscr2(9,nwa))
allocate (iscr3(7,nwa))
allocate (iscr4(9,nwa))

allocate (rscr1 (3,nwa))
allocate (rscr2 (7,nwa))
allocate (rscr3 (7,nwa))
allocate (rscr4 (7,nwa))
allocate (rscr5 (9,nwa))
allocate (rscr6 (3,nwa))
allocate (rscr7 (3,nwa))
allocate (rscr8 (3,nwa))
allocate (rscr9 (3,nwa))
allocate (rscr10(3,nwa))
allocate (rscr11(7,nwa))
allocate (rscr12(7,nwa))
allocate (rscr13(7,nwa))
allocate (rscr14(7,nwa))
allocate (rscr15(7,nwa))
allocate (rscr16(7,nwa))

do iw = 1,nwa
   lscr(1:mloops_w,iw) = itab_w(iw)%loop(1:mloops_w)

   iscr1(1:7,iw) = itab_w(iw)%im(1:7)
   iscr2(1:9,iw) = itab_w(iw)%iu(1:9)
   iscr3(1:7,iw) = itab_w(iw)%iv(1:7)
   iscr4(1:9,iw) = itab_w(iw)%iw(1:9)

   rscr1 (1:3,iw) = itab_w(iw)%diru (1:3)
   rscr2 (1:7,iw) = itab_w(iw)%dirv (1:7)
   rscr3 (1:7,iw) = itab_w(iw)%fwv  (1:7)
   rscr4 (1:7,iw) = itab_w(iw)%fww  (1:7)

   rscr5 (1:9,iw) = itab_w(iw)%fwu  (1:9)
   rscr6 (1:3,iw) = itab_w(iw)%vxu  (1:3)
   rscr7 (1:3,iw) = itab_w(iw)%vyu  (1:3)
   rscr8 (1:3,iw) = itab_w(iw)%vzu  (1:3)
   rscr9 (1:3,iw) = itab_w(iw)%vxu_w(1:3)
   rscr10(1:3,iw) = itab_w(iw)%vyu_w(1:3)

   rscr11(1:7,iw) = itab_w(iw)%farm (1:7)
   rscr12(1:7,iw) = itab_w(iw)%farv (1:7)

   rscr13(1:7,iw) = itab_w(iw)%gxps1(1:7)
   rscr14(1:7,iw) = itab_w(iw)%gyps1(1:7)
   rscr15(1:7,iw) = itab_w(iw)%gxps2(1:7)
   rscr16(1:7,iw) = itab_w(iw)%gyps2(1:7)

enddo

ndims = 2
idims(1) = mloops_w
idims(2) = nwa

call shdf5_orec(ndims,idims,'itab_w%loop',lvara=lscr)

idims(1) = 3

call shdf5_orec(ndims,idims,'itab_w%diru' ,rvara=rscr1)
call shdf5_orec(ndims,idims,'itab_w%vxu'  ,rvara=rscr6)
call shdf5_orec(ndims,idims,'itab_w%vyu'  ,rvara=rscr7)
call shdf5_orec(ndims,idims,'itab_w%vzu'  ,rvara=rscr8)
call shdf5_orec(ndims,idims,'itab_w%vxu_w',rvara=rscr9)
call shdf5_orec(ndims,idims,'itab_w%vyu_w',rvara=rscr10)

idims(1) = 7

call shdf5_orec(ndims,idims,'itab_w%im'  ,ivara=iscr1)
call shdf5_orec(ndims,idims,'itab_w%iv'  ,ivara=iscr3)

call shdf5_orec(ndims,idims,'itab_w%dirv',rvara=rscr2)
call shdf5_orec(ndims,idims,'itab_w%fwv' ,rvara=rscr3)
call shdf5_orec(ndims,idims,'itab_w%fww' ,rvara=rscr4)
call shdf5_orec(ndims,idims,'itab_w%farm',rvara=rscr11)
call shdf5_orec(ndims,idims,'itab_w%farv',rvara=rscr12)

call shdf5_orec(ndims,idims,'itab_w%gxps1',rvara=rscr13)
call shdf5_orec(ndims,idims,'itab_w%gyps1',rvara=rscr14)
call shdf5_orec(ndims,idims,'itab_w%gxps2',rvara=rscr15)
call shdf5_orec(ndims,idims,'itab_w%gyps2',rvara=rscr16)

idims(1) = 9

call shdf5_orec(ndims,idims,'itab_w%iu'  ,ivara=iscr2)
call shdf5_orec(ndims,idims,'itab_w%iw'  ,ivara=iscr4)
call shdf5_orec(ndims,idims,'itab_w%fwu' ,rvara=rscr5)

deallocate(lscr,iscr1,iscr2,iscr3, &
    rscr1,rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9,rscr10, &
    rscr11,rscr12,rscr13,rscr14,rscr15,rscr16)

! Check whether LAND/SEA models are used

if (isfcl == 1) then

! Write SEAFLUX VALUES

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_orec(ndims, idims, 'NSEAFLUX' ,ivars=nseaflux)
   call shdf5_orec(ndims, idims, 'NSFPATS',ivars=nsfpats)

   idims(1) = nseaflux

   call shdf5_orec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%isfglobe)
   call shdf5_orec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
   call shdf5_orec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
   call shdf5_orec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iws)
   call shdf5_orec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
   call shdf5_orec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
   call shdf5_orec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
   call shdf5_orec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
   call shdf5_orec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
   call shdf5_orec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
   call shdf5_orec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
   call shdf5_orec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sea)

   if (nsfpats > 0) then

      idims(1) = nsfpats

      call shdf5_orec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

      ndims = 2
      idims(1) = 5
      idims(2) = nsfpats

      call shdf5_orec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
      call shdf5_orec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
      call shdf5_orec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

   endif

! Write LANDFLUX VALUES

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_orec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
   call shdf5_orec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

   idims(1) = nlandflux

   call shdf5_orec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ilfglobe)
   call shdf5_orec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
   call shdf5_orec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
   call shdf5_orec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwl)
   call shdf5_orec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
   call shdf5_orec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
   call shdf5_orec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
   call shdf5_orec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
   call shdf5_orec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
   call shdf5_orec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
   call shdf5_orec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
   call shdf5_orec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_land)

   if (nlfpats > 0) then

      idims(1) = nlfpats

      call shdf5_orec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

      ndims = 2
      idims(1) = 5
      idims(2) = nlfpats

      call shdf5_orec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
      call shdf5_orec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
      call shdf5_orec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

   endif

endif

! Close GRIDFILE

call shdf5_close()

return
end subroutine gridfile_write

!===============================================================================

subroutine gridfile_read()

use max_dims,   only: maxngrdll
use misc_coms,  only: io6, ngrids, gridfile, mdomain, meshtype, nzp, nxp,  &
                      itopoflg,  &
                      deltax, deltaz, dzmax, dzrat, zbase,  &
                      ngrdll, grdrad, grdlat, grdlon, meshtype
use mem_ijtabs, only: mloops_m, mloops_u, mloops_v, mloops_w, mrls,  &
                      itab_m, itab_u, itab_v, itab_w, alloc_itabs
use mem_grid,   only: nza, nma, nua, nva, nwa,  &
                      mza, mma, mua, mva, mwa, nsw_max,  &
                      zm, zt, dzm, dzt, dzim, dzit,  &
                      zfacm, zfact, zfacim, zfacit, &
                      lpm, lpu, lcu, lpv, lcv, lpw, lsw,  &
                      topm, topw, xem, yem, zem, xeu, yeu, zeu,  &
                      xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz,  &
                      dnu, dniu, dnv, dniv, arw0, arm0,  &
                      glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
                      aru, arv, volui, volvi, arw, volwi, volt, volti, &
                      alloc_gridz, alloc_xyzem, alloc_xyzew,  &
                      alloc_grid1, alloc_grid2
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, mseaflux, mlandflux,  &
                      nsfpats, nlfpats, msfpats, mlfpats, nsfpatm, nlfpatm,  &
                      seaflux, landflux,  &
                      xemsfpat, yemsfpat, zemsfpat,  &
                      xemlfpat, yemlfpat, zemlfpat

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_para,   only: myrank

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iu, iv, iw

integer :: ierr

integer :: ngr, i
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, meshtype0, nxp0, nzp0, itopoflg0, isfcl0

real :: deltax0, deltaz0, dzrat0, dzmax0, zbase0

character(len=128) :: flnm
character(len=10) :: number

logical  :: exans

integer :: ngrdll0(ngrids)

real :: grdrad0(ngrids)
real :: grdlat0(maxngrdll,ngrids)
real :: grdlon0(maxngrdll,ngrids)

! Scratch arrays for copying input

logical, allocatable :: lscr(:,:)
integer, allocatable :: iscr1(:,:),iscr2(:,:),iscr3(:,:),iscr4(:,:)

real, allocatable :: rscr1(:,:),rscr2(:,:),rscr3(:,:),rscr4(:,:),rscr5(:,:),  &
                     rscr6(:,:),rscr7(:,:),rscr8(:,:),rscr9(:,:),rscr10(:,:), &
                     rscr11(:,:),rscr12(:,:),rscr13(:,:),rscr14(:,:), &
                     rscr15(:,:),rscr16(:,:)

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
   call shdf5_irec(ndims, idims, 'DELTAZ'  , rvars=deltaz0)
   call shdf5_irec(ndims, idims, 'DZRAT'   , rvars=dzrat0)
   call shdf5_irec(ndims, idims, 'DZMAX'   , rvars=dzmax0)
   call shdf5_irec(ndims, idims, 'ZBASE'   , rvars=zbase0)

   idims(1) = ngrids0

   call shdf5_irec(ndims, idims, 'NGRDLL' , ivara=ngrdll0)
   call shdf5_irec(ndims, idims, 'GRDRAD' , rvara=grdrad0)

   ndims = 2
   idims(1) = maxngrdll
   idims(2) = ngrids0

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

   if (abs(deltax0 - deltax) > 1.e-3) ierr = 1 
   if (abs(deltaz0 - deltaz) > 1.e-3) ierr = 1 
   if (abs(dzrat0  - dzrat ) > 1.e-3) ierr = 1 
   if (abs(dzmax0  - dzmax ) > 1.e-3) ierr = 1 
   if (abs(zbase0  - zbase ) > 1.e-3) ierr = 1 

   do ngr = 1,ngrids0
      if (abs(ngrdll0 (ngr) - ngrdll (ngr)) > 1.e1 ) ierr = 1
      if (abs(grdrad0 (ngr) - grdrad (ngr)) > 1.e1 ) ierr = 1

      do i = 1,ngrdll0(ngr)
         if (abs(grdlat0(i,ngr) - grdlat(i,ngr)) > 1.e-3) ierr = 1
         if (abs(grdlon0(i,ngr) - grdlon(i,ngr)) > 1.e-3) ierr = 1
      enddo
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
      write(io6,*)              'deltaz:   ',deltaz0  ,deltaz
      write(io6,*)              'dzrat:    ',dzrat0   ,dzrat
      write(io6,*)              'dzmax:    ',dzmax0   ,dzmax
      write(io6,*)              'zbase:    ',zbase0   ,zbase
      write(io6,*) ' '
      write(io6, '(a,20i12)')   'ngrdll0:  ',ngrdll0 (1:ngrids)
      write(io6, '(a,20i12)')   'ngrdll:   ',ngrdll  (1:ngrids)
      write(io6,*) ' '
      write(io6, '(a,20f12.1)') 'grdrad0:  ',grdrad0 (1:ngrids)
      write(io6, '(a,20f12.1)') 'grdrad:   ',grdrad  (1:ngrids)
      write(io6,*) ' '

      do ngr = 1,ngrids
         write(io6, '(a,i5)') 'ngr: ',ngr
         write(io6,*) ' '
         write(io6, '(a,20f10.3)') 'grdlat0: ',grdlat0(1:ngrdll(ngr),ngr)
         write(io6, '(a,20f10.3)') 'grdlat:  ',grdlat (1:ngrdll(ngr),ngr)
         write(io6,*) ' '
         write(io6, '(a,20f10.3)') 'grdlon0: ',grdlon0(1:ngrdll(ngr),ngr)
         write(io6, '(a,20f10.3)') 'grdlon:  ',grdlon (1:ngrdll(ngr),ngr)
         write(io6,*) ' '
      enddo

      write(io6,*) '-----------------------------------------------'

      stop 'stop - gridfile mismatch'
   
   endif

! Read the grid dimensions

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

   call alloc_gridz()
   call alloc_itabs(meshtype,nma,nua,nva,nwa)
   call alloc_xyzem(nma)
   call alloc_xyzew(nwa)
   call alloc_grid1(meshtype)
   call alloc_grid2(meshtype)

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

   idims(1) = nma

   call shdf5_irec(ndims, idims, 'LPM'  , ivara=lpm)
   call shdf5_irec(ndims, idims, 'TOPM' , rvara=topm)
   call shdf5_irec(ndims, idims, 'XEM'  , rvara=xem)
   call shdf5_irec(ndims, idims, 'YEM'  , rvara=yem)
   call shdf5_irec(ndims, idims, 'ZEM'  , rvara=zem)
   call shdf5_irec(ndims, idims, 'ARM0' , rvara=arm0)
   call shdf5_irec(ndims, idims, 'GLATM', rvara=glatm)
   call shdf5_irec(ndims, idims, 'GLONM', rvara=glonm)

   if (meshtype == 1) then

      idims(1) = nua

      call shdf5_irec(ndims, idims, 'LPU'  , ivara=lpu)
      call shdf5_irec(ndims, idims, 'LCU'  , ivara=lcu)
      call shdf5_irec(ndims, idims, 'XEU'  , rvara=xeu)
      call shdf5_irec(ndims, idims, 'YEU'  , rvara=yeu)
      call shdf5_irec(ndims, idims, 'ZEU'  , rvara=zeu)
      call shdf5_irec(ndims, idims, 'GLATU', rvara=glatu)
      call shdf5_irec(ndims, idims, 'GLONU', rvara=glonu)

   else

      idims(1) = nva

      call shdf5_irec(ndims, idims, 'LPV' , ivara=lpv)
      call shdf5_irec(ndims, idims, 'LCV' , ivara=lcv)
      call shdf5_irec(ndims, idims, 'XEV' , rvara=xev)
      call shdf5_irec(ndims, idims, 'YEV' , rvara=yev)
      call shdf5_irec(ndims, idims, 'ZEV' , rvara=zev)
      call shdf5_irec(ndims, idims, 'GLATV', rvara=glatv)
      call shdf5_irec(ndims, idims, 'GLONV', rvara=glonv)

   endif

   call shdf5_irec(ndims, idims, 'DNU' , rvara=dnu)
   call shdf5_irec(ndims, idims, 'DNV' , rvara=dnv)
   call shdf5_irec(ndims, idims, 'DNIU', rvara=dniu)
   call shdf5_irec(ndims, idims, 'DNIV', rvara=dniv)
   call shdf5_irec(ndims, idims, 'UNX' , rvara=unx)
   call shdf5_irec(ndims, idims, 'UNY' , rvara=uny)
   call shdf5_irec(ndims, idims, 'UNZ' , rvara=unz)
   call shdf5_irec(ndims, idims, 'VNX' , rvara=vnx)
   call shdf5_irec(ndims, idims, 'VNY' , rvara=vny)
   call shdf5_irec(ndims, idims, 'VNZ' , rvara=vnz)

   idims(1) = nwa

   call shdf5_irec(ndims, idims, 'LPW'  , ivara=lpw)
   call shdf5_irec(ndims, idims, 'LSW'  , ivara=lsw)
   call shdf5_irec(ndims, idims, 'XEW'  , rvara=xew)
   call shdf5_irec(ndims, idims, 'YEW'  , rvara=yew)
   call shdf5_irec(ndims, idims, 'ZEW'  , rvara=zew)
   call shdf5_irec(ndims, idims, 'TOPW' , rvara=topw)
   call shdf5_irec(ndims, idims, 'ARW0' , rvara=arw0)
   call shdf5_irec(ndims, idims, 'GLATW', rvara=glatw)
   call shdf5_irec(ndims, idims, 'GLONW', rvara=glonw)
   call shdf5_irec(ndims, idims, 'WNX'  , rvara=wnx)
   call shdf5_irec(ndims, idims, 'WNY'  , rvara=wny)
   call shdf5_irec(ndims, idims, 'WNZ'  , rvara=wnz)

   ndims = 2
   idims(1) = nza

   if (meshtype == 1) then

      idims(2) = nua

      call shdf5_irec(ndims, idims, 'VOLUI', rvara=volui)

   else

      idims(2) = nva

      call shdf5_irec(ndims, idims, 'ARV'  , rvara=arv)
      call shdf5_irec(ndims, idims, 'VOLVI', rvara=volvi)
   
   endif

   call shdf5_irec(ndims, idims, 'ARU'  , rvara=aru)

   idims(2) = nwa

   call shdf5_irec(ndims, idims, 'ARW'  , rvara=arw)
   call shdf5_irec(ndims, idims, 'VOLWI', rvara=volwi)
   call shdf5_irec(ndims, idims, 'VOLT' , dvara=volt)
   call shdf5_irec(ndims, idims, 'VOLTI', dvara=volti)
   
! Read ITAB_M SCALARS

   ndims = 1
   idims(1) = nma
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
   call shdf5_irec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm)
   call shdf5_irec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe)
   call shdf5_irec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
   call shdf5_irec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
   call shdf5_irec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
   call shdf5_irec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh)

! Read ITAB_M ARRAYS

   allocate (lscr(mloops_m,nma))
   allocate (iscr1(7,nma))
   allocate (iscr2(7,nma))
   allocate (iscr3(7,nma))

   allocate (rscr1(7,nma))

   ndims = 2
   idims(1) = mloops_m
   idims(2) = nma

   call shdf5_irec(ndims,idims,'itab_m%loop',lvara=lscr)

   idims(1) = 7

   if (meshtype == 1) then
      call shdf5_irec(ndims,idims,'itab_m%iu',ivara=iscr1)
   else
      call shdf5_irec(ndims,idims,'itab_m%iv' ,ivara=iscr2)
   endif

   call shdf5_irec(ndims,idims,'itab_m%iw' ,ivara=iscr3)
   call shdf5_irec(ndims,idims,'itab_m%fmw',rvara=rscr1)

   do im = 1,nma
      itab_m(im)%loop(1:mloops_m) = lscr(1:mloops_m,im)
      itab_m(im)%iu(1:7) = iscr1(1:7,im)
      itab_m(im)%iv(1:7) = iscr2(1:7,im)
      itab_m(im)%iw(1:7) = iscr3(1:7,im)

      itab_m(im)%fmw(1:7) = rscr1(1:7,im)
   enddo

   deallocate (lscr,iscr1,iscr2,iscr3,rscr1)

   if (meshtype == 1) then

! Read ITAB_U SCALARS

      ndims = 1
      idims(1) = nua
      idims(2) = 1

      call shdf5_irec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)
      call shdf5_irec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)
      call shdf5_irec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)

      call shdf5_irec(ndims,idims,'itab_u%gcf36'  ,rvara=itab_u(:)%gcf36)
      call shdf5_irec(ndims,idims,'itab_u%gcf45'  ,rvara=itab_u(:)%gcf45)

      call shdf5_irec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
      call shdf5_irec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
      call shdf5_irec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
      call shdf5_irec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
      call shdf5_irec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
      call shdf5_irec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
      call shdf5_irec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
      call shdf5_irec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)
      call shdf5_irec(ndims,idims,'itab_u%crossmm',rvara=itab_u(:)%crossmm)
      call shdf5_irec(ndims,idims,'itab_u%crossww',rvara=itab_u(:)%crossww)

! Read ITAB_U ARRAYS

      allocate (lscr(mloops_u,nua))

      allocate (iscr1( 2,nua))
      allocate (iscr2(12,nua))
      allocate (iscr3( 6,nua))

      allocate (rscr1 ( 4,nua))
      allocate (rscr2 (12,nua))
      allocate (rscr3 ( 6,nua))
      allocate (rscr4 ( 4,nua))
      allocate (rscr5 ( 4,nua))
      allocate (rscr6 ( 4,nua))
      allocate (rscr7 ( 2,nua))
      allocate (rscr8 ( 2,nua))
      allocate (rscr9 ( 4,nua))

      ndims = 2
      idims(1) = mloops_u
      idims(2) = nua

      call shdf5_irec(ndims,idims,'itab_u%loop',lvara=lscr)

      idims(1) = 2

      call shdf5_irec(ndims,idims,'itab_u%im'   ,ivara=iscr1)
      call shdf5_irec(ndims,idims,'itab_u%vxw_u',rvara=rscr7)
      call shdf5_irec(ndims,idims,'itab_u%vyw_u',rvara=rscr8)

      idims(1) = 2

      call shdf5_irec(ndims,idims,'itab_u%diru' ,rvara=rscr1)
      call shdf5_irec(ndims,idims,'itab_u%tuu'  ,rvara=rscr4)
      call shdf5_irec(ndims,idims,'itab_u%guw'  ,rvara=rscr9)
      call shdf5_irec(ndims,idims,'itab_u%vxu_u',rvara=rscr5)
      call shdf5_irec(ndims,idims,'itab_u%vyu_u',rvara=rscr6)

      idims(1) = 6

      call shdf5_irec(ndims,idims,'itab_u%fuw' ,rvara=rscr3)
      call shdf5_irec(ndims,idims,'itab_u%iw'  ,ivara=iscr3)

      idims(1) = 12

      call shdf5_irec(ndims,idims,'itab_u%iu' ,ivara=iscr2)
      call shdf5_irec(ndims,idims,'itab_u%fuu',rvara=rscr2)

      do iu = 1,nua
         itab_u(iu)%loop(1:mloops_u) = lscr(1:mloops_u,iu)

         itab_u(iu)%im(1: 2) = iscr1(1: 2,iu)
         itab_u(iu)%iu(1:12) = iscr2(1:12,iu)
         itab_u(iu)%iw(1: 6) = iscr3(1: 6,iu)

         itab_u(iu)%diru (1: 4) = rscr1 (1: 4,iu)
         itab_u(iu)%fuu  (1:12) = rscr2 (1:12,iu)
         itab_u(iu)%fuw  (1: 6) = rscr3 (1: 6,iu)
         itab_u(iu)%tuu  (1: 4) = rscr4 (1: 4,iu)
         itab_u(iu)%vxu_u(1: 4) = rscr5 (1: 4,iu)
         itab_u(iu)%vyu_u(1: 4) = rscr6 (1: 4,iu)
         itab_u(iu)%vxw_u(1: 2) = rscr7 (1: 2,iu)
         itab_u(iu)%vyw_u(1: 2) = rscr8 (1: 2,iu)
         itab_u(iu)%guw  (1: 4) = rscr9 (1: 4,iu)
      enddo

      deallocate (lscr,iscr1,iscr2,iscr3,  &
                  rscr1,rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9)

   endif

   if (meshtype == 2) then

! Read ITAB_V SCALARS

      ndims = 1
      idims(1) = nva
      idims(2) = 1

      call shdf5_irec(ndims,idims,'itab_v%ivp'      ,ivara=itab_v(:)%ivp)
      call shdf5_irec(ndims,idims,'itab_v%mrlv'     ,ivara=itab_v(:)%mrlv)
      call shdf5_irec(ndims,idims,'itab_v%ivglobe'  ,ivara=itab_v(:)%ivglobe)

! Read ITAB_V ARRAYS

      allocate (lscr(mloops_v,nva))

      allocate (iscr1( 6,nva))
      allocate (iscr2(16,nva))
      allocate (iscr3( 4,nva))

      allocate (rscr2(12,nva))
      allocate (rscr3( 4,nva))
      allocate (rscr4(16,nva))
      allocate (rscr5( 2,nva))
      allocate (rscr6( 2,nva))
      allocate (rscr7( 2,nva))
      allocate (rscr8( 2,nva))
      allocate (rscr9( 2,nva))

      ndims = 2
      idims(1) = mloops_v
      idims(2) = nva

      call shdf5_irec(ndims,idims,'itab_v%loop',lvara=lscr)

      idims(1) = 2

      call shdf5_irec(ndims,idims,'itab_v%farw',rvara=rscr5)
      call shdf5_irec(ndims,idims,'itab_v%cosv',rvara=rscr6)
      call shdf5_irec(ndims,idims,'itab_v%sinv',rvara=rscr7)
      call shdf5_irec(ndims,idims,'itab_v%dxps',rvara=rscr8)
      call shdf5_irec(ndims,idims,'itab_v%dyps',rvara=rscr9)

      idims(1) = 4

      call shdf5_irec(ndims,idims,'itab_v%iw'  ,ivara=iscr3)
      call shdf5_irec(ndims,idims,'itab_v%fvw' ,rvara=rscr3)

      idims(1) = 6

      call shdf5_irec(ndims,idims,'itab_v%im'  ,ivara=iscr1)

      idims(1) = 12

      call shdf5_irec(ndims,idims,'itab_v%fvv' ,rvara=rscr2)

      idims(1) = 16

      call shdf5_irec(ndims,idims,'itab_v%iv'  ,ivara=iscr2)
      call shdf5_irec(ndims,idims,'itab_v%fuv' ,rvara=rscr4)

      do iv = 1,nva
         itab_v(iv)%loop(1:mloops_v) = lscr(1:mloops_v,iv)

         itab_v(iv)%im(1: 6) = iscr1(1: 6,iv)
         itab_v(iv)%iv(1:16) = iscr2(1:16,iv)
         itab_v(iv)%iw(1: 4) = iscr3(1: 4,iv)

         itab_v(iv)%fvv (1:12) = rscr2(1:12,iv)
         itab_v(iv)%fvw (1: 4) = rscr3(1: 4,iv)
         itab_v(iv)%fuv (1:16) = rscr4(1:16,iv)
         itab_v(iv)%farw(1: 2) = rscr5(1: 2,iv)
         itab_v(iv)%cosv(1: 2) = rscr6(1: 2,iv)
         itab_v(iv)%sinv(1: 2) = rscr7(1: 2,iv)
         itab_v(iv)%dxps(1: 2) = rscr8(1: 2,iv)
         itab_v(iv)%dyps(1: 2) = rscr9(1: 2,iv)
      enddo

      deallocate (lscr,iscr1,iscr2,iscr3, &
                  rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9)

   endif

! Read ITAB_W SCALARS

   ndims = 1
   idims(1) = nwa
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly)
   call shdf5_irec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)
   call shdf5_irec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
   call shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
   call shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
   call shdf5_irec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
   call shdf5_irec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

   call shdf5_irec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)
   call shdf5_irec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)
   call shdf5_irec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

   call shdf5_irec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w)
   call shdf5_irec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w)

   call shdf5_irec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w)
   call shdf5_irec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w)
   call shdf5_irec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w)

! Read ITAB_W ARRAYS

   allocate (lscr(mloops_w,nwa))

   allocate (iscr1(7,nwa))
   allocate (iscr2(9,nwa))
   allocate (iscr3(7,nwa))
   allocate (iscr4(9,nwa))

   allocate (rscr1 (3,nwa))
   allocate (rscr2 (7,nwa))
   allocate (rscr3 (7,nwa))
   allocate (rscr4 (7,nwa))
   allocate (rscr5 (9,nwa))
   allocate (rscr6 (3,nwa))
   allocate (rscr7 (3,nwa))
   allocate (rscr8 (3,nwa))
   allocate (rscr9 (3,nwa))
   allocate (rscr10(3,nwa))
   allocate (rscr11(7,nwa))
   allocate (rscr12(7,nwa))
   allocate (rscr13(7,nwa))
   allocate (rscr14(7,nwa))
   allocate (rscr15(7,nwa))
   allocate (rscr16(7,nwa))

   ndims = 2
   idims(1) = mloops_w
   idims(2) = nwa

   call shdf5_irec(ndims,idims,'itab_w%loop',lvara=lscr)

   idims(1) = 3

   call shdf5_irec(ndims,idims,'itab_w%diru' ,rvara=rscr1)
   call shdf5_irec(ndims,idims,'itab_w%vxu'  ,rvara=rscr6)
   call shdf5_irec(ndims,idims,'itab_w%vyu'  ,rvara=rscr7)
   call shdf5_irec(ndims,idims,'itab_w%vzu'  ,rvara=rscr8)
   call shdf5_irec(ndims,idims,'itab_w%vxu_w',rvara=rscr9)
   call shdf5_irec(ndims,idims,'itab_w%vyu_w',rvara=rscr10)

   idims(1) = 7

   call shdf5_irec(ndims,idims,'itab_w%im'  ,ivara=iscr1)
   call shdf5_irec(ndims,idims,'itab_w%iv'  ,ivara=iscr3)

   call shdf5_irec(ndims,idims,'itab_w%dirv',rvara=rscr2)
   call shdf5_irec(ndims,idims,'itab_w%fwv' ,rvara=rscr3)
   call shdf5_irec(ndims,idims,'itab_w%fww' ,rvara=rscr4)
   call shdf5_irec(ndims,idims,'itab_w%farm',rvara=rscr11)
   call shdf5_irec(ndims,idims,'itab_w%farv',rvara=rscr12)

   call shdf5_irec(ndims,idims,'itab_w%gxps1',rvara=rscr13)
   call shdf5_irec(ndims,idims,'itab_w%gyps1',rvara=rscr14)
   call shdf5_irec(ndims,idims,'itab_w%gxps2',rvara=rscr15)
   call shdf5_irec(ndims,idims,'itab_w%gyps2',rvara=rscr16)

   idims(1) = 9

   call shdf5_irec(ndims,idims,'itab_w%iu'  ,ivara=iscr2)
   call shdf5_irec(ndims,idims,'itab_w%iw'  ,ivara=iscr4)
   call shdf5_irec(ndims,idims,'itab_w%fwu' ,rvara=rscr5)

   do iw = 1,nwa
      itab_w(iw)%loop(1:mloops_w) = lscr(1:mloops_w,iw)

      itab_w(iw)%im(1:7) = iscr1(1:7,iw)
      itab_w(iw)%iu(1:9) = iscr2(1:9,iw)
      itab_w(iw)%iv(1:7) = iscr3(1:7,iw)
      itab_w(iw)%iw(1:9) = iscr4(1:9,iw)

      itab_w(iw)%diru  (1:3) = rscr1 (1:3,iw)
      itab_w(iw)%dirv  (1:7) = rscr2 (1:7,iw)
      itab_w(iw)%fwv   (1:7) = rscr3 (1:7,iw)
      itab_w(iw)%fww   (1:7) = rscr4 (1:7,iw)
      itab_w(iw)%fwu   (1:9) = rscr5 (1:9,iw)
      itab_w(iw)%vxu   (1:3) = rscr6 (1:3,iw)
      itab_w(iw)%vyu   (1:3) = rscr7 (1:3,iw)
      itab_w(iw)%vzu   (1:3) = rscr8 (1:3,iw)
      itab_w(iw)%vxu_w (1:3) = rscr9 (1:3,iw)
      itab_w(iw)%vyu_w (1:3) = rscr10(1:3,iw)
      itab_w(iw)%farm  (1:7) = rscr11(1:7,iw)
      itab_w(iw)%farv  (1:7) = rscr12(1:7,iw)
      itab_w(iw)%gxps1 (1:7) = rscr13(1:7,iw)
      itab_w(iw)%gyps1 (1:7) = rscr14(1:7,iw)
      itab_w(iw)%gxps2 (1:7) = rscr15(1:7,iw)
      itab_w(iw)%gyps2 (1:7) = rscr16(1:7,iw)
   enddo

   deallocate(lscr,iscr1,iscr2,iscr3,  &
      rscr1,rscr2,rscr3,rscr4,rscr5,rscr6,rscr7,rscr8,rscr9,rscr10, &
      rscr11,rscr12,rscr13,rscr14,rscr15,rscr16)

! Check whether LAND/SEA models are used

   if (isfcl == 1) then

! Read SEAFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      call shdf5_irec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
      call shdf5_irec(ndims, idims, 'NSFPATS' ,ivars=nsfpats)

      mseaflux = nseaflux
      msfpats = nsfpats

      allocate (seaflux(nseaflux))

      allocate (nsfpatm(nsfpats))

      allocate (xemsfpat(5,nsfpats))
      allocate (yemsfpat(5,nsfpats))
      allocate (zemsfpat(5,nsfpats))


      idims(1) = nseaflux

      call shdf5_irec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%isfglobe)
      call shdf5_irec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
      call shdf5_irec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
      call shdf5_irec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iws)
      call shdf5_irec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
      call shdf5_irec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
      call shdf5_irec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
      call shdf5_irec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
      call shdf5_irec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
      call shdf5_irec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
      call shdf5_irec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sea)

      if (nsfpats > 0) then

         idims(1) = nsfpats

         call shdf5_irec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

         ndims = 2
         idims(1) = 5
         idims(2) = nsfpats

         call shdf5_irec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
         call shdf5_irec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
         call shdf5_irec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

      endif

! Read LANDFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      call shdf5_irec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
      call shdf5_irec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

      mlandflux = nlandflux
      mlfpats = nlfpats

      allocate (landflux(nlandflux))

      allocate (nlfpatm(nlfpats))

      allocate (xemlfpat(5,nlfpats))
      allocate (yemlfpat(5,nlfpats))
      allocate (zemlfpat(5,nlfpats))

      idims(1) = nlandflux

      call shdf5_irec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ilfglobe)
      call shdf5_irec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
      call shdf5_irec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
      call shdf5_irec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwl)
      call shdf5_irec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
      call shdf5_irec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
      call shdf5_irec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
      call shdf5_irec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
      call shdf5_irec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
      call shdf5_irec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
      call shdf5_irec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_land)

      if (nlfpats > 0) then

         idims(1) = nlfpats

         call shdf5_irec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

         ndims = 2
         idims(1) = 5
         idims(2) = nlfpats

         call shdf5_irec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
         call shdf5_irec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
         call shdf5_irec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

      endif

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



