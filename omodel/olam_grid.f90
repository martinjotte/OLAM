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
SUBROUTINE gridinit()

USE misc_coms,   ONLY: io6, runtype, mdomain, ngrids, initial, nxp, nzp, &
                       nzaux, timmax8, alloc_misc, iparallel, meshtype, &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim

USE leaf_coms,   ONLY: nzg, nzs, isfcl, nwl
USE sea_coms,    ONLY: nws

USE mem_ijtabs,  ONLY: istp, mrls, fill_jtabs, itab_md, itab_ud, itab_wd, itab_w
USE oplot_coms,  ONLY: op
USE mem_grid,    ONLY: nza, nma, nua, nva, nwa, &
                       zm, zt, dzt, xem, yem, zem, xew, yew, zew, glatw, glonw, &
                       alloc_grid1, alloc_grid2
USE mem_sflux,   ONLY: nlandflux, nseaflux
USE mem_nudge,   ONLY: nudflag, nudnxp, nwnud, itab_wnud, alloc_nudge1, &
                       xewnud, yewnud, zewnud

IMPLICIT NONE

REAL, ALLOCATABLE :: quarter_kite(:,:)

INTEGER :: npoly
INTEGER :: j, jmaxneg, jminpos, k
INTEGER :: imd,imd1,imd2,iud,iwd,iwnud,iwnudn,iw

REAL :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
REAL :: b11,b21,b31,b12,b22,b32,b13,b23,b33
REAL :: dist_wnud,dist,dist_min,xi,yi
REAL :: xin(6),yin(6)

! Read LAND and SEA files

IF (runtype /= 'MAKESFC' .AND. isfcl == 1) THEN

   WRITE(io6,'(/,a)') 'gridinit calling landfile_read'
   CALL landfile_read()

   WRITE(io6,'(/,a)') 'gridinit calling seafile_read'
   CALL seafile_read()

ENDIF

! Generate OLAM grid structure for 'MAKESFC' or 'MAKEGRID' runtype

IF (runtype == 'MAKESFC' .OR. runtype == 'MAKEGRID') THEN

! Vertical grid coordinate setup

   IF (nzaux < 0) THEN
      WRITE(io6,'(/,a)') 'gridinit calling gridset'
      CALL gridset()
   ELSE
      WRITE(io6,'(/,a)') 'gridinit calling gridset2'
      CALL gridset2()
   ENDIF

   WRITE(io6,'(a,f8.1)') ' Model top height = ',zm(nza-1)

! Print out vertical grid structure 

   CALL gridset_print()

! Horizontal grid setup

   IF (mdomain == 0) THEN

! If doing 'MAKEGRID' run and using nudging on global domain, generate
! nudging grid here

      if (runtype == 'MAKEGRID' .and. nudflag > 0) then

         write(io6,'(/,a)') 'gridinit calling icosahedron for nudging grid'

         CALL icosahedron(nudnxp)  ! global spherical domain; calls 2 allocs

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

         DEALLOCATE (itab_md, itab_ud, itab_wd)
         DEALLOCATE (xem, yem, zem)

         WRITE(io6,'(/,a)') 'gridinit after nudging grid construction'
         WRITE(io6,'(a,i8)')    ' nwnud = ',nwnud

      ENDIF

! Now generate global atmospheric grid

      WRITE(io6,'(/,a)') 'gridinit calling icosahedron'
      CALL icosahedron(nxp)  ! global spherical domain; calls 2 allocs

   ELSEIF (mdomain == 2 .OR. mdomain == 3) THEN

      WRITE(io6,'(/,a)') 'gridinit calling cartesian_2d'
      CALL cartesian_2d()    ! 2D cartesian channel domain; calls 2 allocs

   ELSEIF (mdomain == 4) THEN

      WRITE(io6,'(/,a)') 'gridinit calling cartesian_3d'
      CALL cartesian_3d()    ! 3D cartesian channel domain; calls 2 allocs

   ENDIF

   WRITE(io6,'(/,a)') 'gridinit after icosahedron or cartesian'
   WRITE(io6,'(a,i8)')    ' nma = ',nma
   WRITE(io6,'(a,i8)')    ' nua = ',nua
   WRITE(io6,'(a,i8)')    ' nwa = ',nwa

   IF (ngrids > 1 .AND. mdomain /= 2 .AND. mdomain /= 3) THEN

      WRITE(io6,'(/,a,i5)') 'gridinit calling spawn_nest; ngrids = ',ngrids

      CALL spawn_nest()  ! calls 2 allocs

      WRITE(io6,'(/,a)') 'gridinit after spawn_nest'
      WRITE(io6,'(a,i8)')   ' nma = ',nma
      WRITE(io6,'(a,i8)')   ' nua = ',nua
      WRITE(io6,'(a,i8)')   ' nwa = ',nwa

   ENDIF

! 'MAKESFC' run uses OLAM triangular atmos grid configuration
! to generate land/sea cells

   IF (runtype == 'MAKESFC') THEN
      WRITE(io6,'(/,a)') 'gridinit calling makesfc'
      CALL makesfc()
      RETURN
   ENDIF

! CHECK WHETHER THIS IS TRIANGLE OR HEXAGON MESH

   IF (meshtype == 1) THEN
      CALL delaunay()

      WRITE(io6,'(/,a)') 'gridinit after delaunay'
      WRITE(io6,'(a,i8)')    ' nma = ',nma
      WRITE(io6,'(a,i8)')    ' nua = ',nua
      WRITE(io6,'(a,i8)')    ' nva = ',nva
      WRITE(io6,'(a,i8)')    ' nwa = ',nwa
   ELSE
      CALL voronoi()
      CALL pcvt()

      WRITE(io6,'(/,a)') 'gridinit after voronoi and pcvt'
      WRITE(io6,'(a,i8)')    ' nma = ',nma
      WRITE(io6,'(a,i8)')    ' nua = ',nua
      WRITE(io6,'(a,i8)')    ' nva = ',nva
      WRITE(io6,'(a,i8)')    ' nwa = ',nwa
   ENDIF

! Allocate remaining GRID FOOTPRINT arrays for full domain

   WRITE(io6,'(/,a)') 'gridinit calling alloc_grid1 for full domain'

   CALL alloc_grid1(meshtype)

! Fill remaining GRID FOOTPRINT geometry for full domain

   WRITE(io6,'(/,a)') 'gridinit calling grid_geometry'

   IF (meshtype == 1) THEN
      CALL grid_geometry_tri()
   ELSE
      ALLOCATE (quarter_kite(2,nva))
      CALL grid_geometry_hex(quarter_kite)
   ENDIF

! Initialize dtlm, dtsm, ndtrat, and nacoust,
! and compute the timestep schedule for all grid operations.

   WRITE(io6,'(/,a)') 'gridinit calling modsched'

   CALL modsched()

   WRITE(io6,'(/,a)') 'gridinit calling fill_jtabs'

   CALL fill_jtabs(nma,nua,nva,nwa)

! Allocate remaining unstructured grid geometry arrays

   WRITE(io6,'(/,a)') 'gridinit calling alloc_grid2'

   CALL alloc_grid2(meshtype)

! Set up control volumes

   WRITE(io6,'(/,a)') 'gridinit calling ctrlvols'

   IF (meshtype == 1) THEN
      CALL ctrlvols_tri()
   ELSE
      CALL ctrlvols_hex(quarter_kite)
      DEALLOCATE (quarter_kite)
   ENDIF

! If doing 'MAKEGRID' run and using nudging on global grid, compute
! nudging indices and coefficients of W points

   IF (runtype == 'MAKEGRID' .AND. mdomain == 0 .AND. nudflag > 0) THEN

! Compute a mean distance between adjacent nudging points.  It is safe to assume
! that any ATM grid IW point will be closer than this distance to at least one
! nuding point

      dist_wnud = 7100.e3 / REAL(nudnxp)

! Loop over all W points

      DO iw = 2,nwa

! Initialize minimum distance variable to large number

         dist_min = 1.e12

! Loop over all nudging points

         DO iwnud = 2,nwnud

! Skip interaction if ATM grid point and nudging point are more than
! dist_wnud apart in any earth-frame direction

            IF (ABS(xew(iw) - xewnud(iwnud)) > dist_wnud) CYCLE
            IF (ABS(yew(iw) - yewnud(iwnud)) > dist_wnud) CYCLE
            IF (ABs(zew(iw) - zewnud(iwnud)) > dist_wnud) CYCLE

! Compute distance between current ATM IW grid point and nudging point,
! and reset minimum distance and nudging point index if new minimum is found

            dist = sqrt((xew(iw) - xewnud(iwnud))**2 &
                      + (yew(iw) - yewnud(iwnud))**2 &
                      + (zew(iw) - zewnud(iwnud))**2)

            IF (dist_min > dist) THEN
               dist_min = dist
               itab_w(iw)%iwnud(1) = iwnud
            ENDIF

         ENDDO

! Now that primary nudging point has been found for current IW point,
! compute its x,y components on a polar stereographic plane tangent
! at W point (W point is at 0,0)

         iwnud = itab_w(iw)%iwnud(1)

         CALL e_ps(xewnud(iwnud),yewnud(iwnud),zewnud(iwnud), &
                   glatw(iw),glonw(iw),xi,yi)

! Initialize vecprodz_minpos and vecprodz_maxneg

         vecprodz_maxneg = -1.e15
         vecprodz_minpos =  1.e15

! Loop through nearest polygon neighbors (j, iwnudn) of nudging point IWNUD

         npoly = itab_wnud(iwnud)%npoly

         DO j = 1,npoly

! Get nudging point index (iwnudn) for current polygon neighbor of iwnud.

            iwnudn = itab_wnud(iwnud)%iwnud(j)

! Compute x,y components of iwnudn polygon center on a polar stereographic 
! plane tangent at W point

!write(6,'(a,3i7,8e13.3)') 'og1 ',iw,j,iwnudn, &
!          xewnud(iwnudn),yewnud(iwnudn),zewnud(iwnudn), &
!          glatw(iw),glonw(iw)

            CALL e_ps(xewnud(iwnudn),yewnud(iwnudn),zewnud(iwnudn), &
                      glatw(iw),glonw(iw),xin(j),yin(j))

! Compute z component (in polar stereographic space) of vector product of
! the vector from iwnud to iw (at 0,0) and the vector from iwnud to iwnudn.

            vecprodz = -xi * (yin(j) - yi) + yi * (xin(j) - xi)

! Compute scalar product of the vector from iwnud to iw (at 0,0) and the
! vector from iwnud to iwnudn in polar stereographic space.

            scalprod = -xi * (xin(j) - xi) - yi * (yin(j) - yi)

! Check whether scalar product is positive for current J point.  If so,
! J point is a candidate for the nudging triad for IW.

            IF (scalprod > 0.) THEN

! Identify maximum negative vecprodz among all iwnudn polygon neighbors of
! iwnud.  This iwnudn will be second point of nudging triad for IW

               IF (vecprodz < 0. .AND. vecprodz > vecprodz_maxneg) THEN
                  vecprodz_maxneg = vecprodz
                  jmaxneg = j
                  itab_w(iw)%iwnud(2) = iwnudn
               ENDIF

! Identify minimum positive vecprodz among all iwnudn polygon neighbors of
! iwnud.  This iwnudn will be third point of nudging triad for IW

               IF (vecprodz >= 0. .AND. vecprodz < vecprodz_minpos) THEN
                  vecprodz_minpos = vecprodz
                  jminpos = j
                  itab_w(iw)%iwnud(3) = iwnudn
               ENDIF
            ENDIF

         ENDDO

! Lastly, fill 3 nudging weight coefficients for this W point.
! Weights are computed in 2_d polar stereographic space.

! Invert matrix of coordinates

         CALL matinv3x3(1.,xi,yi, &
                        1.,xin(jmaxneg),yin(jmaxneg), &
                        1.,xin(jminpos),yin(jminpos), &
                        b11,b21,b31,b12,b22,b32,b13,b23,b33)

! Assign coefficients

         itab_w(iw)%fnud(1) = b11
         itab_w(iw)%fnud(2) = b21
         itab_w(iw)%fnud(3) = b31

      ENDDO

   ENDIF

! Write GRIDFILE

   WRITE(io6,'(/,a)') 'gridinit calling gridfile_write'
   CALL gridfile_write()

ELSE

! Read atmos grid for INITIAL/HISTORY/PLOTONLY/PARCOMBINE run

   WRITE(io6,'(/,a)') 'gridinit calling gridfile_read'
   CALL gridfile_read()
   
   CALL gridset_print()

ENDIF

! Print out land, sea, landflux, and seaflux grid sizes

IF (isfcl == 1) THEN
   WRITE(io6,'(/,a)') 'gridinit before return to olam_run'
   WRITE(io6,'(a,i8)')   ' nwl       = ',nwl
   WRITE(io6,'(a,i8)')   ' nws       = ',nws
   WRITE(io6,'(a,i8)')   ' nlandflux = ',nlandflux
   WRITE(io6,'(a,i8)')   ' nseaflux  = ',nseaflux
ENDIF

RETURN
END SUBROUTINE gridinit

!===============================================================================

SUBROUTINE gridset()

USE mem_grid,    ONLY: nza, mza, &
                       zm, zt, dzm, dzt, dzim, dzit, &
                       zfacm, zfact, zfacim, zfacit, &
                       alloc_gridz
USE misc_coms,   ONLY: io6, nzp, deltaz, dzrat, dzmax, ztop, zbase, mdomain
USE consts_coms, ONLY: erad
USE oname_coms,  ONLY: nl

IMPLICIT NONE

INTEGER :: k
REAL :: dzr,dzrati
REAL, ALLOCATABLE :: zmvec(:),ztvec(:)

nza = nzp
mza = nzp

CALL alloc_gridz()
ALLOCATE (zmvec(-1:mza+1),ztvec(-1:mza+1))

! calculate zm

IF ( deltaz < SPACING(0.) ) THEN
   zmvec(1:nzp) = nl%zz(1:nzp)
   zmvec(nzp+1) = 2. * zmvec(nzp) - zmvec(nzp-1)
ELSE
   zmvec(1) = zbase
   zmvec(2) = zbase + deltaz
   DO k = 3,nzp+1
      zmvec(k) = zmvec(k-1) &
               + MIN(dzrat * (zmvec(k-1) - zmvec(k-2)),MAX(deltaz,dzmax))
   ENDDO
ENDIF
dzrati = (zmvec(2) - zmvec(1)) / (zmvec(3) - zmvec(2))
zmvec(0) = zmvec(1) - (zmvec(2) - zmvec(1)) * dzrati
zmvec(-1) = zmvec(0) - (zmvec(1) - zmvec(0)) * dzrati

! compute zt values by geometric interpolation.

DO k = 1,nza
   dzr = SQRT(SQRT((zmvec(k+1) - zmvec(k)) / (zmvec(k-1) - zmvec(k-2))))
   ztvec(k) = zmvec(k-1) + (zmvec(k) - zmvec(k-1)) / (1. + dzr)
ENDDO
ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

! Other vertical coordinate values

DO k = 1,nza
   zm(k) = zmvec(k)
   zt(k) = ztvec(k)

   dzm(k) = ztvec(k+1) - ztvec(k)
   dzt(k) = zmvec(k) - zmvec(k-1)

   dzim(k) = 1. / dzm(k)
   dzit(k) = 1. / dzt(k)

   IF (mdomain < 2) THEN
      zfacm(k) = (erad + zm(k)) / erad
      zfact(k) = (erad + zt(k)) / erad
   ELSE
      zfacm(k) = 1.
      zfact(k) = 1.
   ENDIF

   zfacim(k) = 1. / zfacm(k)
   zfacit(k) = 1. / zfact(k)
ENDDO

DEALLOCATE (zmvec,ztvec)

RETURN
END SUBROUTINE gridset

!===============================================================================

SUBROUTINE gridset2()

USE mem_grid,    ONLY: nza, mza, &
                       zm, zt, dzm, dzt, dzim, dzit, &
                       zfacm, zfact, zfacim, zfacit, &
                       alloc_gridz
USE misc_coms,   ONLY: io6, nzp, nzaux, mdomain, &
                       zbase, dzbase, ztop, dztop, zaux, dzaux
USE consts_coms, ONLY: erad

IMPLICIT NONE

INTEGER :: izaux, kvec, nseries, iseries, k

REAL :: zend, dzend, dzbeg, ztarg, dzr, rdzr

REAL, ALLOCATABLE :: zmvec(:),ztvec(:)

! This subroutine will compute number of vertical levels.
! First, allocate zmvec and ztvec arrays to large size.

ALLOCATE (zmvec(0:500),ztvec(500))

! Fill starting values of zmvec and ztvec

zmvec(0) = zbase - dzbase
zmvec(1) = zbase
zmvec(2) = zbase + dzbase

ztvec(1) = zbase - dzbase * .5
ztvec(2) = zbase + dzbase * .5

kvec = 2

! Loop through any (zaux,dzaux) values, plus one extra cycle for (ztop,dztop)

DO izaux = 1,nzaux+1

   IF (izaux < nzaux + 1) THEN
      zend = zaux(izaux)
      dzend = dzaux(izaux)
   ELSE
      zend = ztop
      dzend = dztop
   ENDIF

   dzbeg = zmvec(kvec) - zmvec(kvec-1)

! Target height for next series

   ztarg = zend + .5 * dzend

! Stretch ratio based on geometric series sum

   dzr = (ztarg - zmvec(kvec)) / (ztarg - zmvec(kvec-1) - dzend)

! NSERIES from stretch ratio

   IF (abs(dzr - 1.) > 1.e-3) THEN
      nseries = max(1,nint(log(dzend/dzbeg) / log(dzr)))
   ELSE
      nseries = max(1,nint((ztarg - zmvec(kvec)) / dzbeg))
   ENDIF

! Stretch ratio based on N

   dzr = (dzend / dzbeg) ** (1./REAL(nseries))
   rdzr = sqrt(dzr)

! Loop over series and compute new levels

   DO iseries = 1,nseries

      if (kvec >= nzp - 1) GO TO 5

      kvec = kvec + 1
      dzbeg = dzbeg * dzr

      zmvec(kvec) = zmvec(kvec-1) + dzbeg
      ztvec(kvec) = zmvec(kvec-1) + (zmvec(kvec) - zmvec(kvec-1)) / (1. + rdzr)

   ENDDO

ENDDO

5 CONTINUE

! Fill NZA and MZA values

nza = kvec + 1
mza = nza

IF (nza < 3) THEN
   WRITE(io6,'(a)')      'OLAM requires that NZA >= 3.'
   WRITE(io6,'(a,i5,a)') 'NZA = ',nza,' in subroutine gridset2.  Stopping model.'
   STOP 'stop nza in gridset2'
ENDIF

! Fill top 2 ZMVEC and ZTVEC values

zmvec(nza)   = zmvec(nza-1) * 2. - zmvec(nza-2)
zmvec(nza+1) = zmvec(nza)   * 2. - zmvec(nza-1)

ztvec(nza)   = .5 * (zmvec(nza-1) + zmvec(nza))
ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

! Allocate main grid arrays

CALL alloc_gridz()

! Other vertical coordinate values

DO k = 1,nza
   zm(k) = zmvec(k)
   zt(k) = ztvec(k)

   dzm(k) = ztvec(k+1) - ztvec(k)
   dzt(k) = zmvec(k) - zmvec(k-1)

   dzim(k) = 1. / dzm(k)
   dzit(k) = 1. / dzt(k)

   IF (mdomain < 2) THEN
      zfacm(k) = (erad + zm(k)) / erad
      zfact(k) = (erad + zt(k)) / erad
   ELSE
      zfacm(k) = 1.
      zfact(k) = 1.
   ENDIF

   zfacim(k) = 1. / zfacm(k)
   zfacit(k) = 1. / zfact(k)
ENDDO

DEALLOCATE (zmvec,ztvec)

RETURN
END SUBROUTINE gridset2

!===============================================================================

SUBROUTINE gridset_print()

USE misc_coms, ONLY: io6
USE mem_grid,  ONLY: nza, zm, zt, dzt

implicit none

INTEGER :: k

! Print vertical grid structure 

WRITE(io6,'(/,a)') '================================================='
WRITE(io6,'(a)'  ) '         OLAM VERTICAL GRID STRUCTURE'
WRITE(io6,'(a)'  ) '================================================='
WRITE(io6,'(a)'  ) '   k    zm(m)     k     zt(m)    dzm(m)    dzrat'
WRITE(io6,'(a,/)') '================================================='

DO k = nza-1,1,-1
   IF (k == nza-1 .or. k == 1) THEN
      WRITE(io6,11) k, zm(k), dzt(k+1)/dzt(k)
   ELSE
      WRITE(io6,12) k, zm(k), dzt(k+1)/dzt(k)
   ENDIF
   IF (k > 1) WRITE(io6,13) k,zt(k),dzt(k)
ENDDO
   
11 format (i4,f10.2,1x,3('========='),f6.3)
12 format (i4,f10.2,1x,3('---------'),f6.3)
13 format (15x,i4,2f10.2)

END

!===============================================================================

SUBROUTINE topo_init(nqa,topq,glatq,glonq,xeq,yeq,zeq)

USE misc_coms,   ONLY: io6, deltax
USE consts_coms, ONLY: pi1, pio180

IMPLICIT NONE

INTEGER, INTENT(in) :: nqa
REAL, INTENT(in) :: glatq(nqa),glonq(nqa),xeq(nqa),yeq(nqa),zeq(nqa)

REAL, INTENT(out) :: topq(nqa)

INTEGER :: iq

REAL :: hfwid
REAL :: hgt
REAL :: hfwid2

REAL :: r, r0

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

DO iq = 2,nqa
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

ENDDO

RETURN
END SUBROUTINE topo_init

!===============================================================================

SUBROUTINE gridfile_write()

  USE max_dims,   ONLY: maxngrdll
  USE misc_coms,  ONLY: io6, ngrids, gridfile, mdomain, meshtype, nzp, nxp, &
       iclobber, itopoflg, &
       deltax, deltaz, dzmax, dzrat, zbase, dzbase, ztop, dztop, &
       nzaux, zaux, dzaux, ngrdll, grdrad, grdlat, grdlon, meshtype
  USE mem_ijtabs, ONLY: mloops_m, mloops_u, mloops_v, mloops_w, mrls, &
       itab_m, itab_u, itab_v, itab_w
  USE mem_grid,   ONLY: nza, nma, nua, nva, nwa, nsw_max, &
       zm, zt, dzm, dzt, dzim, dzit, &
       zfacm, zfact, zfacim, zfacit, &
       lpm, lpu, lcu, lpv, lcv, lpw, lsw, &
       topm, topw, xem, yem, zem, xeu, yeu, zeu, &
       xev, yev, zev, xew, yew, zew, &
       unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
       dnu, dniu, dnv, dniv, arw0, arm0, &
       glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
       aru, arv, volui, volvi, arw, volwi, volt, volti
  USE leaf_coms,  ONLY: isfcl
  USE mem_sflux,  ONLY: nseaflux, nlandflux, seaflux, landflux, &
       nsfpats, nlfpats, nsfpatm, nlfpatm, &
       xemsfpat, yemsfpat, zemsfpat, &
       xemlfpat, yemlfpat, zemlfpat

  USE hdf5_utils, ONLY: shdf5_orec, shdf5_open, shdf5_close

  USE mem_nudge,  ONLY: nudflag, nudnxp, nwnud, itab_wnud, &
                        xewnud, yewnud, zewnud

 IMPLICIT NONE

  ! This routine writes the grid variables to the grid file.

  INTEGER :: im, iu, iv, iw, iwnud

  INTEGER :: ndims, idims(2)

  ! Scratch arrays for copying output

  LOGICAL, ALLOCATABLE :: lscr(:,:)
  INTEGER, ALLOCATABLE :: iscr(:,:)

  REAL, ALLOCATABLE :: rscr(:,:)

  WRITE(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(io6,*) 'grid_write: opening file:', TRIM(gridfile)
  WRITE(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  CALL shdf5_open(gridfile,'W',iclobber)

  ! Write the gridfile information that exists in namelist

  ndims = 1
  idims(1) = 1
  idims(2) = 1

  CALL shdf5_orec(ndims, idims, 'NZP'     , ivars=nzp)
  CALL shdf5_orec(ndims, idims, 'NXP'     , ivars=nxp)
  CALL shdf5_orec(ndims, idims, 'MDOMAIN' , ivars=mdomain)
  CALL shdf5_orec(ndims, idims, 'MESHTYPE', ivars=meshtype)
  CALL shdf5_orec(ndims, idims, 'NGRIDS'  , ivars=ngrids)
  CALL shdf5_orec(ndims, idims, 'ISFCL'   , ivars=isfcl)
  CALL shdf5_orec(ndims, idims, 'ITOPOFLG', ivars=itopoflg)
  CALL shdf5_orec(ndims, idims, 'DELTAX'  , rvars=deltax)
  CALL shdf5_orec(ndims, idims, 'DELTAZ'  , rvars=deltaz)
  CALL shdf5_orec(ndims, idims, 'DZRAT'   , rvars=dzrat)
  CALL shdf5_orec(ndims, idims, 'DZMAX'   , rvars=dzmax)
  CALL shdf5_orec(ndims, idims, 'ZBASE'   , rvars=zbase)
  CALL shdf5_orec(ndims, idims, 'DZBASE'  , rvars=dzbase)
  CALL shdf5_orec(ndims, idims, 'ZTOP'    , rvars=ztop)
  CALL shdf5_orec(ndims, idims, 'DZTOP'   , rvars=dztop)

  CALL shdf5_orec(ndims, idims, 'NZAUX'   , ivars=nzaux)

  IF (nzaux > 0) THEN
     idims(1) = nzaux

     CALL shdf5_orec(ndims, idims, 'ZAUX'    , rvara=zaux)
     CALL shdf5_orec(ndims, idims, 'DZAUX'   , rvara=dzaux)
  ENDIF

  idims(1) = ngrids

  CALL shdf5_orec(ndims, idims, 'NGRDLL' , ivara=ngrdll)
  CALL shdf5_orec(ndims, idims, 'GRDRAD' , rvara=grdrad)

  ndims = 2
  idims(1) = ngrids
  idims(2) = maxngrdll

  CALL shdf5_orec(ndims, idims, 'GRDLAT', rvara=grdlat(1:ngrids,:))
  CALL shdf5_orec(ndims, idims, 'GRDLON', rvara=grdlon(1:ngrids,:))

  ! Write the grid dimensions

  ndims = 1
  idims(1) = 1

  CALL shdf5_orec(ndims, idims, 'NZA'    , ivars=nza)
  CALL shdf5_orec(ndims, idims, 'NMA'    , ivars=nma)
  CALL shdf5_orec(ndims, idims, 'NUA'    , ivars=nua)
  CALL shdf5_orec(ndims, idims, 'NVA'    , ivars=nva)
  CALL shdf5_orec(ndims, idims, 'NWA'    , ivars=nwa)
  CALL shdf5_orec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
  CALL shdf5_orec(ndims, idims, 'MRLS'   , ivars=mrls)

  ! Write grid structure variables

  idims(1) = nza

  CALL shdf5_orec(ndims, idims, 'ZM'    , rvara=zm)
  CALL shdf5_orec(ndims, idims, 'ZT'    , rvara=zt)
  CALL shdf5_orec(ndims, idims, 'DZM'   , rvara=dzm)
  CALL shdf5_orec(ndims, idims, 'DZT'   , rvara=dzt)
  CALL shdf5_orec(ndims, idims, 'DZIM'  , rvara=dzim)
  CALL shdf5_orec(ndims, idims, 'DZIT'  , rvara=dzit)
  CALL shdf5_orec(ndims, idims, 'ZFACM' , rvara=zfacm)
  CALL shdf5_orec(ndims, idims, 'ZFACT' , rvara=zfact)
  CALL shdf5_orec(ndims, idims, 'ZFACIM', rvara=zfacim)
  CALL shdf5_orec(ndims, idims, 'ZFACIT', rvara=zfacit)

  idims(1) = nma

  CALL shdf5_orec(ndims, idims, 'LPM'  , ivara=lpm)
  CALL shdf5_orec(ndims, idims, 'TOPM' , rvara=topm)
  CALL shdf5_orec(ndims, idims, 'XEM'  , rvara=xem)
  CALL shdf5_orec(ndims, idims, 'YEM'  , rvara=yem)
  CALL shdf5_orec(ndims, idims, 'ZEM'  , rvara=zem)
  CALL shdf5_orec(ndims, idims, 'ARM0' , rvara=arm0)
  CALL shdf5_orec(ndims, idims, 'GLATM', rvara=glatm)
  CALL shdf5_orec(ndims, idims, 'GLONM', rvara=glonm)

  IF (meshtype == 1) THEN

     idims(1) = nua

     CALL shdf5_orec(ndims, idims, 'LPU'  , ivara=lpu)
     CALL shdf5_orec(ndims, idims, 'LCU'  , ivara=lcu)
     CALL shdf5_orec(ndims, idims, 'XEU'  , rvara=xeu)
     CALL shdf5_orec(ndims, idims, 'YEU'  , rvara=yeu)
     CALL shdf5_orec(ndims, idims, 'ZEU'  , rvara=zeu)
     CALL shdf5_orec(ndims, idims, 'GLATU', rvara=glatu)
     CALL shdf5_orec(ndims, idims, 'GLONU', rvara=glonu)

  ELSE

     idims(1) = nva

     CALL shdf5_orec(ndims, idims, 'LPV'  , ivara=lpv)
     CALL shdf5_orec(ndims, idims, 'LCV'  , ivara=lcv)
     CALL shdf5_orec(ndims, idims, 'XEV'  , rvara=xev)
     CALL shdf5_orec(ndims, idims, 'YEV'  , rvara=yev)
     CALL shdf5_orec(ndims, idims, 'ZEV'  , rvara=zev)
     CALL shdf5_orec(ndims, idims, 'GLATV', rvara=glatv)
     CALL shdf5_orec(ndims, idims, 'GLONV', rvara=glonv)

  ENDIF

  CALL shdf5_orec(ndims, idims, 'DNU' , rvara=dnu)
  CALL shdf5_orec(ndims, idims, 'DNV' , rvara=dnv)
  CALL shdf5_orec(ndims, idims, 'DNIU', rvara=dniu)
  CALL shdf5_orec(ndims, idims, 'DNIV', rvara=dniv)
  CALL shdf5_orec(ndims, idims, 'UNX' , rvara=unx)
  CALL shdf5_orec(ndims, idims, 'UNY' , rvara=uny)
  CALL shdf5_orec(ndims, idims, 'UNZ' , rvara=unz)
  CALL shdf5_orec(ndims, idims, 'VNX' , rvara=vnx)
  CALL shdf5_orec(ndims, idims, 'VNY' , rvara=vny)
  CALL shdf5_orec(ndims, idims, 'VNZ' , rvara=vnz)

  idims(1) = nwa

  CALL shdf5_orec(ndims, idims, 'LPW'  , ivara=lpw)
  CALL shdf5_orec(ndims, idims, 'LSW'  , ivara=lsw)
  CALL shdf5_orec(ndims, idims, 'XEW'  , rvara=xew)
  CALL shdf5_orec(ndims, idims, 'YEW'  , rvara=yew)
  CALL shdf5_orec(ndims, idims, 'ZEW'  , rvara=zew)
  CALL shdf5_orec(ndims, idims, 'TOPW' , rvara=topw)
  CALL shdf5_orec(ndims, idims, 'ARW0' , rvara=arw0)
  CALL shdf5_orec(ndims, idims, 'GLATW', rvara=glatw)
  CALL shdf5_orec(ndims, idims, 'GLONW', rvara=glonw)
  CALL shdf5_orec(ndims, idims, 'WNX'  , rvara=wnx)
  CALL shdf5_orec(ndims, idims, 'WNY'  , rvara=wny)
  CALL shdf5_orec(ndims, idims, 'WNZ'  , rvara=wnz)

  ndims = 2
  idims(1) = nza

  IF (meshtype == 1) THEN

     idims(2) = nua

     CALL shdf5_orec(ndims, idims, 'VOLUI', rvara=volui)

  ELSE

     idims(2) = nva

     CALL shdf5_orec(ndims, idims, 'ARV'  , rvara=arv)
     CALL shdf5_orec(ndims, idims, 'VOLVI', rvara=volvi)

  ENDIF

  CALL shdf5_orec(ndims, idims, 'ARU'  , rvara=aru)

  idims(2) = nwa

  CALL shdf5_orec(ndims, idims, 'ARW'  , rvara=arw)
  CALL shdf5_orec(ndims, idims, 'VOLWI', rvara=volwi)
  CALL shdf5_orec(ndims, idims, 'VOLT' , dvara=volt)
  CALL shdf5_orec(ndims, idims, 'VOLTI', dvara=volti)

  ! Write ITAB_M SCALARS

  ndims = 1
  idims(1) = nma
  idims(2) = 1

  CALL shdf5_orec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
  CALL shdf5_orec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm)
  CALL shdf5_orec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe)
  CALL shdf5_orec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
  CALL shdf5_orec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
  CALL shdf5_orec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
  CALL shdf5_orec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh)

  ! Write ITAB_M ARRAYS

  ndims = 2
  idims(1) = mloops_m
  idims(2) = nma

  ALLOCATE (lscr(mloops_m,nma)) ; lscr = .FALSE.
  DO im = 1,nma
     lscr(1:mloops_m,im) = itab_m(im)%loop(1:mloops_m)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_m%loop',lvara=lscr)
  DEALLOCATE (lscr)

  idims(1) = 7

  ALLOCATE (iscr(7,nma)) ; iscr = 0
  IF (meshtype == 1) THEN
     DO im = 1,nma
        iscr(1:7,im) = itab_m(im)%iu(1:7)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_m%iu',ivara=iscr)
  ELSE
     DO im = 1,nma
        iscr(1:7,im) = itab_m(im)%iv(1:7)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_m%iv',ivara=iscr)
  ENDIF
  DEALLOCATE (iscr)

  ALLOCATE (iscr(7,nma)) ; iscr = 0
  DO im = 1,nma
     iscr(1:7,im) = itab_m(im)%iw(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_m%iw',ivara=iscr)
  DEALLOCATE (iscr)

  ALLOCATE (rscr(7,nma)) ; rscr = 0.0
  DO im = 1,nma
     rscr(1:7,im) = itab_m(im)%fmw(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_m%fmw',rvara=rscr)
  DEALLOCATE (rscr)


  IF (meshtype == 1) THEN

     ! Write ITAB_U SCALARS

     ndims = 1
     idims(1) = nua
     idims(2) = 1

     CALL shdf5_orec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)
     CALL shdf5_orec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)
     CALL shdf5_orec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)

     CALL shdf5_orec(ndims,idims,'itab_u%gcf36'  ,rvara=itab_u(:)%gcf36)
     CALL shdf5_orec(ndims,idims,'itab_u%gcf45'  ,rvara=itab_u(:)%gcf45)

     CALL shdf5_orec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
     CALL shdf5_orec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)
     CALL shdf5_orec(ndims,idims,'itab_u%crossmm',rvara=itab_u(:)%crossmm)
     CALL shdf5_orec(ndims,idims,'itab_u%crossww',rvara=itab_u(:)%crossww)

     ! Write ITAB_U ARRAYS

     ndims = 2
     idims(1) = mloops_u
     idims(2) = nua

     ALLOCATE (lscr(mloops_u,nua))
     DO iu = 1,nua
        lscr(1:mloops_u,iu) = itab_u(iu)%loop(1:mloops_u)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%loop',lvara=lscr)
     DEALLOCATE (lscr)

     idims(1) = 2

     ALLOCATE (iscr(2,nua))
     DO iu = 1,nua
        iscr(1:2,iu) = itab_u(iu)%im(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%im',ivara=iscr)
     DEALLOCATE (iscr)

     ALLOCATE (rscr(2,nua))
     DO iu = 1,nua
        rscr(1:2,iu) = itab_u(iu)%vxw_u(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%vxw_u',rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (rscr(2,nua))
     DO iu = 1,nua
        rscr(1:2,iu) = itab_u(iu)%vyw_u(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%vyw_u',rvara=rscr)
     DEALLOCATE (rscr)

     idims(1) = 4

     ALLOCATE (rscr(4,nua))
     DO iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%diru(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%diru' ,rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (rscr(4,nua))
     DO iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%tuu(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%tuu',rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (rscr(4,nua))
     DO iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%guw(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%guw',rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (rscr(4,nua))
     DO iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%vxu_u(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%vxu_u',rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (rscr(4,nua))
     DO iu = 1,nua
        rscr(1:4,iu) = itab_u(iu)%vyu_u(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%vyu_u',rvara=rscr)
     DEALLOCATE (rscr)

     idims(1) = 6

     ALLOCATE (rscr(6,nua))
     DO iu = 1,nua
        rscr(1:6,iu) = itab_u(iu)%fuw(1:6)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%fuw',rvara=rscr)
     DEALLOCATE (rscr)

     ALLOCATE (iscr(6,nua))
     DO iu = 1,nua
        iscr(1:6,iu) = itab_u(iu)%iw(1:6)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%iw' ,ivara=iscr)
     DEALLOCATE (iscr)

     idims(1) = 12

     ALLOCATE (iscr(12,nua))
     DO iu = 1,nua
        iscr(1:12,iu) = itab_u(iu)%iu(1:12)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%iu',ivara=iscr)
     DEALLOCATE (iscr)

     ALLOCATE (rscr(12,nua))
     DO iu = 1,nua
        rscr(1:12,iu) = itab_u(iu)%fuu(1:12)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_u%fuu',rvara=rscr)
     DEALLOCATE (rscr)

  ENDIF

  IF (meshtype == 2) THEN

     ! Write ITAB_V SCALARS

     ndims = 1
     idims(1) = nva
     idims(2) = 1

     CALL shdf5_orec(ndims,idims,'itab_v%ivp'    ,ivara=itab_v(:)%ivp)
     CALL shdf5_orec(ndims,idims,'itab_v%mrlv'   ,ivara=itab_v(:)%mrlv)
     CALL shdf5_orec(ndims,idims,'itab_v%ivglobe',ivara=itab_v(:)%ivglobe)

     ! Write ITAB_V ARRAYS

     ndims = 2
     idims(1) = mloops_v
     idims(2) = nva

     ALLOCATE (lscr(mloops_v,nva))
     DO iv = 1,nva
        lscr(1:mloops_v,iv) = itab_v(iv)%loop(1:mloops_v)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%loop',lvara=lscr)
     DEALLOCATE(lscr)

     idims(1) = 2

     ALLOCATE (rscr(2,nva))
     DO iv = 1,nva
        rscr(1:2,iv) = itab_v(iv)%farw(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%farw',rvara=rscr)
     DEALLOCATE(rscr)

     ALLOCATE (rscr(2,nva))
     DO iv = 1,nva
        rscr(1:2,iv) = itab_v(iv)%cosv(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%cosv',rvara=rscr)
     DEALLOCATE(rscr)

     ALLOCATE (rscr(2,nva))
     DO iv = 1,nva
        rscr(1:2,iv) = itab_v(iv)%sinv(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%sinv',rvara=rscr)
     DEALLOCATE(rscr)

     ALLOCATE (rscr(2,nva))
     DO iv = 1,nva
        rscr(1:2,iv) = itab_v(iv)%dxps(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%dxps',rvara=rscr)
     DEALLOCATE(rscr)

     ALLOCATE (rscr(2,nva))
     DO iv = 1,nva
        rscr(1:2,iv) = itab_v(iv)%dyps(1:2)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%dyps',rvara=rscr)
     DEALLOCATE(rscr)

     idims(1) = 4

     ALLOCATE (iscr(4,nva))
     DO iv = 1,nva
        iscr(1:4,iv) = itab_v(iv)%iw(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%iw',ivara=iscr)
     DEALLOCATE(iscr)

     ALLOCATE (rscr(4,nva))
     DO iv = 1,nva
        rscr(1:4,iv) = itab_v(iv)%fvw(1:4)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%fvw',rvara=rscr)
     DEALLOCATE(rscr)

     idims(1) = 6

     ALLOCATE (iscr(6,nva))
     DO iv = 1,nva
        iscr(1:6,iv) = itab_v(iv)%im(1:6)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%im',ivara=iscr)
     DEALLOCATE(iscr)

     idims(1) = 12

     ALLOCATE (rscr(12,nva))
     DO iv = 1,nva
        rscr(1:12,iv) = itab_v(iv)%fvv(1:12)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%fvv',rvara=rscr)
     DEALLOCATE(rscr)

     idims(1) = 16

     ALLOCATE (iscr(16,nva))
     DO iv = 1,nva
        iscr(1:16,iv) = itab_v(iv)%iv(1:16)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%iv',ivara=iscr)
     DEALLOCATE(iscr)

     ALLOCATE (rscr(16,nva))
     DO iv = 1,nva
        rscr(1:16,iv) = itab_v(iv)%fuv(1:16)
     ENDDO
     CALL shdf5_orec(ndims,idims,'itab_v%fuv',rvara=rscr)
     DEALLOCATE(rscr)

  ENDIF

  ! Write ITAB_W SCALARS

  ndims = 1
  idims(1) = nwa
  idims(2) = 1

  CALL shdf5_orec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly)
  CALL shdf5_orec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)
  CALL shdf5_orec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
  CALL shdf5_orec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
  CALL shdf5_orec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
  CALL shdf5_orec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
  CALL shdf5_orec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

  CALL shdf5_orec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)
  CALL shdf5_orec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)
  CALL shdf5_orec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

  CALL shdf5_orec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w)
  CALL shdf5_orec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w)

  CALL shdf5_orec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w)
  CALL shdf5_orec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w)
  CALL shdf5_orec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w)

  ! Write ITAB_W ARRAYS

  ndims = 2
  idims(1) = mloops_w
  idims(2) = nwa

  ALLOCATE (lscr(mloops_w,nwa))
  DO iw = 1,nwa
     lscr(1:mloops_w,iw) = itab_w(iw)%loop(1:mloops_w)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%loop',lvara=lscr)
  DEALLOCATE(lscr)

  idims(1) = 3

  ALLOCATE (iscr(3,nwa))
  DO iw = 1,nwa
     iscr(1:3,iw) = itab_w(iw)%iwnud(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%iwnud',ivara=iscr)
  DEALLOCATE(iscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%fnud(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%fnud',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr (3,nwa))
  DO iw = 1,nwa
     rscr (1:3,iw) = itab_w(iw)%diru (1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%diru',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vxu(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%vxu',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vyu(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%vyu',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vzu(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%vzu',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vxu_w(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%vxu_w',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(3,nwa))
  DO iw = 1,nwa
     rscr(1:3,iw) = itab_w(iw)%vyu_w(1:3)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%vyu_w',rvara=rscr)
  DEALLOCATE(rscr)

  idims(1) = 7

  ALLOCATE (iscr(7,nwa))
  DO iw = 1,nwa
     iscr(1:7,iw) = itab_w(iw)%im(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%im',ivara=iscr)
  DEALLOCATE(iscr)

  ALLOCATE (iscr(7,nwa))
  DO iw = 1,nwa
     iscr(1:7,iw) = itab_w(iw)%iv(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%iv',ivara=iscr)
  DEALLOCATE(iscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%dirv(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%dirv',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%fwv(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%fwv',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%fww(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%fww',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%farm(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%farm',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%farv(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%farv',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gxps1(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%gxps1',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gyps1(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%gyps1',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gxps2(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%gxps2',rvara=rscr)
  DEALLOCATE(rscr)

  ALLOCATE (rscr(7,nwa))
  DO iw = 1,nwa
     rscr(1:7,iw) = itab_w(iw)%gyps2(1:7)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%gyps2',rvara=rscr)
  DEALLOCATE(rscr)

  idims(1) = 9

  ALLOCATE (iscr(9,nwa))
  DO iw = 1,nwa
     iscr(1:9,iw) = itab_w(iw)%iu(1:9)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%iu',ivara=iscr)
  DEALLOCATE(iscr)

  ALLOCATE (iscr(9,nwa))
  DO iw = 1,nwa
     iscr(1:9,iw) = itab_w(iw)%iw(1:9)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%iw',ivara=iscr)
  DEALLOCATE(iscr)

  ALLOCATE (rscr(9,nwa))
  DO iw = 1,nwa
     rscr(1:9,iw) = itab_w(iw)%fwu(1:9)
  ENDDO
  CALL shdf5_orec(ndims,idims,'itab_w%fwu',rvara=rscr)
  DEALLOCATE(rscr)

  ! Check whether LAND/SEA models are used

  IF (isfcl == 1) THEN

     ! Write SEAFLUX VALUES

     ndims = 1
     idims(1) = 1
     idims(2) = 1

     CALL shdf5_orec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
     CALL shdf5_orec(ndims, idims, 'NSFPATS',ivars=nsfpats)

     idims(1) = nseaflux

     CALL shdf5_orec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%ifglobe)
     CALL shdf5_orec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
     CALL shdf5_orec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
     CALL shdf5_orec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iwls)
     CALL shdf5_orec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
     CALL shdf5_orec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
     CALL shdf5_orec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
     CALL shdf5_orec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
     CALL shdf5_orec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
     CALL shdf5_orec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
     CALL shdf5_orec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
     CALL shdf5_orec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sfc)

     IF (nsfpats > 0) THEN

        idims(1) = nsfpats

        CALL shdf5_orec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

        ndims = 2
        idims(1) = 5
        idims(2) = nsfpats

        CALL shdf5_orec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
        CALL shdf5_orec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
        CALL shdf5_orec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

     ENDIF

     ! Write LANDFLUX VALUES

     ndims = 1
     idims(1) = 1
     idims(2) = 1

     CALL shdf5_orec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
     CALL shdf5_orec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

     idims(1) = nlandflux

     CALL shdf5_orec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ifglobe)
     CALL shdf5_orec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
     CALL shdf5_orec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
     CALL shdf5_orec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwls)
     CALL shdf5_orec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
     CALL shdf5_orec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
     CALL shdf5_orec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
     CALL shdf5_orec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
     CALL shdf5_orec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
     CALL shdf5_orec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
     CALL shdf5_orec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
     CALL shdf5_orec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_sfc)

     IF (nlfpats > 0) THEN

        idims(1) = nlfpats

        CALL shdf5_orec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

        ndims = 2
        idims(1) = 5
        idims(2) = nlfpats

        CALL shdf5_orec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
        CALL shdf5_orec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
        CALL shdf5_orec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

     ENDIF

  ENDIF

  ! Check whether NUDGING arrays are used

  IF (mdomain == 0 .AND. nudflag > 0) THEN

     ndims = 1
     idims(1) = 1
     idims(2) = 1

     CALL shdf5_orec(ndims, idims, 'NUDNXP' , ivars=nudnxp)
     CALL shdf5_orec(ndims, idims, 'NWNUD'  , ivars=nwnud)

     idims(1) = nwnud

     CALL shdf5_orec(ndims, idims, 'XEWNUD'  , rvara=xewnud)
     CALL shdf5_orec(ndims, idims, 'YEWNUD'  , rvara=yewnud)
     CALL shdf5_orec(ndims, idims, 'ZEWNUD'  , rvara=zewnud)

     CALL shdf5_orec(ndims,idims,'itab_wnud%npoly',ivara=itab_wnud(:)%npoly)

     ALLOCATE (iscr(6,nwnud))

     DO iwnud = 1,nwnud
        iscr(1:6,iwnud) = itab_wnud(iwnud)%iwnud(1:6)
     ENDDO

     ndims = 2
     idims(1) = 6
     idims(2) = nwnud

     CALL shdf5_orec(ndims,idims,'itab_wnud%iwnud' ,ivara=iscr)

     DEALLOCATE(iscr)

  ENDIF

  ! Close GRIDFILE

  CALL shdf5_close()

  RETURN
END SUBROUTINE gridfile_write

!===============================================================================

SUBROUTINE gridfile_read()

USE max_dims,   ONLY: maxngrdll
USE misc_coms,  ONLY: io6, ngrids, gridfile, mdomain, meshtype, nzp, nxp, &
                      itopoflg, deltax, deltaz, dzmax, dzrat, &
                      zbase, dzbase, ztop, dztop, nzaux, zaux, dzaux, &
                      ngrdll, grdrad, grdlat, grdlon, meshtype
USE mem_ijtabs, ONLY: mloops_m, mloops_u, mloops_v, mloops_w, mrls, &
                      itab_m, itab_u, itab_v, itab_w, alloc_itabs
USE mem_grid,   ONLY: nza, nma, nua, nva, nwa, &
                      mza, mma, mua, mva, mwa, nsw_max, &
                      zm, zt, dzm, dzt, dzim, dzit, &
                      zfacm, zfact, zfacim, zfacit, &
                      lpm, lpu, lcu, lpv, lcv, lpw, lsw, &
                      topm, topw, xem, yem, zem, xeu, yeu, zeu, &
                      xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                      dnu, dniu, dnv, dniv, arw0, arm0, &
                      glatw, glonw, glatm, glonm, glatu, glonu, glatv, glonv, &
                      aru, arv, volui, volvi, arw, volwi, volt, volti, &
                      alloc_gridz, alloc_xyzem, alloc_xyzew, &
                      alloc_grid1, alloc_grid2
USE leaf_coms,  ONLY: isfcl
USE mem_sflux,  ONLY: nseaflux, nlandflux, mseaflux, mlandflux, &
                      nsfpats, nlfpats, msfpats, mlfpats, nsfpatm, nlfpatm, &
                      seaflux, landflux, &
                      xemsfpat, yemsfpat, zemsfpat, &
                      xemlfpat, yemlfpat, zemlfpat

USE hdf5_utils, ONLY: shdf5_irec, shdf5_open, shdf5_close
USE mem_para,   ONLY: myrank

USE mem_nudge,  ONLY: nudflag, nudnxp, nwnud, mwnud, itab_wnud, &
                       xewnud, yewnud, zewnud, alloc_nudge1

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

IMPLICIT NONE

INTEGER :: im, iu, iv, iw, iwnud, izaux

INTEGER :: ierr

INTEGER :: ngr, i
INTEGER :: ndims, idims(2)

INTEGER :: ngrids0, mdomain0, meshtype0, nxp0, nzp0, itopoflg0, isfcl0, nzaux0

REAL    :: deltax0, deltaz0, dzrat0, dzmax0, zbase0, dzbase0, ztop0, dztop0

LOGICAL :: exans

INTEGER, ALLOCATABLE :: ngrdll0(:)
REAL,    ALLOCATABLE :: grdrad0(:)
REAL,    ALLOCATABLE :: grdlat0(:,:)
REAL,    ALLOCATABLE :: grdlon0(:,:)
REAL,    ALLOCATABLE :: zaux0(:)
REAL,    ALLOCATABLE :: dzaux0(:)

! Scratch arrays for copying input

LOGICAL, ALLOCATABLE :: lscr(:,:)
INTEGER, ALLOCATABLE :: iscr(:,:)

REAL, ALLOCATABLE :: rscr(:,:)

! Check if grid file exists

INQUIRE(file=gridfile, exist=exans)

IF (exans) THEN

! Grid file exists.  Open, read, and close file.

   WRITE(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(io6,*) 'Opening grid file ', TRIM(gridfile)
   WRITE(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   CALL shdf5_open(TRIM(gridfile),'R')

! Read the grid information that exists in namelist

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   CALL shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp0)
   CALL shdf5_irec(ndims, idims, 'NXP'     , ivars=nxp0)
   CALL shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain0)
   CALL shdf5_irec(ndims, idims, 'MESHTYPE', ivars=meshtype0)
   CALL shdf5_irec(ndims, idims, 'NGRIDS'  , ivars=ngrids0)
   CALL shdf5_irec(ndims, idims, 'ISFCL'   , ivars=isfcl0)
   CALL shdf5_irec(ndims, idims, 'ITOPOFLG', ivars=itopoflg0)
   CALL shdf5_irec(ndims, idims, 'DELTAX'  , rvars=deltax0)
   CALL shdf5_irec(ndims, idims, 'DELTAZ'  , rvars=deltaz0)
   CALL shdf5_irec(ndims, idims, 'DZRAT'   , rvars=dzrat0)
   CALL shdf5_irec(ndims, idims, 'DZMAX'   , rvars=dzmax0)
   CALL shdf5_irec(ndims, idims, 'ZBASE'   , rvars=zbase0)
   CALL shdf5_irec(ndims, idims, 'DZBASE'  , rvars=dzbase0)
   CALL shdf5_irec(ndims, idims, 'ZTOP'    , rvars=ztop0)
   CALL shdf5_irec(ndims, idims, 'DZTOP'   , rvars=dztop0)

   CALL shdf5_irec(ndims, idims, 'NZAUX'   , ivars=nzaux0)

   ALLOCATE( ngrdll0 (ngrids0) )
   ALLOCATE( grdrad0 (ngrids0) )
   ALLOCATE( grdlat0 (ngrids0, maxngrdll) )
   ALLOCATE( grdlon0 (ngrids0, maxngrdll) )
   ALLOCATE( zaux0   (nzaux0))
   ALLOCATE( dzaux0  (nzaux0))

   IF (nzaux > 0) THEN
      idims(1) = nzaux

      CALL shdf5_irec(ndims, idims, 'ZAUX'   , rvara=zaux0)
      CALL shdf5_irec(ndims, idims, 'DZAUX'  , rvara=dzaux0)
   ENDIF

   idims(1) = ngrids0

   CALL shdf5_irec(ndims, idims, 'NGRDLL' , ivara=ngrdll0)
   CALL shdf5_irec(ndims, idims, 'GRDRAD' , rvara=grdrad0)

   ndims = 2
   idims(1) = ngrids0
   idims(2) = maxngrdll

   CALL shdf5_irec(ndims, idims, 'GRDLAT', rvara=grdlat0)
   CALL shdf5_irec(ndims, idims, 'GRDLON', rvara=grdlon0)

! Check equality between grid file information and namelist variables

   ierr = 0

   IF (nzp0      /= nzp     ) ierr = 1 
   IF (nxp0      /= nxp     ) ierr = 1 
   IF (mdomain0  /= mdomain ) ierr = 1 
   IF (meshtype0 /= meshtype) ierr = 1 
   IF (ngrids0   /= ngrids  ) ierr = 1 
   IF (isfcl0    /= isfcl   ) ierr = 1 
   IF (itopoflg0 /= itopoflg) ierr = 1 
   IF (nzaux0    /= nzaux   ) ierr = 1 

   IF (ABS(deltax0 - deltax) > 1.e-3) ierr = 1 
   IF (ABS(deltaz0 - deltaz) > 1.e-3) ierr = 1 
   IF (ABS(dzrat0  - dzrat ) > 1.e-3) ierr = 1 
   IF (ABS(dzmax0  - dzmax ) > 1.e-3) ierr = 1 
   IF (ABS(zbase0  - zbase ) > 1.e-3) ierr = 1 
   IF (ABS(dzbase0 - dzbase) > 1.e-3) ierr = 1 
   IF (ABS(ztop0   - ztop  ) > 1.e-3) ierr = 1 
   IF (ABS(dztop0  - dztop ) > 1.e-3) ierr = 1 

   DO ngr = 2, MIN(ngrids0,ngrids)
      IF (ABS(ngrdll0 (ngr) - ngrdll (ngr)) > 1.e1 ) ierr = 1
      IF (ABS(grdrad0 (ngr) - grdrad (ngr)) > 1.e1 ) ierr = 1

      DO i = 1,ngrdll0(ngr)
         IF (ABS(grdlat0(ngr,i) - grdlat(ngr,i)) > 1.e-3) ierr = 1
         IF (ABS(grdlon0(ngr,i) - grdlon(ngr,i)) > 1.e-3) ierr = 1
      ENDDO
   ENDDO

   DO izaux = 1, MIN(nzaux0,nzaux)
      IF (ABS(zaux0 (izaux) - zaux (izaux)) > 1.e1 ) ierr = 1
      IF (ABS(dzaux0(izaux) - dzaux(izaux)) > 1.e1 ) ierr = 1
   ENDDO

   IF (ierr == 1) THEN

      WRITE(io6,*) 'GRIDFILE mismatch with OLAMIN namelist: Stopping model run'
      WRITE(io6,*) 'Values: gridfile, namelist'
      WRITE(io6,*) '-----------------------------------------------'
      WRITE(io6,*)              'nzp:      ',nzp0     ,nzp
      WRITE(io6,*)              'nxp:      ',nxp0     ,nxp
      WRITE(io6,*)              'mdomain:  ',mdomain0 ,mdomain
      WRITE(io6,*)              'meshtype: ',meshtype0,meshtype
      WRITE(io6,*)              'ngrids:   ',ngrids0  ,ngrids
      WRITE(io6,*)              'isfcl:    ',isfcl0   ,isfcl
      WRITE(io6,*)              'itopoflg: ',itopoflg0,itopoflg
      WRITE(io6,*)              'deltax:   ',deltax0  ,deltax
      WRITE(io6,*)              'deltaz:   ',deltaz0  ,deltaz
      WRITE(io6,*)              'dzrat:    ',dzrat0   ,dzrat
      WRITE(io6,*)              'dzmax:    ',dzmax0   ,dzmax
      WRITE(io6,*)              'zbase:    ',zbase0   ,zbase
      WRITE(io6,*)              'dzbase:   ',dzbase0  ,dzbase
      WRITE(io6,*)              'ztop:     ',ztop0    ,ztop
      WRITE(io6,*)              'dztop:    ',dztop0   ,dztop
      WRITE(io6,*)              'nzaux:    ',nzaux0   ,nzaux
      WRITE(io6,*) ' '
      WRITE(io6, '(a,20f10.1)') 'zaux0:    ',zaux0 (1:nzaux)
      WRITE(io6, '(a,20f10.1)') 'zaux:     ',zaux  (1:nzaux)
      WRITE(io6,*) ' '
      WRITE(io6, '(a,20f10.1)') 'dzaux0:   ',dzaux0 (1:nzaux)
      WRITE(io6, '(a,20f10.1)') 'dzaux:    ',dzaux  (1:nzaux)
      WRITE(io6,*) ' '
      WRITE(io6, '(a,20i12)')   'ngrdll0:  ',ngrdll0 (1:ngrids)
      WRITE(io6, '(a,20i12)')   'ngrdll:   ',ngrdll  (1:ngrids)
      WRITE(io6,*) ' '
      WRITE(io6, '(a,20f12.1)') 'grdrad0:  ',grdrad0 (1:ngrids)
      WRITE(io6, '(a,20f12.1)') 'grdrad:   ',grdrad  (1:ngrids)
      WRITE(io6,*) ' '

      DO ngr = 2, MIN(ngrids0,ngrids)
         WRITE(io6, '(a,i5)') 'ngr: ',ngr
         WRITE(io6,*) ' '
         WRITE(io6, '(a,20f10.3)') 'grdlat0: ',grdlat0(ngr,1:ngrdll(ngr))
         WRITE(io6, '(a,20f10.3)') 'grdlat:  ',grdlat (ngr,1:ngrdll(ngr))
         WRITE(io6,*) ' '
         WRITE(io6, '(a,20f10.3)') 'grdlon0: ',grdlon0(ngr,1:ngrdll(ngr))
         WRITE(io6, '(a,20f10.3)') 'grdlon:  ',grdlon (ngr,1:ngrdll(ngr))
         WRITE(io6,*) ' '
      ENDDO

      WRITE(io6,*) '-----------------------------------------------'

      STOP 'stop - gridfile mismatch'
   
   ENDIF

   DEALLOCATE (ngrdll0, grdrad0, grdlat0, grdlon0)

! Read the grid dimensions

   CALL shdf5_irec(ndims, idims, 'NZA'    , ivars=nza)
   CALL shdf5_irec(ndims, idims, 'NMA'    , ivars=nma)
   CALL shdf5_irec(ndims, idims, 'NUA'    , ivars=nua)
   CALL shdf5_irec(ndims, idims, 'NVA'    , ivars=nva)
   CALL shdf5_irec(ndims, idims, 'NWA'    , ivars=nwa)
   CALL shdf5_irec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
   CALL shdf5_irec(ndims, idims, 'MRLS'   , ivars=mrls)

! Copy grid dimensions

   mza = nza
   mma = nma
   mva = nva
   mua = nua
   mwa = nwa

! Allocate and read grid structure variables

   CALL alloc_gridz()
   CALL alloc_itabs(meshtype,nma,nua,nva,nwa)
   CALL alloc_xyzem(nma)
   CALL alloc_xyzew(nwa)
   CALL alloc_grid1(meshtype)
   CALL alloc_grid2(meshtype)

   idims(1) = nza

   CALL shdf5_irec(ndims, idims, 'ZM'    , rvara=zm)
   CALL shdf5_irec(ndims, idims, 'ZT'    , rvara=zt)
   CALL shdf5_irec(ndims, idims, 'DZM'   , rvara=dzm)
   CALL shdf5_irec(ndims, idims, 'DZT'   , rvara=dzt)
   CALL shdf5_irec(ndims, idims, 'DZIM'  , rvara=dzim)
   CALL shdf5_irec(ndims, idims, 'DZIT'  , rvara=dzit)
   CALL shdf5_irec(ndims, idims, 'ZFACM' , rvara=zfacm)
   CALL shdf5_irec(ndims, idims, 'ZFACT' , rvara=zfact)
   CALL shdf5_irec(ndims, idims, 'ZFACIM', rvara=zfacim)
   CALL shdf5_irec(ndims, idims, 'ZFACIT', rvara=zfacit)

   idims(1) = nma

   CALL shdf5_irec(ndims, idims, 'LPM'  , ivara=lpm)
   CALL shdf5_irec(ndims, idims, 'TOPM' , rvara=topm)
   CALL shdf5_irec(ndims, idims, 'XEM'  , rvara=xem)
   CALL shdf5_irec(ndims, idims, 'YEM'  , rvara=yem)
   CALL shdf5_irec(ndims, idims, 'ZEM'  , rvara=zem)
   CALL shdf5_irec(ndims, idims, 'ARM0' , rvara=arm0)
   CALL shdf5_irec(ndims, idims, 'GLATM', rvara=glatm)
   CALL shdf5_irec(ndims, idims, 'GLONM', rvara=glonm)

   IF (meshtype == 1) THEN

      idims(1) = nua

      CALL shdf5_irec(ndims, idims, 'LPU'  , ivara=lpu)
      CALL shdf5_irec(ndims, idims, 'LCU'  , ivara=lcu)
      CALL shdf5_irec(ndims, idims, 'XEU'  , rvara=xeu)
      CALL shdf5_irec(ndims, idims, 'YEU'  , rvara=yeu)
      CALL shdf5_irec(ndims, idims, 'ZEU'  , rvara=zeu)
      CALL shdf5_irec(ndims, idims, 'GLATU', rvara=glatu)
      CALL shdf5_irec(ndims, idims, 'GLONU', rvara=glonu)

   ELSE

      idims(1) = nva

      CALL shdf5_irec(ndims, idims, 'LPV' , ivara=lpv)
      CALL shdf5_irec(ndims, idims, 'LCV' , ivara=lcv)
      CALL shdf5_irec(ndims, idims, 'XEV' , rvara=xev)
      CALL shdf5_irec(ndims, idims, 'YEV' , rvara=yev)
      CALL shdf5_irec(ndims, idims, 'ZEV' , rvara=zev)
      CALL shdf5_irec(ndims, idims, 'GLATV', rvara=glatv)
      CALL shdf5_irec(ndims, idims, 'GLONV', rvara=glonv)

   ENDIF

   CALL shdf5_irec(ndims, idims, 'DNU' , rvara=dnu)
   CALL shdf5_irec(ndims, idims, 'DNV' , rvara=dnv)
   CALL shdf5_irec(ndims, idims, 'DNIU', rvara=dniu)
   CALL shdf5_irec(ndims, idims, 'DNIV', rvara=dniv)
   CALL shdf5_irec(ndims, idims, 'UNX' , rvara=unx)
   CALL shdf5_irec(ndims, idims, 'UNY' , rvara=uny)
   CALL shdf5_irec(ndims, idims, 'UNZ' , rvara=unz)
   CALL shdf5_irec(ndims, idims, 'VNX' , rvara=vnx)
   CALL shdf5_irec(ndims, idims, 'VNY' , rvara=vny)
   CALL shdf5_irec(ndims, idims, 'VNZ' , rvara=vnz)

   idims(1) = nwa

   CALL shdf5_irec(ndims, idims, 'LPW'  , ivara=lpw)
   CALL shdf5_irec(ndims, idims, 'LSW'  , ivara=lsw)
   CALL shdf5_irec(ndims, idims, 'XEW'  , rvara=xew)
   CALL shdf5_irec(ndims, idims, 'YEW'  , rvara=yew)
   CALL shdf5_irec(ndims, idims, 'ZEW'  , rvara=zew)
   CALL shdf5_irec(ndims, idims, 'TOPW' , rvara=topw)
   CALL shdf5_irec(ndims, idims, 'ARW0' , rvara=arw0)
   CALL shdf5_irec(ndims, idims, 'GLATW', rvara=glatw)
   CALL shdf5_irec(ndims, idims, 'GLONW', rvara=glonw)
   CALL shdf5_irec(ndims, idims, 'WNX'  , rvara=wnx)
   CALL shdf5_irec(ndims, idims, 'WNY'  , rvara=wny)
   CALL shdf5_irec(ndims, idims, 'WNZ'  , rvara=wnz)

   ndims = 2
   idims(1) = nza

   IF (meshtype == 1) THEN

      idims(2) = nua

      CALL shdf5_irec(ndims, idims, 'VOLUI', rvara=volui)

   ELSE

      idims(2) = nva

      CALL shdf5_irec(ndims, idims, 'ARV'  , rvara=arv)
      CALL shdf5_irec(ndims, idims, 'VOLVI', rvara=volvi)
   
   ENDIF

   CALL shdf5_irec(ndims, idims, 'ARU'  , rvara=aru)

   idims(2) = nwa

   CALL shdf5_irec(ndims, idims, 'ARW'  , rvara=arw)
   CALL shdf5_irec(ndims, idims, 'VOLWI', rvara=volwi)
   CALL shdf5_irec(ndims, idims, 'VOLT' , dvara=volt)
   CALL shdf5_irec(ndims, idims, 'VOLTI', dvara=volti)
   
! Read ITAB_M SCALARS

   ndims = 1
   idims(1) = nma
   idims(2) = 1

   CALL shdf5_irec(ndims,idims,'itab_m%npoly'    ,ivara=itab_m(:)%npoly)
   CALL shdf5_irec(ndims,idims,'itab_m%itopm'    ,ivara=itab_m(:)%itopm)
   CALL shdf5_irec(ndims,idims,'itab_m%imglobe'  ,ivara=itab_m(:)%imglobe)
   CALL shdf5_irec(ndims,idims,'itab_m%mrlm'     ,ivara=itab_m(:)%mrlm)
   CALL shdf5_irec(ndims,idims,'itab_m%mrlm_orig',ivara=itab_m(:)%mrlm_orig)
   CALL shdf5_irec(ndims,idims,'itab_m%mrow'     ,ivara=itab_m(:)%mrow)
   CALL shdf5_irec(ndims,idims,'itab_m%mrowh'    ,ivara=itab_m(:)%mrowh)

! Read ITAB_M ARRAYS

   ndims = 2
   idims(1) = mloops_m
   idims(2) = nma

   ALLOCATE (lscr(mloops_m,nma))
   CALL shdf5_irec(ndims,idims,'itab_m%loop',lvara=lscr)
   DO im = 1,nma
      itab_m(im)%loop(1:mloops_m) = lscr(1:mloops_m,im)
   ENDDO
   DEALLOCATE (lscr)

   idims(1) = 7

   ALLOCATE (iscr(7,nma))
   IF (meshtype == 1) THEN
      CALL shdf5_irec(ndims,idims,'itab_m%iu',ivara=iscr)
      DO im = 1,nma
         itab_m(im)%iu(1:7) = iscr(1:7,im)
      ENDDO
   ELSE
      CALL shdf5_irec(ndims,idims,'itab_m%iv',ivara=iscr)
      DO im = 1,nma
         itab_m(im)%iv(1:7) = iscr(1:7,im)
      ENDDO
   ENDIF
   DEALLOCATE (iscr)

   ALLOCATE (iscr(7,nma))
   CALL shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr)
   DO im = 1,nma
      itab_m(im)%iw(1:7) = iscr(1:7,im)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (rscr(7,nma))
   CALL shdf5_irec(ndims,idims,'itab_m%fmw',rvara=rscr)
   DO im = 1,nma
      itab_m(im)%fmw(1:7) = rscr(1:7,im)
   ENDDO
   DEALLOCATE (rscr)

   IF (meshtype == 1) THEN

! Read ITAB_U SCALARS

      ndims = 1
      idims(1) = nua
      idims(2) = 1

      CALL shdf5_irec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)
      CALL shdf5_irec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)
      CALL shdf5_irec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)

      CALL shdf5_irec(ndims,idims,'itab_u%gcf36'  ,rvara=itab_u(:)%gcf36)
      CALL shdf5_irec(ndims,idims,'itab_u%gcf45'  ,rvara=itab_u(:)%gcf45)

      CALL shdf5_irec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
      CALL shdf5_irec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)
      CALL shdf5_irec(ndims,idims,'itab_u%crossmm',rvara=itab_u(:)%crossmm)
      CALL shdf5_irec(ndims,idims,'itab_u%crossww',rvara=itab_u(:)%crossww)

! Read ITAB_U ARRAYS

      ndims = 2
      idims(1) = mloops_u
      idims(2) = nua

      ALLOCATE (lscr(mloops_u,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%loop',lvara=lscr)
      DO iu = 1,nua
         itab_u(iu)%loop(1:mloops_u) = lscr(1:mloops_u,iu)
      ENDDO
      DEALLOCATE (lscr)

      idims(1) = 2

      ALLOCATE (iscr(2,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%im',ivara=iscr)
      DO iu = 1,nua
         itab_u(iu)%im(1:2) = iscr(1:2,iu)
      ENDDO
      DEALLOCATE (iscr)

      ALLOCATE (rscr(2,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%vxw_u',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%vxw_u(1:2) = rscr(1:2,iu)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(2,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%vyw_u',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%vyw_u(1:2) = rscr(1:2,iu)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 4

      ALLOCATE (rscr(4,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%diru',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%diru(1:4) = rscr(1:4,iu)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(4,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%tuu',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%tuu(1:4) = rscr(1:4,iu)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(4,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%guw',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%guw(1:4) = rscr(1:4,iu)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(4,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%vxu_u',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%vxu_u(1:4) = rscr(1:4,iu)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(4,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%vyu_u',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%vyu_u(1:4) = rscr(1:4,iu)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 6

      ALLOCATE (iscr(6,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%iw',ivara=iscr)
      DO iu = 1,nua
         itab_u(iu)%iw(1:6) = iscr(1:6,iu)
      ENDDO
      DEALLOCATE (iscr)

      ALLOCATE (rscr(6,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%fuw',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%fuw(1:6) = rscr(1:6,iu)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 12

      ALLOCATE (iscr(12,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%iu',ivara=iscr)
      DO iu = 1,nua
         itab_u(iu)%iu(1:12) = iscr(1:12,iu)
      ENDDO
      DEALLOCATE (iscr)

      ALLOCATE (rscr(12,nua))
      CALL shdf5_irec(ndims,idims,'itab_u%fuu',rvara=rscr)
      DO iu = 1,nua
         itab_u(iu)%fuu(1:12) = rscr(1:12,iu)
      ENDDO
      DEALLOCATE (rscr)

   ENDIF

   IF (meshtype == 2) THEN

! Read ITAB_V SCALARS

      ndims = 1
      idims(1) = nva
      idims(2) = 1

      CALL shdf5_irec(ndims,idims,'itab_v%ivp'      ,ivara=itab_v(:)%ivp)
      CALL shdf5_irec(ndims,idims,'itab_v%mrlv'     ,ivara=itab_v(:)%mrlv)
      CALL shdf5_irec(ndims,idims,'itab_v%ivglobe'  ,ivara=itab_v(:)%ivglobe)

! Read ITAB_V ARRAYS

      ndims = 2
      idims(1) = mloops_v
      idims(2) = nva

      ALLOCATE (lscr(mloops_v,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%loop',lvara=lscr)
      DO iv = 1,nva
         itab_v(iv)%loop(1:mloops_v) = lscr(1:mloops_v,iv)
      ENDDO
      DEALLOCATE (lscr)

      idims(1) = 2

      ALLOCATE (rscr(2,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%farw',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%farw(1:2) = rscr(1:2,iv)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(2,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%cosv',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%cosv(1:2) = rscr(1:2,iv)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(2,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%sinv',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%sinv(1:2) = rscr(1:2,iv)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(2,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%dxps',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%dxps(1:2) = rscr(1:2,iv)
      ENDDO
      DEALLOCATE (rscr)

      ALLOCATE (rscr(2,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%dyps',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%dyps(1:2) = rscr(1:2,iv)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 4

      ALLOCATE (iscr(4,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%iw',ivara=iscr)
      DO iv = 1,nva
         itab_v(iv)%iw(1:4) = iscr(1:4,iv)
      ENDDO
      DEALLOCATE (iscr)

      ALLOCATE (rscr(4,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%fvw',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%fvw (1:4) = rscr(1:4,iv)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 6

      ALLOCATE (iscr(6,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%im',ivara=iscr)
      DO iv = 1,nva
         itab_v(iv)%im(1:6) = iscr(1:6,iv)
      ENDDO
      DEALLOCATE (iscr)

      idims(1) = 12

      ALLOCATE (rscr(12,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%fvv',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%fvv(1:12) = rscr(1:12,iv)
      ENDDO
      DEALLOCATE (rscr)

      idims(1) = 16

      ALLOCATE (iscr(16,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%iv',ivara=iscr)
      DO iv = 1,nva
         itab_v(iv)%iv(1:16) = iscr(1:16,iv)
      ENDDO
      DEALLOCATE (iscr)

      ALLOCATE (rscr(16,nva))
      CALL shdf5_irec(ndims,idims,'itab_v%fuv',rvara=rscr)
      DO iv = 1,nva
         itab_v(iv)%fuv(1:16) = rscr(1:16,iv)
      ENDDO
      DEALLOCATE (rscr)

   ENDIF

! Read ITAB_W SCALARS

   ndims = 1
   idims(1) = nwa
   idims(2) = 1

   CALL shdf5_irec(ndims,idims,'itab_w%npoly'    ,ivara=itab_w(:)%npoly)
   CALL shdf5_irec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)
   CALL shdf5_irec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
   CALL shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
   CALL shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
   CALL shdf5_irec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
   CALL shdf5_irec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

   CALL shdf5_irec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)
   CALL shdf5_irec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)
   CALL shdf5_irec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

   CALL shdf5_irec(ndims,idims,'itab_w%unx_w' ,rvara=itab_w(:)%unx_w)
   CALL shdf5_irec(ndims,idims,'itab_w%uny_w' ,rvara=itab_w(:)%uny_w)

   CALL shdf5_irec(ndims,idims,'itab_w%vnx_w' ,rvara=itab_w(:)%vnx_w)
   CALL shdf5_irec(ndims,idims,'itab_w%vny_w' ,rvara=itab_w(:)%vny_w)
   CALL shdf5_irec(ndims,idims,'itab_w%vnz_w' ,rvara=itab_w(:)%vnz_w)

! Read ITAB_W ARRAYS

   ndims = 2
   idims(1) = mloops_w
   idims(2) = nwa

   ALLOCATE (lscr(mloops_w,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%loop',lvara=lscr)
   DO iw = 1,nwa
      itab_w(iw)%loop(1:mloops_w) = lscr(1:mloops_w,iw)
   ENDDO
   DEALLOCATE (lscr)

   idims(1) = 3

   ALLOCATE (iscr(3,nwa))
   call shdf5_irec(ndims,idims,'itab_w%iwnud',ivara=iscr)
   DO iw = 1,nwa
      itab_w(iw)%iwnud(1:3) = iscr(1:3,iw)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (rscr(3,nwa))
   call shdf5_irec(ndims,idims,'itab_w%fnud',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%fnud(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%diru',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%diru(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%vxu',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%vxu(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%vyu',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%vyu(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%vzu',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%vzu(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%vxu_w',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%vxu_w(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(3,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%vyu_w',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%vyu_w(1:3) = rscr(1:3,iw)
   ENDDO
   DEALLOCATE (rscr)

   idims(1) = 7

   ALLOCATE (iscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%im',ivara=iscr)
   DO iw = 1,nwa
      itab_w(iw)%im(1:7) = iscr(1:7,iw)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (iscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%iv',ivara=iscr)
   DO iw = 1,nwa
      itab_w(iw)%iv(1:7) = iscr(1:7,iw)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%dirv',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%dirv(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%fwv',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%fwv(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%fww',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%fww(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%farm',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%farm(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%farv',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%farv(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%gxps1',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%gxps1(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%gyps1',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%gyps1(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%gxps2',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%gxps2(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   ALLOCATE (rscr(7,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%gyps2',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%gyps2(1:7) = rscr(1:7,iw)
   ENDDO
   DEALLOCATE (rscr)

   idims(1) = 9

   ALLOCATE (iscr(9,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%iu',ivara=iscr)
   DO iw = 1,nwa
      itab_w(iw)%iu(1:9) = iscr(1:9,iw)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (iscr(9,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%iw',ivara=iscr)
   DO iw = 1,nwa
      itab_w(iw)%iw(1:9) = iscr(1:9,iw)
   ENDDO
   DEALLOCATE (iscr)

   ALLOCATE (rscr(9,nwa))
   CALL shdf5_irec(ndims,idims,'itab_w%fwu',rvara=rscr)
   DO iw = 1,nwa
      itab_w(iw)%fwu(1:9) = rscr(1:9,iw)
   ENDDO
   DEALLOCATE (rscr)

! Check whether LAND/SEA models are used

   IF (isfcl == 1) THEN

! Read SEAFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      CALL shdf5_irec(ndims, idims, 'NSEAFLUX',ivars=nseaflux)
      CALL shdf5_irec(ndims, idims, 'NSFPATS' ,ivars=nsfpats)

      mseaflux = nseaflux
      msfpats = nsfpats

      ALLOCATE (seaflux(nseaflux))

      ALLOCATE (nsfpatm(nsfpats))

      ALLOCATE (xemsfpat(5,nsfpats))
      ALLOCATE (yemsfpat(5,nsfpats))
      ALLOCATE (zemsfpat(5,nsfpats))

      idims(1) = nseaflux

      CALL shdf5_irec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%ifglobe)
      CALL shdf5_irec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
      CALL shdf5_irec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
      CALL shdf5_irec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iwls)
      CALL shdf5_irec(ndims,idims,'seaflux%jpats'   ,ivara=seaflux(:)%jpats)
      CALL shdf5_irec(ndims,idims,'seaflux%ipat'    ,ivara=seaflux(:)%ipat)
      CALL shdf5_irec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
      CALL shdf5_irec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
      CALL shdf5_irec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
      CALL shdf5_irec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
      CALL shdf5_irec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
      CALL shdf5_irec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sfc)

      IF (nsfpats > 0) THEN

         idims(1) = nsfpats

         CALL shdf5_irec(ndims,idims,'nsfpatm' ,ivara=nsfpatm)

         ndims = 2
         idims(1) = 5
         idims(2) = nsfpats

         CALL shdf5_irec(ndims,idims,'xemsfpat' ,rvara=xemsfpat)
         CALL shdf5_irec(ndims,idims,'yemsfpat' ,rvara=yemsfpat)
         CALL shdf5_irec(ndims,idims,'zemsfpat' ,rvara=zemsfpat)

      ENDIF

! Read LANDFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      CALL shdf5_irec(ndims, idims, 'NLANDFLUX',ivars=nlandflux)
      CALL shdf5_irec(ndims, idims, 'NLFPATS'  ,ivars=nlfpats)

      mlandflux = nlandflux
      mlfpats = nlfpats

      ALLOCATE (landflux(nlandflux))

      ALLOCATE (nlfpatm(nlfpats))

      ALLOCATE (xemlfpat(5,nlfpats))
      ALLOCATE (yemlfpat(5,nlfpats))
      ALLOCATE (zemlfpat(5,nlfpats))

      idims(1) = nlandflux

      CALL shdf5_irec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ifglobe)
      CALL shdf5_irec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
      CALL shdf5_irec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
      CALL shdf5_irec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwls)
      CALL shdf5_irec(ndims,idims,'landflux%jpats'   ,ivara=landflux(:)%jpats)
      CALL shdf5_irec(ndims,idims,'landflux%ipat'    ,ivara=landflux(:)%ipat)
      CALL shdf5_irec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
      CALL shdf5_irec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
      CALL shdf5_irec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
      CALL shdf5_irec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
      CALL shdf5_irec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
      CALL shdf5_irec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_sfc)

      IF (nlfpats > 0) THEN

         idims(1) = nlfpats

         CALL shdf5_irec(ndims,idims,'nlfpatm' ,ivara=nlfpatm)

         ndims = 2
         idims(1) = 5
         idims(2) = nlfpats

         CALL shdf5_irec(ndims,idims,'xemlfpat' ,rvara=xemlfpat)
         CALL shdf5_irec(ndims,idims,'yemlfpat' ,rvara=yemlfpat)
         CALL shdf5_irec(ndims,idims,'zemlfpat' ,rvara=zemlfpat)

      ENDIF

   ENDIF

! Check whether NUDGING arrays are used

   IF (mdomain == 0 .AND. nudflag > 0) THEN

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      CALL shdf5_irec(ndims, idims, 'NUDNXP' , ivars=nudnxp)
      CALL shdf5_irec(ndims, idims, 'NWNUD'  , ivars=nwnud)

      mwnud = nwnud

      idims(1) = nwnud

      CALL alloc_nudge1(nwnud)

      CALL shdf5_irec(ndims, idims, 'XEWNUD'  , rvara=xewnud)
      CALL shdf5_irec(ndims, idims, 'YEWNUD'  , rvara=yewnud)
      CALL shdf5_irec(ndims, idims, 'ZEWNUD'  , rvara=zewnud)

      CALL shdf5_irec(ndims,idims,'itab_wnud%npoly' ,ivara=itab_wnud(:)%npoly)

      ALLOCATE (iscr(6,nwnud))

      ndims = 2
      idims(1) = 6
      idims(2) = nwnud

      CALL shdf5_irec(ndims,idims,'itab_wnud%iwnud',ivara=iscr)

      DO iwnud = 1,nwnud
         itab_wnud(iwnud)%iwnud(1:6) = iscr(1:6,iwnud)
      ENDDO

      DEALLOCATE(iscr)

   ENDIF

! Close the GRIDFILE

   CALL shdf5_close()

ELSE

! Grid file does not exist.

   WRITE(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   WRITE(io6,*) '!!!  Gridfile does not exist:'
   WRITE(io6,*) '!!!  '//TRIM(gridfile)
   WRITE(io6,*) '!!!  Stopping run'
   WRITE(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   
   STOP 'stop - no gridfile'
   
ENDIF

WRITE(io6,*) 'end of gridfile_read '

RETURN
END SUBROUTINE gridfile_read



