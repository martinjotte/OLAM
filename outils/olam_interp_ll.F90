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

subroutine interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,nlevin,nlevout,alon,alat,field,field_ll)

  use mem_grid,   only: mma, mza, mwa, zm, zt, lpw, &
                        xem, yem, zem, xew, yew, zew, glatm, glonm

  use mem_ijtabs, only: itab_m, jtab_m, jtm_vadj, itabg_m, itab_v
  use misc_coms,  only: io6
  use consts_coms, only: eradi,piu180,pio180

  implicit none

  integer, intent(in)    :: nlon, nlat, npts, nlevin, nlevout
  integer, intent(in)    :: ims_loc(npts), lls_loc(npts), ijs_loc(npts)
  real,    intent(in)    :: alon(nlon),alat(nlat)
  real,    intent(in)    :: field(nlevin,mwa)
  real,    intent(inout) :: field_ll(npts,nlevout)

  integer :: k,im,iw,j,jnext,iwnext,ilat,ilon,kll
  integer :: k1,k2,koff

  real :: x(3),y(3),z(3)

  real :: a(nlevin),b(nlevin),c(nlevin),pts(nlevin)
  real :: field_avg(nlevin)

  real :: qx,qy
  real :: dxe, dye, dze
  real :: cosmlon, sinmlon
  real :: cosmlat, sinmlat

  integer :: n, nn, img

  if (npts == 0) return

  if (nlevin == 1) then
     k1 = 1
     k2 = 1
     koff = 0
  else
     k2 = mza
     koff = 1
  endif

  x(1) = 0.
  y(1) = 0.

do nn = 1, npts

   n   = lls_loc(nn)            ! 1d lat/lon index
   img = ims_loc(nn)            ! global index
   im  = itabg_m(img)%im_myrank ! local index

   ! convert 1d lat/lon index to 2d

   ilat = (n-1)/nlon + 1
   ilon = mod(n-1,nlon) + 1

   ! Initialize maximum lpw and field average

   if (nlevin > 1) k1 = maxval( lpw( itab_m(im)%iw(1:3) ) )

   field_avg(1:nlevin) = 0.

   ! store sin/cos of this M point

   sinmlat = sin(glatm(im) * pio180)
   cosmlat = cos(glatm(im) * pio180)
   sinmlon = sin(glonm(im) * pio180)
   cosmlon = cos(glonm(im) * pio180)

   ! Transform current lat-lon point to PS space

   call ll_xy(alat(ilat),alon(ilon),glatm(im),glonm(im),qx,qy)

   ! Loop over all W points that surround current M point to compute average

   pts(:) = 0.0

   do j = 1, 3

      ! Current W point index

      iw = itab_m(im)%iw(j)

      ! Skip this IM point if iw < 2
      ! if (iw < 2) then
      !   go to 9
      ! endif

      do k = k1, k2
         if (nlevin > 1 .and. k < lpw(iw)) cycle
         pts(k) = pts(k) + 1.0
         field_avg(k) = field_avg(k) + field(k,iw)
      enddo
   enddo

   field_avg(k1:k2) = field_avg(k1:k2) / pts(k1:k2)

   ! fill field values for triangle that contains this lat/lon point 

   j = ijs_loc(nn)

   jnext = j + 1
   if (j == 3) jnext = 1

   iw     = itab_m(im)%iw(j)
   iwnext = itab_m(im)%iw(jnext)

   ! Transform iw point to PS coordinates tangent at M point

   dxe = xew(iw) - xem(im)
   dye = yew(iw) - yem(im)
   dze = zew(iw) - zem(im)
   call de_ps(dxe,dye,dze,cosmlat,sinmlat,cosmlon,sinmlon,x(2),y(2))

   ! Transform iwnext point to PS coordinates tangent at M point

   dxe = xew(iwnext) - xem(im)
   dye = yew(iwnext) - yem(im)
   dze = zew(iwnext) - zem(im)
   call de_ps(dxe,dye,dze,cosmlat,sinmlat,cosmlon,sinmlon,x(3),y(3))

   ! Loop over vertical levels

   do k = k1, k2

      z(1) = field_avg(k)

      if (nlevin > 1 .and. k < lpw(iw)) then
         z(2) = field_avg(k)
      else
         z(2) = field(k,iw)
      endif

      if (nlevin > 1 .and. k < lpw(iwnext)) then
         z(3) = field_avg(k)
      else
         z(3) = field(k,iwnext)
      endif

      ! Evaluate interpolation coefficients for current trio of points

      call matrix_3x3(1.  , x(1), y(1),  &
                      1.  , x(2), y(2),  &
                      1.  , x(3), y(3),  &
                      z(1), z(2), z(3),  &
                      a(k), b(k), c(k)   )
   enddo

   ! Interpolate to current field point

   do k = k1, k2
      kll = k - koff
      field_ll(nn,kll) = a(k) + b(k) * qx + c(k) * qy
   enddo

enddo

end subroutine interp_htw_ll

!================================================================================

subroutine find_closest_m_ll(nlon,nlat,alon,alat,ims,ijs)

  use mem_grid,   only: glatm, glonm, mma, nma, xew, yew, zew, arw0, xem, yem, zem
  use mem_ijtabs, only: itab_m, jtab_m, jtm_vadj, itabg_m
  use misc_coms,  only: io6, iparallel, nxp
  use mem_para,   only: mgroupsize, myrank
  use consts_coms,only: pio180, r8, erad

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in)  :: nlon, nlat
  real,    intent(in)  :: alon(nlon), alat(nlat)
  integer, intent(out) :: ims(nlon,nlat)
  integer, intent(out) :: ijs(nlon,nlat)

  real :: alon360(nlon)
  real :: glonm360(mma)
  real :: dllolam(mma)

  integer :: img(nlon), imgs(nlon,mgroupsize)
  integer :: ntri(nlon), ntris(nlon,mgroupsize)

  integer :: im, ier, ilat, ilon, jm, j, iw, np, jnext

  real :: deglon, qx, qy, xw(3), yw(3), fact
  real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v

  real :: sinalat, er_sinalon(nlon)
  real :: cosalat, er_cosalon(nlon)
  real :: xea, yea, zea
  real :: dxe, dye, dze

  real :: cosmlon(mma), sinmlon(mma)
  real :: cosmlat(mma), sinmlat(mma)

  real, parameter :: fuzz = 0.01
  real, parameter :: xo = 0.0
  real, parameter :: yo = 0.0

  do ilon = 1, nlon
     if (alon(ilon) < 0.0) then
        alon360(ilon) = 360.0 + alon(ilon)
     else
        alon360(ilon) = alon(ilon)
     endif

     er_sinalon(ilon) = erad * sin(alon(ilon) * pio180)
     er_cosalon(ilon) = erad * cos(alon(ilon) * pio180)
  enddo

  do jm = 1, jtab_m(jtm_vadj)%jend(1)
     im = jtab_m(jtm_vadj)%im(jm)

     dllolam(im) = maxval( sqrt( arw0( itab_m(im)%iw(1:3) ) ) ) / 108.e3

     if (glonm(im) < 0.0) then
        glonm360(im) = 360.0 + glonm(im)
     else
        glonm360(im) = glonm(im)
     endif

     sinmlat(im) = sin(glatm(im) * pio180)
     cosmlat(im) = cos(glatm(im) * pio180)
     sinmlon(im) = sin(glonm(im) * pio180)
     cosmlon(im) = cos(glonm(im) * pio180)
  enddo

  do ilat = 1, nlat

     fact = 1.0 / cos( min( abs(alat(ilat)), 89.0 ) * pio180 )

     sinalat = sin(alat(ilat) * pio180)
     cosalat = cos(alat(ilat) * pio180)

     zea = erad * sinalat

     do ilon = 1, nlon

        xea = cosalat * er_cosalon(ilon)
        yea = cosalat * er_sinalon(ilon)

        img(ilon) = 1
        ntri(ilon) = 0

        mloop: do jm = 1, jtab_m(jtm_vadj)%jend(1)
                  im = jtab_m(jtm_vadj)%im(jm)

           ! Skip if latitudes are far from current M point

           if ( (glatm(im) > alat(ilat) + dllolam(im)) .or. &
                (glatm(im) < alat(ilat) - dllolam(im)) ) cycle

           ! Skip if longitudes are far from current M point,
           ! being careful at the left or right cyclical boundary

           deglon = dllolam(im) * fact

           if (alon360(ilon) > 90.0 .and. alon360(ilon) < 270.0) then

              if ( (glonm360(im) > alon360(ilon) + deglon) .or. &
                   (glonm360(im) < alon360(ilon) - deglon) ) cycle

           else

              if ( (glonm(im) > alon(ilon) + deglon) .or. &
                   (glonm(im) < alon(ilon) - deglon) ) cycle

           endif

           ! Transform current lat/lon point to PS coordinates tangent at M point

           dxe = xea - xem(im)
           dye = yea - yem(im)
           dze = zea - zem(im)

           call de_ps(dxe,dye,dze,cosmlat(im),sinmlat(im),cosmlon(im),sinmlon(im),qx,qy)

           ! Transform surrounding W points to PS coordinates tangent at M point

           do j = 1, 3
              iw = itab_m(im)%iw(j)

              dxe = xew(iw) - xem(im)
              dye = yew(iw) - yem(im)
              dze = zew(iw) - zem(im)

              call de_ps(dxe,dye,dze,cosmlat(im),sinmlat(im),cosmlon(im),sinmlon(im),xw(j),yw(j))
           enddo
           
           ! Determine if the current lat/lon point is in/near one of the 
           ! 3 triangles adjactent to this M point

           v2x = qx !- x0
           v2y = qy !- yo

           do j = 1, 3
              jnext = j + 1
              if (j == 3) jnext = 1

              ! Set up some triangle-check coefficients

              v0x = xw(j) !- xo
              v0y = yw(j) !- yo

              v1x = xw(jnext) !- xo
              v1y = yw(jnext) !- yo

              dot00 = v0x * v0x + v0y * v0y
              dot01 = v0x * v1x + v0y * v1y
              dot11 = v1x * v1x + v1y * v1y

              denomi = 1. / (dot00 * dot11 - dot01 * dot01)              

              ! Set up remaining triangle_check coefficients

              dot02 = v0x * v2x + v0y * v2y
              dot12 = v1x * v2x + v1y * v2y

              u = (dot11 * dot02 - dot01 * dot12) * denomi
              v = (dot00 * dot12 - dot01 * dot02) * denomi

              if (u > -fuzz .and. v > -fuzz .and. u + v < 1.0+fuzz) then
                 img(ilon) = itab_m(im)%imglobe
                 ntri(ilon) = j
                 exit mloop
              endif

           enddo  ! j loop

        enddo mloop

     enddo  ! ilon loop

     if (iparallel == 0) then

        ims(1:nlon,ilat) = img(1:nlon)
        ijs(1:nlon,ilat) = ntri(1:nlon)

     else
         
#ifdef OLAM_MPI
        call MPI_Allgather(img,  nlon, MPI_INTEGER, imgs,  nlon, MPI_INTEGER, MPI_COMM_WORLD, ier)
        call MPI_Allgather(ntri, nlon, MPI_INTEGER, ntris, nlon, MPI_INTEGER, MPI_COMM_WORLD, ier)
#endif

        ! In case of multiple matches in parallel, select the cell with the
        ! lowest global index to match the single-processor result

        do ilon = 1, nlon
           np = minloc(imgs(ilon,:), mask=imgs(ilon,:)>1, dim=1)
           ims(ilon,ilat) = imgs(ilon,np)
           ijs(ilon,ilat) = ntris(ilon,np)
        enddo
        
     endif

  enddo  ! ilat loop

end subroutine find_closest_m_ll
