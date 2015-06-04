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

subroutine vectslab_horiz_v(iplt)

  use oplot_coms,  only: op
  use mem_grid,    only: mza, mva, mwa, zm, lpw, unx, uny, unz, vnx, vny, vnz, &
                         xev, yev, zev
  use mem_ijtabs,  only: itab_m, itab_v, itab_w, jtab_v, jtv_wadj
  use consts_coms, only: eradi
  use misc_coms,   only: io6, mdomain, iparallel
  use mem_para,    only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: jv,iv,iw1,iw2,notavail,im1,im2,iw1_v1,iw2_v1,jvmax

  real :: fldval_v,pointx,pointy,tailx,taily
  real :: headlen,head1x,head1y,head2x,head2y
  real :: tailxe,tailye,tailze,stemlen
  real :: stemx,stemy,stemz,snx,sny,snz,rnx,rny,rnz
  real :: head1xe,head1ye,head1ze,head2xe,head2ye,head2ze

  integer :: ktf(mwa),kv(mva)
  real :: wtbot(mva),wttop(mva)

  ! Wind barbs (not used for V plotting)
  ! real, parameter :: pf = .3  ! full tic length as a fraction of shaft length
  ! real, parameter :: pi = .15 ! distance between tics as a fraction of shaft length
  ! real, parameter :: ba = 50. ! value for which triangle is drawn
  ! real, parameter :: bb = 10. ! value for which tic is drawn - half tic for half of bb

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, j, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  op%stagpt = 'V'
  call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

  jvmax = jtab_v(jtv_wadj)%jend(1)

  nu   = 0
  ipos = 0

  base = 8 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(jvmax) / 5. )
  else
     inc = jvmax
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  do jv = 1, jvmax

     iv  = jtab_v(jtv_wadj)%iv(jv)

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     iw1_v1 = itab_w(iw1)%iv(1)
     iw2_v1 = itab_w(iw2)%iv(1)
      
     ! Skip this point if we want to plot vectors only on a coarser mesh level

     if ( itab_m(im1)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im2)%mrlm_orig > op%vec_maxmrl ) cycle

     if ( itab_m(im1)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im2)%iv(1) /= iv) cycle

     if ( itab_m(im2)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im1)%iv(1) /= iv) cycle

     ! Transform IV coordinates

     call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),pointx,pointy)

     ! Jump out of loop if vector head is outside plot window. 

     if ( pointx < op%xmin .or. pointx > op%xmax .or.  &
          pointy < op%ymin .or. pointy > op%ymax ) cycle

     ! Check if both neighboring W cells are above ground
     ! (Is this what we want to do here?)

     if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

        ! Cell is above ground 

        call oplot_lib(kv(iv),iv,'VALUV','VC',wtbot(iv),wttop(iv), &
                       fldval_v,notavail)

        ! 3D vector displacement (in time interval op%dtvec)...
        ! Plot normal component to V face

        stemx = vnx(iv) * fldval_v * op%dtvec
        stemy = vny(iv) * fldval_v * op%dtvec
        stemz = vnz(iv) * fldval_v * op%dtvec

        ! Vector length and unit components

        stemlen = max(1.e-6,sqrt(stemx**2 + stemy**2 + stemz**2))      

        snx = stemx / stemlen
        sny = stemy / stemlen
        snz = stemz / stemlen

        ! "Right" unit components

        if (mdomain <= 1) then  ! Spherical geometry case
           rnx = (sny * zev(iv) - snz * yev(iv)) * eradi
           rny = (snz * xev(iv) - snx * zev(iv)) * eradi
           rnz = (snx * yev(iv) - sny * xev(iv)) * eradi
        else                    ! Cartesian case
           rnx = sny
           rny = - snx
           rnz = 0.
        endif

        ! Earth coordinates of tail

        tailxe = xev(iv) - stemx
        tailye = yev(iv) - stemy
        tailze = zev(iv) - stemz

        ! Earth coordinates of left and right head tips

        headlen = op%headspeed * op%dtvec

        head1xe = xev(iv) + rnx * .42 * headlen - snx * .91 * headlen
        head1ye = yev(iv) + rny * .42 * headlen - sny * .91 * headlen
        head1ze = zev(iv) + rnz * .42 * headlen - snz * .91 * headlen

        head2xe = xev(iv) - rnx * .42 * headlen - snx * .91 * headlen
        head2ye = yev(iv) - rny * .42 * headlen - sny * .91 * headlen
        head2ze = zev(iv) - rnz * .42 * headlen - snz * .91 * headlen

        ! Transform other tail and coordinates

        call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)
        call oplot_transform(iplt,head1xe,head1ye,head1ze,head1x,head1y)
        call oplot_transform(iplt,head2xe,head2ye,head2ze,head2x,head2y)

        ! Avoid wrap-around

        if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,tailx)
        if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head1x)
        if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head2x)

        ! Jump out of loop if tail or sides of head are outside plot window. 

        if ( tailx  < op%xmin .or. tailx  > op%xmax .or.  &
             taily  < op%ymin .or. taily  > op%ymax .or.  &
             head1x < op%xmin .or. head1x > op%xmax .or.  &
             head1y < op%ymin .or. head1y > op%ymax .or.  &
             head2x < op%xmin .or. head2x > op%xmax .or.  &
             head2y < op%ymin .or. head2y > op%ymax ) cycle

        ! Draw vector

        if (myrank == 0) then

           call o_frstpt(tailx,taily)
           call o_vector(pointx,pointy)
           call o_frstpt(head1x,head1y)
           call o_vector(pointx,pointy)
           call o_vector(head2x,head2y)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(tailx,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(taily,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(pointx, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(pointy, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head1x, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head1y, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head2x, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head2y, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif

  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)
               
                 call MPI_Unpack(buffer, buffsize, ipos, tailx,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, taily,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, pointx, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, pointy, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, head1x, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, head1y, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, head2x, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, head2y, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call o_frstpt(tailx,taily)
                 call o_vector(pointx,pointy)
                 call o_frstpt(head1x,head1y)
                 call o_vector(pointx,pointy)
                 call o_vector(head2x,head2y)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine vectslab_horiz_v

!===============================================================================

subroutine vectslab_horiz_w(iplt)

  use oplot_coms,  only: op
  use mem_grid,    only: mza, mwa, zm, lpw, xew, yew, zew, wnx, wny, wnz
  use mem_ijtabs,  only: itab_m, itab_v, itab_w, jtab_w, jtw_prog
  use mem_basic,   only: vxe, vye, vze
  use consts_coms, only: eradi
  use misc_coms,   only: io6, mdomain, iparallel
  use mem_para,    only: myrank, mgroupsize, nbytes_int, nbytes_real

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

  integer, parameter :: maxtris = 10 ! max wind 545 m/s
  integer, parameter :: maxbarb =  4
  integer, parameter :: maxhalf =  1
  integer            :: ntris, nbarb, nhalf

real :: xtris(3,maxtris), ytris(3,maxtris)
real :: xbarb(2,maxbarb), ybarb(2,maxbarb)
real :: xhalf(2,maxhalf), yhalf(2,maxhalf)

real :: speed, pc, xt, xea, yea, zea, xeb, yeb, zeb, xec, yec, zec
real :: xa, ya, xb, yb, xc, yc



  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, j, n, is
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

! Find cell K indices on the given plot surface

op%stagpt = 'T'
call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

  nu   = 0
  ipos = 0

  if (op%vectbarb(iplt) /= 'B') then
     base = 8 * nbytes_real
  else
     base = (4 + maxtris*6 + maxbarb*4 + maxhalf*4) * nbytes_real &
          + 3 * nbytes_int
  endif

  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mwa) / 5. )
  else
     inc = mwa
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

do jw = 1, jtab_w(jtw_prog)%jend(1)
   iw = jtab_w(jtw_prog)%iw(jw)

! Skip this point if it is underground

   if (ktf(iw) /= 0) cycle

! Skip this point if we want to plot vectors only on a coarser mesh level

   if (itab_w(iw)%mrlw_orig > op%vec_maxmrl) cycle

! Transform IV coordinates

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),pointx,pointy)

! Jump out of loop if vector head is outside plot window. 

   if (pointx < op%xmin .or. pointx > op%xmax .or.  &
       pointy < op%ymin .or. pointy > op%ymax) cycle

! 3D earth-coordinate wind components

   vx = wtbot(iw) * vxe(kw(iw),iw) + wttop(iw) * vxe(kw(iw)+1,iw)
   vy = wtbot(iw) * vye(kw(iw),iw) + wttop(iw) * vye(kw(iw)+1,iw)
   vz = wtbot(iw) * vze(kw(iw),iw) + wttop(iw) * vze(kw(iw)+1,iw)

! Projection onto local vertical vector

   vert_vel = vx * wnx(iw) + vy * wny(iw) + vz * wnz(iw)

! Subtract local vertical wind vector to give horizontal wind

   vx = vx - vert_vel * wnx(iw)
   vy = vy - vert_vel * wny(iw)
   vz = vz - vert_vel * wnz(iw)

   speed = sqrt(vx**2 + vy**2 + vz**2)

   if (speed < 1.e-9) cycle

! 3D vector displacement (in time interval op%dtvec)...

   if (op%vectbarb(iplt) == 'w') then

! Plot total horizontal vector at W point

      stemx = vx * op%dtvec
      stemy = vy * op%dtvec
      stemz = vz * op%dtvec

   else    ! case for op%vectbarb(iplt) == 'B'

! Plot wind barb at V point

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
      rnx = (sny * zew(iw) - snz * yew(iw)) * eradi
      rny = (snz * xew(iw) - snx * zew(iw)) * eradi
      rnz = (snx * yew(iw) - sny * xew(iw)) * eradi
   else                    ! Cartesian case
      rnx = sny
      rny = - snx
      rnz = 0.
   endif

! Earth coordinates of tail

   tailxe = xew(iw) - stemx
   tailye = yew(iw) - stemy
   tailze = zew(iw) - stemz

   if (op%vectbarb(iplt) /= 'B') then

! Earth coordinates of left and right head tips

      headlen = op%headspeed * op%dtvec

      head1xe = xew(iw) + rnx * .42 * headlen - snx * .91 * headlen
      head1ye = yew(iw) + rny * .42 * headlen - sny * .91 * headlen
      head1ze = zew(iw) + rnz * .42 * headlen - snz * .91 * headlen

      head2xe = xew(iw) - rnx * .42 * headlen - snx * .91 * headlen
      head2ye = yew(iw) - rny * .42 * headlen - sny * .91 * headlen
      head2ze = zew(iw) - rnz * .42 * headlen - snz * .91 * headlen

! Transform other tail and coordinates

      call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)
      call oplot_transform(iplt,head1xe,head1ye,head1ze,head1x,head1y)
      call oplot_transform(iplt,head2xe,head2ye,head2ze,head2x,head2y)

! Avoid wrap-around

      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,tailx)
      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head1x)
      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head2x)

! Jump out of loop if tail or sides of head are outside plot window. 

      if (tailx  < op%xmin .or. tailx  > op%xmax .or.  &
          taily  < op%ymin .or. taily  > op%ymax .or.  &
          head1x < op%xmin .or. head1x > op%xmax .or.  &
          head1y < op%ymin .or. head1y > op%ymax .or.  &
          head2x < op%xmin .or. head2x > op%xmax .or.  &
          head2y < op%ymin .or. head2y > op%ymax) cycle

! Draw vector

      if (myrank == 0) then
         call o_frstpt(tailx,taily)
         call o_vector(pointx,pointy)
         call o_frstpt(head1x,head1y)
         call o_vector(pointx,pointy)
         call o_vector(head2x,head2y)
      else
#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif
           call MPI_Pack(tailx,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(taily,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(pointx, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(pointy, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head1x, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head1y, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head2x, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(head2y, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
        endif

   else    ! case for op%vectbarb(iplt) == 'B'

! Transform stem tail

      call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)

! Draw stem

      if (myrank == 0) then
         call o_frstpt(tailx,taily)
         call o_vector(pointx,pointy)
      endif

      pc = 1.
      xt = speed + .25 * bb

! Draw triangles (if any)

      ntris = 0
      do while (xt >= ba)

         xea = xew(iw) + (tailxe - xew(iw)) * pc
         yea = yew(iw) + (tailye - yew(iw)) * pc
         zea = zew(iw) + (tailze - zew(iw)) * pc

         xeb = xea - rnx * pf * op%stemlength
         yeb = yea - rny * pf * op%stemlength
         zeb = zea - rnz * pf * op%stemlength

         pc = pc - pi

         xec = xew(iw) + (tailxe - xew(iw)) * pc
         yec = yew(iw) + (tailye - yew(iw)) * pc
         zec = zew(iw) + (tailze - zew(iw)) * pc

! Transform triangle coordinates

         call oplot_transform(iplt,xea,yea,zea,xa,ya)
         call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)
         call oplot_transform(iplt,xec,yec,zec,xc,yc)

! Avoid wrap-around

         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)
         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xc)

         if (myrank == 0) then
            call o_frstpt(xa,ya)
            call o_vector(xb,yb)
            call o_vector(xc,yc)
         elseif (ntris < maxtris) then
            ntris = ntris + 1
            xtris(1,ntris) = xa
            ytris(1,ntris) = ya
            xtris(2,ntris) = xb
            ytris(2,ntris) = yb
            xtris(3,ntris) = xc
            ytris(3,ntris) = yc
         endif

         pc = pc - pi
         xt = xt - ba

      enddo

! Draw barbs (if any)

      nbarb = 0
      do while (xt >= bb)

         xea = xew(iw) + (tailxe - xew(iw)) * pc
         yea = yew(iw) + (tailye - yew(iw)) * pc
         zea = zew(iw) + (tailze - zew(iw)) * pc

         xeb = xea - rnx * pf * op%stemlength
         yeb = yea - rny * pf * op%stemlength
         zeb = zea - rnz * pf * op%stemlength

! Transform triangle coordinates

         call oplot_transform(iplt,xea,yea,zea,xa,ya)
         call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)

! Avoid wrap-around

         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

         if (myrank == 0) then
            call o_frstpt(xa,ya)
            call o_vector(xb,yb)
         elseif (nbarb < maxbarb) then
            nbarb = nbarb + 1
            xbarb(1,nbarb) = xa
            ybarb(1,nbarb) = ya
            xbarb(2,nbarb) = xb
            ybarb(2,nbarb) = yb
         endif

         pc = pc - pi
         xt = xt - bb

      enddo

! Draw half barb (if any)

      nhalf = 0
      if (xt >= .5*bb) then

         xea = xew(iw) + (tailxe - xew(iw)) * pc
         yea = yew(iw) + (tailye - yew(iw)) * pc
         zea = zew(iw) + (tailze - zew(iw)) * pc

         xeb = xea - rnx * .5 * pf * op%stemlength
         yeb = yea - rny * .5 * pf * op%stemlength
         zeb = zea - rnz * .5 * pf * op%stemlength

! Transform triangle coordinates

         call oplot_transform(iplt,xea,yea,zea,xa,ya)
         call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)

! Avoid wrap-around

         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
         if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

         if (myrank == 0) then
            call o_frstpt(xa,ya)
            call o_vector(xb,yb)
         else
            nhalf = nhalf + 1
            xhalf(1,nhalf) = xa
            yhalf(1,nhalf) = ya
            xhalf(2,nhalf) = xb
            yhalf(2,nhalf) = yb
         endif

      endif

#ifdef OLAM_MPI
      if (myrank /= 0) then

         nu = nu + 1
         if (buffsize < ipos + base) then
            allocate( bcopy (buffsize + inc * base) )
            bcopy(1:buffsize) = buffer
            call move_alloc(bcopy, buffer)
            buffsize = size(buffer)
         endif

         call MPI_Pack(tailx,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         call MPI_Pack(taily,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         call MPI_Pack(pointx, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         call MPI_Pack(pointy, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

         call MPI_Pack(ntris, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         call MPI_Pack(nbarb, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         call MPI_Pack(nhalf, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

         if (ntris > 0) then
            call MPI_Pack(xtris, 3*ntris, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
            call MPI_Pack(ytris, 3*ntris, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         endif

         if (nbarb > 0) then
            call MPI_Pack(xbarb, 2*nbarb, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
            call MPI_Pack(ybarb, 2*nbarb, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         endif

         if (nhalf > 0) then
            call MPI_Pack(xhalf, 2*nhalf, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
            call MPI_Pack(yhalf, 2*nhalf, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
         endif
      endif
#endif

   endif

enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)

                 if (op%vectbarb(iplt) /= 'B') then

                    call MPI_Unpack(buffer, buffsize, ipos, tailx,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, taily,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, pointx, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, pointy, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, head1x, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, head1y, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, head2x, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, head2y, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                    call o_frstpt(tailx,taily)
                    call o_vector(pointx,pointy)
                    call o_frstpt(head1x,head1y)
                    call o_vector(pointx,pointy)
                    call o_vector(head2x,head2y)

                 else

                    call MPI_Unpack(buffer, buffsize, ipos, tailx,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, taily,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, pointx, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, pointy, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                    call o_frstpt(tailx,taily)
                    call o_vector(pointx,pointy)

                    call MPI_Unpack(buffer, buffsize, ipos, ntris, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, nbarb, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, nhalf, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)

                    if (ntris > 0) then

                       call MPI_Unpack(buffer, buffsize, ipos, xtris, ntris*3, MPI_REAL, MPI_COMM_WORLD, ier)
                       call MPI_Unpack(buffer, buffsize, ipos, ytris, ntris*3, MPI_REAL, MPI_COMM_WORLD, ier)

                       do is = 1, ntris
                          call o_frstpt(xtris(1,is),ytris(1,is))
                          call o_vector(xtris(2,is),ytris(2,is))
                          call o_vector(xtris(3,is),ytris(3,is))
                       enddo

                    endif

                    if (nbarb > 0) then

                       call MPI_Unpack(buffer, buffsize, ipos, xbarb, nbarb*2, MPI_REAL, MPI_COMM_WORLD, ier)
                       call MPI_Unpack(buffer, buffsize, ipos, ybarb, nbarb*2, MPI_REAL, MPI_COMM_WORLD, ier)

                       do is = 1, nbarb
                          call o_frstpt(xbarb(1,is),ybarb(1,is))
                          call o_vector(xbarb(2,is),ybarb(2,is))
                       enddo

                    endif

                    if (nhalf > 0) then

                       call MPI_Unpack(buffer, buffsize, ipos, xhalf, nhalf*2, MPI_REAL, MPI_COMM_WORLD, ier)
                       call MPI_Unpack(buffer, buffsize, ipos, yhalf, nhalf*2, MPI_REAL, MPI_COMM_WORLD, ier)

                       do is = 1, nhalf
                          call o_frstpt(xhalf(1,is),yhalf(1,is))
                          call o_vector(xhalf(2,is),yhalf(2,is))
                       enddo

                    endif

                 endif

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine vectslab_horiz_w
