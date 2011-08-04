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
subroutine vectslab_horiz_v(iplt)

use oplot_coms,  only: op
use mem_grid,    only: mza, mva, mwa, zm, lpw, unx, uny, unz, vnx, vny, vnz, &
                       xeu, yeu, zeu, xev, yev, zev
use mem_ijtabs,  only: itab_m, itab_u, itab_v, itab_w
use misc_coms,   only: io6, mdomain, meshtype
use consts_coms, only: erad

implicit none

integer, intent(in) :: iplt

integer :: iv,iw1,iw2,notavail,im1,im2,iw1_v1,iw2_v1

real :: fldval_u,fldval_v,uavg,vavg,pointx  &
   ,pointy,tailx,taily,stemangle,headangle1,headangle2   &
   ,headlen,head1x  &
   ,head1y,head2x,head2y,ptx,pty,tlx,tly,h1x,h1y,h2x,h2y  &
   ,fac,headlent  &
   ,tailxe,tailye,tailze,stemlen
real :: stemx,stemy,stemz,snx,sny,snz,rnx,rny,rnz
real :: head1xe,head1ye,head1ze,head2xe,head2ye,head2ze
real :: xeuv,yeuv,zeuv

integer :: ktf(mwa),kv(mva)
real :: wtbot(mva),wttop(mva)

! Wind barbs

real, parameter :: pf = .3  ! full tic length as a fraction of shaft length
real, parameter :: pi = .15 ! distance between tics as a fraction of shaft length
real, parameter :: ba = 50. ! value for which triangle is drawn
real, parameter :: bb = 10. ! value for which tic is drawn - half tic for half of bb

real :: speed, pc, xt, xea, yea, zea, xeb, yeb, zeb, xec, yec, zec
real :: xa, ya, xb, yb, xc, yc

! Find cell K indices on the given plot surface

op%stagpt = 'V'
call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

do iv = 2,mva

   if (meshtype == 1) then

      im1 = itab_u(iv)%im(1)
      im2 = itab_u(iv)%im(2)

      iw1 = itab_u(iv)%iw(1)
      iw2 = itab_u(iv)%iw(2)

      iw1_v1 = itab_w(iw1)%iu(1)
      iw2_v1 = itab_w(iw2)%iu(1)
      
      xeuv = xeu(iv)
      yeuv = yeu(iv)
      zeuv = zeu(iv)

! Skip this point if we want to plot vectors only on a coarser mesh level

      if (itab_w(iw1)%mrlw_orig > op%vec_maxmrl .and. &
          itab_w(iw2)%mrlw_orig > op%vec_maxmrl) cycle

      if (itab_w(iw1)%mrlw_orig > op%vec_maxmrl .and. &
          itab_w(iw2)%iu(1) /= iv) cycle

      if (itab_w(iw2)%mrlw_orig > op%vec_maxmrl .and. &
          itab_w(iw1)%iu(1) /= iv) cycle

   else

      im1 = itab_v(iv)%im(1)
      im2 = itab_v(iv)%im(2)

      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

      iw1_v1 = itab_w(iw1)%iv(1)
      iw2_v1 = itab_w(iw2)%iv(1)
      
      xeuv = xev(iv)
      yeuv = yev(iv)
      zeuv = zev(iv)

! Skip this point if we want to plot vectors only on a coarser mesh level

      if (itab_m(im1)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im2)%mrlm_orig > op%vec_maxmrl) cycle

      if (itab_m(im1)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im2)%iv(1) /= iv) cycle

      if (itab_m(im2)%mrlm_orig > op%vec_maxmrl .and. &
          itab_m(im1)%iv(1) /= iv) cycle

   endif

! Transform IV coordinates

   call oplot_transform(iplt,xeuv,yeuv,zeuv,pointx,pointy)

! Jump out of loop if vector head is outside plot window. 

   if (pointx < op%xmin .or. pointx > op%xmax .or.  &
       pointy < op%ymin .or. pointy > op%ymax) cycle

! Check if both neighboring W cells are above ground
! (Is this what we want to do here?)

   if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

! Cell is above ground 

      call oplot_lib(kv(iv),iv,'VALUV','UC',wtbot(iv),wttop(iv), &
                     fldval_u,notavail)
      call oplot_lib(kv(iv),iv,'VALUV','VC',wtbot(iv),wttop(iv), &
                     fldval_v,notavail)

! 3D vector displacement (in time interval op%dtvec)...

      if (op%vectbarb(iplt) == 'U') then

! Plot normal component to U face

         stemx = unx(iv) * fldval_u * op%dtvec
         stemy = uny(iv) * fldval_u * op%dtvec
         stemz = unz(iv) * fldval_u * op%dtvec

      elseif (op%vectbarb(iplt) == 'V') then

! Plot normal component to V face

         stemx = vnx(iv) * fldval_v * op%dtvec
         stemy = vny(iv) * fldval_v * op%dtvec
         stemz = vnz(iv) * fldval_v * op%dtvec

      elseif (op%vectbarb(iplt) == 'v') then

! Plot total horizontal vector at V point

         stemx = (unx(iv) * fldval_u + vnx(iv) * fldval_v) * op%dtvec
         stemy = (uny(iv) * fldval_u + vny(iv) * fldval_v) * op%dtvec
         stemz = (unz(iv) * fldval_u + vnz(iv) * fldval_v) * op%dtvec

      else    ! case for op%vectbarb(iplt) == 'B'

! Plot wind barb at V point

         speed = sqrt(fldval_u**2 + fldval_v**2)
         
         if (speed < 1.e-9) cycle

         stemx = (unx(iv) * fldval_u + vnx(iv) * fldval_v) * op%stemlength / speed
         stemy = (uny(iv) * fldval_u + vny(iv) * fldval_v) * op%stemlength / speed
         stemz = (unz(iv) * fldval_u + vnz(iv) * fldval_v) * op%stemlength / speed

      endif

! Vector length and unit components

      stemlen = max(1.e-6,sqrt(stemx**2 + stemy**2 + stemz**2))      

      snx = stemx / stemlen
      sny = stemy / stemlen
      snz = stemz / stemlen

! "Right" unit components

      if (mdomain <= 1) then  ! Spherical geometry case
         rnx = (sny * zeuv - snz * yeuv) / erad
         rny = (snz * xeuv - snx * zeuv) / erad
         rnz = (snx * yeuv - sny * xeuv) / erad
      else                    ! Cartesian case
         rnx = sny
         rny = - snx
         rnz = 0.
      endif

! Earth coordinates of tail

      tailxe = xeuv - stemx
      tailye = yeuv - stemy
      tailze = zeuv - stemz

      if (op%vectbarb(iplt) /= 'B') then

! Earth coordinates of left and right head tips

         headlen = op%headspeed * op%dtvec

         head1xe = xeuv + rnx * .42 * headlen - snx * .91 * headlen
         head1ye = yeuv + rny * .42 * headlen - sny * .91 * headlen
         head1ze = zeuv + rnz * .42 * headlen - snz * .91 * headlen

         head2xe = xeuv - rnx * .42 * headlen - snx * .91 * headlen
         head2ye = yeuv - rny * .42 * headlen - sny * .91 * headlen
         head2ze = zeuv - rnz * .42 * headlen - snz * .91 * headlen

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

! Compute re-scaling factor from transformation- STILL NEED THIS??????????

         fac = 1.0

! Draw vector

         call o_frstpt(pointx+fac*(tailx-pointx),pointy+fac*(taily-pointy))
         call o_vector(pointx,pointy)
         call o_frstpt(pointx+fac*(head1x-pointx),pointy+fac*(head1y-pointy))
         call o_vector(pointx,pointy)
         call o_vector(pointx+fac*(head2x-pointx),pointy+fac*(head2y-pointy))

      else    ! case for op%vectbarb(iplt) == 'B'

! Transform stem tail

         call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)

! Draw stem

         call o_frstpt(tailx,taily)
         call o_vector(pointx,pointy)

         pc = 1.
         xt = speed + .25 * bb

! Draw triangles (if any)

         do while (xt >= ba)

            xea = xeuv + (tailxe - xeuv) * pc
            yea = yeuv + (tailye - yeuv) * pc
            zea = zeuv + (tailze - zeuv) * pc

            xeb = xea - rnx * pf * op%stemlength
            yeb = yea - rny * pf * op%stemlength
            zeb = zea - rnz * pf * op%stemlength

            pc = pc - pi

            xec = xeuv + (tailxe - xeuv) * pc
            yec = yeuv + (tailye - yeuv) * pc
            zec = zeuv + (tailze - zeuv) * pc

! Transform triangle coordinates

            call oplot_transform(iplt,xea,yea,zea,xa,ya)
            call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)
            call oplot_transform(iplt,xec,yec,zec,xc,yc)

! Avoid wrap-around

            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)
            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xc)

            call frstpt(xa,ya)
            call vector(xb,yb)
            call vector(xc,yc)

            pc = pc - pi
            xt = xt - ba

         enddo

! Draw barbs (if any)

         do while (xt >= bb)

            xea = xeuv + (tailxe - xeuv) * pc
            yea = yeuv + (tailye - yeuv) * pc
            zea = zeuv + (tailze - zeuv) * pc

            xeb = xea - rnx * pf * op%stemlength
            yeb = yea - rny * pf * op%stemlength
            zeb = zea - rnz * pf * op%stemlength

! Transform triangle coordinates

            call oplot_transform(iplt,xea,yea,zea,xa,ya)
            call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)

! Avoid wrap-around

            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

            call frstpt(xa,ya)
            call vector(xb,yb)

            pc = pc - pi
            xt = xt - bb

         enddo

! Draw half barb (if any)

         if (xt >= .5*bb) then

            xea = xeuv + (tailxe - xeuv) * pc
            yea = yeuv + (tailye - yeuv) * pc
            zea = zeuv + (tailze - zeuv) * pc

            xeb = xea - rnx * .5 * pf * op%stemlength
            yeb = yea - rny * .5 * pf * op%stemlength
            zeb = zea - rnz * .5 * pf * op%stemlength

! Transform triangle coordinates

            call oplot_transform(iplt,xea,yea,zea,xa,ya)
            call oplot_transform(iplt,xeb,yeb,zeb,xb,yb)

! Avoid wrap-around

            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xa)
            if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,xb)

            call frstpt(xa,ya)
            call vector(xb,yb)

         endif

      endif

   endif

enddo

return
end subroutine vectslab_horiz_v
