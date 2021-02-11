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
subroutine spring_dynamics(mrows, moveint, ngr, nma, nua, nwa, xem, yem, zem, &
                           itab_md, itab_ud, itab_wd)

  ! Subroutine spring_dynamics is used only for mesh refinements (on both the
  ! atm and sfc grids).  Call subroutine spring_dynamics1 to adjust global grid 1.

  use mem_ijtabs,  only: itab_md_vars, itab_ud_vars, itab_wd_vars
  use mem_grid,    only: impent
  use consts_coms, only: pi2, erad, erador5, r8
  use misc_coms,   only: io6, nxp, mdomain, deltax
  use oplot_coms,  only: op

  implicit none

  integer, intent(in) :: mrows, moveint, ngr, nma, nua, nwa

  real, intent(inout) :: xem(nma), yem(nma), zem(nma)

  type (itab_md_vars), intent(inout) :: itab_md(nma)
  type (itab_ud_vars), intent(inout) :: itab_ud(nua)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwa)

  integer, parameter :: niter = 5000
  integer, parameter :: nprnt = 50
  real,    parameter :: relax = .05
  real,    parameter :: beta  = 1.24
  real,    parameter :: onethird = 1./3.

  ! Automatic arrays

  real :: dist (nua)
  real :: dist0(nua)

  real :: dx(nua)
  real :: dy(nua)
  real :: dz(nua)

  real :: dirs(7,nma)
  real :: dsm(nma)

  integer :: iu,iu1,iu2,iu3,iu4
  integer :: im,im1,im2,im3,im4
  integer :: iw,iw1,iw2,j
  integer :: iter,ipent,mrow1,mrow2,npoly,ngrw,mrmax,mrmin

  integer :: iskip
  real :: xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2

  real :: twocosphi3, twocosphi4, ratio

  real :: dist00, distm, frac_change

  integer :: ipic

  ! Array MOVEM is used to flag selected MD points that will be allowed to move
  ! in the spring dynamics procedure.  If MOVEM = .true. the point is allowed to
  ! move, and if MOVEM = 0 the point remains stationary.

  ! Starting with OLAM-SOIL and carrying over to merged OLAM version 6.0, MD
  ! points OUTSIDE the current newly refined mesh are NEVER allowed to move.
  ! This is because ATM and SURFACE grids can be independently refined, and
  ! neither must be allowed to impact the coarser grids that are common to
  ! both ATM and SURFACE.

  logical  :: movem(nma)
  logical  :: moveu(nua)
  logical  :: compu(nua)

  real(r8) :: xem8(nma),yem8(nma),zem8(nma)
  real(r8) :: xem0(nma),yem0(nma),zem0(nma)
  real(r8) :: expansion, erad8

  character(10) :: string

  erad8 = real(erad,r8)

  dsm(:) = 0.0

  xem8(:) = real(xem(:),r8)
  yem8(:) = real(yem(:),r8)
  zem8(:) = real(zem(:),r8)

  ! First, turn off movement of ALL M points

  movem(:) = .false.

  do im = 2,nma

     ! Check grid number of M point.  M point will move in spring dynamics only
     ! if it is of the current grid number (ngr) that was just spawned.  On the
     ! atmosphere grid, the row containing 7-way vertices (future centers of
     ! heptagon Voronoi cells) is the outermost row of M points that are
     ! designated to be on ngr.  On the surface grid, however, the 7-way
     ! vertices are NOT on ngr but rather on the lower-numbered parent grid.
     ! This prevents these M points from moving in spring dynamics.

     if (itab_md(im)%ngr /= ngr) cycle

     ! If interior M points are flagged to move, then all M points on ngr will
     ! be moved, so set movem flag to 1.

     if (moveint == 1) then
        movem(im) = .true.
     else 

        ! Even if interior points are not to be moved, M points in the
        ! transition row must always be moveable, so set their movem flag to 1. 

        npoly = itab_md(im)%npoly

        do j = 1,npoly
           iw = itab_md(im)%iw(j)

           ngrw = itab_wd(iw)%ngr

           if (ngrw /= ngr) cycle

           if (itab_wd(iw)%mrow == -3 .or. &
               itab_wd(iw)%mrow == -2 .or. &
               itab_wd(iw)%mrow == -1 .or. &
               itab_wd(iw)%mrow ==  1) then
              movem(im) = .true.
              cycle ! once this point has been flagged, we can move on
           endif

        enddo

     endif

  enddo      

  ! Compute mean length of coarse mesh U segments

  if (mdomain < 2) then
     dist00 = beta * pi2 * erad / (5. * real(nxp))
  else
     dist00 = deltax * sqrt(2.0) / sqrt(sqrt(3.0))
  endif

  ! Compute target length of each triangle U segment
  ! Loop over all U points

  !$omp parallel private(iter)
  !$omp do private(iu1,iu3,im1,im2,im3,im4,iw1,iw2,mrow1,mrow2,mrmax,mrmin)
  do iu = 2, nua

     iu1 = itab_ud(iu)%iu(1)
     iu3 = itab_ud(iu)%iu(3)

     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     if (itab_ud(iu1)%im(1) == im1) then
        im3 = itab_ud(iu1)%im(2)
     else
        im3 = itab_ud(iu1)%im(1)
     endif

     if (itab_ud(iu3)%im(1) == im1) then
        im4 = itab_ud(iu3)%im(2)
     else
        im4 = itab_ud(iu3)%im(1)
     endif

     moveu(iu) = movem(im1) .or. movem(im2)
     compu(iu) = movem(im1) .or. movem(im2) .or. movem(im3) .or. movem(im4)

     if (.not. compu(iu)) cycle

     ! Compute target distance for any MRL value

     dist0(iu) = dist00 / real( 2**(itab_ud(iu)%mrlu - 1) )

     ! Modified distance in MRL border zone

     if (ngr > 1) then

        iw1 = itab_ud(iu)%iw(1)
        iw2 = itab_ud(iu)%iw(2)

        mrow1 = itab_wd(iw1)%mrow
        mrow2 = itab_wd(iw2)%mrow

        mrmax = max(mrow1,mrow2)
        mrmin = min(mrow1,mrow2)
      
        if (mrows >= 3) then

           if (.false.) then  ! Martin's adjustments to expand smallest cells

              if (mrmax == -5 .and. mrmin == -5) then
                 dist0(iu) = dist0(iu) * 12.12 / 12.
              elseif (mrmax == -4 .and. mrmin == -5) then
                 dist0(iu) = dist0(iu) * 12.25 / 12.
              elseif (mrmax == -4 .and. mrmin == -4) then
                 dist0(iu) = dist0(iu) * 12.65 / 12.
              elseif (mrmax == -3 .and. mrmin == -4) then
                 dist0(iu) = dist0(iu) * 12.75 / 12.
              elseif (mrmax == -3 .and. mrmin == -3) then
                 dist0(iu) = dist0(iu) * 13.25 / 12.
              elseif (mrmax == -2 .and. mrmin == -3) then
                 dist0(iu) = dist0(iu) * 13.40 / 12.
              elseif (mrmax == -2 .and. mrmin == -2) then
                dist0(iu) = dist0(iu) * 13.80 / 12.
              elseif (mrmax == -1 .and. mrmin == -2) then
                 dist0(iu) = dist0(iu) * 14.00 / 12.
              elseif (mrmax == -1 .and. mrmin == -1) then
                 dist0(iu) = dist0(iu) * 15.60 / 12.
              elseif (mrmax == 1 .and. mrmin == -1) then
                 dist0(iu) = dist0(iu) * 15.80 / 12.
              elseif (mrmax == 1 .and. mrmin == 1) then
                 dist0(iu) = dist0(iu) * 10.35 / 12.
              elseif (mrmax == 1 .and. mrmin == 2) then
                 dist0(iu) = dist0(iu) * 11.10 / 12.
              elseif (mrmax == 2 .and. mrmin == 2) then
                 dist0(iu) = dist0(iu) * 11.25 / 12.
              elseif (mrmin >= 1) then
                 dist0(iu) = dist0(iu) * 11.45 / 12.
              endif

           else

              if     (mrmax == -2 .and. mrmin == -2) then
                 dist0(iu) = dist0(iu) *  7. / 6.  !* .90
              elseif (mrmax == -1 .and. mrmin == -2) then
                 dist0(iu) = dist0(iu) *  8. / 6.  !* .90
              elseif (mrmax == -1 .and. mrmin == -1) then
                 dist0(iu) = dist0(iu) *  9. / 6.  !* .90
              elseif (mrmax == 1 .and. mrmin == -1) then
                 dist0(iu) = dist0(iu) * 10. / 6.  !* .90
              elseif (mrmax == 1 .and. mrmin == 1) then
                 dist0(iu) = dist0(iu) * 11. / 12. !* .90
              endif

           endif

        else ! (mrows = 1) Can be used for sfc grid but not atm grid

           if (mrmax == 1 .and. mrmin == 1) then
              dist0(iu) = dist0(iu) * 0.75
           endif

        endif

     endif ! ngr

  enddo
  !$omp end do

  !$omp do private(j,iu) schedule(guided)
  do im = 2, nma
     if (movem(im)) then

        do j = 1, itab_md(im)%npoly
           iu = itab_md(im)%iu(j)

           if (itab_ud(iu)%im(2) == im) then
              dirs(j,im) =  relax
           else
              dirs(j,im) = -relax
           endif

        enddo

     endif
  enddo
  !$omp end do

  !$omp single
  write(*,'(4(A,I0))') "In spring dynamics: ngr = ", ngr, &
       ", nma = ", nma, ", nmovem = ", count(movem), ", niter = ", niter
  !$omp end single

  ! Main iteration loop

  do iter = 1, niter

     if (iter == 1 .or. mod(iter,nprnt) == 0) then
        !$omp sections
        !$omp section
        xem0 = xem8
        !$omp section
        yem0 = yem8
        !$omp section
        zem0 = zem8
        !$omp end sections
     endif

     ! Compute length of each U segment

     !$omp do private(im1,im2) schedule(guided)
     do iu = 2, nua

        if (compu(iu)) then
           im1 = itab_ud(iu)%im(1)
           im2 = itab_ud(iu)%im(2)

           dx(iu) = real( xem8(im2) - xem8(im1) )
           dy(iu) = real( yem8(im2) - yem8(im1) )
           dz(iu) = real( zem8(im2) - zem8(im1) )

           dist(iu) = sqrt( dx(iu) * dx(iu) &
                          + dy(iu) * dy(iu) &
                          + dz(iu) * dz(iu) )
        endif
     enddo
     !$omp end do

     ! Adjustment of dist0 based on opposite angles of triangles

     !$omp do private(iu1,iu2,iu3,iu4,twocosphi3,twocosphi4,&
     !$omp            ratio,distm,frac_change) schedule(guided)
     do iu = 2, nua
        if (moveu(iu)) then

           iu1 = itab_ud(iu)%iu(1)
           iu2 = itab_ud(iu)%iu(2)
           iu3 = itab_ud(iu)%iu(3)
           iu4 = itab_ud(iu)%iu(4)

           ! Compute two*cosine of angles at IM3 and IM4

           twocosphi3 = (dist(iu1)**2 + dist(iu2)**2 - dist(iu)**2) / (dist(iu1) * dist(iu2))
           twocosphi4 = (dist(iu3)**2 + dist(iu4)**2 - dist(iu)**2) / (dist(iu3) * dist(iu4))

           ! Decrease dist0 if smaller twocosine is less than limiting value of two*cos(72 deg)

           ratio = min(twocosphi3,twocosphi4)

           if (ratio < .618) then
              distm = dist0(iu) * max(ratio / .618, 0.01)               ! Bob's
             !distm = dist0(iu) * max(ratio / .618, 0.01) ** onethird   ! Martin's (why the onethird?)
           else
              distm = dist0(iu)
           endif

           ! Fractional change to dist that would make it equal dist0

           frac_change = (distm - dist(iu)) / dist(iu)

           ! Compute components of displacement that gives dist0

           dx(iu) = dx(iu) * frac_change
           dy(iu) = dy(iu) * frac_change
           dz(iu) = dz(iu) * frac_change
        endif

     enddo
     !$omp end do

     !$omp do private(j,iu,expansion) schedule(guided)
     do im = 2, nma

        ! For preventing either polar M point from moving:
        ! if (im == impent(1 )) cycle
        ! if (im == impent(12)) cycle

        ! For preventing all pentagonal points from moving:
        ! if (any(im == impent(1:12))) cycle

        ! Apply the displacement components to each M point

        if (movem(im)) then

           do j = 1, itab_md(im)%npoly
              iu = itab_md(im)%iu(j)
              xem8(im) = xem8(im) + dirs(j,im) * dx(iu)
              yem8(im) = yem8(im) + dirs(j,im) * dy(iu)
              zem8(im) = zem8(im) + dirs(j,im) * dz(iu)
           enddo

           ! Push M point coordinates out to earth radius

           if (mdomain < 2) then

              expansion = erad8 / sqrt( xem8(im) ** 2 + yem8(im) ** 2 + zem8(im) ** 2 )
              xem8(im) = xem8(im) * expansion
              yem8(im) = yem8(im) * expansion
              zem8(im) = zem8(im) * expansion

           endif

        endif

     enddo
     !$omp end do

     ! Print iteration status

     if (iter == 1 .or. mod(iter,nprnt) == 0) then

        !$omp do schedule(guided)
        do im = 2, nma
           if (movem(im)) then

              dsm(im) = sqrt( (xem8(im) - xem0(im))**2 &
                            + (yem8(im) - yem0(im))**2 &
                            + (zem8(im) - zem0(im))**2 )

           endif
        enddo
        !$omp end do

        !$omp single
        write(*,'(3x,A,I5,A,I5,A,f0.4)') &
             "Iteration ", iter, " of ", niter, ",  Max DS = ", maxval(dsm)
        !$omp end single

     endif

     ! Section for plotting grid at intermediate stages of spring dynamics adjustment
     ! Remove .false. to activate

     if (.false. .and. (iter == 1 .or. mod(iter,100) == 0)) then

        !$omp single
        xem(:) = real(xem8(:))
        yem(:) = real(yem8(:))
        zem(:) = real(zem8(:))

        ! Plot grid lines

        call o_reopnwk()

        do ipic = 1,3

        call plotback()
        call oplot_set(ipic)
 
        do iu = 2,nua
           im1 = itab_ud(iu)%im(1)
           im2 = itab_ud(iu)%im(2)

           call oplot_transform(1,xem(im1),yem(im1),zem(im1),xp1,yp1)
           call oplot_transform(1,xem(im2),yem(im2),zem(im2),xp2,yp2)

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           if (iskip == 1) cycle

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)
        enddo

        if (ipic == 1) then

           ! print mrow values

           do iw = 2, nwa
              im1 = itab_wd(iw)%im(1)
              im2 = itab_wd(iw)%im(2)
              im3 = itab_wd(iw)%im(3) 

              call oplot_transform(1, (xem(im1)+xem(im2)+xem(im3))/3., &
                                      (yem(im1)+yem(im2)+yem(im3))/3., &
                                      (zem(im1)+zem(im2)+zem(im3))/3., &
                                      xp1, yp1                         )

              if ( xp1 < op%xmin .or.  &
                   xp1 > op%xmax .or.  &
                   yp1 < op%ymin .or.  &
                   yp1 > op%ymax ) cycle

              write(string,'(I0)') itab_wd(iw)%mrow
              call o_plchlq (xp1,yp1,trim(adjustl(string)),0.003,0.,0.)
           enddo

        elseif (ipic == 2) then

           ! print ngr values for M points

           do im = 2, nma

              call oplot_transform(1, xem(im), yem(im), zem(im), xp1, yp1)

              if ( xp1 < op%xmin .or.  &
                   xp1 > op%xmax .or.  &
                   yp1 < op%ymin .or.  &
                   yp1 > op%ymax ) cycle

              write(string,'(I0)') itab_md(im)%ngr
              call o_plchlq (xp1,yp1,trim(adjustl(string)),0.003,0.,0.)
           enddo

        elseif (ipic == 3) then

           ! print im values for M points

           do im = 2, nma

              call oplot_transform(1, xem(im), yem(im), zem(im), xp1, yp1)

              if ( xp1 < op%xmin .or.  &
                   xp1 > op%xmax .or.  &
                   yp1 < op%ymin .or.  &
                   yp1 > op%ymax ) cycle

              write(string,'(I0)') im
              call o_plchlq (xp1,yp1,trim(adjustl(string)),0.003,0.,0.)
           enddo

        endif

        call o_frame()

        enddo

        call o_clswk()
        !$omp end single

     endif ! mod(iter,*)

  enddo ! iter
  !$omp end parallel

  xem(:) = real(xem8(:))
  yem(:) = real(yem8(:))
  zem(:) = real(zem8(:))

end subroutine spring_dynamics

!===============================================================================

subroutine spring_dynamics1()

  ! Subroutine spring_dynamics1 is used only for adjusting grid 1 (i.e.,
  ! the quasi-uniform global atm grid) prior to any mesh refinements.
  ! Call subroutine spring_dynamics to adjust mesh refinements of either
  ! the atm or surface grids.

  use mem_ijtabs,  only: itab_md, itab_ud, itab_wd
  use mem_grid,    only: nma, nua, nwa, xem, yem, zem, impent
  use consts_coms, only: pi2, erad, r8, piu180
  use misc_coms,   only: io6, nxp, mdomain, deltax
  use oplot_coms,  only: op

  implicit none

  integer, parameter :: ngr   = 1
  integer, parameter :: niter = 5000
  integer, parameter :: nprnt = 50
  real,    parameter :: relax = .08
  real,    parameter :: beta  = 1.24
  real,    parameter :: onethird = 1./3.

! Automatic arrays

  real     :: dist(nua), distm(nua)
  real     :: ratio, frac_change
  real     :: dx(nua), dy(nua), dz(nua)
  real     :: dirs(7,nma)

  integer  :: iu,iu1,iu2,iu3,iu4
  integer  :: im,im1,im2,im3
  integer  :: iw,j
  integer  :: iter

  integer  :: iskip
  real     :: xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2
  real     :: twocosphi3, twocosphi4
  real     :: dist00
  real     :: dsm(nma)

  real(r8) :: xem8(nma),yem8(nma),zem8(nma)
  real(r8) :: xem0(nma),yem0(nma),zem0(nma)
  real(r8) :: expansion, erad8

  character(10) :: string

! special
! RETURN
! end special

  erad8 = real(erad,r8)

  dsm(1) = 0.0

  xem8(:) = real(xem(:),r8)
  yem8(:) = real(yem(:),r8)
  zem8(:) = real(zem(:),r8)

! Compute mean length of coarse mesh U segments

  if (mdomain < 2) then
     dist00 = beta * pi2 * erad / (5. * real(nxp))
  else
     dist00 = deltax * sqrt( 2.0 / sqrt(3.0) )
  endif

  write(io6,'(a,4i9)') "In spring dynamics: ngr,nma,niter = ",1,nma,niter

  !$omp parallel private(iter)
  !$omp do private(j,iu)
  do im = 2, nma

     do j = 1, itab_md(im)%npoly
        iu = itab_md(im)%iu(j)

        if (itab_ud(iu)%im(2) == im) then
           dirs(j,im) =  relax
        else
           dirs(j,im) = -relax
        endif

     enddo
  enddo
  !$omp end do

! Main iteration loop

  do iter = 1, niter

     if (iter == 1 .or. mod(iter,nprnt) == 0) then
        !$omp sections
        !$omp section
        xem0 = xem8
        !$omp section
        yem0 = yem8
        !$omp section
        zem0 = zem8
        !$omp end sections
     endif

! Compute length of each U segment

     !$omp do private(im1,im2)
     do iu = 2, nua
        im1 = itab_ud(iu)%im(1)
        im2 = itab_ud(iu)%im(2)

        dx(iu) = real( xem8(im2) - xem8(im1) )
        dy(iu) = real( yem8(im2) - yem8(im1) )
        dz(iu) = real( zem8(im2) - zem8(im1) )
     enddo
     !$omp end do

     !$omp do
     do iu = 2, nua
        dist(iu) = sqrt( dx(iu) * dx(iu) &
                       + dy(iu) * dy(iu) &
                       + dz(iu) * dz(iu) )
     enddo
     !$omp end do

! Adjustment of dist0 based on opposite angles of triangles

     !$omp do private(iu1,iu2,iu3,iu4,twocosphi3,twocosphi4,ratio)
     do iu = 2, nua

        iu1 = itab_ud(iu)%iu(1)
        iu2 = itab_ud(iu)%iu(2)
        iu3 = itab_ud(iu)%iu(3)
        iu4 = itab_ud(iu)%iu(4)

        ! Compute cosine of angles at IM3 and IM4

        twocosphi3 = (dist(iu1)**2 + dist(iu2)**2 - dist(iu)**2) / (dist(iu1) * dist(iu2))
        twocosphi4 = (dist(iu3)**2 + dist(iu4)**2 - dist(iu)**2) / (dist(iu3) * dist(iu4))

        ! Ratio of smaller cosine to limiting value of cos(72 deg)

        ratio = min(twocosphi3,twocosphi4)

        if (ratio < .618) then
           distm(iu) = dist00 * max(ratio/.618,0.01) ** onethird
        else
           distm(iu) = dist00
        endif

     enddo
     !$omp end do

     !$omp do private(frac_change)
     do iu = 2, nua

        ! Fractional change to dist that would make it equal dist0

        frac_change = (distm(iu) - dist(iu)) / dist(iu)

        ! Compute components of displacement that gives dist0

        dx(iu) = dx(iu) * frac_change
        dy(iu) = dy(iu) * frac_change
        dz(iu) = dz(iu) * frac_change

     enddo
     !$omp end do

     !$omp do private(j,iu)
     do im = 2, nma

        ! For preventing either polar M point from moving:
        ! if (im == impent(1 )) cycle
        ! if (im == impent(12)) cycle

        ! For preventing all pentagonal points from moving:
        if (any(im == impent(1:12))) cycle

        ! Apply the displacement components to each M point
        do j = 1, itab_md(im)%npoly
           iu = itab_md(im)%iu(j)
           xem8(im) = xem8(im) + dirs(j,im) * dx(iu)
           yem8(im) = yem8(im) + dirs(j,im) * dy(iu)
           zem8(im) = zem8(im) + dirs(j,im) * dz(iu)
        enddo

     enddo
     !$omp end do

! Push M point coordinates out to earth radius

     if (mdomain < 2) then

        !$omp do private(expansion)
        do im = 2, nma
           expansion = erad8 / sqrt( xem8(im) ** 2 + yem8(im) ** 2 + zem8(im) ** 2 )
           xem8(im) = xem8(im) * expansion
           yem8(im) = yem8(im) * expansion
           zem8(im) = zem8(im) * expansion
        enddo
        !$omp end do

     endif

! Print iteration status

     if (iter == 1 .or. mod(iter,nprnt) == 0) then

        !$omp do
        do im = 2, nma
           dsm(im) = sqrt( (xem8(im) - xem0(im))**2 &
                         + (yem8(im) - yem0(im))**2 &
                         + (zem8(im) - zem0(im))**2 )
        enddo
        !$omp end do

        !$omp single
        write(*,'(3x,A,I5,A,I5,A,f0.4)') &
             "Iteration ", iter, " of ", niter, ",  Max DS = ", maxval(dsm)
        !$omp end single

     endif

! Section for plotting grid at intermediate stages of spring dynamics adjustment
! Remove .false. to activate

     if (.false. .and. (iter == 1 .or. mod(iter,100) == 0)) then

        !$omp single
        xem(:) = real(xem8(:))
        yem(:) = real(yem8(:))
        zem(:) = real(zem8(:))

        ! Plot grid lines

        call o_reopnwk()
        call plotback()

        call oplot_set(1)

        do iu = 2,nua
           im1 = itab_ud(iu)%im(1)
           im2 = itab_ud(iu)%im(2)

           call oplot_transform(1,xem(im1),yem(im1),zem(im1),xp1,yp1)
           call oplot_transform(1,xem(im2),yem(im2),zem(im2),xp2,yp2)

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           if (iskip == 1) cycle

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)
        enddo

        ! Print mrow values

        do iw = 2, nwa
           im1 = itab_wd(iw)%im(1)
           im2 = itab_wd(iw)%im(2)
           im3 = itab_wd(iw)%im(3)

           call oplot_transform(1, (xem(im1)+xem(im2)+xem(im3))/3., &
                                   (yem(im1)+yem(im2)+yem(im3))/3., &
                                   (zem(im1)+zem(im2)+zem(im3))/3., &
                                   xp1, yp1                         )

           if ( xp1 < op%xmin .or.  &
                xp1 > op%xmax .or.  &
                yp1 < op%ymin .or.  &
                yp1 > op%ymax ) cycle

           write(string,'(I0)') itab_wd(iw)%mrow
           call o_plchlq (xp1,yp1,trim(adjustl(string)),0.0025,0.,0.)
        enddo

        call o_frame()
        call o_clswk()
        !$omp end single

     endif ! mod(iter,*)

  enddo ! iter
  !$omp end parallel

  xem(:) = real(xem8(:))
  yem(:) = real(yem8(:))
  zem(:) = real(zem8(:))

end subroutine spring_dynamics1
