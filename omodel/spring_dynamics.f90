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

subroutine spring_dynamics( mrows, moveint, ngr, nxp, nma, nua, nwa, &
                            xem, yem, zem, itab_md, itab_ud, itab_wd )

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars

  implicit none

  integer, intent(in) :: mrows, moveint, ngr, nxp, nma, nua, nwa

  real, intent(inout) :: xem(nma), yem(nma), zem(nma)

  type (itab_md_vars), intent(inout) :: itab_md(nma)
  type (itab_ud_vars), intent(inout) :: itab_ud(nua)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwa)

  if (ngr <= 1) then

     call spring_dynamics1( mrows, moveint, ngr, nxp, nma, nua, nwa, &
                            xem, yem, zem, itab_md, itab_ud, itab_wd )

  else

     call spring_dynamics_nest( mrows, moveint, ngr, nxp, nma, nua, nwa, &
                                xem, yem, zem, itab_md, itab_ud, itab_wd )

  endif

end subroutine spring_dynamics

!===============================================================================

subroutine spring_dynamics_nest( mrows, moveint, ngr, nxp, nma, nua, nwa, &
                                 xem, yem, zem, itab_md, itab_ud, itab_wd )

  ! Subroutine spring_dynamics_nest is used only for mesh refinements (on both the
  ! atm and sfc grids).  Call subroutine spring_dynamics1 to adjust base global grid

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars
  use mem_grid,     only: impent
  use consts_coms,  only: pi2, erad, r8
  use misc_coms,    only: io6, mdomain, deltax, runtype
  use oplot_coms,   only: op

  implicit none

  integer, intent(in) :: mrows, moveint, ngr, nxp, nma, nua, nwa

  real, intent(inout) :: xem(nma), yem(nma), zem(nma)

  type (itab_md_vars), intent(inout) :: itab_md(nma)
  type (itab_ud_vars), intent(inout) :: itab_ud(nua)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwa)

  integer            :: niter
  integer, parameter :: nprnt = 50
  real,    parameter :: relax = .04
! real,    parameter :: beta  = 1.242
  real               :: beta

  real, allocatable :: dist (:)
  real, allocatable :: dist0(:)
  real, allocatable :: dx   (:)
  real, allocatable :: dy   (:)
  real, allocatable :: dz   (:)
  real, allocatable :: dsm  (:)
  real, allocatable :: dirs (:,:)

  integer :: iu,iu1,iu3
  integer :: im,im1,im2,im3,im4
  integer :: iw,iw1,iw2,j
  integer :: iter,mrow1,mrow2,ngrw,mrmax,mrmin

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

  logical, allocatable :: moveu(:)
  logical, allocatable :: compu(:)

  integer :: nmovem, mend
  integer :: nmoveu
  integer :: ncompu
  integer :: nm, nc, jm, jm1, jm2
  integer :: ju, ju1, ju2, ju3, ju4

  real(r8), allocatable :: xem8(:), yem8(:), zem8(:)
  real(r8), allocatable :: xem0(:), yem0(:), zem0(:)

  real(r8) :: expansion, erad8

  character(10) :: string

  integer, allocatable :: im_jm(:)
  integer, allocatable :: jm_im(:)
  logical, allocatable :: movem(:)

  integer, allocatable :: jm_ju(:,:)
  integer, allocatable :: ju_ju(:,:)
  integer, allocatable :: ju_jm(:,:)

  integer, allocatable :: iu_ju(:)
  integer, allocatable :: ju_iu(:)

  integer, allocatable :: npoly_jm(:)

  if (runtype == 'MAKEGRID_PLOT') then
     niter = 50
  elseif (mrows == 3) then
     niter = 5000
     beta  = 1.242
  else
     niter = 3000
     beta  = 1.1
  endif

  erad8 = real(erad,r8)

  allocate(movem(nma))
  allocate(im_jm(nma))
  allocate(jm_im(nma))

  ! First, turn off movement of ALL M points

  nmovem   = 0
  movem(:) = .false.

  do im = 2, nma

     ! Check grid number of M point.  M point will move in spring dynamics only
     ! if it is of the current grid number (ngr) that was just spawned.  On the
     ! atmosphere grid, the row containing 7-way vertices (future centers of
     ! heptagon Voronoi cells) is the outermost row of M points that are
     ! designated to be on ngr.  On the surface grid, however, the 7-way
     ! vertices are NOT on ngr but rather on the lower-numbered parent grid.
     ! This prevents these M points from moving in spring dynamics.

     if (itab_md(im)%ngr /= ngr) cycle

     ! For preventing either polar M point from moving:
     ! if (im == impent(1 )) cycle
     ! if (im == impent(12)) cycle

     ! For preventing all pentagonal points from moving:
     ! if (any(im == impent(1:12))) cycle

     ! If interior M points are flagged to move, then all M points on ngr will
     ! be moved, so set movem flag to 1.

     if (moveint == 1) then

        nmovem = nmovem + 1
        im_jm( nmovem ) = im
        movem(im) = .true.
        jm_im(im) = nmovem

     else

        ! Even if interior points are not to be moved, M points in the
        ! transition row must always be moveable, so set their movem flag to 1.

        do j = 1, itab_md(im)%npoly
           iw = itab_md(im)%iw(j)

!          ngrw = itab_wd(iw)%ngr
!
!          if (ngrw /= ngr) cycle
!
!          if (itab_wd(iw)%mrow == -3 .or. &
!              itab_wd(iw)%mrow == -2 .or. &
!              itab_wd(iw)%mrow == -1 .or. &
!              itab_wd(iw)%mrow ==  1) then
           if (itab_wd(iw)%ngr == ngr .and. itab_wd(iw)%mrow /= 0) then

              nmovem = nmovem + 1
              im_jm( nmovem ) = im
              movem(im) = .true.
              jm_im(im) = nmovem

              exit  ! once this point has been flagged, we can move on
           endif
        enddo

     endif
  enddo

! Find U points adjacent to our flagged M points

  allocate(moveu(nua))
  allocate(compu(nua))

  nmoveu = 0
  ncompu = 0

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

     if (moveu(iu)) nmoveu = nmoveu + 1
     if (compu(iu)) ncompu = ncompu + 1
  enddo

  allocate( iu_ju(ncompu) )
  allocate( ju_iu(nua) )

! Arrays to convert between IU and our computational stencil

  nm = 0
  nc = nmoveu

  do iu = 2, nua
     if (moveu(iu)) then

        nm = nm + 1
        iu_ju(nm) = iu
        ju_iu(iu) = nm

     elseif (compu(iu)) then

        nc = nc + 1
        iu_ju(nc) = iu
        ju_iu(iu) = nc

     endif
  enddo

  deallocate(moveu)
  deallocate(compu)

! Add M points adjacent to flagged U points for our computational stencil

  mend = nmovem

  do ju = 1, ncompu
     iu = iu_ju(ju)

     do j = 1, 2
        im = itab_ud(iu)%im(j)

        if (.not. movem(im)) then
           mend = mend + 1
           im_jm( mend ) = im
           jm_im( im ) = mend
        endif
     enddo
  enddo

  deallocate(movem)

  allocate( dirs (7,nmovem) )
  allocate( ju_jm(7,nmovem) )
  allocate( npoly_jm(nmovem) )

  allocate( jm_ju(nmoveu,2) )
  allocate( ju_ju(nmoveu,4) )
  allocate( dist0(nmoveu) )

  write(*,'(4(A,I0))') "In spring dynamics: ngr = ", ngr, &
       ", nma = ", nma, ", nmovem = ", nmovem, ", niter = ", niter

  ! Compute mean length of coarse mesh U segments

  if (mdomain < 2) then
     dist00 = beta * pi2 * erad / (5. * real(nxp))
  else
     dist00 = deltax * sqrt( 2.0 / sqrt(3.0) )
  endif

  ! Compute target length of each triangle U segment
  ! Loop over all U points

  !$omp parallel private(iter)
  !$omp do private(iu,iw1,iw2,mrow1,mrow2,mrmax,mrmin)
  do ju = 1, nmoveu
     iu = iu_ju(ju)

     ! Store derived type values in arrays

     jm_ju(ju,1) = jm_im( itab_ud(iu)%im(1) )
     jm_ju(ju,2) = jm_im( itab_ud(iu)%im(2) )

     ju_ju(ju,1) = ju_iu( itab_ud(iu)%iu(1) )
     ju_ju(ju,2) = ju_iu( itab_ud(iu)%iu(2) )
     ju_ju(ju,3) = ju_iu( itab_ud(iu)%iu(3) )
     ju_ju(ju,4) = ju_iu( itab_ud(iu)%iu(4) )

     ! Compute target distance for any MRL value

     dist0(ju) = dist00 / real( 2**(itab_ud(iu)%mrlu - 1) )

     ! Modified distance in MRL border zone

     if (ngr > 1) then

        iw1 = itab_ud(iu)%iw(1)
        iw2 = itab_ud(iu)%iw(2)

        mrow1 = itab_wd(iw1)%mrow
        mrow2 = itab_wd(iw2)%mrow

        mrmax = max(mrow1,mrow2)
        mrmin = min(mrow1,mrow2)

        if (mrows >= 3) then

           if (.true.) then  ! Martin's adjustments to expand smallest cells

              if (mrmax == -4 .and. mrmin == -5) then
                 dist0(ju) = dist0(ju) * 12.1 / 12.
              elseif (mrmax == -4 .and. mrmin == -4) then
                 dist0(ju) = dist0(ju) * 12.4 / 12.
              elseif (mrmax == -3 .and. mrmin == -4) then
                 dist0(ju) = dist0(ju) * 12.8 / 12.
              elseif (mrmax == -3 .and. mrmin == -3) then
                 dist0(ju) = dist0(ju) * 13.2 / 12.
              elseif (mrmax == -2 .and. mrmin == -3) then
                 dist0(ju) = dist0(ju) * 13.6 / 12.
              elseif  (mrmax == -2 .and. mrmin == -2) then
                 dist0(ju) = dist0(ju) * 14.0 / 12.
              elseif (mrmax == -1 .and. mrmin == -2) then
                 dist0(ju) = dist0(ju) * 16.0 / 12.
              elseif (mrmax == -1 .and. mrmin == -1) then
                 dist0(ju) = dist0(ju) * 17.0 / 12.
              elseif (mrmax == 1 .and. mrmin == -1) then
                 dist0(ju) = dist0(ju) * 18.0 / 12.
              elseif (mrmax == 1 .and. mrmin == 1) then
                 dist0(ju) = dist0(ju) * 11.0 / 12.
              elseif (mrmax == 2 .and. mrmin == 1) then
                 dist0(ju) = dist0(ju) * 11.3 / 12.
              elseif (mrmax == 2 .and. mrmin == 2) then
                 dist0(ju) = dist0(ju) * 11.5 / 12.
              elseif (mrmax == 3 .and. mrmin == 2) then
                 dist0(ju) = dist0(ju) * 11.6 / 12.
              elseif (mrmax == 3 .and. mrmin == 3) then
                 dist0(ju) = dist0(ju) * 11.7 / 12.
              elseif (mrmax == 4 .and. mrmin == 3) then
                 dist0(ju) = dist0(ju) * 11.8 / 12.
              elseif (mrmax == 4 .and. mrmin == 4) then
                 dist0(ju) = dist0(ju) * 11.9 / 12.
              elseif (mrmax == 5 .and. mrmin == 4) then
                 dist0(ju) = dist0(ju) * 11.95 / 12.
              endif

           else

              if     (mrmax == -2 .and. mrmin == -2) then
                 dist0(ju) = dist0(ju) *  7. / 6.  !* .90
              elseif (mrmax == -1 .and. mrmin == -2) then
                 dist0(ju) = dist0(ju) *  8. / 6.  !* .90
              elseif (mrmax == -1 .and. mrmin == -1) then
                 dist0(ju) = dist0(ju) *  9. / 6.  !* .90
              elseif (mrmax == 1 .and. mrmin == -1) then
                 dist0(ju) = dist0(ju) * 10. / 6.  !* .90
              elseif (mrmax == 1 .and. mrmin == 1) then
                 dist0(ju) = dist0(ju) * 11. / 12. !* .90
              endif

           endif

        else ! (mrows = 1) Can be used for sfc grid but not atm grid

           if (mrmax == 1 .and. mrmin == 1) then
              dist0(ju) = dist0(ju) * 0.75
           endif

        endif

     endif ! ngr

  enddo
  !$omp end do

  !$omp do private(im,j,iu)
  do jm = 1, nmovem
     im = im_jm(jm)

     npoly_jm(jm) = itab_md(im)%npoly

     do j = 1, itab_md(im)%npoly
        iu = itab_md(im)%iu(j)

        ju_jm(j,jm) = ju_iu(iu)

        if (itab_ud(iu)%im(2) == im) then
           dirs(j,jm) =  relax
        else
           dirs(j,jm) = -relax
        endif

     enddo
  enddo
  !$omp end do

  !$omp single
  deallocate( ju_iu )

  allocate( dist(ncompu)) ; dist = 0.
  allocate( dx  (nmoveu)) ; dx = 0.
  allocate( dy  (nmoveu)) ; dy = 0.
  allocate( dz  (nmoveu)) ; dz = 0.

  allocate(xem8(mend))
  allocate(yem8(mend))
  allocate(zem8(mend))

  allocate(xem0(nmovem))
  allocate(yem0(nmovem))
  allocate(zem0(nmovem))
  allocate(dsm (nmovem))
  !$omp end single

! Copy Cartesian coordinates to double precision arrays

  !$omp do private(im)
  do jm = 1, mend
     im = im_jm(jm)

     xem8(jm) = xem(im)
     yem8(jm) = yem(im)
     zem8(jm) = zem(im)
  enddo
  !$omp end do

! These points only need to be computed and stored once; they don't change
! during the iteration

  !$omp do private(iu,jm1,jm2)
  do ju = nmoveu+1, ncompu
     iu = iu_ju(ju)

     jm1 = jm_im( itab_ud(iu)%im(1) )
     jm2 = jm_im( itab_ud(iu)%im(2) )

     dist(ju) = sqrt( real( xem8(jm2) - xem8(jm1) ) ** 2 &
                    + real( yem8(jm2) - yem8(jm1) ) ** 2 &
                    + real( zem8(jm2) - zem8(jm1) ) ** 2 )
  enddo
  !$omp end do

  !$omp single
  deallocate(jm_im)
  !$omp end single

! Main iteration loop

  do iter = 1, niter

     if (iter == 1 .or. mod(iter,nprnt) == 0) then

        !$omp do
        do jm = 1, nmovem
           xem0(jm) = xem8(jm)
           yem0(jm) = yem8(jm)
           zem0(jm) = zem8(jm)
        enddo
        !$omp end do nowait

     endif

     ! Compute length of each U segment

     !$omp do private(jm1,jm2)
     do ju = 1, nmoveu

        jm1 = jm_ju(ju,1)
        jm2 = jm_ju(ju,2)

        dx(ju) = real( xem8(jm2) - xem8(jm1) )
        dy(ju) = real( yem8(jm2) - yem8(jm1) )
        dz(ju) = real( zem8(jm2) - zem8(jm1) )

        dist(ju) = sqrt( dx(ju) * dx(ju) &
                       + dy(ju) * dy(ju) &
                       + dz(ju) * dz(ju) )
     enddo
     !$omp end do

     ! Adjustment of dist0 based on opposite angles of triangles

     !$omp do private(ju1,ju2,ju3,ju4,twocosphi3,twocosphi4,ratio, &
     !$omp            distm,frac_change)
     do ju = 1, nmoveu

        ju1 = ju_ju(ju,1)
        ju2 = ju_ju(ju,2)
        ju3 = ju_ju(ju,3)
        ju4 = ju_ju(ju,4)

        ! Compute two*cosine of angles at IM3 and IM4

        twocosphi3 = (dist(ju1)**2 + dist(ju2)**2 - dist(ju)**2) / (dist(ju1) * dist(ju2))
        twocosphi4 = (dist(ju3)**2 + dist(ju4)**2 - dist(ju)**2) / (dist(ju3) * dist(ju4))

        ! Decrease dist0 if smaller twocosine is less than limiting value of two*cos(72 deg)

        ratio = min(twocosphi3,twocosphi4)

        if (ratio < .61) then
           distm = dist0(ju) * max(ratio / .61, 0.1)
        else
           distm = dist0(ju)
        endif

        ! Fractional change to dist that would make it equal dist0

        frac_change = (distm - dist(ju)) / dist(ju)

        ! Compute components of displacement that gives dist0

        dx(ju) = dx(ju) * frac_change
        dy(ju) = dy(ju) * frac_change
        dz(ju) = dz(ju) * frac_change

     enddo
     !$omp end do

     ! Apply the displacement components to each M point

     !$omp do private(j,ju)
     do jm = 1, nmovem

        do j = 1, npoly_jm(jm)
           ju = ju_jm(j,jm)
           xem8(jm) = xem8(jm) + dirs(j,jm) * dx(ju)
           yem8(jm) = yem8(jm) + dirs(j,jm) * dy(ju)
           zem8(jm) = zem8(jm) + dirs(j,jm) * dz(ju)
        enddo

     enddo
     !$omp end do

     ! Push M point coordinates out to earth radius

     if (mdomain < 2) then

        !$omp do private(expansion)
        do jm = 1, nmovem

           expansion = erad8 / sqrt( xem8(jm) ** 2 + yem8(jm) ** 2 + zem8(jm) ** 2 )
           xem8(jm) = xem8(jm) * expansion
           yem8(jm) = yem8(jm) * expansion
           zem8(jm) = zem8(jm) * expansion

        enddo
        !$omp end do

     endif

     ! Print iteration status

     if (iter == 1 .or. mod(iter,nprnt) == 0) then

        !$omp do
        do jm = 1, nmovem

           dsm(jm) = real( (xem8(jm) - xem0(jm))**2 &
                         + (yem8(jm) - yem0(jm))**2 &
                         + (zem8(jm) - zem0(jm))**2 )
        enddo
        !$omp end do

        !$omp single
        write(*,'(3x,A,I5,A,I5,A,f0.4,A)') &
             "Iteration ", iter, " of ", niter, ",  Max DS = ", &
             sqrt( maxval( dsm ) ), " meters."
        !$omp end single

     endif

     ! Section for plotting grid at intermediate stages of spring dynamics adjustment
     ! Remove .false. to activate

     if (.false. .and. (iter == 1 .or. mod(iter,100) == 0)) then

        !$omp do private(im)
        do jm = 1, nmovem
           im = im_jm(jm)

           xem(im) = real(xem8(jm))
           yem(im) = real(yem8(jm))
           zem(im) = real(zem8(jm))
        enddo
        !$omp end do

        ! Plot grid lines

        !$omp single
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

  !$omp do private(im)
  do jm = 1, nmovem
     im = im_jm(jm)

     xem(im) = real(xem8(jm))
     yem(im) = real(yem8(jm))
     zem(im) = real(zem8(jm))
  enddo
  !$omp end do
  !$omp end parallel

end subroutine spring_dynamics_nest

!===============================================================================

subroutine spring_dynamics1( mrows, moveint, ngr, nxp, nma, nua, nwa, &
                             xem, yem, zem, itab_md, itab_ud, itab_wd )

  ! Subroutine spring_dynamics1 is used only for adjusting grid 1 (i.e.,
  ! the quasi-uniform global atm grid) prior to any mesh refinements.
  ! Call subroutine spring_dynamics to adjust mesh refinements of either
  ! the atm or surface grids.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars
  use mem_grid,     only: impent
  use consts_coms,  only: pi2, erad, r8, piu180
  use misc_coms,    only: io6, mdomain, deltax, runtype
  use oplot_coms,   only: op

  implicit none

  integer, intent(in) :: mrows, moveint, ngr, nxp, nma, nua, nwa

  real, intent(inout) :: xem(nma), yem(nma), zem(nma)

  type (itab_md_vars), intent(inout) :: itab_md(nma)
  type (itab_ud_vars), intent(inout) :: itab_ud(nua)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwa)

  integer            :: niter
  integer, parameter :: nprnt = 50
  real,    parameter :: relax = .035
  real,    parameter :: beta  = 1.242

! Automatic arrays

  real     :: dist(nua), distm
  real     :: ratio(nua), frac_change
  real     :: dx(nua), dy(nua), dz(nua)
  real     :: dirs(7,nma)

  integer  :: iu,iu1,iu2,iu3,iu4
  integer  :: im,im1,im2,im3
  integer  :: iw,j
  integer  :: iter

  integer  :: iskip
  real     :: xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2
  real     :: twocosphi3, twocosphi4
  real     :: dist00, disto61
  real     :: dsm(nma)

  real(r8) :: xem8(nma),yem8(nma),zem8(nma)
  real(r8) :: xem0(nma),yem0(nma),zem0(nma)
  real(r8) :: expansion, erad8

  character(10) :: string

  integer :: iumn(nua,2)
  integer :: iuun(nua,4)
  integer :: imnp(nma), imiu(7,nma)

  niter = 5000
  if (runtype == 'MAKEGRID_PLOT') niter = 50

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

  disto61 = dist00 / .61

  write(io6,'(a,4i9)') "In spring dynamics: ngr,nma,niter = ",1,nma,niter

  !$omp parallel private(iter)
  !$omp do
  do iu = 2, nua
     iumn(iu,1) = itab_ud(iu)%im(1)
     iumn(iu,2) = itab_ud(iu)%im(2)

     iuun(iu,1) = itab_ud(iu)%iu(1)
     iuun(iu,2) = itab_ud(iu)%iu(2)
     iuun(iu,3) = itab_ud(iu)%iu(3)
     iuun(iu,4) = itab_ud(iu)%iu(4)
  enddo
  !$omp end do nowait

  !$omp do private(j,iu)
  do im = 2, nma
     imnp(im) = itab_md(im)%npoly

     do j = 1, itab_md(im)%npoly
        iu = itab_md(im)%iu(j)

        imiu(j,im) = iu

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

        !$omp do
        do im = 2, nma
           xem0(im) = xem8(im)
           yem0(im) = yem8(im)
           zem0(im) = zem8(im)
        enddo
        !$omp end do nowait

     endif

! Compute length of each U segment

     !$omp do private(im1,im2)
     do iu = 2, nua
        im1 = iumn(iu,1)
        im2 = iumn(iu,2)

        dx(iu) = real( xem8(im2) - xem8(im1) )
        dy(iu) = real( yem8(im2) - yem8(im1) )
        dz(iu) = real( zem8(im2) - zem8(im1) )
     enddo
     !$omp end do nowait

     !$omp do
     do iu = 2, nua
        dist(iu) = sqrt( dx(iu) * dx(iu) &
                       + dy(iu) * dy(iu) &
                       + dz(iu) * dz(iu) )
     enddo
     !$omp end do

! Adjustment of dist0 based on opposite angles of triangles

     !$omp do private(iu1,iu2,iu3,iu4,twocosphi3,twocosphi4)
     do iu = 2, nua
        iu1 = iuun(iu,1)
        iu2 = iuun(iu,2)
        iu3 = iuun(iu,3)
        iu4 = iuun(iu,4)

        ! Compute cosine of angles at IM3 and IM4

        twocosphi3 = (dist(iu1)**2 + dist(iu2)**2 - dist(iu)**2) / (dist(iu1) * dist(iu2))
        twocosphi4 = (dist(iu3)**2 + dist(iu4)**2 - dist(iu)**2) / (dist(iu3) * dist(iu4))

        ! Ratio of smaller cosine to limiting value of cos(72 deg)

        ratio(iu) = min(twocosphi3,twocosphi4)
     enddo
     !$omp end do nowait

     !$omp do private(distm,frac_change)
     do iu = 2, nua

        if (ratio(iu) < .61) then
           distm = disto61 * ratio(iu)
        else
           distm = dist00
        endif

        ! Fractional change to dist that would make it equal dist0

        frac_change = (distm - dist(iu)) / dist(iu)

        ! Compute components of displacement that gives dist0

        dx(iu) = dx(iu) * frac_change
        dy(iu) = dy(iu) * frac_change
        dz(iu) = dz(iu) * frac_change

     enddo
     !$omp end do

     !$omp do private(j,iu)
     do im = 2, nma

        ! For preventing either polar M point from moving:
        if (im == impent(1 )) cycle
        if (im == impent(12)) cycle

        ! For preventing all pentagonal points from moving:
        ! if (any(im == impent(1:12))) cycle

        ! Apply the displacement components to each M point
        do j = 1, imnp(im)
           iu = imiu(j,im)
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
           dsm(im) = real( (xem8(im) - xem0(im))**2 &
                         + (yem8(im) - yem0(im))**2 &
                         + (zem8(im) - zem0(im))**2 )
        enddo
        !$omp end do

        !$omp single
        write(*,'(3x,A,I5,A,I5,A,f0.4,A)') &
             "Iteration ", iter, " of ", niter, ",  Max DS = ", &
             sqrt( maxval(dsm) ), " meters."
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
