Module mem_regrid

  use ll_bins, only: itab_grid_vars
  implicit none

  integer :: nzp_og
  integer :: mdomain_og
  integer :: ndz_og
  integer :: mrls_og
  integer :: nxp_og

  real, allocatable :: hdz_og(:)
  real, allocatable :: dz_og(:)

  integer :: nza_og
  integer :: nwa_og

  integer, allocatable :: lpw_og(:) ! (nwa_og)

  real, allocatable :: wnx_og(:) ! (nwa_og)
  real, allocatable :: wny_og(:) ! (nwa_og)
  real, allocatable :: wnz_og(:) ! (nwa_og)

  real, allocatable :: xew_og(:) ! (nwa_og)
  real, allocatable :: yew_og(:) ! (nwa_og)
  real, allocatable :: zew_og(:) ! (nwa_og)

  integer, allocatable :: iwn_og(:,:) ! itab_w(:)%iw(:)

  Type itab_wadd_vars         ! data structure for WADD pts (individual rank)
     integer :: iw_og(3) = 1  ! intrp pts
     real    :: f_og(3) = 0.  ! intrp coeffs
  End Type itab_wadd_vars

  type (itab_wadd_vars), allocatable, target :: itab_wadd(:)  ! (dimension mwa)
  type (itab_grid_vars), allocatable         :: itab_wog (:)

  private :: itab_grid_vars

Contains

!==============================================================================

subroutine init_regrid()

  use misc_coms,  only: io6
  use mem_grid,   only: mwa, xew, yew, zew, glatw, glonw
  use ll_bins,    only: latlon_bins, gridcells_from_latlon_bins, binset_vars

  implicit none

  real :: scalprod, vecprodz, vecprodz_maxneg, vecprodz_minpos
  real :: b11,b21,b31,b12,b22,b32,b13,b23,b33
  real :: dist, dist_min, hlatw, hlonw, xi, yi, res
  real :: xin(7),yin(7)

  integer :: iw, jw_og, iw_og, iwnog, npoly, jmaxneg, jminpos
  integer :: iset, i, j
  type(binset_vars) :: bset

  integer, pointer             :: nwbin
  integer, pointer, contiguous :: iwbin(:)

  ! Allocate and fill latlon bins at different resolutions to compartmentalize
  ! all IW points of the OLD ATM grid

  res = 7150.e3 / real( nxp_og * 2**(mrls_og-1) ) * 10.

  call latlon_bins(nwa_og, res, itab_wog, bset)

  ! Allocate data structure that contains interpolation indices and weights

  allocate(itab_wadd(mwa))

  ! Loop over all local-domain IW points in NEW grid

  do iw = 2,mwa

     if (mod(iw,100000) == 0) write(io6,*) 'init_regrid iw, mwa ',iw,mwa

     ! Get ATM bin indices for this IW point of NEW grid

     hlonw = max(-179.9999,min(179.9999,glonw(iw)))
     hlatw = max( -89.9999,min( 89.9999,glatw(iw)))

     call gridcells_from_latlon_bins( hlatw, hlonw, bset, nwbin, iwbin )

     ! j = int(delat(iset) * (hlatw +  90.)) + 1
     ! i = int(bset_og(iset)%blat(j)%delon * (hlonw + 180.)) + 1

     !   ! Check if bins are allocated, starting with smallest
     !   if (allocated(bset_og(iset)%blat(j)%bins(i)%iw)) then
     !      nwbin = bset_og(iset)%blat(j)%bins(i)%nw
     !      allocate(iwbin(nwbin))
     !      iwbin(:) = bset_og(iset)%blat(j)%bins(i)%iw(:)
     !      exit
     !   endif

     ! If iset loop has not branched to statement 10, search over all IW points
     !
     !print*, 'No allocated ATM bin found for iw = ', iw, hlatw, hlonw
     !print*, 'Conducting search over all IW points '

     !nwbin = nwa_og - 1
     !allocate(iwbin(nwbin))
     !do iw_og = 2,nwa_og
     !   iwbin(iw_og-1) = iw_og
     !enddo

     ! Loop over all OLD GRID IW points in bin

     dist_min = 1.e9

     do jw_og = 1, nwbin
        iw_og = iwbin(jw_og)

        ! Compute distance between NEW and OLD grid IW points, and reset minimum
        ! distance and OLD grid point index if new minimum is found

        dist = sqrt( (xew(iw) - xew_og(iw_og))**2 &
                   + (yew(iw) - yew_og(iw_og))**2 &
                   + (zew(iw) - zew_og(iw_og))**2 )

        if (dist_min > dist) then
           dist_min = dist
           itab_wadd(iw)%iw_og(1) = iw_og
        endif
     enddo

     deallocate(iwbin)

     ! Now that primary OLD grid point has been found for NEW grid IW point,
     ! compute its x,y components on a polar stereographic plane tangent
     ! at NEW grid point (NEW grid W point is at 0,0)

     iw_og = itab_wadd(iw)%iw_og(1)

     call e_ps(xew_og(iw_og),yew_og(iw_og),zew_og(iw_og), &
               glatw(iw),glonw(iw),xi,yi)

     ! Initialize vecprodz_minpos and vecprodz_maxneg

     vecprodz_maxneg = -1.e15
     vecprodz_minpos =  1.e15

     ! Loop through nearest polygon neighbors (j, iwn_og) of OLD grid point IW_og

     npoly = itab_wog(iw_og)%np

     do j = 1,npoly

        ! Get nudging point index (iwnudn) for current polygon neighbor of iwnud.

        iwnog = iwn_og(j,iw_og)

        ! Compute x,y components of iwn_og polygon center on a polar stereographic
        ! plane tangent at IW point

        call e_ps(xew_og(iwnog),yew_og(iwnog),zew_og(iwnog), &
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

           ! Identify maximum negative vecprodz among all iwn_og polygon neighbors of
           ! iw_og.  This iwn_og will be second point of HISTORY grid IW_og triad for IW

           if (vecprodz < 0. .and. vecprodz > vecprodz_maxneg) then
              vecprodz_maxneg = vecprodz
              jmaxneg = j
              itab_wadd(iw)%iw_og(2) = iwnog
           endif

           ! Identify minimum positive vecprodz among all iwn_og polygon neighbors of
           ! iw_og.  This iwn_og will be third point of HISTORY grid IW_og triad for IW

           if (vecprodz >= 0. .and. vecprodz < vecprodz_minpos) then
              vecprodz_minpos = vecprodz
              jminpos = j
              itab_wadd(iw)%iw_og(3) = iwnog
           endif
        endif

     enddo

     ! Lastly, fill 3 HISTORY grid weight coefficients for this W point.
     ! Weights are computed in 2_d polar stereographic space.

     ! Invert matrix of coordinates

     call matinv3x3(1.,xi,yi, &
                    1.,xin(jmaxneg),yin(jmaxneg), &
                    1.,xin(jminpos),yin(jminpos), &
                    b11,b21,b31,b12,b22,b32,b13,b23,b33)

     ! Assign coefficients

     itab_wadd(iw)%f_og(1) = b11
     itab_wadd(iw)%f_og(2) = b21
     itab_wadd(iw)%f_og(3) = b31

  enddo

end subroutine init_regrid

!=========================================================================

subroutine interp_regrid(ndims, idims, jdims, varn, stagpt,  &
                         ivara1,rvara1,rvara2,dvara1,dvara2, &
                         ivarb1,rvarb1,rvarb2,dvarb1,dvarb2  )

  use mem_grid,    only: mwa, zt, lpw
  use consts_coms, only: r8

  implicit none

  integer,      intent(in) :: ndims
  integer,      intent(in) :: idims(3), jdims(3)
  character(*), intent(in) :: varn ! Variable label
  character(2), intent(in) :: stagpt

  integer,  optional, intent(inout) :: ivara1(idims(1))
  real,     optional, intent(inout) :: rvara1(idims(1))
  real(r8), optional, intent(inout) :: dvara1(idims(1))
  real,     optional, intent(inout) :: rvara2(idims(1),idims(2))
  real(r8), optional, intent(inout) :: dvara2(idims(1),idims(2))

  integer,  optional, intent(inout) :: ivarb1(jdims(1))
  real,     optional, intent(inout) :: rvarb1(jdims(1))
  real(r8), optional, intent(inout) :: dvarb1(jdims(1))
  real,     optional, intent(inout) :: rvarb2(jdims(1),jdims(2))
  real(r8), optional, intent(inout) :: dvarb2(jdims(1),jdims(2))

  integer :: iw, k
  integer :: iw_og, iw_og1, iw_og2, iw_og3, ka_og
  real :: f_og1, f_og2, f_og3, cof

  ! Apply lower boundary condition to OLD grid fields in case they have not been
  ! defined below LPW level

  do iw_og = 2,nwa_og
     ka_og = lpw_og(iw_og)

     if     (present(rvara2)) then
        rvara2(1:ka_og-1,iw_og) = rvara2(ka_og,iw_og)
     elseif (present(dvara2)) then
        if (trim(varn) == 'RHO'        .or. &
            trim(varn) == 'PRESS'      .or. &
            trim(varn) == 'PRESS_ACCUM') then

           ! Extrapolate density and pressure below ground level to give
           ! values roughly in accord with hydrostatic balance

           do k = ka_og-1,2,-1
              cof = (zt(k) - zt(k+1)) / (zt(k+1) - zt(k+2))
              dvara2(k,iw_og) = dvara2(k+1,iw_og) * (1. + cof) &
                   - dvara2(k+2,iw_og) * cof
           enddo
        else
           dvara2(1:ka_og-1,iw_og) = dvara2(ka_og,iw_og)
        endif
     endif
  enddo

  ! Loop over all IW points on NEW grid

  do iw = 2,mwa

     ! Interpolate values for this IW point from a trio of OLD grid points

     iw_og1 = itab_wadd(iw)%iw_og(1); f_og1 = itab_wadd(iw)%f_og(1)
     iw_og2 = itab_wadd(iw)%iw_og(2); f_og2 = itab_wadd(iw)%f_og(2)
     iw_og3 = itab_wadd(iw)%iw_og(3); f_og3 = itab_wadd(iw)%f_og(3)

     if     (present(ivara1)) then
        ivarb1(iw) = nint(f_og1 * real(ivara1(iw_og1)) &
                   +      f_og2 * real(ivara1(iw_og2)) &
                   +      f_og3 * real(ivara1(iw_og3)))

        if (trim(varn) == 'KPBLH') ivarb1(iw) = max(ivarb1(iw),lpw(iw))

     elseif (present(rvara1)) then
        rvarb1(iw) = f_og1 * rvara1(iw_og1) &
                   + f_og2 * rvara1(iw_og2) &
                   + f_og3 * rvara1(iw_og3)
     elseif (present(rvara2)) then
        rvarb2(:,iw) = f_og1 * rvara2(:,iw_og1) &
                     + f_og2 * rvara2(:,iw_og2) &
                     + f_og3 * rvara2(:,iw_og3)
     elseif (present(dvara1)) then
        dvarb1(iw) = f_og1 * dvara1(iw_og1) &
                   + f_og2 * dvara1(iw_og2) &
                   + f_og3 * dvara1(iw_og3)
     elseif (present(dvara2)) then
        dvarb2(:,iw) = f_og1 * dvara2(:,iw_og1) &
                     + f_og2 * dvara2(:,iw_og2) &
                     + f_og3 * dvara2(:,iw_og3)
     endif

  enddo

end subroutine interp_regrid

!==============================================================================

subroutine close_regrid()

  implicit none

  if ( allocated(lpw_og) ) deallocate(lpw_og)

  if ( allocated(wnx_og) ) deallocate( wnx_og )
  if ( allocated(wny_og) ) deallocate( wny_og )
  if ( allocated(wnz_og) ) deallocate( wnz_og )

  if ( allocated(xew_og) ) deallocate( xew_og )
  if ( allocated(yew_og) ) deallocate( yew_og )
  if ( allocated(zew_og) ) deallocate( zew_og )

  if ( allocated(itab_wog)  ) deallocate( itab_wog  )
  if ( allocated(itab_wadd) ) deallocate( itab_wadd )

  if ( allocated(iwn_og) ) deallocate( iwn_og )

end subroutine close_regrid

!==============================================================================

End Module mem_regrid
