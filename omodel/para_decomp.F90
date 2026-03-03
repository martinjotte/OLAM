subroutine para_decomp()

! Decompose global grid into multiple subdomains for parallel computation

  use mem_para,   only: mgroupsize, myrank
  use misc_coms,  only: mdomain, iparallel
  use mem_ijtabs, only: itab_v_pd, itab_m_pd, itabg_m, itabg_v, itabg_w
  use mem_nudge,  only: nudflag, nudnxp, nwnud, itabg_wnud, &
                        xewnud, yewnud, zewnud
  use mem_grid,   only: nma, nva, nwa, xew, yew, zew, arw0
  use mem_sfcg,   only: nmsfc, itabg_msfc, itab_msfc_pd, &
                        nvsfc, itabg_vsfc, itab_wsfc_pd, &
                        nwsfc, itabg_wsfc, itab_vsfc_pd
  use consts_coms,only: i8
  !$ use omp_lib, only: omp_get_max_threads

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real,    parameter :: scal = 1.0
  integer, parameter :: nbins = 3000

  integer :: im,iv,iw,iwnud
  integer :: iw1,iw2,im1,im2
  integer :: igp,jgp
  integer :: i, j, ii, jj, iinud, jjnud

  integer :: iter,ibin,ierr,irk
  integer :: ngroups
  integer :: numtot, numcut, numcent

  real :: xmin,ymin,zmin
  real :: xmax,ymax,zmax
  real :: cmin, cmax, cmin0, cut

  integer :: igsize(mgroupsize)
  integer :: nwg   (mgroupsize)
  integer :: nwgnud(mgroupsize)

  integer :: irank(mgroupsize)
  real    :: cpufac(mgroupsize)
  real    :: fi, fj
!$ real, allocatable :: ompfac(:)

  integer, allocatable :: iwtemp(:), jwtemp(:)
  integer, allocatable :: iwnudtemp(:), jwnudtemp(:)

  integer :: num(nbins + 2)

  real :: val (nwa)
  real :: valnud(nwnud)

  Type grp_var
     integer, allocatable :: iw(:)
     integer, allocatable :: iwnud(:)
  End Type grp_var

  type (grp_var) :: grp(mgroupsize)

  integer :: iwsfc,ivsfc,imsfc,iwn
  integer :: irankw(3), irankv(3)

  real :: xsca, ysca, zsca, xo, yo, zo

  real, allocatable :: ds(:)

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
     allocate (itabg_wnud(nwnud))
  endif

  if (iparallel == 1) then

     ! cpufac is to take into account the relative cpu speed difference between
     ! nodes. For example, with 4 nodes and nodes 3 and 4 are half the spped of
     ! the others, then set cpufac = [1.0, 1.0, 0.5, 0.5].

     cpufac(:) = 1.0

     ! ompfac is to take into account differences in the number of cpu cores
     ! available between nodes when running in hybrid MPI/OpenMP mode on
     ! heterogeneous clusters.

#ifdef OLAM_MPI
     !$ allocate( ompfac( mgroupsize ) )
     !$ ompfac(myrank+1) = real( omp_get_max_threads() )
     !$ call MPI_Allgather( MPI_IN_PLACE, 1, MPI_REAL,  &
     !$                     ompfac, 1, MPI_REAL, MPI_COMM_WORLD, ierr )
     !$ ompfac(:) = ompfac(:) / maxval( ompfac )
     !$ cpufac(:) = cpufac(:) * ompfac(:)
     !$ deallocate( ompfac )
#endif

     allocate(iwtemp(nwa-1))
     allocate(jwtemp(nwa-1))

     ! Allocate and fill grp%iw for group 1

     allocate (grp(1)%iw(nwa-1))
     allocate (ds(nwa))

     do iw = 2,nwa
        grp(1)%iw(iw-1) = iw
        ds(iw) = scal * sqrt(arw0(iw))
     enddo

     ! Check if NUDGING GRID is being used

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

        allocate(iwnudtemp(nwnud-1))
        allocate(jwnudtemp(nwnud-1))

        allocate (grp(1)%iwnud(nwnud-1))

        do iw = 2,nwnud
           grp(1)%iwnud(iw-1) = iw
        enddo
     endif

     ngroups = 1
     jgp = 1
     igsize(1) = mgroupsize

     irank(1) = 1

     nwg (1) = nwa - 1
     nwgnud(1) = nwnud - 1

     do while (ngroups < mgroupsize)

        do igp = 1,ngroups

           if (igsize(igp) > 1) then

              jgp = jgp + 1

              xmin = 1.e9
              ymin = 1.e9
              zmin = 1.e9

              xmax = -1.e9
              ymax = -1.e9
              zmax = -1.e9

! Find max and min x,y,z of current group

              do i = 1,nwg(igp)
                 iw = grp(igp)%iw(i)

                 if (xmin > xew(iw) - ds(iw)) xmin = xew(iw) - ds(iw)
                 if (ymin > yew(iw) - ds(iw)) ymin = yew(iw) - ds(iw)
                 if (zmin > zew(iw) - ds(iw)) zmin = zew(iw) - ds(iw)

                 if (xmax < xew(iw) + ds(iw)) xmax = xew(iw) + ds(iw)
                 if (ymax < yew(iw) + ds(iw)) ymax = yew(iw) + ds(iw)
                 if (zmax < zew(iw) + ds(iw)) zmax = zew(iw) + ds(iw)
              enddo

! Determine whether to cut in x, y, or z direction

              if ( 1.01 * (zmax - zmin) > xmax - xmin  .and.  &
                   1.01 * (zmax - zmin) > ymax - ymin ) then

! ATM cells - z direction

                 if (xmax - xmin > ymax - ymin) then

                    xsca = 2.0 / (xmax - xmin)
                    xo   = 0.5 * (xmax + xmin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = zew(iw) + ds(iw) * (xew(iw) - xo) * xsca
                    enddo

                 else

                    ysca = 2.0 / (ymax - ymin)
                    yo   = 0.5 * (ymax + ymin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = zew(iw) + ds(iw) * (yew(iw) - yo) * ysca
                    enddo

                 endif

                 cmin = zmin
                 cmax = zmax

! WNUD cells - z direction

                 if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
                    do i = 1,nwgnud(igp)
                       iw = grp(igp)%iwnud(i)
                       valnud(iw) = zewnud(iw)
                    enddo
                 endif

              elseif (1.01 * (xmax - xmin) > ymax - ymin) then

! ATM cells - x direction

                 if (ymax - ymin > zmax - zmin) then

                    ysca = 2.0 / (ymax - ymin)
                    yo   = 0.5 * (ymax + ymin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = xew(iw) + ds(iw) * (yew(iw) - yo) * ysca
                    enddo

                 else

                    zsca = 2.0 / (zmax - zmin)
                    zo   = 0.5 * (zmax + zmin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = xew(iw) + ds(iw) * (zew(iw) - zo) * zsca
                    enddo

                 endif

                 cmin = xmin
                 cmax = xmax

! WNUD cells - x direction

                 if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
                    do i = 1,nwgnud(igp)
                       iw = grp(igp)%iwnud(i)
                       valnud(iw) = xewnud(iw)
                    enddo
                 endif

              else

! ATM cells - y direction

                 if (xmax - xmin > zmax - zmin) then

                    xsca = 2.0 / (xmax - xmin)
                    xo   = 0.5 * (xmax + xmin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = yew(iw) + ds(iw) * (xew(iw) - xo) * xsca
                    enddo

                 else

                    zsca = 2.0 / (zmax - zmin)
                    zo   = 0.5 * (zmax + zmin)

                    do i = 1,nwg(igp)
                       iw = grp(igp)%iw(i)
                       val(iw) = yew(iw) + ds(iw) * (zew(iw) - zo) * zsca
                    enddo

                 endif

                 cmin = ymin
                 cmax = ymax

! WNUD cells - y direction

                 if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
                    do i = 1,nwgnud(igp)
                       iw = grp(igp)%iwnud(i)
                       valnud(iw) = yewnud(iw)
                    enddo
                 endif

              endif

! Determine cut value, iterating 3 times

              do iter = 1,3

                 num(:) = 0
                 numcent = 0

                 do i = 1,nwg(igp)
                    iw = grp(igp)%iw(i)
                    if (val(iw) <= cmin) then
                       ibin = 1
                    elseif (val(iw) >= cmax) then
                       ibin = nbins + 2
                    else
                       ibin = int(real(nbins) * (val(iw) - cmin) / (cmax - cmin)) + 2
                    endif
                    num(ibin) = num(ibin) + 1

                    ! Sum points that are beyond geometric center of group
                    if (ibin > nbins/2 + 2) numcent = numcent + 1
                 enddo

! Set igsize(jgp) based on numcent

                 if (iter == 1) then

!                   igsize(jgp) = (numcent * igsize(igp)) / nwg(igp)
                    igsize(jgp) = real(numcent) * real(igsize(igp)) / real(nwg(igp))

                    igsize(jgp) = max(1, min(igsize(igp)-1, igsize(jgp)))

                    igsize(igp) = igsize(igp) - igsize(jgp)

                    irank(jgp) = irank(igp) + igsize(igp)

!                   numcut = nwg(igp) * igsize(igp) / (igsize(igp) + igsize(jgp))
                    fi = sum( cpufac( irank(igp):irank(igp)+igsize(igp)-1 ) )
                    fj = sum( cpufac( irank(jgp):irank(jgp)+igsize(jgp)-1 ) )
                    numcut = real(nwg(igp)) * fi / (fi + fj)

                 endif

! Sum number in each bin until reaching half of total

                 numtot = num(1)
                 ibin = 1
                 do while (numtot + num(ibin+1) < numcut)
                    numtot = numtot + num(ibin+1)
                    ibin = ibin + 1
                 enddo
                 cmin0 = cmin + (cmax - cmin) * (real(ibin) - 1.1) / real(nbins)
                 cmax = cmin + (cmax - cmin) * (real(ibin) +  .1) / real(nbins)
                 cmin = cmin0

              enddo
              cut = .5 * (cmin + cmax)

! Transfer a number of IW points from igp to jgp

              jj = 0
              ii = 0
              do i = 1,nwg(igp)
                 iw = grp(igp)%iw(i)

                 if (val(iw) > cut) then
                    jj = jj + 1
                    jwtemp(jj) = iw
                 else
                    ii = ii + 1
                    iwtemp(ii) = iw
                 endif
              enddo

              nwg(igp) = ii
              nwg(jgp) = jj

! Deallocate 1 old group and allocate 2 new groups

              deallocate(grp(igp)%iw)
              allocate(grp(igp)%iw(ii))
              allocate(grp(jgp)%iw(jj))

! Fill 2 new groups of ATM cellsfrom temporary arrays

              do j = 1,jj
                 grp(jgp)%iw(j)=jwtemp(j)
              enddo

              do i = 1,ii
                 grp(igp)%iw(i)=iwtemp(i)
              enddo

! Transfer a number of WNUD points from igp to jgp

              if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
                 jjnud = 0
                 iinud = 0
                 do i = 1,nwgnud(igp)
                    iw = grp(igp)%iwnud(i)

                    if (valnud(iw) > cut) then
                       jjnud = jjnud + 1
                       jwnudtemp(jjnud) = iw
                    else
                       iinud = iinud + 1
                       iwnudtemp(iinud) = iw
                    endif
                 enddo

                 nwgnud(igp) = iinud
                 nwgnud(jgp) = jjnud

! Deallocate 1 old group and allocate 2 new groups

                 deallocate(grp(igp)%iwnud)
                 allocate(grp(igp)%iwnud(iinud))
                 allocate(grp(jgp)%iwnud(jjnud))

! Fill 2 new groups of WNUD cells from temporary arrays

                 do j = 1,jjnud
                    grp(jgp)%iwnud(j)=jwnudtemp(j)
                 enddo

                 do i = 1,iinud
                    grp(igp)%iwnud(i)=iwnudtemp(i)
                 enddo
              endif

           endif
        enddo

        ngroups = jgp

     enddo

   ! Fill irank for each IW point from group array

     do igp = 1,ngroups
        irk = irank(igp) - 1

        ! ATM cells

        do i = 1,nwg(igp)
           iw = grp(igp)%iw(i)
           itabg_w(iw)%irank = irk
        enddo

        ! WNUD cells

        if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
           do i = 1,nwgnud(igp)
              iw = grp(igp)%iwnud(i)
              itabg_wnud(iw)%irank = irk
           enddo
        endif

     enddo

     ! Set V rank based on its nearest W neighbors

     do iv = 2,nva
        itabg_v(iv)%irank = itabg_w( itab_v_pd(iv)%iw( mod(iv,2)+1 ) )%irank
     enddo

     ! Set M rank based on its nearest W and V neighbors

     do im = 2, nma
        if (itab_m_pd(im)%npoly < 3) then

           itabg_m(im)%irank = itabg_w( itab_m_pd(im)%iw(1) )%irank

        else

           irankw(1:3) = itabg_w( itab_m_pd(im)%iw(1:3) )%irank
           irankv(1:3) = itabg_v( itab_m_pd(im)%iv(1:3) )%irank

           if (irankw(1) == irankw(2) .or. irankw(1) == irankw(3)) then
              itabg_m(im)%irank = irankw(1)
           elseif (irankw(2) == irankw(3)) then
              itabg_m(im)%irank = irankw(2)

           elseif (irankv(1) == irankv(2) .or. irankv(1) == irankv(3)) then
              itabg_m(im)%irank = irankv(1)
           elseif (irankv(2) == irankv(3)) then
              itabg_m(im)%irank = irankv(2)

           else
              itabg_m(im)%irank = irankw( mod(im,3)+1 )
           endif

        endif
     enddo

     ! Re-set V rank if nearest M neighbors have the same rank

     do iv = 2, nva
        iw1 = itab_v_pd(iv)%iw(1)
        iw2 = itab_v_pd(iv)%iw(2)
        if (itabg_w(iw1)%irank /= itabg_w(iw2)%irank) then
           im1 = itab_v_pd(iv)%im(1)
           im2 = itab_v_pd(iv)%im(2)
           if (itabg_m(im1)%irank == itabg_m(im2)%irank) itabg_v(iv)%irank = itabg_m(im1)%irank
        endif
     enddo

     if (mdomain <= 1) then

        ! Set the rank of the WSFC point to that of the ATM cell to which
        ! it is coupled with the largest coupling area

        do iwsfc = 2, nwsfc
           iw = itab_wsfc_pd(iwsfc)%iwatm(1) ! global index
           itabg_wsfc(iwsfc)%irank = itabg_w(iw)%irank
        enddo

        ! Set the rank of the VSFC point to that of either of its nearest
        ! WSFC neighbors

        do ivsfc = 2,nvsfc
           iwn = itab_vsfc_pd(ivsfc)%iwn( 1+mod(ivsfc,2) )
           itabg_vsfc(ivsfc)%irank = itabg_wsfc(iwn)%irank
        enddo

        ! Set the rank of the MSFC point to the dominant rank of its immediate
        ! neighbors

        do imsfc = 2,nmsfc
           irankw(1:3) = itabg_wsfc( itab_msfc_pd(imsfc)%iwn(1:3) )%irank
           irankv(1:3) = itabg_vsfc( itab_msfc_pd(imsfc)%ivn(1:3) )%irank

           ! Also ensure at least one adjacent W and one adjacent V neighbor
           ! share the same rank as this M point

           if (irankw(1) == irankw(2) .or. irankw(1) == irankw(3)) then
              itabg_msfc(imsfc)%irank = irankw(1)
           elseif (irankw(2) == irankw(3)) then
              itabg_msfc(imsfc)%irank = irankw(2)

           elseif (irankv(1) == irankv(2) .or. irankv(1) == irankv(3)) then
              itabg_msfc(imsfc)%irank = irankv(1)
           elseif (irankv(2) == irankv(3)) then
              itabg_msfc(imsfc)%irank = irankv(2)

           else
              itabg_msfc(imsfc)%irank = irankw( mod(imsfc,3)+1 )
           endif
        enddo

     endif

  else   ! If not running in parallel:

     do iw = 2, nwa
        itabg_w(iw)%irank = 0
     enddo

     do iv = 2, nva
        itabg_v(iv)%irank = 0
     enddo

     do im = 2, nma
        itabg_m(im)%irank = 0
     enddo

     if (mdomain <= 1) then
        do iwsfc = 2, nwsfc
           itabg_wsfc(iwsfc)%irank = 0
        enddo

        do ivsfc = 2, nvsfc
           itabg_vsfc(ivsfc)%irank = 0
        enddo

        do imsfc = 2, nmsfc
           itabg_msfc(imsfc)%irank = 0
        enddo
     endif

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        do iwnud = 2, nwnud
           itabg_wnud(iwnud)%irank = 0
        enddo
     endif

  endif

end subroutine para_decomp
