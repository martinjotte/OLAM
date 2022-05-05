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

subroutine tileslab_horiz_mp(iplt,action)

  use oplot_coms, only: op
  use mem_grid,   only: mza, mma, mwa, lpw, zm, zt, xem, yem, zem, &
                        xev, yev, zev, xew, yew, zew, arm0
  use mem_ijtabs, only: itab_m, itabg_m, jtab_m, jtm_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: j, jn, jnn, jj, jm, im, iw, npoly, iv1, iv2
  integer :: notavail, navail

  real :: hpt, vpt, psiz, vsprd
  real :: fldval

  real :: htpn(7), vtpn(7)
  real :: hqpn(4), vqpn(4)

  integer :: ktf(mwa), km(mma)
  real :: wtbot(mma), wttop(mma)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (action == 'T' .and. op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 19 * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mma) / 5. )
  else
     inc = mma
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over M points

  mloop: do jm = 1, jtab_m(jtm_prog)%jend(1)
            im = jtab_m(jtm_prog)%im(jm)

     ! For now, skip pts that don't read in topm
     ! if (.not. itab_m(im)%loop(1)) cycle mloop

     ! Get tile plot coordinates.

     call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

     npoly = itab_m(im)%npoly

     ! Initialize navail counter and iflag180

     navail   = 0
     iflag180 = 0

     do j = 1,npoly

        iw = itab_m(im)%iw(j)

        ! Skip this M point if current IW point index < 2 (which occurs at
        ! lateral boundary of limited-area domain or parallel subdomain)      

        if (iw < 2) cycle mloop

        call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(j),vtpn(j))

        ! Avoid wrap-around for lat-lon plot and set iflag180

        if (op%projectn(iplt) == 'L') then
           call ll_unwrap(hpt,htpn(j))
           if (htpn(j) < -180.001) iflag180 =  1
           if (htpn(j) >  180.001) iflag180 = -1
        endif

        ! Jump out of loop if any cell corner is on other side of earth

        if (htpn(j) > 1.e11) cycle mloop

        ! If any IW point for this IM point is above ground and below model top,
        ! set notavail flag to zero

        if (ktf(iw) == 0) navail = navail + 1

     enddo

     if (navail == 0) cycle mloop

     ! Jump out of loop if entire cell is outside plot window. 

     if ( all(htpn(1:npoly) < op%xmin) .or. &
          all(htpn(1:npoly) > op%xmax) .or. &
          all(vtpn(1:npoly) < op%ymin) .or. &
          all(vtpn(1:npoly) > op%ymax) ) cycle mloop

     ! Get field value

     call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                    fldval,notavail)

     if (notavail > 0) cycle mloop

     ! If ALL surrounding T cells are available (e.g., above ground), plot entire
     ! M cell as single polygon.

     if (navail == npoly) then

        call get_psiz(iplt,sqrt(arm0(im)),psiz,vsprd)

        if (myrank == 0) then

           call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(npoly,  1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(htpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif
        
        ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
        ! at other end

        if (iflag180 /= 0) then

           do j = 1,npoly
              htpn(j) = htpn(j) + 360. * iflag180
           enddo

           if (myrank == 0) then

              call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(npoly,  1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(htpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vtpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(hpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(psiz,   1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vsprd,  1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldval, 1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif ! iflag180

        ! If any surrounding T cell is unavailable, plot M cell by available sectors

     else  ! navail /= npoly

        hqpn(1) = hpt
        vqpn(1) = vpt

        iflag180 = 0

        do j = 1,npoly

           iw = itab_m(im)%iw(j)

           if (ktf(iw) == 0) then
              jn = j + 1
              if (jn > npoly) jn = 1
              jnn = jn + 1
              if (jnn > npoly) jnn = 1

              hqpn(3) = htpn(j)
              vqpn(3) = vtpn(j)

              ! Specific way to get IV1 and IV2 since ordering of W and U/V 
              ! neighbors of M is not identical for both grid systems

              iv1 = itab_m(im)%iv(jn)
              iv2 = itab_m(im)%iv(jnn)

              call oplot_transform(iplt,xev(iv1),yev(iv1),zev(iv1),hqpn(2),vqpn(2))
              call oplot_transform(iplt,xev(iv2),yev(iv2),zev(iv2),hqpn(4),vqpn(4))

              ! Avoid wrap-around for lat-lon plot and set iflag180

              if (op%projectn(iplt) == 'L') then
                 do jj = 1, 4
                    call ll_unwrap(hpt,hqpn(jj))
                    if (hqpn(jj) < -180.001) iflag180 =  1
                    if (hqpn(jj) >  180.001) iflag180 = -1
                 enddo
              endif

              if (myrank == 0) then

                 call celltile(iplt,4,hqpn,vqpn,hpt,vpt,psiz,vsprd,fldval,action)

              else

#ifdef OLAM_MPI
                 nu = nu + 1
                 if (buffsize < ipos + base) then
                    allocate( bcopy (buffsize + inc * base) )
                    bcopy(1:buffsize) = buffer
                    call move_alloc(bcopy, buffer)
                    buffsize = size(buffer)
                 endif

                 call MPI_Pack(4,      1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(hqpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vqpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(hpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(psiz,   1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(vsprd,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(fldval, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

              endif

              ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
              ! at other end

              if (iflag180 /= 0) then

                 do jj = 1, 4
                    hqpn(jj) = hqpn(jj) + 360. * iflag180
                 enddo

                 if (myrank == 0) then

                    call celltile(iplt,4,hqpn,vqpn,hpt,vpt,psiz,vsprd,fldval,action)

                 else

#ifdef OLAM_MPI
                    nu = nu + 1
                    if (buffsize < ipos + base) then
                       allocate( bcopy (buffsize + inc * base) )
                       bcopy(1:buffsize) = buffer
                       call move_alloc(bcopy, buffer)
                       buffsize = size(buffer)
                    endif

                    call MPI_Pack(4,      1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(hqpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(vqpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(hpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(vpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(psiz,   1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(vsprd,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                    call MPI_Pack(fldval, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

                 endif

              endif ! iflag180
           endif
        enddo
     endif ! navail == npoly

  enddo mloop

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, npoly,  1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,     MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine tileslab_horiz_mp

!===============================================================================

subroutine tileslab_horiz_tw(iplt,action)

  use oplot_coms, only: op, xepc, yepc, zepc
  use mem_grid,   only: mza, mwa, zm, zt, xew, yew, zew, xem, yem, zem, lpw, arw0
  use mem_ijtabs, only: itab_w, jtab_w, jtw_wadj
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: npoly, j, jw, im
  integer :: iw
  integer :: iv1,iv2
  integer :: iok
  integer :: notavail

  real :: hpt, vpt, psiz, vsprd
  real :: fldval

  real :: topo1, topo2

  real :: htpn(7), vtpn(7)
  real :: xq(2), yq(2)

  integer :: ktf(mwa), kw(mwa)
  real :: wtbot(mwa), wttop(mwa)

  real :: arw0_tot, field_tot
  real :: arws(mgroupsize), fields(mgroupsize)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (action == 'T' .and. op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  iok = 0
  xq = 0.
  yq = 0.

  arw0_tot  = 0.
  field_tot = 0.

  nu   = 0
  ipos = 0

  base = 24 * nbytes_real + 2 * nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mwa) / 5. )
  else
     inc = mwa
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over W points

  wloop: do jw = 1, jtab_w(jtw_wadj)%jend(1)
     iw = jtab_w(jtw_wadj)%iw(jw)

     if (itab_w(iw)%irank /= myrank) cycle wloop

     npoly = itab_w(iw)%npoly

     ! For hexagon grid (in limited-area domain), do not tile plot
     ! incomplete boundary cells

     if (action == 'T' .and. npoly < 5) cycle wloop

     ! Skip this point if it is underground

     if (ktf(iw) /= 0) cycle wloop

     ! Initialize iflag180
     
     iflag180 = 0

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

     ! If only printing value, skip polygon section

     if (.not. (action == 'P' .and. npoly < 5) ) then

        do j = 1,npoly

           im = itab_w(iw)%im(j)

           call oplot_transform(iplt,xem(im),yem(im),zem(im),htpn(j),vtpn(j))

           ! Avoid wrap-around for lat-lon plot and set iflag180

           if (op%projectn(iplt) == 'L') then
              call ll_unwrap(hpt,htpn(j))
              if (htpn(j) < -180.001) iflag180 =  1
              if (htpn(j) >  180.001) iflag180 = -1
           endif

           ! Jump out of loop if cell corner is on other side of earth

           if (htpn(j) > 1.e11) cycle wloop

        enddo

        ! Jump out of loop if entire cell is outside plot window

        if ( all(htpn(1:npoly) < op%xmin) .or. all(htpn(1:npoly) > op%xmax) .or.  &
             all(vtpn(1:npoly) < op%ymin) .or. all(vtpn(1:npoly) > op%ymax) ) cycle wloop

     endif  ! (.not. (action == 'P' .and. npoly < 5))

     ! Get cell value and plot if 'available'

     call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                    fldval,notavail)

     if (notavail > 0) cycle wloop

     arw0_tot  = arw0_tot  + arw0(iw)
     field_tot = field_tot + fldval * arw0(iw)

     call get_psiz(iplt,sqrt(arw0(iw)),psiz,vsprd)

     if (myrank == 0) then
        call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)
     endif

     ! Plot cone circle if so specified

     if (op%pltcone(iplt) == 'C' .and. action == 'T') then

        call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
  
        if (iok == 1) then

           call oplot_transform(iplt,xepc(1),yepc(1),zepc(1),xq(1),yq(1))
           call oplot_transform(iplt,xepc(2),yepc(2),zepc(2),xq(2),yq(2))

           if (myrank == 0) then
              call o_frstpt(xq(1),yq(1))
              call o_vector(xq(2),yq(2))
           endif

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

        call MPI_Pack(npoly,  1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(htpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(psiz,   1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vsprd,  1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

        if (op%pltcone(iplt) == 'C' .and. action == 'T') then
           call MPI_Pack(iok, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

           if (iok == 1) then
              call MPI_Pack(xq, 2, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq, 2, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           endif
        endif
     endif
#endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
     ! at other end

     if (iflag180 /= 0) then

        do j = 1,npoly
           htpn(j) = htpn(j) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(npoly,  1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(htpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

           if (op%pltcone(iplt) == 'C' .and. action == 'T') then
              call MPI_Pack(iok, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

              if (iok == 1) then
                 call MPI_Pack(xq, 2, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
                 call MPI_Pack(yq, 2, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              endif
           endif
#endif

        endif
            
     endif ! iflag180

  enddo wloop

#ifdef OLAM_MPI
  if (iparallel == 1) then

     ! first collect the areas and totals to compute the average

     call MPI_Gather(arw0_tot,  1, MPI_REAL, arws,   1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Gather(field_tot, 1, MPI_REAL, fields, 1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     if (myrank == 0 ) then
        arw0_tot  = sum(arws)
        field_tot = sum(fields)
     endif

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, npoly,  1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1,     MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,     MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

                 if (op%pltcone(iplt) == 'C' .and. action == 'T') then
                    call MPI_Unpack(buffer, buffsize, ipos, iok, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    if (iok == 1) then
                       call MPI_Unpack(buffer, buffsize, ipos, xq, 2, MPI_REAL, MPI_COMM_WORLD, ier)
                       call MPI_Unpack(buffer, buffsize, ipos, yq, 2, MPI_REAL, MPI_COMM_WORLD, ier)
                       call o_frstpt(xq(1),yq(1))
                       call o_vector(xq(2),yq(2))
                    endif
                 endif

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

  if (myrank == 0 .and. action == 'T') then
     if (arw0_tot > 1.) print*, trim(op%fldname(iplt)),' field_avg ', &
                                field_tot, arw0_tot, field_tot/arw0_tot
  endif

end subroutine tileslab_horiz_tw

!===============================================================================

subroutine tileslab_horiz_vn(iplt,action)

  use oplot_coms, only: op
  use mem_grid,   only: mza, mva, mwa, zm, xev, yev, zev, &
                        xem, yem, zem, xew, yew, zew, lpw, dnv
  use mem_ijtabs, only: itab_v, jtab_v, jtv_wadj
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: iv
  integer :: iw1, iw2
  integer :: itpn
  integer :: im1, im2
  integer :: notavail

  real :: fldval
  real :: hpt, vpt, psiz, vsprd
  real :: htpn(4), vtpn(4)

  integer :: ktf(mwa), kv(mva)
  real :: wtbot(mva), wttop(mva)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (action == 'T' .and. op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 13 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mva) / 5. )
  else
     inc = mva
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over V points

! do jv = 1, jtab_v(jtv_wadj)%jend(1)
  do iv = 2,mva

     if (itab_v(iv)%irank /= myrank) cycle 

     ! Transform tile plot X and Y coordinates.  

!     iv  = jtab_v(jtv_wadj)%iv(jv)
     call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hpt,vpt)

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! Skip this V point if both its W neighbors are below ground

     if (ktf(iw1) /= 0 .and. ktf(iw2) /= 0) cycle

     ! Initialize iflag180

     iflag180 = 0

     if (im1 > 1) then
        call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),htpn(1),vtpn(1))
     else
        htpn(1) = hpt
        vtpn(1) = vpt
     endif

     if (im2 > 1) then
        call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),htpn(3),vtpn(3))
     else
        htpn(3) = hpt
        vtpn(3) = vpt
     endif

     if (iw2 > 1) then
        call oplot_transform(iplt,xew(iw2),yew(iw2),zew(iw2),htpn(2),vtpn(2))
     else
        htpn(2) = htpn(3)
        vtpn(2) = vtpn(3)
     endif

     if (iw1 > 1) then
        call oplot_transform(iplt,xew(iw1),yew(iw1),zew(iw1),htpn(4),vtpn(4))
     else
        htpn(4) = htpn(3)
        vtpn(4) = vtpn(3)
     endif

     !  Avoid wrap-around for lat-lon plot

     if (op%projectn(iplt) == 'L') then
        do itpn = 1,4
           call ll_unwrap(hpt,htpn(itpn))
        enddo
     endif

     ! Jump out of loop if any cell corner is on other side of earth

     if (any(htpn(1:4) > 1.e11)) cycle

     ! Jump out of loop if entire cell is outside plot window. 

     if ( all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
          all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax)) cycle

     ! Get cell value and 'available' flag

     call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                    fldval,notavail)    

     if (notavail > 0) cycle 

     ! If both W neighbors are above ground, plot full cell

     if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

        ! Else, check if IW1 is above ground

     elseif (ktf(iw1) == 0) then

        ! Plot IW1 half of cell with tile color

        htpn(2) = htpn(3)
        vtpn(2) = vtpn(3)

     else

        ! Plot IW2 half of cell with tile color

        htpn(4) = htpn(3)
        vtpn(4) = vtpn(3)

     endif

     ! Set iflag180
        
     if (op%projectn(iplt) == 'L') then
        do itpn = 1,4
           if (htpn(itpn) < -180.001) iflag180 =  1
           if (htpn(itpn) >  180.001) iflag180 = -1
        enddo
     endif

     call get_psiz(iplt,dnv(iv),psiz,vsprd)

     if (myrank == 0) then

        call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)
        
     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(psiz,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vsprd,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
     ! at other end

     if (iflag180 /= 0) then

        do itpn = 1, 4
           htpn(itpn) = htpn(itpn) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif ! iflag180

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)
           
              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine tileslab_horiz_vn

!===============================================================================

subroutine tileslab_horiz_wsfc(iplt,action)

  use max_dims,   only: maxnlspoly
  use oplot_coms, only: op
  use mem_sfcg,   only: mwsfc, itab_wsfc, sfcg
  use mem_land,   only: nzg
  use leaf_coms,  only: nzs
  use misc_coms,  only: io6, isubdomain, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: k
  integer :: notavail
  integer :: iwsfc, imsfc
  integer :: npoly

  real :: hpt, vpt, psiz, vsprd
  real :: fldval

  real :: htpn(maxnlspoly), vtpn(maxnlspoly)
  real :: wtbot = 1., wttop = 0.

  real :: area_tot, field_tot
  real :: areas(mgroupsize), fields(mgroupsize)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j, nr
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  area_tot = 0.
  field_tot   = 0.

  nu   = 0
  ipos = 0

  nr   = 2 * maxnlspoly + 5
  base = nr * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mwsfc) / 5. )
  else
     inc = mwsfc
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Find K level to plot if field is 3d

  if (trim(op%dimens) == '3G') then
     k = min(nzg,max(1,nint(op%slabloc(iplt))))
  elseif (trim(op%dimens) == '3S') then
     k = min(nzs,max(1,nint(op%slabloc(iplt))))
  else
     k = 1
  endif

  do iwsfc = 2, mwsfc

     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     if (op%stagpt == 'L' .and. sfcg%leaf_class(iwsfc) <  2) cycle
     if (op%stagpt == 'R' .and. sfcg%leaf_class(iwsfc) /= 1) cycle
     if (op%stagpt == 'S' .and. sfcg%leaf_class(iwsfc) /= 0) cycle

     npoly = itab_wsfc(iwsfc)%npoly

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,sfcg%xew(iwsfc),sfcg%yew(iwsfc),sfcg%zew(iwsfc),hpt,vpt)

     ! Initialize iflag180

     iflag180 = 0

     do j = 1,npoly      
        imsfc = itab_wsfc(iwsfc)%imn(j)

        call oplot_transform(iplt, sfcg%xem(imsfc), sfcg%yem(imsfc), sfcg%zem(imsfc), &
                             htpn(j), vtpn(j))

        ! Avoid wrap-around for lat-lon plot and set iflag180

        if (op%projectn(iplt) == 'L') then
           call ll_unwrap(hpt,htpn(j))
           if (htpn(j) < -180.001) iflag180 =  1
           if (htpn(j) >  180.001) iflag180 = -1
        endif
     enddo

     ! Jump out of loop if any cell corner is on other side of earth

     if (any(htpn(1:npoly) > 1.e11)) cycle

     ! Jump out of loop if entire cell is outside plot window.

     if ( all(htpn(1:npoly) < op%xmin) .or. all(htpn(1:npoly) > op%xmax) .or. &
          all(vtpn(1:npoly) < op%ymin) .or. all(vtpn(1:npoly) > op%ymax)) cycle

     ! Plot cell

     call oplot_lib(k,iwsfc,'VALUE',op%fldname(iplt),wtbot,wttop, &
                    fldval,notavail)    
     if (notavail > 0) cycle 

     area_tot  = area_tot  + sfcg%area(iwsfc)
     field_tot = field_tot + fldval * sfcg%area(iwsfc)

     call get_psiz(iplt,sqrt(sfcg%area(iwsfc)),psiz,vsprd)

     if (myrank == 0) then

        call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(npoly,  1,      MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(htpn,   npoly,  MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   npoly,  MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(psiz,   1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vsprd,  1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
     ! at other end

     if (iflag180 /= 0) then

        do j = 1,npoly
           htpn(j) = htpn(j) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(npoly,  1,      MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(htpn,   npoly,  MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   npoly,  MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif ! iflag180

  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then

     ! first collect the areas and totals to compute the average

     call MPI_Gather(area_tot,  1, MPI_REAL, areas,  1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Gather(field_tot, 1, MPI_REAL, fields, 1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     if (myrank == 0 ) then
        area_tot  = sum(areas)
        field_tot = sum(fields)
     endif

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, npoly,  1,      MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   npoly,  MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   npoly,  MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,      MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

  if (myrank == 0 .and. action == 'T') then
     if (area_tot > 1.) print*, trim(op%fldname(iplt)),' field_Savg ', &
                                field_tot, area_tot, field_tot/area_tot
  endif

end subroutine tileslab_horiz_wsfc

!===============================================================================

subroutine tileslab_horiz_vsfc(iplt,action)

  use oplot_coms, only: op
  use mem_sfcg,   only: mvsfc, itab_vsfc, sfcg
  use mem_land,   only: nzg
  use leaf_coms,  only: nzs
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: ivsfc
  integer :: iw1, iw2
  integer :: itpn
  integer :: im1, im2
  integer :: notavail

  real :: fldval
  real :: hpt, vpt, psiz, vsprd
  real :: htpn(4), vtpn(4)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  nu   = 0
  ipos = 0

  base = 13 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mvsfc) / 5. )
  else
     inc = mvsfc
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  do ivsfc = 2,mvsfc

     if (itab_vsfc(ivsfc)%irank /= myrank) cycle 

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,sfcg%xev(ivsfc),sfcg%yev(ivsfc),sfcg%zev(ivsfc),hpt,vpt)

     im1 = itab_vsfc(ivsfc)%imn(1)
     im2 = itab_vsfc(ivsfc)%imn(2)
     iw1 = itab_vsfc(ivsfc)%iwn(1)
     iw2 = itab_vsfc(ivsfc)%iwn(2)

     ! Initialize iflag180

     iflag180 = 0

     if (im1 > 1) then
        call oplot_transform(iplt,sfcg%xem(im1),sfcg%yem(im1),sfcg%zem(im1),htpn(1),vtpn(1))
     else
        htpn(1) = hpt
        vtpn(1) = vpt
     endif

     if (im2 > 1) then
        call oplot_transform(iplt,sfcg%xem(im2),sfcg%yem(im2),sfcg%zem(im2),htpn(3),vtpn(3))
     else
        htpn(3) = hpt
        vtpn(3) = vpt
     endif

     if (iw2 > 1) then
        call oplot_transform(iplt,sfcg%xew(iw2),sfcg%yew(iw2),sfcg%zew(iw2),htpn(2),vtpn(2))
     else
        htpn(2) = htpn(3)
        vtpn(2) = vtpn(3)
     endif

     if (iw1 > 1) then
        call oplot_transform(iplt,sfcg%xew(iw1),sfcg%yew(iw1),sfcg%zew(iw1),htpn(4),vtpn(4))
     else
        htpn(4) = htpn(3)
        vtpn(4) = vtpn(3)
     endif

     !  Avoid wrap-around for lat-lon plot

     if (op%projectn(iplt) == 'L') then
        do itpn = 1,4
           call ll_unwrap(hpt,htpn(itpn))
        enddo
     endif

     ! Jump out of loop if any cell corner is on other side of earth

     if (any(htpn(1:4) > 1.e11)) cycle

     ! Jump out of loop if entire cell is outside plot window. 

     if ( all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
          all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax)) cycle

     ! Get cell value and 'available' flag

     call oplot_lib(1,ivsfc,'VALUE',op%fldname(iplt),1.,0., &
                    fldval,notavail)    

     if (notavail > 0) cycle 

     ! Set iflag180
        
     if (op%projectn(iplt) == 'L') then
        do itpn = 1,4
           if (htpn(itpn) < -180.001) iflag180 =  1
           if (htpn(itpn) >  180.001) iflag180 = -1
        enddo
     endif

     call get_psiz(iplt,sfcg%dnv(ivsfc),psiz,vsprd)

     if (myrank == 0) then

        call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)
        
     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(psiz,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vsprd,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
     ! at other end

     if (iflag180 /= 0) then

        do itpn = 1, 4
           htpn(itpn) = htpn(itpn) + 360. * iflag180
        enddo

        if (myrank == 0) then

           call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     endif ! iflag180

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)
           
              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine tileslab_horiz_vsfc

!===============================================================================

subroutine tileslab_vert_tw(iplt,action)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mza, zm, zt, lpw
  use mem_ijtabs, only: jtab_w, jtw_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: k
  integer :: iw, jw
  integer :: iv1,iv2
  integer :: iok
  integer :: notavail

  real :: hpt,  vpt, psiz, vsprd
  real :: fldval

  real :: htpn(4), vtpn(4)
  real :: topo1,topo2
  real, parameter :: wtbot = 1., wttop = 0.

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  nu   = 0
  ipos = 0

  base = 13 * nbytes_real
  inc = ceiling( real(mwa) / 10. ) * mza

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over W points

  do jw = 1, jtab_w(jtw_prog)%jend(1)
     iw = jtab_w(jtw_prog)%iw(jw)

     ! Get horizontal plot coordinates for this W point

     if (op%projectn(iplt) == 'C') then
        call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
     elseif (op%projectn(iplt)(1:1) == 'V') then
        call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)
     endif

     hpt = .5 * (htpn(1) + htpn(2))

     ! Jump out of loop if this W point does not intersect plot cone

     if (iok /= 1) cycle

     ! Skip if entire cell is outside plot window

     if ( (htpn(1) < op%xmin .and. htpn(2) < op%xmin) .or.  &
          (htpn(1) > op%xmax .and. htpn(2) > op%xmax) ) cycle         

     do k = lpw(iw)-1, mza

        if (op%stagpt == 'T') then

           if (k == lpw(iw)-1) cycle

           vtpn(1) = zm(k-1)
           vtpn(2) = vtpn(1)
           vtpn(3) = zm(k)
           vtpn(4) = vtpn(3)
           vpt = zt(k)

        else ! stagpt == 'W'

           if (k == lpw(iw)-1) then
              vtpn(1) = zm(k)
           else
              vtpn(1) = zt(k)
           endif
           vtpn(2) = vtpn(1)
           if (k == mza) then
              vtpn(3) = zm(mza)
           else
              vtpn(3) = zt(k+1)
           endif
           vtpn(4) = vtpn(3)
           vpt = zm(k)

        endif

        ! Skip if entire cell is above or below plot window

        if ( all(vtpn(1:4) < op%ymin) .or. &
             all(vtpn(1:4) > op%ymax) ) cycle

        call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                       fldval,notavail)
        if (notavail > 0) cycle 
 
        psiz = op%psiz
        vsprd = op%vsprd

        if (myrank == 0) then

           call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(htpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(psiz,   1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vsprd,  1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(fldval, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
        endif
     enddo
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
               
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif
  
  ! Now plot underground T cells with underground color
  if (action == 'T') call plot_underground_w(iplt,(/0/))

end subroutine tileslab_vert_tw

!===============================================================================

subroutine tileslab_vert_v(iplt,action)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mza, zm, zt, lpw, lpv
  use mem_ijtabs, only: itab_v, jtab_w, jtw_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: k
  integer :: iw,jw
  integer :: iv1,iv2
  integer :: iok
  integer :: notavail

  real :: hptl,hptr,hpt,vpt, psiz, vsprd
  real :: fldvall, fldvalr, fldval
  real :: topo1, topo2

  real :: hcpn(4)
  real :: htpnl(4), htpnr(4), htpn(4), vtpn(4)

  real, parameter :: wtbot = 1., wttop = 0.

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  nu   = 0
  ipos = 0

  base = 13 * nbytes_real
  inc = ceiling( real(mwa) / 10. ) * mza

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Loop over W points

  do jw = 1, jtab_w(jtw_prog)%jend(1)
     iw = jtab_w(jtw_prog)%iw(jw)

     ! Get horizontal plot coordinates for IW point

     if (op%projectn(iplt) == 'C') then
        call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,hcpn)
     elseif (op%projectn(iplt)(1:1) == 'V') then
        call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,hcpn)
     endif

     ! Jump out of loop if this IW point does not intersect plot cone

     if (iok /= 1) cycle

     ! Avoid wrap_around

     if (op%projectn(iplt) == 'C') then
        call ll_unwrap(hcpn(1),hcpn(2))
     endif

     ! Skip current IV point if entire cell is outside plot window. 

     if ( (hcpn(1) < op%xmin .and. hcpn(2) < op%xmin) .or.  &
          (hcpn(1) > op%xmax .and. hcpn(2) > op%xmax) ) cycle         
   
     hcpn(3) = hcpn(2)
     hcpn(4) = hcpn(1)
   
     htpnl(1) = hcpn(1)
     htpnl(2) = .5 * (hcpn(1) + hcpn(2))
     htpnl(3) = htpnl(2)
     htpnl(4) = htpnl(1)
     hptl = hcpn(1)

     htpnr(1) = .5 * (hcpn(1) + hcpn(2))
     htpnr(2) = hcpn(2)
     htpnr(3) = htpnr(2)
     htpnr(4) = htpnr(1)
     hptr = op%xmax + 1 ! disable printing redundant right-side values

     do k = lpw(iw),mza  ! Loop is over T levels

        ! Skip this K point if entire cell center is above or below plot window

        if (zm(k) < op%ymin .or. zm(k-1) > op%ymax) cycle
   
        ! Get vertical coordinates

        vtpn(1) = zm(k-1)
        vtpn(2) = vtpn(1)
        vtpn(3) = zm(k)
        vtpn(4) = vtpn(3)

        if (k < lpv(iv1)) then
           vpt = op%ymin - 1. ! disable printing values below terrain
        else
           vpt = zt(k)
        endif

        ! Get value for left half of cell
   
        if (k >= lpv(iv1)) then

           ! extend lowest cell down slightly to avoid gaps in plotted terain
           if (k == lpv(iv1)) then
              vtpn(1) = zm(max(k-3,1))
              vtpn(2) = vtpn(1)
           else
              vtpn(1) = zm(k-1)
              vtpn(2) = vtpn(1)
           endif

           call oplot_lib(k,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                          fldvall,notavail)
           if (notavail > 0) cycle 

           psiz = op%psiz
           vsprd = op%vsprd

           if (myrank ==0) then
              call celltile(iplt,4,htpnl,vtpn,hptl,vpt,psiz,vsprd,fldvall,action)
           else
              nu = nu + 1
#ifdef OLAM_MPI
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(htpnl,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vtpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(hptl,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vpt,     1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(psiz,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vsprd,   1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvall, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
           endif

        endif
           
        ! Get value for right half of cell
        
        if (k >= lpv(iv2)) then

           ! extend lowest cell down slightly to avoid gaps in plotted terain
           if (k == lpv(iv2)) then
              vtpn(1) = zm(max(k-3,1))
              vtpn(2) = vtpn(1)
           else
              vtpn(1) = zm(k-1)
              vtpn(2) = vtpn(1)
           endif

           call oplot_lib(k,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                          fldvalr,notavail)
           if (notavail > 0) cycle 

           if (myrank ==0) then
              call celltile(iplt,4,htpnr,vtpn,hptr,vpt,psiz,vsprd,fldvalr,action)
           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(htpnr,   4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vtpn,    4, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(hptr,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vpt,     1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(psiz,    1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vsprd,   1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(fldvalr, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif
           
        endif
     enddo
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
               
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   4, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, psiz,   1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vsprd,  1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  ! Plot underground cells with underground color
  if (action == 'T') call plot_underground_w(iplt,(/0/))

end subroutine tileslab_vert_v

!===============================================================================

!!subroutine tileslab_vert_l(iplt,action)
!!
!!use oplot_coms, only: op
!!use leaf_coms,  only: nzs
!!use mem_land,   only: mland, land, nzg
!!use misc_coms,  only: io6
!!
!!implicit none
!!
!!integer,      intent(in) :: iplt
!!character(1), intent(in) :: action
!!
!!integer :: k
!!integer :: iw
!!integer :: ng
!!integer :: iv1, iv2
!!integer :: iu1, iu2
!!integer :: ip
!!integer :: itpn
!!integer :: notavail
!!
!!real :: hpt, vpt
!!real :: fldval
!!real :: patwidth
!!real :: botk
!!real :: delzk
!!
!!real :: htpn(4), vtpn(4)
!!real :: wtbot = 1., wttop = 0.
!!
!!!!!! VERTICAL XSECTION NEEDS WORK
!!
!!do ip = 2,mland
!!!   iw = land%iw(ip)
!!
!!! Skip iw column if not intersected by slabloc or if outside window bounds
!!
!!
!!! Get horizontal plot coordinates for cells in this column
!!
!!   op%stagpt = 'LA'  ! Get land cell fractional area
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!
!!!   patwidth = (xu(iu2) - xu(iu1)) * fldval
!!
!!!   htpn(1) = vt2da(iw)
!!!   htpn(2) = htpn(1) + patwidth
!!!   htpn(3) = htpn(2)
!!!   htpn(4) = htpn(1)
!!
!!!   vt2da(iw) = htpn(2)
!!
!!!   hpt = .5 * (htpn(1) + htpn(2))
!!
!!   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
!!   delzk = op%ymin / botk  ! This is a positive number
!!
!!   do k = 1,nzg  ! Loop over soil layers
!!
!!! Get vertical coordinates for soil layers
!!
!!      vtpn(1) = delzk * (float(k-1) + botk) 
!!      vtpn(2) = vtpn(1)
!!      vtpn(3) = delzk * (float(k)   + botk)
!!      vtpn(4) = vtpn(3)
!!      vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!! plot soil layers
!!
!!      op%stagpt = 'L'  ! Get soil value
!!      call oplot_lib(k,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                     fldval,notavail)
!!      call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!      if (op%pltgrid(iplt) == 'G') then
!!         call o_frstpt(htpn(4),vtpn(4))
!!         do itpn = 1,4
!!            call o_vector(htpn(itpn),vtpn(itpn))
!!         enddo
!!      endif
!!
!!   enddo
!!
!!   do k = 1,nzs  ! Loop over sfcwater layers
!!
!!! check for existence of sfcwater in current level k
!!
!!      call oplot_lib(k,ip,'VALUE','SFCWATER_MASS',wtbot,wttop, &
!!                     fldval,notavail)
!!
!!      if (fldval > 0.) then
!!      
!!! Get vertical coordinates for sfcwater layer k
!!
!!         vtpn(1) = delzk * (float(nzg+k-1) + .5 + botk) 
!!         vtpn(2) = vtpn(1)
!!         vtpn(3) = delzk * (float(nzg+k)   + .5 + botk) 
!!         vtpn(4) = vtpn(3)
!!         vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!! plot sfcwater layer k
!!
!!         op%stagpt = 'LW'  ! Get sfcwater value
!!         call oplot_lib(k,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                        fldval,notavail)
!!         call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!         if (op%pltgrid(iplt) == 'G') then
!!            call o_frstpt(htpn(4),vtpn(4))
!!            do itpn = 1,4
!!               call o_vector(htpn(itpn),vtpn(itpn))
!!            enddo
!!         endif
!!         
!!      endif
!!      
!!   enddo
!!
!!! plot vegetation layer
!!
!!   vtpn(1) = delzk * (float(nzg+nzs) + 1. + botk) 
!!   vtpn(2) = vtpn(1)
!!   vtpn(3) = delzk * (float(nzg+nzs+1) + 1. + botk) 
!!   vtpn(4) = vtpn(3)
!!   vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!   op%stagpt = 'LV'  ! Get vegetation value
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!   call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!   if (op%pltgrid(iplt) == 'G') then
!!      call o_frstpt(htpn(4),vtpn(4))
!!      do itpn = 1,4
!!         call o_vector(htpn(itpn),vtpn(itpn))
!!      enddo
!!   endif
!!
!!! plot canopy air layer
!!
!!   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
!!   vtpn(2) = vtpn(1)
!!   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
!!   vtpn(4) = vtpn(3)
!!   vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!   op%stagpt = 'LC'  ! Get vegetation value
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!   call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!   if (op%pltgrid(iplt) == 'G') then
!!      call o_frstpt(htpn(4),vtpn(4))
!!      do itpn = 1,4
!!         call o_vector(htpn(itpn),vtpn(itpn))
!!      enddo
!!   endif
!!
!!enddo
!!   
!!return
!!end subroutine tileslab_vert_l

!===============================================================================

!!subroutine tileslab_vert_s(iplt,action)
!!
!!use oplot_coms, only: op
!!use leaf_coms,  only: nzs
!!use mem_land,   only: nzg
!!use mem_sea,    only: msea, sea
!!use misc_coms,  only: io6
!!
!!implicit none
!!
!!integer,      intent(in) :: iplt
!!character(1), intent(in) :: action
!!
!!integer :: k
!!integer :: iw
!!integer :: ng
!!integer :: iv1, iv2
!!integer :: ip
!!integer :: itpn
!!integer :: notavail
!!
!!real :: hpt, vpt
!!real :: fldval
!!real :: patwidth
!!real :: botk
!!real :: delzk
!!
!!real :: htpn(4), vtpn(4)
!!real :: wtbot = 1., wttop = 0.
!!
!!!!!! VERTICAL XSECTION NEEDS WORK
!!
!!do ip = 2,msea
!!!   iw = sea%iw(ip)
!!
!!! Skip iw column if not intersected by slabloc or if outside window bounds
!!
!!! Get horizontal plot coordinates for cells in this column
!!
!!   op%stagpt = 'SA'  ! Get sea cell fractional area
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!
!!   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
!!   delzk = op%ymin / botk  ! This is a positive number
!!
!!! plot (top) sea layer
!!
!!   vtpn(1) = delzk * (float(nzg-1) + botk) 
!!   vtpn(2) = vtpn(1)
!!   vtpn(3) = delzk * (float(nzg) + botk) 
!!   vtpn(4) = vtpn(3)
!!   vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!   op%stagpt = 'S'  ! Get vegetation value
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!   call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!   if (op%pltgrid(iplt) == 'G') then
!!      call o_frstpt(htpn(4),vtpn(4))
!!      do itpn = 1,4
!!         call o_vector(htpn(itpn),vtpn(itpn))
!!      enddo
!!   endif
!!
!!! plot canopy air layer
!!
!!   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
!!   vtpn(2) = vtpn(1)
!!   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
!!   vtpn(4) = vtpn(3)
!!   vpt = .5 * (vtpn(1) + vtpn(3))
!!
!!   op%stagpt = 'SC'  ! Get canopy air value
!!   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
!!                  fldval,notavail)
!!   call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'
!!
!!   if (op%pltgrid(iplt) == 'G') then
!!      call o_frstpt(htpn(4),vtpn(4))
!!      do itpn = 1,4
!!         call o_vector(htpn(itpn),vtpn(itpn))
!!      enddo
!!   endif
!!
!!enddo
!!   
!!return
!!end subroutine tileslab_vert_s

!===============================================================================

subroutine celltile(iplt,npoly,htpn,vtpn,hpt,vpt,psiz,vsprd,fldval,action)

  use oplot_coms, only: op
  use plotcolors, only: clrtab
  use misc_coms,  only: io6

  implicit none

  integer, intent(in) :: iplt
  integer, intent(in) :: npoly

  real, intent(in) :: htpn(*)
  real, intent(in) :: vtpn(*)
  real, intent(in) :: hpt
  real, intent(in) :: vpt
  real, intent(in) :: psiz
  real, intent(in) :: vsprd
  real, intent(in) :: fldval

  character(1), intent(in) :: action

  real    :: fldval1
  integer :: icolor
  integer :: itab
  integer :: ival

  if (action == 'T') then

     itab = op%icolortab(iplt)

     ! Cyclic treatment of color palette (used with integer-type data)

     fldval1 = fldval

     if (clrtab(itab)%ifmt(1) == 20) &
!         fldval1 = mod(fldval-1.,real(clrtab(itab)%nvals-2)) - 1. ! table 124 has 2 neg vals
          fldval1 = mod(fldval,real(clrtab(itab)%nvals-2)) - 1. ! table 124 has 2 neg vals

     ! Case for color table 130

     if (clrtab(itab)%ifmt(1) == 30) &
          fldval1 = mod(fldval,100.)

     ! Extract contour color from color table

     ival = 1
     do while (fldval1 > clrtab(itab)%vals(ival) .and. &
               ival < clrtab(itab)%nvals               )
        ival = ival + 1
     enddo
     icolor = clrtab(itab)%ipal(ival)
   
     call fillpolyg(npoly,htpn,vtpn,icolor)

     if (op%fldval_min > fldval) op%fldval_min = fldval
     if (op%fldval_max < fldval) op%fldval_max = fldval

  elseif (action == 'P') then
     if ( (hpt > op%xmin) .and. (hpt < op%xmax) .and. &
          (vpt > op%ymin) .and. (vpt < op%ymax) ) then
        call oplot_prtvalue(fldval,hpt,vpt,vsprd,psiz,op%icolortab(iplt))
     endif
  endif

end subroutine celltile

!===============================================================================

subroutine xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)

  ! This subroutine finds points of intersection between vertical plot slab and
  ! hexagonal grid cells

  use oplot_coms,  only: op, xepc, yepc
  use mem_grid,    only: xem, yem, topm
  use mem_ijtabs,  only: itab_w
  use consts_coms, only: pio180
  use oname_coms,  only: nl
  use misc_coms,  only: io6

  implicit none

  integer, intent(in)  :: iplt
  integer, intent(in)  :: iw
  integer, intent(out) :: iv1,iv2
  integer, intent(out) :: iok

  real, intent(out) :: htpn(4)
  real, intent(out) :: topo1,topo2

  integer :: npoly,j,jm1,jm2,jmin,jv,im,im1,im2,iv

  real :: smin,smax
  real :: wt1,wt2
  real :: sinvaz,cosvaz

  real :: sm(7),sm1,sm2
  real :: topm1,topm2 ! Topo height of 2 M points
  real :: xem1, xem2, yem1, yem2

  sinvaz = sin((90. - op%viewazim) * pio180)
  cosvaz = cos((90. - op%viewazim) * pio180)

  smin =  1.e9
  smax = -1.e9

  npoly = itab_w(iw)%npoly

! Loop over neighbor M points for this IW cell

  do j = 1,npoly
     im = itab_w(iw)%im(j)

    sm(j) = (xem(im) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
          + (yem(im) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

    if (smin > sm(j)) then
       jmin = j
       smin = sm(j) 
    endif

    if (smax < sm(j)) then
       smax = sm(j)
    endif
  enddo

  ! Return with iok = 0 if iw column is not intersected by plot slab

  iok = 0

  if (smin > op%slabloc(iplt) .or. smax < op%slabloc(iplt)) return

  iok = 1  ! Since we got here, IW column is intersected by plot cone

! Fill arrays of values at M points in cyclic order around IW column, beginning
! with M point that is closest to cone axis

  do j = 1,npoly
     jm1 = j + jmin - 1
     jm2 = jm1 + 1

     if (jm1 > npoly) jm1 = jm1 - npoly
     if (jm2 > npoly) jm2 = jm2 - npoly
   
     im1 = itab_w(iw)%im(jm1)
     im2 = itab_w(iw)%im(jm2)

     jv = jm2
     iv  = itab_w(iw)%iv(jv)
      
     xem1 = xem(im1)
     yem1 = yem(im1)

     xem2 = xem(im2)
     yem2 = yem(im2)
   
     topm1 = topm(im1)
     topm2 = topm(im2)
   
     sm1 = sm(jm1)
     sm2 = sm(jm2)

! Find two points of intersection between current IW polygon and slab

     if (sm1 <= op%slabloc(iplt) .and. sm2 >= op%slabloc(iplt)) then

! This interval touches slab

        if (sm1 == sm2) then

! This interval is tangent to slab      

           if (jm1 == 1 .or. jm2 == 1) then

! This interval is on minimum side of polygon

              xepc(1) = xem1
              yepc(1) = yem1
              topo1   = topm1
              iv1     = iv

              xepc(2) = xem2
              yepc(2) = yem2
              topo2   = topm2
              iv2     = iv

           else

! This interval is on maximum side of polygon

              xepc(1) = xem2
              yepc(1) = yem2
              topo1   = topm2
              iv1     = iv

              xepc(2) = xem1
              yepc(2) = yem1
              topo2   = topm1
              iv2     = iv

           endif

           exit

        else
         
! This interval touches slab at 1 point

           wt2 = (op%slabloc(iplt) - sm1) / (sm2 - sm1)
           wt1 = 1. - wt2

           xepc(2) = wt1 * xem1  + wt2 * xem2
           yepc(2) = wt1 * yem1  + wt2 * yem2
           topo2   = wt1 * topm1 + wt2 * topm2
           iv2     = iv

        endif
      
     elseif (sm1 > op%slabloc(iplt) .and. sm2 <= op%slabloc(iplt)) then

! This interval touches slab at 1 point

        wt2 = (op%slabloc(iplt) - sm2) / (sm1 - sm2)
        wt1 = 1. - wt2

        xepc(1) = wt1 * xem2  + wt2 * xem1
        yepc(1) = wt1 * yem2  + wt2 * yem1
        topo1   = wt1 * topm2 + wt2 * topm1
        iv1     = iv

     endif

  enddo

! Transform horizontal point coordinates

  htpn(1) = (xepc(1) - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
          - (yepc(1) - nl%plotspecs(iplt)%plotcoord2) * cosvaz

  htpn(2) = (xepc(2) - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
          - (yepc(2) - nl%plotspecs(iplt)%plotcoord2) * cosvaz
           
  htpn(3) = htpn(2)
  htpn(4) = htpn(1)

end subroutine xyplot_w
