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
                        xev, yev, zev, xew, yew, zew
  use mem_ijtabs, only: itab_m, jtab_m, jtm_vadj
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: j, jm, jn, jnn, im, iw, npoly, iv1, iv2
  integer :: notavail, navail

  real :: hpt, vpt
  real :: fldval

  real :: htpn(7), vtpn(7)
  real :: hqpn(4), vqpn(4)

  integer :: ktf(mwa), km(mma)
  real :: wtbot(mma), wttop(mma)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (action == 'T' .and. op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 17 * nbytes_real + nbytes_int
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

  mloop: do jm = 1, jtab_m(jtm_vadj)%jend(1)
            im = jtab_m(jtm_vadj)%im(jm)

     ! For now, skip pts that don't read in topm

     if (.not. itab_m(im)%loop(1)) cycle mloop

     ! Get tile plot coordinates.

     call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

     npoly = itab_m(im)%npoly

     ! Initialize navail counter

     navail = 0

     do j = 1,npoly

        iw = itab_m(im)%iw(j)

        ! Skip this M point if current IW point index < 2 (which occurs at
        ! lateral boundary of limited-area domain or parallel subdomain)      

        if (iw < 2) cycle mloop

        call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(j),vtpn(j))

        ! Avoid wrap-around for lat-lon plot

        if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))

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

        if (myrank == 0) then

           call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,fldval,action)

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
           call MPI_Pack(fldval, 1,     MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If any surrounding T cell is unavailable, plot M cell by avalable sectors

     else

        hqpn(1) = hpt
        vqpn(1) = vpt

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

              if (myrank == 0) then

                 call celltile(iplt,4,hqpn,vqpn,hpt,vpt,fldval,action)

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
                 call MPI_Pack(fldval, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

              endif

           endif
        enddo
     endif

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
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,     MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,fldval,action)

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
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
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

  real :: hpt, vpt
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

  base = 22 * nbytes_real + 2 * nbytes_int
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

  wloop: do jw = 1, jtab_w(jtw_prog)%jend(1)
     iw = jtab_w(jtw_prog)%iw(jw)

     npoly = itab_w(iw)%npoly

     ! For hexagon grid (in limited-area domain), do not tile plot
     ! incomplete boundary cells

     if (action == 'T' .and. npoly < 5) cycle wloop

     ! Skip this point if it is underground

     if (ktf(iw) /= 0) cycle wloop

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

     do j = 1,npoly

        im = itab_w(iw)%im(j)
      
        call oplot_transform(iplt,xem(im),yem(im),zem(im),htpn(j),vtpn(j))

        ! Avoid wrap-around for lat-lon plot

        if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))

        ! Jump out of loop if cell corner is on other side of earth

        if (htpn(j) > 1.e11) cycle wloop

     enddo

     ! Jump out of loop if entire cell is outside plot window

     if ( all(htpn(1:npoly) < op%xmin) .or. all(htpn(1:npoly) > op%xmax) .or.  &
          all(vtpn(1:npoly) < op%ymin) .or. all(vtpn(1:npoly) > op%ymax) ) cycle wloop

     ! Get cell value and plot if 'available'

     call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                    fldval,notavail)

     if (notavail > 0) cycle wloop

     arw0_tot  = arw0_tot  + arw0(iw)
     field_tot = field_tot + fldval * arw0(iw)

     if (myrank == 0) then
        call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,fldval,action)
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
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,     MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,npoly,htpn,vtpn,hpt,vpt,fldval,action)

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
                        xem, yem, zem, xew, yew, zew, lpw
  use mem_ijtabs, only: itab_v, jtab_v, jtv_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: iv, jv
  integer :: iw1, iw2
  integer :: itpn
  integer :: im1, im2
  integer :: notavail

  real :: fldval
  real :: hpt, vpt
  real :: htpn(4), vtpn(4)
  real :: htpn2, htpn4
  real :: vtpn2, vtpn4

  integer :: ktf(mwa), kv(mva)
  real :: wtbot(mva), wttop(mva)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Find cell K indices on the given plot surface

  call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

  ! If field is 3d, first plot underground points with underground color

  if (action == 'T' .and. op%dimens == '3') then
     call plot_underground_w(iplt,ktf)
  endif

  nu   = 0
  ipos = 0

  base = 11 * nbytes_real
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

  do jv = 1, jtab_v(jtv_prog)%jend(1)

     ! Transform tile plot X and Y coordinates.  

     iv  = jtab_v(jtv_prog)%iv(jv)
     call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hpt,vpt)

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! Skip this V point if both its W neighbors are below ground

     if (ktf(iw1) /= 0 .and. ktf(iw2) /= 0) cycle

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

     ! Save copy of 2nd and 4th htpn,vtpn pts

     htpn2 = htpn(2)
     htpn4 = htpn(4)
     vtpn2 = vtpn(2)
     vtpn4 = vtpn(4)
   
     ! Get cell value and 'available' flag

     call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                    fldval,notavail)    

     if (notavail > 0) cycle 

     ! If both W neighbors are above ground, plot full cell

     if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

        if (myrank == 0) then
           call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)
        endif
      
        ! Else, check if IW1 is above ground

     elseif (ktf(iw1) == 0) then

        ! Plot IW1 half of cell with tile color

        htpn(4) = htpn4
        vtpn(4) = vtpn4

        htpn(2) = htpn(3)
        vtpn(2) = vtpn(3)

        if (myrank == 0) then
           call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)
        endif

     else

        ! Plot IW2 half of cell with tile color

        htpn(2) = htpn2
        vtpn(2) = vtpn2

        htpn(4) = htpn(3)
        vtpn(4) = vtpn(3)

        if (myrank == 0) then
           call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)
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

        call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)

     endif
#endif

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
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)
           
              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

end subroutine tileslab_horiz_vn

!===============================================================================

subroutine tileslab_horiz_s(iplt,action)

  use max_dims,   only: maxnlspoly
  use oplot_coms, only: op
  use mem_sea,    only: sea, itab_ws, itabg_ws
  use sea_coms,   only: mws, iseagrid
  use misc_coms,  only: io6, isubdomain, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer,      intent(in) :: iplt
  character(1), intent(in) :: action

  integer :: notavail
  integer :: iws
  integer :: nspoly, jms

  real :: hpt, vpt
  real :: fldval

  real :: htpn(maxnlspoly), vtpn(maxnlspoly)
  real :: wtbot = 1., wttop = 0.

  real :: areasea_tot, field_tot
  real :: areas(mgroupsize), fields(mgroupsize)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j, nr
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  areasea_tot = 0.
  field_tot   = 0.

  nu   = 0
  ipos = 0

  nr   = 2 * maxnlspoly + 3
  base = nr * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mws) / 5. )
  else
     inc = mws
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  do iws = 2, mws

     ! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

     if (isubdomain == 1 .and. itab_ws(iws)%irank /= myrank) cycle

     nspoly = itab_ws(iws)%npoly

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,sea%xew(iws),sea%yew(iws),sea%zew(iws),hpt,vpt)

     do jms = 1,nspoly      
        call oplot_transform(iplt, itab_ws(iws)%xem(jms), itab_ws(iws)%yem(jms), &
                                   itab_ws(iws)%zem(jms), htpn(jms), vtpn(jms))

        ! Avoid wrap-around for lat-lon plot

        if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(jms))
     enddo

     ! Jump out of loop if any cell corner is on other side of earth

     if (any(htpn(1:nspoly) > 1.e11)) cycle

     ! Jump out of loop if entire cell is outside plot window.

     if ( all(htpn(1:nspoly) < op%xmin) .or. all(htpn(1:nspoly) > op%xmax) .or. &
          all(vtpn(1:nspoly) < op%ymin) .or. all(vtpn(1:nspoly) > op%ymax)) cycle

     ! Plot cell

     call oplot_lib(1,iws,'VALUE',op%fldname(iplt),wtbot,wttop, &
                    fldval,notavail)    
     if (notavail > 0) cycle 

     areasea_tot = areasea_tot + sea%area(iws)
     field_tot   = field_tot   + fldval * sea%area(iws)

     if (myrank == 0) then

        call celltile(iplt,nspoly,htpn,vtpn,hpt,vpt,fldval,action)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(nspoly, 1,      MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(htpn,   nspoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   nspoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then

     ! first collect the areas and totals to compute the average

     call MPI_Gather(areasea_tot, 1, MPI_REAL, areas,  1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Gather(field_tot,   1, MPI_REAL, fields, 1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     if (myrank == 0 ) then
        areasea_tot = sum(areas)
        field_tot   = sum(fields)
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
               
                 call MPI_Unpack(buffer, buffsize, ipos, nspoly, 1,      MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   nspoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   nspoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,      MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,nspoly,htpn,vtpn,hpt,vpt,fldval,action)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

  if (myrank == 0 .and. action == 'T') then
     if (areasea_tot > 1.) print*, trim(op%fldname(iplt)),' field_Savg ', &
                                   field_tot, areasea_tot, field_tot/areasea_tot
  endif

end subroutine tileslab_horiz_s

!===============================================================================

subroutine tileslab_horiz_l(iplt,action)

  use max_dims,   only: maxnlspoly
  use oplot_coms, only: op
  use mem_leaf,   only: land, itab_wl, itabg_wl
  use leaf_coms,  only: mwl, nzg, nzs
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
  integer :: iwl
  integer :: nlpoly, jml

  real :: hpt, vpt
  real :: fldval

  real :: htpn(maxnlspoly), vtpn(maxnlspoly)
  real :: wtbot = 1., wttop = 0.

  real :: arealand_tot, field_tot
  real :: areas(mgroupsize), fields(mgroupsize)

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j, nr
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  arealand_tot = 0.
  field_tot    = 0.

  nu   = 0
  ipos = 0

  nr   = 2 * maxnlspoly + 3
  base = nr * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mwl) / 5. )
  else
     inc = mwl
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  ! Find K level to plot if field is 3d

  if (op%dimens == '3G') then
     k = min(nzg,max(1,nint(op%slabloc(iplt))))
  elseif (op%dimens == '3S') then
     k = min(nzs,max(1,nint(op%slabloc(iplt))))
  else
     k = 1
  endif

  do iwl = 2, mwl

     ! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

     if (isubdomain == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

     nlpoly = itab_wl(iwl)%npoly

     ! Get tile plot coordinates.  

     call oplot_transform(iplt,land%xew(iwl),land%yew(iwl),land%zew(iwl),hpt,vpt)

     do jml = 1,nlpoly
        call oplot_transform(iplt, itab_wl(iwl)%xem(jml), itab_wl(iwl)%yem(jml), &
                                   itab_wl(iwl)%zem(jml), htpn(jml), vtpn(jml))

        ! Avoid wrap-around for lat-lon plots

        if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(jml))
     enddo

     ! Jump out of loop if any cell corner is on other side of earth

     if (any(htpn(1:nlpoly) > 1.e11)) cycle

     ! Jump out of loop if entire cell is outside plot window. 

     if ( all(htpn(1:nlpoly) < op%xmin) .or. all(htpn(1:nlpoly) > op%xmax) .or. &
          all(vtpn(1:nlpoly) < op%ymin) .or. all(vtpn(1:nlpoly) > op%ymax) ) cycle

     ! Plot cell

     call oplot_lib(k,iwl,'VALUE',op%fldname(iplt),wtbot,wttop, &
                    fldval,notavail)    
     if (notavail > 0) cycle 


     arealand_tot = arealand_tot + land%area(iwl)
     field_tot    = field_tot + fldval * land%area(iwl)

     if (myrank == 0) then

        call celltile(iplt,nlpoly,htpn,vtpn,hpt,vpt,fldval,action)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + base) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(nlpoly, 1,      MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(htpn,   nlpoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   nlpoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(hpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vpt,    1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(fldval, 1,      MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then

     ! first collect the areas and totals to compute the average

     call MPI_Gather(arealand_tot, 1, MPI_REAL, areas,  1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Gather(field_tot,    1, MPI_REAL, fields, 1, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     if (myrank == 0 ) then
        arealand_tot = sum(areas)
        field_tot    = sum(fields)
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
               
                 call MPI_Unpack(buffer, buffsize, ipos, nlpoly, 1,      MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   nlpoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   nlpoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt,    1,      MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1,      MPI_REAL,    MPI_COMM_WORLD, ier)

                 call celltile(iplt,nlpoly,htpn,vtpn,hpt,vpt,fldval,action)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

  if (myrank == 0 .and. action == 'T') then
     if (arealand_tot > 1.) print*, trim(op%fldname(iplt)),' field_Lavg ', &
                                    field_tot, arealand_tot, field_tot/arealand_tot
  endif

end subroutine tileslab_horiz_l

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

  real :: hpt,  vpt
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

  base = 11 * nbytes_real
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

        ! The following call was designed for triangle grid; 
        ! need to make xyplot_w version for hexagons   
        ! call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)

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
 
        if (myrank == 0) then

           call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)

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
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)

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

  real :: hptl,hptr,hpt,vpt
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

  base = 11 * nbytes_real
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
        call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,hcpn) ! need to fix for hex??
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

           if (myrank ==0) then
              call celltile(iplt,4,htpnl,vtpn,hptl,vpt,fldvall,action)
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
              call celltile(iplt,4,htpnr,vtpn,hptr,vpt,fldvalr,action)
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
                 call MPI_Unpack(buffer, buffsize, ipos, fldval, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call celltile(iplt,4,htpn,vtpn,hpt,vpt,fldval,action)

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
!!use leaf_coms,  only: mwl, nzg, nzs
!!use mem_leaf,   only: land
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
!!do ip = 2,mwl
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
!!use leaf_coms,  only: nzg, nzs
!!use sea_coms,   only: mws
!!use mem_sea,    only: sea
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
!!do ip = 2,mws
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

subroutine celltile(iplt,npoly,htpn,vtpn,hpt,vpt,fldval,action)

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

     ! Case for color table 130; used with itab_w_mrowh

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
        call oplot_prtvalue(fldval,hpt,vpt,op%vsprd,.7*op%psiz,op%icolortab(iplt))
     endif
  endif

end subroutine celltile

!===============================================================================

subroutine xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)

  ! This subroutine finds points of intersection between vertical plot slab and
  ! triangular (prism) grid cells

  use oplot_coms,  only: op
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

  integer :: im1,im2,im3
  integer :: iu1,iu2,iu3
  integer :: iuc1,iuc2,iuc3

  real :: s1,s2,s3
  real :: sc1,sc2,sc3
  real :: smin,smax
  real :: wt1,wt2
  real :: x1,x2
  real :: y1,y2
  real :: sinvaz,cosvaz

  real :: xemc(3),yemc(3)  ! Coords of 3 M points
  real :: topmc(3) ! Topo height of 3 M points

  im1 = itab_w(iw)%im(1)
  im2 = itab_w(iw)%im(2)
  im3 = itab_w(iw)%im(3)

! iu1 = itab_w(iw)%iu(1)  [%iu not defined in hex grid]
! iu2 = itab_w(iw)%iu(2)
! iu3 = itab_w(iw)%iu(3)

  sinvaz = sin((90. - op%viewazim) * pio180)
  cosvaz = cos((90. - op%viewazim) * pio180)

  ! Location of 3 M points along line perpendicular to plot slab

  s1 = (xem(im1) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
     + (yem(im1) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

  s2 = (xem(im2) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
     + (yem(im2) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

  s3 = (xem(im3) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
     + (yem(im3) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

  smin = min(s1,s2,s3)
  smax = max(s1,s2,s3)

  ! Return with iok = 0 if iw column is not intersected by plot slab

  iok = 0

  if (smin > op%slabloc(iplt) .or. smax < op%slabloc(iplt)) return

  iok = 1  ! Since we got here, IW column is intersected by plot cone

  ! Find IM point that has lowest S coordinate and copy to temporary M point #1
  ! Fill other 2 points in cyclic order

  if (s1 <= s2 .and. s1 <= s3) then  ! m1 has lowest S value

     xemc(1)  = xem(im1)
     yemc(1)  = yem(im1)
     topmc(1) = topm(im1)

     xemc(2)  = xem(im2)
     yemc(2)  = yem(im2)
     topmc(2) = topm(im2)

     xemc(3)  = xem(im3)
     yemc(3)  = yem(im3)
     topmc(3) = topm(im3)

     sc1 = s1
     sc2 = s2
     sc3 = s3

     iuc1 = iu1
     iuc2 = iu2
     iuc3 = iu3

  elseif (s2 <= s1 .and. s2 <= s3) then  ! m2 has lowest S value

     xemc(1)  = xem(im2)
     yemc(1)  = yem(im2)
     topmc(1) = topm(im2)

     xemc(2)  = xem(im3)
     yemc(2)  = yem(im3)
     topmc(2) = topm(im3)

     xemc(3)  = xem(im1)
     yemc(3)  = yem(im1)
     topmc(3) = topm(im1)

     sc1 = s2
     sc2 = s3
     sc3 = s1

     iuc1 = iu2
     iuc2 = iu3
     iuc3 = iu1

  elseif (s3 <= s1 .and. s3 <= s2) then  ! m3 has lowest S value

     xemc(1)  = xem(im3)
     yemc(1)  = yem(im3)
     topmc(1) = topm(im3)

     xemc(2)  = xem(im1)
     yemc(2)  = yem(im1)
     topmc(2) = topm(im1)

     xemc(3)  = xem(im2)
     yemc(3)  = yem(im2)
     topmc(3) = topm(im2)

     sc1 = s3
     sc2 = s1
     sc3 = s2

     iuc1 = iu3
     iuc2 = iu1
     iuc3 = iu2

  endif
   
  ! Find two points of intersection between current IW triangle and X slab

  if (sc2 > op%slabloc(iplt)) then

     wt2 = (op%slabloc(iplt) - sc1) / (sc2 - sc1)
     wt1 = 1. - wt2

     x2    = wt1 * xemc(1)  + wt2 * xemc(2)  ! x coord of htpn(2)
     y2    = wt1 * yemc(1)  + wt2 * yemc(2)  ! y coord of htpn(2)
     topo2 = wt1 * topmc(1) + wt2 * topmc(2) ! topo height(2)

     htpn(2) = (x2 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
             - (y2 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
           
     iv2 = iuc3
   
     if (sc3 > op%slabloc(iplt)) then

        wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
        wt1 = 1. - wt2

        x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
        y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
        topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

        htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
                - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
        iv1 = iuc2
      
     else

        wt2 = (op%slabloc(iplt) - sc3) / (sc2 - sc3)
        wt1 = 1. - wt2

        x1    = wt1 * xemc(3)  + wt2 * xemc(2)  ! x coord of htpn(1)
        y1    = wt1 * yemc(3)  + wt2 * yemc(2)  ! y coord of htpn(1)
        topo1 = wt1 * topmc(3) + wt2 * topmc(2) ! topo height(1)

        htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
                - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
        iv1 = iuc1
      
     endif

  else

     wt2 = (op%slabloc(iplt) - sc2) / (sc3 - sc2)
     wt1 = 1. - wt2

     x2    = wt1 * xemc(2)  + wt2 * xemc(3)  ! x coord of htpn(2)
     y2    = wt1 * yemc(2)  + wt2 * yemc(3)  ! y coord of htpn(2)
     topo2 = wt1 * topmc(2) + wt2 * topmc(3) ! topo height(2)

     htpn(2) = (x2 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
            - (y2 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
     iv2 = iuc1

     wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
     wt1 = 1. - wt2

     x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
     y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
     topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

     htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
             - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz

     iv1 = iuc2

  endif

  htpn(3) = htpn(2)
  htpn(4) = htpn(1)

end subroutine xyplot_w
