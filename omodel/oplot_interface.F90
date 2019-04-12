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
subroutine oplot_init()

  use oplot_coms, only: op
  use plotcolors, only: gks_colors
  use misc_coms,  only: io6
  use mem_para,   only: myrank

  implicit none

  if (myrank == 0) then
     call o_opngks()
     call gks_colors(1)
  endif

  ! Initialize some oplot parameters (not set in namelist).

  op%loopplot  = 0
  op%icigrnd   = 13
  op%iplotback = 0

end subroutine oplot_init

!===============================================================================

subroutine plot_fields(id)

  use misc_coms,  only: io6, iyear1, imonth1, idate1, itime1, time_istp8
  use oplot_coms, only: op
  use mem_grid,   only: mza, mwa, lpw
  use misc_coms,  only: runtype, iparallel
  use mem_ijtabs, only: jtab_w, jtw_prog
  use mem_micro,  only: pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                        accpd, accpr, accpp, accps, accpa, accpg, accph
  use mem_cuparm, only: conprr, aconpr
  use obnd,       only: lbcopy_w1d
  use mem_para,   only: myrank, mgroupsize

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: id

  integer :: iplt,labincx,labincy,notavail,k,iw,iok
  real    :: fldval,bsize,dum1(1),dum2(1)
  integer :: outyear,outmonth,outdate,outhour
  real    :: xinc,yinc
  character(len=30) :: ylabel,title

  integer :: notavails(mgroupsize)
  dum1 = 0.
  dum2 = 0.

  ! Lateral boundary copy of surface precipitation quantities 
  ! (only needed for plotting)

  if (allocated(pcprd)) call lbcopy_w1d(1,a1=pcprd,d1=accpd)
  if (allocated(pcprr)) call lbcopy_w1d(1,a1=pcprr,d1=accpr)
  if (allocated(pcprp)) call lbcopy_w1d(1,a1=pcprp,d1=accpp)
  if (allocated(pcprs)) call lbcopy_w1d(1,a1=pcprs,d1=accps)
  if (allocated(pcpra)) call lbcopy_w1d(1,a1=pcpra,d1=accpa)
  if (allocated(pcprg)) call lbcopy_w1d(1,a1=pcprg,d1=accpg)
  if (allocated(pcprh)) call lbcopy_w1d(1,a1=pcprh,d1=accph)

  if (allocated(conprr)) call lbcopy_w1d(1,a1=conprr,d1=aconpr)

  ! Reopen the current graphics output workstation if it is closed

  if (myrank == 0) call o_reopnwk()

  do iplt = 1,op%nplt

     ! If external plot field is allocated, plot only the member(s) of iplt loop
     ! that match the external field name and are listed as external 'e' fields

     if (allocated(op%extfld) .and. op%fldname(iplt) /= op%extfldname) cycle

     if (allocated(op%extfld) .and. op%ext(iplt) /= 'e') cycle

     ! If external plot field is not allocated, plot only the member(s) of iplt
     ! loop that are NOT listed as external 'e' fields

     if (.not. allocated(op%extfld) .and. op%ext(iplt) == 'e') cycle

     ! Plot a white background for this frame

     if (op%iplotback /= 1 .and. myrank == 0) then
        call plotback()
     endif

     ! Draw filled-in map if so specified

!    if (( op%projectn(iplt) == 'L'  .or.   &
!          op%projectn(iplt) == 'P'  .or.   &
!          op%projectn(iplt) == 'G'  .or.   &
!          op%projectn(iplt) == 'O' ).and.  &
!         (op%maptyp(iplt)   == 'M'  .and.  &
!          op%contrtyp(iplt) /= 'T'  .and.  &
!          op%contrtyp(iplt) /= 'F'  .and.  &
!          op%contrtyp(iplt) /= 'O' )) then
!       call mkmap(iplt,'FILL')
!    endif

     ! Get units and stagpoint information for this field

     iw = 2
     k  = 2

     call oplot_lib(k,iw,'UNITS',trim(op%fldname(iplt)),1.,0.,fldval,notavail)

     ! If notavail = 3 is returned from above call, current plot field is unavailable

#ifdef OLAM_MPI
     if (iparallel == 1) then
        call MPI_Allgather(notavail, 1, MPI_INTEGER, notavails, 1, MPI_INTEGER, MPI_COMM_WORLD, iok)
        notavail = minval(notavails)
     endif
#endif

     if (notavail == 3) then
        if (myrank == 0) then
           write(io6,*) 'FIELD ',trim(op%fldname(iplt)),' NOT AVAILABLE IN THIS RUN.'
           call oplot_panel(op%panel(iplt),op%frameoff(iplt),op%colorbar(iplt),1.,op%projectn(iplt))
           call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)
           call o_plchhq(op%fx1,op%fnamey  &
                ,trim(op%fldname(iplt))//' NOT AVAILABLE IN THIS RUN'  &
                ,.012 * (op%hp2 - op%hp1),0.,-1.)
           call o_frame()
        endif
        cycle
     endif

     ! Plot current field values with tiles or filled contours

     if ( op%contrtyp(iplt) == 'T' .or. op%contrtyp(iplt) == 'F' .or. &
          op%contrtyp(iplt) == 'O' ) then
        op%ifill = 1
        call slab(iplt)
     endif

     ! Plot current field values with contour lines

     if ( op%contrtyp(iplt) == 'L' .or. op%contrtyp(iplt) == 'O' ) then
        op%ifill = 0
        call slab(iplt)
     endif

     ! Plot grid cell indices if so specified

     if ( op%pltindx(iplt) == 'I' .or. op%pltindx(iplt) == 'J' ) then
        call plot_index(iplt)
     endif

     ! Plot current field values with printed numbers if so specified

     if ( op%prtval(iplt) == 'P' ) then
        call slab_val(iplt)
     endif

     ! Plot T-point or V-point vectors or wind barbs if so specified

     call vectslab(iplt)

     ! SPECIAL HURRICANE FRANCES LOCATION PLOT
!    call plot_frances(iplt)

     ! Draw land/sea cell boundaries if so specified

     if ( op%pltgrid_landsea(iplt) == 'g' ) then
        call plot_grid_landsea(iplt)
     endif

     ! Draw dual grid cell boundaries if so specified

     if ( op%pltdualgrid(iplt) == 'D' .or. op%pltdualgrid(iplt) == 'd' ) then
        call plot_dualgrid(iplt)
     endif

     ! Draw grid cell boundaries if so specified

     if ( op%pltgrid(iplt) == 'G' ) then
        call plot_grid(iplt)
     endif

     ! Draw outer grid (frame) boundary if so specified

     if ( op%pltborder(iplt) == 'b' .and. myrank == 0 ) then
        call plot_grid_frame()
     endif

     ! Draw plot frame, tick marks, X and Y tick labels, X and Y axis labels
     ! if so specified.

     if ( op%pltborder(iplt) == 't' .and. myrank == 0 ) then

        if (op%projectn(iplt) == 'L') then
           call niceinc20(op%xmin,op%xmax,xinc,labincx)
           call niceinc20(op%ymin,op%ymax,yinc,labincy)

           if (op%windowin(iplt) /= 'W') then
              xinc = 10.
              yinc = 10.
              labincx = 3
              labincy = 3
           endif

           call oplot_xy2(op%panel(iplt),op%frameoff(iplt),op%colorbar(iplt),0.,.016,10,0, &
                          1,dum1,dum2,                                   &
                          'LONGITUDE (deg)','LATITUDE (deg)',            &
                          op%xmin,op%xmax,xinc,labincx,                  &
                          op%ymin,op%ymax,yinc,labincy                   )
        else
           call niceinc20(.001*op%xmin,.001*op%xmax,xinc,labincx)
           call niceinc20(.001*op%ymin,.001*op%ymax,yinc,labincy)

           if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') then
              ylabel = 'Z (km)'
           else
              ylabel = 'Y (km)'
           endif

           call oplot_xy2(op%panel(iplt),op%frameoff(iplt),op%colorbar(iplt),1.0,.016,10,0, &
                          1,dum1,dum2,                                    &
                          'X (km)',ylabel,                                &
                          .001*op%xmin,.001*op%xmax,xinc,labincx,         &
                          .001*op%ymin,.001*op%ymax,yinc,labincy          )
        endif

     endif

     ! Draw line map if so specified

     if ( (op%projectn(iplt) == 'L'  .or.   &
           op%projectn(iplt) == 'P'  .or.   & 
           op%projectn(iplt) == 'G'  .or.   & 
           op%projectn(iplt) == 'O') .and.  & 
          (op%maptyp(iplt)   /= 'N') .and.  &
          (op%maptyp(iplt)   == 'm'  .or.   &
           op%contrtyp(iplt) == 'T'  .or.   &
           op%contrtyp(iplt) == 'F'  .or.   &
           op%contrtyp(iplt) == 'O') ) then
        call mkmap(iplt,'LINE')
     endif

     ! Draw colorbar if so specified

     if ( op%colorbar(iplt) == 'c' .and. myrank == 0 ) then
        call plot_colorbar(op%icolortab(iplt))
     endif

     ! Draw label bar if so specified

     if ( op%labelbar(iplt) == 'n' .or. op%labelbar(iplt) == 'i' .and. &
          myrank == 0 ) then
        call date_add_to8 (iyear1,imonth1,idate1,itime1*100  &
             ,time_istp8,'s',outyear,outmonth,outdate,outhour)

        call plot_labelbar (iplt,op%label  &
             ,op%units,op%projectn(iplt),op%slabloc(iplt) &
             ,op%fldval_min,op%fldval_max  &
             ,time_istp8,outhour,outdate,outmonth,outyear)
     endif

     ! Plot ID if > 0

     if (id > 0 .and. myrank == 0) then
        write (title, '(i6)') id
        bsize = .010 * (op%hp2 - op%hp1)
        call o_plchhq(op%fx1,op%smaxy,'ID = '//trim(adjustl(title)),bsize,0.,-1.)
     endif

     ! Complete current frame if so specified

     if (op%frameoff(iplt) /= 'f' .and. op%frameoff(iplt) /= 'h' .and. myrank == 0) then
        call o_frame()
     endif

  enddo

  ! Close the current workstation if not a plotonly run and if output
  ! is to a NCAR graphics meta file. This allows viewing the complete
  ! meta file (including the last frame) during a run and in case the
  ! simulation crashes.

  if ( trim(runtype) /= 'PLOTONLY' .and. &
       op%plttype == 0 .and. myrank == 0 ) then
     call o_clswk()
  endif

end subroutine plot_fields

!===============================================================================

subroutine slab(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa
  use misc_coms,  only: io6
  use mem_para,   only: myrank

  implicit none

  integer, intent(in) :: iplt

  call oplot_set(iplt)

  ! Set plot line color (black)

  if (myrank == 0) then
     call o_gsplci(10)
     call o_gsfaci(10)
     call o_gstxci(10)
     call o_sflush()
  endif

  if ( op%projectn(iplt) == 'L' .or.  &
       op%projectn(iplt) == 'P' .or.  &
       op%projectn(iplt) == 'G' .or.  &
       op%projectn(iplt) == 'O' .or.  &
       op%projectn(iplt) == 'Z' ) then  ! For horizontal plotting...

     if ( op%contrtyp(iplt) /= 'F' .and. op%contrtyp(iplt) /= 'L' .and.  &
          op%contrtyp(iplt) /= 'O' ) then

        if (op%stagpt == 'T' .or. op%stagpt == 'W') then
           call tileslab_horiz_tw(iplt,'T')
        elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
           call tileslab_horiz_mp(iplt,'T')
        elseif (op%stagpt == 'V' .or. op%stagpt == 'N') then
           call tileslab_horiz_vn(iplt,'T')
        elseif (op%stagpt == 'S') then  ! for sea cells
           call tileslab_horiz_s(iplt,'T')
        elseif (op%stagpt == 'L') then  ! for land cells
           call tileslab_horiz_l(iplt,'T') 
        elseif (op%stagpt == 'B') then  ! for both sea and land cells
           op%stagpt = 'S'
           call tileslab_horiz_s(iplt,'T')
           op%stagpt = 'L'
           call tileslab_horiz_l(iplt,'T') 
           op%stagpt = 'B'
        endif

     else

        if     (op%stagpt == 'T' .or. op%stagpt == 'W') then
           call contslab_horiz_tw(iplt)
        elseif (op%stagpt == 'V' .or. op%stagpt == 'N') then
           call contslab_horiz_vn(iplt)
        elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
           call contslab_horiz_mp(iplt)
        elseif (op%stagpt == 'H') then
           call contslab_topmw(iplt)
        endif

     endif

  else  ! Vertical plots

     if ( op%contrtyp(iplt) /= 'F' .and. op%contrtyp(iplt) /= 'L' .and.  &
          op%contrtyp(iplt) /= 'O' ) then

        if     (op%stagpt == 'T' .or. op%stagpt == 'W') then
           call tileslab_vert_tw(iplt,'T')
        elseif (op%stagpt == 'V') then
           call tileslab_vert_v(iplt,'T')
        endif

     else

        if     (op%stagpt == 'T') then
           call contslab_vert_t(iplt)
        elseif (op%stagpt == 'W') then
           call contslab_vert_w(iplt)
        elseif (op%stagpt == 'V') then
           call contslab_vert_v(iplt)
        endif

     endif

  endif

end subroutine slab

!===============================================================================

subroutine slab_val(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa
  use misc_coms,  only: io6

  implicit none

  integer, intent(in) :: iplt

  call oplot_set(iplt)

  if ( op%projectn(iplt) == 'L' .or.  &
       op%projectn(iplt) == 'P' .or.  &
       op%projectn(iplt) == 'G' .or.  &
       op%projectn(iplt) == 'O' .or.  &
       op%projectn(iplt) == 'Z') then  ! For horizontal plotting...

     if     (op%stagpt == 'T' .or. op%stagpt == 'W') then
        call tileslab_horiz_tw(iplt,'P')
     elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
        call tileslab_horiz_mp(iplt,'P')
     elseif (op%stagpt == 'V' .or. op%stagpt == 'N') then
        call tileslab_horiz_vn(iplt,'P')
     elseif (op%stagpt == 'S') then  ! for sea cells
        call tileslab_horiz_s(iplt,'P')
     elseif (op%stagpt == 'L') then  ! for land cells
        call tileslab_horiz_l(iplt,'P') 
     elseif (op%stagpt == 'B') then  ! for both sea and land cells
        op%stagpt = 'S'
        call tileslab_horiz_s(iplt,'P')
        op%stagpt = 'L'
        call tileslab_horiz_l(iplt,'P') 
        op%stagpt = 'B'
     endif

  else  ! Vertical plots

     if     (op%stagpt == 'T' .or. op%stagpt == 'W') then
        call tileslab_vert_tw(iplt,'P')
     elseif (op%stagpt == 'V') then
        call tileslab_vert_v(iplt,'P')
     endif

  endif

end subroutine slab_val

!===============================================================================

subroutine plot_index(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mma, mva, mwa, xem, yem, zem, &
                        xev, yev, zev, xew, yew, zew
  use mem_ijtabs, only: itab_w, jtab_w, jtw_wadj, &
                        itab_v, jtab_v, jtv_wadj, &
                        itab_m, jtab_m, jtm_vadj
  use sea_coms,   only: mws, mus
  use mem_sea,    only: sea, itab_ws
  use leaf_coms,  only: mwl, mul
  use mem_leaf,   only: land, itab_wl
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: im,jm,img,iv,jv,ivg,iw,jw,iwg,iwl,iwlg,iws,iwsg
  real :: hpt,vpt

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j, i
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! This subroutine plots grid point indices for staggered points
  ! not covered by the current field

  call oplot_set(iplt)  ! not needed here?

  if ( op%projectn(iplt) /= 'L' .and.  &
       op%projectn(iplt) /= 'P' .and.  &
       op%projectn(iplt) /= 'G' .and.  &
       op%projectn(iplt) /= 'O' .and.  &
       op%projectn(iplt) /= 'Z' ) return

  ! Plot M point indices

  if (myrank == 0) then
     call o_sflush()
     call o_gsplci(10)
     call o_gstxci(10)
     call o_sflush()
  endif

  nu   = 0
  ipos = 0

  base = 2 * nbytes_real + nbytes_int
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mva) / 5. )
  else
     inc = mva
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  if ( op%pltindx(iplt) == 'J' .or.  &
       trim(op%stagpt)  == 'M' .or.  &
       trim(op%stagpt)  == 'P' ) then

!     do jm = 1, jtab_m(jtm_vadj)%jend(1)
!        im  = jtab_m(jtm_vadj)%im(jm)

        do im = 2,mma

        img = itab_m(im)%imglobe

        call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)
        if ( hpt < op%xmin .or. hpt > op%xmax .or.  &
             vpt < op%ymin .or. vpt > op%ymax ) cycle

        if (myrank == 0) then

           call oplot_locindex(img,hpt,vpt,op%psiz)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(img, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     enddo
  endif

  ! Plot U/V point indices

  if ( op%pltindx(iplt) == 'J' .or.  &
       trim(op%stagpt)  == 'V' .or.  &
       trim(op%stagpt)  == 'N' ) then

!     do jv = 1, jtab_v(jtv_wadj)%jend(1)
!        iv  = jtab_v(jtv_wadj)%iv(jv)
        do iv = 2,mva

        ivg = itab_v(iv)%ivglobe

        call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hpt,vpt)
        if ( hpt < op%xmin .or. hpt > op%xmax .or.  &
             vpt < op%ymin .or. vpt > op%ymax ) cycle

        if (myrank == 0) then

           call oplot_locindex(ivg,hpt,vpt,op%psiz)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(ivg, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     enddo
  endif

  ! Plot W point indices

  if ( op%pltindx(iplt) == 'J' .or.  &
       trim(op%stagpt)  == 'W' .or.  &
       trim(op%stagpt)  == 'T' ) then

!     do jw = 1, jtab_w(jtw_wadj)%jend(1)
!        iw  = jtab_w(jtw_wadj)%iw(jw)
        do iw = 2,mwa

        iwg = itab_w(iw)%iwglobe

        call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)
        if ( hpt < op%xmin .or. hpt > op%xmax .or.  &
             vpt < op%ymin .or. vpt > op%ymax ) cycle
        
        if (myrank == 0) then

           call oplot_locindex(iwg,hpt,vpt,op%psiz)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(iwg, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt, 1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

     enddo
  endif

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
               
                 call MPI_Unpack(buffer, buffsize, ipos, i,   1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, hpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 
                 call oplot_locindex(i,hpt,vpt,op%psiz)
               
              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  ! Plot WL point indices

  if ( trim(op%stagpt) == 'L' .or.  &
       trim(op%stagpt) == 'B' ) then

     if (myrank == 0) then
        call o_sflush()
        call o_gsplci(6)
        call o_gstxci(6)
        call o_sflush()
     endif

     nu   = 0
     ipos = 0

     base = 2 * nbytes_real + nbytes_int
     if (op%windowin(iplt) == 'W') then
        inc = ceiling( real(mwl) / 5. )
     else
        inc = mwl
     endif

     if (myrank > 0) then
        buffsize = inc * base
        allocate( buffer( buffsize ) )
     endif

     do iwl = 2, mwl
        iwlg = itab_wl(iwl)%iwglobe

        call oplot_transform(iplt,land%xew(iwl),land%yew(iwl),land%zew(iwl),hpt,vpt)

        if ( hpt < op%xmin .or. hpt > op%xmax .or.  &
             vpt < op%ymin .or. vpt > op%ymax ) cycle

        if (myrank == 0) then

           call oplot_locindex(iwlg,hpt,vpt,.7 * op%psiz)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(iwlg, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
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
               
                    call MPI_Unpack(buffer, buffsize, ipos, i,   1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, hpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, vpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 
                    call oplot_locindex(i,hpt,vpt,.7 * op%psiz)
               
                 enddo

              endif
           enddo
        endif

        deallocate(buffer)
     endif
#endif

     if (myrank == 0) call o_sflush()
  endif

  ! Plot WS point indices

  if ( trim(op%stagpt) == 'S' .or.  &
       trim(op%stagpt) == 'B' ) then

     if (myrank == 0) then
        call o_sflush()
        call o_gsplci(3)
        call o_gstxci(3)
        call o_sflush()
     endif

     nu   = 0
     ipos = 0

     base = 2 * nbytes_real + nbytes_int
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
        iwsg = itab_ws(iws)%iwglobe

        call oplot_transform(iplt, sea%xew(iws), sea%yew(iws), sea%zew(iws), hpt, vpt)

        if ( hpt < op%xmin .or. hpt > op%xmax .or.  &
             vpt < op%ymin .or. vpt > op%ymax ) cycle

        if (myrank == 0) then

           call oplot_locindex(iwsg,hpt,vpt,.7 * op%psiz) 

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(iwsg, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(hpt,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vpt,  1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
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
               
                    call MPI_Unpack(buffer, buffsize, ipos, i,   1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, hpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, vpt, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                 
                    call oplot_locindex(i,hpt,vpt,.7 * op%psiz)
               
                 enddo

              endif
           enddo
        endif

        deallocate(buffer)
     endif
#endif

     if (myrank == 0) call o_sflush()

  endif

end subroutine plot_index

!===============================================================================

subroutine vectslab(iplt)

  use oplot_coms, only: op
  use misc_coms,  only: io6
  use mem_para,   only: myrank

  implicit none

  integer, intent(in) :: iplt

  call oplot_set(iplt)

  if (myrank == 0) then
     call o_gsplci(10)
     call o_gstxci(10)
  endif

  if ( op%projectn(iplt) == 'L' .or.  &
       op%projectn(iplt) == 'P' .or.  &
       op%projectn(iplt) == 'G' .or.  &
       op%projectn(iplt) == 'O' .or.  &
       op%projectn(iplt) == 'Z' ) then

     if ( op%vectbarb(iplt) == 'V' ) call vectslab_horiz_v(iplt)

     if ( op%vectbarb(iplt) == 'w' .or. &
          op%vectbarb(iplt) == 'B' ) call vectslab_horiz_w(iplt)

  elseif (op%projectn(iplt) == 'C') then  ! etc for 'V'?

  endif

end subroutine vectslab

!===============================================================================

subroutine horizplot_k(iplt,mha,ktf,kc,wtbot,wttop)

  use mem_basic,  only: press
  use mem_grid,   only: mza, mma, mva, mwa, lpw, zm, zt
  use mem_ijtabs, only: itab_m, itab_v
  use oplot_coms, only: op
  use misc_coms,  only: io6

  implicit none

  integer, intent(in)  :: iplt,mha ! mha is horizontal index for given stagpt
  integer, intent(out) :: ktf(mwa),kc(mha)
  real,    intent(out) :: wtbot(mha),wttop(mha)

  integer :: iw,iv,ip,iw1,iw2,j,k,kp,kpf,npoly
  real :: plev,pressp1,pressp2,count2,pressv1,pressv2

  ! This subroutine determines the following quantities for horizontal plots:
  ! (1) ktf flag for T cells indicating 0-in plot range, 1-below ground, 
  !     2-above model top
  ! (2) vertical level kc and vertical weights wtbot, wttop for given stagpt

  if (op%pltlev(iplt) == 'p') plev = op%slabloc(iplt) * 100.  ! convert from mb to Pa

  do iw = 2, mwa

     if (op%dimens == '2') then

        ktf(iw) = 0

        if (op%stagpt == 'T' .or. op%stagpt == 'W') then
           kc(iw) = lpw(iw) ! not usually needed (except for 'PSFC', so that 's' is not needed)
           wtbot(iw) = 1.
           wttop(iw) = 0.
        endif

     elseif (op%pltlev(iplt) == 's') then

        ktf(iw) = 0

        if (op%stagpt == 'T' .or. op%stagpt == 'W') then
           kc(iw) = lpw(iw)
           wtbot(iw) = 1.
           wttop(iw) = 0.
        endif

     elseif (op%pltlev(iplt) == 'p') then

        k = lpw(iw) - 1

        do while (press(k+1,iw) > plev .and. k < mza)
           k = k + 1
        enddo

        if (k < lpw(iw)) then
           ktf(iw) = 1
        elseif (k >= mza) then
           ktf(iw) = 2
        else
           ktf(iw) = 0

           if (op%stagpt == 'T' .or. op%stagpt == 'W') then
              wtbot(iw) = (plev - press(k+1,iw)) / (press(k,iw) - press(k+1,iw))
              wttop(iw) = 1. - wtbot(iw)

              if (op%stagpt == 'W') then
                 wtbot(iw) = wtbot(iw) + .5
                 wttop(iw) = wttop(iw) - .5

                 if (wtbot(iw) > 1.) then
                    wtbot(iw) = wtbot(iw) - 1.
                    wttop(iw) = wttop(iw) + 1.
               
                    k = k - 1
                 endif

              endif

              kc(iw) = k

           endif

        endif

     else ! plot at specified height op%slabloc(iplt)

        k = 2

        do while (k < mza .and. zm(k) < op%slabloc(iplt))
           k = k + 1
        enddo

        if (k < lpw(iw)) then
           ktf(iw) = 1
        elseif (k > mza) then
           ktf(iw) = 2
        else
           ktf(iw) = 0
           if (op%stagpt == 'T' .or. op%stagpt == 'W') then
              wtbot(iw) = 1.
              wttop(iw) = 0.
              if (op%stagpt == 'W' .and. op%slabloc(iplt) < zt(k)) k = k - 1
              kc(iw) = k
           endif
        endif

     endif

  enddo

  ! 'V' or 'N' stagpt

  if (op%stagpt == 'V' .or. op%stagpt == 'N') then

     do iv = 1,mva

        if (op%dimens == '2') then

           kc(iv) = 2 ! not actually needed
           wtbot(iv) = 1.
           wttop(iv) = 0.

        elseif (op%pltlev(iplt) == 's') then

           iw1 = itab_v(iv)%iw(1)
           iw2 = itab_v(iv)%iw(2)

           kc(iv) = max(lpw(iw1),lpw(iw2))
           wtbot(iv) = 1.
           wttop(iv) = 0.

        elseif (op%pltlev(iplt) == 'p') then

           iw1 = itab_v(iv)%iw(1)
           iw2 = itab_v(iv)%iw(2)

           if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

              do k = max(lpw(iw1),lpw(iw2)),mza-1
                 kc(iv) = k
                 pressv1 = .5 * (press(k  ,iw1) + press(k  ,iw2))
                 pressv2 = .5 * (press(k+1,iw1) + press(k+1,iw2))
                 if (pressv2 < plev) exit
              enddo

              wtbot(iv) = (plev - pressv2) / (pressv1 - pressv2)
              wttop(iv) = 1. - wtbot(iv)

           elseif (ktf(iw1) == 0) then

              do k = lpw(iw1),mza-1
                 kc(iv) = k
                 pressv1 = press(k  ,iw1)
                 pressv2 = press(k+1,iw1)
                 if (pressv2 < plev) exit
              enddo

              wtbot(iv) = (plev - pressv2) / (pressv1 - pressv2)
              wttop(iv) = 1. - wtbot(iv)

           elseif (ktf(iw2) == 0) then

              do k = lpw(iw2),mza-1
                 kc(iv) = k
                 pressv1 = press(k  ,iw2)
                 pressv2 = press(k+1,iw2)
                 if (pressv2 < plev) exit
              enddo

              wtbot(iv) = (plev - pressv2) / (pressv1 - pressv2)
              wttop(iv) = 1. - wtbot(iv)

           else

              ! defaults, not used

              kc(iv) = 2
              wtbot(iv) = 1.
              wttop(iv) = 0.

           endif

        else ! plot at specified height op%slabloc(iplt)

           k = 2

           do while (k < mza .and. zm(k) < op%slabloc(iplt))
              k = k + 1
           enddo

           wtbot(iv) = 1.
           wttop(iv) = 0.

           if (op%stagpt == 'N' .and. op%slabloc(iplt) < zt(k)) k = k - 1

           kc(iv) = k

        endif

     enddo

  endif

  ! 'P' or 'M" stagpt

  if (op%stagpt == 'P' .or. op%stagpt == 'M') then

     do ip = 1,mma

        if (op%dimens == '2') then

           kc(ip) = 2 ! not actually needed
           wtbot(ip) = 1.
           wttop(ip) = 0.

        elseif (op%pltlev(iplt) == 's') then

           kc(ip) = 2
           npoly = itab_m(ip)%npoly
         
           do j = 1,npoly
              iw = itab_m(ip)%iw(j)
              if (kc(ip) < lpw(iw)) kc(ip) = lpw(iw)
           enddo

           wtbot(ip) = 1.
           wttop(ip) = 0.

        elseif (op%pltlev(iplt) == 'p') then

           kp = mza ! initial value
           kpf = 1  ! initial value

           npoly = itab_m(ip)%npoly
         
           do j = 1,npoly
              iw = itab_m(ip)%iw(j)
              if (ktf(iw) == 0) then
                 kpf = 0
                 if (kp > lpw(iw)) kp = lpw(iw)
              endif
           enddo
         
           if (kpf /= 0) then
         
              ! defaults, not used

              kc(ip) = 2
              wtbot(ip) = 1.
              wttop(ip) = 0.

           else

              ! Now, kp is minimum of surrounding lpw(iw) points, at least some
              ! of which are below constant-pressure surface

              do k = kp,mza-1

                 pressp2 = 0.
                 count2 = 0.

                 do j = 1,npoly
                    iw = itab_m(ip)%iw(j)
                    if (ktf(iw) == 0 .and. k >= lpw(iw)) then
                       count2 = count2 + 1.
                       pressp2 = pressp2 + press(k,iw)
                    endif
                 enddo
         
                 pressp2 = pressp2 / count2
         
                 if (k > kp .and. pressp2 < plev) then
                    wtbot(ip) = (plev - pressp2) / (pressp1 - pressp2)
                    wttop(ip) = 1. - wtbot(ip)
                    kc(ip) = k
                    exit
                 endif

                 pressp1 = pressp2  

              enddo

           endif
                            
        else ! plot at specified height op%slabloc(iplt)

           k = 2

           do while (k < mza .and. zm(k) < op%slabloc(iplt))
              k = k + 1
           enddo

           wtbot(ip) = 1.
           wttop(ip) = 0.

           if (op%stagpt == 'M' .and. op%slabloc(iplt) < zt(k)) k = k - 1

           kc(ip) = k

        endif

     enddo

  endif

end subroutine horizplot_k

!===============================================================================

subroutine plot_underground_w(iplt,ktzone)

  use oplot_coms, only: op
  use mem_grid,   only: mza, mwa, zm, zt, lpw, xem, yem, zem, xew, yew, zew
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt
  integer, intent(in) :: ktzone(mwa)

  real :: xpt,ypt
  integer :: k, n
  integer :: npoly,j,iw,jw,im,iok,iv1,iv2
  real :: htpn(7), vtpn(7)
  real :: topo1,topo2

  integer, allocatable :: buffer(:)
  integer :: nu, ier, buffsize, ipos, lpwiw
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40
  integer :: iflag180

  if (op%noundrg(iplt) == 'u') return

! This subroutine plots the W control volumes that are underground, whether
! they be triangles or hexagons.

  if (op%projectn(iplt) == 'L' .or.  &
      op%projectn(iplt) == 'P' .or.  &
      op%projectn(iplt) == 'G' .or.  &
      op%projectn(iplt) == 'O' .or.  &
      op%projectn(iplt) == 'Z') then  ! Horizontal cross section

     nu   = 0
     ipos = 0

     if (myrank > 0) then
        buffsize = jtab_w(jtw_prog)%jend(1) * (nbytes_int + 14 * nbytes_real)
        allocate( buffer( buffsize ) )
     endif

     do jw = 1,jtab_w(jtw_prog)%jend(1)
        iw = jtab_w(jtw_prog)%iw(jw)

        ! Skip this iw point if it is not underground

        if (ktzone(iw) /= 1) cycle

        ! Initialize iflag180
     
        iflag180 = 0

        ! Get tile plot coordinates

        if (op%projectn(iplt) == 'L') then   ! For wrap-around direction
           call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),xpt,ypt)
        endif

        npoly = itab_w(iw)%npoly

        do j = 1,npoly
           im = itab_w(iw)%im(j)

           call oplot_transform(iplt,xem(im),yem(im),zem(im),htpn(j),vtpn(j))

           ! Avoid wrap-around and set iflag180

           if (op%projectn(iplt) == 'L') then
              call ll_unwrap(xpt,htpn(j))
              if (htpn(j) < -180.001) iflag180 =  1
              if (htpn(j) >  180.001) iflag180 = -1
           endif

        enddo

        ! Jump out of loop if any cell corner is on other side of earth

        if (any(htpn(1:npoly) > 1.e11)) cycle

        ! Skip current iw cell if entire cell is outside plot window.

        if ( all(htpn(1:npoly) < op%xmin) .or. &
             all(htpn(1:npoly) > op%xmax) .or. &
             all(vtpn(1:npoly) < op%ymin) .or. &
             all(vtpn(1:npoly) > op%ymax) ) cycle

        if (myrank == 0) then

           ! Tile-plot cell with underground color
           call fillpolyg(npoly,htpn,vtpn,op%icigrnd)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           call MPI_Pack(npoly, 1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(htpn,  npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(vtpn,  npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
        ! at other end

        if (iflag180 /= 0) then

           do j = 1,npoly
              htpn(j) = htpn(j) + 360. * iflag180
           enddo

           if (myrank == 0) then

              call fillpolyg(npoly,htpn,vtpn,op%icigrnd)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              call MPI_Pack(npoly, 1,     MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(htpn,  npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(vtpn,  npoly, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
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

           buffsize = maxval(nus(2:mgroupsize)) * (nbytes_int + 14 * nbytes_real)
           allocate( buffer( buffsize ) )

           do n = 2, mgroupsize

              if (nus(n) > 0) then

                 call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

                 ipos = 0

                 do j = 1, nus(n)

                    call MPI_Unpack(buffer, buffsize, ipos, npoly, 1,     MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, htpn,  npoly, MPI_REAL,    MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, vtpn,  npoly, MPI_REAL,    MPI_COMM_WORLD, ier)

                    ! Tile-plot cell with underground color
                    call fillpolyg(npoly,htpn,vtpn,op%icigrnd)

                 enddo

              endif
           enddo
        endif

        deallocate(buffer)
     endif
#endif

  else  ! Vertical cross section

     nu   = 0
     ipos = 0

     if (myrank > 0) then
        buffsize = jtab_w(jtw_prog)%jend(1) * (nbytes_int + 6 * nbytes_real)
        allocate( buffer( buffsize ) )
     endif

     do jw = 1,jtab_w(jtw_prog)%jend(1)
        iw = jtab_w(jtw_prog)%iw(jw)

        ! Get horizontal plot coordinates for cells in this column

        if (op%projectn(iplt) == 'C') then
           call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
        elseif (op%projectn(iplt) == 'V') then
           call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn) ! need to fix for hex?
        endif

        ! Skip this cell if it does not intersect plot cone

        if (iok /= 1) cycle

        ! Skip if entire cell is outside plot window

        if ( (htpn(1) < op%xmin .and. htpn(2) < op%xmin)  .or. &
             (htpn(1) > op%xmax .and. htpn(2) > op%xmax) ) cycle

        if (myrank == 0) then

           do k = 2, lpw(iw)

              if (k < lpw(iw)) then 

                 ! Get T-cell vertical coordinates for fully underground cells
                 vtpn(1) = zm(k-1)
                 vtpn(2) = vtpn(1)
                 vtpn(3) = zm(k)
                 vtpn(4) = vtpn(3)

              else

                 ! Get T-cell vertical coordinates
                 vtpn(1) = zm(k-1)
                 vtpn(2) = vtpn(1)
                 vtpn(3) = topo2
                 vtpn(4) = topo1

              endif

              call fillpolyg(4,htpn,vtpn,op%icigrnd)

           enddo

        else

#ifdef OLAM_MPI
           nu = nu + 1
           call MPI_Pack(lpw(iw), 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(htpn,    4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(topo1,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(topo2,   1, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
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

           buffsize = maxval(nus(2:mgroupsize)) * (nbytes_int + 6 * nbytes_real)
           allocate( buffer( buffsize ) )

           do n = 2, mgroupsize

              if (nus(n) > 0) then

                 call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

                 ipos = 0

                 do j = 1, nus(n)
               
                    call MPI_Unpack(buffer, buffsize, ipos, lpwiw, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, htpn,  4, MPI_REAL,    MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, topo1, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
                    call MPI_Unpack(buffer, buffsize, ipos, topo2, 1, MPI_REAL,    MPI_COMM_WORLD, ier)
      
                    do k = 2, lpwiw

                       if (k < lpwiw) then 

                          ! Get T-cell vertical coordinates for fully underground cells
                          vtpn(1) = zm(k-1)
                          vtpn(2) = vtpn(1)
                          vtpn(3) = zm(k)
                          vtpn(4) = vtpn(3)

                       else

                          ! Get T-cell vertical coordinates
                          vtpn(1) = zm(k-1)
                          vtpn(2) = vtpn(1)
                          vtpn(3) = topo2
                          vtpn(4) = topo1

                       endif

                       call fillpolyg(4,htpn,vtpn,op%icigrnd)

                    enddo
                 enddo
              endif
           enddo
        endif

        deallocate(buffer)
     endif
#endif

  endif

end subroutine plot_underground_w

!===============================================================================

!!subroutine plot_underground_m(iplt,kt)
!!
!!use oplot_coms, only: op
!!use mem_grid,   only: mma, mza, zm, zt, lpw, xem, yem, zem, xew, yew, zew
!!use mem_ijtabs, only: itab_m, jtab_m, jtm_vadj
!!use misc_coms,  only: io6
!!use mem_basic,  only: press
!!
!!implicit none
!!
!!integer, intent(in) :: iplt,kt
!!
!!real :: xpt,ypt,plev
!!integer :: npoly,j,iw,k,im,jm,iok,iv1,iv2
!!real :: htpn(7),vtpn(7)
!!integer :: kw(7)
!!
!!kw(1:7) = kt
!!
!!if (op%projectn(iplt) == 'L' .or.  &
!!    op%projectn(iplt) == 'P' .or.  &
!!    op%projectn(iplt) == 'G' .or.  &
!!    op%projectn(iplt) == 'O' .or.  &
!!    op%projectn(iplt) == 'Z') then  ! Horizontal cross section
!!
!!   do jm = 1, jtab_m(jtm_vadj)%jend(1)
!!      im = jtab_m(jtm_vadj)%im(jm)
!!
!!!     Get tile plot coordinates
!!
!!      if (op%projectn(iplt) == 'L') then    ! For determining wrap-around direction
!!         call oplot_transform(iplt,xem(im),yem(im),zem(im),xpt,ypt)
!!      endif
!!      
!!      npoly = itab_m(im)%npoly
!!      do j = 1,npoly
!!
!!         iw = itab_m(im)%iw(j)
!!
!!! Reset kw values if plotting on pressure level
!!
!!         if (op%projectn(iplt) /= 'Z' .and.  &
!!             op%pltlev(iplt) == 'p') then
!!
!!            plev = op%slabloc(iplt) * 100.  ! convert from mb to Pa
!!
!!            kw(j) = lpw(iw) - 1
!!
!!            do while (press(kw(j)+1,iw) > plev .and. kw(j) < mza-1)
!!               kw(j) = kw(j) + 1
!!            enddo
!!         
!!         endif
!!
!!         call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(j),vtpn(j))
!!
!!         ! Avoid wrap-around
!!         if (op%projectn(iplt) == 'L') call ll_unwrap(xpt,htpn(j))
!!      enddo
!!
!!! Jump out of loop if any cell corner is on other side of earth
!!
!!      if (any(htpn(1:npoly) > 1.e11)) cycle
!!
!!!     Jump out of loop if entire cell is outside plot window.
!!
!!      if ( all(htpn(1:npoly) < op%xmin) .or. all(htpn(1:npoly) > op%xmax) .or.  &
!!           all(vtpn(1:npoly) < op%ymin) .or. all(vtpn(1:npoly) > op%ymax) ) cycle
!!      
!!!     If cell is underground, tile-plot it with underground color
!!
!!      if (any(kw(1:npoly) < lpw(itab_m(im)%iw(1:npoly)))) then
!!         call fillpolyg(npoly,htpn,vtpn,op%icigrnd)
!!      endif
!!
!!   enddo
!!
!!endif
!!return
!!end subroutine plot_underground_m

!===============================================================================

subroutine plot_grid(iplt)

  use oplot_coms, only: op
  use mem_grid,   only: mwa, mva, mza, xem, yem, zem, xew, yew, zew, zm
  use mem_ijtabs, only: itab_w, itab_v, jtab_v, jtab_w, jtv_wadj, jtw_prog
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: iflag180,iskip,jvmax
  integer :: j,k,im1,im2,iok,iv,jv

! integer :: iw,jw,iv1,iv2
! real :: topo1,topo2
! real :: vtpn(4)

  real :: htpn(4)
  real :: wta1,wta2,wta3,wtb1,wtb2,wtb3
  real :: xp1,xp2,yp1,yp2
  real :: xq1,xq2,yq1,yq2

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  if (myrank == 0) then
     call o_sflush()
     call o_gsplci(10)
     call o_gstxci(10)
  endif

  nu   = 0
  ipos = 0

  jvmax = jtab_v(jtv_wadj)%jend(1)

  base = (4 * nbytes_real) * mza * 2
  if ( op%windowin(iplt) == 'W' .or. &
       op%projectn(iplt) == 'C' .or. &
       op%projectn(iplt) == 'V' ) then
     inc = ceiling( real(jvmax) / 5. )
  else
     inc = jvmax
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  if ( op%projectn(iplt) == 'L'  .or.  &
       op%projectn(iplt) == 'P'  .or.  &
       op%projectn(iplt) == 'G'  .or.  &
       op%projectn(iplt) == 'O'  .or.  &
       op%projectn(iplt) == 'Z') then   ! Horizontal cross section

!     do jv = 1, jvmax
     do iv = 2, mva

        ! Get tile plot coordinates.

!        iv  = jtab_v(jtv_wadj)%iv(jv)
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        if (im1 < 2 .or. im2 < 2) cycle

        call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),xp1,yp1)
        call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),xp2,yp2)

        ! Avoid wrap-around and set iflag180

        iflag180 = 0

        if (op%projectn(iplt) == 'L') then
           call ll_unwrap(xp1,xp2)

           if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
           if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
        endif

        ! Truncate segment if it crosses plot window boundary or skip if
        ! both endpoints are outside plot window

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        ! Plot line segment
      
        if (myrank == 0) then

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + 4 * nbytes_real) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If this segment crosses +/- 180 degrees longitude in lat/lon plot,
        ! re-plot at other end

        if (iflag180 /= 0) then

           xp1 = xp1 + 360. * iflag180
           xp2 = xp2 + 360. * iflag180

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           ! Plot line segment

           if (myrank == 0) then
              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)
           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + 4 * nbytes_real) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
           endif
      
        endif

     enddo

  else  ! Vertical cross section

     ! First, draw vertical lines

!     do jv = 1, jvmax
         do iv = 1,mva

!        iv  = jtab_v(jtv_wadj)%iv(jv)
        
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        if (im1 < 2 .or. im2 < 2) cycle

        ! Subroutine coneplot_tri does not use iv index (which it calls iw).  It
        ! works on triangle concept but here is given limiting case of triangle
        ! with zero angle.

        call coneplot_tri(iplt,iv,xem(im1),yem(im1),zem(im1), &
             xem(im2),yem(im2),zem(im2),xem(im2),yem(im2),zem(im2), &
             wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,htpn)

!q      do jw = 1, jtab_w(jtw_prog)%jend(1)
!q         iw = jtab_w(jtw_prog)%iw(jw)
!q
!q         ! Get horizontal plot coordinates for cells in this column
!q
!q         if (op%projectn(iplt) == 'C') then
!q             call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
!q         elseif (op%projectn(iplt) == 'V') then
!q             call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn) ! need to fix for hex??
!q         endif

! Jump out of loop if this W point does not intersect plot cone or plane

        if (iok /= 1) cycle
 
        call trunc_segment(htpn(1),htpn(1),op%ymin,op%ymax,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 0) then

           if (myrank == 0) then

              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + 4 * nbytes_real) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif
           endif

        endif
         
        call trunc_segment(htpn(2),htpn(2),op%ymin,op%ymax,xq1,xq2,yq1,yq2,iskip)
      
        if (iskip == 0) then

           if (myrank == 0) then

              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + 4 * nbytes_real) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif

! Jump out of loop if either cell side is outside plot window. 

!t      if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
!t          htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle

!q      do k = 2,mza

! Get T-cell vertical coordinates

!q         vtpn(1:2) = zm(k-1)
!q         vtpn(3:4) = zm(k)

! Plot 4 sides of (k,iw) grid box...

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

!q         call trunc_segment(htpn(1),htpn(2),vtpn(1),vtpn(2),xq1,xq2,yq1,yq2,iskip)

!q         if (iskip == 0) then
!q            call o_frstpt (xq1,yq1)
!q            call o_vector (xq2,yq2)
!q         endif

!q         call trunc_segment(htpn(2),htpn(3),vtpn(2),vtpn(3),xq1,xq2,yq1,yq2,iskip)

!q         if (iskip == 0) then
!q            call o_frstpt (xq1,yq1)
!q            call o_vector (xq2,yq2)
!q         endif

!q         call trunc_segment(htpn(3),htpn(4),vtpn(3),vtpn(4),xq1,xq2,yq1,yq2,iskip)

!q         if (iskip == 0) then
!q            call o_frstpt (xq1,yq1)
!q            call o_vector (xq2,yq2)
!q         endif

!q         call trunc_segment(htpn(4),htpn(1),vtpn(4),vtpn(1),xq1,xq2,yq1,yq2,iskip)

!q         if (iskip == 0) then
!q            call o_frstpt (xq1,yq1)
!q            call o_vector (xq2,yq2)
!q         endif

!q      enddo

     enddo

     ! Second, draw horizontal lines (no parallel communication necessary)

     if (myrank == 0) then
        do k = 1,mza

           call trunc_segment(op%xmin,op%xmax,zm(k),zm(k),xq1,xq2,yq1,yq2,iskip)

           if (iskip == 0) then
              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)
           endif

        enddo
     endif

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * (4 * nbytes_real)
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0
              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, xq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, xq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call o_frstpt (xq1,yq1)
                 call o_vector (xq2,yq2)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

end subroutine plot_grid

!===============================================================================

subroutine plot_dualgrid(iplt)

  use oplot_coms,  only: op
  use mem_grid,    only: mwa, mva, mza, xem, yem, zem, xew, yew, zew, zm,  &
                         wnx, wny, wnz
  use mem_ijtabs,  only: itab_w, itab_v, jtab_v, jtv_wadj
  use consts_coms, only: erad
  use misc_coms,   only: io6, iparallel
  use mem_para,    only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: iflag180,iskip,jvmax
  integer :: iv,jv,iw1,iw2

  real :: xp1,xp2,yp1,yp2
  real :: xq1,xq2,yq1,yq2

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Do not plot anything for vertical cross sections

  if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') return

  if (myrank == 0) then
     call o_sflush()
     call o_gsplci(6)
     call o_gstxci(6)
     call o_sflush()
  endif

  nu   = 0
  ipos = 0

  jvmax = jtab_v(jtv_wadj)%jend(1)

  base = (4 * nbytes_real) * 2
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

     ! Get tile plot coordinates.

     iv = jtab_v(jtv_wadj)%iv(jv)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     if (iw1 < 2 .or. iw2 < 2) cycle

     if (op%pltdualgrid(iplt) == 'D') then
        call oplot_transform(iplt,xew(iw1),yew(iw1),zew(iw1),xp1,yp1)
        call oplot_transform(iplt,xew(iw2),yew(iw2),zew(iw2),xp2,yp2)
     endif
   
     ! Avoid wrap-around and set iflag180

     iflag180 = 0

     if (op%projectn(iplt) == 'L') then
        call ll_unwrap(xp1,xp2)
        if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
        if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
     endif

     ! Truncate segment if it crosses plot window boundary or skip if
     ! both endpoints are outside plot window

     call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

     if (iskip == 1) cycle

     ! Plot line segment

     if (myrank == 0) then

        call o_frstpt (xq1,yq1)
        call o_vector (xq2,yq2)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        if (buffsize < ipos + 4 * nbytes_real) then
           allocate( bcopy (buffsize + inc * base) )
           bcopy(1:buffsize) = buffer
           call move_alloc(bcopy, buffer)
           buffsize = size(buffer)
        endif

        call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

     ! If this segment crosses +/- 180 degrees longitude in lat/lon plot,
     ! re-plot at other end

     if (iflag180 /= 0) then

        xp1 = xp1 + 360. * iflag180
        xp2 = xp2 + 360. * iflag180

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        ! Plot line segment

        if (myrank == 0) then

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + 4 * nbytes_real) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
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

        buffsize = maxval(nus(2:mgroupsize)) * (4 * nbytes_real)
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0
              do j = 1, nus(n)

                 call MPI_Unpack(buffer, buffsize, ipos, xq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, xq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call o_frstpt (xq1,yq1)
                 call o_vector (xq2,yq2)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  if (myrank == 0) call o_sflush()

end subroutine plot_dualgrid

!===============================================================================

subroutine plot_grid_landsea(iplt)

  use oplot_coms, only: op
  use sea_coms,   only: mws, iseagrid
  use mem_sea,    only: sea, itab_ws
  use leaf_coms,  only: mwl
  use mem_leaf,   only: land, itab_wl
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real
  use max_dims,   only: maxnlspoly

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: iflag180,iskip

  integer :: iws,iwl,jm1,jm2,npoly

  real :: xp1,xp2,yp1,yp2
  real :: xq1,xq2,yq1,yq2

  integer, allocatable :: buffer(:), bcopy(:)
  integer :: nu, ier, buffsize, ipos, base, inc, n, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  ! Do not plot anything for vertical cross sections

  if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') then
     write(io6,*) 'Land/sea grid plot not available in vertical cross section.'
     return
  endif

!-------------------------------------------------
! Plot sea grid
!-------------------------------------------------

  if (myrank == 0) then
     call o_sflush()
     call o_gsplci(3)
     call o_gstxci(3)
  endif

  nu   = 0
  ipos = 0

  base = 4 * nbytes_real
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

     npoly = itab_ws(iws)%npoly

     do jm1 = 1,npoly
        jm2 = jm1 + 1
        if (jm2 > npoly) jm2 = 1

        ! Get tile plot coordinates.

        call oplot_transform(iplt, itab_ws(iws)%xem(jm1), &
                                   itab_ws(iws)%yem(jm1), &
                                   itab_ws(iws)%zem(jm1), xp1, yp1)
        call oplot_transform(iplt, itab_ws(iws)%xem(jm2), &
                                   itab_ws(iws)%yem(jm2), &
                                   itab_ws(iws)%zem(jm2), xp2, yp2)

        ! Avoid wrap-around and set iflag180

        iflag180 = 0

        if (op%projectn(iplt) == 'L') then
           call ll_unwrap(xp1,xp2)

           if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
           if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
        endif

        ! Truncate segment if it crosses plot window boundary or skip
        ! if both endpoints are outside plot window

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        ! Plot line segment

        if (myrank == 0) then

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If this segment crosses +/- 180 degrees longitude in lat/lon plot, re-plot
        ! re-plot at other end

        if (iflag180 /= 0) then

           xp1 = xp1 + 360. * iflag180
           xp2 = xp2 + 360. * iflag180

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           ! Plot line segment

           if (myrank == 0) then

              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif

     enddo  ! jm1

  enddo  ! iws

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

                 call MPI_Unpack(buffer, buffsize, ipos, xq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, xq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call o_frstpt (xq1,yq1)
                 call o_vector (xq2,yq2)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

!-------------------------------------------------
! Plot land grid
!-------------------------------------------------

  if (myrank == 0) then
     call o_sflush()
     call o_gsplci(6)
     call o_gstxci(6)
  endif
     
  nu   = 0
  ipos = 0

  base = 4 * nbytes_real
  if (op%windowin(iplt) == 'W') then
     inc = ceiling( real(mws) / 5. )
  else
     inc = mws
  endif

  if (myrank > 0) then
     buffsize = inc * base
     allocate( buffer( buffsize ) )
  endif

  do iwl = 2, mwl

     npoly = itab_wl(iwl)%npoly

     do jm1 = 1,npoly
        jm2 = jm1 + 1
        if (jm2 > npoly) jm2 = 1

        ! Get tile plot coordinates.

        call oplot_transform(iplt, itab_wl(iwl)%xem(jm1), &
                                   itab_wl(iwl)%yem(jm1), &
                                   itab_wl(iwl)%zem(jm1), xp1, yp1)
        call oplot_transform(iplt, itab_wl(iwl)%xem(jm2), &
                                   itab_wl(iwl)%yem(jm2), &
                                   itab_wl(iwl)%zem(jm2), xp2, yp2)

        ! Avoid wrap-around and set iflag180

        iflag180 = 0

        if (op%projectn(iplt) == 'L') then
           call ll_unwrap(xp1,xp2)

           if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
           if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
        endif

        ! Truncate segment if it crosses plot window boundary or skip
        ! if both endpoints are outside plot window

        call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        ! Plot line segment

        if (myrank == 0) then

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)

        else

#ifdef OLAM_MPI
           nu = nu + 1
           if (buffsize < ipos + base) then
              allocate( bcopy (buffsize + inc * base) )
              bcopy(1:buffsize) = buffer
              call move_alloc(bcopy, buffer)
              buffsize = size(buffer)
           endif

           call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
           call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

        endif

        ! If this segment crosses +/- 180 degrees longitude in lat/lon plot,
        ! re-plot at other end

        if (iflag180 /= 0) then

           xp1 = xp1 + 360. * iflag180
           xp2 = xp2 + 360. * iflag180

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           ! Plot line segment

           if (myrank == 0) then

              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)

           else

#ifdef OLAM_MPI
              nu = nu + 1
              if (buffsize < ipos + base) then
                 allocate( bcopy (buffsize + inc * base) )
                 bcopy(1:buffsize) = buffer
                 call move_alloc(bcopy, buffer)
                 buffsize = size(buffer)
              endif

              call MPI_Pack(xq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(xq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq1, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
              call MPI_Pack(yq2, 1, MPI_REAL, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

           endif

        endif

     enddo  ! jm1

  enddo  ! iwl

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

                 call MPI_Unpack(buffer, buffsize, ipos, xq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, xq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq1, 1, MPI_REAL, MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, yq2, 1, MPI_REAL, MPI_COMM_WORLD, ier)

                 call o_frstpt (xq1,yq1)
                 call o_vector (xq2,yq2)

              enddo

           endif
        enddo
     endif

     deallocate(buffer)
  endif
#endif

  if (myrank == 0) call o_sflush()

end subroutine plot_grid_landsea

!===============================================================================

subroutine plot_grid_frame()

  use oplot_coms, only: op
  use misc_coms,  only: io6
  use mem_para,   only: myrank

  implicit none

  if (myrank /= 0) return

  call o_sflush()
  call o_gsplci(10)
  call o_gstxci(10)

  call o_frstpt(op%xmin,op%ymin)
  call o_vector(op%xmax,op%ymin)
  call o_vector(op%xmax,op%ymax)
  call o_vector(op%xmin,op%ymax)
  call o_vector(op%xmin,op%ymin)
  call o_sflush()

end subroutine plot_grid_frame

!===============================================================================

subroutine mkmap(iplt,filltype)

  use oplot_coms,  only: op
  use consts_coms, only: erad, erad2
  use misc_coms,   only: io6
  use mem_para,    only: myrank

  implicit none

  integer, intent(in) :: iplt
  character(len=*), intent(in) :: filltype

  call oplot_set(iplt)

  if (myrank /= 0) return

  ! Set plot color (black)

  call o_gsplci(10)
  call o_gsfaci(10)
  call o_gstxci(10)
  call o_sflush()

  if (trim(filltype) == 'LINE') then

     call o_mapint()
     call o_mappos(op%h1,op%h2,op%v1,op%v2)
     call o_mapsti('GR',10)    ! For 10-degree spacing of plotted parallels and meridians
     call o_mapsti('DA',65535) ! To plot parallels and meridians with solid lines
!    call o_mapstc('OU','PS')  ! To plot bnds of continents, countries, and US  states
!    call o_mapstc('OU','NO')  ! To plot NO geographic information
     call o_mapstc('OU','CO')  ! To plot continental outlines
!    call o_mapstc('OU','US')  ! To plot US state outlines
!    call o_mapstc('OU','PO')  ! To plot continental outlines + international outlines

     if (op%projectn(iplt) == 'L') then
        call o_maproj('CE',0.,op%plon3,0.)  ! Force plat to be zero for CE projection
        call o_mapset('LI',op%xmin,op%xmax,op%ymin,op%ymax)
     elseif (op%projectn(iplt) == 'P') then
        call o_maproj('ST',op%plat3,op%plon3,0.)
        call o_mapset('LI',op%xmin/erad2,op%xmax/erad2  &
                          ,op%ymin/erad2,op%ymax/erad2)
     elseif (op%projectn(iplt) == 'O') then
        call o_maproj('OR',op%plat3,op%plon3,0.)
        call o_mapset('LI',op%xmin/erad,op%xmax/erad  &
                          ,op%ymin/erad,op%ymax/erad)
     endif

     call o_mapint()             ! Initialize the above parameters
!    call mapgrd()
     call o_maplot()
     call o_sflush()

  endif

end subroutine mkmap

!===============================================================================

subroutine oplot_transform(iplt,xeq,yeq,zeq,xout,yout)

  use oplot_coms, only: op
  use misc_coms,  only: io6

  implicit none

  integer, intent(in)  :: iplt
  real,    intent(in)  :: xeq,yeq,zeq
  real,    intent(out) :: xout,yout

  ! Transform to plot frame coordinates

  if (op%projectn(iplt) == 'L') then

     call e_ec(xeq,yeq,zeq,xout,yout)
     ! xout is now true geographic longitude; shift according to specified window center
     xout = xout - op%plon3
     if (xout > 180.) then
        xout = xout - 360.
     elseif (xout < -180.) then
        xout = xout + 360.
     endif

  elseif (op%projectn(iplt) == 'P') then
     call e_ps(xeq,yeq,zeq,op%plat3,op%plon3,xout,yout)
  elseif (op%projectn(iplt) == 'G') then
     call e_gn(xeq,yeq,zeq,op%plat3,op%plon3,xout,yout)
  elseif (op%projectn(iplt) == 'O') then
     call e_or(xeq,yeq,zeq,op%plat3,op%plon3,xout,yout)
  else  ! Cartesian domain
     xout = xeq - op%plon3
     yout = yeq - op%plat3
  endif

end subroutine oplot_transform

!===============================================================================

subroutine oplot_panel(panel,frameoff,colorbar0,aspect,projectn)

  use oplot_coms, only: op, offlabs9
  use misc_coms,  only: io6

  implicit none

  character(len=1), intent(in) :: panel,frameoff,colorbar0,projectn
  real,             intent(in) :: aspect

  real :: pscale2x2, pscale3x3, poffset2x2, poffset3x3, asp

 ! offlabs9 = .true. ! Uncomment this to remove some axis labels in 3x3 panel plots
                     ! and to slightly expand size of each panel

  ! Relative positions within panel (if colorbar is used)

  op%fx1 = .12  ! plot frame left side x coord
  op%fx2 = .86  ! plot frame right side x coord
  op%fy1 = .20  ! plot frame bottom y coord

  op%fy2 = .94  ! plot frame top y coord

  op%cbx1 = .88 ! colorbar left side x coord
  op%cbx2 = .91 ! colorbar right side x coord
  op%cblx = .92 ! colorbar label left end x coord

  op%xlaby = .14   ! x-axis label y coord
  op%xtlaby = .175 ! x-axis tick label y coord

  op%ylabx = .03  ! y-axis label x coord
  op%ytlabx = .11 ! y-axis tick label right end x coord

  op%timex = .78 ! elapsed time (sec) left end x coord

  ! Information block coordinates

  op%slabx = .78 ! slab left end x coord

  op%timsy = .11 ! elapsed time (sec) y coord
  op%timdy = .09 ! elapsed time (day) y coord
  op%slaby = .07 ! slab y coord
  op%sminy = .05 ! slab min y coord
  op%smaxy = .03 ! slab max y coord

  if (offlabs9) then
     if (frameoff == 'h' .or. &
            panel == '5' .or. &
            panel == '6' .or. &
            panel == '7' .or. &
            panel == '8' .or. &
            panel == '9') then

        ! Expanded panels in 3x3 configuration (with some axis labels suppressed)

        op%fx1    = .09 ! plot frame left side x coord
        op%fx2    = .93 ! plot frame right side x coord
        op%ytlabx = .08 ! y-axis tick label right end x coord
        op%cbx1   = .95 ! colorbar left side x coord
        op%cbx2   = .98 ! colorbar right side x coord
        op%cblx   = .99 ! colorbar label left end x coord

        op%fy1    = .10  ! plot frame bottom y coord
        op%xtlaby = .075 ! x-axis tick label y coord
        op%xlaby  = .04  ! x-axis label y coord
     endif
  endif

  ! Expand plot area if no colorbar

  if (colorbar0 /= 'c') then

     op%fx1 = .14  ! plot frame left side x coord
     op%fx2 = .97  ! plot frame right side x coord
     op%fy1 = .11  ! plot frame bottom y coord
     op%fy2 = .94  ! plot frame top y coord

!    ! special for hurricane
!    if (projectn == 'C' .or. projectn == 'V') then
!       op%fy2 = .22  ! plot frame top y coord
!    endif

     op%xlaby = .04 ! x-axis label y coord
     op%xtlaby = .08 ! x-axis tick label y coord

     op%ylabx = .03 ! y-axis label x coord
     op%ytlabx = .13 ! y-axis tick label right end x coord

     op%timex = .03  ! elapsed time (sec) left end x coord
     op%timsy = .035 ! elapsed time (sec) y coord
     op%timdy = .015 ! elapsed time (day) y coord

     op%slabx = .78  ! slab left end x coord
     op%slaby = .055 ! slab y coord
     op%sminy = .035 ! slab min y coord
     op%smaxy = .015 ! slab max y coord

  endif

  ! Reduce plot height under some circumstances (make plot window rectangular
  ! instead of square)

  if (aspect < 1.e-3) then

     ! Aspect set to 0 scales plot window to physical coordinate limits

     asp = (op%ymax - op%ymin) / (op%xmax - op%xmin)
     op%fy2 = op%fy1 + (op%fx2 - op%fx1) * asp

  elseif (aspect < .999) then

     ! Aspect between 0 and 1 is used as specified plot window height/width ratio

     op%fy2 = op%fy1 + (op%fx2 - op%fx1) * aspect
   
  endif

  op%fnamey = op%fy2 + .025 ! field name y coord

  ! Panel borders

  op%hp1 = 0.0
  op%hp2 = 1.0
  op%vp1 = 0.0
  op%vp2 = 1.0

  pscale2x2  = 0.5
  poffset2x2 = 0.5

  pscale3x3  = 0.33333
  poffset3x3 = 0.325

!!!!!!!!!!!!!!!!!!!!!!!!!! special - modify plot coordinates

! op%fx2 = .97  ! plot frame right side x coord
! op%fy1 = .14  ! plot frame bottom y coord

! op%xlaby = .06 ! x-axis label y coord
! op%xtlaby = .10 ! x-axis tick label y coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if     (panel == '1' .and. frameoff /= 'h') then
     op%hp1 = pscale2x2 * op%hp1
     op%hp2 = pscale2x2 * op%hp2
     op%vp1 = pscale2x2 * op%vp1 + poffset2x2 * aspect
     op%vp2 = pscale2x2 * op%vp2 + poffset2x2 * aspect
  elseif (panel == '2' .and. frameoff /= 'h') then
     op%hp1 = pscale2x2 * op%hp1 + poffset2x2
     op%hp2 = pscale2x2 * op%hp2 + poffset2x2
     op%vp1 = pscale2x2 * op%vp1 + poffset2x2 * aspect
     op%vp2 = pscale2x2 * op%vp2 + poffset2x2 * aspect
  elseif (panel == '3' .and. frameoff /= 'h') then
     op%hp1 = pscale2x2 * op%hp1
     op%hp2 = pscale2x2 * op%hp2
     op%vp1 = pscale2x2 * op%vp1
     op%vp2 = pscale2x2 * op%vp2
  elseif (panel == '4' .and. frameoff /= 'h') then
     op%hp1 = pscale2x2 * op%hp1 + poffset2x2
     op%hp2 = pscale2x2 * op%hp2 + poffset2x2
     op%vp1 = pscale2x2 * op%vp1
     op%vp2 = pscale2x2 * op%vp2
  elseif (panel == '1') then
     op%hp1 = pscale3x3 * op%hp1
     op%hp2 = pscale3x3 * op%hp2
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * 2.0 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * 2.0 * aspect
  elseif (panel == '2') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * 2.0 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * 2.0 * aspect
  elseif (panel == '3') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3 * 2.0
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3 * 2.0
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * 2.0 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * 2.0 * aspect
  elseif (panel == '4') then
     op%hp1 = pscale3x3 * op%hp1
     op%hp2 = pscale3x3 * op%hp2
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * aspect
  elseif (panel == '5') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * aspect
  elseif (panel == '6') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3 * 2.0
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3 * 2.0
     op%vp1 = pscale3x3 * op%vp1 + poffset3x3 * aspect
     op%vp2 = pscale3x3 * op%vp2 + poffset3x3 * aspect
  elseif (panel == '7') then
     op%hp1 = pscale3x3 * op%hp1
     op%hp2 = pscale3x3 * op%hp2
     op%vp1 = pscale3x3 * op%vp1
     op%vp2 = pscale3x3 * op%vp2
  elseif (panel == '8') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3
     op%vp1 = pscale3x3 * op%vp1
     op%vp2 = pscale3x3 * op%vp2
  elseif (panel == '9') then
     op%hp1 = pscale3x3 * op%hp1 + poffset3x3 * 2.0
     op%hp2 = pscale3x3 * op%hp2 + poffset3x3 * 2.0
     op%vp1 = pscale3x3 * op%vp1
     op%vp2 = pscale3x3 * op%vp2
  endif

  ! Frame borders

  op%h1 = op%hp1 + op%fx1 * (op%hp2 - op%hp1)
  op%h2 = op%hp1 + op%fx2 * (op%hp2 - op%hp1)
  op%v1 = op%vp1 + op%fy1 * (op%vp2 - op%vp1)
  op%v2 = op%vp1 + op%fy2 * (op%vp2 - op%vp1)

end subroutine oplot_panel

!===============================================================================

subroutine oplot_set(iplt)

  use misc_coms,   only: io6, mdomain, deltax, iparallel
  use oplot_coms,  only: op
  use oname_coms,  only: nl
  use mem_grid,    only: xem, yem, xew, yew, zew, zm, mma, mwa, arw0, nza, mza
  use consts_coms, only: erad, pio180
  use mem_para,    only: myrank, mgroupsize

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iplt

  integer :: iw,iok,im,iv1,iv2
  real :: dist,radcone,xmid,ymid
  real :: delxmin,arwmin,xpt,ypt
  real :: aspect
  real :: topo1,topo2
  real :: htpn(8)
  real, allocatable :: buffer(:,:)

  ! Set plot window center or cone center coordinates (whichever is relevant),
  ! plus cone angle

  op%plon3    = nl%plotspecs(iplt)%plotcoord1
  op%plat3    = nl%plotspecs(iplt)%plotcoord2
  op%coneang  = nl%plotspecs(iplt)%slabloc
  op%viewazim = nl%plotspecs(iplt)%viewazim

  ! Set limits of plot window in model domain coordinates

  if (op%windowin(iplt) == 'W') then

     ! If windowing in, use specified window size

     op%xmin = -.5 * nl%plotspecs(iplt)%plotwid
     op%xmax =  .5 * nl%plotspecs(iplt)%plotwid
     op%ymin = -.5 * nl%plotspecs(iplt)%plotwid
     op%ymax =  .5 * nl%plotspecs(iplt)%plotwid

     ! Special treatment for EC (lat/lon) plot: use plot center lat as window
     ! center (Plot center longitude is already accounted for in transformation)

     if (op%projectn(iplt) == 'L') then
        op%xmin = op%plon3 - .5 * nl%plotspecs(iplt)%plotwid
        op%xmax = op%plon3 + .5 * nl%plotspecs(iplt)%plotwid
        op%ymin = op%plat3 - .5 * nl%plotspecs(iplt)%plotwid
        op%ymax = op%plat3 + .5 * nl%plotspecs(iplt)%plotwid
     endif

     ! Special treatment for vertical plot: ymin/ymax related to model Z coordinate

     if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') then

        if (abs(nl%zplot_min + 1.0) < 1.e-3) then
           op%ymin = zm(1)
        else
           op%ymin = max(nl%zplot_min, zm(1))
        endif

        if (abs(nl%zplot_max + 1.0) < 1.e-3) then
           op%ymax = zm(nza)
        else
           op%ymax = min(nl%zplot_max, zm(nza))
        endif

        ! special for hurricane
!       op%xmin = 0.
!       op%ymax = 6.e3
      
        ! special for dudhia expts
!       op%ymax = 10.e3

        if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax
     endif

  else

     ! If not windowing in, use following defaults

     if (op%projectn(iplt) == 'L') then

        op%xmin = -180.
        op%xmax =  180.
        op%ymin =  -90.
        op%ymax =   90.

     elseif (op%projectn(iplt) == 'P') then

        op%xmin = - 5. * erad
        op%xmax =   5. * erad
        op%ymin = - 5. * erad
        op%ymax =   5. * erad

     elseif (op%projectn(iplt) == 'G') then

        op%xmin = - 5. * erad
        op%xmax =   5. * erad
        op%ymin = - 5. * erad
        op%ymax =   5. * erad

     elseif (op%projectn(iplt) == 'O') then

        op%xmin = - 1.1 * erad
        op%xmax =   1.1 * erad
        op%ymin = - 1.1 * erad
        op%ymax =   1.1 * erad

     elseif (op%projectn(iplt) == 'V' .or. &
            (op%projectn(iplt) == 'C' .and. mdomain >= 2)) then

        op%xmin =  1.e9
        op%xmax = -1.e9
        do im = 2, mma
           op%xmin = min(op%xmin,xem(im),yem(im))
           op%xmax = max(op%xmax,xem(im),yem(im))
        enddo
        op%xmin = 1.01 * op%xmin
        op%xmax = 1.01 * op%xmax

        if (abs(nl%zplot_min + 1.0) < 1.e-3) then
           op%ymin = zm(1)
        else
           op%ymin = max(nl%zplot_min, zm(1))
        endif

        if (abs(nl%zplot_max + 1.0) < 1.e-3) then
           op%ymax = zm(nza)
        else
           op%ymax = min(nl%zplot_max, zm(nza))
        endif

        if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

     elseif (op%projectn(iplt) == 'C') then

        op%xmin = -3.4 * erad * sin(op%coneang * pio180)
        op%xmax =  3.4 * erad * sin(op%coneang * pio180)

        if (abs(nl%zplot_min + 1.0) < 1.e-3) then
           op%ymin = zm(1)
        else
           op%ymin = max(nl%zplot_min, zm(1))
        endif

        if (abs(nl%zplot_max + 1.0) < 1.e-3) then
           op%ymax = zm(nza)
        else
           op%ymax = min(nl%zplot_max, zm(nza))
        endif

        if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

     elseif (op%projectn(iplt) == 'Z') then

        op%xmin = minval(xem(2:mma))
        op%xmax = maxval(xem(2:mma))
        
        op%ymin = minval(yem(2:mma))
        op%ymax = maxval(yem(2:mma))

#ifdef OLAM_MPI
        if (iparallel == 1) then

           allocate(buffer(4,mgroupsize))
           buffer(1,myrank) = op%xmin
           buffer(2,myrank) = op%xmax
           buffer(3,myrank) = op%ymin
           buffer(4,myrank) = op%ymax

           call MPI_Allgather(MPI_IN_PLACE, 4, MPI_REAL, buffer, 4, MPI_REAL, MPI_COMM_WORLD, iok)

           op%xmin = minval(buffer(1,:))
           op%xmax = maxval(buffer(2,:))
           op%ymin = minval(buffer(3,:))
           op%ymax = maxval(buffer(4,:))
           deallocate(buffer)

        endif
#endif

        dist = max(op%xmax - op%xmin,op%ymax - op%ymin)
        xmid = .5 * (op%xmin + op%xmax)
        ymid = .5 * (op%ymin + op%ymax)

        op%xmin = xmid - .52 * dist
        op%xmax = xmid + .52 * dist
        op%ymin = ymid - .52 * dist
        op%ymax = ymid + .52 * dist

     endif

  endif

  if ( op%projectn(iplt) == 'L' .or.  &
       op%projectn(iplt) == 'P' .or.  &
       op%projectn(iplt) == 'G' .or.  &
       op%projectn(iplt) == 'O' .or.  &
       op%projectn(iplt) == 'Z' ) then   ! Horizontal cross section

     ! Find deltax of finest grid cell IN CURRENT PLOT WINDOW

     arwmin = 1.e13
     do iw = 2,mwa
        call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),xpt,ypt)

        if ( xpt < op%xmin .or. xpt > op%xmax .or.  &
             ypt < op%ymin .or. ypt > op%ymax ) cycle

        arwmin = min(arwmin,arw0(iw))
     enddo

#ifdef OLAM_MPI
     if (iparallel == 1) then
        allocate(buffer(1,mgroupsize))
        call MPI_Allgather(arwmin, 1, MPI_REAL, buffer, 1, MPI_REAL, MPI_COMM_WORLD, iok)
        arwmin = minval(buffer)
        deallocate(buffer)
     endif
#endif

     delxmin = sqrt(arwmin)

     if (op%projectn(iplt) == 'L') then

        if     (op%prtval_size == 'small' .or. op%prtval_size == 'SMALL') then
           op%psiz = 0.04 * (delxmin / (6. * erad)) * (400. / (op%xmax - op%xmin))
        elseif (op%prtval_size == 'large' .or. op%prtval_size == 'LARGE') then
           op%psiz = 0.12 * (delxmin / (6. * erad)) * (400. / (op%xmax - op%xmin))
        else
           op%psiz = 0.08 * (delxmin / (6. * erad)) * (400. / (op%xmax - op%xmin))
        endif

        op%vsprd = .3 * (delxmin / (6. * erad)) * 400.

     else

        if     (op%prtval_size == 'small' .or. op%prtval_size == 'SMALL') then
           op%psiz = 0.04 * delxmin / (op%xmax - op%xmin)
        elseif (op%prtval_size == 'large' .or. op%prtval_size == 'LARGE') then
           op%psiz = 0.12 * delxmin / (op%xmax - op%xmin)
        else
           op%psiz = 0.08 * delxmin / (op%xmax - op%xmin)
        endif

        op%vsprd = .3 * delxmin

     endif

  else                           ! Vertical cross section

     ! Get horizontal plot coordinates for cells in this column

     arwmin = 1.e13
     do iw = 2, mwa

        if (op%projectn(iplt) == 'C') then
           call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
        elseif (op%projectn(iplt) == 'V') then
           call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn) ! Need to fix for hex?
        endif

        if (iok /= 1) cycle ! If this W pt does not intersect plot cone

        ! Jump out of loop if either cell side is outside plot window. 

        if ( htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
             htpn(2) < op%xmin .or. htpn(2) > op%xmax ) cycle

        arwmin = min(arwmin,arw0(iw))
     enddo

#ifdef OLAM_MPI
     if (iparallel == 1) then
        allocate(buffer(1,mgroupsize))
        call MPI_Allgather(arwmin, 1, MPI_REAL, buffer, 1, MPI_REAL, MPI_COMM_WORLD, iok)
        arwmin = minval(buffer)
        deallocate(buffer)
     endif
#endif

     delxmin = sqrt(arwmin)

!    if (mdomain <= 1) then
!       radcone = erad * sin(op%coneang * pio180)
!       op%psiz = .05 * delxmin / (radcone * pio180 * (op%xmax - op%xmin))
!    else
        op%psiz = .07 * delxmin / (op%xmax - op%xmin)
!    endif
     op%vsprd = 50.  ! for vert xsectn

  endif

  ! Set limits of panel and frame window in plotter coordinates

  if ( op%projectn(iplt) == 'L' .and.  &
       op%windowin(iplt) /= 'W' ) then
     aspect = 0.
  else
     aspect = 1.
  endif
    
  call oplot_panel(op%panel(iplt),op%frameoff(iplt),op%colorbar(iplt),aspect,op%projectn(iplt))

  ! Scale plot window to selected model domain

  if (myrank == 0) then
     call o_set(op%h1,op%h2,op%v1,op%v2,op%xmin,op%xmax,op%ymin,op%ymax,1)
  endif

end subroutine oplot_set
