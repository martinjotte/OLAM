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

subroutine init_nudgrid()

  use misc_coms, only: iparallel
  implicit none

  if (iparallel == 1) then

     call para_init_nud()

  else

     call serial_init_nud()

  endif

end subroutine init_nudgrid

!===============================================================================


subroutine para_init_nud()

  use misc_coms,  only: io6, mdomain, isubdomain, iparallel

  use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                        itab_m_vars, itab_v_vars, itab_w_vars, &
                        itabg_m,     itabg_v,     itabg_w,     &
                        itab_w_pd,   &
                        alloc_itabs, mrls, &
                                  jtm_prog, jtm_wadj, jtm_vadj,           jtm_lbcp, &
                        jtv_init, jtv_prog, jtv_wadj,           jtv_wstn, jtv_lbcp, &
                        jtw_init, jtw_prog, jtw_wadj, jtw_vadj, jtw_wstn, jtw_lbcp

  use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, mza, &
                        alloc_gridz, alloc_xyzem, alloc_xyzew, &
                        alloc_grid1, alloc_grid2
  use max_dims,     only: maxremote
  use mem_para,     only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8
  use olam_mpi_nud, only: nsends_wnud, nrecvs_wnud, send_wnud, recv_wnud
  use mem_nudge,    only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, itabg_wnud, &
                          alloc_nudge1, xewnud, yewnud, zewnud

  implicit none

  integer :: j,imn,ivn,iwn,jnud
  integer :: im,iv,iw,iw1,iw2,iwnud,iwnud1,iwnud2,iwnud3
  integer :: imp,ivp,iwp,ips
  integer :: jsend, jrecv, jend
  integer :: npoly,nv
  integer :: iwnud_myrank ! Counter for WNUD points to be included on this rank
  integer :: nbytes_per_iwnud

  ! Automatic arrays

  logical :: myrankflag_wnud(nwnud)  ! Flag that ITABW(IW)%IWNUD(1:3) exist on
                                     ! this rank (for IW primary on this rank)
  logical :: myrankflag_wnud1(nwnud) ! Flag that ITABW(IW)%IWNUD(1) exists on
                                     ! this rank (for IW primary on this rank)

  integer :: send_wnud_iremotes(maxremote) ! holds ranks we will send to
  integer :: recv_wnud_iremotes(maxremote) ! holds ranks we will recv from

  ! Allocatable arrays

  logical, allocatable :: iwnud_sends(:,:) ! temporary send table
  integer, allocatable :: iwnud_recvs(:)   ! temporary recv table

  real, allocatable :: xewnudg(:), yewnudg(:), zewnudg(:)

  ! Skip if we are not doing spectral-like nudging

  if (.not. (mdomain == 0 .and. nudflag == 0 .and. nudnxp == 0)) return

  ! Initialize send & recv counters to zero

  nsends_wnud = 0
  nrecvs_wnud = 0

  ! Initialize myrank flag arrays to .false.

  myrankflag_wnud (:) = .false.
  myrankflag_wnud1(:) = .false.

  ! Loop over all W points, and for each whose assigned irank is equal to myrank,
  ! flag all WNUD points in its computational stencil for inclusion on this rank

  do iw = 2,nwa
     if (itabg_w(iw)%irank == myrank) then
        myrankflag_wnud ( itab_w_pd(iw)%iwnud(1:3) ) = .true.
        myrankflag_wnud1( itab_w_pd(iw)%iwnud(1)   ) = .true.
     endif
  enddo

  ! Ignore the "dummy" 1st point in case it was activated

  myrankflag_wnud (1) = .false.
  myrankflag_wnud1(1) = .false.

  ! Loop over all WNUD points and count the ones that have been flagged
  ! for inclusion on this rank.

  iwnud_myrank = 1

  do iwnud = 2,nwnud
     if (myrankflag_wnud(iwnud)) iwnud_myrank = iwnud_myrank + 1
  enddo

  mwnud = iwnud_myrank

  ! Allocate grid structure variables

  call move_alloc(xewnud, xewnudg)
  call move_alloc(yewnud, yewnudg)
  call move_alloc(zewnud, zewnudg)

  call alloc_nudge1(mwnud,2)

  ! Reset point counts to 1

  iwnud_myrank = 1

  ! Store new WNUD indices in itabg data structures

  do iwnud = 2, nwnud
     if (myrankflag_wnud(iwnud)) then
        iwnud_myrank = iwnud_myrank + 1

        itabg_wnud(iwnud)%iwnud_myrank = iwnud_myrank
        itab_wnud(iwnud_myrank)%iwnudglobe = iwnud
        itab_wnud(iwnud_myrank)%irank = itabg_wnud(iwnud)%irank

        xewnud(iwnud_myrank) = xewnudg(iwnud)
        yewnud(iwnud_myrank) = yewnudg(iwnud)
        zewnud(iwnud_myrank) = zewnudg(iwnud)
     endif
  enddo

  deallocate (xewnudg, yewnudg, zewnudg)

  ! Read the grid structure for all points in local parallel subdomain
  ! for this rank (or for all points in domain if run is sequential)

  call gridfile_read_nudge()

  ! Convert global itab_w nudging index information to local rank

  do iw = 2, mwa ! local index

     ! Set local indices of WNUD points in W stencil
     do jnud = 1,3
        iwnud = itab_w(iw)%iwnud(jnud) ! Global index
        itab_w(iw)%iwnud(jnud) = itabg_wnud(iwnud)%iwnud_myrank
     enddo

  enddo

  ! Initialize temporary send/recv tables

  ips = min(maxremote, mgroupsize-1)
  allocate(iwnud_sends(mwnud,ips)) ; iwnud_sends = .false.

  allocate(iwnud_recvs(mwnud)) ; iwnud_recvs = -1

  ! Fill temporary MPI SEND/RECV tables

  do iw = 2, nwa ! global index
     if (itabg_w(iw)%irank /= myrank) then

        ! IW point is primary on remote rank.

        ! Although a primary rank is assigned to each nudging point, the assignment is
        ! used only to designate which rank writes values for that nudging point to
        ! history files.  Parallel MPI communication, on the other hand, must take into
        ! account that computation of nudging point values is based on sums over model
        ! grid points on multiple ranks.  Each of the ranks that contribute to a given
        ! nudging point accumulate a partial sum and send the partial sum to all other
        ! ranks that need that nudging point value.  After receiving the partial sums,
        ! each rank completes the sum (a small duplication of effort) over the partial
        ! sums to obtain its own copy of the full sum.  Thus, for building send and
        ! recv tables, we check itab_w_pd(iw)%iwnud(1:3) and itab_w_pd(iw)%iwnud(1)
        ! instead of itab_wnud()%irank.

        iwnud1 = itab_w_pd(iw)%iwnud(1)
        iwnud2 = itab_w_pd(iw)%iwnud(2)
        iwnud3 = itab_w_pd(iw)%iwnud(3)

        if (myrankflag_wnud1(iwnud1)) call send_table_wnud(iwnud1,itabg_w(iw)%irank)
        if (myrankflag_wnud1(iwnud2)) call send_table_wnud(iwnud2,itabg_w(iw)%irank)
        if (myrankflag_wnud1(iwnud3)) call send_table_wnud(iwnud3,itabg_w(iw)%irank)

        if (myrankflag_wnud(iwnud1)) call recv_table_wnud(iwnud1,itabg_w(iw)%irank)

     endif
  enddo

  ! Determine number of bytes to send per IW column

  nbytes_per_iwnud = 6 * mza * nbytes_real8

  ! Copy data from temp arrays to w send table

  allocate( send_wnud(nsends_wnud) )

  do jsend = 1, nsends_wnud
     jend = count( iwnud_sends(2:mwnud,jsend) )

     send_wnud(jsend)%jend    = jend
     send_wnud(jsend)%iremote = send_wnud_iremotes(jsend)
     send_wnud(jsend)%nbytes  = jend * nbytes_per_iwnud

     allocate( send_wnud(jsend)%ipts( send_wnud(jsend)%jend   ) )
     allocate( send_wnud(jsend)%buff( send_wnud(jsend)%nbytes ) )

     j = 0
     do iwnud = 2, mwnud
        if (iwnud_sends(iwnud,jsend)) then
           j = j + 1
           send_wnud(jsend)%ipts(j) = iwnud
        endif
     enddo
  enddo

  ! Copy data from temp arrays to wnud recv table

  allocate( recv_wnud(nrecvs_wnud) )

  do jrecv = 1, nrecvs_wnud
     jend = count( iwnud_recvs(2:mwnud) == jrecv )

     recv_wnud(jrecv)%jend    = jend
     recv_wnud(jrecv)%iremote = recv_wnud_iremotes(jrecv)
     recv_wnud(jrecv)%nbytes  = jend * nbytes_per_iwnud

     allocate( recv_wnud(jrecv)%ipts( recv_wnud(jrecv)%jend   ) )
     allocate( recv_wnud(jrecv)%buff( recv_wnud(jrecv)%nbytes ) )

     j = 0
     do iwnud = 2, mwnud
        if (iwnud_recvs(iwnud) == jrecv) then
           j = j + 1
           recv_wnud(jrecv)%ipts(j) = iwnud
        endif
     enddo
  enddo

  ! Deallocate para_decomp _pd arrays

  if (allocated(itab_w_pd)) deallocate(itab_w_pd)


Contains


  subroutine recv_table_wnud(iwnud, iremote)

    use mem_nudge,    only: itabg_wnud
    use olam_mpi_nud, only: nrecvs_wnud

    implicit none

    integer, intent(in) :: iwnud   ! global index
    integer, intent(in) :: iremote ! remote rank
    integer             :: jrecv

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_wnud
       if (recv_wnud_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_wnud, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_wnud.

    if (jrecv > nrecvs_wnud) then
       nrecvs_wnud = jrecv
       recv_wnud_iremotes(jrecv) = iremote
    endif

    ! Enter remote rank in recv-remote-rank table.

    iwnud_recvs( itabg_wnud(iwnud)%iwnud_myrank ) = jrecv

  end subroutine recv_table_wnud


  subroutine send_table_wnud(iwnud, iremote)

    use mem_nudge,    only: itabg_wnud
    use olam_mpi_nud, only: nsends_wnud

    implicit none

    integer,  intent(in) :: iwnud    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_w is already in table of ranks to send to

    do jsend = 1, nsends_wnud
       if (send_wnud_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_wnud, jsend represents a rank not yet entered in
    ! the table, so increase nsends_w.

    if (jsend > nsends_wnud) then
       nsends_wnud = jsend
       send_wnud_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table

    iwnud_sends( itabg_wnud(iwnud)%iwnud_myrank, jsend ) = .true.

  end subroutine send_table_wnud


end subroutine para_init_nud

!===============================================================================

subroutine serial_init_nud()

  use misc_coms,  only: io6, mdomain, isubdomain, iparallel

  use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                        itab_m_vars, itab_v_vars, itab_w_vars, &
                        itabg_m,     itabg_v,     itabg_w,     &
                        itab_w_pd,   &
                        alloc_itabs, mrls, &
                                  jtm_prog, jtm_wadj, jtm_vadj,           jtm_lbcp, &
                        jtv_init, jtv_prog, jtv_wadj,           jtv_wstn, jtv_lbcp, &
                        jtw_init, jtw_prog, jtw_wadj, jtw_vadj, jtw_wstn, jtw_lbcp

  use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, mza, &
                        alloc_gridz, alloc_xyzem, alloc_xyzew, &
                        alloc_grid1, alloc_grid2
  use max_dims,     only: maxremote
  use mem_para,     only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8
  use olam_mpi_nud, only: nsends_wnud, nrecvs_wnud, send_wnud, recv_wnud
  use mem_nudge,    only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, itabg_wnud, &
                          alloc_nudge1, xewnud, yewnud, zewnud

  implicit none

  integer :: j,imn,ivn,iwn,jnud
  integer :: im,iv,iw,iw1,iw2,iwnud,iwnud1,iwnud2,iwnud3
  integer :: imp,ivp,iwp,ips
  integer :: jsend, jrecv, jend
  integer :: npoly,nv
  integer :: iwnud_myrank ! Counter for WNUD points to be included on this rank
  integer :: nbytes_per_iwnud

  ! Skip if we are not doing spectral-like nudging

  if (.not. (mdomain == 0 .and. nudflag == 0 .and. nudnxp == 0)) return

  ! Allocate grid structure variables

  call alloc_nudge1(mwnud,2)

  ! Store new WNUD indices in itabg data structures

  do iwnud = 1, nwnud
     itabg_wnud(iwnud)%iwnud_myrank = iwnud
     itab_wnud (iwnud)%iwnudglobe   = iwnud
     itab_wnud (iwnud)%irank        = 0
  enddo

  ! Read the grid structure for all points in local parallel subdomain
  ! for this rank (or for all points in domain if run is sequential)

  call gridfile_read_nudge()

end subroutine serial_init_nud
