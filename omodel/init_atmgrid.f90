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

subroutine init_atmgrid()

  use misc_coms, only: iparallel
  implicit none

  if (iparallel == 1) then

     call para_init_atm()

  else

     call serial_init_atm()

  endif

end subroutine init_atmgrid

!===============================================================================

subroutine para_init_atm()

  use misc_coms,  only: mdomain

  use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                        itab_m_vars, itab_v_vars, itab_w_vars, &
                        itabg_m,     itabg_v,     itabg_w,     &
                        itab_m_pd,   itab_v_pd,   itab_w_pd,   &
                        alloc_itabs, &
                                  jtm_prog, jtm_wadj, jtm_vadj,           jtm_lbcp, &
                        jtv_init, jtv_prog, jtv_wadj,           jtv_wstn, jtv_lbcp, &
                        jtw_init, jtw_prog, jtw_wadj, jtw_vadj, jtw_wstn, jtw_lbcp

  use mem_grid,   only: nma, nva, nwa, mma, mua, mva, mwa, mza, &
                        alloc_gridz, alloc_xyzem, alloc_xyzew, &
                        alloc_grid1, alloc_grid2, xew, yew, zew, arw0

  use mem_para,     only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8
  use olam_mpi_atm, only: nsends_v, nrecvs_v, send_v, recv_v, &
                          nsends_m, nrecvs_m, send_m, recv_m, &
                          nsends_w, nrecvs_w, send_w, recv_w
  use max_dims,   only: maxremote
  use var_tables, only: nvar_par
  use mem_nudge,  only: nudflag, nudnxp

  implicit none

  integer :: j,imn,ivn,iwn
  integer :: im,iv,iw,iw1,iw2
  integer :: imp,ivp,iwp,ips
  integer :: jsend, jrecv, jend
  integer :: npoly,nv

  integer :: im_myrank ! Counter for M points to be included on this rank
  integer :: iv_myrank ! Counter for V points to be included on this rank
  integer :: iw_myrank ! Counter for W points to be included on this rank

  integer :: nbytes_per_iw
  integer :: nbytes_per_iv
  integer :: nbytes_per_im

  ! Automatic arrays

  logical :: myrankflag_m(nma) ! Flag for M points existing on this rank
  logical :: myrankflag_v(nva) ! Flag for V points existing on this rank
  logical :: myrankflag_w(nwa) ! Flag for W points existing on this rank

  integer :: send_v_iremotes(maxremote) ! holds ranks we will send to
  integer :: send_w_iremotes(maxremote) ! holds ranks we will send to
  integer :: send_m_iremotes(maxremote) ! holds ranks we will send to

  integer :: recv_v_iremotes(maxremote) ! holds ranks we will recv from
  integer :: recv_w_iremotes(maxremote) ! holds ranks we will recv from
  integer :: recv_m_iremotes(maxremote) ! holds ranks we will recv from

  ! Allocatable arrays

  logical, allocatable :: iv_sends(:,:) ! temporary send table
  logical, allocatable :: iw_sends(:,:) ! temporary send table
  logical, allocatable :: im_sends(:,:) ! temporary send table

  integer :: ivp_recv_jrecv(nva) ! temporary recv table
  integer :: ivp_recv_iv   (nva) ! temporary recv table

  integer :: iwp_recv_jrecv(nwa) ! temporary recv table
  integer :: iwp_recv_iw   (nwa) ! temporary recv table

  integer :: imp_recv_jrecv(nma) ! temporary recv table
  integer :: imp_recv_im   (nma) ! temporary recv table

  real, allocatable :: xewg(:), yewg(:), zewg(:), arw0g(:)

  ! Initialize send & recv counters to zero

  nsends_v = 0
  nsends_w = 0
  nsends_m = 0

  nrecvs_v = 0
  nrecvs_w = 0
  nrecvs_m = 0

  ! Initialize myrank flag arrays to .false.

  myrankflag_m(:) = .false.
  myrankflag_v(:) = .false.
  myrankflag_w(:) = .false.

  ! Loop over all V points, and for each whose assigned irank is equal to myrank,
  ! flag all M, V, and W points in its computational stencil for inclusion on this
  ! rank, excluding IVP.

  do iv = 2,nva
     if (itabg_v(iv)%irank == myrank) then

        myrankflag_v( iv ) = .true.
        myrankflag_m( itab_v_pd(iv)%im(1:2) ) = .true.
        myrankflag_w( itab_v_pd(iv)%iw(1:2) ) = .true.

     endif
  enddo

  ! Loop over all W points, and for each whose assigned irank is equal to myrank,
  ! flag all M, V, and W points in its computational stencil for
  ! inclusion on this rank, excluding IWP.

  do iw = 2,nwa
     if (itabg_w(iw)%irank == myrank) then

        npoly = itab_w_pd(iw)%npoly

        myrankflag_w( iw ) = .true.
        myrankflag_m( itab_w_pd(iw)%im(1:npoly) ) = .true.
        myrankflag_v( itab_w_pd(iw)%iv(1:npoly) ) = .true.
        myrankflag_w( itab_w_pd(iw)%iw(1:npoly) ) = .true.

     endif
  enddo

  ! Loop over all M points, and for each whose assigned irank is equal to myrank,
  ! flag all M, V, and W points in its computational stencil for inclusion on this
  ! rank, excluding IMP.

  do im = 2,nma
     if (itabg_m(im)%irank == myrank) then

        npoly = itab_m_pd(im)%npoly

        myrankflag_m( im ) = .true.
        myrankflag_m( itab_m_pd(im)%im(1:npoly) ) = .true.
        myrankflag_v( itab_m_pd(im)%iv(1:npoly) ) = .true.
        myrankflag_w( itab_m_pd(im)%iw(1:npoly) ) = .true.

     endif
  enddo

  ! Ignore the "dummy" 1st point in case it was activated

  myrankflag_m(1) = .false.
  myrankflag_v(1) = .false.
  myrankflag_w(1) = .false.

  ! Loop over all M, V, and W points and count the ones that have
  ! been flagged for inclusion on this rank.

  im_myrank = 1
  iv_myrank = 1
  iw_myrank = 1

  do im = 2,nma
     if (myrankflag_m(im)) im_myrank = im_myrank + 1
  enddo

  do iv = 2,nva
     if (myrankflag_v(iv)) iv_myrank = iv_myrank + 1
  enddo

  do iw = 2,nwa
     if (myrankflag_w(iw)) iw_myrank = iw_myrank + 1
  enddo

  ! Set mma, mva, mwa values for this rank

  mma = im_myrank
  mva = iv_myrank
  mua = mva
  mwa = iw_myrank

  ! Allocate grid structure variables

  call move_alloc(xew,xewg)
  call move_alloc(yew,yewg)
  call move_alloc(zew,zewg)
  call move_alloc(arw0,arw0g)

  call alloc_gridz()
  call alloc_itabs(mma, mva, mwa, 1)
  call alloc_xyzem(mma)
  call alloc_xyzew(mwa)
  call alloc_grid1(mma, mva, mwa)
  call alloc_grid2(mma, mva, mwa)

  ! Reset point counts to 1

  im_myrank = 1
  iv_myrank = 1
  iw_myrank = 1

  ! Store new myrank M, V, and W indices in itabg data structures

  do im = 2, nma
     if (myrankflag_m(im)) then
        im_myrank = im_myrank + 1

        itabg_m(im)%im_myrank = im_myrank
        itab_m(im_myrank)%imglobe = im
        itab_m(im_myrank)%irank = itabg_m(im)%irank
     endif
  enddo

  do iv = 2, nva
     if (myrankflag_v(iv)) then
        iv_myrank = iv_myrank + 1

        itabg_v(iv)%iv_myrank = iv_myrank
        itab_v(iv_myrank)%ivglobe = iv
        itab_v(iv_myrank)%irank = itabg_v(iv)%irank
     endif
  enddo

  do iw = 2, nwa
     if (myrankflag_w(iw)) then
        iw_myrank = iw_myrank + 1

        itabg_w(iw)%iw_myrank = iw_myrank
        itab_w(iw_myrank)%iwglobe = iw
        itab_w(iw_myrank)%irank = itabg_w(iw)%irank

        xew (iw_myrank) = xewg (iw)
        yew (iw_myrank) = yewg (iw)
        zew (iw_myrank) = zewg (iw)
        arw0(iw_myrank) = arw0g(iw)
     endif
  enddo

  deallocate (xewg, yewg, zewg, arw0g)

  ! Read the grid structure for all points in local parallel subdomain
  ! for this rank (or for all points in domain if run is sequential)

  call gridfile_read()

  ! Copy and convert global itab_m_pd index information to local rank

  do im = 2, nma ! global index

     if (myrankflag_m(im)) then ! M point is in memory on local subdomain

        im_myrank = itabg_m(im)%im_myrank ! Local index of M point
        npoly     = itab_m_pd(im)%npoly

        itab_m(im_myrank)%npoly = npoly

        ! Set indices of IM neighbors that are present on this rank

        do j = 1, npoly
           ivn = itab_m_pd(im)%iv(j) ! Global index
           iwn = itab_m_pd(im)%iw(j) ! Global index
           imn = itab_m_pd(im)%im(j) ! Global index

           if (myrankflag_v(ivn)) &
                itab_m(im_myrank)%iv(j) = itabg_v(ivn)%iv_myrank

           if (myrankflag_w(iwn)) &
                itab_m(im_myrank)%iw(j) = itabg_w(iwn)%iw_myrank

           if (myrankflag_m(imn)) &
                itab_m(im_myrank)%im(j) = itabg_m(imn)%im_myrank
        enddo

        if (itabg_m(im)%irank /= myrank) then

           ! Turn off jtm_prog loop flag if M point is not primary
           ! on this subdomain

           call mloopf('n', im_myrank, -jtm_prog, 0, 0, 0, 0, 0)

           ! Turn off jtm_vadj loop flag if no adjacent V faces are primary
           ! on this subdomain

           if ( all( itabg_v( itab_m_pd(im)%iv(1:npoly) )%irank /= myrank ) ) then
              call mloopf('n', im_myrank, -jtm_vadj, 0, 0, 0, 0, 0)
           endif

           ! Turn off jtm_wadj loop flag if no adjacent W faces are primary
           ! on this subdomain

           if ( all( itabg_w( itab_m_pd(im)%iw(1:npoly) )%irank /= myrank ) ) then
              call mloopf('n', im_myrank, -jtm_wadj, 0, 0, 0, 0, 0)
           endif

        endif

        ! If IMP point exists on local parallel subdomain, assign its local
        ! index to IM stencil; otherwise set index to 1 and turn off lbcp flag

        imp = itab_m_pd(im)%imp ! Global index

        if (myrankflag_m(imp)) then
           itab_m(im_myrank)%imp = itabg_m(imp)%im_myrank ! Local index
        else
           itab_m(im_myrank)%imp = 1
           call mloopf('n', im_myrank, -jtm_lbcp, 0, 0, 0, 0, 0)
        endif

     endif
  enddo

  ! Copy and convert global itab_v_pd index information to local rank

  do iv = 2, nva ! global index

     if (myrankflag_v(iv)) then ! V point is in memory on local subdomain

        iv_myrank = itabg_v(iv)%iv_myrank ! Local index of V point

        ! Set indices of IV neighbors that are present on this rank

        do j = 1, 2
           imn = itab_v_pd(iv)%im(j) ! Global index
           iwn = itab_v_pd(iv)%iw(j) ! Global index

           if (myrankflag_m(imn)) &
                itab_v(iv_myrank)%im(j) = itabg_m(imn)%im_myrank

           if (myrankflag_w(iwn)) &
                itab_v(iv_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
        enddo

        if (itabg_v(iv)%irank /= myrank) then

           ! If IV point is not primary on local parallel subdomain,
           ! turn off jtv_init and jtv_prog flags

           call vloopf('n',iv_myrank, -jtv_init, -jtv_prog, 0, 0, 0, 0)

           ! If adjacent W points are both not primary on local parallel
           !  subdomain, turn off jtv_wadj and jtv_wstn flags

           iw1 = itab_v_pd(iv)%iw(1) ! Global index
           iw2 = itab_v_pd(iv)%iw(2) ! Global index

           if (itabg_w(iw1)%irank /= myrank .and. itabg_w(iw2)%irank /= myrank) then
              call vloopf('n',iv_myrank, -jtv_wadj, -jtv_wstn, 0, 0, 0, 0)
           endif

        endif

        ! If IVP point exists on local parallel subdomain, assign its local
        ! index to IV stencil; otherwise set index to 1 and turn off lbcp flag

        ivp = itab_v_pd(iv)%ivp ! Global index

        if (myrankflag_v(ivp)) then
           itab_v(iv_myrank)%ivp = itabg_v(ivp)%iv_myrank
        else
           itab_v(iv_myrank)%ivp = 1
           call vloopf('n',iv_myrank, -jtv_lbcp, 0, 0, 0, 0, 0)
        endif

     endif
  enddo

  ! Copy and convert global itab_w_pd index information to local rank

  do iw = 2, nwa ! global index

     if (myrankflag_w(iw)) then ! W point is in memory on local subdomain

        iw_myrank = itabg_w(iw)%iw_myrank ! Local index of W point
        npoly     = itab_w_pd(iw)%npoly

        itab_w(iw_myrank)%npoly = npoly

        ! Set local indices of M, V, W points in W stencil

        do j = 1, npoly
           imn = itab_w_pd(iw)%im(j) ! Global index
           ivn = itab_w_pd(iw)%iv(j) ! Global index
           iwn = itab_w_pd(iw)%iw(j) ! Global index

           if (myrankflag_m(imn)) &
                itab_w(iw_myrank)%im(j) = itabg_m(imn)%im_myrank

           if (myrankflag_v(ivn)) &
                itab_w(iw_myrank)%iv(j) = itabg_v(ivn)%iv_myrank

           if (myrankflag_w(iwn)) &
                itab_w(iw_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
        enddo

        if (itabg_w(iw)%irank /= myrank) then

           ! If IW point is not primary on local parallel subdomain,
           ! turn off jtw_init and jtw_prog flags

           call wloopf('n',iw_myrank, -jtw_init, -jtw_prog, 0, 0, 0, 0)

           ! If none of the adjacent W points is primary on local parallel
           ! subdomain, turn off jtw_wadj and jtw_wstn flags

           if ( all( itabg_w( itab_w_pd(iw)%iw(1:npoly) )%irank /= myrank ) ) then
              call wloopf('n',iw_myrank, -jtw_wadj, -jtw_wstn, 0, 0, 0, 0)
           endif

           ! If none of the adjacent V points is primary on local parallel
           ! subdomain, turn off jtw_vadj flag

           if ( all( itabg_v( itab_w_pd(iw)%iv(1:npoly) )%irank /= myrank ) ) then
              call wloopf('n',iw_myrank, -jtw_vadj, 0, 0, 0, 0, 0)
           endif

        endif

        ! If IWP point exists on local parallel subdomain, assign its local
        ! index to IW stencil; otherwise set index to 1 and turn off lbcp flag

        iwp = itab_w_pd(iw)%iwp ! Global index

        if (myrankflag_w(iwp)) then
           itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank
        else
           itab_w(iw_myrank)%iwp = 1
           call wloopf('n',iw_myrank, -jtw_lbcp, 0, 0, 0, 0, 0)
        endif

     endif
  enddo

  ! Special handling of boprder copy points for periodic limited-area runs

!!  do iw = 2, nwa
!!     if (myrankflag_w(iw)) then
!!        iwp = itab_w_pd(iw)%iwp
!!        if (.not. myrankflag_w(iwp)) then
!!           if (itab_w_pd(iwp)%iw_myrank_iwp == -1) then
!!              ! never been set, we will MPI copy to this point
!!              itab_w_pd(iwp)%iw_myrank_iwp = itabg_w(iw)%iw_myrank
!!           else
!!              ! We are already MPI copying this point, so do an LBC copy instead
!!              iw_myrank  = itabg_w(iw)%iw_myrank
!!              itab_w(iw_myrank)%iwp = itab_w_pd(iwp)%iw_myrank_iwp
!!              call wloopf('n', iw_myrank, jtw_lbcp, 0, 0, 0, 0, 0)
!!           endif
!!        endif
!!     endif
!!  enddo
!!
!!  do iv = 2, nva
!!     if (myrankflag_v(iv)) then
!!        ivp = itab_v_pd(iv)%ivp
!!        if (.not. myrankflag_v(ivp)) then
!!           if (itab_v_pd(ivp)%iv_myrank_ivp == -1) then
!!              ! never been set, we will MPI copy to this point
!!              itab_v_pd(ivp)%iv_myrank_ivp = itabg_v(iv)%iv_myrank
!!           else
!!              ! We are already MPI copying this point, so do an LBC copy instead
!!              iv_myrank  = itabg_v(iv)%iv_myrank
!!              itab_v(iv_myrank)%ivp = itab_v_pd(ivp)%iv_myrank_ivp
!!              call vloopf('n', iv_myrank, jtv_lbcp, 0, 0, 0, 0, 0)
!!           endif
!!        endif
!!     endif
!!  enddo
!!
!!  do im = 2, nma
!!     if (myrankflag_m(im)) then
!!        imp = itab_m_pd(im)%imp
!!        if (.not. myrankflag_m(imp)) then
!!           if (itab_m_pd(imp)%im_myrank_imp == -1) then
!!              ! never been set, we will MPI copy to this point
!!              itab_m_pd(imp)%im_myrank_imp = itabg_m(im)%im_myrank
!!           else
!!              ! We are already MPI copying this point, so do an LBC copy instead
!!              im_myrank  = itabg_m(im)%im_myrank
!!              itab_m(im_myrank)%imp = itab_m_pd(imp)%im_myrank_imp
!!              call mloopf('n', im_myrank, jtm_lbcp, 0, 0, 0, 0, 0)
!!           endif
!!        endif
!!     endif
!!  enddo


  ! Allocate/initialize temporary send/recv tables

  ips = min(maxremote, mgroupsize)

  allocate(iv_sends(mva,ips)) ; iv_sends = .false.
  allocate(im_sends(mma,ips)) ; im_sends = .false.
  allocate(iw_sends(mwa,ips)) ; iw_sends = .false.

  ivp_recv_jrecv = 0
  ivp_recv_iv    = 0

  iwp_recv_jrecv = 0
  iwp_recv_iw    = 0

  imp_recv_jrecv = 0
  imp_recv_im    = 0

  ! Fill temporary MPI send/receive tables

  do iv = 2, nva ! global index
     if (itabg_v(iv)%irank /= myrank) then

        ! IV point is primary on remote rank.

        ! Add to send table any V or W point that is primary on myrank and is
        ! in the stencil of IV.

        ivp = itab_v_pd(iv)%ivp
        if (itabg_v(ivp)%irank == myrank) call send_table_v(ivp,itabg_v(iv)%irank)

        do j = 1,2
           iwn = itab_v_pd(iv)%iw(j)
           iwn = itab_w_pd(iwn)%iwp
           if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_v(iv)%irank)

           imn = itab_v_pd(iv)%im(j)
           imn = itab_m_pd(imn)%imp
           if (itabg_m(imn)%irank == myrank) call send_table_m(imn,itabg_v(iv)%irank)
        enddo

     endif

     ! If IV point is in memory of myrank but its IVP is primary on another rank,
     ! add the remote rank to the receive table

     if (myrankflag_v(iv)) then
        ivp = itab_v_pd(iv)%ivp

        if (itabg_v(ivp)%irank /= myrank) then
           iv_myrank = itabg_v(iv)%iv_myrank

           if (itab_v_pd(ivp)%iv_myrank_ivp < 1) then

              ! Never been set, MPI copy to this point and turn off any LBC copy
              itab_v_pd(ivp)%iv_myrank_ivp = iv_myrank
              call recv_table_v(iv_myrank, ivp)
              call vloopf('n', iv_myrank, -jtv_lbcp, 0, 0, 0, 0, 0)

           else

              ! We are already MPI copying this point, so do an LBC copy instead
              itab_v(iv_myrank)%ivp = itab_v_pd(ivp)%iv_myrank_ivp
              call vloopf('n', iv_myrank, jtv_lbcp, 0, 0, 0, 0, 0)

           endif

        endif
     endif

  enddo

  do iw = 2,nwa ! global index
     if (itabg_w(iw)%irank /= myrank) then

        ! IW point is primary on remote rank.

        ! Add to send table any V, W, or M point that is primary on myrank and is
        ! in the stencil of IW.

        iwp = itab_w_pd(iw)%iwp
        if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

        npoly = itab_w_pd(iw)%npoly

        do j = 1, npoly

           ivn = itab_w_pd(iw)%iv(j)
           ivn = itab_v_pd(ivn)%ivp
           if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_w(iw)%irank)

           iwn = itab_w_pd(iw)%iw(j)
           iwn = itab_w_pd(iwn)%iwp
           if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_w(iw)%irank)

           imn = itab_w_pd(iw)%im(j)
           imn = itab_m_pd(imn)%imp
           if (itabg_m(imn)%irank == myrank) call send_table_m(imn,itabg_w(iw)%irank)

        enddo

     endif

     ! If IW point is in memory of myrank but its IWP is primary on another rank,
     ! add the remote rank to the receive table

     if (myrankflag_w(iw)) then
        iwp = itab_w_pd(iw)%iwp

        if (itabg_w(iwp)%irank /= myrank) then
           iw_myrank = itabg_w(iw)%iw_myrank

           if (itab_w_pd(iwp)%iw_myrank_iwp < 1) then

              ! Never been set, MPI copy to this point and turn off any LBC copy
              itab_w_pd(iwp)%iw_myrank_iwp = iw_myrank
              call recv_table_w(iw_myrank, iwp)
              call wloopf('n', iw_myrank, -jtw_lbcp, 0, 0, 0, 0, 0)

           else

              ! We are already MPI copying this point, so do an LBC copy instead
              itab_w(iw_myrank)%iwp = itab_w_pd(iwp)%iw_myrank_iwp
              call wloopf('n', iw_myrank, jtw_lbcp, 0, 0, 0, 0, 0)

           endif

        endif
     endif

  enddo

  do im = 2, nma ! global index
     if (itabg_m(im)%irank /= myrank) then

        ! IM point is primary on remote rank.

        ! Add to send table any V or M point that is primary on myrank and is
        ! in the stencil of IM.

        imp = itab_m_pd(im)%imp
        if (itabg_m(imp)%irank == myrank) call send_table_m(imp,itabg_m(im)%irank)

        npoly = itab_m_pd(im)%npoly

        do j = 1, npoly
           ivn = itab_m_pd(im)%iv(j)
           ivn = itab_v_pd(ivn)%ivp
           if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_m(im)%irank)

           iwn = itab_m_pd(im)%iw(j)
           iwn = itab_w_pd(iwn)%iwp
           if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_m(im)%irank)

           imn = itab_m_pd(im)%im(j)
           imn = itab_m_pd(imn)%imp
           if (itabg_m(imn)%irank == myrank) call send_table_m(imn,itabg_m(im)%irank)
        enddo

     endif

     ! If IM point is in memory of myrank but its IMP is primary on another rank,
     ! add the remote rank to the receive table

     if (myrankflag_m(im)) then
        imp = itab_m_pd(im)%imp

        if (itabg_m(imp)%irank /= myrank) then
           im_myrank = itabg_m(im)%im_myrank

           if (itab_m_pd(imp)%im_myrank_imp < 1) then

              ! Never been set, MPI copy to this point and turn off any LBC copy
              itab_m_pd(imp)%im_myrank_imp = im_myrank
              call recv_table_m(im_myrank, imp)
              call mloopf('n', im_myrank, -jtm_lbcp, 0, 0, 0, 0, 0)

           else

              ! We are already MPI copying this point, so do an LBC copy instead
              itab_m(im_myrank)%imp = itab_m_pd(imp)%im_myrank_imp
              call mloopf('n', im_myrank, jtm_lbcp, 0, 0, 0, 0, 0)

           endif

        endif
     endif

  enddo

  ! Determine number of bytes to send per IW column

  nv = max(nvar_par, 20)

  nbytes_per_iw =  2 * mza * nbytes_real8 &
                + nv * mza * nbytes_real

  ! Copy data from temp arrays to w send table

  allocate( send_w(nsends_w) )

  do jsend = 1, nsends_w
     jend = count( iw_sends(2:mwa,jsend) )

     send_w(jsend)%jend    = jend
     send_w(jsend)%iremote = send_w_iremotes(jsend)
     send_w(jsend)%nbytes  = jend * nbytes_per_iw

     allocate( send_w(jsend)%ipts( send_w(jsend)%jend   ) )
     allocate( send_w(jsend)%buff( send_w(jsend)%nbytes ) )

     j = 0
     do iw = 2, mwa
        if (iw_sends(iw,jsend)) then
           j = j + 1
           send_w(jsend)%ipts(j) = iw
        endif
     enddo
  enddo

  ! Determine number of bytes to send per IV column

  nbytes_per_iv = 3       * nbytes_int &
                + 4 * mza * nbytes_real

  ! Copy data from temp arrays to v send table

  allocate( send_v(nsends_v) )

  do jsend = 1, nsends_v
     jend = count( iv_sends(2:mva,jsend) )

     send_v(jsend)%jend    = jend
     send_v(jsend)%iremote = send_v_iremotes(jsend)
     send_v(jsend)%nbytes  = jend * nbytes_per_iv

     allocate( send_v(jsend)%ipts( send_v(jsend)%jend   ) )
     allocate( send_v(jsend)%buff( send_v(jsend)%nbytes ) )

     j = 0
     do iv = 2, mva
        if (iv_sends(iv,jsend)) then
           j = j + 1
           send_v(jsend)%ipts(j) = iv
        endif
     enddo
  enddo

  ! Determine number of bytes to send per IM column

  nbytes_per_im = 2 * mza * nbytes_real

  ! Copy data from temp arrays to m send table

  allocate( send_m(nsends_m) )

  do jsend = 1, nsends_m
     jend = count( im_sends(2:mma,jsend) )

     send_m(jsend)%jend    = jend
     send_m(jsend)%iremote = send_m_iremotes(jsend)
     send_m(jsend)%nbytes  = jend * nbytes_per_im

     allocate( send_m(jsend)%ipts( send_m(jsend)%jend   ) )
     allocate( send_m(jsend)%buff( send_m(jsend)%nbytes ) )

     j = 0
     do im = 2, mma
        if (im_sends(im,jsend)) then
           j = j + 1
           send_m(jsend)%ipts(j) = im
        endif
     enddo
  enddo

  ! Copy data from temp arrays to w recv table

  allocate( recv_w(nrecvs_w) )

  do jrecv = 1, nrecvs_w
     jend = count( iwp_recv_jrecv(2:nwa) == jrecv )

     recv_w(jrecv)%jend    = jend
     recv_w(jrecv)%iremote = recv_w_iremotes(jrecv)
     recv_w(jrecv)%nbytes  = jend * nbytes_per_iw

     allocate( recv_w(jrecv)%ipts( recv_w(jrecv)%jend   ) )
     allocate( recv_w(jrecv)%buff( recv_w(jrecv)%nbytes ) )

     j = 0
     do iw = 2, nwa
        if (iwp_recv_jrecv(iw) == jrecv) then
           j = j + 1
           recv_w(jrecv)%ipts(j) = iwp_recv_iw(iw)
        endif
     enddo
  enddo

  ! Copy data from temp arrays to v recv table

  allocate( recv_v(nrecvs_v) )

  do jrecv = 1, nrecvs_v
     jend = count( ivp_recv_jrecv(2:nva) == jrecv )

     recv_v(jrecv)%jend    = jend
     recv_v(jrecv)%iremote = recv_v_iremotes(jrecv)
     recv_v(jrecv)%nbytes  = jend * nbytes_per_iv

     allocate( recv_v(jrecv)%ipts( recv_v(jrecv)%jend   ) )
     allocate( recv_v(jrecv)%buff( recv_v(jrecv)%nbytes ) )

     j = 0
     do iv = 2, nva
        if (ivp_recv_jrecv(iv) == jrecv) then
           j = j + 1
           recv_v(jrecv)%ipts(j) = ivp_recv_iv(iv)
        endif
     enddo
  enddo

  ! Copy data from temp arrays to m recv table

  allocate( recv_m(nrecvs_m) )

  do jrecv = 1, nrecvs_m
     jend = count( imp_recv_jrecv(2:nma) == jrecv )

     recv_m(jrecv)%jend    = jend
     recv_m(jrecv)%iremote = recv_m_iremotes(jrecv)
     recv_m(jrecv)%nbytes  = jend * nbytes_per_im

     allocate( recv_m(jrecv)%ipts( recv_m(jrecv)%jend   ) )
     allocate( recv_m(jrecv)%buff( recv_m(jrecv)%nbytes ) )

     j = 0
     do im = 2, nma
        if (imp_recv_jrecv(im) == jrecv) then
           j = j + 1
           recv_m(jrecv)%ipts(j) = imp_recv_im(im)
        endif
     enddo
  enddo

  ! Deallocate para_decomp _pd arrays

  deallocate( itab_m_pd )
  deallocate( itab_v_pd )
  if (.not. (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0)) then
     deallocate( itab_w_pd )
  endif


Contains


  subroutine recv_table_w(iw, iwp_globe)

    use mem_ijtabs,   only: itabg_w
    use olam_mpi_atm, only: nrecvs_w

    implicit none

    integer, intent(in)  :: iw, iwp_globe
    integer              :: jrecv, iremote

    iremote = itabg_w(iwp_globe)%irank

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_w
       if (recv_w_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_w, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_w.

    if (jrecv > nrecvs_w) then
       nrecvs_w = jrecv
       recv_w_iremotes(jrecv) = iremote
    endif

    ! Enter remote and local ranks in recv-remote-rank tables

    iwp_recv_jrecv(iwp_globe) = jrecv
    iwp_recv_iw   (iwp_globe) = iw

  end subroutine recv_table_w


  subroutine recv_table_v(iv, ivp_globe)

    use mem_ijtabs,   only: itabg_v
    use olam_mpi_atm, only: nrecvs_v

    implicit none

    integer, intent(in)  :: iv, ivp_globe
    integer              :: jrecv, iremote

    iremote = itabg_v(ivp_globe)%irank

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_v
       if (recv_v_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_v, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_v.

    if (jrecv > nrecvs_v) then
       nrecvs_v = jrecv
       recv_v_iremotes(jrecv) = iremote
    endif

    ! Enter remote and local ranks in recv-remote-rank tables

    ivp_recv_jrecv(ivp_globe) = jrecv
    ivp_recv_iv   (ivp_globe) = iv

  end subroutine recv_table_v


  subroutine recv_table_m(im, imp_globe)

    use mem_ijtabs,   only: itabg_m
    use olam_mpi_atm, only: nrecvs_m

    implicit none

    integer, intent(in)  :: im, imp_globe
    integer              :: jrecv, iremote

    iremote = itabg_m(imp_globe)%irank

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_m
       if (recv_m_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_m, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_m.

    if (jrecv > nrecvs_m) then
       nrecvs_m = jrecv
       recv_m_iremotes(jrecv) = iremote
    endif

    ! Enter remote and local ranks in recv-remote-rank tables

    imp_recv_jrecv(imp_globe) = jrecv
    imp_recv_im   (imp_globe) = im

  end subroutine recv_table_m


  subroutine send_table_w(iw,iremote)

    use mem_ijtabs,   only: itabg_w
    use olam_mpi_atm, only: nsends_w

    implicit none

    integer,  intent(in) :: iw       ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_w is already in table of ranks to send to

    do jsend = 1, nsends_w
       if (send_w_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_w, jsend represents a rank not yet entered in
    ! the table, so increase nsends_w.

    if (jsend > nsends_w) then
       nsends_w = jsend
       send_w_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table

    iw_sends( itabg_w(iw)%iw_myrank, jsend ) = .true.

  end subroutine send_table_w


  subroutine send_table_v(iv,iremote)

    use mem_ijtabs,   only: itabg_v
    use olam_mpi_atm, only: nsends_v

    implicit none

    integer,  intent(in) :: iv       ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_v is already in table of ranks to send to

    do jsend = 1, nsends_v
       if (send_v_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_v, jsend represents a rank not yet entered in
    ! the table, so increase nsends_v.

    if (jsend > nsends_v) then
       nsends_v = jsend
       send_v_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table

    iv_sends( itabg_v(iv)%iv_myrank, jsend ) = .true.

  end subroutine send_table_v


  subroutine send_table_m(im,iremote)

    use mem_ijtabs,   only: itabg_m
    use olam_mpi_atm, only: nsends_m

    implicit none

    integer,  intent(in) :: im       ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_m is already in table of ranks to send to

    do jsend = 1, nsends_m
       if (send_m_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_m, jsend represents a rank not yet entered in the
    ! table, so increase nsends_m.

    if (jsend > nsends_m) then
       nsends_m = jsend
       send_m_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table

    im_sends( itabg_m(im)%im_myrank, jsend ) = .true.

  end subroutine send_table_m

end subroutine para_init_atm

!===============================================================================

subroutine serial_init_atm()

  use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                        itab_m_vars, itab_v_vars, itab_w_vars, &
                        itabg_m,     itabg_v,     itabg_w,     &
                        itab_m_pd,   itab_v_pd,   itab_w_pd,   &
                        alloc_itabs

  use mem_grid,   only: nma, nva, nwa, mma, mva, mwa, &
                        alloc_gridz, alloc_xyzem, alloc_xyzew, &
                        alloc_grid1, alloc_grid2
  implicit none

  integer :: im, iv, iw

  ! Allocate grid structure variables

  call alloc_gridz()
  call alloc_itabs(mma, mva, mwa, 1)
  call alloc_xyzem(mma)
  call alloc_grid1(mma, mva, mwa)
  call alloc_grid2(mma, mva, mwa)

  ! Store new myrank M, V, and W indices in itabg data structures

  do im = 1, nma
     itabg_m(im)%im_myrank = im
     itab_m (im)%imglobe   = im
     itab_m (im)%irank     = 0
  enddo

  do iv = 1, nva
     itabg_v(iv)%iv_myrank = iv
     itab_v (iv)%ivglobe   = iv
     itab_v (iv)%irank     = 0
  enddo

  do iw = 1, nwa
     itabg_w(iw)%iw_myrank = iw
     itab_w (iw)%iwglobe   = iw
     itab_w (iw)%irank     = 0
  enddo

  ! Read the grid structure for all points in local parallel subdomain
  ! for this rank (or for all points in domain if run is sequential)

  call gridfile_read()

  ! Copy and convert global itab_m_pd index information to itab_m

  do im = 1, mma
     itab_m(im)%npoly   = itab_m_pd(im)%npoly
     itab_m(im)%im(1:3) = itab_m_pd(im)%im(1:3)
     itab_m(im)%iv(1:3) = itab_m_pd(im)%iv(1:3)
     itab_m(im)%iw(1:3) = itab_m_pd(im)%iw(1:3)
     itab_m(im)%imp     = itab_m_pd(im)%imp
  enddo

  ! Copy and convert global itab_v_pd index information to itab_v

  do iv = 1, mva
     itab_v(iv)%im(1:2) = itab_v_pd(iv)%im(1:2)
     itab_v(iv)%iw(1:2) = itab_v_pd(iv)%iw(1:2)
     itab_v(iv)%ivp     = itab_v_pd(iv)%ivp
  enddo

  ! Copy and convert global itab_w_pd index information to itab_w

  do iw = 1, mwa ! global index
     itab_w(iw)%npoly   = itab_w_pd(iw)%npoly
     itab_w(iw)%im(1:7) = itab_w_pd(iw)%im(1:7)
     itab_w(iw)%iv(1:7) = itab_w_pd(iw)%iv(1:7)
     itab_w(iw)%iw(1:7) = itab_w_pd(iw)%iw(1:7)
     itab_w(iw)%iwp     = itab_w_pd(iw)%iwp
  enddo

  ! Deallocate para_decomp _pd arrays

  deallocate( itab_m_pd )
  deallocate( itab_v_pd )
  deallocate( itab_w_pd )

end subroutine serial_init_atm
