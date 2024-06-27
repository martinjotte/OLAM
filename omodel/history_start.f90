subroutine history_start(action)

  ! This routine initializes the model from the history file

  use misc_coms,  only: io6, hfilin, time8, time_istp8
  use hdf5_utils, only: shdf5_exists, shdf5_irec, shdf5_open, shdf5_close
  use max_dims,   only: pathlen
  use mem_para,   only: olam_mpi_finalize

  implicit none

  character(*), intent(in) :: action
  logical                  :: exists

  ! Check if history files exist

  call shdf5_exists(hfilin, exists)

  if (exists) then

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening history file '//trim(hfilin)//' for '//trim(action)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(hfilin, 'R')

     if (action == 'COMMIO') then

        ! Read the common variables
        call commio('READ')

     else

        ! Read the model fields
        call shdf5_irec(1, (/1/), 'time8', dvars=time8)

        write(io6,*) 'calling hist_read'
        call hist_read()

        write(io6,*) 'finished history read'

     endif

     time_istp8 = time8  ! time_istp8 is used for plotting time
     call shdf5_close()

  else

     ! History files do not exist, stop model.

     write(io6,*)
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open history file file:'
     write(io6,*) '!!!   '//trim(hfilin)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) 'in history_start'

     call olam_mpi_finalize()
     stop

  endif

end subroutine history_start

!=========================================================================

subroutine hist_read()

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza, lpv, &
                         vnxo2, vnyo2, vnzo2, lpw
  use mem_ijtabs,  only: itab_w, itab_v, itab_m, jtab_v, jtv_prog
  use misc_coms,   only: io6, runtype, iparallel
  use mem_basic,   only: vc, vmc, vmp, wc, wmc, rho, vxe, vye, vze
  use var_tables,  only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils,  only: shdf5_info, shdf5_irec
  use mem_sfcg,    only: nwsfc, nvsfc, nmsfc, mwsfc, mvsfc, mmsfc, &
                         itab_wsfc, itab_vsfc, itab_msfc
  use mem_land,    only: itab_land, nland, mland
  use mem_lake,    only: itab_lake, nlake, mlake
  use mem_sea,     only: itab_sea, nsea, msea
  use mem_nudge,   only: nwnud, mwnud, itab_wnud
  use consts_coms, only: r8
  use mem_regrid, only: interp_regrid
  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

  implicit none

  integer       :: j, iw, iw1, iw2, iv, kb, k
  integer       :: nv, nvcnt, ns, ndims, idims(3), jdims(3)
  character(32) :: varn
  character (2) :: stagpt
  integer       :: ilocal(max(mwa,mva,mma,mwsfc,mvsfc,mland,mlake,msea,mwnud))

! Scratch arrays for copying input

  integer,  allocatable :: iscr1(:)
  real,     allocatable :: rscr1(:)
  real,     allocatable :: rscr2(:,:)
  real(r8), allocatable :: dscr1(:)
  real(r8), allocatable :: dscr2(:,:)

  nvcnt =  0

  ! Loop through all variables in the model vtables

  do nv = 1, num_var

     ndims = -1
     idims =  0

     ! Skip to next variable if we don't want the current one

     if (.not. vtab_r(nv)%ihist) cycle
     if (vtab_r(nv)%nread) cycle

     varn    = trim(vtab_r(nv)%name)
     stagpt  = vtab_r(nv)%stagpt

     call shdf5_info(varn, ndims, idims)

     if (runtype == 'HISTREGRID') then

        ! In a HISTREGRID restart, the following need not or should not be read
        ! from the OLD grid history file.  Most are diagnostic or auxiliary
        ! variables and are initialized or computed elsewhere, no later than on
        ! the first timestep of the history restart.  The remaining ones are
        ! horizontal velocity quantities defined at IV stagger points and cannot
        ! be interpolated.  Instead, they are re-initialized by projecting
        ! interpolated VXE, VYE, VZE components onto grid cell faces at IV points.

        if (trim(varn) == 'VMC'              .or. &
            trim(varn) == 'VC'               .or. &
            trim(varn) == 'VMP'              .or. &
            trim(varn) == 'VC_ACCUM'         .or. &
            trim(varn) == 'THSRC'            .or. &
            trim(varn) == 'RTSRC'            .or. &
            trim(varn) == 'UMSRC'            .or. &
            trim(varn) == 'VMSRC'            .or. &
            trim(varn) == 'QWCON'            .or. &
            trim(varn) == 'CBMF'             .or. &
            trim(varn) == 'CDDF'             .or. &
            trim(varn) == 'KCUTOP'           .or. &
            trim(varn) == 'KCUBOT'           .or. &
            trim(varn) == 'KUDBOT'           .or. &
            trim(varn) == 'KDDTOP'           .or. &
            trim(varn) == 'KDDMAX'           .or. &
            trim(varn) == 'KDDBOT'           .or. &
            trim(varn) == 'KSTABI'           .or. &
            trim(varn) == 'CU_PWA'           .or. &
            trim(varn) == 'CU_PEV'           .or. &
            trim(varn) == 'IACTCU'           .or. &
            trim(varn) == 'FTHRD_SW'         .or. &
            trim(varn) == 'FTHRD_LW'         .or. &
            trim(varn) == 'CLOUD_FRAC'       .or. &
            trim(varn) == 'RSHORT'           .or. &
            trim(varn) == 'RLONG'            .or. &
            trim(varn) == 'RLONGUP'          .or. &
            trim(varn) == 'RSHORT_TOP'       .or. &
            trim(varn) == 'RSHORTUP_TOP'     .or. &
            trim(varn) == 'RLONGUP_TOP'      .or. &
            trim(varn) == 'ALBEDT'           .or. &
            trim(varn) == 'COSZ'             .or. &
            trim(varn) == 'RSHORT_CLR'       .or. &
            trim(varn) == 'RSHORTUP_CLR'     .or. &
            trim(varn) == 'RSHORT_TOP_CLR'   .or. &
            trim(varn) == 'RSHORTUP_TOP_CLR' .or. &
            trim(varn) == 'RLONG_CLR'        .or. &
            trim(varn) == 'RLONGUP_CLR'      .or. &
            trim(varn) == 'RLONGUP_TOP_CLR'  .or. &
            trim(varn) == 'PBL_CLD_FORC'     .or. &
            trim(varn) == 'PAR_TOT'          .or. &
            trim(varn) == 'PAR_DIFFUSE'      .or. &
            trim(varn) == 'PPFD_TOT'         .or. &
            trim(varn) == 'PPFD_DIFFUSE'     .or. &
            trim(varn) == 'MCICA_SEED'       .or. &
            trim(varn) == 'VKM'              .or. &
            trim(varn) == 'VKH'              .or. &
            trim(varn) == 'VKM_SFC'          .or. &
            trim(varn) == 'SFLUXT'           .or. &
            trim(varn) == 'SFLUXR'           .or. &
            trim(varn) == 'USTAR'            .or. &
            trim(varn) == 'WSTAR'            .or. &
            trim(varn) == 'WTV0')  cycle

     endif

     ! We want it...read it if it's in the history file

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        if (varn == 'VMP') then

           write(io6,*)
           write(io6,*) 'Variable VMP is not in the history file.'
           write(io6,*) 'Setting VMP to VMC.'
           write(io6,*)
           vmp = vmc
           cycle

        else

           write(io6,*)
           write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
           write(io6,*)
           cycle

        endif
     endif

     ! Identify the points we want to read from the history file

     if     (stagpt == 'AW') then
        if (runtype /= 'HISTREGRID' .and. idims(ndims) == nwa) then
           ilocal(1:mwa) = itab_w(1:mwa)%iwglobe
           idims(ndims) = mwa
        elseif (runtype == 'HISTREGRID') then
           jdims = idims
           jdims(ndims) = nwa
        else
           stop "invalid AW array size in history_read"
        endif
     elseif (stagpt == 'AV' .and. idims(ndims) == nva) then
        ilocal(1:mva) = itab_v(1:mva)%ivglobe
        idims(ndims) = mva
     elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
        ilocal(1:mma) = itab_m(1:mma)%imglobe
        idims(ndims) = mma
     elseif (stagpt == 'CW' .and. idims(ndims) == nwsfc) then
        ilocal(1:mwsfc) = itab_wsfc(1:mwsfc)%iwglobe
        idims(ndims) = mwsfc
     elseif (stagpt == 'CV' .and. idims(ndims) == nvsfc) then
        ilocal(1:mvsfc) = itab_vsfc(1:mvsfc)%ivglobe
        idims(ndims) = mvsfc
     elseif (stagpt == 'CM' .and. idims(ndims) == nmsfc) then
        ilocal(1:mmsfc) = itab_msfc(1:mmsfc)%imglobe
        idims(ndims) = mmsfc
     elseif (stagpt == 'LW' .and. idims(ndims) == nland) then
        ilocal(1:mland) = itab_land(1:mland)%iwglobe
        idims(ndims) = mland
     elseif (stagpt == 'RW' .and. idims(ndims) == nlake) then
        ilocal(1:mlake) = itab_lake(1:mlake)%iwglobe
        idims(ndims) = mlake
     elseif (stagpt == 'SW' .and. idims(ndims) == nsea) then
        ilocal(1:msea) = itab_sea(1:msea)%iwglobe
        idims(ndims) = msea
     elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
        ilocal(1:mwnud) = itab_wnud(1:mwnud)%iwnudglobe
        idims(ndims) = mwnud
     else
        stop "invalid array size in history_read"
     endif

     ns = idims(ndims)

     if     (associated(vtab_r(nv)%ivar1_p)) then

        if (stagpt == 'AW' .and. runtype == 'HISTREGRID') then
           allocate(iscr1(idims(1)))
           call shdf5_irec(ndims, idims, varn, ivar1=iscr1)
           call interp_regrid(ndims, idims, jdims, varn, stagpt, &
                               ivara1=iscr1, ivarb1=vtab_r(nv)%ivar1_p)
           deallocate(iscr1)
        else
           call shdf5_irec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, &
                           points=ilocal(1:ns))
        endif

     elseif (associated(vtab_r(nv)%ivar2_p)) then

        call shdf5_irec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, &
                        points=ilocal(1:ns))

     elseif (associated(vtab_r(nv)%ivar3_p)) then

        call shdf5_irec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, &
                        points=ilocal(1:ns))


     elseif (associated(vtab_r(nv)%rvar1_p)) then

        if (stagpt == 'AW' .and. runtype == 'HISTREGRID') then
           allocate(rscr1(idims(1)))
           call shdf5_irec(ndims, idims, varn, rvar1=rscr1)
           call interp_regrid(ndims, idims, jdims, varn, stagpt, &
                               rvara1=rscr1, rvarb1=vtab_r(nv)%rvar1_p)
           deallocate(rscr1)
        else
           call shdf5_irec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, &
                           points=ilocal(1:ns))
        endif

     elseif (associated(vtab_r(nv)%rvar2_p)) then

        if (stagpt == 'AW' .and. runtype == 'HISTREGRID') then
           allocate(rscr2(idims(1),idims(2)))
           call shdf5_irec(ndims, idims, varn, rvar2=rscr2)
           call interp_regrid(ndims, idims, jdims, varn, stagpt, &
                               rvara2=rscr2, rvarb2=vtab_r(nv)%rvar2_p)
           deallocate(rscr2)
        else
           call shdf5_irec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, &
                           points=ilocal(1:ns))
        endif

     elseif (associated(vtab_r(nv)%rvar3_p)) then

        call shdf5_irec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, &
                        points=ilocal(1:ns))

     elseif (associated(vtab_r(nv)%dvar1_p)) then

        if (stagpt == 'AW' .and. runtype == 'HISTREGRID') then
           allocate(dscr1(idims(1)))
           call shdf5_irec(ndims, idims, varn, dvar1=dscr1)
           call interp_regrid(ndims, idims, jdims, varn, stagpt, &
                               dvara1=dscr1, dvarb1=vtab_r(nv)%dvar1_p)
           deallocate(dscr1)
        else
           call shdf5_irec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, &
                           points=ilocal(1:ns))
        endif

     elseif (associated(vtab_r(nv)%dvar2_p)) then

        if (stagpt == 'AW' .and. runtype == 'HISTREGRID') then
           allocate(dscr2(idims(1),idims(2)))
           call shdf5_irec(ndims, idims, varn, dvar2=dscr2)
           call interp_regrid(ndims, idims, jdims, varn, stagpt, &
                               dvara2=dscr2, dvarb2=vtab_r(nv)%dvar2_p)
           deallocate(dscr2)
        else
           call shdf5_irec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, &
                           points=ilocal(1:ns))
        endif

     elseif (associated(vtab_r(nv)%dvar3_p)) then

        call shdf5_irec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, &
                        points=ilocal(1:ns))
     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

  ! For HISTREGRID history start, diagnose VC and VMC from VXE, VYE, VZE

  if (runtype == 'HISTREGRID') then

     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

        ! Compute VC by projecting vxe, vye, vze onto IV face

        kb = lpv(iv)

        do k = kb,mza
           vc(k,iv) = vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                    + vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                    + vnzo2(iv) * (vze(k,iw1) + vze(k,iw2))

           vmc(k,iv) = vc(k,iv) * .5 * real(rho(k,iw1) + rho(k,iw2))
        enddo

        vc (1:kb-1,iv) = 0.
        vmc(1:kb-1,iv) = 0.
     enddo

     ! MPI parallel send/recv of V group

     if (iparallel == 1) then
        call mpi_send_v(rvara1=vmc, rvara2=vc)
        call mpi_recv_v(rvara1=vmc, rvara2=vc)
     endif

     ! Now set underground W momentum/velocity to 0

     do iw = 2, mwa
        wmc(1:lpw(iw)-1,iw) = 0.
        wc (1:lpw(iw)-1,iw) = 0.
     enddo

  endif

end subroutine hist_read

