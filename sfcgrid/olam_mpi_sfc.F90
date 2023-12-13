Module olam_mpi_sfc

  integer, parameter :: itagwsfc = 11
  integer, parameter :: itagvsfc = 12
  integer, parameter :: itagmsfc = 13

  integer :: nsends_wsfc, nrecvs_wsfc
  integer :: nsends_vsfc, nrecvs_vsfc
  integer :: nsends_msfc, nrecvs_msfc

  Type nodebuffs
     character, allocatable :: buff(:)
     integer,   allocatable :: ipts(:)
     integer                :: nbytes = 0
     integer                :: jend = 0
     integer                :: iremote = -1
  End Type nodebuffs

  type(nodebuffs), allocatable :: send_wsfc(:), recv_wsfc(:)
  type(nodebuffs), allocatable :: send_vsfc(:), recv_vsfc(:)
  type(nodebuffs), allocatable :: send_msfc(:), recv_msfc(:)

  integer              :: icurr_wsfc = 1
  integer              :: inext_wsfc = 2
  integer, allocatable :: ireqr_wsfc(:,:)
  integer, allocatable :: ireqs_wsfc(:,:)

  integer              :: icurr_vsfc = 1
  integer              :: inext_vsfc = 2
  integer, allocatable :: ireqr_vsfc(:,:)
  integer, allocatable :: ireqs_vsfc(:,:)

  integer              :: icurr_msfc = 1
  integer              :: inext_msfc = 2
  integer, allocatable :: ireqr_msfc(:,:)
  integer, allocatable :: ireqs_msfc(:,:)

Contains

subroutine olam_mpi_sfc_start()

#ifdef OLAM_MPI
  use mem_para,  only: nbytes_int, nbytes_real
  use mem_land,  only: nzg
  use leaf_coms, only: nzs
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI
  integer :: ierr, jrecv, jsend
  integer :: nbytes_per_iwsfc
  integer :: nbytes_per_ivsfc
  integer :: nbytes_per_imsfc

  allocate( ireqs_wsfc(nsends_wsfc,2) ) ; ireqs_wsfc = MPI_REQUEST_NULL
  allocate( ireqr_wsfc(nrecvs_wsfc,2) ) ; ireqr_wsfc = MPI_REQUEST_NULL

  allocate( ireqs_vsfc(nsends_vsfc,2) ) ; ireqs_vsfc = MPI_REQUEST_NULL
  allocate( ireqr_vsfc(nrecvs_vsfc,2) ) ; ireqr_vsfc = MPI_REQUEST_NULL

  allocate( ireqs_msfc(nsends_msfc,2) ) ; ireqs_msfc = MPI_REQUEST_NULL
  allocate( ireqr_msfc(nrecvs_msfc,2) ) ; ireqr_msfc = MPI_REQUEST_NULL

  ! Determine number of bytes to send per IWSFC column

  nbytes_per_iwsfc =  3       * nbytes_int  &
                   + 20       * nbytes_real &
                   +  2 * nzg * nbytes_real &
                   +  3 * nzs * nbytes_real

  do jrecv = 1, nrecvs_wsfc
     recv_wsfc(jrecv)%nbytes  = recv_wsfc(jrecv)%jend * nbytes_per_iwsfc
     allocate( recv_wsfc(jrecv)%buff( recv_wsfc(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_wsfc
     send_wsfc(jsend)%nbytes  = send_wsfc(jsend)%jend * nbytes_per_iwsfc
     allocate( send_wsfc(jsend)%buff( send_wsfc(jsend)%nbytes ) )
  enddo

  ! Determine number of bytes to send per IVSFC column

  nbytes_per_ivsfc = 8       * nbytes_real &
                   + 2 * nzg * nbytes_real

  do jrecv = 1, nrecvs_vsfc
     recv_vsfc(jrecv)%nbytes  = recv_vsfc(jrecv)%jend * nbytes_per_ivsfc
     allocate( recv_vsfc(jrecv)%buff( recv_vsfc(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_vsfc
     send_vsfc(jsend)%nbytes  = send_vsfc(jsend)%jend * nbytes_per_ivsfc
     allocate( send_vsfc(jsend)%buff( send_vsfc(jsend)%nbytes ) )
  enddo

  ! Determine number of bytes to send per IMSFC column

  nbytes_per_imsfc = 1 * nbytes_real

  do jrecv = 1, nrecvs_msfc
     recv_msfc(jrecv)%nbytes  = recv_msfc(jrecv)%jend * nbytes_per_imsfc
     allocate( recv_msfc(jrecv)%buff( recv_msfc(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_msfc
     send_msfc(jsend)%nbytes = send_msfc(jsend)%jend * nbytes_per_imsfc
     allocate( send_msfc(jsend)%buff( send_msfc(jsend)%nbytes ) )
  enddo

  ! Pre-post non-blocking receives

  do jrecv = 1, nrecvs_wsfc
     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr                          )
  enddo

  do jrecv = 1, nrecvs_vsfc
     call MPI_Irecv(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_vsfc(jrecv)%iremote, itagvsfc, MPI_COMM_WORLD,         &
                    ireqr_vsfc(jrecv,inext_vsfc), ierr                          )
  enddo

  do jrecv = 1, nrecvs_msfc
     call MPI_Irecv(recv_msfc(jrecv)%buff, recv_msfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_msfc(jrecv)%iremote, itagmsfc, MPI_COMM_WORLD,         &
                    ireqr_msfc(jrecv,inext_msfc), ierr                          )
  enddo

#endif

end subroutine olam_mpi_sfc_start

!===============================================================================

subroutine mpi_send_wsfc(set, soil_watfrac, div2d_ex)

  ! Subroutine to perform a parallel MPI send of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,  only: sfcg
  use mem_land,  only: nzg, land, mland, omland
  use mem_lake,  only: lake, omlake
  use mem_sea,   only: sea, msea, omsea
  use leaf_coms, only: nzs
  use sea_coms,  only: nzi
  use misc_coms, only: do_chem

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: soil_watfrac(nzg,mland)
  real, optional, intent(inout) :: div2d_ex(msea)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp, jj, nb
  integer :: iwsfc, iland, ilake, isea

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_wsfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_wsfc, ireqs_wsfc(:,icurr_wsfc), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,jj,iwsfc,nb,iland,ilake,isea,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0
     nb   = send_wsfc(jsend)%nbytes

     ! Loop over number of columns for this jsend

     do jj = 1, send_wsfc(jsend)%jend
        iwsfc = send_wsfc(jsend)%ipts(jj)

        ! Pack the messages into send buffers

        if (set == 'head') then

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(land%head   (1,iland),nzg,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(soil_watfrac(1,iland),nzg,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_grad') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%gxps_vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gxps_vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gxps_vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_progw') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%vmxet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmyet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmzet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmxet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmyet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmzet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%swmdepth  (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%head1   (iwsfc),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_div2d_ex') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(div2d_ex(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           endif

        elseif (set == 'swm_diagvel') then  ! [call from swm_driver after swm_diagvel]

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then

              isea = iwsfc - omsea

              call MPI_Pack(sea%vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           endif

        elseif (set == 'sfc_driv_end') then

           call MPI_Pack(sfcg%cantemp(iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%canrrv (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%wthv   (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rough  (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%nlev_seaice   (isea),  1,MPI_INTEGER, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_cantemp   (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_cantemp   (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_canrrv    (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_canrrv    (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_wthv      (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_wthv      (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_rough     (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_rough     (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaicec       (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seatc         (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaice_tempk(1,isea),nzi,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              call MPI_Pack(lake%lake_energy(ilake),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(lake%depth      (ilake),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(land%sfcwater_mass  (1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_energy(1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_depth (1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%nlev_sfcwater    (iland),1,MPI_INTEGER, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_fracarea     (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_albedo       (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_rough        (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_temp         (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_water     (1,iland),nzg, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_energy    (1,iland),nzg, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

              if (do_chem ==1) then
              call MPI_Pack(land%ppfd             (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd_diffuse     (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              endif

              call MPI_Pack(land%cosz             (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'avgatm') then

           call MPI_Pack(sfcg%vels    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%prss    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rhos    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtemp (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtheta(iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airrrv  (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%windxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%windye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%windze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_wsfc(jsend)%buff, ipos, MPI_PACKED,            &
                    send_wsfc(jsend)%iremote, itagwsfc, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,inext_wsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_wsfc = mod(icurr_wsfc,2) + 1
  inext_wsfc = mod(inext_wsfc,2) + 1

#endif

end subroutine mpi_send_wsfc

!===============================================================================

subroutine mpi_send_vsfc(set, watflux, energyflux, vc_ex)

  ! Subroutine to perform a parallel MPI send of a "VSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,   only: sfcg, itab_vsfc, mvsfc
  use mem_land,   only: nzg
  use oname_coms, only: nl

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: watflux   (nzg,mvsfc)
  real, optional, intent(inout) :: energyflux(nzg,mvsfc)
  real, optional, intent(inout) :: vc_ex         (mvsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp, jj, nb
  integer :: ivsfc, iw1, iw2

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_vsfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_vsfc, ireqs_vsfc(:,icurr_vsfc), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,jj,ivsfc,nb,iw1,iw2,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0
     nb = send_vsfc(jsend)%nbytes

     ! Loop over number of VSFC points for this jsend

     do jj = 1, send_vsfc(jsend)%jend
        ivsfc = send_vsfc(jsend)%ipts(jj)

        iw1 = itab_vsfc(ivsfc)%iwn(1)
        iw2 = itab_vsfc(ivsfc)%iwn(2)

        if (present(watflux)) then  ! follows ivsfc horiz flux loop in surface_driver

           call MPI_Pack(watflux   (1,ivsfc),nzg,MPI_REAL,send_vsfc(jsend)%buff, &
                                                 nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(energyflux(1,ivsfc),nzg,MPI_REAL,send_vsfc(jsend)%buff, &
                                                 nb,ipos,MPI_COMM_WORLD,ierr)

        elseif (set == 'swm_hflux') then

           ! Updates from swm_hflux; limited to swm_active points

           if ((sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) .and. &
                nl%igw_spinup /= 1) then

              call MPI_Pack(sfcg%hflux_wat(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_enr(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vxe(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vye(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vze(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vmp      (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vmc      (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vc       (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (present(vc_ex)) then ! follows vc_ex computation in swm_progv

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Pack(vc_ex(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                           nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        else  ! follows swm_progv call in swm_driver

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Pack(sfcg%vmc(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                              nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vc (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                              nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_vsfc(jsend)%buff, ipos, MPI_PACKED,            &
                    send_vsfc(jsend)%iremote, itagvsfc, MPI_COMM_WORLD, &
                    ireqs_vsfc(jsend,inext_vsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_vsfc = mod(icurr_vsfc,2) + 1
  inext_vsfc = mod(inext_vsfc,2) + 1

#endif

end subroutine mpi_send_vsfc

!===============================================================================

subroutine mpi_send_msfc(vort)

  ! Subroutine to perform a parallel MPI send of a "MSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,   only: sfcg, itab_msfc, mmsfc

  implicit none

  real, optional, intent(inout) :: vort(mmsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp, jj, nb
  integer :: imsfc, iw1, iw2

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_msfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_msfc, ireqs_msfc(:,icurr_msfc), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,nb,jj,imsfc,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0
     nb = send_msfc(jsend)%nbytes

     ! Loop over number of MSFC points for this jsend

     do jj = 1, send_msfc(jsend)%jend
        imsfc = send_msfc(jsend)%ipts(jj)

        if (present(vort)) then

           call MPI_Pack(vort(imsfc),1,MPI_REAL,send_msfc(jsend)%buff, &
                                       nb,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_msfc(jsend)%buff, ipos, MPI_PACKED,            &
                    send_msfc(jsend)%iremote, itagmsfc, MPI_COMM_WORLD, &
                    ireqs_msfc(jsend,inext_msfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_msfc = mod(icurr_msfc,2) + 1
  inext_msfc = mod(inext_msfc,2) + 1

#endif

end subroutine mpi_send_msfc

!=============================================================================

subroutine mpi_recv_wsfc(set, soil_watfrac, div2d_ex)

  ! Subroutine to perform a parallel MPI receive of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,  only: sfcg
  use leaf_coms, only: nzs
  use sea_coms,  only: nzi
  use mem_land,  only: nzg, mland, land, omland
  use mem_lake,  only: lake, omlake
  use mem_sea,   only: sea, msea, omsea
  use misc_coms, only: do_chem

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: soil_watfrac(nzg,mland)
  real, optional, intent(inout) :: div2d_ex(msea)

#ifdef OLAM_MPI

  integer :: ierr, ipos, count
  integer :: jrecv, jtmp, jj
  integer :: iwsfc, iland, ilake, isea
  integer :: status(mpi_status_size)

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_wsfc

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, status, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(ipos,jj,iwsfc,iland,ilake,isea,ierr,count) &
     !$omp      firstprivate(jrecv,status) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

     call MPI_Get_count(status, MPI_PACKED, count, ierr)

     if (count /= 0) then
     do jj = 1, recv_wsfc(jrecv)%jend
        iwsfc = recv_wsfc(jrecv)%ipts(jj)

        if (set == 'head') then

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%head   (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              soil_watfrac(1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_grad') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_progw') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmxet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmyet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmzet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmxet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmyet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmzet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%swmdepth  (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%head1   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_div2d_ex') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              div2d_ex     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_diagvel') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vxe    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vye    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vze    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'sfc_driv_end') then

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%cantemp  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%canrrv   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%wthv      (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%rough     (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%nlev_seaice (isea),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_cantemp (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_cantemp (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_canrrv  (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_canrrv  (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_wthv    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_wthv    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_rough   (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_rough   (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaicec     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seatc       (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaice_tempk(1,isea),nzi,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              lake%lake_energy(ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              lake%depth      (ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_mass  (1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_energy(1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_depth (1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%nlev_sfcwater    (iland),  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_fracarea     (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_albedo       (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_rough        (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_temp         (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_water     (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_energy    (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

              if (do_chem ==1) then
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd             (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd_diffuse     (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              endif

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%cosz             (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'avgatm') then

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%vels    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%prss    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%rhos    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%airtemp (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%airtheta(iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%airrrv  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo
     endif

     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_wsfc

!=============================================================================

subroutine mpi_recv_vsfc(set, watflux, energyflux, vc_ex)

  ! Subroutine to perform a parallel MPI receive of a "VSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,   only: sfcg, mvsfc, itab_vsfc
  use mem_land,   only: nzg
  use oname_coms, only: nl

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: watflux   (nzg,mvsfc)
  real, optional, intent(inout) :: energyflux(nzg,mvsfc)
  real, optional, intent(inout) :: vc_ex         (mvsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos, count
  integer :: jrecv, jtmp
  integer :: jj, ivsfc, iw1, iw2
  integer :: status(mpi_status_size)

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_vsfc

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_vsfc, ireqr_vsfc(:,icurr_vsfc), jrecv, status, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(ipos,jj,ivsfc,iw1,iw2,ierr,count) &
     !$omp      firstprivate(jrecv,status) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

     call MPI_Get_count(status, MPI_PACKED, count, ierr)

     if (count /= 0) then
     do jj = 1, recv_vsfc(jrecv)%jend
        ivsfc = recv_vsfc(jrecv)%ipts(jj)

        iw1 = itab_vsfc(ivsfc)%iwn(1)
        iw2 = itab_vsfc(ivsfc)%iwn(2)

        if (present(watflux)) then  ! follows ivsfc horiz flux loop in surface_driver

           call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                           watflux   (1,ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                           energyflux(1,ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

        elseif (set == 'swm_hflux') then

           ! Updates from swm_hflux; limited to swm_active points

           if ((sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) .and. &
                nl%igw_spinup /= 1) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_wat(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_enr(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vxe(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vye(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vze(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmp      (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmc      (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vc       (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (present(vc_ex)) then ! follows vc_ex computation in swm_progv

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              vc_ex(ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        else  ! follows swm_progv call in swm_driver

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmc(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vc (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo
     endif

     call MPI_Irecv(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_vsfc(jrecv)%iremote, itagvsfc, MPI_COMM_WORLD,         &
                    ireqr_vsfc(jrecv,inext_vsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_vsfc

!=============================================================================

subroutine mpi_recv_msfc(vort)

  ! Subroutine to perform a parallel MPI receive of a "VSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,   only: sfcg, mmsfc, itab_msfc
  use mem_land,   only: nzg
  use oname_coms, only: nl

  implicit none

  real, optional, intent(inout) :: vort(mmsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos, count
  integer :: jrecv, jtmp
  integer :: jj, imsfc
  integer :: status(mpi_status_size)

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_msfc

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_msfc, ireqr_msfc(:,icurr_msfc), jrecv, status, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(ipos,jj,imsfc,ierr,count) &
     !$omp      firstprivate(jrecv,status) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

     call MPI_Get_count(status, MPI_PACKED, count, ierr)

     if (count /= 0) then
     do jj = 1, recv_msfc(jrecv)%jend
        imsfc = recv_msfc(jrecv)%ipts(jj)

        if (present(vort)) then

           call MPI_Unpack(recv_msfc(jrecv)%buff,recv_msfc(jrecv)%nbytes,ipos, &
                           vort(imsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

     enddo
     endif

     call MPI_Irecv(recv_msfc(jrecv)%buff, recv_msfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_msfc(jrecv)%iremote, itagmsfc, MPI_COMM_WORLD,         &
                    ireqr_msfc(jrecv,inext_msfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_msfc

!===============================================================================

subroutine olam_mpi_sfc_stop()

#ifdef OLAM_MPI
  use mpi
#endif

  use misc_coms, only: iparallel
  use mem_para,  only: olam_mpi_barrier

  implicit none

#ifdef OLAM_MPI

  integer :: ierr, jrecv, ii, jsend
  logical :: flags

  if (iparallel == 1) then

     ! First make sure all processes are here so that all sends are finished

     call olam_mpi_barrier()

     ! Cancel any pending communication

     do ii = 1, 2

        do jrecv = 1, nrecvs_wsfc
           if (ireqr_wsfc(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_wsfc(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_wsfc
           if (ireqs_wsfc(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_wsfc(jsend,ii), ierr)
           endif
        enddo

        do jrecv = 1, nrecvs_vsfc
           if (ireqr_vsfc(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_vsfc(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_vsfc
           if (ireqs_vsfc(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_vsfc(jsend,ii), ierr)
           endif
        enddo

        do jrecv = 1, nrecvs_msfc
           if (ireqr_msfc(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_msfc(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_msfc
           if (ireqs_msfc(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_msfc(jsend,ii), ierr)
           endif
        enddo

     enddo

     ! Test that all communication requests have been completed or cancelled

     call olam_mpi_barrier()

     call MPI_Testall(nrecvs_wsfc, ireqr_wsfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_wsfc, ireqr_wsfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_wsfc, ireqs_wsfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_wsfc, ireqs_wsfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nrecvs_vsfc, ireqr_vsfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_vsfc, ireqr_vsfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_vsfc, ireqs_vsfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_vsfc, ireqs_vsfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nrecvs_msfc, ireqr_msfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_msfc, ireqr_msfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_msfc, ireqs_msfc(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_msfc, ireqs_msfc(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call olam_mpi_barrier()

  endif

  ! Deallocate unused arrays

  if (allocated(send_wsfc)) deallocate(send_wsfc)
  if (allocated(recv_wsfc)) deallocate(recv_wsfc)

  if (allocated(send_vsfc)) deallocate(send_vsfc)
  if (allocated(recv_vsfc)) deallocate(recv_vsfc)

  if (allocated(send_msfc)) deallocate(send_msfc)
  if (allocated(recv_msfc)) deallocate(recv_msfc)

  if (allocated(ireqs_wsfc)) deallocate(ireqs_wsfc)
  if (allocated(ireqr_wsfc)) deallocate(ireqr_wsfc)

  if (allocated(ireqs_vsfc)) deallocate(ireqs_vsfc)
  if (allocated(ireqr_vsfc)) deallocate(ireqr_vsfc)

  if (allocated(ireqs_msfc)) deallocate(ireqs_msfc)
  if (allocated(ireqr_msfc)) deallocate(ireqr_msfc)

#endif

end subroutine olam_mpi_sfc_stop

!==================================================================

End Module olam_mpi_sfc
