module clouds_gno

  real :: rh_start
  real :: rh_end
  real :: rh_delta

  real :: qc_start
  real :: qc_end
  real :: qc_delta

  integer :: n_rh
  integer :: n_qc

  real, allocatable :: cldf_tab(:,:)


contains


  subroutine cu_cldfrac(ks, ke, rhin, qcin, cldfrac)

    use mem_grid, only: mza
    implicit none

!===============================================================================
!
! Purpose:
! --------
!
! Parameterization of the cloudiness associated with cumulus convection.
!
! This routine reads in a lookup table of parameterized values of cloud fraction
! as a function of resolved relative humidity (defined as total resolved water
! specific humidity devided by saturation specific humidity) and the square root
! of subgrid convective cloud water normalized by the total resolved water in
! the grid box. The routine to create the lookup table is in the etc/ directory
! in the bony_emanuel_table subdirectory, and uses the routine described in
! Bony and Emanuel (JAS, 2001, pp. 3158 - 3183).
!
!===============================================================================

    integer, intent(in)  :: ks, ke
    real,    intent(in)  :: rhin(mza)  ! q_total / q_sat
    real,    intent(in)  :: qcin(mza)  ! sqrt( qc_cu / q_total )

    real,    intent(out) :: cldfrac(mza)

    integer :: k
    real    :: rh, qc

    real    :: rscale, rint1, rint0
    integer :: rindex

    real    :: qscale, qint1, qint0
    integer :: qindex

    cldfrac(   1:ks-1) = 0.0
    cldfrac(ke+1:mza ) = 0.0

    do k = ks, ke

       rh = rhin(k)
       rh = max( rh_start, min(rh_end, rh) )

       qc = qcin(k)
       qc = max( qc_start, min(qc_end, qc) )

       rscale = (rh - rh_start) / rh_delta
       rindex = max(1, min(n_rh-1, int(rscale) + 1))
       rint1  = rscale - real(rindex-1)
       rint0  = 1.0 - rint1

       qscale = (qc - qc_start) / qc_delta
       qindex  = max(1, min(n_qc-1, int(qscale) + 1))
       qint1  = qscale - real(qindex-1)
       qint0  = 1.0 - qint1

       cldfrac(k) = qint0 * ( rint0 * cldf_tab(qindex  , rindex  )   &
                            + rint1 * cldf_tab(qindex  , rindex+1) ) &
                  + qint1 * ( rint0 * cldf_tab(qindex+1, rindex  )   &
                            + rint1 * cldf_tab(qindex+1, rindex+1) )

    enddo

  end subroutine cu_cldfrac


  subroutine gno_lookup_init()

    use hdf5_utils, only: shdf5_exists, shdf5_open, shdf5_close, shdf5_irec
    use max_dims,   only: pathlen
    use oname_coms, only: nl
    use mem_para,   only: olam_mpi_finalize

    implicit none

    logical :: exists

    character(20),     parameter :: gno_file = "bony_eman_table.h5"
    character(pathlen)           :: inputfile

    integer :: ndims, idims(3)

    inputfile =  trim(nl%rrtmg_datadir) // "/bony_emanuel_table/" &
              // trim(gno_file)

    call shdf5_exists(inputfile, exists)

    if (.not. exists) then
       write(*,*) "gno_clouds_init: Error opening data file " // trim(inputfile)
       write(*,*) "Bony-Emanuel look table datafile cannot be found."
       write(*,*) "Stopping model"
       call olam_mpi_finalize()
       stop
    endif

    call shdf5_open(inputfile, 'R')

    ndims=1 ; idims(1)=1
    call shdf5_irec(ndims, idims, "rh_start", rvars=rh_start)
    call shdf5_irec(ndims, idims, "rh_end"  , rvars=rh_end)
    call shdf5_irec(ndims, idims, "rh_size" , ivars=n_rh)

    rh_delta = (rh_end - rh_start) / real(n_rh - 1)

    call shdf5_irec(ndims, idims, "qc_start", rvars=qc_start)
    call shdf5_irec(ndims, idims, "qc_end"  , rvars=qc_end)
    call shdf5_irec(ndims, idims, "qc_size" , ivars=n_qc)

    qc_delta = (qc_end - qc_start) / real(n_qc - 1)

    allocate(cldf_tab(n_qc, n_rh))

    ndims=2 ; idims(1:2) = (/ n_qc, n_rh /)
    call shdf5_irec(ndims, idims, "cldfrac" , rvar2=cldf_tab)

    call shdf5_close()

  end subroutine gno_lookup_init

end module clouds_gno
