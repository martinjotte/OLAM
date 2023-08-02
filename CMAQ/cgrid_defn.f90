module cgrid_defn

  real, allocatable :: cgrid(:,:,:)

  real, allocatable :: gc_tend(:,:,:)
  real, allocatable :: ae_tend(:,:,:)
  real, allocatable :: nr_tend(:,:,:)

  integer :: ns_co  = 0
  integer :: ns_no  = 0
  integer :: ns_no2 = 0
  integer :: ns_o3  = 0

  integer, allocatable :: gc_spc_to_trns_map(:)
  integer, allocatable :: ae_spc_to_trns_map(:)
  integer, allocatable :: nr_spc_to_trns_map(:)

  real, allocatable :: vdemis_gc(:,:,:)
  real, allocatable :: vdemis_ae(:,:,:)
  real, allocatable :: vdemis_nr(:,:,:)

  character(16), allocatable :: cgrid_names(:)

contains

  subroutine alloc_cgrid(mza, mwa)

    use mem_grid,   only: nsw_max
    use cgrid_spcs, only: n_gc_trns, n_ae_trns, n_nr_trns, &
                          n_gc_emis, n_ae_emis, n_nr_emis, &
                          n_gc_depv, n_ae_depv, n_nr_depv, &
                          n_gc_spc,  n_ae_spc,  n_nr_spc,  &
                          gc_trns_map, ae_trns_map, nr_trns_map, &
                          nspcsd
    use misc_coms,  only: rinit
    use oname_coms, only: nl

    implicit none

    integer, intent(in) :: mza, mwa

    real,    parameter  :: cmin = 1.e-30
    integer             :: n

    allocate( cgrid( mza, mwa, nspcsd ) )
    cgrid = cmin

    allocate( gc_tend( mza, mwa, n_gc_trns) )
    gc_tend = rinit

    allocate( ae_tend( mza, mwa, n_ae_trns) )
    ae_tend = rinit

    allocate( nr_tend( mza, mwa, n_nr_trns) )
    nr_tend = rinit

!!    if (nl%do_emis /= 0) then
!!       allocate( vdemis_gc( mza, mwa, n_gc_emis) )
!!       vdemis_gc = rinit
!!
!!       allocate( vdemis_ae( mza, mwa, n_ae_emis) )
!!       vdemis_ae = rinit
!!
!!       allocate( vdemis_nr( mza, mwa, n_nr_emis) )
!!       vdemis_nr = rinit
!!    endif

!!    if (nl%do_drydep /= 0) then
!!       allocate( sxfer_gc( nsw_max, mwa, n_gc_depv ) )
!!       sxfer_gc = 0.0
!!
!!       allocate( sxfer_ae( nsw_max, mwa, n_ae_depv ) )
!!       sxfer_ae = 0.0
!!
!!       allocate( sxfer_nr( nsw_max, mwa, n_nr_depv ) )
!!       sxfer_nr = 0.0
!!    endif

    allocate(cgrid_names(nspcsd))
    cgrid_names = ' '

    allocate(gc_spc_to_trns_map(n_gc_spc))
    allocate(ae_spc_to_trns_map(n_ae_spc))
    allocate(nr_spc_to_trns_map(n_nr_spc))

    gc_spc_to_trns_map = 0
    ae_spc_to_trns_map = 0
    nr_spc_to_trns_map = 0

    do n = 1, n_gc_trns
       gc_spc_to_trns_map( gc_trns_map(n) ) = n
    enddo

    do n = 1, n_ae_trns
       ae_spc_to_trns_map( ae_trns_map(n) ) = n
    enddo

    do n = 1, n_nr_trns
       nr_spc_to_trns_map( nr_trns_map(n) ) = n
    enddo

  end subroutine alloc_cgrid


  subroutine filltab_cgrid()

    use var_tables, only: increment_vtable
    use misc_coms,  only: io6
    use cgrid_spcs
    implicit none

    integer :: n, ns, j
    integer :: n_spc_adv

    n_spc_adv = n_gc_trns + n_ae_trns + n_nr_trns
    write(io6,'(/,a,i0,a)') " There are ", n_spc_adv, " transported chemical species."

    do n = 1, n_gc_spc
       ns = gc_strt - 1 + n

       cgrid_names(ns) = gc_spc(n)

       if (any( n == gc_trns_map(1:n_gc_trns) )) then

          ! transported species need to be communicated and saved
          call increment_vtable( gc_spc(n), 'AW', rvar2=cgrid(:,:,ns), mpt1=.true.)

       else

          ! skip RXN counter species
          j = max(len_trim(gc_spc(n)), 3)
          if (gc_spc(n)(j-2:j) /= "RXN") then

             ! non-transported species only need to be saved
             call increment_vtable( gc_spc(n), 'AW', rvar2=cgrid(:,:,ns))

          endif

       endif
    enddo

    do n = 1, n_ae_spc
       ns = ae_strt - 1 + n

       cgrid_names(ns) = ae_spc(n)

       if (any( n == ae_trns_map(1:n_ae_trns) )) then

          ! transported species need to be communicated and saved
          call increment_vtable( ae_spc(n), 'AW', rvar2=cgrid(:,:,ns), mpt1=.true.)

       else

          ! non-transported species only need to be saved
          call increment_vtable( ae_spc(n), 'AW', rvar2=cgrid(:,:,ns))

       endif
    enddo

    do n = 1, n_nr_spc
       ns = nr_strt - 1 + n

       cgrid_names(ns) = nr_spc(n)

       if (any( n == nr_trns_map(1:n_nr_trns) )) then

          ! transported species need to be communicated and saved
          call increment_vtable( nr_spc(n), 'AW', rvar2=cgrid(:,:,ns), mpt1=.true.)

       else

          ! non-transported species only need to be saved
          call increment_vtable( nr_spc(n), 'AW', rvar2=cgrid(:,:,ns))

       endif
    enddo

  end subroutine filltab_cgrid


  subroutine cgrid_scalar_tabs()

    use cgrid_spcs
    use var_tables, only: vtables_scalar, scalar_tab, num_scalar
    use oname_coms, only: nl

    implicit none

    integer           :: n, nc, indxe, indxd, ng, nr, na
    integer, external :: findex

    ! Gas chemistry species

    do n = 1, n_gc_trns
       nc = gc_trns_map(n)

       call vtables_scalar( cgrid(:,:,nc), gc_tend(:,:,n), gc_trns(n) )
    enddo

    ! Aerosol species

    do n = 1, n_ae_trns
       nc = ae_trns_map(n) + ae_strt - 1

       call vtables_scalar( cgrid(:,:,nc), ae_tend(:,:,n), ae_trns(n) )
    enddo

    ! Non-reacting species

    do n = 1, n_nr_trns
       nc = nr_trns_map(n) + nr_strt - 1

       call vtables_scalar( cgrid(:,:,nc), nr_tend(:,:,n), nr_trns(n) )
    enddo

    ! Save some scalar indices for common species

    do n = 1, num_scalar
       if ( scalar_tab(n)%name == 'CO'  ) ns_co  = n
       if ( scalar_tab(n)%name == 'NO'  ) ns_no  = n
       if ( scalar_tab(n)%name == 'NO2' ) ns_no2 = n
       if ( scalar_tab(n)%name == 'O3'  ) ns_o3  = n
    enddo

  end subroutine cgrid_scalar_tabs

end module cgrid_defn
