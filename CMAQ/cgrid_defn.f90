module cgrid_defn

  real, allocatable :: cgrid(:,:,:)

  real, allocatable :: gc_tend(:,:,:)
  real, allocatable :: ae_tend(:,:,:)
  real, allocatable :: nr_tend(:,:,:)

  real, allocatable :: vdemis_gc(:,:,:)
  real, allocatable :: vdemis_ae(:,:,:)
  real, allocatable :: vdemis_nr(:,:,:)

  real, allocatable :: sxfer_gc(:,:,:)
  real, allocatable :: sxfer_ae(:,:,:)
  real, allocatable :: sxfer_nr(:,:,:)

  integer :: ns_co  = 0
  integer :: ns_no  = 0
  integer :: ns_no2 = 0
  integer :: ns_o3  = 0

  character(16), allocatable :: cgrid_names(:)

contains

  subroutine alloc_cgrid(mza, mwa)

    use mem_grid,   only: nsw_max
    use cgrid_spcs, only: n_gc_trns, n_ae_trns, n_nr_trns, &
                          n_gc_emis, n_ae_emis, n_nr_emis, &
                          n_gc_depv, n_ae_depv, n_nr_depv, &
                          nspcsd
    use misc_coms,  only: rinit
    implicit none

    integer, intent(in) :: mza, mwa

    real,    parameter  :: cmin = 1.e-30

    allocate( cgrid( mza, mwa, nspcsd ) )
    cgrid = cmin

    allocate( gc_tend( mza, mwa, n_gc_trns) )
    gc_tend = rinit

    allocate( ae_tend( mza, mwa, n_ae_trns) )
    ae_tend = rinit

    allocate( nr_tend( mza, mwa, n_nr_trns) )
    nr_tend = rinit

    allocate( vdemis_gc( mza, mwa, n_gc_emis) )
    vdemis_gc = rinit

    allocate( vdemis_ae( mza, mwa, n_ae_emis) )
    vdemis_ae = rinit

    allocate( vdemis_nr( mza, mwa, n_nr_emis) )
    vdemis_nr = rinit

    allocate( sxfer_gc( nsw_max, mwa, n_gc_depv ) )
    sxfer_gc = 0.0

    allocate( sxfer_ae( nsw_max, mwa, n_ae_depv ) )
    sxfer_ae = 0.0

    allocate( sxfer_nr( nsw_max, mwa, n_nr_depv ) )
    sxfer_nr = 0.0

    allocate(cgrid_names(nspcsd))
    cgrid_names = ''

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
    implicit none

    integer           :: n, nc, indxe, indxd, ng, nr, na
    integer, external :: findex

    ! Gas chemistry species

    do n = 1, n_gc_trns
       ng = gc_trns_map(n)
       nc = gc_strt - 1 + ng

       indxe = findex( ng, n_gc_emis, gc_emis_map )
       indxd = findex( ng, n_gc_depv, gc_depv_map )

       if (indxe > 0 .and. indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), gc_tend(:,:,n), gc_trns(n), &
                               emis   = vdemis_gc(:,:,indxe),             &
                               sxfer  = sxfer_gc(:,:,indxd),              &
                               cu_mix = .true. )

       elseif (indxe > 0) then

          call vtables_scalar( cgrid(:,:,nc), gc_tend(:,:,n), gc_trns(n), &
                               emis   = vdemis_gc(:,:,indxe),             &
                               cu_mix = .true. )

       elseif (indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), gc_tend(:,:,n), gc_trns(n), &
                               sxfer  = sxfer_gc(:,:,indxd),              &
                               cu_mix = .true. )

       else

          call vtables_scalar( cgrid(:,:,nc), gc_tend(:,:,n), gc_trns(n), &
                               cu_mix = .true. )

       endif

    enddo

    ! Aerosol species

    do n = 1, n_ae_trns
       na = ae_trns_map(n)
       nc = ae_strt - 1 + na

       indxe = findex( na, n_ae_emis, ae_emis_map )
       indxd = findex( na, n_ae_depv, ae_depv_map )

       if (indxe > 0 .and. indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), ae_tend(:,:,n), ae_trns(n), &
                               emis   = vdemis_ae(:,:,indxe),             &
                               sxfer  = sxfer_ae(:,:,indxd),              &
                               cu_mix = .true. )

       elseif (indxe > 0) then

          call vtables_scalar( cgrid(:,:,nc), ae_tend(:,:,n), ae_trns(n), &
                               emis   = vdemis_ae(:,:,indxe),             &
                               cu_mix = .true. )

       elseif (indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), ae_tend(:,:,n), ae_trns(n), &
                               sxfer  = sxfer_ae(:,:,indxd),              &
                               cu_mix = .true. )

       else

          call vtables_scalar( cgrid(:,:,nc), ae_tend(:,:,n), ae_trns(n), &
                               cu_mix = .true. )

       endif

    enddo

    ! Non-reacting species

    do n = 1, n_nr_trns
       nr = nr_trns_map(n)
       nc = nr_strt - 1 + nr

       indxe = findex( nr, n_nr_emis, nr_emis_map )
       indxd = findex( nr, n_nr_depv, nr_depv_map )

       if (indxe > 0 .and. indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), nr_tend(:,:,n), nr_trns(n), &
                               emis   = vdemis_nr(:,:,indxe),             &
                               sxfer  = sxfer_nr(:,:,indxd),              &
                               cu_mix = .true. )

       elseif (indxe > 0) then

          call vtables_scalar( cgrid(:,:,nc), nr_tend(:,:,n), nr_trns(n), &
                               emis   = vdemis_nr(:,:,indxe),             &
                               cu_mix = .true. )

       elseif (indxd > 0) then

          call vtables_scalar( cgrid(:,:,nc), nr_tend(:,:,n), nr_trns(n), &
                               sxfer  = sxfer_nr(:,:,indxd),              &
                               cu_mix = .true. )

       else

          call vtables_scalar( cgrid(:,:,nc), nr_tend(:,:,n), nr_trns(n), &
                               cu_mix = .true. )

       endif

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
