subroutine cmaq_driver( mrl )

  use cgrid_conv, only: conv_cgrid, rev_cgrid
  use mem_grid,   only: lpw
  use mem_ijtabs, only: jtab_w, jtw_prog
  use cgrid_defn, only: cldfrac_3d

  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iw

  ! Diagnose cloud fraction for photolysis and aqueous

  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
     call get_cloud_frac(iw, lpw(iw), cldfrac_3d(:,iw))
  enddo

  ! Convert aerosol species to densities

  call rev_cgrid( mrl )

  ! Cloud/aqueous processes

  call cldproc( mrl )

  ! Gas chemistry

  call chem( mrl )

  ! Aerosol microphysics

  call aero( mrl )

  ! Convert aerosol species back to concentrations

  call conv_cgrid ( mrl )

end subroutine cmaq_driver

  
