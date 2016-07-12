subroutine cmaq_driver( mrl )

  use cgrid_conv, only: conv_cgrid, rev_cgrid
  use mem_grid,   only: lpw
  use mem_ijtabs, only: jtab_w, jtw_prog

  implicit none

  integer, intent(in) :: mrl

  ! Cloud/aqueous processes

  call rescld( mrl )

  ! Convert aerosol species to densities

  call rev_cgrid( mrl )

  ! Gas chemistry

  call chem( mrl )

  ! Aerosol microphysics

  call aero( mrl )

  ! Convert aerosol species back to concentrations

  call conv_cgrid ( mrl )

end subroutine cmaq_driver

  
