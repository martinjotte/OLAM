subroutine cmaq_driver( mrl )

  use cgrid_conv, only: conv_cgrid, rev_cgrid

  implicit none

  integer, intent(in) :: mrl

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

  
