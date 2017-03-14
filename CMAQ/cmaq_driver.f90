subroutine cmaq_driver( mrl )

  use cgrid_conv, only: conv_cgrid, rev_cgrid
  use sedv_defn,  only: aero_sedi

  implicit none

  integer, intent(in) :: mrl

  ! Cloud/aqueous processes

  call rescld( mrl )

  ! Convert aerosol species to densities

  call rev_cgrid( mrl )

  ! Aerosol sedimentation

  call aero_sedi( mrl )

  ! Gas chemistry

  call chem( mrl )

  ! Aerosol microphysics

  call aero( mrl )

  ! Convert aerosol species back to concentrations

  call conv_cgrid ( mrl )

end subroutine cmaq_driver

  
