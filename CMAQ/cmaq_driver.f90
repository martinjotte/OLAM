subroutine cmaq_driver( )

  use cgrid_conv, only: conv_cgrid, rev_cgrid
  use sedv_defn,  only: aero_sedi
  use mem_ijtabs, only: istp, mrl_endl

  implicit none

  integer :: mrl

  mrl = mrl_endl(istp)

  ! Photolysis rates

  if (istp == 1) call phot( )

  ! Cloud/aqueous processes

  if (mrl > 0) then
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
  endif

end subroutine cmaq_driver
