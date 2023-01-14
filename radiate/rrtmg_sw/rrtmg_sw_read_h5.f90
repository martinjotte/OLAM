subroutine sw_kgb16_h5
  use rrsw_kg16, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no16
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 1, numGPoints = no16
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo16', rvar1=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl16',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao16',       rvar4=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo16',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo16',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo16',   rvar2=forrefo)

end subroutine sw_kgb16_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb17_h5
  use rrsw_kg17, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no17
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 2
  integer(kind=im), parameter :: numGPoints = no17
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keyupper/),              'sfluxrefo17', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl17',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao17',       rvar4=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo17',       rvar4=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo17',  rvar2=selfrefo)
  call shdf5_irec(2, (/4,numGPoints/),                     'forrefo17',   rvar2=forrefo)

end subroutine sw_kgb17_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb18_h5
  use rrsw_kg18, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no18
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 3
  integer(kind=im), parameter :: numGPoints = no18
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'sfluxrefo18', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl18',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao18',       rvar4=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo18',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo18',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo18',   rvar2=forrefo)

end subroutine sw_kgb18_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb19_h5
  use rrsw_kg19, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no19
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 4
  integer(kind=im), parameter :: numGPoints = no19
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'sfluxrefo19', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl19',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao19',       rvar4=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo19',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo19',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo19',   rvar2=forrefo)

end subroutine sw_kgb19_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb20_h5
  use rrsw_kg20, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, absch4o, no20
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 5
  integer(kind=im), parameter :: numGPoints = no20
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo20', rvar1=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl20',      rvars=rayl)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao20',       rvar3=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo20',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo20',  rvar2=selfrefo)
  call shdf5_irec(2, (/4,numGPoints/),                     'forrefo20',   rvar2=forrefo)
  call shdf5_irec(1, (/numGPoints/),                       'absch4o20',   rvar1=absch4o)

end subroutine sw_kgb20_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb21_h5
  use rrsw_kg21, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no21
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 6
  integer(kind=im), parameter :: numGPoints = no21
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'sfluxrefo21', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl21',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao21',       rvar4=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo21',       rvar4=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo21',  rvar2=selfrefo)
  call shdf5_irec(2, (/4,numGPoints/),                     'forrefo21',   rvar2=forrefo)

end subroutine sw_kgb21_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb22_h5
  use rrsw_kg22, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no22
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 7
  integer(kind=im), parameter :: numGPoints = no22
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'sfluxrefo22', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl22',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao22',       rvar4=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo22',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo22',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo22',   rvar2=forrefo)

end subroutine sw_kgb22_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb23_h5
  use rrsw_kg23, only: sfluxrefo, kao, selfrefo, forrefo, raylo, no23
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 8
  integer(kind=im), parameter :: numGPoints = no23
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo23', rvar1=sfluxrefo)
  call shdf5_irec(1, (/numGPoints/),                       'raylo23',     rvar1=raylo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao23',       rvar3=kao)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo23',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo23',   rvar2=forrefo)

end subroutine sw_kgb23_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb24_h5
  use rrsw_kg24, only: sfluxrefo, kao, kbo, selfrefo, forrefo, &
                       raylao, raylbo, abso3ao, abso3bo, no24
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 9
  integer(kind=im), parameter :: numGPoints = no24
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'sfluxrefo24', rvar2=sfluxrefo)
  call shdf5_irec(2, (/numGPoints,keylower/),              'raylao24',    rvar2=raylao)
  call shdf5_irec(1, (/numGPoints/),                       'raylbo24',    rvar1=raylbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao24',       rvar4=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo24',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo24',  rvar2=selfrefo)
  call shdf5_irec(2, (/Tforeignlower,numGPoints/),         'forrefo24',   rvar2=forrefo)
  call shdf5_irec(1, (/numGPoints/),                       'abso3ao24',   rvar1=abso3ao)
  call shdf5_irec(1, (/numGPoints/),                       'abso3bo24',   rvar1=abso3bo)

end subroutine sw_kgb24_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb25_h5
  use rrsw_kg25, only: sfluxrefo, kao, raylo, abso3ao, abso3bo, no25
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 10
  integer(kind=im), parameter :: numGPoints = no25
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo25', rvar1=sfluxrefo)
  call shdf5_irec(1, (/numGPoints/),                       'raylo25',     rvar1=raylo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao25',       rvar3=kao)
  call shdf5_irec(1, (/numGPoints/),                       'abso3ao25',   rvar1=abso3ao)
  call shdf5_irec(1, (/numGPoints/),                       'abso3bo25',   rvar1=abso3bo)

end subroutine sw_kgb25_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb26_h5
  use rrsw_kg26, only: sfluxrefo, raylo, no26
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 11
  integer(kind=im), parameter :: numGPoints = no26
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo26', rvar1=sfluxrefo)
  call shdf5_irec(1, (/numGPoints/),                       'raylo26',     rvar1=raylo)

end subroutine sw_kgb26_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb27_h5
  use rrsw_kg27, only: sfluxrefo, kao, kbo, raylo, no27
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 12
  integer(kind=im), parameter :: numGPoints = no27
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo27', rvar1=sfluxrefo)
  call shdf5_irec(1, (/numGPoints/),                       'raylo27',     rvar1=raylo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao27',       rvar3=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo27',       rvar3=kbo)

end subroutine sw_kgb27_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb28_h5
  use rrsw_kg28, only: sfluxrefo, kao, kbo, rayl, no28
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 13
  integer(kind=im), parameter :: numGPoints = no28
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keyupper/),              'sfluxrefo28', rvar2=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl28',      rvars=rayl)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao28',       rvar4=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo28',       rvar4=kbo)

end subroutine sw_kgb28_h5
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb29_h5
  use rrsw_kg29, only: sfluxrefo, kao, kbo, selfrefo, forrefo, &
                       absh2oo, absco2o, rayl, no29
  use rrsw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 14
  integer(kind=im), parameter :: numGPoints = no29
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'sfluxrefo29', rvar1=sfluxrefo)
  call shdf5_irec(1, (/1/),                                'rayl29',      rvars=rayl)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao29',       rvar3=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo29',       rvar3=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo29',  rvar2=selfrefo)
  call shdf5_irec(2, (/4,numGPoints/),                     'forrefo29',   rvar2=forrefo)
  call shdf5_irec(1, (/numGPoints/),                       'absh2oo29',   rvar1=absh2oo)
  call shdf5_irec(1, (/numGPoints/),                       'absco2o29',   rvar1=absco2o)

end subroutine sw_kgb29_h5
!*******************************************************************************
