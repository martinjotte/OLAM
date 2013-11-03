subroutine lw_kgb01_h5
  use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                        selfrefo, forrefo, no1
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 1, numGPoints = no1
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao01', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo01', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao01',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo01',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo01',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo01',   rvara=forrefo)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mn201',   rvara=kao_mn2)
  call shdf5_irec(2, (/T,numGPoints/),                     'kbo_mn201',   rvara=kbo_mn2)

end subroutine lw_kgb01_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb02_h5
  use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no2
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 2
  integer(kind=im), parameter :: numGPoints = no2
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao02', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo02', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao02',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo02',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo02',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo02',   rvara=forrefo)

end subroutine lw_kgb02_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb03_h5
  use rrlw_kg03, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                        kbo_mn2o, selfrefo, forrefo, no3
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 3
  integer(kind=im), parameter :: numGPoints = no3
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao03', rvara=fracrefao)
  call shdf5_irec(2, (/numGPoints,keyupper/),              'fracrefbo03', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao03',       rvara=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo03',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo03',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo03',   rvara=forrefo)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mn2o03',  rvara=kao_mn2o)
  call shdf5_irec(3, (/keyupper,T,numGPoints/),            'kbo_mn2o03',  rvara=kbo_mn2o)

end subroutine lw_kgb03_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb04_h5
  use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no4
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 4
  integer(kind=im), parameter :: numGPoints = no4
  integer(kind=im), parameter :: gPointSetxNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao04', rvara=fracrefao)
  call shdf5_irec(2, (/numGPoints,keyupper/),              'fracrefbo04', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao04',       rvara=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo04',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo04',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo04',   rvara=forrefo)

end subroutine lw_kgb04_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb05_h5
  use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
                        selfrefo, forrefo, ccl4o, no5
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 5
  integer(kind=im), parameter :: numGPoints = no5
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao05', rvara=fracrefao)
  call shdf5_irec(2, (/numGPoints,keyupper/),              'fracrefbo05', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao05',       rvara=kao)
  call shdf5_irec(4, (/keyupper,Tdiff,pupper,numGPoints/), 'kbo05',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo05',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo05',   rvara=forrefo)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mo305',   rvara=kao_mo3)
  call shdf5_irec(1, (/numGPoints/),                       'ccl4o05',     rvara=ccl4o)

end subroutine lw_kgb05_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb06_h5
  use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
                        cfc11adjo, cfc12o, no6
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 6
  integer(kind=im), parameter :: numGPoints = no6
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao06', rvara=fracrefao)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao06',       rvara=kao)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo06',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo06',   rvara=forrefo)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mco206',  rvara=kao_mco2)
  call shdf5_irec(1, (/numGPoints/),                       'cfc11adjo06', rvara=cfc11adjo)
  call shdf5_irec(1, (/numGPoints/),                       'cfc12o06',    rvara=cfc12o)

end subroutine lw_kgb06_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb07_h5
  use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
                        kbo_mco2, selfrefo, forrefo, no7
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 7
  integer(kind=im), parameter :: numGPoints = no7
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao07', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo07', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao07',       rvara=kao)
  call shdf5_irec(3, (/         Tdiff,pupper,numGPoints/), 'kbo07',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo07',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo07',   rvara=forrefo)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mco207',  rvara=kao_mco2)
  call shdf5_irec(2, (/         T,numGPoints/),            'kbo_mco207',  rvara=kbo_mco2)

end subroutine lw_kgb07_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb08_h5
  use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                        kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                        cfc12o, cfc22adjo, no8
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 8
  integer(kind=im), parameter :: numGPoints = no8
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao08', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo08', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao08',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo08',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo08',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo08',   rvara=forrefo)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mo308',   rvara=kao_mo3)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mco208',  rvara=kao_mco2)
  call shdf5_irec(2, (/T,numGPoints/),                     'kbo_mco208',  rvara=kbo_mco2)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mn2o08',  rvara=kao_mn2o)
  call shdf5_irec(2, (/T,numGPoints/),                     'kbo_mn2o08',  rvara=kbo_mn2o)
  call shdf5_irec(1, (/numGPoints/),                       'cfc12o08',    rvara=cfc12o)
  call shdf5_irec(1, (/numGPoints/),                       'cfc22adjo08', rvara=cfc22adjo)

end subroutine lw_kgb08_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb09_h5
  use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                        kbo_mn2o, selfrefo, forrefo, no9
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 9
  integer(kind=im), parameter :: numGPoints = no9
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao09', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo09', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao09',       rvara=kao)
  call shdf5_irec(3, (/         Tdiff,pupper,numGPoints/), 'kbo09',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo09',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo09',   rvara=forrefo)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mn2o09',  rvara=kao_mn2o)
  call shdf5_irec(2, (/         T,numGPoints/),            'kbo_mn2o09',  rvara=kbo_mn2o)

end subroutine lw_kgb09_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb10_h5
  use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no10
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 10
  integer(kind=im), parameter :: numGPoints = no10
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao10', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo10', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao10',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo10',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo10',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo10',   rvara=forrefo)

end subroutine lw_kgb10_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb11_h5
  use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
                        kbo_mo2, selfrefo, forrefo, no11
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 11
  integer(kind=im), parameter :: numGPoints = no11
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao11', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo11', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao11',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo11',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo11',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo11',   rvara=forrefo)
  call shdf5_irec(2, (/T,numGPoints/),                     'kao_mo211',   rvara=kao_mo2)
  call shdf5_irec(2, (/T,numGPoints/),                     'kbo_mo211',   rvara=kbo_mo2)

end subroutine lw_kgb11_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb12_h5
  use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo, no12
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 12
  integer(kind=im), parameter :: numGPoints = no12
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao12', rvara=fracrefao)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao12',       rvara=kao)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo12',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo12',   rvara=forrefo)

end subroutine lw_kgb12_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb13_h5
  use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                        kbo_mo3, selfrefo, forrefo, no13
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 13
  integer(kind=im), parameter :: numGPoints = no13  
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao13', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo13', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao13',       rvara=kao)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo13',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo13',   rvara=forrefo)
  call shdf5_irec(2, (/         T,numGPoints/),            'kbo_mo313',   rvara=kbo_mo3)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mco213',  rvara=kao_mco2)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mco13',   rvara=kao_mco)

end subroutine lw_kgb13_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb14_h5
  use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no14
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 14
  integer(kind=im), parameter :: numGPoints = no14
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(1, (/numGPoints/),                       'fracrefao14', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo14', rvara=fracrefbo)
  call shdf5_irec(3, (/Tdiff,plower,numGPoints/),          'kao14',       rvara=kao)
  call shdf5_irec(3, (/Tdiff,pupper,numGPoints/),          'kbo14',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo14',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo14',   rvara=forrefo)

end subroutine lw_kgb14_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb15_h5
  use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo, no15
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 15
  integer(kind=im), parameter :: numGPoints = no15
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao15', rvara=fracrefao)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao15',       rvara=kao)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo15',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo15',   rvara=forrefo)
  call shdf5_irec(3, (/keylower,T,numGPoints/),            'kao_mn215',   rvara=kao_mn2)

end subroutine lw_kgb15_h5
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb16_h5
  use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no16
  use rrlw_ncpar
  use hdf5_utils

  implicit none

  integer(kind=im), parameter :: bandNumber = 16
  integer(kind=im), parameter :: numGPoints = no16
  integer(kind=im), parameter :: gPointSetNumber = 1

  call shdf5_irec(2, (/numGPoints,keylower/),              'fracrefao16', rvara=fracrefao)
  call shdf5_irec(1, (/numGPoints/),                       'fracrefbo16', rvara=fracrefbo)
  call shdf5_irec(4, (/keylower,Tdiff,plower,numGPoints/), 'kao16',       rvara=kao)
  call shdf5_irec(3, (/         Tdiff,pupper,numGPoints/), 'kbo16',       rvara=kbo)
  call shdf5_irec(2, (/Tself,numGPoints/),                 'selfrefo16',  rvara=selfrefo)
  call shdf5_irec(2, (/Tforeign,numGPoints/),              'forrefo16',   rvara=forrefo)

end subroutine lw_kgb16_h5
!*******************************************************************************
