subroutine oplot_chem_lib(kk,ii,infotyp,fldname0,wtbot,wttop,fldval,notavail)

  use oplot_coms,   only: op
  use mem_para,     only: myrank
  use mem_grid,     only: mza, lpm, lpv, lpw, arw0
  use soil_nox,     only: soilnox, pfactor, drytime
  use misc_coms,    only: do_chem, io6
  use mem_ijtabs,   only: itab_w, jtab_w, jtab_v, jtab_m, jtw_prog, jtv_prog, jtm_vadj
  use mem_sfcg,     only: itab_wsfc, sfcg
  use mem_land,     only: mland, omland
  use mem_sea,      only: msea, omsea
  use soil_nox,     only: soilnox, pfactor, drytime
  use ltng_defn,    only: column_ltng_no, vdemis_lt
  use depv_defn,    only: depvel_gas, n_gas_depv, gas_depv_sur
  use cgrid_spcs,   only: gc_spc, gc_depv_map


  implicit none

  integer,      intent(in)  :: kk, ii
  character(*), intent(in)  :: infotyp, fldname0
  real,         intent(in)  :: wtbot, wttop
  real,         intent(out) :: fldval
  integer,      intent(out) :: notavail ! 0 - field available
                                        ! 1 - field below ground
                                        ! 2 - field above model top
                                        ! 3 - field not available in this run
                                        ! 4 - field not available in current grid cell
                                        !     (e.g., no sfcwater fracliq when no sfcwater)
  integer, parameter :: nfields = 9
  character(40)      :: fldlib(4,nfields)
  character(40)      :: fldname
  integer            :: i, k, ifield, kp, n, jsfc, iwsfc, jasfc, isea
  integer, save      :: icase

  integer, save      :: io3_depv = -1

!  fldname    stagpt/dimens      field description & units                 field #
!---------------------------------------------------------------------------------

  ! Soil NOx fields

  data fldlib(1:4,  1:6) / &
    'SOILNOX'       ,'L2' ,'SOIL NOX EMISSIONS', ' [g(N) km:S:-2:N: s:S:-1:N:]'    ,& !  1
    'SOILNOXT'      ,'T2' ,'SOIL NOX EMISSIONS', ' [g(N) km:S:-2:N: s:S:-1:N:]'    ,& !  2
    'SNPULSE'       ,'L2' ,'SOIL NOX PULSE FACTOR', ' [ ]'                         ,& !  3
    'SNPULSET'      ,'T2' ,'SOIL NOX PULSE FACTOR', ' [ ]'                         ,& !  4
    'SNDRYTIME'     ,'L2' ,'SOIL NOX DRY TIME', ' [hr]'                            ,& !  5
    'SNDRYTIMET'    ,'T2' ,'SOIL NOX DRY TIME', ' [hr]'                             / !  6

  ! Lightning NOx fields
  data fldlib(1:4,  7:8) / &
    'LTNGNOX_TOT'   ,'T2' ,'LTNG NOX COLUMN TOTAL', ' [g(N) km:S:-2:N: s:S:-1:N:]' ,& !  7
    'LTNGNOX'       ,'T3' ,'LTNG NOX', ' [mol (N) s:S:-1:N:]'                       / !  8

  ! Deposition velocities
  data fldlib(1:4,  9) / &
    'O3_DEPV'       ,'T2' ,'OZONE DEPOSITION VELOCITY', ' [cm s:S2:-1:]'  / ! 10

  if (do_chem /= 1) then
     write(*,*) 'Plot field ',trim(fldname0),' not available; chemistry not selected.'
     notavail = 3
     return
  endif

  if (io3_depv < 0) then
     io3_depv = 0
     do n = 1, n_gas_depv
        if (gc_spc(gc_depv_map(n)) == 'O3') then
           io3_depv = gas_depv_sur(n)
           exit
        endif
     enddo
  endif

  k = kk
  i = ii

  if (infotyp == 'UNITS') then

! Print field name to be plotted

     if (myrank == 0) write(io6,*) 'oplib ',trim(fldname0)

! Default values

     fldname = fldname0(6:)
     op%label = ' '
     op%units = ' '

     do ifield = 1, nfields
        if (fldname == fldlib(1,ifield)) then
           icase = ifield
           exit
        endif

        if (ifield == nfields) then
           write(io6,*) 'Plot field ',trim(fldname),' not available in oplot_chem_lib.'
           notavail = 3
           return
        endif
     enddo

     op%stagpt = fldlib(2,icase)(1:1)
     op%dimens = fldlib(2,icase)(2:3)
     op%label  = fldlib(3,icase)
     op%units  = fldlib(4,icase)

     op%fldval_min  =  1.e30
     op%fldval_max  = -1.e30
     op%fldvalv_min =  1.e30
     op%fldvalv_max = -1.e30

     if (op%stagpt == 'T' .or. op%stagpt == 'W') then
        i = jtab_w(jtw_prog)%iw(1)
        k = lpw(i)
     elseif (op%stagpt == 'V') then
        i = jtab_v(jtv_prog)%iv(1)
        k = lpv(i)
     elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
        i = jtab_m(jtm_vadj)%im(1)
        k = lpm(i)
     elseif (op%stagpt == 'L' .or. op%stagpt == 'S' .or. op%stagpt == 'B') then
        i = 2
        k = 1
     endif
  endif

  if (op%stagpt == 'L' .and. mland < 2) then
     notavail = 3
     return
  endif

  if (op%stagpt == 'S' .and. msea < 2) then
     notavail = 3
     return
  endif

  notavail = 0
  kp = min(k+1,mza)

! Execute IF block below even when infotyp == 'UNITS'
! in order to check whether current plot field is available in this model run.
! For this type of call to this subroutine, (k,i) are passed in as (1,2).

  select case(icase)

  case(1) ! SOILNOX

     if (.not. allocated(soilnox)) goto 1000
     iwsfc = i
     fldval = soilnox(i) * 14.0 / sfcg%area(iwsfc) * 1.e6  ! convert moles N to grams N

  case(2) ! SOILNOXT

     if (.not. allocated(soilnox)) goto 1000
     fldval = average_land(i,soilnox) * 14.0 / arw0(i) * 1.e6 ! convert moles N to grams N

  case(3) ! SNPULSE

     if (.not. allocated(pfactor)) goto 1000
     fldval = pfactor(i)

  case(4) ! SNPULSET

     if (.not. allocated(pfactor)) goto 1000
     fldval = average_land(i,pfactor)

  case(5) ! SNDRYTIME

     if (.not. allocated(drytime)) goto 1000
     fldval = drytime(i)

  case(6) ! SNDRYTIMET

     if (.not. allocated(drytime)) goto 1000
     fldval = average_land(i,drytime)

  case(7) ! LTNG NOX_TOT

     if (.not. allocated(column_ltng_no)) goto 1000
     fldval = column_ltng_no(i) * 14.0 / arw0(i) * 1.e6

  case(8) ! LTNG NOX [mol / sec]

     if (.not. allocated(vdemis_lt)) goto 1000
     fldval = vdemis_lt(k,i)

  case(9)

     fldval = DEPVEL_GAS(i,io3_depv) * 100.0

  end select

  return

1000 continue

  notavail = 3


contains


    real function average_land(iw,field)
      implicit none

      real,    intent(in) :: field(mland)
      integer, intent(in) :: iw
      integer             :: iland

      average_land = 0.0

      if (itab_w(iw)%jland2 > 0) then

         do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
            iwsfc = itab_w(iw)%iwsfc(jsfc)
            jasfc = itab_w(iw)%jasfc(jsfc)

            iland = iwsfc - omland

            average_land = average_land + field(iland) * itab_wsfc(iwsfc)%arcoariw(jasfc)
         enddo

      endif

    end function average_land


    real function average_sea(iw,field)
      implicit none

      real,    intent(in) :: field(msea)
      integer, intent(in) :: iw

      average_sea = 0.0

      if (itab_w(iw)%jsea2 > 0) then

         do jsfc = itab_w(iw)%jsea1, itab_w(iw)%jsea2
            iwsfc = itab_w(iw)%iwsfc(jsfc)
            jasfc = itab_w(iw)%jasfc(jsfc)

            isea = iwsfc - omsea

            average_sea = average_sea + field(isea) * itab_wsfc(iwsfc)%arcoariw(jasfc)
         enddo

      endif

    end function average_sea

end subroutine oplot_chem_lib
