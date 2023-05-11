subroutine inc_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, nz_avg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     sfluxt_davg, sfluxr_davg, airtempk_davg, cantempk_davg, &
     vegtempk_davg, soiltempk_davg

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: lpw
  use mem_basic,   only: tair, press, vxe, vye, vze
  use consts_coms, only: alvl, cp
  use misc_coms,   only: dtlong, iparallel
  use mem_cuparm,  only: conprr
  use mem_micro,   only: pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh
  use mem_radiate, only: rshort
  use therm_lib,   only: qwtk
  use mem_sfcg,    only: sfcg, mwsfc, itab_wsfc
  use mem_land,    only: land, mland, omland, nzg
  use mem_para,    only: myrank

  implicit none

  integer :: j, iw, k, iland, iwsfc
  real :: soiltempk, fracliq
  integer, dimension(nz_avg-1) :: level_list = (/ 8, 13, 18, 26, 30, 34 /)

  npoints_davg = npoints_davg + 1

  !$omp parallel

! Horizontal loop over all prognostic W/T points

  !$omp do private (iw,k)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

! Accumulate of sums for 2D surface quantities

     k = lpw(iw)

      press_davg(iw) =  press_davg(iw) + real(press(k,iw))
        vxe_davg(iw) =    vxe_davg(iw) + vxe(k,iw)
        vye_davg(iw) =    vye_davg(iw) + vye(k,iw)
        vze_davg(iw) =    vze_davg(iw) + vze(k,iw)
     rshort_davg(iw) = rshort_davg(iw) + rshort(iw)
      tempk_davg(iw) =  tempk_davg(iw) + tair(k,iw)

     if (tempk_dmin(iw) > tair(k,iw)) tempk_dmin(iw) = tair(k,iw)
     if (tempk_dmax(iw) < tair(k,iw)) tempk_dmax(iw) = tair(k,iw)

     if (allocated(pcprd)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprd(iw) * dtlong
     if (allocated(pcprr)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprr(iw) * dtlong
     if (allocated(pcprp)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprp(iw) * dtlong
     if (allocated(pcprs)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprs(iw) * dtlong
     if (allocated(pcpra)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcpra(iw) * dtlong
     if (allocated(pcprg)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprg(iw) * dtlong
     if (allocated(pcprh)) accpmic_dtot(iw) = accpmic_dtot(iw) &
                                            + pcprh(iw) * dtlong

     if (allocated(conprr)) then
        accpcon_dtot(iw) = accpcon_dtot(iw) + conprr(iw) * dtlong
     endif

     press_ul_davg(iw) = press_ul_davg(iw) + real(press(level_list(5),iw))
       vxe_ul_davg(iw) =   vxe_ul_davg(iw) + vxe(level_list(5),iw)
       vye_ul_davg(iw) =   vye_ul_davg(iw) + vye(level_list(5),iw)
       vze_ul_davg(iw) =   vze_ul_davg(iw) + vze(level_list(5),iw)

  enddo
  !$omp end do

! Horizontal loop over all SFC grid cells

  !$omp do
  do iwsfc = 2, mwsfc

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     airtempk_davg(iwsfc) = airtempk_davg(iwsfc) + sfcg%airtemp(iwsfc)
     cantempk_davg(iwsfc) = cantempk_davg(iwsfc) + sfcg%cantemp(iwsfc)
       sfluxt_davg(iwsfc) =   sfluxt_davg(iwsfc) + sfcg%sfluxt (iwsfc)
       sfluxr_davg(iwsfc) =   sfluxr_davg(iwsfc) + sfcg%sfluxr (iwsfc)
  enddo
  !$omp end do nowait

! Horizontal loop over all land cells

  !$omp do private(soiltempk,fracliq)
  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     call qwtk(land%soil_energy(nzg,iland),        &
               land%soil_water(nzg,iland)*1.e3,    &
               land%specifheat_drysoil(nzg,iland), &
               soiltempk, fracliq)

       vegtempk_davg(iland) =   vegtempk_davg(iland) + land%veg_temp(iland)
      soiltempk_davg(iland) =  soiltempk_davg(iland) + soiltempk
  enddo
  !$omp end do nowait

  !$omp end parallel
end subroutine inc_davg_vars

!========================================================================

subroutine norm_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, tempk_davg, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
         vegtempk_davg, soiltempk_davg, &
     cantempk_davg, airtempk_davg, sfluxt_davg, sfluxr_davg

  use mem_sfcg,  only: mwsfc, itab_wsfc
  use mem_land,  only: mland, omland
  use mem_ijtabs,only: jtab_w, jtw_prog
  use misc_coms, only: iparallel
  use mem_para,  only: myrank

  implicit none

  integer :: iw, iland, j, iwsfc
  real :: ni

  ni = 1. / npoints_davg

  !$omp parallel

  !$omp do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        press_davg(iw) =    press_davg(iw) * ni
          vxe_davg(iw) =      vxe_davg(iw) * ni
          vye_davg(iw) =      vye_davg(iw) * ni
          vze_davg(iw) =      vze_davg(iw) * ni
       rshort_davg(iw) =   rshort_davg(iw) * ni
        tempk_davg(iw) =    tempk_davg(iw) * ni
     press_ul_davg(iw) = press_ul_davg(iw) * ni
       vxe_ul_davg(iw) =   vxe_ul_davg(iw) * ni
       vye_ul_davg(iw) =   vye_ul_davg(iw) * ni
       vze_ul_davg(iw) =   vze_ul_davg(iw) * ni
  enddo
  !$omp end do nowait

  !$omp do
  do iwsfc = 2, mwsfc

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     airtempk_davg(iwsfc) = airtempk_davg(iwsfc) * ni
     cantempk_davg(iwsfc) = cantempk_davg(iwsfc) * ni
       sfluxt_davg(iwsfc) =   sfluxt_davg(iwsfc) * ni
       sfluxr_davg(iwsfc) =   sfluxr_davg(iwsfc) * ni
  enddo
  !$omp end do nowait

  !$omp do
  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

       vegtempk_davg(iland) =   vegtempk_davg(iland) * ni
      soiltempk_davg(iland) =  soiltempk_davg(iland) * ni
  enddo
  !$omp end do nowait

  !$omp end parallel
end subroutine norm_davg_vars

!===============================================================

subroutine write_davg_vars(outyear,outmonth,outdate)

  use mem_average_vars, only: &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     airtempk_davg, airtempk_dmin, airtempk_dmax, &
     cantempk_davg, cantempk_dmin, cantempk_dmax, &
     vegtempk_davg,   vegtempk_dmin,   vegtempk_dmax, &
     soiltempk_davg,  soiltempk_dmin,  soiltempk_dmax, &
     sfluxt_davg,     sfluxr_davg

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use misc_coms,  only: hfilepref, iclobber, io6
  use mem_para,   only: iwa_globe_primary, iwa_local_primary, &
                        iwsfc_globe_primary, iwsfc_local_primary, &
                        iland_globe_primary, iland_local_primary, &
                        ilake_globe_primary, ilake_local_primary, &
                        isea_globe_primary, isea_local_primary

  use mem_sfcg,   only: mwsfc, nwsfc
  use mem_land,   only: mland, nland

  implicit none

  integer, intent(in) :: outyear, outmonth, outdate

  character(len=128) :: hnamel
  integer :: lenl
  logical exans
  integer :: ndims,idims(2), month_use, year_use, date_use, nglobe

  integer, pointer, contiguous :: ilpts(:), igpts(:)

! Determine output file name and open file

  date_use = outdate - 1
  month_use = outmonth
  year_use = outyear
  if(date_use == 0)then
     month_use = outmonth - 1
     if(month_use == 0)then
        month_use = 12
        year_use = year_use - 1
     endif
     date_use = 31
     if(month_use == 4 .or. month_use == 6 .or. month_use == 9   &
          .or. month_use == 11)then
        date_use = 30
     elseif(month_use == 2)then
        date_use = 28
        if(mod(year_use,4) == 0)date_use = 29
     endif
  endif

  call makefnam8(hnamel,hfilepref,0.d0,year_use,month_use,date_use,  &
       0,'D','$','h5')  ! '$' argument suppresses grid# encoding

  lenl = len_trim(hnamel)

  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif

  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  write(io6,*)'history_write:open file:',trim(hnamel)
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber,trypario=.true.)

  !  Write day-average variables

  ndims    = 1
  idims(1) = mwa

  ilpts => iwa_local_primary
  igpts => iwa_globe_primary
  nglobe = nwa

  call shdf5_orec(ndims,idims,  'PRESS_DAVG',rvar1=  press_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VXE_DAVG',rvar1=    vxe_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VYE_DAVG',rvar1=    vye_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VZE_DAVG',rvar1=    vze_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'RSHORT_DAVG',rvar1= rshort_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DAVG',rvar1=  tempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DMIN',  rvar1=tempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DMAX',  rvar1=tempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'ACCPMIC_DTOT',rvar1=accpmic_dtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'ACCPCON_DTOT',rvar1=accpcon_dtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  call shdf5_orec(ndims,idims,'PRESS_UL_DAVG',rvar1=press_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VXE_UL_DAVG',rvar1=  vxe_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VYE_UL_DAVG',rvar1=  vye_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VZE_UL_DAVG',rvar1=  vze_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  idims(1) = mwsfc

  ilpts => iwsfc_local_primary
  igpts => iwsfc_globe_primary
  nglobe = nwsfc

  call shdf5_orec(ndims,idims,'AIRTEMPK_DAVG',rvar1=airtempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_DMIN',rvar1=airtempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_DMAX',rvar1=airtempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_DAVG',rvar1=cantempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_DMIN',rvar1=cantempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_DMAX',rvar1=cantempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXT_DAVG',rvar1=  sfluxt_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXR_DAVG',rvar1=  sfluxr_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  idims(1) = mland

  ilpts => iland_local_primary
  igpts => iland_globe_primary
  nglobe = nland

  call shdf5_orec(ndims,idims,  'VEGTEMPK_DAVG',rvar1=  vegtempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VEGTEMPK_DMIN',rvar1=  vegtempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VEGTEMPK_DMAX',rvar1=  vegtempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DAVG',rvar1= soiltempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DMIN',rvar1= soiltempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DMAX',rvar1= soiltempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  call shdf5_close()

end subroutine write_davg_vars

!===============================================================

subroutine read_davg_vars(davgfile)

  use mem_average_vars, only: &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
      airtempk_davg,  airtempk_dmin,  airtempk_dmax, &
      cantempk_davg,  cantempk_dmin,  cantempk_dmax, &
      vegtempk_davg,  vegtempk_dmin,  vegtempk_dmax, &
     soiltempk_davg, soiltempk_dmin, soiltempk_dmax, &
        sfluxt_davg,    sfluxr_davg

  use misc_coms,  only: io6
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa
  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: itab_wsfc, mwsfc
  use mem_land,   only: itab_land, mland

  implicit none

  character(*), intent(in) :: davgfile

  logical :: exans
  integer :: ndims, idims(2)

  integer :: ilocalw(mwa)
  integer :: ilocall(mland)
  integer :: ilocalc(mwsfc)

  inquire(file=davgfile,exist=exans)

  if (exans) then

! Day-average file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening day-average file ', trim(davgfile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(davgfile,'R',trypario=.true.)

 ! Read day-average variables

     ndims    = 1
     idims(1) = mwa
     ilocalw(1:mwa) = itab_w(1:mwa)%iwglobe

     call shdf5_irec(ndims,idims,  'PRESS_DAVG',rvar1=  press_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,    'VXE_DAVG',rvar1=    vxe_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,    'VYE_DAVG',rvar1=    vye_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,    'VZE_DAVG',rvar1=    vze_davg, points=ilocalw)
     call shdf5_irec(ndims,idims, 'RSHORT_DAVG',rvar1= rshort_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,  'TEMPK_DAVG',rvar1=  tempk_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,  'TEMPK_DMIN',rvar1=  tempk_dmin, points=ilocalw)
     call shdf5_irec(ndims,idims,  'TEMPK_DMAX',rvar1=  tempk_dmax, points=ilocalw)
     call shdf5_irec(ndims,idims,'ACCPMIC_DTOT',rvar1=accpmic_dtot, points=ilocalw)
     call shdf5_irec(ndims,idims,'ACCPCON_DTOT',rvar1=accpcon_dtot, points=ilocalw)

     call shdf5_irec(ndims,idims,'PRESS_UL_DAVG',rvar1=press_ul_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,  'VXE_UL_DAVG',rvar1=  vxe_ul_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,  'VYE_UL_DAVG',rvar1=  vye_ul_davg, points=ilocalw)
     call shdf5_irec(ndims,idims,  'VZE_UL_DAVG',rvar1=  vze_ul_davg, points=ilocalw)

     idims(1) = mwsfc
     ilocalc(1:mwsfc) = itab_wsfc(1:mwsfc)%iwglobe

     call shdf5_irec(ndims,idims,'AIRTEMPK_DAVG',rvar1=airtempk_davg, points=ilocalc)
     call shdf5_irec(ndims,idims,'AIRTEMPK_DMIN',rvar1=airtempk_dmin, points=ilocalc)
     call shdf5_irec(ndims,idims,'AIRTEMPK_DMAX',rvar1=airtempk_dmax, points=ilocalc)
     call shdf5_irec(ndims,idims,'CANTEMPK_DAVG',rvar1=cantempk_davg, points=ilocalc)
     call shdf5_irec(ndims,idims,'CANTEMPK_DMIN',rvar1=cantempk_dmin, points=ilocalc)
     call shdf5_irec(ndims,idims,'CANTEMPK_DMAX',rvar1=cantempk_dmax, points=ilocalc)
     call shdf5_irec(ndims,idims,  'SFLUXT_DAVG',rvar1=  sfluxt_davg, points=ilocalc)
     call shdf5_irec(ndims,idims,  'SFLUXR_DAVG',rvar1=  sfluxr_davg, points=ilocalc)

     idims(1) = mland
     ilocall(1:mland) = itab_land(1:mland)%iwglobe

     call shdf5_irec(ndims,idims,  'VEGTEMPK_DAVG',rvar1=  vegtempk_davg, points=ilocall)
     call shdf5_irec(ndims,idims,  'VEGTEMPK_DMIN',rvar1=  vegtempk_dmin, points=ilocall)
     call shdf5_irec(ndims,idims,  'VEGTEMPK_DMAX',rvar1=  vegtempk_dmax, points=ilocall)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DAVG',rvar1= soiltempk_davg, points=ilocall)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DMIN',rvar1= soiltempk_dmin, points=ilocall)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DMAX',rvar1= soiltempk_dmax, points=ilocall)

     call shdf5_close()
  else

   ! Day-average file does not exist, stop model.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open day-average file:'
     write(io6,*) '!!!   '//trim(davgfile)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in read_davg_vars'

  endif

end subroutine read_davg_vars

!=============================================================================

subroutine makefnam8(fname,prefix,tinc8,iyr,imn,idy,itm,ftype,post,fmt)

! creates standard timestamped filename

  implicit none

  character(*), intent(out) :: fname
  character(*), intent(in) :: prefix,post,fmt
  real(kind=8), intent(in) :: tinc8
  integer, intent(in) :: iyr, imn, idy, itm
  character(*), intent(in) :: ftype

  integer :: oyr, omn, ody, otm

  character(len=40) :: dstring

  if (tinc8 == 0.) then
     oyr = iyr ; omn = imn ; ody = idy ; otm = itm
  else
     call date_add_to8(iyr,imn,idy,itm,tinc8,'s',oyr,omn,ody,otm)
  endif

  write(dstring,100) '-',trim(ftype),'-',oyr,'-',omn,'-',ody,'-',otm
  100 format(3a,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

  fname = trim(prefix)//trim(dstring)
  if (post(1:1) /= '$') then
     fname = trim(fname)//'_'//trim(post)
  endif
  fname = trim(fname)//'.'//trim(fmt)

end subroutine makefnam8
