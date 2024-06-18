subroutine sea_init_atm()

  use mem_sea,     only: sea, msea, omsea
  use sea_coms,    only: iupdsst, s1900_sst, isstfile, nsstfiles, dt_sea,  &
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles, nzi
  use misc_coms,   only: s1900_sim, iparallel, runtype
  use mem_sfcg,    only: itab_wsfc, sfcg, nswmzons
  use consts_coms, only: t00
  use therm_lib,   only: rhovsl, rhovsil
  use mem_para,    only: myrank
  use sea_swm,     only: depthmax_swe

  implicit none

  integer :: isea, iwsfc

  real :: timefac_sst
  real :: timefac_seaice
  real :: dum1, dum2

  timefac_sst    = 0.0
  timefac_seaice = 0.0

  if (iupdsst == 1 .and. nsstfiles > 1) then
     timefac_sst = (s1900_sim           - s1900_sst(isstfile-1))  &
                 / (s1900_sst(isstfile) - s1900_sst(isstfile-1))
  endif

  if (iupdseaice == 1 .and. nseaicefiles > 1) then
     timefac_seaice = (s1900_sim                 - s1900_seaice(iseaicefile-1)) &
                    / (s1900_seaice(iseaicefile) - s1900_seaice(iseaicefile-1))
  endif

  !$omp parallel do private (iwsfc)
  do isea = 2,msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Initialize sea temperature, sea ice, and canopy depth

     sea%seatc(isea) = sea%seatp(isea)  &
                     + (sea%seatf(isea) - sea%seatp(isea)) * timefac_sst

     sea%seaicec(isea) = sea%seaicep(isea)  &
                       + (sea%seaicef(isea) - sea%seaicep(isea)) * timefac_seaice

     sfcg%can_depth(iwsfc) = 20. * max(1.,.025 * dt_sea)

     if (allocated(sea%spraytemp )) sea%spraytemp (isea) = sea%seatc(isea)
     if (allocated(sea%spray2temp)) sea%spray2temp(isea) = sea%seatc(isea)

  enddo
  !$omp end parallel do

  ! End of initialization that does not depend on atmospheric conditions

  if (runtype /= "INITIAL") return

  !$omp parallel do private (iwsfc,dum1,dum2)
  do isea = 2,msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply initial atmospheric properties to sea "canopy"

     sfcg%cantemp(iwsfc) = sfcg%airtheta(iwsfc) * sfcg%canexner(iwsfc)
     sfcg%canrrv (iwsfc) = sfcg%airrrv(iwsfc)
     sfcg%ustar  (iwsfc) = 0.1
     sfcg%ggaer  (iwsfc) = 0.0
     sfcg%wthv   (iwsfc) = 0.0
     sfcg%sfluxt (iwsfc) = 0.0
     sfcg%sfluxr (iwsfc) = 0.0

     ! Seawater quantities

     sea%sea_rough   (isea) = .001
     sea%sea_cantemp (isea) = sfcg%cantemp(iwsfc)
     sea%sea_canrrv  (isea) = sfcg%canrrv(iwsfc)
     sea%sea_ustar   (isea) = 0.1
     sea%sea_ggaer   (isea) = 0.0
     sea%sea_wthv    (isea) = 0.0
     sea%sea_sfluxt  (isea) = 0.0
     sea%sea_sfluxr  (isea) = 0.0
     sea%nlev_seaice (isea) = 0

     ! Seaice quantities

     sea%ice_sfluxt  (isea) = 0.0
     sea%ice_sfluxr  (isea) = 0.0
     sea%nlev_seaice (isea) = 0

     call prep_seaice(sea%seatc              (isea), &
                      sea%seaicec            (isea), &
                      sea%sea_cantemp        (isea), &
                      sea%ice_cantemp        (isea), &
                      sea%seaice_energy(1:nzi,isea), &
                      sea%seaice_tempk (1:nzi,isea), &
                      sea%nlev_seaice        (isea), &
                      sea%ice_albedo         (isea), &
                      sea%ice_rlongup        (isea), &
                      dum1                         , &
                      dum2                         , &
                      0.0                          , &
                      0.0                          , &
                      sea%ice_rough          (isea), &
                      sea%sea_canrrv         (isea), &
                      sea%ice_canrrv         (isea), &
                      sea%sea_ustar          (isea), &
                      sea%ice_ustar          (isea), &
                      sea%sea_ggaer          (isea), &
                      sea%ice_ggaer          (isea), &
                      sea%sea_wthv           (isea), &
                      sea%ice_wthv           (isea)  )

     if (sea%nlev_seaice(isea) > 0) then

        sfcg%rough    (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_rough  (isea) + &
                                       sea%seaicec(isea)  * sea%ice_rough  (isea)

        sfcg%cantemp  (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_cantemp(isea) + &
                                       sea%seaicec(isea)  * sea%ice_cantemp(isea)

        sfcg%canrrv   (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_canrrv (isea) + &
                                       sea%seaicec(isea)  * sea%ice_canrrv (isea)
     else

        sfcg%rough      (iwsfc) = sea%sea_rough   (isea)
        sfcg%cantemp    (iwsfc) = sea%sea_cantemp (isea)
        sfcg%canrrv     (iwsfc) = sea%sea_canrrv  (isea)

     endif

  enddo
  !$omp end parallel do

end subroutine sea_init_atm
