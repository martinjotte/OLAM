subroutine lake_init_atm()

  use mem_lake,    only: lake, mlake, omlake
  use lake_coms,   only: dt_lake
  use misc_coms,   only: iparallel, runtype
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_para,    only: myrank
  use consts_coms, only: t00, cliq, alli, p00i, grav, rocp
  use therm_lib,   only: rhovsl, rhovsil

  implicit none

  integer :: ilake, iwsfc
  real    :: laketemp

  ! Initialize lake quantities that do not depend on atmospheric conditions

  !$omp parallel do private (iwsfc)
  do ilake = 2,mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Initialize lake depth and canopy depth

     lake%depth(ilake) = sfcg%topw(iwsfc) - sfcg%bathym(iwsfc)
     sfcg%can_depth(iwsfc) = 20. * max(1.,.025 * dt_lake)

     sfcg%head1(iwsfc) = lake%depth(ilake) + sfcg%bathym(iwsfc) - sfcg%topw(iwsfc)

  enddo  ! ilake
  !$omp end parallel do

  ! End of initialization that does not depend on atmospheric conditions

  if (runtype /= "INITIAL") return

  ! Sea quantities that are initialized only on 'INITIAL' run

  !$omp parallel do private (iwsfc,laketemp)
  do ilake = 2,mlake
     iwsfc = ilake + omlake

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! Apply initial atmospheric properties to lake "canopy"

     sfcg%cantemp(iwsfc) = sfcg%airtheta(iwsfc) * (sfcg%prss(iwsfc) * p00i)**rocp
     sfcg%canrrv (iwsfc) = sfcg%airrrv(iwsfc)
     sfcg%ustar  (iwsfc) = 0.1
     sfcg%ggaer  (iwsfc) = 0.0
     sfcg%wthv   (iwsfc) = 0.0

     ! Lake water quantities (eventually enable lake_energy to be initialized from spin-up simulation)

     laketemp                 = sfcg%cantemp(iwsfc)
     lake%lake_energy (ilake) = (laketemp - t00) * cliq + alli
     lake%surface_srrv(ilake) = rhovsl(laketemp - t00) / sfcg%rhos(iwsfc)
     sfcg%rough       (iwsfc) = .001
     sfcg%ustar       (iwsfc) = 0.1
     sfcg%ggaer       (iwsfc) = 0.0
     sfcg%wthv        (iwsfc) = 0.0

  enddo
  !$omp end parallel do

end subroutine lake_init_atm
