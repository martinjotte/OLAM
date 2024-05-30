module umwm_oforcing

  use umwm_module

  implicit none

contains

!===============================================================================

subroutine oforcing()

  use mem_sea,        only: sea, msea, omsea
  use mem_sfcg,       only: sfcg
  use consts_coms,    only: erad

  integer :: i, iwsfc
  real :: gustu, gustv, raxis

  !$omp parallel do private(iwsfc, raxis, gustu, gustv)
  do i = 2,msea
     iwsfc = i + omsea

     !Bob: STD form of umwm%wspd (not first going to uw,vw); does not accept gustiness modification
     ! umwm%wspd(i) = sqrt(sea%windxe(i)**2 + sea%windye(i)**2 + sea%windze(i)**2)

     raxis = sqrt(sfcg%xew(iwsfc)**2 + sfcg%yew(iwsfc)**2)  ! dist from earth axis

     if (raxis > 1.e3) then

        ! Wind components
        umwm%uwind(i) = (sea%windye(i) * sfcg%xew(iwsfc) - sea%windxe(i) * sfcg%yew(iwsfc)) / raxis
        umwm%vwind(i) =  sea%windze(i) * raxis / erad &
              - (sea%windxe(i) * sfcg%xew(iwsfc) + sea%windye(i) * sfcg%yew(iwsfc)) * sfcg%zew(iwsfc) / (raxis * erad)

        ! Ocean current components
        ucurr(i) = (sea%vye   (i) * sfcg%xew(iwsfc) - sea%vxe   (i) * sfcg%yew(iwsfc)) / raxis
        vcurr(i) =  sea%vze   (i) * raxis / erad &
                 - (sea%vxe   (i) * sfcg%xew(iwsfc) + sea%vye   (i) * sfcg%yew(iwsfc)) * sfcg%zew(iwsfc) / (raxis * erad)

     else
        umwm%uwind(i) = 0.
        umwm%vwind(i) = 0.

        ucurr(i) = 0.
        vcurr(i) = 0.
     endif

     ! add wind gustiness if requested in the namelist;
     ! uw, vw will get a uniformly distributed gust component
     if (gustiness > 0.) then

        ! get random numbers [0, 1]
        call random_number(gustu)
        call random_number(gustv)

        gustu = gustiness * (2. * gustu - 1.)
        gustv = gustiness * (2. * gustv - 1.)

        ! add gustiness to the wind fields
        umwm%uwind(i) = umwm%uwind(i) * (1. + gustu)
        umwm%vwind(i) = umwm%vwind(i) * (1. + gustv)

     endif

     umwm%wspd(i) = sqrt(umwm%uwind(i)**2 + umwm%vwind(i)**2)
     umwm%wdir(i) = atan2(umwm%vwind(i),umwm%uwind(i))

  enddo
  !$omp end parallel do

end subroutine oforcing

end module umwm_oforcing
