module umwm_oforcing

  implicit none
  real, parameter :: alog10 = log(10.)

contains

!===============================================================================

subroutine oforcing( i, iwsfc )

  use mem_sea,     only: sea
  use mem_sfcg,    only: sfcg, nswmzons
  use consts_coms, only: erad, r8, vonk
  use umwm_module, only: umwm, ucurr, vcurr, gustiness

  implicit none

  integer, intent(in) :: i, iwsfc
  real                :: gustu, gustv, raxis

  !Bob: STD form of umwm%wspd (not first going to uw,vw); does not accept gustiness modification
  ! umwm%wspd(i) = sqrt(sea%windxe(i)**2 + sea%windye(i)**2 + sea%windze(i)**2)

  raxis = sqrt(sfcg%xew(iwsfc)**2 + sfcg%yew(iwsfc)**2)  ! dist from earth axis

  if (raxis > 1.e3) then

     ! Wind components
     umwm%uwind(i) = (sea%windye(i) * sfcg%xew(iwsfc) - sea%windxe(i) * sfcg%yew(iwsfc)) / raxis
     umwm%vwind(i) =  sea%windze(i) * raxis / erad &
          - (sea%windxe(i) * sfcg%xew(iwsfc) + sea%windye(i) * sfcg%yew(iwsfc)) * sfcg%zew(iwsfc) / (raxis * erad)

     if (nswmzons > 0) then
        ! Ocean current components
        ucurr(i) = (sea%vye(i) * sfcg%xew(iwsfc) - sea%vxe(i) * sfcg%yew(iwsfc)) / raxis
        vcurr(i) =  sea%vze(i) * raxis / erad &
                 - (sea%vxe(i) * sfcg%xew(iwsfc) + sea%vye(i) * sfcg%yew(iwsfc)) * sfcg%zew(iwsfc) / (raxis * erad)
     else
        ucurr(i) = 0.
        vcurr(i) = 0.
     endif

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
  umwm%wspd(i) = max(umwm%wspd(i), 1.e-1)
  umwm%wdir(i) = atan2(umwm%vwind(i),umwm%uwind(i))

  if (sfcg%dzt_bot(iwsfc) <= 10.) then

     umwm%wspd10m(i) = umwm%wspd(i)

  else

     if (umwm%iactive(i)) then
        umwm%wspd10m(i) = umwm%wspd(i) * (alog10 - umwm%alogzo(i)) / (umwm%alogzs(i) - umwm%alogzo(i))
     else
        call init_u10m(umwm%wspd(i), umwm%alogzs(i), umwm%wspd10m(i), umwm%alogzo(i))
     endif

  endif

  umwm%ustar(i) = umwm%wspd(i) * vonk / (umwm%alogzs(i) - umwm%alogzo(i))

end subroutine oforcing

!===============================================================================

subroutine zo_andreas(w10m, alogzo)

  use consts_coms, only: vonk
  implicit none

  real, intent(in)  :: w10m   ! 10m wind speed
  real, intent(out) :: alogzo ! ln( z0 )
  real              :: ww, xx

  ww = w10m - 8.271
  xx = 0.239 + 0.0433 * ww + sqrt( .339e-3 + .225e-3 * ww**2 )

  alogzo = alog10 - vonk * w10m / xx

end subroutine zo_andreas

!===============================================================================

subroutine init_u10m(wspd, alogzbot, u10m, alogzo)

  implicit none

  real, intent(in)  :: wspd     ! wind speed observation
  real, intent(in)  :: alogzbot ! log(zbot)
  real, intent(out) :: u10m     ! wind speed at 10m
  real, intent(out) :: alogzo   ! log(z0)
  integer           :: i

  u10m = wspd
  do i = 1, 5
     call zo_andreas(u10m, alogzo)
     u10m = wspd * (alog10 - alogzo) / (alogzbot - alogzo)
  enddo

end subroutine init_u10m

!===============================================================================

end module umwm_oforcing
