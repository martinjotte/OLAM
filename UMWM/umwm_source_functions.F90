module umwm_source_functions
  use umwm_module
  use umwm_sheltering, only: sheltering_coare35, sheltering_reynolds
  use mem_sea, only: sea, msea, omsea
  use mem_sfcg, only: sfcg

  implicit none

  private

  public :: sin_d12, sds_d12, snl_d12, s_ice

contains

!===============================================================================

  subroutine sin_d12()

    use consts_coms, only: pi2, grav

    ! Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).
    integer :: i, o, p, iwsfc
    real    :: c6, oneovsin_fac, c7

    c6 = pi2 * rhoswi 
    oneovsin_fac = 1. / sin_fac

    !$omp parallel do private(o,p,iwsfc,c7)
    do i = 2,msea

       umwm%wspd(i)    = max(umwm%wspd(i), 1.e-2)               ! protection against low wind speed values
       fcutoff(i) = min(fprog, 0.53 * grav / umwm%wspd(i)) ! cut-off frequency (4*pierson-moskowitz peak frequency)
       shelt(i)   = sheltering_coare35(umwm%wspd(i))       ! compute variable sheltering coefficient

       ! search for the cut-off frequency bin:
       do o = om-2, 1, -1
          oc(i) = o
          !oc(i) = om-2
          if (fcutoff(i) > umwm%f(o)) exit
       enddo

       do p = 1,pm
          do o = 1,om
             ! compute wind input at height of half wavelength:
             ssin(o,p,i) = (umwm%wspd(i) + 2.5 * umwm%ustar(i) * logl2overz(o,i))&
                         * cos(umwm%wdir(i) - th(p)) - cp0(o,i)&
                        - ucurr(i) * cth(p) - vcurr(i) * sth(p)

             ssin(o,p,i) = sin_fac * abs(ssin(o,p,i)) * ssin(o,p,i)
          enddo
       enddo

       iwsfc = i + omsea
       c7 = c6 * (1 - sea%seaicec(i)) * sfcg%rhos(iwsfc)
       do p = 1,pm
          do o = 1,om
             ! apply variable sheltering coefficient
             if (ssin(o,p,i) > 0.) then
                ssin(o,p,i) = ssin(o,p,i) * shelt(i) * oneovsin_fac
             else
                ! adjust input for opposing winds
                ssin(o,p,i) = ssin(o,p,i) * fieldscale1
                if (cos(umwm%wdir(i)-th(p)) > 0) then
                   ! further reduce for swell that overruns the wind
                   ssin(o,p,i) = ssin(o,p,i) * fieldscale2
                endif
             endif
             ssin(o,p,i) = c7 * ssin(o,p,i) * fkovg(o,i)
          enddo
       enddo

       do p = 1,pm
          do o = oc(i)+1,om
             ! prevent negative ssin for diagnostic tail
             ssin(o,p,i) = max(ssin(o,p,i), 0.)
          enddo
       enddo

    enddo
    !$omp end parallel do

  end subroutine sin_d12

!===============================================================================

  subroutine sds_d12()

    ! Wave dissipation function described by Donelan et al. (2012).
    integer :: i, o, p

    if (mss_fac > 0) then

       !$omp parallel do private(p,o)
       do i = 2,msea
          dummy(:,:,i) = 0.
          do p = 1,pm
             do o = 2, om
                dummy(o,p,i) = dummy(o-1,p,i) + sum(evs(o-1,:,i) * cth2pp(:,p)) * k3dk(o-1,i)
             enddo
          enddo
          dummy(:,:,i) = (1. + mss_fac * dummy(:,:,i))**2

          do p = 1,pm
             do o = 1,om
                sds(o,p,i) = twopisds_fac * umwm%f(o) * dummy(o,p,i) * (evs(o,p,i) * k4(o,i))**sds_power
             enddo
          enddo

       enddo
       !$omp end parallel do

    else

       !$omp parallel do private(p,o)
       do i = 2,msea
          do p = 1,pm
             do o = 1,om
                sds(o,p,i) = twopisds_fac * umwm%f(o) * (evs(o,p,i) * k4(o,i))**sds_power
             enddo
          enddo
       enddo
       !$omp end parallel do

    endif

  end subroutine sds_d12

!===============================================================================

  subroutine snl_d12()

    integer :: o, p, i, iwsfc
    real    :: c6

    c6 = sdt_fac * sqrt(rhoswi)

    ! spread wave energy to 2 next longer wavenumbers exponentially decaying
    ! as distance from donating wavenumber, and remove the energy from
    ! donating wavenumbers:
    !$omp parallel do private(iwsfc,p,o)
    do i = 2,msea

       iwsfc = i + omsea
       snl(:,:,i) = 0.
       do concurrent(o = 1:oc(i), p = 1:pm)
          snl(o,p,i) = bf1_renorm(o,i) * sds(o+1,p,i) * evs(o+1,p,i)&
                     + bf2_renorm(o,i) * sds(o+2,p,i) * evs(o+2,p,i)&
                     - snl_fac * sds(o,p,i) * evs(o,p,i)
       enddo

       do concurrent(o = 1:om, p = 1:pm)
          ! account for plunging breakers
          sds(o,p,i) = sds(o,p,i) * cothkd(o,i)
       enddo

       do concurrent(o = 1:om)
          ! compute dissipation due to turbulence
          sdt(o,i) = c6 * sqrt(sfcg%rhos(iwsfc)) * umwm%ustar(i) * k(o,i)
       enddo

    enddo
    !$omp end parallel do

  end subroutine snl_d12

!===============================================================================

  subroutine s_ice()

    ! Wave attenuation by sea ice, following Kohout et al. (2014).

    integer :: i, o, p

    ! parameters from Kohout et al. 2014
    real, parameter :: H_th = 3.0       ! [m]
    real, parameter :: C1   = -5.35e-6  ! [m^-1]
    real, parameter :: C2   = C1 * H_th ! []

    real, dimension(om,pm) :: spectrumbin

    real :: ht_
  
    !$omp parallel do private(p,o,ht_)
    do i = 2,msea

       sice(:,i) = 0.0

       if (sea%seaicec(i) > fice_lth) then
 
          ht_ = 0.0
 
          do p = 1, pm
             do o = 1, om
                spectrumbin(o,p) = evs(o,p,i) * kdk(o,i)
 
                ht_ = ht_ + spectrumbin(o,p)
             enddo
          enddo

          ht_ = 4 * sqrt(ht_ * dth) ! significant wave height

          ! wave attenuation from sea ice in the two SWH regimes
          if (ht_ < H_th) then
             sice(:,i) = C1 * ht_
          else
             sice(:,i) = C2
          endif

          sice(:,i) = 2 * cg0(:,i) * sice(:,i)
         
       endif
 
    enddo
    !$omp end parallel do

  end subroutine s_ice

end module umwm_source_functions
