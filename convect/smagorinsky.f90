!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================

module smagorinsky

contains

  ! THIS IS SMAGORINSKY-LILLY-HILL TURBULENCE PARAMETERIZATION

  ! EDDY DIFFUSIVITIES AND TENDENCIES DUE TO VERTICAL MIXING
  ! (INCLUDING THE CONTRIBUTION OF SURFACE FLUXES) ARE COMPUTED.

  subroutine turb_k(iw, mrlw)

    use mem_turb,    only: vkm_sfc, vkm, vkh, kpblh
    use mem_ijtabs,  only: itab_w
    use mem_grid,    only: mza, lpw, zm, zm, arw0, dzm, dzimsq, dzit_bot
    use misc_coms,   only: idiffk, csx, csz, dtlm
    use mem_basic,   only: rho, vxe, vye, vze, theta, tair
    use consts_coms, only: vonk, grav, grav2, eps_virt, cpio2
    use mem_tend,    only: thilt
    use tridiag,     only: tridv
    use oname_coms,  only: nl
    use buoyancy,    only: comp_buoy
    use grad_lib,    only: comp_vel_grads_ec

    implicit none

    integer, intent(in) :: iw, mrlw

    integer :: k, ka

    real :: richnum,ambda,ambda2,hill_term,richnum_term
    real :: scalen_asympt,scalen_vert,scalen_horiz
    real :: prinv
    real :: strain2(mza), strain_t(mza)
    real :: dissipation(mza) ! Dissipation rate (conversion of KE to heat) [kg/m^3 m^2/s^3]
    real :: bvfreq2
    real :: vkz2
    real :: dtl
    real :: vels2
    real :: zkm(mza), zkh(mza)

    real :: dudx(mza), dudy(mza), dudz(mza)
    real :: dvdx(mza), dvdy(mza), dvdz(mza)
    real :: dwdx(mza), dwdy(mza), dwdz(mza)
    real :: buoy(mza)

    real, parameter :: asym_len =  400.
    real, parameter :: rchmax   =    3.   ! Test with new asympt vert scale length
    real, parameter :: rmin     = -100.
    real, parameter :: rmax     = 1. / 3. ! Critical Ri number
    real, parameter :: rmaxi    = 1. / rmax
    real, parameter :: pri0     = 1.43    ! Inverse Prandtl number at neutral
    real, parameter :: fact     = (pri0 - 1.0) / rmax
    real, parameter :: onethird = 1. / 3.

    ka  = lpw(iw)
    dtl = dtlm(itab_w(iw)%mrlw)

    scalen_horiz  = csx(mrlw) * sqrt(arw0(iw))  ! change this later?
    scalen_asympt = csx(mrlw) * asym_len

    if (idiffk(mrlw) == 3) then

       ! Compute earth-cartesian velocity gradients at cell centers

       call comp_vel_grads_ec(iw, dudx, dudy, dudz, &
                                  dvdx, dvdy, dvdz, &
                                  dwdx, dwdy, dwdz  )

       ! Strain rate squared at cell centers (T levels)

       do k = ka, mza
          strain_t(k) = 2.0 * (dudx(k)**2 + dvdy(k)**2 + dwdz(k)**2) &
                      + (dudy(k) + dvdx(k))**2 &
                      + (dudz(k) + dwdx(k))**2 &
                      + (dvdz(k) + dwdy(k))**2
       enddo

       ! Strain rate squared at W levels

       do k = ka, mza-1
          strain2(k) = 0.5 * (strain_t(k) + strain_t(k+1))
       enddo

    else

       ! Vertical strain rate squared (based on dV/dz only)

!       strain2(ka:mza-1) = shearh(ka:mza-1)

       do k = ka, mza-1
          strain2(k) = dzimsq(k) * ((vxe(k+1,iw) - vxe(k,iw))**2  &
                                 +  (vye(k+1,iw) - vye(k,iw))**2  &
                                 +  (vze(k+1,iw) - vze(k,iw))**2)
       enddo

    endif ! (idiffk(mrlw) == 3)

    ! Compute buoyancy variable

    call comp_buoy(iw, buoy=buoy)

    ! Loop over W levels: Compute diffusivities for momentum and scalars

    do k = ka, mza-1

       ! Compute buoyancy frequency

       bvfreq2 = grav2 * buoy(k) / (theta(k+1,iw) + theta(k,iw))

       !  Compute Richardson number and Lilly Richardson-number term

       richnum = max(rmin, min(rmax, bvfreq2 / max(strain2(k),1.e-15)))
       richnum_term = min(rchmax, sqrt(max(0.,(1.-rmaxi*richnum))))

       ! Compute vertical and net scale lengths: scalen_vert & ambda

       scalen_vert = csz(mrlw) * dzm(k)

!      ambda = max(scalen_vert,min(scalen_asympt,scalen_horiz))
       ambda = (scalen_vert * min(scalen_asympt,scalen_horiz)**2)**onethird

       ambda2 = ambda ** 2

       vkz2 = (vonk * min( zm(k) - zm(ka-1), 1.e5)) ** 2

       if (bvfreq2 < -1.e-12 .and. k <= kpblh(iw)) then
          hill_term = sqrt(-bvfreq2)
       else
          hill_term = 0.
       endif

       ! Eddy diffusivity for momentum

       zkm(k) = .5 * real(rho(k,iw) + rho(k+1,iw))  & ! density factor
              * vkz2 * ambda2 / (vkz2 + ambda2)     & ! lengthscale^2 factor
              * (sqrt(strain2(k)) + hill_term)      & ! strain rate + Hill term
              * richnum_term                          ! Lilly Richnum term

       ! Subgrid Prandl number from Mason and Brown (1999, JAS)

       if (richnum >= 0.0) then
          prinv = pri0 - fact * richnum
       else
          prinv = pri0 * sqrt((1.-40.*richnum) / (1.-16.*richnum))
       endif

       ! Eddy diffusivity for heat/scalars

       zkh(k) = zkm(k) * prinv

       ! Dissipation due to mechanical strain (applies at W points)

       dissipation(k) = strain2(k) * zkm(k)

    enddo

    ! Compute surface dissipation

    vels2             = vxe(ka,iw)**2 + vye(ka,iw)**2 + vze(ka,iw)**2
    dissipation(ka-1) = vkm_sfc(iw) * vels2 * dzit_bot(ka)**2
    dissipation(mza)  = 0.0

    ! For DCMIP 2016 supercell test case, specify spatially constant K's

    if (nl%test_case == 131) then
       do k = ka, mza-1
          zkm(k) =  500. * .5 * real(rho(k,iw) + rho(k+1,iw))
          zkh(k) = 1500. * .5 * real(rho(k,iw) + rho(k+1,iw))
       enddo
    endif

    ! Apply dissipative heating to THIL
    do k = ka, mza
       thilt(k,iw) = thilt(k,iw) + (dissipation(k) + dissipation(k-1)) * cpio2 * theta(k,iw) / tair(k,iw)
    enddo

    vkm(ka:mza-1,iw) = zkm(ka:mza-1)
    vkh(ka:mza-1,iw) = zkh(ka:mza-1)

  end subroutine turb_k

end module smagorinsky
