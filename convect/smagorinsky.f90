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

    use mem_turb,    only: vkm_sfc, vkm, vkh, wtv0_k, ustar_k, frac_sfc
    use mem_ijtabs,  only: itab_w
    use mem_grid,    only: mza, lpw, dzim, dzit, zm, arw, lsw, volti, zm, &
                           arw0, dzm, dzimsq, dzt_bot
    use misc_coms,   only: idiffk, csx, csz, dtlm
    use mem_basic,   only: rho, vxe, vye, vze, thil, theta, tair, rr_w, rr_v
    use consts_coms, only: vonk, grav, grav2, eps_virt, cpio2
    use mem_tend,    only: vmxet, vmyet, vmzet, thilt
    use var_tables,  only: num_scalar, scalar_tab
    use tridiag,     only: tridv
    use oname_coms,  only: nl
    use buoyancy,    only: comp_buoy
    use grad_lib,    only: comp_vel_grads_ec
    use supercell_testm, only: vxe_init, vye_init, vze_init, thil_init, rr_w_init

    implicit none

    integer, intent(in) :: iw, mrlw

    integer :: n, k, ka, ks

    real :: richnum,ambda,ambda2,hill_term,richnum_term
    real :: scalen_asympt,scalen_vert,scalen_horiz
    real :: prinv
    real :: strain2(mza), strain_t(mza)
    real :: dissipation(mza) ! Dissipation rate (conversion of KE to heat) [kg/m^3 m^2/s^3]
    real :: bvfreq2
    real :: vkz2
    real :: dtl
    real :: vels2, usfc, ufree

    real :: zkm(mza), zkh(mza)
    real :: akodz(mza), dtomass(mza)
    real :: vctr3(mza),vctr5(mza),vctr6(mza),vctr7(mza)

    real :: rhs(mza,max(3,num_scalar+1)), soln(mza,max(3,num_scalar+1))
    real :: varp(mza)
    real :: vctr2(mza,3)
    real :: ql
    real :: dudx(mza), dudy(mza), dudz(mza)
    real :: dvdx(mza), dvdy(mza), dvdz(mza)
    real :: dwdx(mza), dwdy(mza), dwdz(mza)
    real :: buoy(mza), thetav(mza)

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

       do k = ka, mza-1
          strain2(k) = dzimsq(k) * ((vxe(k+1,iw) - vxe(k,iw))**2  &
                                 +  (vye(k+1,iw) - vye(k,iw))**2  &
                                 +  (vze(k+1,iw) - vze(k,iw))**2)
       enddo

    endif ! (idiffk(mrlw) == 3)

    ! Virtual potential temperature

    do k = ka, mza
       ql = rr_w(k,iw) - rr_v(k,iw)
       thetav(k) = theta(k,iw) * (1.0 + eps_virt * rr_v(k,iw) - ql)
    enddo

    ! Include surface shear at first level?

    do ks = 1, lsw(iw)
       k  = ks + lpw(iw) - 1

       ufree = (grav * max(wtv0_k(ks,iw),0.) * dzt_bot(k) / thetav(k)) ** onethird
       usfc  = max(ustar_k(ks,iw), ufree)

       strain2(k) = (1. - frac_sfc(ks,iw)) * strain2(k) &
                  +       frac_sfc(ks,iw)  * max(strain2(k), (usfc * dzit(k))**2)
    enddo

    ! Compute buoyancy variable

    call comp_buoy(iw, buoy=buoy)

    ! Loop over W levels: Compute diffusivities for momentum and scalars

    do k = ka, mza-1

       ! Compute buoyancy frequency

       bvfreq2 = grav2 * buoy(k) / (thetav(k+1) + thetav(k))

       !  Compute Richardson number and Lilly Richardson-number term

       richnum = max(rmin, min(rmax, bvfreq2 / max(strain2(k),1.e-15)))
       richnum_term = min(rchmax, sqrt(max(0.,(1.-rmaxi*richnum))))

       ! Compute vertical and net scale lengths: scalen_vert & ambda

       scalen_vert = csz(mrlw) * dzm(k)

!      ambda = max(scalen_vert,min(scalen_asympt,scalen_horiz))
       ambda = (scalen_vert * min(scalen_asympt,scalen_horiz)**2)**onethird

       ambda2 = ambda ** 2

       vkz2 = (vonk * min( zm(k) - zm(ka-1), 1.e5)) ** 2

       if (bvfreq2 < -1.e-12) then
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
    dissipation(ka-1) = vkm_sfc(1,iw) * vels2 / dzt_bot(ka)**2
    dissipation(mza)  = 0.0

    ! Zero values for top and bottom boundaries

    zkm(ka-1) = 0.
    zkh(ka-1) = 0.
    zkm(mza)  = 0.
    zkh(mza)  = 0.

    ! For DCMIP 2015 baroclinic and tropical cyclone tests, return here so that
    ! vertical mixing is not done by standard OLAM code; instead it will
    ! (optionally) be done in dcmip_physics routine.

    if (nl%test_case == 110 .or. &
        nl%test_case == 111 .or. &
        nl%test_case == 112 .or. &
        nl%test_case == 113 .or. &
        nl%test_case == 114 .or. &
        nl%test_case == 121 .or. &
        nl%test_case == 122) return

    ! For DCMIP 2016 supercell test case, specify spatially constant K's

    if (nl%test_case == 131) then
       do k = ka, mza-1
          zkm(k) =  500. * .5 * real(rho(k,iw) + rho(k+1,iw))
          zkh(k) = 1500. * .5 * real(rho(k,iw) + rho(k+1,iw))
       enddo
    endif

    ! Vertical loop over W levels

    do k = ka,mza-1
       akodz(k) = arw(k,iw) * (zkm(k) + vkm(k,iw)) * dzim(k) * 0.5
    enddo

    akodz(ka-1) = 0.
    akodz(mza)  = 0.

    ! Distribution of surface drag over multiple levels

    vctr3 = 0.

    if (nl%test_case /= 131) then
       do k = ka, ka + lsw(iw) - 1
          ks = k - ka + 1
          vctr3(k) = 2.0 * dzim(k-1) * (arw(k,iw) - arw(k-1,iw)) * vkm_sfc(ks,iw)
       enddo
    endif

    ! COMPUTE MOMENTUM TENDENCIES DUE TO VERTICAL MIXING

    ! Vertical loop over T levels
    do k = ka, mza

       ! Mass in t control volume and its inverse times dtl

       dtomass(k) = dtl * volti(k,iw) / real(rho(k,iw))

       ! Fill tri-diagonal matrix coefficients

       vctr5(k) = -dtomass(k) * akodz(k-1)
       vctr7(k) = -dtomass(k) * akodz(k)
       vctr6(k) = 1. - vctr5(k) - vctr7(k) + dtomass(k) * vctr3(k)

       ! Fill r.h.s. vectors

       rhs(k,1) = vxe(k,iw)
       rhs(k,2) = vye(k,iw)
       rhs(k,3) = vze(k,iw)

    enddo

    ! Supercell test case only

    if (nl%test_case == 131) then
       do k = ka, mza
          rhs(k,1) = rhs(k,1) - vxe_init(k,iw)
          rhs(k,2) = rhs(k,2) - vye_init(k,iw)
          rhs(k,3) = rhs(k,3) - vze_init(k,iw)
       enddo
    endif

    ! Solve tri-diagonal matrix for each component

    if (ka <= mza) then
       call tridv( vctr5, vctr6, vctr7, rhs, soln, ka, mza, mza, 3 )
    endif

    ! Now, soln contains velocity(t+1) values

    ! Compute internal vertical turbulent fluxes

    do k = ka, mza-1
       vctr2(k,1) = akodz(k) * (soln(k,1) - soln(k+1,1))
       vctr2(k,2) = akodz(k) * (soln(k,2) - soln(k+1,2))
       vctr2(k,3) = akodz(k) * (soln(k,3) - soln(k+1,3))
    enddo

    ! Set bottom and top internal fluxes to zero

    vctr2(ka-1, 1:3) = 0.
    vctr2(mza , 1:3) = 0.

    ! Compute momentum tendencies

    do k = ka, mza
       vmxet(k,iw) = vmxet(k,iw) + volti(k,iw)  &
                   * (vctr2(k-1,1) - vctr2(k,1) - vctr3(k) * soln(k,1))

       vmyet(k,iw) = vmyet(k,iw) + volti(k,iw)  &
                   * (vctr2(k-1,2) - vctr2(k,2) - vctr3(k) * soln(k,2))

       vmzet(k,iw) = vmzet(k,iw) + volti(k,iw)  &
                   * (vctr2(k-1,3) - vctr2(k,3) - vctr3(k) * soln(k,3))
    enddo

    ! COMPUTE TEMPERATURE AND SCALAR TENDENCIES DUE TO VERTICAL MIXING

    ! Vertical loop over W levels: fill tri-diagonal matrix coefficients and r.h.s.

    do k = ka, mza-1
       akodz(k) = arw(k,iw) * (zkh(k) + vkh(k,iw)) * dzim(k) * 0.5
       vctr5(k) = -akodz(k) * dtomass(k)
       vctr7(k) = -akodz(k) * dtomass(k+1)
       vctr6(k) = 1. - vctr5(k) - vctr7(k)
    enddo

    ! Load scalars into tridiagonal arrays

    do n = 1, num_scalar

       do k = ka, mza
          varp(k) = scalar_tab(n)%var_p(k,iw) + dtl * scalar_tab(n)%var_t(k,iw) / real(rho(k,iw))
       enddo

       do k = ka, mza-1
          rhs(k,n) = akodz(k) * (varp(k) - varp(k+1))
       enddo

       ! For supercell test case, for rr_w, subtract gradient of initial field

       if (nl%test_case == 131) then
          if (scalar_tab(n)%name == 'RR_W') then
             do k = ka, mza-1
                rhs(k,n) = rhs(k,n) - akodz(k) * (rr_w_init(k,iw) - rr_w_init(k+1,iw))
             enddo
          endif
       endif
    enddo

    ! Load temperature into tridiagional arrays

    n = num_scalar + 1

    do k = ka, mza
       varp(k) = thil(k,iw) + dtl * thilt(k,iw) / real(rho(k,iw))
    enddo

    do k = ka, mza-1
       rhs(k,n) = akodz(k) * (varp(k) - varp(k+1))
    enddo

    ! For supercell test case, for thil, subtract gradient of initial field

    if (nl%test_case == 131) then
       do k = ka, mza-1
          rhs(k,n) = rhs(k,n) - akodz(k) * (thil_init(k,iw) - thil_init(k+1,iw))
       enddo
    endif

    ! Solve tri-diagonal matrix equation

    if (ka <= mza-1) then
       call tridv( vctr5, vctr6, vctr7, rhs, soln, ka, mza-1, mza, num_scalar+1 )
    endif

    ! Set bottom and top vertical internal turbulent fluxes to zero

    soln(ka-1, 1:num_scalar+1) = 0.0
    soln(mza , 1:num_scalar+1) = 0.0

    ! Apply fluxes to conserved scalars

    do n = 1, num_scalar
       ! Vertical loop over T levels
       do k = ka, mza
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + volti(k,iw) * (soln(k-1,n) - soln(k,n))
       enddo
    enddo

    ! Apply fluxes and dissipative heating to THIL

    n = num_scalar + 1

    ! Vertical loop over T levels
    do k = ka, mza
       thilt(k,iw) = thilt(k,iw) + volti(k,iw) * (soln(k-1,n) - soln(k,n)) &
                   + (dissipation(k) + dissipation(k-1)) * cpio2 * theta(k,iw) / tair(k,iw)
    enddo

    vkm(ka-1:mza,iw) = zkm(ka-1:mza)
    vkh(ka-1:mza,iw) = zkh(ka-1:mza)

  end subroutine turb_k

end module smagorinsky
