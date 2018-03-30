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

!!!!!
!!
!! THIS IS SMAGORINSKY-LILLY-HILL TURBULENCE PARAMETERIZATION
!!
!! TENDENCIES DUE TO VERTICAL MIXING (INCLUDING THE CONTRIBUTION OF SURFACE
!! FLUXES) ARE COMPUTED. HORIZONTAL EDDY DIFFUSIVITIES ARE ALSO COMPUTED
!! BUT APPLIED ELSWHERE
!!
!!!!!

  subroutine turb_k(iw, mrlw, thetav, vkh, vkm)

    use mem_turb,    only: vkm_sfc
    use mem_ijtabs,  only: itab_w
    use mem_grid,    only: mza, lpv, lpw, dzim, zm, volt, arw, lsw, volti
    use misc_coms,   only: io6, idiffk, csx, csz, dtlm
    use mem_basic,   only: rho, vxe, vye, vze, thil
    use consts_coms, only: vonk, grav2
    use mem_grid,    only: mza, lpw, arw0, dzm, dzim, zm
    use mem_tend,    only: vmxet, vmyet, vmzet, thilt
    use var_tables,  only: num_scalar, scalar_tab
    use tridiag,     only: tridv
    use oname_coms,  only: nl

    use supercell_testm, only: vxe_init, vye_init, vze_init, thil_init, sh_w_init

    implicit none

    integer, intent(in)    :: iw, mrlw
    real,    intent(in)    :: thetav(mza)
    real,    intent(inout) :: vkh(mza)
    real,    intent(inout) :: vkm(mza)

    integer :: npoly, jw1, jw2, iw1, iw2, iv1, iv2, n, ns, k, ka, ks

    real :: richnum,ambda,ambda2,hill_term,richnum_term
    real :: scalen_asympt,scalen_vert,scalen_horiz
    real :: prinv
    real :: strain2(mza), strain2h(mza)
    real :: bvfreq2(mza)
    real :: vkz2(mza)
    real :: dzim2 (mza)
    real :: dtl, dtli
    real :: gxps1, gyps1, gxps2, gyps2

    real :: akodz(mza), dtomass(mza)
    real :: vctr3(mza),vctr5(mza),vctr6(mza),vctr7(mza)

    real :: rhs(mza,max(3,num_scalar+1)), soln(mza,max(3,num_scalar+1))
    real :: del_scp(mza), varp(mza)
    real :: vctr2(mza,3)

    real :: gxps_vxe(mza), gxps_vye(mza), gxps_vze(mza)
    real :: gyps_vxe(mza), gyps_vye(mza), gyps_vze(mza)

    real, parameter :: rchmax =    3.   ! Test with new asympt vert scale length
    real, parameter :: rmin   = -100.
    real, parameter :: rmax   = 1. / 3. ! Critical Ri number
    real, parameter :: rmaxi  = 1. / rmax

    real, parameter :: pri0   = 1.43    ! Inverse Prandtl number at neutral
    real, parameter :: fact   = (pri0 - 1.0) / rmax

!!!!!
!!
!!  FIRST COMPUTE THE EDDY DIFFUSIVITIES, AND SAVE
!!  THE HORIZONTAL DIFFUSIVITIES HKM FOR USE ELSEWHERE
!!
!!!!!

    do k = 2,mza
       dzim2(k) = dzim(k) * dzim(k)

       strain2h(k) = 0.
    enddo

    ka = lpw(iw)

    scalen_horiz = csx(mrlw) * sqrt(arw0(iw))  ! change this later?
    scalen_asympt = csx(mrlw) * 300.

    if (idiffk(mrlw) == 3) then

    ! Compute horizontal gradients of vxe, vye, vze for full 3D strain rate

       npoly = itab_w(iw)%npoly
   
       gxps_vxe(:) = 0.
       gyps_vxe(:) = 0.

       gxps_vye(:) = 0.
       gyps_vye(:) = 0.

       gxps_vze(:) = 0.
       gyps_vze(:) = 0.

! Loop over W neighbors of this W cell

       do jw1 = 1, npoly

          jw2 = mod(jw1,npoly) + 1

          iw1 = itab_w(iw)%iw(jw1)
          iw2 = itab_w(iw)%iw(jw2)

          iv1 = itab_w(iw)%iv(jw1)
          iv2 = itab_w(iw)%iv(jw2)

          gxps1 = itab_w(iw)%gxps1(jw1)
          gyps1 = itab_w(iw)%gyps1(jw1)

          gxps2 = itab_w(iw)%gxps2(jw1)
          gyps2 = itab_w(iw)%gyps2(jw1)

! Vertical loop over T levels
! Zero-gradient lateral B.C. below lpv(iv1)

          do k = lpv(iv1), mza
             gxps_vxe(k) = gxps_vxe(k) + gxps1 * (vxe(k,iw1) - vxe(k,iw))
             gyps_vxe(k) = gyps_vxe(k) + gyps1 * (vxe(k,iw1) - vxe(k,iw))

             gxps_vye(k) = gxps_vye(k) + gxps1 * (vye(k,iw1) - vye(k,iw))
             gyps_vye(k) = gyps_vye(k) + gyps1 * (vye(k,iw1) - vye(k,iw))

             gxps_vze(k) = gxps_vze(k) + gxps1 * (vze(k,iw1) - vze(k,iw))
             gyps_vze(k) = gyps_vze(k) + gyps1 * (vze(k,iw1) - vze(k,iw))
          enddo

! Vertical loop over T levels
! Zero-gradient lateral B.C. below lpv(iv2)

          do k = lpv(iv2), mza
             gxps_vxe(k) = gxps_vxe(k) + gxps2 * (vxe(k,iw2) - vxe(k,iw))
             gyps_vxe(k) = gyps_vxe(k) + gyps2 * (vxe(k,iw2) - vxe(k,iw))

             gxps_vye(k) = gxps_vye(k) + gxps2 * (vye(k,iw2) - vye(k,iw))
             gyps_vye(k) = gyps_vye(k) + gyps2 * (vye(k,iw2) - vye(k,iw))

             gxps_vze(k) = gxps_vze(k) + gxps2 * (vze(k,iw2) - vze(k,iw))
             gyps_vze(k) = gyps_vze(k) + gyps2 * (vze(k,iw2) - vze(k,iw))
          enddo

       enddo

    ! Loop over T levels: horizontal gradient contributions to 3D strain rate

       do k = ka,mza
          strain2h(k) = gxps_vxe(k)**2 + gyps_vxe(k)**2 &
                      + gxps_vye(k)**2 + gyps_vye(k)**2 &
                      + gxps_vze(k)**2 + gyps_vze(k)**2
       enddo

    endif ! (idiffk(mrlw) == 3)

    ! Loop over W levels

    do k = ka,mza-1

    ! Vertical strain rate squared (based on dV/dz only)

       strain2(k) = dzim2(k) * ((vxe(k+1,iw) - vxe(k,iw))**2  &
                             +  (vye(k+1,iw) - vye(k,iw))**2  &
                             +  (vze(k+1,iw) - vze(k,iw))**2)

       ! Add horizontal-gradient contribution to strain rate squared

       strain2(k) = strain2(k) + .5 * (strain2h(k) + strain2h(k+1))

       bvfreq2(k) = grav2 * dzim(k)  &
                  * (thetav(k+1) - thetav(k)) / (thetav(k+1) + thetav(k))

    enddo

    do k = ka,mza-1

       !  Compute Richardson number and Lilly Richardson-number term

       richnum = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-15)))
       richnum_term = min(rchmax,sqrt(max(0.,(1.-rmaxi*richnum))))

       ! Compute vertical and net scale lengths: scalen_vert & ambda

       scalen_vert = csz(mrlw) * dzm(k)
       ambda = max(scalen_vert,min(scalen_asympt,scalen_horiz))
       ambda2 = ambda ** 2

       vkz2(k) = (vonk * min( zm(k) - zm(ka-1), 1.e5)) ** 2

       if (bvfreq2(k) < -1.e-12) then
          hill_term = sqrt(-bvfreq2(k))
       else
          hill_term = 0.
       endif

       vkm(k) = .5 * (rho(k,iw) + rho(k+1,iw))        & ! density factor
              * vkz2(k) * ambda2 / (vkz2(k) + ambda2) & ! lengthscale^2 factor
              * (sqrt(strain2(k)) + hill_term)        & ! strain rate + Hill term
              * richnum_term                            ! Lilly Richnum term

       ! Subgrid Prandl number from Mason and Brown (1999, JAS)

       if (richnum >= 0.0) then
          prinv = pri0 - fact * richnum
       else
          prinv = pri0 * sqrt((1.-40.*richnum) / (1.-16.*richnum))
       endif

       vkh(k) = vkm(k) * prinv
      
    enddo
      
    ! Zero values for top and bottom boundaries   
   
    vkm(ka-1) = 0.
    vkh(ka-1) = 0.
    vkm(mza)  = 0.
    vkh(mza)  = 0.

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
       do k = ka,mza-1
          vkm(k) =  500. * .5 * real(rho(k,iw) + rho(k+1,iw))
          vkh(k) = 1500. * .5 * real(rho(k,iw) + rho(k+1,iw))
       enddo
    endif

!!!!!
!!
!! COMPUTE THE MOMENTUM TENDENCIES DUE TO VERTICAL MIXING
!!
!!!!!

    dtl  = dtlm(itab_w(iw)%mrlw)
    dtli = 1. / dtl

    ! Vertical loop over W levels

    do k = ka,mza-1
       akodz(k) = arw(k,iw) * vkm(k) * dzim(k)
    enddo

    akodz(ka-1)  = 0.
    akodz(mza) = 0.

    ! Distribution of surface flux over multiple levels
    
    vctr3(1:mza) = 0.

    if (nl%test_case /= 131) then ! No surface flux for supercell test case
       do k = ka, ka + lsw(iw) - 1
          ks = k - ka + 1
          vctr3(k) = (arw(k,iw) - arw(k-1,iw)) * vkm_sfc(ks,iw) * dzim(k-1) * 2.0
       enddo
    endif

    ! Vertical loop over T levels

    do k = ka,mza
       
       ! Mass in t control volume and its inverse times dtl
       
       dtomass(k) = dtl / ( rho(k,iw) * volt(k,iw) )

       ! Fill tri-diagonal matrix coefficients

       vctr5(k) = -dtomass(k) * akodz(k-1)
       vctr7(k) = -dtomass(k) * akodz(k)
       vctr6(k) = 1. - vctr5(k) - vctr7(k) + dtomass(k) * vctr3(k)

       ! Fill r.h.s. vectors

       if (nl%test_case /= 131) then  ! Standard form
          rhs(k,1) = vxe(k,iw)
          rhs(k,2) = vye(k,iw)
          rhs(k,3) = vze(k,iw)
       else                           ! Supercell test case only
          rhs(k,1) = vxe(k,iw) - vxe_init(k,iw)
          rhs(k,2) = vye(k,iw) - vye_init(k,iw)
          rhs(k,3) = vze(k,iw) - vze_init(k,iw)
       endif

    enddo

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

    ! Compute tendencies

    do k = ka, mza
       vmxet(k,iw) = vmxet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,1) - vctr2(k,1) - vctr3(k) * soln(k,1))
       vmyet(k,iw) = vmyet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,2) - vctr2(k,2) - vctr3(k) * soln(k,2))
       vmzet(k,iw) = vmzet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,3) - vctr2(k,3) - vctr3(k) * soln(k,3))
    enddo

!!!!!
!!
!! COMPUTE TEMPERATURE AND SCALAR TENDENCIES DUE TO VERTICAL MIXING
!! THE SOLVER COMPUTES THE FUTURE FLUXES, NOT THE FUTURE SCALAR VALUES
!!
!!!!!

    ! Vertical loop over W levels: fill tri-diagonal matrix coefficients and r.h.s.

    do k = ka, mza-1
       akodz(k) = arw(k,iw) * vkh(k) * dzim(k)
       vctr5(k) = -akodz(k) * dtomass(k)
       vctr7(k) = -akodz(k) * dtomass(k+1)
       vctr6(k) = 1. - vctr5(k) - vctr7(k)
    enddo

    do n = 1, num_scalar
       do k = ka, mza
          varp(k) = scalar_tab(n)%var_p(k,iw) + dtl * scalar_tab(n)%var_t(k,iw) / rho(k,iw)
       enddo
       do k = ka, mza-1
          rhs(k,n) = akodz(k) * (varp(k) - varp(k+1))

          ! For supercell test case, for sh_w, subtract gradient of initial field

          if (nl%test_case == 131 .and. scalar_tab(n)%name == 'SH_W') then
             rhs(k,n) = rhs(k,n) - akodz(k) * (sh_w_init(k,iw) - sh_w_init(k+1,iw))
          endif
       enddo
    enddo

    n = num_scalar + 1
    do k = ka, mza
       varp(k) = thil(k,iw) + dtl * thilt(k,iw) / rho(k,iw)
    enddo
    do k = ka, mza-1
       rhs(k,n) = akodz(k) * (varp(k) - varp(k+1))

       ! For supercell test case, for thil, subtract gradient of initial field

       if (nl%test_case == 131) then
          rhs(k,n) = rhs(k,n) - akodz(k) * (thil_init(k,iw) - thil_init(k+1,iw))
       endif
    enddo

    ! Solve tri-diagonal matrix equation

    if (ka <= mza-1) then
       call tridv( vctr5, vctr6, vctr7, rhs, soln, ka, mza-1, mza, num_scalar+1 )
    endif

    ! Set bottom and top vertical internal turbulent fluxes to zero

    soln( ka-1, 1:num_scalar+1) = 0.0
    soln(mza  , 1:num_scalar+1) = 0.0

    ! Vertical loop over T levels

    do n = 1, num_scalar
       do k = ka, mza
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + volti(k,iw) * (soln(k-1,n) - soln(k,n))
       enddo
    enddo

    n = num_scalar + 1
    do k = ka, mza
       thilt(k,iw) = thilt(k,iw) + volti(k,iw) * (soln(k-1,n) - soln(k,n))
    enddo

  end subroutine turb_k

end module smagorinsky
