!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
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

  subroutine turb_k(iw, mrlw, rhot)

    use mem_turb,    only: hkm, pblh, kpblh, sflux_t, sflux_r, ustar, vkm_sfc, &
                           sxfer_tk, sxfer_rk, fthpbl, fqtpbl
    use mem_ijtabs,  only: istp, jtab_w, itab_w, mrl_begl
    use mem_grid,    only: mza, mwa, lpw, dzim, zt, zm, volt, arw, lsw, volti
    use misc_coms,   only: io6, idiffk, csx, csz, zkhkm, akmin, meshtype, &
                           dtlm
    use mem_basic,   only: rho, thil, theta, sh_v, vxe, vye, vze, sh_w
    use consts_coms, only: vonk, alvl, grav, grav2, cp, rvap
    use mem_grid,    only: mza, lpw, arw0, dzm, dzim, zm, lpu, lpv,  &
                           unx, uny, unz, vnx, vny, vnz
    use micro_coms,  only: level
    use mem_tend,    only: vmxet, vmyet, vmzet, thilt, sh_wt
    use var_tables,  only: num_scalar, scalar_tab

 !$ use omp_lib

    implicit none

    integer, intent(in)    :: iw, mrlw
    real,    intent(inout) :: rhot(mza,mwa)

    integer :: j,k,ka,n,ks

    real :: richnum,ambda,ambda2,hill_term,richnum_term,sbf
    real :: scalen_asympt,scalen_vert,scalen_horiz,bkmin

    real :: thetav(mza)
    real :: strain2(mza)
    real :: bvfreq2(mza)
    real :: vkz2(mza)
    real :: dzim2 (mza)
    real :: vkh(mza)
    real :: vkm(mza)

    real :: dtl, dtl2, dtli

    real :: s1, s2

    ! Automatic arrays:

    real :: akodz(mza), dtomass(mza)
    real :: vctr3(mza),vctr5(mza),vctr6(mza),vctr7(mza)
    real :: vctr8a(mza),vctr8b(mza),vctr8c(mza)
    real :: vctr9a(mza),vctr9b(mza),vctr9c(mza)

    real :: rhs(mza,max(3,num_scalar)), soln(mza,max(3,num_scalar))
    real :: del_scp(mza,2)
    real :: vctr2(mza,3)

    real, parameter :: rchmax =    3.   ! Test with new asympt vert scale length
    real, parameter :: rmin   = -100.
    real, parameter :: rmax   = 1. / 3. ! 1. / zkhkm(mrl)

!!!!!
!!
!!  FIRST COMPUTE THE EDDY DIFFUSIVITIES, AND SAVE
!!  THE HORIZONTAL DIFFUSIVITIES HKM FOR USE ELSEWHERE
!!
!!!!!

    do k = 2,mza-2
       dzim2(k) = dzim(k) * dzim(k)
    enddo

    ka = lpw(iw)

    scalen_horiz = csx(mrlw) * sqrt(arw0(iw))  ! change this later?
    scalen_asympt = csx(mrlw) * 300.

    ! Loop over T levels

    do k = ka,mza-1
       thetav(k) = theta(k,iw) * (1. + .61 * sh_v(k,iw))
    enddo

    ! Loop over W levels

    do k = ka,mza-2

       ! Vertical strain rate squared (based on dV/dz only)

       strain2(k) = dzim2(k) * ((vxe(k+1,iw) - vxe(k,iw))**2  &
                             + (vye(k+1,iw) - vye(k,iw))**2  &
                             + (vze(k+1,iw) - vze(k,iw))**2)

       bvfreq2(k) = grav2 * dzim(k)  &
                  * (thetav(k+1) - thetav(k)) / (thetav(k+1) + thetav(k))

    enddo

    do k = ka,mza-2

       !  Compute Richardson number and Lilly Richardson-number term

       richnum = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-15)))
       richnum_term = min(rchmax,sqrt(max(0.,(1.-zkhkm(mrlw)*richnum))))

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

       vkh(k) = vkm(k) * zkhkm(mrlw)  ! Need to change this factor later
      
    enddo
      
    ! Zero values for top and bottom boundaries   
   
    vkm(ka-1)  = 0.
    vkh(ka-1)  = 0.
    vkm(mza-1) = 0.
    vkh(mza-1) = 0.

    ! Horizontal diffusion coefficient (current version)

    bkmin = akmin(1) * .075 * arw0(iw) ** .666667

    ! akmin hardwired for "grid 1" value for now

    do k = ka,mza-1
       hkm(k,iw) = max(vkm(k-1), vkm(k), bkmin * real(rho(k,iw)))
    enddo

!!!!!
!!
!! DIAGNOSE PBL HEIGHT, WHICH IS OFTEN NEEDED BY OTHER PARAMETERIZATIONS
!!
!!!!!

    ! Compute surface buoyancy flux

    sbf = grav / thetav(ka)  &      
         * (sflux_t(iw) * (1. + .61 * sh_v(ka,iw))  &
         +  sflux_r(iw) * .61 * theta(ka,iw))

    if (sbf < 0.0) then
       
       ! Stable case (from Zilitinkevich 1972)      

       pblh(iw) = 2.4e3 * ustar(iw) ** 1.5
      
    else
      
       ! Unstable case - using gradient Richardon number...

       ! Initialize boundary layer height (minimum value is 1/2 deltaz)

       pblh(iw) = zt(ka) - zm(ka-1)
   
       ! Vertical loop over W levels

       do k = ka, mza-2
      
          !  Compute Richardson number and Lilly Richardson-number term

          richnum      = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-7)))
          richnum_term = min(10.,sqrt(max(0.,(1.-zkhkm(mrlw)*richnum))))

          ! Find the lowest level where richnum_term < 1.e-5 (the lowest neutral
          ! or stable layer).  Take zt(k) to be top of unstable layer.

          if (richnum_term <= 1.e-5) then
             pblh(iw) = zt(k) - zm(ka-1)
             exit
          endif

       enddo
      
    endif

    ! Find the highest model level zt(k) that is at or below
    ! boundary layer height pblh

    kpblh(iw) = ka
    do while (zm(kpblh(iw)) < zm(ka-1) + pblh(iw))
       kpblh(iw) = kpblh(iw) + 1
    enddo

!!!!!
!!
!! COMPUTE THE MOMENTUM TENDENCIES DUE TO VERTICAL MIXING
!!
!!!!!

    dtl  = dtlm(itab_w(iw)%mrlw)
    dtli = 1. / dtl

    ! Vertical loop over W levels

    do k = ka,mza-2
       akodz(k) = arw(k,iw) * vkm(k) * dzim(k)
    enddo

    akodz(ka-1)  = 0.
    akodz(mza-1) = 0.

    ! Distribution of surface flux over multiple levels
    
    vctr3(1:mza) = 0.

    do k = ka, ka + lsw(iw) - 1
       vctr3(k) = (arw(k,iw) - arw(k-1,iw)) * vkm_sfc(iw) * dzim(k-1) * 2.0
    enddo

    ! Vertical loop over T levels

    do k = ka,mza-1
       
       ! Mass in t control volume and its inverse times dtl
       
       dtomass(k) = dtl / ( rho(k,iw) * volt(k,iw) )

       ! Fill tri-diagonal matrix coefficients

       vctr5(k) = -dtomass(k) * akodz(k-1)
       vctr7(k) = -dtomass(k) * akodz(k)
       vctr6(k) = 1. - vctr5(k) - vctr7(k) + dtomass(k) * vctr3(k)

       ! Fill r.h.s. vectors

       rhs(k,1) = vxe(k,iw)
       rhs(k,2) = vye(k,iw)
       rhs(k,3) = vze(k,iw)

    enddo

    ! Solve tri-diagonal matrix for each component

    if (ka <= mza-1) then
       call tridv( vctr5, vctr6, vctr7, rhs, soln, ka, mza-1, mza, 3 )
    endif

    ! Now, soln contains velocity(t+1) values

    ! Compute internal vertical turbulent fluxes

    do k = ka, mza-2
       vctr2(k,1) = akodz(k) * (soln(k,1) - soln(k+1,1))
       vctr2(k,2) = akodz(k) * (soln(k,2) - soln(k+1,2))
       vctr2(k,3) = akodz(k) * (soln(k,3) - soln(k+1,3))
    enddo

    ! Set bottom and top internal fluxes to zero
    
    vctr2(ka-1, 1:3) = 0.
    vctr2(mza-1,1:3) = 0.

    ! Compute tendencies

    do k = ka, mza-1
       vmxet(k,iw) = vmxet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,1) - vctr2(k,1) - vctr3(k) * soln(k,1))
       vmyet(k,iw) = vmyet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,2) - vctr2(k,2) - vctr3(k) * soln(k,2))
       vmzet(k,iw) = vmzet(k,iw) + volti(k,iw) &
                   * (vctr2(k-1,3) - vctr2(k,3) - vctr3(k) * soln(k,3))
    enddo

!!!!!
!!
!! APPLY SURFACE HEAT AND MOISTURE FLUXES DIRECTLY TO TENDENCY ARRAYS
!!
!!!!!
    
    ! Vertical loop over T levels that are adjacent to surface

    do ks = 1,lsw(iw)
       k = ka + ks - 1

       ! Apply surface heat xfer [kg_a K] directly to thilt [kg_a K / s]

       thilt(k,iw) = thilt(k,iw) + dtli * volti(k,iw) * sxfer_tk(ks,iw)

       ! Apply surface vapor xfer [kg_vap] directly to sh_wt [kg_vap / (m^3 s)]

       sh_wt(k,iw) = sh_wt(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

       ! Apply surface vapor xfer [kg_vap] directly to rhot [kg_air / s]

       rhot(k,iw) = rhot(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

       if (allocated(fthpbl)) then
          fthpbl(k,iw) = dtli * volti(k,iw) * sxfer_tk(ks,iw) / rho(k,iw)
       endif
       
       if (allocated(fqtpbl)) then
          fqtpbl(k,iw)  = dtli * volti(k,iw) * sxfer_rk(ks,iw) / rho(k,iw)
       endif
    enddo

!!!!!
!!
!! COMPUTE TEMPERATURE AND SCALAR TENDENCIES DUE TO VERTICAL MIXING
!! THE SOLVER COMPUTES THE FUTURE FLUXES, NOT THE FUTURE SCALAR VALUES
!!
!!!!!

    ! Vertical loop over W levels: fill tri-diagonal matrix coefficients and r.h.s.

    do k = ka, mza-2
       akodz(k) = arw(k,iw) * vkh(k) * dzim(k)
       vctr5(k) = -akodz(k) * dtomass(k)
       vctr7(k) = -akodz(k) * dtomass(k+1)
       vctr6(k) = 1. - vctr5(k) - vctr7(k)

       do n = 1, num_scalar
          rhs(k,n) =  akodz(k) * (scalar_tab(n)%var_p(k,iw) - scalar_tab(n)%var_p(k+1,iw))
       enddo
    enddo

    ! Include the changes in the future fluxes due to surface heat and vapor transfer

    do ks = 1,lsw(iw)
       k = ka + ks - 1
       del_scp(k,1) = sxfer_tk(ks,iw) / (rho(k,iw) * volt(k,iw))
       del_scp(k,2) = sxfer_rk(ks,iw) / (rho(k,iw) * volt(k,iw))
    enddo

    del_scp(ka+lsw(iw),:) = 0.0
    
    do ks = 1,lsw(iw)
       k = ka + ks - 1
       rhs(k,1) = rhs(k,1) + akodz(k) * (del_scp(k,1) - del_scp(k+1,1)) ! temperature
       rhs(k,2) = rhs(k,2) + akodz(k) * (del_scp(k,2) - del_scp(k+1,2)) ! water vapor
    enddo

    ! Solve tri-diagonal matrix equation

    if (ka <= mza-2) then
       call tridv( vctr5, vctr6, vctr7, rhs, soln, ka, mza-2, mza, num_scalar )
    endif

    ! Set bottom and top vertical internal turbulent fluxes to zero

    soln( ka-1, 1:num_scalar) = 0.0
    soln(mza-1, 1:num_scalar) = 0.0

    ! Vertical loop over T levels

    do n = 1, num_scalar
       do k = ka, mza-1
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + volti(k,iw) * (soln(k-1,n) - soln(k,n))
       enddo
    enddo

    ! special for Grell scheme inputs:
    do k = ka, mza-1
       if (allocated(fthpbl)) then
          fthpbl(k,iw) = fthpbl(k,iw) + dtli * volti(k,iw) * (soln(k-1,1) - soln(k,1))
       endif
       
       if (allocated(fqtpbl)) then
          fqtpbl(k,iw) = fqtpbl(k,iw) + dtli * volti(k,iw) * (soln(k-1,2) - soln(k,2))
       endif
    enddo
    
  end subroutine turb_k


  subroutine tridv ( l, d, u, b, x, ka, kz, nz, nsp )

!   Solves tridiagonal system by Thomas algorithm. 
!   The associated tri-diagonal system is stored in 3 arrays
!   D : diagonal
!   L : sub-diagonal
!   U : super-diagonal
!   B : right hand side function
!   X : return solution from tridiagonal solver

!     [ D(1) U(1) 0    0    0 ...       0     ]
!     [ L(2) D(2) U(2) 0    0 ...       .     ]
!     [ 0    L(3) D(3) U(3) 0 ...       .     ]
!     [ .       .     .     .           .     ] X(i) = B(i)
!     [ .             .     .     .     0     ]
!     [ .                   .     .     .     ]
!     [ 0                           L(n) D(n) ]

!-----------------------------------------------------------------------

    implicit none
      
! Arguments:
    
    integer, intent(in)  :: ka, kz, nz
    integer, intent(in)  :: nsp

    real,    intent(in)  :: l(nz)      ! subdiagonal
    real,    intent(in)  :: d(nz)      ! diagonal
    real,    intent(in)  :: u(nz)      ! superdiagonal
    real,    intent(in)  :: b(nz,nsp)  ! r.h. side
    real,    intent(out) :: x(nz,nsp)  ! solution

! Local Variables:

    real    ::  gam(kz)
    real    ::  bet
    integer ::  v, k

! Decomposition and forward substitution:

    bet = 1.0 / d( ka )
    do v = 1, nsp
       x( ka,v ) = bet * b(ka,v)
    enddo

    do k = ka+1, kz
       gam(k) = bet * u( k-1 )
       bet = 1.0 / ( d( k ) - l( k ) * gam( k ) )
       do v = 1, nsp
          x( k,v ) = bet * ( b( k,v ) - l( k ) * x( k-1,v ) )
       enddo
    enddo

! Back-substitution:

    do v = 1, nsp
       do k = kz - 1, ka, -1
          x( k,v ) = x( k,v ) - gam( k+1 ) * x( k+1,v )
       enddo
    enddo
     
  end subroutine tridv


end module smagorinsky
