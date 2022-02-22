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

Module mem_rayf

  implicit none

! Memory for Rayleigh friction layer

  real :: rayf_zmin
  real :: rayf_expon
  real :: rayf_fact

  real :: rayfw_zmin
  real :: rayfw_expon
  real :: rayfw_fact

  real :: rayfdiv_zmin
  real :: rayfdiv_expon
  real :: rayfdiv_fact

  real :: rayfmix_zmin
  real :: rayfmix_expon
  real :: rayfmix_fact

  real, allocatable :: rayf_cof   (:)
  real, allocatable :: rayf_cofw  (:)
  real, allocatable :: rayf_cofdiv(:)
  real, allocatable :: rayf_cofmix(:)

  integer :: krayf_bot
  integer :: krayfw_bot
  integer :: krayfdiv_bot
  integer :: krayfmix_bot

  logical :: dorayf    = .false.
  logical :: dorayfw   = .false.
  logical :: dorayfdiv = .false.
  logical :: dorayfmix = .false.

Contains

  subroutine rayf_init()

! Initialize Rayleigh friction vertical profile coefficients

    use misc_coms,  only: dtsm, initial
    use mem_grid,   only: mza, zm, zt
    use oname_coms, only: nl

    implicit none

    integer :: k
    real    :: distimi, distim0, dti

    dti = 1.0 / real(dtsm(1))

    allocate( rayf_cof   (mza) )
    allocate( rayf_cofw  (mza) )
    allocate( rayf_cofdiv(mza) )
    allocate( rayf_cofmix(mza) )

    rayf_cof    = 0.0
    rayf_cofw   = 0.0
    rayf_cofdiv = 0.0
    rayf_cofmix = 0.0

    krayf_bot    = mza + 1
    krayfw_bot   = mza + 1
    krayfdiv_bot = mza + 1
    krayfmix_bot = mza + 1

! RAYF COEFFICIENT FOR THIL AND VMC (ONLY FOR HORIZ. HOMOG. INITIALIZATION)

    if (rayf_fact > 1.e-7 .and. initial == 1) then

       do k = 2, mza
          if (zt(k) > rayf_zmin) then
             krayf_bot = k
             exit
          endif
       enddo

       if (krayf_bot <= mza) then

          dorayf  = .true.
          distimi = rayf_fact * dti

          do k = krayf_bot, mza
             rayf_cof(k) = distimi   &
                  * ((zt(k) - rayf_zmin) / (zm(mza) - rayf_zmin)) ** rayf_expon
          enddo

       endif
    endif

! RAYF coefficient for WMC damping

    if (rayfw_fact > 1.e-7) then

       do k = 2, mza-1
          if (zm(k) > rayfw_zmin) then
             krayfw_bot = k
             exit
          endif
       enddo

       if (krayfw_bot < mza) then

          dorayfw = .true.
          distimi = rayfw_fact * dti

          do k = krayfw_bot, mza-1
             rayf_cofw(k) = distimi   &
                  * ((zm(k) - rayfw_zmin) / (zm(mza) - rayfw_zmin)) ** rayfw_expon
          enddo

       endif
    endif

! RAYF coefficient for Horiz Divergence

    if (rayfdiv_fact > 1.e-7) then

       do k = 2, mza
          if (zt(k) > rayfdiv_zmin) then
             krayfdiv_bot = k
             exit
          endif
       enddo

       if (krayfdiv_bot <= mza) then

          dorayfdiv = .true.
          distim0   = max(nl%divh_damp_fact, 0.) * dti
          distimi   = max(rayfdiv_fact,nl%divh_damp_fact) * dti - distim0

          do k = krayfdiv_bot, mza
             rayf_cofdiv(k) = distim0 + distimi   &
                  * ((zt(k) - rayfdiv_zmin) / (zm(mza) - rayfdiv_zmin)) ** rayfdiv_expon
          enddo

       endif
    endif

! RAYF coefficient for horizontal velocity (VC) mixing

    if (rayfmix_fact > 1.e-7) then

       do k = 2, mza-1
          if (zm(k) > rayfmix_zmin) then
             krayfmix_bot = k
             exit
          endif
       enddo

       if (krayfmix_bot < mza) then

          dorayfmix = .true.

          do k = krayfmix_bot, mza-1
             rayf_cofmix(k) = rayfmix_fact &
                  * ((zm(k) - rayfmix_zmin) / (zm(mza) - rayfw_zmin)) ** rayfmix_expon
          enddo

       endif
    endif

  end subroutine rayf_init


  subroutine rayf_mix_top_vc( iv, vmt )

    use mem_ijtabs,  only: itab_v
    use mem_basic,   only: vc, rho
    use mem_grid,    only: mza, arw, volvi

    implicit none

    integer, intent(in   ) :: iv
    real,    intent(inout) :: vmt(mza)

    real    :: vflux(mza)
    integer :: k, iw1, iw2

    iw1 = itab_v(iv)%iw(1)
    iw2 = itab_v(iv)%iw(2)

    vflux (mza)           = 0.0
    vflux(krayfmix_bot-1) = 0.0

    do k = krayfmix_bot, mza-1
       vflux(k) = 0.25 * (arw(k,iw1) + arw(k,iw2)) * rayf_cofmix(k) * (vc(k,iv) - vc(k+1,iv)) &
            * real(rho(k+1,iw1) + rho(k+1,iw2) + rho(k,iw1) + rho(k,iw2))
    enddo

    do k = krayfmix_bot, mza
       vmt(k) = vmt(k) + (vflux(k-1) - vflux(k)) * volvi(k,iv)
    enddo

  end subroutine rayf_mix_top_vc


end module mem_rayf
