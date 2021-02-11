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

  real, allocatable :: rayf_cof   (:)
  real, allocatable :: rayf_cofw  (:)
  real, allocatable :: rayf_cofdiv(:)
  real, allocatable :: rayf_cofdif(:)

  integer :: krayf_bot
  integer :: krayfw_bot
  integer :: krayfdiv_bot
  integer :: krayfdif_bot

  logical :: dorayf    = .false.
  logical :: dorayfw   = .false.
  logical :: dorayfdiv = .false.

Contains

  subroutine rayf_init()

! Initialize Rayleigh friction vertical profile coefficients

    use misc_coms,    only: dtsm
    use mem_grid,     only: mza, zm, zt
    use oname_coms,   only: nl

    implicit none

    integer             :: k
    real                :: distimi, distim0

    allocate( rayf_cof   (mza) )
    allocate( rayf_cofw  (mza) )
    allocate( rayf_cofdiv(mza) )
    allocate( rayf_cofdif(mza) )

    krayf_bot    = mza + 1
    krayfw_bot   = mza + 1
    krayfdiv_bot = mza + 1
    krayfdif_bot = mza + 1

! RAYF coefficient for THIL and VMC

    rayf_cof(1:mza) = 0.

    if (rayf_fact > 1.e-7) then

       dorayf = .true.
       distimi = rayf_fact / real(dtsm(1))

       do k = 2, mza
          if (zt(k) > rayf_zmin) then
             krayf_bot = k
             exit
          endif
       enddo

       do k = krayf_bot, mza
          rayf_cof(k) = distimi   &
               * ((zt(k) - rayf_zmin) / (zm(mza) - rayf_zmin)) ** rayf_expon
       enddo

    endif

! RAYF coefficient for WMC

    rayf_cofw(1:mza) = 0.

    if (rayfw_fact > 1.e-7) then

       dorayfw = .true.
       distimi = rayfw_fact / real(dtsm(1))

       do k = 2, mza
          if (zm(k) > rayfw_zmin) then
             krayfw_bot = k
             exit
          endif
       enddo

       do k = krayfw_bot, mza-1
          rayf_cofw(k) = distimi   &
               * ((zm(k) - rayfw_zmin) / (zm(mza) - rayfw_zmin)) ** rayfw_expon
       enddo

    endif

! RAYF coefficient for Horiz Divergence

    rayf_cofdiv(1:mza) = 0.

    if (rayfdiv_fact > 1.e-7) then

       dorayfdiv = .true.
       distim0   = max(nl%divh_damp_fact, 0.) / real(dtsm(1))
       distimi   = (max(rayfdiv_fact,nl%divh_damp_fact)  - max(nl%divh_damp_fact,0.)) / real(dtsm(1))

       do k = 2, mza
          if (zt(k) > rayfdiv_zmin) then
             krayfdiv_bot = k
             exit
          endif
       enddo

       do k = krayfdiv_bot, mza
          rayf_cofdiv(k) = distim0 + distimi   &
               * ((zt(k) - rayfdiv_zmin) / (zm(mza) - rayfdiv_zmin)) ** rayfdiv_expon
       enddo

    endif

! RAYF coefficient for velocity mixing

    rayf_cofdif(1:mza) = 0.

    do k = 2, mza
       if (zt(k) > 27000.) then
          krayfdif_bot = k
          exit
       endif
    enddo

    do k = krayfdif_bot, mza
       rayf_cofdif(k) = 0.35 * (zt(k) - 27000.) / (zm(mza) - 27000.)
    enddo

  end subroutine rayf_init

end module mem_rayf
