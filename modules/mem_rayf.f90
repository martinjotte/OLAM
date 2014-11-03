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
  real :: rayf_distim
  real :: rayf_expon

  real :: rayfw_zmin
  real :: rayfw_distim
  real :: rayfw_expon  

  real :: rayfdiv_zmin
  real :: rayfdiv_distim
  real :: rayfdiv_expon  

  real, allocatable :: rayf_cof(:)
  real, allocatable :: rayf_cofw(:)
  real, allocatable :: rayf_cofdiv(:)

  real, allocatable :: vc03d(:,:)
  real, allocatable :: dn03d(:,:)

  integer :: krayf_bot
  integer :: krayfw_bot
  integer :: krayfdiv_bot

  logical :: dorayf    = .false.
  logical :: dorayfw   = .false.
  logical :: dorayfdiv = .false.

Contains

  subroutine rayf_init(mva,mwa,mza)

! Initialize Rayleigh friction vertical profile coefficients

    use misc_coms,    only: initial, rinit, iparallel
    use mem_grid,     only: lpv, zm, zt
    use mem_ijtabs,   only: jtab_v, jtv_init, itab_v
    use mem_basic,    only: rho, vc
    use olam_mpi_atm, only: mpi_send_v, mpi_recv_v

    implicit none
    
    integer, intent(in) :: mva, mwa, mza

    integer :: k, j, iv, iu, iw1, iw2, mrl
    real    :: distimi

    allocate( rayf_cof   (mza) )
    allocate( rayf_cofw  (mza) )
    allocate( rayf_cofdiv(mza) )

    krayf_bot    = mza
    krayfw_bot   = mza
    krayfdiv_bot = mza

! RAYF coefficient for THIL and VMC  

    rayf_cof(1:mza) = 0.

    if (rayf_distim > 1.e-6) then

       dorayf = .true.
       distimi = 1. / rayf_distim
       
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

    if (rayfw_distim > 1.e-6) then

       dorayfw = .true.
       distimi = 1. / rayfw_distim

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

    if (rayfdiv_distim > 1.e-6) then

       dorayfdiv = .true.
       distimi = 1. / rayfdiv_distim

       do k = 2, mza
          if (zt(k) > rayfdiv_zmin) then
             krayfdiv_bot = k
             exit
          endif
       enddo

       do k = krayfdiv_bot, mza
          rayf_cofdiv(k) = distimi   &
               * ((zt(k) - rayfdiv_zmin) / (zm(mza) - rayfdiv_zmin)) ** rayfdiv_expon
       enddo

    endif

! For a horizontally or latitudinally homogeneous run, allocate the arrays
! to store the momentum values that the model will relax towards

    if (dorayf .and. (initial == 1 .or. initial == 3)) then
       allocate( vc03d(mza,mva) ) ; vc03d = rinit
       allocate( dn03d(mza,mva) ) ; dn03d = rinit
    endif

! For a horizontally homogeneous run, the initial momentum that the model will
! relax towards can be saved here, after inithh() and before history_read()
    
    if (dorayf .and. initial == 1) then

       do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
          iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
          do k = lpv(iv), mza
             dn03d(k,iv) = 0.5 * (rho(k,iw1) + rho(k,iw2))
          enddo
       enddo

       mrl = 1

       if (iparallel == 1) call mpi_send_v(mrl, rvara1=dn03d)
          
       do iv = 2, mva
          do k = lpv(iv), mza
             vc03d(k,iv) = vc(k,iv)
          enddo
       enddo
          
       if (iparallel == 1) call mpi_recv_v(mrl, rvara1=dn03d)

    endif ! (dorayf .and. initial)

    return
  end subroutine rayf_init

end module mem_rayf
