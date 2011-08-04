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
subroutine cup_env(k2, t, p, tcrit, itest, z1, psur, q, he, hes, qes, z)

  use consts_coms, only: alvl, cp, alvi, rdry, grav, cpi

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: t
  real, dimension(k2), intent(in) :: p
  real, intent(in) :: tcrit
  integer, intent(in) :: itest
  real, intent(in) :: z1
  real, intent(in) :: psur

  real, dimension(k2), intent(inout) :: q
  real, dimension(k2), intent(inout) :: he
  real, dimension(k2), intent(out) :: hes

  real, dimension(k2), intent(out) :: qes
  real, dimension(k2), intent(out) :: z

  real, dimension(2) :: ht
  real, dimension(2) :: be
  real, dimension(2) :: ae
  integer :: k
  integer :: iph
  real :: e
  real, dimension(k2) :: tv
  real :: tvbar

  ht(1) = alvl * cpi
  ht(2) = alvi * cpi
  be(1) = .622 * ht(1) / .286
  be(2) = .622 * ht(2) / .286
  ae(1) = be(1) / 273. + log(610.71)
  ae(2) = be(2) / 273. + log(610.71)

  do k=1,k2
     !sgb - iph is for phase, dependent on tcrit (water or ice)
     iph = 1
     if(t(k) <= tcrit)  iph = 2
     e = exp(ae(iph) - be(iph) / t(k))
     e = min(e, 90.*p(k))
     qes(k) = .622 * e / (100. * p(k) - e)
     if(qes(k) <= 1.0e-8) qes(k) = 1.0e-8
     if(q(k) > qes(k)) q(k) = qes(k)
     tv(k) = t(k) * (1.0 + 0.608 * q(k))
  enddo

!--- z's are calculated with changed h's and q's and t's
!--- if itest=2

  if(itest /= 2)then
     z(1) = z1 - (log(p(1)) - log(psur)) * rdry * tv(1) / grav

!--- calculate heights
     do k=2,k2
        tvbar = .5 * tv(k) + .5 * tv(k-1)
        z(k) = z(k-1) - (log(p(k)) - log(p(k-1))) * rdry * tvbar / grav
     enddo

  else

     do k=1,k2
        z(k) = max(1.0e-3, (he(k) - cp * t(k) - alvl * q(k)) / grav)
     enddo

  endif

!--- calculate moist static energy - he
!--- saturated moist static energy - hes
  
  do k=1,k2
     if(itest == 0) he(k) = grav * z(k) + cp * t(k) + alvl * q(k)
     hes(k) = grav * z(k) + cp * t(k) + alvl * qes(k)
     if(he(k) >= hes(k)) he(k) = hes(k)
  enddo

  return
end subroutine cup_env

!======================================================================

subroutine cup_env_clev(k2, psur, z1, q, q_cup, z, z_cup, p, p_cup,  &
     t, t_cup, qes, qes_cup, hes, hes_cup, he, he_cup, gamma_cup)

  use consts_coms, only: alvl, cp, rvap

  implicit none

  integer, intent(in) :: k2
  real, intent(in) :: psur
  real, intent(in) :: z1
  real, dimension(k2), intent(in) :: q
  real, dimension(k2), intent(out) :: q_cup
  real, dimension(k2), intent(in) :: z
  real, dimension(k2), intent(out) :: z_cup
  real, dimension(k2), intent(in) :: p
  real, dimension(k2), intent(out) :: p_cup
  real, dimension(k2), intent(in) :: t
  real, dimension(k2), intent(out) :: t_cup
  real, dimension(k2), intent(in) :: qes
  real, dimension(k2), intent(out) :: qes_cup
  real, dimension(k2), intent(in) :: hes
  real, dimension(k2), intent(out) :: hes_cup
  real, dimension(k2), intent(in) :: he
  real, dimension(k2), intent(out) :: he_cup
  real, dimension(k2), intent(out) :: gamma_cup
  integer :: k

  do k=2,k2
     qes_cup(k) = .5 * (qes(k-1) + qes(k))
     q_cup(k) = .5 * (q(k-1) + q(k))
     hes_cup(k) = .5 * (hes(k-1) + hes(k))
     he_cup(k) = .5 * (he(k-1) +  he(k))
     if(he_cup(k) > hes_cup(k)) he_cup(k) = hes_cup(k)
     z_cup(k) = .5*(z(k-1) + z(k))
     p_cup(k) = .5*(p(k-1) + p(k))
     t_cup(k) = .5*(t(k-1) + t(k))
     gamma_cup(k) = (alvl / cp) * (alvl / (rvap * t_cup(k)**2)) * qes_cup(k)
  enddo

  qes_cup(1) = qes(1)
  q_cup(1) = q(1)
  hes_cup(1) = hes(1)
  he_cup(1) = he(1)
  z_cup(1) =  z1
  p_cup(1) = psur
  t_cup(1) = t(1)
  
  gamma_cup(1) = alvl / cp * (alvl / (rvap * t_cup(1)**2)) * qes_cup(1)

  return
end subroutine cup_env_clev

!=================================================================

subroutine get_zi(k2, z, rcpg, tkeg, ztop, tkmin, kzi)

  implicit none
      
  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: z
  real, dimension(k2), intent(in) :: rcpg
  real, dimension(k2), intent(in) :: tkeg
  real, intent(in) :: ztop
  real, intent(in) :: tkmin
  integer, intent(out) :: kzi

  real, parameter :: rcpmin = 1.0e-5 
  real, parameter :: pblhmax = 3000.0
  real :: tke_tmp
  integer :: ktke_max
  integer :: kzimax
  integer :: k

  kzi  = 2
  tke_tmp = 0.
  ktke_max= 1
  kzimax = 0

!---  max level for kzi
  find_maxlev: do k=1,k2
     if(z(k) >= (pblhmax+ztop)) then
        kzimax = k
        exit find_maxlev
     endif
  enddo find_maxlev

!level of max tke  below kzimax and w/out clouds
  find_maxtke: do k = 1, kzimax
     if(rcpg(k) < rcpmin) then
        if( tkeg(k) >= tke_tmp) then
           tke_tmp = tkeg(k)
           ktke_max = k
        else
           exit find_maxtke
        endif
     endif
  enddo find_maxtke

  find_kzi: do k = ktke_max, kzimax + 1
     if(rcpg(k) < rcpmin) then
        if(tkeg(k) > (1.1 * tkmin))then
           kzi = k
        endif
     else
        kzi = k
        exit find_kzi
     endif
  enddo find_kzi

  kzi = min(kzimax, max(2, kzi))

  return
end subroutine get_zi

!==========================================================================

subroutine cup_kbcon(k22, kbcon, he_cup, hes_cup, ierr,   &
     kbmax, p_cup, cap_max, k2)

  implicit none

  integer, intent(inout) :: k22
  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: he_cup
  real, dimension(k2), intent(in) :: p_cup
  real, dimension(k2), intent(in) :: hes_cup
  integer, intent(in) :: kbmax
  real, intent(in) :: cap_max
  integer, intent(inout) :: ierr
  integer, intent(out) :: kbcon
  real :: pbcdif
  real :: plus

!--- Determine the level of convective cloud base  - KBCON

  kbcon = 1
  if(ierr /= 0)return

  kbcon = k22

  pbcdiff_loop: do

     kbcon_cycle: do
        if(he_cup(k22) < hes_cup(kbcon))then
           kbcon = kbcon + 1
           if(kbcon > (kbmax + 2))then
              ierr = 3
              return
           endif
        else
           exit kbcon_cycle
        endif
     enddo kbcon_cycle
     
     if( (kbcon - k22) == 1)return
     
     ! Cloud base pressure and max moist static energy pressure
     ! i.e., the depth (in mb) of the layer of negative buoyancy
     pbcdif = p_cup(k22) - p_cup(kbcon)
     if(pbcdif <= cap_max)return
     k22 = k22 + 1
     kbcon = k22
  enddo pbcdiff_loop

  return
end subroutine cup_kbcon

!==========================================================================

subroutine cup_kbcon_catt(iloop, k22, kbcon, he_cup, hes_cup, ierr,   &
     kbmax, p_cup, cap_max, k2)

  implicit none

  integer, intent(inout) :: k22
  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: he_cup
  real, dimension(k2), intent(in) :: p_cup
  real, dimension(k2), intent(in) :: hes_cup
  integer, intent(in) :: kbmax
  integer, intent(in) :: iloop
  real, intent(in) :: cap_max
  integer, intent(inout) :: ierr
  integer, intent(out) :: kbcon
  real :: pbcdif
  real :: plus

!--- Determine the level of convective cloud base  - KBCON

  kbcon = 1
  if(ierr /= 0)return

  kbcon = k22

  pbcdiff_loop: do

     kbcon_cycle: do
        if(he_cup(k22) < hes_cup(kbcon))then
           kbcon = kbcon + 1
           if(kbcon > (kbmax + 2))then
              if(iloop < 4) ierr = 3
              if(iloop == 4) ierr = 997
              return
           endif
        else
           exit kbcon_cycle
        endif
     enddo kbcon_cycle
     
     if( (kbcon - k22) == 1)return
     
     ! Cloud base pressure and max moist static energy pressure
     ! i.e., the depth (in mb) of the layer of negative buoyancy
     pbcdif = p_cup(k22) - p_cup(kbcon)
     plus = max(25.0, cap_max)
     if(pbcdif <= plus)return
     k22 = k22 + 1
     kbcon = k22
  enddo pbcdiff_loop

  return
end subroutine cup_kbcon_catt

!==========================================================================

subroutine cup_ktop(ilo, dby, kbcon, ktop, k2, ierr)
  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: ilo
  real, dimension(k2), intent(inout) :: dby
  integer, intent(in) :: kbcon
  integer, intent(out) :: ktop
  integer, intent(inout) :: ierr
  integer :: k
  
  ktop = 1
  if(ierr == 0)then
     do k = kbcon + 1, k2 - 2
        if(dby(k) <= 0.0)then
           ktop = k - 1
           dby((ktop+1):k2) = 0.0
           return
        endif
     enddo
     if(ilo == 1) ierr = 5
     if(ilo == 2) ierr = 998
  endif

  return
end subroutine cup_ktop

!======================================================================

subroutine cup_direction2(imass, massflx, massflx_upwind, id, id_upwind,   &
     iresult, massfld)

  implicit none

  integer, intent(in) :: imass
  integer, intent(in) :: id
  integer, intent(in) :: id_upwind
  real, intent(in) :: massflx
  real, intent(in) :: massflx_upwind
  real, intent(out) :: massfld
  integer, intent(out) :: iresult

  if(imass == 1 .and. id_upwind == 1)then
     massfld = max(massflx_upwind, massflx)
  elseif(imass == 1)then
     massfld = massflx
  else
     massfld = 0.0
  endif

  if(id == 1 .or. id_upwind == 1)then
     iresult = 1
  else
     iresult = 0
  endif

  return
end subroutine cup_direction2

!=========================================================================

subroutine maximi(k2, array, ks, ke, ierr, maxx)

  implicit none
  
  integer, intent(in) :: k2
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: array
  integer, intent(in) :: ks
  integer, intent(in) :: ke
  integer, intent(out) :: maxx
  real :: x
  integer :: k

  maxx = ks
  if(ierr == 0)then
     x = array(ks)
     do k = ks, ke
        if(array(k) >= x)then
           x = array(k)
           maxx = k
        endif
     enddo
  endif

  return
end subroutine maximi

!========================================================================

subroutine minimi(k2, array, ks, kend, kt, ierr)

  implicit none
  
  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: array
  integer, intent(in) :: ks
  integer, intent(in) :: kend
  integer, intent(out) :: kt
  integer, intent(in) :: ierr
  real :: x
  integer :: kstop
  integer :: k

  kt = ks
  if(ierr == 0)then
     x = array(ks)
     kstop = max(ks + 1, kend)
     do k = ks + 1, kstop
        if(array(k) < x) then
           x = array(k)
           kt = k
        endif
     enddo
  endif

  return
end subroutine minimi
