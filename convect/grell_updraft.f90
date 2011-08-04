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
subroutine cup_up_he(k2, k22, hkb, z_cup, cd, entr, he_cup, hc,  &
     kbcon, dby, he, hes_cup)

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: k22
  real, intent(out) :: hkb
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: cd
  real, intent(in) :: entr
  real, dimension(k2), intent(in) :: he_cup
  real, dimension(k2), intent(out) :: hc
  integer, intent(in) :: kbcon
  real, dimension(k2), intent(out) :: dby
  real, dimension(k2), intent(in) :: he
  real, dimension(k2), intent(in) :: hes_cup
  real :: dz
  integer :: k
  

!--- Moist static energy inside cloud
      
  hkb = he_cup(k22)
  hc(1:k2) = 0.0
  dby(1:k2) = 0.0

  do k = 1, k22
     hc(k) = he_cup(k)
  enddo

  do k = k22, kbcon - 1
     hc(k) = hkb
  enddo

  k = kbcon
  hc(k) = hkb
  dby(kbcon) = hkb - hes_cup(k)
  
  do k = 2, k2 - 1
     if(k > kbcon)then
        dz = z_cup(k) - z_cup(k-1)
        hc(k) = (hc(k-1) * (1.0 - 0.5 * cd(k) * dz) + entr *    &
             dz * he(k-1)) / (1.0 + entr * dz - 0.5 * cd(k) * dz)
        dby(k) = hc(k) - hes_cup(k)
     endif
  enddo

  return
end subroutine cup_up_he

!===============================================================================

subroutine cup_up_nms(zu, z_cup, entr, cd, kbcon, ktop, k2, ierr, k22)

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(out) :: zu
  real, dimension(k2), intent(in) :: z_cup
  real, intent(in) :: entr
  real, dimension(k2), intent(in) :: cd
  integer, intent(in) :: kbcon
  integer, intent(in) :: ktop
  integer, intent(in) :: ierr
  integer, intent(in) :: k22
  integer :: k
  real :: dz

  zu(1:k2) = 0.0
  if(ierr == 0)then
     zu(k22:kbcon) = 1.0
     do k = kbcon + 1, ktop
        dz = z_cup(k) - z_cup(k-1)
        zu(k) = zu(k-1) * (1.0 + (entr - cd(k)) * dz)
     enddo
  endif

  return
end subroutine cup_up_nms

!===============================================================================

subroutine cup_up_moisture(k2, ierr, z_cup, qc, qrc, pw, pwav, &
     kbcon, ktop, cd, dby, mentr_rate, q, gamma_cup, zu, qes_cup,   &
     k22, qe_cup)

  use consts_coms, only: alvl

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(out) :: qc
  real, dimension(k2), intent(out) :: qrc
  real, dimension(k2), intent(out) :: pw
  real, intent(out) :: pwav
  integer, intent(in) :: kbcon
  integer, intent(in) :: ktop
  real, dimension(k2), intent(in) :: cd
  real, dimension(k2), intent(in) :: dby
  real, intent(in) :: mentr_rate
  real, dimension(k2), intent(in) :: q
  real, dimension(k2), intent(in) :: gamma_cup
  real, dimension(k2), intent(in) :: zu
  real, dimension(k2), intent(in) :: qes_cup
  integer, intent(in) :: k22
  real, dimension(k2), intent(in) :: qe_cup
  real :: c0
  real :: radius
  integer :: k
  real :: dz
  real :: qrch

!---  No precip for small clouds

  c0 = 0.002
  if(mentr_rate > 0.)then
     radius = 0.2 / mentr_rate
     if(radius < 900.)c0 = 0.0
  endif
  
  pwav = 0.
  
  pw(1:k2) = 0.
  qc(1:k2) = qe_cup(1:k2)
  qrc(1:k2) = 0.
  
  if(ierr == 0)qc(k22:(kbcon-1)) = qe_cup(k22)
  
  do k = 2, k2 - 1
     if(ierr /= 0 .or. k < kbcon .or. k > ktop) cycle
     dz = z_cup(k) - z_cup(k-1)
            
!---  1. steady state plume equation, for what could
!---  be in cloud without condensation
            
     qc(k) = (qc(k-1) * (1.0 - 0.5 * cd(k) * dz) + mentr_rate *    &
          dz * q(k-1)) / (1.0 + mentr_rate * dz - 0.5 * cd(k) * dz)
            
!---  saturation  in cloud, this is what is allowed to be in it

     qrch = qes_cup(k) + (1.0 / alvl) * (gamma_cup(k) /   &
          (1. + gamma_cup(k))) * dby(k)

!---  liquid water content in cloud after rainout

     qrc(k) = (qc(k) - qrch) / (1.0 + c0 * dz)
     if(qrc(k) < 0.0) qrc(k) = 0.0

!---  3.Condensation

     pw(k) = c0 * dz * qrc(k) * zu(k) !unit: kg[liq water]/kg[air]
                                !unit of c0 is m^-1

!---  Set next level

     qc(k) = qrc(k) + qrch

!---  Integrated normalized ondensate

     pwav = pwav + pw(k)

  enddo

  return
end subroutine cup_up_moisture

!===============================================================================

subroutine cup_up_aa0(k2, aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr)

  use consts_coms, only: grav, cp

  implicit none

  integer, intent(in) :: k2
  real, intent(out) :: aa0
  real, dimension(k2), intent(in) :: z
  real, dimension(k2), intent(in) :: zu
  real, dimension(k2), intent(in) :: dby
  real, dimension(k2), intent(in) :: gamma_cup
  real, dimension(k2), intent(in) :: t_cup
  integer, intent(in) :: kbcon
  integer, intent(in) :: ktop
  integer, intent(in) :: ierr
  integer :: k
  real :: dz
  real :: da


  aa0 = 0.0

  do k = 2, k2 - 1
     if(ierr /= 0 .or. k <= kbcon .or. k > ktop) cycle
     dz = z(k) - z(k-1)
     da = zu(k) * dz * (grav / (cp * t_cup(k))) * dby(k-1) /    &
          (1.0 + gamma_cup(k))
     if(k == ktop .and. da <= 0.0) cycle
     aa0 = aa0 + da
     if(aa0 < 0.0)aa0 = 0.0
  enddo

  return
end subroutine cup_up_aa0
