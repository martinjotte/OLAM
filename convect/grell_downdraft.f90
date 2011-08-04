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
subroutine cup_dd_nms(k2, zd, z_cup, cdd, entr, jmin, ierr, itest, kdet, z1)

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(out) :: zd
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(inout) :: cdd
  real, intent(in) :: entr
  integer, intent(in) :: jmin
  integer, intent(in) :: ierr
  integer, intent(in) :: itest
  integer, intent(in) :: kdet
  real, intent(in) :: z1

  real, parameter :: perc = 0.03
  integer :: k
  real :: a
  integer :: ki
  real :: dz

!--- Perc is the percentage of mass left when hitting the ground
  do k = 1, k2
     zd(k) = 0.0
     if(itest == 0) cdd(k) = 0.0
  enddo
  a = 1.0 - perc

  if(ierr == 0)then
     zd(jmin) = 1.0

!--- Integrate downward, specify detrainment(cdd)!
     do ki = jmin - 1, 1, -1
        dz = z_cup(ki+1) - z_cup(ki)
        if(ki <= kdet .and. itest == 0)then
           cdd(ki) = entr + (1.0 - (a * (z_cup(ki) - z1) +   &
                perc * (z_cup(kdet) - z1) ) /   &
                (a * (z_cup(ki+1) - z1) +   &
                perc * (z_cup(kdet) - z1))) / dz
        endif
        zd(ki) = zd(ki+1) * (1.0 + (entr - cdd(ki)) * dz)
     enddo
     
  endif

  return
end subroutine cup_dd_nms

!==================================================================

subroutine cup_dd_he(k2, hes_cup, hcd, z_cup, cdd, entr, jmin, ierr, he, dby)

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: hes_cup
  real, dimension(k2), intent(inout) :: hcd
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: cdd
  real, intent(in) :: entr
  integer, intent(in) :: jmin
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: he
  real, dimension(k2), intent(out) :: dby
  integer :: ki
  real :: dz

  dby(2:k2) = 0.0
  if(ierr == 0)then
     hcd(2:k2) = hes_cup(2:k2)
     hcd(jmin) = hes_cup(jmin)
     dby(jmin) = 0.0
     
     do ki = jmin - 1, 1, -1
        dz = z_cup(ki+1) - z_cup(ki)
        hcd(ki) = (hcd(ki+1) * (1.0 - 0.5 * cdd(ki) * dz) +    &
             entr * dz * he(ki) ) / (1.0 + entr * dz - 0.5 * cdd(ki) * dz)
        dby(ki) = hcd(ki) - hes_cup(ki)
     enddo
  endif

  return
end subroutine cup_dd_he

!=========================================================================

subroutine cup_dd_moisture(k2, zd, hcd, hes_cup, qcd, qes_cup, &
     pwd, q_cup, z_cup, cdd, entr, jmin, ierr, gamma_cup, pwev, bu, qrcd, &
     q, iloop)

  use consts_coms, only: alvl

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: zd
  real, dimension(k2), intent(in) :: hcd
  real, dimension(k2), intent(in) :: hes_cup
  real, dimension(k2), intent(out) :: qcd
  real, dimension(k2), intent(in) :: qes_cup
  real, dimension(k2), intent(out) :: pwd
  real, dimension(k2), intent(in) :: q_cup
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: cdd
  real, intent(in) :: entr
  integer, intent(in) :: jmin
  integer, intent(inout) :: ierr
  real, dimension(k2), intent(in) :: gamma_cup
  real, intent(out) :: pwev
  real, intent(out) :: bu
  real, dimension(k2), intent(out) :: qrcd
  real, dimension(k2), intent(in) :: q
  integer, intent(in) :: iloop
  integer :: k
  real :: dz
  real :: dh
  integer :: ki
  real :: dqeva

  bu = 0.0
  pwev = 0.0
  do k = 1, k2
     qcd(k) = 0.0
     qrcd(k) = 0.0
     pwd(k) = 0.0
  enddo

  if(ierr == 0)then
     k = jmin
     dz = z_cup(k+1) - z_cup(k)
     qcd(k) = q_cup(k)
     qrcd(k) = qes_cup(k)
     pwd(jmin) = min(0.0, qcd(k) - qrcd(k))
     pwev = pwev + pwd(jmin)
     qcd(k) = qes_cup(k)
     dh = hcd(k) - hes_cup(k)
     bu = dz * dh
        
     do ki = jmin - 1, 1, -1
        dz = z_cup(ki+1) - z_cup(ki)
        qcd(ki) = (qcd(ki+1) * (1.0 - 0.5 * cdd(ki) * dz) +   &
             entr * dz * q(ki)) / (1.0 + entr * dz - 0.5 * cdd(ki) * dz)
!--- to be negatively buoyant, hcd should be smaller than hes!
        dh = hcd(ki) - hes_cup(ki)
        bu = bu + dz * dh
        qrcd(ki) = qes_cup(ki) + gamma_cup(ki) / (alvl *  &
             (1.0 + gamma_cup(ki))) * dh
        dqeva = qcd(ki) - qrcd(ki)
        if(dqeva > 0.0) dqeva = 0.0
        pwd(ki) = zd(ki) * dqeva
        qcd(ki) = qrcd(ki)
        pwev = pwev + pwd(ki)
     enddo
     
     if(bu >= 0 .and. iloop == 1) ierr = 7

  endif

  return
end subroutine cup_dd_moisture

!===========================================================================

subroutine cup_dd_edt(k2, ierr, us, vs, z, ktop, kbcon, edt, p, pwav,   &
     pwev, edtmax, edtmin, edtc, vshear, sdp, vws)

  use grell_deep_coms, only: maxens2

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: us
  real, dimension(k2), intent(in) :: vs
  real, dimension(k2), intent(in) :: z
  integer, intent(in) :: ktop
  integer, intent(in) :: kbcon
  real, intent(out) :: edt
  real, dimension(k2), intent(in) :: p
  real, intent(in) :: pwav
  real, intent(in) :: pwev
  real, intent(in) :: edtmax
  real, intent(in) :: edtmin
  real, dimension(k2), intent(inout) :: edtc
  real, intent(out) :: vshear
  real, intent(out) :: sdp
  real, intent(out) :: vws

  integer :: kk
  real :: pef
  real :: pefb
  real :: zkbc
  real :: prezk
  real :: einc
  integer :: k

!--- determine downdraft strength in terms of windshear

! */ calculate an average wind shear over the depth of the cloud

  edt = 0.0
  vws = 0.0
  sdp = 0.0
  vshear = 0.0

  if(ierr == 0)then
     do kk = 1, k2 - 1
        if (kk <= min(ktop,k2-1) .and. kk >= kbcon) then
           vws = vws + (abs((us(kk+1) - us(kk)) / (z(kk+1) - z(kk))) +   &
                abs((vs(kk+1) - vs(kk)) / (z(kk+1) - z(kk)))) *    &
                (p(kk) - p(kk+1))
           sdp = sdp + p(kk) - p(kk+1)
        endif
        if (kk == (k2-1)) vshear = 1000.0 * vws / sdp
     enddo
     
     pef = ((-0.00496 * vshear + 0.0953) * vshear - 0.639) * vshear + 1.591 
     if(pef > edtmax) pef = edtmax
     if(pef < edtmin) pef = edtmin

!--- cloud base precip efficiency
     zkbc = z(kbcon) * 0.003281
     if(zkbc > 25.0)then
        prezk = 2.4
     elseif(zkbc > 3.0)then
        prezk = 0.96729352 + zkbc * (-0.70034167 + zkbc *   &
             (0.162179896 + zkbc * (-1.2569798e-2 + zkbc *   &
             (4.2772e-4 - zkbc * 5.44e-6))))
     else
        prezk = 0.02
     endif
     pefb = 1.0 / (1.0 + prezk)
     if(pefb > edtmax)pefb = edtmax
     if(pefb < edtmin)pefb = edtmin
     edt = 1.0 - 0.5 * (pefb + pef)
     einc = edt / real(maxens2 + 1)
     do k = 1, maxens2
        edtc(k) = edt - real(k) * einc
        edtc(k) = -edtc(k) * pwav / pwev
        if(edtc(k) > edtmax) edtc(k) = edtmax
        if(edtc(k) < edtmin) edtc(k) = edtmin
     enddo
  endif

  return
end subroutine cup_dd_edt

!======================================================================

subroutine cup_dd_aa0(mza, ka, edt, ierr, aa0, jmin, gamma_cup, t_cup, &
     hcd, hes_cup, z, zd)

  use consts_coms, only: grav, cp

  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: ka
  real, intent(in) :: edt
  integer, intent(in) :: ierr
  real, intent(inout) :: aa0
  integer, intent(in) :: jmin
  real, dimension(mza), intent(in) :: gamma_cup
  real, dimension(mza), intent(in) :: t_cup
  real, dimension(mza), intent(in) :: hcd
  real, dimension(mza), intent(in) :: hes_cup
  real, dimension(mza), intent(in) :: z
  real, dimension(mza), intent(in) :: zd
  integer :: k
  integer :: kk
  real :: dz


  do k = 1, ka - 1
     if(ierr == 0 .and. k < jmin)then
        kk = jmin - k
        dz = z(kk) - z(kk+1)
        aa0 = aa0 + zd(kk) * edt * dz * (grav / (cp * t_cup(kk))) *   &
             ((hcd(kk) - hes_cup(kk)) / (1.0 + gamma_cup(kk)))
     endif
  enddo

  return
end subroutine cup_dd_aa0
